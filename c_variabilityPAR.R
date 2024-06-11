## Three packages are used.

#install.packages("MASS")
#install.package("tidyverse")
#install.packages("parallel")

library(tidyverse)
library(MASS)  # mvrnorm
library(parallel)

rm(list=ls())

# Assume at least maximum age of 2 because it is an age structured model (annual species is excluded)
# If simulation suggests Lm > Loo (very rare occasions), Loo is assumed to be 10% larger than Lm

#options(warn = 2)  # This will give error to warning, which can be used to stop using tryCatch()
options(warn = 0) # This brings back to warning as warning. 

cl <- makeCluster(8) # Ran on 16 Core Machine, but using only 8 cores

load ("MPM_DATA.Rdata")

sci_name<-sapply(PREDICTION, function(x){x$SciName})

set.seed(1234)

#########################################

FINAL_RESULTS <- list()
MPM_Results <- list()

for (k in 1:length(LW$Species)){ # Run over species
  ## k is an index for fishbase data
  
  sci_name_fishbase <- LW$Species[k]
  
  ## Obtain an index for FishLife data
  ID <- which(as.vector(sci_name_fishlife==LW$Species[k]))
  
  ## Extract the means and covariance of FishLife data
  COV_PRED <- as.data.frame(PREDICTION[[ID]]$Cov_pred)
  MEAN_PRED <- as.data.frame(t(PREDICTION[[ID]]$Mean_pred))
  rownames(MEAN_PRED)<-NULL
  
  commonName <- PREDICTION[[ID]]$commonName
  
  ## Simulate the FishLife data and FishBase data
  NS <- 1000 #Number of simulations
  
  SIM_PRED <- as.data.frame(mvrnorm(n=NS,mu=as.vector(t(MEAN_PRED)),Sigma=COV_PRED)) %>%
    mutate(tmax=ceiling(exp(tmax)), tm=exp(tm), Loo=exp(Loo), Kappa=exp(K), Lm=exp(Lm),MASPS=exp(ln_MASPS)) %>% #FishLife
    mutate (a=LW$a[k], b=LW$b[k], alpha=LW$alpha[k], beta=LW$beta[k],B_L=B_L,M_u=M_u) %>%
    mutate(tmax=if_else(tmax<2,2,tmax),Loo=if_else(Loo<Lm,Lm*1.1,Loo)) %>% # Here the maximum age is set to be at least 2. Sometimes although very rare, Loo<Lm. If so Loo=1.1*Lm
    mutate(A0 = min(-0.1,tm+log(1-Lm/Loo)/Kappa))
  
  SIM_PRED$SIM_ID <- c(1:nrow(SIM_PRED))
  
  SIM_PRED_LIST <- split(SIM_PRED,as.vector(SIM_PRED$SIM_ID )) # Change to list
  
  demog_calc <- function(data) {
    with(as.list(data),{
      library(dplyr)
      Age<-c(1:tmax)
      L <- Loo*(1-exp(-Kappa*(Age-A0)))  # Length at given age
      WT <- a*(alpha+beta*L)^b
      
      ## Calculate the survival rate
      mu <- function(Age){
        L <- Loo*(1-exp(-Kappa*(Age-A0)))  # Length at given age
        WT <- a*(alpha+beta*L)^b  # Biomass at a given age - based on FishBase
        return(-M_u * WT ^ B_L)   # Instantaneous mortality at given age - based on Lorenzen (1996)
      }
      
      SHC<-exp(sapply(1:tmax,function(x){integrate(mu,1,x)$value[1]})) # Integrate instantaneous mortality to get the natural log of survivorship and then take exponential
      S <- (SHC/lag(SHC))[-1] # Calculate the survival rate from the survivorship curve
      
      ###### Assumes recruitment = maturity
      
      ## Create Leslie Matrix
      tmD<-min(ceiling(tm),tmax)
      M<-rep(0,tmax^2) # all elements of matrix
      ID<- c(1:(tmax-1)) # Get an index for the location of age-specific survival rate
      M[ID*(1+tmax)]<-S[1:(length(S))] # Assign the survival rate
      M[tmD:tmax]<-WT[tmD:tmax] # Assign the weight for reproduction to be scaled later
      M <- matrix(M,nrow = tmax, ncol = tmax, byrow = TRUE) # Create a matrix

      ## First find the multiplier that gives labmda >1 (it is to be used as end point in optimization)
      maxint<-1
      M2<-M
      Ltest <- Re(eigen(M2)[[1]][1])
     
      while (!Ltest>1.00){
        M2[1,]<-M[1,]*maxint
        Ltest<- Re(eigen(M2)[[1]][1])
        maxint<-maxint*2
      }
      
      ## The function to be minimized to find the multiplier
      FGR<-function (C,M,LAM){
        M[1,]<- C*M[1,]
        FGR  <- (Re(eigen(M)[[1]][1])-LAM)^2
      }
      
      ## Find the multiplier when lambda is 1
      RESULT1<- optimize(f = FGR,
                         interval =c(0,maxint),
                         M = M,
                         LAM = 1,
                         tol = 10^-12)
      C1<-(RESULT1$minimum)
      
      ## The final linear matrix
      M[1,] <- C1*M[1,]
      
      ### DD Process
      ## The maximum annual spawner biomass per spawner in excess of replacement
      pi <- as.numeric(1+MASPS/(1-exp(mean(log(S[which(Age[-tmax]>=(tmD-1))])))) ) #
      SPR <- sum(ifelse(tmax==tmD,WT[tmD],cumprod(c(1,S[tmD:(tmax-1)]))*WT[tmD:tmax])) # Spawner biomass per age 1 individual (Density Independent)
      S2M <- if_else(tmD==1,1,prod(S[1:tmD-1]))  # Survivorship to maturity (recruitment)
      DD_a <- pi/SPR/S2M  # Maximum recruit per spawner biomass
      DD_a2 <- if_else(M[1,]==0,0, DD_a*WT)
      
      MD <- M
      
      DDF<-function (DD_b,MD,WT,NT,tmax,tmD,DD_a2){
        DD_b2 <- DD_b 
        n <- matrix(1,nrow=tmax,ncol=1)
        D <- 10
        while (D > 10^(-12)){
          TB<-t(WT[tmD:tmax])%*%n[tmD:tmax]/1000  # Weight in Kg
          MD[1,] <- DD_a2/c(1+DD_b2*TB)
          n2<-MD%*%n  # Population is projected
          D <- sum((n2-n)^2)
          n<-n2
        }
        NSum<-sum(n[tmD:tmax,1])
        return((NSum-NT)^2)
      }
      
      NT <- 10000 # Total number of mature fish
      # When lambda is 1
      RESULT2<- optimize(f = DDF,
                         interval =c(0,1),
                         MD = MD,
                         WT = WT,
                         NT = NT,
                         tmax=tmax,
                         tmD=tmD,
                         DD_a2=DD_a2,
                         tol = 10^-12)
     
      DD_b <- RESULT2$minimum/1000
      
      ######Demographic Calculations
      
      lambda <- eigen(M)[[1]]
      w<-Re(eigen(M)$vectors[,1])
      v<-Re(eigen(t(M))$vectors[,1])
      rho<- Re(lambda[1])/abs(lambda[2]) # Damping ratio
      w<-w/sum(w)       # Normalize w so that it sums to 1
      v<-v/v[1]         # Normalize v so that the first component is 1
      S<-((v%*%t(w))/as.vector(t(v)%*%w)) # Sensitivity Matrix
      E<-M*S/Re(lambda[1])      # Elasticity Matrix
      S <- S*(M>0) # Eliminate the elements of Sensitivity matrix where corresponding population matrix has 0
      
      ### Generation time
      
      ## Fertility matrix
      F<-matrix(0,nrow=nrow(M),ncol=ncol(M))
      F[1,1]<-1
      
      G <-c(Re(lambda[1])*(t(as.vector(v))%*%as.vector(w))/(t(as.vector(v))%*% (F%*%M) %*%as.vector(w)))
      
      ### Resilience
      MD <- M

      n <- matrix(1,nrow=tmax,ncol=1)
      
      index <-1
      D <- 10
      while (D>10^(-12)){
        TB<-t(WT[tmD:tmax])%*%n[tmD:tmax]  # Weight in g
        MD[1,] <- DD_a2/c(1+DD_b*TB)
        n2 <- MD%*%n  # Population is projected
        D <- sum((n2-n)^2)
        n<-n2
        index <- index +1
      }
      #     #
      J <- MD
      J[1,] <- DD_a2/c(1+DD_b*TB)-DD_a2*DD_b*c(TB)/(c(1+DD_b*TB))^2
      EXP_RES <- -log(abs(eigen(J)[[1]][1]))
      results <- list(M=M,MD=MD,DD_a=DD_a,DD_b=DD_b,WT=WT,
                      EXP_RES = EXP_RES,G=G,S=S,E=E,v=v,w=w,rho=rho,SIM_ID=SIM_ID)
    })
  }
  
  results<-parLapply(cl=cl,SIM_PRED_LIST,demog_calc)
  
  TEMP_RESULTS <-as.data.frame(t(do.call(cbind, lapply(results,function(x)c(EXP_RES=x$EXP_RES,GEN=x$G,DPR=x$rho,SIM_ID=x$SIM_ID)))))
  MPM_Results[[k]] <- lapply(results,function(x)list(SIM_ID=x$SIM_ID,w=x$w,v=x$v,S=x$S,E=x$E,sci_name = sci_name_fishbase, commonName = commonName))
  
  
  FINAL_RESULTS[[k]] <- SIM_PRED %>%
    left_join(TEMP_RESULTS, by ="SIM_ID") %>%
    mutate(commonName=commonName,sci_name_fishbase=sci_name_fishbase)
}

stopCluster(cl)
save(FINAL_RESULTS,MPM_Results, file="Demo_Results_Variable_V2.Rdata")

