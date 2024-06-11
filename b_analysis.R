library(tidyverse)

rm(list=ls())
#options(warn = 2)  # This will give error to warning, which can be used to stop using tryCatch()
options(warn = 0) # This brings back to warning as warning. 

load ("MPM_DATA.Rdata")

#########################################
### MPM for each species
dem_param <- list()

for (k in 1:length(LW$Species)){
  sci_name_fishbase <- LW$Species[k]
  
  a <- LW$a[k]
  b <- LW$b[k]
  alpha <- LW$alpha[k]
  beta <- LW$beta[k]
  
  ID <- which((sci_name_fishlife==LW$Species[k]))
  COV_PRED <- as.data.frame(PREDICTION[[ID]]$Cov_pred)
  MEAN_PRED <- as.data.frame(t(PREDICTION[[ID]]$Mean_pred))
  rownames(MEAN_PRED)<-NULL
  commonName <- PREDICTION[[ID]]$commonName
  
  Amax <- ceiling(as.numeric(exp(MEAN_PRED["tmax"])))  # Maximum Age
  if (Amax <2) stop(paste ("Maximum age < 2. ",k,". ",commonName))
  Amat <- min(ceiling(as.numeric(exp(MEAN_PRED["tm"]))),Amax)    # Age of Maturation
  if (Amax < Amat) stop(paste ("Age of maturity is smaller than the maximum age. ",k,". ",commonName))
  Linf <- as.numeric(exp(MEAN_PRED["Loo"]))   # L_inf in the VB curve
  Kapp <- as.numeric(exp(MEAN_PRED["K"]) )    # kappa in the VB curve
  Lmat <- as.numeric(exp(MEAN_PRED["Lm"]))     # Length at maturity
  if (Lmat > Linf) stop(paste("The asymptotic length is smaller than the length of maturity. ",k,". ",commonName))
 
  ## a V-B Parameter (estimated based on length and age of maturation from FishLife Package)
  A0 <- min(-0.1,Amat+log(1-Lmat/Linf)/Kapp ) # age at zero length (Thorson 2019 fixed it at -0.1)
  
  ## Age, Length, biomass relationship
  Age<-c(1:Amax)
  L <- Linf*(1-exp(-Kapp*(Age-A0)))  # Length at given age
  W <- a*(alpha+beta*L)^b
  
  ## Calculate the survival rate  
  mu <- function(Age){
    L <- Linf*(1-exp(-Kapp*(Age-A0)))  # Length at given age
    W <- a*(alpha+beta*L)^b  # Biomass at a given age - based on FishBase
    return(-M_u * W ^ B_L)   # Instantaneous mortality at given age - based on Lorenzen (1996)
  }
  
  SHC<-exp(sapply(1:(Amax),function(x){integrate(mu,1,x)$value[1]})) # Integrate instantaneous mortality to get the natural log of survivorship and then take exponential
  S <- (SHC/lag(SHC))[-1]# Calculate the survival rate from survivorship curve
  ###### Assumes recruitment = maturity
  
  ## Create Leslie Matrix
  M<-rep(0,Amax^2) # all elements of matrix
  ID2<- c(1:(Amax-1)) # Get an index for the location of age-specific survival rate
  M[ID2*(1+Amax)]<-S[1:(length(S))] 
  M[Amat:Amax]<-W[Amat:Amax] # Assign the weight for reproduction to be scaled later
  M<- matrix(M,nrow = Amax, ncol = Amax, byrow = TRUE) # Create a matrix
  
  ## Find the multiplier on fertility rate to make lambda =1
  maxint<-1
  M2<-M
  Ltest<-Re(eigen(M2)[[1]][1])
  
  while (Ltest<1.1){
    M2[1,]<-M[1,]*maxint
    Ltest<-Re(eigen(M2)[[1]][1])
    maxint<-maxint*2
  }
  
  ## The function to be minimized to find the coefficient
  FGR<-function (C,M,LAM){
    M[1,]<- C^2*M[1,]
    FGR  <- (Re(eigen(M)[[1]][1])-LAM)^2
  } 
  
  ## When lambda is 1
  RESULT1<- optimize(f = FGR,
                     interval =c(0,maxint),
                     M = M,
                     LAM = 1,
                     tol = 10^-12)
  C1<-(RESULT1$minimum)^2
  
  M[1,] <- C1*M[1,]
  
  ## DD Process
  MASPS <- as.numeric(exp(MEAN_PRED["ln_MASPS"])) # The maximum annual spawner biomass per spawner biomass in excess of replacement
  pi <- as.numeric(1+MASPS/(1-exp(mean(log(S[which(Age[-Amax]>=(Amat-1))])))) ) # The maximum life-time spawner biomass per spawner biomass (not in excess)
  SPR <- sum(ifelse(Amax==Amat,W[Amat],cumprod(c(1,S[Amat:(Amax-1)]))*W[Amat:Amax])) # Spawner biomass per age 1 individual (Density Independent)
  S2M <- if_else(Amat==1,1,prod(S[1:(Amat-1)]))  # Survivorship to maturity (recruitment)
  DD_a <- pi/SPR/S2M  # Maximum recruit per spawner biomass
  DD_a2 <- if_else(M[1,]==0,0, DD_a*W)
  
  MD <- M
  
  DDF<-function (DD_b,MD,W,NT){
    DD_b2 <- DD_b #sqrt(DD_b^2)
    n <- matrix(1,nrow=Amax,ncol=1)
    D <- 10
    tryCatch({ while (D>10^(-12)){
      TB<-t(W[Amat:Amax])%*%n[Amat:Amax]/1000  # Weight in Kg
      MD[1,] <- DD_a2/c(1+DD_b2*TB)
      n2<-MD%*%n  # Population is projected
      
      D <- sum((n2-n)^2)
     
      n<-n2
    }
    },error=function(e){browser()})
    NSum<-sum(n[Amat:Amax,1])
    return((NSum-NT)^2)
  }
  
  NT <- 10000 # Total number of mature fish
  # When lambda is 1
  RESULT2<- optimize(f = DDF,
                     interval =c(0,1),
                     MD = MD,
                     W = W,
                     NT = NT,
                     tol = 10^-12)
  
  dem_param[[k]] <- list(sci_name_fishbase=sci_name_fishbase,M=M,MD=MD,DD_a=DD_a2,DD_b=RESULT2$minimum/1000, NT=NT, COV_PRED=COV_PRED,MEAN_PRED=MEAN_PRED, a=a,b=b, alpha=alpha, beta=beta, WT=W, commonName=commonName) 
  
}

dem_results <- list()

#length(dem_param)
for (k in 1:length(dem_param)){
  list2env(dem_param[[k]],.GlobalEnv)
  lambda <- eigen(M)[[1]]
  w<-Re(eigen(M)$vectors[,1])  
  v<-Re(eigen(t(M))$vectors[,1])   
  rho<- Re(lambda[1])/abs(lambda[2]) # Damping ratio
  w<-w/sum(w)       # Normalize w so that it sums to 1
  v<-v/v[1]         # Normalize v so that the first component is 1
  S<-((v%*%t(w))/as.vector(t(v)%*%w)) # Sensitivity Matrix
  E<-M*S/Re(lambda[1])      # Elasticity Matrix
  S <- S*(M>0) # Eliminate the elements of Sensitivity matrix where corresponding population matrix has 0
  
  ## Fertility matrix
  F<-matrix(0,nrow=nrow(M),ncol=ncol(M))
  F[1,1]<-1
  
  ## Generation time
  G <-c(Re(lambda[1])*(t(as.vector(v))%*%as.vector(w))/(t(as.vector(v))%*% (F%*%M) %*%as.vector(w)))
  
  ## Resilience
  MD <- M
  
  Amax <- ceiling(as.numeric(exp(MEAN_PRED["tmax"]))) 
  Amat <- ceiling(as.numeric(exp(MEAN_PRED["tm"])))   # Age of Maturation
  n <- matrix(1,nrow=Amax,ncol=1)
  
  index <-1
  D <- 10
  while (D > 10^(-12)) {
    TB<-t(WT[Amat:Amax])%*%n[Amat:Amax]/1000  # Weight in Kg
    MD[1,] <- DD_a/c(1+DD_b*TB)
    n2<-MD%*%n  # Population is projected
    D <- sum((n2-n)^2)
    n<-n2
    index <- index +1
  }
  
  J <- MD
  J[1,] <- DD_a/c(1+DD_b*TB)-DD_a*DD_b*sum(WT[Amat:Amax])*sum(n[Amat:Amax])/(c(1+DD_b*TB))^2
  RESILIENCE <- -log(Re(eigen(M)[[1]][1]))
  
  dem_results[[k]] <- list(sci_name=sci_name_fishbase, commonName=commonName,lambda = lambda, w=w, v=v, E=E, S=S, G=G,
                           RESILIENCE = RESILIENCE, rho=rho, Amax = Amax, Amat=Amat, M=M,MD=MD,DD_a=DD_a,DD_b=DD_b, 
                           NT=NT, COV_PRED=COV_PRED,MEAN_PRED=MEAN_PRED, a=a,b=b, alpha=alpha, beta=beta, WT=WT)
}

save(dem_results, file="Demo_Results.Rdata")
