## Sensitivity of the LM test based on different parameters, SENS.R #new
## to sample from the null hypothesis#new
## author: Yukai Yang
## CORE, Université catholique de Louvain
## November 2014 #new

rm(list=ls(all=TRUE))
source("PSTR.R")
require(snowfall)

## initialization part

# the seed for random number generator
seed = 1

# number of cpu
cpus = 40

# number of replications
iM = 10000

# number of nonlinear components
ir = 1
ikk = ir+2

vT = c(5,20)
vN = c(20,80)

# 
msigma = 10

## parameters

kappa = c(0.2, 0.2, rep(2.45,ir))
Theta = diag(c(0.5, 0.4, rep(0.3,ir)), ir+2)
mR = diag(1, ir+2); mR[mR==0] = 1/3
mD = diag(sqrt(0.3), ir+2)
mSigma = mD%*%mR%*%t(mD)
chS = chol(mSigma)
# B01 = c(1,1) #new
B011 = seq(-1,1,.2); B012 = B011

set.seed(seed)

## experiment

# B011, B012, (lm, wb, wcb), order, T, N, M
aLM1 = array(0, dim=c(length(B011), length(B012), 3, 3,length(vT),length(vN), iM)) #hom
aLM2 = array(0, dim=c(length(B011), length(B012), 3, 3,length(vT),length(vN), iM)) #het 

cat("Experiment SENS starts!\n")
ptm <- proc.time()

for(tter in 1:length(vT)) for(nter in 1:length(vN)) for(b1 in 1:length(B011)) for(b2 in 1:length(B012)){
  
  iT = vT[tter]; iN = vN[nter]; B01 = c(B011[b1],B012[b2]) # new
  
  ftmp <- function(mter){
    ## some temporary functions
    
    ftmp1 <- function(inn){
      mXX = matrix(0,iT+1,ikk); vY = NULL
      for(itt in 2:(iT+1)){
        mXX[itt,] = Theta%*%mXX[itt-1,] + kappa + c(rnorm(ikk) %*% chS)
        vY = c(vY, mu[inn] + B01%*%c(mXX[itt,1:2]) + rnorm(1) )
      }
      return(cbind(vY,mXX[2:(iT+1),]))
    }
    
    ftmp2 <- function(inn){
      mXX = matrix(0,iT+1,ikk); vY = NULL
      for(itt in 2:(iT+1)){
        mXX[itt,] = Theta%*%mXX[itt-1,] + kappa + c(rnorm(ikk) %*% chS)
        vY = c(vY, mu[inn] + B02[,inn]%*%c(mXX[itt,1:2]) + rnorm(1) )
      }
      return(cbind(vY,mXX[2:(iT+1),]))
    }
    
    mu = rnorm(iN)*msigma
    B02 = matrix(rnorm(iN*2),2,iN) + B01
    
    ret = sapply(1:iN,ftmp1)
    mX1 = array(ret,c(iT,ikk+1,iN))
    
    ret = sapply(1:iN,ftmp2)
    mX2 = array(ret,c(iT,ikk+1,iN))
    
    RET1 = WCBLinTest(aDat=mX1,im=3)
    RET2 = WCBLinTest(aDat=mX2,im=3)
    
    return(list(RET1,RET2))
  }
  
  # parallel computation
  sfInit(parallel=T,cpus=cpus)
  sfExport(list=c("WCBLinTest","sLMTEST","iT","iN","msigma","B01","ikk","Theta","kappa","chS"))
  
  RET = sfLapply(1:iM,ftmp)
  sfStop()
  
  for(mter in 1:iM) for(lter in 1:3) for(oter in 1:3){
    aLM1[b1,b2,lter,oter,tter,nter,mter] = RET[[mter]][[1]][oter,lter]
    aLM2[b1,b2,lter,oter,tter,nter,mter] = RET[[mter]][[2]][oter,lter]
  }
  cat(">> T = ",tter," N = ",nter," B01 = (",b1,", ",b2,") done!",sep='')
  
}

save.image("RESULTS/SENS.RData") #new

proc.time() - ptm
cat("Done!\n")
