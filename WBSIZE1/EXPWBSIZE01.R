## Experiment01WB.R #new
## Size of homogeneity test, table 1, wild bootstrap #new
## author: Yukai Yang
## CORE, Université catholique de Louvain
## October 2014 #new

rm(list=ls(all=TRUE))
source("NewPSTR.R")
require(snowfall)

## initialization part

# the seed for random number generator
seed = 1

# number of cpu
cpus = 10

# number of replications
iM = 1000

# size of wb
iB = 500 #new

# number of nonlinear components
ir = 1
ikk = ir+2

vT = c(5,10,20)
vN = c(20,40,80,160)

# 
msigma = 10

## parameters

kappa = c(0.2, 0.2, rep(2.45,ir))
Theta = diag(c(0.5, 0.4, rep(0.3,ir)), ir+2)
mR = diag(1, ir+2); mR[mR==0] = 1/3
mD = diag(sqrt(0.3), ir+2)
mSigma = mD%*%mR%*%t(mD)
chS = chol(mSigma)
B01 = c(1,1) #new

set.seed(seed)

## experiment

# order, T, N
aRet1 = array(0, dim=c(3,length(vT),length(vN))) #new
aRet2 = array(0, dim=c(3,length(vT),length(vN))) #new

cat("Experiment 01WB starts!\n") #new
ptm <- proc.time()

for(tter in 1:length(vT)) for(nter in 1:length(vN)){
  
  iT = vT[tter]; iN = vN[nter]
  
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
        vY = c(vY, mu[inn] + c(B02[,inn])%*%c(mXX[itt,1:2]) + rnorm(1) )
      }
      return(cbind(vY,mXX[2:(iT+1),]))
    }
    
    mu = rnorm(iN)*msigma
    B02 = matrix(rnorm(iN*2),2,iN) + B01
    
    ret = sapply(1:iN,ftmp1)
    mX1 = array(ret,c(iT,ikk+1,iN))
    
    ret = sapply(1:iN,ftmp2)
    mX2 = array(ret,c(iT,ikk+1,iN))
    
    RET1 = WbLinTest(aDat=mX1,im=3,iB=iB)
    RET2 = WbLinTest(aDat=mX2,im=3,iB=iB)
    
    return(c(RET1,RET2))
  }
  
  # parallel computation
  sfInit(parallel=T,cpus=cpus)
  sfExport(list=c("WbLinTest","sLMTEST","iT","iN","msigma","B01","ikk","Theta","kappa","chS","iB"))
  
  RET = sfLapply(1:iM,ftmp)
  sfStop()
  
  mret = NULL
  for(mter in 1:iM) mret = rbind(mret,RET[[mter]])
  
  #new
  aRet1[1,tter,nter] = sum(mret[,1])/iM
  aRet1[2,tter,nter] = sum(mret[,2])/iM
  aRet1[3,tter,nter] = sum(mret[,3])/iM
  
  aRet2[1,tter,nter] = sum(mret[,4])/iM
  aRet2[2,tter,nter] = sum(mret[,5])/iM
  aRet2[3,tter,nter] = sum(mret[,6])/iM
  
  cat(">")
  
}

save.image("RESULTS/EXPWBSIZE01.RData") #new

proc.time() - ptm
cat("Done!\n")
