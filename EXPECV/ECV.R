## Empirical critical values, ECV.R #new
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
cpus = 20

# number of replications
iM = 10000

# number of nonlinear components
ir = 1
ikk = ir+2

vT = c(5,10,20)
vN = c(20,40,80,160)
vL = c(.99,.95,.90)

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

# parameters in the nonlinear part
B1 = c(0.7,0.7)
GA = 4; C1 = 3.5; C2 = c(3,4)
# different from 02_2
vc = C1

set.seed(seed)

## experiment

# (lm, wb, wcb), order, T, N, M
aLM1 = array(0, dim=c(3, 3,length(vT),length(vN), iM)) #hom
aLM2 = array(0, dim=c(3, 3,length(vT),length(vN), iM)) #het 

cat("Experiment ECV starts!\n")
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
        #tmp = fTF(vx=rep(mXX[itt,3],length(vc)),gamma=GA,vc=vc)
        #vY = c(vY, mu[inn] + (B01 + B1*tmp)%*%c(mXX[itt,1:2]) + rnorm(1) )
      }
      return(cbind(vY,mXX[2:(iT+1),]))
    }
    
    ftmp2 <- function(inn){
      mXX = matrix(0,iT+1,ikk); vY = NULL
      for(itt in 2:(iT+1)){
        mXX[itt,] = Theta%*%mXX[itt-1,] + kappa + c(rnorm(ikk) %*% chS)
        vY = c(vY, mu[inn] + B02[,inn]%*%c(mXX[itt,1:2]) + rnorm(1) )
        #tmp = fTF(vx=rep(mXX[itt,3],length(vc)),gamma=GA,vc=vc)
        #vY = c(vY, mu[inn] + (B02[,inn] + B1*tmp)%*%c(mXX[itt,1:2]) + rnorm(1) )
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
  sfExport(list=c("WCBLinTest","sLMTEST","iT","iN","msigma","B01","B1","ikk","Theta","kappa","chS"))
  
  RET = sfLapply(1:iM,ftmp)
  sfStop()
  
  for(mter in 1:iM) for(lter in 1:3) for(oter in 1:3){
    aLM1[lter,oter,tter,nter,mter] = RET[[mter]][[1]][oter,lter]
    aLM2[lter,oter,tter,nter,mter] = RET[[mter]][[2]][oter,lter]
  }
  cat(">")
  
}

save.image("RESULTS/ECV.RData") #new

proc.time() - ptm
cat("Done!\n")
