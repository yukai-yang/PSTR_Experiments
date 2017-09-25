## EXPWARP02_2.R #new
## Power of homogeneity test, table 2, wild bootstrap #new
## vc = C2
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
cpus = 10

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
vc = C2

set.seed(seed)

## experiment

# wb, level, order, T, N
aRet1 = array(0, dim=c(2, 3, 3,length(vT),length(vN))) #new
aRet2 = array(0, dim=c(2, 3, 3,length(vT),length(vN))) #new

cat("Experiment EXPWARP02_2 starts!\n")
ptm <- proc.time()

for(tter in 1:length(vT)) for(nter in 1:length(vN)){
  
  iT = vT[tter]; iN = vN[nter]
  
  ftmp <- function(mter){
    ## some temporary functions
    
    ftmp1 <- function(inn){
      mXX = matrix(0,iT+1,ikk); vY = NULL
      for(itt in 2:(iT+1)){
        mXX[itt,] = Theta%*%mXX[itt-1,] + kappa + c(rnorm(ikk) %*% chS)
        tmp = fTF(vx=rep(mXX[itt,3],length(vc)),gamma=GA,vc=vc)
        vY = c(vY, mu[inn] + (B01 + B1*tmp)%*%c(mXX[itt,1:2]) + rnorm(1) )
      }
      return(cbind(vY,mXX[2:(iT+1),]))
    }
    
    ftmp2 <- function(inn){
      mXX = matrix(0,iT+1,ikk); vY = NULL
      for(itt in 2:(iT+1)){
        mXX[itt,] = Theta%*%mXX[itt-1,] + kappa + c(rnorm(ikk) %*% chS)
        tmp = fTF(vx=rep(mXX[itt,3],length(vc)),gamma=GA,vc=vc)
        vY = c(vY, mu[inn] + (B02[,inn] + B1*tmp)%*%c(mXX[itt,1:2]) + rnorm(1) )
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
  sfExport(list=c("fTF","WCBLinTest","sLMTEST","iT","iN","GA","vc","msigma","B01","B1","ikk","Theta","kappa","chS"))
  
  RET = sfLapply(1:iM,ftmp)
  sfStop()
  
  mLM1 = NULL; qLM11 = NULL; qLM12 = NULL
  mLM2 = NULL; qLM21 = NULL; qLM22 = NULL
  for(mter in 1:iM){
    mLM1 = cbind(mLM1,RET[[mter]][[1]][,1])
    qLM11 = cbind(qLM11,RET[[mter]][[1]][,2])
    qLM12 = cbind(qLM12,RET[[mter]][[1]][,3])

    mLM2 = cbind(mLM2,RET[[mter]][[2]][,1])
    qLM21 = cbind(qLM21,RET[[mter]][[2]][,2])
    qLM22 = cbind(qLM22,RET[[mter]][[2]][,3])
  }
  
  for(oter in 1:3){
    tmp11 = quantile(qLM11[oter,],prob=vL) # hom wb
    tmp12 = quantile(qLM12[oter,],prob=vL) # hom wcb
    tmp21 = quantile(qLM21[oter,],prob=vL) # het wb
    tmp22 = quantile(qLM22[oter,],prob=vL) # het wcb
    for(lter in 1:3){
      # wb
      aRet1[1,lter,oter,tter,nter] = sum(mLM1[oter,]<=tmp11[lter])/iM
      aRet2[1,lter,oter,tter,nter] = sum(mLM2[oter,]<=tmp21[lter])/iM

      # wcb
      aRet1[2,lter,oter,tter,nter] = sum(mLM1[oter,]<=tmp12[lter])/iM
      aRet2[2,lter,oter,tter,nter] = sum(mLM2[oter,]<=tmp22[lter])/iM 
    }
  }
  cat(">")
  
}

save.image("RESULTS/EXPWARP02_2.RData") #new

proc.time() - ptm
cat("Done!\n")


