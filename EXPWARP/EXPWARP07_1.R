## EXPWARP07_1.R
## vc = C1, B2 = B1
## Power of test of no remaining heterogeneity, table 7
## author: Yukai Yang
## CORE, Université catholique de Louvain
## November 2014

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
B01 = c(1,1)

# parameters in the nonlinear part
B1 = c(.7,.7)
GA = 8; C1 = 3; C2 = c(3,4) 
# parameters in the nonlinear part 2
B2 = B1
#B2 = -B1
GB =8; CC = 4
vc = C1 

set.seed(seed)

## experiment

# wb, level, order, T, N
aRet1 = array(0, dim=c(2,3,3,length(vT),length(vN)))
aRet2 = array(0, dim=c(2,3,3,length(vT),length(vN)))

cat("Experiment EXPWARP07_1 starts!\n")
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
        tmpp = fTF(vx=rep(mXX[itt,3],length(vc)),gamma=GB,vc=CC)
        vY = c(vY, mu[inn] + (B01 + B1*tmp + B2*tmpp)%*%c(mXX[itt,1:2]) + rnorm(1) )
      }
      return(cbind(vY,mXX[2:(iT+1),]))
    }

    ftmp2 <- function(inn){
      mXX = matrix(0,iT+1,ikk); vY = NULL
      for(itt in 2:(iT+1)){
        mXX[itt,] = Theta%*%mXX[itt-1,] + kappa + c(rnorm(ikk) %*% chS)
        tmp = fTF(vx=rep(mXX[itt,3],length(vc)),gamma=GA,vc=vc)
        tmpp = fTF(vx=rep(mXX[itt,3],length(vc)),gamma=GB,vc=CC)
        vY = c(vY, mu[inn] + (B02[,inn] + B1*tmp + B2*tmpp)%*%c(mXX[itt,1:2]) + rnorm(1) )
      }
      return(cbind(vY,mXX[2:(iT+1),]))
    }

    while(T){
      mu = rnorm(iN)*msigma
      B02 = matrix(rnorm(iN*2),2,iN) + B01

      ret = sapply(1:iN,ftmp1)
      mX1 = array(ret,c(iT,ikk+1,iN))

      ret = sapply(1:iN,ftmp2)
      mX2 = array(ret,c(iT,ikk+1,iN))

      EST1 = try(EstPSTR(aDat=mX1,im=1,par=c(log(GA),vc)), T)
      if(class(EST1)=='try-error') next
      RET1 = try(WCBHETest(pstrest=EST1,im=3,vq=c(mX1[,ikk+1,])), T)
      if(class(RET1)=='try-error') next

      EST2 = try(EstPSTR(aDat=mX2,im=1,par=c(log(GA),vc)), T)
      if(class(EST2)=='try-error') next
      RET2 = try(WCBHETest(pstrest=EST2,im=3,vq=c(mX2[,ikk+1,])), T)
      if(class(RET2)=='try-error') next
      
      #EST1 = EstPSTR(aDat=mX1,im=1,par=c(log(GA),vc))
      #RET1 = WCBHETest(pstrest=EST1,im=3,vq=c(mX1[,ikk+1,]))
      #EST2 = EstPSTR(aDat=mX2,im=1,par=c(log(GA),vc))
      #RET2 = WCBHETest(pstrest=EST2,im=3,vq=c(mX2[,ikk+1,]))

      break
    }
    return(list(RET1,RET2))

  }

  # parallel computation
  sfInit(parallel=T,cpus=cpus)
  sfExport(list=c("fTF","WCBHETest","sLMTEST","EstPSTR","DerGFunc","iT","iN","GA","vc","msigma","B01","B1","ikk","Theta","kappa","chS","B2","GB","CC"))
  
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

save.image("RESULTS/EXPWARP07_1.RData")

proc.time() - ptm
cat("Done!\n")

