## EXPWARP05_2.R
## vc = C1 and CC a vector, multiple switches
## t, unnormalized
## Power of test of parameter constancy, table 5
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
cpus = 2

# number of replications
iM = 2

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
B1 = c(1,1)
GA = 3; C1 = 3.5; C2 = c(3,4) 
# parameters in the nonlinear part 2, time-varying part
B2 = B01*.7; B3 = B1*.7
GT =4; vc = C1 
# different from 05_1
CC = c(.3, .7) 

set.seed(seed)

## experiment

# wb, level, order, T, N
aRet1 = array(0, dim=c(2,3,3,length(vT),length(vN)))
aRet2 = array(0, dim=c(2,3,3,length(vT),length(vN)))

cat("Experiment EXPWARP05_2 starts!\n")
ptm <- proc.time()

for(tter in 1:length(vT)) for(nter in 1:length(vN)){

  iT = vT[tter]; iN = vN[nter]
  vt = 1:iT; ct = CC*iT

  ftmp <- function(mter){
    ## some temporary functions

    ftmp1 <- function(inn){
      mXX = matrix(0,iT+1,ikk); vY = NULL
      for(itt in 2:(iT+1)){
        mXX[itt,] = Theta%*%mXX[itt-1,] + kappa + c(rnorm(ikk) %*% chS)
        tmp = fTF(vx=rep(mXX[itt,3],length(vc)),gamma=GA,vc=vc)
        tmpp = fTF(vx=rep(vt[itt-1],length(ct)),gamma=GT,vc=ct)
        vY = c(vY, mu[inn] + (B01 + B1*tmp + B2*tmpp + B3*tmp*tmpp)%*%c(mXX[itt,1:2]) + rnorm(1) )
      }
      return(cbind(vY,mXX[2:(iT+1),]))
    }

    ftmp2 <- function(inn){
      mXX = matrix(0,iT+1,ikk); vY = NULL
      for(itt in 2:(iT+1)){
        mXX[itt,] = Theta%*%mXX[itt-1,] + kappa + c(rnorm(ikk) %*% chS)
        tmp = fTF(vx=rep(mXX[itt,3],length(vc)),gamma=GA,vc=vc)
        tmpp = fTF(vx=rep(vt[itt-1],length(ct)),gamma=GT,vc=ct)
        vY = c(vY, mu[inn] + (B02[,inn] + B1*tmp + B2*tmpp + B3*tmp*tmpp)%*%c(mXX[itt,1:2]) + rnorm(1) )
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
      RET1 = try(WCBTVTest(pstrest=EST1,im=3), T)
      if(class(RET1)=='try-error') next

      EST2 = try(EstPSTR(aDat=mX2,im=1,par=c(log(GA),vc)), T)
      if(class(EST2)=='try-error') next
      RET2 = try(WCBTVTest(pstrest=EST2,im=3), T)
      if(class(RET2)=='try-error') next

      break
    }
    return(list(RET1,RET2))

  }

  # parallel computation
  sfInit(parallel=T,cpus=cpus)
  sfExport(list=c("fTF","WCBTVTest","sLMTEST","EstPSTR","DerGFunc","iT","iN","GA","vc","msigma","B01","B1","ikk","Theta","kappa","chS","vt","ct","GT","B2","B3"))
  
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

save.image("RESULTS/EXPWARP05_2.RData")

proc.time() - ptm
cat("Done!\n")

