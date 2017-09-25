## Experiment04_1.R
## vc = C1
## Size of test of parameter constancy, table 4
## author: Yukai Yang
## CORE, Université catholique de Louvain
## April 2014

rm(list=ls(all=TRUE))
source("NewPSTR.R")

## initialization part

# number of nonlinear components
ir = 1
ikk = ir+2

# the seed for random number generator
seed = 1

vT = c(5,10,20)
vN = c(20,40,80,160)

# number of replications
iM = 10000

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

# different from 04_2
vc = C1

set.seed(seed)

## experiment

# sig_level, order, T, N, type(ST_X, ST_F, HAC_X, HAC_F)
aRet1 = array(0, dim=c(3,3,length(vT),length(vN),4))
aRet2 = array(0, dim=c(3,3,length(vT),length(vN),4))

iErrEst = 0; iErrEva = 0
iFinish = 0

cat("Experiment 04_1 starts!\n")
ptm <- proc.time()

for(tter in 1:length(vT)) for(nter in 1:length(vN)){

  iT = vT[tter]; iN = vN[nter]

  R1_PV1_M1_X=NULL; R1_PV1_M1_F=NULL; R1_PV2_M1_X=NULL; R1_PV2_M1_F=NULL
  R1_PV1_M2_X=NULL; R1_PV1_M2_F=NULL; R1_PV2_M2_X=NULL; R1_PV2_M2_F=NULL
  R1_PV1_M3_X=NULL; R1_PV1_M3_F=NULL; R1_PV2_M3_X=NULL; R1_PV2_M3_F=NULL

  R2_PV1_M1_X=NULL; R2_PV1_M1_F=NULL; R2_PV2_M1_X=NULL; R2_PV2_M1_F=NULL
  R2_PV1_M2_X=NULL; R2_PV1_M2_F=NULL; R2_PV2_M2_X=NULL; R2_PV2_M2_F=NULL
  R2_PV1_M3_X=NULL; R2_PV1_M3_F=NULL; R2_PV2_M3_X=NULL; R2_PV2_M3_F=NULL

  for(mter in 1:iM){

    while(T){
    mu = rnorm(iN)*msigma
    B02 = matrix(rnorm(iN*2),2,iN) + B01
    ret = sapply(1:iN,ftmp1)
    mX1 = array(ret,c(iT,ikk+1,iN))
    ret = sapply(1:iN,ftmp2)
    mX2 = array(ret,c(iT,ikk+1,iN))

    EST1 = try( EstPSTR(aDat=mX1,im=1,par=c(log(GA),vc)), T)
    if(class(EST1)=='try-error') {iErrEst=iErrEst+1; next}
    RET1 = try( EvalTest(pstrest=EST1,type="time-varying",im=3), T)
    if(class(RET1)=='try-error') {iErrEva=iErrEva+1; next}
    RET1 = RET1$TV

    R1_PV1_M1_X = c(R1_PV1_M1_X, RET1[[1]]$PV1_X); R1_PV1_M1_F = c(R1_PV1_M1_F, RET1[[1]]$PV1_F)
    R1_PV1_M2_X = c(R1_PV1_M2_X, RET1[[2]]$PV1_X); R1_PV1_M2_F = c(R1_PV1_M2_F, RET1[[2]]$PV1_F)
    R1_PV1_M3_X = c(R1_PV1_M3_X, RET1[[3]]$PV1_X); R1_PV1_M3_F = c(R1_PV1_M3_F, RET1[[3]]$PV1_F)
    R1_PV2_M1_X = c(R1_PV2_M1_X, RET1[[1]]$PV2_X); R1_PV2_M1_F = c(R1_PV2_M1_F, RET1[[1]]$PV2_F)
    R1_PV2_M2_X = c(R1_PV2_M2_X, RET1[[2]]$PV2_X); R1_PV2_M2_F = c(R1_PV2_M2_F, RET1[[2]]$PV2_F)
    R1_PV2_M3_X = c(R1_PV2_M3_X, RET1[[3]]$PV2_X); R1_PV2_M3_F = c(R1_PV2_M3_F, RET1[[3]]$PV2_F)

    EST2 = try( EstPSTR(aDat=mX2,im=1,par=c(log(GA),vc)), T)
    if(class(EST2)=='try-error') {iErrEst=iErrEst+1; next}
    RET2 = try( EvalTest(pstrest=EST2,type="time-varying",im=3), T)
    if(class(RET2)=='try-error') {iErrEva=iErrEva+1; next}
    RET2= RET2$TV

    R2_PV1_M1_X = c(R2_PV1_M1_X, RET2[[1]]$PV1_X); R2_PV1_M1_F = c(R2_PV1_M1_F, RET2[[1]]$PV1_F)
    R2_PV1_M2_X = c(R2_PV1_M2_X, RET2[[2]]$PV1_X); R2_PV1_M2_F = c(R2_PV1_M2_F, RET2[[2]]$PV1_F)
    R2_PV1_M3_X = c(R2_PV1_M3_X, RET2[[3]]$PV1_X); R2_PV1_M3_F = c(R2_PV1_M3_F, RET2[[3]]$PV1_F)
    R2_PV2_M1_X = c(R2_PV2_M1_X, RET2[[1]]$PV2_X); R2_PV2_M1_F = c(R2_PV2_M1_F, RET2[[1]]$PV2_F)
    R2_PV2_M2_X = c(R2_PV2_M2_X, RET2[[2]]$PV2_X); R2_PV2_M2_F = c(R2_PV2_M2_F, RET2[[2]]$PV2_F)
    R2_PV2_M3_X = c(R2_PV2_M3_X, RET2[[3]]$PV2_X); R2_PV2_M3_F = c(R2_PV2_M3_F, RET2[[3]]$PV2_F)

    break
    }
    iFinish = iFinish+1

  }

  aRet1[1,1,tter,nter,1] = sum(R1_PV1_M1_X < 0.01)/iM
  aRet1[1,1,tter,nter,2] = sum(R1_PV1_M1_F < 0.01)/iM
  aRet1[1,1,tter,nter,3] = sum(R1_PV2_M1_X < 0.01)/iM
  aRet1[1,1,tter,nter,4] = sum(R1_PV2_M1_F < 0.01)/iM
  aRet1[1,2,tter,nter,1] = sum(R1_PV1_M2_X < 0.01)/iM
  aRet1[1,2,tter,nter,2] = sum(R1_PV1_M2_F < 0.01)/iM
  aRet1[1,2,tter,nter,3] = sum(R1_PV2_M2_X < 0.01)/iM
  aRet1[1,2,tter,nter,4] = sum(R1_PV2_M2_F < 0.01)/iM
  aRet1[1,3,tter,nter,1] = sum(R1_PV1_M3_X < 0.01)/iM
  aRet1[1,3,tter,nter,2] = sum(R1_PV1_M3_F < 0.01)/iM
  aRet1[1,3,tter,nter,3] = sum(R1_PV2_M3_X < 0.01)/iM
  aRet1[1,3,tter,nter,4] = sum(R1_PV2_M3_F < 0.01)/iM

  aRet1[2,1,tter,nter,1] = sum(R1_PV1_M1_X < 0.05)/iM
  aRet1[2,1,tter,nter,2] = sum(R1_PV1_M1_F < 0.05)/iM
  aRet1[2,1,tter,nter,3] = sum(R1_PV2_M1_X < 0.05)/iM
  aRet1[2,1,tter,nter,4] = sum(R1_PV2_M1_F < 0.05)/iM
  aRet1[2,2,tter,nter,1] = sum(R1_PV1_M2_X < 0.05)/iM
  aRet1[2,2,tter,nter,2] = sum(R1_PV1_M2_F < 0.05)/iM
  aRet1[2,2,tter,nter,3] = sum(R1_PV2_M2_X < 0.05)/iM
  aRet1[2,2,tter,nter,4] = sum(R1_PV2_M2_F < 0.05)/iM
  aRet1[2,3,tter,nter,1] = sum(R1_PV1_M3_X < 0.05)/iM
  aRet1[2,3,tter,nter,2] = sum(R1_PV1_M3_F < 0.05)/iM
  aRet1[2,3,tter,nter,3] = sum(R1_PV2_M3_X < 0.05)/iM
  aRet1[2,3,tter,nter,4] = sum(R1_PV2_M3_F < 0.05)/iM

  aRet1[3,1,tter,nter,1] = sum(R1_PV1_M1_X < 0.1)/iM
  aRet1[3,1,tter,nter,2] = sum(R1_PV1_M1_F < 0.1)/iM
  aRet1[3,1,tter,nter,3] = sum(R1_PV2_M1_X < 0.1)/iM
  aRet1[3,1,tter,nter,4] = sum(R1_PV2_M1_F < 0.1)/iM
  aRet1[3,2,tter,nter,1] = sum(R1_PV1_M2_X < 0.1)/iM
  aRet1[3,2,tter,nter,2] = sum(R1_PV1_M2_F < 0.1)/iM
  aRet1[3,2,tter,nter,3] = sum(R1_PV2_M2_X < 0.1)/iM
  aRet1[3,2,tter,nter,4] = sum(R1_PV2_M2_F < 0.1)/iM
  aRet1[3,3,tter,nter,1] = sum(R1_PV1_M3_X < 0.1)/iM
  aRet1[3,3,tter,nter,2] = sum(R1_PV1_M3_F < 0.1)/iM
  aRet1[3,3,tter,nter,3] = sum(R1_PV2_M3_X < 0.1)/iM
  aRet1[3,3,tter,nter,4] = sum(R1_PV2_M3_F < 0.1)/iM


  aRet2[1,1,tter,nter,1] = sum(R2_PV1_M1_X < 0.01)/iM
  aRet2[1,1,tter,nter,2] = sum(R2_PV1_M1_F < 0.01)/iM
  aRet2[1,1,tter,nter,3] = sum(R2_PV2_M1_X < 0.01)/iM
  aRet2[1,1,tter,nter,4] = sum(R2_PV2_M1_F < 0.01)/iM
  aRet2[1,2,tter,nter,1] = sum(R2_PV1_M2_X < 0.01)/iM
  aRet2[1,2,tter,nter,2] = sum(R2_PV1_M2_F < 0.01)/iM
  aRet2[1,2,tter,nter,3] = sum(R2_PV2_M2_X < 0.01)/iM
  aRet2[1,2,tter,nter,4] = sum(R2_PV2_M2_F < 0.01)/iM
  aRet2[1,3,tter,nter,1] = sum(R2_PV1_M3_X < 0.01)/iM
  aRet2[1,3,tter,nter,2] = sum(R2_PV1_M3_F < 0.01)/iM
  aRet2[1,3,tter,nter,3] = sum(R2_PV2_M3_X < 0.01)/iM
  aRet2[1,3,tter,nter,4] = sum(R2_PV2_M3_F < 0.01)/iM

  aRet2[2,1,tter,nter,1] = sum(R2_PV1_M1_X < 0.05)/iM
  aRet2[2,1,tter,nter,2] = sum(R2_PV1_M1_F < 0.05)/iM
  aRet2[2,1,tter,nter,3] = sum(R2_PV2_M1_X < 0.05)/iM
  aRet2[2,1,tter,nter,4] = sum(R2_PV2_M1_F < 0.05)/iM
  aRet2[2,2,tter,nter,1] = sum(R2_PV1_M2_X < 0.05)/iM
  aRet2[2,2,tter,nter,2] = sum(R2_PV1_M2_F < 0.05)/iM
  aRet2[2,2,tter,nter,3] = sum(R2_PV2_M2_X < 0.05)/iM
  aRet2[2,2,tter,nter,4] = sum(R2_PV2_M2_F < 0.05)/iM
  aRet2[2,3,tter,nter,1] = sum(R2_PV1_M3_X < 0.05)/iM
  aRet2[2,3,tter,nter,2] = sum(R2_PV1_M3_F < 0.05)/iM
  aRet2[2,3,tter,nter,3] = sum(R2_PV2_M3_X < 0.05)/iM
  aRet2[2,3,tter,nter,4] = sum(R2_PV2_M3_F < 0.05)/iM

  aRet2[3,1,tter,nter,1] = sum(R2_PV1_M1_X < 0.1)/iM
  aRet2[3,1,tter,nter,2] = sum(R2_PV1_M1_F < 0.1)/iM
  aRet2[3,1,tter,nter,3] = sum(R2_PV2_M1_X < 0.1)/iM
  aRet2[3,1,tter,nter,4] = sum(R2_PV2_M1_F < 0.1)/iM
  aRet2[3,2,tter,nter,1] = sum(R2_PV1_M2_X < 0.1)/iM
  aRet2[3,2,tter,nter,2] = sum(R2_PV1_M2_F < 0.1)/iM
  aRet2[3,2,tter,nter,3] = sum(R2_PV2_M2_X < 0.1)/iM
  aRet2[3,2,tter,nter,4] = sum(R2_PV2_M2_F < 0.1)/iM
  aRet2[3,3,tter,nter,1] = sum(R2_PV1_M3_X < 0.1)/iM
  aRet2[3,3,tter,nter,2] = sum(R2_PV1_M3_F < 0.1)/iM
  aRet2[3,3,tter,nter,3] = sum(R2_PV2_M3_X < 0.1)/iM
  aRet2[3,3,tter,nter,4] = sum(R2_PV2_M3_F < 0.1)/iM

  cat(">")

}

save.image("RESULTS/EXPERIMENT04_1.RData")

proc.time() - ptm
cat("Done!\n")