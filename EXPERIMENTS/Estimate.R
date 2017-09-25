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


iT = vT[3]; iN = vN[4]
vc = C2

mu = rnorm(iN)*msigma
B02 = matrix(rnorm(iN*2),2,iN) + B01

ret = sapply(1:iN,ftmp1)
mX1 = array(ret,c(iT,ikk+1,iN))

ret = sapply(1:iN,ftmp2)
mX2 = array(ret,c(iT,ikk+1,iN))

ptm <- proc.time()
EST1 = EstPSTR(aDat=mX1, im=2, par=c(log(GA),vc))
EST2 = EstPSTR(aDat=mX2, im=2, par=c(log(GA),vc))

c(GA,vc)
c(EST1$gamma, EST1$c)
c(EST2$gamma, EST2$c)

EvalTest(pstrest=EST1,type="time-varying",im=3)$TV
EvalTest(pstrest=EST2,type="time-varying",im=3)$TV
proc.time() - ptm