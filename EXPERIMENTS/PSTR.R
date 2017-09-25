## PSTR.R
## Utility functions for PSTR
## author: Yukai Yang
## CORE, Université catholique de Louvain
## April 2014

fTF <- function(vx, gamma, vc)
# to evaluate the transition function
# input: vx, gamma, vc
# the length of vc is the number of switches m or h
# vx can be a vector or a matrix
# if vx is a matrix, its ROW number must be equal to the length of vc!
# output: if vx is a vector, then a scalor retured, otherwise a vector.
{
  tmp = matrix(vx-vc, nrow=length(vc))
  tmp = -apply(tmp,2,prod)*gamma
  return(1 / ( exp(tmp) + 1 ))
}

LMTEST <- function(iT, iN, vU, mX, mW, mM, s2, mX2, invXX)
# to compute the LM tests and p-values
{
  df1 = ncol(mW)
  df2 = iT*iN - df1 - iN - ncol(mX)

  mW2 = mM %*% mW
  mXW2 = t(mX2) %*% mW2

  S1 = ( t(mW2)%*%mW2 - t(mXW2) %*% invXX %*% mXW2 ) * s2

  LM1_X = c(t(vU) %*% mW2 %*% chol2inv(chol(S1)) %*% t(mW2) %*% vU)
  PV1_X = 1-pchisq(LM1_X,df=df1)
  LM1_F = LM1_X * df2 / (iT*iN*df1)
  PV1_F = 1-pf(LM1_F,df1=df1,df2=df2)

  mZ = cbind(mX2, mW2)
  Delta = 0
  for(nter in  0:(iN-1)*iT){
    itmp = (nter+1):(nter+iT)
    tmp = t(mZ[itmp,]) %*% vU[itmp] %*% t(vU[itmp]) %*% mZ[itmp,]
    Delta = Delta + tmp
  }
  tmp = cbind( - t(mXW2) %*% invXX, diag(1, df1) )
  S2 = tmp %*% Delta %*% t(tmp)

  LM2_X = c(t(vU) %*% mW2 %*% chol2inv(chol(S2)) %*% t(mW2) %*% vU)
  PV2_X = 1-pchisq(LM2_X,df=df1)
  LM2_F = LM2_X * df2 / (iT*iN*df1)
  PV2_F = 1-pf(LM2_F,df1=df1,df2=df2)

  return(list(LM1_X=LM1_X, PV1_X=PV1_X, LM1_F=LM1_F, PV1_F=PV1_F,
    LM2_X=LM2_X, PV2_X=PV2_X, LM2_F=LM2_F, PV2_F=PV2_F))
}

LinTest <- function(aDat, ik=2, im=1, iq=NULL)
# linearity test
# input: aDat, an array of the data, T which N
# ik, NOT the lag, number of x
# im, the biggest number of switches, or better say, the biggest order of the test
{
  tmp = dim(aDat)
  iT = tmp[1]
  iN = tmp[3]
  if(is.null(iq)) iq = ik+2

  vY = NULL; vYb = NULL
  mX = NULL; mXb = NULL
  vQ = NULL
  for(nter in 1:iN){
    vY = c(vY, aDat[,1,nter])
    vYb = c(vYb, aDat[,1,nter] - mean(aDat[,1,nter]))
    tmp = aDat[,2:(ik+1),nter]
    mX = rbind(mX, tmp); mXb = rbind(mXb, t(t(tmp)-apply(t(tmp),1,mean)))
    vQ = c(vQ, aDat[,iq,nter])
  }
  tmp = chol2inv(chol(t(mXb)%*%mXb)) %*% t(mXb) %*% vYb
  vU = matrix(c(vY-mX%*%tmp), iT, iN)
  vU = c(t(t(vU)-apply(t(vU), 1, mean)))
  s2 = sum((vU-mean(vU))**2)/(iT*iN) # sigma^2

  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - mD%*%t(mD)/iT
  mX2 = mM %*% mX
  invXX = chol2inv(chol(t(mX2)%*%mX2))

  RET = list()
  length(RET) = im

  mW = NULL
  for(mter in 1:im){
    mW = cbind(mW, mX*(vQ**mter))
    RET[[mter]] = LMTEST(iT=iT,iN=iN,vU=vU,mX=mX,mW=mW,mM=mM,s2=s2,mX2=mX2,invXX=invXX)
  }

  return(RET)
  
}
