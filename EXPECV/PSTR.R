## PSTR.R
## Utility functions for PSTR
## author: Yukai Yang
## CORE, Université catholique de Louvain
## November 2014

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

sLMTEST <- function(iT, iN, vU, mX, mW, mM, s2, mX2, invXX)
# to compute the LM tests and p-values
{
  df1 = ncol(mW)
  df2 = iT*iN - df1 - iN - ncol(mX)

  mW2 = mM %*% mW
  mXW2 = t(mX2) %*% mW2

  S1 = ( t(mW2)%*%mW2 - t(mXW2) %*% invXX %*% mXW2 ) * s2

  LM1_X = c(t(vU) %*% mW2 %*% chol2inv(chol(S1)) %*% t(mW2) %*% vU)

  return(LM1_X)
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

WCBLinTest <- function(aDat, ik=2, im=1, iq=NULL)
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
  beta = chol2inv(chol(t(mXb)%*%mXb)) %*% t(mXb) %*% vYb
  vU = matrix(vY-mX%*%beta, iT, iN)
  mu = apply(vU, 2, mean)
  vU = c(t(t(vU) - mu))
  s2 = sum((vU-mean(vU))**2)/(iT*iN) # sigma^2

  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - mD%*%t(mD)/iT
  mX2 = mM %*% mX
  invXX = chol2inv(chol(t(mX2)%*%mX2)) 
  eY = mD%*%mu + mX%*%beta

  # WB
  ve1 = sample(c(1,-1),iT*iN,replace=T)*vU
  my1 = matrix(eY + ve1, iT, iN)
  vyb1 = c(t(t(my1) - apply(my1, 2, mean)))
  tmpb1 = chol2inv(chol(t(mXb)%*%mXb)) %*% t(mXb) %*% vyb1
  vu1 = matrix(c(c(my1)-mX%*%tmpb1), iT, iN)
  vu1 = c(t(t(vu1)-apply(vu1, 2, mean)))
  ss1 = sum((vu1-mean(vu1))**2)/(iT*iN) # sigma^2 

  # WCB
  ve2 = c(t(matrix(sample(c(1,-1),iN,replace=T), iN, iT))) * vU
  my2 = matrix(eY + ve2, iT, iN)
  vyb2 = c(t(t(my2) - apply(my2, 2, mean))) 
  tmpb2 = chol2inv(chol(t(mXb)%*%mXb)) %*% t(mXb) %*% vyb2
  vu2 = matrix(c(c(my2)-mX%*%tmpb2), iT, iN)
  vu2 = c(t(t(vu2)-apply(vu2, 2, mean)))
  ss2 = sum((vu2-mean(vu2))**2)/(iT*iN) # sigma^2 

  RET = NULL; mW = NULL
  for(mter in 1:im){
    mW = cbind(mW, mX*(vQ**mter)) 
    LM = sLMTEST(iT=iT,iN=iN,vU=vU,mX=mX,mW=mW,mM=mM,s2=s2,mX2=mX2,invXX=invXX)

    # WB
    qLM1 = sLMTEST(iT=iT,iN=iN,vU=vu1,mX=mX,mW=mW,mM=mM,s2=ss1,mX2=mX2,invXX=invXX) 
    # WCB
    qLM2 = sLMTEST(iT=iT,iN=iN,vU=vu2,mX=mX,mW=mW,mM=mM,s2=ss2,mX2=mX2,invXX=invXX)

    RET = rbind(RET, c(LM,qLM1,qLM2))
  }

  return(RET) 
}

WCBLinTest2 <- function(aDat, ik=2, im=1, iq=NULL)
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
  beta = chol2inv(chol(t(mXb)%*%mXb)) %*% t(mXb) %*% vYb
  vU = matrix(vY-mX%*%beta, iT, iN)
  mu = apply(vU, 2, mean)
  vU = c(t(t(vU) - mu))
  s2 = sum((vU-mean(vU))**2)/(iT*iN) # sigma^2

  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - mD%*%t(mD)/iT
  mX2 = mM %*% mX
  invXX = chol2inv(chol(t(mX2)%*%mX2)) 
  eY = mD%*%mu + mX%*%beta

  RET = NULL; mW = NULL
  for(mter in 1:im){
    mW = cbind(mW, mX*(vQ**mter)) 
    LM = sLMTEST(iT=iT,iN=iN,vU=vU,mX=mX,mW=mW,mM=mM,s2=s2,mX2=mX2,invXX=invXX)

    # WB
    ve = sample(c(1,-1),iT*iN,replace=T)*vU
    my = matrix(eY + ve, iT, iN)
    vyb = c(t(t(my) - apply(my, 2, mean)))
    tmp = chol2inv(chol(t(mXb)%*%mXb)) %*% t(mXb) %*% vyb 
    vu = matrix(c(c(my)-mX%*%tmp), iT, iN)
    vu = c(t(t(vu)-apply(vu, 2, mean)))
    ss = sum((vu-mean(vu))**2)/(iT*iN) # sigma^2 
    qLM1 = sLMTEST(iT=iT,iN=iN,vU=vu,mX=mX,mW=mW,mM=mM,s2=ss,mX2=mX2,invXX=invXX) 

    # WCB
    ve = c(t(matrix(sample(c(1,-1),iN,replace=T), iN, iT))) * vU
    my = matrix(eY + ve, iT, iN)
    vyb = c(t(t(my) - apply(my, 2, mean)))
    tmp = chol2inv(chol(t(mXb)%*%mXb)) %*% t(mXb) %*% vyb 
    vu = matrix(c(c(my)-mX%*%tmp), iT, iN)
    vu = c(t(t(vu)-apply(vu, 2, mean)))
    ss = sum((vu-mean(vu))**2)/(iT*iN) # sigma^2 
    qLM2 = sLMTEST(iT=iT,iN=iN,vU=vu,mX=mX,mW=mW,mM=mM,s2=ss,mX2=mX2,invXX=invXX)

    RET = rbind(RET, c(LM,qLM1,qLM2))
  }

  return(RET) 
}

DerGFunc <- function(vg,vs,vp)
# compute the derivative of dg/dgamma, dg/dc
# input:
#	vg, vector of the transition functions
#	vs, vector of the transition variables
#	vp, est. of nonlinear parameters, vp[1] gamma, otherwise c
# output: matrix of derivatives, row=length
{
  gamma = vp[1]; cc = vp[2:length(vp)]
  tmp1 = vg * (1-vg)
  tmp2 = matrix(vs, length(vs), length(cc))
  tmp2 = t(tmp2) - cc

  ret = tmp1 * apply(tmp2,2,prod)

  ftmp <- function(iter){
    tmp3 = tmp2; tmp3[iter,] = 1
    return(- tmp1 * apply(tmp3,2,prod) * gamma)
  }

  ret = cbind(ret, sapply(1:length(cc),ftmp))
  return(ret)
}

EstPSTR <- function(aDat, ik=2, im=1, iq=NULL, par)
# to estimate the PSTR model
# input: aDat, an array of the data, T which N
# ik, NOT the lag, number of x
# im, the biggest number of switches
# par, initial values of delta and c, vec 1+im
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
  mQ = t(matrix(vQ,iT*iN,im))

  ftmp <- function(vx) return(vx - mean(vx))

  residuleSumSquare <- function(vp){
    # vp[1] = log(gamma) or delta
    vg = fTF(vx=mQ,gamma=exp(vp[1]),vc=vp[2:length(vp)]) 
    mXX = mX * vg
    aXX = array(c(mXX), dim=c(iT,iN,ik))
    mXXb = cbind(mXb, matrix(c(apply(aXX,c(2,3),ftmp)), iT*iN, ik))
    tmp = chol2inv(chol(t(mXXb)%*%mXXb)) %*% t(mXXb) %*% vYb 
    mXX = cbind(mX, mXX)
    vU = c(apply(matrix(c(vY-mXX%*%tmp),iT,iN),2,ftmp))
    return(sum(vU*vU)/iT/iN) 
  }

  opt = optim(par=par,fn=residuleSumSquare,method="L-BFGS-B",
    lower=par-2,upper=par+2)

  # return value
  ret = list(); ret$iT=iT; ret$iN=iN; ret$iq=iq
  ret$aDat = aDat
  ret$delta = opt$par[1]; ret$gamma = exp(ret$delta)
  ret$c = opt$par[2:length(opt$par)]

  ret$mX = mX

  vg = fTF(vx=mQ,gamma=ret$gamma,vc=ret$c) # g_it
  ret$vg = vg
  mXX = mX * vg # x_it * g_it

  aXX = array(c(mXX), dim=c(iT,iN,ik))
  mXXb = cbind(mXb, matrix(c(apply(aXX,c(2,3),ftmp)), iT*iN, ik)) # mean adjust
  tmp = chol2inv(chol(t(mXXb)%*%mXXb)) %*% t(mXXb) %*% vYb # beta
  ret$beta = c(tmp)

  mXX = cbind(mX, mXX) # (x_it, x_it*g_it)
  ret$mXX = mXX

  ret$vU = c(apply(matrix(c(vY-mXX%*%tmp),iT,iN),2,ftmp))
  ret$eY = vY - ret$vU
  ret$s2 = c(ret$vU %*% ret$vU) / (iT*iN)

  ret$mD = DerGFunc(vg=vg,vs=vQ,vp=c(ret$gamma,ret$c))

  return(ret)
}

EvalTest <- function(pstrest, type=c("time-varying","heterogeneity"), im=1, vq=NULL)
# to implement the evaluation tests
# input:
#   pstrest, a return value from EstPSTR
#   im, the biggest number of switches, or better say, the biggest order of the test
#   vq, a new transition variable in vector
{
  RET = list()

  mD = diag(1,pstrest$iN) %x% rep(1,pstrest$iT)
  mM = diag(1, pstrest$iN*pstrest$iT) - mD%*%t(mD)/pstrest$iT

  tmp = c(pstrest$mX %*% pstrest$beta[(ncol(pstrest$mX)+1):length(pstrest$beta)])
  tmp = pstrest$mD * tmp
  mV = cbind(pstrest$mXX, tmp)
  mV2 = mM %*% mV
  invVV = chol2inv(chol(t(mV2)%*%mV2))

  if(length(grep("time-varying",type))>0){
    RET$TV = list(); length(RET$TV) = im

    vt = 1:pstrest$iT/pstrest$iT

    mW = NULL
    for(mter in 1:im){
      mW = cbind(mW, pstrest$mXX*(vt**mter))
      RET$TV[[mter]] = LMTEST(iT=pstrest$iT,iN=pstrest$iN,vU=pstrest$vU,
        mX=mV,mW=mW,mM=mM,s2=pstrest$s2,mX2=mV2,invXX=invVV)
    }
  }

  if(length(grep("heterogeneity",type))>0){
    RET$HT = list(); length(RET$HT) = im

    mW = NULL
    for(mter in 1:im){
      mW = cbind(mW, pstrest$mXX*(vq**mter))
      RET$HT[[mter]] = LMTEST(iT=pstrest$iT,iN=pstrest$iN,vU=pstrest$vU,
        mX=mV,mW=mW,mM=mM,s2=pstrest$s2,mX2=mV2,invXX=invVV)
    }
  }

  return(RET)
}


WCBTVTest <- function(pstrest, im=1)
# to implement the evaluation tests
# input:
#   pstrest, a return value from EstPSTR
#   im, the biggest number of switches, or better say, the biggest order of the test
{
  iT = pstrest$iT; iN = pstrest$iN
  aDat = pstrest$aDat
  vU = pstrest$vU; eY = pstrest$eY

  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - mD%*%t(mD)/iT

  tmp = c(pstrest$mX %*% pstrest$beta[(ncol(pstrest$mX)+1):length(pstrest$beta)])
  tmp = pstrest$mD * tmp
  mV = cbind(pstrest$mXX, tmp)
  mV2 = mM %*% mV
  invVV = chol2inv(chol(t(mV2)%*%mV2))

  vt = 1:iT/iT

  # WB
  ve1 = sample(c(1,-1),iT*iN,replace=T)*vU
  my1 = matrix(eY + ve1, iT, iN) 
  for(nter in 1:iN) aDat[,1,nter] = my1[,nter]
  EST = try(EstPSTR(aDat=aDat,im=1,par=c(pstrest$delta,pstrest$c)),T)
  if(class(EST)=='try-error') return(EST) 
  vu1 = EST$vU; ss1 = EST$s2 # sigma^2
  tmp = c(EST$mX%*%EST$beta[(ncol(EST$mX)+1):length(EST$beta)])
  tmp = EST$mD * tmp
  mV11 = cbind(EST$mXX, tmp)
  mV12 = mM %*% mV11
  invVV1 = chol2inv(chol(t(mV12)%*%mV12))

  # WCB
  ve2 = c(t(matrix(sample(c(1,-1),iN,replace=T),iN,iT)))*vU
  my2 = matrix(eY + ve2, iT, iN) 
  for(nter in 1:iN) aDat[,1,nter] = my2[,nter]
  EST = try(EstPSTR(aDat=aDat,im=1,par=c(pstrest$delta,pstrest$c)),T)
  if(class(EST)=='try-error') return(EST) 
  vu2 = EST$vU; ss2 = EST$s2 # sigma^2
  tmp = c(EST$mX%*%EST$beta[(ncol(EST$mX)+1):length(EST$beta)])
  tmp = EST$mD * tmp
  mV21 = cbind(EST$mXX, tmp)
  mV22 = mM %*% mV21
  invVV2 = chol2inv(chol(t(mV22)%*%mV22))

  RET= NULL; mW = NULL
  for(mter in 1:im){
    mW = cbind(mW, pstrest$mXX*(vt**mter))
    LM = sLMTEST(iT=iT,iN=iN,vU=vU,
      mX=mV,mW=mW,mM=mM,s2=pstrest$s2,mX2=mV2,invXX=invVV)

    # WB
    qLM1 = sLMTEST(iT=iT,iN=iN,vU=vu1,mX=mV11,mW=mW,mM=mM,s2=ss1,mX2=mV12,invXX=invVV1)
    # WCB
    qLM2 = sLMTEST(iT=iT,iN=iN,vU=vu2,mX=mV21,mW=mW,mM=mM,s2=ss2,mX2=mV22,invXX=invVV2) 

    RET = rbind(RET, c(LM,qLM1,qLM2))
  }

  return(RET)
}

WCBHETest <- function(pstrest, im=1, vq)
# to implement the evaluation tests
# input:
#   pstrest, a return value from EstPSTR
#   im, the biggest number of switches, or better say, the biggest order of the test
#   vq, a new transition variable in vector
{
  iT = pstrest$iT; iN = pstrest$iN
  aDat = pstrest$aDat
  vU = pstrest$vU; eY = pstrest$eY

  mD = diag(1,iN) %x% rep(1,iT)
  mM = diag(1, iN*iT) - mD%*%t(mD)/iT

  tmp = c(pstrest$mX %*% pstrest$beta[(ncol(pstrest$mX)+1):length(pstrest$beta)])
  tmp = pstrest$mD * tmp
  mV = cbind(pstrest$mXX, tmp)
  mV2 = mM %*% mV
  invVV = chol2inv(chol(t(mV2)%*%mV2))

  # WB
  ve1 = sample(c(1,-1),iT*iN,replace=T)*vU
  my1 = matrix(eY + ve1, iT, iN) 
  for(nter in 1:iN) aDat[,1,nter] = my1[,nter]
  EST = try(EstPSTR(aDat=aDat,im=1,par=c(pstrest$delta,pstrest$c)),T)
  if(class(EST)=='try-error') return(EST) 
  vu1 = EST$vU; ss1 = EST$s2 # sigma^2
  tmp = c(EST$mX%*%EST$beta[(ncol(EST$mX)+1):length(EST$beta)])
  tmp = EST$mD * tmp
  mV11 = cbind(EST$mXX, tmp)
  mV12 = mM %*% mV11
  invVV1 = chol2inv(chol(t(mV12)%*%mV12))

  # WCB
  ve2 = c(t(matrix(sample(c(1,-1),iN,replace=T),iN,iT)))*vU
  my2 = matrix(eY + ve2, iT, iN) 
  for(nter in 1:iN) aDat[,1,nter] = my2[,nter]
  EST = try(EstPSTR(aDat=aDat,im=1,par=c(pstrest$delta,pstrest$c)),T)
  if(class(EST)=='try-error') return(EST) 
  vu2 = EST$vU; ss2 = EST$s2 # sigma^2
  tmp = c(EST$mX%*%EST$beta[(ncol(EST$mX)+1):length(EST$beta)])
  tmp = EST$mD * tmp
  mV21 = cbind(EST$mXX, tmp)
  mV22 = mM %*% mV21
  invVV2 = chol2inv(chol(t(mV22)%*%mV22))

  RET = NULL; mW = NULL
  for(mter in 1:im){
    mW = cbind(mW, pstrest$mXX*(vq**mter))
    LM = sLMTEST(iT=iT,iN=iN,vU=vU,
      mX=mV,mW=mW,mM=mM,s2=pstrest$s2,mX2=mV2,invXX=invVV)

    # WB
    qLM1 = sLMTEST(iT=iT,iN=iN,vU=vu1,mX=mV11,mW=mW,mM=mM,s2=ss1,mX2=mV12,invXX=invVV1)
    # WCB
    qLM2 = sLMTEST(iT=iT,iN=iN,vU=vu2,mX=mV21,mW=mW,mM=mM,s2=ss2,mX2=mV22,invXX=invVV2) 

    RET = rbind(RET, c(LM,qLM1,qLM2))
  }

  return(RET)
}

