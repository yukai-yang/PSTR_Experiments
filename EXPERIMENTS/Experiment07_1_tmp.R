## Experiment07_1_tmp.R
## vc = C1, B2 = B1
## Power of test of no remaining heterogeneity, table 7
## author: Yukai Yang
## CORE, Université catholique de Louvain
## May 2014

rm(list=ls(all=TRUE))

load("TMP/EXPERIMENT07_1.RData")

## experiment

cat("Experiment 07_1_tmp starts!\n")
ptm <- proc.time()

for(tter in 3:length(vT)) for(nter in 4:length(vN)){

  #if(tter == 3 && nter < 3) next

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
    RET1 = try( EvalTest(pstrest=EST1,type="heterogeneity",im=3,vq=c(mX1[,ikk+1,])), T)
    if(class(RET1)=='try-error') {iErrEva=iErrEva+1; next}
    RET1 = RET1$HT#; print(RET1)

    R1_PV1_M1_X = c(R1_PV1_M1_X, RET1[[1]]$PV1_X); R1_PV1_M1_F = c(R1_PV1_M1_F, RET1[[1]]$PV1_F)
    R1_PV1_M2_X = c(R1_PV1_M2_X, RET1[[2]]$PV1_X); R1_PV1_M2_F = c(R1_PV1_M2_F, RET1[[2]]$PV1_F)
    R1_PV1_M3_X = c(R1_PV1_M3_X, RET1[[3]]$PV1_X); R1_PV1_M3_F = c(R1_PV1_M3_F, RET1[[3]]$PV1_F)
    R1_PV2_M1_X = c(R1_PV2_M1_X, RET1[[1]]$PV2_X); R1_PV2_M1_F = c(R1_PV2_M1_F, RET1[[1]]$PV2_F)
    R1_PV2_M2_X = c(R1_PV2_M2_X, RET1[[2]]$PV2_X); R1_PV2_M2_F = c(R1_PV2_M2_F, RET1[[2]]$PV2_F)
    R1_PV2_M3_X = c(R1_PV2_M3_X, RET1[[3]]$PV2_X); R1_PV2_M3_F = c(R1_PV2_M3_F, RET1[[3]]$PV2_F)

    EST2 = try( EstPSTR(aDat=mX2,im=1,par=c(log(GA),vc)), T)
    if(class(EST2)=='try-error') {iErrEst=iErrEst+1; next}
    RET2 = try( EvalTest(pstrest=EST2,type="heterogeneity",im=3,vq=c(mX2[,ikk+1,])), T)
    if(class(RET2)=='try-error') {iErrEva=iErrEva+1; next}
    RET2= RET2$HT#; print(RET2)

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

save.image("RESULTS/EXPERIMENT07_1.RData")

proc.time() - ptm
cat("Done!\n")