################################################################################
# This file includes all steps to Generate Annual Emulations                   #
################################################################################
library(R.matlab)
library(nloptr)
library(fpp2)
library(irlba)
library(fdaoutlier)
library(QZ)
library(pracma)
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(tseries)
library(approxOT)
source("Functions.R")
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



##### Inputs
###################################################################################
R=7
YN=86
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
id.land=read.csv("LENS2_Data/landid.csv")$x
thetas=read.csv("LENS2_Data/thetas.csv")$x
phis=read.csv("LENS2_Data/phis.csv")$x
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
forcing=forcing$total[266:351] 
# Load rho
Res.hatrho=read.csv("Annual/Outputs/Res_hatrho.csv")$x
# Load beta_0, beta_1, beta_2
Beta=as.matrix(read.csv("Annual/Outputs/Res_hatBeta.csv")[,-1])
# Load sigma
Sig=read.csv("Annual/Outputs/Sig_Annual.csv")$x
# Load v
v2hat=read.csv("Annual/Outputs/v2hat.csv")$x
# Load phi
TPhi.hat=read.csv("Annual/Outputs/Phihat.csv")$x
# Load U (Please uncompressed it first)
TU.axial=readMat("Annual/Outputs/U.mat")$U
# Load Q_l and Q_o
Ll=35
Lo=L=69
# Load P
P=1
# Give A
lseq=mseq=0
for(i in 1:143){
  lseq=c(lseq,rep(i,2*i+1))
  mseq=c(mseq,seq(-i,i,by=1))
}
A=matrix(0,L^2,L^2)
id1=which(mseq[1:L^2]==0)
id2=which(mseq[1:L^2]>0)
for(i in id1){
  A[i,i]=1
}
for(i in id2){
  l=lseq[i]
  m=mseq[i]
  A[c(l^2+l-m+1,i),c(l^2+l-m+1,i)]=matrix(c((-1)^(m+1)*1i,1i,(-1)^m,1),2,2)
}  



##### Stage 1: Preliminary
###################################################################################
getX=function(rho){
  X=matrix(1,YN,3)
  X[,2]=forcing
  X[1,3]=0
  for(i in 2:YN){
    X[i,3]=(1-rho)*(t(rho^seq(i-2,0,by=-1))%*%forcing[1:(i-1)])
  }
  return(X)
}
Dathat=function(i){
  X=getX(Res.hatrho[i])
  return(X%*%Beta[,i])
}
cl<- makeCluster(4) 
registerDoParallel(cl) 
Dat.hat=foreach(i=1:nrow(Dat.loc),
                   .combine=cbind,
                   .packages=c("nloptr")) %dopar% Dathat(i)
stopCluster(cl)
Dat.hat=t(Dat.hat)




##### Stage 2: Emulation
###################################################################################
LL.axial=t(chol(TU.axial))

set.seed(888)
RNM=array(0,c(R,L^2,YN+1))
for(r in 1:R){
  for(t in 1:(YN+1)){
    RNM[r,,t]=rnorm(L^2)
  }
}
set.seed(666)
EPS=array(0,c(R,nrow(Dat.loc),YN))
for(i in 1:nrow(Dat.loc)){
  EPS[,i,]=matrix(rnorm(R*YN,mean=0,sd=sqrt(v2hat[i])),R,YN)
}

LL=LL.axial
Gen.Dat.Y=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  TGen.Dat.SHT=matrix(0,L^2,YN+1)
  TGen.Dat.SHT[,1]=LL%*%RNM[r,,1]
  for(t in 2:(YN+1)){
    TGen.Dat.SHT[,t]=LL%*%RNM[r,,t]+TPhi.hat*TGen.Dat.SHT[,t-1]
  }
  TGen.Dat.SHT=TGen.Dat.SHT[,-1]
  Gen.Dat.SHT=A%*%TGen.Dat.SHT
  funcT=function(t){
    flm=Gen.Dat.SHT[,t]
    fs=rep(0,nrow(Dat.loc))
    fs[-id.land]=(Re(c(t(emsht_inverse(flm,thetas,phis,L)))[-id.land])+EPS[r,-id.land,t])*Sig[-id.land]
    fs[id.land]=(Re(c(t(emsht_inverse(flm,thetas,phis,Ll)))[id.land])+EPS[r,id.land,t])*Sig[id.land]
    return(fs)
  }
  cl=makeCluster(4)
  registerDoParallel(cl)
  Gen.Dat=foreach(i=1:YN,
                  .combine = cbind,
                  .packages = c("QZ","pracma")) %dopar% funcT(i)
  stopCluster(cl)
  Gen.Dat.Y[r,,]=Gen.Dat+Dat.hat
}








