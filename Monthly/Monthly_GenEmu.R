################################################################################
# This file includes all steps to Generate Monthly Emulations                  #
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
library(moments)
library(tseries)
library(approxOT)
source("Functions.R")
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



##### Inputs
###################################################################################
R=7
YN=86*12
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
id.land=read.csv("LENS2_Data/landid.csv")$x
thetas=read.csv("LENS2_Data/thetas.csv")$x
phis=read.csv("LENS2_Data/phis.csv")$x
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
# Load rho
Res.hatrho=read.csv("Monthly/Outputs/Res_Hatrho.csv")$x
# Load beta_0, beta_1, beta_2, {a_k,b_k},k=1,2,3
BetaAB=readMat("Monthly/Outputs/BetaAB.mat")$BetaAB
# Load sigma
Sig=read.csv("Monthly/Outputs/Monthly_Sig.csv")$x
# Load v
v2hat=read.csv("Monthly/Outputs/v2hat.csv")$x
# Load parameters in TGH Auto-regressive model
Tukeyres=as.matrix(read.csv("Monthly/Tukeyres.csv")[,-1])
rtc=read.csv("Monthly/Outputs/rtc.csv")$x
TPhi.hat=read.csv("Monthly/Outputs/Phihat_Tukey.csv")$x
idBera.new=which(Tukeyres[,1]!=0)
# Load U (Please uncompress it)
CU.axial=readMat("Monthly/Outputs/U.mat")$U    
# Load Q_l and Q_o
Ll=36
Lo=L=70
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
  M=3
  X=matrix(0,YN,3+2*M)
  X[,1]=rep(1,YN)
  for(i in 1:YN){
    ty=ceiling(i/12)
    X[i,2]=forcing$total[265+ty] # forcing$X0[265]=2014
    X[i,3]=(1-rho)*(t(rho^seq(ty+263,0,by=-1))%*%forcing$total[1:(264+ty)])
    for(j in 1:M){
      X[i,(2*j+2):(3+2*j)]=c(cos(2*pi*i*j/12),sin(2*pi*i*j/12))
    }
  }
  return(X)
}
Dathat=function(i){
  X=getX(Res.hatrho[i])
  return(X%*%BetaAB[,i])
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
ares=svd(CU.axial)
LL.axial=ares$u%*%diag(sqrt(ares$d))

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
  CGen.Dat.SHT=matrix(0,L^2,YN+1)
  CGen.Dat.SHT[,1]=LL%*%RNM[r,,1]
  for(t in 2:(YN+1)){
    CGen.Dat.SHT[,t]=LL%*%RNM[r,,t]+TPhi.hat*CGen.Dat.SHT[,t-1]
  }
  CGen.Dat.SHT=CGen.Dat.SHT[,-1]
  
  TGen.Dat.SHT=CGen.Dat.SHT
  for(j in idBera.new){
    taugh=function(z){return(Tukeyres[j,1]*(exp(Tukeyres[j,2]*z)-1)/Tukeyres[j,2]*exp(Tukeyres[j,3]*z^2/2))}
    TGen.Dat.SHT[j,]=taugh(CGen.Dat.SHT[j,]/rtc[j])
  }
  
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














