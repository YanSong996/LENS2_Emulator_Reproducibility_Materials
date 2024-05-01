################################################################################
# This file includes all steps to Generate Daily Emulations                    #
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
Tseq=c((5*365+1):(6*365),(25*365+1):(26*365),(45*365+1):(46*365),
       (65*365+1):(66*365),(85*365+1):(86*365))
YN=length(Tseq)
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
thetas=read.csv("LENS2_Data/thetas.csv")$x
phis=read.csv("LENS2_Data/phis.csv")$x
id.land=read.csv("LENS2_Data/landid.csv")$x
# Load rho
Res.hatrho=read.csv("Daily/Outputs/Res_Hatrho.csv")$x
# Load beta_0, beta_1, beta_2, {a_k,b_k},k=1,2,3
BetaAB=readMat("Daily/Outputs/BetaAB.mat")$BetaAB
# Load sigma
Sig=read.csv("Daily/Outputs/Daily_Sig.csv")$x
# Load v
v2hat=read.csv("Daily/Outputs/v2hat_do.csv")$x
# Load parameters in TGH Auto-regressive model
Tukeyres=as.matrix(read.csv("Daily/Outputs/Tukeyres.csv")[,-1])
rtc=read.csv("Daily/Outputs/rtc.csv")$x
TPhi.hat=read.csv("Daily/Outputs/Phihat_Tukey.csv")$x
idBera.new=which(Tukeyres[,1]!=0)
# Load U (Please uncompress it)
CU.axial=readMat("Daily/Outputs/U.mat")$U
# Load Q_l and Q_o
Ll=36
Lo=L=68
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
  M=4
  X=matrix(0,YN,3+2*M)
  X[,1]=rep(1,YN)
  for(i in 1:YN){
    ty=ceiling(Tseq[i]/365)
    X[i,2]=forcing$total[265+ty] # forcing$X0[265]=2014
    X[i,3]=(1-rho)*(t(rho^seq(ty+263,0,by=-1))%*%(forcing$total[1:(264+ty)]))
    for(j in 1:M){
      X[i,(2+2*j):(3+2*j)]=c(cos(2*pi*Tseq[i]*j/365),sin(2*pi*Tseq[i]*j/365))
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
RNM=array(0,c(R,L^2,YN+5))
for(r in 1:R){
  for(t in 1:(YN+5)){
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
  for(m in 1:5){
    CGen.Dat.SHT=matrix(0,L^2,365+1)
    CGen.Dat.SHT[,1]=LL%*%RNM[r,,365*(m-1)+1]
    for(t in 2:(365+1)){
      CGen.Dat.SHT[,t]=LL%*%RNM[r,,365*(m-1)+t]+TPhi.hat*CGen.Dat.SHT[,t-1]
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
    Gen.Dat=foreach(i=1:365,
                    .combine = cbind,
                    .packages = c("QZ","pracma")) %dopar% funcT(i)
    stopCluster(cl)
    Gen.Dat.Y[r,,(365*(m-1)+1):(365*m)]=Gen.Dat+Dat.hat[,(365*(m-1)+1):(365*m)]
  }
}




