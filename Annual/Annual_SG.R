################################################################################
# This file includes all steps to develop an SG for Annual data                #
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


##### Inputs: R annual ensembles, K=0, Q_l(Ll)=35, Q_o(Lo)=69, P=1
###################################################################################
### Load data information
R=7
YN=86
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
forcing=forcing$total[266:351]   # Year 2015--2100
thetas=read.csv("LENS2_Data/thetas.csv")$x
phis=read.csv("LENS2_Data/phis.csv")$x
id.land=read.csv("LENS2_Data/landid.csv")$x
Dat=array(0,c(R,nrow(Dat.loc),YN))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Annual/dat_em1_year.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Annual/dat_em2_year.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Annual/dat_em3_year.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Annual/dat_em4_year.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Annual/dat_em5_year.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Annual/dat_em6_year.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Annual/dat_em7_year.mat")$Dat

### Load model information
Ll=35    #  Ll and Lo represent Q_l and Q_o in paper, respectively
L=Lo=69
P=1
lseq=mseq=0
for(i in 1:143){
  lseq=c(lseq,rep(i,2*i+1))
  mseq=c(mseq,seq(-i,i,by=1))
}




##### Stage 1: Deterministic component m_t and sigma
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
profile_negllh=function(y,rho){
  part1=t(y)%*%y
  X=getX(rho)
  part2=t(X)%*%y
  part2=t(part2)%*%solve(t(X)%*%X,part2)
  value=log(part1-part2)
  return(value)
}

##### Steps 1 and 2
### Calculate and Store parameters rho(L_i,l_j) 
hat.rho=function(i){
  obj=function(rho){
    value=0
    for(r in 1:R){value=value+profile_negllh(Dat[r,i,],rho)}
    return(value)
  }
  res=bobyqa(x0=c(0.9),fn=obj,lower=c(0.01),upper=c(0.99))
  return(res$par)
}
cl<- makeCluster(4) 
registerDoParallel(cl) 
Res.hatrho=foreach(i=1:nrow(Dat.loc),
                   .combine=cbind,
                   .packages=c("nloptr")) %dopar% hat.rho(i)
stopCluster(cl)
write.csv(c(Res.hatrho),"Annual/Outputs/Res_hatrho.csv") 

### Calculate and Store parameters beta_0, beta_1, beta_2, and sigma. 
### Calculate m_t at the same time
BetaSighat=function(i){
  X=getX(Res.hatrho[i])
  beta=solve(t(X)%*%X)%*%t(X)%*%t(Dat[,i,])
  Beta=apply(beta,1,mean)
  meanhat=X%*%Beta
  Sigmahat=Dat[,i,]-rep(1,R)%*%t(meanhat)
  sigmahat=mean(Sigmahat^2)
  return(c(Beta,sigmahat,meanhat))
}
cl<- makeCluster(4)
registerDoParallel(cl)
Res.BetaSig= foreach(i=1:nrow(Dat.loc),
                 .combine=cbind,
                 .packages=c("nloptr")) %dopar% BetaSighat(i)
stopCluster(cl)
write.csv(as.matrix(Res.BetaSig[1:3,]),"Annual/Outputs/Res_hatBeta.csv") 
Sig=sqrt(c(Res.BetaSig[4,]))
write.csv(Sig,"Annual/Outputs/Sig_Annual.csv") 
Dat.hat=t(as.matrix(Res.BetaSig[-(1:4),]))

##### Step 3
### Calculate stochastic component Z by detrending m_t and rescaling sigma
Dat.srsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.srsd[r,,]=(Dat[r,,]-Dat.hat)/Sig
}




##### Stage 2: Stochastic component Z_t 
###################################################################################
##### Step 1
### Do SHT 
Dat.rsd.SHT=array(0,c(R,144^2,YN))
for(r in 1:R){
  Dat.rsdd=Dat.srsd[r,,]
  
  cl<- makeCluster(4) 
  registerDoParallel(cl) 
  Dat.rsd.SHT[r,,]=foreach(t=1:YN,
                           .combine=cbind,
                           .packages=c("nloptr","QZ","pracma")) %dopar% emsht_forward(t(matrix(Dat.rsdd[,t],288,192)),thetas,phis,144)
  stopCluster(cl)
}

### Calculate and Store v
Dat.rsdd=Dat.rsd.hat=matrix(0,nrow(Dat.loc),YN*R)
for(r in 1:R){
  Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.srsd[r,,]
  
  cl<- makeCluster(4) 
  registerDoParallel(cl) 
  Dat.rsd.hat[,((r-1)*YN+1):(YN*r)]=foreach(t=1:YN,
                                            .combine=cbind,
                                            .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:Ll^2,t],thetas,phis,Ll)))
  stopCluster(cl)
}
v2hat=apply((Dat.rsdd-Re(Dat.rsd.hat))^2,1,mean)
Dat.rsdd=Dat.rsd.hat=matrix(0,nrow(Dat.loc),YN*R)
for(r in 1:R){
  Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.srsd[r,,]
  
  cl<- makeCluster(4) 
  registerDoParallel(cl) 
  Dat.rsd.hat[,((r-1)*YN+1):(YN*r)]=foreach(t=1:YN,
                                            .combine=cbind,
                                            .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:Lo^2,t],thetas,phis,Lo)))
  stopCluster(cl)
}
v2hat[-id.land]=apply((Dat.rsdd[-id.land,]-Re(Dat.rsd.hat[-id.land,]))^2,1,mean)
write.csv(v2hat,"Annual/Outputs/v2hat.csv")


##### Step 2
### Real-valued transformation
TDat.rsd.SHT=array(0,c(R,L^2,YN))
id1=which(mseq[1:L^2]==0)
id2=which(mseq[1:L^2]>0)
TDat.rsd.SHT[,id1,]=Re(Dat.rsd.SHT[,id1,])
TDat.rsd.SHT[,id2,]=Re(Dat.rsd.SHT[,id2,])
for(i in id2){
  l=lseq[i]
  m=mseq[i]
  TDat.rsd.SHT[,l^2+l-m+1,]=Im(Dat.rsd.SHT[,i,])
}  


##### Step 3
### Auto-regressive model, Calculate and Store phi
TPhi.hat=rep(0,L^2)
for(i in 1:L^2){
  yy=c(TDat.rsd.SHT[,i,2:YN])
  yx=c(TDat.rsd.SHT[,i,-YN])
  TPhi.hat[i]=t(yx)%*%yy/(t(yx)%*%yx)
}
write.csv(TPhi.hat,"Annual/Outputs/Phihat.csv")


##### Step 4
### Calculate tilde K_0
TK.axial=matrix(0,L^2,L^2)
for(m in 0:(L-1)){
  idm=which(mseq[1:L^2]==m)
  idM=which(mseq[1:L^2]==-m)
  for(r in 1:R){
    for(t in 1:YN){
      TK.axial[idm,idm]=TK.axial[idm,idm]+crossprod(t(TDat.rsd.SHT[r,idm,t]))+
        crossprod(t(TDat.rsd.SHT[r,idM,t]))
    }
  }
  TK.axial[idm,idm]=TK.axial[idm,idm]/R/YN/2
  TK.axial[idM,idM]=TK.axial[idm,idm]
}


##### Step 5
### Calculate and Store U
TU.axial=TK.axial-outer(TPhi.hat,TPhi.hat,"*")*TK.axial
writeMat("Annual/Outputs/U.mat",U=TU.axial)  # This file is compressed to reduce size






