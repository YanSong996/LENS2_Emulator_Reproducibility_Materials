################################################################################
# This file includes all steps to develop an SG for Daily data                 #
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


##### Inputs: R monthly ensembles, KM=3, Q_l(Ll)=36, Q_o(Lo)=70, P=1
###################################################################################
### Load data information
R=7
Tseq=c((5*365+1):(6*365),(25*365+1):(26*365),(45*365+1):(46*365),
       (65*365+1):(66*365),(85*365+1):(86*365))
YN=length(Tseq)
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
thetas=read.csv("LENS2_Data/thetas.csv")$x
phis=read.csv("LENS2_Data/phis.csv")$x
id.land=read.csv("LENS2_Data/landid.csv")$x
Dat=array(0,c(R,nrow(Dat.loc),length(Tseq)))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Daily/dat_em1_day.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Daily/dat_em2_day.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Daily/dat_em3_day.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Daily/dat_em4_day.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Daily/dat_em5_day.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Daily/dat_em6_day.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Daily/dat_em7_day.mat")$Dat

### Load model information
Ll=36    #  Ll and Lo represent Q_l and Q_o in paper, respectively
L=Lo=68
P=1
lseq=mseq=0
for(i in 1:143){
  lseq=c(lseq,rep(i,2*i+1))
  mseq=c(mseq,seq(-i,i,by=1))
}



##### Stage 1: Deterministic component m_t and sigma
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
cl<- makeCluster(6)
registerDoParallel(cl)
Res.hatrho= foreach(i=1:nrow(Dat.loc),
                    .combine=cbind,
                    .packages=c("nloptr")) %dopar% hat.rho(i)
stopCluster(cl)
write.csv(Res.hatrho,"Daily/Outputs/Res_Hatrho.csv")


### Calculate and Store parameters beta_0, beta_1, beta_2, {a_k,b_k},k=1,2,3,4, and sigma 
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
cl<- makeCluster(6)
registerDoParallel(cl)
Res.BetaSig= foreach(i=1:nrow(Dat.loc),
                     .combine=cbind,
                     .packages=c("nloptr")) %dopar% BetaSighat(i)
stopCluster(cl)
writeMat("Daily/Outputs/BetaAB.mat",BetaAB=as.matrix(Res.BetaSig[1:11,]))
Sig=sqrt(c(Res.BetaSig[12,]))
write.csv(Sig,"Daily/Outputs/Daily_Sig.csv") 
Dat.hat=t(as.matrix(Res.BetaSig[-(1:12),]))


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
Dat.rsd.SHT=array(0,c(R,144^2,length(Tseq)))
for(r in 1:R){
  Dat.rsdd=Dat.srsd[r,,]
  
  cl<- makeCluster(15) 
  registerDoParallel(cl) 
  Dat.rsd.SHT[r,,]=foreach(t=1:length(Tseq),
                           .combine=cbind,
                           .packages=c("nloptr","QZ","pracma")) %dopar% emsht_forward(t(matrix(Dat.rsdd[,t],288,192)),thetas,phis,144)
  stopCluster(cl)
}

# Calculate and Store v
Dat.err=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  cl<- makeCluster(8) 
  registerDoParallel(cl) 
  Dat.rsd.hat=foreach(t=1:YN,
                      .combine=cbind,
                      .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:Ll^2,t],thetas,phis,Ll)))[id.land]
  stopCluster(cl)
  Dat.rsd.hat=Re(Dat.rsd.hat)
  Dat.err[r,id.land,]=(Dat.srsd[r,id.land,]-Dat.rsd.hat)^2
  
  cl<- makeCluster(4) 
  registerDoParallel(cl) 
  Dat.rsd.hat=foreach(t=1:YN,
                      .combine=cbind,
                      .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:L^2,t],thetas,phis,L)))[-id.land]
  stopCluster(cl)
  Dat.rsd.hat=Re(Dat.rsd.hat)
  Dat.err[r,-id.land,]=(Dat.srsd[r,-id.land,]-Dat.rsd.hat)^2
}
Dat.err=cbind(Dat.err[1,,],Dat.err[2,,],Dat.err[3,,],Dat.err[4,,],Dat.err[5,,],Dat.err[6,,],Dat.err[7,,])
v2hat=apply(Dat.err,1,mean)
write.csv(v2hat,"Daily/Outputs/v2hat_do.csv")


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
### TGH Auto-regressive model, Calculate and Store omega, g, h, lambda, and phi
### Test the normality of coefficients
cSke=cKur=cBera=rep(0,L^2)
for(i in 1:L^2){
  cSke[i]=skewness(c(TDat.rsd.SHT[,i,]))
  cKur[i]=kurtosis(c(TDat.rsd.SHT[,i,]))
  cBera[i]=jarque.bera.test(c(TDat.rsd.SHT[,i,]))$p.value
}
idBera=which(cBera<=0.05)  # non-normal
length(idBera)/L/L   # 0.5668253

### Calculate and Store parameters in TGH
p=P=1
Tukeyini=function(y){
  xx=c(y)
  p=seq(0.001,0.999,by=0.001)
  z=qnorm(p,median(xx),1)
  zz=quantile(xx,p)
  Tgh=function(para){
    b=para[1]
    g=para[2]
    h=para[3]
    res=b*(exp(g*z)-1)/g*exp(h*z^2/2)
    return(sum((zz-res)^2))
  }
  xxres=optim(par=c(0.5,-0.1,0.01),Tgh,lower=c(0.01,-0.5,0),method = "L-BFGS-B",upper=c(4,0.5,0.5))$par
  
  tildez=rep(0,length(xx))
  for(j in 1:length(xx)){
    taugh=function(x){
      if(abs(xxres[2])<1e-3){
        return(xxres[1]*x*exp(xxres[3]*x^2/2)-xx[j])
      }
      if(abs(xxres[2])>=1e-3){
        return(xxres[1]*(exp(xxres[2]*x)-1)*exp(xxres[3]*x^2/2)/xxres[2]-xx[j])
      }
    }
    tildez[j]=uniroot(taugh,c(-10,10))$root
  }
  return(c(xxres,jarque.bera.test(tildez)$p.value))
}

TukeyPara=function(para,y){
  b=para[1]
  g=para[2]
  h=para[3]
  phi=para[-(1:3)]
  taugh=function(z){
    if(g==0){
      return(b*z*exp(h*z^2/2))
    }
    else{return(b*(exp(g*z)-1)/g*exp(h*z^2/2))}
  }
  K=max(1000,ncol(y))
  Z=seq(-60,60,length.out=6*K)
  Y=taugh(Z)
  kid=matrix(.bincode(y,Y,right=FALSE),dim(y))
  tildez=Phi=matrix(0,nrow(y),ncol(y))
  u2=rep(0,nrow(y))
  for(r in 1:nrow(y)){
    for(i in 1:ncol(y)){
      tildez[r,i]=Z[kid[r,i]]+(y[r,i]-Y[kid[r,i]])/(Y[kid[r,i]+1]-Y[kid[r,i]])*(Z[2]-Z[1])
      Phi[r,i]=-h/2*tildez[r,i]^2-log(exp(g*tildez[r,i])+(exp(g*tildez[r,i])-1)/g*h*tildez[r,i])-log(b)
    }
    xx=rep(0,(365-p)*5)
    for(tt in p:1){
      xx=cbind(xx,tildez[r,c(tt:(tt+365-p-1),(tt+365):(tt+365*2-p-1),(tt+365*2):(tt+365*3-p-1),
                             (tt+365*3):(tt+365*4-p-1),(tt+365*4):(tt+365*5-p-1))])
    }
    xx=as.matrix(xx[,-1])
    u2[r]=mean((tildez[r,c((p+1):365,(p+366):(365*2),(p+365*2+1):(365*3),(p+365*3+1):(365*4),(p+365*4+1):(365*5))]-xx%*%phi)^2)
  }
  u2=mean(u2)
  res=sum(Phi)-nrow(y)*(365-p)*5/2*log(2*pi)-nrow(y)*(365-p)*5/2*log(u2)
  return(-res)
}

TukeyAuto=function(y){
  Beray=jarque.bera.test(c(y))$p.value
  skey=skewness(c(y))
  Tukeyinip=Tukeyini(y)
  
  if(skey<0){
    TukeyP=function(para){return(TukeyPara(para,-y))}
    res=optim(c(Tukeyinip[1],-Tukeyinip[2],max(Tukeyinip[3],0.01),rep(0.8,p)),TukeyP,lower=c(0.1,0.01,0.01,rep(0.1,p)),method="L-BFGS-B",upper=c(4,0.5,0.5,rep(1,p)))
    Paraest=res$par
    Paraest[2]=-Paraest[2]
  }
  else{
    TukeyP=function(para){return(TukeyPara(para,y))}
    res=optim(c(Tukeyinip[1],Tukeyinip[2],max(Tukeyinip[3],0.01),rep(0.8,p)),TukeyP,lower=c(0.1,0.01,0.01,rep(0.1,p)),method="L-BFGS-B",upper=c(4,0.5,0.5,rep(1,p)))
    Paraest=res$par
  }
  
  taugh=function(z){return(Paraest[1]*(exp(Paraest[2]*z)-1)/Paraest[2]*exp(Paraest[3]*z^2/2))}
  K=max(1000,ncol(y))
  Z=seq(-60,60,length.out=6*K)
  Y=taugh(Z)
  kid=matrix(.bincode(y,Y,right=FALSE),dim(y))
  tildez=matrix(0,nrow(y),ncol(y))
  for(r in 1:nrow(y)){
    for(i in 1:ncol(y)){
      tildez[r,i]=Z[kid[r,i]]+(y[r,i]-Y[kid[r,i]])/(Y[kid[r,i]+1]-Y[kid[r,i]])*(Z[2]-Z[1])
    }
  }
  Beraz=jarque.bera.test(c(tildez))$p.value
  
  Berares=c(Beray,Tukeyinip[4],Beraz)
  if(which.max(Berares)==1){
    return(c(0,0,0,rep(0,p)))
  }
  if(which.max(Berares)==2){
    return(c(Tukeyinip[1:3],rep(0,p)))
  }
  if(which.max(Berares)==3){
    return(Paraest)
  }
}

Tukeyres=matrix(0,L^2,4)
for(j in idBera){
  Tukeyres[j,]=TukeyAuto(TDat.rsd.SHT[,j,])
}
write.csv(Tukeyres,"Daily/Outputs/Tukeyres.csv")

idBera.new=which(Tukeyres[,1]!=0)    # need Tukey, i.e., set S_{gh}
idphi1=which(Tukeyres[,1]==0)     #    no need Tukey   
idphi2=intersect(which(Tukeyres[,1]!=0),which(Tukeyres[,4]==0))

### TGH transformation for idBera.new=S_{gh}
rtc=rep(1,L^2)
CDat.rsd.SHT=TDat.rsd.SHT
for(j in idBera.new){
  y=TDat.rsd.SHT[,j,]
  taugh=function(z){return(Tukeyres[j,1]*(exp(Tukeyres[j,2]*z)-1)/Tukeyres[j,2]*exp(Tukeyres[j,3]*z^2/2))}
  K=max(1000,ncol(y))
  Z=seq(-10,10,length.out=K)
  Y=taugh(Z)
  kid=matrix(.bincode(y,Y,right=FALSE),dim(y))
  tildez=matrix(0,nrow(y),ncol(y))
  for(r in 1:nrow(y)){
    for(i in 1:ncol(y)){
      tildez[r,i]=Z[kid[r,i]]+(y[r,i]-Y[kid[r,i]])/(Y[kid[r,i]+1]-Y[kid[r,i]])*(Z[2]-Z[1])
    }
  }
  rtc[j]=sd(y)/sd(tildez)
  CDat.rsd.SHT[,j,]=tildez*rtc[j]  
}
write.csv(rtc,"Daily/Outputs/rtc.csv")    # omega

### Calculate \phi for idphi1 and idphi2
TPhi.hat=rep(0,L^2)
for(i in c(idphi1,idphi2)){
  yy=c(t(CDat.rsd.SHT[,i,c((1+1):365,(1+366):(365*2),(1+365*2+1):(365*3),(1+365*3+1):(365*4),(1+365*4+1):(365*5))]))
  yx=matrix(0,(365-1)*5*R,1)
  for(r in 1:R){
    xx=rep(0,(365-1)*5)
    for(tt in 1:1){
      xx=cbind(xx,CDat.rsd.SHT[r,i,c(tt:(tt+365-1-1),(tt+365):(tt+365*2-1-1),(tt+365*2):(tt+365*3-1-1),
                                     (tt+365*3):(tt+365*4-1-1),(tt+365*4):(tt+365*5-1-1))])
    }
    yx[((365-1)*5*(r-1)+1):((365-1)*5*r),]=as.matrix(xx[,-1])
  }
  TPhi.hat[i]=solve(t(yx)%*%yx)%*%t(yx)%*%yy
}
TPhi.hat[-c(idphi1,idphi2)]=Tukeyres[-c(idphi1,idphi2),4]
write.csv(TPhi.hat,"Daily/Outputs/Phihat_Tukey.csv")


##### Step 4
### Calculate check K_0
CK.axial=matrix(0,L^2,L^2)
for(m in 0:(L-1)){
  idm=which(mseq[1:L^2]==m)
  idM=which(mseq[1:L^2]==-m)
  for(r in 1:R){
    for(t in 1:YN){
      CK.axial[idm,idm]=CK.axial[idm,idm]+crossprod(t(CDat.rsd.SHT[r,idm,t]))+
        crossprod(t(CDat.rsd.SHT[r,idM,t]))
    }
  }
  CK.axial[idm,idm]=CK.axial[idm,idm]/R/YN/2
  CK.axial[idM,idM]=CK.axial[idm,idm]
}


##### Step 5
### Calculate and Store U
CU.axial=CK.axial-outer(TPhi.hat,TPhi.hat,"*")*CK.axial
writeMat("Daily/Outputs/U.mat",U=CU.axial)  # Need to be compressed











