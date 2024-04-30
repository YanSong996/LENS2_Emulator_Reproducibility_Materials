################################################################################################
# This file includes all implementations about Monthly Case Study (Fig.S9--Fig.S12)            #
################################################################################################
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


###############################################################################################
##### load the data
R=7
YN=86*12
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
thetas=read.csv("LENS2_Data/thetas.csv")$x
phis=read.csv("LENS2_Data/phis.csv")$x
id.land=read.csv("LENS2_Data/landid.csv")$x
Dat=array(0,c(R,nrow(Dat.loc),YN))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Monthly/dat_em1_month.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Monthly/dat_em2_month.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Monthly/dat_em3_month.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Monthly/dat_em4_month.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Monthly/dat_em5_month.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Monthly/dat_em6_month.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Monthly/dat_em7_month.mat")$Dat



############################################################################################
######################      Model the deterministic component       ######################
getX=function(rho,M){
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
profile_negllh=function(y,rho,M){
  part1=t(y)%*%y
  X=getX(rho,M)
  part2=t(X)%*%y
  part2=t(part2)%*%solve(t(X)%*%X,part2)
  value=log(part1-part2)
  return(value)
}
hat.rho=function(i){
  obj=function(rho){
    value=0
    for(r in 1:R){value=value+profile_negllh(Dat[r,i,],rho,M)}
    return(value)
  }
  res=bobyqa(x0=c(0.9),fn=obj,lower=c(0.01),upper=c(0.99))
  return(res$par)
}


################ Choose M
M=1  # M=1--5
getX=function(rho,M){
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
profile_negllh=function(y,rho,M){
  part1=t(y)%*%y
  X=getX(rho,M)
  part2=t(X)%*%y
  part2=t(part2)%*%solve(t(X)%*%X,part2)
  value=log(part1-part2)
  return(value)
}
hat.rho=function(i){
  obj=function(rho){
    value=0
    for(r in 1:R){value=value+profile_negllh(Dat[r,i,],rho,M)}
    return(value)
  }
  res=bobyqa(x0=c(0.9),fn=obj,lower=c(0.01),upper=c(0.99))
  return(res$par)
}
cl<- makeCluster(4) 
registerDoParallel(cl) 
Res.hatrho1= foreach(i=1:nrow(Dat.loc),
                     .combine=cbind,
                     .packages=c("nloptr")) %dopar% hat.rho(i)
stopCluster(cl)
# write.csv(c(Res.hatrho1),"Monthly/Res_hatrho_M1.csv")

Res.hatRho=matrix(0,5,nrow(Dat.loc))
Res.hatRho[1,]=read.csv("Monthly/Res_hatrho_M1.csv")$x
Res.hatRho[2,]=read.csv("Monthly/Res_hatrho_M2.csv")$x
Res.hatRho[3,]=read.csv("Monthly/Res_Hatrho.csv")$x
Res.hatRho[4,]=read.csv("Monthly/Res_hatrho_M4.csv")$x
Res.hatRho[5,]=read.csv("Monthly/Res_hatrho_M5.csv")$x
Bic=matrix(0,nrow(Dat.loc),5)
for(M in 1:5){
  Dathat=function(i){
    X=getX(Res.hatRho[M,i],M)
    beta=solve(t(X)%*%X)%*%t(X)%*%t(Dat[,i,])
    Beta=apply(beta,1,mean)
    meanhat=X%*%Beta
    Sigmahat=Dat[,i,]-rep(1,R)%*%t(meanhat)
    sigmahat=mean(Sigmahat^2)
    bicvalue=(2*M+3)*log(R*YN)+R*YN*log(sigmahat)
    return(bicvalue)
  }
  cl<- makeCluster(4) 
  registerDoParallel(cl) 
  Bic[,M]= foreach(i=1:nrow(Dat.loc),
                   .combine=cbind,
                   .packages=c("nloptr")) %dopar% Dathat(i)
  stopCluster(cl)
}
# write.csv(as.matrix(Bic),"Monthly/BIC_KM.csv")

a=apply(Bic,1,which.min)
length(which(a==1))/nrow(Dat.loc)  # 0.005009404
length(which(a==2))/nrow(Dat.loc)  # 0.1556894
length(which(a==3))/nrow(Dat.loc)  # 0.3378364
length(which(a==4))/nrow(Dat.loc)  # 0.340115
length(which(a==5))/nrow(Dat.loc)  # 0.1613498

### Plot Fig.S13(f)
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degree=as.factor(a))
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degree),size=0.8,data=dataF)+
  scale_color_manual(values = Col[c(1,2,3,5,6)])+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position ="right",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color=expression(K[M]))
PT

##### Choose M=3
M=3
cl<- makeCluster(4) 
registerDoParallel(cl) 
Res.hatrho= foreach(i=1:nrow(Dat.loc),
                    .combine=cbind,
                    .packages=c("nloptr")) %dopar% hat.rho(i)
stopCluster(cl)
write.csv(c(Res.hatrho),"Monthly/Outputs/Res_Hatrho.csv")


Dathat=function(i){
  X=getX(Res.hatrho[i],M=3)
  beta=solve(t(X)%*%X)%*%t(X)%*%t(Dat[,i,])
  Beta=apply(beta,1,mean)
  meanhat=X%*%Beta
  Sigmahat=Dat[,i,]-rep(1,R)%*%t(meanhat)
  sigmahat=mean(Sigmahat^2)
  return(c(meanhat,sigmahat))
}
cl<- makeCluster(4)
registerDoParallel(cl)
Dat.hat= foreach(i=1:nrow(Dat.loc),
                 .combine=cbind,
                 .packages=c()) %dopar% Dathat(i)
stopCluster(cl)
Sig=sqrt(Dat.hat[1033,])
Dat.hat=t(Dat.hat[1:1032,])
# writeMat("Monthly/Dathat.mat",Dathat=Dat.hat)
write.csv(Sig,"Monthly/Outputs/Monthly_Sig.csv")

### Plot Fig.S9(a)
Dat.mean=(Dat[1,,105]+Dat[2,,105]+Dat[3,,105]+Dat[4,,105]+Dat[5,,105]+Dat[6,,105]+Dat[7,,105])/7
Dat.hat=readMat("Monthly/Dathat.mat")$Dat.hat
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),z=abs(Dat.mean-Dat.hat[,105]))
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradient(low = "#3288BD",high ="#FEE08B")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color="\u00b0C")
PT

### Plot Fig.S9(b)
# Dat.hat=readMat("Monthly/Dathat.mat")$Dathat
# Sig=read.csv("Monthly/Monthly_Sig.csv")$x
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degree=Sig)
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degree),size=0.8,data=dataF)+
  scale_color_gradientn(values = seq(0,1,0.125),colours = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position ="right",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color="\u00b0C")
PT



############################################################################################
################            Model the stochastic component              ###############
Dat.srsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.srsd[r,,]=(Dat[r,,]-Dat.hat)/Sig
}

##### do SHT with Q=144
lseq=mseq=0
for(i in 1:143){
  lseq=c(lseq,rep(i,2*i+1))
  mseq=c(mseq,seq(-i,i,by=1))
}
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
# writeMat("Monthly/DatrsdshtDesig7.mat",Datrsdsht=Dat.rsd.SHT[7,,])
# Dat.rsd.SHT[7,,]=readMat("Monthly/DatrsdshtDesig7.mat")$Datrsdsht


##### Choose Q
llseq=c(20,30,40,50,60,70,80,90)
BIC.land=BIC.ocean=matrix(0,YN*R,length(llseq))
for(j in 1:length(llseq)){
  L=llseq[j]
  Dat.err=array(0,c(R,nrow(Dat.loc),YN))
  for(r in 1:R){
    cl<- makeCluster(5) 
    registerDoParallel(cl) 
    Dat.rsd.hat=foreach(t=1:YN,
                        .combine=cbind,
                        .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:L^2,t],thetas,phis,L)))
    stopCluster(cl)
    Dat.rsd.hat=Re(Dat.rsd.hat)
    Dat.err[r,,]=(Dat.srsd[r,,]-Dat.rsd.hat)^2
  }
  Dat.err=cbind(Dat.err[1,,],Dat.err[2,,],Dat.err[3,,],Dat.err[4,,],Dat.err[5,,],Dat.err[6,,],Dat.err[7,,])
  v2hat=apply(Dat.err,1,mean)
  BIC.land[,j]=apply((Dat.err/v2hat)[id.land,],2,sum)+sum(log(v2hat[id.land]))+log(2*pi)*length(id.land)+log(length(id.land))*(L^2)
  BIC.ocean[,j]=apply((Dat.err/v2hat)[-id.land,],2,sum)+sum(log(v2hat[-id.land]))+log(2*pi)*(nrow(Dat.loc)-length(id.land))+log(nrow(Dat.loc)-length(id.land))*(L^2)
}
write.csv(BIC.land,"Monthly/Outputs/BIC_land.csv")
write.csv(BIC.ocean,"Monthly/Outputs/BIC_ocean.csv")
llseq1=65:75
BICd.ocean=matrix(0,YN*R,length(llseq1))
for(j in 1:length(llseq1)){
  L=llseq1[j]
  Dat.err=array(0,c(R,nrow(Dat.loc),YN))
  for(r in 1:R){
    cl<- makeCluster(5) 
    registerDoParallel(cl) 
    Dat.rsd.hat=foreach(t=1:YN,
                        .combine=cbind,
                        .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:L^2,t],thetas,phis,L)))
    stopCluster(cl)
    Dat.rsd.hat=Re(Dat.rsd.hat)
    Dat.err[r,,]=(Dat.srsd[r,,]-Dat.rsd.hat)^2
  }
  Dat.err=cbind(Dat.err[1,,],Dat.err[2,,],Dat.err[3,,],Dat.err[4,,],Dat.err[5,,],Dat.err[6,,],Dat.err[7,,])
  v2hat=apply(Dat.err,1,mean)
  BICd.ocean[,j]=apply((Dat.err/v2hat)[-id.land,],2,sum)+sum(log(v2hat[-id.land]))+log(2*pi)*(nrow(Dat.loc)-length(id.land))+log(nrow(Dat.loc)-length(id.land))*(L^2)
}
write.csv(BICd.ocean,"Monthly/Outputs/BICd_ocean.csv")
llseq2=30:40
BICd.land=matrix(0,YN*R,length(llseq2))
for(j in 1:length(llseq2)){
  L=llseq2[j]
  Dat.err=array(0,c(R,nrow(Dat.loc),YN))
  for(r in 1:R){
    cl<- makeCluster(4) 
    registerDoParallel(cl) 
    Dat.rsd.hat=foreach(t=1:YN,
                        .combine=cbind,
                        .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:L^2,t],thetas,phis,L)))
    stopCluster(cl)
    Dat.rsd.hat=Re(Dat.rsd.hat)
    Dat.err[r,,]=(Dat.srsd[r,,]-Dat.rsd.hat)^2
  }
  Dat.err=cbind(Dat.err[1,,],Dat.err[2,,],Dat.err[3,,],Dat.err[4,,],Dat.err[5,,],Dat.err[6,,],Dat.err[7,,])
  v2hat=apply(Dat.err,1,mean)
  BICd.land[,j]=apply((Dat.err/v2hat)[id.land,],2,sum)+sum(log(v2hat[id.land]))+log(2*pi)*length(id.land)+log(length(id.land))*(L^2)
}
write.csv(BICd.land,"Monthly/Outputs/BICd_land.csv")

### Plot Fig.S10(a)
dataF=data.frame(bic=c(c(BIC.land),c(BIC.ocean)),
                 qseq=as.factor(c(rep(llseq,each=YN*R),rep(llseq,each=YN*R))),
                 group=as.factor(c(rep("Land",YN*R*length(llseq)),rep("Ocean",YN*R*length(llseq)))))
p=ggplot(data = dataF,aes(x=qseq,y=bic,color=group))+geom_boxplot()+
  scale_color_manual(values=c("#009E73", "#0072B2"))+
  scale_y_continuous(n.breaks = 3)+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(0,0),
                   legend.position =c(0,0),
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(0.8,"line"))+
  xlab(expression(Q))+ylab("BIC")
p
dataf1=data.frame(bic=apply(BICd.land,2,median),qseq=llseq2)
dataf2=data.frame(bic=apply(BICd.ocean,2,median),qseq=llseq1)
p+geom_point(mapping = aes(x=qseq/10-1.2,y=bic),data=dataf1,shape=4,col="black",size=0.7)+
  geom_point(mapping = aes(x=qseq/10-1.2,y=bic),data=dataf1[7,],shape=4,col="red",size=0.7)+
  geom_point(mapping = aes(x=qseq/10-0.8,y=bic),data=dataf2,shape=4,col="black",size=0.7)+
  geom_point(mapping = aes(x=qseq/10-0.8,y=bic),data=dataf2[6,],shape=4,col="red",size=0.7)



Ll=36
L=Lo=70
Dat.err=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  cl<- makeCluster(4) 
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
write.csv(v2hat,"Monthly/Outputs/v2hat.csv")

### Plot Fig.S10(b)
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z=sqrt(v2hat))
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),size=0.8,data=dataF)+
  scale_color_gradientn(values = seq(0,1,0.125),colours = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position ="right",
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  
  labs(color="\u00b0C")
PT

### Real-valued transformation
TDat.rsd.SHT=array(0,c(R,L^2,YN))
A=InvA=matrix(0,L^2,L^2)
id1=which(mseq[1:L^2]==0)
id2=which(mseq[1:L^2]>0)
for(i in id1){
  A[i,i]=1
  InvA[i,i]=1
}
TDat.rsd.SHT[,id1,]=Re(Dat.rsd.SHT[,id1,])
TDat.rsd.SHT[,id2,]=Re(Dat.rsd.SHT[,id2,])
for(i in id2){
  l=lseq[i]
  m=mseq[i]
  TDat.rsd.SHT[,l^2+l-m+1,]=Im(Dat.rsd.SHT[,i,])
  A[c(l^2+l-m+1,i),c(l^2+l-m+1,i)]=matrix(c((-1)^(m+1)*1i,1i,(-1)^m,1),2,2)
  InvA[c(l^2+l-m+1,i),c(l^2+l-m+1,i)]=matrix(c(-0.5*1i,-0.5,-0.5*1i,0.5),2,2)
}  

### Test the normality of coefficients
cSke=cKur=cBera=rep(0,L^2)
for(i in 1:L^2){
  cSke[i]=skewness(c(TDat.rsd.SHT[,i,]))
  cKur[i]=kurtosis(c(TDat.rsd.SHT[,i,]))
  cBera[i]=jarque.bera.test(c(TDat.rsd.SHT[,i,]))$p.value
}
idBera=which(cBera<=0.05)  # non-normal
length(idBera)/L/L   # 0.4263265



################################################################################
##########   Case 1: no Tukey g-and-h
### choose P
bicp=matrix(0,L^2,5)
for(i in 1:(L^2)){
  ts=TDat.rsd.SHT[,i,]
  for(p in 1:5){
    u2=rep(0,R)
    for(r in 1:R){
      xx=rep(0,YN-p)
      for(tt in p:1){
        xx=cbind(xx,ts[r,tt:(tt+YN-p-1)])
      }
      xx=as.matrix(xx[,-1])
      yy=as.matrix(ts[r,(p+1):YN])
      u2[r]=t(yy)%*%(diag(1,YN-p)-xx%*%solve(t(xx)%*%xx)%*%t(xx))%*%yy/(YN-p)
    }
    u2=mean(u2)
    bicp[i,p]=p*log((YN-p)*R)+R*(YN-p)*log(2*pi)+R*(YN-p)*log(u2)
  }
}
write.csv(bicp,"Monthly/Outputs/bicp_noTukey.csv")
length(which(apply(bicp,1,which.min)==1))/L/L   # 0.875102
length(which(apply(bicp,1,which.min)==2))/L/L   # 0.07632653
length(which(apply(bicp,1,which.min)==3))/L/L   # 0.03755102
length(which(apply(bicp,1,which.min)==4))/L/L   # 0.006734694
length(which(apply(bicp,1,which.min)==5))/L/L   # 0.004285714

### Choose P=1 
TPhi.hat=rep(0,L^2)
for(i in 1:L^2){
  yy=c(TDat.rsd.SHT[,i,2:YN])
  yx=c(TDat.rsd.SHT[,i,-YN])
  TPhi.hat[i]=t(yx)%*%yy/(t(yx)%*%yx)
}
write.csv(TPhi.hat,"Monthly/Outputs/Phihat_noTukey.csv")

### Model the spatial dependence K
# axial symmetric
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

### Generate emulations 
TU.axial=TK.axial-outer(TPhi.hat,TPhi.hat,"*")*TK.axial
#LL.axial=t(chol(TU.axial))
ares=svd(TU.axial)
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
# writeMat("Monthly/GenDatnoTukey7.mat",GenDat=Gen.Dat.Y[7,,])

### Calculate I.uq
central.region.area=function(y){
  yn=ncol(y)
  R=nrow(y)
  res=modified_band_depth(y)
  id.res=which(rank(-res,ties.method = "first")<=round(R/2))
  centralRegion=y[id.res,]
  return(sum(apply(centralRegion,2,max)-apply(centralRegion,2,min)))
}
get.IUQ=function(haty,y){
  return(central.region.area(haty)/central.region.area(y))
}
findIuq=function(i){
  return(get.IUQ(Gen.Dat.Y[,i,],Dat[1:R,i,]))
}
cl=makeCluster(4)
registerDoParallel(cl)
I.uq=foreach(i=1:nrow(Dat.loc),
             .combine = cbind,
             .packages = c("fdaoutlier")) %dopar% findIuq(i)
stopCluster(cl)
write.csv(c(I.uq),"Monthly/Outputs/Iuq_noTukey.csv")

### Calculate Wasserstein
WD.time=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  WD.time[i]=wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat.Y[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(WD.time,"Monthly/Outputs/WD_time_noTukey.csv")

WD.space=rep(0,YN)
for(t in 1:YN){
  WD.space[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat.Y[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
# write.csv(WD.space,"Monthly/WD_space_noTukey.csv")


################################################################################
##########   Case 2: Tukey g-and-h
### model the temporal dependence: choose P
bicpfunc1=function(j){
  bicp=rep(0,5)
  y=TDat.rsd.SHT[,j,]
  for(p in 1:5){
    skey=skewness(c(y))
    y1=y*sign(skey)
    TukeyPara=function(para){
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
      Z=seq(-10,10,length.out=K)
      Y=taugh(Z)
      kid=matrix(.bincode(y1,Y,right=FALSE),dim(y))
      tildez=Phi=matrix(0,nrow(y),ncol(y))
      u2=rep(0,nrow(y))
      for(r in 1:nrow(y)){
        for(i in 1:ncol(y)){
          tildez[r,i]=Z[kid[r,i]]+(y1[r,i]-Y[kid[r,i]])/(Y[kid[r,i]+1]-Y[kid[r,i]])*(Z[2]-Z[1])
          Phi[r,i]=-h/2*tildez[r,i]^2-log(exp(g*tildez[r,i])+(exp(g*tildez[r,i])-1)/g*h*tildez[r,i])-log(b)
        }
        xx=rep(0,ncol(y)-p)
        for(tt in p:1){
          xx=cbind(xx,tildez[r,tt:(tt+ncol(y)-p-1)])
        }
        xx=as.matrix(xx[,-1])
        u2[r]=mean((tildez[r,-(1:p)]-xx%*%phi)^2)
      }
      u2=mean(u2)
      res=sum(Phi)-nrow(y)*(ncol(y)-p)/2*log(2*pi)-nrow(y)*(ncol(y)-p)/2*log(u2)
      return(-res)
    }
    res=optim(c(0.5,0.2,0.1,rep(0.8,p)),TukeyPara,lower=c(0.1,0.05,0.05,rep(0.1,p)),method="L-BFGS-B",upper=c(4,0.5,0.5,rep(1,p)))
    Paraest=res$par
    Paraest[2]=sign(skey)*Paraest[2]
    
    taugh=function(z){return(Paraest[1]*(exp(Paraest[2]*z)-1)/Paraest[2]*exp(Paraest[3]*z^2/2))}
    
    K=max(1000,ncol(y))
    Z=seq(-10,10,length.out=K)
    Y=taugh(Z)
    kid=matrix(.bincode(y,Y,right=FALSE),dim(y))
    
    tildez=Phi=matrix(0,nrow(y),ncol(y))
    u2=rep(0,nrow(y))
    for(r in 1:nrow(y)){
      for(i in 1:ncol(y)){
        tildez[r,i]=Z[kid[r,i]]+(y[r,i]-Y[kid[r,i]])/(Y[kid[r,i]+1]-Y[kid[r,i]])*(Z[2]-Z[1])
        Phi[r,i]=-Paraest[3]/2*tildez[r,i]^2-log(exp(Paraest[2]*tildez[r,i])+
                                                   (exp(Paraest[2]*tildez[r,i])-1)/Paraest[2]*Paraest[3]*tildez[r,i])-log(Paraest[1])
      }
      xx=rep(0,ncol(y)-p)
      for(tt in p:1){
        xx=cbind(xx,tildez[r,tt:(tt+ncol(y)-p-1)])
      }
      xx=as.matrix(xx[,-1])
      u2[r]=mean((tildez[r,-(1:p)]-xx%*%Paraest[-(1:3)])^2)
    }
    u2=mean(u2)
    bicp[p]=p*log((ncol(y)-p)*nrow(y))+nrow(y)*(ncol(y)-p)*log(2*pi)+nrow(y)*(ncol(y)-p)*log(u2)
  }
  return(bicp)
}
cl<- makeCluster(4) 
registerDoParallel(cl) 
bicp1=foreach(j=idBera,
              .combine=rbind,
              .packages=c("moments","nloptr")) %dopar% bicpfunc1(j)
stopCluster(cl)
write.csv(bicp1,"Monthly/Outputs/bicp_Tukey.csv")

bicp=as.matrix(read.csv("Monthly/Outputs/bicp_noTukey.csv")[,-1])
bicp[idBera,]=as.matrix(bicp1)
length(which(apply(bicp,1,which.min)==1))/L/L   # 0.8728571
length(which(apply(bicp,1,which.min)==2))/L/L   # 0.08816327
length(which(apply(bicp,1,which.min)==3))/L/L   # 0.03408163
length(which(apply(bicp,1,which.min)==4))/L/L   # 0.03801722
length(which(apply(bicp,1,which.min)==5))/L/L   # 0.001836735


# Choose P=1 and calculate the parameters in Tukey g-and-h autoregressive model
p=1
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
    u2[r]=mean((tildez[r,-1]-phi*tildez[r,-ncol(y)])^2)
  }
  u2=mean(u2)
  res=sum(Phi)-nrow(y)*(ncol(y)-1)/2*log(2*pi)-nrow(y)*(ncol(y)-1)/2*log(u2)
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
write.csv(Tukeyres,"Monthly/Outputs/Tukeyres.csv")
idBera.new=which(Tukeyres[,1]!=0)    # need Tukey 
idphi1=which(Tukeyres[,1]==0)     #    no need Tukey   
idphi2=intersect(which(Tukeyres[,1]!=0),which(Tukeyres[,4]==0))


### Tukey g-and-h transformation for idBera.new
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
  CDat.rsd.SHT[,j,]=tildez*rtc[j]  # the sd of tildez will be the same as the sd of y
}
write.csv(rtc,"Monthly/Outputs/rtc.csv")


### Calculate \phi for idphi1 and idphi2
TPhi.hat=rep(0,L^2)
for(i in c(idphi1,idphi2)){
  yy=c(t(CDat.rsd.SHT[,i,2:YN]))
  yx=c(t(CDat.rsd.SHT[,i,-YN]))
  TPhi.hat[i]=solve(t(yx)%*%yx)%*%t(yx)%*%yy
}
TPhi.hat[-c(idphi1,idphi2)]=Tukeyres[-c(idphi1,idphi2),4]
write.csv(TPhi.hat,"Monthly/Outputs/Phihat_Tukey.csv")


### Plot Fig.S10(c)
dataF=data.frame(lon=c(lseq[1:L^2],0.1,0.1),lat=c(mseq[1:L^2],0.1,-0.1),z=c(TPhi.hat,0,1))
PT=ggplot()+xlab(expression(q))+ylab(expression(m))+
  geom_point(aes(x=lon,y=lat,colour=z),size=1.8,shape=15,data=dataF[1:L^2,])+
  geom_point(aes(x=lon,y=lat,colour=z),size=0.1,alpha=0.1,shape=15,data=dataF[-(1:L^2),])+
  #scale_color_gradientn(values = seq(0,1,0.125),colours = Col)+
  scale_color_gradient2(low="white",high="#0072B2")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.position ="right",
                   legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT


### Model the spatial dependence K
# axial symmetric
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

### Generate emulations 
CU.axial=CK.axial-outer(TPhi.hat,TPhi.hat,"*")*CK.axial
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
# writeMat("Monthly/GenDatTukey1.mat",GenDat=Gen.Dat.Y[1,,])
# Gen.Dat.Y[3,,]=readMat("Monthly/GenDatTukey3.mat")$GenDat

### Calculate I.uq
central.region.area=function(y){
  yn=ncol(y)
  R=nrow(y)
  res=modified_band_depth(y)
  id.res=which(rank(-res,ties.method = "first")<=round(R/2))
  centralRegion=y[id.res,]
  return(sum(apply(centralRegion,2,max)-apply(centralRegion,2,min)))
}
get.IUQ=function(haty,y){
  return(central.region.area(haty)/central.region.area(y))
}
findIuq=function(i){
  return(get.IUQ(Gen.Dat.Y[,i,],Dat[1:R,i,]))
}
cl=makeCluster(4)
registerDoParallel(cl)
I.uq.Tukey=foreach(i=1:nrow(Dat.loc),
             .combine = cbind,
             .packages = c("fdaoutlier")) %dopar% findIuq(i)
stopCluster(cl)
write.csv(c(I.uq.Tukey),"Monthly/Outputs/Iuq_Tukey.csv")

### Calculate Wasserstein
WD.time.Tukey=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  WD.time.Tukey[i]=wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat.Y[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(WD.time.Tukey,"Monthly/Outputs/WD_time_Tukey.csv")

WD.space.Tukey=rep(0,YN)
for(t in 1:YN){
  WD.space.Tukey[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat.Y[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(WD.space.Tukey,"Monthly/Outputs/WD_space_Tukey.csv")


### Plot Figs.S11(a) and S11(b)
I.uq=read.csv("Monthly/Outputs/Iuq_Tukey.csv")$x
I.uq.Huang=read.csv("Suppelment/Huang/Monthly/Outputs/Iuq_Huang.csv")$x
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z.H=I.uq.Huang,z=I.uq)
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradient2(low ="#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 1,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",limit=c(min(dataF$z.H),max(dataF$z)))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size = 12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color=expression(I[uq]))
PT
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z.H),data=dataF[which(dataF$z.H<=max(dataF$z)),],shape=15)+
  scale_colour_gradient2(low ="#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 1,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",limit=c(min(dataF$z.H),max(dataF$z)))+
  #geom_point(aes(x=lon,y=lat),data=dataF[c(7785,10377,18153,26793),],shape=4)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size = 12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color=expression(I[uq]))
PT

### Plot Figs.S11(c) and S11(d)
WD.time=read.csv("Monthly/Outputs/WD_time_Tukey.csv")$x
WD.Huang=read.csv("Supplement/Huang/Monthly/Outputs/WD_Huang.csv")$x
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z=c(WD.time),z.H=c(WD.Huang))
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradient(low = "#3288BD",high ="#FEE08B",limit=c(0,max(WD.time)))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color=expression(WD[S]))
PT
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z.H),data=dataF[which(dataF$z.H<=max(dataF$z)),],shape=15)+
  scale_colour_gradient(low = "#3288BD",high ="#FEE08B",limit=c(0,max(WD.time)))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color=expression(WD[S]))
PT

### Plot Figs.S11(e) and S11(f)
I.uq.noTukey=read.csv("Monthly/Outputs/Iuq_noTukey.csv")$x
WD.noTukey=read.csv("Monthly/Outputs/WD_time_noTukey.csv")$x
dataF=data.frame(Iuq=c(I.uq,I.uq.noTukey),
                 WD=c(WD.time,WD.noTukey),
                 Type=as.factor(rep(c("With TGH","Without TGH"),each=55296)))
PT=ggplot()+geom_boxplot(aes(x=Type,y=Iuq),data=dataF)+
  xlab(" ")+ylab(expression(I[uq]))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=14),
                   axis.title = element_text(size=12),
                   legend.justification = c(1,1),
                   legend.position =c(1,1),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  geom_hline(yintercept = 1,color="red")
PT
PT=ggplot()+geom_boxplot(aes(x=Type,y=WD),data=dataF)+
  xlab(" ")+ylab(expression(WD[S]))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=14),
                   axis.title = element_text(size=12),
                   legend.justification = c(1,1),
                   legend.position =c(1,1),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  geom_hline(yintercept = 0,color="red")
PT

### Plot Figs.S12(a) and S12(b)
DatP=Dat[,18281,]
DatQ=Gen.Dat.Y[,18281,]
dataF=data.frame(x=c(c(DatQ),c(DatP)),
                 group=as.factor(c(rep("Emulations",YN*R),rep("Simulations",YN*R))))
DatP=Dat[,,555]
DatQ=Gen.Dat.Y[,,555]
dataF=data.frame(x=c(c(DatQ),c(DatP)),
                 group=as.factor(c(rep("Emulations",nrow(Dat.loc)*R),rep("Simulations",nrow(Dat.loc)*R))))
PT=ggplot(dataF,aes(x,fill=group,color=group))+xlab(" ")+ylab("Density")+
  geom_histogram(aes(y=after_stat(density)),position = "identity",binwidth = 1,alpha=.5)+
  scale_fill_manual(values = c("#3288BD","#FEE08B"))+
  scale_color_manual(values = c("#3288BD","#FEE08B"))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.justification = c(0,1),
                   #legend.position =c(0.01,0.99),
                   legend.position = "none",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))
PT

### Plot Figs.S12(c) and S12(d)
Mbe=((2021-2015)*12+1):((2025-2014)*12)
dataF1=data.frame(Time=rep(Mbe,times=2),
                  TempO=c(Dat[1,18281,Mbe]-5,Gen.Dat.Y[1,18281,Mbe]),
                  TempL=c(Dat[1,39780,Mbe]-15,Gen.Dat.Y[1,39780,Mbe]),
                  Type=as.factor(rep(c("Simulations-5\u00b0C","Emulations"),each=60)))
dataF2=data.frame(Time=rep(Mbe,times=2),
                  TempO=c(Dat[2,18281,Mbe]-5,Gen.Dat.Y[2,18281,Mbe]),
                  TempL=c(Dat[2,39780,Mbe]-15,Gen.Dat.Y[2,39780,Mbe]),
                  Type=as.factor(rep(c("Simulations-5\u00b0C","Emulations"),each=60)))
dataF3=data.frame(Time=rep(Mbe,times=2),
                  TempO=c(Dat[3,18281,Mbe]-5,Gen.Dat.Y[3,18281,Mbe]),
                  TempL=c(Dat[3,39780,Mbe]-15,Gen.Dat.Y[3,39780,Mbe]),
                  Type=as.factor(rep(c("Simulations-5\u00b0C","Emulations"),each=60)))
dataF4=data.frame(Time=rep(Mbe,times=2),
                  TempO=c(Dat[4,18281,Mbe]-5,Gen.Dat.Y[4,18281,Mbe]),
                  TempL=c(Dat[4,39780,Mbe]-15,Gen.Dat.Y[4,39780,Mbe]),
                  Type=as.factor(rep(c("Simulations-5\u00b0C","Emulations"),each=60)))
dataF5=data.frame(Time=rep(Mbe,times=2),
                  TempO=c(Dat[5,18281,Mbe]-5,Gen.Dat.Y[5,18281,Mbe]),
                  TempL=c(Dat[5,39780,Mbe]-15,Gen.Dat.Y[5,39780,Mbe]),
                  Type=as.factor(rep(c("Simulations-5\u00b0C","Emulations"),each=60)))
dataF6=data.frame(Time=rep(Mbe,times=2),
                  TempO=c(Dat[6,18281,Mbe]-5,Gen.Dat.Y[6,18281,Mbe]),
                  TempL=c(Dat[6,39780,Mbe]-15,Gen.Dat.Y[6,39780,Mbe]),
                  Type=as.factor(rep(c("Simulations-5\u00b0C","Emulations"),each=60)))
dataF7=data.frame(Time=rep(Mbe,times=2),
                  TempO=c(Dat[7,18281,Mbe]-5,Gen.Dat.Y[7,18281,Mbe]),
                  TempL=c(Dat[7,39780,Mbe]-15,Gen.Dat.Y[7,39780,Mbe]),
                  Type=as.factor(rep(c("Simulations-5\u00b0C","Emulations"),each=60)))
dataFM=data.frame(Time=rep(Mbe,times=2),
                  TempO=c(apply(Dat[,18281,Mbe],2,mean)-5,apply(Gen.Dat.Y[,18281,Mbe],2,mean)),
                  TempL=c(apply(Dat[,39780,Mbe],2,mean)-15,apply(Gen.Dat.Y[,39780,Mbe],2,mean)),
                  Type=as.factor(rep(c("Simulations-5\u00b0C","Emulations"),each=60)))
dataL=data.frame(Time=rep(Mbe,times=2),
                 TempO=c(Dat.hat[18281,Mbe]-5,Dat.hat[18281,Mbe]),
                 TempL=c(Dat.hat[39780,Mbe]-15,Dat.hat[39780,Mbe]),
                 Type=as.factor(rep(c("Simulations-5\u00b0C","Emulations"),each=60)))
PT=ggplot()+
  geom_line(aes(x=Time,y=TempO,colour=Type),alpha=0.5,data=dataF1)+
  geom_line(aes(x=Time,y=TempO,colour=Type),alpha=0.5,data=dataF2)+
  geom_line(aes(x=Time,y=TempO,colour=Type),alpha=0.5,data=dataF3)+
  geom_line(aes(x=Time,y=TempO,colour=Type),alpha=0.5,data=dataF4)+
  geom_line(aes(x=Time,y=TempO,colour=Type),alpha=0.5,data=dataF5)+
  geom_line(aes(x=Time,y=TempO,colour=Type),alpha=0.5,data=dataF6)+
  geom_line(aes(x=Time,y=TempO,colour=Type),alpha=0.5,data=dataF7)+
  scale_color_manual(values = c(cbPalette[4],cbPalette[6]))+
  geom_line(aes(x=Time,y=TempO,linetype=Type),data=dataFM)+
  scale_linetype_manual(values=c(1,1),guide="none")+
  geom_line(aes(x=Time,y=TempO),color="#E41A1C",linetype=2,data=dataL[1:60,])+
  geom_line(aes(x=Time,y=TempO),color="#E41A1C",linetype=2,data=dataL[61:(60*2),])+
  xlab("Year")+ylab("Temperature(\u00b0C)")+
  scale_x_continuous(breaks=c(79,91,103,115,127),labels=2021:2025)+
  # scale_x_continuous(breaks=c(183,548,913,1278,1643),labels=c(2020,2040,2060,2080,2100))+
  scale_y_continuous(limits = c(10.5,26))+ # for GO
  # scale_y_continuous(limits = c(-43,30))+  # for GL
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.justification = c(0.99,0.01),
                   legend.position =c(0.99,0.01),
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(0.8,"line"))
PT
# 5.55*3.00

