##################################################################################################
# This file includes all implementations about Daily Case Study (Figs.6 and 7, Figs.S14 and S15) #
##################################################################################################
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




############################################################################################
######################      Model the deterministic component       ######################
getX=function(rho,M){
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

#########################################################
##### Choose M
set.seed(100)
locrandid=sample(1:55296,round((55296/5)))
write.csv("Daily/Outputs/locrandid.csv")

M=1
cl<- makeCluster(8)
registerDoParallel(cl)
Res.hatrho1= foreach(i=locrandid,
                     .combine=cbind,
                     .packages=c("nloptr")) %dopar% hat.rho(i)
stopCluster(cl)
write.csv(Res.hatrho1,"Res_hatrho1.csv")

Res.hatrho=rep(0,55296)
Res.hatrho[locrandid]=Res.hatrho1
Dathat=function(i){
  X=getX(Res.hatrho[i],M)
  beta=solve(t(X)%*%X)%*%t(X)%*%t(Dat[,i,])
  Beta=apply(beta,1,mean)
  meanhat=X%*%Beta
  Sigmahat=Dat[,i,]-rep(1,R)%*%t(meanhat)
  sigmahat=mean(Sigmahat^2)
  bicvalue=(2*M+3)*log(R*YN)+R*YN*log(sigmahat)
  return(bicvalue)
}
cl<- makeCluster(8)
registerDoParallel(cl)
Bic1= foreach(i=locrandid,
              .combine=cbind,
              .packages=c("nloptr")) %dopar% Dathat(i)
stopCluster(cl)
Bic=cbind(Bic,c(Bic5))
write.csv(Bic,"Daily/Outputs/BicforM.csv")



#########################################################
##### Evaluate the deterministic component with M=4
M=4
cl<- makeCluster(6)
registerDoParallel(cl)
Res.hatrho= foreach(i=1:nrow(Dat.loc),
                    .combine=cbind,
                    .packages=c("nloptr")) %dopar% hat.rho(i)
stopCluster(cl)
write.csv(Res.hatrho,"Daily/Outputs/Res_Hatrho.csv")

Dathat=function(i){
  X=getX(Res.hatrho[i],M=4)
  beta=solve(t(X)%*%X)%*%t(X)%*%t(Dat[,i,])
  Beta=apply(beta,1,mean)
  meanhat=X%*%Beta
  Sigmahat=Dat[,i,]-rep(1,R)%*%t(meanhat)
  sigmahat=mean(Sigmahat^2)
  return(c(meanhat,sigmahat))
}
cl<- makeCluster(8)
registerDoParallel(cl)
Dat.hat= foreach(i=1:nrow(Dat.loc),
                 .combine=cbind,
                 .packages=c("nloptr")) %dopar% Dathat(i)
stopCluster(cl)
Sig=sqrt(Dat.hat[1826,])
Dat.hat=t(Dat.hat[1:1825,])
# writeMat("Daily/Dathat.mat",Dathat=Dat.hat)
write.csv(Sig,"Daily/Outputs/Daily_Sig.csv")

### Plot Fig.S14(a)
Dat.mean=(Dat[1,,365*1+12*8+13]+Dat[2,,365*1+12*8+13]+Dat[3,,365*1+12*8+13]+Dat[4,,365*1+12*8+13]+
            Dat[5,,365*1+12*8+13]+Dat[6,,365*1+12*8+13]+Dat[7,,365*1+12*8+13])/7
Dat.hat=readMat("Daily/Dathat.mat")$Dat.hat
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),z=abs(Dat.mean-Dat.hat[,365*1+12*8+13]))
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


### Plot Fig.S14(b)
# Dat.hat=readMat("Daily/Dathat.mat")$Dathat
# Sig=read.csv("Daily/Outputs/Daily_Sig.csv")$x
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

### do SHT with Q=144
lseq=mseq=0
for(i in 1:143){
  lseq=c(lseq,rep(i,2*i+1))
  mseq=c(mseq,seq(-i,i,by=1))
}
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
writeMat("DatrsdshtDesig7.mat",Datrsdsht=Dat.rsd.SHT[7,,])
Dat.rsd.SHT[1,,]=readMat("Daily/DatrsdshtDesig1.mat")$Datrsdsht


### Choose Q
llseq=c(20,30,40,50,60,70,80,90)
BIC.land=BIC.ocean=matrix(0,YN*R,length(llseq))
for(j in 1:length(llseq)){
  L=llseq[j]
  Dat.err=array(0,c(R,nrow(Dat.loc),YN))
  for(r in 1:R){
    cl<- makeCluster(8) 
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
write.csv(BIC.land,"Daily/Outputs/BIC_land.csv")
write.csv(BIC.ocean,"Daily/Outputs/BIC_ocean.csv")
llseq1=65:75
BICd.ocean=matrix(0,YN*R,length(llseq1))
for(j in 2:length(llseq1)){
  L=llseq1[j]
  Dat.err=array(0,c(R,nrow(Dat.loc),YN))
  for(r in 1:R){
    cl<- makeCluster(8) 
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
write.csv(BICd.ocean,"Daily/Outputs/BICd_ocean.csv")
llseq2=30:40
BICd.land=matrix(0,YN*R,length(llseq2))
for(j in 1:length(llseq2)){
  L=llseq2[j]
  Dat.err=array(0,c(R,nrow(Dat.loc),YN))
  for(r in 1:R){
    cl<- makeCluster(8) 
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
write.csv(BICd.land,"Daily/Outputs/BICd_land.csv")

### Plot Fig.S14(c)
dataF=data.frame(bic=c(c(BIC.land),c(BIC.ocean)),
                 qseq=as.factor(c(rep(llseq,each=1825*R),rep(llseq,each=1825*R))),
                 group=as.factor(c(rep("Land",YN*R*length(llseq)),rep("Ocean",YN*R*length(llseq)))))
p=ggplot(data = dataF,aes(x=qseq,y=bic,color=group))+geom_boxplot()+
  scale_color_manual(values=c("#009E73", "#0072B2"))+
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
  geom_point(mapping = aes(x=qseq/10-0.8,y=bic),data=dataf2[4,],shape=4,col="red",size=0.7)


Ll=36
L=Lo=68
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
  
  cl<- makeCluster(8) 
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

### Plot Fig.S14(d)
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


### real-valued transformation
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
length(idBera)/L/L   # 0.5668253



################################################################################
##########   Case 1: no Tukey g-and-h
### choose P
bicpfunc=function(i){
  bicp=rep(0,5)
  ts=TDat.rsd.SHT[,i,]
  for(p in 1:5){
    u2=rep(0,R)
    for(r in 1:R){
      xx=rep(0,(365-p)*5)
      for(tt in p:1){
        xx=cbind(xx,TDat.rsd.SHT[r,i,c(tt:(tt+365-p-1),(tt+365):(tt+365*2-p-1),(tt+365*2):(tt+365*3-p-1),
                                       (tt+365*3):(tt+365*4-p-1),(tt+365*4):(tt+365*5-p-1))])
      }
      xx=as.matrix(xx[,-1])
      yy=as.matrix(ts[r,c((p+1):365,(p+366):(365*2),(p+365*2+1):(365*3),(p+365*3+1):(365*4),(p+365*4+1):(365*5))])
      u2[r]=t(yy)%*%(diag(1,(365-p)*5)-xx%*%solve(t(xx)%*%xx)%*%t(xx))%*%yy/((365-p)*5)
    }
    u2=mean(u2)
    bicp[p]=p*log((365-p)*5*R)+R*(365-p)*5*log(2*pi)+R*(365-p)*5*log(u2)
  }
  return(bicp)
}
cl<- makeCluster(8) 
registerDoParallel(cl) 
bicp=foreach(i=1:(L^2),
             .combine=rbind,
             .packages=c("moments","nloptr")) %dopar% bicpfunc(i)
stopCluster(cl)
write.csv(bicp,"Daily/Outputs/bicp_noTukey.csv")


### Choose P=1 because we need to compare with the Tukey case 
p=1
TPhi.hat=matrix(0,L^2,p)
for(i in 1:L^2){
  yy=c(t((TDat.rsd.SHT[,i,c((p+1):365,(p+366):(365*2),(p+365*2+1):(365*3),(p+365*3+1):(365*4),(p+365*4+1):(365*5))])))
  yx=matrix(0,(365-p)*5*R,p)
  for(r in 1:R){
    xx=rep(0,(365-p)*5)
    for(tt in p:1){
      xx=cbind(xx,TDat.rsd.SHT[r,i,c(tt:(tt+365-p-1),(tt+365):(tt+365*2-p-1),(tt+365*2):(tt+365*3-p-1),
                                     (tt+365*3):(tt+365*4-p-1),(tt+365*4):(tt+365*5-p-1))])
    }
    yx[((365-p)*5*(r-1)+1):((365-p)*5*r),]=as.matrix(xx[,-1])
  }
  TPhi.hat[i,]=solve(t(yx)%*%yx)%*%t(yx)%*%yy
}
write.csv(TPhi.hat,"Daily/Outputs/Phihat_noTukey.csv")

### Model the spatial dependence K under the axial symmetry
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
TU.axial=TK.axial-outer(TPhi.hat[,1],TPhi.hat[,1],"*")*TK.axial
ares=svd(TU.axial)
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
    TGen.Dat.SHT=matrix(0,L^2,365+1)
    TGen.Dat.SHT[,1]=LL%*%RNM[r,,365*(m-1)+1]
    for(t in 2:(365+1)){
      TGen.Dat.SHT[,t]=LL%*%RNM[r,,365*(m-1)+t]+TPhi.hat[,1]*TGen.Dat.SHT[,t-1]
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
    Gen.Dat=foreach(i=1:365,
                    .combine = cbind,
                    .packages = c("QZ","pracma")) %dopar% funcT(i)
    stopCluster(cl)
    Gen.Dat.Y[r,,(365*(m-1)+1):(365*m)]=Gen.Dat+Dat.hat[,(365*(m-1)+1):(365*m)]
  }
}
# writeMat("Daily/GenDatnoTukey7.mat",GenDat=Gen.Dat.Y[7,,])
# Gen.Dat.Y[7,,]=readMat("Daily/GenDatnoTukey7.mat")$GenDat

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
write.csv(c(I.uq),"Daily/Outputs/Iuq_noTukey.csv")


### Calculate Wasserstein
WDtime=function(i){
  return(wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat.Y[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate"))
}
cl=makeCluster(4)
registerDoParallel(cl)
WD.time.Tukey=foreach(i=1:nrow(Dat.loc),
                      .combine = cbind,
                      .packages = c("approxOT")) %dopar% WDtime(i)
stopCluster(cl)
write.csv(c(WD.time.Tukey),"Daily/Outputs/WD_time_noTukey.csv")


WD.space.Tukey=rep(0,YN)
for(t in 1:YN){
  WD.space.Tukey[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat.Y[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(WD.space.Tukey,"Daily/Outputs/WD_space_noTukey.csv")



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
        xx=rep(0,(365-p)*5)
        for(tt in p:1){
          xx=cbind(xx,tildez[r,c(tt:(tt+365-p-1),(tt+365):(tt+365*2-p-1),(tt+365*2):(tt+365*3-p-1),
                                 (tt+365*3):(tt+365*4-p-1),(tt+365*4):(tt+365*5-p-1))])
        }
        xx=as.matrix(xx[,-1])
        u2[r]=mean((tildez[r,c((p+1):365,(p+366):(365*2),(p+365*2+1):(365*3),(p+365*3+1):(365*4),(p+365*4+1):(365*5))]-xx%*%phi)^2)
      }
      u2=mean(u2)
      res=sum(Phi)-nrow(y)*((365-p)*5)/2*log(2*pi)-nrow(y)*((365-p)*5)/2*log(u2)
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
      xx=rep(0,(365-p)*5)
      for(tt in p:1){
        xx=cbind(xx,tildez[r,c(tt:(tt+365-p-1),(tt+365):(tt+365*2-p-1),(tt+365*2):(tt+365*3-p-1),
                               (tt+365*3):(tt+365*4-p-1),(tt+365*4):(tt+365*5-p-1))])
      }
      xx=as.matrix(xx[,-1])
      u2[r]=mean((tildez[r,c((p+1):365,(p+366):(365*2),(p+365*2+1):(365*3),(p+365*3+1):(365*4),(p+365*4+1):(365*5))]-xx%*%Paraest[-(1:3)])^2)
    }
    u2=mean(u2)
    bicp[p]=-2*(sum(Phi)-nrow(y)*((365-p)*5)/2*log(2*pi)-nrow(y)*((365-p)*5)/2*log(u2))+(p+3)*log(nrow(y)*(365-p)*5)
  }
  return(bicp)
}
cl<- makeCluster(8) 
registerDoParallel(cl) 
bicp1=foreach(j=idBera,
              .combine=rbind,
              .packages=c("moments","nloptr")) %dopar% bicpfunc1(j)
stopCluster(cl)
write.csv(bicp1,"Daily/Outputs/bicp_Tukey.csv")
bicp=as.matrix(read.csv("Daily/Outputs/bicp_noTukey.csv")[,-1])
bicp[idBera,]=as.matrix(bicp1)
length(which(apply(bicp,1,which.min)==1))/L/L  # 0.5549308
length(which(apply(bicp,1,which.min)==2))/L/L  # 0.07352941
length(which(apply(bicp,1,which.min)==3))/L/L  # 0.3551038
length(which(apply(bicp,1,which.min)==4))/L/L  # 0.01643599
length(which(apply(bicp,1,which.min)==5))/L/L  # 0

### Plot Fig.S14(e)
dataF=data.frame(lon=c(lseq[1:L^2],5),lat=c(mseq[1:L^2],0),z=c(apply(bicp,1,which.min),5))
PT=ggplot()+xlab(expression(q))+ylab(expression(m))+
  geom_point(aes(x=lon,y=lat,colour=z),size=1.8,shape=15,data=dataF[1:4624,])+
  geom_point(aes(x=lon,y=lat,colour=z),size=0.1,alpha=0.1,shape=15,data=dataF[-(1:4624),])+
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

### Plot Fig.S14(f)
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
# writeMat("Daily/GenDatTukey7.mat",GenDat=Gen.Dat.Y[7,,])
# Gen.Dat.Y[1,,]=readMat("Daily/GenDatTukey1.mat")$GenDat


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
write.csv(c(I.uq.Tukey),"Daily/Outputs/Iuq_Tukey.csv")


### Calculate Wasserstein
WDtime=function(i){
  return(wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat.Y[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate"))
}
cl=makeCluster(4)
registerDoParallel(cl)
WD.time=foreach(i=1:nrow(Dat.loc),
                .combine = cbind,
                .packages = c("approxOT")) %dopar% WDtime(i)
stopCluster(cl)
write.csv(c(WD.time),"Daily/Outputs/WD_time_Tukey.csv")

WD.space=rep(0,YN)
for(t in 1:YN){
  WD.space[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat.Y[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(WD.space,"Daily/Outputs/WD_space_Tukey.csv")


### Plot Figs.6(a) and 6(b)
I.uq=read.csv("Daily/Outputs/Iuq_Tukey.csv")$x
I.uq.Huang=read.csv("Supplement/Huang/Daily/Outputs/Iuq_Huang.csv")$x
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z.H=I.uq.Huang,z=I.uq)
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradient2(low ="#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 1,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",limit=c(min(dataF$z.H),max(dataF$z)))+
  geom_point(aes(x=lon,y=lat),data=dataF[c(7785,10377,18153,26793),],shape=4)+
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

### Plot Figs.6(c) and 6(d)
WD.time=read.csv("Daily/Outputs/WD_time_Tukey.csv")$x
WD.Huang=read.csv("Supplement/Huang/Daily/Outputs/WD_Huang.csv")$x
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z=c(WD.time),z.H=c(WD.Huang))
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradient(low = "#3288BD",high ="#FEE08B",limit=c(0,max(WD.time)))+
  geom_point(aes(x=lon,y=lat),data=dataF[c(7785,10377,18153,26793),],shape=4)+
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

### Plot Figs.6(e) and 6(f)
I.uq.noTukey=read.csv("Daily/Outputs/Iuq_noTukey.csv")$x
WD.noTukey=read.csv("Daily/Outputs/WD_time_noTukey.csv")$x
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


### Plot Fig.7
writeMat("Daily/A.mat",A=Dat[,c(7785,10377,18153,26793),])
writeMat("Daily/B.mat",B=Gen.Dat.Y[,c(7785,10377,18153,26793),])   
A=readMat("Daily/A.mat")$A
B=readMat("Daily/B.mat")$B
dataF1=data.frame(Time=rep(1:1825,times=2),
                  TempB=c(A[1,1,]-20,B[1,1,]),TempN=c(A[1,2,]-10,B[1,2,]),
                  TempO=c(A[1,3,]-10,B[1,3,]),TempL=c(A[1,4,]-10,B[1,4,]),
                  Type=as.factor(rep(c("Simulations-10\u00b0C","Emulations"),each=1825)))
dataF2=data.frame(Time=rep(1:1825,times=2),
                  TempB=c(A[2,1,]-20,B[2,1,]),TempN=c(A[2,2,]-10,B[2,2,]),
                  TempO=c(A[2,3,]-10,B[2,3,]),TempL=c(A[2,4,]-10,B[2,4,]),
                  Type=as.factor(rep(c("Simulations-10\u00b0C","Emulations"),each=1825)))
dataF3=data.frame(Time=rep(1:1825,times=2),
                  TempB=c(A[3,1,]-20,B[3,1,]),TempN=c(A[3,2,]-10,B[3,2,]),
                  TempO=c(A[3,3,]-10,B[3,3,]),TempL=c(A[3,4,]-10,B[3,4,]),
                  Type=as.factor(rep(c("Simulations-10\u00b0C","Emulations"),each=1825)))
dataF4=data.frame(Time=rep(1:1825,times=2),
                  TempB=c(A[4,1,]-20,B[4,1,]),TempN=c(A[4,2,]-10,B[4,2,]),
                  TempO=c(A[4,3,]-10,B[4,3,]),TempL=c(A[4,4,]-10,B[4,4,]),
                  Type=as.factor(rep(c("Simulations-10\u00b0C","Emulations"),each=1825)))
dataF5=data.frame(Time=rep(1:1825,times=2),
                  TempB=c(A[5,1,]-20,B[5,1,]),TempN=c(A[5,2,]-10,B[5,2,]),
                  TempO=c(A[5,3,]-10,B[5,3,]),TempL=c(A[5,4,]-10,B[5,4,]),
                  Type=as.factor(rep(c("Simulations-10\u00b0C","Emulations"),each=1825)))
dataF6=data.frame(Time=rep(1:1825,times=2),
                  TempB=c(A[6,1,]-20,B[6,1,]),TempN=c(A[6,2,]-10,B[6,2,]),
                  TempO=c(A[6,3,]-10,B[6,3,]),TempL=c(A[6,4,]-10,B[6,4,]),
                  Type=as.factor(rep(c("Simulations-10\u00b0C","Emulations"),each=1825)))
dataF7=data.frame(Time=rep(1:1825,times=2),
                  TempB=c(A[7,1,]-20,B[7,1,]),TempN=c(A[7,2,]-10,B[7,2,]),
                  TempO=c(A[7,3,]-10,B[7,3,]),TempL=c(A[7,4,]-10,B[7,4,]),
                  Type=as.factor(rep(c("Simulations-10\u00b0C","Emulations"),each=1825)))
dataFM=data.frame(Time=rep(1:1825,times=2),
                  TempB=c(apply(A[,1,],2,mean)-20,apply(B[,1,],2,mean)),
                  TempN=c(apply(A[,2,],2,mean)-10,apply(B[,2,],2,mean)),
                  TempO=c(apply(A[,3,],2,mean)-10,apply(B[,3,],2,mean)),
                  TempL=c(apply(A[,4,],2,mean)-10,apply(B[,4,],2,mean)),
                  Type=as.factor(rep(c("Simulations-10\u00b0C","Emulations"),each=1825)))
dataL=data.frame(Time=rep(1:1825,times=2),
                 TempB=c(Dat.hat[7785,]-20,Dat.hat[7785,]),
                 TempN=c(Dat.hat[10377,]-10,Dat.hat[10377,]),
                 TempO=c(Dat.hat[18153,]-10,Dat.hat[18153,]),
                 TempL=c(Dat.hat[27081,]-10,Dat.hat[27081,]),
                 Type=as.factor(rep(c("Simulations-10\u00b0C","Emulations"),each=1825)))

PT=ggplot()+
  geom_line(aes(x=Time,y=TempL,colour=Type),alpha=0.5,data=dataF1)+
  geom_line(aes(x=Time,y=TempL,colour=Type),alpha=0.5,data=dataF2)+
  geom_line(aes(x=Time,y=TempL,colour=Type),alpha=0.5,data=dataF3)+
  geom_line(aes(x=Time,y=TempL,colour=Type),alpha=0.5,data=dataF4)+
  geom_line(aes(x=Time,y=TempL,colour=Type),alpha=0.5,data=dataF5)+
  geom_line(aes(x=Time,y=TempL,colour=Type),alpha=0.5,data=dataF6)+
  geom_line(aes(x=Time,y=TempL,colour=Type),alpha=0.5,data=dataF7)+
  scale_color_manual(values = c(cbPalette[4],cbPalette[6]))+
  geom_line(aes(x=Time,y=TempL,linetype=Type),data=dataFM)+
  scale_linetype_manual(values=c(1,1),guide="none")+
  geom_line(aes(x=Time,y=TempL),color="#E41A1C",linetype=2,data=dataL[1:1825,])+
  geom_line(aes(x=Time,y=TempL),color="#E41A1C",linetype=2,data=dataL[1826:(1825*2),])+
  xlab("Year")+ylab("Temperature(\u00b0C)")+
  geom_vline(xintercept = c(365,365*2,365*3,365*4),linetype=3)+
  scale_x_continuous(breaks=c(183,548,913,1278,1643),labels=c(2020,2040,2060,2080,2100))+
  # scale_y_continuous(limits = c(1,25))+ # for GO
  scale_y_continuous(limits = c(8,35))+  # for GL
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.justification = c(0.99,0.01),
                   legend.position =c(0.99,0.01),
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))
PT
# 5.55*3.00

### Plot Fig.S15
dataF=data.frame(Time=rep(1:1825,times=7),
                 TempL=c(t(Dat[,10377,])),
                 TempO=c(t(Dat[,10477,])),
                 LEN=as.factor(rep(c("1","2","3","4","5","6","7"),each=1825)))
PT=ggplot()+
  xlab("Year")+ylab("Temperature(\u00b0C)")+
  geom_line(aes(x=Time,y=TempL,colour=LEN),alpha=0.5,data=dataF)+
  geom_line(aes(x=Time,y=TempO,colour=LEN),alpha=0.5,data=dataF)+
  scale_color_manual(values = c(cbPalette[1:7]),guide="none")+
  geom_vline(xintercept = c(365,365*2,365*3,365*4),linetype=3)+
  scale_x_continuous(breaks=c(183,548,913,1278,1643),labels=c(2020,2040,2060,2080,2100))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   # legend.justification = c(0.99,0.01),
                   # legend.position =c(0.99,0.01),
                   legend.position = "none",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))
PT+annotate("text",x=1800,y=9,label="TGN'")+annotate("text",x=1800,y=0.5,label="TGN")


dataF=data.frame(Time=rep(1:1825,times=7),
                 TempL=c(t(Dat[,7785,]))-20,
                 TempO=c(t(Dat[,51085,]))+20,
                 LEN=as.factor(rep(c("1","2","3","4","5","6","7"),each=1825)))
PT=ggplot()+
  xlab("Year")+ylab("Temperature(\u00b0C)")+
  geom_line(aes(x=Time,y=TempL,colour=LEN),alpha=0.5,data=dataF)+
  geom_line(aes(x=Time,y=TempO,colour=LEN),alpha=0.5,data=dataF)+
  scale_color_manual(values = c(cbPalette[1:7]),guide="none")+
  geom_vline(xintercept = c(365,365*2,365*3,365*4),linetype=3)+
  scale_x_continuous(breaks=c(183,548,913,1278,1643),labels=c(2020,2040,2060,2080,2100))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   # legend.justification = c(0.99,0.01),
                   # legend.position =c(0.99,0.01),
                   legend.position = "none",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))
PT+annotate("text",x=1700,y=-25,label="TGB'")+annotate("text",x=1700,y=5,label="TGNP")




