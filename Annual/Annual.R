################################################################################################
# This file includes all implementations about Annual Case Study (Figs.3,4,5 and Fig.S5,S7,S8) #
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
library(tseries)
library(approxOT)
source("Functions.R")
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###############################################################################################
##### Load the data
YN=86
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
forcing=forcing$total[266:351]   # Year 2015--2100
thetas=read.csv("LENS2_Data/thetas.csv")$x
phis=read.csv("LENS2_Data/phis.csv")$x
id.land=read.csv("LENS2_Data/landid.csv")$x
Dat=array(0,c(10,nrow(Dat.loc),YN))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Annual/dat_em1_year.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Annual/dat_em2_year.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Annual/dat_em3_year.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Annual/dat_em4_year.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Annual/dat_em5_year.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Annual/dat_em6_year.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Annual/dat_em7_year.mat")$Dat
Dat[8,,]=readMat("LENS2_Data/Annual/dat_em8_year.mat")$Dat
Dat[9,,]=readMat("LENS2_Data/Annual/dat_em9_year.mat")$Dat
Dat[10,,]=readMat("LENS2_Data/Annual/dat_em10_year.mat")$Dat
lseq=mseq=0
for(i in 1:143){
  lseq=c(lseq,rep(i,2*i+1))
  mseq=c(mseq,seq(-i,i,by=1))
}



############################################################################################
######################      Model the deterministic component       ######################
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


##### Choose R using I.fit
Dat.hat=matrix(0,nrow(Dat.loc),YN)
I.fit=matrix(0,nrow(Dat.loc),9)
for(R in 2:10){
  ### evaluate the mean trend
  mthat=function(i){
    obj=function(rho){
      value=0
      for(r in 1:R){value=value+profile_negllh(Dat[r,i,],rho)}
      return(value)
    }
    rho.hat=bobyqa(x0=c(0.9),fn=obj,lower=c(0.01),upper=c(0.999))$par
    X=getX(rho.hat)
    beta=solve(t(X)%*%X)%*%t(X)%*%t(Dat[1:R,i,])
    Beta=apply(beta,1,mean)
    return(t(X%*%Beta))
  }
  cl<- makeCluster(10) 
  registerDoParallel(cl) 
  Dat.hat=foreach(i=1:nrow(Dat.loc),
                  .combine=rbind,
                  .packages=c("nloptr")) %dopar% mthat(i)
  stopCluster(cl)
  # writeMat("Annual/Dathat_R5.mat",Dathat=Dat.hat)
  
  ### calculate the sample mean
  Dat.mean=matrix(0,nrow(Dat.loc),YN)
  for(r in 1:R){
    Dat.mean=Dat.mean+Dat[r,,]
  }
  Dat.mean=Dat.mean/R
  
  ### calculate I.fit for each location
  I.fit.up=I.fit.down=rep(0,nrow(Dat.hat))
  for(r in 1:R){
    I.fit.up=I.fit.up+apply((Dat[r,,]-Dat.hat)^2,1,sum)
    I.fit.down=I.fit.down+apply((Dat[r,,]-Dat.mean)^2,1,sum)
  }
  I.fit[,R-1]=I.fit.up/I.fit.down*(R-1)/R
}
write.csv(I.fit,"Annual/Outputs/IfitwithRs.csv")


### Plot Fig.3(a)
dataF=data.frame(Ifit=c(I.fit),
                 Method=as.factor(rep(c(2:10),each=nrow(Dat.hat))))
p=ggplot(data = dataF,aes(x=Method,y=Ifit))+geom_boxplot()+
  geom_hline(yintercept = 1.0,colour="red",linetype=2)+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(1,1),
                   legend.position ="right",
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1,"line"))+
  xlab(expression(R))+ylab(expression(I[fit]))
p
# 5.55*3.00 IfitwithRs.pdf


### Choose R=7
R=7
Dat=Dat[1:7,,]
Dat.hat=readMat("Annual/Dathat_R7.mat")$Dathat
Dat.rsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.rsd[r,,]=Dat[r,,]-Dat.hat
}
Sig=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  Sig[i]=sqrt(sum(Dat.rsd[,i,]^2)/R/YN)
}
write.csv(Sig,"Annual/Outputs/Sig_Annual.csv")

### Plot Fig.3(b)
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
# 5.55*3.00 sigma_demo.pdf

### Test the normality
# Bera=rep(0,nrow(Dat.loc))
# for(j in 1:nrow(Dat.loc)){
#   Bera[j]=jarque.bera.test(c(Dat.rsd[,j,]))$p.value
# }
# idBera=which(Bera<=0.05)
# length(idBera)/nrow(Dat.loc)
# 0.3013419



############################################################################################
################            Model the stochastic component              ###############
### Do SHT with Q=144
Dat.rsd.SHT=array(0,c(R,144^2,YN))
for(r in 1:R){
  Dat.rsdd=Dat.rsd[r,,]/Sig
  
  cl<- makeCluster(4) 
  registerDoParallel(cl) 
  Dat.rsd.SHT[r,,]=foreach(t=1:YN,
                           .combine=cbind,
                           .packages=c("nloptr","QZ","pracma")) %dopar% emsht_forward(t(matrix(Dat.rsdd[,t],288,192)),thetas,phis,144)
  stopCluster(cl)
}
writeMat("Annual/DatrsdshtDesig.mat",Datrsdsht=Dat.rsd.SHT)
# Dat.rsd.SHT=readMat("Annual/DatrsdshtDesig.mat")$Datrsdsht


##### choose Qs with BIC for land and ocean
llseq=c(20,30,40,50,60,70,80,90,100)
BIC.land=BIC.ocean=matrix(0,YN*R,length(llseq))
for(j in 1:length(llseq)){
  L=llseq[j]
  Dat.rsdd=Dat.rsd.hat=matrix(0,nrow(Dat.loc),YN*R)
  for(r in 1:R){
    Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.rsd[r,,]/(Sig%*%t(rep(1,YN)))
    
    cl<- makeCluster(4) 
    registerDoParallel(cl) 
    Dat.rsd.hat[,((r-1)*YN+1):(YN*r)]=foreach(t=1:YN,
                                              .combine=cbind,
                                              .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:L^2,t],thetas,phis,L)))
    stopCluster(cl)
  }
  v2hat=apply((Dat.rsdd-Re(Dat.rsd.hat))^2,1,mean)
  BIC.land[,j]=apply(((Dat.rsdd-Re(Dat.rsd.hat))^2/v2hat)[id.land,],2,sum)+sum(log(v2hat[id.land]))+log(2*pi)*length(id.land)+log(length(id.land))*(L^2) # should be L^2, but no influence
  BIC.ocean[,j]=apply(((Dat.rsdd-Re(Dat.rsd.hat))^2/v2hat)[-id.land,],2,sum)+sum(log(v2hat[-id.land]))+log(2*pi)*(nrow(Dat.loc)-length(id.land))+log(nrow(Dat.loc)-length(id.land))*(L^2)
}
write.csv(BIC.land,"Annual/Outputs/BIC_land.csv")
write.csv(BIC.ocean,"Annual/Outputs/BIC_ocean.csv")

llseq1=60:80
BICd.ocean=matrix(0,YN*R,length(llseq1))
for(j in 1:length(llseq1)){
  L=llseq1[j]
  Dat.rsdd=Dat.rsd.hat=matrix(0,nrow(Dat.loc),YN*R)
  for(r in 1:R){
    Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.rsd[r,,]/(Sig%*%t(rep(1,YN)))
    
    cl<- makeCluster(4) 
    registerDoParallel(cl) 
    Dat.rsd.hat[,((r-1)*YN+1):(YN*r)]=foreach(t=1:YN,
                                              .combine=cbind,
                                              .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:L^2,t],thetas,phis,L)))
    stopCluster(cl)
  }
  v2hat=apply((Dat.rsdd-Re(Dat.rsd.hat))^2,1,mean)
  BICd.ocean[,j]=apply(((Dat.rsdd-Re(Dat.rsd.hat))^2/v2hat)[-id.land,],2,sum)+sum(log(v2hat[-id.land]))+log(2*pi)*(nrow(Dat.loc)-length(id.land))+log(nrow(Dat.loc)-length(id.land))*(L^2)
}
write.csv(BICd.ocean,"Annual/Outputs/BICd_ocean.csv")

llseq2=30:50
BICd.land=matrix(0,YN*R,length(llseq2))
for(j in 1:length(llseq2)){
  L=llseq2[j]
  Dat.rsdd=Dat.rsd.hat=matrix(0,nrow(Dat.loc),YN*R)
  for(r in 1:R){
    Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.rsd[r,,]/(Sig%*%t(rep(1,YN)))
    
    cl<- makeCluster(4) 
    registerDoParallel(cl) 
    Dat.rsd.hat[,((r-1)*YN+1):(YN*r)]=foreach(t=1:YN,
                                              .combine=cbind,
                                              .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:L^2,t],thetas,phis,L)))
    stopCluster(cl)
  }
  v2hat=apply((Dat.rsdd-Re(Dat.rsd.hat))^2,1,mean)
  BICd.land[,j]=apply(((Dat.rsdd-Re(Dat.rsd.hat))^2/v2hat)[id.land,],2,sum)+sum(log(v2hat[id.land]))+log(2*pi)*(length(id.land))+log(length(id.land))*(L^2)
}
write.csv(BICd.ocean,"Annual/Outputs/BICd_land.csv")

### Plot Fig.3(c)
dataF=data.frame(bic=c(c(BIC.land),c(BIC.ocean)),
                 qseq=as.factor(c(rep(llseq,each=R*YN),rep(llseq,each=R*YN))),
                 group=as.factor(c(rep("Land",R*YN*length(llseq)),rep("Ocean",R*YN*length(llseq)))))
p=ggplot(data = dataF,aes(x=qseq,y=bic,color=group))+geom_boxplot()+
  scale_color_manual(values=c("#009E73", "#0072B2"))+
  theme_bw()+theme(panel.border = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.text.y = element_text(angle=90,size=12),
                   axis.title=element_text(size=14),
                   legend.justification=c(0,1),
                   legend.position =c(0,1),
                   legend.title = element_blank(),
                   legend.background = element_rect(colour = "black"),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(0.8,"line"))+
  xlab(expression(Q))+ylab("BIC")
p
dataf1=data.frame(bic=apply(BICd.land,2,median),qseq=llseq2)
dataf2=data.frame(bic=apply(BICd.ocean,2,median),qseq=llseq1)
p+geom_point(mapping = aes(x=qseq/10-1.2,y=bic),data=dataf1,shape=4,col="black",size=0.7)+
  geom_point(mapping = aes(x=qseq/10-1.2,y=bic),data=dataf1[6,],shape=4,col="red",size=0.7)+
  geom_point(mapping = aes(x=qseq/10-0.8,y=bic),data=dataf2,shape=4,col="black",size=0.7)+
  geom_point(mapping = aes(x=qseq/10-0.8,y=bic),data=dataf2[10,],shape=4,col="red",size=0.7)


Ll=35
Dat.rsdd=Dat.rsd.hat=matrix(0,nrow(Dat.loc),YN*R)
for(r in 1:R){
  Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.rsd[r,,]/(Sig%*%t(rep(1,YN)))
  
  cl<- makeCluster(4) 
  registerDoParallel(cl) 
  Dat.rsd.hat[,((r-1)*YN+1):(YN*r)]=foreach(t=1:YN,
                                            .combine=cbind,
                                            .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:Ll^2,t],thetas,phis,Ll)))
  stopCluster(cl)
}
v2hat=apply((Dat.rsdd-Re(Dat.rsd.hat))^2,1,mean)
Lo=69
Dat.rsdd=Dat.rsd.hat=matrix(0,nrow(Dat.loc),YN*R)
for(r in 1:R){
  Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.rsd[r,,]/(Sig%*%t(rep(1,YN)))
  
  cl<- makeCluster(4) 
  registerDoParallel(cl) 
  Dat.rsd.hat[,((r-1)*YN+1):(YN*r)]=foreach(t=1:YN,
                                            .combine=cbind,
                                            .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:Lo^2,t],thetas,phis,Lo)))
  stopCluster(cl)
}
v2hat[-id.land]=apply((Dat.rsdd[-id.land,]-Re(Dat.rsd.hat[-id.land,]))^2,1,mean)
write.csv(v2hat,"Annual/Outputs/v2hat.csv")

### Plot Fig.3(d)
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z=sqrt(v2hat))
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
# 5.55*3.00


### real-valued transformation
L=Lo
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


##### model the temporal dependence: choose P
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
length(which(apply(bicp,1,which.min)==1))/L/L   # 0.9930687
length(which(apply(bicp,1,which.min)==2))/L/L   # 0.002100399
length(which(apply(bicp,1,which.min)==3))/L/L   # 0.001890359
length(which(apply(bicp,1,which.min)==4))/L/L   # 0.0006301197
length(which(apply(bicp,1,which.min)==5))/L/L   # 0.002310439


### Choose P=1
TPhi.hat=rep(0,L^2)
for(i in 1:L^2){
  yy=c(TDat.rsd.SHT[,i,2:YN])
  yx=c(TDat.rsd.SHT[,i,-YN])
  TPhi.hat[i]=t(yx)%*%yy/(t(yx)%*%yx)
}
write.csv(TPhi.hat,"Annual/Outputs/Phihat.csv")

### Plot Fig.4(a)
#TPhi.hat=as.matrix(read.csv("Annual/Phihat.csv")[,-1])
dataF=data.frame(lon=lseq[1:L^2],lat=mseq[1:L^2],z=TPhi.hat)
PT=ggplot()+xlab(expression(q))+ylab(expression(m))+
  geom_point(aes(x=lon,y=lat,colour=z),size=1.8,shape=15,data=dataF)+
  scale_colour_gradient2(low = "#0072B2",mid = "white",high = "#E41A1C",
                         midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
  # scale_color_gradient2(low="white",high="#0072B2")+
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
# sample covariance
TK.sample=matrix(0,L^2,L^2)
for(r in 1:R){
  for(t in 1:YN){
    TK.sample=TK.sample+crossprod(t(TDat.rsd.SHT[r,,t]))}
}
TK.sample=TK.sample/R/YN

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

### Plot Figs.4(b) and (c)
dataF=data.frame(V1=rep(1:L,each=L),V2=rep(seq(-1,-L,by=-1),times=L),V3=c(TK.sample[1:L,1:L]),V4=c(TK.axial[1:L,1:L]))
PT=ggplot()+xlab(" ")+ylab(" ")+
  geom_point(aes(x=V1,y=V2,colour=V3),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#0072B2", mid = "white",high ="#E41A1C",
                         midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_blank(),
                   axis.title = element_blank(),
                   legend.title = element_blank(),
                   legend.position = "right",
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT
# 3.7*3
PT=ggplot()+xlab(" ")+ylab(" ")+
  geom_point(aes(x=V1,y=V2,colour=V4),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#0072B2", mid = "white",high ="#E41A1C",
                         midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_blank(),
                   axis.title = element_blank(),
                   legend.title = element_blank(),
                   legend.position = "right",
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT



##### Calculate auto-covariance
Ang=cbind(rep(thetas,each=288),rep(phis,times=192))
Yfunc=function(para){
  theta=para[1]
  phi=para[2]
  Yseq=rep(0,L^2)
  for(l in 0:(L-1)){
    Leg=legendre(l,cos(theta))
    Reg=rep(0,l+1)
    for(m in 0:l){
      Reg[m+1]=sqrt(factorial(l-m)/factorial(l+m)*(2*l+1)/4/pi)*exp(m*phi*1i)
    }
    Yseq[(l^2+l+1):(l^2+l+l+1)]=Leg*Reg
    for(m in 1:l){
      Yseq[l^2+l-m+1]=(-1)^m*Conj(Yseq[l^2+l+m+1])
    }
  }
  return(Yseq)
}
# id.test=which(round(Dat.loc[,2],digit=1)==36.3)
id.test=which(round(Dat.loc[,2],digit=1)==-11.8)
Ang1=Ang[id.test,]
Y=apply(Ang1,1,Yfunc)

# Empirical
AC.emp=rep(0,287)
for(r in 1:R){
  for(t in 1:YN){
    AC.emp=AC.emp+(Dat.rsd[r,id.test[1:287],t]/Sig[id.test[1:287]])*(Dat.rsd[r,id.test[2:288],t]/Sig[id.test[2:288]])
  }
}
AC.emp=AC.emp/R/YN

# sample
# k.sample=A%*%TK.sample%*%t(A)
# AC.sample=rep(0,287)
# for(i in 1:287){
#   AC.sample[i]=Re(t(Y[,i])%*%k.sample%*%Y[,i+1])
# }

# Axial-nonstationary
k.axial=A%*%TK.axial%*%t(A)
id.test.land=intersect(id.test,id.land)
aa=nabor::knn(id.test.land,id.test,k=1)
aid=which(aa$nn.dists==0)
AC.axialn=rep(0,287)
for(i in 1:287){
  if(length(intersect(c(i,i+1),aid))==0){
    AC.axialn[i]=Re(t(Y[,i])%*%k.axial%*%Y[,i+1])
  }
  if(length(intersect(c(i,i+1),aid))==2){
    AC.axialn[i]=Re(t(Y[1:Ll^2,i])%*%k.axial[1:Ll^2,1:Ll^2]%*%Y[1:Ll^2,i+1])
  }
  if(length(intersect(c(i,i+1),aid))==1){
    bid=intersect(c(i,i+1),aid)
    if(bid==i){
      AC.axialn[i]=Re(t(Y[1:Ll^2,i])%*%k.axial[1:Ll^2,]%*%Y[,i+1])
    }
    if(bid==i+1){
      AC.axialn[i]=Re(t(Y[,i])%*%k.axial[,1:Ll^2]%*%Y[1:Ll^2,i+1])
    }
  }
}


# Axial (Codes from line 527 to line 606 is used to calculate the auto-covariance under the Axial symmetry)
###################
BIC=matrix(0,YN*R,length(llseq1))
for(j in 1:length(llseq1)){
  L=llseq1[j]
  Dat.rsdd=Dat.rsd.hat=matrix(0,nrow(Dat.loc),YN*R)
  for(r in 1:R){
    Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.rsd[r,,]/(Sig%*%t(rep(1,YN)))
    
    cl<- makeCluster(4) 
    registerDoParallel(cl) 
    Dat.rsd.hat[,((r-1)*YN+1):(YN*r)]=foreach(t=1:YN,
                                              .combine=cbind,
                                              .packages=c("nloptr","QZ","pracma")) %dopar% c(t(emsht_inverse(Dat.rsd.SHT[r,1:L^2,t],thetas,phis,L)))
    stopCluster(cl)
  }
  v2hat=apply((Dat.rsdd-Re(Dat.rsd.hat))^2,1,mean)
  BIC[,j]=apply((Dat.rsdd-Re(Dat.rsd.hat))^2/v2hat,2,sum)+sum(log(v2hat))+log(2*pi)*nrow(Dat.loc)+log(nrow(Dat.loc))*L^2
}
write.csv(BIC,"Annual/Outputs/BIC_axial.csv")
#which.min(apply(BIC,2,median))
La=77
TDat.rsd.SHT1=array(0,c(R,La^2,YN))
A1=InvA1=matrix(0,La^2,La^2)
id1=which(mseq[1:La^2]==0)
id2=which(mseq[1:La^2]>0)
for(i in id1){
  A1[i,i]=1
  InvA1[i,i]=1
}
TDat.rsd.SHT1[,id1,]=Re(Dat.rsd.SHT[,id1,])
TDat.rsd.SHT1[,id2,]=Re(Dat.rsd.SHT[,id2,])
for(i in id2){
  l=lseq[i]
  m=mseq[i]
  TDat.rsd.SHT1[,l^2+l-m+1,]=Im(Dat.rsd.SHT[,i,])
  A1[c(l^2+l-m+1,i),c(l^2+l-m+1,i)]=matrix(c((-1)^(m+1)*1i,1i,(-1)^m,1),2,2)
  InvA1[c(l^2+l-m+1,i),c(l^2+l-m+1,i)]=matrix(c(-0.5*1i,-0.5,-0.5*1i,0.5),2,2)
}  

TK.axial1=matrix(0,La^2,La^2)
for(m in 0:(La-1)){
  idm=which(mseq[1:La^2]==m)
  idM=which(mseq[1:La^2]==-m)
  for(r in 1:R){
    for(t in 1:YN){
      TK.axial1[idm,idm]=TK.axial1[idm,idm]+crossprod(t(TDat.rsd.SHT1[r,idm,t]))+
        crossprod(t(TDat.rsd.SHT1[r,idM,t]))
    }
  }
  TK.axial1[idm,idm]=TK.axial1[idm,idm]/R/YN/2
  TK.axial1[idM,idM]=TK.axial1[idm,idm]
}
k.axial1=A1%*%TK.axial1%*%t(A1)

Yfunc1=function(para){
  theta=para[1]
  phi=para[2]
  Yseq=rep(0,La^2)
  for(l in 0:(La-1)){
    Leg=legendre(l,cos(theta))
    Reg=rep(0,l+1)
    for(m in 0:l){
      Reg[m+1]=sqrt(factorial(l-m)/factorial(l+m)*(2*l+1)/4/pi)*exp(m*phi*1i)
    }
    Yseq[(l^2+l+1):(l^2+l+l+1)]=Leg*Reg
    for(m in 1:l){
      Yseq[l^2+l-m+1]=(-1)^m*Conj(Yseq[l^2+l+m+1])
    }
  }
  return(Yseq)
}
Y1=apply(Ang1,1,Yfunc1)

AC.axial=rep(0,287)
for(i in 1:287){
  AC.axial[i]=Re(t(Y1[,i])%*%k.axial1%*%Y1[,i+1])
}
 
####################################################
### Plot Fig4(d) and (e)
dataF=data.frame(lon=rep(phis[-288],times=3),
                 cov=c(AC.emp,AC.axialn,AC.axial),
                 Method=as.factor(rep(c("Empirical","Axial-land/ocean","Axial"),each=287)))
PT=ggplot(dataF,aes(x=180*lon/pi,y=cov,colour=Method,shape=Method))+
  xlab("Longitude")+ylab(" ")+
  geom_point(size=1)+
  scale_color_manual(values=c("#0072B2","#E41A1C","black"))+
  scale_shape_manual(values=c(18,16,17))+
  scale_y_continuous(limits = c(0.3,1.01))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   # legend.justification = c(1,1),
                   # legend.position =c(0.95,0.4),
                   legend.position ="none",
                   legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(1,"line"))

PT



############ Generate emulations 
TU.axial=TK.axial-outer(TPhi.hat,TPhi.hat,"*")*TK.axial
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
# writeMat("Annual/GenDatY.mat",GenDat=Gen.Dat.Y)


##### Calculate I.uq
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
write.csv(c(I.uq),"Annual/Outputs/Iuq_axialnon.csv")

# Plot Figs. 5(a) and 5(b)
I.uq.Huang=read.csv("Supplement/Huang/Annual/Outputs/Iuq.csv")$x
I.uq=read.csv("Annual/Outputs/Iuq_axialnon.csv")$x
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),z=I.uq,z.H=I.uq.Huang)
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 1,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",
                         limits=range(I.uq.Huang))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color=expression(I[uq]))
PT


##### Calculate Wasserstein
WD.time=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  WD.time[i]=wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat.Y[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(WD.time,"Annual/Outputs/WD_time.csv")

WD.space=rep(0,YN)
for(t in 1:YN){
  WD.space[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat.Y[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(WD.space,"Annual/Outputs/WD_space.csv")

### Plot Fig.6(d)
WD.Huang=read.csv("Supplement/Huang/Output/Annual/WD_time.csv")$x
WD.time=read.csv("Annual/Outputs/WD_time.csv")$x
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),z=WD.time,Z.H=WD.Huang)
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradient(low = "#3288BD",high ="#FEE08B",limit=c(0,max(WD.Huang)))+
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

### Computational time
# t1=proc.time()[[3]]
# TGen.Dat.SHT=matrix(0,L^2,YN+1)
# TGen.Dat.SHT[,1]=LL%*%RNM[r,,1]
# for(t in 2:(YN+1)){
#   TGen.Dat.SHT[,t]=LL%*%RNM[r,,t]+TPhi.hat*TGen.Dat.SHT[,t-1]
# }
# TGen.Dat.SHT=TGen.Dat.SHT[,-1]
# Gen.Dat.SHT=A%*%TGen.Dat.SHT
# funcT=function(t){
#   flm=Gen.Dat.SHT[,t]
#   fs=rep(0,nrow(Dat.loc))
#   fs[-id.land]=(Re(c(t(emsht_inverse(flm,thetas,phis,L)))[-id.land])+EPS[r,-id.land])*Sig[-id.land]
#   fs[id.land]=(Re(c(t(emsht_inverse(flm,thetas,phis,Ll)))[id.land])+EPS[r,id.land])*Sig[id.land]
#   return(fs)
# }
# Gen.Dat=sapply(1:YN, funcT)
# a=Gen.Dat+Dat.hat
# t2=proc.time()[[3]]
# 
# # t2-t1=31.558



##### Plots in the Supplementary Materials
#################################
### Plot Fig.S5(a) 
I.fit.montoan=read.csv("Monthly/Aggregate/Ifit_Monthly_to_Annual.csv")$x
I.fit=as.matrix(read.csv("Annual/Outputs/IfitwithRs.csv")[,-1])
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),z=c(I.fit[,R-1]))
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  # geom_point(aes(x=lon,y=lat,colour=z),data=dataF[-ida,],shape=15)+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  geom_point(aes(x=lon,y=lat),data=dataF[c(45461,28236),],shape=4)+
  scale_colour_gradient2(low = "#0072B2",mid = "white",high = "#E41A1C",
                         midpoint = 1,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",
                         limit=c(min(I.fit[,R-1],I.fit.montoan),max(I.fit[,R-1],I.fit.montoan)))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color=expression(I[fit]))
PT+annotate("text",x=13.75,y=10.5,label="G1",size=3)+annotate("text",x=305.00,y=50,label="G2",size=3)


### Plot Fig.S5(b)
dataF=data.frame(Time=rep(2015:2100,times=R),
                 TempL=c(t(Dat[,28236,]))-15,
                 TempO=c(t(Dat[,45461,])),
                 LEN=as.factor(rep(c("1","2","3","4","5","6","7"),each=YN)))
dataL=data.frame(Time=rep(2015:2100,times=4),
                 Temp=c(apply(Dat[,28236,],2,mean)-15,apply(Dat[,45461,],2,mean),Dat.hat[28236,]-15,Dat.hat[45461,]),
                 Type=as.factor(rep(c("Mean on G1-15\u00b0C","Mean on G2",
                                      "Evaluated mean on G1-15\u00b0C","Evaluated mean on G2"),each=YN)))
PT=ggplot()+
  xlab("Year")+ylab("Temperature(\u00b0C)")+
  geom_line(aes(x=Time,y=TempL),colour=cbPalette[1],alpha=0.5,data=dataF[which(dataF$LEN=="1"),])+
  geom_line(aes(x=Time,y=TempO),colour=cbPalette[1],alpha=0.5,data=dataF[which(dataF$LEN=="1"),])+
  geom_line(aes(x=Time,y=TempL),colour=cbPalette[2],alpha=0.5,data=dataF[which(dataF$LEN=="2"),])+
  geom_line(aes(x=Time,y=TempO),colour=cbPalette[2],alpha=0.5,data=dataF[which(dataF$LEN=="2"),])+
  geom_line(aes(x=Time,y=TempL),colour=cbPalette[3],alpha=0.5,data=dataF[which(dataF$LEN=="3"),])+
  geom_line(aes(x=Time,y=TempO),colour=cbPalette[3],alpha=0.5,data=dataF[which(dataF$LEN=="3"),])+
  geom_line(aes(x=Time,y=TempL),colour=cbPalette[4],alpha=0.5,data=dataF[which(dataF$LEN=="4"),])+
  geom_line(aes(x=Time,y=TempO),colour=cbPalette[4],alpha=0.5,data=dataF[which(dataF$LEN=="4"),])+
  geom_line(aes(x=Time,y=TempL),colour=cbPalette[5],alpha=0.5,data=dataF[which(dataF$LEN=="5"),])+
  geom_line(aes(x=Time,y=TempO),colour=cbPalette[5],alpha=0.5,data=dataF[which(dataF$LEN=="5"),])+
  geom_line(aes(x=Time,y=TempL),colour=cbPalette[6],alpha=0.5,data=dataF[which(dataF$LEN=="6"),])+
  geom_line(aes(x=Time,y=TempO),colour=cbPalette[6],alpha=0.5,data=dataF[which(dataF$LEN=="6"),])+
  geom_line(aes(x=Time,y=TempL),colour=cbPalette[7],alpha=0.5,data=dataF[which(dataF$LEN=="7"),])+
  geom_line(aes(x=Time,y=TempO),colour=cbPalette[7],alpha=0.5,data=dataF[which(dataF$LEN=="7"),])+
  geom_line(aes(x=Time,y=Temp,linetype=Type,colour=Type),data=dataL)+
  scale_linetype_manual(values=c(2,1,2,1))+
  scale_colour_manual(values=c("red","red","black","black"))+
  #scale_y_continuous(limits = c(0,28))+
  # scale_color_manual(values = c(cbPalette[1:7]),guide="none")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.justification = c(0.005,0.995),
                   legend.position =c(0.005,0.995),
                   # legend.position = "bottom",
                   legend.title = element_blank(),
                   legend.key.width=unit(0.8,"line"),
                   legend.key.height=unit(0.7,"line"))

PT+annotate("text",x=2100,y=14,label="G1")+annotate("text",x=2100,y=3.5,label="G2")


### Compare with Figure 1, plot Figs. S7(a)--S7(d)
Gen.Dat.Y=readMat("Annual/GenDatY.mat")$GenDat
Gen.Dat.mean=apply(Gen.Dat.Y[,,2083-2014],2,mean) 
Gen.Dat.sd=apply(Gen.Dat.Y[,,2083-2014],2,sd)
Dat.mean=apply(Dat[,,2083-2014],2,mean) 
Dat.sd=apply(Dat[,,2083-2014],2,sd)
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degree=abs(Gen.Dat.mean-Dat.mean))
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degree),size=0.8,data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=4,data=dataF[c(39780,18281),])+
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

dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degree=abs(Gen.Dat.sd-Dat.sd))
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degree),size=0.8,data=dataF)+
  geom_point(aes(x=lon,y=lat),shape=4,data=dataF[c(39780,18281),])+
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


### Plot Figs. S7(e) and S7(f)
Dat.hat=readMat("Annual/Dathat_R7.mat")$Dathat
dataF1=data.frame(Time=rep(2015:2100,times=2),
                  TempO=c(Dat[1,18281,]-3,Gen.Dat.Y[1,18281,]),
                  TempL=c(Dat[1,39780,]-5,Gen.Dat.Y[1,39780,]),
                  Type=as.factor(rep(c("Simulations-3\u00b0C","Emulations"),each=86)))
dataF2=data.frame(Time=rep(2015:2100,times=2),
                  TempO=c(Dat[2,18281,]-3,Gen.Dat.Y[2,18281,]),
                  TempL=c(Dat[2,39780,]-5,Gen.Dat.Y[2,39780,]),
                  Type=as.factor(rep(c("Simulations-3\u00b0C","Emulations"),each=86)))
dataF3=data.frame(Time=rep(2015:2100,times=2),
                  TempO=c(Dat[3,18281,]-3,Gen.Dat.Y[3,18281,]),
                  TempL=c(Dat[3,39780,]-5,Gen.Dat.Y[3,39780,]),
                  Type=as.factor(rep(c("Simulations-3\u00b0C","Emulations"),each=86)))
dataF4=data.frame(Time=rep(2015:2100,times=2),
                  TempO=c(Dat[4,18281,]-3,Gen.Dat.Y[4,18281,]),
                  TempL=c(Dat[4,39780,]-5,Gen.Dat.Y[4,39780,]),
                  Type=as.factor(rep(c("Simulations-3\u00b0C","Emulations"),each=86)))
dataF5=data.frame(Time=rep(2015:2100,times=2),
                  TempO=c(Dat[5,18281,]-3,Gen.Dat.Y[5,18281,]),
                  TempL=c(Dat[5,39780,]-5,Gen.Dat.Y[5,39780,]),
                  Type=as.factor(rep(c("Simulations-3\u00b0C","Emulations"),each=86)))
dataF6=data.frame(Time=rep(2015:2100,times=2),
                  TempO=c(Dat[6,18281,]-3,Gen.Dat.Y[6,18281,]),
                  TempL=c(Dat[6,39780,]-5,Gen.Dat.Y[6,39780,]),
                  Type=as.factor(rep(c("Simulations-3\u00b0C","Emulations"),each=86)))
dataF7=data.frame(Time=rep(2015:2100,times=2),
                  TempO=c(Dat[7,18281,]-3,Gen.Dat.Y[7,18281,]),
                  TempL=c(Dat[7,39780,]-5,Gen.Dat.Y[7,39780,]),
                  Type=as.factor(rep(c("Simulations-3\u00b0C","Emulations"),each=86)))
dataFM=data.frame(Time=rep(2015:2100,times=2),
                  TempO=c(apply(Dat[,18281,],2,mean)-3,apply(Gen.Dat.Y[,18281,],2,mean)),
                  TempL=c(apply(Dat[,39780,],2,mean)-5,apply(Gen.Dat.Y[,39780,],2,mean)),
                  Type=as.factor(rep(c("Simulations-3\u00b0C","Emulations"),each=86)))
dataL=data.frame(Time=rep(2015:2100,times=2),
                 TempO=c(Dat.hat[18281,]-3,Dat.hat[18281,]),
                 TempL=c(Dat.hat[39780,]-5,Dat.hat[39780,]),
                 Type=as.factor(rep(c("Simulations-3\u00b0C","Emulations"),each=86)))
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
  geom_line(aes(x=Time,y=TempO),color="#E41A1C",linetype=2,data=dataL[1:86,])+
  geom_line(aes(x=Time,y=TempO),color="#E41A1C",linetype=2,data=dataL[87:(86*2),])+
  xlab("Year")+ylab("Temperature(\u00b0C)")+
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


#### Comparison of periodograms, plot Fig. S8
Dat.srsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.srsd[r,,]=(Dat[r,,]-Dat.hat)/Sig
}
Dat.Gen.srsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.Gen.srsd[r,,]=(Gen.Dat.Y[r,,]-Dat.hat)/Sig
}

id.test=which(round(Dat.loc[,2],digit=1)==-11.8)
# id.test=which(round(Dat.loc[,2],digit=1)==36.3)
Peridogram.sim=rep(0,144)
for(r in 1:R){
  for(t in 1:YN){
    Peridogram.sim=Peridogram.sim+(abs(fft(Dat.srsd[r,id.test,t]))[1:144])^2
  }
}
Peridogram.emu=rep(0,144)
for(r in 1:R){
  for(t in 1:YN){
    Peridogram.emu=Peridogram.emu+(abs(fft(Dat.Gen.srsd[r,id.test,t]))[1:144])^2
  }
}
dataF=data.frame(Wavenumber=rep(1:144,times=2),
                 Periodogram=c(log(Peridogram.sim/288/YN/R),log(Peridogram.emu/288/YN/R)),
                 Type=as.factor(rep(c("Simulations","Emulations"),each=144)))
PT=ggplot()+xlab("Wave number")+ylab("log(Periodogram)")+
  geom_point(aes(x=Wavenumber,y=Periodogram,colour=Type,shape=Type),size=0.8,data=dataF)+
  geom_line(aes(x=Wavenumber,y=Periodogram,colour=Type),data=dataF)+
  scale_color_manual(values = c("#E41A1C","black"))+
  scale_shape_manual(values=c(16,17))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.justification = c(0.99,0.99),
                   legend.position ="none",
                   legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT













