##############  Wrapper file for sequentially implementing each figure in the Supplementary Materials ##############
###### Load necessary R packages and functions
# R version 3.6.3
library(R.matlab)
library(nloptr)
library(fpp2)
library(irlba)
library(fdaoutlier)
library(QZ)
library(pracma)
library(foreach)
library(doParallel)
library(moments)
library(tseries)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(LatticeKrig)
library(approxOT)
library(mvtnorm)
library(reticulate)
source("Functions.R")
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

########     Step 1. Download and process the data in sub-repository "LENS2_Data"   ########
# Please follow the README.md file in the sub-repository "LENS2_Data" to download          #
# the raw data, process and store them using file “Data_treatment.R” in the "LENS2_Data".  # 
#   Processing montly and annual data together for each ensemble takes about 63 seconds.   #
#   Processing daily data for each ensemble takes about 32.5 seconds.                      #
# Assume that we have downloaded the raw data, processed and stored them in sub-repository #
# "LENS2_Data/Annual", "LENS2_Data/Monthly", and "LENS2_Data/Daily".                       #
############################################################################################



#########         Step 2. Reproduce each section and figure sequentially           #########
# Please refer to the README_Supplement.md file for more details                           #
############################################################################################

######   Figure S1 in Section S2   #############################################
# Please refer to the reproduction of Figure 1 in Section 2                    #
################################################################################


######   Figure S2 in Section S3.2.3   #########################################
# Please refer to the reproduction of Figure 2 in Section 3.2.1                #
################################################################################


######   Figure S3 in Section S3.2.4   #########################################
# Figure S3 illustrates the approximation performance of LatticeKrig and SHT.  #
# Figures S3(d) and S3(b) are Figures 2(e) and 2(f), respectively. Therefore,  #
# we would not repeat these two figures.                                       #
# All intermediate outputs are in sub-repository "Supplement/FigureS3/Outputs".#
################################################################################
### Load the annual data and necessary information
R=7
YN=86
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
Dat.loc1=cbind(rep(c(seq(0,180,by=1.25),seq(-178.75,-1.25,by=1.25)),times=192),Dat.loc[,2])
thetas=read.csv("LENS2_Data/thetas.csv")$x
phis=read.csv("LENS2_Data/phis.csv")$x
Dat=array(0,c(7,nrow(Dat.loc),YN))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Annual/dat_em1_year.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Annual/dat_em2_year.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Annual/dat_em3_year.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Annual/dat_em4_year.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Annual/dat_em5_year.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Annual/dat_em6_year.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Annual/dat_em7_year.mat")$Dat
lseq=mseq=0
for(i in 1:(144-1)){
  lseq=c(lseq,rep(i,2*i+1))
  mseq=c(mseq,seq(-i,i,by=1))
}
# 5.182

### Calculate one set of stochastic component Z_9^(1)(L_i,l_j)
Dat.mean=apply(Dat[,,2023-2014],2,mean)
Dat.sd=apply(Dat[,,2023-2014],2,sd)
Dat.S1=(Dat[1,,2023-2014]-Dat.mean)/Dat.sd
# 2.38

### Use LatticeKrig with nlevel=5 to approximate Z_9^(1)(L_i,l_j) and plot Figure S3(a)
# Details of tuning parameters selection can be found in "Supplement/FigureS3/LatticeKrig_Annual.R".
# Related intermediate results "Res_FindNu.csv" are provided.
# t1=proc.time()[[3]]
LKinfol3=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=5,a.wght=5.482,nu=1,LKGeometry="LKSphere") 
LKres3=LatticeKrig(Dat.loc1,Dat.S1,LKinfo=LKinfol3)
# t2=proc.time()[[3]]
# t2-t1=180.13
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(LKres3$residuals))
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degrees),size=0.8,data=dataF)+
  scale_color_gradientn(values = seq(0,1,0.125),colours = Col,limits=c(0,2.974017))+
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

### Use SHT with Q=58 to approximate Z_9^(1)(L_i,l_j) and plot Figure S3(c)
# t1=proc.time()[[3]]
DatSHT=emsht_forward(t(matrix(Dat.S1,288,192)),thetas,phis,144)
DatSHTinv4.E1=emsht_inverse(DatSHT,thetas,phis,58)
# t2=proc.time()[[3]]
# t2-t1=4.233
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(c(t(Re(DatSHTinv4.E1)))-Dat.S1))  # 2.974017
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degrees),size=0.8,data=dataF)+
  scale_color_gradientn(values = seq(0,1,0.125),colours = Col,limit=c(0,2.974017))+
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

### Compare LatticeKrig and SHT using multiple stochastic components and plot Figures S3(e) and S3(f)
# Details of tuning parameters selection can be found in "Supplement/FigureS3/LatticeKrig_Annual.R".
# Related intermediate results "Res_FindLambda.csv" are provided.
# t1=proc.time()[[3]]
LKinfol=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=5,a.wght=5.482,nu=1,lambda=0.000124893,LKGeometry="LKSphere")
LKres=LatticeKrig(Dat.loc1,Dat.srsd[2,,20],LKinfo=LKinfol)      
wX=LKres$wX
wU=LKres$wU
LK.res=function(y){
  LKres=LKrig(Dat.loc1,y,LKinfo=LKinfol,wX=wX,wU=wU)
  return(LKres$residuals)
}
Rest=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  cl<- makeCluster(4)
  registerDoParallel(cl)
  Rest[r,,]=foreach(t=1:YN,
                    .combine=cbind,
                    .packages=c("LatticeKrig")) %dopar% LK.res(Dat.srsd[r,,t])
  stopCluster(cl)
}
SHT.res=function(y){
  DatSHT=emsht_forward(t(matrix(y,288,192)),thetas,phis,144)
  DatSHTinv=emsht_inverse(DatSHT,thetas,phis,58)
  return(y-c(t(Re(DatSHTinv))))
}
Rest.SHT=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  cl<- makeCluster(4)
  registerDoParallel(cl)
  Rest.SHT[r,,]=foreach(t=1:YN,
                        .combine=cbind,
                        .packages=c("pracma","nloptr","QZ")) %dopar% SHT.res(Dat.srsd[r,,t])
  stopCluster(cl)
}
# t2=proc.time()[[3]]
# t2-t1=2377.255 
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),
                 z=apply(apply(abs(Rest),c(2,3),mean),1,mean)-apply(apply(abs(Rest.SHT),c(2,3),mean),1,mean))
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  PT
dataF=data.frame(time=rep(2015:2100,times=2),
                 MSE=c(apply(apply(abs(Rest),c(2,3),mean),2,mean),apply(apply(abs(Rest.SHT),c(2,3),mean),2,mean)),
                 Method=as.factor(rep(c("LatticeKrig","SHT"),each=YN)))
PT=ggplot()+xlab("Year")+ylab("MSE")+
  #geom_point(aes(x=time,y=MSE,colour=type),data=dataF)+
  geom_line(aes(x=time,y=MSE,colour=Method,linetype=Method),data=dataF)+
  scale_color_manual(values=c("#FEE08B","#3288BD"))+
  scale_linetype_manual(values=c(1,1))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.justification = c(1,1),
                   legend.position =c(0.99,0.99),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(1,"line"))
PT



######   Figure S4 in Section S3.4     #########################################
# Figure S4 validates the assumption of diagonal matrix Phi_p.                 #
# Here we take the annual case as an example. Assume that we have got the real-#
# valued SHT coefficients "TDat.rsd.SHT" and their temporal dependence         #
# "TPhi.hat" by following lines 1-192 in "Annual/Annual_SG.R".                 #
################################################################################
### Calculate the residuals of autoregressive model
# t1=proc.time()[[3]]
RTDat.rsd.SHT=array(0,c(R,L^2,YN-1))
for(i in 1:L^2){
  xX=TDat.rsd.SHT[,i,] 
  # xX=CDat.rsd.SHT[,i,] # for monthly data
  RTDat.rsd.SHT[,i,]=xX[,2:YN]-TPhi.hat[i]*xX[,1:(YN-1)]
}
# t2=proc.time()[[3]]
# t2-t1=0.22

### Calculate p-values of the first temporal lag of the cross-correlation and plot Figure S4(a)
# t1=proc.time()[[3]]
crosfunc=function(x){
  i=x[1]
  j=x[2]
  a=0
  for(r in 1:R){
    a=a+ccf(RTDat.rsd.SHT[r,i,],RTDat.rsd.SHT[r,j,],lag.max=1,plot = FALSE)$acf
  }
  return(Box.test(a/R,lag=1,type="Ljung-Box")$p.value)
}
cros1=apply(cbind(rep(1:100,times=100),rep(1:100,each=100)),1,crosfunc)
cros2=apply(cbind(rep(L^2-100:1+1,times=100),rep(L^2-100:1+1,each=100)),1,crosfunc)
cros3=apply(cbind(rep(L^2-100:1+1,times=100),rep(1:100,each=100)),1,crosfunc)
cros4=apply(cbind(rep(L^2-100:1+1,times=100),rep(2380+(-49):50,each=100)),1,crosfunc)
# t2=proc.time()[[3]]
# t2-t1=395.86
dataF=data.frame(loc1=c(rep((1:100),times=100),1.01,1.01),loc2=c(-rep(1:100,each=100),-1.01,-1.02),
                 z1=c(cros1,0,1),z2=c(cros2,0,1),z3=c(cros3,0,1),z4=c(cros4,0,1))
PT1=ggplot()+xlab(" ")+ylab(" ")+
  geom_point(aes(x=loc1,y=loc2,colour=z1),data=dataF[-c(10001,10002),],shape=15,size=0.5)+
  geom_point(aes(x=loc1,y=loc2,colour=z1),data=dataF[c(10001,10002),],alpha=0.1,size=0.01,shape=15)+
  scale_colour_gradient2(low = "#0072B2", mid = "white",high ="#E41A1C",
                         midpoint = 0.05,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(),
                   axis.title = element_blank(),
                   legend.title = element_blank(),
                   legend.position = "right",
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  scale_y_continuous(labels = c("100","75","50","25","1"))+
  scale_x_continuous(labels = c("1","25","50","75","100"))
PT2=ggplot()+xlab(" ")+ylab(" ")+
  geom_point(aes(x=loc1,y=loc2,colour=z2),data=dataF[-c(10001,10002),],shape=15,size=0.5)+
  geom_point(aes(x=loc1,y=loc2,colour=z2),data=dataF[c(10001,10002),],alpha=0.1,size=0.01,shape=15)+
  scale_colour_gradient2(low = "#0072B2", mid = "white",high ="#E41A1C",
                         midpoint = 0.05,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(),
                   axis.title = element_blank(),
                   legend.title = element_blank(),
                   legend.position = "right",
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  scale_y_continuous(labels = c("4761","4737","4712","4687","4662"))+
  scale_x_continuous(labels = c("4662","4687","4712","4737","4761"))
PT3=ggplot()+xlab(" ")+ylab(" ")+
  geom_point(aes(x=loc1,y=loc2,colour=z3),data=dataF[-c(10001,10002),],shape=15,size=0.5)+
  geom_point(aes(x=loc1,y=loc2,colour=z3),data=dataF[c(10001,10002),],alpha=0.1,size=0.01,shape=15)+
  scale_colour_gradient2(low = "#0072B2", mid = "white",high ="#E41A1C",
                         midpoint = 0.05,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(),
                   axis.title = element_blank(),
                   legend.title = element_blank(),
                   legend.position = "right",
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  scale_y_continuous(labels = c("4761","4737","4712","4687","4662"))+
  scale_x_continuous(labels = c("1","25","50","75","100"))
PT4=ggplot()+xlab(" ")+ylab(" ")+
  geom_point(aes(x=loc1,y=loc2,colour=z4),data=dataF[-c(10001,10002),],shape=15,size=0.5)+
  geom_point(aes(x=loc1,y=loc2,colour=z4),data=dataF[c(10001,10002),],alpha=0.1,size=0.01,shape=15)+
  scale_colour_gradient2(low = "#0072B2", mid = "white",high ="#E41A1C",
                         midpoint = 0.05,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(),
                   axis.title = element_blank(),
                   legend.title = element_blank(),
                   legend.position = "right",
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  scale_y_continuous(labels = c("4761","4737","4712","4687","4662"))+
  scale_x_continuous(labels = c("2331","2355","2380","2405","2430"))
ggarrange(PT1,PT3,PT4,PT2,ncol=4,nrow=1,legend="right",common.legend = TRUE)
# 10.29*2.24


######   Figure S5 in Section S4.1.1   #########################################
# Please refer to the reproduction of Figure 3 in Section 4.1                  #
################################################################################


######   Figure S6 in Section S3.4     #########################################
# Figure S6 illustrates the inference results of annual data obtained by HCBG. #
# All intermediate outputs are in sub-repository "Supplement/Huang/Annual/Outputs".#
################################################################################
### Load 7 ensembles of annual data and necessary information for HCBG
# t1=proc.time()[[3]]
R=7
YN=86
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
forcing=forcing$total[266:351]   # Year 2015--2100
id.land=read.csv("LENS2_Data/landid.csv")$x
landFlag=rep(0,nrow(Dat.loc))
landFlag[id.land]=1
landFlag=t(matrix(landFlag,288,192))
Dat=array(0,c(R,nrow(Dat.loc),YN))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Annual/dat_em1_year.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Annual/dat_em2_year.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Annual/dat_em3_year.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Annual/dat_em4_year.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Annual/dat_em5_year.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Annual/dat_em6_year.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Annual/dat_em7_year.mat")$Dat
np=import("numpy")
npdata=np$load('Supplement/Huang/results_model5_step2.npz')
psiLand=npdata$f[["psiLand"]]
alphaLand=npdata$f[["alphaLand"]]
nuLand=npdata$f[["nuLand"]]
psiOcean=npdata$f[["psiOcean"]]
alphaOcean=npdata$f[["alphaOcean"]]
nuOcean=npdata$f[["nuOcean"]]
npdata=np$load('Supplement/Huang/results_model5_step3.npz')
xi=npdata$f[["xi"]]
kappa=npdata$f[["kappa"]]
# t2=proc.time()[[3]]
# t2-t1=4.641

### Model the mean trend and temporal dependence at each grid point and plot Figures S6(a)
getK=function(phi){
  K=matrix(0,YN,YN)
  K[1,1]=1
  for(t in 2:YN){
    K[1,t]=phi*K[1,t-1]
  }
  for(i in 2:YN){
    K[i,i:YN]=K[i-1,(i-1):(YN-1)]
  }
  K=K+t(K)
  diag(K)=1
  return(K)
}
getX=function(rho){
  X=matrix(1,YN,3)
  X[,2]=forcing
  X[1,3]=0
  for(i in 2:YN){
    X[i,3]=(1-rho)*(t(rho^seq(i-2,0,by=-1))%*%forcing[1:(i-1)])
  }
  return(X)
}
profile_negllh=function(y,para.rhophi){
  rho=para.rhophi[1]
  phi=para.rhophi[2]
  K=getK(phi)
  X=getX(rho)
  YN=length(y)
  value=(YN-1)/YN*log(1-phi^2)
  Kiy=solve(K,y)
  part1=t(y)%*%Kiy
  part2=t(X)%*%Kiy
  part2=t(part2)%*%solve(t(X)%*%solve(K,X),part2)
  value=value+log(part1-part2)
  return(value)
}
hat.rhophi=function(i){ # estimate rho and phi
  obj=function(para.rhophi){
    value=0
    for(r in 1:R){value=value+profile_negllh(Dat[r,i,],para.rhophi)}
    return(value)
  }
  res=bobyqa(x0=c(0.9,0.9),fn=obj,lower=c(0.01,-0.99),upper=c(0.99,0.99))
  return(res$par)
}
cl<- makeCluster(4)
registerDoParallel(cl)
Res.hatrhophi= foreach(i=1:nrow(Dat.loc),
                       .combine=cbind,
                       .packages=c("nloptr")) %dopar% hat.rhophi(i)
stopCluster(cl)
write.csv(as.matrix(Res.hatrhophi),"Supplement/Huang/Annual/Outputs/Res_hatrhophi.csv")
# Tim.EhoPhi=rep(0,100)
# for(i in 1:100){
#   t1=proc.time()[[3]]
#   a=hat.rhophi(sample(1:nrow(Dat.loc),1))
#   t2=proc.time()[[3]]
#   Tim.EhoPhi[i]=t2-t1
# }
# mean(Tim.EhoPhi)=2.04137
Res.hatrhophi=as.matrix(read.csv("Supplement/Huang/Annual/Outputs/Res_hatrhophi.csv")[,-1])
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z=Res.hatrhophi[2,])
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradientn(colours = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT

### Evaluate mean trend and sd at each grid point and plot Figure S6(b)
# t1=proc.time()[[3]]
Eta=array(0,c(R,nrow(Dat.loc),YN))
Dat.hat=matrix(0,nrow(Dat.loc),YN)
Sig=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  K=getK(Res.hatrhophi[2,i])
  X=getX(Res.hatrhophi[1,i])
  Sigma=0
  Beta=rep(0,3)
  for(r in 1:R){
    beta=solve(t(X)%*%solve(K,X))%*%t(X)%*%solve(K,Dat[r,i,])
    Beta=Beta+beta
    err=Dat[r,i,]-X%*%beta
    Sigma=Sigma+t(solve(K,err))%*%err/YN
  }
  Beta=Beta/R
  Sigma=sqrt(Sigma/R)
  Sig[i]=Sigma
  Dat.hat[i,]=c(X%*%Beta)
  epsilon=Dat[,i,]-rep(1,R)%*%t(Dat.hat[i,])
  Eta[,i,1]=epsilon[,1]/c(Sigma)
  Eta[,i,2:YN]=(epsilon[,2:YN]-Res.hatrhophi[2,i]*epsilon[,1:(YN-1)])/c(Sigma)
}
# t2=proc.time()[[3]]
# t2-t1=723.85
write.csv(Sig,"Supplement/Huang/Annual/Outputs/Sig_Annual.csv")
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z=Sig)
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  scale_colour_gradientn(colours = Col)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))
PT

### Model longitude dependence at each latitude and plot Figure S6(c)
omega=2*pi/288*(0:287)
negllh_long=function(j,arg){
  psiL=arg[1]
  alphaL=arg[2]
  nuL=arg[3]
  psiO=arg[4]
  alphaO=arg[5]
  nuO=arg[6]
  
  tmpL=sqrt(psiL/((alphaL*alphaL+4*sin(omega/2)*sin(omega/2))^(nuL+0.5)))
  tmpO=sqrt(psiO/((alphaO*alphaO+4*sin(omega/2)*sin(omega/2))^(nuO+0.5)))
  flag=landFlag[j,]
  SigmaFactor=(outer(flag,tmpL)+outer(1-flag,tmpO))*exp(1i*outer(0:287,omega))/sqrt(288)
  eigres=svd(SigmaFactor)
  eig=eigres$d
  
  ldet=2*sum(log(abs(eig)))
  value=ldet*YN*R
  Y=Eta[,(288*(j-1)+1):(288*j),]  # R*288*YN
  for(r in 1:R){
    value=value+sum(apply((t(Y[r,,])%*%solve(SigmaFactor))^2,1,sum))
  }
  return(Re(value))
}
Res.longLO.ini=as.matrix(read.csv("Supplement/Huang/Annual/Outputs/Res_longOL_ini.csv")[,-1])
hat.longLO=function(j){
  paraLand=Res.longLO.ini[1:3,j]
  paraOcean=Res.longLO.ini[4:6,j]
  obj.ocean=function(argOcean){return(negllh_long(j,c(paraLand,argOcean)))}
  res.ocean=bobyqa(x0=paraOcean,fn=obj.ocean,lower=0.75*paraOcean,upper=1.25*paraOcean)$par
  obj.land=function(argLand){return(negllh_long(j,c(argLand,res.ocean)))}
  res.land=bobyqa(x0=paraLand,fn=obj.land,lower=0.75*paraLand,upper=1.25*paraLand)$par
  negllh=negllh_long(j,c(res.land,res.ocean))
  negllh.ini=negllh_long(j,Res.longLO.ini[,j])
  if(negllh<negllh.ini){
    return(c(res.land,res.ocean))
  }
  else{
    return(Res.longLO.ini[,j])
  }
}
cl<- makeCluster(4)
registerDoParallel(cl)
Res.longLO= foreach(j=1:192,
                    .combine=cbind,
                    .packages=c("nloptr")) %dopar% hat.longLO(j)
stopCluster(cl)
write.csv(as.matrix(Res.longLO),"Supplement/Huang/Annual/Outputs/Res_longLO.csv")
# t1=proc.time()[[3]]
# hat.longLO(sample(1:192,1))
# t2=proc.time()[[3]]
# t2-t1=10.207
Res.longLO=as.matrix(read.csv("Supplement/Huang/Annual/Outputs/Res_longLO.csv")[,-1])
dataF=data.frame(psi=c(Res.longLO[1,],Res.longLO[4,]),
                 alpha=c(Res.longLO[2,],Res.longLO[5,]),
                 nu=c(Res.longLO[3,],Res.longLO[6,]),
                 Type=as.factor(c(rep("Land",times=192),rep("Ocean",times=192))),
                 Index=rep(as.numeric(levels(as.factor(Dat.loc[,2]))),times=2))
PT1=ggplot()+geom_line(aes(x=Index,y=psi,color=Type),dataF)+
  scale_color_manual(values = c(cbPalette[4],cbPalette[6]))+
  xlab("Lattitude")+ylab(expression(hat(psi)))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   legend.justification = c(1,1),
                   legend.position = c(0.99,0.99),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(1,"line"))
PT2=ggplot()+geom_line(aes(x=Index,y=alpha,color=Type),dataF)+
  scale_color_manual(values = c(cbPalette[4],cbPalette[6]))+
  xlab("Lattitude")+ylab(expression(hat(alpha)))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   legend.justification = c(1,1),
                   legend.position ="none",
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(1,"line"))
PT3=ggplot()+geom_line(aes(x=Index,y=nu,color=Type),dataF)+
  scale_color_manual(values = c(cbPalette[4],cbPalette[6]))+
  xlab("Lattitude")+ylab(expression(hat(nu)))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   legend.justification = c(1,1),
                   legend.position ="none",
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_blank(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(1,"line"))
PT1+PT2+PT3



######   Figures S7 and S8 in Section S4.1.3   #################################
# Please refer to the reproduction of Figure 5 in Section 4.1                  #
################################################################################


######   Figure S9 in Section S4.2     #########################################
# Figure S9 illustrates inference results for deterministic components of the  #
# monthly data. Note that Figures S9-S13 are all for the monthly data, so      #
# please do not clear the environment after reproducing Figure S9 if you want  #
# to reproduce Figures S10-S13.                                                #
# All intermediate outputs are in sub-repository "Monthly/Outputs".            #
################################################################################
### Load 7 ensembles of monthly data and necessary information
# t1=proc.time()[[3]]
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
lseq=mseq=0
for(i in 1:143){
  lseq=c(lseq,rep(i,2*i+1))
  mseq=c(mseq,seq(-i,i,by=1))
}
# t2=proc.time()[[3]]
# t2-t1=66.098

### Model the deterministic components, evaluate mean trend and sd for each grid point, and plot Figure S9
getX=function(rho){
  M=3    # Please refer lines 77-150 of "Monthly/Monthly.R" to see how to choose M=3 by using BIC, provided with intermediate results "BIC_KM.csv"
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
profile_negllh=function(y,rho){
  part1=t(y)%*%y
  X=getX(rho)
  part2=t(X)%*%y
  part2=t(part2)%*%solve(t(X)%*%X,part2)
  value=log(part1-part2)
  return(value)
}
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
write.csv(c(Res.hatrho),"Monthly/Outputs/Res_Hatrho.csv") 
# t1=proc.time()[[3]]
# hat.rho(sample(1:nrow(Dat.loc),1))
# t2=proc.time()[[3]]
# t2-t1=13.0 
Res.hatrho=read.csv("Monthly/Outputs/Res_Hatrho.csv")$x
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
writeMat("Monthly/Outputs/BetaAB.mat",BetaAB=as.matrix(Res.BetaSig[1:9,]))
Sig=sqrt(c(Res.BetaSig[10,]))
write.csv(Sig,"Monthly/Outputs/Monthly_Sig.csv") 
Dat.hat=t(as.matrix(Res.BetaSig[-(1:10),]))
# t1=proc.time()[[3]]
# BetaSighat(sample(1:nrow(Dat.loc),1))
# t2=proc.time()[[3]]
# t2-t1=0.037
Dat.mean=(Dat[1,,105]+Dat[2,,105]+Dat[3,,105]+Dat[4,,105]+Dat[5,,105]+Dat[6,,105]+Dat[7,,105])/7
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
Sig=read.csv("Monthly/Outputs/Monthly_Sig.csv")$x
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

######   Figure S10 in Section S4.2     ########################################
# Figure S10 illustrates inference results for stochastic components of the    #
# monthly data. Assume that we keep all the intermediate results of Figure S9. #
# All intermediate outputs are in sub-repository "Monthly/Outputs".            #
################################################################################
### Calculate stochastic components Z_t^{(r)}(L_i,l_j)  by detrending and rescaling 
t1=proc.time()[[3]]
Dat.srsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.srsd[r,,]=(Dat[r,,]-Dat.hat)/Sig
}
t2=proc.time()[[3]]
# t2-t1=19.792

### Do SHT with Q=144 for the stochastic component at each ensemble and time point 
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
# As we illustrated before, it takes about 4.2 seconds for each ensemble r and each time point t.
# Without parallel, it will take 4.2*7*(86*12)=30340.8 seconds

### Calculate BIC values under different Q values and plot Figure S10(a)
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
# The major computation time is for calculating the inverse SHT. But this step takes time because it calculates 
# the inverse SHT for each ensemble r, each time point t, and all candidates of Q (represented by llseq, llseq1, and llseq2).
# Take Q=90 as an example, which maximizes the computational time of inversing SHT, it takes about
# 0.9 seconds for each ensemble r and time point t. If we run this step without doing parallel,
# we would take at most 0.9*7*(86*12)*(8+11+11)=195048 seconds.
BIC.land=as.matrix(read.csv("Monthly/Outputs/BIC_land.csv")[,-1])
BIC.ocean=as.matrix(read.csv("Monthly/Outputs/BIC_ocean.csv")[,-1])
BICd.land=as.matrix(read.csv("Monthly/Outputs/BICd_land.csv")[,-1])
BICd.ocean=as.matrix(read.csv("Mothly/Outputs/BICd_ocean.csv")[,-1])
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

### Calculate v^2(L_i,l_j) under Ql=36 and Qo=70 and plot Figure S10(b)
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
# As before, the major computational time in this step is to calculate the inverse SHT with Q=36 and 70 for 
# each ensemble r, time point t. It will take about 0.8*7*86*12=5779.2 seconds without the parallel.
v2hat=read.csv("Monthly/Outputs/v2hat.csv")$x
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

### Do the real-valued transformation to SHT coefficients so that they are real values
# t1=proc.time()[[3]]
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
# t2=proc.time()[[3]]
# t2-t1=2.847

### Test the normality of coefficients
# t1=proc.time()[[3]]
cBera=rep(0,L^2)
for(i in 1:L^2){
  cBera[i]=jarque.bera.test(c(TDat.rsd.SHT[,i,]))$p.value
}
idBera=which(cBera<=0.05)  # non-normal
length(idBera)/L/L   # 0.4263265
# t2=proc.time()[[3]]
# t2-t1=6.734

### Model the temporal dependence structure using a Tukey g-and-h autoregressive model with order P=1 and plot Figure S10(c)
# bicpfunc1=function(j){ # choose order p using BIC
#   bicp=rep(0,5)
#   y=TDat.rsd.SHT[,j,]
#   for(p in 1:5){
#     skey=skewness(c(y))
#     y1=y*sign(skey)
#     TukeyPara=function(para){
#       b=para[1]
#       g=para[2]
#       h=para[3]
#       phi=para[-(1:3)]
#       taugh=function(z){
#         if(g==0){
#           return(b*z*exp(h*z^2/2))
#         }
#         else{return(b*(exp(g*z)-1)/g*exp(h*z^2/2))}
#       }
#       K=max(1000,ncol(y))
#       Z=seq(-10,10,length.out=K)
#       Y=taugh(Z)
#       kid=matrix(.bincode(y1,Y,right=FALSE),dim(y))
#       tildez=Phi=matrix(0,nrow(y),ncol(y))
#       u2=rep(0,nrow(y))
#       for(r in 1:nrow(y)){
#         for(i in 1:ncol(y)){
#           tildez[r,i]=Z[kid[r,i]]+(y1[r,i]-Y[kid[r,i]])/(Y[kid[r,i]+1]-Y[kid[r,i]])*(Z[2]-Z[1])
#           Phi[r,i]=-h/2*tildez[r,i]^2-log(exp(g*tildez[r,i])+(exp(g*tildez[r,i])-1)/g*h*tildez[r,i])-log(b)
#         }
#         xx=rep(0,ncol(y)-p)
#         for(tt in p:1){
#           xx=cbind(xx,tildez[r,tt:(tt+ncol(y)-p-1)])
#         }
#         xx=as.matrix(xx[,-1])
#         u2[r]=mean((tildez[r,-(1:p)]-xx%*%phi)^2)
#       }
#       u2=mean(u2)
#       res=sum(Phi)-nrow(y)*(ncol(y)-p)/2*log(2*pi)-nrow(y)*(ncol(y)-p)/2*log(u2)
#       return(-res)
#     }
#     res=optim(c(0.5,0.2,0.1,rep(0.8,p)),TukeyPara,lower=c(0.1,0.05,0.05,rep(0.1,p)),method="L-BFGS-B",upper=c(4,0.5,0.5,rep(1,p)))
#     Paraest=res$par
#     Paraest[2]=sign(skey)*Paraest[2]
#     
#     taugh=function(z){return(Paraest[1]*(exp(Paraest[2]*z)-1)/Paraest[2]*exp(Paraest[3]*z^2/2))}
#     
#     K=max(1000,ncol(y))
#     Z=seq(-10,10,length.out=K)
#     Y=taugh(Z)
#     kid=matrix(.bincode(y,Y,right=FALSE),dim(y))
#     
#     tildez=Phi=matrix(0,nrow(y),ncol(y))
#     u2=rep(0,nrow(y))
#     for(r in 1:nrow(y)){
#       for(i in 1:ncol(y)){
#         tildez[r,i]=Z[kid[r,i]]+(y[r,i]-Y[kid[r,i]])/(Y[kid[r,i]+1]-Y[kid[r,i]])*(Z[2]-Z[1])
#         Phi[r,i]=-Paraest[3]/2*tildez[r,i]^2-log(exp(Paraest[2]*tildez[r,i])+
#                                                    (exp(Paraest[2]*tildez[r,i])-1)/Paraest[2]*Paraest[3]*tildez[r,i])-log(Paraest[1])
#       }
#       xx=rep(0,ncol(y)-p)
#       for(tt in p:1){
#         xx=cbind(xx,tildez[r,tt:(tt+ncol(y)-p-1)])
#       }
#       xx=as.matrix(xx[,-1])
#       u2[r]=mean((tildez[r,-(1:p)]-xx%*%Paraest[-(1:3)])^2)
#     }
#     u2=mean(u2)
#     bicp[p]=p*log((ncol(y)-p)*nrow(y))+nrow(y)*(ncol(y)-p)*log(2*pi)+nrow(y)*(ncol(y)-p)*log(u2)
#   }
#   return(bicp)
# }
# cl<- makeCluster(4) 
# registerDoParallel(cl) 
# bicp1=foreach(j=idBera,
#               .combine=rbind,
#               .packages=c("moments","nloptr")) %dopar% bicpfunc1(j)
# stopCluster(cl)
# write.csv(bicp1,"Monthly/Outputs/bicp_Tukey.csv")
# bicp=as.matrix(read.csv("Monthly/Outputs/bicp_noTukey.csv")[,-1]) # bicp is an intermediate result
# # for SG without using Tukey g-and-h. Please take a look at lines 436-570 of "Monthly/Monthly.R" for more details.
# bicp1=as.matrix(read.csv("Monthly/Outputs/bicp_Tukey.csv")[,-1])
# bicp[idBera,]=as.matrix(bicp1)
# length(which(apply(bicp,1,which.min)==1))/L/L   # 0.8728571
# length(which(apply(bicp,1,which.min)==2))/L/L   # 0.08816327
# length(which(apply(bicp,1,which.min)==3))/L/L   # 0.03408163
# length(which(apply(bicp,1,which.min)==4))/L/L   # 0.03801722
# length(which(apply(bicp,1,which.min)==5))/L/L   # 0.001836735
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
  # t1=proc.time()[[3]]
  Tukeyres[j,]=TukeyAuto(TDat.rsd.SHT[,j,])
  # t2=proc.time()[[3]]
  # t2-t1=2.823
}  # 2.823*2089=5897.247
write.csv(Tukeyres,"Monthly/Outputs/Tukeyres.csv")
Tukeyres=as.matrix(read.csv("Monthly/Tukeyres.csv")[,-1])
# t1=proc.time()[[3]]
idBera.new=which(Tukeyres[,1]!=0)    # need Tukey 
idphi1=which(Tukeyres[,1]==0)     #    no need Tukey   
idphi2=intersect(which(Tukeyres[,1]!=0),which(Tukeyres[,4]==0))
rtc=rep(1,L^2)  
CDat.rsd.SHT=TDat.rsd.SHT
for(j in idBera.new){ # Tukey g-and-h transformation for idBera.new
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
# t2=proc.time()[[3]]
# t2-t1=10.326
# t1=proc.time()[[3]]
TPhi.hat=rep(0,L^2) # calculate \phi for idphi1 and idphi2
for(i in c(idphi1,idphi2)){
  yy=c(t(CDat.rsd.SHT[,i,2:YN]))
  yx=c(t(CDat.rsd.SHT[,i,-YN]))
  TPhi.hat[i]=solve(t(yx)%*%yx)%*%t(yx)%*%yy
}
TPhi.hat[-c(idphi1,idphi2)]=Tukeyres[-c(idphi1,idphi2),4]
write.csv(TPhi.hat,"Monthly/Outputs/Phihat_Tukey.csv")
# t2=proc.time()[[3]]
# t2-t1=2.329
TPhi.hat=read.csv("Monthly/Outputs/Phihat_Tukey.csv")$x
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


######   Figures S11 and S12 in Section S4.2     ###############################
# Figures S11 and S12 illustrates the performance of the generated monthly     #
# emulations. Assume that we keep all the intermediate results of Figures S9 and S10.#
# All intermediate outputs are in sub-repository "Monthly/Outputs".            #
################################################################################
### Model the spatial dependence by evaluating the covariance matrix of (read-valued and Gaussianized) SHT coefficients
# t1=proc.time()[[3]]
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
# t2=proc.time()[[3]]
# t2-t1=48.262

### Calculate the covaraince matrix \check U 
# t1=proc.time()[[3]]
CU.axial=CK.axial-outer(TPhi.hat,TPhi.hat,"*")*CK.axial
# t2=proc.time()[[3]]
# t2-t1=0.283
writeMat("Monthly/Outputs/U.mat",U=CU.axial)      # This file is compressed to reduce size

### Generate R'=7 ensembles of monthly emulations using 4 cores
# t1=proc.time()[[3]]
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
# t2=proc.time()[[3]]
# t2-t1=85.784
Gen.Dat.Y=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  # t1=proc.time()[[3]]
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
  # t2=proc.time()[[3]]
  # t2-t1=410.501
}
# 410.501*7=2873.507

### Calculate I.uq values and plot Figures S11(a), S11(b), and S11(e)
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
# Tim=0
# for(i in 1:100){
#   t1=proc.time()[[3]]
#   findIuq(sample(1:nrow(Dat.loc),1))
#   t2=proc.time()[[3]]
#   Tim=Tim+t2-t1
# }
# Tim/100=0.08984
# 0.08984*192*288/4=1241.948
write.csv(c(I.uq.Tukey),"Monthly/Outputs/Iuq_Tukey.csv")
I.uq=read.csv("Monthly/Outputs/Iuq_Tukey.csv")$x
I.uq.Huang=read.csv("Supplement/Huang/Monthly/Outputs/Iuq_Huang.csv")$x  # Load the Iuq values of monthly emulations generated by HCBG, 
# details about them can be found in "Supplement/Huang/Monthly/Huang_Monthly.R"
I.uq.noTukey=read.csv("Monthly/Outputs/Iuq_noTukey.csv")$x               # Load the Iuq values of monthly emulations generated by SHT without TGH, 
# details about them can be found in lines 436-570 of "Monthly/Monthly.R"
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
dataF=data.frame(Iuq=c(I.uq,I.uq.noTukey),
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

### Calculate WD_S values and plot Figures S11(c), S11(d), and S11(f)
WD.time.Tukey=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  # t1=proc.time()[[3]]
  WD.time.Tukey[i]=wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat.Y[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
  # t2=proc.time()[[3]]
  # t2-t1=0.019
}
# 0.019*192*288/4=262.656
write.csv(WD.time.Tukey,"Monthly/Outputs/WD_time_Tukey.csv")
# WD.space.Tukey=rep(0,YN)
# for(t in 1:YN){
#   WD.space.Tukey[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat.Y[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
# }
# write.csv(WD.space.Tukey,"Monthly/Outputs/WD_space_Tukey.csv")
WD.time=read.csv("Monthly/Outputs/WD_time_Tukey.csv")$x
WD.Huang=read.csv("Supplement/Huang/Monthly/Outputs/WD_Huang.csv")$x  # Load the WD values of monthly emulations generated by HCBG, 
# details about them can be found in "Supplement/Huang/Monthly/Huang_Monthly.R"
WD.noTukey=read.csv("Monthly/Outputs/WD_time_noTukey.csv")$x          # Load the WD values of monthly emulations generated by SHT without TGH, 
# details about them can be found in lines 436-570 of "Monthly/Monthly.R"
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
dataF=data.frame(WD=c(WD.time,WD.noTukey),
                 Type=as.factor(rep(c("With TGH","Without TGH"),each=55296)))
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

### Plot Figures S12(a) and S12(b)
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

### Plot Figures S12(c) and S12(d)
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


######   Figure S13 in Section S4.2     ########################################
# Figure S13 illustrates the performance of the generated monthly emulations by#
# aggregating them to be annual emulations. Assume that we keep all the inter- #
# mediate results of Figures S9-S12, especially the generated monthly emulations.#
# All intermediate outputs are in sub-repository "Monthly/Aggregate".          #
################################################################################
### Aggregate monthly emulations to annual emulations
# t1=proc.time()[[3]]
Gen.Dat.year=array(0,c(R,nrow(Dat.loc),86))
for(r in 1:R){
  for(i in 1:86){
    Gen.Dat.year[r,,i]=apply(Gen.Dat.Y[r,,(i-1)*12+(1:12)],1,mean)
  }
}
# t2=proc.time()[[3]]
# t2-t1=429.646

### Load the annual simulations
# t1=proc.time()[[3]]
Dat=array(0,c(R,nrow(Dat.loc),86))
Dat[1,,]=readMat("LENS2_Data/Annual/dat_em1_year.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Annual/dat_em2_year.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Annual/dat_em3_year.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Annual/dat_em4_year.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Annual/dat_em5_year.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Annual/dat_em6_year.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Annual/dat_em7_year.mat")$Dat
# t2=proc.time()[[3]]
# t2-t1=5.196

### Calculate I.uq values and plot Figures S13(a) and S13(b)
# t1=proc.time()[[3]]
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
  return(get.IUQ(Gen.Dat.year[,i,],Dat[1:R,i,]))
}
cl=makeCluster(4)
registerDoParallel(cl)
I.uq.montoan=foreach(i=1:nrow(Dat.loc),
                     .combine = cbind,
                     .packages = c("fdaoutlier")) %dopar% findIuq(i)
stopCluster(cl)
# t2=proc.time()[[3]]
# t2-t1=139.211
write.csv(c(I.uq.montoan),"Monthly/Aggregate/Iuq_Monthly_to_Annual.csv")
I.uq=read.csv("Annual/Outputs/Iuq_axialnon.csv")$x
I.uq.montoan=read.csv("Monthly/Aggregate/Iuq_Monthly_to_Annual.csv")$x
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z1=c(I.uq),z2=c(I.uq.montoan))
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z2),data=dataF,shape=15)+
  geom_point(aes(x=lon,y=lat),data=dataF[c(45461,28236),],shape=4)+
  scale_colour_gradient2(low ="#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 1,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",
                         limit=c(min(dataF$z2),max(dataF$z1)))+
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
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z1),data=dataF,shape=15)+
  geom_point(aes(x=lon,y=lat),data=dataF[c(45461,28236),],shape=4)+
  scale_colour_gradient2(low ="#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 1,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",
                         limit=c(min(dataF$z2),max(dataF$z1)))+
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

### Calculate WD values and plot Figures S13(b) and S13(c)
# t1=proc.time()[[3]]
WD.montoan=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  WD.montoan[i]=wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat.year[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
# t2=proc.time()[[3]]
# t2-t1=42.303
write.csv(c(WD.montoan),"Monthly/Aggregate/WD_Monthly_to_Annual.csv")
WD.time=read.csv("Annual/Outputs/WD_time.csv")$x
WD.montoan=read.csv("Monthly/Aggregate/WD_Monthly_to_Annual.csv")$x
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z1=c(WD.time),z2=c(WD.montoan))
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z2),data=dataF,shape=15)+
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
  geom_point(aes(x=lon,y=lat,colour=z1),data=dataF,shape=15)+
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

### Calculate I.fit values and plot Figure S13(e)
# t1=proc.time()[[3]]
Gen.Dat.mean=apply(Gen.Dat.year,c(2,3),mean)
Dat.mean=apply(Dat,c(2,3),mean)
I.fit.up=I.fit.down=rep(0,nrow(Dat.loc))
for(r in 1:R){
  I.fit.up=I.fit.up+apply((Dat[r,,]-Gen.Dat.mean)^2,1,sum)
  I.fit.down=I.fit.down+apply((Dat[r,,]-Dat.mean)^2,1,sum)
}
I.fit.montoan=I.fit.up/I.fit.down*(R-1)/R
# t2=proc.time()[[3]]
# t2-t1=129.05
write.csv(I.fit.montoan,"Monthly/Aggregate/Ifit_Monthly_to_Annual.csv")
I.fit=as.matrix(read.csv("Annual/Outputs/IfitwithRs.csv")[,-1])
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),z1=c(I.fit[,R-1]),z2=I.fit.montoan)
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z1),data=dataF,shape=15)+
  geom_point(aes(x=lon,y=lat),data=dataF[c(45461,28236),],shape=4)+
  scale_colour_gradient2(low = "#0072B2",mid = "white",high = "#E41A1C",
                         midpoint = 1,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",limit=c(min(I.fit),max(I.fit.montoan)))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color=expression(I[fit]))
PT

### Plot Figure S13(f) using "BIC_KM.csv" in "Monthly/Outputs"





