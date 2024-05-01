################  Wrapper file for sequentially implementing each figure in the main manuscript  ################
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
# Please refer to the README.md file for more details                                      #
############################################################################################

######   Figure 1 in Section 2 (and Figure S1 in Section S2)  ##################
# Figures 1 and S1 are used to illustrate the characteristics of surface       #
# temperature simulations in different scales.                                 #
# All the intermediate outputs can be found in sub-repository "Figure1".       #
################################################################################
### Load the annual data and necessary information
# t1=proc.time()[[3]]
R=7
YN=86
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
forcing=forcing$total[266:351] 
Dat=array(0,c(R,nrow(Dat.loc),YN))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Annual/dat_em1_year.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Annual/dat_em2_year.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Annual/dat_em3_year.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Annual/dat_em4_year.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Annual/dat_em5_year.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Annual/dat_em6_year.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Annual/dat_em7_year.mat")$Dat
# t2=proc.time()[[3]]
# t2-t1=5.182

### Calculate annual ensemble mean to plot Figure 1(a) 
# t1=proc.time()[[3]]
Dat.mean=apply(Dat[,,2023-2014],2,mean) 
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degree=Dat.mean)
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
# t2=proc.time()[[3]]
# t2-t1=1.939

### Calculate annual ensemble sd to plot Figure 1(b) 
# t1=proc.time()[[3]]
Dat.sd=apply(Dat[,,2023-2014],2,sd) 
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),Degree=Dat.sd)
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
# t2=proc.time()[[3]]
# t2-t1=3.771

### Plot Figure 1(c)
# t1=proc.time()[[3]]
dataF=data.frame(Time=rep(2015:2100,times=R),
                 TempL=c(t(Dat[,39780,])),
                 TempO=c(t(Dat[,18281,])),
                 LEN=as.factor(rep(c("1","2","3","4","5","6","7"),each=YN)))
dataFf=data.frame(Time=2015:2100,Forc=forcing+13,Type=as.factor(rep("Forcing",YN)))
dataL=data.frame(Time=rep(2015:2100,times=3),
                 Temp=c(apply(Dat[,39780,],2,mean),apply(Dat[,18281,],2,mean),forcing+13),
                 Type=as.factor(rep(c("Mean on GL","Mean on GO","Forcing"),each=YN)))
PT=ggplot()+
  geom_line(aes(x=Time,y=TempL,colour=LEN),alpha=0.5,data=dataF)+
  geom_line(aes(x=Time,y=TempO,colour=LEN),alpha=0.5,data=dataF)+
  xlab("Year")+ylab("Temperature(\u00b0C)")+
  geom_line(aes(x=Time,y=Temp,linetype=Type,size=Type),data=dataL)+
  scale_linetype_manual(values=c(4,2,1))+
  scale_size_manual(values=c(0.8,0.6,0.6))+
  scale_y_continuous(limits = c(0,28),
                     sec.axis = sec_axis(~.-13,name=expression(RF(Wm^{-2}))))+
  scale_color_manual(values = c(cbPalette[1:7]),guide="none")+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.justification = c(0.99,0.01),
                   legend.position =c(0.99,0.01),
                   # legend.position = "right",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1.2,"line"))
PT+annotate("text",x=2100,y=24,label="GO")+annotate("text",x=2100,y=11.7,label="GL")
# t2=proc.time()[[3]]
# t2-t1=0.615

### Save R annual time series at two grid points for ploting Figure 1(e)
# t1=proc.time()[[3]]
DatA1=Dat[,39780,]
writeMat("Figure1/DatA1.mat",Dat=DatA1)
DatA2=Dat[,18281,]
writeMat("Figure1/DatA2.mat",Dat=DatA2)
# t2=proc.time()[[3]]
# t2-t1=0.402

### Calculate skewness and kurtosis for residuals of annaul data to plot Figures S1(a) and S1(d)
# t1=proc.time()[[3]]
Dat.mean=(Dat[1,,]+Dat[2,,]+Dat[3,,]+Dat[4,,]+Dat[5,,]+Dat[6,,]+Dat[7,,])/7
Skewfunc=function(i){return(skewness(c(Dat[,i,]-rep(1,R)%*%t(Dat.mean[i,]))))}
Kurtfunc=function(i){return(kurtosis(c(Dat[,i,]-rep(1,R)%*%t(Dat.mean[i,]))))}
SkewA=sapply(1:nrow(Dat.loc),Skewfunc)
write.csv(SkewA,"Figure1/SkewA.csv")
KurtA=sapply(1:nrow(Dat.loc),Kurtfunc)
write.csv(KurtA,"Figure1/KurtA.csv")
# t2=proc.time()[[3]]
# t2-t1=13.882

### Load the monthly data 
# t1=proc.time()[[3]]
YN=86*12
Dat=array(0,c(R,nrow(Dat.loc),YN))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Monthly/dat_em1_month.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Monthly/dat_em2_month.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Monthly/dat_em3_month.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Monthly/dat_em4_month.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Monthly/dat_em5_month.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Monthly/dat_em6_month.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Monthly/dat_em7_month.mat")$Dat
# t2=proc.time()[[3]]
# t2-t1=71.644

### Plot Figure 1(d) 
# t1=proc.time()[[3]]
Mbe=((2021-2015)*12+1):((2025-2014)*12) # 2021--2025
dataF=data.frame(Time=rep(Mbe,times=7),
                 TempL=c(t(Dat[,39780,Mbe])),
                 TempO=c(t(Dat[,18281,Mbe])),
                 LEN=as.factor(rep(c("1","2","3","4","5","6","7"),each=60)))
dataL=data.frame(Time=rep(Mbe,times=2),
                 Temp=c(apply(Dat[,39780,Mbe],2,mean),apply(Dat[,18281,Mbe],2,mean)),
                 Type=as.factor(rep(c("Mean on GL","Mean on GO"),each=60)))
PT=ggplot()+
  xlab("Year")+ylab("Temperature(\u00b0C)")+
  geom_line(aes(x=Time,y=TempL,colour=LEN),alpha=0.5,data=dataF)+
  geom_line(aes(x=Time,y=TempO,colour=LEN),alpha=0.5,data=dataF)+
  scale_color_manual(values = c(cbPalette[1:7]),guide="none")+
  scale_x_continuous(breaks=c(79,91,103,115,127),labels=2021:2025)+
  geom_line(aes(x=Time,y=Temp,linetype=Type),data=dataL)+
  scale_linetype_manual(values=c(2,1))+
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
PT+annotate("text",x=134,y=2,label="GL")+annotate("text",x=134,y=20,label="GO")
# t2=proc.time()[[3]]
# t2-t1=0.366

### Save R=7 monthly time series at two grid points for ploting Figure 1(f) 
# t1=proc.time()[[3]]
DatM1=Dat[,39780,]
writeMat("Figure1/DatM1.mat",Dat=DatM1)
DatM2=Dat[,18281,]
writeMat("Figure1/DatM2.mat",Dat=DatM2)
# t2=proc.time()[[3]]
# t2-t1=0.492

### Calculate skewness and kurtosis for residuals of monthly data to plot Figures S1(b) and S1(e)
# t1=proc.time()[[3]]
Dat.mean=(Dat[1,,]+Dat[2,,]+Dat[3,,]+Dat[4,,]+Dat[5,,]+Dat[6,,]+Dat[7,,])/7
Skewfunc=function(i){return(skewness(c(Dat[,i,]-rep(1,R)%*%t(Dat.mean[i,]))))}
Kurtfunc=function(i){return(kurtosis(c(Dat[,i,]-rep(1,R)%*%t(Dat.mean[i,]))))}
SkewM=sapply(1:nrow(Dat.loc),Skewfunc)
write.csv(SkewM,"Figure1/SkewM.csv")
KurtM=sapply(1:nrow(Dat.loc),Kurtfunc)
write.csv(KurtM,"Figure1/KurtM.csv")
# t2=proc.time()[[3]]
# t2-t1=100.835

### Load the daily data
# t1=proc.time()[[3]]
YN=365*5
Dat=array(0,c(R,nrow(Dat.loc),YN))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Daily/dat_em1_day.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Daily/dat_em2_day.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Daily/dat_em3_day.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Daily/dat_em4_day.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Daily/dat_em5_day.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Daily/dat_em6_day.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Daily/dat_em7_day.mat")$Dat
# t2=proc.time()[[3]]
# t2-t1=104.142

### Save R=7 daily time series at two grid points for ploting Figure 1(g)
# t1=proc.time()[[3]]
DatD1=Dat[,39780,]
writeMat("Figure1/DatD1.mat",Dat=DatD1)
DatD2=Dat[,18281,]
writeMat("Figure1/DatD2.mat",Dat=DatD2)
# t2=proc.time()[[3]]
# t2-t1=0.37

### Calculate skewness and kurtosis for residuals of daily data to plot Figures S1(c) and S1(f)
# t1=proc.time()[[3]]
Dat.mean=(Dat[1,,]+Dat[2,,]+Dat[3,,]+Dat[4,,]+Dat[5,,]+Dat[6,,]+Dat[7,,])/7
Skewfunc=function(i){return(skewness(c(Dat[,i,]-rep(1,R)%*%t(Dat.mean[i,]))))}
Kurtfunc=function(i){return(kurtosis(c(Dat[,i,]-rep(1,R)%*%t(Dat.mean[i,]))))}
SkewD=sapply(1:nrow(Dat.loc),Skewfunc)
write.csv(SkewD,"Figure1/SkewD.csv")
KurtD=sapply(1:nrow(Dat.loc),Kurtfunc)
write.csv(KurtD,"Figure1/KurtD.csv")
# t2=proc.time()[[3]]
# t2-t1=180.749

### Plot Figures 1(e)--1(g)
# t1=proc.time()[[3]]
DatA1.mean=apply(DatA1,2,mean)
DatA2.mean=apply(DatA2,2,mean)
DatM1.mean=apply(DatM1,2,mean)
DatM2.mean=apply(DatM2,2,mean)
DatD1.mean=apply(DatD1,2,mean)
DatD2.mean=apply(DatD2,2,mean)
dataF1=data.frame(Temp=c(c(DatA1-rep(1,R)%*%t(DatA1.mean)),c(DatM1-rep(1,R)%*%t(DatM1.mean)),c(DatD1-rep(1,R)%*%t(DatD1.mean))),
                  Type=as.factor(c(rep("Annual",86*R),rep("Monthly",86*12*R),rep("Daily",1825*R))))
dataF2=data.frame(Temp=c(c(DatA2-rep(1,R)%*%t(DatA2.mean)),c(DatM2-rep(1,R)%*%t(DatM2.mean)),c(DatD2-rep(1,R)%*%t(DatD2.mean))),
                  Type=as.factor(c(rep("Annual",86*R),rep("Monthly",86*12*R),rep("Daily",1825*R))))
PT1=ggplot(dataF1[which(dataF1$Type=="Annual"),],aes(x=Temp))+
  geom_histogram(aes(y=after_stat(density)),
                 bins=10,fill=cbPalette[4],alpha=0.8,color="black")+
  geom_density(aes(y=..density..))+
  scale_x_continuous(limits = c(-2.5,2.5))+
  xlab(" ")+ylab("Density")+
  annotate("text",x=-2.5,y=0.67,label="GL")+
  annotate("text",x=1.76,y=0.67,label="Skewness: -0.11",size=3.3)+
  annotate("text",x=1.9,y=0.57,label="Kurtosis: 2.96",size=3.3)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   #legend.justification = c(0.99,0.01),
                   #legend.position =c(0.99,0.01),
                   legend.position = "right",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1.3,"line"))
PT2=ggplot(dataF2[which(dataF2$Type=="Annual"),],aes(x=Temp))+
  geom_histogram(aes(y=after_stat(density)),
                 bins=10,fill=cbPalette[6],alpha=0.8,color="black")+
  geom_density(aes(y=..density..))+
  scale_x_continuous(limits = c(-1.7,1.7))+
  xlab(" ")+ylab("Density")+
  annotate("text",x=-1.67,y=1.17,label="GO")+
  annotate("text",x=1.22,y=1.17,label="Skewness: 0.00",size=3.3)+
  annotate("text",x=1.295,y=0.99,label="Kurtosis: 3.00",size=3.3)+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   #legend.justification = c(0.99,0.01),
                   #legend.position =c(0.99,0.01),
                   legend.position = "right",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1.3,"line"))

PT1/PT2
PT1=ggplot(dataF1[which(dataF1$Type=="Monthly"),],aes(x=Temp))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth=0.5,fill=cbPalette[4],alpha=0.8,color="black")+
  geom_density(aes(y=..density..))+
  xlab(" ")+ylab("Density")+
  annotate("text",x=-15,y=0.245,label="GL")+
  annotate("text",x=10.5,y=0.245,label="Skewness: -0.28",size=3.3)+
  annotate("text",x=11.4,y=0.209,label="Kurtosis: 3.97",size=3.3)+
  scale_x_continuous(limits = c(-15,15))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   #legend.justification = c(0.99,0.01),
                   #legend.position =c(0.99,0.01),
                   legend.position = "right",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1.3,"line"))
PT2=ggplot(dataF2[which(dataF2$Type=="Monthly"),],aes(x=Temp))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth=0.12,fill=cbPalette[6],alpha=0.8,color="black")+
  geom_density(aes(y=..density..))+
  xlab(" ")+ylab("Density")+
  annotate("text",x=-3.96,y=0.78,label="GO")+
  annotate("text",x=2.91,y=0.78,label="Skewness: 0.21",size=3.3)+
  annotate("text",x=3.06,y=0.665,label="Kurtosis: 4.04",size=3.3)+
  scale_x_continuous(limits = c(-4,4))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   #legend.justification = c(0.99,0.01),
                   #legend.position =c(0.99,0.01),
                   legend.position = "right",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1.3,"line"))
PT1/PT2
PT1=ggplot(dataF1[which(dataF1$Type=="Daily"),],aes(x=Temp))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth=0.6,fill=cbPalette[4],alpha=0.8,color="black")+
  geom_density(aes(y=..density..))+
  xlab(" ")+ylab("Density")+
  annotate("text",x=-21,y=0.152,label="GL")+
  annotate("text",x=14.8,y=0.152,label="Skewness: -0.61",size=3.3)+
  annotate("text",x=16,y=0.129,label="Kurtosis: 4.27",size=3.3)+
  scale_x_continuous(limits = c(-21,21))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   #legend.justification = c(0.99,0.01),
                   #legend.position =c(0.99,0.01),
                   legend.position = "right",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1.3,"line"))
PT2=ggplot(dataF2[which(dataF2$Type=="Daily"),],aes(x=Temp))+
  geom_histogram(aes(y=after_stat(density)),
                 binwidth=0.08,fill=cbPalette[6],alpha=0.8,color="black")+
  geom_density(aes(y=..density..))+
  xlab(" ")+ylab("Density")+
  annotate("text",x=-3.15,y=0.81,label="GO")+
  annotate("text",x=2.3,y=0.81,label="Skewness: 0.83",size=3.3)+
  annotate("text",x=2.45,y=0.68,label="Kurtosis: 5.57",size=3.3)+
  scale_x_continuous(limits = c(-3.2,3.2))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   #legend.justification = c(0.99,0.01),
                   #legend.position =c(0.99,0.01),
                   legend.position = "right",
                   legend.title = element_blank(),
                   legend.key.width=unit(2,"line"),
                   legend.key.height=unit(1.3,"line"))
PT1/PT2
# t2=proc.time()[[3]]
# t2-t1=3.213

### Plot Figures S1(a)--S1(c)
# t1=proc.time()[[3]]
SkewA=read.csv("Figure1/SkewA.csv")$x
SkewM=read.csv("Figure1/SkewM.csv")$x
SkewD=read.csv("Figure1/SkewD.csv")$x
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),
                 SkewA=SkewA,SkewM=SkewM,SkewD=SkewD)
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=SkewA),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 0,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",
                         limits=c(min(c(SkewA,SkewM,SkewD)),max(c(SkewA,SkewM,SkewD))))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color="Skewness")
PT
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=SkewM),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 0,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",
                         limits=c(min(c(SkewA,SkewM,SkewD)),max(c(SkewA,SkewM,SkewD))))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color="Skewness")
PT
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=SkewD),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 0,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",
                         limits=c(min(c(SkewA,SkewM,SkewD)),max(c(SkewA,SkewM,SkewD))))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color="Skewness")
PT
# t2=proc.time()[[3]]
# t2-t1=1.817

### Plot Figures S1(d)--S1(f)
t1=proc.time()[[3]]
KurtA=read.csv("Figure1/KurtA.csv")$x
KurtM=read.csv("Figure1/KurtM.csv")$x
KurtD=read.csv("Figure1/KurtD.csv")$x
idM=which(KurtM>15)
KurtM[idM]=15
idD=which(KurtD>15)
KurtD[idD]=15
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),
                 KurtA=KurtA,KurtM=KurtM,KurtD=KurtD)
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=KurtA),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 3,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",limits=c(min(c(KurtA,KurtM,KurtD)),15))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color="Kurtosis")
PT
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=KurtM),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 3,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",limits=c(min(c(KurtA,KurtM,KurtD)),15))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color="Kurtosis")
PT
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=KurtD),data=dataF,shape=15)+
  scale_colour_gradient2(low = "#377EB8",mid = "white",high = "#E41A1C",
                         midpoint = 3,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",limits=c(min(c(KurtA,KurtM,KurtD)),15))+
  theme_bw()+theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(colour = "black"),
                   axis.text=element_text(size=12),
                   axis.title = element_text(size=14),
                   legend.title = element_text(),
                   legend.key.width=unit(1,"line"),
                   legend.key.height=unit(2,"line"))+
  labs(color="Kurtosis")
PT
t2=proc.time()[[3]]
# t2-t1=1.736



######  Figure 2 in Section 3.2.1 (and Figure S2 in Section S3.2.3)  ###########
# Figures 2 and S2 are used to illustrate the performance of SHT.              #
# All the intermediate outputs can be found in sub-repository "Figure2".       #
################################################################################
### Load the annual data and necessary information
R=7
YN=86
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
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

### Calculate one set of stochastic component Z_9^(1)(L_i,l_j) and plot Figure 2(a)
# t1=proc.time()[[3]]
Dat.mean=apply(Dat[,,2023-2014],2,mean)
Dat.sd=apply(Dat[,,2023-2014],2,sd)
Dat.S1=(Dat[1,,2023-2014]-Dat.mean)/Dat.sd
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degree=Dat.S1)
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
# t2=proc.time()[[3]]
# t2-t1=2.38

### Do SHT for Z_9^(1)(L_i,l_j) and plot Figure 2(b)
# t1=proc.time()[[3]]
DatSHT=emsht_forward(t(matrix(Dat.S1,288,192)),thetas,phis,144)
# t2=proc.time()[[3]]
# t2-t1=4.158
dataF=data.frame(lon=lseq,lat=mseq,Degrees=log10(abs(c(DatSHT))))
PT=ggplot()+xlab("Degree")+ylab("Order")+
  geom_point(aes(x=lon,y=lat,colour=Degrees),size=0.8,data=dataF)+
  scale_x_continuous(limits = c(0,144),breaks=c(0,36*(1:4)))+
  scale_color_gradientn(values = seq(0,1,0.125),colours = Col)+
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
# t2=proc.time()[[3]]
# t2-t1=4.729

### Do inverse SHT with Q=36 and plot Figure 2(c)
# t1=proc.time()[[3]]
DatSHTinv.E1=emsht_inverse(DatSHT,thetas,phis,36)
# t2=proc.time()[[3]]
# t2-t1=0.27
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(c(t(Re(DatSHTinv.E1)))-Dat.S1))  # 2.974017
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degrees),size=0.8,data=dataF)+
  scale_color_gradientn(values = seq(0,1,0.125),colours = Col,limit=c(0,max(dataF$Degrees)))+
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

### Do inverse SHT with Q=72 and plot Figure 2(d)
# t1=proc.time()[[3]]
DatSHTinv2.E1=emsht_inverse(DatSHT,thetas,phis,72)
# t2=proc.time()[[3]]
# t2-t1=0.423
dataF2=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(c(t(Re(DatSHTinv2.E1)))-Dat.S1))
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degrees),size=0.8,data=dataF2)+
  scale_color_gradientn(values = seq(0,1,0.125),colours = Col,limit=c(0,max(dataF$Degrees)))+
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

### Do inverse SHT with Q=116 and plot Figure 2(e)
# t1=proc.time()[[3]]
DatSHTinv3.E1=emsht_inverse(DatSHT,thetas,phis,116)
# t2=proc.time()[[3]]
# t2-t1=1.316
dataF3=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(c(t(Re(DatSHTinv3.E1)))-Dat.S1))
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degrees),size=0.8,data=dataF3)+
  scale_color_gradientn(values = seq(0,1,0.125),colours = Col,limit=c(0,max(dataF$Degrees)))+
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

### Use LatticeKrig to approximate Z_9^(1)(L_i,l_j) and plot Figure 2(f)
# t1=proc.time()[[3]]
Dat.loc1=cbind(rep(c(seq(0,180,by=1.25),seq(-178.75,-1.25,by=1.25)),times=192),Dat.loc[,2])
LKinfol2=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=6,a.wght=5.482,nu=1,LKGeometry="LKSphere")
LKres2=LatticeKrig(Dat.loc1,(Dat[1,,2023-2014]-Dat.mean)/Dat.sd,LKinfo=LKinfol2)
# t2=proc.time()[[3]]
# t2-t1=814.479
write.csv(LKres2$residuals,"Figure2/LatticeKrig_level6_res.csv")
dataF.Lattice=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(LKres2$residuals))
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degrees),size=0.8,data=dataF.Lattice)+
  scale_color_gradientn(values = seq(0,1,0.125),colours = Col,limits=c(0,max(dataF$Degrees)))+
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

### Calculate one set of stochastic component Z_9^(3)(L_i,l_j), do SHT and inverse SHT, and plot Figure S2(a) and S2(c)
# t1=proc.time()[[3]]
Dat.S3=(Dat[3,,2023-2014]-Dat.mean)/Dat.sd
DatSHT=emsht_forward(t(matrix(Dat.S3,288,192)),thetas,phis,144)
DatSHTinv.E3=emsht_inverse(DatSHT,thetas,phis,36)
DatSHTinv2.E3=emsht_inverse(DatSHT,thetas,phis,72)
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(c(t(Re(DatSHTinv.E3)))-Dat.S3))  # 2.203444
dataF2=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(c(t(Re(DatSHTinv2.E3)))-Dat.S3))
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
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degrees),size=0.8,data=dataF2)+
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
# t2=proc.time()[[3]]
# t2-t1=6.024

### Calculate one set of stochastic component Z_69^(1)(L_i,l_j), do SHT and inverse SHT, and plot Figure S2(b) and S2(d)
# t1=proc.time()[[3]]
Dat.mean2=apply(Dat[,,2083-2014],2,mean)
Dat.sd2=apply(Dat[,,2083-2014],2,sd)
Dat.S2=(Dat[1,,2083-2014]-Dat.mean2)/Dat.sd2
DatSHT=emsht_forward(t(matrix(Dat.S2,288,192)),thetas,phis,144)
DatSHTinv.S2=emsht_inverse(DatSHT,thetas,phis,36)
DatSHTinv2.S2=emsht_inverse(DatSHT,thetas,phis,72)
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(c(t(Re(DatSHTinv.S2)))-Dat.S2)) # 2.057448e+00
dataF2=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(c(t(Re(DatSHTinv2.S2)))-Dat.S2))
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
PT=ggplot()+xlab("Longitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=Degrees),size=0.8,data=dataF2)+
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
# t2=proc.time()[[3]]
# t2-t1=7.362



######  Figure 3 in Section 4.1  (and Figure S5 in Section S4.1.1) #############
# Figure 3 and S5 illustrates inference results of annual data.                #
# All the intermediate outputs can be found in sub-repository "Annual/Outputs" #
# Note that Figures 3-5 are all for annual data, so please do not clear the    #
# environment after reproducing Figure 3 if you want to reproduce Figures 4, 5.#
################################################################################
### Load ten ensembles of annual data and necessary information
# t1=proc.time()[[3]]
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
# t2=proc.time()[[3]]
# t2-t1=11.039

### Calculate I.fit values under different values of R and plot Figures 3(a) and S5(a)
# t1=proc.time()[[2]]
R=10   # take R=10 as an example, which maximizes the computational time
Dat.mean=apply(Dat[1:R,,],c(2,3),mean)   # calculate the ensemble mean
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
} # necessary functions to evaluate the mean trend
cl<- makeCluster(10) 
registerDoParallel(cl) 
Dat.hat=foreach(i=1:nrow(Dat.loc),
                .combine=rbind,
                .packages=c("nloptr")) %dopar% mthat(i)
stopCluster(cl)
I.fit.up=I.fit.down=rep(0,nrow(Dat.hat))
for(r in 1:R){
  I.fit.up=I.fit.up+apply((Dat[r,,]-Dat.hat)^2,1,sum)
  I.fit.down=I.fit.down+apply((Dat[r,,]-Dat.mean)^2,1,sum)
}
I.fit=I.fit.up/I.fit.down*(R-1)/R   # calculate I.fit values
# t2=proc.time()[[3]]
# t2-t1=10917.91
I.fit=as.matrix(read.csv("Annual/Outputs/IfitwithRs.csv")[,-1]) # Load the I.fit values for all R
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
dataF=data.frame(lon=c(Dat.loc[,1]),lat=c(Dat.loc[,2]),z=c(I.fit[,6]))
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  # geom_point(aes(x=lon,y=lat,colour=z),data=dataF[-ida,],shape=15)+
  geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
  geom_point(aes(x=lon,y=lat),data=dataF[c(45461,28236),],shape=4)+
  scale_colour_gradient2(low = "#0072B2",mid = "white",high = "#E41A1C",
                         midpoint = 1,space = "Lab",na.value = "grey50",
                         guide = "colourbar",aesthetics = "colour",
                         limit=c(0.8848317,3.6843153))+
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

### Evaluate the deterministic components m_t(L_i,l_j) and sigma(L_i,l_j) using R=7 ensembles and plot Figures 3(b) and S5(b)
R=7  
Dat=Dat[1:R,,]  # Use R=7 ensembles
hat.rho=function(i){
  obj=function(rho){
    value=0
    for(r in 1:R){value=value+profile_negllh(Dat[r,i,],rho)}
    return(value)
  }
  res=bobyqa(x0=c(0.9),fn=obj,lower=c(0.01),upper=c(0.99))
  return(res$par)
} # Estimate rho at each grid point
cl<- makeCluster(4) 
registerDoParallel(cl) 
Res.hatrho=foreach(i=1:nrow(Dat.loc),
                   .combine=cbind,
                   .packages=c("nloptr")) %dopar% hat.rho(i)
stopCluster(cl)
write.csv(c(Res.hatrho),"Annual/Outputs/Res_hatrho.csv") 
# Tim.rho=rep(0,100)
# for(i in 1:100){
#   t1=proc.time()[[3]]
#   hat.rho(sample(1:nrow(Dat.loc),1))
#   t2=proc.time()[[3]]
#   Tim.rho[i]=t2-t1
# }
# mean(Tim.rho)=0.61765
Res.hatrho=read.csv("Annual/Res_hatrho.csv")$x
BetaSighat=function(i){
  X=getX(Res.hatrho[i])
  beta=solve(t(X)%*%X)%*%t(X)%*%t(Dat[,i,])
  Beta=apply(beta,1,mean)
  meanhat=X%*%Beta
  Sigmahat=Dat[,i,]-rep(1,R)%*%t(meanhat)
  sigmahat=mean(Sigmahat^2)
  return(c(Beta,sigmahat,meanhat))
} # Estimate the mean trend m_t(L_i,l_j) and sd sigma(L_i,l_j) at each grid point
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
# Tim.BetaSig=rep(0,100)
# for(i in 1:100){
#   t1=proc.time()[[3]]
#   BetaSighat(sample(1:nrow(Dat.loc),1))
#   t2=proc.time()[[3]]
#   Tim.BetaSig[i]=t2-t1
# }
# mean(Tim.BetaSig)=0.00494
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

### Calculate stochastic component Z_t^{(r)}(L_i,l_j) by detrending m_t(L_i,l_j) and rescaling sigma(L_i,l_j)
# t1=proc.time()[[3]]
Dat.srsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.srsd[r,,]=(Dat[r,,]-Dat.hat)/Sig
}
# t2=proc.time()[[3]]
# t2-t1=1.41

### Do SHT with Q=144 for the stochastic component at each ensemble r and time point t
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
# Without parallel, it will take 4.2*7*86=2528.4 seconds

### Calculate BIC values under different Q values and plot Figure 3(c)
llseq=c(20,30,40,50,60,70,80,90,100)
BIC.land=BIC.ocean=matrix(0,YN*R,length(llseq))
for(j in 1:length(llseq)){
  L=llseq[j]
  Dat.rsdd=Dat.rsd.hat=matrix(0,nrow(Dat.loc),YN*R)
  for(r in 1:R){
    Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.srsd[r,,]
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
    Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.srsd[r,,]
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
    Dat.rsdd[,((r-1)*YN+1):(YN*r)]=Dat.srsd[r,,]
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
# The major computation time is for calculating the inverse SHT. But this step takes time because it calculates 
# the inverse SHT for each ensemble r, each time point t, and all candidates of Q (represented by llseq, llseq1, and llseq2).
# Take Q=100 as an example, which maximizes the computational time of inversing SHT, it takes about
# 0.9 seconds for each ensemble r and time point t. If we run this step without doing parallel,
# we would take at most 0.9*7*86*(9+21+21)=27631.8 seconds.
BIC.land=as.matrix(read.csv("Annual/Outputs/BIC_land.csv")[,-1])
BIC.ocean=as.matrix(read.csv("Annual/Outputs/BIC_ocean.csv")[,-1])
BICd.land=as.matrix(read.csv("Annual/Outputs/BICd_land.csv")[,-1])
BICd.ocean=as.matrix(read.csv("Annual/Outputs/BICd_ocean.csv")[,-1])
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

### Calculate v^2(L_i,l_j) under Q_l=35 and Q_o=69 and plot Figure 3(d)
Ll=35
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
Lo=69
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
# As before, the major computational time in this step is to calculate the inverse SHT with Q=35 and 69 for 
# each ensemble r, time point t. It will take about (0.27+0.42)*7*86=415.38 seconds.
v2hat=read.csv("Annual/Outputs/v2hat.csv")$x
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



######  Figure 4 in Section 4.1   ##############################################
# Figure 4 illustrates the temporal and spatial dependence structures of annual#
# data in the spectral domain.                                                 #
# Assume that we keep all intermediate results of reproducing Figure 3.        #
# All the intermediate outputs can be found in sub-repository "Annual/Outputs".#
################################################################################
### Do the real-valued transformation to SHT coefficients so that they are real values
# t1=proc.time()[[3]]
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
# t2=proc.time()[[3]]
# t2-t1=1.372

### Model the temporal dependence structure using an autoregressive model with order P=1 and plot Figure 4(a)
# bicp=matrix(0,L^2,5)  # Choose the order P by using BIC
# for(i in 1:(L^2)){
#   ts=TDat.rsd.SHT[,i,]
#   for(p in 1:5){
#     u2=rep(0,R)
#     for(r in 1:R){
#       xx=rep(0,YN-p)
#       for(tt in p:1){
#         xx=cbind(xx,ts[r,tt:(tt+YN-p-1)])
#       }
#       xx=as.matrix(xx[,-1])
#       yy=as.matrix(ts[r,(p+1):YN])
#       u2[r]=t(yy)%*%(diag(1,YN-p)-xx%*%solve(t(xx)%*%xx)%*%t(xx))%*%yy/(YN-p)
#     }
#     u2=mean(u2)
#     bicp[i,p]=p*log((YN-p)*R)+R*(YN-p)*log(2*pi)+R*(YN-p)*log(u2)
#   }
# }
# length(which(apply(bicp,1,which.min)==1))/L/L   # 0.9930687
# length(which(apply(bicp,1,which.min)==2))/L/L   # 0.002100399
# length(which(apply(bicp,1,which.min)==3))/L/L   # 0.001890359
# length(which(apply(bicp,1,which.min)==4))/L/L   # 0.0006301197
# length(which(apply(bicp,1,which.min)==5))/L/L   # 0.002310439
# t1=proc.time()[[3]]
TPhi.hat=rep(0,L^2)
for(i in 1:L^2){
  yy=c(TDat.rsd.SHT[,i,2:YN])
  yx=c(TDat.rsd.SHT[,i,-YN])
  TPhi.hat[i]=t(yx)%*%yy/(t(yx)%*%yx)
}
# t2=proc.time()[[3]]
# t2-t1=0.357
write.csv(TPhi.hat,"Annual/Outputs/Phihat.csv")
TPhi.hat=as.matrix(read.csv("Annual/Outputs/Phihat.csv")[,-1])
dataF=data.frame(lon=lseq[1:L^2],lat=mseq[1:L^2],z=TPhi.hat)
PT=ggplot()+xlab(expression(q))+ylab(expression(m))+
  geom_point(aes(x=lon,y=lat,colour=z),size=1.8,shape=15,data=dataF)+
  scale_colour_gradient2(low = "#0072B2",mid = "white",high = "#E41A1C",
                         midpoint = 0,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
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

### Model the spatial dependence by evaluating the covariance matrix of (read-valued) SHT coefficients and plot Figures 4(b) and 4(c)
# t1=proc.time()[[3]]
TK.sample=matrix(0,L^2,L^2)  # evaluate the sample covariance matrix
for(r in 1:R){
  for(t in 1:YN){
    TK.sample=TK.sample+crossprod(t(TDat.rsd.SHT[r,,t]))}
}
TK.sample=TK.sample/R/YN
TK.axial=matrix(0,L^2,L^2) # evaluate the covariance matrix under axial symmetry
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
# t2=proc.time()[[3]]
# t2-t1=127.059
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
PT # 3.7*3
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

### Assess spatial models by selecting auto-covariance at selected latitides and plot Figures 4(d) (or 4(e))
id.test=which(round(Dat.loc[,2],digit=1)==-11.8) 
# id.test=which(round(Dat.loc[,2],digit=1)==36.3)
# t1=proc.time()[[3]]
AC.emp=rep(0,287)   # evaluate the empirical auto-covariance
for(r in 1:R){
  for(t in 1:YN){
    AC.emp=AC.emp+(Dat.srsd[r,id.test[1:287],t])*(Dat.srsd[r,id.test[2:288],t])
  }
}
AC.emp=AC.emp/R/YN
# t2=proc.time()[[3]] 
# t2-t1=0.055
# evaluate auto-covariance for Axial-land/ocean
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
Ang1=Ang[id.test,]
Y=apply(Ang1,1,Yfunc)
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
# t3=proc.time()[[3]]
# t3-t2=59.713
# evaluate auto-covariance for Axial
La=77    # Q=77 is selected by BICs, which are saved in "Outputs/BIC_axial.csv"
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
# t4=proc.time()[[3]]
# t4-t3=95.523, t4-t1=155.291
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



####  Figure 5 in Section 4.1 (and Figures S7 and S8 in Section S4.1.3)  #######
# Figures 5, S7, and S8 illustrate the performance of annual emulations.       #
# Assume that we keep all intermediate results of reproducing Figures 3 and 4. #
# All the intermediate outputs can be found in sub-repository "Annual/Outputs".#
################################################################################
### Caculate the covaraince matrix \tilde U 
# t1=proc.time()[[3]]
TU.axial=TK.axial-outer(TPhi.hat,TPhi.hat,"*")*TK.axial
# t2=proc.time()[[3]]
# t2-t1=0.365
writeMat("Annual/Outputs/U.mat",U=TU.axial) # This file is compressed to reduce size

### Do Cholescky decompostion on \tilde U 
# t1=proc.time()[[3]]
LL.axial=t(chol(TU.axial))
# t2=proc.time()[[3]]
# t2-t1=0.614

### Generate R'=7 ensembles of annual emulations using 4 cores
# t1=proc.time()[[3]]
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
# t2=proc.time()[[3]]
# t2-t1=207.79

### Calculate I.uq values using 4 cores and plot Figures 5(a) and 5(b)
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
  return(get.IUQ(Gen.Dat.Y[,i,],Dat[1:R,i,]))
}
cl=makeCluster(4)
registerDoParallel(cl)
I.uq=foreach(i=1:nrow(Dat.loc),
             .combine = cbind,
             .packages = c("fdaoutlier")) %dopar% findIuq(i)
stopCluster(cl)
# t2=proc.time()[[3]]
# t2-t1=119.055
write.csv(c(I.uq),"Annual/Outputs/Iuq_axialnon.csv")
I.uq.Huang=read.csv("Supplement/Huang/Annual/Outputs/Iuq.csv")$x    # Load the Iuq values of annual emulations generated by HCBG, 
# details about them can be found in "Supplement/Huang/Annual/Huang_Annual.R"
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
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z.H),data=dataF,shape=15)+
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

### Calculate WD_S values using 4 cores and plot Figures 5(c) and 5(d)
# t1=proc.time()[[3]]
WD.time=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  WD.time[i]=wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat.Y[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
# t2=proc.time()[[3]]
# t2-t1=31.127
write.csv(WD.time,"Annual/Outputs/WD_time.csv")
# WD.space=rep(0,YN)
# for(t in 1:YN){
#   WD.space[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat.Y[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
# }
# write.csv(WD.space,"Annual/WD_space.csv")
WD.Huang=read.csv("Supplement/Huang/Annual/Outputs/WD_time.csv")$x  # Load the WD_S values of annual emulations generated by HCBG, 
# details about them can be found in "Supplement/Huang/Annual/Huang_Annual.R"
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
PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
  geom_point(aes(x=lon,y=lat,colour=z.H),data=dataF,shape=15)+
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

### Compare ensemble means and sds of annual emulations with those of simulations and plot Figures S7(a) and S7(c) (or Figures S7(b) and S7(d))
# t1=proc.time()[[3]]
Gen.Dat.mean=apply(Gen.Dat.Y[,,2083-2014],2,mean) 
Gen.Dat.sd=apply(Gen.Dat.Y[,,2083-2014],2,sd)
Dat.mean=apply(Dat[,,2083-2014],2,mean) 
Dat.sd=apply(Dat[,,2083-2014],2,sd)
# t2=proc.time()[[3]]
# t2-t1=2.364
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

### Compare time series of annual emulations with those of simulations and plot Figures S7(e) and S7(f)
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

#### Compare periodograms of emulations with those of simulations and plot Figure S8(a) (or S8(b))
# t1=proc.time()[[3]]
id.test=which(round(Dat.loc[,2],digit=1)==-11.8)
# id.test=which(round(Dat.loc[,2],digit=1)==36.3)
Dat.Gen.srsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.Gen.srsd[r,,]=(Gen.Dat.Y[r,,]-Dat.hat)/Sig
}
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
t2=proc.time()[[3]]
# t2-t1=1.428
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

#####Figures 6 and 7 in Section 4.2 (and Figures S14 and S15 in Section S4.3)  #
# Figures 6 and 7 illustrate the performance of the generated daily emulations.#
# The inference process, which should be done before generating emulations, is #
# illustrated in Figure S14. Figure S15 just shows daily time series.          #
# All these four plots are for daily data. The order should be Figure S14, 6, 7#
# and S15. Therefore, we give the entire procedure here together.              #
# All intermediate outputs are in sub-repository "Daily/Outputs".              #
################################################################################
### Load 7 ensembles of daily data and necessary information
# t1=proc.time()[[3]]
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
lseq=mseq=0
for(i in 1:143){
  lseq=c(lseq,rep(i,2*i+1))
  mseq=c(mseq,seq(-i,i,by=1))
}
# t2=proc.time()[[3]]
# t2-t1=116.522

### Model the deterministic components, evaluate mean trend and sd for each grid point, and plot Figures S14(a) and S14(b)
getX=function(rho){
  M=4 # Please refer lines 79-113 of "Daily/Daily.R" to see how to choose M=4 by using BIC, provided with intermediate results "BICforM.csv"
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
Res.hatrho= foreach(i=1:nrow(Dat.loc),
                    .combine=cbind,
                    .packages=c("nloptr")) %dopar% hat.rho(i)
stopCluster(cl)
write.csv(Res.hatrho,"Daily/Outputs/Res_Hatrho.csv")
# t1=proc.time()[[3]]
# hat.rho(sample(1:nrow(Dat.loc),1))
# t2=proc.time()[[3]]
# t2-t1=39.328
Res.hatrho=read.csv("Daily/Outputs/Res_Hatrho.csv")$x
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
writeMat("Daily/Outputs/BetaAB.mat",BetaAB=as.matrix(Res.BetaSig[1:11,]))
Sig=sqrt(c(Res.BetaSig[12,]))
write.csv(Sig,"Daily/Outputs/Daily_Sig.csv") 
Dat.hat=t(as.matrix(Res.BetaSig[-(1:12),]))
# t1=proc.time()[[3]]
# BetaSighat(sample(1:nrow(Dat.loc),1))
# t2=proc.time()[[3]]
# t2-t1=0.311
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
Sig=read.csv("Daily/Outputs/Daily_Sig.csv")$x
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

### Calculate stochastic components Z_t^{(r)}(L_i,l_j)  by detrending and rescaling 
# t1=proc.time()[[3]]
Dat.srsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.srsd[r,,]=(Dat[r,,]-Dat.hat)/Sig
}
# t2=proc.time()[[3]]
# t2-t1=32.219

### Do SHT with Q=144 for the stochastic component at each ensemble and time point 
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
# As we illustrated before, it takes about 4.2 seconds for each ensemble r and each time point t.
# Without parallel, it will take 4.2*7*(5*365)=53655 seconds

### Calculate BIC values under different Q values and plot Figure S14(c)
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
# The major computation time is for calculating the inverse SHT. But this step takes time because it calculates 
# the inverse SHT for each ensemble r, each time point t, and all candidates of Q (represented by llseq, llseq1, and llseq2).
# Take Q=90 as an example, which maximizes the computational time of inversing SHT, it takes about
# 0.9 seconds for each ensemble r and time point t. If we run this step without doing parallel,
# we would take at most 0.9*7*(365*5)*(8+11+11)=344925 seconds.
BIC.land=as.matrix(read.csv("Daily/Outputs/BIC_land.csv")[,-1])
BIC.ocean=as.matrix(read.csv("Daily/Outputs/BIC_ocean.csv")[,-1])
BICd.land=as.matrix(read.csv("Daily/Outputs/BICd_land.csv")[,-1])
BICd.ocean=as.matrix(read.csv("Daily/Outputs/BICd_ocean.csv")[,-1])
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

### Calculate v^2(L_i,l_j) under Ql=36 and Qo=68 and plot Figure S14(d)
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
# As before, the major computational time in this step is to calculate the inverse SHT with Q=36 and 68 for 
# each ensemble r, time point t. It will take about 0.8*7*365*5=10220 seconds without the parallel.
v2hat=read.csv("Daily/Outputs/v2hat_do.csv")$x
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
# t2-t1=4.248

### Test the normality of coefficients
# t1=proc.time()[[3]]
cBera=rep(0,L^2)
for(i in 1:L^2){
  cBera[i]=jarque.bera.test(c(TDat.rsd.SHT[,i,]))$p.value
}
idBera=which(cBera<=0.05)  # non-normal
length(idBera)/L/L    # 0.5668253
# t2=proc.time()[[3]]
# t2-t1=11.445

### Choose the order of TGH autoregressive model using BIC and plot Figure S14(e)
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
# t1=proc.time()[[3]]
# bicpfunc1(1)
# t2=proc.time()[[3]]
# (t2-t1)*length(idBera)=166635.3
write.csv(bicp1,"Daily/Outputs/bicp_Tukey.csv")
bicp=as.matrix(read.csv("Daily/Outputs/bicp_noTukey.csv")[,-1])     # bicp is an intermediate result
# for SG without using Tukey g-and-h. Please take a look at lines 388-547 of "Daily/Daily.R" for more details.
bicp1=as.matrix(read.csv("Daily/Outputs/bicp_Tukey.csv")[,-1])
bicp[idBera,]=as.matrix(bicp1)
length(which(apply(bicp,1,which.min)==1))/L/L  # 0.5549308
length(which(apply(bicp,1,which.min)==2))/L/L  # 0.07352941
length(which(apply(bicp,1,which.min)==3))/L/L  # 0.3551038
length(which(apply(bicp,1,which.min)==4))/L/L  # 0.01643599
length(which(apply(bicp,1,which.min)==5))/L/L  # 0
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

### Model the temporal dependence structure using a Tukey g-and-h autoregressive model with order P=1 and plot Figure S14(f)
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
  # t1=proc.time()[[3]]
  Tukeyres[j,]=TukeyAuto(TDat.rsd.SHT[,j,])
  # t2=proc.time()[[3]]
  # t2-t1=6.395
}  # 6.395*2621=16761.29
write.csv(Tukeyres,"Daily/Outputs/Tukeyres.csv")
Tukeyres=as.matrix(read.csv("Daily/Outputs/Tukeyres.csv")[,-1])
# t1=proc.time()[[3]]
idBera.new=which(Tukeyres[,1]!=0)    # need Tukey, i.e., set S_{gh}
idphi1=which(Tukeyres[,1]==0)     #    no need Tukey   
idphi2=intersect(which(Tukeyres[,1]!=0),which(Tukeyres[,4]==0))
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
# t2=proc.time()[[3]]
# t2-t1=21.061
# t1=proc.time()[[3]]
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
# t2=proc.time()[[3]]
# t2-t1=8.978
TPhi.hat=read.csv("Daily/Outputs/Phihat_Tukey.csv")$x
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
# t2-t1=80.436

### Calculate the covaraince matrix \check U 
# t1=proc.time()[[3]]
CU.axial=CK.axial-outer(TPhi.hat,TPhi.hat,"*")*CK.axial
# t2=proc.time()[[3]]
# t2-t1=0.255
writeMat("Daily/Outputs/U.mat",U=CU.axial)      # This file is compressed to reduce size

### Generate R'=7 ensembles of daily emulations using 4 cores
# t1=proc.time()[[3]]
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
# t2=proc.time()[[3]]
# t2-t1=115.752
Gen.Dat.Y=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  # t1=proc.time()[[3]]
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
  # t2=proc.time()[[3]]
  # t2-t1=502.918
}
# 502.918*7=3520.426

### Calculate I.uq values and plot Figures 6(a), 6(b), and 6(e)
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
# Tim=0
# for(i in 1:100){
#   t1=proc.time()[[3]]
#   findIuq(sample(1:nrow(Dat.loc),1))
#   t2=proc.time()[[3]]
#   Tim=Tim+t2-t1
# }
# Tim/100=0.1702
# 0.1702*192*288/4=2352.845
write.csv(c(I.uq.Tukey),"Daily/Outputs/Iuq_Tukey.csv")
I.uq=read.csv("Daily/Outputs/Iuq_Tukey.csv")$x
I.uq.Huang=read.csv("Supplement/Huang/Daily/Outputs/Iuq_Huang.csv")$x  # Load the Iuq values of daily emulations generated by HCBG, 
# details about them can be found in "Supplement/Huang/Daily/Huang_Daily.R"
I.uq.noTukey=read.csv("Daily/Outputs/Iuq_noTukey.csv")$x               # Load the Iuq values of daily emulations generated by SHT without TGH, 
# details about them can be found in lines 388-547 of "Daily/Daily.R"
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
I.uq.noTukey=read.csv("Daily/Outputs/Iuq_noTukey.csv")$x
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

### Calculate WD_S values and plot Figures 6(c), 6(d), and 6(f)
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
# Tim=0
# for(i in 1:100){
#   t1=proc.time()[[3]]
#   WDtime(sample(1:nrow(Dat.loc),1))
#   t2=proc.time()[[3]]
#   Tim=Tim+t2-t1
# }
# Tim/100=0.0094
# 0.0094*192*288/4=129.9456
# WD.space=rep(0,YN)
# for(t in 1:YN){
#   WD.space[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat.Y[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
# }
# write.csv(WD.space,"Daily/Outputs/WD_space_Tukey.csv")
WD.time=read.csv("Daily/Outputs/WD_time_Tukey.csv")$x
WD.Huang=read.csv("Supplement/Huang/Daily/Outputs/WD_Huang.csv")$x    # Load the WD values of daily emulations generated by HCBG, 
# details about them can be found in "Supplement/Huang/Daily/Huang_Daily.R"
WD.noTukey=read.csv("Daily/Outputs/WD_time_noTukey.csv")$x            # Load the WD values of daily emulations generated by SHT without TGH, 
# details about them can be found in lines 388-547 of "Daily/Daily.R"
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

### Plot Figures 7 and S15
A=Dat[,c(7785,10377,18153,26793),]
B=Gen.Dat.Y[,c(7785,10377,18153,26793),]
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




