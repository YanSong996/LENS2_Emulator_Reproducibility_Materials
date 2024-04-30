###################################################################################
# This file includes aggregation Monthly data to Annual data (Fig.S13)            #
###################################################################################
library(R.matlab)
library(fdaoutlier)
library(approxOT)
library(foreach)
library(doParallel)
R=7
YN=86
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])

### Load the generated monthly data (Assume that we have got monthly emulations)
Gen.Dat.Y=array(0,c(R,nrow(Dat.loc),YN*12))
Gen.Dat.Y[1,,]=readMat("Monthly/GenDatTukey1.mat")$GenDat
Gen.Dat.Y[2,,]=readMat("Monthly/GenDatTukey2.mat")$GenDat
Gen.Dat.Y[3,,]=readMat("Monthly/GenDatTukey3.mat")$GenDat
Gen.Dat.Y[4,,]=readMat("Monthly/GenDatTukey4.mat")$GenDat
Gen.Dat.Y[5,,]=readMat("Monthly/GenDatTukey5.mat")$GenDat
Gen.Dat.Y[6,,]=readMat("Monthly/GenDatTukey6.mat")$GenDat
Gen.Dat.Y[7,,]=readMat("Monthly/GenDatTukey7.mat")$GenDat

### Aggregate monthly data to annual data
Gen.Dat.year=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  for(i in 1:YN){
    Gen.Dat.year[r,,i]=apply(Gen.Dat.Y[r,,(i-1)*12+(1:12)],1,mean)
  }
}



##### Test the performance of the annually aggregated generated data
### Load the simulated annual data
Dat=array(0,c(R,nrow(Dat.loc),YN))
Dat[1,,]=readMat("LENS2_Data/Annual/dat_em1_year.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Annual/dat_em2_year.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Annual/dat_em3_year.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Annual/dat_em4_year.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Annual/dat_em5_year.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Annual/dat_em6_year.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Annual/dat_em7_year.mat")$Dat

### Check I.fit
Gen.Dat.mean=apply(Gen.Dat.year,c(2,3),mean)
Dat.mean=apply(Dat,c(2,3),mean)
I.fit.up=I.fit.down=rep(0,nrow(Dat.loc))
for(r in 1:R){
  I.fit.up=I.fit.up+apply((Dat[r,,]-Gen.Dat.mean)^2,1,sum)
  I.fit.down=I.fit.down+apply((Dat[r,,]-Dat.mean)^2,1,sum)
}
I.fit.montoan=I.fit.up/I.fit.down*(R-1)/R
write.csv(I.fit.montoan,"Monthly/Aggregate/Ifit_Monthly_to_Annual.csv")

### Plot Fig.S13(e)
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

### Check Iuq
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
write.csv(c(I.uq.montoan),"Monthly/Aggregate/Iuq_Monthly_to_Annual.csv")

### Plot Figs.S13(a) and S13(b)
I.uq=read.csv("Annual/Outputs/Iuq_axialnon.csv")$x
I.uq.montoan=read.csv("Monthly/Aggregate/Iuq_Monthly_to_Annual.csv")$x
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z1=c(I.uq),z2=c(I.uq.montoan))
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

### Check WD
WD.montoan=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  WD.montoan[i]=wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat.year[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(c(WD.montoan),"Monthly/Aggregate/WD_Monthly_to_Annual.csv")

### Plot Figs.S13(c) and S13(d)
WD.time=read.csv("Annual/Outputs/WD_time.csv")$x
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


