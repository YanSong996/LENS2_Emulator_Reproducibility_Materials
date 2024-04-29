################################################################################################
# This file includes all implementations about LatticeKrig (Figure S3)                         #
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
library(LatticeKrig)
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###############################################################################################
##### Load 7 ensembles of annual data and necessary information
R=7
YN=86
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
forcing=forcing$total[266:351]   # Year 2015--2100
Dat=array(0,c(7,nrow(Dat.loc),YN))   # Original data on sphere
Dat[1,,]=readMat("LENS2_Data/Annual/dat_em1_year.mat")$Dat
Dat[2,,]=readMat("LENS2_Data/Annual/dat_em2_year.mat")$Dat
Dat[3,,]=readMat("LENS2_Data/Annual/dat_em3_year.mat")$Dat
Dat[4,,]=readMat("LENS2_Data/Annual/dat_em4_year.mat")$Dat
Dat[5,,]=readMat("LENS2_Data/Annual/dat_em5_year.mat")$Dat
Dat[6,,]=readMat("LENS2_Data/Annual/dat_em6_year.mat")$Dat
Dat[7,,]=readMat("LENS2_Data/Annual/dat_em7_year.mat")$Dat


############################################################################################
################  Model the stochastic component using LatticeKrig
##### When we use LatticeKrig, we have several parameters to choose.
# Given locations\in[-180,180], the allocation of the center points are fixed.
# Choose startinglevel and nlevel to determine the number of bases, here we choose
# startinglevel=1 and nlevel=5 so that the number of bases is 3420, between 36^2=1296 and 4900
# There are 10238 bases in the sixed level. So we cannot choose such a big number.
# a.wght is always 5.482
# Choose the value of nu. Here given different values of nu, we choose nu which
# achieves a minimum MSE, i.e., nu=1
Dat.loc1=cbind(rep(c(seq(0,180,by=1.25),seq(-178.75,-1.25,by=1.25)),times=192),Dat.loc[,2])
LKinfol1=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=5,a.wght=5.482,nu=0.5,LKGeometry="LKSphere")
LKinfol2=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=5,a.wght=5.482,nu=1,LKGeometry="LKSphere")
LKinfol3=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=5,a.wght=5.482,nu=1.5,LKGeometry="LKSphere")
LKinfol4=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=5,a.wght=5.482,nu=2,LKGeometry="LKSphere")
set.seed(100)
id.test=sample(R*YN,20)
Findnu=function(i){
  r=ceiling(id.test[i]/YN)
  t=id.test[i]-(r-1)*YN
  
  LKres1=LatticeKrig(Dat.loc1,Dat.srsd[r,,t],LKinfo=LKinfol1)
  LKres2=LatticeKrig(Dat.loc1,Dat.srsd[r,,t],LKinfo=LKinfol2)
  LKres3=LatticeKrig(Dat.loc1,Dat.srsd[r,,t],LKinfo=LKinfol3)
  LKres4=LatticeKrig(Dat.loc1,Dat.srsd[r,,t],LKinfo=LKinfol4)
  return(c(mean((LKres1$residuals)^2),mean((LKres2$residuals)^2),mean((LKres3$residuals)^2),
           mean((LKres4$residuals)^2),LKres1$lambda.fixed,LKres2$lambda.fixed,
           LKres3$lambda.fixed,LKres4$lambda.fixed))
}
cl<- makeCluster(2)
registerDoParallel(cl)
Res.FindNu=foreach(i=1:length(id.test),
                   .combine=cbind,
                   .packages=c("LatticeKrig")) %dopar% Findnu(i)
stopCluster(cl)
write.csv(as.matrix(Res.FindNu),"Low_Rank_Approximations/LatticeKrig/Res_FindNu.csv")


##### Replicate Figure 2 and compare LatticeKrig with SHT
# Here we choose nlevel=5 and 6 to see the influence of numbers of bases
Dat.mean=apply(Dat[,,2023-2014],2,mean)
Dat.sd=apply(Dat[,,2023-2014],2,sd)
# Results of SHT in Fig.~2
t1=proc.time()[[3]]
DatSHT=emsht_forward(t(matrix((Dat[1,,2023-2014]-Dat.mean)/Dat.sd,288,192)),thetas,phis,144)
t2=proc.time()[[3]]
DatSHTinv=emsht_inverse(DatSHT,thetas,phis,58)
t3=proc.time()[[3]]
SHT1=abs(c(t(Re(DatSHTinv)))-(Dat[1,,2023-2014]-Dat.mean)/Dat.sd)
t4=proc.time()[[3]]
DatSHTinv2=emsht_inverse(DatSHT,thetas,phis,116)
t5=proc.time()[[3]]
SHT2=abs(c(t(Re(DatSHTinv2)))-(Dat[1,,2023-2014]-Dat.mean)/Dat.sd)
# Results of LatticeKrig in Fig.~2
t1=proc.time()[[3]]
LKinfol3=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=5,a.wght=5.482,nu=1,LKGeometry="LKSphere")
LKres3=LatticeKrig(Dat.loc1,(Dat[1,,2023-2014]-Dat.mean)/Dat.sd,LKinfo=LKinfol3)
t2=proc.time()[[3]]
LKinfol2=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=6,a.wght=5.482,nu=1,LKGeometry="LKSphere")
LKres2=LatticeKrig(Dat.loc1,(Dat[1,,2023-2014]-Dat.mean)/Dat.sd,LKinfo=LKinfol2)
t3=proc.time()[[3]]

### Plot Figs.S3(a)--S3(d)
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(LKres2$residuals))
dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],Degrees=abs(SHT2))
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


##### Then, we still need to evaluate parameters in bases, i.e., rho, sigma, and lambda(key)
# These are used to describe the spatial dependence structure, so we assume that they don't 
# change with time t. 
# Note that we now assume a stationary Gaussian for the data. There are
# options in package which can help use to build non-stationary model in LatticeKrig, e.g., 
# parameters changing with locations. However, we keeps using the stationary one, since 
# the construction of axial-symmetry is one of our contribution and one advantage compared to LatticeKrig.
# We can check that after giving the same lambda, although using different datasets, we still
# get the same wX and wU, which are bases for random and fixed part.
# LKinfol=LKrigSetup(Dat.loc,startingLevel=2,nlevel=5,a.wght=5.482,nu=1.5,lambda=0.01205,LKGeometry="LKSphere")
# LKres1=LatticeKrig(Dat.loc,Dat.srsd[2,,20],LKinfo=LKinfol)
# LKres2=LatticeKrig(Dat.loc,Dat.srsd[3,,20],LKinfo=LKinfol)
# sum(LKres1$wX-LKres2$wX)=0; sum(LKres1$wU-LKres2$wU)=0
# Check the 7th row of Res.FindNu, we can see similar lambda values. More importantly, they are 
# all in 1e-4 scale. Their average is 0.0002906839
# We randomly choose datasets to evaluate the value of lambda
LKinfol=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=5,a.wght=5.482,nu=1,LKGeometry="LKSphere")
set.seed(1)
id.test=sample(R*YN,66)
Findlambda=function(i){
  r=ceiling(id.test[i]/YN)
  t=id.test[i]-(r-1)*YN
  LKres=LatticeKrig(Dat.loc1,Dat.srsd[r,,t],LKinfo=LKinfol)
  return(LKres$lambda.fixed)
}
cl<- makeCluster(2)
registerDoParallel(cl)
Res.FindLambda=foreach(i=1:length(id.test),
                   .combine=cbind,
                   .packages=c("LatticeKrig")) %dopar% Findlambda(i)
stopCluster(cl)
write.csv(c(Res.FindLambda),"Supplement/FigureS3/Outputs/Res_FindLambda.csv")
median(c(Res.FindLambda,Res.FindNu[7,])) # 0.0001236995
mean(c(Res.FindLambda,Res.FindNu[7,]))   # 0.000124893
# Finally, we choose lambda=0.000124893

##### Calculate coefficients 
# Given parameters in covariance and get the bases
LKinfol=LKrigSetup(Dat.loc1,startingLevel=1,nlevel=5,a.wght=5.482,nu=1,lambda=0.000124893,LKGeometry="LKSphere")
t1=proc.time()[[3]]
LKres=LatticeKrig(Dat.loc1,Dat.srsd[2,,20],LKinfo=LKinfol)
t2=proc.time()[[3]]
t2-t1        # 114
wX=LKres$wX
wU=LKres$wU
t1=proc.time()[[3]]
LKres=LKrig(Dat.loc1,Dat.srsd[3,,20],LKinfo=LKinfol,wX=wX,wU=wU)
t2=proc.time()[[3]]
t2-t1  # 6.7

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

### Plot Figs.S3(e) and S3(f)
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




# Use LatticeKrig as an example, if other low-rank approximation methods want to
# have a better approximation, they need thousands of bases. When we model the 
# spatial dependence, without the stationary or axial-symmetric assumption, we 
# need to store the whole covariance matrix, which need a huge storage.
# Even in Huang's paper, they model the data on each latitude, and do expansion
# by latitudes.

