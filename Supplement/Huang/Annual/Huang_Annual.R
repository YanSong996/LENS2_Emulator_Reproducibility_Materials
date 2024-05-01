################################################################################################
# This file includes all implementations about Annual Case Study using Huang's method (Fig.S6) #
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
library(mvtnorm)
library(patchwork)
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


###############################################################################################
##### Load 7 ensembles of annual data and necessary information
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

### Load Huang's parameter estimation as initial values
library(reticulate)
np=import("numpy")
npdata=np$load('Supplement/Huang/results_model5_step2.npz')
psiLand=npdata$f[["psiLand"]]
alphaLand=npdata$f[["alphaLand"]]
nuLand=npdata$f[["nuLand"]]
psiOcean=npdata$f[["psiOcean"]]
alphaOcean=npdata$f[["alphaOcean"]]
nuOcean=npdata$f[["nuOcean"]]
plot(1:192,psiLand,"l")
lines(1:192,psiOcean)
plot(1:192,alphaLand,"l")
lines(1:192,alphaOcean)
plot(1:192,nuLand,"l")
lines(1:192,nuOcean)
npdata=np$load('Supplement/Huang/results_model5_step3.npz')
xi=npdata$f[["xi"]]
kappa=npdata$f[["kappa"]]


############################################################################################
######################    Model mean trend and temporal dependence    ######################
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

##### Estimate rho and phi
hat.rhophi=function(i){
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

### Plot Fig.S6(a) 
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


##### Evaluate mean trend and get the residual eta
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
# writeMat("Huang/Annual/Dathat.mat",Dathat=Dat.hat)
write.csv(Sig,"Supplement/Huang/Annual/Outputs/Sig_Annual.csv")

### Plot Fig.S6(b)
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




############################################################################################
######################  Model longitude dependence at each latitude   ######################
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

###### Estimate psi, alpha, and nu
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

### Plot Fig.S6(c)
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
  



############################################################################################
######################       Model the latitude dependence       ######################
# t1=proc.time()[[3]]
get_SigmaFactor=function(psiL,alphaL,nuL,psiO,alphaO,nuO,j){
  tmpL=sqrt(psiL/((alphaL*alphaL+4*sin(omega/2)*sin(omega/2))^(nuL+0.5)))
  tmpO=sqrt(psiO/((alphaO*alphaO+4*sin(omega/2)*sin(omega/2))^(nuO+0.5)))
  
  flag=landFlag[j,]
  SigmaFactor=(outer(flag,tmpL)+outer(1-flag,tmpO))*exp(1i*outer(0:287,omega))/sqrt(288)

  return(SigmaFactor)
}
get_coherence=function(xi,kappa,latDiff){
  return((xi/((1+4*sin(omega/2)*sin(omega/2))^kappa))^(abs(latDiff)))
}
etaFFTstandardized=array(0,c(YN*R,192,288))
for(m in 1:192){
  SigmaFactor=get_SigmaFactor(Res.longLO[1,m],Res.longLO[2,m],Res.longLO[3,m],
                              Res.longLO[4,m],Res.longLO[5,m],Res.longLO[6,m],m)
  for(r in 1:R){
    etaFFTstandardized[(YN*(r-1)+1):(YN*r),m,]=t(Re(solve(SigmaFactor,Eta[r,(288*(m-1)+1):(288*m),])))
  }
}
negllh_lat=function(arg){
  Mlist=1:192
  Msmall=192
  xi=arg[1]
  kappa=arg[2]
  
  # Compute the log determinant
  tmp=get_coherence(xi,kappa,1)
  ldet=YN*R*(Msmall-1)*sum(log(1-tmp*tmp))
  
  value=0
  for(c in 1:(288/2+1)){
    Finv=matrix(0,Msmall,Msmall)
    diag(Finv)=1
    diag(Finv[-1,-1])=diag(Finv[-1,-1])+tmp[c]*tmp[c]
    diag(Finv[-Msmall,-1])=diag(Finv[-Msmall,-1])-tmp[c]
    diag(Finv[-1,-Msmall])=diag(Finv[-1,-Msmall])-tmp[c]
    Finv=Finv/(1-tmp[c]*tmp[c])
    
    FinvSparse=Matrix(Finv,sparse=TRUE)
    if(c==1 || c==(288/2+1)){
      value=value+sum(etaFFTstandardized[,Mlist,c]%*%FinvSparse*etaFFTstandardized[,Mlist,c])
    }
    else{
      value=value+2*sum(etaFFTstandardized[,Mlist,c]%*%FinvSparse*etaFFTstandardized[,Mlist,c])
    }
  }
  value=value+ldet
  return(value)
}
Res.lat=bobyqa(x0=c(0.9,2),fn=negllh_lat,lower=c(1e-2,0.5),upper=c(0.99999,10))$par
# 0.9287769 0.5000000
xi=Res.lat[1]
kappa=Res.lat[2]
# t2=proc.time()[[3]]
# t2-t1=44.264


############################################################################################
######################       Generate Emulations       ######################
# t1=proc.time()[[3]]
tmp=get_coherence(c(xi),c(kappa),1)
# Generate e
set.seed(666)
EE=array(0,c(YN*R,192,288))
for(t in 1:(YN*R)){
  for(c in 1:(288/2+1)){
    EE[t,,c]=rmvnorm(1,rep(0,192),diag(c(1,rep(sqrt(1-tmp[c]*tmp[c]),191))))
    if(c>1 && c<(288/2+1)){
      EE[t,,2*145-c]=EE[t,,c]
    }
  }
}
# Generate \tilde\eta
etaFFT=array(0,c(YN*R,192,288))
for(t in 1:(YN*R)){
  etaFFT[t,1,]=EE[t,1,]
  for(m in 2:192){
    etaFFT[t,m,]=tmp*etaFFT[t,m-1,]+EE[t,m,]
  }
}
# Calculate \eta
GenEta=array(0,c(R,nrow(Dat.loc),YN))
for(m in 1:192){
  SigmaFactor=get_SigmaFactor(Res.longLO[1,m],Res.longLO[2,m],Res.longLO[3,m],
                              Res.longLO[4,m],Res.longLO[5,m],Res.longLO[6,m],m)
  for(r in 1:R){
    GenEta[r,(288*(m-1)+1):(288*m),]=Re(SigmaFactor%*%t(etaFFT[(YN*(r-1)+1):(YN*r),m,]))
  }
}
# Generate emulations
Gen.Dat=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  GenEps=matrix(0,nrow(Dat.loc),YN)
  GenEps[,1]=Sig*GenEta[r,,1]
  for(t in 2:YN){
    GenEps[,t]=Res.hatrhophi[2,]*GenEps[,t-1]+Sig*GenEta[r,,t]
  }
  Gen.Dat[r,,]=Dat.hat+GenEps
}
# t2=proc.time()[[3]]
# t2-t1=822.321
# writeMat("Huang/Annual/GenDat.Mat",GenDat=Gen.Dat)

### Calculate I.uq
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
  return(get.IUQ(Gen.Dat[,i,],Dat[1:R,i,]))
}
cl=makeCluster(4)
registerDoParallel(cl)
I.uq=foreach(i=1:nrow(Dat.loc),
                   .combine = cbind,
                   .packages = c("fdaoutlier")) %dopar% findIuq(i)
stopCluster(cl)
# t2=proc.time()[[3]]
# t2-t1=118.181
write.csv(c(I.uq),"Supplement/Huang/Annual/Outputs/Iuq.csv") 

### Calculate WD
# t1=proc.time()[[3]]
WDtime=function(i){
  return(wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate"))
}
cl=makeCluster(4)
registerDoParallel(cl)
WD.time=foreach(i=1:nrow(Dat.loc),
                .combine = cbind,
                .packages = c("approxOT")) %dopar% WDtime(i)
stopCluster(cl)
write.csv(c(WD.time),"Supplement/Huang/Annual/Outputs/WD_time.csv") # 0.07640789
# t2=proc.time()[[3]]
# t2-t1=71.064
WD.space.Huang=rep(0,YN)
for(t in 1:YN){
  WD.space.Huang[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(WD.space.Huang,"Supplement/Huang/Annual/Outputs/WD_space_Huang.csv")
# t3=proc.time()[[3]]
# t3-t2=11.637

# dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z=c(I.uq.Huang))
# PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
#   geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
#   scale_colour_gradient2(low ="#377EB8",mid = "white",high = "#E41A1C",
#                          midpoint = 1,space = "Lab",na.value = "grey50",guide = "colourbar",aesthetics = "colour")+
#   #geom_point(aes(x=lon,y=lat),data=dataF[c(7785,10377,18153,26793),],shape=4)+
#   theme_bw()+theme(panel.grid.major = element_blank(),
#                    panel.grid.minor = element_blank(),
#                    panel.background = element_rect(colour = "black"),
#                    axis.text=element_text(size = 12),
#                    axis.title = element_text(size=14),
#                    legend.title = element_text(),
#                    legend.key.width=unit(1,"line"),
#                    legend.key.height=unit(2,"line"))+
#   labs(color=expression(I[uq]))
# PT
# 
# dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z=c(WD.time))
# PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
#   geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
#   scale_colour_gradient(low = "#3288BD",high ="#FEE08B")+
#   theme_bw()+theme(panel.grid.major = element_blank(),
#                    panel.grid.minor = element_blank(),
#                    panel.background = element_rect(colour = "black"),
#                    axis.text=element_text(size=12),
#                    axis.title = element_text(size=14),
#                    legend.title = element_text(),
#                    legend.key.width=unit(1,"line"),
#                    legend.key.height=unit(2,"line"))+
#   labs(color=expression(WD[S]))
# PT

