################################################################################################
# This file includes all implementations about Monthly Case Study using HCBG method            #  
# The procedures are similar to those in "Supplement/Huang/Annual/Huang_Annual.R" but have     #
# additional steps for TGH. We will not illustrate the computational time here, since it can   #
# be evaluated from the annual case.                                                           #
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
library(moments)
library(approxOT)
library(mvtnorm)
Col=c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FEE08B","#FDAE61","#F46D43","#D53E4F")
cbPalette=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

###############################################################################################
##### Load the data
R=7
YN=86*12
Dat.loc=as.matrix(read.csv("LENS2_Data/Dat_loc.csv")[,-1])
forcing=read.csv("LENS2_Data/ERF_ssp370_1750-2500.csv")
id.land=read.csv("LENS2_Data/landid.csv")$x
landFlag=rep(0,nrow(Dat.loc))
landFlag[id.land]=1
landFlag=t(matrix(landFlag,288,192))
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
getX=function(rho){
  X=matrix(0,YN,3+2*3)
  X[,1]=rep(1,YN)
  for(i in 1:YN){
    ty=ceiling(i/12)
    X[i,2]=forcing$total[265+ty] # forcing$X0[265]=2014
    X[i,3]=(1-rho)*(t(rho^seq(ty+263,0,by=-1))%*%forcing$total[1:(264+ty)])
    for(j in 1:3){
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
Res.hatrho= foreach(i=1:nrow(Dat.loc),
                    .combine=cbind,
                    .packages=c("nloptr")) %dopar% hat.rho(i)
stopCluster(cl)
# Res.hatrho=read.csv("Monthly/Outputs/Res_Hatrho.csv")

Dathat=function(i){
  X=getX(Res.hatrho[i])
  beta=solve(t(X)%*%X)%*%t(X)%*%t(Dat[,i,])
  Beta=apply(beta,1,mean)
  meanhat=X%*%Beta
  # Sigmahat=Dat[,i,]-rep(1,R)%*%t(meanhat)
  # sigmahat=mean(Sigmahat^2)
  # return(c(meanhat,sigmahat))
  return(meanhat)
}
cl<- makeCluster(4)
registerDoParallel(cl)
Dat.hat= foreach(i=1:nrow(Dat.loc),
                 .combine=cbind,
                 .packages=c("nloptr")) %dopar% Dathat(i)
stopCluster(cl)
Dat.hat=t(Dat.hat)


######## Calculate the residuals
Dat.rsd=array(0,c(R,nrow(Dat.loc),YN))
for(r in 1:R){
  Dat.rsd[r,,]=Dat[r,,]-Dat.hat
}

### Test the normality
Ske=Kur=Bera=rep(0,nrow(Dat.loc))
for(i in 1:nrow(Dat.loc)){
  Ske[i]=skewness(c(Dat.rsd[,i,]))
  Kur[i]=kurtosis(c(Dat.rsd[,i,]))
  Bera[i]=jarque.bera.test(c(Dat.rsd[,i,]))$p.value
}
idBera=which(Bera<=0.05)
length(idBera)/nrow(Dat.loc)   # 0.9532697
write.csv(Bera,"Supplement/Huang/Monthly/Outputs/Bera.csv")



######### TGH Auto-regressive model for the residuals
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
    res=optim(c(Tukeyinip[1],-Tukeyinip[2],max(Tukeyinip[3],0.01),rep(0.8,1)),TukeyP,lower=c(0.1,0.01,0.01,rep(0.1,1)),method="L-BFGS-B",upper=c(4,0.5,0.5,rep(1,1)))
    Paraest=res$par
    Paraest[2]=-Paraest[2]
  }
  else{
    TukeyP=function(para){return(TukeyPara(para,y))}
    res=optim(c(Tukeyinip[1],Tukeyinip[2],max(Tukeyinip[3],0.01),rep(0.8,1)),TukeyP,lower=c(0.1,0.01,0.01,rep(0.1,1)),method="L-BFGS-B",upper=c(4,0.5,0.5,rep(1,1)))
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
    return(c(0,0,0,rep(0,1)))
  }
  if(which.max(Berares)==2){
    return(c(Tukeyinip[1:3],rep(0,1)))
  }
  if(which.max(Berares)==3){
    return(Paraest)
  }
}
cl<- makeCluster(6)
registerDoParallel(cl)
Tukeyres= foreach(i=idBera,
                  .combine=cbind,
                  .packages=c("nloptr","tseries","moments")) %dopar% TukeyAuto(Dat.rsd[,i,])
stopCluster(cl)
TukeyRes=matrix(0,4,nrow(Dat.loc))
TukeyRes[,idBera]=Tukeyres
write.csv(TukeyRes,"Supplement/Huang/Monthly/Outputs/TukeyRes.csv")
idBera.new=which(TukeyRes[1,]!=0)    # need Tukey 
idphi1=which(TukeyRes[1,]==0)     #    no need Tukey   
idphi2=intersect(which(TukeyRes[1,]!=0),which(TukeyRes[4,]==0))


### Tukey g-and-h transformation for idBera.new
TGH.Dat.rsd=Dat.rsd
for(j in idBera.new){
  y=Dat.rsd[,j,]
  taugh=function(z){return(TukeyRes[1,j]*(exp(TukeyRes[2,j]*z)-1)/TukeyRes[2,j]*exp(TukeyRes[3,j]*z^2/2))}
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
  TGH.Dat.rsd[,j,]=tildez  
}

### Calculate phi and sigma for TGH.Dat.rsd
Temp.Func=function(i){
  yy=c(t(TGH.Dat.rsd[,i,2:YN]))
  yx=c(t(TGH.Dat.rsd[,i,-YN]))
  phihat=c(solve(t(yx)%*%yx)%*%t(yx)%*%yy)
  sighat=c(t(yy-phihat*yx)%*%(yy-phihat*yx)/R/(YN-1))
  return(c(phihat,sighat))
}
cl<- makeCluster(6)
registerDoParallel(cl)
Tempres= foreach(i=1:nrow(Dat.loc),
                  .combine=cbind,
                  .packages=c()) %dopar% Temp.Func(i)
stopCluster(cl)
TPhi.hat=c(Tempres[1,])
TSig.hat=sqrt(c(Tempres[2,]))
write.csv(TPhi.hat,"Supplement/Huang/Monthly/Outputs/TPhihat.csv")
write.csv(TSig.hat,"Supplement/Huang/Monthly/Outputs/TSighat.csv")


############################################################################################
######################  Model longitude dependence at each latitude   ######################
####### Get the residual eta
Eta=array(0,c(R,nrow(Dat.loc),YN))
for(i in 1:nrow(Dat.loc)){
  epsilon=TGH.Dat.rsd[,i,]
  Eta[,i,1]=epsilon[,1]/TSig.hat[i]
  Eta[,i,2:YN]=(epsilon[,2:YN]-TPhi.hat[i]*epsilon[,1:(YN-1)])/(TSig.hat[i])
}


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

# Estimate psi, alpha, and nu
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

Res.longLO=matrix(0,6,192)
for(j in 18:192){
  Res.longLO[,j]=hat.longLO(j)
  print(j)
}
write.csv(Res.longLO,"Supplement/Huang/Monthly/Outputs/Res_longLO.csv")
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
# 0.9165393 0.5
xi=Res.lat[1]
kappa=Res.lat[2]


############################################################################################
######################       Generate Emulations       ######################
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
  TGenEps=matrix(0,nrow(Dat.loc),YN)
  TGenEps[,1]=TSig.hat*GenEta[r,,1]
  for(t in 2:YN){
    TGenEps[,t]=TPhi.hat*TGenEps[,t-1]+TSig.hat*GenEta[r,,t]
  }
  
  GenEps=TGenEps
  for(j in idBera.new){
    taugh=function(z){return(TukeyRes[1,j]*(exp(TukeyRes[2,j]*z)-1)/TukeyRes[2,j]*exp(TukeyRes[3,j]*z^2/2))}
    GenEps[j,]=taugh(TGenEps[j,])
  }
  
  Gen.Dat[r,,]=Dat.hat+GenEps
}
# writeMat("Huang/Monthly/GenDat7.Mat",GenDat=Gen.Dat[7,,])

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
  return(get.IUQ(Gen.Dat[,i,],Dat[1:R,i,]))
}
cl=makeCluster(4)
registerDoParallel(cl)
I.uq.Huang=foreach(i=1:nrow(Dat.loc),
             .combine = cbind,
             .packages = c("fdaoutlier")) %dopar% findIuq(i)
stopCluster(cl)
write.csv(c(I.uq.Huang),"Supplement/Huang/Monthly/Outputs/Iuq_Huang.csv") 

WDtime=function(i){
  return(wasserstein(X=c(Dat[,i,]),Y=c(Gen.Dat[,i,]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate"))
}
cl=makeCluster(4)
registerDoParallel(cl)
WD.Huang=foreach(i=1:nrow(Dat.loc),
                .combine = cbind,
                .packages = c("approxOT")) %dopar% WDtime(i)
stopCluster(cl)
write.csv(c(WD.Huang),"Supplement/Huang/Monthly/Outputs/WD_Huang.csv") 

WD.space.Huang=rep(0,YN)
for(t in 1:YN){
  WD.space.Huang[t]=wasserstein(X=c(Dat[,,t]),Y=c(Gen.Dat[,,t]),p=1,ground_p = 1,observation.orientation="colwise",method="univariate")
}
write.csv(WD.space.Huang,"Supplement/Huang/Monthly/Outputs/WD_space_Huang.csv")


# I.uq=read.csv("Monthly/Iuq_Tukey.csv")$x
# dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z.H=I.uq.Huang,z=I.uq)
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
# PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
#   geom_point(aes(x=lon,y=lat,colour=z.H),data=dataF[which(dataF$z.H<=max(dataF$z)),],shape=15)+
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
# WD.time=read.csv("Monthly/WD_time_Tukey.csv")$x
# dataF=data.frame(lon=Dat.loc[,1],lat=Dat.loc[,2],z=c(WD.time),z.H=c(WD.Huang))
# PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
#   geom_point(aes(x=lon,y=lat,colour=z),data=dataF,shape=15)+
#   scale_colour_gradient(low = "#3288BD",high ="#FEE08B",limit=c(0,max(WD.time)))+
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
# 
# PT=ggplot()+xlab("Longtitude")+ylab("Latitude")+
#   geom_point(aes(x=lon,y=lat,colour=z.H),data=dataF[which(dataF$z.H<=max(dataF$z)),],shape=15)+
#   scale_colour_gradient(low = "#3288BD",high ="#FEE08B",limit=c(0,max(WD.time)))+
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
# 
