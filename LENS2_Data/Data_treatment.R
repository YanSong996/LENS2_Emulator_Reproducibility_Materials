################################################################################################
# This file includes all data processing procedures                                            #
################################################################################################
library(ncdf4)
library(abind)
library(R.matlab)

###### Daily data : choose data of year 2020, 2040, 2060, 2080, 2100
# Ensemble 1, 2020
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h1.TS.20150101-20241231.nc")
# lon=ncvar_get(daT,varid="lon")
# lat=ncvar_get(daT,varid="lat")
# Dat.loc=cbind(rep(lon,times=length(lat)),rep(lat,each=length(lon)))
# write.csv(LENS2_Data//Dat_loc.csv")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=dat.em.day
# Ensemble 1, 2040
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h1.TS.20350101-20441231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 1, 2060
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h1.TS.20550101-20641231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 1, 2080
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h1.TS.20750101-20841231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 1, 2100
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h1.TS.20950101-21001231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
plot(1:(365*5),a[39969,],"l")
writeMat("LENS2_Data/Daily/dat_em1_day.mat",Dat=a)

# Ensemble 2, 2020
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h1.TS.20150101-20241231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=dat.em.day
# Ensemble 2, 2040
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h1.TS.20350101-20441231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 2, 2060
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h1.TS.20550101-20641231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 2, 2080
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h1.TS.20750101-20841231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 2, 2100
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h1.TS.20950101-21001231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
plot(1:(365*5),a[39969,],"l")
writeMat("LENS2_Data/Daily/dat_em2_day.mat",Dat=a)

# Ensemble 3, 2020
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h1.TS.20150101-20241231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=dat.em.day
# Ensemble 3, 2040
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h1.TS.20350101-20441231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 3, 2060
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h1.TS.20550101-20641231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 3, 2080
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h1.TS.20750101-20841231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 3, 2100
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h1.TS.20950101-21001231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
plot(1:(365*5),a[39969,],"l")
writeMat("LENS2_Data/Daily/dat_em3_day.mat",Dat=a)

# Ensemble 4, 2020
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h1.TS.20150101-20241231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=dat.em.day
# Ensemble 4, 2040
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h1.TS.20350101-20441231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 4, 2060
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h1.TS.20550101-20641231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 4, 2080
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h1.TS.20750101-20841231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 4, 2100
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h1.TS.20950101-21001231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
plot(1:(365*5),a[39969,],"l")
writeMat("LENS2_Data/Daily/dat_em4_day.mat",Dat=a)

# Ensemble 5, 2020
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h1.TS.20150101-20241231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=dat.em.day
# Ensemble 5, 2040
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h1.TS.20350101-20441231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 5, 2060
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h1.TS.20550101-20641231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 5, 2080
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h1.TS.20750101-20841231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 5, 2100
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h1.TS.20950101-21001231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
plot(1:(365*5),a[39969,],"l")
writeMat("LENS2_Data/Daily/dat_em5_day.mat",Dat=a)

# Ensemble 6, 2020
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h1.TS.20150101-20241231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=dat.em.day
# Ensemble 6, 2040
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h1.TS.20350101-20441231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 6, 2060
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h1.TS.20550101-20641231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 6, 2080
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h1.TS.20750101-20841231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 6, 2100
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h1.TS.20950101-21001231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
plot(1:(365*5),a[39969,],"l")
writeMat("LENS2_Data/Daily/dat_em6_day.mat",Dat=a)

# Ensemble 7, 2020
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h1.TS.20150101-20241231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=dat.em.day
# Ensemble 7, 2040
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h1.TS.20350101-20441231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 7, 2060
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h1.TS.20550101-20641231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 7, 2080
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h1.TS.20750101-20841231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
# Ensemble 7, 2100
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h1.TS.20950101-21001231.nc")
Dat.em=ncvar_get(daT)
Dat.em=Dat.em[,,(5*365+1):(6*365)]
dat.em.day=matrix(0,55296,365)
for(i in 1:365){
  dat.em.day[,i]=c(Dat.em[,,i])-273.15
}
a=cbind(a,dat.em.day)
plot(1:(365*5),a[39969,],"l")
writeMat("LENS2_Data/Daily/dat_em7_day.mat",Dat=a)




####### Monthly and annual data
# Ensemble 1
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.011.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em1_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em1_year.mat",Dat=dat.em.year)

# Ensemble 2
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.012.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em2_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em2_year.mat",Dat=dat.em.year)

# Ensemble 3
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.013.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em3_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em3_year.mat",Dat=dat.em.year)

# Ensemble 4
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.014.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em4_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em4_year.mat",Dat=dat.em.year)

# Ensemble 5
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.015.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em5_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em5_year.mat",Dat=dat.em.year)

# Ensemble 6
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.016.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em6_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em6_year.mat",Dat=dat.em.year)

# Ensemble 7
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.017.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em7_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em7_year.mat",Dat=dat.em.year)

# Ensemble 8
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.018.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.018.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.018.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.018.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.018.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.018.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.018.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.018.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.018.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em8_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em8_year.mat",Dat=dat.em.year)

# Ensemble 9
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.019.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.019.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.019.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.019.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.019.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.019.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.019.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.019.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.019.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em9_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em9_year.mat",Dat=dat.em.year)

# Ensemble 10
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.020.cam.h0.TS.201501-202412.nc")
Dat.em1=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.020.cam.h0.TS.202501-203412.nc")
Dat.em2=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.020.cam.h0.TS.203501-204412.nc")
Dat.em3=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.020.cam.h0.TS.204501-205412.nc")
Dat.em4=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.020.cam.h0.TS.205501-206412.nc")
Dat.em5=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.020.cam.h0.TS.206501-207412.nc")
Dat.em6=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.020.cam.h0.TS.207501-208412.nc")
Dat.em7=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.020.cam.h0.TS.208501-209412.nc")
Dat.em8=ncvar_get(daT)
daT=nc_open("LENS2_Data/Raw/b.e21.BSSP370smbb.f09_g17.LE2-1251.020.cam.h0.TS.209501-210012.nc")
Dat.em9=ncvar_get(daT)
Dat.em.month=abind(Dat.em1,Dat.em2,Dat.em3,Dat.em4,Dat.em5,Dat.em6,Dat.em7,Dat.em8,Dat.em9,along=3)
dat.em.month=matrix(0,55296,1032)
for(i in 1:1032){
  dat.em.month[,i]=c(Dat.em.month[,,i])-273.15
}
plot(1:1032,dat.em.month[39969,],"l")
writeMat("LENS2_Data/Monthly/dat_em10_month.mat",Dat=dat.em.month)
dat.em.year=matrix(0,55296,86)
for(i in 1:86){
  dat.em.year[,i]=apply(dat.em.month[,(i-1)*12+(1:12)],1,mean)
}
plot(1:86,dat.em.year[39969,],"l")
writeMat("LENS2_Data/Annual/dat_em10_year.mat",Dat=dat.em.year)




# library(sf)
# library(spData)
# points=data.frame(Var1=Dat.loc[,1],Var2=Dat.loc[,2])
# pts=st_as_sf(points, coords=1:2, crs=4326)
# ## Find which points fall over land
# ii=!is.na(as.numeric(st_intersects(pts, world)))
# il=(1:55296)[ii]
# plot(Dat.loc[il,])
# write.csv(il,'LENS2_Data/landid.csv')


