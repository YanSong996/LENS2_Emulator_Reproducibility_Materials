library(R.matlab)
###############################################################
###                    reflected_wts                        ###
###############################################################
reflected_wts=function(mm){
  if(mod(mm,2)==0)
    w=2/(1-mm^2)
  else
    w=0
  if(mm==1)
    w=1i*pi/2
  if(mm==-1)
    w=-1i*pi/2
  return(w)
}


###############################################################
###             dl = dl_trapani_eighth(dl, L, el)           ###
# % Calculates *eighth* (for m = 0:el and mm = 0:m) of        #
# el'th plane of a d-matrix for PI/2 using Trapani & Navaza's #
# recursion method. For el>0, require the dl plane to be      #
# computed already with values for el-1.                      #
# L: harmonic band-limit,                                     #
# el: current harmonic index,                                 #
# dl: the Wigner plane computed for el-1,                     #
# the returned dl contains the Wigner plane computed for el.  #
# % Author: Jason McEwen (www.jasonmcewen.org)                #
###############################################################
dl_trapani_eighth=function(dl,L,el){
  if(el==0)
    dl[el+L,el+L]=1
  else{
    dmm=matrix(0,1,el+1)
    dmm[1]=-sqrt((2*el-1)/(2*el))*dl[el-1+L,0+L]   # Eqn (9) of T&N (2006)
    # Eqn (10) of T&N (2006)
    for(mm in 1:el){
      dmm[mm+1]=sqrt(el)/sqrt(2)*sqrt(2*el-1)/sqrt(el+mm)/sqrt(el+mm-1)*dl[el-1+L,mm-1+L]
    }
    # Initialise dl for next el
    for(mm in 0:el){dl[el+L,mm+L]=dmm[mm+1]}
    # Eqn (11) of T&N (2006)
    for(mm in 0:el){
      m=el-1
      dl[m+L,mm+L]=2*mm/sqrt(el-m)/sqrt(el+m+1)*dl[m+1+L,mm+L]
      # Remaining m cases
      if(el-2>=mm){
        for(m in seq(el-2,mm,by=-1)){
          t1=2*mm/sqrt(el-m)/sqrt(el+m+1)*dl[m+1+L,mm+L]
          t2=sqrt(el-m-1)*sqrt(el+m+2)/sqrt(el-m)/sqrt(el+m+1)*dl[m+2+L,mm+L]
          dl[m+L,mm+L]=t1-t2
        }
      }
    }
  }
  return(dl)
}


###############################################################
###             dl = dl_trapani_full(dl, L, el)             ###
# Calculates el'th plane of a d-matrix for PI/2 using         #
# Trapani & Navaza's recursion method. For el>0, require the  #
# dl plane to be computed already with values for el-1.       #
# L: harmonic band-limit                                      #
# el: current harmonic index                                  # 
# dl: the Wigner plane computed for el-1                      #
# the returned dl contains the Wigner plane computed for el.  #
#   Author: Jason McEwen (www.jasonmcewen.org)                #
###############################################################
dl_trapani_full=function(dl,L,el){
  signs=matrix(0,L+1,1)      # Precompute signs
  for(m in seq(0,(L-1),by=2)){
    signs[m+1]=1.0
    signs[m+2]=-1.0
  }
  dl=dl_trapani_eighth(dl,L,el) # Compute eighth plane
  if(el>0){
    # Diagonal symmetry to fill in quarter
    for(m in 0:(el-1)){
      for(mm in (m+1):el){
        dl[m+L,mm+L]=signs[m+1]*signs[mm+1]*dl[mm+L,m+L]
      }
    }
    # Symmetry in m to fill in half
    for(mm in 0:el){
      for(m in -el:(-1)){
        dl[m+L,mm+L]=signs[el+1]*signs[mm+1]*dl[-m+L,mm+L]
      }
    }
    # Symmetry in mm to fill in remaining plane
    for(mm in -el:(-1)){
      for(m in -el:el){
        dl[m+L,mm+L]=signs[el+1]*signs[abs(m)+1]*dl[m+L,-mm+L]
      }
    }
  }
  return(dl)
}


###############################################################
###        flmn=emsht_forward_sub(fr_spatial_ext,L)         ###
###############################################################
emsht_forward_sub=function(fr_spatial_ext,L){
  flmn=matrix(0,1,L^2)
  reflected_wts_vec=complex(matrix(0,1,length(-(2*L-2):(2*L-2))))
  for(m2 in (-2*L*2):(2*L-2)){
    reflected_wts_vec[(2*L-1)+m2]=reflected_wts(m2)
  }
  # Now compute Fmnm; m' is the first one
  Fmnm=fr_spatial_ext    #just to initalize and make sure it is complex
  for(m in -(L-1):(L-1)){
    for(m1 in -(L-1):(L-1)){
      Fmnm[L+m1,L+m]=t(reflected_wts_vec[((L-1)+m1+1):((2*L-1)+m1+(L-1))])%*%fr_spatial_ext[,L+m]
    }
  }
  # Compute Trapani Initialize
  dl_base=matrix(0,2*L-1,2*L-1)
  dl1=dl_trapani_full(dl_base,L,0)
  # This gives \delta_{m,n}^\ell for each \ell. 
  # m varies along row (first index) and n rows along column (second index)
  
  # exp_mat for use in the computation of Wigner-d functions
  for(ell in 0:(L-1)){
    for(m in -ell:ell){
      flmn[ell^2+ell+m+1]=2*pi*1i^(-m)*sqrt((2*ell+1)/(4*pi))*t(dl1[,L+m]*dl1[,L])%*%(Fmnm[,L+m])
    }
    # next Trapani iteration
    if(ell<L-1){
      dl1=dl_trapani_full(dl1,L,ell+1)
    }
  }
  return(flmn)
}


###########################################################################
###                   Spatial to spectral transform                     ###
# flm=emsht_forward(f_spatial,thetas,phis,L)                              #
# f_spatial: a matrix of of length(thetas)*length(phi)                    #
# L: band-limit                                                           #
# flm: a vector of length L^2 containing spherical harmonic coefficients  #
#                                                                         #
# Author: Zubair Khalid                                                   #
# Copyright (C) 2022  Zubair Khalid                                       #
###########################################################################
library(QZ)
library(pracma)
emsht_forward=function(f_spatial,thetas,phis,L){
  Exp_mat_phis1=exp(-1i*outer(phis,(-(L-1):(L-1)),"*"))
  Gmtheta_r=f_spatial%*%Exp_mat_phis1/length(phis)
  thetas_ext=c(thetas,2*pi-flipud(as.matrix(thetas[-c(1,length(thetas))]))) # Extension of domain
  Gmtheta_r_ext=matrix(0,length(thetas_ext),2*L-1)  # Extension of signal
  for(m in -(L-1):(L-1)){
    Gmtheta_r_ext[1:length(thetas),L+m]=Gmtheta_r[,L+m]
    Gmtheta_r_ext[(length(thetas)+1):length(thetas_ext),L+m]=(-1)^(m)*flipud(as.matrix(Gmtheta_r[-c(1,nrow(Gmtheta_r)),L+m]))
  }
  Exp_mat_thetas_ext=exp(1i*outer(thetas_ext,-(L-1):(L-1),"*"))   # Take FFT of Gmtheta_r_ext
  Bmnm=H(Exp_mat_thetas_ext)%*%Gmtheta_r_ext/length(thetas_ext)
  flm=emsht_forward_sub(Bmnm,L)
}

# Test
# FLM=matrix(0,95,L^2)
# for(t in 1:95){
#   fs = t(Dat.1$values[,,t])
#   FLM[t,]= emsht_forward(fs,thetas,phis,L)
# }



############################################################################
###         Gmtheta = emsht_inverse_trapani_sub(flmn,L,thetas )          ###
############################################################################
emsht_inverse_trapani_sub=function(flmn,L,thetas){
  Gmtheta=matrix(0,length(thetas),2*L-1)
  # Trapani Initialize
  dl_base=matrix(0,2*L-1,2*L-1)
  dl1=dl_trapani_full(dl_base,L,0)
  # This gives \delta_{m,n}^\ell for each \ell. 
  # m varies along row (first index) and n rowas along column (second index)
  # exp_mat for use in the computation of Wigner-d functions
  Exp_mat=exp(1i*outer(thetas,-(L-1):(L-1),"*"))
  # Gmnm=matrix(runif((2*L-1)^2),2*L-1,2*L-1)+1i  #first dimension is m'
  Gmnm=matrix(0,2*L-1,2*L-1)
  for(ell in 0:(L-1)){
    for(m in -ell:ell){
      n=0
      delm=dl1[,L+m]
      deln=dl1[,L+n]
      Gmnm[,L+m]=Gmnm[,L+m]+1i^(n-m)*sqrt((2*ell+1)/(4*pi))*(delm*deln)*(flmn[ell^2+ell+m+1])
    }
    # next Trapani iteration
    if(ell<L-1){
      dl1=dl_trapani_full(dl1,L,ell+1)
    }
  }
  for(m in -(L-1):(L-1)){Gmtheta[,L+m]=Exp_mat%*%Gmnm[,L+m]}
  return(Gmtheta)
}


###########################   Inverse DHT   ################################
###         f_spatial  = emsht_inverse( flm,thetas,phis,L )              ###
############################################################################
emsht_inverse=function(flm,thetas,phis,L){
  Gmtheta=emsht_inverse_trapani_sub(flm,L,thetas)
  Exp_mat_phis=exp(1i*outer(phis,-(L-1):(L-1),"*"))
  f_spatial=Gmtheta%*%t(Exp_mat_phis)
  return(f_spatial)
}

