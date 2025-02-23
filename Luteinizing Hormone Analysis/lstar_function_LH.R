library("readxl")
library("tsDyn")

LH1o=read_excel("LH1.xlsx",col_names = FALSE)
LH2o=read_excel("LH2.xlsx",col_names = FALSE)

LH1=LH1o
LH2=LH2o

LH1=unname(as.matrix(LH1o))
LH1
head(LH1)

J=dim(LH1)[2]



LH2=unname(as.matrix(LH2o))
LH2
head(LH2)


zNJ1=LH1
zNJ2=LH2



##############analysis


backlstar=function(xt3,xt2,coe,qt)
{
  gamma=coe[7]
  c=coe[8]
  Gvalue=plogis(xt2,location=c,scale=1/gamma)
  xt1=(xt3-( coe[1]+coe[4]*Gvalue+(coe[2]+coe[5]*Gvalue)*xt2 ) ) / (coe[3]+coe[6]*Gvalue)
  return(unname(xt1))
}


qua_lstar=function(k,xt3,xt2,coe)
{
  backvec=rep(0,40)
  vv3=xt3
  vv2=xt2
  for (i in 1:k)
  {
    vv1=backlstar(vv3,vv2,coe,vv2 )
    vv3=vv2
    vv2=vv1
    backvec[i]=vv1
  }
  return(backvec)
}



findtau_lstar=function(vec)
{
  len=length(vec)
  for (i in 1:len)
  {
    tau=len-i
    if (vec[len-i]<=0){ break}
  }
  mo_lstar=vec[(tau+2):len]/vec[(tau+1):(len-1)] 
  for (j in 1:length(mo_lstar))
  {
    tau_final=tau+j
    if ( mo_lstar[j]<=2){ break}
  }
  return(tau_final)
}


findtau_lstar2_LH=function(vec)
{
  len=length(vec)
  for (i in 1:len)
  {
    tau=len-i
    if ((vec[len-i]<=0) | (vec[len-i]>vec[len-i+1]) && (vec[len-i]<0.3 ) ){ break}
  }
  mo_lstar=vec[(tau+2):len]/vec[(tau+1):(len-1)] 
  for (j in 1:length(mo_lstar))
  {
    tau_final=tau+j
    if ( mo_lstar[j]<=2){ break}
  }
  return(tau_final)
}



lstar_z0_1=rep(0,16)

lstar_z0_2=rep(0,16)

mo1=rep(0,16)

mo2=rep(0,16)

##############

for (j in c(1:7,9:16) )
{
  # k=findtau_lstar(zNJ1[,j])
  
  k=findtau_lstar2_LH(zNJ1[,j])
  
  sl1=lstar(ts(zNJ1[k:40,j]),m=2,mTh=c(1,0),d=1,)
  
  coe= sl1$coefficients
  
  sl1$coefficients
  
  xt3=sl1$fitted.values[1]
  
  xt2=sl1$fitted.values[2]
  
  zhat=qua_lstar(k+2,xt3,xt2,coe)
  
  lstar_z0_1[j]=zhat[k+2]
  
  length_fitted=length(sl1$fitted.values)
  
  mo1[j]=mean(sl1$fitted.values[-1]/sl1$fitted.values[-length_fitted])
}



for (j in c(1:8,12:15))
{
  # k=findtau_lstar(zNJ2[,j])
  
  k=findtau_lstar2_LH(zNJ2[,j])
  
  sl1=lstar(ts(zNJ2[k:40,j]),m=2,mTh=c(1,0),d=1)
  
  coe= sl1$coefficients
  
  sl1$coefficients
  
  xt3=sl1$fitted.values[1]
  
  xt2=sl1$fitted.values[2]
  
  zhat=qua_lstar(k+2,xt3,xt2,coe)
  
  lstar_z0_2[j]=zhat[k+2]
  
  length_fitted=length(sl1$fitted.values)
  
  mo2[j]=mean(sl1$fitted.values[-1]/sl1$fitted.values[-length_fitted])
}



z01_final=lstar_z0_1[which( lstar_z0_1<1 & lstar_z0_1>0) ]

z02_final=lstar_z0_2[which(lstar_z0_2<1 & lstar_z0_2>0)]

mean(z01_final)/mean(z02_final)



mo1_final=mo1[which(mo1>0)]
mean(mo1_final)
sqrt(var(mo1_final))

mo2_final=mo2[which(mo2>0)]
mean(mo2_final)
sqrt(var(mo2_final))

################### method

tau_vec1=matrix(0,nrow=2,ncol=15)

zNJ1_hat=matrix(0,nrow=40 ,ncol=15)
kk=1
for (j in c(1:7,9:16) )
{
  # k=findtau_lstar(zNJ1[,j])
  
  k=findtau_lstar2_LH(zNJ1[,j])
  
  sl1=lstar(ts(zNJ1[k:40,j]),m=2,mTh=c(1,0),d=1,)
  
  length_fitted=length(sl1$fitted.values)
  
  tauhat=40-length_fitted+1
  
  zNJ1_hat[tauhat:40,kk]=sl1$fitted.values
  
  tau_vec1[1,kk]=tauhat
  
  tau_vec1[2,kk]=40
  
  kk=kk+1

}


J=15
mo1_vec=rep(0,J)
for (j in 1:J)
{
  tau_1=tau_vec1[1,j]
  tau_2=tau_vec1[2,j]
  mo1_vec[j]= mean(zNJ1_hat[(tau_1+1):tau_2,j ]/zNJ1_hat[tau_1:(tau_2-1),j])
}
mo1_vec
mo1_hat=mean(mo1_vec)
ma1_vec=rep(0,J)
for (j in 1:J)
{
  tau_1=tau_vec1[1,j]
  tau_2=tau_vec1[2,j]
  zsum=sum(zNJ1_hat[(tau_1:tau_2),j])
  ma1_vec[j]=zsum*(  (mo1_hat-1)/(mo1_hat^(tau_2+1) - mo1_hat^tau_1 ) )
}
ma1_hat=mean(ma1_vec)





tau_vec2=matrix(0,nrow=2,ncol=12)

zNJ2_hat=matrix(0,nrow=40 ,ncol=12 )
kk=1
for (j in c(1:8,12:15))
{
  # k=findtau_lstar(zNJ2[,j])
  
  k=findtau_lstar2_LH(zNJ2[,j])
  
  sl1=lstar(ts(zNJ2[k:40,j]),m=2,mTh=c(1,0),d=1,)
  
  length_fitted=length(sl1$fitted.values)
  
  tauhat=40-length_fitted+1
  
  zNJ2_hat[tauhat:40,kk]=sl1$fitted.values
  
  tau_vec2[1,kk]=tauhat
  
  tau_vec2[2,kk]=40
  
  kk=kk+1
  
}

J=12
mo2_vec=rep(0,J)
for (j in 1:J)
{
  tau_1=tau_vec2[1,j]
  tau_2=tau_vec2[2,j]
  mo2_vec[j]= mean(zNJ2_hat[(tau_1+1):tau_2,j ]/zNJ2_hat[tau_1:(tau_2-1),j])
}
mo2_vec
mo2_hat=mean(mo2_vec)
ma2_vec=rep(0,J)
for (j in 1:J)
{
  tau_1=tau_vec2[1,j]
  tau_2=tau_vec2[2,j]
  zsum=sum(zNJ2_hat[(tau_1:tau_2),j])
  ma2_vec[j]=zsum*(  (mo2_hat-1)/(mo2_hat^(tau_2+1) - mo2_hat^tau_1 ) )
}
ma2_hat=mean(ma2_vec)




sig1=var(ma1_vec)/ma1_hat^2/J
sig2= var(ma2_vec)/ma2_hat^2/J
R_hat=ma1_hat/ma2_hat
v_hat=R_hat^2*(sig1+sig2)


R_hat
c(R_hat-qnorm(0.975)*sqrt(v_hat),R_hat+qnorm(0.975)*sqrt(v_hat))

2*qnorm(0.975)*sqrt(v_hat)





