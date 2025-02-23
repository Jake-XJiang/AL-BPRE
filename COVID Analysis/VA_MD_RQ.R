##### BPRE script for VA and MA data


### include data
library(readxl)
library(writexl)
library(dplyr)


# fun. for imputation
inc_ensure <- function(vec)
{
  l=length(vec)
  output=vec
  for (i in l:2)
  {if (output[i] < output[i-1]) 
  {output[i-1] = output[i]}
  }
  return(output)
}



# fun. to find the day range for:
# (counts > tau1) and (growth rate > tau2) and (growth rate < tau3)
find_range=function(vec,pat1,tau2,tau3)
{
  l=length(vec)
  min1=which(vec>=pat1)
  if (length(min1)==0){ return(c(0,0))}
  else
  {
    t1=min(min1)
    t_vec=vec[(t1+1):l]/vec[t1:(l-1)]
    t2=t1+min(which(t_vec<tau2))-1
    t3=t1+min(which(t_vec<tau3))-1
    return(c(t2,t3))
  }
}


##########


vapopo =read_excel("vapop.xlsx",col_names = FALSE)
va_pop=vapopo[-1,c(2,4)]
colnames(va_pop)=c("County Name","pop")
county_name=va_pop$`County Name`
summary(va_pop$pop)

va_pop_vec=(va_pop$pop)
va_pop_vec[order(va_pop_vec)]



mdpopo <- read_excel("md_pop.xlsx")
md_pop=mdpopo[,c(2,4)]
colnames(md_pop)=c("County Name","pop")
county_name=md_pop$`County Name`

md_pop_vec=(md_pop$pop)
md_pop_vec[order(md_pop_vec)]


va_pop$pop[order(va_pop$pop)]
md_pop$pop[order(md_pop$pop)]


vapop3 <- va_pop %>% 
  filter( pop > 90000 & pop < 300000)
dim(vapop3)
mean(vapop3$pop)

mdpop3 <- md_pop %>% 
  filter( pop > 90000 & pop < 300000)
dim(mdpop3)
mean(mdpop3$pop)




#####


va=read_excel("va_daily.xlsx")
va_puredata_day=as.matrix(va[,-1])
va_county_name=va_pop$`County Name`

md<- read_excel("md_daily.xlsx")
md_puredata_day=as.matrix(md[,-1])
md_county_name=md_pop$`County Name`


#####



week_day_vec=45+(0:80)*7 

va_puredata_week=va_puredata_day[,week_day_vec] # get pure week data 
colnames(va_puredata_week)=paste0("w",1:81)
va_puredata_week=t(apply(va_puredata_week,1,FUN=inc_ensure))
colnames(va_puredata_week)=paste0("w",1:81)
va_week=data.frame( va_county_name , va_puredata_week)


md_puredata_week=md_puredata_day[,week_day_vec] # get pure week data 
colnames(md_puredata_week)=paste0("w",1:81)
md_puredata_week=t(apply(md_puredata_week,1,FUN=inc_ensure))
colnames(md_puredata_week)=paste0("w",1:81)
md_week=data.frame( md_county_name , md_puredata_week)

dim(va_week)
dim(md_week)

#####

va_tau_b_e=apply(va_puredata_week,1,FUN=find_range,pat1=100,tau2=2,tau3=1.05) # find start gen. and end gen.
va_tau_b_e
va_tau_gap=va_tau_b_e[2,]-va_tau_b_e[1,]  # generations width in use
va_tau_mat=data.frame(va_county_name ,t(va_tau_b_e))

md_tau_b_e=apply(md_puredata_week,1,FUN=find_range,pat1=100,tau2=2,tau3=1.05) # find start gen. and end gen.
md_tau_b_e
md_tau_gap=md_tau_b_e[2,]-md_tau_b_e[1,]  # generations width in use
md_tau_mat=data.frame(md_county_name ,t(md_tau_b_e))


va_week_3 <- va_week %>%
  filter ( va_county_name  %in% vapop3$`County Name`)
va_tau_3 <- va_tau_mat %>%
  filter ( va_county_name  %in% vapop3$`County Name`)
va_tau_3$gen=va_tau_3$X2-va_tau_3$X1+1


md_week_3 <- md_week %>%
  filter ( md_county_name  %in% mdpop3$`County Name`)
md_tau_3 <- md_tau_mat %>%
  filter ( md_county_name  %in% mdpop3$`County Name`)
md_tau_3$gen=md_tau_3$X2-md_tau_3$X1+1



#####


zNJ1=t(unname(as.matrix(md_week_3[,-1] )))

zNJ2=t(unname(as.matrix(va_week_3[,-1] )))

dim(zNJ1)

dim(zNJ2)

tau_vec1 = unname(t(md_tau_3[,c(2,3)]) )

tau_vec2 = unname(t(va_tau_3[,c(2,3)]) )


##SV1 
J=dim(zNJ1)[2]

### new method 
nJ1=0
mo1_vec=rep(0,J)
mos1_vec=rep(0,J)
for (j in 1:J)
{
  tau_1=tau_vec1[1,j]
  tau_2=tau_vec1[2,j]
  mo1_vec[j]= mean(zNJ1[(tau_1+1):tau_2,j ]/zNJ1[tau_1:(tau_2-1),j])
  mos1_vec[j]= mean((zNJ1[(tau_1+1):tau_2,j ]/zNJ1[tau_1:(tau_2-1),j])^2)
  nJ1=nJ1+length(zNJ1[(tau_1+1):tau_2,j ]/zNJ1[tau_1:(tau_2-1),j])
}
mo1_vec
mo1_hat=mean(mo1_vec)
mo1_sig=(mean(mos1_vec)-(mo1_hat^2))/nJ1
mo1_sd=sqrt(mo1_sig)

ma1_vec=rep(0,J)
#

### old method
# psum=0
# csum=0
# for (j in 1:J)
# {
#   tau_1=tau_vec1[1,j]
#   tau_2=tau_vec1[2,j]
#   psum=psum+ sum(zNJ1[(tau_1+1):tau_2,j ])
#   csum=csum+sum(zNJ1[tau_1:(tau_2-1),j])
# }
# mo1_hat=psum/csum
#

for (j in 1:J)
{
  tau_1=tau_vec1[1,j]
  tau_2=tau_vec1[2,j]
  zsum=sum(zNJ1[(tau_1:tau_2),j])
  ma1_vec[j]=zsum*(  (mo1_hat-1)/(mo1_hat^(tau_2+1) - mo1_hat^tau_1 ) )
}

ma1_hat=mean(ma1_vec)





##SV2  
J=dim(zNJ2)[2]
### new method
nJ2=0
mo2_vec=rep(0,J)
mos2_vec=rep(0,J)
for (j in 1:J)
{
  tau_1=tau_vec2[1,j]
  tau_2=tau_vec2[2,j]
  mo2_vec[j]= mean(zNJ2[(tau_1+1):tau_2,j ]/zNJ2[tau_1:(tau_2-1),j])
  mos2_vec[j]= mean((zNJ2[(tau_1+1):tau_2,j ]/zNJ2[tau_1:(tau_2-1),j])^2)
  nJ2=nJ2+length(zNJ2[(tau_1+1):tau_2,j ]/zNJ2[tau_1:(tau_2-1),j])
}
mo2_vec
mo2_hat=mean(mo2_vec)
mo2_sig=(mean(mos2_vec)-(mo2_hat^2))/nJ2
mo2_sd=sqrt(mo2_sig)

ma2_vec=rep(0,J)
#

### old method
# psum=0
# csum=0
# for (j in 1:J)
# {
#   tau_1=tau_vec2[1,j]
#   tau_2=tau_vec2[2,j]
#   psum=psum+ sum(zNJ2[(tau_1+1):tau_2,j ])
#   csum=csum+sum(zNJ2[tau_1:(tau_2-1),j])
# }
# mo2_hat=psum/csum
#

for (j in 1:J)
{
  tau_1=tau_vec2[1,j]
  tau_2=tau_vec2[2,j]
  zsum=sum(zNJ2[(tau_1:tau_2),j])
  ma2_vec[j]=zsum*(  (mo2_hat-1)/(mo2_hat^(tau_2+1) - mo2_hat^tau_1 ) )
}

ma2_hat=mean(ma2_vec)

sig1=var(ma1_vec)/ma1_hat^2/J
sig2= var(ma2_vec)/ma2_hat^2/J
R_hat=ma1_hat/ma2_hat
v_hat=R_hat^2*(sig1+sig2)

# end relative ------

mo1_hat
mo2_hat
ma1_hat
ma2_hat
R_hat
c(R_hat-qnorm(0.975)*sqrt(v_hat),R_hat+qnorm(0.975)*sqrt(v_hat))

c(mo1_hat,mo1_sd, ma1_hat,sqrt(sig1))
c(mo2_hat,mo2_sd, ma2_hat,sqrt(sig2))

mean(mdpop3$pop)/mean(vapop3$pop)


