##### BPRE script for md data
library(readxl)
library(writexl)
library(dplyr)

### include data
mdpopo <- read_excel("md_pop.xlsx")
md_pop=mdpopo[,c(2,4)]
colnames(md_pop)=c("County Name","pop")
county_name=md_pop$`County Name`

pop_vec=(md_pop$pop)
pop_vec[order(pop_vec)]

length(pop_vec)

mdpop1 <- md_pop %>% 
  filter(  pop > 120000)
dim(mdpop1)
mean(mdpop1$pop)

mdpop2 <- md_pop %>%
  filter(pop > 30000 & pop < 120000) 
dim(mdpop2)
mean(mdpop2$pop)

mean(mdpop1$pop)/mean(mdpop2$pop)
### data

md<- read_excel("md_daily.xlsx")
puredata_day=as.matrix(md[,-1])
dim(puredata_day)


## week data prepration 

# fun. to find the first day larger than k patients
find1=function(vec,k)  
{
  for (i in 1:(length(vec)))
  {
    if (vec[i]>=k){ break}
  }
  return (i)
}


startday=apply(puredata_day,1,FUN=find1,k=1) # find the day of the 1 patient
length(startday)
min(startday)
max(startday)


week_day_vec=min(startday)+(0:80)*7 
puredata_week=puredata_day[,week_day_vec] # get pure week data 
colnames(puredata_week)=paste0("w",1:81)

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

puredata_week=t(apply(puredata_week,1,FUN=inc_ensure))


colnames(puredata_week)=paste0("w",1:81)
md_week=data.frame( county_name , puredata_week)

# write_xlsx(md_week, "md_week.xlsx")

ratio_chart=data.frame(puredata_week[,-1]/puredata_week[,-81])
apply(ratio_chart,2,min)



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

tau_b_e=apply(puredata_week,1,FUN=find_range,pat1=100,tau2=2,tau3=1.05) # find start gen. and end gen.
tau_b_e
tau_gap=tau_b_e[2,]-tau_b_e[1,]  # generations width in use

md_tau_mat=data.frame(county_name ,t(tau_b_e))


# # analysis 2 --------------------------------------------------------------

# data modification

mdpop1
mdpop2
dim(md_week)

md_week_1 <- md_week %>%
  filter ( county_name  %in% mdpop1$`County Name`)
md_tau_1 <- md_tau_mat %>%
  filter ( county_name  %in% mdpop1$`County Name`)
md_tau_1$gen=md_tau_1$X2-md_tau_1$X1+1

md_week_2 <- md_week %>%
  filter ( county_name  %in% mdpop2$`County Name`)
md_tau_2 <- md_tau_mat %>%
  filter ( county_name  %in% mdpop2$`County Name`)
md_tau_2$gen=md_tau_2$X2-md_tau_2$X1+1

cbind(mdpop1$pop,md_tau_1)

cbind(mdpop2$pop,md_tau_2)


cbind(md_week_2[,1],md_week_2[,(10:31)]/md_week_2[,(9:30)])


# sample_county="Talbot County"
# mdpop1=mdpop1 %>%
#   filter( `County Name` != sample_county )
# md_week_1 <- md_week %>%
#   filter ( county_name  %in% mdpop1$`County Name`)
# md_tau_1 <- md_tau_mat %>%
#   filter ( county_name  %in% mdpop1$`County Name`)


sample_county="Talbot County"
md_tau_2[md_tau_2==sample_county,]$X1=15

md_tau_2[md_tau_2==sample_county,]$X2=27

####


# end of tau modification-----------------------

# relative quantitation 


zNJ1=t(unname(as.matrix(md_week_1[,-1] )))

zNJ2=t(unname(as.matrix(md_week_2[,-1] )))

dim(zNJ1)

dim(zNJ2)

tau_vec1 = unname(t(md_tau_1[,c(2,3)]) )

tau_vec2 = unname(t(md_tau_2[,c(2,3)]) )


### start relative --------
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

ma1_hat
ma2_hat
R_hat
c(R_hat-qnorm(0.975)*sqrt(v_hat),R_hat+qnorm(0.975)*sqrt(v_hat))

pop_r=mean(mdpop1$pop)/mean(mdpop2$pop)

2*qnorm(0.975)*sqrt(v_hat)


c("pat1=100","tau2=50","tau3=1.05")
c(mo1_hat,mo1_sd, ma1_hat,sqrt(sig1))
c(mo2_hat,mo2_sd, ma2_hat,sqrt(sig2))
c(mean(mdpop1$pop),mean(mdpop2$pop),pop_r)
c(R_hat,c(R_hat-qnorm(0.975)*sqrt(v_hat),R_hat+qnorm(0.975)*sqrt(v_hat)))
cbind(mdpop1,ma1_vec,t(tau_vec1) )
cbind(mdpop2,ma2_vec,t(tau_vec2))
