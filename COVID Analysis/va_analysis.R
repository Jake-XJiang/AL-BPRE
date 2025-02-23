##### BPRE script for VA data


### include data
library(readxl)
library(writexl)
library(dplyr)

# population research
vapopo =read_excel("vapop.xlsx",col_names = FALSE)
va_pop=vapopo[-1,c(2,4)]
colnames(va_pop)=c("County Name","pop")
county_name=va_pop$`County Name`
summary(va_pop$pop)

pop_vec=(va_pop$pop)
pop_vec[order(pop_vec)]



vapop1 <- va_pop %>% 
filter( pop < 110000 & pop > 89000)
dim(vapop1)
mean(vapop1$pop)

vapop2 <- va_pop %>%
  filter(pop > 37000 & pop < 41000) 
dim(vapop2)
mean(vapop2$pop)

# va <- read_excel("va.xlsx")
# View(va)
# 

### data preparation
 
# before change city
va=read_excel("va_daily.xlsx")
puredata_day=as.matrix(va[,-1])

# after change city
# va=read_excel("va_daily_combined.xlsx")
# puredata_day=as.matrix(va_data[,-1])


#
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
va_week=data.frame( county_name , puredata_week)



ratio_chart=data.frame(puredata_week[,-1]/puredata_week[,-81])
apply(ratio_chart,2,min)


# write_xlsx(va_week, "va_week.xlsx")

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

va_tau_mat=data.frame(county_name ,t(tau_b_e))




# ################# data analysis 1
# 
# # data analysis1 ----------------------------------------------------------
# 
# 
# ## find candidate county based on day of 50 patients
# 
# tau1_st=apply(puredata_week,1,FUN=find1,k=50) # find the day of the 50 patient
# 
# can_con_vec=which((tau_gap>5)& (tau1_st<=10) ) # choose candidate county
# 
# can_con_vec
# 
# length(can_con_vec)
# 
# can_week=puredata_week[can_con_vec,] # get candidate data
# 
# can_pop=purepop[can_con_vec] # get candidate population
# 
# 
# ##
# 
# zNJ1=t(can_week)
# 
# tau_vec1=tau_b_e[,can_con_vec]
# 
# J=dim(can_week)[1]
# 
# mo1_vec=rep(0,J)
# for (j in 1:J)
# {
#   tau_1=tau_vec1[1,j]
#   tau_2=tau_vec1[2,j]
#   mo1_vec[j]= mean(zNJ1[(tau_1+1):tau_2,j ]/zNJ1[tau_1:(tau_2-1),j])
# }
# mo1_vec
# mo1_hat=mean(mo1_vec)
# 
# ma1_vec=rep(0,J)
# for (j in 1:J)
# {
#   tau_1=tau_vec1[1,j]
#   tau_2=tau_vec1[2,j]
#   zsum=sum(zNJ1[(tau_1:tau_2),j])
#   ma1_vec[j]=zsum*(  (mo1_hat-1)/(mo1_hat^(tau_2+1) - mo1_hat^tau_1 ) )
# }
# 
# ma1_vec
# 
# ma1_hat=mean(ma1_vec)
# 
# cbind(ma1_vec,can_pop)
# 
# ra_pop_vec=ma1_vec/can_pop*100
# 
# ra_pop_vec
# 
# mean(ra_pop_vec)
# var(ra_pop_vec)
# 
# mean(ra_pop_vec)/sqrt(var(ra_pop_vec))
# 
# boxplot(ma1_vec)
# 
# boxplot(ra_pop_vec)
# 
# # ----------------
# 
# ################### Analysis for 2 Relative "Quantitation"
# 
# # analysis 2 --------------------------------------------------------------

# data modification
vapop1
vapop2
dim(va_week)

va_week_1 <- va_week %>%
  filter ( county_name  %in% vapop1$`County Name`)
va_tau_1 <- va_tau_mat %>%
  filter ( county_name  %in% vapop1$`County Name`)
va_tau_1$gen=va_tau_1$X2-va_tau_1$X1+1

va_week_2 <- va_week %>%
  filter ( county_name  %in% vapop2$`County Name`)
va_tau_2 <- va_tau_mat %>%
  filter ( county_name  %in% vapop2$`County Name`)
va_tau_2$gen=va_tau_2$X2-va_tau_2$X1+1

cbind(vapop1$pop,va_tau_1)

cbind(vapop2$pop,va_tau_2)

sample_county="Isle of Wight County"
vapop2=vapop2 %>%
  filter( `County Name` != sample_county )
va_week_2 <- va_week %>%
  filter ( county_name  %in% vapop2$`County Name`)
va_tau_2 <- va_tau_mat %>%
  filter ( county_name  %in% vapop2$`County Name`)
va_tau_2$gen=va_tau_2$X2-va_tau_2$X1+1

# # ------------- plot 0 to 40 days
# library(tidyverse)
# 
# va_gp1_plot=cbind("Group1",va_week_1[,1],log(va_week_1[,2:41]))
# 
# colnames(va_gp1_plot)[c(1,2)]=c("Group","County")
# 
# va_gp2_plot=cbind("Group2",va_week_2[,1],log(va_week_2[,2:41]))
# 
# colnames(va_gp2_plot)[c(1,2)]=c("Group","County")
# 
# va_gp_df= data.frame(rbind(va_gp1_plot,va_gp2_plot)) 
# 
# library(tidyr)
# library(dplyr)
# 
# # Convert to long format
# long_df <- pivot_longer(va_gp_df, cols = starts_with("w"), names_to = "Time", values_to = "Data")
# 
# # Ensure group information is included correctly
# long_df$Time <- factor(long_df$Time, levels = paste0("w", 1:40))
# 
# library(ggplot2)
# 
# # Create the line plot
# gpp_plot <- ggplot(long_df, aes(x = Time, y = Data, group = County, color = Group)) +
#   geom_line() +
#   scale_color_manual(values = c("Group1" = "red", "Group2" = "blue"),  # Colors for groups
#                      labels = c("Group1" = "Large Pop", "Group2" = "Small Pop") ) +   # legend label names for groups
#   labs(title = "Log of Weekly Patients of Counties by Grouped by Population",
#        x = "Week",
#        y = "Log patients") +
#   theme_minimal() +
#   theme( plot.caption = element_text(size = 19),              # Change caption text size
#          legend.text = element_text(size = 16) )+ 
#   scale_x_discrete(breaks = c("w1", "w10","w20", "w30","w40"), 
#                    labels = c("0", "10", "20", "30", "40")) 
# 
# # Print the plot
# print(gpp_plot)


# ggsave("gpp_plot.pdf",gpp_plot,width =14,height = 8 )



# -----------------------

# tau modification ----------------------------
# sample_county="Hanover County"
# 
# sample_county="Isle of Wight County"
# 
# sample_county="Warren County"
# 
# sample_ratio_c=va_week[va_week$county_name==sample_county, ]
# sample_ratio_cn=unname(unlist(sample_ratio_c[-1]))
# ISLE_ratio=sample_ratio_cn[-1]/sample_ratio_cn[-length(sample_ratio_cn) ]
# ISLE_ratio[20:40]
# va_week[va_week$county_name==sample_county, 20:40]
# 
# 
# vapop1=vapop1 %>%
#   filter( `County Name` != sample_county )
# va_week_1 <- va_week %>%
#   filter ( county_name  %in% vapop1$`County Name`)
# va_tau_1 <- va_tau_mat %>%
#   filter ( county_name  %in% vapop1$`County Name`)
# 
# va_tau_1$gen=va_tau_1$X2-va_tau_1$X1+1
# 
# va_tau_2[va_tau_2==sample_county,]$X2=29
# 
# va_tau_1$gen=va_tau_1$X2-va_tau_1$X1+1
# 
cbind(vapop1$pop,va_tau_1)
cbind(vapop2$pop,va_tau_2)



va_week_1

va_week_2

va_tau_1

va_tau_2
# end of tau modification-----------------------

# relative quantitation 


zNJ1=t(unname(as.matrix(va_week_1[,-1] )))

zNJ2=t(unname(as.matrix(va_week_2[,-1] )))

dim(zNJ1)

dim(zNJ2)

tau_vec1 = unname(t(va_tau_1[,c(2,3)]) )

tau_vec2 = unname(t(va_tau_2[,c(2,3)]) )

# fixed_window=20
# 
# tau_vec1[2,]=tau_vec1[1,]+fixed_window
# 
# tau_vec2[2,]=tau_vec2[1,]+fixed_window


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

pop_r=mean(vapop1$pop)/mean(vapop2$pop)

2*qnorm(0.975)*sqrt(v_hat)


c("pat1=100","tau2=50","tau3=1.05")
c(mo1_hat,mo1_sd, ma1_hat,sqrt(sig1))
c(mo2_hat,mo2_sd, ma2_hat,sqrt(sig2))
c(mean(vapop1$pop),mean(vapop2$pop),pop_r)
c(R_hat,c(R_hat-qnorm(0.975)*sqrt(v_hat),R_hat+qnorm(0.975)*sqrt(v_hat)))
cbind(vapop1,ma1_vec,t(tau_vec1) )
cbind(vapop2,ma2_vec,t(tau_vec2))


zNJ1

zNJ2

for ( i in 1: length(tau_vec1[1,]))
{
  cat(vapop2$`County Name`[i]," ")
  cat(zNJ1[(tau_vec1[1,i] :tau_vec1[2,i] ),i])
  cat("\n")
}


for ( i in 1: length(tau_vec1[1,]))
{
  cat(vapop2$`County Name`[i]," ")
  cat(zNJ1[(tau_vec1[1,i] :tau_vec1[2,i] ),i])
  cat("\n")
}
ma2_vec


