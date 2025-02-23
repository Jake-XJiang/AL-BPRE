
library(readxl)

Ten=read_excel("Ten_replicates_10to7_vs_10to6.xlsx")

sv1p=Ten[1:10,5:44]

sv2p=Ten[11:20,5:44]

# sv1p=Ten_replicates_10to7_vs_10to6[1:10,5:44]
# 
# sv2p=Ten_replicates_10to7_vs_10to6[11:20,5:44]

zNJ1=unname(as.matrix(t(sv1p)))
  
zNJ2=unname(as.matrix(t(sv2p)))


zNJ1

zNJ2

J=dim(zNJ1)[2]

############function
find_tau2=function(vec,tau1,tau2,tau3)
{
  l=length(vec)
  min1=which(vec>=tau1)
  if (length(min1)==0){ return(c(0,0))}
  else
  {
    t1=min(min1)
    t_vec=vec[(t1+1):l]/vec[t1:(l-1)]
    t2=t1+min(which(t_vec>tau2))-1
    t3=t1+max(which(t_vec>tau3))
    return(c(t2,t3))
  }
}
#########

tau_vec1=apply(zNJ1,2,FUN=find_tau2,tau1=.2,tau2=1.55,tau3=1.55)

tau_vec2=apply(zNJ2,2,FUN=find_tau2,tau1=.2,tau2=1.55,tau3=1.55)

tau_vec1

tau_vec2


##SV1  
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
########

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
#########

for (j in 1:J)
{
  tau_1=tau_vec1[1,j]
  tau_2=tau_vec1[2,j]
  zsum=sum(zNJ1[(tau_1:tau_2),j])
  ma1_vec[j]=zsum*(  (mo1_hat-1)/(mo1_hat^(tau_2+1) - mo1_hat^tau_1 ) )
}

ma1_hat=mean(ma1_vec)





##SV2  
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
#########

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
#########

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


R_hat
c(R_hat-qnorm(0.975)*sqrt(v_hat),R_hat+qnorm(0.975)*sqrt(v_hat))


2*qnorm(0.975)*sqrt(v_hat)



mo1_hat
mo1_sd

mo2_hat
mo2_sd






#######################
# library(tidyverse)
# 
# library(latex2exp)
# 
# ma1_plot=ma1_vec*10^8
# 
# ma2_plot=ma2_vec*10^8
# 
# par(mar = c(5, 5.5, 4, 2) + 0.1)
# 
# plot(c(1:10),ma1_plot,ylim=c(0,5),ylab=TeX("$\\hat{m}_A$"),xlab="replicate number")
# 
# points(c(1:10),ma2_plot,pch=4)
# 
# J_vec=1:10
# 
# lh_plot_data=data.frame(J_vec,ma1_plot,ma2_plot)
# lh_plot_data_long <- pivot_longer(lh_plot_data, cols = c(ma1_plot, ma2_plot), names_to = "series", values_to = "ph12")
# 
# finalplot2=ggplot(lh_plot_data_long, aes(x = J_vec, y = ph12 ,color=series)) + 
#   geom_point() +  # Add points
#   scale_y_continuous(limits = c(0, 5)) + # Set y-axis limits
#   theme_minimal() +  # Optional: Use a minimal theme for aesthetics
#   labs(x = "Replicate Number", y = TeX("$\\hat{Z}_{0,j}$"), title = TeX("Estimators of $Z_{0,j}$")) + 
#   theme_light() + # Use the light theme
#   scale_x_continuous(breaks = 1:15) + # Ensure x-axis has marks for each replicate +
#   scale_color_manual(values = c("ma1_plot" = "red", "ma2_plot" = "blue"),
#                      labels = c("SV1", "SV2")  ) + # Manual color assignment
#   labs( ) +
#   theme(panel.background = element_rect(fill = "#e0e0e0"), 
#         plot.background = element_rect(fill = "#f0f0f0"),
#         panel.grid.major.x = element_line(color = "black",  size = 0.05,linetype = "dotted"),  # Change major x grid lines to red
#         panel.grid.minor.x = element_line(color = "black",  size = 0.05,linetype = "dotted"), 
#         panel.grid.major.y = element_line(color = "black",  size = 0.05,linetype = "dotted"),  # Change major x grid lines to red
#         panel.grid.minor.y = element_line(color = "black", size = 0.05,linetype = "dotted")
#   )
# 
# finalplot2
# 
# # ggsave("SV_z0j.pdf", plot = finalplot2, width = 6, height = 4)








