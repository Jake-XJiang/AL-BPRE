# This is a BPRE(i.i.d) simulation with a poisson offspring distribution 
# lambda is determined by gamma distribution

################# data size setting
#####
set.seed(200)

B=5000

## test Setting



N=40

tau1=5

tau2=20

tau3=35


# poisson ancestor
m_A=10
sigma2_A=m_A


b=1
ii=1
jj=1


## Table Setting

Jvec=c(100,200,300,400)
L_Jvec=length(Jvec)
#####


################# parameter setting
#######

############### gamma poisson
# alpha=10
# beta=.03
# mo=1+alpha*beta
# sigma2=alpha*beta^2
# mo2=sigma2+mo^2
# r_o= sqrt(mo2/(mo^2)  )
# gamma2=alpha*beta



###################    beta binomial
alpha=90
beta=10
mo=1+alpha/(alpha+beta)
sigma2=alpha*beta/ ( (alpha+beta)^2*(alpha+beta+1) )
mo2=sigma2+mo^2
r_o= sqrt(mo2/(mo^2)  )
gamma2= (alpha*beta)/((alpha+beta)*(alpha+beta+1) )

###

eta=(mo^2)/mo2
kap=mo/mo2

mfD=m_A*gamma2*((1)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1)/(1-eta))
mfDn=m_A*gamma2*((1-kap^N)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1-eta^N)/(1-eta))

#####


################# B pre setting 
#####
mo1_est_vec=rep(0,B)
ma1_est_vec=rep(0,B)
varma1_est_vec=rep(0,B)
ma1b_vec=rep(0,B)
ma1b_var_vec=rep(0,B)
cilength_phi1=rep(0,B)
cilength_t1=rep(0,B)
cilength_boot1=rep(0,B)


mo2_est_vec=rep(0,B)
ma2_est_vec=rep(0,B)
varma2_est_vec=rep(0,B)
ma2b_vec=rep(0,B)
ma2b_var_vec=rep(0,B)
cilength_phi2=rep(0,B)
cilength_t2=rep(0,B)
cilength_boot2=rep(0,B)


mo3_est_vec=rep(0,B)
ma3_est_vec=rep(0,B)
varma3_est_vec=rep(0,B)
ma3b_vec=rep(0,B)
ma3b_var_vec=rep(0,B)
cilength_phi3=rep(0,B)
cilength_t3=rep(0,B)
cilength_boot3=rep(0,B)

#####


################# Simulation: jj loop
#####

for (jj in 1:L_Jvec)
{
  J=Jvec[jj]
  
  ################# jj loop reset
  ##### 
  # z0J=matrix(0,B,J)
  zNJ=matrix(0,N,J)
  
  cover_num_phi1=0
  cover_num_t1=0
  cover_num_boot1=0
  
  cover_num_phi2=0
  cover_num_t2=0
  cover_num_boot2=0
  
  cover_num_phi3=0
  cover_num_t3=0
  cover_num_boot3=0
  
  
  #####

################# Simulation: B loop
#####


 for (b in (1:B))
 {
  
  #### create data  
  for (j in (1:J) )
  {
    ###poisson ancestor  
    
    # ## Generate  poisson - gamma  distribution
    # test1=rgamma(N,shape=alpha,scale=beta)   # poisson - gamma
    # repeat
    # {
    # z0=rpois(1,m_A)
    # if (z0>0){break}
    # }
    # # Generate Generations
    #  zNJ[1,j]=z0+rpois(1,z0*test1[1])
    # for (n in 2:(N))
    # {
    #    zNJ[n,j]=zNJ[n-1,j]+rpois(1,zNJ[n-1,j]*test1[n])
    # }
    
    
    # # ## Generate binomial - beta distribution
    test1=rbeta(N,shape1=alpha,shape2=beta)   #  binomial - beta
    repeat
    {
      z0=rpois(1,m_A)
      if (z0>0){break}
    }
    # Generate Generations
    zNJ[1,j]=z0+rbinom(1,z0,test1[1])   #
    for (n in 2:N)
    {
      zNJ[n,j]=zNJ[n-1,j]+rbinom(1,zNJ[n-1,j],test1[n])
    }
    
    
    
    # z0J[b,j]=z0
    
  }

  ######## ma1tau1
  ##### 
 
  mo1_hat= mean( zNJ[(tau1+1):(N),]/zNJ[tau1:(N-1),] )
  mo1_est_vec[b]=mo1_hat
  N_J=( ( mo1_hat-1)/( mo1_hat^tau1*( mo1_hat^(N-tau1+1)-1)) )
  ma1_hat= sum( zNJ[(tau1):(N),] )/J * N_J
  ma1_est_vec[b]=ma1_hat
  varma1_hat= 1/(J-1)* sum ( ( apply(zNJ[tau1:N , ],2,sum) *N_J -ma1_hat )^2 )
  varma1_est_vec[b]=varma1_hat/J
  
  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varma1_hat)/sqrt(J)
  cilength_phi1[b]=2*boundl_phi
  if ( ( ma1_hat-boundl_phi <m_A)&( ma1_hat+boundl_phi > m_A) ) { cover_num_phi1=cover_num_phi1+1 }
  
  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varma1_hat)/sqrt(J)
  if ( ( ma1_hat-boundl_t <m_A)&( ma1_hat+boundl_t > m_A) ) { cover_num_t1=cover_num_t1+1 }
  cilength_t1[b]=2*boundl_t
  
  ## bootstrap CI 
  BT=200
  mo1b_hat_vec=rep(0,BT)
  ma1b_hat_vec=rep(0,BT)
  
  for (bt in 1:BT)
  {  
    samp=sample(1:J,size=J,replace=TRUE)
    if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
    ## ma1b(2)
    mo1b_hat_vec[bt]=mean(zNJ[  (tau1+1):(N),samp]/zNJ[ (tau1):(N-1),samp ])
    N_J_b=( ( mo1b_hat_vec[bt]-1)/( mo1b_hat_vec[bt]^tau1*( mo1b_hat_vec[bt]^(N-tau1+1)-1)) )
    ma1b_hat_vec[bt]=sum( zNJ[(tau1):(N),samp] )/J * N_J_b
    if (is.nan(ma1b_hat_vec[bt])==TRUE)  {stop()}
  } # end for bt
  
  ## mo1b_hat=mean(mo1b_hat_vec)
  ma1b_hat=mean(ma1b_hat_vec)
  ma1b_vec[b]=ma1b_hat
  ma1b_var_vec[b]=var(ma1b_hat_vec)
  if ( ( quantile(ma1b_hat_vec,0.025)<m_A)&(quantile(ma1b_hat_vec,0.975)>m_A ) ) { cover_num_boot1=cover_num_boot1+1 }
  cilength_boot1[b]=quantile(ma1b_hat_vec,0.975)-quantile(ma1b_hat_vec,0.025)
  
  #####
    
  
  ######## ma2tau2  
  ##### 
  
  mo2_hat= mean( zNJ[(tau2+1):(N),]/zNJ[tau2:(N-1),] )
  mo2_est_vec[b]=mo2_hat
  N_J=( ( mo2_hat-1)/( mo2_hat^tau2*( mo2_hat^(N-tau2+1)-1)) )
  ma2_hat= sum( zNJ[(tau2):(N),] )/J * N_J
  ma2_est_vec[b]=ma2_hat
  varma2_hat= 1/(J-1)* sum ( ( apply(zNJ[tau2:N , ],2,sum) *N_J -ma2_hat )^2 )
  varma2_est_vec[b]=varma2_hat/J
  

  ##phi CI
  boundl_phi= qnorm(0.975) * sqrt(varma2_hat)/sqrt(J)
  cilength_phi2[b]=2*boundl_phi
  if ( ( ma2_hat-boundl_phi <m_A)&( ma2_hat+boundl_phi > m_A) ) { cover_num_phi2=cover_num_phi2+1 }
  
  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varma2_hat)/sqrt(J)
  if ( ( ma2_hat-boundl_t <m_A)&( ma2_hat+boundl_t > m_A) ) { cover_num_t2=cover_num_t2+1 }
  cilength_t2[b]=2*boundl_t
  
  ## bootstrap CI
  BT=200
  mo2b_hat_vec=rep(0,BT)
  ma2b_hat_vec=rep(0,BT)

  for (bt in 1:BT)
  {  
    samp=sample(1:J,size=J,replace=TRUE)
    if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
    ## ma2b(2)
    mo2b_hat_vec[bt]=mean(zNJ[  (tau2+1):(N),samp]/zNJ[ (tau2):(N-1),samp ])
    N_J_b=( ( mo2b_hat_vec[bt]-1)/( mo2b_hat_vec[bt]^tau2*( mo2b_hat_vec[bt]^(N-tau2+1)-1)) )
    ma2b_hat_vec[bt]=sum( zNJ[(tau2):(N),samp] )/J * N_J_b
    if (is.nan(ma2b_hat_vec[bt])==TRUE)  {stop()}
  } # end for bt
  
  ## mo2b_hat=mean(mo2b_hat_vec)
  ma2b_hat=mean(ma2b_hat_vec)
  ma2b_vec[b]=ma2b_hat
  ma2b_var_vec[b]=var(ma2b_hat_vec)
  if ( ( quantile(ma2b_hat_vec,0.025)<m_A)&(quantile(ma2b_hat_vec,0.975)>m_A ) ) { cover_num_boot2=cover_num_boot2+1 }
  cilength_boot2[b]=quantile(ma2b_hat_vec,0.975)-quantile(ma2b_hat_vec,0.025)
  #####
  
  
  ######## ma3tau3  
  ##### 
  mo3_hat= mean( zNJ[(tau3+1):(N),]/zNJ[tau3:(N-1),] )
  mo3_est_vec[b]=mo3_hat
  N_J=( ( mo3_hat-1)/( mo3_hat^tau3*( mo3_hat^(N-tau3+1)-1)) )
  ma3_hat= sum( zNJ[(tau3):(N),] )/J * N_J
  ma3_est_vec[b]=ma3_hat
  varma3_hat= 1/(J-1)* sum ( ( apply(zNJ[tau3:N , ],2,sum) *N_J -ma3_hat )^2 )
  varma3_est_vec[b]=varma3_hat/J
  
  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varma3_hat)/sqrt(J)
  cilength_phi3[b]=2*boundl_phi
  if ( ( ma3_hat-boundl_phi <m_A)&( ma3_hat+boundl_phi > m_A) ) { cover_num_phi3=cover_num_phi3+1 }
  
  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varma3_hat)/sqrt(J)
  if ( ( ma3_hat-boundl_t <m_A)&( ma3_hat+boundl_t > m_A) ) { cover_num_t3=cover_num_t3+1 }
  cilength_t3[b]=2*boundl_t
  
  ## bootstrap CI 
  BT=200
  mo3b_hat_vec=rep(0,BT)
  ma3b_hat_vec=rep(0,BT)
  
  for (bt in 1:BT)
  {  
    samp=sample(1:J,size=J,replace=TRUE)
    if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
    ## ma3b(2)
    mo3b_hat_vec[bt]=mean(zNJ[  (tau3+1):(N),samp]/zNJ[ (tau3):(N-1),samp ])
    N_J_b=( ( mo3b_hat_vec[bt]-1)/( mo3b_hat_vec[bt]^tau3*( mo3b_hat_vec[bt]^(N-tau3+1)-1)) )
    ma3b_hat_vec[bt]=sum( zNJ[(tau3):(N),samp] )/J * N_J_b
    if (is.nan(ma3b_hat_vec[bt])==TRUE)  {stop()}
  } # end for bt
  
  ## mo3b_hat=mean(mo3b_hat_vec)
  ma3b_hat=mean(ma3b_hat_vec)
  ma3b_vec[b]=ma3b_hat
  ma3b_var_vec[b]=var(ma3b_hat_vec)
  if ( ( quantile(ma3b_hat_vec,0.025)<m_A)&(quantile(ma3b_hat_vec,0.975)>m_A ) ) { cover_num_boot3=cover_num_boot3+1 }
  cilength_boot3[b]=quantile(ma3b_hat_vec,0.975)-quantile(ma3b_hat_vec,0.025)
  
  #####
  
  
    ## process counter
    if ( b%%(B/10)==1 ) 
    {cat(c( (b-1)/(B/10))*10 ) 
      cat("% \n")}
     } # end loop for b

######## result table 
#####
  if (jj==1)
  {
    output_tab1=c(mean(ma1_est_vec),mean(varma1_est_vec),cover_num_phi1/B,cover_num_t1/B,
                  mean(ma1b_vec),mean(ma1b_var_vec),cover_num_boot1/B)
    output_tab2=c(mean(ma2_est_vec),mean(varma2_est_vec),cover_num_phi2/B,cover_num_t2/B,
                  mean(ma2b_vec),mean(ma2b_var_vec),cover_num_boot2/B)
    output_tab3=c(mean(ma3_est_vec),mean(varma3_est_vec),cover_num_phi3/B,cover_num_t3/B,
                  mean(ma3b_vec),mean(ma3b_var_vec),cover_num_boot3/B)
  } else
  {
    output_tab1=rbind(output_tab1,c(mean(ma1_est_vec),mean(varma1_est_vec),cover_num_phi1/B,cover_num_t1/B,
                                    mean(ma1b_vec),mean(ma1b_var_vec),cover_num_boot1/B))
    output_tab2=rbind(output_tab2,c(mean(ma2_est_vec),mean(varma2_est_vec),cover_num_phi2/B,cover_num_t2/B,
                                    mean(ma2b_vec),mean(ma2b_var_vec),cover_num_boot2/B))
    output_tab3=rbind(output_tab3,c(mean(ma3_est_vec),mean(varma3_est_vec),cover_num_phi3/B,cover_num_t3/B,
                                    mean(ma3b_vec),mean(ma3b_var_vec),cover_num_boot3/B))
  }
#####

} # end loop for jj




output_tab1
output_tab2
output_tab3



save(B,N,tau1,tau2,tau3,Jvec,output_tab1, output_tab2, output_tab3, file = "T21result.RData")






