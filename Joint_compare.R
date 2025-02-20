# This is a BPRE(i.i.d) simulation with a poisson offspring distribution 
# lambda is determined by gamma distribution

################# data size setting
#####
set.seed(200)

B=500
J=400
N=30

tau1_1=6
tau1_2=18

tau2_1=16
tau2_2=23

tau3_1=24
tau3_2=27


# poisson ancestor
m_A=10
sigma2_A=m_A

b=1
ii=1
jj=1


## Table Setting




Jvec=c(50,100,150,200)
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
mo_s1_est_vec=rep(0,B)
ma_s1_est_vec=rep(0,B)
varma_s1_est_vec=rep(0,B)
ma_s1b_vec=rep(0,B)
ma_s1b_var_vec=rep(0,B)
cilength_phi_s1=rep(0,B)
cilength_t_s1=rep(0,B)
cilength_boot_s1=rep(0,B)

mo_e1_est_vec=rep(0,B)
ma_e1_est_vec=rep(0,B)
varma_e1_est_vec=rep(0,B)
ma_e1b_vec=rep(0,B)
ma_e1b_var_vec=rep(0,B)
cilength_phi_e1=rep(0,B)
cilength_t_e1=rep(0,B)
cilength_boot_e1=rep(0,B)


mo_s2_est_vec=rep(0,B)
ma_s2_est_vec=rep(0,B)
varma_s2_est_vec=rep(0,B)
ma_s2b_vec=rep(0,B)
ma_s2b_var_vec=rep(0,B)
cilength_phi_s2=rep(0,B)
cilength_t_s2=rep(0,B)
cilength_boot_s2=rep(0,B)

mo_e2_est_vec=rep(0,B)
ma_e2_est_vec=rep(0,B)
varma_e2_est_vec=rep(0,B)
ma_e2b_vec=rep(0,B)
ma_e2b_var_vec=rep(0,B)
cilength_phi_e2=rep(0,B)
cilength_t_e2=rep(0,B)
cilength_boot_e2=rep(0,B)


mo_s3_est_vec=rep(0,B)
ma_s3_est_vec=rep(0,B)
varma_s3_est_vec=rep(0,B)
ma_s3b_vec=rep(0,B)
ma_s3b_var_vec=rep(0,B)
cilength_phi_s3=rep(0,B)
cilength_t_s3=rep(0,B)
cilength_boot_s3=rep(0,B)

mo_e3_est_vec=rep(0,B)
ma_e3_est_vec=rep(0,B)
varma_e3_est_vec=rep(0,B)
ma_e3b_vec=rep(0,B)
ma_e3b_var_vec=rep(0,B)
cilength_phi_e3=rep(0,B)
cilength_t_e3=rep(0,B)
cilength_boot_e3=rep(0,B)

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
  
  cover_num_phi_s1=0
  cover_num_t_s1=0
  cover_num_boot_s1=0
  cover_num_phi_e1=0
  cover_num_t_e1=0
  cover_num_boot_e1=0
  
  cover_num_phi_s2=0
  cover_num_t_s2=0
  cover_num_boot_s2=0
  cover_num_phi_e2=0
  cover_num_t_e2=0
  cover_num_boot_e2=0
  
  cover_num_phi_s3=0
  cover_num_t_s3=0
  cover_num_boot_s3=0
  cover_num_phi_e3=0
  cover_num_t_e3=0
  cover_num_boot_e3=0   
  
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
  
  ######## ma_s1tau1_s1
  ##### 
  
  mo_s1_hat= mean( zNJ[(tau1_2+1):(N),]/zNJ[tau1_2:(N-1),] )
  mo_s1_est_vec[b]=mo_s1_hat
  N_J=( ( mo_s1_hat-1)/( mo_s1_hat^tau1_1*( mo_s1_hat^(tau1_2-tau1_1+1)-1)) )
  ma_s1_hat= sum( zNJ[(tau1_1):(tau1_2),] )/J * N_J
  ma_s1_est_vec[b]=ma_s1_hat
  varma_s1_hat= 1/(J-1)* sum ( ( apply(zNJ[(tau1_1):(tau1_2) , ],2,sum) *N_J -ma_s1_hat )^2 )
  varma_s1_est_vec[b]=varma_s1_hat/J
  
  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varma_s1_hat)/sqrt(J)
  cilength_phi_s1[b]=2*boundl_phi
  if ( ( ma_s1_hat-boundl_phi <m_A)&( ma_s1_hat+boundl_phi > m_A) ) { cover_num_phi_s1=cover_num_phi_s1+1 }
  
  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varma_s1_hat)/sqrt(J)
  if ( ( ma_s1_hat-boundl_t <m_A)&( ma_s1_hat+boundl_t > m_A) ) { cover_num_t_s1=cover_num_t_s1+1 }
  cilength_t_s1[b]=2*boundl_t
  
  # ## bootstrap CI
  # BT=200
  # mo_s1b_hat_vec=rep(0,BT)
  # ma_s1b_hat_vec=rep(0,BT)
  # 
  # for (bt in 1:BT)
  # {
  #   samp=sample(1:J,size=J,replace=TRUE)
  #   if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
  #   ## ma_s1b(2)
  #   mo_s1b_hat_vec[bt]=mean(zNJ[  (tau1_2+1):(N),samp]/zNJ[ (tau1_2):(N-1),samp ])
  #   N_J_b=( ( mo_s1b_hat_vec[bt]-1)/( mo_s1b_hat_vec[bt]^tau1_1*( mo_s1b_hat_vec[bt]^(tau1_2-tau1_1+1)-1)) )
  #   ma_s1b_hat_vec[bt]=sum( zNJ[(tau1_1):(tau1_2),samp] )/J * N_J_b
  #   if (is.nan(ma_s1b_hat_vec[bt])==TRUE)  {stop()}
  # } # end for bt
  # 
  # ## mo_s1b_hat=mean(mo_s1b_hat_vec)
  # ma_s1b_hat=mean(ma_s1b_hat_vec)
  # ma_s1b_vec[b]=ma_s1b_hat
  # ma_s1b_var_vec[b]=var(ma_s1b_hat_vec)
  # if ( ( quantile(ma_s1b_hat_vec,0.025)<m_A)&(quantile(ma_s1b_hat_vec,0.975)>m_A ) ) { cover_num_boot_s1=cover_num_boot_s1+1 }
  # cilength_boot_s1[b]=quantile(ma_s1b_hat_vec,0.975)-quantile(ma_s1b_hat_vec,0.025)
  
  #####
  
  ######## ma_e1tau1_e1
  ##### 
  
  mo_e1_hat= mean( zNJ[(tau1_1+1):(tau1_2),]/zNJ[tau1_1:(tau1_2-1),] )
  mo_e1_est_vec[b]=mo_e1_hat
  N_J=( ( mo_e1_hat-1)/( mo_e1_hat^tau1_2*( mo_e1_hat^(n-tau1_2+1)-1)) )
  ma_e1_hat= sum( zNJ[(tau1_2):(N),] )/J * N_J
  ma_e1_est_vec[b]=ma_e1_hat
  varma_e1_hat= 1/(J-1)* sum ( ( apply(zNJ[(tau1_2):(N) , ],2,sum) *N_J -ma_e1_hat )^2 )
  varma_e1_est_vec[b]=varma_e1_hat/J
  
  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varma_e1_hat)/sqrt(J)
  cilength_phi_e1[b]=2*boundl_phi
  if ( ( ma_e1_hat-boundl_phi <m_A)&( ma_e1_hat+boundl_phi > m_A) ) { cover_num_phi_e1=cover_num_phi_e1+1 }
  
  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varma_e1_hat)/sqrt(J)
  if ( ( ma_e1_hat-boundl_t <m_A)&( ma_e1_hat+boundl_t > m_A) ) { cover_num_t_e1=cover_num_t_e1+1 }
  cilength_t_e1[b]=2*boundl_t
  
  # ## bootstrap CI
  # BT=200
  # mo_e1b_hat_vec=rep(0,BT)
  # ma_e1b_hat_vec=rep(0,BT)
  # 
  # for (bt in 1:BT)
  # {
  #   samp=sample(1:J,size=J,replace=TRUE)
  #   if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
  #   ## ma_e1b(2)
  #   mo_e1b_hat_vec[bt]=mean(zNJ[  (tau1_1+1):(tau1_2),samp]/zNJ[ tau1_1:(tau1_2-1),samp ])
  #   N_J_b=( ( mo_e1b_hat_vec[bt]-1)/( mo_e1b_hat_vec[bt]^tau1_2*( mo_e1b_hat_vec[bt]^(n-tau1_2+1)-1)) )
  #   ma_e1b_hat_vec[bt]=sum( zNJ[(tau1_2):(N),samp] )/J * N_J_b
  #   if (is.nan(ma_e1b_hat_vec[bt])==TRUE)  {stop()}
  # # } # end for bt
  # 
  # ## mo_e1b_hat=mean(mo_e1b_hat_vec)
  # ma_e1b_hat=mean(ma_e1b_hat_vec)
  # ma_e1b_vec[b]=ma_e1b_hat
  # ma_e1b_var_vec[b]=var(ma_e1b_hat_vec)
  # if ( ( quantile(ma_e1b_hat_vec,0.025)<m_A)&(quantile(ma_e1b_hat_vec,0.975)>m_A ) ) { cover_num_boot_e1=cover_num_boot_e1+1 }
  # cilength_boot_e1[b]=quantile(ma_e1b_hat_vec,0.975)-quantile(ma_e1b_hat_vec,0.025)
  
  #####
  
  
  ######## ma_s2tau2_s2
  #####

  mo_s2_hat= mean( zNJ[(tau2_2+1):(N),]/zNJ[tau2_2:(N-1),] )
  mo_s2_est_vec[b]=mo_s2_hat
  N_J=( ( mo_s2_hat-1)/( mo_s2_hat^tau2_1*( mo_s2_hat^(tau2_2-tau2_1+1)-1)) )
  ma_s2_hat= sum( zNJ[(tau2_1):(tau2_2),] )/J * N_J
  ma_s2_est_vec[b]=ma_s2_hat
  varma_s2_hat= 1/(J-1)* sum ( ( apply(zNJ[(tau2_1):(tau2_2) , ],2,sum) *N_J -ma_s2_hat )^2 )
  varma_s2_est_vec[b]=varma_s2_hat/J

  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varma_s2_hat)/sqrt(J)
  cilength_phi_s2[b]=2*boundl_phi
  if ( ( ma_s2_hat-boundl_phi <m_A)&( ma_s2_hat+boundl_phi > m_A) ) { cover_num_phi_s2=cover_num_phi_s2+1 }

  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varma_s2_hat)/sqrt(J)
  if ( ( ma_s2_hat-boundl_t <m_A)&( ma_s2_hat+boundl_t > m_A) ) { cover_num_t_s2=cover_num_t_s2+1 }
  cilength_t_s2[b]=2*boundl_t

  # ## bootstrap CI
  # BT=200
  # mo_s2b_hat_vec=rep(0,BT)
  # ma_s2b_hat_vec=rep(0,BT)
  # 
  # for (bt in 1:BT)
  # {
  #   samp=sample(1:J,size=J,replace=TRUE)
  #   if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
  #   ## ma_s2b(2)
  #   mo_s2b_hat_vec[bt]=mean(zNJ[  (tau2_2+1):(N),samp]/zNJ[ (tau2_2):(N-1),samp ])
  #   N_J_b=( ( mo_s2b_hat_vec[bt]-1)/( mo_s2b_hat_vec[bt]^tau2_1*( mo_s2b_hat_vec[bt]^(tau2_2-tau2_1+1)-1)) )
  #   ma_s2b_hat_vec[bt]=sum( zNJ[(tau2_1):(tau2_2),samp] )/J * N_J_b
  #   if (is.nan(ma_s2b_hat_vec[bt])==TRUE)  {stop()}
  # } # end for bt
  # 
  # ## mo_s2b_hat=mean(mo_s2b_hat_vec)
  # ma_s2b_hat=mean(ma_s2b_hat_vec)
  # ma_s2b_vec[b]=ma_s2b_hat
  # ma_s2b_var_vec[b]=var(ma_s2b_hat_vec)
  # if ( ( quantile(ma_s2b_hat_vec,0.025)<m_A)&(quantile(ma_s2b_hat_vec,0.975)>m_A ) ) { cover_num_boot_s2=cover_num_boot_s2+1 }
  # cilength_boot_s2[b]=quantile(ma_s2b_hat_vec,0.975)-quantile(ma_s2b_hat_vec,0.025)

  #####

  ######## ma_e2tau2_e2
  #####

  mo_e2_hat= mean( zNJ[(tau2_1+1):(tau2_2),]/zNJ[tau2_1:(tau2_2-1),] )
  mo_e2_est_vec[b]=mo_e2_hat
  N_J=( ( mo_e2_hat-1)/( mo_e2_hat^tau2_2*( mo_e2_hat^(n-tau2_2+1)-1)) )
  ma_e2_hat= sum( zNJ[(tau2_2):(N),] )/J * N_J
  ma_e2_est_vec[b]=ma_e2_hat
  varma_e2_hat= 1/(J-1)* sum ( ( apply(zNJ[(tau2_2):(N) , ],2,sum) *N_J -ma_e2_hat )^2 )
  varma_e2_est_vec[b]=varma_e2_hat/J

  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varma_e2_hat)/sqrt(J)
  cilength_phi_e2[b]=2*boundl_phi
  if ( ( ma_e2_hat-boundl_phi <m_A)&( ma_e2_hat+boundl_phi > m_A) ) { cover_num_phi_e2=cover_num_phi_e2+1 }

  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varma_e2_hat)/sqrt(J)
  if ( ( ma_e2_hat-boundl_t <m_A)&( ma_e2_hat+boundl_t > m_A) ) { cover_num_t_e2=cover_num_t_e2+1 }
  cilength_t_e2[b]=2*boundl_t

  # ## bootstrap CI
  # BT=200
  # mo_e2b_hat_vec=rep(0,BT)
  # ma_e2b_hat_vec=rep(0,BT)
  # 
  # for (bt in 1:BT)
  # {
  #   samp=sample(1:J,size=J,replace=TRUE)
  #   if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
  #   ## ma_e2b(2)
  #   mo_e2b_hat_vec[bt]=mean(zNJ[  (tau2_1+1):(tau2_2),samp]/zNJ[ tau2_1:(tau2_2-1),samp ])
  #   N_J_b=( ( mo_e2b_hat_vec[bt]-1)/( mo_e2b_hat_vec[bt]^tau2_2*( mo_e2b_hat_vec[bt]^(n-tau2_2+1)-1)) )
  #   ma_e2b_hat_vec[bt]=sum( zNJ[(tau2_2):(N),samp] )/J * N_J_b
  #   if (is.nan(ma_e2b_hat_vec[bt])==TRUE)  {stop()}
  # } # end for bt
  # 
  # ## mo_e2b_hat=mean(mo_e2b_hat_vec)
  # ma_e2b_hat=mean(ma_e2b_hat_vec)
  # ma_e2b_vec[b]=ma_e2b_hat
  # ma_e2b_var_vec[b]=var(ma_e2b_hat_vec)
  # if ( ( quantile(ma_e2b_hat_vec,0.025)<m_A)&(quantile(ma_e2b_hat_vec,0.975)>m_A ) ) { cover_num_boot_e2=cover_num_boot_e2+1 }
  # cilength_boot_e2[b]=quantile(ma_e2b_hat_vec,0.975)-quantile(ma_e2b_hat_vec,0.025)

  #####


  ######## ma_s3tau1_s3
  #####

  mo_s3_hat= mean( zNJ[(tau3_2+1):(N),]/zNJ[tau3_2:(N-1),] )
  mo_s3_est_vec[b]=mo_s3_hat
  N_J=( ( mo_s3_hat-1)/( mo_s3_hat^tau3_1*( mo_s3_hat^(tau3_2-tau3_1+1)-1)) )
  ma_s3_hat= sum( zNJ[(tau3_1):(tau3_2),] )/J * N_J
  ma_s3_est_vec[b]=ma_s3_hat
  varma_s3_hat= 1/(J-1)* sum ( ( apply(zNJ[(tau3_1):(tau3_2) , ],2,sum) *N_J -ma_s3_hat )^2 )
  varma_s3_est_vec[b]=varma_s3_hat/J

  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varma_s3_hat)/sqrt(J)
  cilength_phi_s3[b]=2*boundl_phi
  if ( ( ma_s3_hat-boundl_phi <m_A)&( ma_s3_hat+boundl_phi > m_A) ) { cover_num_phi_s3=cover_num_phi_s3+1 }

  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varma_s3_hat)/sqrt(J)
  if ( ( ma_s3_hat-boundl_t <m_A)&( ma_s3_hat+boundl_t > m_A) ) { cover_num_t_s3=cover_num_t_s3+1 }
  cilength_t_s3[b]=2*boundl_t

  # ## bootstrap CI
  # BT=200
  # mo_s3b_hat_vec=rep(0,BT)
  # ma_s3b_hat_vec=rep(0,BT)
  # 
  # for (bt in 1:BT)
  # {
  #   samp=sample(1:J,size=J,replace=TRUE)
  #   if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
  #   ## ma_s3b(2)
  #   mo_s3b_hat_vec[bt]=mean(zNJ[  (tau3_2+1):(N),samp]/zNJ[ (tau3_2):(N-1),samp ])
  #   N_J_b=( ( mo_s3b_hat_vec[bt]-1)/( mo_s3b_hat_vec[bt]^tau3_1*( mo_s3b_hat_vec[bt]^(tau3_2-tau3_1+1)-1)) )
  #   ma_s3b_hat_vec[bt]=sum( zNJ[(tau3_1):(tau3_2),samp] )/J * N_J_b
  #   if (is.nan(ma_s3b_hat_vec[bt])==TRUE)  {stop()}
  # } # end for bt
  # 
  # ## mo_s3b_hat=mean(mo_s3b_hat_vec)
  # ma_s3b_hat=mean(ma_s3b_hat_vec)
  # ma_s3b_vec[b]=ma_s3b_hat
  # ma_s3b_var_vec[b]=var(ma_s3b_hat_vec)
  # if ( ( quantile(ma_s3b_hat_vec,0.025)<m_A)&(quantile(ma_s3b_hat_vec,0.975)>m_A ) ) { cover_num_boot_s3=cover_num_boot_s3+1 }
  # cilength_boot_s3[b]=quantile(ma_s3b_hat_vec,0.975)-quantile(ma_s3b_hat_vec,0.025)

  #####

  ######## ma_e3tau1_e3
  #####

  mo_e3_hat= mean( zNJ[(tau3_1+1):(tau3_2),]/zNJ[tau3_1:(tau3_2-1),] )
  mo_e3_est_vec[b]=mo_e3_hat
  N_J=( ( mo_e3_hat-1)/( mo_e3_hat^tau3_2*( mo_e3_hat^(n-tau3_2+1)-1)) )
  ma_e3_hat= sum( zNJ[(tau3_2):(N),] )/J * N_J
  ma_e3_est_vec[b]=ma_e3_hat
  varma_e3_hat= 1/(J-1)* sum ( ( apply(zNJ[(tau3_2):(N) , ],2,sum) *N_J -ma_e3_hat )^2 )
  varma_e3_est_vec[b]=varma_e3_hat/J

  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varma_e3_hat)/sqrt(J)
  cilength_phi_e3[b]=2*boundl_phi
  if ( ( ma_e3_hat-boundl_phi <m_A)&( ma_e3_hat+boundl_phi > m_A) ) { cover_num_phi_e3=cover_num_phi_e3+1 }

  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varma_e3_hat)/sqrt(J)
  if ( ( ma_e3_hat-boundl_t <m_A)&( ma_e3_hat+boundl_t > m_A) ) { cover_num_t_e3=cover_num_t_e3+1 }
  cilength_t_e3[b]=2*boundl_t

  # ## bootstrap CI
  # BT=200
  # mo_e3b_hat_vec=rep(0,BT)
  # ma_e3b_hat_vec=rep(0,BT)
  # 
  # for (bt in 1:BT)
  # {
  #   samp=sample(1:J,size=J,replace=TRUE)
  #   if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
  #   ## ma_e3b(2)
  #   mo_e3b_hat_vec[bt]=mean(zNJ[  (tau3_1+1):(tau3_2),samp]/zNJ[ tau3_1:(tau3_2-1),samp ])
  #   N_J_b=( ( mo_e3b_hat_vec[bt]-1)/( mo_e3b_hat_vec[bt]^tau3_2*( mo_e3b_hat_vec[bt]^(n-tau3_2+1)-1)) )
  #   ma_e3b_hat_vec[bt]=sum( zNJ[(tau3_2):(N),samp] )/J * N_J_b
  #   if (is.nan(ma_e3b_hat_vec[bt])==TRUE)  {stop()}
  # } # end for bt
  # 
  # ## mo_e3b_hat=mean(mo_e3b_hat_vec)
  # ma_e3b_hat=mean(ma_e3b_hat_vec)
  # ma_e3b_vec[b]=ma_e3b_hat
  # ma_e3b_var_vec[b]=var(ma_e3b_hat_vec)
  # if ( ( quantile(ma_e3b_hat_vec,0.025)<m_A)&(quantile(ma_e3b_hat_vec,0.975)>m_A ) ) { cover_num_boot_e3=cover_num_boot_e3+1 }
  # cilength_boot_e3[b]=quantile(ma_e3b_hat_vec,0.975)-quantile(ma_e3b_hat_vec,0.025)

  #####
  
  
  ###### process counter
  #####
  if ( b%%(B/10)==1 ) 
  {cat(c( (b-1)/(B/10))*10 ) 
    cat("% \n")}
  #####
  
} # end loop for b
  
  ######## result table 
  #####
  if (jj==1)
  {
    output_tab_s1=c(mean(ma_s1_est_vec),mean(varma_s1_est_vec),cover_num_phi_s1/B,cover_num_t_s1/B)
    output_tab_e1=c(mean(ma_e1_est_vec),mean(varma_e1_est_vec),cover_num_phi_e1/B,cover_num_t_e1/B)
    output_tab_s2=c(mean(ma_s2_est_vec),mean(varma_s2_est_vec),cover_num_phi_s2/B,cover_num_t_s2/B)
    output_tab_e2=c(mean(ma_e2_est_vec),mean(varma_e2_est_vec),cover_num_phi_e2/B,cover_num_t_e2/B)
    output_tab_s3=c(mean(ma_s3_est_vec),mean(varma_s3_est_vec),cover_num_phi_s3/B,cover_num_t_s3/B)
    output_tab_e3=c(mean(ma_e3_est_vec),mean(varma_e3_est_vec),cover_num_phi_e3/B,cover_num_t_e3/B)
  }
  else 
  {
    output_tab_s1=rbind(output_tab_s1,c(mean(ma_s1_est_vec),mean(varma_s1_est_vec),cover_num_phi_s1/B,cover_num_t_s1/B))
    output_tab_e1=rbind(output_tab_e1,c(mean(ma_e1_est_vec),mean(varma_e1_est_vec),cover_num_phi_e1/B,cover_num_t_e1/B))
    output_tab_s2=rbind(output_tab_s2,c(mean(ma_s2_est_vec),mean(varma_s2_est_vec),cover_num_phi_s2/B,cover_num_t_s2/B))
    output_tab_e2=rbind(output_tab_e2,c(mean(ma_e2_est_vec),mean(varma_e2_est_vec),cover_num_phi_e2/B,cover_num_t_e2/B))
    output_tab_s3=rbind(output_tab_s3,c(mean(ma_s3_est_vec),mean(varma_s3_est_vec),cover_num_phi_s3/B,cover_num_t_s3/B))
    output_tab_e3=rbind(output_tab_e3,c(mean(ma_e3_est_vec),mean(varma_e3_est_vec),cover_num_phi_e3/B,cover_num_t_e3/B))
    
  }
  #####
  
}


#####


# mean( ma_s1_est_vec )
# mean(varma_s1_est_vec)
# mean( ma_e1_est_vec )
# mean(varma_e1_est_vec)
# 
# 
# mean( ma_s2_est_vec )
# mean(varma_s2_est_vec)
# mean( ma_e2_est_vec )
# mean(varma_e2_est_vec)
# 
# 
# mean( ma_s3_est_vec )
# mean(varma_s3_est_vec)
# mean( ma_e3_est_vec )
# mean(varma_e3_est_vec)

#####






output_tab_s1
output_tab_e1
output_tab_s2
output_tab_e2
output_tab_s3
output_tab_e3





