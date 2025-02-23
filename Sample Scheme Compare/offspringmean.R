# This is a BPRE(i.i.d) simulation with a poisson offspring distribution 
# lambda is determined by gamma distribution

################# data size setting
#####
set.seed(200)

B=500

## test Setting

J=40

N=20

tau1=5

tau2=12

tau3=18


# poisson ancestor
m_A=10
sigma2_A=m_A




b=1
ii=1
jj=1


## Table Setting
mavec=c(10,20,30,40)

# Jvec=c(5,6,8,10)

Jvec=c(20,30,40,50)

# Nvec=c(10,20,30,40)

mafinalvec=matrix(0,4,4)

mavarfinalvec=matrix(0,4,4)

testmatrix=matrix(0,nrow=5,5)
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
NJ1_est_vec=rep(0,B)
varmo1_est_vec=rep(0,B)
mo1b_vec=rep(0,B)
mo1b_var_vec=rep(0,B)
cover_num_phi1=0
cilength_phi1=rep(0,B)
cover_num_t1=0
cilength_t1=rep(0,B)
cover_num_boot1=0
cilength_boot1=rep(0,B)


mo2_est_vec=rep(0,B)
NJ2_est_vec=rep(0,B)
varmo2_est_vec=rep(0,B)
mo2b_vec=rep(0,B)
mo2b_var_vec=rep(0,B)
cover_num_phi2=0
cilength_phi2=rep(0,B)
cover_num_t2=0
cilength_t2=rep(0,B)
cover_num_boot2=0
cilength_boot2=rep(0,B)


mo3_est_vec=rep(0,B)
NJ3_est_vec=rep(0,B)
varmo3_est_vec=rep(0,B)
mo3b_vec=rep(0,B)
mo3b_var_vec=rep(0,B)
cover_num_phi3=0
cilength_phi3=rep(0,B)
cover_num_t3=0
cilength_t3=rep(0,B)
cover_num_boot3=0
cilength_boot3=rep(0,B)

#####


################# Simulation: B loop
#####
# z0J=matrix(0,B,J)
zNJ=matrix(0,N,J)

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
  
  ######## mo1tau1
  ##### 
  
  mo1_hat= mean( zNJ[(tau1+1):(N),]/zNJ[tau1:(N-1),] )
  mo1_est_vec[b]=mo1_hat
  N_J=( ( mo1_hat-1)/( mo1_hat^tau1*( mo1_hat^(N-tau1+1)-1)) )
  varmo1_hat= mean (( zNJ[(tau1+1):(N),]/zNJ[tau1:(N-1),] - mo1_hat )^2)  
  varmo1_est_vec[b]=varmo1_hat /(N-tau1-1)
  
  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varmo1_hat)/sqrt(J)
  cilength_phi1[b]=2*boundl_phi
  if ( ( mo1_hat-boundl_phi <mo)&( mo1_hat+boundl_phi > mo) ) { cover_num_phi1=cover_num_phi1+1 }
  
  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varmo1_hat)/sqrt(J)
  if ( ( mo1_hat-boundl_t <mo)&( mo1_hat+boundl_t > mo) ) { cover_num_t1=cover_num_t1+1 }
  cilength_t1[b]=2*boundl_t
  
  ## bootstrap CI 
  BT=200
  mo1b_hat_vec=rep(0,BT)
  mo1b_hat_vec=rep(0,BT)
  
  for (bt in 1:BT)
  {  
    samp=sample(1:J,size=J,replace=TRUE)
    if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
    ## mo1b(2)
    mo1b_hat_vec[bt]=mean(zNJ[  (tau1+1):(N),samp]/zNJ[ (tau1):(N-1),samp ])
    if (is.nan(mo1b_hat_vec[bt])==TRUE)  {stop()}
  } # end for bt
  
  mo1b_hat=mean(mo1b_hat_vec)
  mo1b_vec[b]=mo1b_hat
  mo1b_var_vec[b]=var(mo1b_hat_vec)
  if ( ( quantile(mo1b_hat_vec,0.025)<mo)&(quantile(mo1b_hat_vec,0.975)>mo ) ) { cover_num_boot1=cover_num_boot1+1 }
  cilength_boot1[b]=quantile(mo1b_hat_vec,0.975)-quantile(mo1b_hat_vec,0.025)
  
  #####
  
  
  ######## mo2tau2
  ##### 
  
  mo2_hat= mean( zNJ[(tau2+1):(N),]/zNJ[tau2:(N-1),] )
  mo2_est_vec[b]=mo2_hat
  N_J=( ( mo2_hat-1)/( mo2_hat^tau2*( mo2_hat^(N-tau2+1)-1)) )
  varmo2_hat= mean (( zNJ[(tau2+1):(N),]/zNJ[tau2:(N-1),] - mo2_hat )^2)  
  varmo2_est_vec[b]=varmo2_hat /(N-tau2-1)
  
  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varmo2_hat)/sqrt(J)
  cilength_phi2[b]=2*boundl_phi
  if ( ( mo2_hat-boundl_phi <mo)&( mo2_hat+boundl_phi > mo) ) { cover_num_phi2=cover_num_phi2+1 }
  
  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varmo2_hat)/sqrt(J)
  if ( ( mo2_hat-boundl_t <mo)&( mo2_hat+boundl_t > mo) ) { cover_num_t2=cover_num_t2+1 }
  cilength_t2[b]=2*boundl_t
  
  ## bootstrap CI 
  BT=200
  mo2b_hat_vec=rep(0,BT)
  mo2b_hat_vec=rep(0,BT)
  
  for (bt in 1:BT)
  {  
    samp=sample(1:J,size=J,replace=TRUE)
    if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
    ## mo2b(2)
    mo2b_hat_vec[bt]=mean(zNJ[  (tau2+1):(N),samp]/zNJ[ (tau2):(N-1),samp ])
    if (is.nan(mo2b_hat_vec[bt])==TRUE)  {stop()}
  } # end for bt
  
  mo2b_hat=mean(mo2b_hat_vec)
  mo2b_vec[b]=mo2b_hat
  mo2b_var_vec[b]=var(mo2b_hat_vec)
  if ( ( quantile(mo2b_hat_vec,0.025)<mo)&(quantile(mo2b_hat_vec,0.975)>mo ) ) { cover_num_boot2=cover_num_boot2+1 }
  cilength_boot2[b]=quantile(mo2b_hat_vec,0.975)-quantile(mo2b_hat_vec,0.025)
  
  #####
  
  
  ######## mo3tau3
  ##### 
  
  mo3_hat= mean( zNJ[(tau3+1):(N),]/zNJ[tau3:(N-1),] )
  mo3_est_vec[b]=mo3_hat
  N_J=( ( mo3_hat-1)/( mo3_hat^tau3*( mo3_hat^(N-tau3+1)-1)) )
  varmo3_hat= mean (( zNJ[(tau3+1):(N),]/zNJ[tau3:(N-1),] - mo3_hat )^2)  
  varmo3_est_vec[b]=varmo3_hat /(N-tau3-1)
  
  ## phi CI
  boundl_phi= qnorm(0.975) * sqrt(varmo3_hat)/sqrt(J)
  cilength_phi3[b]=2*boundl_phi
  if ( ( mo3_hat-boundl_phi <mo)&( mo3_hat+boundl_phi > mo) ) { cover_num_phi3=cover_num_phi3+1 }
  
  ## t CI
  boundl_t=qt(0.975,J-1)* sqrt(varmo3_hat)/sqrt(J)
  if ( ( mo3_hat-boundl_t <mo)&( mo3_hat+boundl_t > mo) ) { cover_num_t3=cover_num_t3+1 }
  cilength_t3[b]=2*boundl_t
  
  ## bootstrap CI 
  BT=200
  mo3b_hat_vec=rep(0,BT)
  mo3b_hat_vec=rep(0,BT)
  
  for (bt in 1:BT)
  {  
    samp=sample(1:J,size=J,replace=TRUE)
    if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
    ## mo3b(2)
    mo3b_hat_vec[bt]=mean(zNJ[  (tau3+1):(N),samp]/zNJ[ (tau3):(N-1),samp ])
    if (is.nan(mo3b_hat_vec[bt])==TRUE)  {stop()}
  } # end for bt
  
  mo3b_hat=mean(mo3b_hat_vec)
  mo3b_vec[b]=mo3b_hat
  mo3b_var_vec[b]=var(mo3b_hat_vec)
  if ( ( quantile(mo3b_hat_vec,0.025)<mo)&(quantile(mo3b_hat_vec,0.975)>mo ) ) { cover_num_boot3=cover_num_boot3+1 }
  cilength_boot3[b]=quantile(mo3b_hat_vec,0.975)-quantile(mo3b_hat_vec,0.025)
  
  #####
  
  
  ## process counter
  if ( b%%(B/10)==1 ) 
  {cat(c( (b-1)/(B/10))*10 ) 
    cat("% \n")}
} # end loop for b



mean( mo1_est_vec )
mean(varmo1_est_vec)

mean( mo2_est_vec )
mean(varmo2_est_vec)

mean( mo3_est_vec )
mean(varmo3_est_vec)








