
# This is a BPRE(i.i.d) simulation with a poisson offspring distribution 
# lambda is determined by gamma distribution
set.seed(100)

### Simulation Setting
B=5000
## test Setting

J=20
N=30


tau=N*0.5
m_AT=3000
m_A=m_AT
m_AC=100

R=m_AT/m_AC


sigma2_AT=m_AT

sigma2_A=sigma2_AT

sigma2_AC=m_AC

b=1



## true mean EX= m_A*e^(m_A)/(e^(m_A)-1)  
## VarX=EX*(1+m_A-EX) 




result_matrix=matrix(0,nrow=10,ncol=4)




## Table Setting
mavec=c(10,20,30,40)

Jvec=c(20,30,40,50)

Nvec=2^(5:8)



# for (ii in 1:4)
# {
for (jj in 1:4)
  {

# Offspring Distribution Setting
J=Jvec[jj]
# N=Nvec[jj]
# m_A=mavec[ii]




#############
##         gamma poisson Target
# alpha=10
# beta=0.03
# mo=1+alpha*beta
# sigma2=alpha*beta^2
# mo2=sigma2+mo^2
# r_o= sqrt(mo2/(mo^2)  )
# gamma2=alpha*beta
# 

################# beta binomial  Target


alpha=90
beta=10
mo=1+alpha/(alpha+beta)
sigma2=alpha*beta/ ( (alpha+beta)^2*(alpha+beta+1) )
mo2=sigma2+mo^2
r_o= sqrt(mo2/(mo^2)  )
gamma2= (alpha*beta)/((alpha+beta)*(alpha+beta+1) )

###
###################

############# 
##         gamma poisson    Calibrator
# alphaC=10
# betaC=0.03
# moC=1+alpha*beta
# sigma2C=alpha*beta^2
# mo2C=sigma2+mo^2
# r_oC= sqrt(mo2/(mo^2)  )
# gamma2C=alpha*beta
# 

################# beta binomial  Calibrator


alphaC=90
betaC=10
moC=1+alpha/(alpha+beta)
sigma2C=alpha*beta/ ( (alpha+beta)^2*(alpha+beta+1) )
mo2C=sigma2+mo^2
r_oC= sqrt(mo2/(mo^2)  )
gamma2C= (alpha*beta)/((alpha+beta)*(alpha+beta+1) )

###
################### 


eta=(mo^2)/mo2
kap=mo/mo2

etaC=(moC^2)/mo2C
kapC=moC/mo2C


# 
# mfD=m_A*gamma2*((1)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1)/(1-eta))
# mfDn=m_A*gamma2*((1-kap^N)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1-eta^N)/(1-eta))
# 
# varzn=mo2^(N-1)*mfD+mo^(2*N)*sigma2_A
# 
# sigma2_A/r_o^(2*N)
# 
# mfDn/mo2
# 
# mfD
# 
# testmatrix[,5]=c(mo,mo2,r_o,m_A,mfD)
# 

#setting for estimator
moT_est_vec=rep(0,B)

maT_est_vec=rep(0,B)

moC_est_vec=rep(0,B)

maC_est_vec=rep(0,B)

R_est_vec=rep(0,B)

#setting for variance
varmaT_est_vec=rep(0,B)

varmaC_est_vec=rep(0,B)

#setting for 


Rb_var_vec=rep(0,B)

varR_est_vec=rep(0,B)

#setting for CI length and Coverage rate
cover_num_phi=0
cilength_phi=rep(0,B)

cover_num_t=0
cilength_t=rep(0,B)

cover_num_boot=0
cilength_boot=rep(0,B)


for (b in (1:B))
{
  
  zNJT=matrix(0,N,J)
  zNJC=matrix(0,N,J)
  
  
  
  for (j in (1:J) )
  {
    
    
    ## Generate  poisson - gamma  distribution  Target
    # test1T=rgamma(N,shape=alpha,scale=beta)   # poisson - gamma
    # repeat
    # {
    # z0T=rpois(1,m_AT)
    # if (z0T>0){break}
    # }
    # # Generate Generations
    #  zNJT[1,j]=z0+rpois(1,z0T*test1T[1])     
    # for (n in 2:(N))
    # {
    #    zNJT[n,j]=zNJT[n,j]+rpois(1,zNJT[n-1,j]*test1T[n])     
    # }
    
    
    ## Generate binomial - beta distribution Target
    test1T=rbeta(N,shape1=alpha,shape2=beta)   #  binomial - beta
    
    repeat
    {
      z0T=rpois(1,m_AT)
      if (z0T>0){break}
    } 

    # Generate Generations
    zNJT[1,j]=z0T+rbinom(1,z0T,test1T[1])   # 
    for (n in 2:N)
    {
      zNJT[n,j]=zNJT[n-1,j]+rbinom(1,zNJT[n-1,j],test1T[n])  
    }  

    ## Generate  poisson - gamma  distribution Calibrator
    # test1C=rgamma(N,shape=alpha,scale=beta)   # poisson - gamma
    # repeat
    # {
    # z0C=rpois(1,m_AC)
    # if (z0C>0){break}
    # }
    # # Generate Generations
    #  zNJC[1,j]=z0+rpois(1,z0C*test1C[1])     
    # for (n in 2:(N))
    # {
    #    zNJC[n,j]=zNJC[n,j]+rpois(1,zNJC[n-1,j]*test1C[n])     
    # }
    
    
    ## Generate binomial - beta distribution Calibrator
    test1C=rbeta(N,shape1=alpha,shape2=beta)
    
    repeat
    {
      z0C=rpois(1,m_AC)
      if (z0C>0){break}
    } 
    
    zNJC[1,j]=z0C+rbinom(1,z0C,test1C[1])   # 
    for (n in 2:N)
    {
      zNJC[n,j]=zNJC[n-1,j]+rbinom(1,zNJC[n-1,j],test1C[n])  
    }  
  }
  
  
  
  ## ma(2)
  moT_hat= mean( zNJT[(tau+1):(N),]/zNJT[tau:(N-1),] )
  
  moC_hat= mean( zNJC[(tau+1):(N),]/zNJC[tau:(N-1),] )
  
  moT_est_vec[b]=moT_hat
  
  moC_est_vec[b]=moC_hat
  
  N_JT=( ( moT_hat-1)/( moT_hat^tau*( moT_hat^(N-tau+1)-1)) )
  
  maT_hat= sum( zNJT[(tau):(N),] )/J * N_JT
  
  N_JC=( ( moC_hat-1)/( moC_hat^tau*( moC_hat^(N-tau+1)-1)) )
  
  maC_hat= sum( zNJC[(tau):(N),] )/J * N_JC
  
  maT_est_vec[b]=maT_hat
  
  maC_est_vec[b]=maC_hat
  
  varmaT_hat=1/(J-1)* sum ( ( apply(zNJT[tau:N , ],2,sum) *N_JT -maT_hat )^2 )
  
  varmaC_hat=1/(J-1)* sum ( ( apply(zNJC[tau:N , ],2,sum) *N_JC -maC_hat )^2 )
  
  
  
  # ma(1)
  # moT_hat= mean( zNJT[2:(N-1),]/zNJT[1:(N-2),] )
  # 
  # moC_hat= mean( zNJC[2:(N-1),]/zNJC[1:(N-2,] )
  # 
  # moT_est_vec[b]=moT_hat 
  # 
  # moC_est_vec[b]=moC_hat
  # 
  # maT_hat= mean( zNJT[N,] / moT_hat^(N) )
  # 
  # maC_hat= mean( zNJC[N,] / moC_hat^(N) )
  # 
  # maT_est_vec[b]=maT_hat
  # 
  # maC_est_vec[b]=maC_hat
  #  
  # varmaT_hat=mean( (zNJT[N,]/moT_hat^N-maT_hat)^2     )*(J/(J-1))
  # 
  # varmaC_hat=mean( (zNJC[N,]/moC_hat^N-maC_hat)^2     )*(J/(J-1))
  
  


  # var R
  
  R_hat=maT_hat/maC_hat
  
  R_est_vec[b]=R_hat
  
  varmaT_est_vec[b]=varmaT_hat 
  
  varmaC_est_vec[b]=varmaC_hat 
  
  varR_est= R_hat^2*(  varmaT_hat/maT_hat^2   +  varmaC_hat/maC_hat^2  )
  
  varR_est_vec[b]=varR_est/J
  
  ################# Inference
  ##### phi CI
  
  boundl_phi= qnorm(0.975) * sqrt(varR_est)/sqrt(J)
  
  cilength_phi[b]=2*boundl_phi
  
  if ( ( R_hat-boundl_phi <R)&( R_hat+boundl_phi > R) ) { cover_num_phi=cover_num_phi+1 }
  
  ##### t CI
  
  boundl_t=qt(0.975,J-1)* sqrt(varR_est)/sqrt(J)
  
  if ( ( R_hat-boundl_t <R)&( R_hat+boundl_t > R) ) { cover_num_t=cover_num_t+1 }
  
  cilength_t[b]=2*boundl_t
  
  
  ##### bootstrap CI 
  
  BT=200
  mobT_hat_vec=rep(0,BT)
  mabT_hat_vec=rep(0,BT)
  mobC_hat_vec=rep(0,BT)
  mabC_hat_vec=rep(0,BT)
  

  
  for (bt in 1:BT)
  {
    samp=sample(1:J,size=J,replace=TRUE)
    
    ## mab(2)
    
    mobT_hat_vec[bt]=mean(zNJT[  (tau+1):(N),samp]/zNJT[ (tau):(N-1),samp ])
    
    N_JT_b=( ( mobT_hat_vec[bt]-1)/( mobT_hat_vec[bt]^tau*( mobT_hat_vec[bt]^(N-tau+1)-1)) )
    
    mabT_hat_vec[bt]=sum( zNJT[(tau):(N),samp] )/J * N_JT_b
    
    mobC_hat_vec[bt]=mean(zNJC[  (tau+1):(N),samp]/zNJC[ (tau):(N-1),samp ])
    
    N_JC_b=( ( mobC_hat_vec[bt]-1)/( mobC_hat_vec[bt]^tau*( mobC_hat_vec[bt]^(N-tau+1)-1)) )
    
    mabC_hat_vec[bt]=sum( zNJC[(tau):(N),samp] )/J * N_JC_b
    
    
    ## mab(1)
    
    # mobT_hat_vec[bt]=mean(zNJT[  2:(N-1),samp]/zNJT[ 1:(N-2),samp ])
    # 
    # mabT_hat_vec[bt]=mean(zNJT[N,samp ] )/ mobT_hat_vec[r]^N
    # 
    # mobC_hat_vec[bt]=mean(zNJC[  2:(N-1),samp]/zNJC[ 1:(N-2),samp ])
    # 
    # mabC_hat_vec[bt]=mean(zNJC[N,samp ] )/ mobC_hat_vec[r]^N

  }

  

  Rb_hat_vec= mabT_hat_vec/mabC_hat_vec
  
  Rb_var_vec[b]=var(Rb_hat_vec)
  
  if ( ( quantile(Rb_hat_vec,0.025)<R)&(quantile(Rb_hat_vec,0.975)>R ) ) { cover_num_boot=cover_num_boot+1 }
  
  cilength_boot[b]=quantile(Rb_hat_vec,0.975)-quantile(Rb_hat_vec,0.025)
  
  ## process counter
  if ( b%%(B/10)==1 ) 
    {cat(c( (b-1)/(B/10))*10 ) 
    cat("% \n")}
} # end loop for b




result_matrix[1,jj]=mean( R_est_vec )

result_matrix[2,jj]=mean(varR_est_vec)

cover_rate_phi= cover_num_phi /B

result_matrix[3,jj]=cover_rate_phi

result_matrix[4,jj]=mean( cilength_phi)

cover_rate_t= cover_num_t /B

result_matrix[5,jj]=cover_rate_t

result_matrix[6,jj]=mean( cilength_t)



result_matrix[7,jj]=mean(Rb_hat_vec)

result_matrix[8,jj]=mean( Rb_var_vec)

cover_rate_boot= cover_num_boot /B

result_matrix[9,jj]=cover_rate_boot

result_matrix[10,jj]=mean( cilength_boot)


} # end loop for jj



### theoretical asympototic variance

mfD=m_A*gamma2*((1)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1)/(1-eta))

mfDn=m_A*gamma2*((1-kap^N)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1-eta^N)/(1-eta))

Ntau=c(tau:N)

mfDtau=m_A*gamma2*((1-kap^Ntau)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1-eta^Ntau)/(1-eta))



mfDC=m_AC*gamma2C*((1)/(1-kapC))+ (m_AC^2+sigma2_AC)*sigma2C*((1)/(1-etaC))

mfDnC=m_AC*gamma2C*((1-kapC^N)/(1-kapC))+ (m_AC^2+sigma2_AC)*sigma2C*((1-etaC^N)/(1-etaC))

Ntau=c(tau:N)

mfDtauC=m_AC*gamma2C*((1-kapC^Ntau)/(1-kapC))+ (m_AC^2+sigma2_AC)*sigma2C*((1-etaC^Ntau)/(1-etaC))




## with sigma_A 

# ma(1)

# var_the_whole=(mfDn/mo2*r_o^(2*N)+sigma2_A) /Jvec
# 
# var_the_wholeC=(mfDnC/mo2C*r_oC^(2*N)+sigma2_AC) /Jvec
# 


## ma(2)

N_J_t=( ( mo-1)/( mo^tau*( mo^(N-tau+1)-1)) )

tauvec=tau:N

var_znvec=mfDtau*mo2^(tauvec-1)+mo^(2*tauvec)*sigma2_A

var_zntau=0

for (l in tau:(N-1) )
{
  var_zntau= var_zntau +( 1+2*  (sum( mo^(1:(N-l)) ))  ) *var_znvec[l-tau+1]
 # print(var_zntau)
}
var_zntau=var_zntau+ var_znvec[N-tau+1]

var_the_whole=N_J_t^2*var_zntau/Jvec




N_J_tC=( ( moC-1)/( moC^tau*( moC^(N-tau+1)-1)) )

tauvec=tau:N

var_znvecC=mfDtauC*mo2C^(tauvec-1)+moC^(2*tauvec)*sigma2_AC

var_zntauC=0

for (l in tau:(N-1) )
{
  var_zntauC= var_zntauC +( 1+2*  (sum( moC^(1:(N-l)) ))  ) *var_znvecC[l-tau+1]
  # print(var_zntauC)
}
var_zntauC=var_zntauC + var_znvecC[N-tau+1]

var_the_wholeC=N_J_tC^2*var_zntauC/Jvec






var_the_wholeR=R^2*(  var_the_whole/m_AT^2   +  var_the_wholeC/m_AC^2  )

var_the_asy=var_the_wholeR



result_matrix=round(result_matrix,3)


#without asy
# table_mat=result_matrix
# 
# table_name=c("J", " $\\hat{R}$" ,"$\\hat{ \\Lambda}_{R,2,n,J}$ ", "G CR" ,"G ML" ,"t CR" ,"t ML" ,"$\\hat{R}_b$" ,"Var$_b$","b CR","b ML")

##with asy

table_mat=matrix(0,nrow=nrow(result_matrix)+1 ,ncol=ncol(result_matrix))

table_mat[1:2,]=result_matrix[1:2,]

table_mat[3,]=round(var_the_asy,3)

table_mat[4:nrow(table_mat), ]=result_matrix[3:nrow(result_matrix),  ]

table_name=c("J", " $\\hat{m}_A$" ,"$\\hat{ \\Lambda}_{1,n,J}$ ", "Asy Var " ,"G CR" ,"G ML" ,"t CR" ,"t ML" ,"$\\hat{m}_{A,b}$" ,"Var$_b$","b CR","b ML")






for (jj in 1:(dim(table_mat)[1]+1) )
{
  cat(table_name[jj])
  for (ii in 1:4)
  {
    if (jj==1)
    {
      cat(" & ")
      cat(Jvec[ii])
    }
    else
    {
      cat(" & ")
      cat(table_mat[jj-1,ii])   
    }
  }
  cat(" \\\\ \n")
  cat("\\hline \n")
}

c(R,N ,tau  )
