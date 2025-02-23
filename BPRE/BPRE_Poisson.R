
# This is a BPRE(i.i.d) simulation with a poisson offspring distribution 
# lambda is determined by gamma distribution
set.seed(100)

################## Simulation Setting
B=5000

## test Setting


J=10

N=20

tau=12


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



result_matrix=matrix(0,nrow=10,ncol=4)

result_matrix_2=matrix(0,nrow=6,ncol=4)



# for (ii in 1:4)
# {
for (jj in 1:4)
  {
  
# Offspring Distribution Setting
J=Jvec[jj]
# N=Nvec[jj]
# m_A=mavec[ii]




#############
## gamma poisson
# alpha=10
# beta=.03
# mo=1+alpha*beta
# sigma2=alpha*beta^2
# mo2=sigma2+mo^2
# r_o= sqrt(mo2/(mo^2)  )
# gamma2=alpha*beta


################# 
##   beta binomial
alpha=90
beta=10
mo=1+alpha/(alpha+beta)
sigma2=alpha*beta/ ( (alpha+beta)^2*(alpha+beta+1) )
mo2=sigma2+mo^2
r_o= sqrt(mo2/(mo^2)  )
gamma2= (alpha*beta)/((alpha+beta)*(alpha+beta+1) )

###
###################


eta=(mo^2)/mo2
kap=mo/mo2

mfD=m_A*gamma2*((1)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1)/(1-eta))
mfDn=m_A*gamma2*((1-kap^N)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1-eta^N)/(1-eta))



################# Simulation part

mo_est_vec=rep(0,B)

ma_est_vec=rep(0,B)

varma_est_vec=rep(0,B)

# mo2_est_vec=rep(0,B)

# ro_est_vec=rep(0,B)

# mfD_est_vec=rep(0,B)

mab_vec=rep(0,B)
mab_var_vec=rep(0,B)


zNJ=matrix(0,N,J)

cover_num_phi=0
cilength_phi=rep(0,B)

cover_num_t=0
cilength_t=rep(0,B)

cover_num_boot=0
cilength_boot=rep(0,B)

z0J=matrix(0,B,J)

for (b in (1:B))
{

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
    
       

      z0J[b,j]=z0
    
    }
  
   ###### store all estimators for simulations
  
   
   #### ma(2)
   
   mo_hat= mean( zNJ[(tau+1):(N),]/zNJ[tau:(N-1),] )
   mo_est_vec[b]=mo_hat

   N_J=( ( mo_hat-1)/( mo_hat^tau*( mo_hat^(N-tau+1)-1)) )

   ma_hat= sum( zNJ[(tau):(N),] )/J * N_J
   ma_est_vec[b]=ma_hat

   varma_hat= 1/(J-1)* sum ( ( apply(zNJ[tau:N , ],2,sum) *N_J -ma_hat )^2 )
   varma_est_vec[b]=varma_hat/J
   
   ## ma(1)
   
   # mo_hat= mean( zNJ[2:(N-1),]/zNJ[1:(N-2),] ) ## calculate m_{o,n,J}
   # mo_est_vec[b]=mo_hat
   # 
   # ma_hat= mean (zNJ[N,]) / ( mo_hat^N   )   ## calculate m_{a,n,J} (1)
   # ma_est_vec[b]=ma_hat
   # 
   # varma_hat=mean( (zNJ[N,]/mo_hat^N-ma_hat)^2     )*(J/(J-1))
   # varma_est_vec[b]=varma_hat/J
   

   
   ################# Inference
   ##### phi CI
   
   boundl_phi= qnorm(0.975) * sqrt(varma_hat)/sqrt(J)
   
   cilength_phi[b]=2*boundl_phi
   
   if ( ( ma_hat-boundl_phi <m_A)&( ma_hat+boundl_phi > m_A) ) { cover_num_phi=cover_num_phi+1 }
   
   ##### t CI
   
   boundl_t=qt(0.975,J-1)* sqrt(varma_hat)/sqrt(J)
   
   if ( ( ma_hat-boundl_t <m_A)&( ma_hat+boundl_t > m_A) ) { cover_num_t=cover_num_t+1 }
   
   cilength_t[b]=2*boundl_t
   
   ##### bootstrap CI 
   
   BT=200
   mob_hat_vec=rep(0,BT)
   mab_hat_vec=rep(0,BT)
   
 
   
   for (bt in 1:BT)
   {  
     samp=sample(1:J,size=J,replace=TRUE)
     
     if (min(samp)==max(samp)) {samp[J]=sample( (1:J)[-(samp[1])],size=1)} ## avoid same value for all replicates
     
     ## mab(2)
     
     mob_hat_vec[bt]=mean(zNJ[  (tau+1):(N),samp]/zNJ[ (tau):(N-1),samp ])

     N_J_b=( ( mob_hat_vec[bt]-1)/( mob_hat_vec[bt]^tau*( mob_hat_vec[bt]^(N-tau+1)-1)) )
      
     mab_hat_vec[bt]=sum( zNJ[(tau):(N),samp] )/J * N_J_b
     
     if (is.nan(mab_hat_vec[bt])==TRUE)  {stop()}
     
     ## mab(1)
     
     # mob_hat_vec[bt]=mean(zNJ[  2:(N-1),samp]/zNJ[ 1:(N-2),samp ])
     # 
     # mab_hat_vec[bt]=mean(zNJ[N,samp ] )/ mob_hat_vec[bt]^N

   } # end for bt
   
   
   # mob_hat=mean(mob_hat_vec)
   
   mab_hat=mean(mab_hat_vec)
   
   mab_vec[b]=mab_hat
   
   mab_var_vec[b]=var(mab_hat_vec)
   
   mab_vec[b]
   
   if ( ( quantile(mab_hat_vec,0.025)<m_A)&(quantile(mab_hat_vec,0.975)>m_A ) ) { cover_num_boot=cover_num_boot+1 }
   
   cilength_boot[b]=quantile(mab_hat_vec,0.975)-quantile(mab_hat_vec,0.025)
   
   ## process counter
   if ( b%%(B/10)==1 ) 
   {cat(c( (b-1)/(B/10))*10 ) 
     cat("% \n")}
} # end loop for b




result_matrix[1,jj]=mean( ma_est_vec )

result_matrix[2,jj]=mean(varma_est_vec)

cover_rate_phi= cover_num_phi /B

result_matrix[3,jj]=cover_rate_phi

result_matrix[4,jj]=mean( cilength_phi)

cover_rate_t= cover_num_t /B

result_matrix[5,jj]=cover_rate_t

result_matrix[6,jj]=mean( cilength_t)



result_matrix[7,jj]=mean(mab_vec)

result_matrix[8,jj]=mean( mab_var_vec)

cover_rate_boot= cover_num_boot /B

result_matrix[9,jj]=cover_rate_boot

result_matrix[10,jj]=mean( cilength_boot)


## estimate based on z0

z0_e_mean=apply(z0J,1,mean)
result_matrix_2[1,jj]=mean(z0_e_mean)

z0_e_var=apply(z0J,1,var)/J
result_matrix_2[2,jj]=mean(z0_e_var)

z0_l=z0_e_mean-qnorm(0.975)*sqrt(z0_e_var)
z0_u=z0_e_mean+qnorm(0.975)*sqrt(z0_e_var)
z0_ml=mean(z0_u-z0_l)

z0_CR= sum( (z0_l<m_A) & (z0_u>m_A)  )/B
result_matrix_2[3,jj]=z0_CR

result_matrix_2[4,jj]=z0_ml


z0_hl_t=qt(0.975,J-1)* sqrt(z0_e_var)

z0_l_t=z0_e_mean-z0_hl_t

z0_u_t=z0_e_mean+z0_hl_t

z0_CR_t=sum( (z0_l_t<m_A) & (z0_u_t>m_A)  )/B
result_matrix_2[5,jj]=z0_CR_t

result_matrix_2[6,jj]=mean(2*z0_hl_t)
   
# mafinalvec[ii,jj]=mean( ma_est_vec)
# mavarfinalvec[ii,jj]=mean( varma_est_vec)
 

# mean( varma_est_vec)
# r_o^(2*N)/J*mfDn/mo2

} # end loop for jj
# } # end loop for ii





mfD=m_A*gamma2*((1)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1)/(1-eta))

mfDn=m_A*gamma2*((1-kap^N)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1-eta^N)/(1-eta))

Ntau=c(tau:N)

mfDtau=m_A*gamma2*((1-kap^Ntau)/(1-kap))+ (m_A^2+sigma2_A)*sigma2*((1-eta^Ntau)/(1-eta))


## without sigma_A
# var_the_asy=(mfDn/mo2*r_o^(2*N)) /Jvec



## with sigma_A 

# ma(1)

# var_the_whole=(mfDn/mo2*r_o^(2*N)+sigma2_A) /Jvec

## ma(2)

N_J_t=( ( mo-1)/( mo^tau*( mo^(N-tau+1)-1)) )

tauvec=tau:N

var_znvec=mfDtau*mo2^(tauvec-1)+mo^(2*tauvec)*sigma2_A

var_zntau=0

for (l in tau:(N-1) )
{
  var_zntau= var_zntau +( 1+2*  (sum( mo^(1:(N-l)) ))  ) *var_znvec[l-tau+1]
  print(var_zntau)
}
var_zntau=var_zntau+ var_znvec[N-tau+1]

var_the_whole=N_J_t^2*var_zntau/Jvec

var_the_whole




##

var_the_asy=var_the_whole


result_matrix=round(result_matrix,3)

##without asy
# table_mat=result_matrix
# 
# table_name=c("J", " $\\hat{m}_A$" ,"$\\hat{ \\Lambda}_{2,n,J}$ ", "G CR" ,"G ML" ,"t CR" ,"t ML" ,"$\\hat{m}_{A,b}$" ,"Var$_b$","b CR","b ML")

##with asy 

table_mat=matrix(0,nrow=nrow(result_matrix)+1 ,ncol=ncol(result_matrix))

table_mat[1:2,]=result_matrix[1:2,]

table_mat[3,]=round(var_the_asy,3)

table_mat[4:nrow(table_mat), ]=result_matrix[3:nrow(result_matrix),  ]

table_name=c("J", " $\\hat{m}_A$" ,"$\\hat{ \\Lambda}_{1,n,J}$ ", "Asy Var " ,"G CR" ,"G ML" ,"t CR" ,"t ML" ,"$\\hat{m}_{A,b}$" ,"Var$_b$","b CR","b ML")

{
  cat( "\\begin{center}\n")
  cat("\\begin{tabular}[H]{|c| c | c  | c | c | } \n")
  cat("\\hline\n")
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
 cat("\\end{tabular}\n")
 cat("\\end{center}\n")
 
}



c(m_A,N,tau)



round(result_matrix_2,3)

table_mat2=round(result_matrix_2,3)

table_name2=c("$\\bar{Z_0}$","var$(Z_0)$ ","$(Z_0)$ G CR ","$(Z_0)$ G ML","$(Z_0)$ t CR","$(Z_0)$ t ML")

for (jj in 1:(dim(table_mat2)[1]) )
{
  cat(table_name2[jj])
  for (ii in 1:4)
  {
 
      cat(" & ")
      cat(table_mat2[jj,ii])   
    
  }
  cat(" \\\\ \n")
  cat("\\hline \n")
}


cat(paste0("$m_A=",m_A,"$, $n=",N, "$, $\\tau=",tau,"$: \n"))





