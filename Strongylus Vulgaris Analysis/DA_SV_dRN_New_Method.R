
sv1 <- Copy_number_3_dRn[1:11, 5:44]    # copy number 10^7
sv2 <- Copy_number_3_dRn[13:23, 5:44]   # copy number 10^6

aa1 <-  as.matrix(sapply(sv1, as.numeric))
aa2 <-  as.matrix(sapply(sv2, as.numeric))

# c1 = aa1[c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11), ]
# c2 = aa2[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), ]
# 
# J = dim(c1)[1]
# F = abs(c2)
# # > t(ctau1)
# # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# # [1,]   34   34   33   33   33   33   34   34   33    34
# # > t(ctau2)
# # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# # [1,]   35   35   35   35   35   35   36   36   35    36
# ctau1 = 33
# ctau2 = 35
# cc1 = 0
# for (j in 1:J){ cc1 = cc1 + sum(F[j, ctau1:ctau2])}
# F_hat = cc1/J
# m1 = matrix(0, J, 1)
# for (j in 1:J){m1[j] = sum(F[j, (ctau1+1):ctau2])}
# m2 = matrix(0, J, 1)
# for (j in 1:J){m2[j] = sum(F[j, ctau1:(ctau2-1)])}
# m_star = sum(m1)/sum(m2)
# ma_hat2 = F_hat*(m_star-1)/(m_star^ctau1*(m_star^(ctau2-ctau1+1)-1))
# ma_hat2
# 
# F = abs(c1)
# # > t(ctau1)
# # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# # [1,]   31   31   30   30   31   31   31   30   30    31
# # > t( ctau2)
# # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# # [1,]   33   32   32   32   32   33   32   32   32    33
# ctau1 = 30
# ctau2 = 32
# cc1 = 0
# for (j in 1:J){ cc1 = cc1 + sum(F[j, ctau1:ctau2])}
# F_hat = cc1/J
# 
# m1 = matrix(0, J, 1)
# for (j in 1:J){m1[j] = sum(F[j, (ctau1+1):ctau2])}
# m2 = matrix(0, J, 1)
# for (j in 1:J){m2[j] = sum(F[j, ctau1:(ctau2-1)])}
# m_star = sum(m1)/sum(m2)
# 
# ma_hat1 = F_hat*(m_star-1)/(m_star^ctau1*(m_star^(ctau2-ctau1+1)-1))
# ma_hat1
# ma_hat1/ma_hat2


# another estimator
### est for c2 (10^6) ###
c1 = aa1[c(1, 2, 4, 5, 6, 7, 8, 9, 10, 11), ]
c2 = aa2[c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11), ]

J = dim(c1)[1]
F = c2
# determine the starting cycles and ending cycles
r = dim(F)[1]
n = dim(F)[2]
tau1 <- 1.55
tau2 <- 1.55
ctau1 <- matrix(0, r, 1)
ctau2 <- matrix(0, r, 1)
cgam = matrix(0, r, 1)
for ( i in 1:r){ cgam[i] = min(which(F[i, ]>=.2))}
for (i in 1:r){ctau1[i] = min(which(F[i, (cgam[i]+1):n]/F[i, cgam[i]:(n-1)] >= tau1 ) + cgam[i] ) }
ctau1 = ctau1 - 1
for (i in 1:r){ctau2[i] = max(which(F[i, (cgam[i]+1):n]/F[i, cgam[i]:(n-1)] >= tau2 ) + cgam[i] ) }

# estimation
m1 = matrix(0, J, 1)
for (j in 1:J){m1[j] = sum(F[j, (ctau1[j] +1):ctau2[j] ])}
m2 = matrix(0, J, 1)
for (j in 1:J){m2[j] = sum(F[j, ctau1[j] :(ctau2[j]-1)])}
m_star = sum(m1)/sum(m2)
F1_hat = matrix(0, J, 1)
for (j in 1:J){F1_hat[j] = sum(F[j, ctau1[j]:ctau2[j]] ) }
m_hat = matrix(0, J, 1)
for (j in 1:J){m_hat[j] = F1_hat[j]*(m_star-1)/(m_star^ctau1[j]*(m_star^(ctau2[j] -ctau1[j]+1)-1)) }
ma_hat2 = mean(m_hat)

### -----############
# variance est
var_mu = matrix(0, J, 1)
for (j in 1:J){var_mu[j] = (m_hat[j] - mean(m_hat))^2/(J-1) }
sigma11 = mean(var_mu)/ma_hat2^2



### est for c1 (10^7) #####
F = c1
# determine the starting cycles and ending cycles
  r = dim(F)[1]
  n = dim(F)[2]
  ctau1 <- matrix(0, r, 1)
  ctau2 <- matrix(0, r, 1)
  cgam = matrix(0, r, 1)
  for ( i in 1:r){ cgam[i] = min(which(F[i, ]>=.2))}
  for (i in 1:r){ctau1[i] = min(which(F[i, (cgam[i]+1):n]/F[i, cgam[i]:(n-1)] >= tau1 ) + cgam[i] ) }
  ctau1 = ctau1 - 1
  for (i in 1:r){ctau2[i] = max(which(F[i, (cgam[i]+1):n]/F[i, cgam[i]:(n-1)] >= tau2 ) + cgam[i] ) }
# estimation 
m1 = matrix(0, J, 1)
for (j in 1:J){m1[j] = sum(F[j, (ctau1[j] +1):ctau2[j] ])}
m2 = matrix(0, J, 1)
for (j in 1:J){m2[j] = sum(F[j, ctau1[j] :(ctau2[j]-1)])}
m_star = sum(m1)/sum(m2)
F1_hat = matrix(0, J, 1)
for (j in 1:J){F1_hat[j] = sum(F[j, ctau1[j]:ctau2[j]] ) }
m_hat = matrix(0, J, 1)
for (j in 1:J){m_hat[j] = F1_hat[j]*(m_star-1)/(m_star^ctau1[j]*(m_star^(ctau2[j] -ctau1[j]+1)-1)) }
ma_hat1 = mean(m_hat)
### -----############

ma_hat1/ma_hat2


# variance est
var_mu1 = matrix(0, J, 1)
for (j in 1:J){var_mu1[j] = (m_hat[j] -mean(m_hat))^2/(J-1) }
sigma21 = mean(var_mu1)/ma_hat1^2



variance1 = 9.5972^2*(sigma11 + sigma21)

# 95% CI
9.5972-1.96*sqrt(variance1) 
9.5972+1.96*sqrt(variance1) 




