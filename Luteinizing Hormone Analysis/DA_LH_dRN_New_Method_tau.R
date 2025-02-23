# J = dim(F)[1]
# k = dim(F)[2]
#
# b1 <- matrix(0, 16, 40)
# for (i in 1:16){
#   b1[i, ] = a1[ ,i]
# }
#
# b2 <- matrix(0, 16, 40)
# for (i in 1:16){
#   b2[i, ] = a2[ ,i]
# }
#
# ldata1 <- b1[c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), ]
# ldata2 <- b2[c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16 ), ]
#
# F = ldata2
#
# # ctau1 = c(21, 21, 20, 20, 20, 20, 20, 20, 20, 19, 20, 20, 20, 20, 20)
# # ctau2 = c(24, 23, 23, 23, 23, 23, 24, 23, 24, 23, 24, 24, 24, 24, 24)
#
# ctau1 = matrix(20, 15, 1)
# ctau2 = matrix(24, 15, 1)
#
# cc1 = 0
# for (j in 1:J){ cc1 = cc1 + sum(F[j, ctau1[j]:ctau2[j]])}
# F_hat = cc1/J
#
# m1 = matrix(0, J, 1)
# for (j in 1:J){m1[j] = sum(F[j, (ctau1[j]+1):ctau2[j]])}
# m2 = matrix(0, J, 1)
# for (j in 1:J){m2[j] = sum(F[j, ctau1[j]:(ctau2[j]-1)])}
# m_star = sum(m1)/sum(m2)
# ma_hat = F_hat*(m_star-1)/(m_star^mean(ctau1)*(m_star^(mean(ctau2)-mean(ctau1)+1)-1))
# ma_hat
#
#
# F = ldata1
# # ctau1 =  c(18, 19, 19, 19, 18, 19, 18, 18, 18, 19, 18, 19, 18, 19, 19)
# # ctau2 =  c(21, 22, 22, 22, 22, 22, 22, 22, 21, 22, 22, 22, 22, 22, 22)
#
# ctau1 = matrix(19, 15, 1)
# ctau2 = matrix(22, 15, 1)
#
# cc1 = 0
# for (j in 1:J){ cc1 = cc1 + sum(F[j, ctau1[j]:ctau2[j]])}
# F_hat = cc1/J
#
# m1 = matrix(0, J, 1)
# for (j in 1:J){m1[j] = sum(F[j, (ctau1[j]+1):ctau2[j]])}
# m2 = matrix(0, J, 1)
# for (j in 1:J){m2[j] = sum(F[j, ctau1[j]:(ctau2[j]-1)])}
# m_star = sum(m1)/sum(m2)
# ma_hat = F_hat*(m_star-1)/(m_star^mean(ctau1)*(m_star^(mean(ctau2)-mean(ctau1)+1)-1))
# ma_hat


#### another method tau is not a vector, which works fine

a1 <- as.matrix(sapply(LH1, as.numeric))
a2 <- as.matrix(sapply(LH2, as.numeric))

b1 <- matrix(0, 16, 40)
for (i in 1:16){
  b1[i, ] = a1[ ,i]
}

b2 <- matrix(0, 16, 40)
for (i in 1:16){
  b2[i, ] = a2[ ,i]
}

ldata1 <- b1[c(1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16), ]
ldata2 <- b2[c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16 ), ]

J = dim(ldata1)[1]

F = ldata2
ctau1 = 19
ctau2 = 24
cc1 = 0
for (j in 1:J){ cc1 = cc1 + sum(F[j, ctau1:ctau2])}
F_hat = cc1/J
m1 = matrix(0, J, 1)
for (j in 1:J){m1[j] = sum(F[j, (ctau1+1):ctau2])}
m2 = matrix(0, J, 1)
for (j in 1:J){m2[j] = sum(F[j, ctau1:(ctau2-1)])}
m_star = sum(m1)/sum(m2)
ma_hat1 = F_hat*(m_star-1)/(m_star^ctau1*(m_star^(ctau2-ctau1+1)-1))
ma_hat1
# get 1.3151e-05
var_mu = matrix(0, J, 1)
for (j in 1:J){var_mu[j] = ((m_star-1)/(m_star^ctau1*(m_star^(ctau2-ctau1+1)-1))*sum(F[j, ctau1:ctau2]) -ma_hat1)^2/(J-1)}
sigma11 = sum(var_mu)/ma_hat1^2


F = ldata1
ctau1 = 19
ctau2 = 22
cc1 = 0
for (j in 1:J){ cc1 = cc1 + sum(F[j, ctau1:ctau2])}
F_hat = cc1/J
m1 = matrix(0, J, 1)
for (j in 1:J){m1[j] = sum(F[j, (ctau1+1):ctau2])}
m2 = matrix(0, J, 1)
for (j in 1:J){m2[j] = sum(F[j, ctau1:(ctau2-1)])}
m_star = sum(m1)/sum(m2)
ma_hat2 = F_hat*(m_star-1)/(m_star^ctau1*(m_star^(ctau2-ctau1+1)-1))
ma_hat2
# get 3.7431e-05

var_mu2 = matrix(0, J, 1)
for (j in 1:J){var_mu2[j] = ((m_star-1)/(m_star^ctau1*(m_star^(ctau2-ctau1+1)-1))*sum(F[j, ctau1:ctau2]) -ma_hat2)^2/(J-1) }
sigma21 = sum(var_mu2)/ma_hat2^2

ma_hat2/ma_hat1
# the ratio is 2.8462

variance1 = 2.8462^2*(sigma11 + sigma21)

# 95% CI
2.8462-1.96*sqrt(variance1)
2.8462+1.96*sqrt(variance1)







