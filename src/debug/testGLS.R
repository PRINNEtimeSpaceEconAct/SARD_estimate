library(nlme)
library(MASS)

N = 1000
alpha = 0.1
beta = c(0.5,0.7)
sigma = 0.1

X = matrix(rnorm(2*N),ncol=2)

# easy errors ----
eps = sigma*rnorm(N)
Y = alpha + X %*% beta + eps
data = data.frame(Y,X)

## OLS ----
summary(lm(Y ~ X, data))

## GLS ----
summary(gls(Y ~ X, data))

# bad errors ----
A = matrix(rnorm(N*N),nrow=N)
Omega = t(A) %*% A 
eps = sigma * mvrnorm(n=1,mu=rep(0,N),Omega)

Y = alpha + X %*% beta + eps
data = data.frame(Y,X)

## OLS ----
summary(lm(Y ~ X, data))

## GLS ----
summary(gls(Y ~ X, data))
