

using GLM
using Distributions

N01 = Normal(0,1)
N = 1000
X = hcat(ones(N,1),reshape(rand(N01,2*N),N,2))
beta = reshape([0.1, 0.5, 0.7],3,1)

# easy errors 
eps = 0.1*rand(N01,N)
Y = (X*beta + eps)[:]

## OLS 
lm(X,Y)

## GLS
glm(X,Y,Normal())

# spatially correlated errors
A = rand(N01,N,N)
Σ = A'*A/N
Xepsilon = MvNormal(zeros(N),Σ)

eps = 0.1*rand(Xepsilon,1)[:]
Y = (X*beta + eps)[:]

## OLS
lm(X,Y)

## GLS
gls(X,Y,Normal())
