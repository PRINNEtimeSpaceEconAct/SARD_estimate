using CovarianceEstimation
using Distributions
using LinearAlgebra

N = 1000
X = Normal(0,1)
A = rand(X,N,N)
Σ = A'*A

Xepsilon = MvNormal(zeros(N),Σ)
epsilon = reshape(Matrix(rand(Xepsilon,1)),1,N)
# epsilon = vcat(epsilon,epsilon)
Σ̂ = cov(LinearShrinkage(DiagonalUnequalVariance()),epsilon)
Σ̂ = cov(LinearShrinkage(ConstantCorrelation()),epsilon)
Σ̂ = cov(SimpleCovariance(corrected=true),epsilon)
Σ̂ = cov(SimpleCovariance(corrected=false),epsilon)

