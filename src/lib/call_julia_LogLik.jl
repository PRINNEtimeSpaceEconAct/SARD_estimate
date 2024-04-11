
# using LinearAlgebra
# using Optimization
# using OptimizationOptimJL
# using SparseArrays
# using Distributions
# using FiniteDiff

function julia_LogLik_4Mat(X,Y,MS,MA,MR,MD,Weps,initialCondition)
    """
    Estimate SARD model by using LogLikelyhood
    """

    MS = sparse(Matrix{Float64}(MS)); dropzeros!(MS)
    MA = sparse(Matrix{Float64}(MA)); dropzeros!(MA)
    MR = sparse(Matrix{Float64}(MR)); dropzeros!(MR)
    MD = sparse(Matrix{Float64}(MD)); dropzeros!(MD)
    Weps = sparse(Matrix{Float64}(Weps)); dropzeros!(Weps)

    lb = [-1,-200,-200,-1,-1]
    ub = [1,200,200,1,1]
    p_start = initialCondition
    
    # starting optimization
    f_opt = OptimizationFunction(negLogLik_4MatSpatErr,Optimization.AutoFiniteDiff())
    prob = OptimizationProblem(f_opt,p_start,(Y,X,MS,MA,MR,MD,Weps),
    lb=lb,ub=ub)
    opt = Optimization.solve(prob,Optim.NelderMead(),show_trace=true)
    optResult = opt.u

    XNormal = Normal(0,1)

    # computing se, pvalue of beta
    N = size(MS,1)
    A = (I(N) - optResult[1]*MS - optResult[2]*MA +
                - optResult[3]*MR - optResult[4]*MD)
    B = (I(N)-optResult[5]*Weps)
    spatY = A * Y
    beta = (transpose(X)*X) \ ( transpose(X) * spatY )
    resid = spatY - X*beta
    mu = B * resid
    sigma2est = sqrt(dot(mu,mu))/N
    covBeta = sigma2est * inv(transpose(X)*transpose(B)*B*X)
    seBeta = sqrt.(diag(covBeta))
    tStatBeta = beta./seBeta
    pValueBeta = 2*cdf.(XNormal,-abs.(tStatBeta))

    # computing se, pvalue of thetas (lags)
    LL(p) =  LogLik_4MatSpatErr(p,(Y,X,MS,MA,MR,MD,Weps))
    Hess = -FiniteDiff.finite_difference_hessian(LL,optResult)
    covTheta = inv(Hess)
    seTheta = sqrt.(diag(covTheta))
    tStatTheta = optResult./seTheta
    pValueTheta = 2*cdf.(XNormal,-abs.(tStatTheta))

    # collecting output
    coef = [beta..., optResult...]
    se_coef = [seBeta...,seTheta...]
    pvalue_coef = [pValueBeta...,pValueTheta...]
    residuals = resid

    return(coef, se_coef, pvalue_coef, residuals)
end

function LogLik_4MatSpatErr(p,param)  
    # LogLikelyhood of spatial lag model with 4 matrices and no spatial error

    y,x,w1,w2,w3,w4,weps = param

    N = size(w1,1)
    A = I(N) - p[1]*w1 - p[2]*w2 - p[3]*w3 - p[4]*w4
    B = I(N) - p[5]*weps
    spatY = A*y
    beta = (transpose(x)*x) \ ( transpose(x) * spatY )
    u = spatY - x * beta
    mu = B * u
    sigma2 = dot(mu,mu) / N
    nu = mu / sqrt(sigma2)

    logAbsDetA = logabsdet(A)
    logAbsDetB = logabsdet(B)

    LogLik = -(N/2)*(log(2*pi)) - (N/2)*log(sigma2) + logAbsDetB[1] + logAbsDetA[1] - 1/2*dot(nu,nu)
    
    return LogLik
end

function negLogLik_4MatSpatErr(p,param)
    # negative LogLikelyhood to minimize
    return -LogLik_4MatSpatErr(p,param)
end

function julia_LogLik_WN_4Mat(X,Y,MS,MA,MR,MD,initialCondition)
    """
    Estimate SARD WN model by using LogLikelyhood
    """

    MS = sparse(Matrix{Float64}(MS)); dropzeros!(MS)
    MA = sparse(Matrix{Float64}(MA)); dropzeros!(MA)
    MR = sparse(Matrix{Float64}(MR)); dropzeros!(MR)
    MD = sparse(Matrix{Float64}(MD)); dropzeros!(MD)

    lb = [-1,-200,-200,-1]
    ub = [1,200,200,1]
    p_start = initialCondition
    
   
    # starting optimization
    f_opt = OptimizationFunction(negLogLik_4Mat,Optimization.AutoFiniteDiff())
    prob = OptimizationProblem(f_opt,p_start,(Y,X,MS,MA,MR,MD),
    lb=lb,ub=ub)
    opt = Optimization.solve(prob,Optim.NelderMead(),show_trace=true)
    optResult = opt.u

    XNormal = Normal(0,1)

    # computing se, pvalue of beta
    N = size(MS,1)
    A = (I(N) - optResult[1]*MS - optResult[2]*MA +
                - optResult[3]*MR - optResult[4]*MD)
    spatY = A * Y
    beta = (transpose(X)*X) \ ( transpose(X) * spatY )
    resid = spatY - X*beta
    sigma2est = sqrt(dot(resid,resid))/N
    covBeta = sigma2est * inv(transpose(X)*X)
    seBeta = sqrt.(diag(covBeta))
    tStatBeta = beta./seBeta
    pValueBeta = 2*cdf.(XNormal,-abs.(tStatBeta))

    # computing se, pvalue of thetas (lags)
    LL(p) =  LogLik_4Mat(p,(Y,X,MS,MA,MR,MD))
    Hess = -FiniteDiff.finite_difference_hessian(LL,optResult)
    covTheta = inv(Hess)
    seTheta = sqrt.(diag(covTheta))
    tStatTheta = optResult./seTheta
    pValueTheta = 2*cdf.(XNormal,-abs.(tStatTheta))

    # collecting output
    coef = [beta..., optResult...]
    se_coef = [seBeta...,seTheta...]
    pvalue_coef = [pValueBeta...,pValueTheta...]
    residuals = resid

    return(coef, se_coef, pvalue_coef, residuals)
end

function LogLik_4Mat(p,param)  
    # LogLikelyhood of spatial lag model with 4 matrices and no spatial error

    y,x,w1,w2,w3,w4 = param

    N = size(w1,1)
    A = I(N) - p[1]*w1 - p[2]*w2 - p[3]*w3 - p[4]*w4
    spatY = A*y
    beta = (transpose(x)*x) \ ( transpose(x) * spatY )
    u = spatY - x * beta
    sigma2 = dot(u,u) / N

    logAbsDet = logabsdet(A)
    LogLik = -(N/2)*(log(2*pi)) - (N/2)*log(sigma2) + logAbsDet[1] - N/2
    return LogLik
end

function negLogLik_4Mat(p,param)
    # negative LogLikelyhood to minimize
    return -LogLik_4Mat(p,param)
end

function julia_LogLik_WN_3Mat(X,Y,MA,MR,MD,initialCondition)
    """
    Estimate SARD WN model by using LogLikelyhood
    """

    MA = sparse(Matrix{Float64}(MA)); dropzeros!(MA)
    MR = sparse(Matrix{Float64}(MR)); dropzeros!(MR)
    MD = sparse(Matrix{Float64}(MD)); dropzeros!(MD)

    lb = [-1,-1,-1]
    ub = [1,1,1]
    p_start = initialCondition
    
   
    # starting optimization
    f_opt = OptimizationFunction(negLogLik_3Mat,Optimization.AutoFiniteDiff())
    prob = OptimizationProblem(f_opt,p_start,(Y,X,MA,MR,MD),
    lb=lb,ub=ub)
    opt = Optimization.solve(prob,Optim.NelderMead(),show_trace=false)
    optResult = opt.u

    XNormal = Normal(0,1)

    # computing se, pvalue of beta
    N = size(MA,1)
    A = (I(N) - optResult[1]*MA +
                - optResult[2]*MR - optResult[3]*MD)
    spatY = A * Y
    beta = (transpose(X)*X) \ ( transpose(X) * spatY )
    resid = spatY - X*beta
    sigma2est = sqrt(dot(resid,resid))/N

    covBeta = sigma2est * inv(transpose(X)*X)
    seBeta = sqrt.(diag(covBeta))
    tStatBeta = beta./seBeta
    pValueBeta = 2*cdf.(XNormal,-abs.(tStatBeta))

    # computing se, pvalue of thetas (lags)
    LL(p) =  LogLik_3Mat(p,(Y,X,MA,MR,MD))
    Hess = -FiniteDiff.finite_difference_hessian(LL,optResult)
    covTheta = inv(Hess)
    seTheta = sqrt.(diag(covTheta))
    tStatTheta = optResult./seTheta
    pValueTheta = 2*cdf.(XNormal,-abs.(tStatTheta))

    # collecting output
    coef = [beta..., optResult...]
    se_coef = [seBeta...,seTheta...]
    pvalue_coef = [pValueBeta...,pValueTheta...]
    residuals = resid

    return(coef, se_coef, pvalue_coef, residuals)
end

function LogLik_3Mat(p,param)  
    # LogLikelyhood of spatial lag model with 4 matrices and no spatial error

    y,x,w1,w2,w3 = param

    N = size(w1,1)
    A = I(N) - p[1]*w1 - p[2]*w2 - p[3]*w3 
    spatY = A*y
    beta = (transpose(x)*x) \ ( transpose(x) * spatY )
    u = spatY - x * beta
    sigma2 = dot(u,u) / N

    logAbsDet = logabsdet(A)
    LogLik = -(N/2)*(log(2*pi)) - (N/2)*log(sigma2) + logAbsDet[1] - N/2
    return LogLik
end

function negLogLik_3Mat(p,param)
    # negative LogLikelyhood to minimize
    return -LogLik_3Mat(p,param)
end

function julia_LogLik_3Mat(X,Y,MA,MR,MD,Weps,initialCondition)
    """
    Estimate SARD model by using LogLikelyhood
    """

    MA = sparse(Matrix{Float64}(MA)); dropzeros!(MA)
    MR = sparse(Matrix{Float64}(MR)); dropzeros!(MR)
    MD = sparse(Matrix{Float64}(MD)); dropzeros!(MD)
    Weps = sparse(Matrix{Float64}(Weps)); dropzeros!(Weps)

    lb = [-1,-1,-1,-1]
    ub = [1,1,1,1]
    p_start = initialCondition
    
    # starting optimization
    f_opt = OptimizationFunction(negLogLik_3MatSpatErr,Optimization.AutoFiniteDiff())
    prob = OptimizationProblem(f_opt,p_start,(Y,X,MA,MR,MD,Weps),
    lb=lb,ub=ub)
    opt = Optimization.solve(prob,Optim.NelderMead(),show_trace=false)
    optResult = opt.u

    XNormal = Normal(0,1)

    # computing se, pvalue of beta
    N = size(MA,1)
    A = (I(N) - optResult[1]*MA +
                - optResult[2]*MR - optResult[3]*MD)
    B = (I(N)-optResult[4]*Weps)
    spatY = A * Y
    beta = (transpose(X)*X) \ ( transpose(X) * spatY )
    resid = spatY - X*beta
    mu = B * resid
    sigma2est = sqrt(dot(mu,mu))/N
    covBeta = sigma2est * inv(transpose(X)*transpose(B)*B*X)
    seBeta = sqrt.(diag(covBeta))
    tStatBeta = beta./seBeta
    pValueBeta = 2*cdf.(XNormal,-abs.(tStatBeta))

    # computing se, pvalue of thetas (lags)
    LL(p) =  LogLik_3MatSpatErr(p,(Y,X,MA,MR,MD,Weps))
    Hess = -FiniteDiff.finite_difference_hessian(LL,optResult)
    covTheta = inv(Hess)
    seTheta = sqrt.(diag(covTheta))
    tStatTheta = optResult./seTheta
    pValueTheta = 2*cdf.(XNormal,-abs.(tStatTheta))

    # collecting output
    coef = [beta..., optResult...]
    se_coef = [seBeta...,seTheta...]
    pvalue_coef = [pValueBeta...,pValueTheta...]
    residuals = resid

    return(coef, se_coef, pvalue_coef, residuals, covBeta)
end

function LogLik_3MatSpatErr(p,param)  
    # LogLikelyhood of spatial lag model with 4 matrices and no spatial error

    y,x,w1,w2,w3,weps = param

    N = size(w1,1)
    A = I(N) - p[1]*w1 - p[2]*w2 - p[3]*w3
    B = I(N) - p[4]*weps
    spatY = A*y
    beta = (transpose(x)*x) \ ( transpose(x) * spatY )
    u = spatY - x * beta
    mu = B * u
    sigma2 = dot(mu,mu) / N
    nu = mu / sqrt(sigma2)

    logAbsDetA = logabsdet(A)
    logAbsDetB = logabsdet(B)

    LogLik = -(N/2)*(log(2*pi)) - (N/2)*log(sigma2) + logAbsDetB[1] + logAbsDetA[1] - 1/2*dot(nu,nu)
    
    return LogLik
end

function negLogLik_3MatSpatErr(p,param)
    # negative LogLikelyhood to minimize
    return -LogLik_3MatSpatErr(p,param)
end

