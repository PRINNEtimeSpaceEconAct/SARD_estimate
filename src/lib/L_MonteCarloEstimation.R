createModelVariables <- function(SARDp){
    regressors = as.formula(delta ~ 0)
    instruments = as.formula(~ 0)
    variables = c()
    if (SARDp$gammaS != 0){
        regressors=add.terms(regressors,c("xS","MSDelta"))
        instruments=add.terms(instruments,c("xS","MS2X"))
        variables[length(variables)+1] = "xS"
    }
    if (SARDp$gammaA != 0){
        regressors=add.terms(regressors,c("xA","MADelta"))
        instruments=add.terms(instruments,c("xA","MA2X"))
        variables[length(variables)+1] = "xA"
    }
    if (SARDp$gammaR != 0){
        regressors=add.terms(regressors,c("xR","MRDelta"))
        instruments=add.terms(instruments,c("xR","MR2X"))
        variables[length(variables)+1] = "xR"
    }
    if (SARDp$gammaD != 0){
        regressors=add.terms(regressors,c("xD","MDDelta"))
        instruments=add.terms(instruments,c("xD","MD2X"))
        variables[length(variables)+1] = "xD"
    }
    
    model = listN(regressors,instruments)    
    return(listN(model,variables))
}

estimate_ARD_MC_LM_NAIVE <- function(df,shp,xA,xR,xD,MA,MR,MD){
    
    # lag regressors
    MADelta = as.numeric(MA %*% matrix(df$delta))
    MRDelta = as.numeric(MR %*% matrix(df$delta))
    MDDelta = as.numeric(MD %*% matrix(df$delta))
    
    df = df %>% mutate(xA = xA, xR = xR, xD = xD)
    
    LM_est = lm(delta ~ y0 + xA  + xR  + xD, data = df)
    AICc = AICc(LM_est)[1]
    R2 = AICc2R2Nagelkerke(AICc,df$delta,4)
    
    return(listN(LM_est,AICc,R2))
}

estimate_ARD_MC_LM_NAIVEBoot <- function(df,shp,xA,xR,xD,MADelta,MRDelta,MDDelta){
    
    df = df %>% mutate(xA = xA, xR = xR, xD = xD)
    
    LM_est = lm(delta ~ y0 + xA  + xR  + xD, data = df)
    return(LM_est)
}


estimate_ARD_MC_LM <- function(df,shp,xA,xR,xD,MA,MR,MD,weights=rep(1,nrow(df))){
    
    # lag regressors
    MADelta = as.numeric(MA %*% matrix(df$delta))
    MRDelta = as.numeric(MR %*% matrix(df$delta))
    MDDelta = as.numeric(MD %*% matrix(df$delta))
    
    df = df %>% mutate(xA = xA, xR = xR, xD = xD,
                       MADelta=MADelta, MRDelta=MRDelta,MDDelta=MDDelta)
    
    LM_est = lm(delta ~ y0 + xA + MADelta + xR + MRDelta + xD + MDDelta, data = df, weights = weights)
    
    AICc = AICc(LM_est)[1]
    R2 = AICc2R2Nagelkerke(AICc,df$delta,8)
        
    return(listN(LM_est,AICc,R2))
}

estimate_ARD_MC_LMBoot <- function(df,shp,xA,xR,xD,MADelta,MRDelta,MDDelta,weights=rep(1,nrow(df))){
    
    
    df = df %>% mutate(xA = xA, xR = xR, xD = xD,
                       MADelta=MADelta, MRDelta=MRDelta,MDDelta=MDDelta)
    
    LM_est = lm(delta ~ y0 + xA + MADelta + xR + MRDelta + xD + MDDelta, data = df, weights = weights)
    return(LM_est)
}

estimate_ARD_MC_IV <- function(df,shp,xA,xR,xD,MA,MR,MD){
    
    # lag regressors
    MADelta = as.numeric(MA %*% matrix(df$delta))
    MRDelta = as.numeric(MR %*% matrix(df$delta))
    MDDelta = as.numeric(MD %*% matrix(df$delta))
    
    X = data.frame(xA=xA,xR=xR,xD=xD)
    X = as.matrix(X)
    
    # instruments for IV
    MA2X=as.matrix(MA %*% MA %*% X)
    MR2X=as.matrix(MR %*% MR %*% X) 
    MD2X=as.matrix(MD %*% MD %*% X)
    
    df = df %>% mutate(xA = xA, xR = xR, xD = xD,
                       MADelta=MADelta, MRDelta=MRDelta,MDDelta=MDDelta,
                       MA2X=MA2X[,],MR2X=MR2X[,],MD2X=MD2X[,])
    
    IV_est = ivreg(delta ~ y0 + xA + xR + xD + 
                       + MADelta + MRDelta + MDDelta  | 
                       xA + xR + xD + MA2X + MR2X + MD2X , data=df)
    
    LAR = LogLikAICcR2_MC(df,IV_est$coefficients,8,xA, xR, xD, MA, MR, MD)
    AICc = LAR$AICc
    R2 = LAR$R2Nagelkerke
    
    return(listN(IV_est,AICc,R2))
}

estimate_ARD_MC_IVBoot <- function(df,shp,xA,xR,xD,MADelta,MRDelta,MDDelta,MA2X,MR2X,MD2X){
    

    df = df %>% mutate(xA = xA, xR = xR, xD = xD,
                       MADelta=MADelta, MRDelta=MRDelta,MDDelta=MDDelta,
                       MA2X=MA2X[,],MR2X=MR2X[,],MD2X=MD2X[,])
    
    IV_est = ivreg(delta ~ y0 + xA + xR + xD + 
                       + MADelta + MRDelta + MDDelta  | 
                       xA + xR + xD + MA2X + MR2X + MD2X , data=df)
    
    return(IV_est)
}

estimate_ARD_MC_LL <- function(df,shp,xA,xR,xD,MA,MR,MD){
    
    Y = as.matrix(df$delta)
    X = data.frame(ones=df$ones,y0=df$y0,xA=xA,xR=xR,xD=xD)
    X = as.matrix(X)
    
    coefLM = coef(estimate_ARD_MC_LM(df,shp,xA,xR,xD,MA,MR,MD)$LM_est)
    initialOptim = coefLM[c("MADelta","MRDelta","MDDelta")]
    outARD_WN3MatEstimate = call_julia_LogLik_WN_3Mat(X,Y,MA,MR,MD,initialOptim)
    
    SpatError = compute_spatial_error_mat(outARD_WN3MatEstimate$residuals,shp,maxLag = 10)
    Werr = SpatError$Werr
    initialOptim = c(outARD_WN3MatEstimate$coef[c(6,7,8)],0.5)  
    outARD_3MatEstimate = call_julia_LogLik_3Mat(X,Y,MA,MR,MD,Werr,initialOptim)
    
    LAR = LogLikAICcR2_MC_SpatError(df,outARD_3MatEstimate$coef,9,xA, xR, xD, MA, MR, MD, Werr)
    AICc = LAR$AICc
    R2 = LAR$R2Nagelkerke
    
    
    return(listN(outARD_WN3MatEstimate,outARD_3MatEstimate,SpatError,AICc,R2))
}

LogLikAICcR2_MC_SpatError <- function(df, coef, k, xA, xR, xD, MA, MR, MD, W_eps){
    # compute LogLik and AICc of the SARD Model, with dof degree of freedom
    # coefs is of length 9, in order
    # a, phi, gammA, gammaR, gammaD, rhoA, rhoR, rhoD, lambda.
        
    if (DEBUG == TRUE){ print("computing LogLik") }
    
    # thankyou R for being so modern    
    a=coef[1]; phi=coef[2]; gammA=coef[3]; gammaR = coef[4];
    gammaD=coef[5]; rhoA=coef[6]; rhoR=coef[7]; rhoD=coef[8]; 
    lambda=coef[9]
    
    N = nrow(MA)
    Y = df$delta
    X = cbind(df$ones,df$y0,xA,xR,xD)
    colnames(X) <- c("ones","y0","xA","xR","xD")
    
    A = diag(N) - rhoA*MA - rhoR*MR - rhoD*MD
    B = diag(N) - lambda*W_eps
    spatY = A %*% as.numeric(Y)
    
    if (det(t(X) %*% X) != 0) {
        beta = solve( t(X) %*% X, t(X) %*% spatY)}
    else {     return(list(LogLiKelyhood=NA,AICc=NA,R2Nagelkerke=NA)) }
    # beta = solve( t(X) %*% X, t(X) %*% spatY)
    
    epsilon = spatY - X %*% beta
    
    mu = B %*% epsilon
    sigma2 = as.numeric(t(mu) %*% mu) / N
    Omega = sigma2 * diag(N)
    nu = mu / sqrt(sigma2)
    
    logAbsDetA = determinant(A)$modulus   # log(|det(A)|)
    logAbsDetB = determinant(B)$modulus
    
    LogLiKelyhood = -(N/2)*(log(2*pi)) - (N/2)*log(sigma2) + 
        + logAbsDetB + logAbsDetA - 1/2 * as.numeric(t(nu) %*% nu)
    AIC = 2*k - 2*LogLiKelyhood
    AICc = AIC + (2*k^2+2*k)/(N-k-1)
    
    LL0 = logLik(lm(Y ~ 1))
    R2Nagelkerke =  c(1 - exp(-(2/N)*(LogLiKelyhood - LL0)))
    
    return(listN(LogLiKelyhood,AICc,R2Nagelkerke))
}

LogLikAICcR2_MC <- function(df, coef, k, xA, xR, xD, MA, MR, MD){
    # compute LogLik and AICc of the SARD Model, with dof degree of freedom
    # coefs is of length 8, in order
    # a, phi, gammA, gammaR, gammaD, rhoA, rhoR, rhoD
    
    
    if (DEBUG == TRUE){ print("computing LogLik") }
    
    # thankyou R for being so modern    
    a=coef[1]; phi=coef[2]; gammA=coef[3]; gammaR = coef[4];
    gammaD=coef[5]; rhoA=coef[6]; rhoR=coef[7]; rhoD=coef[8]; 

    N = nrow(MA)
    Y = df$delta
    X = cbind(df$ones,df$y0,xA,xR,xD)
    colnames(X) <- c("ones","y0","xA","xR","xD")
    
    A = diag(N) - rhoA*MA - rhoR*MR - rhoD*MD
    B = diag(N) 
    spatY = A %*% as.numeric(Y)
    
    if (det(t(X) %*% X) != 0) {
        beta = solve( t(X) %*% X, t(X) %*% spatY)}
    else {     return(list(LogLiKelyhood=NA,AICc=NA,R2Nagelkerke=NA)) }
    # beta = solve( t(X) %*% X, t(X) %*% spatY)
    
    epsilon = spatY - X %*% beta
    
    mu = B %*% epsilon
    sigma2 = as.numeric(t(mu) %*% mu) / N
    Omega = sigma2 * diag(N)
    nu = mu / sqrt(sigma2)
    
    logAbsDetA = determinant(A)$modulus   # log(|det(A)|)
    logAbsDetB = determinant(B)$modulus
    
    LogLiKelyhood = -(N/2)*(log(2*pi)) - (N/2)*log(sigma2) + 
        + logAbsDetB + logAbsDetA - 1/2 * as.numeric(t(nu) %*% nu)
    AIC = 2*k - 2*LogLiKelyhood
    AICc = AIC + (2*k^2+2*k)/(N-k-1)
    
    LL0 = logLik(lm(Y ~ 1))
    R2Nagelkerke =  c(1 - exp(-(2/N)*(LogLiKelyhood - LL0)))
    
    return(listN(LogLiKelyhood,AICc,R2Nagelkerke))
}
