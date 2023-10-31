require(ivreg)

estimate_IV_SARD <- function(data,MsDeriv,xS,xD,MS,MD,D,hA,hR){
    # Estimate SARD WN via IV with distances hA, hR
    # return AICc
    
    WhA = compute_WhA(D,hA)
    WhR = compute_WhA(D,hR)
    xA = compute_xAR(data,MsDeriv, WhA)
    xR = compute_xAR(data,MsDeriv, WhR)
    MA = compute_MARLag(data,MsDeriv,WhA)
    MR = compute_MARLag(data,MsDeriv,WhR)
    
    MS2X = MA2X = MR2X = MD2X = None # instruments for IV
    
    IV_est = None # estimate via IV
    
    coef = nparam = W_eps = None # coef = IV_est$coef + 0 lambda, nparam = 6, W_eps = I
    LAR_IV = LogLikAICcR2(data, coef, nparam, xS, xA, xR, xD, 
                                                        MS, MA, MR, MD, W_eps)
    LogLik = LAR_IV$LogLik
    AICc = LAR_IV$AICc
    R2N = LAR_IV$R2Nagelkerke
    
    return(listN(IV_est, LogLik, AICc,R2N))
}

estimate_IV_SARD_auto <- function(df,hA,hR){
    # Estimate SARD WN via IV with distances hA, hR
    # return AICc
    
    MsDeriv = GFDM(df)
    D = compute_D(df)
    
    WhA = compute_WhAR(D,df,hA)
    WhR = compute_WhAR(D,df,hR)
    xS = compute_xS(df,MsDeriv)
    xA = compute_xAR(df,MsDeriv, WhA)
    xR = compute_xAR(df,MsDeriv, WhR)
    xD = compute_xD(df,MsDeriv)
    MS = compute_MSLag(df,MsDeriv)
    MA = compute_MARLag(df,MsDeriv,WhA,parallel=TRUE)
    MR = compute_MARLag(df,MsDeriv,WhR,parallel=TRUE)
    MD = compute_MDLag(MsDeriv)
    
    X = as.matrix(cbind(df$ones,df$y0,xS,xA,xR,xD))
    
    # lag regressors
    MSDelta = as.numeric(MS %*% matrix(df$delta))
    MADelta = as.numeric(MA %*% matrix(df$delta))
    MRDelta = as.numeric(MR %*% matrix(df$delta))
    MDDelta = as.numeric(MD %*% matrix(df$delta))
    
    # instruments for IV
    MS2X=as.matrix(MS %*% MS %*% X)
    MA2X=as.matrix(MA %*% MA %*% X)
    MR2X=as.matrix(MR %*% MR %*% X) 
    MD2X=as.matrix(MD %*% MD %*% X)
    
    IV_est = ivreg(delta ~ y0 + xS + xA + xR + xD + 
               + MSDelta + MADelta + MRDelta + MDDelta  | 
               y0 + xS + xA + xR + xD + MS2X + MA2X + MR2X + MD2X , data=df)
    
    summary(IV_est)
    
    LAR_IV = LogLikAICcR2(df, c(coef(IV_est),0), 10, xS, xA, xR, xD, 
                                                MS, MA, MR, MD, diag(nrow(df)))
    LogLik = LAR_IV$LogLik
    AICc = LAR_IV$AICc
    R2N = LAR_IV$R2Nagelkerke
    
    return(listN(IV_est, LogLik, AICc,R2N))
    
}

