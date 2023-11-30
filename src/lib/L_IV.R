require(AER)

chose_hAhR <- function(df,hA_range,hR_range,longlat=TRUE) {
    # select best hA and hR based on the AICc score
    
    print("            selecting best hA hR via IV")
    
    MsDeriv = GFDM(df)
    D = compute_D(df,longlat=longlat)
    
    xS = compute_xS(df,MsDeriv)
    xD = compute_xD(df,MsDeriv)
    MS = compute_MSLag(df,MsDeriv)
    MD = compute_MDLag(MsDeriv)
    
    allPairs = expand.grid(hA=hA_range,hR=hR_range)
    allAICc = list()
    allIV_est = list()
    for (i in 1:nrow(allPairs)) {
        
        hA_i = allPairs[i,]$hA
        hR_i = allPairs[i,]$hR
        
        print(paste("    Estimating IV (",i,"/",nrow(allPairs),"), with hA = ",
                                            hA_i,", and hR = ", hR_i ,sep=""))
        
        outIVEstimate_i = estimate_IV_SARD_hAhRGiven(df,MsDeriv,xS,xD,MS,MD,D,hA_i,hR_i)
        allAICc[[i]] = outIVEstimate_i$AICc
        allIV_est[[i]] = outIVEstimate_i$IV_est
    }
    
    iBest = which.min(allAICc)
    
    hABest = allPairs[iBest,]$hA
    hRBest = allPairs[iBest,]$hR
    IV_est = allIV_est[[iBest]]

    return (listN(hABest,hRBest,IV_est,allPairs,allAICc,allIV_est))

}

estimate_IV_SARD_allGiven <- function(df,xS,xA,xR,xD,MS,MA,MR,MD) {
    # Estimate SARD WN via IV with distances hA, hR
    
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
    

    LAR_IV = LogLikAICcR2(df, c(coef(IV_est),0), 10, xS, xA, xR, xD, 
                          MS, MA, MR, MD, diag(nrow(df)))
    LogLik = LAR_IV$LogLiKelyhood
    AICc = LAR_IV$AICc
    R2N = LAR_IV$R2Nagelkerke
    
    return(listN(IV_est, LogLik, AICc,R2N))
    
}

estimate_IV_SARD_hAhRGiven <- function(df,MsDeriv,xS,xD,MS,MD,D,hA,hR) {
    # Estimate SARD WN via IV with distances hA, hR

    
    WhA = compute_WhAR(D,df,hA)
    WhR = compute_WhAR(D,df,hR)
    xA = compute_xAR(df,MsDeriv, WhA)
    xR = compute_xAR(df,MsDeriv, WhR)
    MA = compute_MARLag(df,MsDeriv,WhA)
    MR = compute_MARLag(df,MsDeriv,WhR)

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
    

    LAR_IV = LogLikAICcR2(df, c(coef(IV_est),0), 10, xS, xA, xR, xD, 
                          MS, MA, MR, MD, diag(nrow(df)))
    LogLik = LAR_IV$LogLiKelyhood
    AICc = LAR_IV$AICc
    R2N = LAR_IV$R2Nagelkerke
    
    return(listN(IV_est, LogLik, AICc,R2N))
    
}

estimate_IV_SARD_auto <- function(df,hA,hR,longlat=TRUE){
    # Estimate SARD WN via IV with distances hA, hR

    print("Estimating IV for given hA hR") 
    
    MsDeriv = GFDM(df)
    D = compute_D(df,longlat=longlat)
    
    WhA = compute_WhAR(D,df,hA)
    WhR = compute_WhAR(D,df,hR)
    xS = compute_xS(df,MsDeriv)
    xA = compute_xAR(df,MsDeriv, WhA)
    xR = compute_xAR(df,MsDeriv, WhR)
    xD = compute_xD(df,MsDeriv)
    MS = compute_MSLag(df,MsDeriv)
    MA = compute_MARLag(df,MsDeriv,WhA)
    MR = compute_MARLag(df,MsDeriv,WhR)
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
    

    LAR_IV = LogLikAICcR2(df, c(coef(IV_est),0), 10, xS, xA, xR, xD, 
                                                MS, MA, MR, MD, diag(nrow(df)))
    LogLik = LAR_IV$LogLiKelyhood
    AICc = LAR_IV$AICc
    R2N = LAR_IV$R2Nagelkerke
    
    return(listN(IV_est, LogLik, AICc,R2N))
    
}

