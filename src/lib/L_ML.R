estimate_SARD_auto <- function(df,shp,hA,hR,longlat=TRUE){
    # Estimate SARD via ML with distances hA, hR
    
    if (DEBUG == TRUE){ print("Estimating SARD for given hA hR") }
    
    initJulia()
    
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
    Y = as.matrix(df$delta)
    
    IV_estimator = estimate_IV_SARD_allGiven(df,xS,xA,xR,xD,MS,MA,MR,MD)
    
    WN_SARD_est = estimate_WN_SARD_allGiven(df,xS,xA,xR,xD,MS,MA,MR,MD)
    residualsSARD_WN = WN_SARD_est$residualsSARD_WN
    
    SpatError = compute_spatial_error_mat(residualsSARD_WN,shp)
    Werr = SpatError$Werr
    
    # initial condition of optimizer via WN SARD
    # initial spatial error to 0.5 
    allCoef = WN_SARD_est$coefSARD_WN
    initialOptim = c(allCoef[7:10],0.5)
    
    # call julia optimization
    if (DEBUG == TRUE){ print("Starting Julia optimization for SARD") }
    outSARD_Estimate = call_julia_LogLik(X,Y,MS,MA,MR,MD,Werr,initialOptim)
    coefSARD = outSARD_Estimate$coef
    se_coefSARD = outSARD_Estimate$se_coef
    pvalue_coefSARD = outSARD_Estimate$pvalue_coef
    residualsSARD = outSARD_Estimate$residuals - coefSARD[11]*Werr %*% outSARD_Estimate$residuals
    
    
    LAR_LL = LogLikAICcR2(df, coefSARD, 11, xS, xA, xR, xD, 
                          MS, MA, MR, MD, Werr)
    
    LogLik = LAR_LL$LogLiKelyhood
    AICc = LAR_LL$AICc
    R2N = LAR_LL$R2Nagelkerke
    
    return(listN(coefSARD, se_coefSARD, pvalue_coefSARD,
                 residualsSARD,SpatError,LogLik,AICc,R2N,IV_estimator,WN_SARD_est))
    
}

estimate_WN_SARD_auto <- function(df,hA,hR,longlat=TRUE){
    # Estimate SARD WN via ML with distances hA, hR
    
    if (DEBUG == TRUE){ print("Estimating SARD WN for given hA hR") }
    
    initJulia()
    
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
    Y = as.matrix(df$delta)
    
    # initial condition of optimizer via IV
    IV_estimator = estimate_IV_SARD_allGiven(df,xS,xA,xR,xD,MS,MA,MR,MD)
    IV_est = IV_estimator$IV_est
    initialOptim = c(coef(IV_est)["MSDelta"],coef(IV_est)["MADelta"],
                     coef(IV_est)["MRDelta"],coef(IV_est)["MDDelta"])
    
    # call julia optimization
    if (DEBUG == TRUE){ print("Starting Julia optimization for WN SARD") }
    outSARD_WNEstimate = call_julia_LogLik_WN(X,Y,MS,MA,MR,MD,initialOptim)
    coefSARD_WN = outSARD_WNEstimate$coef
    se_coefSARD_WN = outSARD_WNEstimate$se_coef
    pvalue_coefSARD_WN = outSARD_WNEstimate$pvalue_coef
    residualsSARD_WN = outSARD_WNEstimate$residuals
    
    LAR_LL = LogLikAICcR2(df, c(coefSARD_WN,0), 10, xS, xA, xR, xD, 
                          MS, MA, MR, MD, diag(nrow(df)))
    
    LogLik = LAR_LL$LogLiKelyhood
    AICc = LAR_LL$AICc
    R2N = LAR_LL$R2Nagelkerke
    
    return(listN(coefSARD_WN, se_coefSARD_WN, pvalue_coefSARD_WN,
                 residualsSARD_WN,LogLik,AICc,R2N,IV_est))
    
}


estimate_WN_SARD_allGiven <- function(df,xS,xA,xR,xD,MS,MA,MR,MD){
    # Estimate SARD WN via ML with distances hA, hR
    
    if (DEBUG == TRUE){ print("Estimating SARD WN for given hA hR") }
    
    
    X = as.matrix(cbind(df$ones,df$y0,xS,xA,xR,xD))
    Y = as.matrix(df$delta)
    
    # initial condition of optimizer via IV
    IV_estimator = estimate_IV_SARD_allGiven(df,xS,xA,xR,xD,MS,MA,MR,MD)
    IV_est = IV_estimator$IV_est
    initialOptim = c(coef(IV_est)["MSDelta"],coef(IV_est)["MADelta"],
                     coef(IV_est)["MRDelta"],coef(IV_est)["MDDelta"])
    
    # call julia optimization
    outSARD_WNEstimate = call_julia_LogLik_WN(X,Y,MS,MA,MR,MD,initialOptim)
    coefSARD_WN = outSARD_WNEstimate$coef
    se_coefSARD_WN = outSARD_WNEstimate$se_coef
    pvalue_coefSARD_WN = outSARD_WNEstimate$pvalue_coef
    residualsSARD_WN = outSARD_WNEstimate$residuals
    
    LAR_LL = LogLikAICcR2(df, c(coefSARD_WN,0), 10, xS, xA, xR, xD, 
                          MS, MA, MR, MD, diag(nrow(df)))
    
    LogLik = LAR_LL$LogLiKelyhood
    AICc = LAR_LL$AICc
    R2N = LAR_LL$R2Nagelkerke
    
    return(listN(coefSARD_WN, se_coefSARD_WN, pvalue_coefSARD_WN,
                 residualsSARD_WN,LogLik,AICc,R2N,IV_est))
    
}

compute_spatial_error_mat <- function(residWN,shp,
                                      maxLag = 20, pThreshold = 0.1){
    # decompose the residual into spatial components with different contiguity 
    # rule of thumb: select as maximum contiguity for the spatial error
    # the first contiguity order that is followed by two consecutive 
    # nonsignificant coefficient in the decomposition
    maxLag=min(maxLag, floor(sqrt(length(residWN))))
    
    sink <- capture.output(sf::sf_use_s2(FALSE))
    suppressMessages(spatialNeighbors <- poly2nb(shp,queen=TRUE))
    suppressWarnings(spatialNeighbors.lag <- nblag(spatialNeighbors, maxLag))
    
    maxLag = length(spatialNeighbors.lag)
    
    WAll = list(); 
    WAll[[1]] = nb2mat(spatialNeighbors.lag[[1]],style="B",zero.policy = T)
    WAll[[1]] = as(WAll[[1]],"sparseMatrix")
    WPartial = list(); 
    WPartial[[1]] = nb2mat(spatialNeighbors.lag[[1]],style="B",zero.policy = T) 
    WPartial[[1]] = as(WPartial[[1]],"sparseMatrix")
    errCorPartial = list(); 
    errCorPartial[[1]] = WAll[[1]] %*% residWN
    errCorPartial[[1]] = as.numeric(errCorPartial[[1]])
    
    for (i in 2:maxLag){
        WAll[[i]] = nb2mat(spatialNeighbors.lag[[i]],style="B",zero.policy = T)
        WAll[[i]] = as(WAll[[i]],"sparseMatrix")
        WPartial[[i]] = WAll[[i]] - WAll[[i-1]]
        errCorPartial[[i]] = WPartial[[i]] %*% residWN
        errCorPartial[[i]] = as.numeric(errCorPartial[[i]])
    }
    
    errCorPartialMatrix = matrix(unlist(errCorPartial),ncol=maxLag,byrow=FALSE)
    Mdf = errCorPartialMatrix
    df = data.frame(residWN,Mdf)
    lm.errDecompose = lm(residWN ~ . -1, data = df)
    s.errDecompose = summary(lm.errDecompose)
    pValues = s.errDecompose$coefficients[,"Pr(>|t|)"]
    
    # rule of thumb:
    # magic: index of the first followed by 2 consecutive nonsignificant coefs
    pValuesLarge = pValues > pThreshold
    maxSignifLag = max(Position(function(x) x==TRUE,
                                colSums(rbind(pValuesLarge,c(pValuesLarge[2:maxLag],
                                                             FALSE))) == 2) - 1,1)
    if (is.na(maxSignifLag)) { maxSignifLag = maxLag}    
    
    Werr = s.errDecompose$coefficients[1,"Estimate"]*WPartial[[1]]
    
    if (maxSignifLag == 1){ return(listN(maxSignifLag,Werr,lm.errDecompose,spatialNeighbors.lag)) }
    
    for (i in 2:maxSignifLag){
        Werr = Werr + s.errDecompose$coefficients[i,"Estimate"]*WPartial[[i]]}
    
    return(listN(maxSignifLag,Werr,lm.errDecompose,spatialNeighbors))
}

estimate_SARD_auto_NOCONSTANT <- function(df,shp,hA,hR,longlat=TRUE){
    # Estimate SARD via ML with distances hA, hR
    
    if (DEBUG == TRUE){ print("Estimating SARD for given hA hR") }
    
    initJulia()
    
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
    
    
    X = as.matrix(cbind(df$y0,xS,xA,xR,xD))
    Y = as.matrix(df$delta)
    
    # initial condition of optimizer via IV
    IV_estimator = estimate_IV_SARD_allGiven(df,xS,xA,xR,xD,MS,MA,MR,MD)
    IV_est = IV_estimator$IV_est
    initialOptim = c(coef(IV_est)["MSDelta"],coef(IV_est)["MADelta"],
                     coef(IV_est)["MRDelta"],coef(IV_est)["MDDelta"])
    
    # call julia optimization
    outSARD_WNEstimate = call_julia_LogLik_WN(X,Y,MS,MA,MR,MD,initialOptim)
    coefSARD_WN = outSARD_WNEstimate$coef
    se_coefSARD_WN = outSARD_WNEstimate$se_coef
    pvalue_coefSARD_WN = outSARD_WNEstimate$pvalue_coef
    residualsSARD_WN = outSARD_WNEstimate$residuals
    
    
    SpatError = compute_spatial_error_mat(residualsSARD_WN,shp)
    Werr = SpatError$Werr
    
    # initial condition of optimizer via WN SARD
    # initial spatial error to 0.5 
    allCoef = coefSARD_WN
    initialOptim = c(allCoef[6:9],0.5)
    
    # call julia optimization
    if (DEBUG == TRUE){ print("Starting Julia optimization for SARD") }
    outSARD_Estimate = call_julia_LogLik(X,Y,MS,MA,MR,MD,Werr,initialOptim)
    coefSARD = outSARD_Estimate$coef
    se_coefSARD = outSARD_Estimate$se_coef
    pvalue_coefSARD = outSARD_Estimate$pvalue_coef
    residualsSARD = outSARD_Estimate$residuals - coefSARD[10]*Werr %*% outSARD_Estimate$residuals
    
    
    LAR_LL = LogLikAICcR2(df, coefSARD, 10, xS, xA, xR, xD, 
                          MS, MA, MR, MD, Werr)
    
    LogLik = LAR_LL$LogLiKelyhood
    AICc = LAR_LL$AICc
    R2N = LAR_LL$R2Nagelkerke
    
    return(listN(coefSARD, se_coefSARD, pvalue_coefSARD,
                 residualsSARD,SpatError,LogLik,AICc,R2N,outSARD_Estimate))
    
}
