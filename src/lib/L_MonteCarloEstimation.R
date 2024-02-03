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


estimate_ARD_MC_LM <- function(df,shp,xA,xR,xD,MA,MR,MD){
    
    # lag regressors
    MADelta = as.numeric(MA %*% matrix(df$delta))
    MRDelta = as.numeric(MR %*% matrix(df$delta))
    MDDelta = as.numeric(MD %*% matrix(df$delta))
    
    df = df %>% mutate(xA = xA, xR = xR, xD = xD,
                       MADelta=MADelta, MRDelta=MRDelta,MDDelta=MDDelta)
    
    LM_est = lm(delta ~ 0 + xA + MADelta + xR + MRDelta + xD + MDDelta, data = df)
    return(LM_est)
}


estimate_ARD_MC_LL <- function(df,shp,xA,xR,xD,MA,MR,MD){
    
    Y = as.matrix(df$delta)
    X = data.frame(xA=xA,xR=xR,xD=xD)
    X = as.matrix(X)

    coefLM = coef(estimate_ARD_MC_LM(df,shp,xA,xR,xD,MA,MR,MD))
    initialOptim = coefLM[c(2,4,6)]
    outARD_WN3MatEstimate = call_julia_LogLik_WN_3Mat(X,Y,MA,MR,MD,initialOptim)
    
    SpatError = compute_spatial_error_mat(outARD_WN3MatEstimate$residuals,shp,maxLag = 10)
    Werr = SpatError$Werr
    initialOptim = c(outARD_WN3MatEstimate$coef[c(4,5,6)],0.5)  
    outARD_3MatEstimate = call_julia_LogLik_3Mat(X,Y,MA,MR,MD,Werr,initialOptim)
    
    return(listN(outARD_WN3MatEstimate,outARD_3MatEstimate,SpatError))
}

