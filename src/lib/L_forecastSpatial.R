forcastSpatial <- function(NPeriods,DURBIN,hDurbin,df,tau=11){
    # df in the same format as data
    # SARDCoeff = tilde(alpha,phi,gammaS,gammaA,gammaR,gammaD,rhoS,rhoA,rhoR,rhoD,lambda)
    # correction = 1/(1-tau*rhoPhi/2)
    
    
    D = compute_D(df,longlat=TRUE)
    W = compute_Wh(D,df,hDurbin) 
    Id <- as(diag(nrow(df)),"sparseMatrix")
    
    alpha = DURBIN$coefficients["(Intercept)"]
    phi = DURBIN$coefficients["y0"]
    theta = phi = DURBIN$coefficients["lag.y0"]
    gammaAlt = DURBIN$coefficients["s"]
    thetaAlt = phi = DURBIN$coefficients["lag.s"]
    rho = DURBIN$rho
    
    s = DURBIN$X[,"s"]
    Ws = DURBIN$X[,"lag.s"]
    A = Id - rho*W
    
    y = df$yT
    for (t in 1:NPeriods){
        print(t)
        
        yTmp = alpha + phi*y + gammaAlt * s + theta*(W %*% y) + thetaAlt * Ws 
        y = y + solve(A,yTmp)
        y[y<0] = 0
    }
    
    return(as.numeric(y))
}