forcastSARD <- function(NPeriods,SARDCoeff,hA,hR,df,tau=11){
    # df in the same format as data
    # SARDCoeff = tilde(alpha,phi,gammaS,gammaA,gammaR,gammaD,rhoS,rhoA,rhoR,rhoD,lambda)
    # correction = 1/(1-tau*rhoPhi/2)
    
    yT = sum(df$yT*df$km2)
    y0 = sum(df$y0*df$km2)
    rhoPhi = (2/tau)*(1-1/(tau*SARDCoeff[2])*log((yT+SARDCoeff[1]/SARDCoeff[2])/(y0+SARDCoeff[1]/SARDCoeff[2])))
    
    alpha = SARDCoeff[1]*(1-tau/2*rhoPhi)
    phi = SARDCoeff[2]*(1-tau/2*rhoPhi)
    gammaS = SARDCoeff[3]*(1-tau/2*rhoPhi)
    gammaA = SARDCoeff[4]*(1-tau/2*rhoPhi)
    gammaR = SARDCoeff[5]*(1-tau/2*rhoPhi)
    gammaD = SARDCoeff[6]*(1-tau/2*rhoPhi)
    rhoS = SARDCoeff[7]*(1-tau/2*rhoPhi)*2/tau
    rhoA = SARDCoeff[8]*(1-tau/2*rhoPhi)*2/tau
    rhoR = SARDCoeff[9]*(1-tau/2*rhoPhi)*2/tau
    rhoD = SARDCoeff[10]*(1-tau/2*rhoPhi)*2/tau
    
    MsDeriv = GFDM(df)
    D = compute_D(df)
    WhA = compute_WhAR(D,df,hA)
    WhR = compute_WhAR(D,df,hR)
    
    MxWA = MsDeriv$Mx %*% as(WhA,"sparseMatrix")
    MyWA = MsDeriv$My %*% as(WhA,"sparseMatrix")
    MxWR = MsDeriv$Mx %*% as(WhR,"sparseMatrix")
    MyWR = MsDeriv$My %*% as(WhR,"sparseMatrix")
    
    MS = compute_MSLag(df,MsDeriv)
    MD = compute_MDLag(MsDeriv)
    
    Id <- as(diag(nrow(df)),"sparseMatrix")
    
    y = df$yT
    for (t in 1:NPeriods){
        print(t)
        
        MA = compute_MARLag(df,MsDeriv,WhA)
        MR = compute_MARLag(df,MsDeriv,WhR)
        yPre = y
        yTmp = alpha + phi*y + 
            + gammaS * (MsDeriv$Mx %*% (y * MsDeriv$Mx %*% df$s) + MsDeriv$My %*% (y * MsDeriv$My %*% df$s)) + 
            + gammaA * (MsDeriv$Mx %*% (y * MxWA %*% y) + MsDeriv$My %*% (y * MyWA %*% y)) + 
            + gammaR * (MsDeriv$Mx %*% (y * MxWR %*% y) + MsDeriv$My %*% (y * MyWR %*% y)) +
            + gammaD * (MD %*% y) 
        A = Id - rhoS*MS - rhoA*MA - rhoR*MR - rhoD*MD
        y = y + solve(A,yTmp)
        y[y<0] = yPre[matrix(y)<0]*exp(-(yPre[matrix(y)<0]-y[matrix(y)<0])/yPre[matrix(y)<0])
    }
    
    return(as.numeric(y))
}

forcastSARD_NOCONSTANT <- function(NPeriods,SARDCoeff,hA,hR,df,tau=11){
    # df in the same format as data
    # SARDCoeff = tilde(alpha,phi,gammaS,gammaA,gammaR,gammaD,rhoS,rhoA,rhoR,rhoD,lambda)
    # correction = 1/(1-tau*rhoPhi/2)
    
    yT = sum(df$yT*df$km2)
    y0 = sum(df$y0*df$km2)
    rhoPhi = (2/tau)*(1-1/(tau*SARDCoeff[1])*log((yT+0/SARDCoeff[1])/(y0+0/SARDCoeff[1])))
    
    phi = SARDCoeff[1]*(1-tau/2*rhoPhi)
    gammaS = SARDCoeff[2]*(1-tau/2*rhoPhi)
    gammaA = SARDCoeff[3]*(1-tau/2*rhoPhi)
    gammaR = SARDCoeff[4]*(1-tau/2*rhoPhi)
    gammaD = SARDCoeff[5]*(1-tau/2*rhoPhi)
    rhoS = SARDCoeff[6]*(1-tau/2*rhoPhi)*2/tau
    rhoA = SARDCoeff[7]*(1-tau/2*rhoPhi)*2/tau
    rhoR = SARDCoeff[8]*(1-tau/2*rhoPhi)*2/tau
    rhoD = SARDCoeff[9]*(1-tau/2*rhoPhi)*2/tau
    
    MsDeriv = GFDM(df)
    D = compute_D(df)
    WhA = compute_WhAR(D,df,hA)
    WhR = compute_WhAR(D,df,hR)
    
    MxWA = MsDeriv$Mx %*% as(WhA,"sparseMatrix")
    MyWA = MsDeriv$My %*% as(WhA,"sparseMatrix")
    MxWR = MsDeriv$Mx %*% as(WhR,"sparseMatrix")
    MyWR = MsDeriv$My %*% as(WhR,"sparseMatrix")
    
    MS = compute_MSLag(df,MsDeriv)
    MD = compute_MDLag(MsDeriv)
    
    Id <- as(diag(nrow(df)),"sparseMatrix")
    
    y = df$yT
    yOut = matrix(NA, nrow = nrow(df), ncol = (NPeriods+1))
    yOut[,1] = y
    for (t in 1:NPeriods){
        print(t)
        
        MA = compute_MARLag(df,MsDeriv,WhA)
        MR = compute_MARLag(df,MsDeriv,WhR)
        yPre = y
        yTmp = phi*y + 
            + gammaS * (MsDeriv$Mx %*% (y * MsDeriv$Mx %*% df$s) + MsDeriv$My %*% (y * MsDeriv$My %*% df$s)) + 
            + gammaA * (MsDeriv$Mx %*% (y * MxWA %*% y) + MsDeriv$My %*% (y * MyWA %*% y)) + 
            + gammaR * (MsDeriv$Mx %*% (y * MxWR %*% y) + MsDeriv$My %*% (y * MyWR %*% y)) +
            + gammaD * (MD %*% y) 
        A = Id - rhoS*MS - rhoA*MA - rhoR*MR - rhoD*MD
        y = y + solve(A,yTmp)
        y[y<0] = yPre[matrix(y)<0]*exp(-(yPre[matrix(y)<0]-y[matrix(y)<0])/yPre[matrix(y)<0])
        
        yOut[,t+1] = as.numeric(y)
    }
    
    return(yOut)
}