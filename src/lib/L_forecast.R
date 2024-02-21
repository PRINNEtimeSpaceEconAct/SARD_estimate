forcastSARD <- function(NPeriods,SARDCoeff,hA,hR,df,tau=11){
    # df in the same format as data
    # SARDCoeff = tilde(alpha,phi,gammaS,gammaA,gammaR,gammaD,rhoS,rhoA,rhoR,rhoD,lambda)
    # correction = 1/(1-tau*rhoPhi/2)
    
    
    phiTrue = log(sum(df$yT*df$km2)/sum(df$y0*df$km2))/tau 
    correction = SARDCoeff[2]/phiTrue
    
    alpha = SARDCoeff[1]/correction
    phi = SARDCoeff[2]/correction
    gammaS = SARDCoeff[3]/correction
    gammaA = SARDCoeff[4]/correction
    gammaR = SARDCoeff[5]/correction
    gammaD = SARDCoeff[6]/correction
    rhoS = SARDCoeff[7]/correction*2/tau
    rhoA = SARDCoeff[8]/correction*2/tau
    rhoR = SARDCoeff[9]/correction*2/tau
    rhoD = SARDCoeff[10]/correction*2/tau

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