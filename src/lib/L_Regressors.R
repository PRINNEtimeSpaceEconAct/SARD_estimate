require(Matrix)

compute_xS <- function(df,MsDeriv){
    # compute regressor related to gamma_S
    
    if (DEBUG == TRUE){ print("computing xS") }
    
    Mx = MsDeriv$Mx
    My = MsDeriv$My
    xS = Mx %*% (matrix(df$y0) * (Mx %*% matrix(df$s))) + 
        + My %*% (matrix(df$y0) * (My %*% matrix(df$s)))
    xS = as.numeric(xS)
    
    return(xS)
}

compute_xAR <- function(df,MsDeriv,Wh){
    # compute regressor related to gamma_A or gamma_R, depending on Wh
    
    if (DEBUG == TRUE){ print("computing xAR") }
    
    Mx = MsDeriv$Mx
    My = MsDeriv$My
    xAR = ( Mx %*% (matrix(df$y0) * ((Mx %*% Wh) %*% matrix(df$y0))) + 
        + My %*% (matrix(df$y0) * ((My %*% Wh) %*% matrix(df$y0))) ) 
    xAR = as.numeric(xAR)
    
    return(xAR)
}

compute_xD <- function(df,MsDeriv){
    # compute regressor related to gamma_D
    
    if (DEBUG == TRUE){ print("computing xD") }
    
    Mxx = MsDeriv$Mxx
    Myy = MsDeriv$Myy
    xD = (Mxx + Myy) %*% matrix(df$y0)
    xD = as.numeric(xD)
    
    return(xD)
}

compute_MSLag <- function(df,MsDeriv){
    # compute the lag matrix MS
    
    if (DEBUG == TRUE){ print("computing MSLag") }
    
    Mx = MsDeriv$Mx
    My = MsDeriv$My
    fS <- function(y) as.numeric( Mx %*% (matrix(y) * (Mx %*% matrix(df$s))) + 
                                  + My %*% (matrix(y) * (My %*% matrix(df$s))) )
    MS = apply(diag(nrow(df)), 1, fS); 
    MS = as(MS, "sparseMatrix")
    
    return(MS)
}

compute_MARLag <- function(df,MsDeriv,Wh){
    # compute the lag matrix MA or MR, depending on Wh
    
    if (DEBUG == TRUE){ print("computing MARLag") }
    
    Mx = MsDeriv$Mx
    My = MsDeriv$My
    MxWh = Mx %*% Wh
    MyWh = My %*% Wh
    fAR <- function(delta,y) as.numeric( 
              Mx %*% (matrix(delta) * MxWh %*% matrix(y)) + 
            + Mx %*% (matrix(y) * MxWh %*% matrix(delta)) + 
            + My %*% (matrix(delta) * MyWh %*% matrix(y)) + 
            + My %*% (matrix(y) * MyWh %*% matrix(delta)) )

    if (PARALLEL == FALSE){
        MAR = apply(diag(nrow(df)), 1, fAR, df$y0)
        MAR = as(MAR, "sparseMatrix")}
    if (PARALLEL == TRUE){
        parallelCluster <- snow::makeCluster(NPROCS,type = "SOCK")
        clusterEvalQ(parallelCluster, library(Matrix))
        clusterExport(parallelCluster,c("Mx","My","MxWh","MyWh","df"),
                                                            envir=environment())
        MAR = snow::parApply(parallelCluster,diag(nrow(df)),1,fAR,df$y0);
        stopCluster(parallelCluster)}
    
    MAR = as(MAR, "sparseMatrix")

    return(MAR)
}

compute_MDLag <- function(MsDeriv){
    # compute the lag matrix MD
    
    if (DEBUG == TRUE){ print("computing MDLAg") }
    
    Mxx = MsDeriv$Mxx
    Myy = MsDeriv$Myy
    MD = Mxx + Myy

    return(MD)
}