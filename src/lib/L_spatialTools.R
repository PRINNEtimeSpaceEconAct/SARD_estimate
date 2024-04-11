require(spdep)
require(Matrix)
require(snow)
require(pracma)

compute_D_meanDist <- function(shp,Xi,Yi,fi,NSamplePerPoly = 1000){
    
    # km2 <- st_area(shp)/sum(st_area(shp))
    # y0 = continuous2shp(Xi, Yi, fi, shp)
    
    Np = nrow(shp)
    Xistack = c(Xi)
    Yistack = c(Yi)
    fistack = c(fi)
    
    xSort = sort(unique(Xistack))
    ySort = sort(unique(Xistack))
    fI <- function(x,y) {interp2(xSort,ySort,fi,x,y)}
    
    xy = data.frame(x=Xistack, y=Yistack, f=fistack)
    xys = st_as_sf(xy, coords=c("x","y"))
    st_crs(xys) ="WGS84"
    intersection = st_intersects(shp,xys)

    samplePoly = list()
    for (indexPoly in 1:Np){
        poly = shp$geom[[indexPoly]]
        maxDensity = max(fistack[intersection[[indexPoly]]])
        samplePoly[[indexPoly]] = sampleFromDensity(fI,maxDensity,NSamplePerPoly,shp$geom[[indexPoly]])
    }
    
    # plot(poly)
    # points(Xistack[intersection[[indexPoly]]],Yistack[intersection[[indexPoly]]])
    # points(sampleFromDensity(fI,maxDensity,NSamplePerPoly,shp$geom[[indexPoly]]),col="red")
    
    D = matrix(0,nrow=Np,ncol=Np)
    for (i in 1:Np){
        for (j in 1:Np){
            if (i != j) {
            sample_i = samplePoly[[i]]
            sample_j = samplePoly[[j]]
            D[i,j] = mean((crossdist.default(sample_i[,1],sample_i[,2],sample_j[,1],sample_j[,2],period=c(1,1)))) }
        }
    }

    return(D)
}


sampleFromDensity <- function(densityFun,maxDensity,N,poly){
    samples = suppressWarnings(matrix(data=c(0,0),nrow=0,ncol=2))
    while (nrow(samples) < N){
        Nmissing = N - nrow(samples)
        x = st_coordinates(st_sample(poly,Nmissing))
        y = runif(Nmissing,min=0,max=maxDensity)
        acceptedSamples = x[(y <= densityFun(x[,1],x[,2])),]
        samples = rbind(samples,acceptedSamples)
        samples = samples[!(is.na(samples[,1]) | is.na(samples[,2])),]
    }
    samples = samples[1:N,]
    return(samples)
}


compute_D <- function(df,dMax=100, longlat=TRUE, torus=FALSE){
    # dMax = 100Km by default. Assumes that we will never consider interaction
    # at distance higher than 100Km.
    # longlat = TRUE by default. Assume we want to consider great circle
    # distance.
    # In DataFrame data we assume that Latitude Longitude are in degrees with
    # projection "+proj=longlat  +datum=WGS84  +no_defs"
    
    if (DEBUG == TRUE){ print("computing all distances") }
    
    if (torus == FALSE){
        coord = cbind(df$Longitude,df$Latitude)
        spatialNeighbors <- dnearneigh(coord, 0,dMax, row.names=NULL, longlat=longlat)
        dlist <- nbdists(spatialNeighbors, coord, longlat=longlat)
        dlist1 <- lapply(dlist, function(x) x)          
        spatialNeighbors <- suppressWarnings(nb2listw(spatialNeighbors,
                                                      glist=dlist1, style="B", zero.policy=TRUE))
        D <- listw2mat(spatialNeighbors)
        D <- as(D, "sparseMatrix")
    }
    else {
        # Dx = as.matrix(stats::dist(df$Longitude, diag = TRUE, upper = TRUE))
        # Dy = as.matrix(stats::dist(df$Latitude, diag = TRUE, upper = TRUE))
        # D = sqrt( pmin( Dx,1-Dx )^2 + pmin( Dy,1-Dy )^2 )
        D = crossdist.default(df$Longitude,df$Latitude,df$Longitude,df$Latitude,period=c(1,1))
    }
    return(D)
    
}

GFDM <- function(df,torus=FALSE){
    # starting from coordinates returns a list with all sparse matrices
    # Mx,My,Mxx,Myy,Mxy
    
    if (DEBUG == TRUE){ print("computing derivative matrices") }
    
    coord = cbind(df$Longitude,df$Latitude)
    MsDeriv = compute_MDiff(coord,torus=torus)
    return(MsDeriv)
}

compute_WhAR <- function(D,df,h){
    # weight matrices WhA or WhR
    
    if (DEBUG == TRUE){ print("computing WhAR") }

    dInvSq <- function(d,h){ 1/(2*pi*(log(2)-1/2)) * 1/h^2 * 1/(d/h+1)^2 * ((d <= h) & (d > 0)) }

    Wh = dInvSq(D,h)
    Wh[is.na(Wh)] = 0
    diag(Wh) = 1/(2*pi*(log(2)-1/2)) * 1/h^2
    
    Wh = t(apply(Wh, 1, function(x) x * as.numeric(df$km2)))
    Wh = Wh / sum((Wh%*%df$y0)*df$km2)
    return(Wh)
}

compute_Wh <- function(D,df,h){
    # weight matrices for spatial durbin
    
    if (DEBUG == TRUE){ print("computing W for spatial Durbin") }
    
    dInvSq <- function(d,cut) 1/d^2 * ((d <= cut) & (d > 0))
    
    W = dInvSq(D,h)
    W[is.na(W)] = 0
    W = W/rowSums(W)
    W[is.na(W)] = 0
    W = as(W,"sparseMatrix")

    return(W)
}

LogLikAICcR2 <- function(df, coef, k, xS, xA, xR, xD, MS, MA, MR, MD, W_eps){
    # compute LogLik and AICc of the SARD Model, with dof degree of freedom
    # coefs is of length 11, in order
    # a, phi, gammaS, gammA, gammaR, gammaD, rhoS, rhoA, rhoR, rhoD, lambda.
    # used also for NAIVE, IV WN, SARD WN. 
    
    if (DEBUG == TRUE){ print("computing LogLik") }
    
    # thankyou R for being so modern    
    a=coef[1]; phi=coef[2]; gammaS=coef[3]; gammA=coef[4]; gammaR = coef[5];
    gammaD=coef[6]; rhoS=coef[7]; rhoA=coef[8]; rhoR=coef[9]; rhoD=coef[10]; 
    lambda=coef[11]

    N = nrow(MS)
    Y = df$delta
    X = cbind(df$ones,df$y0,xS,xA,xR,xD)
    colnames(X) <- c("ones","y0","xS","xA","xR","xD")
    
    A = diag(N) - rhoS*MS - rhoA*MA - rhoR*MR - rhoD*MD
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

LogLikelyhood2AICc <- function(LogLiKelyhood,k,N){
    AIC = 2*k - 2*LogLiKelyhood
    AICc = AIC + (2*k^2+2*k)/(N-k-1)
    return(AICc)
}


AICc2R2Nagelkerke <- function(aicc,Y,k){
    # from AIC corrected to R^2 Nagelkerke
    
    N = length(Y)
    
    # null model
    LL0 = logLik(lm(Y ~ 1))
    LogLiKelyhood = (2*k*N/(N-k-1) - aicc)/2
    
    R2Nagelkerke =  c(1 - exp(-(2/N)*(LogLiKelyhood - LL0)))
    
    return(R2Nagelkerke)
}


listN <- function(...){
    # automatically give names to list elements = var name
    anonList <- list(...)
    names(anonList) <- as.character(substitute(list(...)))[-1]
    anonList
}

correlogram <- function(resid,shp,maxLag = 20){
    # make the spatial correlogram with contiguity matrices
    
    sf_use_s2(FALSE)
    spatialNeighbors <- poly2nb(shp)
    correlogram_resid = sp.correlogram(spatialNeighbors,as.numeric(resid),
                    order = maxLag, style="B", method = "I",zero.policy = T)
    
    plot(correlogram_resid,main="")
    
    return(listN(correlogram_resid))
}
