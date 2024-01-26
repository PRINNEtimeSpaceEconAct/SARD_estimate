# MonteCarlo <- function(Np, Na, Nm, tau, typeOfDist = "uniform", SARDp, hA_range, 
#                        hR_range, NeS = 1000){
#     # add documentation of the function
#     # INPUT
#     # OUTPUT
#     
#     #Create shape file
#     shpMC = createShape(Np, typeOfDist = typeOfDist, plot=FALSE)
#     
#     #Create agents
#     AgentsAll = call_julia_computeAgents(Nm,Na,tau,SARDp)
#     agents0_MC = AgentsAll$agents0_MC
#     agentsT_MC = AgentsAll$agentsT_MC
#     
#     #For s
#     SComputed = call_julia_computeS(NeS)
#     X = SComputed$X
#     Y = SComputed$Y
#     S = SComputed$S
#     
#     outSARDEstimateMonteCarlo = list()
#     for (m in 1:Nm){
#         agents0 = agents0_MC[m,,]
#         agentsT = agentsT_MC[m,,]
#         
#         data_shp = createDataframe(agents0, agentsT, shpMC, tau, X, Y, S)
#         # plot(data_shp$shp_sf)
#         data = data_shp$data
#         
#         select_h=chose_hAhR(data,hA_range,hR_range,longlat=FALSE)
#         outSARDEstimate = estimate_SARD_auto(data,shpMC,select_h$hABest,
#                                              select_h$hRBest,longlat=FALSE)
#         
#         outSARDEstimate$hABest = select_h$hABest
#         outSARDEstimate$hRBest = select_h$hRBest
#         
#         outSARDEstimateMonteCarlo[[m]] = outSARDEstimate
#     }
#     
#     return(outSARDEstimateMonteCarlo)
#     
# }
# 
# MonteCarloFixedhAhR_IV <- function(Np, Na, Nm, tau, typeOfDist = "uniform", SARDp,hA,hR,
#                                    NeS = 100,model=list()){
#     # add documentation of the function
#     # INPUT
#     # OUTPUT
#     
#     #Create shape file
#     shpMC = createShape(Np, typeOfDist = typeOfDist, plot=FALSE)
#     
#     #Create agents
#     AgentsAll = call_julia_computeAgents(Nm,Na,tau,SARDp)
#     agents0_MC = AgentsAll$agents0_MC
#     agentsT_MC = AgentsAll$agentsT_MC
#     
#     #For s
#     SComputed = call_julia_computeS(NeS)
#     X = SComputed$X
#     Y = SComputed$Y
#     S = SComputed$S
#     
#     outIVEstimateMonteCarlo = list()
#     for (m in 1:Nm){
#         agents0 = agents0_MC[m,,]
#         agentsT = agentsT_MC[m,,]
#         
#         data_shp = createDataframe(agents0, agentsT, shpMC, tau, X, Y, S)
#         data = data_shp$data
#         
#         outIVEstimate = estimate_IV_SARD_autoMC(data,hA,hR,longlat = FALSE,model=model)
#         
#         outIVEstimateMonteCarlo[[m]] = outIVEstimate
#     }
#     
#     return(outIVEstimateMonteCarlo)
#     
# }

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

MonteCarloOneRun_LM_Agents_fixedhAhR <- function(Np, Na, tau, typeOfDist = "uniform", SARDp,hA,hR,
                                   NeS = 100,torus=TRUE){
    # add documentation of the function
    # INPUT
    # OUTPUT
    allVar = createModelVariables(SARDp)
    model = allVar$model
    variables = allVar$variables
    
    
    #Create shape file
    shpMC = createShape(Np, typeOfDist = typeOfDist, plot=FALSE)
    
    #Create agents
    AgentsAll = call_julia_computeAgents(1,Na,tau,SARDp)
    agents0_MC = AgentsAll$agents0_MC
    agentsT_MC = AgentsAll$agentsT_MC
    
    #For s
    SComputed = call_julia_computeS(NeS)
    Xs = SComputed$X
    Ys = SComputed$Y
    S = SComputed$S
    
        agents0 = agents0_MC[1,,]
        agentsT = agentsT_MC[1,,]
        
        data_shp = createDataframeAgents(agents0, agentsT, shpMC, tau, Xs, Ys, S)
        data = data_shp$data
        shp = data_shp$shp_sf

        outLMEstimate = estimate_LM_SARD_autoMC(data,hA,hR,shp,longlat = FALSE,model=model,variables=variables,torus=torus)
        plot(select(outLMEstimate$shp_regressors,-c("km2")))
        
    return(listN(outLMEstimate,data,shp))
    
}

MonteCarloOneRun_LM_PDE_fixedhAhR <- function(Np, tau, typeOfDist = "VoronoiUniform", typeOfEst = "LM", SARDp,hA,hR,
                                                 NeS = 100,torus=TRUE){
    # typeOfEst = "LM", "WNLL", "SELL"
    # add documentation of the function
    # INPUT
    # OUTPUT
    
    allVar = createModelVariables(SARDp)
    model = allVar$model
    variables = allVar$variables
    
    #Create shape file
    shpMC = createShape(Np, typeOfDist = typeOfDist, plot=FALSE)
    
    #Compute PDE
    PDEAll = call_julia_computePDE(tau,SARDp)
    Xpde = PDEAll$X
    Ypde = PDEAll$Y
    PDE0 = PDEAll$PDE0
    PDET = PDEAll$PDET
    
    #For s
    SComputed = call_julia_computeS(NeS)
    Xs = SComputed$X
    Ys = SComputed$Y
    S = SComputed$S
    
    data_shp = createDataframePDE(PDE0, PDET, Xpde, Ypde, shpMC, tau, Xs, Ys, S)
    data = data_shp$data
    shp = data_shp$shp_sf
    
    if (typeOfEst == "LM"){ 
        outEstimate = estimate_LM_SARD_autoMC(data,hA,hR,shp,longlat = FALSE,model=model,variables=variables,torus=torus) }
    if (typeOfEst == "WNLL"){
        outEstimate = estimate_WN1Mat_SARD_auto_MC(data,hA,hR,shp,longlat = FALSE,model=model,variables=variables,torus=torus)
    }
    if (typeOfEst == "SELL"){
        outEstimate = estimate_1Mat_SARD_auto_MC(data,hA,hR,shp,longlat = FALSE,model=model,variables=variables,torus=torus)
    }
        
    plot(select(outEstimate$shp_regressors,-c("km2")))

    return(listN(outEstimate,data,shp))
    
}

MonteCarloOneRun_LM_PDE_fixedhAhR_Voronoi2Unif <- function(Np, tau, SARDp,hA,hR,
                                       typeOfEst = "LM",NeS = 100,torus=TRUE){
    
    allVar = createModelVariables(SARDp)
    model = allVar$model
    variables = allVar$variables
    
    #Create shape file
    shpMC = createShape(Np, typeOfDist = "VoronoiUniform", plot=FALSE)
    
    #Compute PDE
    PDEAll = call_julia_computePDE(tau,SARDp)
    Xpde = PDEAll$X
    Ypde = PDEAll$Y
    PDE0 = PDEAll$PDE0
    PDET = PDEAll$PDET
    
    #For s
    SComputed = call_julia_computeS(NeS)
    Xs = SComputed$X
    Ys = SComputed$Y
    S = SComputed$S
    
    data_shp = createDataframePDE(PDE0, PDET, Xpde, Ypde, shpMC, tau, Xs, Ys, S)
    data = data_shp$data
    shp = data_shp$shp_sf
    
    square = st_polygon(list(cbind(c(0,1,1,0,0), c(0,0,1,1,0))))
    squareBox = st_bbox(c(xmin = 0, xmax = 1, ymax = 1, ymin = 0))
    squareSfc = st_sfc(square,crs = "WGS84")
    grid_st = st_make_grid(squareSfc,cellsize = 1/sqrt(Np))
    grid_sf = st_sf(geom=grid_st,crs = "WGS84")
    grid_sf = grid_sf %>% mutate(Id = 1:nrow(grid_sf))
    attr(st_geometry(grid_sf), "bbox") = squareBox
    coordinatesGrid = suppressWarnings(st_coordinates(st_centroid(grid_sf,of_largest_polygon = TRUE)))
    
    voronoi2grid_shp = suppressWarnings(st_interpolate_aw(select(shp,c(y0,yT,s,)),grid_sf,extensive=F))
    voronoi2grid_shp = voronoi2grid_shp %>% mutate(Id = as.character(1:nrow(voronoi2grid_shp)))
    voronoi2grid_shp = voronoi2grid_shp %>% mutate(delta = (yT-y0)/tau, ones = 1, 
                                km2 = as.numeric(st_area(grid_sf)/sum(st_area(grid_sf))),
                                Longitude = coordinatesGrid[,1],Latitude = coordinatesGrid[,2]) %>%
                                rename(geom = geometry)
    voronoi2grid_shp = voronoi2grid_shp[,colnames(shp)]
    
    shp = voronoi2grid_shp
    data = select(as.data.frame(voronoi2grid_shp),-c("geom")) %>% rename(geo = Id)
    
    if (typeOfEst == "LM"){ 
        outEstimate = estimate_LM_SARD_autoMC(data,hA,hR,shp,longlat = FALSE,model=model,variables=variables,torus=torus) }
    if (typeOfEst == "WNLL"){
        outEstimate = estimate_WN1Mat_SARD_auto_MC(data,hA,hR,shp,longlat = FALSE,model=model,variables=variables,torus=torus)
    }
    if (typeOfEst == "SELL"){
        outEstimate = estimate_1Mat_SARD_auto_MC(data,hA,hR,shp,longlat = FALSE,model=model,variables=variables,torus=torus)
    }
    
    
    plot(select(outEstimate$shp_regressors,-c("km2")))
    
    return(listN(outEstimate,data,shp))
    
}

estimate_LM_SARD_autoMC <- function(df,hA,hR,shp,longlat=FALSE,model=list(),variables,torus=TRUE){
    # Estimate SARD WN via IV with distances hA, hR
    # model = list(regressors=as.formula(delta ~ y0 + xS + xA + xR + xD + MSDelta + MADelta + MRDelta + MDDelta),
                 # instruments=as.formula(~ y0 + xS + xA + xR + xD + MS2X + MA2X + MR2X + MD2X))
    # variables = c("ones","y0","xS","xA","xR","xD")
    
    if (DEBUG == TRUE){ print("Estimating LM for given hA hR") }
  
    MsDeriv = GFDM(df,torus=torus)
    D = compute_D(df,longlat=longlat,torus=torus)
    
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
    
    X = data.frame(ones=df$ones,y0=df$y0,xS=xS,xA=xA,xR=xR,xD=xD)
    X = X %>% select(all_of(variables))
    X = as.matrix(X)
    
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
    
    df = df %>% mutate(xS=xS, xA = xA, xR = xR, xD = xD,
               MSDelta=MSDelta, MADelta=MADelta, MRDelta=MRDelta,MDDelta=MDDelta
               ,MS2X=MS2X[,],MA2X=MA2X[,],MR2X=MR2X[,],MD2X=MD2X[,])
    
    shp_regressors = shp %>% left_join(df, by=c("Id"="geo"),keep=FALSE)
    
    if (DEBUG == TRUE){ print("estimating with LM") }
    # LM_est = ivreg(model$regressors, model$instruments, data=df)
    
    LM_est = lm(model$regressors, data = df)
    
    LAR_LM = NULL
    LogLik = LAR_LM$LogLiKelyhood
    AICc = LAR_LM$AICc
    R2N = LAR_LM$R2Nagelkerke
    
    shp_regressors = shp_regressors %>% select(c("y0"="y0.x","yT"="yT.x","delta"="delta.x","km2"="km2.x","xS","xA","xR","xD","MSDelta","MADelta","MRDelta","MDDelta")) %>% mutate("WhAy0"=WhA%*%df$y0,"WhRy0"=WhR%*%df$y0)
    
    return(listN(LM_est, LogLik, AICc,R2N,shp_regressors,MsDeriv))
    
}

# estimate_WN_SARD_auto <- function(df,hA,hR,longlat=TRUE){
estimate_WN1Mat_SARD_auto_MC <- function(df,hA,hR,shp,longlat=FALSE,model=list(),variables,torus=TRUE){
    
    # Estimate SARD WN via ML with distances hA, hR
    
    if (DEBUG == TRUE){ print("Estimating SARD WN with one regressor for given hA hR") }
    if (length(variables) > 2){ print("Error, this funcitions admits only one regressor")
                                break }

    MsDeriv = GFDM(df,torus=torus)
    D = compute_D(df,longlat=longlat,torus=torus)    
    
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
    
    Y = as.matrix(df$delta)
    X = data.frame(ones=df$ones,y0=df$y0,xS=xS,xA=xA,xR=xR,xD=xD)
    X = X %>% select(all_of(variables))
    X = as.matrix(X)
    
    # lag regressors
    MSDelta = as.numeric(MS %*% matrix(df$delta))
    MADelta = as.numeric(MA %*% matrix(df$delta))
    MRDelta = as.numeric(MR %*% matrix(df$delta))
    MDDelta = as.numeric(MD %*% matrix(df$delta))
    
 
    df = df %>% mutate(xS=xS, xA = xA, xR = xR, xD = xD,
                       MSDelta=MSDelta, MADelta=MADelta, MRDelta=MRDelta,MDDelta=MDDelta)
    
    shp_regressors = shp %>% left_join(df, by=c("Id"="geo"),keep=FALSE)
    

    LM_est = lm(model$regressors, data = df)
    initialOptim = coef(LM_est)[2]

    # call julia optimization
    if (variables == "xS"){ W1 = MS}
    if (variables == "xA"){ W1 = MA}
    if (variables == "xR"){ W1 = MR}
    if (variables == "xD"){ W1 = MD}
    
    if (DEBUG == TRUE){ print("Starting Julia optimization for WN SARD 1Mat") }
    outSARD_WN1MatEstimate = call_julia_LogLik_WN_1Mat(X,Y,W1,initialOptim)

    LAR_LM = NULL
    LogLik = LAR_LM$LogLiKelyhood
    AICc = LAR_LM$AICc
    R2N = LAR_LM$R2Nagelkerke
    
    shp_regressors = shp_regressors %>% select(c("y0"="y0.x","yT"="yT.x","delta"="delta.x","km2"="km2.x","xS","xA","xR","xD","MSDelta","MADelta","MRDelta","MDDelta")) %>% mutate("WhAy0"=WhA%*%df$y0,"WhRy0"=WhR%*%df$y0)
    
    outEstimate = listN(outSARD_WN1MatEstimate, LogLik, AICc,R2N,shp_regressors,MsDeriv)
    return(outEstimate)
    
}

estimate_1Mat_SARD_auto_MC <- function(df,hA,hR,shp,longlat=FALSE,model=list(),variables,torus=TRUE){
    
    # Estimate SARD via ML with distances hA, hR
    
    if (DEBUG == TRUE){ print("Estimating SARD WN with one regressor for given hA hR") }
    if (length(variables) > 2){ print("Error, this funcitions admits only one regressor")
        break }
    
    MsDeriv = GFDM(df,torus=torus)
    D = compute_D(df,longlat=longlat,torus=torus)    
    
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
    
    Y = as.matrix(df$delta)
    X = data.frame(ones=df$ones,y0=df$y0,xS=xS,xA=xA,xR=xR,xD=xD)
    X = X %>% select(all_of(variables))
    X = as.matrix(X)
    
    # lag regressors
    MSDelta = as.numeric(MS %*% matrix(df$delta))
    MADelta = as.numeric(MA %*% matrix(df$delta))
    MRDelta = as.numeric(MR %*% matrix(df$delta))
    MDDelta = as.numeric(MD %*% matrix(df$delta))
    
    
    df = df %>% mutate(xS=xS, xA = xA, xR = xR, xD = xD,
                       MSDelta=MSDelta, MADelta=MADelta, MRDelta=MRDelta,MDDelta=MDDelta)
    
    shp_regressors = shp %>% left_join(df, by=c("Id"="geo"),keep=FALSE)
    
    
    LM_est = lm(model$regressors, data = df)
    initialOptim = coef(LM_est)[2]
    
    # call julia optimization
    if (variables == "xS"){ W1 = MS}
    if (variables == "xA"){ W1 = MA}
    if (variables == "xR"){ W1 = MR}
    if (variables == "xD"){ W1 = MD}
    
    if (DEBUG == TRUE){ print("Starting Julia optimization for WN SARD 1Mat") }
    outSARD_WN1MatEstimate = call_julia_LogLik_WN_1Mat(X,Y,W1,initialOptim)
    coefSARD_WN1Mat = outSARD_WN1MatEstimate$coef
    se_coefSARD_WN1Mat = outSARD_WN1MatEstimate$se_coef
    pvalue_coefSARD_WN1Mat = outSARD_WN1MatEstimate$pvalue_coef
    residualsSARD_WN1Mat = outSARD_WN1MatEstimate$residuals
    
    SpatError = compute_spatial_error_mat(residualsSARD_WN1Mat,shp,maxLag = 5)
    Werr = SpatError$Werr
    
    # initial condition of optimizer via WN SARD
    # initial spatial error to 0.5 
    allCoef = coefSARD_WN1Mat
    initialOptim = c(allCoef[2],0.5)  ## quiiii
    
    # call julia optimization
    if (DEBUG == TRUE){ print("Starting Julia optimization for SARD") }
    outSARD1_Estimate = call_julia_LogLik_1Mat(X,Y,W1,Werr,initialOptim)
    outSARD1_Estimate$residuals = outSARD1_Estimate$residuals - outSARD1_Estimate$coef[3]*Werr %*% outSARD1_Estimate$residuals
    
    LAR_LM = NULL
    LogLik = LAR_LM$LogLiKelyhood
    AICc = LAR_LM$AICc
    R2N = LAR_LM$R2Nagelkerke
    
    shp_regressors = shp_regressors %>% select(c("y0"="y0.x","yT"="yT.x","delta"="delta.x","km2"="km2.x","xS","xA","xR","xD","MSDelta","MADelta","MRDelta","MDDelta")) %>% mutate("WhAy0"=WhA%*%df$y0,"WhRy0"=WhR%*%df$y0)
    
    outEstimate = listN(outSARD1_Estimate, SpatError, LogLik, AICc,R2N,shp_regressors,MsDeriv)
    return(outEstimate)
    
}


