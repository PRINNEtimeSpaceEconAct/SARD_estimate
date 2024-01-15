# MonteCarlo <- function(Np, Na, Nm, tau, typeOfDist = "uniform", SARDp, hA_range, 
#                        hR_range, NeS = 1000){
#     # add documentation of the function
#     # INPUT
#     # OUTPUT
#     
#     #Create shape file
#     shpMC = createShapeVoronoi(Np, typeOfDist = typeOfDist, plot=FALSE)
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
#     shpMC = createShapeVoronoi(Np, typeOfDist = typeOfDist, plot=FALSE)
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

MonteCarloOneRun_LM_Agents_fixedhAhR <- function(Np, Na, tau, typeOfDist = "uniform", SARDp,hA,hR,
                                   NeS = 100,model=list(),variables=c(""),torus=TRUE){
    # add documentation of the function
    # INPUT
    # OUTPUT
    
    #Create shape file
    shpMC = createShapeVoronoi(Np, typeOfDist = typeOfDist, plot=FALSE)
    
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
        plot(outLMEstimate$shp_regressors)
        
    return(listN(outLMEstimate,data,shp))
    
}

MonteCarloOneRun_LM_PDE_fixedhAhR <- function(Np, tau, typeOfDist = "uniform", SARDp,hA,hR,
                                                 NeS = 100,model=list(),variables=c(""),torus=TRUE){
    # add documentation of the function
    # INPUT
    # OUTPUT
    
    #Create shape file
    shpMC = createShapeVoronoi(Np, typeOfDist = typeOfDist, plot=FALSE)
    
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
    
    outLMEstimate = estimate_LM_SARD_autoMC(data,hA,hR,shp,longlat = FALSE,model=model,variables=variables,torus=torus)
    plot(outLMEstimate$shp_regressors)
    
    return(listN(outLMEstimate,data,shp))
    
}


estimate_LM_SARD_autoMC <- function(df,hA,hR,shp,longlat=TRUE,model=list(),variables,torus=TRUE){
    # Estimate SARD WN via IV with distances hA, hR
    # model = list(regressors=as.formula(delta ~ y0 + xS + xA + xR + xD + MSDelta + MADelta + MRDelta + MDDelta),
                 # instruments=as.formula(~ y0 + xS + xA + xR + xD + MS2X + MA2X + MR2X + MD2X))
    # variables = c("ones","y0","xS","xA","xR","xD")
    
    print("Estimating LM for given hA hR") 
    
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
    # IV_est = ivreg(model$regressors, model$instruments, data=df)
    
    LM_est = lm(model$regressors, data = df)
    
    LAR_LM = NULL
    LogLik = LAR_LM$LogLiKelyhood
    AICc = LAR_LM$AICc
    R2N = LAR_LM$R2Nagelkerke
    
    shp_regressors = shp_regressors %>% select(c("y0"="y0.x","yT"="yT.x","delta"="delta.x","xS","xA","xR","xD")) %>% mutate("WhAy0"=WhA%*%df$y0,"WhRy0"=WhR%*%df$y0)
    
    return(listN(LM_est, LogLik, AICc,R2N,shp_regressors))
    
}




