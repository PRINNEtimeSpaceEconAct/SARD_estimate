MonteCarlo <- function(Np, Na, Nm, tau, typeOfDist, SARDp, hA_range, 
                       hR_range, NeS = 1000){
    # add documentation of the function
    # INPUT
    # OUTPUT
    
    #Create shape file
    shpMC = createShapeVoronoi(Np, plot=FALSE)
    
    #Create agents
    AgentsAll = call_julia_computeAgents(Nm,Na,tau,SARDp)
    agents0_MC = AgentsAll$agents0_MC
    agentsT_MC = AgentsAll$agentsT_MC
    
    #For s
    SComputed = call_julia_computeS(NeS)
    X = SComputed$X
    Y = SComputed$Y
    S = SComputed$S
    
    outSARDEstimateMonteCarlo = list()
    for (m in 1:Nm){
        agents0 = agents0_MC[m,,]
        agentsT = agentsT_MC[m,,]
        
        data_shp = createDataframe(agents0, agentsT, shpMC, tau, X, Y, S)
        data = data_shp$data
        
        select_h=chose_hAhR(data,hA_range,hR_range,longlat=FALSE)
        outSARDEstimate = estimate_SARD_auto(data,shpMC,select_h$hABest,
                                             select_h$hRBest,longlat=FALSE)
        
        outSARDEstimate$hABest = select_h$hABest
        outSARDEstimate$hRBest = select_h$hRBest
        
        outSARDEstimateMonteCarlo[[m]] = outSARDEstimate
    }
    
    return(outSARDEstimateMonteCarlo)
    
}

MonteCarloFixedhAhR_IV <- function(Np, Na, Nm, tau, typeOfDist, SARDp,hA,hR,
                                   NeS = 100){
    # add documentation of the function
    # INPUT
    # OUTPUT
    
    #Create shape file
    shpMC = createShapeVoronoi(Np, plot=FALSE)
    
    #Create agents
    AgentsAll = call_julia_computeAgents(Nm,Na,tau,SARDp)
    agents0_MC = AgentsAll$agents0_MC
    agentsT_MC = AgentsAll$agentsT_MC
    
    #For s
    SComputed = call_julia_computeS(NeS)
    X = SComputed$X
    Y = SComputed$Y
    S = SComputed$S
    
    outIVEstimateMonteCarlo = list()
    for (m in 1:Nm){
        agents0 = agents0_MC[m,,]
        agentsT = agentsT_MC[m,,]
        
        data_shp = createDataframe(agents0, agentsT, shpMC, tau, X, Y, S)
        data = data_shp$data
        
        outIVEstimate = estimate_IV_SARD_auto(data,hA,hR)
        
        outIVEstimateMonteCarlo[[m]] = outIVEstimate
    }
    
    return(outIVEstimateMonteCarlo)
    
}



