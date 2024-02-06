rm(list = ls())
source("lib/L_loadAll.R")
DEBUG = FALSE
PARALLEL = TRUE
NPROCS = 12
initJulia()

Na=100
Nm = 1000
tau= 0.05
NeS = 100
SARDp = list(gammaS = 0.0, gammaA = -0.035, gammaR = 0.05, gammaD = 0.105, hA = 0.15, hR = 0.4)

# PDE once
# PDEAll = call_julia_computePDE(tau,SARDp)
# save(PDEAll,file="../datasets_montecarlo/PDE.RData")
# 
# # Agents for error
NSigma = 100

# Na = 50000
# AgentsAll = call_julia_computeAgents(NSigma,Na,tau,SARDp)
# save(AgentsAll,file="../datasets_montecarlo/AgentsAll50k.RData")

# Na = 100000
# AgentsAll = call_julia_computeAgents(NSigma,Na,tau,SARDp)
# save(AgentsAll,file="../datasets_montecarlo/AgentsAll100k.RData")

# Na = 200000
# AgentsAll = call_julia_computeAgents(NSigma,Na,tau,SARDp)
# save(AgentsAll,file="../datasets_montecarlo/AgentsAll200k.RData")

# Np = 225 ----
Np = 225

## create PDE shape ----
shpMC = createShape(Np, typeOfDist = "Uniform")
load(file="../datasets_montecarlo/PDE.RData")
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
save(data,shp,file="../datasets_montecarlo/PDEdatashp225.RData")

## create regressors ----
load(file="../datasets_montecarlo/PDEdatashp225.RData")

MsDeriv = GFDM(data,torus=TRUE)
D = compute_D(data,longlat=FALSE,torus=TRUE)
WhA = compute_WhAR(D,data,SARDp$hA)
WhR = compute_WhAR(D,data,SARDp$hR)
xA = compute_xAR(data,MsDeriv, WhA)
xR = compute_xAR(data,MsDeriv, WhR)
xD = compute_xD(data,MsDeriv)
MA = compute_MARLag(data,MsDeriv,WhA)
MR = compute_MARLag(data,MsDeriv,WhR)
MD = compute_MDLag(MsDeriv)
save(data,shp,xA,xR,xD,MA,MR,MD,file="../datasets_montecarlo/DataShpRegressors225.RData")

## estimate PDE LL ----
load("../datasets_montecarlo/DataShpRegressors225.RData")
estimatePDE225 = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)

# LM_est = estimatePDE225; s = summary(LM_est)
# print(paste("xA in ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/s$coefficients["xA","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/abs(SARDp$gammaA),digits=2), " relative error.", sep=""))
# print(paste("xR in ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/s$coefficients["xR","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/abs(SARDp$gammaR),digits=2), " relative error.", sep=""))
# print(paste("xD in ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/s$coefficients["xD","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/abs(SARDp$gammaD),digits=2), " relative error.", sep=""))

## read and aggregate agents  ----
load("../datasets_montecarlo/DataShpRegressors225.RData")
### 50k ----
load("../datasets_montecarlo/AgentsAll50k.RData")
agents0_MC = AgentsAll$agents0_MC
agentsT_MC = AgentsAll$agentsT_MC
km2 = shp$km2

error50k225 = matrix(data=NA,nrow=NSigma,ncol=225)
pb = txtProgressBar(min = 1, max = NSigma, initial = 1,style=3) 
for (m in 1:NSigma){
    agents0m = agents0_MC[m,,]
    agentsTm = agentsT_MC[m,,]
    
    agents0Sf = sf_point(agents0m)
    st_crs(agents0Sf) ="WGS84"
    aggrAgents0 = lengths(st_intersects(shp,agents0Sf))
    y0m = aggrAgents0/(km2*sum(aggrAgents0))
    
    agentsTSf = sf_point(agentsTm)
    st_crs(agentsTSf) ="WGS84"
    aggrAgentsT = lengths(st_intersects(shp,agentsTSf))
    yTm = aggrAgentsT/(km2*sum(aggrAgentsT))
    
    deltam = yTm-y0m
    errorm = shp$delta-deltam
    error50k225[m,] = errorm
    
    setTxtProgressBar(pb,m)
}
close(pb)

Sigma50k225 = cvCovEst(error50k)
# save(Sigma50k225,file="../datasets_montecarlo/Sigma50k225.RData")
load(file="../datasets_montecarlo/Sigma50k225.RData")

## montecarlo ----
load("../datasets_montecarlo/DataShpRegressors225.RData")

### montecarlo 50k ----
load(file="../datasets_montecarlo/Sigma50k225.RData")
coefMC_225_50k_LM = matrix(data=NA,nrow=Nm,ncol=3)
coefMC_225_50k_IV = matrix(data=NA,nrow=Nm,ncol=3)
coefMC_225_50k_LL = matrix(data=NA,nrow=Nm,ncol=3)
errorMc = mvrnorm(Nm,mu=rep(0,225),Sigma=Sigma50k225$estimate)

pb = txtProgressBar(min = 1, max = Nm, initial = 1,style=3) 
for (m in 1:Nm){
    shpMcm = mutate(shp,delta = delta + errorMc[m,])
    dataMcm = mutate(data,delta = delta + errorMc[m,])
    
    estimatePDE225LMm = estimate_ARD_MC_LM(dataMcm,shpMcm,xA,xR,xD,MA,MR,MD)
    estimatePDE225IVm = estimate_ARD_MC_IV(dataMcm,shpMcm,xA,xR,xD,MA,MR,MD)
    estimatePDE225LLm = estimate_ARD_MC_LL(dataMcm,shpMcm,xA,xR,xD,MA,MR,MD)
    
    coefMC_225_50k_LM[m,] = coef(estimatePDE225LMm)[c(1,3,5)]
    coefMC_225_50k_IV[m,] = coef(estimatePDE225IVm)[c(1,2,3)]
    coefMC_225_50k_LL[m,] = estimatePDE225LLm$outARD_3MatEstimate$coef[c(1,2,3)]
    
    setTxtProgressBar(pb,m)
}
close(pb)
save(coefMC_225_50k_LM,coefMC_225_50k_IV,coefMC_225_50k_LL,file="../datasets_montecarlo/resultMC225.RData")

# Np = 256 ----
Np = 256

## create PDE shape ----
shpMC = createShape(Np, typeOfDist = "Uniform")
load(file="../datasets_montecarlo/PDE.RData")
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
save(data,shp,file="../datasets_montecarlo/PDEdatashp256.RData")

## create regressors ----
load(file="../datasets_montecarlo/PDEdatashp256.RData")

MsDeriv = GFDM(data,torus=TRUE)
D = compute_D(data,longlat=FALSE,torus=TRUE)
WhA = compute_WhAR(D,data,SARDp$hA)
WhR = compute_WhAR(D,data,SARDp$hR)
xA = compute_xAR(data,MsDeriv, WhA)
xR = compute_xAR(data,MsDeriv, WhR)
xD = compute_xD(data,MsDeriv)
MA = compute_MARLag(data,MsDeriv,WhA)
MR = compute_MARLag(data,MsDeriv,WhR)
MD = compute_MDLag(MsDeriv)
save(data,shp,xA,xR,xD,MA,MR,MD,file="../datasets_montecarlo/DataShpRegressors256.RData")

## estimate PDE LL ----
load("../datasets_montecarlo/DataShpRegressors256.RData")
estimatePDE256 = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)

# LM_est = estimatePDE256; s = summary(LM_est)
# print(paste("xA in ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/s$coefficients["xA","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/abs(SARDp$gammaA),digits=2), " relative error.", sep=""))
# print(paste("xR in ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/s$coefficients["xR","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/abs(SARDp$gammaR),digits=2), " relative error.", sep=""))
# print(paste("xD in ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/s$coefficients["xD","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/abs(SARDp$gammaD),digits=2), " relative error.", sep=""))

## read and aggregate agents  ----
load("../datasets_montecarlo/DataShpRegressors256.RData")
### 50k ----
load("../datasets_montecarlo/AgentsAll50k.RData")
agents0_MC = AgentsAll$agents0_MC
agentsT_MC = AgentsAll$agentsT_MC
km2 = shp$km2

error50k256 = matrix(data=NA,nrow=NSigma,ncol=256)
pb = txtProgressBar(min = 1, max = NSigma, initial = 1,style=3) 
for (m in 1:NSigma){
    agents0m = agents0_MC[m,,]
    agentsTm = agentsT_MC[m,,]
    
    agents0Sf = sf_point(agents0m)
    st_crs(agents0Sf) ="WGS84"
    aggrAgents0 = lengths(st_intersects(shp,agents0Sf))
    y0m = aggrAgents0/(km2*sum(aggrAgents0))
    
    agentsTSf = sf_point(agentsTm)
    st_crs(agentsTSf) ="WGS84"
    aggrAgentsT = lengths(st_intersects(shp,agentsTSf))
    yTm = aggrAgentsT/(km2*sum(aggrAgentsT))
    
    deltam = yTm-y0m
    errorm = shp$delta-deltam
    error50k256[m,] = errorm
    
    setTxtProgressBar(pb,m)
}
close(pb)

Sigma50k256 = cvCovEst(error50k)
# save(Sigma50k256,file="../datasets_montecarlo/Sigma50k256.RData")
load(file="../datasets_montecarlo/Sigma50k256.RData")

## montecarlo ----
load("../datasets_montecarlo/DataShpRegressors256.RData")

### montecarlo 50k ----
load(file="../datasets_montecarlo/Sigma50k256.RData")
coefMC_256_50k_LM = matrix(data=NA,nrow=Nm,ncol=3)
coefMC_256_50k_IV = matrix(data=NA,nrow=Nm,ncol=3)
coefMC_256_50k_LL = matrix(data=NA,nrow=Nm,ncol=3)
errorMc = mvrnorm(Nm,mu=rep(0,256),Sigma=Sigma50k256$estimate)

pb = txtProgressBar(min = 1, max = Nm, initial = 1,style=3) 
for (m in 1:Nm){
    shpMcm = mutate(shp,delta = delta + errorMc[m,])
    dataMcm = mutate(data,delta = delta + errorMc[m,])
    
    estimatePDE256LMm = estimate_ARD_MC_LM(dataMcm,shpMcm,xA,xR,xD,MA,MR,MD)
    estimatePDE256IVm = estimate_ARD_MC_IV(dataMcm,shpMcm,xA,xR,xD,MA,MR,MD)
    estimatePDE256LLm = estimate_ARD_MC_LL(dataMcm,shpMcm,xA,xR,xD,MA,MR,MD)
    
    coefMC_256_50k_LM[m,] = coef(estimatePDE256LMm)[c(1,3,5)]
    coefMC_256_50k_IV[m,] = coef(estimatePDE256IVm)[c(1,2,3)]
    coefMC_256_50k_LL[m,] = estimatePDE256LLm$outARD_3MatEstimate$coef[c(1,2,3)]
    
    setTxtProgressBar(pb,m)
}
close(pb)
save(coefMC_256_50k_LM,coefMC_256_50k_IV,coefMC_256_50k_LL,file="../datasets_montecarlo/resultMC256.RData")

# Np = 289 ----
Np = 289

## create PDE shape ----
shpMC = createShape(Np, typeOfDist = "Uniform")
load(file="../datasets_montecarlo/PDE.RData")
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
save(data,shp,file="../datasets_montecarlo/PDEdatashp289.RData")

## create regressors ----
load(file="../datasets_montecarlo/PDEdatashp289.RData")

MsDeriv = GFDM(data,torus=TRUE)
D = compute_D(data,longlat=FALSE,torus=TRUE)
WhA = compute_WhAR(D,data,SARDp$hA)
WhR = compute_WhAR(D,data,SARDp$hR)
xA = compute_xAR(data,MsDeriv, WhA)
xR = compute_xAR(data,MsDeriv, WhR)
xD = compute_xD(data,MsDeriv)
MA = compute_MARLag(data,MsDeriv,WhA)
MR = compute_MARLag(data,MsDeriv,WhR)
MD = compute_MDLag(MsDeriv)
save(data,shp,xA,xR,xD,MA,MR,MD,file="../datasets_montecarlo/DataShpRegressors289.RData")

## estimate PDE LL ----
load("../datasets_montecarlo/DataShpRegressors289.RData")
estimatePDE289 = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)

# LM_est = estimatePDE289; s = summary(LM_est)
# print(paste("xA in ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/s$coefficients["xA","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/abs(SARDp$gammaA),digits=2), " relative error.", sep=""))
# print(paste("xR in ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/s$coefficients["xR","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/abs(SARDp$gammaR),digits=2), " relative error.", sep=""))
# print(paste("xD in ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/s$coefficients["xD","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/abs(SARDp$gammaD),digits=2), " relative error.", sep=""))

## read and aggregate agents  ----
load("../datasets_montecarlo/DataShpRegressors289.RData")
### 50k ----
load("../datasets_montecarlo/AgentsAll50k.RData")
agents0_MC = AgentsAll$agents0_MC
agentsT_MC = AgentsAll$agentsT_MC
km2 = shp$km2

error50k289 = matrix(data=NA,nrow=NSigma,ncol=289)
pb = txtProgressBar(min = 1, max = NSigma, initial = 1,style=3) 
for (m in 1:NSigma){
    agents0m = agents0_MC[m,,]
    agentsTm = agentsT_MC[m,,]
    
    agents0Sf = sf_point(agents0m)
    st_crs(agents0Sf) ="WGS84"
    aggrAgents0 = lengths(st_intersects(shp,agents0Sf))
    y0m = aggrAgents0/(km2*sum(aggrAgents0))
    
    agentsTSf = sf_point(agentsTm)
    st_crs(agentsTSf) ="WGS84"
    aggrAgentsT = lengths(st_intersects(shp,agentsTSf))
    yTm = aggrAgentsT/(km2*sum(aggrAgentsT))
    
    deltam = yTm-y0m
    errorm = shp$delta-deltam
    error50k289[m,] = errorm
    
    setTxtProgressBar(pb,m)
}
close(pb)

Sigma50k289 = cvCovEst(error50k)
# save(Sigma50k289,file="../datasets_montecarlo/Sigma50k289.RData")
load(file="../datasets_montecarlo/Sigma50k289.RData")

## montecarlo ----
load("../datasets_montecarlo/DataShpRegressors289.RData")

### montecarlo 50k ----
load(file="../datasets_montecarlo/Sigma50k289.RData")
coefMC_289_50k_LM = matrix(data=NA,nrow=Nm,ncol=3)
coefMC_289_50k_IV = matrix(data=NA,nrow=Nm,ncol=3)
coefMC_289_50k_LL = matrix(data=NA,nrow=Nm,ncol=3)
errorMc = mvrnorm(Nm,mu=rep(0,289),Sigma=Sigma50k289$estimate)

pb = txtProgressBar(min = 1, max = Nm, initial = 1,style=3) 
for (m in 1:Nm){
    shpMcm = mutate(shp,delta = delta + errorMc[m,])
    dataMcm = mutate(data,delta = delta + errorMc[m,])
    
    estimatePDE289LMm = estimate_ARD_MC_LM(dataMcm,shpMcm,xA,xR,xD,MA,MR,MD)
    estimatePDE289IVm = estimate_ARD_MC_IV(dataMcm,shpMcm,xA,xR,xD,MA,MR,MD)
    estimatePDE289LLm = estimate_ARD_MC_LL(dataMcm,shpMcm,xA,xR,xD,MA,MR,MD)
    
    coefMC_289_50k_LM[m,] = coef(estimatePDE289LMm)[c(1,3,5)]
    coefMC_289_50k_IV[m,] = coef(estimatePDE289IVm)[c(1,2,3)]
    coefMC_289_50k_LL[m,] = estimatePDE289LLm$outARD_3MatEstimate$coef[c(1,2,3)]
    
    setTxtProgressBar(pb,m)
}
close(pb)
save(coefMC_289_50k_LM,coefMC_289_50k_IV,coefMC_289_50k_LL,file="../datasets_montecarlo/resultMC289.RData")





