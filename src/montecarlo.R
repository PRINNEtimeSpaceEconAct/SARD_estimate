rm(list = ls())
source("lib/L_loadAll.R")
DEBUG = FALSE
PARALLEL = TRUE
NPROCS = 12
initJulia()

Na=100000
Nm = 100
tau= 0.05
NeS = 100
SARDp = list(gammaS = 0.0, gammaA = -0.035, gammaR = 0.05, gammaD = 0.105, hA = 0.15, hR = 0.4)

# PDE once
PDEAll = call_julia_computePDE(tau,SARDp)
save(PDEAll,file="PDE.RData")

# Agents for error
NSigma = 100

# Na = 50000
AgentsAll = call_julia_computeAgents(NSigma,Na,tau,SARDp)
save(AgentsAll,file="AgentsAll50k.RData")

# Na = 100000
AgentsAll = call_julia_computeAgents(NSigma,Na,tau,SARDp)
save(AgentsAll,file="AgentsAll100k.RData")

# Na = 200000
AgentsAll = call_julia_computeAgents(NSigma,Na,tau,SARDp)
save(AgentsAll,file="AgentsAll200k.RData")

# Np = 256 ----
Np = 256

## create PDE shape ----
shpMC = createShape(Np, typeOfDist = "Uniform")
# load(file="PDE.RData")
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
save(data,shp,file="PDEdatashp256.RData")

## create regressors ----
# load(file="PDEdatashp256.RData")

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
save(data,shp,xA,xR,xD,MA,MR,MD,file="DataShpRegressors256.RData")

## estimate PDE LL ----
# load("DataShpRegressors256.RData")
estimatePDE256 = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)

print(SARDp)
print(estimatePDE256$outARD_3MatEstimate$coef[1:3])
print(paste("xA in ",round(abs(SARDp$gammaA-estimatePDE256$outARD_3MatEstimate$coef[1])/estimatePDE256$outARD_3MatEstimate$se_coef[1],digits=2)," std. error, and ",round(abs(SARDp$gammaA-estimatePDE256$outARD_3MatEstimate$coef[1])/abs(SARDp$gammaA),digits=2), " relative error.", sep=""))
print(paste("xR in ",round(abs(SARDp$gammaR-estimatePDE256$outARD_3MatEstimate$coef[2])/estimatePDE256$outARD_3MatEstimate$se_coef[2],digits=2)," std. error, and ",round(abs(SARDp$gammaA-estimatePDE256$outARD_3MatEstimate$coef[2])/abs(SARDp$gammaR),digits=2), " relative error.", sep=""))
print(paste("xD in ",round(abs(SARDp$gammaD-estimatePDE256$outARD_3MatEstimate$coef[3])/estimatePDE256$outARD_3MatEstimate$se_coef[3],digits=2)," std. error, and ",round(abs(SARDp$gammaA-estimatePDE256$outARD_3MatEstimate$coef[3])/abs(SARDp$gammaD),digits=2), " relative error.", sep=""))
library(sm)
sm.density(estimatePDE256$outARD_3MatEstimate$residuals,model=TRUE)

