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

# PDE once ----
PDEAll = call_julia_computePDE(tau,SARDp)
# save(PDEAll,file="../datasets_montecarlo/PDE.RData")
# load("../datasets_montecarlo/PDE.RData")

 
# Agents for error ----
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
estimatePDE225_LM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE225_IV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE225_LL = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
save(estimatePDE225_LM,estimatePDE225_IV,estimatePDE225_LL,file="../datasets_montecarlo/resultPDE225.RData")

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

Sigma50k225 = cvCovEst(error50k225)
save(Sigma50k225,file="../datasets_montecarlo/Sigma50k225.RData")
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
estimatePDE256_LM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE256_IV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE256_LL = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
save(estimatePDE256_LM,estimatePDE256_IV,estimatePDE256_LL,file="../datasets_montecarlo/resultPDE256.RData")

LM_est = estimatePDE256_LM; s = summary(LM_est)
print(paste("xA in ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/s$coefficients["xA","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/abs(SARDp$gammaA),digits=2), " relative error.", sep=""))
print(paste("xR in ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/s$coefficients["xR","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/abs(SARDp$gammaR),digits=2), " relative error.", sep=""))
print(paste("xD in ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/s$coefficients["xD","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/abs(SARDp$gammaD),digits=2), " relative error.", sep=""))

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

Sigma50k256 = cvCovEst(error50k256)
save(Sigma50k256,file="../datasets_montecarlo/Sigma50k256.RData")
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

## plot shp ----
load("../datasets_montecarlo/DataShpRegressors256.RData")

dev.new()
plot(shp["y0"])
dev.copy2pdf(file="../datasets_montecarlo/y0256.pdf")

dev.new()
plot(shp["yT"])
dev.copy2pdf(file="../datasets_montecarlo/yT256.pdf")

dev.new()
plot(shp["delta"])
dev.copy2pdf(file="../datasets_montecarlo/delta256.pdf")

dev.new()
plot(cbind(shp,xA)["xA"])
dev.copy2pdf(file="../datasets_montecarlo/xA256.pdf")

dev.new()
plot(cbind(shp,xR)["xR"])
dev.copy2pdf(file="../datasets_montecarlo/xR256.pdf")

dev.new()
plot(cbind(shp,xD)["xD"])
dev.copy2pdf(file="../datasets_montecarlo/xD256.pdf")



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
estimatePDE289_LM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE289_IV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE289_LL = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
save(estimatePDE289_LM,estimatePDE289_IV,estimatePDE289_LL,file="../datasets_montecarlo/resultPDE289.RData")

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

Sigma50k289 = cvCovEst(error50k289)
save(Sigma50k289,file="../datasets_montecarlo/Sigma50k289.RData")
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





# montecarlo PDE different Np ----
lengthAll = seq(from=10,to=20,by=1)
NpAll = lengthAll^2
coefSARDNpAll = matrix(NA,nrow=length(NpAll),ncol=3)
colnames(coefSARDNpAll) = c("gammaA","gammaR","gammaD")
rownames(coefSARDNpAll) = NpAll
for (i in 1:length(NpAll)){
    Np = NpAll[i]
    print(i)
    print(Np)
    ## create PDE shape ----
    shpMC = createShape(Np, typeOfDist = "Uniform")
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
    
    ## create regressors ----
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
    
    ## estimate PDE LM ----
    estimatePDE_LM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
    
    coefSARDNpAll[i,] = coef(estimatePDE_LM)[c(1,3,5)]
}

gammaA = coefSARDNpAll[,"gammaA"]
gammaR = coefSARDNpAll[,"gammaR"]
gammaD = coefSARDNpAll[,"gammaD"]

Np = NpAll
dev.new()
plot(Np,gammaA,type="b")
abline(h=SARDp$gammaA,col="red")
grid()
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridGammaA.pdf")
dev.new()
plot(Np,gammaR,type="b")
abline(h=SARDp$gammaR,col="red")
grid()
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridGammaR.pdf")
dev.new()
plot(Np,gammaD,type="b")
abline(h=SARDp$gammaD,col="red")
grid()
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridGammaD.pdf")


# montecarlo analysis ----
load("../datasets_montecarlo/resultMC225.RData")
load("../datasets_montecarlo/resultMC256.RData")
load("../datasets_montecarlo/resultMC289.RData")
load("../datasets_montecarlo/resultPDE225.RData")
load("../datasets_montecarlo/resultPDE256.RData")
load("../datasets_montecarlo/resultPDE289.RData")


## 225 Result ----
table225 = matrix(data = NA, nrow = 3, ncol = 4*3)
rownames(table225) = c("gammaA","gammaR","gammaD")
colnames(table225) = c("PDE_LM","meanMC_LM","medianMC_LM","sdMC_LM",
                       "PDE_IV","meanMC_IV","medianMC_IV","sdMC_IV",
                       "PDE_LL","meanMC_LL","medianMC_LL","sdMC_LL")

table225[,"PDE_LM"] = coef(estimatePDE225_LM)[c(1,3,5)]
table225[,"meanMC_LM"] = apply(coefMC_225_50k_LM,FUN=mean,MARGIN=2)
table225[,"medianMC_LM"] = apply(coefMC_225_50k_LM,FUN=median,MARGIN=2)
table225[,"sdMC_LM"] = apply(coefMC_225_50k_LM,FUN=sd,MARGIN=2)

table225[,"PDE_IV"] = coef(estimatePDE225_IV)[c(1,2,3)]
table225[,"meanMC_IV"] = apply(coefMC_225_50k_IV,FUN=mean,MARGIN=2)
table225[,"medianMC_IV"] = apply(coefMC_225_50k_IV,FUN=median,MARGIN=2)
table225[,"sdMC_IV"] = apply(coefMC_225_50k_IV,FUN=sd,MARGIN=2)

table225[,"PDE_LL"] = estimatePDE225_LL$outARD_3MatEstimate$coef[c(1,2,3)]
table225[,"meanMC_LL"] = apply(coefMC_225_50k_LL,FUN=mean,MARGIN=2)
table225[,"medianMC_LL"] = apply(coefMC_225_50k_LL,FUN=median,MARGIN=2)
table225[,"sdMC_LL"] = apply(coefMC_225_50k_LL,FUN=sd,MARGIN=2)

# colnames(table225) = sapply(colnames(table225),FUN=function(x) return(paste(x,"225",sep="")))

## 256 Result ----
table256 = matrix(data = NA, nrow = 3, ncol = 4*3)
rownames(table256) = c("gammaA","gammaR","gammaD")
colnames(table256) = c("PDE_LM","meanMC_LM","medianMC_LM","sdMC_LM",
                       "PDE_IV","meanMC_IV","medianMC_IV","sdMC_IV",
                       "PDE_LL","meanMC_LL","medianMC_LL","sdMC_LL")

table256[,"PDE_LM"] = coef(estimatePDE256_LM)[c(1,3,5)]
table256[,"meanMC_LM"] = apply(coefMC_256_50k_LM,FUN=mean,MARGIN=2)
table256[,"medianMC_LM"] = apply(coefMC_256_50k_LM,FUN=median,MARGIN=2)
table256[,"sdMC_LM"] = apply(coefMC_256_50k_LM,FUN=sd,MARGIN=2)

table256[,"PDE_IV"] = coef(estimatePDE256_IV)[c(1,2,3)]
table256[,"meanMC_IV"] = apply(coefMC_256_50k_IV,FUN=mean,MARGIN=2)
table256[,"medianMC_IV"] = apply(coefMC_256_50k_IV,FUN=median,MARGIN=2)
table256[,"sdMC_IV"] = apply(coefMC_256_50k_IV,FUN=sd,MARGIN=2)

table256[,"PDE_LL"] = estimatePDE256_LL$outARD_3MatEstimate$coef[c(1,2,3)]
table256[,"meanMC_LL"] = apply(coefMC_256_50k_LL,FUN=mean,MARGIN=2)
table256[,"medianMC_LL"] = apply(coefMC_256_50k_LL,FUN=median,MARGIN=2)
table256[,"sdMC_LL"] = apply(coefMC_256_50k_LL,FUN=sd,MARGIN=2)

# colnames(table256) = sapply(colnames(table256),FUN=function(x) return(paste(x,"256",sep="")))

## 289 Result ----
table289 = matrix(data = NA, nrow = 3, ncol = 4*3)
rownames(table289) = c("gammaA","gammaR","gammaD")
colnames(table289) = c("PDE_LM","meanMC_LM","medianMC_LM","sdMC_LM",
                       "PDE_IV","meanMC_IV","medianMC_IV","sdMC_IV",
                       "PDE_LL","meanMC_LL","medianMC_LL","sdMC_LL")

table289[,"PDE_LM"] = coef(estimatePDE289_LM)[c(1,3,5)]
table289[,"meanMC_LM"] = apply(coefMC_289_50k_LM,FUN=mean,MARGIN=2)
table289[,"medianMC_LM"] = apply(coefMC_289_50k_LM,FUN=median,MARGIN=2)
table289[,"sdMC_LM"] = apply(coefMC_289_50k_LM,FUN=sd,MARGIN=2)
 
table289[,"PDE_IV"] = coef(estimatePDE289_IV)[c(1,2,3)]
table289[,"meanMC_IV"] = apply(coefMC_289_50k_IV,FUN=mean,MARGIN=2)
table289[,"medianMC_IV"] = apply(coefMC_289_50k_IV,FUN=median,MARGIN=2)
table289[,"sdMC_IV"] = apply(coefMC_289_50k_IV,FUN=sd,MARGIN=2)

table289[,"PDE_LL"] = estimatePDE289_LL$outARD_3MatEstimate$coef[c(1,2,3)]
table289[,"meanMC_LL"] = apply(coefMC_289_50k_LL,FUN=mean,MARGIN=2)
table289[,"medianMC_LL"] = apply(coefMC_289_50k_LL,FUN=median,MARGIN=2)
table289[,"sdMC_LL"] = apply(coefMC_289_50k_LL,FUN=sd,MARGIN=2)

# colnames(table289) = sapply(colnames(table289),FUN=function(x) return(paste(x,"289",sep="")))

## table All ----
tableAll = cbind(table225,table256,table289)

## output in xtable
table256 = rbind(table256[,1:4],table256[,5:8],table256[,9:12])

xtable(table225,auto=TRUE)
xtable(table256,auto=TRUE)
xtable(table289,auto=TRUE)
xtable(tableAll,auto=TRUE)
xtable(tableAll,display=rep("f",ncol(tableAll)+1))


## table All sorted by estimation method ----
tableLM225 = t(table225[,1:4])
tableLM256 = t(table256[,1:4])
tableLM289 = t(table289[,1:4])

tableLM_gammaA = cbind(tableLM225[,colnames(tableLM225)=="gammaA"],
                       tableLM256[,colnames(tableLM256)=="gammaA"],
                       tableLM289[,colnames(tableLM289)=="gammaA"])
colnames(tableLM_gammaA) = c("225", "256", "289")

tableLM_gammaR = cbind(tableLM225[,colnames(tableLM225)=="gammaR"],
                       tableLM256[,colnames(tableLM256)=="gammaR"],
                       tableLM289[,colnames(tableLM289)=="gammaR"])
colnames(tableLM_gammaR) = c("225", "256", "289")

tableLM_gammaD = cbind(tableLM225[,colnames(tableLM225)=="gammaD"],
                       tableLM256[,colnames(tableLM256)=="gammaD"],
                       tableLM289[,colnames(tableLM289)=="gammaD"])
colnames(tableLM_gammaD) = c("225", "256", "289")

table_LM = cbind(tableLM_gammaA, tableLM_gammaR, tableLM_gammaD)
xtable(table_LM, digits=5)


tableIV225 = t(table225[,5:8])
tableIV256 = t(table256[,5:8])
tableIV289 = t(table289[,5:8])

tableIV_gammaA = cbind(tableIV225[,colnames(tableIV225)=="gammaA"],
                       tableIV256[,colnames(tableIV256)=="gammaA"],
                       tableIV289[,colnames(tableIV289)=="gammaA"])
colnames(tableIV_gammaA) = c("225", "256", "289")

tableIV_gammaR = cbind(tableIV225[,colnames(tableIV225)=="gammaR"],
                       tableIV256[,colnames(tableIV256)=="gammaR"],
                       tableIV289[,colnames(tableIV289)=="gammaR"])
colnames(tableIV_gammaR) = c("225", "256", "289")

tableIV_gammaD = cbind(tableIV225[,colnames(tableIV225)=="gammaD"],
                       tableIV256[,colnames(tableIV256)=="gammaD"],
                       tableIV289[,colnames(tableIV289)=="gammaD"])
colnames(tableIV_gammaD) = c("225", "256", "289")

table_IV = cbind(tableIV_gammaA, tableIV_gammaR, tableIV_gammaD)
xtable(table_IV, digits=5)

tableLL225 = t(table225[,9:12])
tableLL256 = t(table256[,9:12])
tableLL289 = t(table289[,9:12])

tableLL_gammaA = cbind(tableLL225[,colnames(tableLL225)=="gammaA"],
                       tableLL256[,colnames(tableLL256)=="gammaA"],
                       tableLL289[,colnames(tableLL289)=="gammaA"])
colnames(tableLL_gammaA) = c("225", "256", "289")

tableLL_gammaR = cbind(tableLL225[,colnames(tableLL225)=="gammaR"],
                       tableLL256[,colnames(tableLL256)=="gammaR"],
                       tableLL289[,colnames(tableLL289)=="gammaR"])
colnames(tableLL_gammaR) = c("225", "256", "289")

tableLL_gammaD = cbind(tableLL225[,colnames(tableLL225)=="gammaD"],
                       tableLL256[,colnames(tableLL256)=="gammaD"],
                       tableLL289[,colnames(tableLL289)=="gammaD"])
colnames(tableLL_gammaD) = c("225", "256", "289")

table_LL = cbind(tableLL_gammaA, tableLL_gammaR, tableLL_gammaD)
xtable(table_LL, digits=5)

