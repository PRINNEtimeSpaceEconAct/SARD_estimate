rm(list = ls())
source("lib/L_loadAll.R")
DEBUG = FALSE
PARALLEL = TRUE
NPROCS = 12
initJulia()

Na=100
Nm = 1000
tau= 1.0
NeS = 100
SARDp = list(gammaS = 0.0, gammaA = -0.00175, gammaR = 0.0025, gammaD = 0.00525, hA = 0.15, hR = 0.4)
load("../datasets_montecarlo/PDE.RData")


lengthAll = seq(from=10,to=50,by=1)
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
save(coefSARDNpAll,file="../datasets_montecarlo/montecarloPDENpGrid.RData")
     
load(file="../datasets_montecarlo/montecarloPDENpGrid.RData")
gammaA = coefSARDNpAll[,"gammaA"]
gammaR = coefSARDNpAll[,"gammaR"]
gammaD = coefSARDNpAll[,"gammaD"]

errRelgammaA = (SARDp$gammaA-gammaA)/SARDp$gammaA
errRelgammaR = (SARDp$gammaR-gammaR)/SARDp$gammaR
errRelgammaD = (SARDp$gammaD-gammaD)/SARDp$gammaD
Np = NpAll

dev.new()
plot(Np,errRelgammaA,type="b",ylab = "")
lines(Np,errRelgammaR,type="b",col="blue")
lines(Np,errRelgammaD,type="b",col="red")
legend("topright",legend=c("relative error gammaA","relative error gammaR","relative error gammaD"),col=c("black","blue","red"),lty=1,cex=1.0)
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel.pdf")


# dev.new()
# plot(Np,gammaA,type="b")
# abline(h=SARDp$gammaA,col="red")
# grid()
# dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridGammaA.pdf")
# dev.new()
# plot(Np,gammaR,type="b")
# abline(h=SARDp$gammaR,col="red")
# grid()
# dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridGammaR.pdf")
# dev.new()
# plot(Np,gammaD,type="b")
# abline(h=SARDp$gammaD,col="red")
# grid()
# dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridGammaD.pdf")
# 
