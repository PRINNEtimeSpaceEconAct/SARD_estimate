# start ----
rm(list = ls())
source("lib/L_loadAll.R")
library(viridis)
DEBUG = FALSE
PARALLEL = TRUE
NPROCS = 12
initJulia()

Na=100
Nm = 1000
tau= 1.0
NeS = 100
SARDp = list(alpha = 0.01, phi = 0.01, gammaS = 0.0, gammaA = -0.00175, gammaR = 0.0025, gammaD = 0.00525, hA = 0.15, hR = 0.4)

# PDE once ----
PDEAll = call_julia_computePDE(tau,SARDp)
save(PDEAll,file="../datasets_montecarlo/PDE.RData")

# Np = 144 ----
Np = 144

## create PDE shape ----
shpMC = createShape(Np, typeOfDist = "Uniform")
# load(file="../datasets_montecarlo/PDE.RData")
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
# save(data,shp,file="../datasets_montecarlo/PDEdatashp144.RData")


## create regressors ----
# load(file="../datasets_montecarlo/PDEdatashp144.RData")

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
shp = cbind(shp,xA,xR,xD)
plot(shp[c("y0","yT","delta","xA","xR","xD")])
# save(data,shp,xA,xR,xD,MA,MR,MD,file="../datasets_montecarlo/DataShpRegressors144.RData")


## estimate PDE LL ----
# load("../datasets_montecarlo/DataShpRegressors1156.RData")
estimatePDE_LM_NAIVE = estimate_ARD_MC_LM_NAIVE(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE_LM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE_IV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE_LL = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
# save(estimatePDE_LM_NAIVE,estimatePDE_LM,estimatePDE_IV,estimatePDE_LL,file="../datasets_montecarlo/resultPDE.RData")

yT = sum(data$yT)/Np
y0 = sum(data$y0)/Np
coefsLMTilde = estimatePDE_LM$coefficients[c(1,2,3,5,7)]
rhoPhiLM = 2/tau*(1-log( (yT + coefsLMTilde["(Intercept)"]/coefsLMTilde["y0"])/(y0 + coefsLMTilde["(Intercept)"]/coefsLMTilde["y0"]) )/(tau*coefsLMTilde["y0"]) )
coefsLM = coefsLMTilde*(1-tau*rhoPhiLM/2)


# load(file="../datasets_montecarlo/resultPDE.RData")
LM_est = estimatePDE_LM; s = summary(LM_est)
print(paste("xA in ",round((SARDp$gammaA-LM_est$coefficients["xA"])/s$coefficients["xA","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/abs(SARDp$gammaA),digits=2), " relative error.", sep=""))
print(paste("xR in ",round((SARDp$gammaR-LM_est$coefficients["xR"])/s$coefficients["xR","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/abs(SARDp$gammaR),digits=2), " relative error.", sep=""))
print(paste("xD in ",round((SARDp$gammaD-LM_est$coefficients["xD"])/s$coefficients["xD","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/abs(SARDp$gammaD),digits=2), " relative error.", sep=""))

print(paste("xA in ",round((SARDp$gammaA- estimatePDE_LL$outARD_3MatEstimate$coef[1])/abs(SARDp$gammaA),digits=2), " relative error.", sep=""))
print(paste("xR in ",round((SARDp$gammaR- estimatePDE_LL$outARD_3MatEstimate$coef[2])/abs(SARDp$gammaR),digits=2), " relative error.", sep=""))
print(paste("xD in ",round((SARDp$gammaD- estimatePDE_LL$outARD_3MatEstimate$coef[3])/abs(SARDp$gammaD),digits=2), " relative error.", sep=""))



## plot shp ----

dev.new()
plot(shp["y0"],main="")
dev.copy2pdf(file="../datasets_montecarlo/y0144.pdf")

dev.new()
plot(shp["yT"],main="")
dev.copy2pdf(file="../datasets_montecarlo/yT144.pdf")

dev.new()
plot(shp["delta"],main="")
dev.copy2pdf(file="../datasets_montecarlo/delta144.pdf")

dev.new()
plot(cbind(shp,xA)["xA"],main="")
dev.copy2pdf(file="../datasets_montecarlo/xA144.pdf")

dev.new()
plot(cbind(shp,xR)["xR"],main="")
dev.copy2pdf(file="../datasets_montecarlo/xR144.pdf")

dev.new()
plot(cbind(shp,xD)["xD"],main="")
dev.copy2pdf(file="../datasets_montecarlo/xD144.pdf")

## plot contour PDE ----
dev.new()
filled.contour(t(PDE0), color.palette=plasma)
dev.copy2pdf(file="../datasets_montecarlo/y0_continuosSpace.pdf")
dev.new()
filled.contour(t(PDET), color.palette=plasma)
dev.copy2pdf(file="../datasets_montecarlo/yT_continuosSpace.pdf")
dev.new()
filled.contour((t(PDET)-t(PDE0))/tau, color.palette=plasma)
dev.copy2pdf(file="../datasets_montecarlo/delta_continuosSpace.pdf")

# montecarlo PDE different Np ----
load(file="../datasets_montecarlo/PDE.RData")

lengthAll = c(12,20,30,40,50)
NpAll = lengthAll^2
Nm = 1000
coefMNpLMNaive = array(data=NA,dim=c(length(NpAll),Nm,5))
coefMNpLM = array(data=NA,dim=c(length(NpAll),Nm,5))
coefMNpIV = array(data=NA,dim=c(length(NpAll),Nm,5))
coefNpLL = matrix(data=NA,nrow=length(NpAll),ncol=5)
estimateNpLL = list()
AICR2LMNaive = matrix(data = NA, nrow = length(NpAll), ncol = 2); colnames(AICR2LMNaive) = c("AICc","R2")
AICR2LM = matrix(data = NA, nrow = length(NpAll), ncol = 2); colnames(AICR2LM) = c("AICc","R2")
AICR2IV = matrix(data = NA, nrow = length(NpAll), ncol = 2); colnames(AICR2IV) = c("AICc","R2")
AICR2LL = matrix(data = NA, nrow = length(NpAll), ncol = 2); colnames(AICR2LL) = c("AICc","R2")
y0yT = array(data=NA,dim=c(length(NpAll),Nm,2)); 

for (i in 1:length(NpAll)){
    Np = NpAll[i]
    print(i)
    print(Np)
    shpMC = createShape(Np, typeOfDist = "Uniform")
    Xpde = PDEAll$X
    Ypde = PDEAll$Y
    PDE0 = PDEAll$PDE0
    PDET = PDEAll$PDET
    
    #For s
    NeS = 100
    Xs = matrix(0,nrow=NeS,ncol=NeS)
    Ys = matrix(0,nrow=NeS,ncol=NeS)
    S = matrix(0,nrow=NeS,ncol=NeS)
    
    data_shp = createDataframePDE(PDE0, PDET, Xpde, Ypde, shpMC, tau, Xs, Ys, S)
    data = data_shp$data
    shp = data_shp$shp_sf
    
    MsDeriv = GFDM(data,torus=TRUE)
    D = compute_D(data,longlat=FALSE,torus=TRUE)
    WhA = compute_WhAR(D,data,SARDp$hA)
    WhR = compute_WhAR(D,data,SARDp$hR)
    xS = runif(Np)
    xA = compute_xAR(data,MsDeriv, WhA)
    xR = compute_xAR(data,MsDeriv, WhR)
    xD = compute_xD(data,MsDeriv)
    MS = matrix(data=0, nrow=Np,ncol=Np)
    MA = compute_MARLag(data,MsDeriv,WhA)
    MR = compute_MARLag(data,MsDeriv,WhR)
    MD = compute_MDLag(MsDeriv)
    
    MADelta = as.numeric(MA %*% matrix(data$delta))
    MRDelta = as.numeric(MR %*% matrix(data$delta))
    MDDelta = as.numeric(MD %*% matrix(data$delta))
    
    X = data.frame(xA=xA,xR=xR,xD=xD)
    X = as.matrix(X)
    
    # instruments for IV
    MA2X=as.matrix(MA %*% MA %*% X)
    MR2X=as.matrix(MR %*% MR %*% X) 
    MD2X=as.matrix(MD %*% MD %*% X)
    
    
    estimateNpLM_NAIVE = estimate_ARD_MC_LM_NAIVE(data,shp,xA,xR,xD,MA,MR,MD)
    AICR2LMNaive[i,] = c(estimateNpLM_NAIVE$AICc,estimateNpLM_NAIVE$R2)
    
    estimateNpLM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
    coefLM = coef(estimateNpLM$LM_est)
    coefLM = c(coefLM[c(1,2)],0,coefLM[c(3,5,7)],0,coefLM[c(4,6,8)])
    WN_SARD_OLS_BYSARD = LogLikAICcR2(data, c(coefLM,0), 10, xS, xA, xR, xD, 
                                      MS, MA, MR, MD, diag(nrow(data)))
    AICR2LM[i,] = c(WN_SARD_OLS_BYSARD$AICc[[1]],WN_SARD_OLS_BYSARD$R2Nagelkerke)
    
    
    estimateNpIV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)
    AICR2IV[i,] = c(estimateNpIV$AICc,estimateNpIV$R2)
    
    # estimateNpLL[[i]] = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
    # AICR2LL[i,] = c(estimateNpLL[[i]]$AICc,estimateNpLL[[i]]$R2)
    # coefNpLL[i,] = estimateNpLL[[i]]$outARD_3MatEstimate$coef[c(1,2,3,4,5)]
    # 
    # for (m in 1:Nm){
    #     # print(m)
    #     
    #     iBoot = sample(1:nrow(data),replace=T)
    #     datam = data[iBoot,]
    #     xAm = xA[iBoot]
    #     xRm = xR[iBoot]
    #     xDm = xD[iBoot]
    #     MADeltam = MADelta[iBoot]
    #     MRDeltam = MRDelta[iBoot]
    #     MDDeltam = MDDelta[iBoot]
    #     MA2Xm=MA2X[iBoot,]
    #     MR2Xm=MR2X[iBoot,]
    #     MD2Xm=MD2X[iBoot,]
    #     y0yT[i,m,1] = sum(datam$y0)*1/Np
    #     y0yT[i,m,2] = sum(datam$yT)*1/Np
    #     
    # 
    #     estimatePDE_LMNaivem = estimate_ARD_MC_LM_NAIVEBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
    #     coefMNpLMNaive[i,m,] = estimatePDE_LMNaivem$coefficients[c(1,2,3,4,5)]
    #     
    #     estimatePDE_LMm = estimate_ARD_MC_LMBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
    #     coefMNpLM[i,m,] = estimatePDE_LMm$coefficients[c(1,2,3,5,7)]
    #     
    #     estimatePDE_IVm = estimate_ARD_MC_IVBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam,MA2Xm,MR2Xm,MD2Xm)
    #     coefMNpIV[i,m,] = coef(estimatePDE_IVm)[c(1,2,3,4,5)]
    # }
}
save(NpAll,coefMNpLM,coefMNpLMNaive,coefMNpIV,coefNpLL,estimateNpLL,AICR2LMNaive,AICR2LM,AICR2IV,AICR2LL,y0yT,file="../datasets_montecarlo/montecarloPDENpGrid.RData")

## AIC plot ----
load(file="../datasets_montecarlo/montecarloPDENpGrid.RData")
Np = NpAll
DeltaX = 1/sqrt(Np)

dev.new()
plot(DeltaX,AICR2LM[,1],type="p",ylab="",xlab="",xaxt="n",pch=1,lwd=1.5,ylim=range(AICR2LMNaive[,1],AICR2LM[,1],AICR2IV[,1],AICR2LL[,1]))
lines(DeltaX,AICR2LMNaive[,1],type="p",col="red",lwd = 1.5, xaxt="n",pch=2)
lines(DeltaX,AICR2IV[,1],type="p",col="blue",lwd = 1.5,xaxt="n",pch=3)
lines(DeltaX,AICR2LL[,1],type="p",col="purple",lwd = 1.5,xaxt="n",pch=4)
mapply(axis, side = 1, at = DeltaX, labels = paste(c("(COARSEST)\n","","","","(A) "), round(DeltaX,digits=3)," \n (N = ",Np ,")",sep=""),las = 3, hadj=0.8, cex.axis=1)
grid()
title(ylab = "AICc",mgp=c(2,0,0))
title(xlab = TeX(r'(                                                                                                                                                                                                 $\frac{1}{\sqrt{N}}$)'),mgp=c(0.5,0,0))
legend("bottomright",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
sapply(DeltaX,FUN = function(x) return(abline(v = x,lty=3)) );
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEAICc.pdf")

## robust SE  ----
library(latex2exp)
load(file="../datasets_montecarlo/montecarloPDENpGrid.RData")
Np = NpAll
Nm = dim(coefMNpLM)[2]
coefRep = matrix(rep(c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),Nm),nrow=Nm,byrow = T)
coefNpLM = matrix(data=NA,nrow=length(Np),ncol=5)
lowerNpLM = matrix(data=NA,nrow=length(Np),ncol=5)
upperNpLM = matrix(data=NA,nrow=length(Np),ncol=5)
seNpLM = matrix(data=NA,nrow=length(Np),ncol=5)
coefNpLMNaive = matrix(data=NA,nrow=length(Np),ncol=5)
lowerNpLMNaive = matrix(data=NA,nrow=length(Np),ncol=5)
upperNpLMNaive = matrix(data=NA,nrow=length(Np),ncol=5)
seNpLMNaive = matrix(data=NA,nrow=length(Np),ncol=5)
coefNpIV = matrix(data=NA,nrow=length(Np),ncol=5)
lowerNpIV = matrix(data=NA,nrow=length(Np),ncol=5)
upperNpIV = matrix(data=NA,nrow=length(Np),ncol=5)
seNpIV = matrix(data=NA,nrow=length(Np),ncol=5)
lowerNpLL = matrix(data=NA,nrow=length(Np),ncol=5)
upperNpLL = matrix(data=NA,nrow=length(Np),ncol=5)
seNpLL = matrix(data=NA,nrow=length(Np),ncol=5)

# new for nonTilde coefficient LL
coefMNpLLBoot = array(data=NA,dim=c(length(NpAll),Nm,5))

for (i in 1:length(Np)){
    coefMLM = coefMNpLM[i,,]
    coefMLMNaive = coefMNpLMNaive[i,,]    
    coefMIV = coefMNpIV[i,,]   
    
    coefLL = coefNpLL[i,]
    covBeta = estimateNpLL[[i]]$outARD_3MatEstimate$covBeta
    
    for (m in 1:Nm) {
        y0 = y0yT[i,m,1]
        yT = y0yT[i,m,2]
        
        # coefficient notTilde by bootstrap for LM LMNaive and IV
        coefLMTilde = coefMLM[m,]
        alphaTildeLM = coefLMTilde[1]
        phiTildeLM = coefLMTilde[2]
        rhoPhiLM = 2/tau*(1-log( (yT + alphaTildeLM/phiTildeLM)/(y0 + alphaTildeLM/phiTildeLM) )/(tau*phiTildeLM) )
        coefMLM[m,] = coefLMTilde*(1-tau*rhoPhiLM/2)
        
        coefLMNaiveTilde = coefMLMNaive[m,]
        alphaTildeLMNaive = coefLMNaiveTilde[1]
        phiTildeLMNaive = coefLMNaiveTilde[2]
        rhoPhiLMNaive = 2/tau*(1-log( (yT + alphaTildeLMNaive/phiTildeLMNaive)/(y0 + alphaTildeLMNaive/phiTildeLMNaive) )/(tau*phiTildeLMNaive) )
        coefMLMNaive[m,] = coefLMNaiveTilde*(1-tau*rhoPhiLMNaive/2)
        
        coefIVTilde = coefMIV[m,]
        alphaTildeIV = coefIVTilde[1]
        phiTildeIV = coefIVTilde[2]
        rhoPhiIV = 2/tau*(1-log( (yT + alphaTildeIV/phiTildeIV)/(y0 + alphaTildeIV/phiTildeIV) )/(tau*phiTildeIV) )
        coefMIV[m,] = coefIVTilde*(1-tau*rhoPhiIV/2)
        
        # coefficient notTilde by resampling for LL
        y0 = mean(y0yT[i,,1])
        yT = mean(y0yT[i,,2])
        
        coefLLTilde = mvrnorm(n=1,mu=coefLL,Sigma=covBeta)
        alphaTildeLL = coefLLTilde[1]
        phiTildeLL = coefLLTilde[2]
        rhoPhiLL = 2/tau*(1-log( (yT + alphaTildeLL/phiTildeLL)/(y0 + alphaTildeLL/phiTildeLL) )/(tau*phiTildeLL) )
        coefMNpLLBoot[i,m,] = coefLLTilde*(1-tau*rhoPhiLL/2)
    }
    coefMLL = coefMNpLLBoot[i,,]
    
    errRelLM = (coefMLM-coefRep)/coefRep
    errRelLMNaive = (coefMLMNaive-coefRep)/coefRep
    errRelIV = (coefMIV-coefRep)/coefRep
    errRelLL = (coefMLL-coefRep)/coefRep
    
    coefNpLM[i,] = apply(errRelLM,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
    lowerNpLM[i,] = apply(errRelLM,MARGIN=2,FUN=function(x) quantile(x,0.05,na.rm=T))
    upperNpLM[i,] = apply(errRelLM,MARGIN=2,FUN=function(x) quantile(x,0.95,na.rm=T))
    seNpLM[i,] = apply(coefMLM,MARGIN=2,FUN=function(x) sd(x,na.rm=T))
    
    coefNpLMNaive[i,] = apply(errRelLMNaive,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
    lowerNpLMNaive[i,] = apply(errRelLMNaive,MARGIN=2,FUN=function(x) quantile(x,0.05,na.rm=T))
    upperNpLMNaive[i,] = apply(errRelLMNaive,MARGIN=2,FUN=function(x) quantile(x,0.95,na.rm=T))
    seNpLMNaive[i,] = apply(coefMLMNaive,MARGIN=2,FUN=function(x) sd(x,na.rm=T))
    
    coefNpIV[i,] = apply(errRelIV,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
    lowerNpIV[i,] = apply(errRelIV,MARGIN=2,FUN=function(x) quantile(x,0.05,na.rm=T))
    upperNpIV[i,] = apply(errRelIV,MARGIN=2,FUN=function(x) quantile(x,0.95,na.rm=T))
    seNpIV[i,] = apply(coefMIV,MARGIN=2,FUN=function(x) sd(x,na.rm=T))
    
    coefNpLL[i,] = apply(errRelLL,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
    lowerNpLL[i,] = apply(errRelLL,MARGIN=2,FUN=function(x) quantile(x,0.05,na.rm=T))
    upperNpLL[i,] = apply(errRelLL,MARGIN=2,FUN=function(x) quantile(x,0.95,na.rm=T))
    seNpLL[i,] = apply(coefMLL,MARGIN=2,FUN=function(x) sd(x,na.rm=T))
    

}
maxExp = 0
f <- function(x) { return (sign(x)*(log10(1+abs(x)/(10^maxExp)))) }
coefNpLM = apply(coefNpLM,MARGIN=c(1,2),FUN=f)
lowerNpLM = apply(lowerNpLM,MARGIN=c(1,2),FUN=f)
upperNpLM = apply(upperNpLM,MARGIN=c(1,2),FUN=f)
coefNpLMNaive = apply(coefNpLMNaive,MARGIN=c(1,2),FUN=f)
lowerNpLMNaive = apply(lowerNpLMNaive,MARGIN=c(1,2),FUN=f)
upperNpLMNaive = apply(upperNpLMNaive,MARGIN=c(1,2),FUN=f)
coefNpIV = apply(coefNpIV,MARGIN=c(1,2),FUN=f)
lowerNpIV = apply(lowerNpIV,MARGIN=c(1,2),FUN=f)
upperNpIV = apply(upperNpIV,MARGIN=c(1,2),FUN=f)
coefNpLL = apply(coefNpLL,MARGIN=c(1,2),FUN=f)
lowerNpLL = apply(lowerNpLL,MARGIN=c(1,2),FUN=f)
upperNpLL = apply(upperNpLL,MARGIN=c(1,2),FUN=f)

## plot ----
DeltaX = 1/sqrt(Np)
dev.new()
plot(DeltaX,coefNpLM[,1],type="b",pch=1,ylab="",xlab="",xaxt = "n",lwd=1.5,ylim=range(coefNpLM[,1],lowerNpLM[,1],upperNpLM[,1],coefNpLMNaive[,1],lowerNpLMNaive[,1],upperNpLMNaive[,1],coefNpIV[,1],lowerNpIV[,1],upperNpIV[,1],coefNpLL[,1],lowerNpLL[,1],upperNpLL[,1]))
lines(DeltaX,lowerNpLM[,1],lty="dotted",col="black",xaxt="n")
lines(DeltaX,upperNpLM[,1],lty="dotted",col="black",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLM[,1],rev(upperNpLM[,1])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLMNaive[,1],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(DeltaX,upperNpLMNaive[,1],lty="dotted",col="red",xaxt="n")
lines(DeltaX,lowerNpLMNaive[,1],lty="dotted",col="red",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLMNaive[,1],rev(upperNpLMNaive[,1])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpIV[,1],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(DeltaX,upperNpIV[,1],lty="dotted",col="blue",xaxt="n")
lines(DeltaX,lowerNpIV[,1],lty="dotted",col="blue",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpIV[,1],rev(upperNpIV[,1])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLL[,1],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(DeltaX,lowerNpLL[,1],lty="dotted",col="purple",xaxt="n")
lines(DeltaX,upperNpLL[,1],lty="dotted",col="purple",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLL[,1],rev(upperNpLL[,1])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
mapply(axis, side = 1, at = DeltaX, labels = paste(c("(WORST)\n","","","","(A) "), round(DeltaX,digits=3)," \n (N = ",Np ,")",sep=""),las = 3, adj=0.5, cex.axis=1)
abline(h=0)
grid()
title(ylab = TeX(r'(${(\hat{\alpha} - \alpha)/{\alpha}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
title(xlab = TeX(r'(                                                                                                                                       $\frac{1}{\sqrt{N}}$)'),mgp=c(0.5,0,0))
legend("topright",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_alpha.pdf")

dev.new()
plot(DeltaX,coefNpLM[,2],type="b",pch=1,ylab="",xlab="",xaxt = "n",lwd=1.5,ylim=range(coefNpLM[,2],lowerNpLM[,2],upperNpLM[,2],coefNpLMNaive[,2],lowerNpLMNaive[,2],upperNpLMNaive[,2],coefNpIV[,2],lowerNpIV[,2],upperNpIV[,2],coefNpLL[,2],lowerNpLL[,2],upperNpLL[,2]))
lines(DeltaX,lowerNpLM[,2],lty="dotted",col="black",xaxt="n")
lines(DeltaX,upperNpLM[,2],lty="dotted",col="black",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLM[,2],rev(upperNpLM[,2])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLMNaive[,2],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(DeltaX,upperNpLMNaive[,2],lty="dotted",col="red",xaxt="n")
lines(DeltaX,lowerNpLMNaive[,2],lty="dotted",col="red",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLMNaive[,2],rev(upperNpLMNaive[,2])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpIV[,2],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(DeltaX,upperNpIV[,2],lty="dotted",col="blue",xaxt="n")
lines(DeltaX,lowerNpIV[,2],lty="dotted",col="blue",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpIV[,2],rev(upperNpIV[,2])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLL[,2],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(DeltaX,lowerNpLL[,2],lty="dotted",col="purple",xaxt="n")
lines(DeltaX,upperNpLL[,2],lty="dotted",col="purple",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLL[,2],rev(upperNpLL[,2])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
mapply(axis, side = 1, at = DeltaX, labels = paste(c("(WORST)\n","","","","(A) "), round(DeltaX,digits=3)," \n (N = ",Np ,")",sep=""),las = 3, adj=0.5, cex.axis=1)
abline(h=0)
grid()
title(ylab = TeX(r'(${(\hat{\varphi} - \varphi)/{\varphi}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
title(xlab = TeX(r'(                                                                                                                                       $\frac{1}{\sqrt{N}}$)'),mgp=c(0.5,0,0))
legend("topright",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_phi.pdf")

dev.new()
plot(DeltaX,coefNpLM[,3],type="b",pch=1,ylab="",xlab="",xaxt = "n",lwd=1.5,ylim=range(coefNpLM[,3],lowerNpLM[,3],upperNpLM[,3],coefNpLMNaive[,3],lowerNpLMNaive[,3],upperNpLMNaive[,3],coefNpIV[,3],lowerNpIV[,3],upperNpIV[,3],coefNpLL[,3],lowerNpLL[,3],upperNpLL[,3]))
lines(DeltaX,lowerNpLM[,3],lty="dotted",col="black",xaxt="n")
lines(DeltaX,upperNpLM[,3],lty="dotted",col="black",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLM[,3],rev(upperNpLM[,3])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLMNaive[,3],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(DeltaX,upperNpLMNaive[,3],lty="dotted",col="red",xaxt="n")
lines(DeltaX,lowerNpLMNaive[,3],lty="dotted",col="red",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLMNaive[,3],rev(upperNpLMNaive[,3])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpIV[,3],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(DeltaX,upperNpIV[,3],lty="dotted",col="blue",xaxt="n")
lines(DeltaX,lowerNpIV[,3],lty="dotted",col="blue",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpIV[,3],rev(upperNpIV[,3])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLL[,3],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(DeltaX,lowerNpLL[,3],lty="dotted",col="purple",xaxt="n")
lines(DeltaX,upperNpLL[,3],lty="dotted",col="purple",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLL[,3],rev(upperNpLL[,3])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
mapply(axis, side = 1, at = DeltaX, labels = paste(c("(WORST)\n","","","","(A) "), round(DeltaX,digits=3)," \n (N = ",Np ,")",sep=""),las = 3, adj=0.5, cex.axis=1)
abline(h=0)
grid()
title(ylab = TeX(r'(${(\hat{\gamma}_A - \gamma_A)/{\gamma_A}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
title(xlab = TeX(r'(                                                                                                                                       $\frac{1}{\sqrt{N}}$)'),mgp=c(0.5,0,0))
legend("topright",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_gammaA.pdf")

dev.new()
plot(DeltaX,coefNpLM[,4],type="b",pch=1,ylab="",xlab="",xaxt = "n",lwd=1.5,ylim=range(coefNpLM[,4],lowerNpLM[,4],upperNpLM[,4],coefNpLMNaive[,4],lowerNpLMNaive[,4],upperNpLMNaive[,4],coefNpIV[,4],lowerNpIV[,4],upperNpIV[,4],coefNpLL[,4],lowerNpLL[,4],upperNpLL[,4]))
lines(DeltaX,lowerNpLM[,4],lty="dotted",col="black",xaxt="n")
lines(DeltaX,upperNpLM[,4],lty="dotted",col="black",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLM[,4],rev(upperNpLM[,4])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLMNaive[,4],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(DeltaX,upperNpLMNaive[,4],lty="dotted",col="red",xaxt="n")
lines(DeltaX,lowerNpLMNaive[,4],lty="dotted",col="red",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLMNaive[,4],rev(upperNpLMNaive[,4])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpIV[,4],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(DeltaX,upperNpIV[,4],lty="dotted",col="blue",xaxt="n")
lines(DeltaX,lowerNpIV[,4],lty="dotted",col="blue",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpIV[,4],rev(upperNpIV[,4])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLL[,4],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(DeltaX,lowerNpLL[,4],lty="dotted",col="purple",xaxt="n")
lines(DeltaX,upperNpLL[,4],lty="dotted",col="purple",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLL[,4],rev(upperNpLL[,4])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
mapply(axis, side = 1, at = DeltaX, labels = paste(c("(WORST)\n","","","","(A) "), round(DeltaX,digits=3)," \n (N = ",Np ,")",sep=""),las = 3, adj=0.5, cex.axis=1)
abline(h=0)
grid()
title(ylab = TeX(r'(${(\hat{\gamma}_R - \gamma_R)/{\gamma_R}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
title(xlab = TeX(r'(                                                                                                                                       $\frac{1}{\sqrt{N}}$)'),mgp=c(0.5,0,0))
legend("topright",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_gammaR.pdf")

dev.new()
plot(DeltaX,coefNpLM[,5],type="b",pch=1,ylab="",xlab="",xaxt = "n",lwd=1.5,ylim=range(coefNpLM[,5],lowerNpLM[,5],upperNpLM[,5],coefNpLMNaive[,5],lowerNpLMNaive[,5],upperNpLMNaive[,5],coefNpIV[,5],lowerNpIV[,5],upperNpIV[,5],coefNpLL[,5],lowerNpLL[,5],upperNpLL[,5]))
lines(DeltaX,lowerNpLM[,5],lty="dotted",col="black",xaxt="n")
lines(DeltaX,upperNpLM[,5],lty="dotted",col="black",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLM[,5],rev(upperNpLM[,5])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLMNaive[,5],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(DeltaX,upperNpLMNaive[,5],lty="dotted",col="red",xaxt="n")
lines(DeltaX,lowerNpLMNaive[,5],lty="dotted",col="red",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLMNaive[,5],rev(upperNpLMNaive[,5])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpIV[,5],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(DeltaX,upperNpIV[,5],lty="dotted",col="blue",xaxt="n")
lines(DeltaX,lowerNpIV[,5],lty="dotted",col="blue",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpIV[,5],rev(upperNpIV[,5])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(DeltaX,coefNpLL[,5],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(DeltaX,lowerNpLL[,5],lty="dotted",col="purple",xaxt="n")
lines(DeltaX,upperNpLL[,5],lty="dotted",col="purple",xaxt="n")
polygon(c(DeltaX, rev(DeltaX)), c(lowerNpLL[,5],rev(upperNpLL[,5])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
mapply(axis, side = 1, at = DeltaX, labels = paste(c("(WORST)\n","","","","(A) "), round(DeltaX,digits=3)," \n (N = ",Np ,")",sep=""),las = 3, adj=0.5, cex.axis=1)
abline(h=0)
grid()
title(ylab = TeX(r'(${(\hat{\gamma}_D - \gamma_D)/{\gamma_D}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
title(xlab = TeX(r'(                                                                                                                                       $\frac{1}{\sqrt{N}}$)'),mgp=c(0.5,0,0))
legend("topright",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_gammaD.pdf")


# montecarlo PDE different Tau ----
Np = 144
TauAll = c(0.1,0.25,0.5,0.75,1.0)
Nm = 1000
coefMTauLMNaive = array(data=NA,dim=c(length(TauAll),Nm,5))
coefMTauLM = array(data=NA,dim=c(length(TauAll),Nm,5))
coefMTauIV = array(data=NA,dim=c(length(TauAll),Nm,5))
coefTauLL = matrix(data=NA,nrow=length(TauAll),ncol=5)
estimateTauLL = list()
AICR2LMNaive = matrix(data = NA, nrow = length(TauAll), ncol = 2); colnames(AICR2LMNaive) = c("AICc","R2")
AICR2LM = matrix(data = NA, nrow = length(TauAll), ncol = 2); colnames(AICR2LM) = c("AICc","R2")
AICR2IV = matrix(data = NA, nrow = length(TauAll), ncol = 2); colnames(AICR2IV) = c("AICc","R2")
AICR2LL = matrix(data = NA, nrow = length(TauAll), ncol = 2); colnames(AICR2LL) = c("AICc","R2")
y0yT = array(data=NA,dim=c(length(TauAll),Nm,2)); 

#For s
NeS = 100
# SComputed = call_julia_computeS(NeS)
# Xs = SComputed$X
# Ys = SComputed$Y
# S = SComputed$S
Xs = matrix(0,nrow=NeS,ncol=NeS)
Ys = matrix(0,nrow=NeS,ncol=NeS)
S = matrix(0,nrow=NeS,ncol=NeS)

shpMC = createShape(Np, typeOfDist = "Uniform")

for (i in 1:length(TauAll)){
    tau = TauAll[i]
    print(i)
    print(tau)
    
    PDEAll = call_julia_computePDE(tau,SARDp)
    Xpde = PDEAll$X
    Ypde = PDEAll$Y
    PDE0 = PDEAll$PDE0
    PDET = PDEAll$PDET
    
    data_shp = createDataframePDE(PDE0, PDET, Xpde, Ypde, shpMC, tau, Xs, Ys, S)
    data = data_shp$data
    shp = data_shp$shp_sf
    
    MsDeriv = GFDM(data,torus=TRUE)
    D = compute_D(data,longlat=FALSE,torus=TRUE)
    WhA = compute_WhAR(D,data,SARDp$hA)
    WhR = compute_WhAR(D,data,SARDp$hR)
    xS = runif(Np)
    xA = compute_xAR(data,MsDeriv, WhA)
    xR = compute_xAR(data,MsDeriv, WhR)
    xD = compute_xD(data,MsDeriv)
    MS = matrix(data=0, nrow=nrow(MD),ncol=nrow(MD))
    MA = compute_MARLag(data,MsDeriv,WhA)
    MR = compute_MARLag(data,MsDeriv,WhR)
    MD = compute_MDLag(MsDeriv)
    
    MADelta = as.numeric(MA %*% matrix(data$delta))
    MRDelta = as.numeric(MR %*% matrix(data$delta))
    MDDelta = as.numeric(MD %*% matrix(data$delta))
    
    X = data.frame(xA=xA,xR=xR,xD=xD)
    X = as.matrix(X)
    
    # instruments for IV
    MA2X=as.matrix(MA %*% MA %*% X)
    MR2X=as.matrix(MR %*% MR %*% X) 
    MD2X=as.matrix(MD %*% MD %*% X)
    
    estimateTauLM_NAIVE = estimate_ARD_MC_LM_NAIVE(data,shp,xA,xR,xD,MA,MR,MD)
    AICR2LMNaive[i,] = c(estimateTauLM_NAIVE$AICc,estimateTauLM_NAIVE$R2)
    
    estimateTauLM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
    coefLM = coef(estimateTauLM$LM_est)
    coefLM = c(coefLM[c(1,2)],0,coefLM[c(3,5,7)],0,coefLM[c(4,6,8)])
    WN_SARD_OLS_BYSARD = LogLikAICcR2(data, c(coefLM,0), 10, xS, xA, xR, xD, 
                                      MS, MA, MR, MD, diag(nrow(data)))
    AICR2LM[i,] = c(WN_SARD_OLS_BYSARD$AICc[[1]],WN_SARD_OLS_BYSARD$R2Nagelkerke)
    
    estimateTauIV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)
    AICR2IV[i,] = c(estimateTauIV$AICc,estimateTauIV$R2)
    
    estimateTauLL[[i]] = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
    AICR2LL[i,] = c(estimateTauLL[[i]]$AICc,estimateTauLL[[i]]$R2)
    coefTauLL[i,] = estimateTauLL[[i]]$outARD_3MatEstimate$coef[c(1,2,3,4,5)]

    for (m in 1:Nm){
        # print(m)

        iBoot = sample(1:nrow(data),replace=T)
        datam = data[iBoot,]
        xAm = xA[iBoot]
        xRm = xR[iBoot]
        xDm = xD[iBoot]
        MADeltam = MADelta[iBoot]
        MRDeltam = MRDelta[iBoot]
        MDDeltam = MDDelta[iBoot]
        MA2Xm=MA2X[iBoot,]
        MR2Xm=MR2X[iBoot,]
        MD2Xm=MD2X[iBoot,]
        y0yT[i,m,1] = sum(datam$y0)*1/Np
        y0yT[i,m,2] = sum(datam$yT)*1/Np


        estimatePDE_LMNaivem = estimate_ARD_MC_LM_NAIVEBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
        coefMTauLMNaive[i,m,] = estimatePDE_LMNaivem$coefficients[c(1,2,3,4,5)]

        estimatePDE_LMm = estimate_ARD_MC_LMBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
        coefMTauLM[i,m,] = estimatePDE_LMm$coefficients[c(1,2,3,5,7)]

        estimatePDE_IVm = estimate_ARD_MC_IVBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam,MA2Xm,MR2Xm,MD2Xm)
        coefMTauIV[i,m,] = coef(estimatePDE_IVm)[c(1,2,3,4,5)]
    }
}
save(TauAll,coefMTauLM,coefMTauLMNaive,coefMTauIV,coefTauLL,estimateTauLL,AICR2LMNaive,AICR2LM,AICR2IV,AICR2LL,y0yT,file="../datasets_montecarlo/montecarloPDETauGrid.RData")

## AIC plot ----
load(file="../datasets_montecarlo/montecarloPDETauGrid.RData")

library(latex2exp)
dev.new()
plot(TauAll,(AICR2LM[,1]),type="p",ylab="",pch=1,xlab = TeX(r'($\tau$)'),xaxt="n",lwd=1.5,ylim=range(AICR2LMNaive[,1],AICR2LM[,1],AICR2IV[,1],AICR2LL[,1]))
lines(TauAll,(AICR2LMNaive[,1]),type="p",col="red",pch=2,lwd = 1.5, xaxt="n")
lines(TauAll,(AICR2IV[,1]),type="p",col="blue",pch=3,lwd = 1.5,xaxt="n")
lines(TauAll,(AICR2LL[,1]),type="p",col="purple",pch=4,lwd = 1.5,xaxt="n")
axis(1, at = TauAll, label = c("(B) 0.1","0.25","0.50","0.75","(COARSEST) 1"))
grid()
title(ylab = "AICc",mgp=c(2,0,0))
sapply(TauAll,FUN = function(x) return(abline(v = x,lty=3)) );
legend("bottomleft",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEAICcTau.pdf")

## robust SE ----
library(latex2exp)
load(file="../datasets_montecarlo/montecarloPDETauGrid.RData")
Nm = dim(coefMTauLM)[2]
coefRep = matrix(rep(c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),Nm),nrow=Nm,byrow = T)
coefTauLM = matrix(data=NA,nrow=length(TauAll),ncol=5)
lowerTauLM = matrix(data=NA,nrow=length(TauAll),ncol=5)
upperTauLM = matrix(data=NA,nrow=length(TauAll),ncol=5)
seTauLM = matrix(data=NA,nrow=length(TauAll),ncol=5)
coefTauLMNaive = matrix(data=NA,nrow=length(TauAll),ncol=5)
lowerTauLMNaive = matrix(data=NA,nrow=length(TauAll),ncol=5)
upperTauLMNaive = matrix(data=NA,nrow=length(TauAll),ncol=5)
seTauLMNaive = matrix(data=NA,nrow=length(TauAll),ncol=5)
coefTauIV = matrix(data=NA,nrow=length(TauAll),ncol=5)
lowerTauIV = matrix(data=NA,nrow=length(TauAll),ncol=5)
upperTauIV = matrix(data=NA,nrow=length(TauAll),ncol=5)
seTauIV = matrix(data=NA,nrow=length(TauAll),ncol=5)
lowerTauLL = matrix(data=NA,nrow=length(TauAll),ncol=5)
upperTauLL = matrix(data=NA,nrow=length(TauAll),ncol=5)
seTauLL = matrix(data=NA,nrow=length(TauAll),ncol=5)

# new for nonTilde coefficient LL
coefMTauLLBoot = array(data=NA,dim=c(length(TauAll),Nm,5))

for (i in 1:length(TauAll)){
    coefMLM = coefMTauLM[i,,]
    coefMLMNaive = coefMTauLMNaive[i,,]    
    coefMIV = coefMTauIV[i,,]   
    
    coefLL = coefTauLL[i,]
    covBeta = estimateTauLL[[i]]$outARD_3MatEstimate$covBeta
    
    for (m in 1:Nm) {
        y0 = y0yT[i,m,1]
        yT = y0yT[i,m,2]
        
        # coefficient notTilde by bootstrap for LM LMNaive and IV
        coefLMTilde = coefMLM[m,]
        alphaTildeLM = coefLMTilde[1]
        phiTildeLM = coefLMTilde[2]
        rhoPhiLM = 2/tau*(1-log( (yT + alphaTildeLM/phiTildeLM)/(y0 + alphaTildeLM/phiTildeLM) )/(tau*phiTildeLM) )
        coefMLM[m,] = coefLMTilde*(1-tau*rhoPhiLM/2)
        
        coefLMNaiveTilde = coefMLMNaive[m,]
        alphaTildeLMNaive = coefLMNaiveTilde[1]
        phiTildeLMNaive = coefLMNaiveTilde[2]
        rhoPhiLMNaive = 2/tau*(1-log( (yT + alphaTildeLMNaive/phiTildeLMNaive)/(y0 + alphaTildeLMNaive/phiTildeLMNaive) )/(tau*phiTildeLMNaive) )
        coefMLMNaive[m,] = coefLMNaiveTilde*(1-tau*rhoPhiLMNaive/2)
        
        coefIVTilde = coefMIV[m,]
        alphaTildeIV = coefIVTilde[1]
        phiTildeIV = coefIVTilde[2]
        rhoPhiIV = 2/tau*(1-log( (yT + alphaTildeIV/phiTildeIV)/(y0 + alphaTildeIV/phiTildeIV) )/(tau*phiTildeIV) )
        coefMIV[m,] = coefIVTilde*(1-tau*rhoPhiIV/2)
        
        # coefficient notTilde by resampling for LL
        y0 = mean(y0yT[i,,1])
        yT = mean(y0yT[i,,2])
        
        coefLLTilde = mvrnorm(n=1,mu=coefLL,Sigma=covBeta)
        alphaTildeLL = coefLLTilde[1]
        phiTildeLL = coefLLTilde[2]
        rhoPhiLL = 2/tau*(1-log( (yT + alphaTildeLL/phiTildeLL)/(y0 + alphaTildeLL/phiTildeLL) )/(tau*phiTildeLL) )
        coefMTauLLBoot[i,m,] = coefLLTilde*(1-tau*rhoPhiLL/2)
    }
    coefMLL = coefMTauLLBoot[i,,]
    
    errRelLM = (coefMLM-coefRep)/coefRep
    errRelLMNaive = (coefMLMNaive-coefRep)/coefRep
    errRelIV = (coefMIV-coefRep)/coefRep
    errRelLL = (coefMLL-coefRep)/coefRep
    
    coefTauLM[i,] = apply(errRelLM,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
    lowerTauLM[i,] = apply(errRelLM,MARGIN=2,FUN=function(x) quantile(x,0.05,na.rm=T))
    upperTauLM[i,] = apply(errRelLM,MARGIN=2,FUN=function(x) quantile(x,0.95,na.rm=T))
    seTauLM[i,] = apply(coefMLM,MARGIN=2,FUN=function(x) sd(x,na.rm=T))
    
    coefTauLMNaive[i,] = apply(errRelLMNaive,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
    lowerTauLMNaive[i,] = apply(errRelLMNaive,MARGIN=2,FUN=function(x) quantile(x,0.05,na.rm=T))
    upperTauLMNaive[i,] = apply(errRelLMNaive,MARGIN=2,FUN=function(x) quantile(x,0.95,na.rm=T))
    seTauLMNaive[i,] = apply(coefMLMNaive,MARGIN=2,FUN=function(x) sd(x,na.rm=T))
    
    coefTauIV[i,] = apply(errRelIV,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
    lowerTauIV[i,] = apply(errRelIV,MARGIN=2,FUN=function(x) quantile(x,0.05,na.rm=T))
    upperTauIV[i,] = apply(errRelIV,MARGIN=2,FUN=function(x) quantile(x,0.95,na.rm=T))
    seTauIV[i,] = apply(coefMIV,MARGIN=2,FUN=function(x) sd(x,na.rm=T))
    
    coefTauLL[i,] = apply(errRelLL,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
    lowerTauLL[i,] = apply(errRelLL,MARGIN=2,FUN=function(x) quantile(x,0.05,na.rm=T))
    upperTauLL[i,] = apply(errRelLL,MARGIN=2,FUN=function(x) quantile(x,0.95,na.rm=T))
    seTauLL[i,] = apply(coefMLL,MARGIN=2,FUN=function(x) sd(x,na.rm=T))
    
    
}
maxExp = 0
f <- function(x) { return (sign(x)*(log10(1+abs(x)/(10^maxExp)))) }


coefTauLM = apply(coefTauLM,MARGIN=c(1,2),FUN=f)
lowerTauLM = apply(lowerTauLM,MARGIN=c(1,2),FUN=f)
upperTauLM = apply(upperTauLM,MARGIN=c(1,2),FUN=f)
coefTauLMNaive = apply(coefTauLMNaive,MARGIN=c(1,2),FUN=f)
lowerTauLMNaive = apply(lowerTauLMNaive,MARGIN=c(1,2),FUN=f)
upperTauLMNaive = apply(upperTauLMNaive,MARGIN=c(1,2),FUN=f)
coefTauIV = apply(coefTauIV,MARGIN=c(1,2),FUN=f)
lowerTauIV = apply(lowerTauIV,MARGIN=c(1,2),FUN=f)
upperTauIV = apply(upperTauIV,MARGIN=c(1,2),FUN=f)
coefTauLL = apply(coefTauLL,MARGIN=c(1,2),FUN=f)
lowerTauLL = apply(lowerTauLL,MARGIN=c(1,2),FUN=f)
upperTauLL = apply(upperTauLL,MARGIN=c(1,2),FUN=f)

## plot ----

dev.new()
plot(TauAll,coefTauLM[,1],type="b",pch=1,ylab="",xlab = TeX(r'($\tau$)'), xaxt = "n",lwd=1.5,ylim=range(coefTauLM[,1],lowerTauLM[,1],upperTauLM[,1],coefTauLMNaive[,1],lowerTauLMNaive[,1],upperTauLMNaive[,1],coefTauIV[,1],lowerTauIV[,1],upperTauIV[,1],coefTauLL[,1],lowerTauLL[,1],upperTauLL[,1]))
lines(TauAll,lowerTauLM[,1],lty="dotted",col="black",xaxt="n")
lines(TauAll,upperTauLM[,1],lty="dotted",col="black",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLM[,1],rev(upperTauLM[,1])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLMNaive[,1],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(TauAll,upperTauLMNaive[,1],lty="dotted",col="red",xaxt="n")
lines(TauAll,lowerTauLMNaive[,1],lty="dotted",col="red",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLMNaive[,1],rev(upperTauLMNaive[,1])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauIV[,1],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(TauAll,upperTauIV[,1],lty="dotted",col="blue",xaxt="n")
lines(TauAll,lowerTauIV[,1],lty="dotted",col="blue",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauIV[,1],rev(upperTauIV[,1])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLL[,1],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(TauAll,lowerTauLL[,1],lty="dotted",col="purple",xaxt="n")
lines(TauAll,upperTauLL[,1],lty="dotted",col="purple",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLL[,1],rev(upperTauLL[,1])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
axis(1, at = TauAll, label = c("(B) 0.1","0.25","0.50","0.75","(WORST) 1"))
grid()
abline(h=0)
title(ylab = TeX(r'(${(\hat{\alpha} - \alpha)/{\alpha}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
legend("bottomleft",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_alphaTau.pdf")

dev.new()
plot(TauAll,coefTauLM[,2],type="b",pch=1,ylab="",xlab = TeX(r'($\tau$)'), xaxt = "n",lwd=1.5,ylim=range(coefTauLM[,2],lowerTauLM[,2],upperTauLM[,2],coefTauLMNaive[,2],lowerTauLMNaive[,2],upperTauLMNaive[,2],coefTauIV[,2],lowerTauIV[,2],upperTauIV[,2],coefTauLL[,2],lowerTauLL[,2],upperTauLL[,2]))
lines(TauAll,lowerTauLM[,2],lty="dotted",col="black",xaxt="n")
lines(TauAll,upperTauLM[,2],lty="dotted",col="black",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLM[,2],rev(upperTauLM[,2])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLMNaive[,2],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(TauAll,upperTauLMNaive[,2],lty="dotted",col="red",xaxt="n")
lines(TauAll,lowerTauLMNaive[,2],lty="dotted",col="red",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLMNaive[,2],rev(upperTauLMNaive[,2])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauIV[,2],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(TauAll,upperTauIV[,2],lty="dotted",col="blue",xaxt="n")
lines(TauAll,lowerTauIV[,2],lty="dotted",col="blue",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauIV[,2],rev(upperTauIV[,2])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLL[,2],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(TauAll,lowerTauLL[,2],lty="dotted",col="purple",xaxt="n")
lines(TauAll,upperTauLL[,2],lty="dotted",col="purple",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLL[,2],rev(upperTauLL[,2])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
axis(1, at = TauAll, label = c("(B) 0.1","0.25","0.50","0.75","(WORST) 1"))
grid()
abline(h=0)
title(ylab = TeX(r'(${(\hat{\varphi} - \varphi)/{\varphi}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
legend("bottomleft",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_phiTau.pdf")

dev.new()
plot(TauAll,coefTauLM[,3],type="b",pch=1,ylab="",xlab = TeX(r'($\tau$)'), xaxt = "n",lwd=1.5,ylim=range(coefTauLM[,3],lowerTauLM[,3],upperTauLM[,3],coefTauLMNaive[,3],lowerTauLMNaive[,3],upperTauLMNaive[,3],coefTauIV[,3],lowerTauIV[,3],upperTauIV[,3],coefTauLL[,3],lowerTauLL[,3],upperTauLL[,3]))
lines(TauAll,lowerTauLM[,3],lty="dotted",col="black",xaxt="n")
lines(TauAll,upperTauLM[,3],lty="dotted",col="black",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLM[,3],rev(upperTauLM[,3])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLMNaive[,3],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(TauAll,upperTauLMNaive[,3],lty="dotted",col="red",xaxt="n")
lines(TauAll,lowerTauLMNaive[,3],lty="dotted",col="red",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLMNaive[,3],rev(upperTauLMNaive[,3])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauIV[,3],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(TauAll,upperTauIV[,3],lty="dotted",col="blue",xaxt="n")
lines(TauAll,lowerTauIV[,3],lty="dotted",col="blue",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauIV[,3],rev(upperTauIV[,3])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLL[,3],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(TauAll,lowerTauLL[,3],lty="dotted",col="purple",xaxt="n")
lines(TauAll,upperTauLL[,3],lty="dotted",col="purple",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLL[,3],rev(upperTauLL[,3])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
axis(1, at = TauAll, label = c("(B) 0.1","0.25","0.50","0.75","(WORST) 1"))
grid()
abline(h=0)
title(ylab = TeX(r'(${(\hat{\gamma}_A - \gamma_A)/{\gamma_A}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
legend("bottomleft",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_gammaATau.pdf")

dev.new()
plot(TauAll,coefTauLM[,4],type="b",pch=1,ylab="",xlab = TeX(r'($\tau$)'), xaxt = "n",lwd=1.5,ylim=range(coefTauLM[,4],lowerTauLM[,4],upperTauLM[,4],coefTauLMNaive[,4],lowerTauLMNaive[,4],upperTauLMNaive[,4],coefTauIV[,4],lowerTauIV[,4],upperTauIV[,4],coefTauLL[,4],lowerTauLL[,4],upperTauLL[,4]))
lines(TauAll,lowerTauLM[,4],lty="dotted",col="black",xaxt="n")
lines(TauAll,upperTauLM[,4],lty="dotted",col="black",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLM[,4],rev(upperTauLM[,4])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLMNaive[,4],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(TauAll,upperTauLMNaive[,4],lty="dotted",col="red",xaxt="n")
lines(TauAll,lowerTauLMNaive[,4],lty="dotted",col="red",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLMNaive[,4],rev(upperTauLMNaive[,4])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauIV[,4],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(TauAll,upperTauIV[,4],lty="dotted",col="blue",xaxt="n")
lines(TauAll,lowerTauIV[,4],lty="dotted",col="blue",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauIV[,4],rev(upperTauIV[,4])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLL[,4],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(TauAll,lowerTauLL[,4],lty="dotted",col="purple",xaxt="n")
lines(TauAll,upperTauLL[,4],lty="dotted",col="purple",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLL[,4],rev(upperTauLL[,4])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
axis(1, at = TauAll, label = c("(B) 0.1","0.25","0.50","0.75","(WORST) 1"))
grid()
abline(h=0)
title(ylab = TeX(r'(${(\hat{\gamma}_R - \gamma_R)/{\gamma_R}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
legend("bottomleft",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_gammaRTau.pdf")

dev.new()
plot(TauAll,coefTauLM[,5],type="b",pch=1,ylab="",xlab = TeX(r'($\tau$)'), xaxt = "n",lwd=1.5,ylim=range(coefTauLM[,5],lowerTauLM[,5],upperTauLM[,5],coefTauLMNaive[,5],lowerTauLMNaive[,5],upperTauLMNaive[,5],coefTauIV[,5],lowerTauIV[,5],upperTauIV[,5],coefTauLL[,5],lowerTauLL[,5],upperTauLL[,5]))
lines(TauAll,lowerTauLM[,5],lty="dotted",col="black",xaxt="n")
lines(TauAll,upperTauLM[,5],lty="dotted",col="black",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLM[,5],rev(upperTauLM[,5])), col = adjustcolor("black", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLMNaive[,5],col="red",xaxt="n",type="b",pch=2,lwd=1.5)
lines(TauAll,upperTauLMNaive[,5],lty="dotted",col="red",xaxt="n")
lines(TauAll,lowerTauLMNaive[,5],lty="dotted",col="red",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLMNaive[,5],rev(upperTauLMNaive[,5])), col = adjustcolor("red", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauIV[,5],col="blue",xaxt="n",type="b",pch=3,lwd=1.5)
lines(TauAll,upperTauIV[,5],lty="dotted",col="blue",xaxt="n")
lines(TauAll,lowerTauIV[,5],lty="dotted",col="blue",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauIV[,5],rev(upperTauIV[,5])), col = adjustcolor("blue", alpha.f=0.1) , lty = 0)
lines(TauAll,coefTauLL[,5],col="purple",xaxt="n",type="b",pch=4,lwd=1.5)
lines(TauAll,lowerTauLL[,5],lty="dotted",col="purple",xaxt="n")
lines(TauAll,upperTauLL[,5],lty="dotted",col="purple",xaxt="n")
polygon(c(TauAll, rev(TauAll)), c(lowerTauLL[,5],rev(upperTauLL[,5])), col = adjustcolor("purple", alpha.f=0.1) , lty = 0)
axis(1, at = TauAll, label = c("(B) 0.1","0.25","0.50","0.75","(WORST) 1"))
grid()
abline(h=0)
title(ylab = TeX(r'(${(\hat{\gamma}_D - \gamma_D)/{\gamma_D}}$ (Symmetric $\log_{10}$ scale) )'),mgp=c(2,0,0))
legend("bottomleft",legend=c("OLS NAIVE","OLS","IV","ML"),col=c("red","black","blue","purple"),lty=1,lwd=1.5,cex=1.0,pch=c(2,1,3,4))
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRel_gammaDTau.pdf")


# WORST CASE ----
Na=100
Nm = 1000
tau= 1.0
NeS = 100
SARDp = list(alpha = 0.01, phi = 0.01, gammaS = 0.0, gammaA = -0.00175, gammaR = 0.0025, gammaD = 0.00525, hA = 0.15, hR = 0.4)

# PDE once 
PDEAll = call_julia_computePDE(tau,SARDp)

# Np = 144 
Np = 144

## create PDE shape 
shpMC = createShape(Np, typeOfDist = "Uniform")
Xpde = PDEAll$X
Ypde = PDEAll$Y
PDE0 = PDEAll$PDE0
PDET = PDEAll$PDET

#For s
# SComputed = call_julia_computeS(NeS)
# Xs = SComputed$X
# Ys = SComputed$Y
# S = SComputed$S
Xs = matrix(0,nrow=NeS,ncol=NeS)
Ys = matrix(0,nrow=NeS,ncol=NeS)
S = matrix(0,nrow=NeS,ncol=NeS)

data_shp = createDataframePDE(PDE0, PDET, Xpde, Ypde, shpMC, tau, Xs, Ys, S)
data = data_shp$data
shp = data_shp$shp_sf


## create regressors
MsDeriv = GFDM(data,torus=TRUE)
D = compute_D(data,longlat=FALSE,torus=TRUE)
WhA = compute_WhAR(D,data,SARDp$hA)
WhR = compute_WhAR(D,data,SARDp$hR)x
xS = runif(Np)
xA = compute_xAR(data,MsDeriv, WhA)
xR = compute_xAR(data,MsDeriv, WhR)
xD = compute_xD(data,MsDeriv)
MS = matrix(data=0, nrow=nrow(MD),ncol=nrow(MD))
MA = compute_MARLag(data,MsDeriv,WhA)
MR = compute_MARLag(data,MsDeriv,WhR)
MD = compute_MDLag(MsDeriv)
shp = cbind(shp,xA,xR,xD)
plot(shp[c("y0","yT","delta","xA","xR","xD")])


## estimate PDE LL
estimatePDE_LM_NAIVE = estimate_ARD_MC_LM_NAIVE(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE_LM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
coefLM = coef(estimatePDE_LM$LM_est)
coefLM = c(coefLM[c(1,2)],0,coefLM[c(3,5,7)],0,coefLM[c(4,6,8)])
WN_SARD_OLS_BYSARD = LogLikAICcR2(data, c(coefLM,0), 10, xS, xA, xR, xD, MS, MA, MR, MD, diag(nrow(data)))


estimatePDE_IV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE_LL = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
save(estimatePDE_LM_NAIVE,estimatePDE_LM,estimatePDE_IV,estimatePDE_LL, file = "../datasets_montecarlo/EstimateWORST.RData")

load("../datasets_montecarlo/EstimateWORST.RData")

coefLL = estimatePDE_LL$outARD_3MatEstimate$coef[c(1,2,3,4,5)]
covBeta = estimatePDE_LL$outARD_3MatEstimate$covBeta
y0yT = matrix(data = NA, nrow = Nm, ncol = 2); 
coefMLMNaive = matrix(data = NA, nrow = Nm, ncol = 5)
coefMLM = matrix(data = NA, nrow = Nm, ncol = 5)
coefMIV = matrix(data = NA, nrow = Nm, ncol = 5)
coefMLL = matrix(data = NA, nrow = Nm, ncol = 5)

X = data.frame(xA=xA,xR=xR,xD=xD)
X = as.matrix(X)
MA2X=as.matrix(MA %*% MA %*% X)
MR2X=as.matrix(MR %*% MR %*% X) 
MD2X=as.matrix(MD %*% MD %*% X)
MADelta = as.numeric(MA %*% matrix(data$delta))
MRDelta = as.numeric(MR %*% matrix(data$delta))
MDDelta = as.numeric(MD %*% matrix(data$delta))
set.seed(1)
for (m in 1:Nm){
    # print(m)
    
    iBoot = sample(1:nrow(data),replace=T)
    datam = data[iBoot,]
    xAm = xA[iBoot]
    xRm = xR[iBoot]
    xDm = xD[iBoot]
    MADeltam = MADelta[iBoot]
    MRDeltam = MRDelta[iBoot]
    MDDeltam = MDDelta[iBoot]
    MA2Xm=MA2X[iBoot,]
    MR2Xm=MR2X[iBoot,]
    MD2Xm=MD2X[iBoot,]
    y0yT[m,1] = sum(datam$y0)*1/Np
    y0yT[m,2] = sum(datam$yT)*1/Np
    
    estimatePDE_LMm = estimate_ARD_MC_LMBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
    coefMLM[m,] = estimatePDE_LMm$coefficients[c(1,2,3,5,7)]
    
    estimatePDE_LMNaivem = estimate_ARD_MC_LM_NAIVEBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
    coefMLMNaive[m,] = estimatePDE_LMNaivem$coefficients[c(1,2,3,4,5)]
    
    estimatePDE_IVm = estimate_ARD_MC_IVBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam,MA2Xm,MR2Xm,MD2Xm)
    coefMIV[m,] = coef(estimatePDE_IVm)[c(1,2,3,4,5)]
    
    # robust SE
    y0 = y0yT[m,1]
    yT = y0yT[m,2]
    
    coefLMTilde = coefMLM[m,]
    alphaTildeLM = coefLMTilde[1]
    phiTildeLM = coefLMTilde[2]
    rhoPhiLM = 2/tau*(1-log( (yT + alphaTildeLM/phiTildeLM)/(y0 + alphaTildeLM/phiTildeLM) )/(tau*phiTildeLM) )
    coefMLM[m,] = coefLMTilde*(1-tau*rhoPhiLM/2)
    
    coefLMNaiveTilde = coefMLMNaive[m,]
    alphaTildeLMNaive = coefLMNaiveTilde[1]
    phiTildeLMNaive = coefLMNaiveTilde[2]
    rhoPhiLMNaive = 2/tau*(1-log( (yT + alphaTildeLMNaive/phiTildeLMNaive)/(y0 + alphaTildeLMNaive/phiTildeLMNaive) )/(tau*phiTildeLMNaive) )
    coefMLMNaive[m,] = coefLMNaiveTilde*(1-tau*rhoPhiLMNaive/2)
    
    coefIVTilde = coefMIV[m,]
    alphaTildeIV = coefIVTilde[1]
    phiTildeIV = coefIVTilde[2]
    rhoPhiIV = 2/tau*(1-log( (yT + alphaTildeIV/phiTildeIV)/(y0 + alphaTildeIV/phiTildeIV) )/(tau*phiTildeIV) )
    coefMIV[m,] = coefIVTilde*(1-tau*rhoPhiIV/2)
    
    # coefficient notTilde by resampling for LL
    y0 = sum(data$y0)*1/Np
    yT = sum(data$yT)*1/Np
    
    coefLLTilde = mvrnorm(n=1,mu=coefLL,Sigma=covBeta)
    alphaTildeLL = coefLLTilde[1]
    phiTildeLL = coefLLTilde[2]
    rhoPhiLL = 2/tau*(1-log( (yT + alphaTildeLL/phiTildeLL)/(y0 + alphaTildeLL/phiTildeLL) )/(tau*phiTildeLL) )
    coefMLL[m,] = coefLLTilde*(1-tau*rhoPhiLL/2)
}

coefLMNaive = apply(coefMLMNaive,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
seLMNaive = apply(coefMLMNaive,MARGIN=2,FUN=function(x) sd(x,na.rm=T))

coefLM = apply(coefMLM,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
seLM = apply(coefMLM,MARGIN=2,FUN=function(x) sd(x,na.rm=T))

coefIV = apply(coefMIV,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
seIV = apply(coefMIV,MARGIN=2,FUN=function(x) sd(x,na.rm=T))

coefNpLL = apply(coefMLL,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
seLL  = apply(coefMLL,MARGIN=2,FUN=function(x) sd(x,na.rm=T))

round(coefLMNaive,digits = 5)
round(coefLMNaive - c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),digits=5)
round(seLMNaive, digits = 5)
round(sqrt(1/Nm * colSums((coefMLMNaive-c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD))^2,na.rm=T)),digits=5)
round(estimatePDE_LM_NAIVE$AICc,digits = 2) 
round(estimatePDE_LM_NAIVE$R2, digits = 5)

round(coefLM,digits = 5)
round(coefLM - c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),digits=5)
round(seLM, digits = 5)
round(sqrt(1/Nm * colSums((coefMLM-c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD))^2,na.rm=T)),digits=5)
round(estimatePDE_LM$R2, digits = 5)
round(estimatePDE_LM$AICc,digits = 2) 

round(coefIV,digits = 5)
round(coefIV - c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),digits=5)
round(seIV, digits = 5)
round(sqrt(1/Nm * colSums((coefMIV-c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD))^2,na.rm=T)),digits=5)
round(estimatePDE_IV$AICc,digits = 2) 
round(estimatePDE_IV$R2, digits = 5)

round(coefLL,digits = 5)
round(coefLL - c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),digits=5)
round(seLL, digits = 5)
round(sqrt(1/Nm * colSums((coefMLL-c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD))^2,na.rm=T)),digits=5)
round(estimatePDE_LL$AICc,digits = 2) 
round(estimatePDE_LL$R2, digits = 5)


# BEST CASE ----
Na=100
Nm = 1000
tau= 0.1
NeS = 100
SARDp = list(alpha = 0.01, phi = 0.01, gammaS = 0.0, gammaA = -0.00175, gammaR = 0.0025, gammaD = 0.00525, hA = 0.15, hR = 0.4)

# PDE once 
PDEAll = call_julia_computePDE(tau,SARDp)

# Np = 144 
Np = 2500

## create PDE shape 
shpMC = createShape(Np, typeOfDist = "Uniform")
Xpde = PDEAll$X
Ypde = PDEAll$Y
PDE0 = PDEAll$PDE0
PDET = PDEAll$PDET

#For s
# SComputed = call_julia_computeS(NeS)
# Xs = SComputed$X
# Ys = SComputed$Y
# S = SComputed$S
Xs = matrix(0,nrow=NeS,ncol=NeS)
Ys = matrix(0,nrow=NeS,ncol=NeS)
S = matrix(0,nrow=NeS,ncol=NeS)


data_shp = createDataframePDE(PDE0, PDET, Xpde, Ypde, shpMC, tau, Xs, Ys, S)
data = data_shp$data
shp = data_shp$shp_sf


## create regressors
MsDeriv = GFDM(data,torus=TRUE)
D = compute_D(data,longlat=FALSE,torus=TRUE)
WhA = compute_WhAR(D,data,SARDp$hA)
WhR = compute_WhAR(D,data,SARDp$hR)
xA = compute_xAR(data,MsDeriv, WhA)
xR = compute_xAR(data,MsDeriv, WhR)
xD = compute_xD(data,MsDeriv)
xS = runif(Np)
MA = compute_MARLag(data,MsDeriv,WhA)
MR = compute_MARLag(data,MsDeriv,WhR)
MD = compute_MDLag(MsDeriv)
MS = matrix(data=0, nrow=nrow(MD),ncol=nrow(MD))
shp = cbind(shp,xA,xR,xD)
plot(shp[c("y0","yT","delta","xA","xR","xD")])


## estimate PDE LL
estimatePDE_LM_NAIVE = estimate_ARD_MC_LM_NAIVE(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE_LM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
coefLM = coef(estimatePDE_LM$LM_est)
coefLM = c(coefLM[c(1,2)],0,coefLM[c(3,5,7)],0,coefLM[c(4,6,8)])
WN_SARD_OLS_BYSARD = LogLikAICcR2(data, c(coefLM,0), 10, xS, xA, xR, xD, MS, MA, MR, MD, diag(nrow(data)))

estimatePDE_IV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE_LL = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
save(estimatePDE_LM_NAIVE,estimatePDE_LM,estimatePDE_IV,estimatePDE_LL, file = "../datasets_montecarlo/EstimateBEST.RData")

load("../datasets_montecarlo/EstimateBEST.RData")
coefLL = estimatePDE_LL$outARD_3MatEstimate$coef[c(1,2,3,4,5)]
covBeta = estimatePDE_LL$outARD_3MatEstimate$covBeta
y0yT = matrix(data = NA, nrow = Nm, ncol = 2); 
coefMLMNaive = matrix(data = NA, nrow = Nm, ncol = 5)
coefMLM = matrix(data = NA, nrow = Nm, ncol = 5)
coefMIV = matrix(data = NA, nrow = Nm, ncol = 5)
coefMLL = matrix(data = NA, nrow = Nm, ncol = 5)

X = data.frame(xA=xA,xR=xR,xD=xD)
X = as.matrix(X)
MA2X=as.matrix(MA %*% MA %*% X)
MR2X=as.matrix(MR %*% MR %*% X) 
MD2X=as.matrix(MD %*% MD %*% X)
MADelta = as.numeric(MA %*% matrix(data$delta))
MRDelta = as.numeric(MR %*% matrix(data$delta))
MDDelta = as.numeric(MD %*% matrix(data$delta))
set.seed(1)
for (m in 1:Nm){
    print(m)
    
    iBoot = sample(1:nrow(data),replace=T)
    datam = data[iBoot,]
    xAm = xA[iBoot]
    xRm = xR[iBoot]
    xDm = xD[iBoot]
    MADeltam = MADelta[iBoot]
    MRDeltam = MRDelta[iBoot]
    MDDeltam = MDDelta[iBoot]
    MA2Xm=MA2X[iBoot,]
    MR2Xm=MR2X[iBoot,]
    MD2Xm=MD2X[iBoot,]
    y0yT[m,1] = sum(datam$y0)*1/Np
    y0yT[m,2] = sum(datam$yT)*1/Np
    
    estimatePDE_LMm = estimate_ARD_MC_LMBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
    coefMLM[m,] = estimatePDE_LMm$coefficients[c(1,2,3,5,7)]
    
    estimatePDE_LMNaivem = estimate_ARD_MC_LM_NAIVEBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
    coefMLMNaive[m,] = estimatePDE_LMNaivem$coefficients[c(1,2,3,4,5)]
    
    estimatePDE_IVm = estimate_ARD_MC_IVBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam,MA2Xm,MR2Xm,MD2Xm)
    coefMIV[m,] = coef(estimatePDE_IVm)[c(1,2,3,4,5)]
    
    # robust SE
    y0 = y0yT[m,1]
    yT = y0yT[m,2]
    
    coefLMTilde = coefMLM[m,]
    alphaTildeLM = coefLMTilde[1]
    phiTildeLM = coefLMTilde[2]
    rhoPhiLM = 2/tau*(1-log( (yT + alphaTildeLM/phiTildeLM)/(y0 + alphaTildeLM/phiTildeLM) )/(tau*phiTildeLM) )
    coefMLM[m,] = coefLMTilde*(1-tau*rhoPhiLM/2)
    
    coefLMNaiveTilde = coefMLMNaive[m,]
    alphaTildeLMNaive = coefLMNaiveTilde[1]
    phiTildeLMNaive = coefLMNaiveTilde[2]
    rhoPhiLMNaive = 2/tau*(1-log( (yT + alphaTildeLMNaive/phiTildeLMNaive)/(y0 + alphaTildeLMNaive/phiTildeLMNaive) )/(tau*phiTildeLMNaive) )
    coefMLMNaive[m,] = coefLMNaiveTilde*(1-tau*rhoPhiLMNaive/2)
    
    coefIVTilde = coefMIV[m,]
    alphaTildeIV = coefIVTilde[1]
    phiTildeIV = coefIVTilde[2]
    rhoPhiIV = 2/tau*(1-log( (yT + alphaTildeIV/phiTildeIV)/(y0 + alphaTildeIV/phiTildeIV) )/(tau*phiTildeIV) )
    coefMIV[m,] = coefIVTilde*(1-tau*rhoPhiIV/2)
    
    # coefficient notTilde by resampling for LL
    y0 = sum(data$y0)*1/Np
    yT = sum(data$yT)*1/Np
    
    coefLLTilde = mvrnorm(n=1,mu=coefLL,Sigma=covBeta)
    alphaTildeLL = coefLLTilde[1]
    phiTildeLL = coefLLTilde[2]
    rhoPhiLL = 2/tau*(1-log( (yT + alphaTildeLL/phiTildeLL)/(y0 + alphaTildeLL/phiTildeLL) )/(tau*phiTildeLL) )
    coefMLL[m,] = coefLLTilde*(1-tau*rhoPhiLL/2)
}

coefLMNaive = apply(coefMLMNaive,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
seLMNaive = apply(coefMLMNaive,MARGIN=2,FUN=function(x) sd(x,na.rm=T))

coefLM = apply(coefMLM,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
seLM = apply(coefMLM,MARGIN=2,FUN=function(x) sd(x,na.rm=T))

coefIV = apply(coefMIV,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
seIV = apply(coefMIV,MARGIN=2,FUN=function(x) sd(x,na.rm=T))

coefNpLL = apply(coefMLL,MARGIN=2,FUN=function(x) mean(x,na.rm=T))
seLL  = apply(coefMLL,MARGIN=2,FUN=function(x) sd(x,na.rm=T))

round(coefLMNaive,digits = 5)
round(coefLMNaive - c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),digits=5)
round(seLMNaive, digits = 5)
round(sqrt(1/Nm * colSums((coefMLMNaive-c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD))^2)),digits=5)
round(estimatePDE_LM_NAIVE$AICc,digits = 2) 
round(estimatePDE_LM_NAIVE$R2, digits = 5)

round(coefLM,digits = 5)
round(coefLM - c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),digits=5)
round(seLM, digits = 5)
round(sqrt(1/Nm * colSums((coefMLM-c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD))^2)),digits=5)
round(estimatePDE_LM$R2, digits = 5)
round(estimatePDE_LM$AICc,digits = 2) 

round(coefIV,digits = 5)
round(coefIV - c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),digits=5)
round(seIV, digits = 5)
round(sqrt(1/Nm * colSums((coefMIV-c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD))^2)),digits=5)
round(estimatePDE_IV$AICc,digits = 2) 
round(estimatePDE_IV$R2, digits = 5)

round(coefLL,digits = 5)
round(coefLL - c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),digits=5)
round(seLL, digits = 5)
round(sqrt(1/Nm * colSums((coefMLL-c(SARDp$alpha,SARDp$phi,SARDp$gammaA,SARDp$gammaR,SARDp$gammaD))^2)),digits=5)
round(estimatePDE_LL$AICc,digits = 2) 
round(estimatePDE_LL$R2, digits = 5)



# Analysis of Wepsilon for different Np and tau----
load(file="../datasets_montecarlo/EstimateWORST.RData")
worstLL = estimatePDE_LL
load(file="../datasets_montecarlo/EstimateBEST.RData")
besttLL = estimatePDE_LL

worstLambdaHat = worstLL$outARD_3MatEstimate$coef[length(worstLL$outARD_3MatEstimate$coef)]
bestLambdaHat = besttLL$outARD_3MatEstimate$coef[length(besttLL$outARD_3MatEstimate$coef)]

worstLambdaHatSE = worstLL$outARD_3MatEstimate$se_coef[length(worstLL$outARD_3MatEstimate$coef)]
bestLambdaHatSE = besttLL$outARD_3MatEstimate$se_coef[length(besttLL$outARD_3MatEstimate$coef)]

worstLambdas = round(coef(worstLL$SpatError$lm.errDecompose), digits=3)
bestLambdas = round(coef(besttLL$SpatError$lm.errDecompose), digits=3)

worstLambdasSE = round(sqrt(diag(vcov((worstLL$SpatError$lm.errDecompose)))), digits=3)
bestLambdasSE = round(sqrt(diag(vcov((besttLL$SpatError$lm.errDecompose)))),  digits=3)


worstMaxLag= worstLL$SpatError$maxSignifLag
bestMaxLag= besttLL$SpatError$maxSignifLag

## Analysis residuals ----
Na=100
Nm = 1000
tau= 1.0
NeS = 100
SARDp = list(alpha = 0.01, phi = 0.01, gammaS = 0.0, gammaA = -0.00175, gammaR = 0.0025, gammaD = 0.00525, hA = 0.15, hR = 0.4)

# PDE once
PDEAll = call_julia_computePDE(tau,SARDp)
save(PDEAll,file="../datasets_montecarlo/PDE.RData")


Np = 144

## create PDE shape 
shpMC = createShape(Np, typeOfDist = "Uniform")
Xpde = PDEAll$X
Ypde = PDEAll$Y
PDE0 = PDEAll$PDE0
PDET = PDEAll$PDET
#For s
Xs = matrix(0,nrow=NeS,ncol=NeS)
Ys = matrix(0,nrow=NeS,ncol=NeS)
S = matrix(0,nrow=NeS,ncol=NeS)

# SComputed = call_julia_computeS(NeS)
# Xs = SComputed$X
# Ys = SComputed$Y
# S = SComputed$S

data_shp = createDataframePDE(PDE0, PDET, Xpde, Ypde, shpMC, tau, Xs, Ys, S)
data = data_shp$data
shp = data_shp$shp_sf


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
shp = cbind(shp,xA,xR,xD)

errorST = as.numeric(shp$delta - SARDp$alpha - SARDp$phi*shp$y0 - SARDp$gammaA*xA - SARDp$gammaR*xR - SARDp$gammaD*xD 
    - tau/2*SARDp$gammaA*MA %*% shp$delta - tau/2*SARDp$gammaR*MR %*% shp$delta - tau/2*SARDp$gammaD*MD %*% shp$delta)

shp = cbind(shp, errorST)
dev.new()
plot(shp["errorST"],main="")
dev.copy2pdf(file="../datasets_montecarlo/residual2500NoFilter.pdf")

dev.new()
plot(shp["errorSTFilter"],main="")
dev.copy2pdf(file="../datasets_montecarlo/residual2500Filter.pdf")


SpatError = compute_spatial_error_mat(errorST,shp,maxLag = 10, pThreshold = 0.1)
Werr = SpatError$Werr


errorSTFilter = as.numeric(errorST - Werr %*% errorST)
shp = cbind(shp, errorSTFilter)
plot(shp["errorST"],breaks = seq(from=-0.4,to=0.4,by=0.05),main="",key.pos=1)
plot(shp["errorSTFilter"],breaks = seq(from=-0.4,to=0.4,by=0.05),main="",key.pos=1)
plot(shp[c("errorST","errorSTFilter")],main="",key.pos=4)
save(shp,SpatError,file="../datasets_montecarlo/MCerror144.RData")


### plot ----
load("../datasets_montecarlo/MCerror2500.RData")
shp2500 = shp
SpatError2500 = SpatError
load("../datasets_montecarlo/MCerror10000.RData")
shp10000 = shp
SpatError10000 = SpatError



qqPre = qqnorm(shp10000$errorST,plot.it = F)
qqPost = qqnorm(shp10000$errorSTFilter,plot.it = F)
plot(qqPre$x,qqPre$y,main="qq 10000")
points(qqPost$x,qqPost$y,col="red")

qq2500 = qqnorm(shp2500$errorSTFilter,plot.it = F)
qq10000 = qqnorm(shp10000$errorSTFilter,plot.it = F)
plot(qq2500$x,qq2500$y,main="qq 2500 vs 10000 Filter")
points(qq10000$x,qq10000$y,col="red")

sm.density(shp2500$errorSTFilter,model="Normal",ylim=c(0,20))
sm.density(shp10000$errorSTFilter,model="Normal",ylim=c(0,5))

library(kldest)
kld_est(shp2500$errorST,q=function(x) dnorm(x,mean=mean(shp2500$errorST),sd=sd(shp2500$errorST)))
kld_est(shp2500$errorSTFilter,q=function(x) dnorm(x,mean=mean(shp2500$errorSTFilter),sd=sd(shp2500$errorSTFilter)))
kld_est(shp10000$errorST,q=function(x) dnorm(x,mean=mean(shp10000$errorST),sd=sd(shp10000$errorST)))
kld_est(shp10000$errorSTFilter,q=function(x) dnorm(x,mean=mean(shp10000$errorSTFilter),sd=sd(shp10000$errorSTFilter)))


## plot correlogram ----
### 2500 ----
load("../datasets_montecarlo/MCerror2500.RData")
shp2500 = shp
maxLag = 10
correlogramm2500 = correlogram(shp2500$errorST,shp2500,maxLag = maxLag)
correlogramm2500Filter = correlogram(shp2500$errorSTFilter,shp2500,maxLag = maxLag)

coeffCorrelogramm2500 = correlogramm2500$correlogram_resid$res[,1]
SECorrelogramm2500 = sqrt(correlogramm2500$correlogram_resid$res[,3])
coeffCorrelogramm2500Filter = correlogramm2500Filter$correlogram_resid$res[,1]
SECorrelogramm2500Filter = sqrt(correlogramm2500Filter$correlogram_resid$res[,3])


dev.new()
plot(1:maxLag,coeffCorrelogramm2500,xlab="Lag",ylab="Moran's I",pch=18,ylim=range(coeffCorrelogramm2500-1.96*SECorrelogramm2500,coeffCorrelogramm2500+1.96*SECorrelogramm2500,coeffCorrelogramm2500Filter-1.96*SECorrelogramm2500Filter,coeffCorrelogramm2500Filter+1.96*SECorrelogramm2500Filter),xlim=c(1,maxLag))
axis(side = 1, at = 1:maxLag)
arrows(1:maxLag,coeffCorrelogramm2500-1.96*SECorrelogramm2500,1:maxLag,coeffCorrelogramm2500+1.96*SECorrelogramm2500,length=0.1,angle=90,code=3)
points(1:maxLag,coeffCorrelogramm2500Filter,xlab="Lag",ylab="Moran's I",pch=18,col="red")
arrows(1:maxLag,coeffCorrelogramm2500Filter-1.96*SECorrelogramm2500Filter,1:maxLag,coeffCorrelogramm2500Filter+1.96*SECorrelogramm2500Filter,length=0.1,angle=90,code=3,col="red")
abline(h=0)
grid()
legend("topright",legend=c("Reminder before filtering","Reminder after filtering"),col=c("black","red"),lty=1,lwd=1.5,cex=1.0)
dev.copy2pdf(file="../datasets_montecarlo/correlogramm2500.pdf")

### 144 ----
load("../datasets_montecarlo/MCerror144.RData")
shp144 = shp
maxLag = 10
correlogramm144 = correlogram(shp144$errorST,shp144,maxLag = maxLag)
correlogramm144Filter = correlogram(shp144$errorSTFilter,shp144,maxLag = maxLag)

coeffCorrelogramm144 = correlogramm144$correlogram_resid$res[,1]
SECorrelogramm144 = sqrt(correlogramm144$correlogram_resid$res[,3])
coeffCorrelogramm144Filter = correlogramm144Filter$correlogram_resid$res[,1]
SECorrelogramm144Filter = sqrt(correlogramm144Filter$correlogram_resid$res[,3])


dev.new()
plot(1:maxLag,coeffCorrelogramm144,xlab="lags",ylab="Moran's I",pch=18,ylim=range(coeffCorrelogramm144-1.96*SECorrelogramm144,coeffCorrelogramm144+1.96*SECorrelogramm144,coeffCorrelogramm144Filter-1.96*SECorrelogramm144Filter,coeffCorrelogramm144Filter+1.96*SECorrelogramm144Filter),xlim=c(1,maxLag))
axis(side = 1, at = 1:maxLag)
arrows(1:maxLag,coeffCorrelogramm144-1.96*SECorrelogramm144,1:maxLag,coeffCorrelogramm144+1.96*SECorrelogramm144,length=0.1,angle=90,code=3)
points(1:maxLag,coeffCorrelogramm144Filter,xlab="lags",ylab="Moran's I",pch=18,col="red")
arrows(1:maxLag,coeffCorrelogramm144Filter-1.96*SECorrelogramm144Filter,1:maxLag,coeffCorrelogramm144Filter+1.96*SECorrelogramm144Filter,length=0.1,angle=90,code=3,col="red")
abline(h=0)
legend("topright",legend=c("reminder","reminder after filtering"),col=c("black","red"),lty=1,lwd=1.5,cex=1.0)
dev.copy2pdf(file="../datasets_montecarlo/correlogramm144.pdf")

