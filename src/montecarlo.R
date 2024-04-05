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
SARDp = list(gammaS = 0.0, gammaA = -0.00175, gammaR = 0.0025, gammaD = 0.00525, hA = 0.15, hR = 0.4)

# PDE once ----
PDEAll = call_julia_computePDE(tau,SARDp)
save(PDEAll,file="../datasets_montecarlo/PDE.RData")

# Np = 144 ----
Np = 144

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
save(data,shp,file="../datasets_montecarlo/PDEdatashp144.RData")


## create regressors ----
load(file="../datasets_montecarlo/PDEdatashp144.RData")

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
save(data,shp,xA,xR,xD,MA,MR,MD,file="../datasets_montecarlo/DataShpRegressors144.RData")


## estimate PDE LL ----
load("../datasets_montecarlo/DataShpRegressors1156.RData")
estimatePDE1156_LM_NAIVE = estimate_ARD_MC_LM_NAIVE(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE1156_LM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE1156_IV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE1156_LL = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
save(estimatePDE1156_LM_NAIVE,estimatePDE1156_LM,estimatePDE1156_IV,estimatePDE1156_LL,file="../datasets_montecarlo/resultPDE1156.RData")


load(file="../datasets_montecarlo/resultPDE1156.RData")

LM_est = estimatePDE1156_LM; s = summary(LM_est)
print(paste("xA in ",round((SARDp$gammaA-LM_est$coefficients["xA"])/s$coefficients["xA","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/abs(SARDp$gammaA),digits=2), " relative error.", sep=""))
print(paste("xR in ",round((SARDp$gammaR-LM_est$coefficients["xR"])/s$coefficients["xR","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/abs(SARDp$gammaR),digits=2), " relative error.", sep=""))
print(paste("xD in ",round((SARDp$gammaD-LM_est$coefficients["xD"])/s$coefficients["xD","Std. Error"],digits=2)," std. error, and ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/abs(SARDp$gammaD),digits=2), " relative error.", sep=""))

print(paste("xA in ",round((SARDp$gammaA- estimatePDE1156_LL$outARD_3MatEstimate$coef[1])/abs(SARDp$gammaA),digits=2), " relative error.", sep=""))
print(paste("xR in ",round((SARDp$gammaR- estimatePDE1156_LL$outARD_3MatEstimate$coef[2])/abs(SARDp$gammaR),digits=2), " relative error.", sep=""))
print(paste("xD in ",round((SARDp$gammaD- estimatePDE1156_LL$outARD_3MatEstimate$coef[3])/abs(SARDp$gammaD),digits=2), " relative error.", sep=""))



## plot shp ----
load("../datasets_montecarlo/DataShpRegressors144.RData")

dev.new()
plot(shp["y0"])
dev.copy2pdf(file="../datasets_montecarlo/y0144.pdf")

dev.new()
plot(shp["yT"])
dev.copy2pdf(file="../datasets_montecarlo/yT144.pdf")

dev.new()
plot(shp["delta"])
dev.copy2pdf(file="../datasets_montecarlo/delta144.pdf")

dev.new()
plot(cbind(shp,xA)["xA"])
dev.copy2pdf(file="../datasets_montecarlo/xA144.pdf")

dev.new()
plot(cbind(shp,xR)["xR"])
dev.copy2pdf(file="../datasets_montecarlo/xR144.pdf")

dev.new()
plot(cbind(shp,xD)["xD"])
dev.copy2pdf(file="../datasets_montecarlo/xD144.pdf")

## plot contour PDE ----
dev.new()
filled.contour(t(PDE0), color.palette=plasma)
dev.copy2pdf(file="../datasets_montecarlo/y01156_continuosSpace.pdf")
dev.new()
filled.contour(t(PDET), color.palette=plasma)
dev.copy2pdf(file="../datasets_montecarlo/yT1156_continuosSpace.pdf")
dev.new()
filled.contour(t(PDET)-t(PDE0), color.palette=plasma)
dev.copy2pdf(file="../datasets_montecarlo/delta_continuosSpace.pdf")

# montecarlo PDE different Np ----
load(file="../datasets_montecarlo/PDE.RData")

lengthAll = seq(from=12,to=50,by=1)
NpAll = lengthAll^2
Nm = 1000
coefMNpLM = array(data=NA,dim=c(length(NpAll),Nm,3))
coefMNpLMNaive = array(data=NA,dim=c(length(NpAll),Nm,3))
coefMNpIV = array(data=NA,dim=c(length(NpAll),Nm,3))
estimateNpLL = list()
coefNpLL = matrix(data=NA,nrow=length(NpAll),ncol=3)
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
    
    MADelta = as.numeric(MA %*% matrix(data$delta))
    MRDelta = as.numeric(MR %*% matrix(data$delta))
    MDDelta = as.numeric(MD %*% matrix(data$delta))
    
    X = data.frame(xA=xA,xR=xR,xD=xD)
    X = as.matrix(X)
    
    # instruments for IV
    MA2X=as.matrix(MA %*% MA %*% X)
    MR2X=as.matrix(MR %*% MR %*% X) 
    MD2X=as.matrix(MD %*% MD %*% X)
    
    estimateNpLL[[i]] = estimate_ARD_MC_LL(data,shp,xA,xR,xD,MA,MR,MD)
    coefNpLL[i,] = estimateNpLL[[i]]$outARD_3MatEstimate$coef[c(1,2,3)]
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
        
        estimatePDE_LMm = estimate_ARD_MC_LMBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
        coefMNpLM[i,m,] = estimatePDE_LMm$coefficients[c(1,3,5)]
        
        estimatePDE_LMNaivem = estimate_ARD_MC_LM_NAIVEBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
        coefMNpLMNaive[i,m,] = estimatePDE_LMNaivem$coefficients[c(1,2,3)]

        estimatePDE_IVm = estimate_ARD_MC_IVBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam,MA2Xm,MR2Xm,MD2Xm)
        coefMNpIV[i,m,] = coef(estimatePDE_IVm)[c(1,2,3)]
        
    }
}
save(NpAll,coefMNpLM,coefMNpLMNaive,coefMNpIV,coefNpLL,estimateNpLL,file="../datasets_montecarlo/montecarloPDENpGrid.RData")

## with SE ----
load(file="../datasets_montecarlo/montecarloPDENpGrid.RData")
library(latex2exp)
Np = NpAll
coefRep = matrix(rep(c(SARDp$gammaA,SARDp$gammaR,SARDp$gammaD),Nm),nrow=Nm,byrow = T)
coefNpLM = matrix(data=NA,nrow=length(Np),ncol=3)
lowerNpLM = matrix(data=NA,nrow=length(Np),ncol=3)
upperNpLM = matrix(data=NA,nrow=length(Np),ncol=3)
seNpLM = matrix(data=NA,nrow=length(Np),ncol=3)
coefNpLMNaive = matrix(data=NA,nrow=length(Np),ncol=3)
lowerNpLMNaive = matrix(data=NA,nrow=length(Np),ncol=3)
upperNpLMNaive = matrix(data=NA,nrow=length(Np),ncol=3)
seNpLMNaive = matrix(data=NA,nrow=length(Np),ncol=3)
coefNpIV = matrix(data=NA,nrow=length(Np),ncol=3)
lowerNpIV = matrix(data=NA,nrow=length(Np),ncol=3)
upperNpIV = matrix(data=NA,nrow=length(Np),ncol=3)
seNpIV = matrix(data=NA,nrow=length(Np),ncol=3)
coefNpLLRel = matrix(data=NA,nrow=length(Np),ncol=3)
coefNpLLSeRel = matrix(data=NA,nrow=length(Np),ncol=3)
    
for (i in 1:length(Np)){
    coefMLM = coefMNpLM[i,,]
    coefMLMNaive = coefMNpLMNaive[i,,]    
    coefMIV = coefMNpIV[i,,]   
    coefLL = coefNpLL[i,]
    SELL = estimateNpLL[[i]]$outARD_3MatEstimate$se_coef[c(1,2,3)]
    
    errRelLM = (coefRep-coefMLM)/coefRep
    errRelLMNaive = (coefRep-coefMLMNaive)/coefRep
    errRelIV = (coefRep-coefMIV)/coefRep
    coefNpLLRel[i,] =  (c(SARDp$gammaA,SARDp$gammaR,SARDp$gammaD)-coefLL)/c(SARDp$gammaA,SARDp$gammaR,SARDp$gammaD)
    coefNpLLSeRel[i,] = SELL/(abs(c(SARDp$gammaA,SARDp$gammaR,SARDp$gammaD)))
    
    coefNpLM[i,] = apply(errRelLM,MARGIN=2,FUN=function(x) mean(x))
    lowerNpLM[i,] = apply(errRelLM,MARGIN=2,FUN=function(x) quantile(x,0.05))
    upperNpLM[i,] = apply(errRelLM,MARGIN=2,FUN=function(x) quantile(x,0.95))
    seNpLM[i,] = apply(coefMLM,MARGIN=2,FUN=function(x) sd(x))
    
    coefNpLMNaive[i,] = apply(errRelLMNaive,MARGIN=2,FUN=function(x) mean(x))
    lowerNpLMNaive[i,] = apply(errRelLMNaive,MARGIN=2,FUN=function(x) quantile(x,0.05))
    upperNpLMNaive[i,] = apply(errRelLMNaive,MARGIN=2,FUN=function(x) quantile(x,0.95))
    seNpLMNaive[i,] = apply(coefMLMNaive,MARGIN=2,FUN=function(x) sd(x))
    
    coefNpIV[i,] = apply(errRelIV,MARGIN=2,FUN=function(x) mean(x))
    lowerNpIV[i,] = apply(errRelIV,MARGIN=2,FUN=function(x) quantile(x,0.05))
    upperNpIV[i,] = apply(errRelIV,MARGIN=2,FUN=function(x) quantile(x,0.95))
    seNpIV[i,] = apply(coefMIV,MARGIN=2,FUN=function(x) sd(x))
}

### plot ----

coefNpLM = coefNpLM[1:(length(Np)-1),]
lowerNpLM = lowerNpLM[1:(length(Np)-1),]
upperNpLM = upperNpLM[1:(length(Np)-1),]
coefNpLMNaive = coefNpLMNaive[1:(length(Np)-1),]
lowerNpLMNaive = lowerNpLMNaive[1:(length(Np)-1),]
upperNpLMNaive = upperNpLMNaive[1:(length(Np)-1),]
coefNpIV = coefNpIV[1:(length(Np)-1),]
lowerNpIV = lowerNpIV[1:(length(Np)-1),]
upperNpIV = upperNpIV[1:(length(Np)-1),]
coefNpLLRel = coefNpLLRel[1:(length(Np)-1),]
coefNpLLSeRel = coefNpLLSeRel[1:(length(Np)-1),]
Np = Np[1:(length(Np)-1)]

dev.new()
plot(Np,coefNpLM[,1],type="l",ylab = TeX("Relative error $\\gamma_A$"), xaxt = "n",ylim=range(coefNpLM[,1],lowerNpLM[,1],upperNpLM[,1],coefNpLMNaive[,1],lowerNpLMNaive[,1],upperNpLMNaive[,1],coefNpIV[,1],lowerNpIV[,1],upperNpIV[,1],coefNpLLRel[,1]))
lines(Np,lowerNpLM[,1],lty="dashed",col="black",xaxt="n")
lines(Np,upperNpLM[,1],lty="dashed",col="black",xaxt="n")
lines(Np,coefNpLMNaive[,1],col="red",xaxt="n")
lines(Np,upperNpLMNaive[,1],lty="dashed",col="red",xaxt="n")
lines(Np,lowerNpLMNaive[,1],lty="dashed",col="red",xaxt="n")
lines(Np,coefNpIV[,1],col="blue",xaxt="n")
lines(Np,upperNpIV[,1],lty="dashed",col="blue",xaxt="n")
lines(Np,lowerNpIV[,1],lty="dashed",col="blue",xaxt="n")
lines(Np,coefNpLLRel[,1],col="purple",xaxt="n")
lines(Np,coefNpLLRel[,1]+qnorm(0.95)*coefNpLLSeRel[,1],lty="dashed",col="purple",xaxt="n")
lines(Np,coefNpLLRel[,1]-qnorm(0.95)*coefNpLLSeRel[,1],lty="dashed",col="purple",xaxt="n")
axis(1, at = Np)
grid()
abline(h=0)
# legend(locator(1),legend=c("OLS","OLS NAIVE","Instrumental Variable","Maximum Likelihood"),col=c("black","red","blue","purple"),lty=1,cex=1.0)
legend("topleft",legend=c("OLS","OLS NAIVE","Instrumental Variable","Maximum Likelihood"),col=c("black","red","blue","purple"),lty=1,cex=1.0)
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRelGammaA.pdf")

dev.new()
plot(Np,coefNpLM[,2],type="l",ylab = TeX("Relative error $\\gamma_R$"), xaxt = "n",ylim=range(coefNpLM[,2],lowerNpLM[,2],upperNpLM[,2],coefNpLMNaive[,2],lowerNpLMNaive[,2],upperNpLMNaive[,2],coefNpIV[,2],lowerNpIV[,2],upperNpIV[,2],coefNpLLRel[,2]))
lines(Np,lowerNpLM[,2],lty="dashed",col="black",xaxt="n")
lines(Np,upperNpLM[,2],lty="dashed",col="black",xaxt="n")
lines(Np,coefNpLMNaive[,2],col="red",xaxt="n")
lines(Np,lowerNpLMNaive[,2],lty="dashed",col="red",xaxt="n")
lines(Np,upperNpLMNaive[,2],lty="dashed",col="red",xaxt="n")
lines(Np,coefNpIV[,2],col="blue",xaxt="n")
lines(Np,upperNpIV[,2],lty="dashed",col="blue",xaxt="n")
lines(Np,lowerNpIV[,2],lty="dashed",col="blue",xaxt="n")
lines(Np,coefNpLLRel[,2],col="purple",xaxt="n")
lines(Np,coefNpLLRel[,2]+qnorm(0.95)*coefNpLLSeRel[,2],lty="dashed",col="purple",xaxt="n")
lines(Np,coefNpLLRel[,2]-qnorm(0.95)*coefNpLLSeRel[,2],lty="dashed",col="purple",xaxt="n")
axis(1, at = Np)
grid()
abline(h=0)
# legend(locator(1),legend=c("LM","LM Naive","IV","LL"),col=c("black","red","blue","purple"),lty=1,cex=1.0)
legend("topleft",legend=c("OLS","OLS NAIVE","Instrumental Variable","Maximum Likelihood"),col=c("black","red","blue","purple"),lty=1,cex=1.0, ncol=4)
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRelGammaR.pdf")

dev.new()
plot(Np,coefNpLM[,3],type="l",ylab = TeX("Relative error $\\gamma_D$"), xaxt = "n",ylim=range(coefNpLM[,3],lowerNpLM[,3],upperNpLM[,3],coefNpLMNaive[,3],lowerNpLMNaive[,3],upperNpLMNaive[,3],coefNpIV[,3],lowerNpIV[,3],upperNpIV[,3],coefNpLLRel[,3]))
lines(Np,lowerNpLM[,3],lty="dashed",col="black",xaxt="n")
lines(Np,upperNpLM[,3],lty="dashed",col="black",xaxt="n")
lines(Np,coefNpLMNaive[,3],col="red",xaxt="n")
lines(Np,lowerNpLMNaive[,3],lty="dashed",col="red",xaxt="n")
lines(Np,upperNpLMNaive[,3],lty="dashed",col="red",xaxt="n")
lines(Np,coefNpIV[,3],col="blue",xaxt="n")
lines(Np,lowerNpIV[,3],lty="dashed",col="blue",xaxt="n")
lines(Np,upperNpIV[,3],lty="dashed",col="blue",xaxt="n")
lines(Np,coefNpLLRel[,3],col="purple",xaxt="n")
lines(Np,coefNpLLRel[,3]+qnorm(0.95)*coefNpLLSeRel[,3],lty="dashed",col="purple",xaxt="n")
lines(Np,coefNpLLRel[,3]-qnorm(0.95)*coefNpLLSeRel[,3],lty="dashed",col="purple",xaxt="n")
axis(1, at = Np)
grid()
abline(h=0)
# legend(locator(1),legend=c("LM","LM Naive","IV","LL"),col=c("black","red","blue","purple"),lty=1,cex=1.0)
legend("topleft",legend=c("OLS","OLS NAIVE","Instrumental Variable","Maximum Likelihood"),col=c("black","red","blue","purple"),lty=1,cex=1.0, ncol=4)
dev.copy2pdf(file="../datasets_montecarlo/montecarloPDEGridErrRelGammaD.pdf")



# montecarlo analysis ----
load("../datasets_montecarlo/resultMC225.RData")
load("../datasets_montecarlo/resultMC1156.RData")
load("../datasets_montecarlo/resultMC289.RData")
load("../datasets_montecarlo/resultPDE225.RData")
load("../datasets_montecarlo/resultPDE1156.RData")
load("../datasets_montecarlo/resultPDE289.RData")

## 2500 Result ----
load("../datasets_montecarlo/DataShpRegressors2500.RData")
estimatePDE2500_LM_NAIVE = estimate_ARD_MC_LM_NAIVE(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE2500_LM = estimate_ARD_MC_LM(data,shp,xA,xR,xD,MA,MR,MD)
estimatePDE2500_IV = estimate_ARD_MC_IV(data,shp,xA,xR,xD,MA,MR,MD)

table2500 = matrix(data = NA, nrow = 3, ncol = 2*4)
rownames(table2500) = c("gammaA","gammaR","gammaD")
colnames(table2500) = c("PDE_LMNaive","seMC_LMNaive",
                        "PDE_LM","seMC_LM",
                        "PDE_IV","seMC_IV",
                        "PDE_LL","seMC_LL")

table2500[,"PDE_LMNaive"] = coef(estimatePDE2500_LM_NAIVE)[c(1,2,3)]
table2500[,"seMC_LMNaive"] = seNpLMNaive[which(Np==2500),] 

table2500[,"PDE_LM"] = coef(estimatePDE2500_LM)[c(1,3,5)]
table2500[,"seMC_LM"] = seNpLM[which(Np==2500),] 

table2500[,"PDE_IV"] = coef(estimatePDE2500_IV)[c(1,2,3)]
table2500[,"seMC_IV"] = seNpIV[which(Np==2500),] 

table2500[,"PDE_LL"] = estimateNpLL[[which(Np==2500)]]$outARD_3MatEstimate$coef[c(1,2,3)]
table2500[,"seMC_LL"] = estimateNpLL[[which(Np==2500)]]$outARD_3MatEstimate$se_coef[c(1,2,3)]

# colnames(table2500) = sapply(colnames(table2500),FUN=function(x) return(paste(x,"2500",sep="")))

## output in xtable
table2500 = rbind(table2500[,1:4],table2500[,5:8],table2500[,9:12])

xtable(table225,auto=TRUE)
xtable(table576,auto=TRUE)
xtable(table289,auto=TRUE)
xtable(tableAll,auto=TRUE)
xtable(tableAll,display=rep("f",ncol(tableAll)+1))


## table All sorted by estimation method ----
tableLM225 = t(table225[,1:4])
tableLM576 = t(table576[,1:4])
tableLM289 = t(table289[,1:4])

tableLM_gammaA = cbind(tableLM225[,colnames(tableLM225)=="gammaA"],
                       tableLM576[,colnames(tableLM576)=="gammaA"],
                       tableLM289[,colnames(tableLM289)=="gammaA"])
colnames(tableLM_gammaA) = c("225", "576", "289")

tableLM_gammaR = cbind(tableLM225[,colnames(tableLM225)=="gammaR"],
                       tableLM576[,colnames(tableLM576)=="gammaR"],
                       tableLM289[,colnames(tableLM289)=="gammaR"])
colnames(tableLM_gammaR) = c("225", "576", "289")

tableLM_gammaD = cbind(tableLM225[,colnames(tableLM225)=="gammaD"],
                       tableLM576[,colnames(tableLM576)=="gammaD"],
                       tableLM289[,colnames(tableLM289)=="gammaD"])
colnames(tableLM_gammaD) = c("225", "576", "289")

table_LM = cbind(tableLM_gammaA, tableLM_gammaR, tableLM_gammaD)
xtable(table_LM, digits=5)


tableIV225 = t(table225[,5:8])
tableIV576 = t(table576[,5:8])
tableIV289 = t(table289[,5:8])

tableIV_gammaA = cbind(tableIV225[,colnames(tableIV225)=="gammaA"],
                       tableIV576[,colnames(tableIV576)=="gammaA"],
                       tableIV289[,colnames(tableIV289)=="gammaA"])
colnames(tableIV_gammaA) = c("225", "576", "289")

tableIV_gammaR = cbind(tableIV225[,colnames(tableIV225)=="gammaR"],
                       tableIV576[,colnames(tableIV576)=="gammaR"],
                       tableIV289[,colnames(tableIV289)=="gammaR"])
colnames(tableIV_gammaR) = c("225", "576", "289")

tableIV_gammaD = cbind(tableIV225[,colnames(tableIV225)=="gammaD"],
                       tableIV576[,colnames(tableIV576)=="gammaD"],
                       tableIV289[,colnames(tableIV289)=="gammaD"])
colnames(tableIV_gammaD) = c("225", "576", "289")

table_IV = cbind(tableIV_gammaA, tableIV_gammaR, tableIV_gammaD)
xtable(table_IV, digits=5)

tableLL225 = t(table225[,9:12])
tableLL576 = t(table576[,9:12])
tableLL289 = t(table289[,9:12])

tableLL_gammaA = cbind(tableLL225[,colnames(tableLL225)=="gammaA"],
                       tableLL576[,colnames(tableLL576)=="gammaA"],
                       tableLL289[,colnames(tableLL289)=="gammaA"])
colnames(tableLL_gammaA) = c("225", "576", "289")

tableLL_gammaR = cbind(tableLL225[,colnames(tableLL225)=="gammaR"],
                       tableLL576[,colnames(tableLL576)=="gammaR"],
                       tableLL289[,colnames(tableLL289)=="gammaR"])
colnames(tableLL_gammaR) = c("225", "576", "289")

tableLL_gammaD = cbind(tableLL225[,colnames(tableLL225)=="gammaD"],
                       tableLL576[,colnames(tableLL576)=="gammaD"],
                       tableLL289[,colnames(tableLL289)=="gammaD"])
colnames(tableLL_gammaD) = c("225", "576", "289")

table_LL = cbind(tableLL_gammaA, tableLL_gammaR, tableLL_gammaD)
xtable(table_LL, digits=5)


# ---- MC finto with bootstrap
load("../datasets_montecarlo/DataShpRegressors576.RData")
MADelta = as.numeric(MA %*% matrix(data$delta))
MRDelta = as.numeric(MR %*% matrix(data$delta))
MDDelta = as.numeric(MD %*% matrix(data$delta))

Nm = 1000
coefM = matrix(data=NA,nrow=Nm,ncol=3)
for (m in 1:Nm){
    print(m)
    
    iBoot = sample(1:nrow(datam),replace=T)
    datam = data[iBoot,]
    xAm = xA[iBoot]
    xRm = xR[iBoot]
    xDm = xD[iBoot]
    MADeltam = MADelta[iBoot]
    MRDeltam = MRDelta[iBoot]
    MDDeltam = MDDelta[iBoot]
    
    estimatePDE576_LMm = estimate_ARD_MC_LMBoot(datam,shp,xAm,xRm,xDm,MADeltam,MRDeltam,MDDeltam)
    coefM[m,] = estimatePDE576_LMm$coefficients[c(1,3,5)]
}

dev.new()
par(mfrow=c(1,3))
# xA
hist(coefM[,1],50,xlim=range(c(coefM[,1],SARDp$gammaA)),main="gammaA")
abline(v=SARDp$gammaA,col="red")
abline(v=estimatePDE576_LM$coefficients[1],col="blue")

# xR
hist(coefM[,2],50,main="gammaR")
abline(v=SARDp$gammaR,col="red")
abline(v=estimatePDE576_LM$coefficients[3],col="blue")

# xD
hist(coefM[,3],50,xlim=range(SARDp$gammaD,coefM[,3]),main="gammaD")
abline(v=SARDp$gammaD,col="red")
abline(v=estimatePDE576_LM$coefficients[5],col="blue")
legend("topright",legend=c("true parameter","estimated via PDE"),col=c("red","blue"),lty=1,cex=1.0)

dev.copy2pdf(file="bootCoeff.pdf")

# Analysis of Wepsilon for different Np ----

lambdaHat = vector("numeric",length(Np))
lambdas = matrix(NA, nrow=length(Np), ncol=10)
maxLag = vector("numeric",length(Np))
for (i in 1:length(Np)){
    lambdaHat[i] = round(estimateNpLL[[i]]$outARD_3MatEstimate$coef[length(estimateNpLL[[i]]$outARD_3MatEstimate$coef)], digits=3)
    lambdas[i,] = round(coef(estimateNpLL[[i]]$SpatError$lm.errDecompose), digits=3)
    maxLag[i] = estimateNpLL[[i]]$SpatError$maxSignifLag
}

tableW = cbind(round(Np, digits=0), lambdaHat, lambdas, maxLag)
colnames(tableW) = c("Np", "lambda hat", paste("lambdaHat",seq(1:10)), "max lag")

print(xtable(tableW), include.rownames=FALSE)
