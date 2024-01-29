rm(list = ls())
source("lib/L_loadAll.R")
DEBUG = FALSE
PARALLEL = TRUE
NPROCS = 12
initJulia()


# parameters DGP ----
Np=250
Na=100000
Nm = 100
tau= 0.05
NeS = 100
typeOfDist = "Uniform"
typeOfEst = "LM"
# typeOfEst = "WNLL"
# typeOfEst = "SELL"


# funzionano A e D
# Np=250
# tau= 0.1
# SARDp = list(gammaS = 0.0, gammaA = -0.05, gammaR = 0.0, gammaD = 0.03, hA = 0.3, hR = 0.4)

# funzionano A R e D
# Np=250
# tau= 0.05
# SARDp = list(gammaS = 0.0, gammaA = -0.03, gammaR = 0.1, gammaD = 0.09, hA = 0.15, hR = 0.4)

# parameters estimation PDE ----
SARDp = list(gammaS = 0.0, gammaA = -0.03, gammaR = 0.1, gammaD = 0.09, hA = 0.15, hR = 0.4)
resultMonteCarloOneRun_LM_PDE_fixedhAhRUniform = MonteCarloOneRun_LM_PDE_fixedhAhR(Np, tau, typeOfDist = "Uniform",typeOfEst=typeOfEst,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
LM_est = resultMonteCarloOneRun_LM_PDE_fixedhAhRUniform$outEstimate$LM_est; s = summary(LM_est)
summary(LM_est)
print(paste("xA in ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/s$coefficients["xA","Std. Error"],digits=2)," std. error",sep=""))
print(paste("xR in ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/s$coefficients["xR","Std. Error"],digits=2)," std. error",sep=""))
print(paste("xD in ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/s$coefficients["xD","Std. Error"],digits=2)," std. error",sep=""))

# parameters estimation Agents ----
SARDp = list(gammaS = 0.0, gammaA = -0.03, gammaR = 0.1, gammaD = 0.09, hA = 0.15, hR = 0.4)
resultMonteCarloOneRun_LM_Agents_fixedhAhRUniform = MonteCarloOneRun_LM_Agents_fixedhAhR(Np,Na,tau,typeOfDist = "Uniform",typeOfEst=typeOfEst,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
LM_est = resultMonteCarloOneRun_LM_Agents_fixedhAhRUniform$outEstimate$LM_est; s = summary(LM_est)
summary(LM_est)
print(paste("xA in ",round(abs(SARDp$gammaA-LM_est$coefficients["xA"])/s$coefficients["xA","Std. Error"],digits=2)," std. error",sep=""))
print(paste("xR in ",round(abs(SARDp$gammaR-LM_est$coefficients["xR"])/s$coefficients["xR","Std. Error"],digits=2)," std. error",sep=""))
print(paste("xD in ",round(abs(SARDp$gammaD-LM_est$coefficients["xD"])/s$coefficients["xD","Std. Error"],digits=2)," std. error",sep=""))



# montecarlo PDE one run with fixed hA hR ----
# resultMonteCarloOneRun_LM_PDE_fixedhAhRUniform = MonteCarloOneRun_LM_PDE_fixedhAhR(Np, tau, typeOfDist = "Uniform",typeOfEst=typeOfEst,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
# set.seed(1); resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi = MonteCarloOneRun_LM_PDE_fixedhAhR(Np, tau, typeOfDist = "VoronoiUniform",typeOfEst=typeOfEst,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
# set.seed(1); resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi2Uniform = MonteCarloOneRun_LM_PDE_fixedhAhR_Voronoi2Unif(Np, tau,typeOfEst=typeOfEst,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
# summary(resultMonteCarloOneRun_LM_PDE_fixedhAhRUniform$outEstimate$LM_est)
# summary(resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$outEstimate$LM_est)
# summary(resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi2Uniform$outEstimate$LM_est)
# 
# dataWithRegressors = cbind(resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$data[,-c(4)]
#                            , resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$outEstimate$LM_est$model)  
# 
# dataWithRegressors = dataWithRegressors %>% mutate(xA=xA-sum(xA*km2),
#                                                    MADelta=MADelta-sum(MADelta*km2))
# 
# scaledLMVoronoi = lm(delta ~ xA + MADelta, data=dataWithRegressors)
# summary(scaledLMVoronoi)





# weigthedLMVoronoi = lm(delta ~ 0 + xA + MADelta, weights = km2, data=dataWithRegressors)
# summary(weigthedLMVoronoi)
# 
# sqrtW = diag(sqrt(dataWithRegressors$km2))
# dataWeighted = dataWithRegressors
# dataWeighted$delta = sqrtW %*% dataWeighted$delta
# dataWeighted$xA = sqrtW %*% sqrtW %*% dataWeighted$xA
# dataWeighted$MADelta = sqrtW %*% sqrtW %*% dataWeighted$MADelta
# weigthedLMVoronoi = lm(delta ~ 0 + xA + MADelta, data=dataWeighted)
# summary(weigthedLMVoronoi)

# deltaWei=diag(sqrt(dataWithRegressors$km2))%*%as.matrix(dataWithRegressors$delta)
# xAWei=diag(sqrt(dataWithRegressors$km2))%*%as.matrix(dataWithRegressors$xA)
# MADeltaWei=diag(sqrt(dataWithRegressors$km2))%*%as.matrix(dataWithRegressors$MADelta)
# 
# weigthedLMVoronoi = lm(deltaWei ~ 0 + xAWei + MADeltaWei)
# summary(weigthedLMVoronoi)



# montecarlo Agents one run with fixed hA hR ----
# set.seed(1)
# resultMonteCarloOneRun_LM_Agents_fixedhAhR = MonteCarloOneRun_LM_Agents_fixedhAhR(Np, Na, tau, typeOfDist = typeOfDist,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
# summary(resultMonteCarloOneRun_LM_Agents_fixedhAhR$outLMEstimate$LM_est)

# plot shp ----
# shp_regressors = resultMonteCarloOneRun_LM_Agents_fixedhAhR$outLMEstimate$shp_regressors
regression = resultMonteCarloOneRun_LM_PDE_fixedhAhRUniform
regressionLM = regression$outEstimate$LM_est
shp_regressors = regression$outEstimate$shp_regressors
plot(cbind(shp_regressors,resid = scale((resid(regressionLM))))["resid"])


plot(shp_regressors[c("y0","yT")],key.pos=4)
plot(shp_regressors[c("y0","WhAy0")],key.pos=4)
plot(shp_regressors[c("y0","WhRy0")],key.pos=4)
plot(shp_regressors[c("delta")],key.pos=4)
plot(shp_regressors[c("xA")],key.pos=4)
plot(shp_regressors[c("MADelta")],key.pos=4)
plot(mutate(shp_regressors,xA=xA)[c("xA","delta")],key.pos=4)
plot(mutate(shp_regressors,ratio = xA/delta)["ratio"],key.pos=4)
plot(mutate(shp_regressors,diff = -0.1*xA-delta)["diff"],key.pos=4)
plot(shp_regressors["xD"],key.pos=4)

testShp = mutate(shp_regressors, "Mxx" = as.matrix(resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$outEstimate$MsDeriv$Mxx) %*% y0, "MxMx" = as.matrix(resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$outLMEstimate$MsDeriv$Mx) %*%  as.matrix(resultMonteCarloOneRun_LM_PDE_fixedhAhR$outLMEstimate$MsDeriv$Mx) %*% y0)
plot(testShp[c("Mxx","MxMx")])

library(plotly)
plotVar <- function(estimate,name){
    df = data.frame(x = estimate$shp$Longitude, y = estimate$shp$Latitude, z = estimate$outLMEstimate$shp_regressors[[eval(name)]])
    plot_ly(df,x = ~x,y = ~y, z = ~z, color = I("red"), type = "mesh3d",size = 0.1, showlegend = FALSE)
}
df = data.frame(x = resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi2Uniform$shp$Longitude, y = resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi2Uniform$shp$Latitude, z = resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi2Uniform$outLMEstimate$shp_regressors$y0)
plot_ly(df,x = ~x,y = ~y, z = ~z, color = I("red"), type = "mesh3d",size = 0.1, showlegend = FALSE)


# confronto voronoi uniforme ----
resultMonteCarloOneRun_LM_PDE_fixedhAhRUniforme = MonteCarloOneRun_LM_PDE_fixedhAhR(Np, tau, typeOfDist = "Uniform",SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi = MonteCarloOneRun_LM_PDE_fixedhAhR(Np, tau, typeOfDist = "VoronoiUniform",SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
cor(resultMonteCarloOneRun_LM_PDE_fixedhAhRUniforme$outLMEstimate$LM_est$model$delta,resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$outLMEstimate$LM_est$model$delta)
cor(resultMonteCarloOneRun_LM_PDE_fixedhAhRUniforme$outLMEstimate$LM_est$model$MADelta,resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$outLMEstimate$LM_est$model$MADelta)

# spatial test on residuals
library(spdep)
# first-order contiguity matrix
shpVoronoi=resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$shp
contig = poly2nb(shpVoronoi)
lm.morantest(resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$outLMEstimate$LM_est, listw=nb2listw(contig, style="W"))

corSpat=sp.correlogram(contig, resultMonteCarloOneRun_LM_PDE_fixedhAhRVoronoi$outLMEstimate$LM_est$residuals, order=5, method = "I")
plot(corSpat)

# pseudo R2 ----
df = resultMonteCarloOneRun_LM_Agents_fixedhAhR$outLMEstimate$shp_regressors
RSSIV = mean((resultMonteCarloOneRun_LM_Agents_fixedhAhR$outLMEstimate$LM_est$residuals)^2)
delta = df$delta.x
RSSNULL = mean((delta-mean(delta))^2)
PR2 = 1-RSSIV/RSSNULL
print(PR2)


# cross-section ----
library(mgcv)
library(sm)
# df = filter(df,df$xA < 0)
lmModel = lm(delta.x ~ 0 + xD + MDDelta ,data = df)
smreg = sm.regression(df$xD,df$delta.x)
# lines(smreg$eval.points,coef(lmModel)[1]*smreg$eval.points,col="red")
# lines(smreg$eval.points,SARDp$gammaA*smreg$eval.points,col="blue")
# abline(v = 0)
gamreg = gam(df$delta.x ~ s(df$xD)+s(df$MDDelta))
plot(gamreg,select=1)
points(df$xA,df$delta.x,pch=19,cex=0.2)
lines(smreg$eval.points,coef(lmModel)[1]*smreg$eval.points,col="red")
lines(smreg$eval.points,SARDp$gammaD*smreg$eval.points,col="blue")
abline(v = 0)
abline(h=0)


