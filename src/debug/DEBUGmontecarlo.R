rm(list = ls())
source("lib/L_loadAll.R")
DEBUG = FALSE
PARALLEL = TRUE
NPROCS = 12
initJulia()


# parameters DGP ----
Np=200
Na=100000
Nm = 100
tau=0.1
NeS = 100
typeOfDist = "Uniform"
# typeOfDist = "VoronoiUniform"

# parameters estimation ----
# SARDp = list(gammaS = 0.0, gammaA = 0.0, gammaR = 0.0, gammaD = 0.1, hA = 0.15, hR = 0.2)
SARDp = list(gammaS = 0.0, gammaA = -0.1, gammaR = 0.0, gammaD = 0.0, hA = 0.4, hR = 0.4)

# montecarlo PDE one run with fixed hA hR ----
set.seed(1)
resultMonteCarloOneRun_LM_PDE_fixedhAhR = MonteCarloOneRun_LM_PDE_fixedhAhR(Np, tau, typeOfDist = typeOfDist,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
summary(resultMonteCarloOneRun_LM_PDE_fixedhAhR$outLMEstimate$LM_est)
 
# modelloLM = resultMonteCarloOneRun_LM_PDE_fixedhAhR$outLMEstimate$LM_est$model
# summary(lm(modelloLM$delta ~ 0 + modelloLM$xA + modelloLM$MADelta))

# montecarlo Agents one run with fixed hA hR ----
set.seed(1)
resultMonteCarloOneRun_LM_Agents_fixedhAhR = MonteCarloOneRun_LM_Agents_fixedhAhR(Np, Na, tau, typeOfDist = typeOfDist,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
summary(resultMonteCarloOneRun_LM_Agents_fixedhAhR$outLMEstimate$LM_est)

# plot shp ----
# shp_regressors = resultMonteCarloOneRun_LM_Agents_fixedhAhR$outLMEstimate$shp_regressors
shp_regressors = resultMonteCarloOneRun_LM_PDE_fixedhAhR$outLMEstimate$shp_regressors
plot(shp_regressors[c("y0","yT")],key.pos=4)
plot(shp_regressors[c("y0","WhAy0")],key.pos=4)
plot(shp_regressors[c("y0","WhRy0")],key.pos=4)
plot(mutate(shp_regressors,xD = 0.1*xD)[c("xD","delta")],key.pos=4)
plot(mutate(shp_regressors,ratio = 0.1*xD/delta)["ratio"],key.pos=4)
plot(shp_regressors["xD"],key.pos=4)

testShp = mutate(shp_regressors, "Mxx" = as.matrix(resultMonteCarloOneRun_LM_PDE_fixedhAhR$outLMEstimate$MsDeriv$Mxx) %*% y0, "MxMx" = as.matrix(resultMonteCarloOneRun_LM_PDE_fixedhAhR$outLMEstimate$MsDeriv$Mx) %*%  as.matrix(resultMonteCarloOneRun_LM_PDE_fixedhAhR$outLMEstimate$MsDeriv$Mx) %*% y0)
plot(testShp[c("Mxx","MxMx")])

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

