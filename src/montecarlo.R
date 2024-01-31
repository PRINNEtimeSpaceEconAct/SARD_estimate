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


# Np 256 ----
Np=256
resultMonteCarloOneRun_LM_PDE_fixedhAhRUniform = MonteCarloOneRun_LM_PDE_fixedhAhR(Np, tau, typeOfDist = "Uniform",typeOfEst=typeOfEst,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)

resultMonteCarloOneRun_LM_Agents_fixedhAhRUniform = MonteCarloOneRun_LM_Agents_fixedhAhR(Np,Na,tau,typeOfDist = "Uniform",typeOfEst=typeOfEst,SARDp,SARDp$hA,SARDp$hR,torus=TRUE)
