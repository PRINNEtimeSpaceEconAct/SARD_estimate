rm(list = ls())

source("lib/L_loadAll.R")
DEBUG = TRUE
PARALLEL = TRUE
NPROCS = 30

initJulia()


Np=1000
Na=10000
Nm = 2
tau=1.0
SARDp = list(gammaS = -0.05, gammaA = -0.2, gammaR = 0.1, gammaD = 0.001, hA = 0.1, hR = 0.3)

# hA_range=c(0.2,0.5)
# hR_range=c(0.3,0.7)
# resultMonteCarlo = MonteCarlo(Np, Na, Nm, tau, typeOfDist, SARDp, hA_range, hR_range)

resultMonteCarlo_FixedhAhR_IV = MonteCarloFixedhAhR_IV(Np, Na, Nm, tau, typeOfDist,SARDp,SARDp$hA,SARDp$hR)
    