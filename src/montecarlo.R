rm(list = ls())

source("lib/L_loadAll.R")
DEBUG = TRUE
PARALLEL = TRUE
NPROCS = 8
torus = TRUE
initJulia()
set.seed(1)


# parameters DGP ----
Np=225
Na=10000
Nm = 100
tau=0.1
NeS = 100
SARDp = list(gammaS = 0.0, gammaA = -0.2, gammaR = 0.0, gammaD = 0.0, hA = 0.3, hR = 0.1)
# SARDp = list(gammaS = 0.1, gammaA = -0.0, gammaR = 0.0, gammaD = 0.0, hA = 0.3, hR = 0.0)

# parameters estimation ----
model = list(regressors=as.formula(delta ~ y0 + xS + xA + xR + xD + MSDelta + MADelta + MRDelta + MDDelta),
             instruments=as.formula(~ y0 + xS + xA + xR + xD + MS2X + MA2X + MR2X + MD2X))
variables = c("ones","y0","xS","xA","xR","xD")

# montecarlo one run with fixed hA hR ----
set.seed(1)
resultMonteCarloOneRunIV_fixedhAhR = MonteCarloOneRunIV_fixedhAhR(Np, Na, tau, typeOfDist = "uniform",SARDp,SARDp$hA,SARDp$hR,model=model,variables=variables)
summary(resultMonteCarloOneRunIV_fixedhAhR$outIVEstimate$IV_est)
