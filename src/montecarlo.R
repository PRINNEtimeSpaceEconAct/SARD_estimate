rm(list = ls())

source("lib/L_loadAll.R")
DEBUG = TRUE
PARALLEL = TRUE
NPROCS = 8
initJulia()
set.seed(1)


# parameters DGP ----
Np=500
Na=10000
Nm = 100
tau=1.0
NeS = 100
SARDp = list(gammaS = 0.0, gammaA = -0.2, gammaR = 0.0, gammaD = 0.0, hA = 0.3, hR = 0.0)
# SARDp = list(gammaS = 0.0, gammaA = 0.0, gammaR = 0.0, gammaD = 0.01, hA = 0.3, hR = 0.4)

# parameters estimation ----
model = list(regressors=as.formula(delta ~ 0 + xA + MADelta),instruments=as.formula(~ 0 + xA + MA2X))
variables = c("xA")

# model = list(regressors=as.formula(delta ~ y0 + xS + xA + xR + xD + MSDelta + MADelta + MRDelta + MDDelta),
#              instruments=as.formula(~ y0 + xS + xA + xR + xD + MS2X + MA2X + MR2X + MD2X))
# variables = c("ones","y0","xS","xA","xR","xD")

# montecarlo one run with fixed hA hR ----
set.seed(1)
resultMonteCarloOneRunIV_fixedhAhR = MonteCarloOneRunIV_fixedhAhR(Np, Na, tau, typeOfDist = "uniform",SARDp,SARDp$hA,SARDp$hR,model=model,variables=variables)
summary(resultMonteCarloOneRunIV_fixedhAhR$outIVEstimate$IV_est)

# full montecarlo with range search ----
# hA_range=c(0.2,0.5)
# hR_range=c(0.3,0.7)
# resultMonteCarlo = MonteCarlo(Np, Na, Nm, tau, typeOfDist, SARDp, hA_range, hR_range)

# montecarlo with fixed hA hR using IV ----
coef = list()
std_coef = list()
set.seed(1)
for  (i in 1:Nm) {
    resultMonteCarloOneRunIV_fixedhAhR = MonteCarloOneRunIV_fixedhAhR(Np, Na, tau, typeOfDist = "uniform",SARDp,SARDp$hA,SARDp$hR,model=model,variables=variables)
    coef[[i]]  = resultMonteCarloOneRunIV_fixedhAhR$outIVEstimate$IV_est$coefficients
    s = summary(resultMonteCarloOneRunIV_fixedhAhR$outIVEstimate$IV_est)
    std_coef[[i]] = s$coefficients[,'Std. Error'] } 



