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

# PDE once
# PDEAll = call_julia_computePDE(tau,SARDp)
# save(PDEAll,file="PDE.RData")

# Agents for error
NSigma = 100

# Na = 50000
# AgentsAll = call_julia_computeAgents(NSigma,Na,tau,SARDp)
# save(AgentsAll,file="AgentsAll50k.RData")

# Na = 100000
# AgentsAll = call_julia_computeAgents(NSigma,Na,tau,SARDp)
# save(AgentsAll,file="AgentsAll100k.RData")

# Na = 200000
# AgentsAll = call_julia_computeAgents(NSigma,Na,tau,SARDp)
# save(AgentsAll,file="AgentsAll200k.RData")

