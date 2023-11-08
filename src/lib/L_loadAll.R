# libraries ----

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(spdep))
suppressPackageStartupMessages(library(spatialreg))
suppressPackageStartupMessages(library(gamlr))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(snow))
suppressPackageStartupMessages(library(plm))
suppressPackageStartupMessages(library(AER))
suppressPackageStartupMessages(library(JuliaCall))

# sources ----
source("lib/L_data.R")
source("lib/L_spatialTools.R")
source("lib/L_Regressors.R")
source("lib/L_callJulia.R")
source("lib/L_IV.R")
source("lib/L_ML.R")
source("lib/L_spatialModels.R")

# globals
DEBUG = FALSE
PARALLEL = FALSE
NPROCS = 1

