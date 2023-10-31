# libraries ----

library(dplyr)
library(spdep)
library(sf)
library(Matrix)
library(snow)
library(plm)
library(ivreg)

# sources ----
source("lib/L_data.R")
source("lib/L_spatialTools.R")
source("lib/L_Regressors.R")
source("IV.R")
source("SARD.R")

# DEBUG
DEBUG = FALSE

