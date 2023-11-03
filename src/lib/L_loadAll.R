# libraries ----

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(spdep))
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
source("lib/callJulia.R")
source("IV.R")
source("ML.R")

# globals
DEBUG = FALSE

