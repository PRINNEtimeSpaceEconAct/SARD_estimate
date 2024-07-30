# libraries ----

suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(spdep))
suppressPackageStartupMessages(library(spatialreg))
suppressPackageStartupMessages(library(gamlr))
suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(sfheaders))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(snow))
suppressPackageStartupMessages(library(plm))
suppressPackageStartupMessages(library(AER))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(spatstat.geom))
suppressPackageStartupMessages(library(buildmer))
suppressPackageStartupMessages(library(JuliaCall))
suppressPackageStartupMessages(library(cvCovEst))
suppressPackageStartupMessages(library(classInt))
suppressPackageStartupMessages(library(tmap))
suppressPackageStartupMessages(library(polyCub))
suppressPackageStartupMessages(library(RColorBrewer))


# sources ----
source("lib/L_data.R")
source("lib/L_spatialTools.R")
source("lib/L_Regressors.R")
source("lib/L_callJulia.R")
source("lib/L_IV.R")
source("lib/L_ML.R")
source("lib/L_spatialModels.R")
source("lib/L_GFDM.R")
source("lib/L_maps.R")
source("lib/L_forecast.R")
source("lib/L_forecastSpatial.R")
source("lib/L_MonteCarloGeneration.R")
source("lib/L_MonteCarloEstimation.R")


# globals
Sys.setenv(JULIA_NUM_THREADS = "8")
DEBUG = FALSE
PARALLEL = FALSE
NPROCS = 1

# install.packages(c("MASS","dplyr","spdep","spatialreg","gamlr","sf","sfheaders","Matrix","snow","plm","AER","pracma","spatstat.geom","buildmer","JuliaCall","cvCovEst","classInt","tmap"))

