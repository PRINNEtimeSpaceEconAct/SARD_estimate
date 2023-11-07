rm(list = ls())

source("lib/L_loadAll.R")
DEBUG = TRUE
PARALLEL = TRUE
NPROCS = 30

folderName = "../datasets_it/"

initJulia()

rawData = create_df(folderName)
data = rawData$data
shp_data = rawData$shp_data

# old coordinates
# load(file = "df.tot.RData")
# data$Latitude = df.tot$Latitude.old
# data$Longitude = df.tot$Longitude.old

# hA_range = seq(from=5,to=20,by=5)
# hR_range = seq(from=20,to=50,by=5)
# bestIV_est = chose_hAhR(data,hA_range,hR_range)
# hA = bestIV_est$hABest
# hR = bestIV_est$hRBest

hA = 10
hR = 30

# outIVEstimate = estimate_IV_SARD_auto(data,hA,hR)
# outWNEstimate = estimate_WN_SARD_auto(data,hA,hR)
outWNEstimate = estimate_SARD_auto(data,shp_data,hA,hR)

load("../outWNEstimate.RData")





