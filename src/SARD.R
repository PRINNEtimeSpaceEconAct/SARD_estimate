rm(list = ls())

source("lib/L_loadAll.R")
DEBUG = TRUE
PARALLEL = TRUE
NPROCS = 8
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
# hR_range = seq(from=25,to=50,by=5)
# bestIV_est = chose_hAhR(data,hA_range,hR_range)
# hA = bestIV_est$hABest
# hR = bestIV_est$hRBest

hA = 10
hR = 45

# outIVEstimate = estimate_IV_SARD_auto(data,hA,hR)
# outSARDEstimate = estimate_SARD_auto(data,shp_data,hA,hR)

# load results ----
load("../datasets_it/outSARDEstimate.RData")
load("../datasets_it/outDURBINEstimate.RData")

## correlogramm SARD with spatial error ----
dev.new()
correlogramSARD = correlogram(outSARDEstimate$residualsSARD,shp_data,maxLag = 10)
dev.copy2pdf(file="../datasets_montecarlo/correlogramSARD_eta_I.pdf")

## map errors ----
mapErrorsSameScale(outSARDEstimate,DURBINEstimate,data,shp_data)

## spatial error matrix decompose ----
xtable(summary(outSARDEstimate$SpatError$lm.errDecompose))

## counterfactuals ----
names(outSARDEstimate$coefSARD) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD","lambda")

### SARD ----
SARDCoeff = outSARDEstimate$coefSARD
yFutSARD = forcastSARD(50,SARDCoeff,hA,hR,data,tau=11)

### ARD ----
ARDCoeff = outSARDEstimate$coefSARD
ARDCoeff["gammaS"] = 0
ARDCoeff["rhoS"] = 0
yFutARD = forcastSARD(50,ARDCoeff,hA,hR,data,tau=11)

### SRD ----
SRDCoeff = outSARDEstimate$coefSARD
SRDCoeff["gammaA"] = 0
SRDCoeff["rhoA"] = 0
yFutSRD = forcastSARD(50,SRDCoeff,hA,hR,data,tau=11)

### SRD ----
SADCoeff = outSARDEstimate$coefSARD
SADCoeff["gammaR"] = 0
SADCoeff["rhoR"] = 0
yFutSAD = forcastSARD(50,SADCoeff,hA,hR,data,tau=11)

### SAR ----
SARCoeff = outSARDEstimate$coefSARD
SARCoeff["gammaR"] = 0
SARCoeff["rhoR"] = 0
yFutSAR = forcastSAR(50,SARCoeff,hA,hR,data,tau=11)

# save
dfForecast = data.frame(yFutSARD=yFutSARD,yFutARD=yFutARD,yFutSRD=yFutSRD,yFutSAD=yFutSAD,yFutSAR=yFutSAR) 
save(dfForecast,file="../datasets_it/dfForecast.RData")


# df <- df %>% mutate(cities=ifelse(geo=="058091" | geo=="015146" | geo=="063049" |geo=="001272" | geo=="072006" | geo=="082053" | geo=="087015" | geo=="037006" | geo=="048017" | geo=="027042" | geo=="010025" | geo=="083048" | geo=="080063" | geo=="092009", 1, 0))

