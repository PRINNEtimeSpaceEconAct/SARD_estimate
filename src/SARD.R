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
shpAll = get_data(folderName)$shp
mapErrorsSameScale(outSARDEstimate,DURBINEstimate,data,shpAll)

## spatial error matrix decompose ----
xtable(summary(outSARDEstimate$SpatError$lm.errDecompose))

## counterfactuals ----
names(outSARDEstimate$coefSARD) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD","lambda")

### SARD ----
SARDCoeff = outSARDEstimate$coefSARD
yFutSARD = forcastSARD(50,SARDCoeff,hA,hR,data,tau=11)
save(yFutSARD,file="../datasets_it/yFutSARD.RData")

### ARD ----
ARDCoeff = outSARDEstimate$coefSARD
ARDCoeff["gammaS"] = 0
ARDCoeff["rhoS"] = 0
yFutARD = forcastSARD(50,ARDCoeff,hA,hR,data,tau=11)
save(yFutARD,file="../datasets_it/yFutARD.RData")

### SRD ----
SRDCoeff = outSARDEstimate$coefSARD
SRDCoeff["gammaA"] = 0
SRDCoeff["rhoA"] = 0
yFutSRD = forcastSARD(50,SRDCoeff,hA,hR,data,tau=11)
save(yFutSRD,file="../datasets_it/yFutSRD.RData")

### SAD ----
SADCoeff = outSARDEstimate$coefSARD
SADCoeff["gammaR"] = 0
SADCoeff["rhoR"] = 0
yFutSAD = forcastSARD(50,SADCoeff,hA,hR,data,tau=11)
save(yFutSAD,file="../datasets_it/yFutSAD.RData")

### SAR ----
SARCoeff = outSARDEstimate$coefSARD
SARCoeff["gammaR"] = 0
SARCoeff["rhoR"] = 0
yFutSAR = forcastSARD(50,SARCoeff,hA,hR,data,tau=11)
save(yFutSAR,file="../datasets_it/yFutSAR.RData")

# save
load(file="../datasets_it/yFutSARD.RData")
load(file="../datasets_it/yFutARD.RData")
load(file="../datasets_it/yFutSRD.RData")
load(file="../datasets_it/yFutSAD.RData")
load(file="../datasets_it/yFutSAR.RData")
dfForecast = data.frame(yFutSARD=yFutSARD,yFutARD=yFutARD,yFutSRD=yFutSRD,yFutSAD=yFutSAD,yFutSAR=yFutSAR) 
save(dfForecast,file="../datasets_it/dfForecast.RData")

# maps
load("../datasets_it/dfForecast.RData")
shpAll = get_data(folderName)$shp
mapYFinalForecast(dfForecast,data,shpAll,"../datasets_it/map_incomeKm2_2069SARD.pdf")
mapGRForecast(data,shpAll,dfForecast$yFutSARD,"../datasets_it/yFutSARDGR.pdf",counterfactual=FALSE)
mapGRForecast(data,shpAll,dfForecast$yFutARD,"../datasets_it/yFutARDGR.pdf",counterfactual=TRUE)
mapGRForecast(data,shpAll,dfForecast$yFutSRD,"../datasets_it/yFutSRDGR.pdf",counterfactual=TRUE)
mapGRForecast(data,shpAll,dfForecast$yFutSAD,"../datasets_it/yFutSADGR.pdf",counterfactual=TRUE)
mapGRForecast(data,shpAll,dfForecast$yFutSAR,"../datasets_it/yFutSARGR.pdf",counterfactual=TRUE)

# Table results ----
replace_numbers = function(x, low=0.01, high=1e3, digits = 3, scipen=-10, ...) {
    nCharx = nchar(x)
    x = gsub(mark,'.',x)
    x.num = as.numeric(x)
    if ((x.num >= low) & (x.num < high)){
        round(x.num, digits = digits)
    } else {
        xNew = prettyNum(x.num, digits=digits, scientific = scipen) 
        # paste(strrep(" ",nCharx-nchar(xNew)),xNew,sep="")
    }
}

names(outSARDEstimate$coefSARD) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD","lambda")
niceOutputSARDCoef = unlist(lapply(outSARDEstimate$coefSARD, FUN = replace_numbers))  
names(outSARDEstimate$se_coefSARD) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD","lambda")
niceOutputSARDse = unlist(lapply(outSARDEstimate$se_coefSARD, FUN = replace_numbers))  

names(outSARDEstimate$WN_SARD_est$coefSARD_WN) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD")
niceOutputWNSARDCoef = unlist(lapply(outSARDEstimate$WN_SARD_est$coefSARD_WN, FUN = replace_numbers))  
names(outSARDEstimate$WN_SARD_est$se_coefSARD_WN) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD")
niceOutputWNSARDse = unlist(lapply(outSARDEstimate$WN_SARD_est$se_coefSARD_WN, FUN = replace_numbers))  

SARDWNIVCoef = coef(outSARDEstimate$IV_estimator$IV_est)
names(SARDWNIVCoef) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD")
niceOutputWNSARDIVCoef = unlist(lapply(SARDWNIVCoef, FUN = replace_numbers))  
SARDWNIVse = sqrt(diag(vcov(outSARDEstimate$IV_estimator$IV_est)))
names(SARDWNIVse) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD")
niceOutputWNSARDIVse = unlist(lapply(SARDWNIVse, FUN = replace_numbers))  

regressors = outSARDEstimate$IV_estimator$IV_est$model
naiveSARD = lm(delta ~ y0 + xS + xA + xR + xD, data = regressors)

naiveSARDCoef = naiveSARD$coefficients
names(naiveSARDCoef) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD")
niceOutputNAIVESARDCoef = unlist(lapply(naiveSARDCoef, FUN = replace_numbers))  

naiveSARDse = sqrt(diag(vcov(naiveSARD)))
names(naiveSARDse) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD")
niceOutputNAIVESARDse= unlist(lapply(naiveSARDse, FUN = replace_numbers))  

k = length(coef(naiveSARD))
N = nrow(regressors)
AICc = AIC(naiveSARD) + (2*k^2+2*k)/(N-k-1)
LL0 = logLik(lm(delta ~ 1, data = regressors))
R2Nagelkerke =  c(1 - exp(-(2/N)*(logLik(naiveSARD) - LL0)))



