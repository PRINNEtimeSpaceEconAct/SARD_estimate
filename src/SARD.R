rm(list = ls())

source("lib/L_loadAll.R")
DEBUG = TRUE
PARALLEL = TRUE
NPROCS = 30
initJulia()



folderName = "../datasets_it/"
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
yFutSARD = forcastSARD(11,SARDCoeff,hA,hR,data,tau=11)
save(yFutSARD,file="../datasets_it/yFutSARD.RData")

### ARD ----
ARDCoeff = outSARDEstimate$coefSARD
ARDCoeff["gammaS"] = 0
ARDCoeff["rhoS"] = 0
yFutARD = forcastSARD(11,ARDCoeff,hA,hR,data,tau=11)
save(yFutARD,file="../datasets_it/yFutARD.RData")

### SRD ----
SRDCoeff = outSARDEstimate$coefSARD
SRDCoeff["gammaA"] = 0
SRDCoeff["rhoA"] = 0
yFutSRD = forcastSARD(11,SRDCoeff,hA,hR,data,tau=11)
save(yFutSRD,file="../datasets_it/yFutSRD.RData")

### SAD ----
SADCoeff = outSARDEstimate$coefSARD
SADCoeff["gammaR"] = 0
SADCoeff["rhoR"] = 0
yFutSAD = forcastSARD(11,SADCoeff,hA,hR,data,tau=11)
save(yFutSAD,file="../datasets_it/yFutSAD.RData")

### SAR ----
SARCoeff = outSARDEstimate$coefSARD
SARCoeff["gammaD"] = 0
SARCoeff["rhoD"] = 0
yFutSAR = forcastSARD(11,SARCoeff,hA,hR,data,tau=11)
save(yFutSAR,file="../datasets_it/yFutSAR.RData")

# save
load(file="../datasets_it/yFutSARD.RData")
load(file="../datasets_it/yFutARD.RData")
load(file="../datasets_it/yFutSRD.RData")
load(file="../datasets_it/yFutSAD.RData")
load(file="../datasets_it/yFutSAR.RData")
dfForecast = data.frame(yFutSARD=yFutSARD,yFutARD=yFutARD,yFutSRD=yFutSRD,yFutSAD=yFutSAD,yFutSAR=yFutSAR) 
save(dfForecast,file="../datasets_it/dfForecast.RData")

## Maps ----
load("../datasets_it/dfForecast2069.RData")
shpAll = get_data(folderName)$shp
mapYFinalForecast(dfForecast,data,shpAll,"../datasets_it/map_incomeKm2_2069SARD.pdf")
mapGRForecast(data,shpAll,dfForecast,dfForecast$yFutSARD,"../datasets_it/deltaSARDnorm",counterfactual=FALSE)

load("../datasets_it/dfForecast.RData")
shpAll = get_data(folderName)$shp
mapGRForecast(data,shpAll,dfForecast,dfForecast$yFutARD,"../datasets_it/deltaSInSample.pdf",counterfactual=TRUE)
mapGRForecast(data,shpAll,dfForecast,dfForecast$yFutSRD,"../datasets_it/deltaAInSample.pdf",counterfactual=TRUE)
mapGRForecast(data,shpAll,dfForecast,dfForecast$yFutSAD,"../datasets_it/deltaRInSample.pdf",counterfactual=TRUE)
mapGRForecast(data,shpAll,dfForecast,dfForecast$yFutSAR,"../datasets_it/deltaDInSample.pdf",counterfactual=TRUE)


## Table metro ----
load("../datasets_it/dfForecast.RData")

data = cbind(data, dfForecast)
data = data %>% mutate(GR.SARD = ifelse(yFutSARD > 0, (log(yFutSARD/y0)/11*100), 0), 
                       GR.S = ifelse(yFutARD > 0, (log(yFutARD/y0)/11*100), 0), 
                       GR.A = ifelse(yFutSRD > 0, (log(yFutSRD/y0)/11*100), 0),
                       GR.R = ifelse(yFutSAD > 0, (log(yFutSAD/y0)/11*100), 0), 
                       GR.D = ifelse(yFutSAR > 0, (log(yFutSAR/y0)/11*100), 0))

D = compute_D(data,longlat=TRUE,torus=FALSE)
dCut <- function(d,cut) ((d <= cut) & (d > 0))
WhA = dCut(D,hA)
diag(WhA) = 0
WhR = dCut(D,hR)
diag(WhR) = 0

data = data %>% mutate(WhAy0 = as.numeric(WhA%*%y0),
                       WhRy0 = as.numeric(WhR%*%y0),
                       ratioWhA=WhAy0/y0,
                       ratioWhR=WhRy0/y0)

data = data %>% mutate(cities=ifelse(geo=="058091" | geo=="015146" 
                                      | geo=="063049" |geo=="001272" | geo=="072006" 
                                      | geo=="082053" | geo=="087015" | geo=="037006" 
                                      | geo=="048017" | geo=="027042" | geo=="010025" 
                                      | geo=="083048" | geo=="080063" | geo=="092009", 1, 0))  
cities = data %>% filter(cities == 1) %>% select(NameGeo, y0, WhAy0, WhRy0, ratioWhA, ratioWhR, 
                                                 GR.SARD, GR.S, GR.A, GR.R, GR.D) %>% arrange(desc(y0)) %>%
                                         mutate(NameGeo=tolower(NameGeo))


print(xtable(cities, digits=2), include.rownames=FALSE)


# Table results SARD ----
replace_numbers = function(x, low=0.01, high=1e3, digits = 3, scipen=-10, ...) {
    
    mark  = '::::'
    reg = paste0("([0-9.\\-]+", mark, "[0-9.\\-]+)")
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


regressors = outSARDEstimate$IV_estimator$IV_est$model
WN_SARD_OLS = lm(delta ~ y0 + xS + xA + xR + xD + MSDelta + MADelta + MRDelta + MDDelta, data = regressors)

WN_SARD_OLSCoef = WN_SARD_OLS$coefficients
names(WN_SARD_OLSCoef) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD")
niceOutputWN_SARD_OLSCoef = unlist(lapply(WN_SARD_OLSCoef, FUN = replace_numbers))  

WN_SARD_OLSse = sqrt(diag(vcov(WN_SARD_OLS)))
names(WN_SARD_OLSse) <- c("alpha","phi","gammaS","gammaA","gammaR","gammaD","rhoS","rhoA","rhoR","rhoD")
niceOutputWN_SARD_OLSse= unlist(lapply(WN_SARD_OLSse, FUN = replace_numbers))  


k = length(coef(WN_SARD_OLS))
N = nrow(regressors)
AICc = AIC(WN_SARD_OLS) + (2*k^2+2*k)/(N-k-1)
LL0 = logLik(lm(delta ~ 1, data = regressors))
R2Nagelkerke =  c(1 - exp(-(2/N)*(logLik(WN_SARD_OLS) - LL0)))
