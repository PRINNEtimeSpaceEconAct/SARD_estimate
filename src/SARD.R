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

# No constant ---- 
outSARDEstimate_NOCONSTANT = estimate_SARD_auto_NOCONSTANT(data,shp_data,hA,hR)
save(outSARDEstimate_NOCONSTANT,file="../datasets_it/outSARDNOCONSTANTEstimate.RData")

yFutSARD_NOCONSTANT = forcastSARD_NOCONSTANT(50,outSARDEstimate_NOCONSTANT$coefSARD,hA,hR,data,tau=11)
save(yFutSARD_NOCONSTANT,file="../datasets_it/yFutSARD_NOCONSTANT.RData")

load("../datasets_it/yFutSARD_NOCONSTANT50.RData")
data = data.frame(yFutSARD_NOCONSTANT)
data = data[,c(1,ncol(data))]
colnames(data) = c("yT", "VarFinal")
shp_data = cbind(shp_data, yT=data$yT, VarFinal=data$VarFinal)


load("../datasets_it/yFutSARD.RData")
data = data.frame(yFutSARD)
data = data[,c(41,51)]
colnames(data) = c("yT", "VarFinal")
shp_data = cbind(shp_data, yT=data$yT, VarFinal=data$VarFinal)


shp <- shp_data %>% mutate(GR.VarFinal = ifelse(VarFinal > 0, (log(VarFinal/yT)/10*100), 0))
breaks_qt = c(min(shp$VarFinal,na.rm=T), 0, classIntervals((filter(shp, GR.VarFinal>0)$GR.VarFinal), n = 8, style = "quantile")$brks)
cols = c("salmon", brewer.pal(9,"Blues"))
brkslab = c("negative")
for (i in 3:length(breaks_qt)){
    brkslab[i] = paste(round(breaks_qt[i-1],digits=2), " to ", round(breaks_qt[i],digits=2))
}
brkslab=brkslab[-c(2)]
shp <- shp %>% mutate(cities=ifelse(PRO_COM_T=="058091" | PRO_COM_T=="015146" | PRO_COM_T=="063049" |PRO_COM_T=="001272" | PRO_COM_T=="072006" | PRO_COM_T=="082053" | PRO_COM_T=="087015" | PRO_COM_T=="037006" | PRO_COM_T=="048017" | PRO_COM_T=="027042" | PRO_COM_T=="010025" | PRO_COM_T=="083048" | PRO_COM_T=="080063" | PRO_COM_T=="092009", 1, 0))  
cities = shp %>% filter(cities == 1)

tm_obj <- tm_shape(shp) + tm_fill("GR.VarFinal", title="Annual growth rate (%)", breaks=breaks_qt,labels=brkslab,palette=cols) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tmap_options(check.and.fix = TRUE)  + tm_shape(cities) +  tm_borders(alpha = 1, col="gold")
# dev.new()
tm_obj

# level final year 
colnames(data) = c(2019:2069)
shp = shp_data %>% bind_cols(data)
varsIndex = which(colnames(shp)%in%2019:2069)
shp = shp %>% mutate_at(varsIndex, log)
dataLog = shp %>% select(which(colnames(shp)%in%2019:2069)) %>%
    st_drop_geometry() %>% mutate_all(as.numeric)
breaks_qt = classIntervals(c(unlist(dataLog)),n=10,na.rm=T)$brks

stringYear = as.character(2019:2069)
tm_obj <- tm_shape(shp) + tm_fill(stringYear[1], title=expression(Income~per~Km^2~(log)), breaks=breaks_qt, palette ="YlGn", midpoint = NA) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0)  + tm_shape(cities) +  tm_borders(alpha = 1, col="gold") + tm_scale_bar(position=c("left", "bottom"),breaks=c(10,30,50,100),text.size=1); tm_obj; tmap_save(filename = paste(VarName,".pdf", sep=""));



shp <- shp %>% mutate(VarFinal.log = log(VarFinal))
breaks_qt = classIntervals(c(shp$VarFinal.log),n=10,na.rm=T)$brks
breaks_qt[1] = quantile(shp$VarFinal.log, 0.01, na.rm = T)
brkslab = c(1)
for (i in 1:(length(breaks_qt)-1)){
    brkslab[i] = paste(round(breaks_qt[i],digits=2), " to ", round(breaks_qt[i+1],digits=2))
}
tm_obj <- tm_shape(shp) + tm_fill("VarFinal.log", title=expression(Income~per~Km^2~(log)), breaks=breaks_qt, labels=brkslab,palette ="YlGn", midpoint = NA) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0)  + tm_shape(cities) +  tm_borders(alpha = 1, col="gold") + tm_scale_bar(position=c("left", "bottom"),breaks=c(10,30,50,100),text.size=1)
tm_obj


# load results ----
load("../datasets_it/outSARDEstimate.RData")
load("../datasets_it/outDURBINEstimate.RData")

## correlogramm SARD with spatial error ----
maxLag = 10
correlogramSARD_WN = correlogram(outSARDEstimate$WN_SARD_est$residualsSARD_WN,shp_data,maxLag = maxLag)
correlogramSARD = correlogram(outSARDEstimate$residualsSARD,shp_data,maxLag = maxLag)

coeffCorrelogramm = correlogramSARD_WN$correlogram_resid$res[,1]
SECorrelogramm = sqrt(correlogramSARD_WN$correlogram_resid$res[,3])
coeffCorrelogrammFilter = correlogramSARD$correlogram_resid$res[,1]
SECorrelogrammFilter = sqrt(correlogramSARD$correlogram_resid$res[,3])


dev.new()
plot(1:maxLag,coeffCorrelogramm,xlab="lags",ylab="Moran's I",pch=18,ylim=c(-0.02,0.25),xlim=c(1,maxLag))
axis(side = 1, at = 1:maxLag)
arrows(1:maxLag,coeffCorrelogramm-1.96*SECorrelogramm,1:maxLag,coeffCorrelogramm+1.96*SECorrelogramm,length=0.1,angle=90,code=3)
points(1:maxLag,coeffCorrelogrammFilter,xlab="lags",ylab="Moran's I",pch=18,col="red")
arrows(1:maxLag,coeffCorrelogrammFilter-1.96*SECorrelogrammFilter,1:maxLag,coeffCorrelogrammFilter+1.96*SECorrelogrammFilter,length=0.1,angle=90,code=3,col="red")
abline(h=0)
legend("topright",legend=c("residual WN SARD","residual SARD"),col=c("black","red"),lty=1,lwd=1.5,cex=1.0)
dev.copy2pdf(file="../datasets_it/correlogramm.pdf")


load("../datasets_it/outDURBINEstimate.RData")
DURBINCorrelogram = correlogram(DURBINEstimate$DURBIN$residuals,shp_data,maxLag)
dev.new()
plot(DURBINCorrelogram$correlogram_resid,ylim=c(-0.02,0.25),main="")
dev.copy2pdf(file="../datasets_it/correlogramDURBIN.pdf")

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
shpAll = get_data(folderName)$shp
load("../datasets_it/dfForecast2069.RData")
mapYFinalForecast(dfForecast,data,shpAll,"../datasets_it/map_incomeKm2_2069SARD.pdf")
mapGRForecast(data,shpAll,dfForecast,dfForecast$yFutSARD,"../datasets_it/deltaSARDnorm.pdf",counterfactual=FALSE)

# alpha = outSARDEstimate$coefSARD[1]
# phi = outSARDEstimate$coefSARD[2]
# y = dfForecast$yFutSARD
# for (i in 1:50){
#     y = y + alpha + phi*y
# }
# dfForecast$yFutSARD = y

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

# Scatterplot GR vs each  ----

data = data %>% mutate(gObs=log(yT/y0)/11, gS=GR.SARD-GR.S, gA=GR.SARD-GR.A, gR=GR.SARD-GR.R, gD=GR.SARD-GR.D)

plot(data$gS, data$gObs, cex=0.5, pch=19)
abline(v=0)
abline(h=0)

library(sm)
library(latex2exp)

dev.new()
#plot(log(data$y0), data$gS, cex=0.5, pch=19)
sm.regression(log(data$y0), data$gS, col="red", lwd=2, ylim=c(-0.5,0.7), pch=1, cex=0.0, se=T, xlab=TeX(r'(Income per Km$^2$ in 2008 (log))'), ylab="")
title(ylab="Contribution to annual growth rate",mgp=c(2,0,0))
sm.regression(log(data$y0), data$gA, col="green", lwd=2, ylim=c(-1,1), pch=1, cex=0.01, se=T, add=T)
sm.regression(log(data$y0), data$gR, col="blue", lwd=2, ylim=c(-1,1), pch=1, cex=0.01, se=T, add=T)
sm.regression(log(data$y0), data$gD, col="purple", lwd=2, ylim=c(-1,1), pch=1, cex=0.01, se=T, add=T)
grid()
abline(h=0,lty=1)
legend("topright", col=c("red", "green", "blue", "purple"), c(TeX(r'($\hat{g}_S$)'), TeX(r'($\hat{g}_A$)'), TeX(r'($\hat{g}_R$)'), TeX(r'($\hat{g}_D$)')),lty=1,lwd=2,ncol=4)
dev.copy2pdf(file="../datasets_it/gjCrossSectional.pdf")


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

incomeLag = lm(delta ~ y0, data = regressors)
mean((incomeLag$residuals)^2)

SincomeLag = lm(delta ~ y0 + xS, data = regressors)
mean((SincomeLag$residuals)^2)

AltincomeLag = lm(delta ~ y0 + s , data = data)
mean((AltincomeLag$residuals)^2)


k = length(coef(WN_SARD_OLS))
N = nrow(regressors)
AICc = AIC(WN_SARD_OLS) + (2*k^2+2*k)/(N-k-1)
LL0 = logLik(lm(delta ~ 1, data = regressors))
R2Nagelkerke =  c(1 - exp(-(2/N)*(logLik(WN_SARD_OLS) - LL0)))


MsDeriv = GFDM(data)
D = compute_D(data,longlat=TRUE)
WhA = compute_WhAR(D,data,hA)
WhR = compute_WhAR(D,data,hR)
xS = compute_xS(data,MsDeriv)
xA = compute_xAR(data,MsDeriv, WhA)
xR = compute_xAR(data,MsDeriv, WhR)
xD = compute_xD(data,MsDeriv)
MS = compute_MSLag(data,MsDeriv)
MA = compute_MARLag(data,MsDeriv,WhA)
MR = compute_MARLag(data,MsDeriv,WhR)
MD = compute_MDLag(MsDeriv)
WN_SARD_OLS_BYSARD = LogLikAICcR2(data, c(coef(WN_SARD_OLS),0), 10, xS, xA, xR, xD, 
                      MS, MA, MR, MD, diag(nrow(data)))
NAIVE_SARD_OLS_BYSARD = LogLikAICcR2(data, c(coef(naiveSARD),0,0,0,0,0), 6, xS, xA, xR, xD, MS, MA, MR, MD, diag(nrow(data)))


# cross-sectional distribution 
# density estimation 
source("../../../../adaptive density estimation R code latest/bivariateKernelDensity.R")
source("../../../../adaptive density estimation R code latest/bivariateKernelDensityAdaptive.R")
source("../../../../adaptive density estimation R code latest/numericalIntegration.r")
source("../../../../adaptive density estimation R code latest/univariateDensityAdaptive.R")
source("../../../../adaptive density estimation R code latest/weighted.sd.R")

folderName = "../datasets_it/"
rawData = create_df(folderName)
data = rawData$data
load("../datasets_it/dfForecast.RData")
dfForecast = dfForecast %>% filter(yFutSARD > 0)

dens2008 = univariateDensityAdaptive(log((data$y0)/mean(data$y0)),ngrid=5000,method.h="sj")
dens2019 = univariateDensityAdaptive(log(data$yT/mean(data$yT)),ngrid=5000,method.h="sj")
dens2069 = univariateDensityAdaptive(log(dfForecast$yFutSARD/mean(dfForecast$yFutSARD)),ngrid=5000,method.h="sj")
plot(dens2008$eval.points,dens2008$estimate,type="l",xlim=c(-5,5))
lines(dens2019$eval.points,dens2019$estimate,col="red")
lines(dens2069$eval.points,dens2069$estimate,col="blue")


