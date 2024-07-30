rm(list = ls())
source("lib/L_loadAll.R")
DEBUG = TRUE

folderName = "../datasets_it/"

rawData = create_df(folderName)
data = rawData$data
shp_data = rawData$shp_data

# old coordinates
# load(file = "df.tot.RData")
# data$Latitude = df.tot$Latitude.old
# data$Longitude = df.tot$Longitude.old


# h_range = seq(from=5,to=100,by=5)
# bestDURBIN_est = chose_hDurbin(data,h_range)
# hBest = 25
# DURBINEstimate = estimate_DURBIN_auto(data,hBest)
# save(DURBINEstimate,file="../outDURBINEstimate.RData")

# SLX LAG
cut = 25
D = compute_D(data,longlat=TRUE)
dInvSq <- function(d,cut) 1/d^2 * ((d <= cut) & (d > 0))

W = compute_Wh(D,data,cut)        
trMatc <- suppressWarnings(trW(W, type="MC"))
lw <- suppressWarnings(spdep::mat2listw(W, style="W"))

LAG = lagsarlm(delta ~ y0 + s, data = data,
                                   type="lag", listw=lw,method="MC",
                                   trs=trMatc, quiet=T,zero.policy = T)
SLX = lmSLX(delta ~ y0 + s, data = data,listw=lw,zero.policy = T)



load("../datasets_it/outDURBINEstimate.RData")
DURBINCorrelogram = correlogram(DURBINEstimate$DURBIN$residuals,shp_data,10)
plot(DURBINCorrelogram$correlogram_resid,ylim=c(-1,1))


load("../datasets_it/dfForecast.RData")
hBest = 25
yFutDurbin = forcastSpatial(50,DURBINEstimate$DURBIN,hBest,data,tau=11)
shpAll = get_data(folderName)$shp
mapGRForecast(data,shpAll,dfForecast,yFutDurbin,"../datasets_it/deltaDURBINNorm.pdf",counterfactual=FALSE)

