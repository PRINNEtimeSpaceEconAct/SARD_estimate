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

h_range = seq(from=5,to=100,by=5)
bestDURBIN_est = chose_hDurbin(data,h_range)
    