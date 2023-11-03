
# CHECK DF ----

## OLD dataframe ----
library(rgdal)
load("../../df0819Full.RData")
load("../../dbEstimate08-19.RData")

geoId = "053011"

common <- shp_comm$PRO_COM_T %in% df$geo
shp_comm  <- shp_comm[common,]
shp_comm <- spTransform(shp_comm, CRS('+proj=longlat +datum=WGS84 +no_defs'))

geo_shp <- shp_comm@data$PRO_COM_T
geo_df <- unique(df$geo)
iii <- order(match(geo_shp, geo_df))
coord <- coordinates(shp_comm)
coord <- coord[iii,]

df <- df %>% mutate(Longitude = coord[,1],Latitude = coord[,2])
df <- df %>% select(-c(time,income))
df <- df %>% rename(y0.old = income.km2, delta.old = delta, s.old = altim, km2.old = area, NameGeo.old = Comune, 
                    xS.old = gammaS, xA.old = gammaA, xR.old = gammaR, xD.old = gammaD,
                    MSDelta.old = MSDelta, MADelta.old = MADelta, MRDelta.old = MRDelta, 
                    MDDelta.old = MDDelta, Latitude.old = Latitude, Longitude.old = Longitude)
df <- df %>% mutate(geo = as.character(geo))
df.old = df


# OUTPUT -> df.old

## NEW dataframe ---                    
rawData = create_df(folderName)
data = rawData$data
shp_data = rawData$shp_data

MsDeriv = GFDM(data)
D = compute_D(data)

WhA = compute_WhAR(D,data,hA)
WhR = compute_WhAR(D,data,hR)
xS = compute_xS(data,MsDeriv)
xA = compute_xAR(data,MsDeriv, WhA)
xR = compute_xAR(data,MsDeriv, WhR)
xD = compute_xD(data,MsDeriv)
MS = compute_MSLag(data,MsDeriv)
MA = compute_MARLag(data,MsDeriv,WhA,parallel=TRUE)
MR = compute_MARLag(data,MsDeriv,WhR,parallel=TRUE)
MD = compute_MDLag(MsDeriv)

X = as.matrix(cbind(data$ones,data$y0,xS,xA,xR,xD))

# lag regressors
MSDelta = as.numeric(MS %*% matrix(data$delta))
MADelta = as.numeric(MA %*% matrix(data$delta))
MRDelta = as.numeric(MR %*% matrix(data$delta))
MDDelta = as.numeric(MD %*% matrix(data$delta))

data <- data %>% select(-c(yT,ones))
data <- cbind(data, xS,xA,xR,xD,MSDelta,MADelta,MRDelta,MDDelta)
data <- data %>% mutate(geo = as.character(geo))

df.new = data

MS2X=as.matrix(MS %*% MS %*% X)
MA2X=as.matrix(MA %*% MA %*% X)
MR2X=as.matrix(MR %*% MR %*% X)
MD2X=as.matrix(MD %*% MD %*% X)

IV_est = ivreg(delta ~ y0 + xS + xA + xR + xD + MSDelta + MADelta + MRDelta + MDDelta  | y0 + xS + xA + xR + xD + MS2X + MA2X + MR2X + MD2X , data=df.new)
summary(IV_est)




# OUTPUT -> df.new

## MERGE ---

df.tot = df.old %>% inner_join(df.new, by = "geo")
save(df.tot, file = "df.tot.RData")


sum(df.tot$y0 - df.tot$y0.old)
sum(df.tot$delta - df.tot$delta.old)
cbind(df.tot$xD, df.tot$xD.old)
cbind(df.tot$xS, df.tot$xS.old)
cbind(df.tot$xA, df.tot$xA.old)
cbind(df.tot$MSDelta, df.tot$MSDelta.old)
cbind(df.tot$s, df.tot$s.old)
distCoord = sqrt(rowSums((cbind(df.tot$Latitude,df.tot$Longitude) - cbind(df.tot$Latitude.old,df.tot$Longitude.old))^2))


# CHECK SHP LAT LONG ----
library(rgdal)
allFiles = dir(folderName)
shpFile = allFiles[endsWith(allFiles,".shp")]
sf_use_s2(FALSE)


# old
shp.old = readOGR(paste(folderName, shpFile, sep=""))
shp.old <- spTransform(shp.old, "+proj=longlat +datum=WGS84 +no_defs")
coord.old = coordinates(shp.old)

# new
shp.new=st_read(paste(folderName, shpFile, sep=""), promote_to_multi=TRUE, quiet=FALSE)
shp.new <- st_transform(shp.new, "+proj=utm +zone=32 +datum=WGS84 +no_defs")
shp.new <- st_centroid(shp.new, of_largest_polygon = TRUE)
shp.new <- st_transform(shp.new, "+proj=longlat +datum=WGS84 +no_defs")
coord.new = st_coordinates(shp.new)

distCoord = sqrt(rowSums((coord.new - coord.old)^2))
max(distCoord)


