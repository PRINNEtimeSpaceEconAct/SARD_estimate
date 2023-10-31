
require(sf)
require(plm)

get_data <- function(fileLocation) {
    # load and returns every file in fileLocation
    # Expected to find exactly  1 RData containing the dataframe
    #                       and 1 ShpFile containing the shape 
    # dataframe contains
    # time, geo, Comune, income (millions of euro), km2, altim (in km)
    # shpFile has units labeled by variable PRO_COM_T
    
    if (DEBUG == TRUE){ print("reading input files") }
    
    allFiles = dir(fileLocation)
    
    dfFile = allFiles[endsWith(allFiles,".RData")]
    load(paste(fileLocation, dfFile, sep=""))   # df
    
    shpFile = allFiles[endsWith(allFiles,".shp")]
    shp_sf=st_read(paste(fileLocation, shpFile, sep=""),
                   promote_to_multi = FALSE,quiet = TRUE)  
    
    
    return(list(df=df, shp=shp_sf))
}

create_df <- function(fileLocation, minDist=21, initialYear=2008, finalYear=2019){
    # input raw dataframe and shapefile
    # create dataframe starting from shp and initial and final data, 
    # and exogenous variables. data contains also coordinates of municipalities,
    # remark: delta already divided by tau
    # data: geo | y0 | yT | delta | ones | s | km2 | Latitude | Longitude
    # shp_data is the shape of only observation considered, that are dose 
    # that have at least one neighbor within distance minDist 
    
    if (DEBUG == TRUE){ print("creating dataframe") }
    
    rawData = get_data(fileLocation)
    df=rawData$df
    shp=rawData$shp
    
    #balance the panel
    df <- df %>% filter(time %in% c(initialYear,finalYear)) %>% na.omit()
    df.p <- pdata.frame(df, index=c("geo", "time")) 
    df.p.b <- make.pbalanced(df.p, balance.type = "shared.individuals")
    
    
    #construct cross-section from panel
    dfInitial <- df.p.b %>% filter(time==initialYear) %>% 
            rename(NameGeo=Comune) %>% mutate(y0 = income/km2, s=altim/1000)
    dfFinal <- df.p.b %>% filter(time==finalYear) %>% 
            mutate(yT = income/km2) %>% select(c(geo,yT))
    data <- dfInitial %>% left_join(dfFinal, by="geo")
    data <- data %>% mutate(delta=(yT-y0)/(finalYear-initialYear))

    #keep common regions
    common <- shp$PRO_COM_T %in% data$geo
    shpCommon  <- shp[common,]

    shpCommon <- st_transform(shpCommon, 
                              "+proj=utm +zone=32 +datum=WGS84 +no_defs")
    shpCommon <- suppressWarnings(st_centroid(shpCommon,
                                              of_largest_polygon = TRUE))
    shpCommon <- st_transform(shpCommon, "+proj=longlat  +datum=WGS84 +no_defs")
    shpCommon <- shpCommon %>%
        dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                      Latitude = sf::st_coordinates(.)[,2])
   
    distances <- compute_D(shpCommon,dMax=minDist)
    
    shpCommonAllNeig <- shpCommon[rowSums(distances)!=0, ]
    dataCommonAllNeig <- data %>% filter(geo %in% shpCommonAllNeig$PRO_COM_T)
    dataCommonAllNeig <- dataCommonAllNeig %>% mutate(ones=1)
    
    #Check that the order of geo is the same
    names1 <- shpCommonAllNeig$PRO_COM_T
    names2 <- dataCommonAllNeig$geo
    #rearrange the order
    orderData <- order(match(names1, names2))
    shpCommonAllNeig <- shpCommonAllNeig[orderData,]
    
    dataCommonAllNeig$Longitude <- shpCommonAllNeig$Longitude
    dataCommonAllNeig$Latitude <- shpCommonAllNeig$Latitude
    
    dataCommonAllNeig <- dataCommonAllNeig %>% select(c(geo,NameGeo,y0,yT,
                                        delta,ones,s,km2,Latitude,Longitude))
    
    data = data.frame(geo = as.character(dataCommonAllNeig$geo), 
                      NameGeo = as.character(dataCommonAllNeig$NameGeo),
                      y0 = as.numeric(dataCommonAllNeig$y0), 
                      yT = as.numeric(dataCommonAllNeig$yT),
                      delta = as.numeric(dataCommonAllNeig$delta), 
                      ones = as.numeric(dataCommonAllNeig$ones),
                      s = as.numeric(dataCommonAllNeig$s), 
                      km2 = as.numeric(dataCommonAllNeig$km2), 
                      Latitude = as.numeric(dataCommonAllNeig$Latitude), 
                      Longitude = as.numeric(dataCommonAllNeig$Longitude))
    
    return(list(data=data, shp_data=shpCommonAllNeig))
}

