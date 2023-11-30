
createShapeVoronoi <- function(N=100, typeOfDist="uniform", meanNorm=0.5, 
                               sdNorm=0.25, plot=FALSE){
    # Create the shape file from Voronoi sets of random generated points
    # N is the number of polygons
    
    # create square polygon
    square = st_polygon(list(cbind(c(0,1,1,0,0), c(0,0,1,1,0))))
    squareBox = st_bbox(c(xmin = 0, xmax = 1, ymax = 1, ymin = 0))
    squareSfc = st_sfc(square,crs = "WGS84")
    
    # check the number of feature after intersection
    voronoiSf <- matrix(0, nrow = 1, ncol = 1)
    while(nrow(voronoiSf)<N){
        # create points to voronoi
        if (typeOfDist=="uniform"){
            voronoiPoints = st_multipoint(matrix(runif(n = 2*N),nrow = N, ncol = 2))
        }else{
            voronoiPoints = st_multipoint(matrix(rnorm(n = 2*N, mean=meanNorm, 
                                                       sd=sdNorm),nrow = N, ncol=2))
        }
        
        # compute voronoy polygon
        voronoiPoly = st_voronoi(st_union(voronoiPoints), squareSfc)
        
        # create sf object
        voronoiGeom = st_geometrycollection(voronoiPoly)
        voronoiGeomColl = st_collection_extract(voronoiGeom)
        voronoiSf = st_sf(geom=voronoiGeomColl,crs = "WGS84")
        attr(st_geometry(voronoiSf), "bbox") = squareBox
        voronoiSf = st_intersection(voronoiSf,squareSfc)
    }       
     
       # plot
        if (plot==TRUE){
            plot(voronoiSf)
            plot(st_centroid(voronoiSf),add = T,pch=19,cex=0.1)
        }
        
        return(voronoiSf)
}

# Create mountains ----
createMountains <- function(X, Y, S, shp_sf){
    # aggregatre levels of S into polygons provided by the shapefile shp_sf
    
    Xstack = c(X)
    Ystack = c(Y)
    Sstack = c(S)
    
    xy = data.frame(x=Xstack, y=Ystack, s=Sstack)
    xys = st_as_sf(xy, coords=c("x","y"))
    st_crs(xys) ="WGS84"
    p_ag1 = aggregate(xys, shp_sf, mean)
    
    return(p_ag1$s)
}

# Create dataframe with all variables ----
createDataframe <- function(agents0, agentsT, shp_sf, tau, X, Y, S){
    
    # data: geo | y0 | yT | delta | ones | s | km2 | Latitude | Longitude
    # shp_data is the shape of only observation considered, that are dose 
    # that have at least one neighbor within distance minDist 
    
    agents0Sf = sf_point(agents0)
    st_crs(agents0Sf) ="WGS84"
    aggrAgents0 = lengths(st_intersects(shp_sf,agents0Sf))
    y0 = aggrAgents0/sum(aggrAgents0)
    
    agentsTSf = sf_point(agentsT)
    st_crs(agentsTSf) ="WGS84"
    aggrAgentsT = lengths(st_intersects(shp_sf,agentsTSf))
    yT = aggrAgentsT/sum(aggrAgentsT)
    
    #codes for cell
    geo = seq(1:length(y0))
    
    #delta y over the period of length tau
    delta = delta=(yT-y0)/tau
    
    #constant vector
    ones = rep(1, length(y0))
    
    #mountains! 
    s=createMountains(X, Y, S, shp_sf)
    
    #km2 of each cell
    km2 <- st_area(shp_sf)/1000000
    
    #Longitute and Latitude
    shpCentroids <- suppressWarnings(st_centroid(shp_sf,
                                                 of_largest_polygon = TRUE))
    
    shpCentroids <- st_transform(shpCentroids, "+proj=longlat  +datum=WGS84 +no_defs")
    shpCentroids <- cbind(shpCentroids, Longitude = sf::st_coordinates(shpCentroids)[,1],
                       Latitude = sf::st_coordinates(shpCentroids)[,2])
    
    Longitude <- shpCentroids$Longitude
    Latitude <- shpCentroids$Latitude
    
    data = data.frame(geo = as.character(geo), 
                      y0 = as.numeric(y0), 
                      yT = as.numeric(yT),
                      delta = as.numeric(delta), 
                      ones = as.numeric(ones),
                      s = as.numeric(s), 
                      km2 = as.numeric(km2), 
                      Latitude = as.numeric(Latitude), 
                      Longitude = as.numeric(Longitude))
    
    shp_sf$Id = as.character(1:nrow(shp_sf))
    shp_sf = shp_sf %>% left_join(data, by=c("Id"="geo"))
    
    return(listN(data, shp_sf))
    
}




