
createShape <- function(N=100, typeOfDist="VoronoiUniform", meanNorm=0.5, 
                               sdNorm=0.25, plot=FALSE){
    # Create the shape file from Voronoi sets of random generated points
    # N is the number of polygons
    # typeOfDist = "VoronoiUniform", "Uniform", "other" 
    # anything of type "other" will make voronoi sets with irregular grid 
    
    # create square polygon
    square = st_polygon(list(cbind(c(0,1,1,0,0), c(0,0,1,1,0))))
    squareBox = st_bbox(c(xmin = 0, xmax = 1, ymax = 1, ymin = 0))
    squareSfc = st_sfc(square,crs = "WGS84")
    
    
    if (typeOfDist != "Uniform"){
        # check the number of feature after intersection
        voronoiSf <- matrix(0, nrow = 1, ncol = 1)
        while(nrow(voronoiSf)<N){
            if (typeOfDist=="VoronoiUniform"){
                voronoiPoints = st_multipoint(matrix(runif(n = 2*N),nrow = N, ncol = 2))
                # x = seq(from = -1, to = 1, length.out = sqrt(N))
                # voronoiPoints = st_multipoint(as.matrix(expand.grid(x,x)))
            }else{
                voronoiPoints = st_multipoint(gaussianMixBoundary(N,10))
            }
            
            # compute voronoy polygon
            voronoiPoly = st_voronoi(st_union(voronoiPoints), squareSfc)
            
            # create sf object
            voronoiGeom = st_geometrycollection(voronoiPoly)
            voronoiGeomColl = st_collection_extract(voronoiGeom)
            voronoiSf = st_sf(geom=voronoiGeomColl,crs = "WGS84")
            attr(st_geometry(voronoiSf), "bbox") = squareBox
            sf_out = st_sf(st_intersection(voronoiSf,squareSfc))
        }
    }
    else {
        grid_st = st_make_grid(squareSfc,n = round(sqrt(N)))
        grid_sf = st_sf(geom=grid_st,crs = "WGS84")
        attr(st_geometry(grid_sf), "bbox") = squareBox
        sf_out = grid_sf
    }
     
   # plot
    if (plot==TRUE){
        plot(sf_out)
        plot(st_centroid(sf_out),add = T,pch=19,cex=0.1)
    }
    
    return(sf_out)
}

# Continuous to shp ----
continuous2shp <- function(X, Y, S, shp_sf){
    # aggregatre levels of S into polygons provided by the shapefile shp_sf
    
    Xstack = c(X)
    Ystack = c(Y)
    Sstack = c(S)
    
    xy = data.frame(x=Xstack, y=Ystack, s=Sstack)
    xys = st_as_sf(xy, coords=c("x","y"))
    st_crs(xys) ="WGS84"
    p_ag1 = aggregate(xys, shp_sf, mean)
    p_ag1$s[is.na(p_ag1$s)] = 0.0
    return(p_ag1$s)
}

gaussianMixBoundary <- function(N,NGaussiansPside=10){
    sigma = 0.05
    NGPS = NGaussiansPside
    centers = seq(0,1,1/NGPS)
    centers = centers[2:(NGPS-1)]
    mu = cbind(rbind(centers,1),rbind(centers,0),rbind(1,centers),
               rbind(0,centers),rbind(0,0),rbind(0,1),rbind(1,0),rbind(1,1))
    Sigma = sigma * diag(1,nrow=2)
    NsEachG = round(N/(4*(NGPS-1)))
    samples = mvrnorm(n = NsEachG, mu = mu[,1],Sigma = Sigma)
    for (i in 2:ncol(mu)){
        samples = rbind(samples,mvrnorm(n = NsEachG, mu = mu[,i],Sigma = Sigma))
    }
    samples = samples[1:N,]
    
    samples[samples[,1]<0,1] = -samples[samples[,1]<0,1]
    samples[samples[,1]>1,1] = 1-(samples[samples[,1]>1,1]-1)
    samples[samples[,2]<0,2] = -samples[samples[,2]<0,2]
    samples[samples[,2]>1,2] = 1-(samples[samples[,2]>1,2]-1)
    
    return(samples) 
}

# Create dataframe with all variables ----
createDataframeAgents <- function(agents0, agentsT, shp_sf, tau, Xs, Ys, S){
    
    # data: geo | y0 | yT | delta | ones | s | km2 | Latitude | Longitude
    # shp_data is the shape of only observation considered, that are dose 
    # that have at least one neighbor within distance minDist 
    
    #km2 of each cell
    km2 <- st_area(shp_sf)/sum(st_area(shp_sf))
    
    agents0Sf = sf_point(agents0)
    st_crs(agents0Sf) ="WGS84"
    aggrAgents0 = lengths(st_intersects(shp_sf,agents0Sf))
    y0 = aggrAgents0/(km2*sum(aggrAgents0))
    # y0 = aggrAgents0/sum(aggrAgents0)
    
    agentsTSf = sf_point(agentsT)
    st_crs(agentsTSf) ="WGS84"
    aggrAgentsT = lengths(st_intersects(shp_sf,agentsTSf))
    yT = aggrAgentsT/(km2*sum(aggrAgentsT))
    # yT = aggrAgentsT/sum(aggrAgentsT)
    
    #codes for cell
    geo = seq(1:length(y0))
    
    #delta y over the period of length tau
    delta = delta=(yT-y0)/tau
    
    #constant vector
    ones = rep(1, length(y0))
    
    #mountains! 
    s=continuous2shp(Xs, Ys, S, shp_sf)
    
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

createDataframePDE <- function(PDE0, PDET, Xpde, Ypde, shp_sf, tau, Xs, Ys, S){
    
    # data: geo | y0 | yT | delta | ones | s | km2 | Latitude | Longitude
    # shp_data is the shape of only observation considered, that are dose 
    # that have at least one neighbor within distance minDist 
    
    #km2 of each cell
    km2 <- st_area(shp_sf)/sum(st_area(shp_sf))
    
    y0 = continuous2shp(Xpde, Ypde, PDE0, shp_sf)
    yT = continuous2shp(Xpde, Ypde, PDET, shp_sf)

    #codes for cell
    geo = seq(1:length(y0))
    
    #delta y over the period of length tau
    delta = delta=(yT-y0)/tau
    
    #constant vector
    ones = rep(1, length(y0))
    
    #mountains! 
    s=continuous2shp(Xs, Ys, S, shp_sf)
    
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


