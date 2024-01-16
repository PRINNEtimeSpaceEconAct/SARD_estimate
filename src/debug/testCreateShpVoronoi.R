rm(list = ls())
library(sf)


# create square polygon
square = st_polygon(list(cbind(c(0,1,1,0,0), c(0,0,1,1,0))))
squareBox = st_bbox(c(xmin = 0, xmax = 1, ymax = 1, ymin = 0))
squareSfc = st_sfc(square,crs = "WGS84")

# create points to voronoi
N = 100
voronoiPoints = st_multipoint(matrix(runif(n = 2*N),nrow = N, ncol = 2))
# voronoiPoints = st_multipoint(matrix(rnorm(n = 2*N,mean=0.5,sd=0.25),nrow = N, ncol = 2))

# compute voronoy polygon
voronoiPoly = st_voronoi(st_union(voronoiPoints), squareSfc)

# create sf object
voronoiGeom = st_geometrycollection(voronoiPoly)
voronoiGeomColl = st_collection_extract(voronoiGeom)
voronoiSf = st_sf(geom=voronoiGeomColl,crs = "WGS84")
attr(st_geometry(voronoiSf), "bbox") = squareBox
voronoiSf = st_intersection(voronoiSf,squareSfc)

# plot
plot(voronoiSf)
plot(st_centroid(voronoiSf),add = T,pch=19,cex=0.1)

# st_write(voronoiSf,dsn="shpTest.shp")
# shp_sf=st_read("shpTest.shp",promote_to_multi = FALSE,quiet = TRUE)  
