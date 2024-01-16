

source("testCreateShpVoronoi.R")
library(sfheaders)

shpMontecarlo = voronoiSf

Na = 1000
agents = sf_point(matrix(runif(n = 2*Na),nrow = Na, ncol = 2))

st_crs(agents) ="WGS84"
lengths(st_intersects(shpMontecarlo,agents))




