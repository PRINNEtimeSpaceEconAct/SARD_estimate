# test that differential matrix are computed correctly
rm(list=ls())
source("lib/L_GFDM.R")
library(plotly)

DEBUG=TRUE

# coordinates grid
N = 100
x = seq(0,1,length.out=N)
y = seq(0,1,length.out=N)
coord = as.matrix(expand.grid(x,y))

# coordinates rand
# N = 50^2
# coord = matrix(runif(2*N,min=0,max=1),nrow=N,ncol=2)


# computing matrices
ALL_MAT = compute_MDiff(coord,torus=TRUE)
Mx.m = ALL_MAT$Mx
My.m = ALL_MAT$My
Mxx.m = ALL_MAT$Mxx
Myy.m = ALL_MAT$Myy


# values
f <- function(x,y){
    return(0.5*1/sqrt(2*pi)*exp(-((x-0.5)^2+(y-0.5)^2)/0.05))
}

fx <- function(x,y){
    return(0.5*1/sqrt(2*pi)*(-2*(x-0.5)/0.05)*exp(-((x-0.5)^2+(y-0.5)^2)/0.05))
}

z = f(coord[,1],coord[,2])
Mxz = as.matrix(Mx.m)%*%z
fx = fx(coord[,1],coord[,2])

coord[,1] = (coord[,1] + 0.5) %% 1.0
# visualize
df = data.frame(x = coord[,1], y = coord[,2], z = z)
plot_ly(df,x = ~x,y = ~y, z = ~z,type = "mesh3d",size = 0.1, showlegend = FALSE)

df = data.frame(x = coord[,1], y = coord[,2], z = Mxz)
plot_ly(df,x = ~x,y = ~y, z = ~z,type = "scatter3d",size = 0.1, showlegend = FALSE)

df = data.frame(x = coord[,1], y = coord[,2], z = fx)
plot_ly(df,x = ~x,y = ~y, z = ~z,type = "scatter3d",size = 0.1, showlegend = FALSE)

df = data.frame(x = coord[,1], y = coord[,2], z = fx-Mxz)
plot_ly(df,x = ~x,y = ~y, z = ~z,type = "scatter3d",size = 0.1, showlegend = FALSE)
