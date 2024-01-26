rm(list=ls())
library(spatstat.geom)
library(plotly)


sampleFromDensity <- function(densityFun,maxDensity,N,xlim,ylim){
    samples = suppressWarnings(matrix(data=c(0,0),nrow=0,ncol=2))
    while (nrow(samples) < N){
        Nmissing = N - nrow(samples)
        x = cbind(runif(Nmissing,min=xlim[1],max=xlim[2]),runif(Nmissing,min=ylim[1],max=ylim[2]))
        y = runif(Nmissing,min=0,max=maxDensity)
        acceptedSamples = x[(y <= densityFun(x[,1],x[,2])),]
        samples = rbind(samples,acceptedSamples)
    }
    samples = samples[1:N,]
    return(samples)
}

computeIntegral <- function(densityFun,xlim,ylim,Deltax = 1e-2){
    x = seq(from=xlim[1],to=xlim[2],by=Deltax)
    y = seq(from=ylim[1],to=ylim[2],by=Deltax)
    points = as.matrix(expand.grid(x,y))
    Int = Deltax^2 * sum(densityFun(points[,1],points[,2]))
    return(Int)
}


# set A = [3,4] x [5,6]
# set B = [-2,0] x [0,1]
a1 = 3
b1 = 4
c1 = 5
d1 = 6
a2 = 3
b2 = 4
c2 = 4
d2 = 5
xlimA = c(a1,b1); ylimA = c(c1,d1);
xlimB = c(a2,b2); ylimB = c(c2,d2);
centerA = c((a1+b1)/2,(c1+d1)/2)
centerB = c((a2+b2)/2,(c2+d2)/2)
DistCenter = norm(centerA-centerB,"2")


# fANoNorm <- function(x,y){ (((a1 <= x) | (x <= b1)) & ((c1 <= y) | (y <= d1))) * 1/((b1-a1)*(d1-c1))  }
fANoNorm <- function(x,y){ (((a1 <= x) | (x <= b1)) & ((c1 <= y) | (y <= d1))) * exp(-((x-xlimA[2])^2+(y-ylimA[2])^2)/0.1)   }
IntFaNoNorm = computeIntegral(fANoNorm,xlimA,ylimA)
fA <- function(x,y){ fANoNorm(x,y)/IntFaNoNorm }

# fBNoNorm <- function(x,y){ (((a2 <= x) | (x <= b2)) & ((c2 <= y) | (y <= d2))) * 1/((b2-a2)*(d2-c2))  }
fBNoNorm <- function(x,y){ (((a2 <= x) | (x <= b2)) & ((c2 <= y) | (y <= d2))) * exp(-((x-xlimB[1])^2+(y-ylimB[1])^2)/0.1)   }
IntFbNoNorm = computeIntegral(fBNoNorm,xlimB,ylimB)
fB <- function(x,y){ fBNoNorm(x,y)/IntFbNoNorm }

Na = 10000
Nb = 10000
XA = sampleFromDensity(fA,1/IntFaNoNorm,Na,xlimA,ylimA)
XB = sampleFromDensity(fB,1/IntFbNoNorm,Nb,xlimB,ylimB)

meanA = apply(XA,MARGIN = 2,FUN = mean)
meanB = apply(XB,MARGIN = 2,FUN = mean)
DistMean = norm(meanA-meanB,"2")

NgridAx = 100
NgridAy = 100
NgridBx = 100
NgridBy = 100
x_A = seq(from=a1,to=b1,length.out=NgridAx)
y_A = seq(from=c1,to=d1,length.out=NgridAy)
x_B = seq(from=a2,to=b2,length.out=NgridBx)
y_B = seq(from=c2,to=d2,length.out=NgridBy)
xi_A = as.matrix(expand.grid(x_A,y_A))
yj_B = as.matrix(expand.grid(x_B,y_B))
fi_A = fA(xi_A[,1],xi_A[,2])
fj_B = fB(yj_B[,1],yj_B[,2])


df = data.frame(x = xi_A[,1], y = xi_A[,2], z = fi_A)
plot_ly(df,x = ~x,y = ~y, z = ~z,type = "scatter3d",size = 0.1, showlegend = FALSE)
df = data.frame(x = yj_B[,1], y = yj_B[,2], z = fj_B)
plot_ly(df,x = ~x,y = ~y, z = ~z,type = "mesh3d",size = 0.1, showlegend = FALSE)


DistanceSample = crossdist.default(XA[,1],XA[,2],XB[,1],XB[,2])
MeanDistanceSample = mean(DistanceSample)

DistanceEvalPoints = crossdist.default(xi_A[,1],xi_A[,2],yj_B[,1],yj_B[,2])
f_AB = as.matrix(fi_A) %*% as.matrix(t(fj_B))
MeanPdfPoints = abs(diff(xlimA))*abs(diff(ylimA))/((NgridAx-1)*(NgridAy-1))*abs(diff(xlimB))*abs(diff(ylimB))/((NgridBx-1)*(NgridBy-1)) * sum(DistanceEvalPoints*f_AB)

meanAx = sum(xi_A[,1]*fi_A)*abs(diff(xlimA))*abs(diff(ylimA))/((NgridAx-1)*(NgridAy-1))
meanAy = sum(xi_A[,2]*fi_A)*abs(diff(xlimA))*abs(diff(ylimA))/((NgridAx-1)*(NgridAy-1))
meanA
meanAx
meanAy

# MeanPdfPoints = 0.0
# for (i in 1:(NgridAx*NgridAy)){
#     print(round(i/(NgridAx*NgridAy)*100,digits=1))
#     for (j in 1:(NgridBx*NgridBy)){
#         MeanPdfPoints = MeanPdfPoints + abs(diff(xlimA))*abs(diff(ylimA))*abs(diff(xlimB))*abs(diff(ylimB))/(NgridAx*NgridAy*NgridBx*NgridBy)*norm(xi_A[i,]-yj_B[j,],"2")*fi_A[i]*fj_B[j]
#     }
# }
 


DistCenter
DistMean
MeanDistanceSample
MeanPdfPoints

