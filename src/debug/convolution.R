rm(list=ls())
library(Matrix)

N = 5
xcoord = seq(from=0,to=1-1/N,length.out=N)
ycoord = seq(from=0,to=1-1/N,length.out=N)
coords = expand.grid(xcoord,ycoord)
# coords = data.frame(x1=runif(N^2),x2=runif(N^2))

plot(coords[,1],coords[,2],col = "white")
text(coords[,1],coords[,2],labels = rownames(coords),col="red")

y0 = runif(N*N)
# y0 = rep(1,N^2)
y0 = y0 / sum(y0) * N^2
km2 = rep(1/N^2,N^2)
df = data.frame(y0=y0,Longitude=coords[,1],Latitude=coords[,2],km2=km2)

xRep = matrix(rep(df$Longitude,N^2),nrow=N^2)
DxFlat = xRep-t(xRep)
AbsDxFlat = as.matrix(stats::dist(df$Longitude, diag = TRUE, upper = TRUE))
AbsDxTorus = pmin(AbsDxFlat,1-AbsDxFlat)
TorusFlatDists = (AbsDxTorus == AbsDxFlat)
TorusFlatDists = 1*TorusFlatDists + -1*(!TorusFlatDists)
DxTorus = AbsDxTorus*sign(DxFlat)*TorusFlatDists

yRep = matrix(rep(df$Latitude,N^2),nrow=N^2)
DyFlat = yRep-t(yRep)
AbsDyFlat = as.matrix(stats::dist(df$Latitude, diag = TRUE, upper = TRUE))
AbsDyTorus = pmin(AbsDyFlat,1-AbsDyFlat)
TorusFlatDists = (AbsDyTorus == AbsDyFlat)
TorusFlatDists = 1*TorusFlatDists + -1*(!TorusFlatDists)
DyTorus = AbsDyTorus*sign(DyFlat)*TorusFlatDists

D = sqrt(AbsDxTorus^2 + AbsDyTorus^2)


Dx = as.matrix(stats::dist(df$Longitude, diag = TRUE, upper = TRUE))
Dy = as.matrix(stats::dist(df$Latitude, diag = TRUE, upper = TRUE))
D1 = sqrt( pmin( Dx,1-Dx )^2 + pmin( Dy,1-Dy )^2 )

# here same as montecarlo
compute_WhAR <- function(D,df,h){
    dInvSq <- function(d,h){ 1/(2*pi*(log(2)-1/2)) * 1/h^2 * 1/(d/h+1)^2 * ((d <= h) & (d > 0)) }
    Wh = dInvSq(D,h)
    Wh[is.na(Wh)] = 0
    diag(Wh) = 1/(2*pi*(log(2)-1/2)) * 1/h^2
 
    Wh = t(apply(Wh, 1, function(x) x * as.numeric(df$km2)))
    Wh = as(Wh,"sparseMatrix")
    return(Wh)
}
 
hA = 0.2
WhA = compute_WhAR(D,df,hA)
sum((WhA%*%df$y0)*df$km2)

WhA = WhA/sum((WhA%*%df$y0)*df$km2)
sum((WhA%*%df$y0)*df$km2)

hRange = seq(from=0.05,to=0.5,by=0.005)
Int = sapply(hRange, function(hA){ WhA = compute_WhAR(D,df,hA); sum((WhA%*%df$y0)*df$km2)})
plot(hRange,Int,type="b")
abline(h=1)
mean(Int[2:length(Int)])

hRange = seq(from=0.05,to=0.5,by=0.005)
NNeigh = sapply(hRange, function(hA){ sum(D[1,] <= hA)})
plot(hRange,NNeigh,type="b")
