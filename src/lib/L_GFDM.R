# Cristiano Ricci - 11 Oct 2022
# compute differential matrices Mx My Mxx Myy Mxy for coordinates in coord

require(Matrix)

compute_MDiff <- function(coord){
    # select the ns k-neigh nodes to approximate derivatives
    ns = 20
    Npt = dim(coord)[1]

    x = coord[,1]
    y = coord[,2]
    Dx = matrix(rep(x, Npt), nrow=Npt)  - matrix(rep(x, Npt), nrow=Npt,byrow=TRUE)
    Dy = matrix(rep(y, Npt), nrow=Npt)  - matrix(rep(y, Npt), nrow=Npt,byrow=TRUE)
    Dist = sqrt(Dx^2+Dy^2)

    Mx = matrix(0,nrow=Npt,ncol=Npt)
    My = matrix(0,nrow=Npt,ncol=Npt)
    Mxx = matrix(0,nrow=Npt,ncol=Npt)
    Myy = matrix(0,nrow=Npt,ncol=Npt)
    Mxy = matrix(0,nrow=Npt,ncol=Npt)
    
    print("computing differential matrices Mx, My, Mxx, Myy")
    pb = txtProgressBar(min = 1, max = Npt, initial = 1, style = 3) 
    for (i in 1:Npt){
        idx_neighbors = closestStar(ns,i,Dist)
        
        h = Dx[idx_neighbors,i]
        k = Dy[idx_neighbors,i]
        d = Dist[idx_neighbors,i]
        
        Di = compute_Di(h,k,d)
        Mx[i,c(i,idx_neighbors)]  = Di[1,]
        My[i,c(i,idx_neighbors)]  = Di[2,]
        Mxx[i,c(i,idx_neighbors)] = Di[3,]
        Myy[i,c(i,idx_neighbors)] = Di[4,]
        Mxy[i,c(i,idx_neighbors)] = Di[5,]
        
        setTxtProgressBar(pb,i)
    }
    close(pb)
    
    Mxx <- as(Mxx, "sparseMatrix")
    Mx <- as(Mx, "sparseMatrix")
    Myy <- as(Myy, "sparseMatrix")
    My <- as(My, "sparseMatrix")
    
    return(  list(Mx=Mx,My=My,Mxx=Mxx,Myy=Myy)  )
}
    
closestStar <- function(ns,i,Dist){
    d = Dist[,i]
    idx_neighbors = (sort(d,index.return=T)$ix)[2:(ns+1)]
    return(idx_neighbors)
}

compute_Di <- function(h,k,d){
    R = max(d)
    w = array(data = 0, dim = length(d))
    for (i in 1:length(d)){
        w[i] = weight(d[i],R)
    }
    
    A = compute_A(h,k,w)
    B = compute_B(h,k,w)
    
    Di = solve(A) %*% B
    return(Di)
}

compute_A <- function(h,k,w){
    A = matrix(data = 0,nrow = 5, ncol = 5)
    w2 = w^2
    
    for (i in 1:length(h)){
        A[1,1] = A[1,1] + h[i]^2 * w2[i]
        A[1,2] = A[1,2] + h[i] * k[i] * w2[i]
        A[1,3] = A[1,3] + h[i]^3 * w2[i] / 2
        A[1,4] = A[1,4] + h[i] * k[i]^2 * w2[i] / 2
        A[1,5] = A[1,5] + h[i]^2 * k[i] * w2[i]
        A[2,2] = A[2,2] + k[i]^2 * w2[i]
        A[2,3] = A[2,3] + h[i]^2 * k[i] * w2[i] / 2
        A[2,4] = A[2,4] + k[i]^3 * w2[i] / 2
        A[2,5] = A[2,5] + h[i] * k[i]^2 * w2[i] 
        A[3,3] = A[3,3] + h[i]^4 * w2[i] / 4
        A[3,4] = A[3,4] + h[i]^2 * k[i]^2 * w2[i] / 4
        A[3,5] = A[3,5] + h[i]^3 * k[i] * w2[i] / 2
        A[4,4] = A[4,4] + k[i]^4 * w2[i] / 4
        A[4,5] = A[4,5] + h[i] * k[i]^3  * w2[i] / 2
        A[5,5] = A[5,5] + h[i]^2 * k[i]^2 * w2[i]        
    }
    A = A + t(triu(A,1))
    return(A)
}

compute_B <- function(h,k,w){
    B = matrix(data = 0, nrow = 5, ncol = (length(h)+1) )
    w2 = w^2
    
    B[1,1] = -sum(h*w2)
    B[2,1] = -sum(k*w2)
    B[3,1] = -sum( (h^2) * w2) / 2
    B[4,1] = -sum( (k^2) * w2) / 2
    B[5,1] = -sum( h * k * w2)
    for (i in 1:length(h)){
        B[1,i+1] = h[i] * w2[i]
        B[2,i+1] = k[i] * w2[i]
        B[3,i+1] = h[i]^2 * w2[i] / 2
        B[4,i+1] = k[i]^2 * w2[i] / 2
        B[5,i+1] = h[i] * k[i] * w2[i] 
    }
    return(B)
}


weight <- function(d,dm){
    x = d/dm
    if (x <= 1){
        return(1 - 6 * x^2 + 8 * x^3 - 3 * x^4)
    }
    else{
        return(0)
    }
}


