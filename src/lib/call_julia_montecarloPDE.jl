# PDE
# using DifferentialEquations
# using ImageFiltering
# using LinearAlgebra
# using Distributions

function computePDE(tau,SARDp; Δx = 1e-3)    
    Nx = Int(1/Δx)
    T_span = (0.0,tau)
    x = LinRange(0,1-Δx,Nx)
    y = LinRange(0,1-Δx,Nx)
    X = repeat(x',Nx,1)
    Y = repeat(y,1,Nx)
    S = Sfun.(X,Y)
    ∂xS = ∂x(S)
    ∂yS = ∂y(S)
    WhA = make_WDiscrete(Δx,SARDp.hA)
    WhR = make_WDiscrete(Δx,SARDp.hR)
    Whu0 = make_WDiscrete(Δx,0.1)

    p = (SARDp,∂xS,∂yS,WhA,WhR)

    u0 = make_u0(Δx,Whu0)

    prob = ODEProblem(df!, u0, T_span, p)
    sol = solve(prob,Tsit5(),progress=true,progress_steps = 1)
    return X,Y,sol[1],sol[end]
end

function Sfun(x,y)
    return exp(-((x-0.6)^2)/0.02)
end

# function make_u0(Δx,Whu0)
# # Np=256
# # tau= 0.05
# # SARDp = list(gammaS = 0.0, gammaA = -0.035, gammaR = 0.05, gammaD = 0.105, hA = 0.15, hR = 0.4)
#     Nx = Int(1/Δx)
#     x = 0.0:Δx:(1-Δx)
#     y = 0.0:Δx:(1-Δx)
    
#     center1 = [0.45,0.75]
#     center2 = [0.65,0.75]
#     center3 = [0.5,0.3]

    
#     u0 = (0.6*[pdf(MvNormal(center1,0.005*I(2)),[xi,yi]) for xi in x, yi in y] .+ 
#          .+ 0.45*[pdf(MvNormal(center2,0.005*I(2)),[xi,yi]) for xi in x, yi in y] .+
#          .+ 0.55*[pdf(MvNormal(center3,0.005*I(2)),[xi,yi]) for xi in x, yi in y])

#     u0 .= u0 .+ 0.1
#     # u0 .= imfilter(u0,Whu0,Pad(:circular,Nx,Nx))*Δx^2
#     u0 .= u0/sum(u0*Δx^2)

#     return u0
# end

function make_u0(Δx,Whu0)
    Nx = Int(1/Δx)
    x = 0.0:Δx:(1-Δx)
    y = 0.0:Δx:(1-Δx)
    
    center1 = [0.45,0.75]
    center2 = [0.65,0.75]
    center3 = [0.5,0.3]

    
    u0 = (0.6*[pdf(MvNormal(center1,0.005*I(2)),[xi,yi]) for xi in x, yi in y] .+ 
         .+ 0.45*[pdf(MvNormal(center2,0.005*I(2)),[xi,yi]) for xi in x, yi in y] .+
         .+ 0.55*[pdf(MvNormal(center3,0.005*I(2)),[xi,yi]) for xi in x, yi in y])

    u0 .= u0 .+ 0.1
    # u0 .= imfilter(u0,Whu0,Pad(:circular,Nx,Nx))*Δx^2
    u0 .= u0/sum(u0*Δx^2)

    return u0
end

function df!(du,u,p,t)
    SARDp,∂xS,∂yS,WhA,WhR = p
    Nx = size(u,1)
    Δx = 1/Nx
    WhAu = imfilter(u,WhA,Pad(:circular,Nx,Nx))*Δx^2
    WhRu = imfilter(u,WhR,Pad(:circular,Nx,Nx))*Δx^2

    du .= ( SARDp.gammaS * (∂x(u .* ∂xS) +  ∂y(u .* ∂yS)) + 
          + SARDp.gammaA * (∂x(u .* ∂x(WhAu)) +  ∂y(u .* ∂y(WhAu))) +
          + SARDp.gammaR * (∂x(u .* ∂x(WhRu)) +  ∂y(u .* ∂y(WhRu))) +
          + SARDp.gammaD * Δ(u))
    @show t
end


function Δ(u)
    Nx = size(u,1)
    Δx = 1/Nx
    Δu = similar(u)

    for i in 2:Nx-1, j in 2:Nx-1
        Δu[i,j] = 1/(Δx)^2 * (u[i-1,j] + u[i+1,j] + u[i,j-1] + u[i,j+1] - 4*u[i,j])
    end

    # torus
    for i in 2:Nx-1   
        Δu[i,1] = 1/(Δx)^2 * (u[i-1,1] + u[i+1,1] + u[i,Nx] + u[i,2] - 4*u[i,1])
        Δu[i,Nx] = 1/(Δx)^2 * (u[i-1,Nx] + u[i+1,Nx] + u[i,Nx-1] + u[i,1] - 4*u[i,Nx])
    end
    for j = 2:Nx-1
        Δu[1,j] = 1/(Δx)^2 * (u[2,j] + u[1,j-1] + u[1,j+1] + u[Nx,j] - 4*u[1,j])
        Δu[Nx,j] = 1/(Δx)^2 * (u[Nx-1,j] + u[Nx,j-1] + u[Nx,j+1] + u[1,j] - 4*u[Nx,j])
    end
    Δu[1,1] = 1/(Δx)^2 * (u[2,1] + u[1,2] + u[1,Nx] + u[Nx,1] - 4*u[1,1])
    Δu[1,Nx] = 1/(Δx)^2 * (u[2,Nx] + u[1,Nx-1] + u[Nx,Nx] + u[1,1] - 4*u[1,Nx])
    Δu[Nx,1] = 1/(Δx)^2 * (u[Nx-1,1] + u[Nx,2] + u[Nx,Nx] + u[1,1] - 4*u[Nx,1])
    Δu[Nx,Nx] = 1/(Δx)^2 * (u[Nx-1,Nx] + u[Nx,Nx-1] + u[Nx,1] + u[1,Nx] - 4*u[Nx,Nx])

    return Δu
end

function ∂x(u)
    Nx = size(u,1)
    Δx = 1/Nx

    ∂x = similar(u)
    for i in 2:(Nx-1), j in 1:Nx
        ∂x[i,j] = 1/(2*Δx) * (u[i+1,j] - u[i-1,j])
    end

    # torus
    for j in 1:Nx
        ∂x[1,j] =  1/(2*Δx) * (u[2,j] - u[Nx,j])
        ∂x[Nx,j] = 1/(2*Δx) * (u[1,j] - u[Nx-1,j]) 
    end
    return ∂x
end

function ∂y(u)
    Nx = size(u,1)
    Δx = 1/Nx

    ∂y = similar(u)

    for i in 1:Nx, j in 2:(Nx-1)
        ∂y[i,j] = 1/(2*Δx) * (u[i,j+1] - u[i,j-1])
    end

    # torus
    for i in 1:Nx
        ∂y[i,1] =  1/(2*Δx) * (u[i,2] - u[i,Nx])
        ∂y[i,Nx] = 1/(2*Δx) * (u[i,1] - u[i,Nx-1]) 
    end
    return ∂y
end

function W(x)
    """ smoothing kernel """
    return norm(x) <= 1 ? 1/(2*pi*(log(2)-1/2)) * 1/(norm(x)+1)^2 : 0.0
end

function W(x,h)
    """ rescaled kernel """
    return 1/h^2*W(x/h)
end

function make_WDiscrete(Δx,bandwith)
    Npt = ceil(Int,bandwith/Δx)
    Wd = zeros(2Npt+1,2Npt+1)
    x = -Npt*Δx:Δx:Npt*Δx
    y = -Npt*Δx:Δx:Npt*Δx
    for i in 1:length(x), j in 1:length(y)
        Wd[i,j] = W([x[i],y[j]],bandwith)
    end
    Wd .= Wd/sum(Wd*Δx^2)
    return Wd
end
