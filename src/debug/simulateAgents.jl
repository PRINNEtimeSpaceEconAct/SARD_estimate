using StaticArrays
using LinearAlgebra
using ForwardDiff
using Distributions
using ProgressMeter
using KernelDensity
using Plots
using Interpolations

# PDE
using LoopVectorization
using DifferentialEquations
using ImageFiltering

# using CellListMap.PeriodicSystems
# import CellListMap.wrap_relative_to
# import CellListMap.limits, CellListMap.Box


# function computeAgents(Nm,Na,tau,SARDp)

#     Nm = Int(Nm)
#     Na = Int(Na)
#     tau = Float64(tau)

#     Agents0 = zeros(Nm,Na,2)
#     AgentsT = zeros(Nm,Na,2)
#     cutoff = max(SARDp.hA, SARDp.hR)

#     for m in 1:Nm
#         println("("*string(m)*"/"*string(Nm)*")")
#         system = init_system(Na = Na, cutoff = cutoff)
#         Agents0[m,:,:] .= vecvec2mat(system.positions)
#         AgentsT[m,:,:] = sampleAgents(system,Na,tau,SARDp)
#     end 

#     return Agents0,AgentsT
# end

# function init_system(;Na::Int=200,cutoff::Float64=0.1)
#     positions = generatePositions(Na)    
#     unitcell = SVector(1.0,1.0)
#     system = PeriodicSystem(
#         positions=positions,
#         cutoff=cutoff,
#         unitcell = unitcell,
#         output=similar(positions),
#         output_name=:interaction,
#         parallel=true,
#     )
#     return system
# end

# function generatePositions(Na)
#     σ = 0.05
#     X = MvNormal([0.5,0.5],σ*I(2))

#     positions = [SVector{2,Float64}(rand(X)) for _ in 1:Na]
#     while any(particleOut.(positions))
#         positions[particleOut.(positions)] = 
#             [SVector{2,Float64}(rand(X)) for _ in 1:sum(particleOut.(positions))]
#     end

#     return positions
# end

# function sampleAgents(system,Na,tau,SARDp; dt=1e-2)

#     nsteps = round(Int,tau/dt)
#     X = MvNormal([0.0,0.0],I(2))

#     @showprogress for step in 1:nsteps

#         # compute pairwise interacitons interaction at this step
#         # key point where interactions are computed
#         map_pairwise!(
#             (x,y,i,j,d2,interaction) -> update_interaction!(x,y,i,j,d2,interaction,Na,SARDp),
#             system)

#         x = system.positions
#         f = system.interaction 
#         s = -SARDp.gammaS * ∇Sfun.(system.positions)
#         noise = sqrt(2*SARDp.gammaD) * [SVector{2,Float64}(rand(X)) for _ in 1:Na]
#         x = x + dt*f + dt*s + sqrt(dt)*noise
#         x = fixPositions(x)
#         system.positions .= x

#         # if particles escape, break
#         # if particlesOut(system.positions)
#         #     error("particles out of box") 
#         # end

#     end
#     return vecvec2mat(system.positions)

# end

# function sampleAgentsInteractive(Na,tau,SARDp; dt=1e-3,isave=1)

#     cutoff = max(SARDp.hA, SARDp.hR)
#     system = init_system(Na=Na,cutoff = cutoff)

#     nsteps = round(Int,tau/dt)
#     X = MvNormal([0.0,0.0],I(2))
#     trajectory = typeof(system.positions)[]
#     push!(trajectory, copy(system.positions))

#     @showprogress for step in 1:nsteps
#         # @show step
#         # compute pairwise interacitons interaction at this step
#         # key point where interactions are computed

#         map_pairwise!((x,y,i,j,d2,interaction) -> update_interaction!(x,y,i,j,d2,interaction,Na,SARDp),system)

#         x = system.positions
#         f = system.interaction 
#         s = -SARDp.gammaS * ∇Sfun.(system.positions)
#         noise = sqrt(2*SARDp.gammaD) * [SVector{2,Float64}(rand(X)) for _ in 1:Na]
#         x = x + dt*f + dt*s + sqrt(dt)*noise
#         x .= fixPositions(x)
#         system.positions .= x

#         if step % isave == 0
#             push!(trajectory, copy(system.positions))
#         end

#     end
#     return trajectory
# end

# function fixPositions(positions)
#     for i in eachindex(positions)
#         pos = positions[i]
#         positions[i] = SVector{2,Float64}([mod(pos[1],1.0),mod(pos[2],1.0)])
#     end
#     return positions
# end

# function update_interaction!(x,y,i,j,d2,interaction,Na,SARDp)
    
#     dxy = x-y
#     d = sqrt(d2)
#     denhA = d * (SARDp.hA + d)^3
#     denhR = d * (SARDp.hR + d)^3
# #    denhA = (SARDp.hA^2 + d2)^2
# #    denhR = (SARDp.hR^2 + d2)^2

    
#     if d2 == 0 
#         drift = zeros(SVector{2})
#     else
#         drift = (-2 * dxy) * SARDp.gammaA * (d2 ≤ SARDp.hA^2) / denhA + 
#                 (-2 * dxy) * SARDp.gammaR * (d2 ≤ SARDp.hR^2) / denhR
                
#         drift = drift / Na

#     end
#     interaction[i] -= drift
#     interaction[j] += drift

#     return interaction
# end

# # function Sfun(x,y)
# #     # return the level of s for every coordinate x,y ∈ [0,1] × [0,1]
# #     b = 0.2

# #     # center
# #     if (b <= x <= 1-b) && (b <= y <= 1-b)
# #         return 0.0
# #     # sides
# #     elseif (x < b) && (b <= y <= 1-b)
# #         return (x-b)^2
# #     elseif (x > 1-b) && (b <= y <= 1-b)
# #         return (x-(1-b))^2
# #     elseif (y < b) && (b <= x <= 1-b)
# #         return (y-b)^2
# #     elseif (y > 1-b) && (b <= x <= 1-b)
# #         return (y-(1-b))^2
# #     # corners
# #     elseif (x < b) && (y < b)
# #         return (x-b)^2+(y-b)^2 # ok
# #     elseif (x < b) && (y > 1 - b)
# #         return (x-b)^2+(y-(1-b))^2
# #     elseif (x > 1-b) && (y < b)
# #         return (x-(1-b))^2+(y-b)^2
# #     elseif (x > 1-b) && (y > 1-b)
# #         return (x-(1-b))^2+(y-(1-b))^2 # ok
# #     else # (x > 1 || x < 0 || y > 1 || y < 0)
# #         return b^2
# #     end
# # end

function Sfun(x,y)
    return exp(-((x-0.6)^2)/0.02)
end

function Sfun(p)
    return (Sfun(p[1],p[2]))
end

function ∇Sfun(p)
    return ForwardDiff.gradient(Sfun,p)
end

function ∇Sfun(x,y)
    return ForwardDiff.gradient(Sfun,[x,y])
end

# function computeS(Npt,Sfun)
#     # return X,Y,S of size (Ne_x x Ne_y) 
#     # with x coordinates, y coordinates and s level

#     Npt = Int(Npt)
#     x = LinRange(0,1,Npt)
#     y = LinRange(0,1,Npt)
#     X = repeat(y',Npt,1)
#     Y = repeat(x,1,Npt)

#     S = reshape([Float64(Sfun(xi,yi)) for xi in x for yi in y],Npt,Npt)
#     return X,Y,S
# end

# function vecvec2mat(x)
#     return reduce(vcat,transpose.(x))
# end

# function mat2vecvec(x)
#     return 
# end

# function particlesOut(positions)
#     if maximum(vecvec2mat(positions)) > 1.0 || minimum(vecvec2mat(positions)) < 0.0
#         return true
#     else 
#         return false
#     end
# end

# function particleOut(position)
#     if maximum(position) > 1.0 || minimum(position) < 0.0
#         return true
#     else 
#         return false
#     end
# end

# function animate(trajectory)
#     for step in trajectory
#         Plots.display(Plots.scatter(Tuple.(step),
#         label=nothing,aspect_ratio=1,lims=(0.0, 1.0),framestyle=:box,ms=1.0))
#         sleep(0.1)
#     end
# end

# function MakeGif(trajectory)
#     for (i,step) in enumerate(trajectory)
#         Plots.scatter(Tuple.(step),
#         label=nothing,aspect_ratio=1,lims=(0.0, 1.0),framestyle=:box)
#         Plots.savefig("Gif" * lpad(string(i),3,"0"))
#     end 
# end


# function plotPoints(points)
#     Plots.scatter(Tuple.(points),
#     label=nothing,aspect_ratio=1,lims=(0.0, 1.0),framestyle=:box,ms=0.2)
# end

# function plotFunc(f)
#     N = 100
#     x = LinRange(0,1,N)
#     y = LinRange(0,1,N)

#     values = [f(xi,yi) for xi in x for yi in y]
#     surface(x,y,values)
# end

# function showDensity(points::Vector{SVector{2, Float64}}) where T
#     # points as vector of vector
#     N = 100
#     x = LinRange(0,1,N)
#     y = LinRange(0,1,N)
#     KDEest = pdf(kde(vecvec2mat(points),boundary=((0.0,1.0),(0.0,1.0))),x,y)
#     heatmap(x,y,KDEest,size=(1000, 1000),cmap=:redsblues,clims = (-0.85, 0.85) .* maximum(abs, KDEest))
# end 

# function showDelta(points1::Vector{SVector{2, Float64}},points2::Vector{SVector{2, Float64}}) where T
#     # points as vector of vector
#     N = 100
#     x = LinRange(0,1,N)
#     y = LinRange(0,1,N)
#     KDEest1 = pdf(kde(vecvec2mat(points1),boundary=((0.0,1.0),(0.0,1.0))),x,y)
#     KDEest2 = pdf(kde(vecvec2mat(points2),boundary=((0.0,1.0),(0.0,1.0))),x,y)
#     hm = heatmap(x,y,KDEest2-KDEest1,size=(1000, 1000),
#         cmap=:redsblues,clims = (-0.85, 0.85) .* maximum(abs, KDEest2-KDEest1))
# end 

# function showDensity(trajectory::Vector{Vector{SVector{2, T}}}) where T
#     # points as vector of vector
#     N = 100
#     x = LinRange(0,1,N)
#     y = LinRange(0,1,N)
#     # t = 0
#     KDEestT = pdf(kde(vecvec2mat(trajectory[end]),boundary=((0.0,1.0),(0.0,1.0))),x,y)
#     contour(x,y,KDEestT,size=(1000, 1000),cbar=false,clabels=true,color=:blue)
#     # t = T
#     KDEest0 = pdf(kde(vecvec2mat(trajectory[1]),boundary=((0.0,1.0),(0.0,1.0))),x,y)
#     contour!(x,y,KDEest0,size=(1000, 1000),cbar=false,clabels=true,color=:red)
# end 

# # test parameter configurations
# SARDp = (gammaS = -0.5, gammaA = -0.1, gammaR = 0.05, gammaD = 0.001, hA = 0.3, hR = 0.35)
# SARDp = (gammaS = 0.0, gammaA = 0.0, gammaR = 0.0, gammaD = 0.1, hA = 0.3, hR = 0.4)
# trajectory = sampleAgentsInteractive(10000,0.1,SARDp; dt = 1e-3, isave = 10)
# animate(trajectory)

############## PDE #############################################################


function computeStep(tau,SARDp; Δx = 1e-2)
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
    xS = tau*SARDp.gammaS*compute_xS(u0,∂xS,∂yS)
    xA = tau*SARDp.gammaA*compute_xA(p,u0,WhA)
    xR = tau*SARDp.gammaR*compute_xR(p,u0,WhR)
    xD = tau*SARDp.gammaD*compute_xD(u0)
    yT = u0 + tau*(xS+xA+xR+xD)
    return yT,xS,xA,xR,xD
end

function computePDEInteractive(tau,SARDp; Δx = 1e-2)    
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

    xS = tau*SARDp.gammaS*compute_xS(u0,∂xS,∂yS)
    xA = tau*SARDp.gammaA*compute_xA(p,u0,WhA)
    xR = tau*SARDp.gammaR*compute_xR(p,u0,WhR)
    xD = tau*SARDp.gammaD*compute_xD(u0)
    return sol,xS,xA,xR,xD
end

function compute_xS(u,∂xS,∂yS)
    Nx = size(u,1)
    Δx = 1/Nx
    xS = (∂x(u .* ∂xS) +  ∂y(u .* ∂yS))
    return xS
end

function compute_xA(p,u,WhA)
    Nx = size(u,1)
    Δx = 1/Nx
    WhAu = imfilter(u,WhA,Pad(:circular,Nx,Nx))*Δx^2
    xA = (∂x(u .* ∂x(WhAu)) +  ∂y(u .* ∂y(WhAu)))
    return xA
end

function compute_xR(p,u,WhR)
    Nx = size(u,1)
    Δx = 1/Nx
    WhRu = imfilter(u,WhR,Pad(:circular,Nx,Nx))*Δx^2
    xR = (∂x(u .* ∂x(WhRu)) +  ∂y(u .* ∂y(WhRu)))
    return xR
end

function compute_xD(u)
    xD = Δ(u)
    return xD
end


function computePDE(tau,SARDp; Δx = 1e-2)    
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
    Whu0 = make_WDiscrete(Δx,0.01)

    p = (SARDp,∂xS,∂yS,WhA,WhR)

    u0 = make_u0(Δx,Whu0)

    prob = ODEProblem(df!, u0, T_span, p)
    sol = solve(prob,Tsit5(),progress=true,progress_steps = 1)
    return sol[end]
end



# function make_u0(Δx,Whu0)
#     Nx = Int(1/Δx)
#     x = 0.0:Δx:(1-Δx)
#     y = 0.0:Δx:(1-Δx)
    
#     center1 = [0.25,0.75]
#     center2 = [0.75,0.70]
#     center3 = [0.75,0.30]

    
#     u0 = (0.5*[pdf(MvNormal(center1,0.05*I(2)),[xi,yi]) for xi in x, yi in y] .+ 
#          .+ 0.25*[pdf(MvNormal(center2,0.01*I(2)),[xi,yi]) for xi in x, yi in y] .+
#          .+ 0.25*[pdf(MvNormal(center3,0.01*I(2)),[xi,yi]) for xi in x, yi in y])

#     u0 .= imfilter(u0,Whu0,Pad(:circular,Nx,Nx))*Δx^2
#     u0 .= u0/sum(u0*Δx^2)

#     return u0
# end


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
    # return norm(x) <= 1 ? 1/(2*pi*(log(2)-1/2)) * 1/(norm(x)+1)^2 : 0.0
    return norm(x) <= 1 ? 1/(2*pi*(log(2)-1/2)) * 1/(norm(x)+1)^2 - 1/4*1/(2*pi*(log(2)-1/2)) : 0.0
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
    # return normalize(Wd,1)
end

function plotPDE(sol)
    Nt = 100
    tau = sol.t[end]
    maxz = maximum(sol)
    for t in LinRange(0,tau,Nt)
        tit = "t = "*string(round(t,digits=3))
        Plots.display(surface(sol(t)',zlims=[0,maxz],size = [1000,1000],title=tit))
        # sleep(0.01)
    end
end

function plotSurf(f)
    Nt = 100
    surface(f,zlims=[min(minimum(f),0),maximum(f)],size = [1000,1000])
end

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

tau = 0.05
SARDp = (gammaS = 0.0, gammaA = -0.035, gammaR = 0.05, gammaD = 0.105, hA = 0.15, hR = 0.4)
sol,xS,xA,xR,xD  = computePDEInteractive(tau,SARDp)
yT,xS,xA,xR,xD = computeStep(tau,SARDp)
plotPDE(sol)
plotSurf(sol[end]-yT)



function densityInitialCondition(;Δx = 1e-2)
    Nx = Int(1/Δx)
    x = LinRange(0,1-Δx,Nx)
    y = LinRange(0,1-Δx,Nx)
    X = repeat(x',Nx,1)
    Y = repeat(y,1,Nx)
    Whu0 = make_WDiscrete(Δx,0.1)
    
    
    u0 = make_u0(Δx,Whu0)
    maxu0 = maximum(u0)
    nodes = (x,y)
    u0Fun = extrapolate(interpolate(nodes,u0,Gridded(Linear(Periodic()))),Periodic())
    
    return u0Fun, maxu0
end

function sampleFromDensity(densityFun,maxDensity,N)
    samples = zeros(0,2)
    while size(samples,1) < N
        Nmissing = N - size(samples,1)
        x = rand(Nmissing,2)
        values = densityFun.(x[:,1],x[:,2])
        y = maxDensity*rand(Nmissing)
        acceptedSamples = x[findall(y .<= values),:]
        samples = vcat(acceptedSamples,samples)
    end
    samples = samples[1:N,:]
    samplesSvec = [SVector{2,Float64}(samples[i,:]) for i in 1:N]
    return samplesSvec
end

Na = 1000
function generatePositions(Na)
    σ = 0.05
    X = MvNormal([0.5,0.5],σ*I(2))

    positions = [SVector{2,Float64}(rand(X)) for _ in 1:Na]
    while any(particleOut.(positions))
        positions[particleOut.(positions)] = 
            [SVector{2,Float64}(rand(X)) for _ in 1:sum(particleOut.(positions))]
    end

    return positions
end

