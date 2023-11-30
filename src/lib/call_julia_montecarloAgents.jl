# using StaticArrays
# using LinearAlgebra
# using ForwardDiff
# using Distributions
# using ProgressMeter
# using KernelDensity
# using Plots

using CellListMap.PeriodicSystems
import CellListMap.wrap_relative_to
import CellListMap.limits, CellListMap.Box




function computeAgents(Nm,Na,tau,SARDp)

    Nm = Int(Nm)
    Na = Int(Na)
    tau = Float64(tau)

    Agents0 = zeros(Nm,Na,2)
    AgentsT = zeros(Nm,Na,2)
    cutoff = max(SARDp.hA, SARDp.hR)

    @showprogress for m in 1:Nm
        system = init_system(Na = Na, cutoff = cutoff)
        Agents0[m,:,:] .= vecvec2mat(system.positions)
        AgentsT[m,:,:] = sampleAgents(system,Na,tau,SARDp)
    end 

    return Agents0,AgentsT
end

function init_system(;Na::Int=200,cutoff::Float64=0.1)
    positions = generatePositions(Na)    
    unitcell = [1.0, 1.0]
    system = PeriodicSystem(
        positions=positions,
        cutoff=cutoff,
        unitcell = [1.0, 1.0],
        output=similar(positions),
        output_name=:interaction,
        parallel=true,
    )
    return system
end

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

function sampleAgents(system,Na,tau,SARDp; dt=1e-2)

    nsteps = round(Int,tau/dt)
    X = MvNormal([0.0,0.0],I(2))

    for step in 1:nsteps

        # compute pairwise interacitons interaction at this step
        # key point where interactions are computed
        map_pairwise!(
            (x,y,i,j,d2,interaction) -> update_interaction!(x,y,i,j,d2,interaction,Na,SARDp),
            system)

        x = system.positions
        f = system.interaction 
        s = SARDp.gammaS * ∇Sfun.(system.positions)
        noise = sqrt(2*SARDp.gammaD) * [SVector{2,Float64}(rand(X)) for _ in 1:Na]
        x = x + dt*f + dt*s + sqrt(dt)*noise
        system.positions .= x

        # if particles escape, break
        # if particlesOut(system.positions)
        #     error("particles out of box") 
        # end


    end
    return vecvec2mat(system.positions)

end

function sampleAgentsInteractive(Na,tau,SARDp; dt=1e-3,isave=1)

    cutoff = max(SARDp.hA, SARDp.hR)
    system = init_system(Na=Na,cutoff = cutoff)

    nsteps = round(Int,tau/dt)
    X = MvNormal([0.0,0.0],I(2))
    trajectory = typeof(system.positions)[]
    push!(trajectory, copy(system.positions))

    @showprogress for step in 1:nsteps

        # compute pairwise interacitons interaction at this step
        # key point where interactions are computed
        map_pairwise!((x,y,i,j,d2,interaction) -> update_interaction!(x,y,i,j,d2,interaction,Na,SARDp),system)

        x = system.positions
        f = system.interaction 
        s = SARDp.gammaS * ∇Sfun.(system.positions)
        noise = sqrt(2*SARDp.gammaD) * [SVector{2,Float64}(rand(X)) for _ in 1:Na]
        x = x + dt*f + dt*s + sqrt(dt)*noise
        system.positions .= x

            
        if step % isave == 0
            push!(trajectory, copy(system.positions))
        end

        # # if particles escape, break
        # if particlesOut(system.positions)
        #     error("particles out of box") 
        # end

    end
    return trajectory
end


function update_interaction!(x,y,i,j,d2,interaction,Na,SARDp)
    
    dxy = x-y
    # drift = 1/Na * (SARDp.gammaA * ∇h(dxy,d2,SARDp.hA) + SARDp.gammaR * ∇h(dxy,d2,SARDp.hR))
    
    den = (sqrt(d2) * (-1 + sqrt(d2))^3)
    drift = -2 * dxy * (SARDp.gammaA * (d2 ≤ SARDp.hA^2) + 
            + SARDp.gammaR * (d2 ≤ SARDp.hR^2)) / (Na * den)

    interaction[i] += drift
    interaction[j] -= drift

    return interaction
end


# function ∇h(dxy,d2,h)
#     # dxy = X^i - X^j (vector 2 x 1)
#     # d2 = distance squared
#     # h = cutoff squared
#     # h(x,y) = 1/(dxy+1)^2 * (Ind(d2) ≤ h2)
#     x = dxy[1]
#     y = dxy[2]
#     if d2 ≤ h^2
#         return [-(2*x)/(sqrt(d2) * (-1 + sqrt(d2))^3),
#          -(2*y)/(sqrt(d2) * (-1 + sqrt(d2))^3)]
#     else
#         return zeros(2)
#     end
# end


function Sfun(x,y)
    # return the level of s for every coordinate x,y ∈ [0,1] × [0,1]
    b = 0.2

    # center
    if (b <= x <= 1-b) && (b <= y <= 1-b)
        return 0.0
    # sides
    elseif (x < b) && (b < y < 1-b)
        return -1/b * x + 1
    elseif (x > 1-b) && (b < y < 1-b)
        return 1/b * x + (b-1)/b
    elseif (y < b) && (b < x < 1-b)
        return -1/b * y + 1
    elseif (y > 1-b) && (b < x < 1-b)
        return 1/b * y + (b-1)/b
    # corners
    elseif (x < b) && (y < b)
        return 1-(1/b * (x+1-b) + (b-1)/b)*(1/b * (y+1-b) + (b-1)/b) # ok
    elseif (x < b) && (y > 1 - b)
        return 1-(1/b * (x+1-b) + (b-1)/b)*(-1/b * (y-(1-b)) + 1)
    elseif (x > 1-b) && (y < b)
        return 1-(-1/b * (x-(1-b)) + 1)*(1/b * (y+1-b) + (b-1)/b)
    elseif (x > 1-b) && (y > 1-b)
        return 1-(-1/b * (x-(1-b)) + 1)*(-1/b * (y-(1-b)) + 1) # ok
    else # (x > 1 || x < 0 || y > 1 || y < 0)
        return 1.0
    end
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

function computeS(Npt,Sfun)
    # return X,Y,S of size (Ne_x x Ne_y) 
    # with x coordinates, y coordinates and s level

    Npt = Int(Npt)
    x = LinRange(0,1,Npt)
    y = LinRange(0,1,Npt)
    X = repeat(y',Npt,1)
    Y = repeat(x,1,Npt)

    S = reshape([Float64(Sfun(xi,yi)) for xi in x for yi in y],Npt,Npt)
    return X,Y,S
end

function vecvec2mat(x)
    return reduce(vcat,transpose.(x))
end

function mat2vecvec(x)
    return 
end

function particlesOut(positions)
    if maximum(vecvec2mat(positions)) > 1.0 || minimum(vecvec2mat(positions)) < 0.0
        return true
    else 
        return false
    end
end

function particleOut(position)
    if maximum(position) > 1.0 || minimum(position) < 0.0
        return true
    else 
        return false
    end
end


function animate(trajectory)
    for step in trajectory
        Plots.display(Plots.scatter(Tuple.(step),
        label=nothing,aspect_ratio=1,lims=(0.0, 1.0),framestyle=:box,ms=0.1))
        sleep(0.1)
    end
end

function MakeGif(trajectory)
    for (i,step) in enumerate(trajectory)
        Plots.scatter(Tuple.(step),
        label=nothing,aspect_ratio=1,lims=(0.0, 1.0),framestyle=:box)
        Plots.savefig("Gif" * lpad(string(i),3,"0"))
    end 
end


function plotPoints(points)
    Plots.scatter(Tuple.(points),
    label=nothing,aspect_ratio=1,lims=(0.0, 1.0),framestyle=:box,ms=0.1)
end

function plotFunc(f)
    N = 100
    x = LinRange(0,1,N)
    y = LinRange(0,1,N)

    values = [f(xi,yi) for xi in x for yi in y]
    surface(x,y,values)
end

function showDensity(points::Vector{SVector{2, Float64}}) where T
    # points as vector of vector
    N = 100
    x = LinRange(0,1,N)
    y = LinRange(0,1,N)
    KDEest = pdf(kde(vecvec2mat(points),boundary=((0.0,1.0),(0.0,1.0))),x,y)
    surface(x,y,KDEest,size=(1000, 1000))
end 

function showDensity(trajectory::Vector{Vector{SVector{2, T}}}) where T
    # points as vector of vector
    N = 100
    x = LinRange(0,1,N)
    y = LinRange(0,1,N)
    # t = 0
    KDEestT = pdf(kde(vecvec2mat(trajectory[end]),boundary=((0.0,1.0),(0.0,1.0))),x,y)
    contour(x,y,KDEestT,size=(1000, 1000),cbar=false,clabels=true,color=:blue)
    # t = T
    KDEest0 = pdf(kde(vecvec2mat(trajectory[1]),boundary=((0.0,1.0),(0.0,1.0))),x,y)
    contour!(x,y,KDEest0,size=(1000, 1000),cbar=false,clabels=true,color=:red)
end 



# # DEBUG
# # generate montecarlo samples
# SARDp = (gammaS = -0.5, gammaA = -0.1, gammaR = 0.05, gammaD = 0.001, hA = 0.3, hR = 0.35)
# Agents0, AgentsT = computeAgents(10,1000,1.0,SARDp)

# # test parameter configurations
# # SARDp = (gammaS = -0.5, gammaA = -0.1, gammaR = 0.05, gammaD = 0.001, hA = 0.3, hR = 0.35)
# SARDp = (gammaS = -0.05, gammaA = -0.2, gammaR = 0.1, gammaD = 0.001, hA = 0.1, hR = 0.3)
# trajectory = sampleAgentsInteractive(10000,0.5,SARDp; dt = 1e-2, isave = 1)
# animate(trajectory)
# showDensity(trajectory[1])
# showDensity(trajectory[end])
# showDensity(trajectory)

# # generate points for S
# X,Y,S = computeS(1000,Sfun)
