# using StaticArrays
# using LinearAlgebra
# using ForwardDiff
# using Distributions
# using ProgressMeter

using CellListMap.PeriodicSystems
import CellListMap.wrap_relative_to
import CellListMap.limits, CellListMap.Box
using Random

# Random.seed!(1)

function computeAgents(Nm,Na,tau,SARDp)

    Nm = Int(Nm)
    Na = Int(Na)
    tau = Float64(tau)

    Agents0 = zeros(Nm,Na,2)
    AgentsT = zeros(Nm,Na,2)
    cutoff = max(SARDp.hA, SARDp.hR)

    for m in 1:Nm
        println("("*string(m)*"/"*string(Nm)*")")
        system = init_system(Na = Na, cutoff = cutoff)
        Agents0[m,:,:] .= vecvec2mat(system.positions)
        AgentsT[m,:,:] = sampleAgents(system,Na,tau,SARDp)
    end 

    return Agents0,AgentsT
end

function init_system(;Na::Int=200,cutoff::Float64=0.1)
    unitcell = SVector(1.0,1.0)
    # positions = generatePositions(Na) 
    u0Fun, maxu0 = densityInitialCondition()
    positions = sampleFromDensity(u0Fun,maxu0,Na)   
    system = PeriodicSystem(
        positions=positions,
        cutoff=cutoff,
        unitcell = unitcell,
        output=similar(positions),
        output_name=:interaction,
        parallel=true,
    )
    return system
end

function densityInitialCondition(;Δx = 1e-2)
    Nx = Int(1/Δx)
    x = LinRange(0,1-Δx,Nx)
    y = LinRange(0,1-Δx,Nx)
    X = repeat(x',Nx,1)
    Y = repeat(y,1,Nx)
    Whu0 = make_WDiscrete(Δx,0.1)
    
    
    u0 = make_u0(Δx,Whu0)
    u0 = u0'
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

function sampleAgents(system,Na,tau,SARDp; dt=1e-3)

    nsteps = round(Int,tau/dt)
    X = MvNormal([0.0,0.0],I(2))

    @showprogress for step in 1:nsteps

        # compute pairwise interacitons interaction at this step
        # key point where interactions are computed
        map_pairwise!(
            (x,y,i,j,d2,interaction) -> update_interaction!(x,y,i,j,d2,interaction,Na,SARDp),
            system)

        x = system.positions
        f = system.interaction 
        s = -SARDp.gammaS * ∇Sfun.(system.positions)
        noise = sqrt(2*SARDp.gammaD) * [SVector{2,Float64}(rand(X)) for _ in 1:Na]
        x = x + dt*f + dt*s + sqrt(dt)*noise
        x = fixPositions(x)
        system.positions .= x

        # if particles escape, break
        # if particlesOut(system.positions)
        #     error("particles out of box") 
        # end

    end
    return vecvec2mat(system.positions)

end

function fixPositions(positions)
    for i in eachindex(positions)
        pos = positions[i]
        positions[i] = SVector{2,Float64}([mod(pos[1],1.0),mod(pos[2],1.0)])
    end
    return positions
end

function update_interaction!(x,y,i,j,d2,interaction,Na,SARDp)
    
    dxy = x-y
    d = sqrt(d2)
    denhA = d * (SARDp.hA + d)^3
    denhR = d * (SARDp.hR + d)^3

    
    if d2 == 0 
        drift = zeros(SVector{2})
    else
        drift = - (-2 * dxy) * SARDp.gammaA * (d2 ≤ SARDp.hA^2) / denhA + 
                - (-2 * dxy) * SARDp.gammaR * (d2 ≤ SARDp.hR^2) / denhR
                
        drift = drift / Na

    end
    interaction[i] += drift
    interaction[j] -= drift

    return interaction
end


function Sfun(x,y)
    return exp(-((x-0.6)^2)/0.02)
end

# function Sfun(x,y)
#     # return the level of s for every coordinate x,y ∈ [0,1] × [0,1]
#     b = 0.2

#     # center
#     if (b <= x <= 1-b) && (b <= y <= 1-b)
#         return 0.0
#     # sides
#     elseif (x < b) && (b <= y <= 1-b)
#         return (x-b)^2
#     elseif (x > 1-b) && (b <= y <= 1-b)
#         return (x-(1-b))^2
#     elseif (y < b) && (b <= x <= 1-b)
#         return (y-b)^2
#     elseif (y > 1-b) && (b <= x <= 1-b)
#         return (y-(1-b))^2
#     # corners
#     elseif (x < b) && (y < b)
#         return (x-b)^2+(y-b)^2 # ok
#     elseif (x < b) && (y > 1 - b)
#         return (x-b)^2+(y-(1-b))^2
#     elseif (x > 1-b) && (y < b)
#         return (x-(1-b))^2+(y-b)^2
#     elseif (x > 1-b) && (y > 1-b)
#         return (x-(1-b))^2+(y-(1-b))^2 # ok
#     else # (x > 1 || x < 0 || y > 1 || y < 0)
#         return b^2
#     end
# end

# function Sfun(x,y)
    # return the level of s for every coordinate x,y ∈ [0,1] × [0,1]
#     return (x-0.5)^2+(y-0.5)^2
# end

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
    Δx = 1/Npt
    x = LinRange(0,1-Δx,Npt)
    y = LinRange(0,1-Δx,Npt)
    X = repeat(x',Npt,1)
    Y = repeat(y,1,Npt)
    
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



