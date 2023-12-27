# this part of the code is extremely platform dependent. fix it by yourself. 
# after the command initJulia() has been run successfully, everything else will. 
# how to make it work in every OS, only god knows
# sincerely, 
# the authors

# on MAC: 
# if necessary on mac press cancel the first time julia is run
# then go to Settings -> Privacy and Security -> Open Anyway and run again

# on Linux: 
# add on your .Renviron in /home/YOURUSERNAME
# JULIA_BINDIR=/home/YOURUSERNAME/.local/share/R/JuliaCall/julia/1.8.0/julia-1.8.0/bin
# or alternatively run
# Sys.setenv(JULIA_BINDIR = "/home/YOURUSERNAME/.local/share/R/JuliaCall/julia/1.8.0/julia-1.8.0/bin")
# possibly, run also 
# JuliaCall::julia_setup(rebuild = TRUE)


require(JuliaCall)

guideInstallPrint <- function(){
    print("Julia not present. Installing Julia 1.8")
    print("On MAC after installation go to")
    print("Settings -> Privacy and Security -> Open Anyway and run again")
    print("on Linux: ")
    print("add on your .Renviron in /home/YOURUSERNAME")
    print("JULIA_BINDIR=/home/YOURUSERNAME/.local/share/R/JuliaCall/julia/1.8.0/julia-1.8.0/bin")
    print("or alternatively run")
    print('Sys.setenv(JULIA_BINDIR = "/home/YOURUSERNAME/.local/share/R/JuliaCall/julia/1.8.0/julia-1.8.0/bin")')
    print("possibly, run also ")
    print('JuliaCall::julia_setup(rebuild = TRUE)')
}

initJulia <- function(){
    if (DEBUG == TRUE){ print("initializing Julia") }
    
    # check if Julia is installed. If not install it 
    # pinned to 1.8 to avoid change in syntax in future versions
    tryCatch({julia <- julia_setup(verbose=FALSE,installJulia=FALSE)}, 
        warning = function(w){}, 
        error = function(e) { 
        guideInstallPrint()
        install_julia(version="1.8.0")
        julia <- julia_setup(verbose=FALSE,installJulia=FALSE)})

    # load and install packages if needed
    julia_install_package_if_needed("Optimization")
    julia_install_package_if_needed("OptimizationOptimJL")
    julia_install_package_if_needed("SparseArrays")
    julia_install_package_if_needed("Distributions")
    julia_install_package_if_needed("FiniteDiff")
    julia_install_package_if_needed("StaticArrays")
    julia_install_package_if_needed("CellListMap")
    julia_install_package_if_needed("ForwardDiff")
    julia_install_package_if_needed("ProgressMeter")
    
    julia_library("LinearAlgebra")
    julia_library("Optimization")
    julia_library("OptimizationOptimJL")
    julia_library("SparseArrays")
    julia_library("Distributions")
    julia_library("FiniteDiff")
    julia_library("StaticArrays")
    # julia_library("CellListMap")
    julia_library("ForwardDiff")
    julia_library("ProgressMeter")
    
    
    
    if (DEBUG == TRUE){ julia_command('println("Julia is working as intended")') }
}

call_julia_LogLik_WN <- function(X,Y,MS,MA,MR,MD,initialCondition){
    
    # initJulia()
    julia_command('include("lib/call_julia_LogLik.jl")')
    julia_assign("X",X)
    julia_assign("Y",Y)
    julia_assign("MS",as.matrix(MS))
    julia_assign("MA",as.matrix(MA))
    julia_assign("MR",as.matrix(MR))
    julia_assign("MD",as.matrix(MD))
    julia_assign("initialCondition",initialCondition)
    
    outJulia = julia_eval('julia_LogLik_WN(X,Y,MS,MA,MR,MD,initialCondition)')
    coef = outJulia[[1]]
    se_coef = outJulia[[2]]
    pvalue_coef = outJulia[[3]]
    residuals = outJulia[[4]]
    
    return(listN(coef, se_coef, pvalue_coef, residuals))
}

call_julia_LogLik <- function(X,Y,MS,MA,MR,MD,Weps,initialCondition){
    
    # initJulia()
    julia_command('include("lib/call_julia_LogLik.jl")')
    julia_assign("X",X)
    julia_assign("Y",Y)
    julia_assign("MS",as.matrix(MS))
    julia_assign("MA",as.matrix(MA))
    julia_assign("MR",as.matrix(MR))
    julia_assign("MD",as.matrix(MD))
    julia_assign("Weps",as.matrix(Weps))
    julia_assign("initialCondition",initialCondition)
    
    outJulia = julia_eval('julia_LogLik(X,Y,MS,MA,MR,MD,Weps,initialCondition)')
    coef = outJulia[[1]]
    se_coef = outJulia[[2]]
    pvalue_coef = outJulia[[3]]
    residuals = outJulia[[4]]
    
    return(listN(coef, se_coef, pvalue_coef, residuals))
}

call_julia_computeS <- function(NeS){
    
    # initJulia()
    julia_command
    julia_command('include("lib/call_julia_montecarloAgents.jl")')
    julia_assign("NeS",NeS)
    
    outJulia = julia_eval('computeS(NeS,Sfun)')
    X = outJulia[[1]]
    Y = outJulia[[2]]
    S = outJulia[[3]]
    
    return(listN(X, Y, S))
}

call_julia_computeAgents <- function(Nm,Na,tau,SARDp){
    
    # initJulia()
    julia_command('include("lib/call_julia_montecarloAgents.jl")')

    julia_assign("Nm",Nm)
    julia_assign("Na",Na)
    julia_assign("tau",tau)

    julia_assign("gammaS",SARDp$gammaS)
    julia_assign("gammaA",SARDp$gammaA)
    julia_assign("gammaR",SARDp$gammaR)
    julia_assign("gammaD",SARDp$gammaD)
    julia_assign("hA",SARDp$hA)
    julia_assign("hR",SARDp$hR)
    julia_command("SARDp = (gammaS = gammaS,gammaA = gammaA,gammaR = gammaR,gammaD = gammaD,hA = hA,hR = hR)")
    
    outJulia = julia_eval('computeAgents(Nm,Na,tau,SARDp)')
    agents0_MC = outJulia[[1]]
    agentsT_MC = outJulia[[2]]
    
    return(listN(agents0_MC, agentsT_MC))
}


