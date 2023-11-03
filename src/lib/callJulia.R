# this part of the code is extremely platform dependent. fix it by yourself. 
# after the command initJulia() has been run successfully, everything else will. 
# how to make it work in every OS, only god knows
# sincerely, 
# the authors

# on MAC: 
# if necessary on mac press cancel the first time julia is run
# then go to Settings -> Privacy and Security -> Open Anyway and run again


require(JuliaCall)

initJulia <- function(){
    if (DEBUG == TRUE){ print("initializing Julia") }
    
    # check if Julia is installed. If not install it 
    # pinned to 1.8 to avoid change in syntax in future versions
    tryCatch({julia <- julia_setup(verbose=FALSE,installJulia=FALSE)}, 
        warning = function(w){}, 
        error = function(e) { 
        print("Julia not present. Installing Julia 1.8")
        print("On MAC after installation go to")
        print("Settings -> Privacy and Security -> Open Anyway and run again")
        install_julia(version="1.8.0")
        julia <- julia_setup(verbose=FALSE,installJulia=FALSE)})

    # load and install packages if needed
    julia_install_package_if_needed("Optimization")
    julia_install_package_if_needed("OptimizationOptimJL")
    julia_install_package_if_needed("SparseArrays")
    julia_install_package_if_needed("Distributions")
    julia_install_package_if_needed("FiniteDiff")
    
    julia_library("LinearAlgebra")
    julia_library("Optimization")
    julia_library("OptimizationOptimJL")
    julia_library("SparseArrays")
    julia_library("Distributions")
    julia_library("FiniteDiff")
    
    if (DEBUG == TRUE){ julia_command('println("Julia is working as intended")') }
}

call_julia_LogLik_WN <- function(X,Y,MS,MA,MR,MD,initialCondition){
    
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


