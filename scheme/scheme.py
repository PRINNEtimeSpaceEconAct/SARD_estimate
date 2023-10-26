
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput

def get_data(fileLocation):
    """
    read all data from files
    """
    shp = df = None
    return shp, df

def create_df(fileLocation):
    """
    create dataframe starting from shp and initial and final data, 
    and exogenous variables. data contains also coordinates of municipalities,
    df : geo | y0 | yT | delta | ones | s | coordx (°) | coordy (°)
    remark: delta already divided by tau
    data: geo | y0 | yT | delta | ones | s | coordx (°) | coordy (°)
    shp_data is the shape of only observation considered
    """
    shp, df = get_data(fileLocation)
    data = None    
    shp_data = None
    return data, shp_data

def LogLikAICc(data, coef, k, xS, xA, xR, xD, MS, MA, MR, MD, W_eps):
    """
    compute LogLik and AICc of the SARD Model, with dof degree of freedom
    coefs is of length 11, in order
    a, phi, gammaS, gammA, gammaR, gammaD, rhoS, rhoA, rhoR, rhoD, lambda.
    used also for NAIVE, IV WN, SARD WN. 
    """
    LogLik = AICc = None
    return LogLik, AICc

def compute_D(data,dMax):
    """
    compute the matrix D of all the mutual distances (km) from all the points 
    of coordinates coords inside data, whose distance is lower than dMax
    """
    coords = None # data$coords
    D = None
    return D

def GFDM(data):
    """
    compute all the differntial matrices Mx My Mxx Myy
    """
    MsDeriv = None
    return MsDeriv

def compute_WhA(D,hmax):
    """
    compute the weight matrix WhA given the mtarices of all the distances 
    and the cut-off threshold hA
    """
    WhA = None
    return WhA

def compute_WhR(D,hmax):
    """
    compute the weight matrix WhR given the mtarices of all the distances 
    and the cut-off threshold hR
    """
    WhR = None
    return WhR

def compute_xS(data,MsDeriv):
    """
    compute regressor related to gamma_S
    """
    y0 = s = None # data$y0, data$s
    xS = None
    return xS


def compute_xAR(data,MsDeriv, Wh):
    """
    compute regressor related to gamma_A, gamma_R
    """
    y0 = None # data$y0
    xAR = None
    return xAR


def compute_xD(data,MsDeriv):
    """
    compute regressor related to gamma_D
    """
    y0 = None # data$y0
    xD = None
    return xD


def compute_MSLag(data,MsDeriv):
    """
    compute the lag matrix MS
    """
    y0 = s = None # data$y0, data$s
    MS = None
    return MS

def compute_MARLag(data,MsDeriv,Wh):
    """
    compute the lag matrix MA, MR
    """
    y0 = delta = None # data$y0, data$delta
    MAR = None
    return MAR

def compute_MDLag(MsDeriv):
    """
    compute the lag matrix MD = Mxx + Myy
    """
    MD = None
    return MD

def estimate_IV_SARD(data,MsDeriv,xS,xD,MS,MD,D,hA,hR):
    """
    Estimate SARD WN via IV with distances hA, hR
    return AICc
    """

    WhA = compute_WhA(D,hA)
    WhR = compute_WhA(D,hR)
    xA = compute_xAR(data,MsDeriv, WhA)
    xR = compute_xAR(data,MsDeriv, WhR)
    MA = compute_MARLag(data,MsDeriv,WhA)
    MR = compute_MARLag(data,MsDeriv,WhR)

    MS2X = MA2X = MR2X = MD2X = None # instruments for IV

    IV_est = None # estimate via IV
    coef = nparam = W_eps = None # coef = IV_est$coef + 0 lambda, nparam = 6, W_eps = I
    LogLik, AICc = LogLikAICc(data, coef, nparam, xS, xA, xR, xD, MS, MA, MR, MD, W_eps)

    return IV_est, LogLik, AICc

def chose_hAhR(data,hA_range,hR_range):
    MsDeriv = GFDM(data)
    D = compute_D(data,dMax = 1.0)

    xS = compute_xS(data,MsDeriv)
    xD = compute_xD(data,MsDeriv)
    MS = compute_MSLag(data,MsDeriv)
    MD = compute_MDLag(MsDeriv)

    for (hA_i,hR_i) in (hA_range,hR_range):
        IV_est_i, LogLik_i, AICc_i = estimate_IV_SARD(data,MsDeriv,xS,xD,MS,MD,D,hA_i,hR_i)
    
    hABest = hRBest = None  # minimum of (AICc_i)

    return hABest,hRBest


def estimate_WN_SARD(data,hA,hR):
    """
    Estimate WN SARD model (without spatial error) using julia to 
    perform the minimization. 
    """
    MsDeriv = GFDM(data)
    D = compute_D(data,dMax = 1.0)
    WhA = compute_WhA(D,hA)
    WhR = compute_WhA(D,hR)
    xS = compute_xS(data,MsDeriv)
    xA = compute_xAR(data,MsDeriv, WhA)
    xR = compute_xAR(data,MsDeriv, WhR)
    xD = compute_xD(data,MsDeriv)
    MS = compute_MSLag(data,MsDeriv)
    MA = compute_MARLag(data,MsDeriv,WhA)
    MR = compute_MARLag(data,MsDeriv,WhR)
    MD = compute_MDLag(MsDeriv)

    IV_est, LogLik, AICc = estimate_IV_SARD(data,MsDeriv,xS,xD,MS,MD,D,hA,hR)
    coef_initial = None # IV_est$coef
    coef, se_coef, pvalue_coef, residuals = call_julia_LogLik_WN(data,coef_initial,xS,xA,xR,xD,MS,MA,MR,MD)
    nparam = W_eps = None # coef = coef + 0 lambda, nparam = 10, W_eps = I
    LogLik, AICc = LogLikAICc(data, coef, nparam, xS, xA, xR, xD, MS, MA, MR, MD, W_eps)
    allMat = (xS,xA,xR,xD,MS,MA,MR,MD)

    return coef, se_coef, pvalue_coef, residuals, LogLik, AICc, allMat

def call_julia_LogLik_WN(data,coef_initial,xS,xA,xR,xD,MS,MA,MR,MD):
    """
    In Julia minimize -LogLik of the WN SARD (no spatial error)
    return coef, se_coef, pvalue_coef, residuals
    """
    coef = se_coef = pvalue_coef = residuals = None
    return coef, se_coef, pvalue_coef, residuals

def estimate_SARD(data,shp_data,hA,hR):
    """
    Estimate full SARD model with spatial error
    perform first estimate_WN_SARD to get residuals and compute spatial error matrix
    """

    coef_WN, se_coef_WN, pvalue_coef_WN, residuals_WN, LogLik_WN, AICc_WN, allMat_WN = estimate_WN_SARD(data,hA,hR)
    W_eps, lambdas, se_lambdas, pvalue_lambdas = compute_spatial_error_mat(shp_data,residuals_WN,20)

    xS,xA,xR,xD,MS,MA,MR,MD = allMat_WN 
    coef_initial = coef_WN # from estimate_WN_SARD
    
    coef, se_coef, pvalue_coef, residuals = call_julia_LogLik_SARD(data,coef_initial,xS,xA,xR,xD,MS,MA,MR,MD,W_eps)

    nparam = None #  nparam = 11
    LogLik, AICc = LogLikAICc(data, coef, nparam, xS, xA, xR, xD, MS, MA, MR, MD, W_eps)

    return coef, se_coef, pvalue_coef, residuals, LogLik, AICc

def compute_spatial_error_mat(shp_data,residuals,max_cont):
    """
    Compute the matrix of spatial errors by ... da vedere
    rule of thumb to select max lambda is: 
    the first significant followed by two consecutive lambda non significant at 5%
    """
    W_eps = lambdas = se_lambdas = pvalue_lambdas = None
    return W_eps, lambdas, se_lambdas, pvalue_lambdas

def call_julia_LogLik_SARD(data,coef_initial,xS,xA,xR,xD,MS,MA,MR,MD,W_eps):
    """
    In Julia minimize -LogLik of the full SARD model
    return coef, se_coef, pvalue_coef, residuals
    """
    coef = se_coef = pvalue_coef = residuals = None
    return coef, se_coef, pvalue_coef, residuals

# input
fileLocation = "cartella"   # shp, data
hA_range = [1,2]
hR_range = [3,4]


def main(fileLocation,hA_range,hR_range):
    data, shp_data = create_df(fileLocation)
    hABest,hRBest = chose_hAhR(data,hA_range,hR_range)
    coef, se_coef, pvalue_coef, residuals, LogLik, AICc = estimate_SARD(data,shp_data,hABest,hRBest)

    return None

graphviz = GraphvizOutput()
graphviz.output_file = 'basic.svg'
with PyCallGraph(output=graphviz):   
    main(fileLocation,hA_range,hR_range)
