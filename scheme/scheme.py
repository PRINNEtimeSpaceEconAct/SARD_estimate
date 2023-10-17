

def get_data(fileLocation):
    """
    read all data from files
    """
    shp = geo = y0 = yT = s = None
    return shp, geo, y0, yT, s

def create_df(fileLocation):
    """
    create dataframe starting from shp and initial and final data, 
    and exogenous variables. data contains also coordinates of municipalities,
    data: geo | y0 | yT | delta | s | coord (°)
    """
    shp, geo, y0, yT, s = get_data(fileLocation)
    data = None    
    return data

def LogLikAICc(data, coefs, dof, xS, xA, xR, xD, MS, MA, MR, MD, Weps):
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

def compute_WhA(D,hA):
    """
    compute the weight matrix WhA given the mtarices of all the distances 
    and the cut-off threshold hA
    """
    WhA = None
    return WhA

def compute_WhR(D,hR):
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

def compute_xA(data,MsDeriv, WhA):
    """
    compute regressor related to gamma_A
    """
    y0 = None # data$y0
    xA = None
    return xA

def compute_xR(data,MsDeriv, WhR):
    """
    compute regressor related to gamma_R
    """
    y0 = None # data$y0
    xR = None
    return xR


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

def compute_MALag(data,MsDeriv,WhA):
    """
    compute the lag matrix MA
    """
    y0 = delta = None # data$y0, data$delta
    MA = None
    return MA

def compute_MRLag(data,MsDeriv,WhR):
    """
    compute the lag matrix MR
    """
    y0 = delta = None # data$y0, data$delta
    MR = None
    return MR

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
    xA = compute_xA(data,MsDeriv, WhA)
    xR = compute_xR(data,MsDeriv, WhR)
    MA = compute_MALag(data,MsDeriv,WhA)
    MR = compute_MRLag(data,MsDeriv,WhR)

    MS2X = MA2X = MR2X = MD2X = None # instruments for IV

    IV_est = None # estimate via IV
    coefs = dof = Weps = None # coefs = IV_est$coefs + 0 lambda, dof = 6, Weps = I
    LogLik, AICc = LogLikAICc(data, coefs, dof, xS, xA, xR, xD, MS, MA, MR, MD, Weps)

    return IV_est, AICc

def chose_hAhR(data,hA_range,hR_range):
    MsDeriv = GFDM(data)
    D = compute_D(data,dMax = 1.0)

    xS = compute_xS(data,MsDeriv)
    xD = compute_xD(data,MsDeriv)
    MS = compute_MSLag(data,MsDeriv)
    MD = compute_MDLag(MsDeriv)

    for (hA_i,hR_i) in (hA_range,hR_range):
        IV_est_i, AICc_i = estimate_IV_SARD(data,MsDeriv,xS,xD,MS,MD,D,hA_i,hR_i)
    
    hABest = hRBest = None  # minimum of (AICc_i)

    return hABest,hRBest




# input
fileLocation = "cartella"   # shp, data
hA_range = [1,2]
hR_range = [3,4]
data = create_df(fileLocation)

hABest,hRBest = chose_hAhR(data,hA_range,hR_range)
# da qui in poi maximum likelihood



