rm(list = ls())

source("lib/L_loadAll.R")
DEBUG = TRUE
PARALLEL = TRUE
NPROCS = 4

folderName = "../datasets_it/"
rawData = create_df(folderName)
data = rawData$data
shp_data = rawData$shp_data


# 500
N = 500
df = data[1:N,]
MsDeriv = GFDM(df)
D = compute_D(df)
Wh = compute_WhAR(D,df,50)
t0 = Sys.time()
MAR = compute_MARLag(df,MsDeriv,Wh)
elapsed = Sys.time() - t0 
print(elapsed)

# 1000
N = 1000
df = data[1:N,]
MsDeriv = GFDM(df)
D = compute_D(df)
Wh = compute_WhAR(D,df,50)
t0 = Sys.time()
MAR = compute_MARLag(df,MsDeriv,Wh)
elapsed = Sys.time() - t0 
print(elapsed)

# 2000
N = 2000
df = data[1:N,]
MsDeriv = GFDM(df)
D = compute_D(df)
Wh = compute_WhAR(D,df,50)
t0 = Sys.time()
MAR = compute_MARLag(df,MsDeriv,Wh)
elapsed = Sys.time() - t0 
print(elapsed)

# all
df = data
MsDeriv = GFDM(df)
D = compute_D(df)
Wh = compute_WhAR(D,df,50)
t0 = Sys.time()
MAR = compute_MARLag(df,MsDeriv,Wh)
elapsed = Sys.time() - t0 
print(elapsed)

