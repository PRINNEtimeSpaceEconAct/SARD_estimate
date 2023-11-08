chose_hDurbin <- function(df,h_range){
    # select best h for Durbin model based on the AICc score
    
    print("            selecting best h for Durbin model")
    
    D = compute_D(df)
    dInvSq <- function(d,cut) 1/d^2 * ((d <= cut) & (d > 0))
    
    allEstimated_DURBIN = list()
    allAICc = list()
    for (i in 1:length(h_range)){
        
        h = h_range[i]
        print(paste("    Estimating Durbin (",i,"/",length(h_range),"), with h = ",
                    h, sep=""))
        
        W = compute_Wh(D,df,h)        
        trMatc <- suppressWarnings(trW(W, type="MC"))
        lw <- suppressWarnings(spdep::mat2listw(W, style="W"))
        
        
        DURBIN = suppressWarnings(lagsarlm(delta ~ y0 + s, data = df,
                               type="Durbin", listw=lw,method="MC",
                               trs=trMatc, quiet=T,zero.policy = T))
        
        allAICc[[i]] = AICc(DURBIN)[1]
        allEstimated_DURBIN[[i]] = DURBIN
    }
    
    iBest = which.min(allAICc)
    hBest = h_range[iBest]
    DURBIN = allEstimated_DURBIN[[iBest]]
    DURBIN_AICc = allAICc[[iBest]]
    DURBIN_R2N = AICc2R2Nagelkerke(DURBIN_AICc,df$delta,6)
    
    return(listN(hBest,DURBIN,DURBIN_AICc,DURBIN_R2N,allAICc,allEstimated_DURBIN))

}

estimate_DURBIN_auto <- function(df,h){
    # Estimate DURBIN with distances h
    
    print("Estimating Durbin for given h")
    
    D = compute_D(df)
    dInvSq <- function(d,cut) 1/d^2 * ((d <= cut) & (d > 0))

    W = compute_Wh(D,df,h)        
    trMatc <- suppressWarnings(trW(W, type="MC"))
    lw <- suppressWarnings(spdep::mat2listw(W, style="W"))
    
    DURBIN = suppressWarnings(lagsarlm(delta ~ y0 + s, data = df,
                                       type="Durbin", listw=lw,method="MC",
                                       trs=trMatc, quiet=T,zero.policy = T))
    
    DURBIN_AICc = AICc(DURBIN)[1]
    DURBIN_R2N = AICc2R2Nagelkerke(DURBIN_AICc,df$delta,6)
    
    return(listN(DURBIN,DURBIN_AICc,DURBIN_R2N))
    
}