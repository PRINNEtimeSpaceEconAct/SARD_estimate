mapErrorsSameScale <- function(SARD,DURBIN,df,shp){
    
    residSARD = SARD$residualsSARD
    residDurbin = DURBIN$DURBIN$residuals
    df = df %>% bind_cols(residSARD = as.numeric(residSARD), residDurbin = as.numeric(residDurbin))
    quantileArray = c(residDurbin,as.numeric(residSARD))
    
    
    plotRescaledVar <- plotRescaledVarf(shp=shp, quantileArray = quantileArray, data=df, var.name="residSARD", id.name="geo", nqq=5, includeZero=F)
    dev.new()
    tm_obj <- tm_shape(plotRescaledVar$shp) + tm_fill("varDraw", title="Residual", breaks=plotRescaledVar$breaks_qt,  palette="RdBu", midpoint =0, labels=plotRescaledVar$lab) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tmap_options(check.and.fix = TRUE)
    tm_obj
    tmap_save(tm = tm_obj,  filename = "../datasets_montecarlo/residSARD.pdf")

    plotRescaledVar <- plotRescaledVarf(shp=shp, quantileArray = quantileArray, data=df, var.name="residDurbin", id.name="geo", nqq=5, includeZero=F)
    dev.new()
    tm_obj <- tm_shape(plotRescaledVar$shp) + tm_fill("varDraw", title="Residual", breaks=plotRescaledVar$breaks_qt,  palette="RdBu", midpoint =0, labels=plotRescaledVar$lab) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tmap_options(check.and.fix = TRUE)
    tm_obj
    tmap_save(tm = tm_obj,  filename = "../datasets_montecarlo/residDurbin.pdf")
    
    
}

plotRescaledVarf <- function(shpData, quantileArray = NULL, data, var.name, id.name, nqq=5, includeZero=TRUE){

    df = data %>% dplyr::select(value=!!var.name, geo=!!id.name)
    shp = shpData
    shp <- shp %>% left_join(df, by=c("PRO_COM_T"="geo"))
    
    if (is.null(quantileArray)){
        quantileArray = shp@data$value
    }    
    df.quantileArray = data.frame(quantileArray)    
    df.quantileArray = df.quantileArray %>% mutate(quantileArrayDraw = ifelse(quantileArray < 0, -quantileArray/min(df.quantileArray$quantileArray,na.rm = T)*max(df.quantileArray$quantileArray,na.rm = T), quantileArray))
    
    shp = shp %>% mutate(varDraw = ifelse(value < 0, -value/min(shp$value,na.rm = T)*max(shp$value,na.rm = T), value))
    
    breaks_neg = classIntervals(filter(df.quantileArray,quantileArrayDraw < 0)$quantileArrayDraw, n = nqq, style = "quantile")
    breaks_pos = classIntervals(filter(df.quantileArray,quantileArrayDraw > 0)$quantileArrayDraw, n = nqq, style = "quantile")
    break_neg_lab = signif(classIntervals(filter(df.quantileArray,quantileArray < 0)$quantileArray, n = nqq, style = "quantile")$brks,digits=3)
    break_pos_lab = signif(classIntervals(filter(df.quantileArray,quantileArray > 0)$quantileArray, n = nqq, style = "quantile")$brks,digits=3)
    
    if (includeZero==T){        
        breaks_qt = c(breaks_neg$brks,0,0,breaks_pos$brks)
        break_lab = c(break_neg_lab,0,0,break_pos_lab)
    }else{
        breaks_qt = c(breaks_neg$brks,breaks_pos$brks)
        break_lab = c(break_neg_lab,break_pos_lab)
    }
    
    lab = character(length=(length(break_lab)-1))
    for (i in 1:(length(break_lab)-1)){
        lab[i] = paste(break_lab[i]," to ",break_lab[i+1],sep="")
    }
    
    list(shp=shp, breaks_qt=breaks_qt, lab=lab)
}
