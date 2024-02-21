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
        quantileArray = shp$value
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


mapYFinalForecast <- function(df,data,shp,fileName){
    
    data = cbind(data,df)
    shp <- shp %>% left_join(select(data,-c("Latitude","Longitude")), by=c("PRO_COM_T"="geo")) %>% rename("geo" = "PRO_COM_T")

    shp <- shp %>% mutate(yFutSARD=ifelse(yFutSARD==0, min(yFutSARD[yFutSARD>0]), yFutSARD)) %>% mutate(logyFutSARD = log(yFutSARD))
    shp <- shp %>% mutate(cities=ifelse(geo=="058091" | geo=="015146" | geo=="063049" |geo=="001272" | geo=="072006" | geo=="082053" | geo=="087015" | geo=="037006" | geo=="048017" | geo=="027042" | geo=="010025" | geo=="083048" | geo=="080063" | geo=="092009", 1, 0))  
    cities = shp %>% filter(cities == 1)
    
    breaks_qt = classIntervals(shp$logyFutSARD,n=10,na.rm=T)$brks
    breaks_qt[1] = quantile(shp$logyFutSARD, 0.01, na.rm = T)
    brkslab = c(1)
    for (i in 1:(length(breaks_qt)-1)){
        brkslab[i] = paste(round(breaks_qt[i],digits=2), " to ", round(breaks_qt[i+1],digits=2))
    }
    
    dev.new()
    tm_obj <- tm_shape(shp) + tm_fill("logyFutSARD", title=expression(Income~per~Km^2~(log)), breaks=breaks_qt, labels=brkslab,palette ="YlGn", midpoint = NA) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0)  + tm_shape(cities) +  tm_borders(alpha = 1, col="gold") + tm_scale_bar(position=c("left", "bottom"),breaks=c(10,30,50,100),text.size=1)
    
    tm_obj
    tmap_save(tm = tm_obj,  filename = fileName)
       
}


mapGRForecast <- function(data,shp,VarFinal,fileName,counterfactual=TRUE){
    
    data = cbind(data,VarFinal)
    shp <- shp %>% left_join(select(data,-c("Latitude","Longitude")), by=c("PRO_COM_T"="geo"))
    shp <- shp %>% mutate(GR.VarFinal = ifelse(VarFinal > 0, (log(VarFinal/yT)/50*100), 0))
    data <- data %>% mutate(GR.VarFinal = ifelse(VarFinal > 0, (log(VarFinal/yT)/50*100), 0))
    
    shp <- shp %>% mutate(cities=ifelse(PRO_COM_T=="058091" | PRO_COM_T=="015146" | PRO_COM_T=="063049" |PRO_COM_T=="001272" | PRO_COM_T=="072006" | PRO_COM_T=="082053" | PRO_COM_T=="087015" | PRO_COM_T=="037006" | PRO_COM_T=="048017" | PRO_COM_T=="027042" | PRO_COM_T=="010025" | PRO_COM_T=="083048" | PRO_COM_T=="080063" | PRO_COM_T=="092009", 1, 0))  
    cities = shp %>% filter(cities == 1)
    

    dev.new()
    if (counterfactual) {
        plotRescaledVar <- plotRescaledVarf(shp=shp, data=data, var.name="GR.VarFinal", id.name="geo", nqq=5, includeZero=F)
        tm_obj <- tm_shape(plotRescaledVar$shp) + tm_fill("varDraw", title="Annual growth rate (%)", breaks=plotRescaledVar$breaks_qt,  palette="RdBu", midpoint =0, labels=plotRescaledVar$lab) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) + tmap_options(check.and.fix = TRUE)
        tm_obj
    }
    else{
        breaks_qt = c(min(shp$VarFinal,na.rm=T), 0, classIntervals((filter(shp, GR.VarFinal>0)$GR.VarFinal), n = 8, style = "quantile")$brks)
        cols = c("salmon", brewer.pal(9,"Blues"))
        brkslab = c("negative")
        for (i in 3:length(breaks_qt)){
            brkslab[i] = paste(round(breaks_qt[i-1],digits=2), " to ", round(breaks_qt[i],digits=2))
        }
        brkslab=brkslab[-c(2)]
        tm_obj <- tm_shape(shp) + tm_fill("GR.VarFinal", title="Annual growth rate (%)", breaks=breaks_qt,labels=brkslab,palette=cols) + tm_layout(frame = FALSE) + tm_borders("white", alpha=0) +     tmap_options(check.and.fix = TRUE)  + tm_shape(cities) +  tm_borders(alpha = 1, col="gold")
        tm_obj
    }
    
    tmap_save(tm = tm_obj,  filename = fileName)
}
