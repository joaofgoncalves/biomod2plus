
#' Convert raster files from grd to GeoTIFF
#' 
#' A function used to convert grd files in a source folder to 
#' GeoTIFF format

convertToGeoTIFF <- function(inputFolder,outputFolder){
  
  fl <- list.files(inputFolder, pattern= ".grd$", full.names = TRUE)
  
  if(!dir.exists(outFolder)){
    dir.create(outFolder)
  }
  
  for(f in fl){
    
    cat("Converting file:\n",f," ...\n",sep="")
    fn <- file_path_sans_ext(basename(f))
    raster::stack(f) %>%
      raster::writeRaster(filename = paste(outputFolder,"/",fn,".tif",sep=""))
    cat("done!\n\n")
  }
}

countDistinct <- function(x) length(unique(x))

'%!in%' <- function(x,y)!('%in%'(x,y))

abbrevNames <- function(x) paste(str_to_title(unlist(strsplit(gsub("\\.","",x),"\\ "))),collapse="")


projectMultipleScenarios <- function(myBiomodModelOut, myBiomodProj, myBiomodEM, 
                                       modelsToUse, projNames=NULL, rstStackList){
  
  if(is.null(projNames)){
    projNames <- names(rstStackList)
    if(is.null(projNames)){
      stop("projNames is not properly defined!")
    }
  }
  
  for(i in 1:length(rstStackList)){
    
    # Obtain spatiotemporal projections
    myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                      new.env = rstStackList[[i]],
                                      proj.name = projNames[i],
                                      selected.models = modelsToUse,
                                      filtered.meth = NULL,
                                      binary.meth = NULL,
                                      compress = 'gzip',
                                      clamping.mask = TRUE,
                                      output.format = '.grd',
                                      do.stack = TRUE)
    
    # Perform the ensembling of projections
    myBiomodEF <- BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                             binary.meth = c('TSS','ROC','KAPPA'),
                                             EM.output = myBiomodEM,
                                             output.format = '.grd')
    
    # Convert all output raster files to GeoTIFF
    inFolder <- paste(getwd(),"/",sp,"/proj_",projName,sep="")
    outFolder <- paste(inFolder,"/","GeoTIFF", sep="")
    dir.create(outFolder)
    
    # Convert to GeoTIFF raster format
    convertToGeoTIFF(inFolder, outFolder)
    
  } 
}


getQuantileThresh <- function(myBiomodModelOut, qt = 0.75, evalMetric = "TSS"){
  
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  evalDF <- as.data.frame(myBiomodModelEval[evalMetric,"Testing.data",,,])
  quantileThresh <- quantile(evalDF, probs = qt, na.rm=TRUE)
  
  return(quantileThresh)
}


countNumberOfGridRecords <- function(spPoints, refGrid){
  
  
}

#' Intersect raster layers to avoid Maxent projection problems
#' 
#' 
#' @param x A RasterStack object
#' 
#' @return A corrected RasterStack

intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(stack(mask(x, intersect_mask(x))))
}



checkModAlgoRun <- function(biomodModelOut){
  
  allModAlgos <- c("GLM", "GBM", "GAM", "CTA", "ANN", "FDA", "MARS", 
                   "RF", "MAXENT.Phillips", "MAXENT.Tsuruoka")
  builtMods <- get_built_models(biomodModelOut)
  wereBuilt <- vector(mode="logical",length = length(allModAlgos))
  
  for(i in 1:length(allModAlgos)){
    wereBuilt[i] <- ifelse(sum(grepl("GLM",builtMods))>0,TRUE,FALSE)
  }
  names(wereBuilt) <- allModAlgos
  return(wereBuilt)
}


responsePlots <- function(biomodModelOut, Data, modelsToUse, showVars = "all", 
                          fixedVarMetric = 'mean', addMarginalPlot = TRUE,
                          varNames = NULL, marginPlotType = "histogram", plotStdErr = FALSE,
                          plot = FALSE, save = TRUE, filePrefix = "ResponsePlot_", 
                          outFolder = getwd(), width = NA, height = NA, dpi = 300){
  
  if(inherits(Data,"RasterStack")){
    cat("Loading data from RasterStack....")
    Data <- as.data.frame(na.omit(values(Data)))
    cat("done.\n\n")
  }
  
  if(!inherits(Data,"data.frame")){
    stop("Data must be set as a data.frame object!")
  }
  
  if(showVars=="all"){
    showVars <- colnames(Data)
  }
  
  respCurvesList <- list()
  respCurvesList[["DATA"]] <- list()
  respCurvesList[["PLOTS"]] <- list()
  
  k <- 0
  pb <- txtProgressBar(1,length(showVars),style=3)
  
  for(modVar in showVars){
    
    k<-k+1
    respCurves <- as.data.frame(response.plot2(models = modelsToUse,
                                               Data = Data,
                                               show.variables = modVar,
                                               do.bivariate = FALSE,
                                               fixed.var.metric = fixedVarMetric,
                                               plot = FALSE,
                                               save.file = "no"))
    
    ## Make one plot per variable -------------------------------------------------
    ##
    ##
    tmpDF <- data.frame(x = respCurves[,1], 
                        avgResp = apply(respCurves[,-1], 1, mean),
                        stdErr = apply(respCurves[,-1], 1, function(x) sd(x)/sqrt(length(x))))
    
    
    cn <- colnames(respCurves)
    
    ## Create plot with response lines per PA/eval_round
    if(!plotStdErr){
      
      # Start the ggplot object
      g <- ggplot(respCurves, aes_string(x = cn[1]))
      
      # Add one response line per model to the ggplot object
      for(i in 2:(ncol(respCurves))){
        g <- g + geom_line(aes_string(y=cn[i]), alpha=0.2)
      }
      
      # Add the average line and change default theme settings
      g <- g + 
        geom_line(mapping = aes(y = avgResp, x = x),data = tmpDF, size = 1.5) + 
        xlab(ifelse(is.null(varNames),modVar,varNames[k])) + 
        ylab("Response") + 
        theme_bw()
      
      # Create plot with average and std-error curves/bands  
    }else{
      g <- ggplot(tmpDF,aes(x = x)) + 
        #geom_line(aes(y = avgResp + stdErr), linetype="dashed", color="blue") + 
        #geom_line(aes(y = avgResp - stdErr), linetype="dashed", color="blue") + 
        geom_ribbon(aes(ymin = avgResp - 2*stdErr, ymax = avgResp + 2*stdErr), 
                    alpha = 0.3, fill = "dark blue") + 
        geom_line(aes(y = avgResp)) + 
        xlab(modVar) + 
        ylab("Response") + 
        theme_bw()
    }
    
    # Adds a marginal plot to show the distribution of the target variable
    # It can be either a boxplot, density, histogram, ...
    if(addMarginalPlot){
      
      minVal <- min(tmpDF$avgResp)
      g <- g + geom_point(mapping = aes_string(y = minVal, x = modVar), 
                          data = Data, alpha = 0)
      g <- ggMarginal(g, margins = "x", type = marginPlotType)
    }
    
    # Plot response curve?
    if(plot){
      plot(g)
    } 
    
    # Save image per variable?
    if(save){
      suppressMessages(ggsave(filename = paste(outFolder,"/",filePrefix,modVar,".png",sep=""),
                              plot = g, width = width, height = height, dpi = dpi))
    }
    
    # Save data to an object
    respCurvesList[["DATA"]][[modVar]] <- respCurves
    respCurvesList[["PLOTS"]][[modVar]] <- g
    
    setTxtProgressBar(pb, k)
  }
  invisible(respCurvesList)
}




