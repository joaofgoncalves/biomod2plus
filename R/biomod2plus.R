
#' Convert raster files from grd format to GeoTIFF
#' 
#' A function used to convert rasters in grd file format from a source 
#' folder to GeoTIFF format to an output folder
#' 
#' @param inputFolder Input folder with .grd files
#' @param  outputFolder Output folder where GeoTIFF raster 
#' files will be placed
#' 
#' @return NULL
#' 
#' @export

convertToGeoTIFF <- function(inputFolder,outputFolder){
  
  fl <- list.files(inputFolder, pattern= ".grd$", full.names = TRUE)
  
  if(!dir.exists(outFolder)){
    dir.create(outFolder)
  }
  
  for(f in fl){
    
    cat("Converting file:\n",f," ...\n",sep="")
    fn <- tools::file_path_sans_ext(basename(f))
    raster::stack(f) %>%
      raster::writeRaster(filename = paste(outputFolder,"/",fn,".tif",sep=""))
    cat("done!\n\n")
  }
}

#' Count distinct elements
#' 
#' A help function used to count the number of distinct elements in a
#' vector
#' 
#' @param x An input vector
#' 
#' @return An integer with the number of unique elements in the input vector
#' 
#' @export
#' 

countDistinct <- function(x) length(unique(x))

#' Not in function
#' 
#' A simple not in function operator
#' 
#' @param x A vector
#' @param y Another vector
#' 
#' @return A logical vector containing the elements of x that are not in y
#' @export
#' 

'%not_in%' <- function(x,y) !('%in%'(x,y))

#' Abbreviate a string to use as a biomod2 species name
#' 
#' A simple function used to abbreviate a species name, remove spaces and points 
#' to use it in biomod2
#' 
#' @param x A string
#' 
#' @return A simplified string
#' 
#' @export
#' 

abbrevNames <- function(x) paste(str_to_title(unlist(strsplit(gsub("\\.","",x),"\\ "))),
                                 collapse="")


#' Project multiple scenarios
#' 
#' A wrapper function used to project and ensemble multiple scenarios contained in a 
#' list object with one \code{RasterStack} per scenario. Useful for forecasting (for example) 
#' multiple climate change RCP scenarios or future dates (or both).
#' 
#' @param myBiomodModelOut A \link[biomod2]{BIOMOD.models.out} class object.
#' @param myBiomodProj A \link[biomod2]{BIOMOD.projection.out} class object.
#' @param myBiomodEM A \link[biomod2]{BIOMOD.EnsembleModeling.out} class object
#' @param modelsToUse A character vector listing the models to use.
#' @param projNames The character vector with names for each different scenario. 
#' Must equal the length of \code{rstStackList}.
#' @param rstStackList A list object with one \code{RasterStack} object by scenario. 
#' The list may be named and these names will be used as the input of \code{proj.name} 
#' parameter from \link[biomod2]{BIOMOD_Projection} function.
#'
#' @return NULL
#' 
#' @seealso 
#' \link[biomod2]{BIOMOD_Projection}, \link[biomod2]{BIOMOD_EnsembleForecasting}, 
#' \link{convertToGeoTIFF}
#' 
#' @export
#'  

projectMultipleScenarios <- function(myBiomodModelOut, myBiomodProj, myBiomodEM, 
                                       modelsToUse, scenarioNames=NULL, rstStackList, 
                                       convert2GeoTIFF=TRUE){
  
  if(is.null(scenarioNames)){
    scenarioNames <- names(rstStackList)
    if(is.null(scenarioNames)){
      stop("scenarioNames is not properly defined!")
    }
  }
  
  for(i in 1:length(rstStackList)){
    
    # Obtain spatiotemporal projections
    myBiomodProj <- BIOMOD_Projection(modeling.output = myBiomodModelOut,
                                      new.env = rstStackList[[i]],
                                      proj.name = scenarioNames[i],
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
    if(convert2GeoTIFF){
      convertToGeoTIFF(inFolder, outFolder)
    }
    
  } 
}

#' Calculate the quantile of an evaluation metric
#' 
#' An ancillary function to calculate a target quantile of an evaluation metric 
#' calculated by \code{biomod2} hold-out evaluation. This function can be used to 
#' calculate the threshold used in \link[biomod2]{BIOMOD_EnsembleModeling} function 
#' to set the parameter \code{eval.metric.quality.threshold}.
#' 
#' @param myBiomodModelOut A \link[biomod2]{BIOMOD.models.out} class object.
#' @param qt The target quantile between 0 and 1.
#' @param evalMetric A character string defining the evaluation metric. 
#' Available options are: 'KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 
#' 'POD', 'CSI' and 'ETS'.
#' 
#' @return The target quantile value.
#' 
#' @export
#' 

getEvalMetricQuantile <- function(myBiomodModelOut, qt = 0.75, evalMetric = "TSS"){
  
  if(length(qt)>1){
    qt<-qt[1]
    warning("qt parameter length is greater than 1. Using only the first element!")
  }
  
  myBiomodModelEval <- biomod2::get_evaluations(myBiomodModelOut)
  evalDF <- as.data.frame(myBiomodModelEval[evalMetric,"Testing.data",,,])
  quantileThresh <- quantile(evalDF, probs = qt, na.rm=TRUE)
  
  return(quantileThresh)
}


#' Intersect raster layers to avoid Maxent projection problems
#' 
#' A help function provided by \code{biomod2} creators to avoid projection 
#' errors in MAXENT.Phillips
#' 
#' @param x A RasterStack object
#' 
#' @return A corrected RasterStack
#' @export
#' 

intersectMask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(stack(mask(x, intersect_mask(x))))
}


#' Check which model algorithms where run
#' 
#' A simple function used to verify which models were actually run by \code{biomod2} 
#' for a given target species.
#' 
#' @param biomodModelOut A \link[biomod2]{BIOMOD.models.out} class object.
#' 
#' @return A named logical vector specifying which models were run.
#' @export
#' 

checkModAlgoRun <- function(biomodModelOut){
  
  allModAlgos <- c("GLM", "GBM", "GAM", "CTA", "ANN", "FDA", "MARS", 
                   "RF", "MAXENT.Phillips", "MAXENT.Tsuruoka")
  builtMods <- biomod2::get_built_models(biomodModelOut)
  wereBuilt <- vector(mode = "logical", length = length(allModAlgos))
  
  for(i in 1:length(allModAlgos)){
    wereBuilt[i] <- ifelse(sum(grepl("GLM",builtMods)) > 0,TRUE,FALSE)
  }
  names(wereBuilt) <- allModAlgos
  return(wereBuilt)
}


## PLOTS AND GRAPHICS ----------------------------------------------------------------------

#' Make response plots for multiple variables 
#' 
#' A flexiblw function used to produce response plots using \code{ggplot2} package. 
#' 
#' @param biomodModelOut A \link[biomod2]{BIOMOD.models.out} class object.
#' @param Data A \code{data.frame} or \code{RasterStack} containing the data used to make 
#' the response plots. They must have the same names as the ones used to calibrate the model.
#' @param modelsToUse
#' @param showVars
#' @param fixedVarMetric
#' @param addMarginalPlot
#' @param varNames
#' @param marginPlotType
#' @param plotStdErr
#' @param plot
#' @param save
#' @param filePrefix
#' @param fileSuffix
#' @param outFolder
#' @param width
#' @param height
#' @param dpi
#' 
#' @return A list object containing two components:     
#' \itemize{
#'   \item DATA - inside this, a list with one element for each predictor variable with the data 
#'   used to make the response plots
#'   \item PLOTS - the \code{ggplot} objects. One for each variable.
#' }
#' 
#' @details 
#' This function uses Adaptation of the Evaluation Strip method proposed by Elith et al.(2005) 
#' and implemented in the \link[biomod2]{response.plot2}. 
#' This function enables to plot the response curves of a model independently of the algorithm 
#' used for building the model. It therefore permits a direct comparison of predicted responses 
#' from the different statistical approaches on the same data.
#' 
#' @seealso 
#' \link[biomod2]{response.plot2}   
#' \link[biomod2]{BIOMOD_LoadModels}
#' @export
#' 

responsePlots <- function(biomodModelOut, Data, modelsToUse, showVars = "all", 
                          fixedVarMetric = 'mean', addMarginalPlot = TRUE,
                          varNames = NULL, marginPlotType = "histogram", plotStdErr = FALSE,
                          plot = FALSE, save = TRUE, filePrefix = "ResponsePlot_", 
                          fileSuffix = "", outFolder = getwd(), width = NA, 
                          height = NA, dpi = 300){
  
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
    
    ## Make one plot per variable ------------------------------------------------- ##
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
      suppressMessages(ggsave(filename = paste(outFolder,"/",filePrefix,modVar,fileSuffix,".png",sep=""),
                              plot = g, width = width, height = height, dpi = dpi))
    }
    
    # Save data to an object
    respCurvesList[["DATA"]][[modVar]] <- respCurves
    respCurvesList[["PLOTS"]][[modVar]] <- g
    
    setTxtProgressBar(pb, k)
  }
  invisible(respCurvesList)
}


#' Plot variable importance from biomod2
#' 
#' Plot variable average importance scores across different PA sets, algorithms and 
#' evaluation rounds. \code{ggplot2} graphics are used to make the plots.
#' 
#' @export
#' 

varImportancePlot <- function(biomodModelsOut, by="all", sortVars=TRUE, 
                              filterAlgos = NULL, plotType="facet_wrap",
                              save=FALSE, plot=TRUE, outputFolder=getwd(), 
                              filename="VarImportancePlot.png",...){ 
  
  # Calculate variable importance
  varImportance <- get_variables_importance(biomodModelsOut)
  
  # Make plot for each variable and across all PA sets, Algos and eval rounds
  if(by=="all"){
    
    # Calculate average and std-error stats
    varImportanceByVariableAVG <- apply(varImportance, 1, mean)
    varImportanceByVariableSTE <- apply(varImportance, 1, function(x) sd(x)/sqrt(length(x)))
    
    vimpDF <- data.frame(cnames = names(varImportanceByVariableAVG),
                         varImpAVG = varImportanceByVariableAVG, 
                         varImpSTE = varImportanceByVariableSTE) %>% 
      arrange(desc(varImpAVG))
    
    if(sortVars){
      vimpDF[,"cnames"] <- factor(vimpDF$cnames, levels = vimpDF$cnames, 
                                  labels = vimpDF$cnames) 
    }
    
    g <- ggplot(vimpDF,aes(x= cnames, y = varImpAVG)) + 
      geom_bar(stat="identity") + 
      geom_errorbar(aes(ymin=varImpAVG - 2*varImpSTE, ymax=varImpAVG + 2*varImpSTE), 
                    width=0.35) + 
      xlab("Variables") + 
      ylab("Variable importance score") + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  # Make a panel wrap or grid by model algorithm
  else if(by=="algo"){
    
    # Calculate average and std-error stats
    # The c(1,2) means average for each variable (row) and model algo (col)
    varImportanceByVariableAVG <- apply(varImportance, c(1,2), mean) %>% 
      as.data.frame %>% 
      mutate(varNames = as.factor(rownames(.))) %>% 
      gather("modAlgo","varImpAVG",setdiff(colnames(.),"varNames"))
    
    varImportanceByVariableSTE <- apply(varImportance, c(1,2), 
                                        function(x) sd(x)/sqrt(length(x))) %>% 
      as.data.frame %>% 
      mutate(varNames = as.factor(rownames(.))) %>% 
      gather("modAlgo","varImpSTE",setdiff(colnames(.),"varNames"))
    
    # Join data
    vimpDF <- cbind(varImportanceByVariableAVG, 
                    varImpSTE=varImportanceByVariableSTE[,"varImpSTE"])
    
    # Filter rows for specific model algos only
    if(!is.null(filterAlgos)){
      vimpDF <- vimpDF %>% filter(modAlgo %in% filterAlgos)
    }
    
    # Make plot
    g <- ggplot(vimpDF,aes(x= varNames, y = varImpAVG)) + 
      geom_bar(stat="identity") + 
      geom_errorbar(aes(ymin=varImpAVG - 2*varImpSTE, ymax=varImpAVG + 2*varImpSTE), 
                    width=0.35) + 
      xlab("Variables") + 
      ylab("Variable importance score") + 
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) 
      
      # Add a panel grid or wrap
      if(plotType=="facet_wrap"){
        g <- g + facet_wrap(modAlgo~.)
      }else if(plotType=="facet_grid"){
        g <- g + facet_grid(modAlgo~.)
      }else{
        stop("Invalid option in plotType!")
      }
  }
  
  if(plot){
    plot(g)
  }
  
  if(save){
    suppressMessages(ggsave(filename = paste(outputFolder,"/",filename,sep=""), 
                            plot = g, ...))
  }
  
  invisible(list(vimpAll=varImportance, vimpSummary=vimpDF, gplot=g))
  
}

#' A plot for evaluation metrics
#' 
#' Makes a plot using ggplot2 for one selected evaluation metric for each 
#' modelling algorithm
#' 
#' @export
#' 

evalMetricPlot <- function(biomodModelOut, evalMetric = "TSS", sort = TRUE, 
                           removeFull=FALSE){
  
  myBiomodModelEval <- biomod2::get_evaluations(biomodModelOut)
  
  evalDF <- as.data.frame(myBiomodModelEval[evalMetric,"Testing.data",,,])
  
  if(removeFull){
    evalDF <- select(evalDF, -starts_with("FULL"))
  }
  
  evalDFsummary <- data.frame(
    algoName = rownames(evalDF),
    evalScoreAVG = apply(evalDF, 1, mean),
    evalScoreSTE = apply(evalDF, 1, function(x) sd(x)/sqrt(length(x)))) %>% 
    arrange(desc(evalScoreAVG))
  
  if(sort){
    evalDFsummary[,"algoName"] <- factor(evalDFsummary$algoName, evalDFsummary$algoName,
                                         labels = evalDFsummary$algoName)
  }
  
  g <- ggplot(evalDFsummary, aes(y = evalScoreAVG, x = algoName)) + 
    geom_bar(stat="identity") + 
    geom_errorbar(aes(ymin = evalScoreAVG - 2*evalScoreSTE, 
                      ymax = evalScoreAVG + 2*evalScoreSTE), width=0.35) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1)) + 
    ylab(paste(evalMetric,"score")) + 
    xlab("Variables")
  
  plot(g)
  
  invisible(list(biomodEvalObj = myBiomodModelEval,
                 evalDataFrame = evalDF,
                 evalSummary   = evalDFsummary,
                 ggPlotObj     = g))
}


