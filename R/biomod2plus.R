
#' Get distinct raster cell centroids from input presence records
#' 
#' This function is used to remove duplicates and NA's in the input points and 
#' return the coordinates for each unique cell centroid
#' 
#' @param rasterObj A raster stack or layer object
#' @param pointObj A set of points as a SpatialPoints* object
#' @param spatial Return as a spatial object? (default: TRUE) 
#' 
#' @return Either a \code{data.frame} with unique cell centroid coordinates for the input 
#' point dataset or a \code{SpatialPoints} object
#' 
#' @importFrom raster extract
#' @importFrom raster xyFromCell
#' @importFrom raster crs
#' @importFrom dplyr mutate
#' @importFrom dplyr distinct
#' @importFrom dplyr select
#' @importFrom sp SpatialPoints
#' 
#' @export
#' 

getUniqueRasterXYCoords <- function(rasterObj, pointObj, spatial=TRUE){
  
  if(!inherits(rasterObj,c("RasterLayer","RasterStack"))){
    stop("rasterObj must be a RasterLayer or a RasterStack object!")
  }
  
  # Extract data to points and remove NA's
  xyData <- raster::extract(rasterObj, pointObj, cellnumbers=TRUE) %>% 
    na.omit
  
  # Get unique xy coordinates from raster cell centroids
  xyCoordsDF <- raster::xyFromCell(rasterObj, xyData[,"cells"]) %>% 
    as.data.frame %>% 
    dplyr::mutate(ID = paste(x, y, sep = "_")) %>% 
    dplyr::distinct(ID, .keep_all = TRUE) %>% 
    dplyr::select(x, y)
  
  if(!spatial){
    return(xyCoordsDF)  
  }else{
    spObj <- sp::SpatialPoints(xyCoordsDF, proj4string = raster::crs(pointObj))
    return(spObj)
  }
}



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
#' @importFrom raster stack
#' @importFrom raster writeRaster 
#' @importFrom tools file_path_sans_ext
#' 
#' @export

convertToGeoTIFF <- function(inputFolder,outputFolder){
  
  fl <- list.files(inputFolder, pattern= ".grd$", full.names = TRUE)
  
  if(!dir.exists(outputFolder)){
    dir.create(outputFolder)
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
#' 
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
#' @importFrom stringr str_to_title
#' 
#' @export
#' 

abbrevNames <- function(x) paste(stringr::str_to_title(unlist(strsplit(gsub("\\.","",x),"\\ "))),
                                 collapse="")


#' Project multiple scenarios
#' 
#' A wrapper function used to project and ensemble multiple scenarios contained in a 
#' list object with one \code{RasterStack} per scenario. Useful for forecasting (for example) 
#' multiple climate change RCP scenarios or future dates (or both).
#' 
#' @param myBiomodModelOut A \link[biomod2]{BIOMOD.models.out-class} class object.
#' @param myBiomodProj A \link[biomod2]{BIOMOD.projection.out-class} class object.
#' @param myBiomodEM A \link[biomod2]{BIOMOD.EnsembleModeling.out-class} class object
#' @param modelsToUse A character vector listing the models to use.
#' @param scenarioNames The character vector with names for each different scenario. 
#' Must equal the length of \code{rstStackList}.
#' @param rstStackList A list object with one \code{RasterStack} object by scenario. 
#' The list may be named and these names will be used as the input of \code{proj.name} 
#' parameter from \link[biomod2]{BIOMOD_Projection} function.
#' @param convert2GeoTIFF Convert to GeoTIFF file format?
#' 
#' @return NULL
#' 
#' @seealso 
#' \link[biomod2]{BIOMOD_Projection}, \link[biomod2]{BIOMOD_EnsembleForecasting}, 
#' \link{convertToGeoTIFF}
#' 
#' @importFrom biomod2 BIOMOD_Projection
#' @importFrom biomod2 BIOMOD_EnsembleForecasting
#' 
#' @export
#'  

projectMultipleScenarios <- function(myBiomodModelOut, myBiomodProj, myBiomodEM, 
                                       modelsToUse, scenarioNames = NULL, rstStackList, 
                                       convert2GeoTIFF = TRUE){
  
  if(is.null(scenarioNames)){
    scenarioNames <- names(rstStackList)
    if(is.null(scenarioNames)){
      stop("scenarioNames is not properly defined!")
    }
  }
  
  sp <- myBiomodModelOut@sp.name
  
  for(i in 1:length(rstStackList)){
    
    # Obtain spatiotemporal projections
    myBiomodProj <- biomod2::BIOMOD_Projection(modeling.output = myBiomodModelOut,
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
    myBiomodEF <- biomod2::BIOMOD_EnsembleForecasting(projection.output = myBiomodProj,
                                             binary.meth = c('TSS','ROC','KAPPA'),
                                             EM.output = myBiomodEM,
                                             output.format = '.grd')
    
    # Convert all output raster files to GeoTIFF
    inFolder <- paste(getwd(), "/", sp, "/proj_", scenarioNames[i], sep="")
    outFolder <- paste(inFolder, "/", "GeoTIFF", sep="")
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
#' @param biomodModelOut A \link[biomod2]{BIOMOD.models.out-class} class object.
#' @param qt The target quantile between 0 and 1.
#' @param evalMetric A character string defining the evaluation metric. 
#' Available options are: 'KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 
#' 'POD', 'CSI' and 'ETS'.
#' @param na.rm Remove NA's?
#' 
#' @return The target quantile value.
#' 
#' @importFrom biomod2 get_evaluations
#' @importFrom stats quantile
#' 
#' @export
#' 

getEvalMetricQuantile <- function(biomodModelOut, qt = 0.75, evalMetric = "TSS", na.rm=TRUE){
  
  if(length(qt)>1){
    qt <- qt[1]
    warning("qt parameter length is greater than 1. Using only the first element!")
  }
  
  myBiomodModelEval <- biomod2::get_evaluations(biomodModelOut)
  evalDF <- as.data.frame(myBiomodModelEval[evalMetric,"Testing.data",,,])
  quantileThresh <- stats::quantile(evalDF, probs = qt, na.rm = na.rm)
  
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
#' 
#' @importFrom raster getValues
#' @importFrom raster nlayers
#' @importFrom raster setValues
#' @importFrom raster subset
#' @importFrom raster stack
#' @importFrom raster mask
#' 
#' @export
#' 

intersectMask <- function(x){
  
  values_x <- raster::getValues(x)
  inter_x <- values_x %*% rep(1,raster::nlayers(x))
  mask <- raster::setValues(raster::subset(x,1), values = (inter_x>0))
  
  return(raster::stack(raster::mask(x, intersectMask(x))))
}


#' Check which model algorithms where run
#' 
#' A simple function used to verify which models were actually run by \code{biomod2} 
#' for a given target species.
#' 
#' @param biomodModelOut A \link[biomod2]{BIOMOD.models.out-class} class object.
#' 
#' @return A named logical vector specifying which models were run.
#' 
#' @importFrom biomod2 get_built_models
#' 
#' @export
#' 

checkModAlgoRun <- function(biomodModelOut){
  
  allModAlgos <- c("GLM", "GBM", "GAM", "CTA", "ANN", "FDA", "MARS", 
                   "RF", "MAXENT.Phillips", "MAXENT.Tsuruoka")
  builtMods <- biomod2::get_built_models(biomodModelOut)
  wereBuilt <- vector(mode = "logical", length = length(allModAlgos))
  
  for(i in 1:length(allModAlgos)){
    wereBuilt[i] <- ifelse(sum(grepl(allModAlgos[i],builtMods)) > 0,TRUE,FALSE)
  }
  names(wereBuilt) <- allModAlgos
  return(wereBuilt)
}


#' Adjust the percentage of the training sample
#' 
#' An ancillary function used to modify the training sample percentage so that the 
#' total amount of training records (presences + pseudo absences) times the train 
#' fraction) is below maxTrain. The function is used in cases where model fitting is 
#' very slow due to the amount of available records.
#' 
#' @param biomodData A \link[biomod2]{BIOMOD.formated.data-class} object
#' @param percTrain The percentage of training data (0-100)
#' @param maxTrain The maximum number of training points (default: 1000)
#' 
#' @return An adjusted training set percentage to guarantee that the total 
#' number of training points does not exceed the maxTrain value.
#' 
#' @export
#' 

adjustPercTrain <- function(biomodData, percTrain, maxTrain=1000){
  
  # Total number of records
  nTotal <- sum(biomodData@PA[,1])
  
  nTrainTotal <- round(nTotal * (percTrain/100))
  if(nTrainTotal > maxTrain){
    percTrainAdj <- round((maxTrain / nTotal)*100)
    return(percTrainAdj)
  }else{
    return(percTrain)
  }
} 

#' Calculate the number of pseudo-absences using a power function
#' 
#' This function is intended to calculate the number of pseudo-absences based on a 
#' power function thus allowing to have a large number of PA's for small sample 
#' sizes and, it approximates the number presences for very large sample sizes. 
#' The function is defined  as follows: \eqn{n_{PA} = n_{P} \times a \times n_{P}^{-b}}, 
#' with \eqn{n_{PA}} equal to the number of pseudo-absences and \eqn{n_{P}} to the number of 
#' presences. The parameters \eqn{a} and \eqn{b} allow to adjust the amount of PA's.
#' 
#' @param nPresences Number of presence records
#' @param a Multiplication parameter (default: 45)
#' @param b Exponential parameter (default: -0.55)
#'
#' @return An integer defining the number of pseudo-absences
#'
#' @export
#'

powerFunPAnumberCalc <- function(nPresences, a = 45, b = -0.55){
  return(round(nPresences * (a * nPresences^b)))
}


#' Performs a two-stage best model selection
#' 
#' The function is used to make a two-stage selection based, first selecting the best 
#' modelling algorithms and then selecting the top best-performing fraction of models 
#' within those techniques. 
#'
#' @param biomodModelOut A \link[biomod2]{BIOMOD.models.out-class} class object
#' @param evalMetric The evaluation metric to use (available options are: 
#' 'KAPPA', TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI' and 'ETS')
#' @param nrBestAlgos Integer number defining how many top-performing 
#' algorithms (default: 5)
#' @param bestAlgoFun Function used to aggregate model scores and rank each 
#' algorithm (default: median)
#' @param topFraction Fraction of best-performing models for each selected modelling 
#' technique 
#' @param na.rm Remove NA's?
#'
#' @return A character vector with selected models
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr select
#' @importFrom dplyr summarize
#' @importFrom dplyr filter
#' @importFrom tidyr gather
#' @importFrom biomod2 get_evaluations
#' @importFrom stats median
#' 
#' @export
#' 

twoStepBestModelSelection <- function(biomodModelOut, 
                                      evalMetric = "TSS", 
                                      nrBestAlgos = 5, 
                                      bestAlgoFun = stats::median, 
                                      topFraction = 0.25,
                                      na.rm = TRUE){
  
  ## STEP #1 - select the top-performing techniques ------------------------
  ##
  ##
  
  # Get model evaluations
  modEvalArray <- biomod2::get_evaluations(biomodModelOut)
  
  modEvalDF <- as.data.frame(modEvalArray[evalMetric, "Testing.data",,,])
  ncols <- ncol(modEvalDF)
  
  modAlgoRank <- sort(apply(modEvalDF, MARGIN = 1, FUN = bestAlgoFun, na.rm = na.rm), decreasing = TRUE)
  bestAlgos <- names(modAlgoRank)[1:nrBestAlgos]
  
  
  ## STEP #2 - select top-performing models --------------------------------------- 
  ##
  
  modEvalDF2 <- tidyr::gather(modEvalDF, key = "modString", value = "modScore")
  
  # Create the full evaluation data frame
  #
  modEvalDF <- modEvalDF %>% 
    rownames %>% 
    rep(ncols) %>% # Generate a char vector with model names
    paste(modEvalDF2$modString %>% 
            strsplit("\\.") %>% 
            lapply(FUN = function(x) return(x[c(2,1)])) %>% # Invert the order to PA set and run
            lapply(FUN = function(x) paste(x, collapse="_",sep="")) %>% 
            unlist(use.names = FALSE), sep="_") %>% 
    strsplit("_") %>% 
    unlist %>%  
    matrix(ncol = 3, byrow = TRUE) %>%
    data.frame(modEvalScore = modEvalDF2$modScore, stringsAsFactors = FALSE) %>% 
    `colnames<-`(c("modAlgo","PAset","evalRun","modEvalScore"))

  # Check if model evaluation score array only has one PA set
  if(dim(modEvalArray)[5]==1){
    modEvalDF[,"PAset"] <- "PA1"
  }
  
  # if(all(is.na(modEvalDF[,"PAset"])) | all(modEvalDF[,"PAset"]=="NA")){
  #   modEvalDF[,"PAset"] <- "PA1"
  # }
  
  # An ugly hack to correct for model calibrations with a single evaluation run!!!...
  if((dim(modEvalArray)[4]==1)){
    modEvalDF[,"PAset"] <- modEvalDF[,"evalRun"]
    modEvalDF[,"evalRun"] <- "RUN1"
  }

  # Calculate the top fraction quantile values for the best modelling algorithms
  # These will be used as cutoff values to select the best models
  #
  topQtValues <- modEvalDF %>% 
    dplyr::filter(modAlgo %in% bestAlgos) %>% 
    dplyr::group_by(modAlgo) %>% 
    dplyr::summarize(qts = stats::quantile(modEvalScore, probs = 1 - topFraction, na.rm = na.rm))
  
  i <- 0
  selMods <- c()
  
  for(modAlgo_i in as.character(topQtValues$modAlgo)){
    i <- i + 1
    # Make and append model character strings in biomod2 style
    selMods <- c(selMods,
                 modEvalDF %>% 
                   dplyr::filter(modAlgo == modAlgo_i, modEvalScore >= topQtValues$qts[i]) %>% 
                   dplyr::select(2,3,1) %>% 
                   apply(MARGIN = 1, FUN = function(x) paste(x, collapse="_",sep="")))
  }
  
  # Append the target species name to the model id strings
  selMods <- paste(biomodModelOut@sp.name,selMods,sep="_")
  
  return(selMods)
}

## PLOTS AND GRAPHICS ----------------------------------------------------------------------

#' Make response plots for multiple variables 
#' 
#' A flexible function used to produce response plots using \code{ggplot2} package. 
#' 
#' @param biomodModelOut A \link[biomod2]{BIOMOD.models.out-class} class object
#' @param Data A \code{data.frame} or \code{RasterStack} containing the data used to make 
#' the response plots. The data must have the same names as the ones used to fit the model.
#' @param modelsToUse A character vector specifying which models to use. The 
#' \link[biomod2]{BIOMOD_LoadModels} function must be used for this purpose allowing to load 
#' the models into the Global Environment so these can be used by the \link[biomod2]{response.plot2} 
#' function
#' @param showVars A character vector specifying which variables should be used. The option "all" 
#' (default) uses all the variables contained in \code{Data}
#' @param fixedVarMetric Either 'mean' (default), 'median', 'min' or 'max' specifying the 
#' statistic used to fix as constant the remaining variables when the predicted response 
#' is estimated for one of the variables
#' @param addMarginalPlot Logical determining if a marginal plot, showing the distribution of 
#' the target variable, is plotted (default: FALSE)
#' @param marginPlotType Type of marginal plot to add. Valid options are: "histogram" (default), 
#' "boxplot", "density", "violin" or "densigram"
#' @param varNames Character vector with the names for each variable to be plotted
#' @param plotStdErr Plot with std-error bands or simply plot all response lines, one 
#' for each model selected (default: FALSE)
#' @param plot Plot to device?
#' @param save Save plot? (by default the png image format is used)
#' @param filePrefix Filename prefix (default: "ResponsePlot_")
#' @param fileSuffix Filename suffix (default: "")
#' @param outFolder Output folder where the plot image will be placed (default: getwd())
#' @param na.rm Remove NA's?
#' @param ... Other parameters passed to \link[ggplot2]{ggsave}
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
#' To use this function, first the selected models have to loaded to the Global Environment 
#' using the \link[biomod2]{BIOMOD_LoadModels} function. 
#' 
#' @seealso 
#' \link[biomod2]{response.plot2}   
#' \link[biomod2]{BIOMOD_LoadModels}
#' 
#' @importFrom raster values
#' @importFrom biomod2 response.plot2
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_ribbon
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 aes_string
#' @importFrom ggExtra ggMarginal
#' @importFrom stats sd
#' @importFrom stats na.omit
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' 
#' @export
#' 

responsePlots <- function(biomodModelOut, Data, modelsToUse, showVars = "all", 
                          fixedVarMetric = 'mean', addMarginalPlot = TRUE, 
                          marginPlotType = "histogram", varNames = NULL,  
                          plotStdErr = FALSE, plot = FALSE, save = TRUE, 
                          filePrefix = "ResponsePlot_", 
                          fileSuffix = "", outFolder = getwd(), na.rm = TRUE, ...){
  
  if(inherits(Data,"RasterStack")){
    cat("Loading data from RasterStack....")
    Data <- as.data.frame(stats::na.omit(raster::values(Data)))
    cat("done.\n\n")
  }
  
  if(!inherits(Data,"data.frame")){
    stop("Data must be set as a data.frame object!")
  }
  
  if(showVars=="all"){
    showVars <- colnames(Data)
  }
  
  stdErr <- function(x, na.rm = TRUE){
    if(na.rm){
      stats::sd(x, na.rm = na.rm)/sqrt(sum(!is.na(x)))
    }else{
      stats::sd(x)/sqrt(length(x))
    }
  }
  
  respCurvesList <- list()
  respCurvesList[["DATA"]] <- list()
  respCurvesList[["PLOTS"]] <- list()
  
  k <- 0
  pb <- utils::txtProgressBar(1,length(showVars),style=3)
  
  for(modVar in showVars){
    
    k<-k+1
    respCurves <- as.data.frame(biomod2::response.plot2(models = modelsToUse,
                                               Data = Data,
                                               show.variables = modVar,
                                               do.bivariate = FALSE,
                                               fixed.var.metric = fixedVarMetric,
                                               plot = FALSE,
                                               save.file = "no"))
    
    ## Make one plot per variable ------------------------------------------------- ##
    ##
    ##
    tmpDF <- data.frame(targetVar = respCurves[,1], 
                        avgResp = apply(respCurves[,-1], 1, mean, na.rm = na.rm),
                        stdErr = apply(respCurves[,-1], 1, stdErr, na.rm = na.rm))
    
    
    cn <- colnames(respCurves)
    
    ## Create plot with response lines per PA/eval_round
    if(!plotStdErr){
      
      # Start the ggplot object
      g <- ggplot2::ggplot(respCurves, ggplot2::aes_string(x = cn[1]))
      
      # Add one response line per model to the ggplot object
      for(i in 2:(ncol(respCurves))){
        g <- g + ggplot2::geom_line(ggplot2::aes_string(y=cn[i]), alpha=0.2)
      }
      
      # Add the average line and change default theme settings
      g <- g + 
            ggplot2::geom_line(mapping = aes(y = avgResp, x = targetVar),data = tmpDF, size = 1.5) + 
            ggplot2::xlab(ifelse(is.null(varNames),modVar,varNames[k])) + 
            ggplot2::ylab("Response") + 
            ggplot2::theme_bw()
      
      # Create plot with average and std-error curves/bands  
    }else{
      g <- ggplot2::ggplot(tmpDF,aes(x = targetVar)) + 
            ggplot2::geom_ribbon(aes(ymin = avgResp - 2*stdErr, ymax = avgResp + 2*stdErr), 
                        alpha = 0.3, fill = "dark blue") + 
            ggplot2::geom_line(aes(y = avgResp)) + 
            ggplot2::xlab(modVar) + 
            ggplot2::ylab("Response") + 
            ggplot2::theme_bw()
    }
    
    # Adds a marginal plot to show the distribution of the target variable
    # It can be either a boxplot, density, histogram, ...
    if(addMarginalPlot){
      
      minVal <- min(tmpDF$avgResp)
      g <- g + ggplot2::geom_point(mapping = ggplot2::aes_string(y = minVal, x = modVar), 
                          data = Data, alpha = 0)
      g <- ggExtra::ggMarginal(g, margins = "x", type = marginPlotType)
    }
    
    # Plot response curve?
    if(plot){
      plot(g)
    } 
    
    # Save image per variable?
    if(save){
      suppressMessages(ggplot2::ggsave(filename = paste(outFolder,"/",filePrefix,modVar,fileSuffix,".png",sep=""),
                              plot = g, ...))
    }
    
    # Save data to an object
    respCurvesList[["DATA"]][[modVar]] <- respCurves
    respCurvesList[["PLOTS"]][[modVar]] <- g
    
    utils::setTxtProgressBar(pb, k)
  }
  invisible(respCurvesList)
}


#' Plot variable importance from biomod2
#' 
#' Plot variable average importance scores across different PA sets, algorithms and 
#' evaluation rounds. \code{ggplot2} graphics are used to make the plots.
#' 
#' @param biomodModelsOut A \link[biomod2]{BIOMOD.models.out-class} class object.
#' @param by Options are "all" (default) for plotting the average across all PA sets, algorithms 
#' and evaluation rounds or "algo" to plot the average importance scores by modelling 
#' algorithm and across PA sets and evaluation rounds
#' @param sortVars Sort variables (in decreasing order) by importance score? (default: TRUE)
#' @param filterAlgos Character vector used to select specific algorithms. OPtions are: 
#' 'GLM', 'GBM', 'GAM', 'CTA', 'ANN', 'SRE', 'FDA', 'MARS', 'RF', 'MAXENT.Phillips' and 
#' "MAXENT.Tsuruoka" (default: NULL)
#' @param plotType Plot type. Available options are "facet_wrap" or "facet_grid"
#' @param save Save plot? (default: FALSE)
#' @param plot Plot to device? (default: TRUE)
#' @param outputFolder Output folder used to save the plot image (uses the 
#' \link[ggplot2]{ggsave}) function
#' @param filename Name of the image file to save
#' @param na.rm Remove NA's?
#' @param ... Other parameters passed to \link[ggplot2]{ggsave}) function
#' 
#' @return A list object containing:
#' \itemize{
#'    \item vimpAll - All importance scores array object
#'    \item vimpSummary - Variable importance summary
#'    \item gplot - The ggplot object
#' }
#' 
#' @importFrom biomod2 get_variables_importance
#' @importFrom dplyr arrange 
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr desc
#' @importFrom magrittr %>%
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggsave
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 facet_grid
#' @importFrom stats sd
#' 
#' @export
#' 

varImportancePlot <- function(biomodModelsOut, by = "all", sortVars = TRUE, 
                              filterAlgos = NULL, plotType = "facet_wrap",
                              save = FALSE, plot = TRUE, outputFolder = getwd(), 
                              filename = "VarImportancePlot.png", na.rm = TRUE, ...){ 
  
  # Calculate variable importance
  varImportance <- biomod2::get_variables_importance(biomodModelsOut)
  
  stdErr <- function(x, na.rm = TRUE){
    if(na.rm){
      stats::sd(x, na.rm = na.rm)/sqrt(sum(!is.na(x)))
    }else{
      stats::sd(x)/sqrt(length(x))
    }
  }
  
  # Make plot for each variable and across all PA sets, Algos and eval rounds
  if(by=="all"){

    # Calculate average and std-error stats
    varImportanceByVariableAVG <- apply(varImportance, 1, mean, na.rm = na.rm)
    varImportanceByVariableSTE <- apply(varImportance, 1, stdErr, na.rm = na.rm)
    
    vimpDF <- data.frame(cnames = names(varImportanceByVariableAVG),
                         varImpAVG = varImportanceByVariableAVG, 
                         varImpSTE = varImportanceByVariableSTE) %>% 
      dplyr::arrange(dplyr::desc(varImpAVG))
    
    if(sortVars){
      vimpDF[,"cnames"] <- factor(vimpDF$cnames, levels = vimpDF$cnames, 
                                  labels = vimpDF$cnames) 
    }
    
    g <- ggplot2::ggplot(vimpDF,ggplot2::aes(x= cnames, y = varImpAVG)) + 
          ggplot2::geom_bar(stat="identity") + 
          ggplot2::geom_errorbar(aes(ymin=varImpAVG - 2*varImpSTE, ymax=varImpAVG + 2*varImpSTE), width=0.35) + 
          ggplot2::xlab("Variables") + 
          ggplot2::ylab("Variable importance score") + 
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  # Make a panel wrap or grid by model algorithm
  else if(by=="algo"){
    
    # Calculate average and std-error stats
    # The c(1,2) means average for each variable (row) and model algo (col)
    varImportanceByVariableAVG <- apply(varImportance, c(1,2), mean, na.rm = na.rm) %>% 
      as.data.frame %>% 
      dplyr::mutate(varNames = as.factor(rownames(.))) %>% 
      tidyr::gather("modAlgo","varImpAVG",setdiff(colnames(.),"varNames"))
    
    varImportanceByVariableSTE <- apply(varImportance, c(1,2), stdErr, na.rm = na.rm) %>% 
      as.data.frame %>% 
      dplyr::mutate(varNames = as.factor(rownames(.))) %>% 
      tidyr::gather("modAlgo","varImpSTE",setdiff(colnames(.),"varNames"))
    
    # Join data with average values and std-errors
    vimpDF <- cbind(varImportanceByVariableAVG, 
                    varImpSTE = varImportanceByVariableSTE[,"varImpSTE"])
    
    # Filter rows for specific model algos only
    if(!is.null(filterAlgos)){
      vimpDF <- vimpDF %>% dplyr::filter(modAlgo %in% filterAlgos)
    }
    
    # Make ggplot object with var importance scores
    g <- ggplot2::ggplot(vimpDF,ggplot2::aes(x = varNames, y = varImpAVG)) + 
          ggplot2::geom_bar(stat="identity") + 
          ggplot2::geom_errorbar(aes(ymin = varImpAVG - 2*varImpSTE, 
                                     ymax = varImpAVG + 2*varImpSTE), 
                                     width = 0.35) + 
          ggplot2::xlab("Variables") + 
          ggplot2::ylab("Variable importance score") + 
          ggplot2::theme_bw() +
          ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust=0.5)) 
      
      # Add a panel grid or wrap
      if(plotType=="facet_wrap"){
        g <- g + ggplot2::facet_wrap(modAlgo~.)
      }else if(plotType=="facet_grid"){
        g <- g + ggplot2::facet_grid(modAlgo~.)
      }else{
        stop("Invalid option in plotType!")
      }
  }
  
  # Plot to device?
  if(plot){
    plot(g)
  }
  
  # Save plot to image file?
  if(save){
    suppressMessages(ggplot2::ggsave(filename = paste(outputFolder,"/",filename,sep=""), 
                            plot = g, ...))
  }
  
  invisible(list(vimpAll=varImportance, vimpSummary=vimpDF, gplot=g))
  
}

#' Plot evaluation metrics calculated by biomod2
#' 
#' Makes a plot using ggplot2 for one selected evaluation metric for each 
#' modelling algorithm
#' 
#' @param biomodModelOut A \link[biomod2]{BIOMOD.models.out-class} class object
#' @param evalMetric String defining the evaluation metric (e.g., 'ROC', 'TSS', 'KAPPA')
#' @param sort Sort the bar plot by the average evaluation score (best to worst, default: TRUE)
#' @param removeFull Remove full evaluation rounds? (default: FALSE)
#' @param save Save the plot? (default: FALSE)
#' @param plot Plot to device? (default: TRUE)
#' @param outputFolder Output folder where the plot will be save to (default: \code{getwd()})
#' @param filename Filename of the image containing the plot (default:"evalPlot.png")
#' @param na.rm Remove NA's?
#' @param ... Other arguments passed to \code{ggsave} function such as height, width, dpi, etc.
#' 
#' @return A list object with the following components:
#' \itemize{
#'    \item biomodEvalObj - A full array with all biomod2 evaluation metrics
#'    \item evalDataFrame - A data.frame with evaluation scores for the target metric
#'    \item evalSummary - A summary data.frame with average and std-error values for 
#'    the target metric
#'    \item ggPlotObj - The ggplot2 object used for plotting the evaluation scores
#' 
#' }
#' 
#' @importFrom biomod2 get_evaluations
#' @importFrom dplyr select
#' @importFrom dplyr arrange
#' @importFrom dplyr desc
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 geom_errorbar
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 ggsave
#' @importFrom tidyselect starts_with
#' @importFrom stats sd
#' 
#' @export
#' 

evalMetricPlot <- function(biomodModelOut, evalMetric = "TSS", sort = TRUE, 
                           removeFull=FALSE, save=FALSE, plot=TRUE, outputFolder=getwd(),
                           filename="evalPlot.png", na.rm = TRUE, ...){
  
  myBiomodModelEval <- biomod2::get_evaluations(biomodModelOut)
  
  evalDF <- as.data.frame(myBiomodModelEval[evalMetric,"Testing.data",,,])
  
  if(removeFull){
    evalDF <- dplyr::select(evalDF, -tidyselect::starts_with("FULL"))
  }
  
  stdErr <- function(x, na.rm = TRUE){
    if(na.rm){
      stats::sd(x, na.rm = na.rm)/sqrt(sum(!is.na(x)))
    }else{
      stats::sd(x)/sqrt(length(x))
    }
  } 
  
  evalDFsummary <- data.frame(
    algoName = rownames(evalDF),
    evalScoreAVG = apply(evalDF, 1, mean, na.rm = na.rm),
    evalScoreSTE = apply(evalDF, 1, stdErr, na.rm = na.rm)) %>% 
    dplyr::arrange(dplyr::desc(evalScoreAVG))
  
  if(sort){
    evalDFsummary[,"algoName"] <- factor(evalDFsummary$algoName, evalDFsummary$algoName,
                                         labels = evalDFsummary$algoName)
  }
  
  g <- ggplot2::ggplot(evalDFsummary, ggplot2::aes(y = evalScoreAVG, x = algoName)) + 
        ggplot2::geom_bar(stat="identity") + 
        ggplot2::geom_errorbar(ggplot2::aes(ymin = evalScoreAVG - 2*evalScoreSTE, 
                                   ymax = evalScoreAVG + 2*evalScoreSTE), width=0.35) +
        ggplot2::theme_bw() + 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust=1)) + 
        ggplot2::ylab(paste(evalMetric,"score")) + 
        ggplot2::xlab("Variables")
  
  if(plot){
    plot(g)
  }
  
  if(save){
    suppressMessages(ggplot2::ggsave(filename = paste(outputFolder,"/",filename,sep=""), 
                            plot = g, ...))
  }
  
  invisible(list(biomodEvalObj = myBiomodModelEval,
                 evalDataFrame = evalDF,
                 evalSummary   = evalDFsummary,
                 ggPlotObj     = g))
}


#' Check raster stacks for Maxent
#'
#' Checks a \code{RasterStack} object to avoid problems during Maxent (Phillips) projection step
#' 
#' @param x A \code{RasterStack} object to check
#'
#' @return A new and checked \code{RasterStack}
#'
#' @note This is based on Damien Georges suggestion at: 
#' \url{https://r-forge.r-project.org/forum/message.php?msg_id=42888&group_id=302}
#'
checkRasterStack <- function(x){
  
  ## function to define the intersect of rasters
  intersect_mask <- function(x){
    values_x <- getValues(x)
    inter_x <- values_x %*% rep(1,nlayers(x))
    mask <- setValues(subset(x,1),values = (inter_x>0))
    return(mask)
  }
  
  out <- stack(mask(x, intersect_mask(x)))
  names(out) <- names(x)
  return(out)
}


