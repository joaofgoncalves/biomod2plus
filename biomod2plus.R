
#' Convert raster files to GeoTIFF

convertToGeoTIFF <- function(inputFolder,outputFolder){
  
  fl <- list.files(inputFolder, pattern= ".grd$", full.names = TRUE)
  
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


projectEnsemblesMultiStack <- function(myBiomodModelOut, myBiomodProj, myBiomodEM, 
                                       modelsToUse, projNames, rstStackList){
  
  
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



