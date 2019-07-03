#' Optimization of GR2M model parameters with MOPSOCD algorithm.
#'
#' @param Location		 General work directory where data is located.
#' @param Qmodel       Simulated streamflow matrix (time, subbasin) from Run_GR2MSemiDistr
#' @param Shapefile		 Subbasins shapefile.
#' @param FlowDir			 Flow direction raster in GRASS format.
#' @param Mask         Subbasins centroids mask raster. 'Centroids_mask.tif' as default.
#' @return  Routing streamflow for each subbasin.
#' @export
#' @import  rgrass7
#' @import  rgdal
#' @import  raster
#' @import  foreach
#' @import  tictoc
Routing_GR2MSemiDistr <- function(Location, Qmodel, Shapefile, FlowDir, Mask='Centroids_mask.tif'){

  require(rgrass7)
  require(foreach)
  require(rgdal)
  require(raster)
  require(tictoc)
  tic()

  # Load rasters
  path.fdr  <- file.path(Location,'Inputs', FlowDir)
  path.mask <- file.path(Location,'Inputs', Mask)
  flowdir   <- raster(path.fdr)
  qMask     <- raster(path.mask)

  # Load shapefile
  path.shp   <- file.path(Location,'Inputs', Shapefile)
  basin      <- readOGR(path.shp, verbose=F)

  # Load GRASS (require to be installed previously)
  loc <- initGRASS('C:/Program Files/GRASS GIS 7.4.4', home=getwd(), gisDbase="GRASS_TEMP", override=TRUE)

  # Import raster to GRASS
  writeRAST(as(flowdir, 'SpatialGridDataFrame'), "fdr", overwrite=T)

  # Set an study extention
  execGRASS("g.region", Sys_show.output.on.console=F, parameters=list(raster="fdr"))

  # Auxiliary variables
  qRas   <- qMask
  qSub   <- matrix(NA, nrow=nrow(Qmodel), ncol=ncol(Qmodel))          # Streamflow time series
  qArray <- array(NA, dim=c(nrow(qMask), ncol(qMask), nrow(Qmodel)))  # Streamflow matrix for each time-step
  qBrick <- brick(nr=nrow(qMask), nc=ncol(qMask), nl=nrow(Qmodel))    # Streamflow raster for each time-step

  for (i in 1:nrow(Qmodel)){

    # Replace values for each subbasin (from Run_GR2MSemiDistr)
    foreach (j=1:ncol(Qmodel)) %do% {
      qRaster[qMask==j] <- Qmodel[i,j]
    }

    # Import Qsim rasters to GRASS
    writeRAST(as(qRas, 'SpatialGridDataFrame'), "qweight", overwrite=T)

    # Weighted flow accumulation
    execGRASS("r.accumulate", flags=c("overwrite"),  Sys_show.output.on.console=F,
              parameters=list(direction="fdr", weight='qweight', accumulation="qacum"))

    # Load GRASS raster into R
    qAcum <- raster(readRAST('qacum'))
    cat('\f')

    # Extract routing Qsim for each subbasin
    foreach (w=1:ncol(Qmodel)) %do% {
      qSub[i,w] <- maxValue(setMinMax(mask(qAcum, basin[w,])))
    }

    # Store Qacum data
    qArray[,,i] <- round(as.matrix(qAcum),3)
  }

  # Convert an array to a raster brick
  qBrick         <- setValues(qBrick, qArray)
  crs(qBrick)    <- crs(qMask)
  extent(qBrick) <- extent(qMask)
  res(qBrick)    <- res(qMask)
  # qBrick[qBrick==0]<-NA

  # Clean GRASS workspace
  unlink(file.path(getwd(), "GRASS_TEMP"), recursive=T)
  toc()

  # Results
  Ans <- list(Qsub=qSub, Qacum=qBrick)
  return(Ans)

} #End (not run)
