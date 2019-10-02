#' Optimization of GR2M model parameters with MOPSOCD algorithm.
#'
#' @param Location		 General work directory where data is located.
#' @param Qmodel       Simulated streamflow matrix (time, subbasin) from Run_GR2MSemiDistr
#' @param Shapefile		 Subbasins shapefile.
#' @param FlowDir			 Flow direction raster in GRASS format. 'Flow_Direction.tif' as default.
#' @param Mask         Subbasins centroids mask raster. 'Centroids_mask.tif' as default.
#' @return  Routing streamflow for each subbasin.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  foreach
#' @import  tictoc
Routing_GR2MSemiDistr <- function(Location, Qmodel, Shapefile, FlowDir='Flow_Direction.tif', Mask='Centroids_mask.tif'){

# Location=Location
# Qmodel=Ans2$Qsub
# Shapefile=File.Shape
# FlowDir='Flow_Direction.tif'
# Mask='Centroids_mask.tif'

  require(foreach)
  require(rgdal)
  require(raster)
  require(tictoc)
  tic()

  # Load rasters
  setwd(file.path(Location,'Inputs'))
  # path.fdr  <- file.path(Location,'Inputs', FlowDir)
  # path.mask <- file.path(Location,'Inputs', Mask)
  # flowdir   <- raster(path.fdr)
  # qMask     <- raster(path.mask)

  # Load shapefile
  path.shp   <- file.path(Location,'Inputs', Shapefile)
  basin      <- readOGR(path.shp, verbose=F)


  # Auxiliary variables
  qRas   <- qMask
  qSub   <- matrix(NA, nrow=nrow(Qmodel), ncol=ncol(Qmodel))          # Streamflow time series
  qArray <- array(NA, dim=c(nrow(qMask), ncol(qMask), nrow(Qmodel)))  # Streamflow matrix for each time-step
  qBrick <- brick(nr=nrow(qMask), nc=ncol(qMask), nl=nrow(Qmodel))    # Streamflow raster for each time-step

  for (i in 1:nrow(Qmodel)){

    # Show message
    cat('\f')
    message('Routing outputs from Semidistribute GR2M model')
    message(paste0('Timestep: ', i, " from ", nrow(Qmodel)))
    message('Please wait..')

    # Replace values for each subbasin (from Run_GR2MSemiDistr)
    foreach (j=1:ncol(Qmodel)) %do% {
      qRas[qMask==j] <- Qmodel[i,j]
    }

    # Weighted Flow Accumulation
    system("mpiexec -n 8 AreaD8 -p Flow_Direction.tif -wg Centroids_mask.tif -ad8 Flow_Accumulation.tif")
    qAcum=raster("Flow_Accumulation.tif")

    # Extract routing Qsim for each subbasin
    foreach (w=1:ncol(Qmodel)) %do% {
      qSub[i,w] <- maxValue(setMinMax(mask(qAcum, basin[w,])))
    }

    # Store Qacum data
    qArray[,,i] <- round(as.matrix(qAcum),3)
  }

  # Show message
  cat('\f')
  message('Generating a raster brick')
  message('Please wait..')

  # Convert an array to a raster brick
  qBrick         <- setValues(qBrick, qArray)
  crs(qBrick)    <- crs(qMask)
  extent(qBrick) <- extent(qMask)
  res(qBrick)    <- res(qMask)
  # qBrick[qBrick==0]<-NA

  toc()

  # Results
  Ans <- list(Qsub=qSub, Qacum=qBrick)
  return(Ans)

} #End (not run)
