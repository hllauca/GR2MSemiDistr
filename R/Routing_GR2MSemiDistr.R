#' Routing simulated monthly streamflows.
#'
#' @param Location		 Work directory where 'Inputs' folder is located.
#' @param Qmodel       Simulated streamflow matrix (time, subbasins) from Run_GR2MSemiDistr
#' @param Shapefile		 Subbasins shapefile.
#' @param Dem          Raster DEM.
#' @param RunIni       Initial date 'mm/yyyy' of the model simulation period.
#' @param RunEnd       Final date 'mm/yyyy' of the model simulation period.
#' @param Save         Conditional to save Weighted Flow Accumulation rasters. TRUE as default.
#' @return  Routed streamflows for each subbasin.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  foreach
#' @import  tictoc
Routing_GR2MSemiDistr <- function(Location, Qmodel, Shapefile, Dem, RunIni, RunEnd, Save=TRUE){

# Location  <- Location
# Qmodel    <- Mod$Qsub
# Shapefile <- File.Shape
# Dem       <- File.Raster
# RunIni    <- RunModel.Ini
# RunEnd    <- RunModel.End
# Save      <- TRUE

  # Load packages
  require(foreach)
  require(rgdal)
  require(rgeos)
  require(raster)
  require(tictoc)
  require(parallel)
  tic()

  # Set work directory
  setwd(file.path(Location,'Inputs'))

  # Load inputs
  area  <- readOGR(file.path(Location,'Inputs', Shapefile), verbose=F)
  dem   <- raster(file.path(Location,'Inputs', Dem))
  rcrop <- !is.na(mask(dem, area))
  rcrop[rcrop==0] <- NA

  # Create dates
  dates <- seq(as.Date(paste0('01/',RunIni), format='%d/%m/%Y'),
               as.Date(paste0('01/',RunEnd), format='%d/%m/%Y'),
               by='months')

  if(save==TRUE){
    baseName <- readline(prompt="Enter WFAC raster basename: " )
  }

  # Auxiliary function (from https://stackoverflow.com/questions/44327994/calculate-centroid-within-inside-a-spatialpolygon)
  gCentroidWithin <- function(pol) {
    require(rgeos)

    pol$.tmpID <- 1:length(pol)
    # initially create centroid points with gCentroid
    initialCents <- gCentroid(pol, byid = T)

    # add data of the polygons to the centroids
    centsDF <- SpatialPointsDataFrame(initialCents, pol@data)
    centsDF$isCentroid <- TRUE

    # check whether the centroids are actually INSIDE their polygon
    centsInOwnPoly <- sapply(1:length(pol), function(x) {
      gIntersects(pol[x,], centsDF[x, ])
    })

    if(all(centsInOwnPoly) == TRUE){
      return(centsDF)
    } else {
      # substitue outside centroids with points INSIDE the polygon
      newPoints <- SpatialPointsDataFrame(gPointOnSurface(pol[!centsInOwnPoly, ],
                                                          byid = T), pol@data[!centsInOwnPoly,])
      newPoints$isCentroid <- FALSE
      centsDF <- rbind(centsDF[centsInOwnPoly,], newPoints)

      # order the points like their polygon counterpart based on `.tmpID`
      centsDF <- centsDF[order(centsDF$.tmpID),]

      # remove `.tmpID` column
      centsDF@data <- centsDF@data[, - which(names(centsDF@data) == ".tmpID")]

      cat(paste(length(pol), "polygons;", sum(centsInOwnPoly), "actual centroids;",
                sum(!centsInOwnPoly), "Points corrected \n"))

      return(centsDF)
    }}

  # Create mask points
  qMask <- dem
  values(qMask) <- 0
  xycen <- gCentroidWithin(area)
  val   <- extract(qMask, xycen, method='simple', cellnumbers=TRUE, df=TRUE)

  # Pitremove DEM
  system(paste0("mpiexec -n 8 pitremove -z ",Dem," -fel Ras.tif"))

  # Create Flow Direction raster
  system("mpiexec -n 8 D8Flowdir -p Flow_Direction.tif -sd8 X.tif -fel Ras.tif",show.output.on.console=F,invisible=F)
  file.remove('Ras.tif')
  file.remove('X.tif')

  # For each time step
  qSub <- matrix(NA, nrow=nrow(Qmodel), ncol=ncol(Qmodel))  # Streamflow time series
    for (i in 1:nrow(Qmodel)){

      # Show message
      cat('\f')
      message('Routing outputs from Semidistribute GR2M model')
      message(paste0('Timestep: ', format(dates[i],'%b-%Y')))
      message('Please wait..')

      # Replace values for each subbasin (from Run_GR2MSemiDistr)
      qMask[val$cells] <- Qmodel[i,]

      # Save raster
      writeRaster(qMask, filename='Weights.tif', overwrite=T)

      # Weighted Flow Accumulation
      system("mpiexec -n 8 AreaD8 -p Flow_Direction.tif -wg Weights.tif -ad8 Flow_Accumulation.tif")
      qAcum <- raster("Flow_Accumulation.tif")*rcrop

      if(Save==TRUE){
        # Create 'Ouput' folder
        dir.create(file.path(Location,'Outputs','Raster_simulation'))
        name   <- paste0(baseName,'_',format(dates[i],'%Y-%m'),'.tif')
        writeRaster(qAcum, file=file.path(Location,'Outputs','Raster_simulation',name))
      }

      # Extract routing Qsim for each subbasin
      if(i==1){
        # Mean areal rainfall for each subbasin
        cl=makeCluster(detectCores()-1) #detectar y asignar numero de cluster
        clusterEvalQ(cl,c(library(raster))) #cargar paquete para aplicar a cada nodo
        clusterExport(cl,varlist=c("area","qAcum"),envir=environment())

        xycoord <- parLapply(cl, 1:ncol(Qmodel), function(z) {
          ans <- extract(qAcum, area[z,], cellnumbers=TRUE, df=TRUE)$cell
          return(ans)
        })
      }

      # Extract routing Qsim for each subbasin
      clusterExport(cl,varlist=c("area","qAcum","xycoord"),envir=environment())
      qAcumOut <- parLapply(cl, 1:ncol(Qmodel), function(z) {
          # ans <- maxValue(mask(qAcum, area[z,]))
          ans <- round(max(qAcum[xycoord[[z]]], na.rm=T),3)
          return(ans)
        })
      qSub[i,] <- unlist(qAcumOut)
      file.remove(name)
    }
    # Remove auxiliary raster
    file.remove('Weights.tif')
    file.remove("Flow_Accumulation.tif")

    # Close the cluster
    stopCluster(cl)

    toc()

  # Results
  Ans <- data.frame(dates, qSub)
  colnames(Ans) <- c('Dates', paste0('ID_',1:ncol(Qmodel)))
  write.table(Ans, file=file.path(Location,'Outputs','Routing_GR2MSemiDistr.csv'), sep='\t', row.names=FALSE)

  # Show message
  message('Done!')
  return(Ans)

} #End (not run)
