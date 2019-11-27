#' Routing simulated monthly streamflows.
#'
#' @param Location		 Work directory where 'Inputs' folder is located.
#' @param Model        Model results from Run_GR2MSemiDistr
#' @param Shapefile		 Subbasins shapefile.
#' @param Dem          Raster DEM.
#' @param AcumIni      Initial date 'mm/yyyy' for flow accumulation.
#' @param AcumEnd      Final date 'mm/yyyy' for flow accumulation.
#' @param Save         Logical value to save results as rasters. TRUE as default.
#' @return  Routed streamflows for each subbasin.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  foreach
#' @import  tictoc
Routing_GR2MSemiDistr <- function(Location, Model, Shapefile, Dem, AcumIni, AcumEnd, Save=TRUE){

# Location  <- Location
# Model     <- Mod
# Shapefile <- File.Shape
# Dem       <- File.Raster
# AcumIni   <- '09/2019'
# AcumEnd   <- '10/2019'
# Save      <- FALSE

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
  dates <- seq(as.Date(paste0('01/',AcumIni), format='%d/%m/%Y'),
               as.Date(paste0('01/',AcumEnd), format='%d/%m/%Y'),
               by='months')

  if(Save==TRUE){
    baseName <- readline(prompt="Enter a raster basename: " )
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
  Ind    <- seq(which(format(Model$Dates, '%d/%m/%Y')==paste0('01/',AcumIni)),
                which(format(Model$Dates, '%d/%m/%Y')==paste0('01/',AcumEnd)))
  Qmodel <- Model$Qsub[Ind,]
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
      # system("mpiexec -n 8 AreaD8 -p Flow_Direction.tif -ad8 Flow_Accumulation.tif")
      qAcum <- raster("Flow_Accumulation.tif")*rcrop

      if(Save==TRUE){
        # Create 'Ouput' folder
        dir.create(file.path(Location,'Outputs','Raster_simulation'))
        name   <- paste0(baseName,'_',format(dates[i],'%Y-%m'),'.tif')
        qAcum2 <- qAcum
        qAcum2[qAcum2==0] <- NA
        writeRaster(qAcum2, file=file.path(Location,'Outputs','Raster_simulation',name))
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
          ans <- round(max(qAcum[xycoord[[z]]], na.rm=T),5)
          return(ans)
        })
      qSub[i,] <- unlist(qAcumOut)

    }
    # Remove auxiliary raster
    file.remove('Weights.tif')
    file.remove("Flow_Accumulation.tif")

    # Close the cluster
    stopCluster(cl)

  # Results
  Ans <- data.frame(dates, qSub)
  colnames(Ans) <- c('Dates', paste0('ID_',1:ncol(Qmodel)))
  write.table(Ans, file=file.path(Location,'Outputs','Routing_GR2MSemiDistr.csv'), sep=',', row.names=FALSE)

  # Show message
  message('Done!')
  toc()
  return(Ans)

} #End (not run)
