#' Routing discharges for each subbasin.
#' @param Model        List of model results from \code{Run_GR2MSemiDistr}.
#' @param Subbasins		 Subbasins' shapefile. Must contain the following attributes: 'Area' (in km2), 'Region' (in letters), and 'COMID' (identifier number).
#' @param Dem          Digital elevation model raster for the extent of the basin.
#' @param AcumIni      Initial date of the model routing in 'mm/yyyy' format. NULL as default
#' @param AcumEnd      Ending date of the model routing in 'mm/yyyy' format. NULL as default
#' @param Save         Boolean to save results as a text file in the 'Outputs' location. FALSE as default.
#' @param Update       Boolean for the updating mode where only the last month's values will be returned. FALSE as default.
#' @return List of model routing outputs.
#' @return QR: Routed discharge timeseries for all subbasins in [m3/s].
#' @return Dates: Vector of dates of the simulation period.
#' @return COMID: Vector of identifier numbers for each subbasin.
#' @export
#' @examples
#' # Load data
#' require(GR2MSemiDistr)
#' data(dem)
#'
#' # Routing discharges in the streamflow network
#' rou <- Routing_GR2MSemiDistr(Model=model,
#'                              Subbasins=roi,
#'                              Dem=dem,
#'                              AcumIni='01/1981',
#'                              AcumEnd='12/2016')
#' View(rou$QR)
#' plot(rou$Dates, rou$QR[,1], type='l)
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  foreach
#' @import  tictoc
#' @import  parallel
#' @import  lubridate
#' @import  exactextractr
#' @import  sf
Routing_GR2MSemiDistr <- function(Model,
                                  Subbasins,
                                  Dem,
                                  AcumIni=NULL,
                                  AcumEnd=NULL,
                                  Save=FALSE,
                                  Update=FALSE){

  # Require packages
  require(rgdal)
  require(raster)
  require(rgeos)
  require(foreach)
  require(tictoc)
  require(parallel)
  require(lubridate)
  require(exactextractr)
  require(sf)
  tic()
  location <- getwd()

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
                                                          byid = T),
                                          pol@data[!centsInOwnPoly,])
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

  # Subsetting data for routing
  if(is.null(AcumIni)==TRUE & is.null(AcumEnd)==TRUE | nrow(Model$QS)==1){
    Dates <- Model$Dates
    QS    <- as.matrix(Model$QS)
  }else{
    Dates <- seq(as.Date(paste0('01/',AcumIni), format='%d/%m/%Y'),
                 as.Date(paste0('01/',AcumEnd), format='%d/%m/%Y'),
                 by='months')
    Ind   <- seq(which(format(as.Date(Model$Dates),'%d/%m/%Y')==paste0('01/',AcumIni)),
                 which(format(as.Date(Model$Dates),'%d/%m/%Y')==paste0('01/',AcumEnd)))
    QS    <- as.matrix(Model$QS[Ind,])
  }
  if(is.null(ncol(QS))==TRUE){
    nsub  <- length(QS)
    ntime <- 1
  }else{
    nsub  <- ncol(QS)
    ntime <- nrow(QS)
  }

  # Extract cell position for each subbasin (centroid)
  Weight <- Dem
  values(Weight) <- 0
  xycen <- gCentroidWithin(Subbasins)
  index <- extract(Weight, xycen, method='simple', cellnumbers=TRUE, df=TRUE)

  # Pitremove DEM
  dir.create('./Inputs', recursive=T, showWarnings=F)
  setwd('./Inputs')
  writeRaster(Dem, file='dem.tif', overwrite=TRUE)
  system("mpiexec -n 8 pitremove -z dem.tif -fel CONDEM.tif")

  # Flow direction
  system("mpiexec -n 8 D8Flowdir -p FlowDir.tif -sd8 X.tif -fel CONDEM.tif",
         show.output.on.console=F,invisible=F)
  fdr <- raster("FlowDir.tif")
  file.remove('X.tif')
  file.remove('dem.tif')

  # Flow acumulation
  roi_sf <- st_as_sf(Subbasins)
  comid  <- as.vector(roi_sf$COMID)
  ans    <- list()
  for (i in 1:ntime){
    # Show message
    cat('\f')
    message(paste0('Routing streamflows for ',nsub,' sub-basins'))
    message(paste0('Processing...',round(100*i/ntime,3),'%'))
    file.remove('Weights.tif')
    file.remove('FlowAcum.tif')

    # Raster of weights
    if(ntime==1){
      Weight[index$cells] <- QS
    } else{
      Weight[index$cells] <- QS[i,]
    }
    writeRaster(Weight, filename='Weights.tif', overwrite=T)

    # Weighted Flow Accumulation
    system("mpiexec -n 8 AreaD8 -p FlowDir.tif -wg Weights.tif -ad8 FlowAcum.tif")
    fac      <- raster("FlowAcum.tif")
    ans[[i]] <- as.numeric(exact_extract(fac, roi_sf, progress=FALSE,
                                         function(values, coverage_fraction)
                                         max(values[coverage_fraction==1], na.rm=TRUE)))
    rm(fac)
  }
  file.remove('Weights.tif')

  if(Update==TRUE | ntime==1){
    qr <- round(unlist(ans),1)
    qr <- matrix(qr,ncol=length(qr))
  }else{
    qr <- do.call(rbind, ans)
    qr <- round(qr,1)
  }

  # Export results
  Ans <- list(QR=qr,
              Dates=Dates,
              COMID=comid)

  # Save results
  if(Save==TRUE){
    if(Update==TRUE){
      MnYr1     <- format(floor_date(Sys.Date()-months(2), "month"),"%Y%m")
      MnYr2     <- format(floor_date(Sys.Date()-months(1), "month"),"%Y%m")
      QR_name1  <- paste0('QR_GR2MSemiDistr_',MnYr1,'.txt')
      QR_name2  <- paste0('QR_GR2MSemiDistr_',MnYr2,'.txt')
      qr_old    <- read.table(file.path(location,'Outputs',QR_name1), header=TRUE, sep='\t')
      qr_new    <- as.data.frame(rbind(as.matrix(qr_old),as.matrix(qr)))
      colnames(qr_new) <- paste0('QR_',comid)
      rownames(qr_new) <- c(as.Date(rownames(qr_old)),Dates)
      write.table(qr_new, file=file.path(location,'Outputs',QR_name2), sep='\t')
      file.remove(file.path(location,'Outputs',QR_name1))
    }
    if(Update==FALSE){
      qr <- as.data.frame(qr)
      colnames(qr) <- paste0('QR_',comid)
      rownames(qr) <- Dates
      QR_name  <- paste0('QR_GR2MSemiDistr_',format(tail(Dates,1),'%Y%m'),'.txt')
      write.table(qr, file=file.path(location,'Outputs',QR_name), sep='\t')
    }
  }

  message('Done!')
  setwd(location)
  toc()
  return(Ans)

}# End (not run)
