#' Routing simulated monthly streamflows for each subbasin.
#'
#' @param Model        Model results from Run_GR2MSemiDistr.
#' @param Subbasin		 Subbasin shapefile.
#' @param Dem          Raster DEM.
#' @param AcumIni      Initial date for accumulation (in mm/yyyy format).
#' @param AcumEnd      Final date for accumulation (in mm/yyyy format).
#' @param Positions    Cell numbers to extract data faster for each subbasin. NULL as default.
#' @param Save         Boolean to save raster results for each time-step. FALSE as default.
#' @param Update       Boolean to update a previous accumulation file. FALSE as default.
#' @param all          Boolean to consider all the period of model from GR2MSemiDistr. FALSE as default
#' @return  Export and save an accumulation csv file.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  foreach
#' @import  tictoc
#' @import  parallel
#' @import  lubridate
Routing_GR2MSemiDistr <- function(Model,
                                  Subbasin,
                                  Dem,
                                  AcumIni,
                                  AcumEnd,
                                  Positions=NULL,
                                  Save=FALSE,
                                  Update=FALSE,
                                  all=FALSE){

  # Model=Ans3
  # Subbasin=roi
  # Dem='Subbasins.tif'
  # AcumIni='11/2016'
  # AcumEnd='12/2016'
  # Positions=NULL
  # Save=FALSE
  # Update=FALSE
  # all=FALSE

  # Require packages
  require(rgdal)
  require(raster)
  require(rgeos)
  require(foreach)
  require(tictoc)
  require(parallel)
  require(lubridate)
  tic()
  loc <- getwd()

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


  # Extract cell position for each subbasin (centroid)
  qMask <- raster(file.path(getwd(),'Inputs',Dem))
  values(qMask) <- 0
  xycen <- gCentroidWithin(Subbasin)
  index <- extract(qMask, xycen, method='simple', cellnumbers=TRUE, df=TRUE)

  # Pitremove DEM
  setwd(file.path(getwd(),'Inputs'))
  system(paste0("mpiexec -n 8 pitremove -z ",Dem," -fel Ras.tif"))

  # Create flow direction raster
  system("mpiexec -n 8 D8Flowdir -p Flow_Direction.tif -sd8 X.tif -fel Ras.tif",
         show.output.on.console=F,invisible=F)
  fdr <- as.matrix(raster("Flow_Direction.tif"))
  file.remove('Ras.tif')
  file.remove('X.tif')

  # Streamflow accumulation with WFAC
  if(all==TRUE){
    Qmodel <- Model$Qsub
    dates  <- Model$Dates
    if(is.null(ncol(Qmodel))==TRUE){
      nsub  <- length(Qmodel)
      ntime <- 1
    }else{
      nsub   <- ncol(Qmodel)
      ntime  <- nrow(Qmodel)
    }
  }else{
    dates  <- format(seq(as.Date(paste0('01/',AcumIni), format='%d/%m/%Y'),
                         as.Date(paste0('01/',AcumEnd), format='%d/%m/%Y'),
                         by='months'),'%Y-%m-%d')
    Ind    <- seq(which(format(as.Date(Model$Dates),'%d/%m/%Y')==paste0('01/',AcumIni)),
                  which(format(as.Date(Model$Dates),'%d/%m/%Y')==paste0('01/',AcumEnd)))

    Qmodel <- Model$Qsub[Ind,]
    if(is.null(ncol(Qmodel))==TRUE){
      nsub  <- length(Qmodel)
      ntime <- 1
    }else{
      nsub  <- ncol(Qmodel)
      ntime <- nrow(Qmodel)
    }
  }

  # Start loop
  qMon <- list()
  for (i in 1:ntime){
    # Show message
    cat('\f')
    message(paste0('Routing streamflows for ',nsub,' sub-basins'))
    message(paste0('Processing...',round(100*i/ntime,3),'%'))

    # Create a raster of weights (streamflow for each subbasin)
    if(ntime==1){
      qMask[index$cells] <- Qmodel
    } else{
      qMask[index$cells] <- Qmodel[i,]
    }
    writeRaster(qMask, filename='Weights.tif', overwrite=T)

    # Weighted Flow Accumulation
    system("mpiexec -n 8 AreaD8 -p Flow_Direction.tif -wg Weights.tif -ad8 Flow_Accumulation.tif")
    qAcum <- raster("Flow_Accumulation.tif")
    qmat  <- as.matrix(qAcum)

    # Save flow accumulation rasters
    if(Save==TRUE){
      dir.create(file.path(loc,'Outputs','Simulations'))
      NameOut <- paste0('Streamflow_GR2MSemiDistr_',
                        format(as.Date(dates[i]),'%Y%m'),'.tif')
      writeRaster(qAcum, filename=file.path(loc,'Outputs','Simulations',NameOut), overwrite=TRUE)
    }

    # Show message
    message(paste0('Extracting accumulated streamflows for ', format(as.Date(dates[i]),'%b-%Y')))
    message('Please wait..')

    # Positions for extracting accumulated streamflows for each subbasin
    if(is.null(Positions)==TRUE){
      if(i==1){
        cl=makeCluster(detectCores()-1)
        clusterEvalQ(cl,c(library(raster)))
        clusterExport(cl,varlist=c("Subbasin","qAcum","nsub"),envir=environment())
        xycoord <- parLapply(cl, 1:nsub, function(z) {
          ans   <- extract(qAcum, Subbasin[z,], cellnumbers=TRUE, df=TRUE)$cell
          return(ans)
        })
        Positions <- xycoord
        save(Positions, file=file.path(getwd(),'Positions_Routing.Rda'))
      }
    }else{
      if(i==1){
        cl=makeCluster(detectCores()-1)
        clusterEvalQ(cl,c(library(raster)))
      }
      xycoord <- Positions
    }

    # Extracting accumulated streamflows for each subbasin
    clusterExport(cl,varlist=c("qAcum","xycoord","qmat","nsub"),envir=environment())
    fstr <- parLapply(cl, 1:nsub, function(z){
      x   <- rowFromCell(qAcum, xycoord[[z]])
      y   <- colFromCell(qAcum, xycoord[[z]])
      ans <- max(diag(qmat[x,y]), na.rm=T)
      return(ans)
    })
    qMon[[i]] <- unlist(fstr)

  }# End loop
  stopCluster(cl)

  if(ntime==1){
    qSub <- matrix(unlist(qMon), nrow=1, ncol=nsub)
  }else{
    qSub <- do.call(rbind,qMon)
  }

  # Remove auxiliary rasters
  file.remove('Weights.tif')
  file.remove("Flow_Accumulation.tif")
  file.remove("Flow_Direction.tif")

  # Export results
  if(Update==TRUE){
    MnYr1     <- format(floor_date(Sys.Date()-months(2), "month"),"%b%y")
    MnYr2     <- format(floor_date(Sys.Date()-months(1), "month"),"%b%y")
    OldName   <- paste0('Routing_GR2MSemiDistr_',MnYr1,'.txt')
    NewName   <- paste0('Routing_GR2MSemiDistr_',MnYr2,'.txt')
    Data      <- read.table(file.path(loc,'Outputs',OldName), header=T, sep='\t')
    Dates     <- as.Date(Data$Dates, tryFormats=c('%Y-%m-%d','%d/%m/%Y'))
    qSub_Old  <- Data[,-1]
    qSub_New  <- rbind(as.matrix(qSub_Old), qSub)
    Dates_New <- c(Dates,as.Date(dates))
    Database  <- data.frame(Dates_New, qSub_New)
    file.remove(file.path(loc,'Outputs',OldName))
  }
  if(Update==FALSE){
    Database <- data.frame(dates, qSub)
    NewName  <- paste0('Routing_GR2MSemiDistr_',format(tail(as.Date(dates),1),'%b%y'),'.txt')
  }
  colnames(Database) <- c('Dates', paste0('GR2M-ID_',1:nsub))
  write.table(Database, file=file.path(loc,'Outputs',NewName), sep='\t',row.names=FALSE)

  message('Done!')
  setwd(loc)
  toc()
  return(Database)
}
