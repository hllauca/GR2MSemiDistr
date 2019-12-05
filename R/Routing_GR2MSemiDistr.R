#' Routing simulated monthly streamflows.
#'
#' @param Location		 Work directory where 'Inputs' folder is located.
#' @param Model        Model results from Run_GR2MSemiDistr.
#' @param Shapefile		 Subbasins shapefile.
#' @param Dem          Raster DEM.
#' @param AcumIni      Initial date 'mm/yyyy' for flow accumulation.
#' @param AcumEnd      Final date 'mm/yyyy' for flow accumulation.
#' @param Save         Logical value to save results as rasters. FALSE as default.
#' @param Update`      Logical value to update a previous csv file (operative mode). FALSE as default.
#' @return  Accumulated streamflows for each subbasin as a csv file
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  foreach
#' @import  tictoc
Routing_GR2MSemiDistr <- function(Location, Model, Shapefile, Dem, AcumIni, AcumEnd, Save=FALSE, Update=FALSE){

# Location  <- Location
# Model     <- Mod
# Shapefile <- File.Shape
# Dem       <- File.Raster
# AcumIni   <- '11/2019'
# AcumEnd   <- '11/2019'
# Save      <- FALSE
# Update    <- FALSE

  # Load packages
    require(foreach)
    require(rgdal)
    require(rgeos)
    require(raster)
    require(tictoc)
    require(parallel)
    tic()

  # Load shape subbasins
    area  <- readOGR(file.path(Location,'Inputs', Shapefile), verbose=F)
    dem   <- raster(file.path(Location,'Inputs', Dem))

  # Create dates
    dates <- seq(as.Date(paste0('01/',AcumIni), format='%d/%m/%Y'),
                 as.Date(paste0('01/',AcumEnd), format='%d/%m/%Y'),
                 by='months')

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
    index <- extract(qMask, xycen, method='simple', cellnumbers=TRUE, df=TRUE)

  # Pitremove DEM
    setwd(file.path(Location,'Inputs'))
    system(paste0("mpiexec -n 8 pitremove -z ",Dem," -fel Ras.tif"))

  # Create Flow Direction raster
    system("mpiexec -n 8 D8Flowdir -p Flow_Direction.tif -sd8 X.tif -fel Ras.tif",show.output.on.console=F,invisible=F)
    file.remove('Ras.tif')
    file.remove('X.tif')

  # Accumulate for each time step
    Ind <- seq(which(format(Model$Dates, '%d/%m/%Y')==paste0('01/',AcumIni)),
               which(format(Model$Dates, '%d/%m/%Y')==paste0('01/',AcumEnd)))

    # Subsetting data
      Qmodel <- Model$Qsub[Ind,]
      nSub   <- ncol(Model$Qsub)
      if(is.null(ncol(Qmodel))==FALSE){
        qSub  <- matrix(NA, nrow=nrow(Qmodel), ncol=ncol(Qmodel))  # Streamflow time series
        ntime <- nrow(Qmodel)
      } else{
        ntime <- 1
      }

      for (i in 1:ntime){

        # Show message
          cat('\f')
          message('Routing outputs from Semidistribute GR2M model')
          message(paste0('Timestep: ', format(dates[i],'%b-%Y')))
          message('Please wait..')

        # Replace values for each subbasin (from Run_GR2MSemiDistr)
          if(is.null(ncol(Qmodel))==FALSE){
            qMask[index$cells] <- Qmodel[i,]
          } else{
            qMask[index$cells] <- Qmodel
          }

        # Save raster
          writeRaster(qMask, filename='Weights.tif', overwrite=T)

        # Weighted Flow Accumulation
          system("mpiexec -n 8 AreaD8 -p Flow_Direction.tif -wg Weights.tif -ad8 Flow_Accumulation.tif")
          qAcum <- raster("Flow_Accumulation.tif")

        # Save rasters
        if(Save==TRUE){
          dir.create(file.path(Location,'Outputs','Raster_simulations'))
          NameOut <- paste0('GR2MSemiDistr_',format(dates[i],'%Y-%m'),'.tif')
          writeRaster(qAcum, file=file.path(Location,'Outputs','Raster_simulation',NameOut))
        }

      # Extract routing Qsim for each subbasin
      if(i==1){
        cl=makeCluster(detectCores()-1)
        clusterEvalQ(cl,c(library(raster)))
        clusterExport(cl,varlist=c("area","qAcum"),envir=environment())
        xycoord <- parLapply(cl, 1:nSub, function(z) {
                   ans <- extract(qAcum, area[z,], cellnumbers=TRUE, df=TRUE)$cell
                   return(ans)
        })
      }

      # Extract routing Qsim for each subbasin
        clusterExport(cl,varlist=c("area","qAcum","xycoord"),envir=environment())
        qAcumOut <- parLapply(cl, 1:nSub, function(z) {
                    ans <- round(max(qAcum[xycoord[[z]]], na.rm=T),5)
                    return(ans)
          })

      # Save routed data
        if(is.null(ncol(Qmodel))==TRUE){
          qSub     <- unlist(qAcumOut)
        } else{
          qSub[i,] <- unlist(qAcumOut)
        }

    } # End loop

  # Remove auxiliary raster
    file.remove('Weights.tif')
    file.remove("Flow_Accumulation.tif")

  # Close the cluster
    stopCluster(cl)

  # Export results
    if (Update==TRUE){
      month     <- as.numeric(format(as.Date(cut(Sys.Date(), "month"), "%Y-%m-%d"), "%m"))
      year      <- as.numeric(format(as.Date(cut(Sys.Date(), "month"), "%Y-%m-%d"), "%Y"))
      MnYr1     <- format(as.Date(paste('01',month-1, year, sep="/"),"%d/%m/%Y"),"%b%y")
      MnYr2     <- format(as.Date(paste('01',month-2, year, sep="/"),"%d/%m/%Y"),"%b%y")
      OldName   <- paste0('Routing_GR2MSemiDistr_',MnYr1,'.csv')
      NewName   <- paste0('Routing_GR2MSemiDistr_',MnYr2,'.csv')
      Data      <- read.table(file.path(Location,'Outputs',OldName), header=T, sep=',')
      Dates     <- as.vector(Data$Dates)
      qSub_Old  <- Data[,-1]
      qSub_New  <- rbind(qSub_Old, qSub)
      Dates_New <- c(Dates,as.character(dates))
      Database  <- data.frame(Dates_New, qSub_New)
    } else{
      Database  <- data.frame(dates, qSub)
      NewName   <- 'Routing_GR2MSemiDistr.csv'
    }
    colnames(Database) <- c('Dates', paste0('ID_',1:nSub))
    write.table(Database, file=file.path(Location,'Outputs',NewName), sep=',', row.names=FALSE)

  # Show message
  message('Done!')
  toc()
  return(Ans)

} #End (not run)
