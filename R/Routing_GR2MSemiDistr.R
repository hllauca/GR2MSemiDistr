#' Routing simulated monthly streamflows for each subbasin.
#'
#' @param Location		 Directory where 'Inputs' folder is located.
#' @param Model        Model results from Run_GR2MSemiDistr.
#' @param Shapefile		 Subbasin shapefile.
#' @param Dem          Raster DEM.
#' @param AcumIni      Initial date 'mm/yyyy' for accumulation.
#' @param AcumEnd      Final date 'mm/yyyy' for accumulation.
#' @param Save         Logical value to save raster results for each time-step. FALSE as default.
#' @param Update       Logical value to update a previous accumulation csv file. FALSE as default.
#' @param Positions    Cell numbers to extract data faster for each subbasin. NULL as default.
#' @param all          Conditional to consider all the period of model from GR2MSemiDistr. FALSE as default
#' @return  Export and save an accumulation csv file.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  foreach
#' @import  tictoc
#' @import  parallel
Routing_GR2MSemiDistr <- function(Location, Model, Shapefile, Dem, AcumIni, AcumEnd,
                                  Save=FALSE, Update=FALSE, Positions=NULL, all=FALSE){

# Location  <- Location
# Model     <- Qsub_forcast
# Shapefile <- File.Shape
# Dem       <- File.Raster
# AcumIni   <- RunModel.Ini
# AcumEnd   <- RunModel.End
# Save      <- FALSE
# Update    <- FALSE
# Positions <- Positions
# all       <- TRUE

  # Load packages
    require(foreach)
    require(rgdal)
    require(rgeos)
    require(raster)
    require(tictoc)
    require(parallel)
    tic()

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

  # Load shapefiles and raster
    area  <- readOGR(file.path(Location,'Inputs', Shapefile), verbose=F)
    dem   <- raster(file.path(Location,'Inputs', Dem))

  # Position for each subbasin centroid
    qMask <- dem
    values(qMask) <- 0
    xycen <- gCentroidWithin(area)
    index <- extract(qMask, xycen, method='simple', cellnumbers=TRUE, df=TRUE)

  # Pitremove DEM
    setwd(file.path(Location,'Inputs'))
    system(paste0("mpiexec -n 8 pitremove -z ",Dem," -fel Ras.tif"))

  # Create raster of flow direction
    system("mpiexec -n 8 D8Flowdir -p Flow_Direction.tif -sd8 X.tif -fel Ras.tif",show.output.on.console=F,invisible=F)
    fdr <- as.matrix(raster("Flow_Direction.tif"))
    file.remove('Ras.tif')
    file.remove('X.tif')

  # Accumulate streamflows for each time step
      # Conditionals
      if (all==TRUE){
        Qmodel <- Model$Qsub
        dates  <- Model$Dates
        if(is.null(ncol(Qmodel))==TRUE){
          nsub  <- length(Qmodel)
          ntime <- 1
        }else{
          nsub   <- ncol(Qmodel)
          ntime  <- nrow(Qmodel)
        }
      } else{
        # Create a vector of dates
        dates  <- seq(as.Date(paste0('01/',AcumIni), format='%d/%m/%Y'),
                      as.Date(paste0('01/',AcumEnd), format='%d/%m/%Y'),
                      by='months')
        Ind    <- seq(which(format(Model$Dates, '%d/%m/%Y')==paste0('01/',AcumIni)),
                      which(format(Model$Dates, '%d/%m/%Y')==paste0('01/',AcumEnd)))
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
      QACUM <- array(NA, dim=c(nrow(fdr), ncol(fdr), ntime))
      for (i in 1:ntime){
        # Show message
          cat('\f')
          message(paste0('Routing streamflows for ',nsub,' sub-basins'))
          message(paste0('Time-step: ', format(dates[i],'%b-%Y')))
          message('Please wait..')

        # Create a raster of weights (streamflow for each subbasin)
          if(ntime == 1){
            qMask[index$cells] <- Qmodel
          } else{
            qMask[index$cells] <- Qmodel[i,]
          }
          writeRaster(qMask, filename='Weights.tif', overwrite=T)

        # Weighted Flow Accumulation
          system("mpiexec -n 8 AreaD8 -p Flow_Direction.tif -wg Weights.tif -ad8 Flow_Accumulation.tif")
          qAcum <- raster("Flow_Accumulation.tif")
          QACUM[,,i] <- as.matrix(qAcum)

        # Save flow accumulation rasters
          if(Save == TRUE){
            dir.create(file.path(Location,'Outputs','Raster_simulations'))
            NameOut <- paste0('GR2MSemiDistr_',format(dates[i],'%Y-%m'),'.tif')
            writeRaster(qAcum, file=file.path(Location,'Outputs','Raster_simulation',NameOut))
          }
      }# End loop

  # Positions for extracting accumulated streamflows for each subbasin
    if(is.null(Positions)==TRUE){
       cl=makeCluster(detectCores()-1)
       clusterEvalQ(cl,c(library(raster)))
       clusterExport(cl,varlist=c("area","qAcum"),envir=environment())
       xycoord <- parLapply(cl, 1:nsub, function(z) {
          ans <- extract(qAcum, area[z,], cellnumbers=TRUE, df=TRUE)$cell
       return(ans)
       })
       stopCluster(cl)
       Positions <- xycoord
       save(Positions, file=file.path(getwd(),'Positions_Routing.Rda'))
    }else{
       cl=makeCluster(detectCores()-1)
       clusterEvalQ(cl,c(library(raster)))
       xycoord <- Positions
    }

  # Extracting accumulated streamflows for each subbasin
    # Show message
    cat('\f')
    message(paste0('Extracting accumulated streamflows for ',nsub,' subbasins'))
    message('Please wait..')
    clusterExport(cl,varlist=c("qAcum","xycoord","QACUM","ntime","nsub"),envir=environment())
    fstr <- parLapply(cl, 1:nsub, function(z) {
            message(paste0('Processing...',round(100*z/nsub,3),'%'))
            x  <- rowFromCell(qAcum, xycoord[[z]])
            y  <- colFromCell(qAcum, xycoord[[z]])
            xy <- QACUM[x,y,]
            if(ntime==1){
              ans <- max(diag(xy))
            }else{
              ans <- c()
              for (k in 1:ntime){
                ans[k] <- max(diag(xy[,,k]))
              }
            }
            return(ans)
          })
    stopCluster(cl)
    if(ntime==1){
      qSub <- unlist(fstr)
    }else{
      qSub <- do.call(cbind, fstr)
    }


  # Remove auxiliary rasters
    file.remove('Weights.tif')
    file.remove("Flow_Accumulation.tif")
    file.remove("Flow_Direction.tif")

  # Export results
  #===============
    if (Update==TRUE){
      month     <- as.numeric(format(as.Date(cut(Sys.Date(), "month"), "%Y-%m-%d"), "%m"))
      year      <- as.numeric(format(as.Date(cut(Sys.Date(), "month"), "%Y-%m-%d"), "%Y"))
      MnYr1     <- format(as.Date(paste('01',month-2, year, sep="/"),"%d/%m/%Y"),"%b%y")
      MnYr2     <- format(as.Date(paste('01',month-1, year, sep="/"),"%d/%m/%Y"),"%b%y")
      OldName   <- paste0('Routing_GR2MSemiDistr_',MnYr1,'.csv')
      NewName   <- paste0('Routing_GR2MSemiDistr_',MnYr2,'.csv')
      Data      <- read.table(file.path(Location,'Outputs',OldName), header=T, sep=',')
      Dates     <- as.vector(Data$Dates)
      qSub_Old  <- Data[,-1]
      qSub_New  <- rbind(qSub_Old, qSub)
      Dates_New <- c(Dates,as.character(dates))
      Database  <- data.frame(Dates_New, qSub_New)
      file.remove(file.path(Location,'Outputs',OldName))
    } else{
      Database  <- data.frame(dates, qSub)
      NewName   <- 'Routing_GR2MSemiDistr.csv'
    }
    colnames(Database) <- c('Dates', paste0('ID_',1:nsub))
    write.table(Database, file=file.path(Location,'Outputs',NewName), sep=',', row.names=FALSE)

    # Show message
      message('Done!')
      toc()
      return(Database)

} #End (not run)
