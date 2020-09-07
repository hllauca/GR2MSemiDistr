#' Routing simulated monthly streamflows for each subbasin.
#'
#' @param Model        Model results from Run_GR2MSemiDistr.
#' @param Subbasins		 Subbasins shapefile.
#' @param Dem          Raster DEM.
#' @param AcumIni      Initial date for accumulation (in mm/yyyy format). NULL as default
#' @param AcumEnd      Final date for accumulation (in mm/yyyy format). NULL as default
#' @param Positions    Cell numbers to extract data faster for each subbasin. NULL as default.
#' @param Save         Boolean to results as text file. FALSE as default.
#' @param Update       Boolean to update a previous accumulation file. FALSE as default.
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
                                  Subbasins,
                                  Dem,
                                  AcumIni=NULL,
                                  AcumEnd=NULL,
                                  Positions=NULL,
                                  Save=FALSE,
                                  Update=FALSE){

  # Model=Ans3
  # Subbasins=roi
  # Dem='Subbasins.tif'
  # AcumIni='11/2016'
  # AcumEnd='12/2016'
  # Positions=NULL
  # Save=FALSE
  # Update=FALSE

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
  Qmask <- raster(file.path(getwd(),'Inputs',Dem))
  values(Qmask) <- 0
  xycen <- gCentroidWithin(Subbasins)
  index <- extract(Qmask, xycen, method='simple', cellnumbers=TRUE, df=TRUE)

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
  if(is.null(AcumIni)==TRUE & is.null(AcumIni)==TRUE){
    Qsub  <- as.matrix(Model$Qsub)
    Dates <- Model$Dates
  }else{
    Dates  <- seq(as.Date(paste0('01/',AcumIni), format='%d/%m/%Y'),
                  as.Date(paste0('01/',AcumEnd), format='%d/%m/%Y'),
                  by='months')
    Ind  <- seq(which(format(as.Date(Model$Dates),'%d/%m/%Y')==paste0('01/',AcumIni)),
                which(format(as.Date(Model$Dates),'%d/%m/%Y')==paste0('01/',AcumEnd)))

    Qsub <- as.matrix(Model$Qsub[Ind,])
  }
  if(is.null(ncol(Qsub))==TRUE){
    nsub  <- length(Qsub)
    ntime <- 1
  }else{
    nsub  <- ncol(Qsub)
    ntime <- nrow(Qsub)
  }

  # Start loop
  Qmon <- list()
  for (i in 1:ntime){
    # Show message
    cat('\f')
    message(paste0('Routing streamflows for ',nsub,' sub-basins'))
    message(paste0('Processing...',round(100*i/ntime,3),'%'))

    # Create a raster of weights (streamflow for each subbasin)
    if(ntime==1){
      Qmask[index$cells] <- Qsub
    } else{
      Qmask[index$cells] <- Qsub[i,]
    }
    writeRaster(Qmask, filename='Weights.tif', overwrite=T)

    # Weighted Flow Accumulation
    system("mpiexec -n 8 AreaD8 -p Flow_Direction.tif -wg Weights.tif -ad8 Flow_Accumulation.tif")
    Qacum <- raster("Flow_Accumulation.tif")
    Qmat  <- as.matrix(Qacum)

    # Show message
    message(paste0('Extracting accumulated streamflows for ', format(as.Date(Dates[i]),'%b-%Y')))
    message('Please wait..')

    # Positions for extracting accumulated streamflows for each subbasin
    if(is.null(Positions)==TRUE){
      if(i==1){
        cl=makeCluster(detectCores()-1)
        clusterEvalQ(cl,c(library(raster)))
        clusterExport(cl,varlist=c("Subbasins","Qacum","nsub"),envir=environment())
        xycoord <- parLapply(cl, 1:nsub, function(z) {
          ans   <- extract(Qacum, Subbasins[z,], cellnumbers=TRUE, df=TRUE)$cell
          return(ans)
        })
        Positions_Rou <- xycoord
        save(Positions_Rou, file=file.path(getwd(),'Positions_Rou.Rda'))
      }
    }else{
      if(i==1){
        cl=makeCluster(detectCores()-1)
        clusterEvalQ(cl,c(library(raster)))
      }
      xycoord <- Positions
    }

    # Extracting accumulated streamflows for each subbasin
    clusterExport(cl,varlist=c("Qacum","xycoord","Qmat","nsub"),envir=environment())
    fstr <- parLapply(cl, 1:nsub, function(z){
      x   <- rowFromCell(Qacum, xycoord[[z]])
      y   <- colFromCell(Qacum, xycoord[[z]])
      ans <- max(diag(Qmat[x,y]), na.rm=T)
      return(ans)
    })
    Qmon[[i]] <- unlist(fstr)

  }# End loop
  stopCluster(cl)

  if(ntime==1){
    Qrou <- round(matrix(unlist(Qmon), nrow=1, ncol=nsub),3)
  }else{
    Qrou <- round(do.call(rbind,Qmon),3)
  }

  # Remove auxiliary rasters
  file.remove('Weights.tif')
  file.remove("Flow_Accumulation.tif")
  file.remove("Flow_Direction.tif")

  # Export results
  if(Save==TRUE){
    if(Update==TRUE){
      MnYr1     <- format(floor_date(Sys.Date()-months(2), "month"),"%b%y")
      MnYr2     <- format(floor_date(Sys.Date()-months(1), "month"),"%b%y")
      OldName   <- paste0('QR_GR2MSemiDistr_',MnYr1,'.txt')
      NewName   <- paste0('QR_GR2MSemiDistr_',MnYr2,'.txt')
      Qrou_Old <- read.table(file.path(loc,'Outputs',OldName), header=TRUE, sep='\t')
      Qrou_New <- as.data.frame(rbind(as.matrix(Qrou_Old),as.matrix(Qrou)))
      colnames(Qrou_New) <- Model$GR2M_ID
      rownames(Qrou_New) <- c(as.Date(rownames(Qrou_Old)),Dates)
      write.table(Qrou_New, file=file.path(loc,'Outputs',NewName), sep='\t')
      file.remove(file.path(loc,'Outputs',OldName))
    }
    if(Update==FALSE){
      Qrou <- as.data.frame(Qrou)
      colnames(Qrou) <- Model$GR2M_ID
      rownames(Qrou) <- Dates
      NewName  <- paste0('QR_GR2MSemiDistr_',
                         format(tail(Dates,1),'%b%y'),'.txt')
      write.table(Qrou, file=file.path(loc,'Outputs',NewName), sep='\t')
    }
  }

  Ans <- as.data.frame(Qrou)
  colnames(Ans) <- Model$GR2M_ID
  if(length(unique(Dates))!=1){
    rownames(Ans) <- Dates
  }
  message('Done!')
  setwd(loc)
  toc()
  return(Ans)
}
