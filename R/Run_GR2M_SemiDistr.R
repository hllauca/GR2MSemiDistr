#' Run lumped GR2M mode for each sub-basin and timestep.
#'
#' @param Parameters  GR2M (X1 and X2) parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Region      Calibration regions for each subbasin.
#' @param Location    General work directory where data is located.
#' @param FlowDir     Flow direction raster in GRASS format (from 1 to 8). 'Flow_Direction.tif' as default.
#' @param Mask        Subbasins centroids mask raster. 'Centroids_mask.tif' as default.
#' @param Shapefile   Subbasins shapefile.
#' @param Input       Model forcing data in airGR format (DatesR,P,T,Qmm). 'Inputs_Basins.txt' as default.
#' @param WarmIni     Initial date 'mm/yyyy' of the warm-up period.
#' @param WarEnd      Final date 'mm/yyyy' of the warm-up period.
#' @param RunIni      Initial date 'mm/yyyy' of the model evaluation period.
#' @param RunEnd      Final date 'mm/yyyy' of the model evaluation period.
#' @param IdBasin     Subbasin ID number to compute outlet model (from shapefile attribute table).
#' @param Remove      Logical value to remove streamflow generated in the IdBasin. FALSE as default.
#' @param Plot        Logical value to plot observed and simulated streamflow timeseries. TRUE as default.
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @param wfac        Logical value to compute streamflows appliying the Weighted Flow Accumulation algorithm. TRUE as default.
#' @return Semidistribute GR2M model outputs for a subbasin.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  rgrass7
#' @import  rtop
#' @import  hydroGOF
#' @import  foreach
#' @import  tictoc
#' @import  ProgGUIinR
Run_GR2M_SemiDistr <- function(Parameters, Region, Location, FlowDir='Flow_Direction.tif', Mask='Centroids_mask.tif',
                               Shapefile, Input='Inputs_Basins.txt', WarmIni, WarmEnd, RunIni, RunEnd,
                               IdBasin, Remove=FALSE, Plot=TRUE, IniState=NULL, wfac=TRUE){

Parameters=Model.Param
FlowDir='Flow_Direction.tif'
Mask='Centroids_mask.tif'
Input='Inputs_Basins.txt'
Region=Model.Region
Location=Location
Shapefile=File.Shape
WarmIni=WarmUp.Ini
WarmEnd=WarmUp.End
RunIni=RunModel.Ini
RunEnd=RunModel.End
IdBasin=Optim.Basin
Remove=Optim.Remove
Plot=TRUE
IniState=NULL
wfac=WFacum

  # Load packages
    require(ProgGUIinR)
    require(rgdal)
    require(raster)
    require(rgeos)
    require(rgrass7)
    require(rtop)
    require(hydroGOF)
    require(foreach)
    require(tictoc)
    tic()

  # Shapefiles and rasters paths
    path.shp   <- file.path(Location,'Inputs', Shapefile)
    path.rast  <- file.path(Location,'Inputs', FlowDir)
    path.mask  <- file.path(Location,'Inputs', Mask)

  # Load shapefiles and rasters
    area       <- readOGR(path.shp, verbose=F)
    surf       <- area@data$Area
    nsub       <- nrow(area@data)
    rast       <- raster(path.rast)

  # Read input data
    Data        <- read.table(file.path(Location, 'Inputs', Input), sep='\t', header=T)
    Data$DatesR <- as.POSIXct(Data$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))

  # Subset data for the study period
    Subset      <- seq(which(format(Data$DatesR, format="%m/%Y") == WarmIni),
                       which(format(Data$DatesR, format="%m/%Y") == RunEnd))
    Database    <- Data[Subset,]
    time        <- length(Subset)

  # Load GRASS (require to be installed previously)
    if (wfac == TRUE){
    loc <- initGRASS('C:/Program Files/GRASS GIS 7.4.4', home=getwd(),
                     gisDbase="GRASS_TEMP", override=TRUE)
    }

  # Auxiliary variables
    qModel     <- matrix(NA, nrow=time , ncol=nsub)
    if (wfac == TRUE){
    qSub       <- matrix(NA, nrow=time , ncol=nsub)
    qMask      <- raster(path.mask)
    qRaster    <- qMask
    qArray     <- array(NA, dim=c(nrow(qMask), ncol(qMask), time))
    qBrick     <- brick(nr=nrow(qMask), nc=ncol(qMask), nl=time)
    } else{
    qSub       <- vector()
    }
    ParamSub   <- list()
    OutModel   <- list()
    States     <- list()
    EndState   <- list()
    Factor     <- list()
    Inputs     <- list()
    FixInputs  <- list()

  # GR2M model parameters
    Zone  <- sort(unique(Region))
    nreg  <- length(Zone)
    Param <- data.frame(Region=sort(unique(Region)),
                            X1=Parameters[1:nreg],
                            X2=Parameters[(nreg+1):(2*nreg)],
                             f=Parameters[((2*nreg)+1):length(Parameters)])

  # Start loop for each timestep
    for (i in 1:time){

      Date  <- format(Database$DatesR[i], "%m/%Y")
      nDays <- days.in.month(as.numeric(format(Database$DatesR[i],'%Y')),
                             as.numeric(format(Database$DatesR[i],'%m')))

          foreach (j=1:nsub) %do% {
                ParamSub[[j]]  <- c(subset(Param$X1, Param$Region==Region[j]), subset(Param$X2, Param$Region==Region[j]))
                Factor[[j]]    <- subset(Param$f, Param$Region==Region[j])
                Inputs[[j]]    <- Database[,c(1,j+1,j+1+nsub)]
                FixInputs[[j]] <- data.frame(DatesR=Inputs[[j]][,1], Factor[[j]]*Inputs[[j]][,c(2,3)])
                FixInputs[[j]]$DatesR <- as.POSIXct(FixInputs[[j]]$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))
                if (i==1){
                OutModel[[j]]  <- GR2MSemiDistr::run_gr2m_step(FixInputs[[j]], ParamSub[[j]], IniState[[j]], Date)
                }else{
                States[[j]]    <- OutModel[[j]]$StateEnd
                OutModel[[j]]  <- GR2MSemiDistr::run_gr2m_step(FixInputs[[j]], ParamSub[[j]], States[[j]], Date)
                }
                qModel[i,j]    <- round(OutModel[[j]]$Qsim*surf[j]/(86.4*nDays),3)
                if (wfac == TRUE){
                qRaster[qMask==j] <- qModel[i,j]
                }
          }

    # Accumulate streamflow at the basin outlet
      if (wfac == TRUE){
        Result      <- GR2MSemiDistr::run_wfac(rast, qRaster, area, nsub, i)
        qSub[i,]    <- round(Result$Qsub,3)
        qArray[,,i] <- as.matrix(Result$Qacum)
      } else{
        qSub[i]  <- round(sum(qModel[i,], na.rm=T),3)
      }
    # Show message
      cat('\f')
      message('Running Semidistribute GR2M model')
      message(paste0('Timestep: ', format(Database$DatesR[i], "%b-%Y")))
      message('Please wait..')
    }# End loop

   # Compute streamflow raster brick
    if (wfac == TRUE){
      qBrick         <- setValues(qBrick, qArray)
      crs(qBrick)    <- crs(qMask)
      extent(qBrick) <- extent(qMask)
      res(qBrick)    <- res(qMask)

      # Clean GRASS workspace
      unlink(file.path(getwd(), "GRASS_TEMP"), recursive=T)
    }

    # Subset data (without warm-up period)
    Subset2     <- seq(which(format(Database$DatesR, format="%m/%Y") == RunIni),
                       which(format(Database$DatesR, format="%m/%Y") == RunEnd))
    Database2   <- Database[Subset2,]


  # Show comparative figure
    Qobs  <- round(Database2$Qm3s,3)
    if (Plot==TRUE){
      x11()
      if (Remove==FALSE){
        if (wfac == TRUE){
          Qsim <- qSub[Subset2, IdBasin]
        } else{
          Qsim <- qSub[Subset2]
        }
      } else{
        if (wfac == TRUE){
          Qsim <- qSub[Subset2, IdBasin] - qModel[Subset2, IdBasin]
          } else{
          Qsim <- qSub[Subset2] - qModel[Subset2]
          }
      }
      ggof(Qsim, Qobs)
    }

  # Streamflow simulated at the basin outlet and raster streamflows
    if (wfac == TRUE){
      QOUT <- qSub[Subset2,]
      QRAS <- qBrick[[Subset2]]
    } else{
      QOUT <- qSub[Subset2]
      QRAS <- NULL
    }

  # Forcing data multiplying by a factor 'f'
    PP  <- matrix(NA, ncol=nsub, nrow=length(Subset2))
    PET <- matrix(NA, ncol=nsub, nrow=length(Subset2))
    for (w in 1:nsub){
      PP[,w] <- subset(Param$f, Param$Region==Region[w])*Database2[,(w+1)]
      PET[,w]<- subset(Param$f, Param$Region==Region[w])*Database2[,(nsub+w)]
    }

  # Model results
  Ans <- list(Qout=QOUT,
              Qras=QRAS,
              Qobs=Qobs,
              Qsub=qModel[Subset2,],
              Precip=PP,
              Pevap=PET,
              Dates=Database2$DatesR,
              EndState=EndState)

  # Show message
  message('Done!')
  toc()

  # Output
  return(Ans)
}
