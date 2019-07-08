#' Run lumped GR2M mode for each sub-basin and timestep.
#'
#' @param Parameters  GR2M (X1 and X2) parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Location    General work directory where data is located.
#' @param Shapefile   Subbasins shapefile.
#' @param Input       Model forcing data in airGR format (DatesR,P,T,Qmm). 'Inputs_Basins.txt' as default.
#' @param WarmIni     Initial date 'mm/yyyy' of the warm-up period.
#' @param RunIni      Initial date 'mm/yyyy' of the model evaluation period.
#' @param RunEnd      Final date 'mm/yyyy' of the model evaluation period.
#' @param IdBasin     Subbasin ID number to compute outlet model (from shapefile attribute table).
#' @param Remove      Logical value to remove streamflow generated in the IdBasin. FALSE as default.
#' @param Plot        Logical value to plot observed and simulated streamflow timeseries. TRUE as default.
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @return Semidistribute GR2M model outputs for a subbasin.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  rtop
#' @import  hydroGOF
#' @import  foreach
#' @import  tictoc
#' @import  ProgGUIinR
Run_GR2MSemiDistr <- function(Parameters, Location, Shapefile, Input='Inputs_Basins.txt',
                              WarmIni=NULL, RunIni, RunEnd, IdBasin, Remove=FALSE,
                              Plot=TRUE, IniState=NULL){

# Parameters=Model.Param
# Input='Inputs_Basins.txt'
# Location=Location
# Shapefile=File.Shape
# WarmIni=WarmUp.Ini
# RunIni=RunModel.Ini
# RunEnd=RunModel.End
# IdBasin=Optim.Basin
# Remove=Optim.Remove
# Plot=TRUE
# IniState=NULL

  # Load packages
    require(ProgGUIinR)
    require(rgdal)
    require(raster)
    require(rgeos)
    require(rtop)
    require(hydroGOF)
    require(foreach)
    require(tictoc)
    tic()

  # Load shapefile
    path.shp   <- file.path(Location,'Inputs', Shapefile)
    basin      <- readOGR(path.shp, verbose=F)
    area       <- basin@data$Area
    region     <- basin@data$Region
    nsub       <- nrow(basin@data)

  # Read and subset input data for the study period
    Data        <- read.table(file.path(Location, 'Inputs', Input), sep='\t', header=T)
    Data$DatesR <- as.POSIXct(Data$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))
    if(is.null(WarmIni)==TRUE){
      Subset      <- seq(which(format(Data$DatesR, format="%m/%Y") == RunIni),
                         which(format(Data$DatesR, format="%m/%Y") == RunEnd))
    } else{
      Subset      <- seq(which(format(Data$DatesR, format="%m/%Y") == WarmIni),
                         which(format(Data$DatesR, format="%m/%Y") == RunEnd))
    }
    Database    <- Data[Subset,]
    time        <- length(Subset)

  # Auxiliary variables
    qSub       <- matrix(NA, nrow=time, ncol=nsub)
    qOut       <- vector()
    ParamSub   <- list()
    OutModel   <- list()
    States     <- list()
    EndState   <- list()
    FactorPET  <- list()
    Inputs     <- list()
    FixInputs  <- list()

  # GR2M model parameters
    Zone  <- sort(unique(region))
    nreg  <- length(Zone)
    Param <- data.frame(Region=sort(unique(region)),
                            X1=Parameters[1:nreg],
                            X2=Parameters[(nreg+1):(2*nreg)],
                            Fpet=Parameters[((2*nreg)+1):length(Parameters)])

  # Start loop for each timestep
    for (i in 1:time){
      Date  <- format(Database$DatesR[i], "%m/%Y")
      nDays <- days.in.month(as.numeric(format(Database$DatesR[i],'%Y')),
                             as.numeric(format(Database$DatesR[i],'%m')))

          foreach (j=1:nsub) %do% {
                ParamSub[[j]]  <- c(subset(Param$X1, Param$Region==region[j]), subset(Param$X2, Param$Region==region[j]))
                FactorPET[[j]] <- subset(Param$Fpet, Param$Region==region[j])
                Inputs[[j]]    <- Database[,c(1,j+1,j+1+nsub)]
                FixInputs[[j]] <- data.frame(DatesR=Inputs[[j]][,1], P=Inputs[[j]][,2], E=round(FactorPET[[j]]*Inputs[[j]][,3],2))
                FixInputs[[j]]$DatesR <- as.POSIXct(FixInputs[[j]]$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))
                if (i==1){
                OutModel[[j]]  <- GR2MSemiDistr::run_gr2m_step(FixInputs[[j]], ParamSub[[j]], IniState[[j]], Date)
                }else{
                States[[j]]    <- OutModel[[j]]$StateEnd
                OutModel[[j]]  <- GR2MSemiDistr::run_gr2m_step(FixInputs[[j]], ParamSub[[j]], States[[j]], Date)
                }
                qSub[i,j]      <- round(OutModel[[j]]$Qsim*area[j]/(86.4*nDays),3)

                # Save last state from GR2M
                if(i==time){
                EndState[[j]]  <- OutModel[[j]]$StateEnd
                }

          }

    # Accumulate streamflow at the basin outlet
      qOut[i]  <- round(sum(qSub[i,]),3)

    # Show message
      cat('\f')
      message('Running Semidistribute GR2M model')
      message(paste0('Timestep: ', format(Database$DatesR[i], "%b-%Y")))
      message('Please wait..')
    }# End loop

    # Subset data (without warm-up period)
      Subset2     <- seq(which(format(Database$DatesR, format="%m/%Y") == RunIni),
                         which(format(Database$DatesR, format="%m/%Y") == RunEnd))
      Database2   <- Database[Subset2,]

    # Evaluation criteria at the outlet
      Qall <- qOut
      Qobs <- Database2$Qm3s
      Qsim <- qOut[Subset2]
      if (Remove==TRUE){
        Qsim <- Qsim - qSub[Subset2, IdBasin]
        Qall <- Qall - qSub[,IdBasin]
      }
      evaluation <- data.frame(KGE=round(KGE(Qsim, Qobs), 3),
                               NSE=round(NSE(Qsim, Qobs), 3),
                               lnNSE=round(NSE(log(Qsim), log(Qobs)), 3),
                               R=round(rPearson(Qsim, Qobs), 3),
                               RMSE=round(rmse(Qsim, Qobs), 3),
                               PBIAS=round(pbias(Qsim, Qobs), 3))

    # Show comparative figure
      if (Plot==TRUE){
        x11()
        ggof(Qsim, Qobs, main=sub('.shp', '',Shapefile), digits=3, gofs=c("NSE", "KGE", "r", "RMSE", "PBIAS"))
      }

  # Forcing data multiplying by a factor 'Fpet'
    pp  <- matrix(NA, ncol=nsub, nrow=length(Subset2))
    pet <- matrix(NA, ncol=nsub, nrow=length(Subset2))
    for (w in 1:nsub){
      pp[,w] <- Database2[,(w+1)]
      pet[,w]<- subset(Param$Fpet, Param$Region==region[w])*Database2[,(nsub+w)]
    }

  # Model results
    Ans <- list(Qsim=Qsim,
                Qobs=Qobs,
                Qsub=qSub[Subset2,],
                Qall=Qall,
                Precip=pp,
                Evaptr=pet,
                Dates=Database2$DatesR,
                EndState=EndState,
                Eval=evaluation)

  # Show message
    message('Done!')
    toc()

  # Output
    return(Ans)
}
