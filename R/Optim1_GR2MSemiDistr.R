#' Optimization of GR2M model parameters with SCE-UA algorithm.
#'
#' @param Parameters      GR2M (X1 and X2) model parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Parameters.Min  Minimum GR2M (X1, X2 and f) model parameters values.
#' @param Parameters.Max  Maximum GR2M (X1, X2 and f) model parameters values.
#' @param Max.Functions 	Maximum number of functions used in the optimization loop. 5000 as default.
#' @param Optimization    Mono-objective evaluation criteria for GR2M (NSE, lnNSE, KGE, RMSE, R, PBIAS).
#' @param Location    General work directory where data is located.
#' @param Raster      Flow direction raster in GRASS format.
#' @param Shapefile   Subbasins shapefile.
#' @param Input       Model forcing data in airGR format (DatesR,P,T,Qmm). 'Inputs_Basins.txt' as default.
#' @param WarmIni     Initial date 'mm/yyyy' of the warm-up period.
#' @param RunIni      Initial date 'mm/yyyy' of the model evaluation period.
#' @param RunEnd      Final date 'mm/yyyy' of the model evaluation period.
#' @param IdBasin     Subbasin ID number to compute outlet model (from shapefile attribute table).
#' @param Remove      Logical value to remove streamflow generated in the IdBasin. FALSE as default.
#' @param No.Optim    Calibration regions not to optimize.
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @return  Best semidistribute GR2M model parameters.
#' @export
#' @import  ProgGUIinR
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  rtop
#' @import  hydroGOF
#' @import  foreach
#' @import  tictoc
Optim1_GR2MSemiDistr <- function(Parameters, Parameters.Min, Parameters.Max, Max.Functions=5000,
									               Optimization='NSE', Location, Shapefile, Input='Inputs_Basins.txt',
									               WarmIni, RunIni, RunEnd, IdBasin, Remove=FALSE, No.Optim=NULL, IniState=NULL){

# Parameters=Model.Param
# Parameters.Min=Model.ParMin
# Parameters.Max=Model.ParMax
# Optimization=Optim.Eval
# Location=Location
# Shapefile=File.Shape
# Input='Inputs_Basins.txt'
# WarmIni=WarmUp.Ini
# RunIni=RunModel.Ini
# RunEnd=RunModel.End
# IdBasin=Optim.Basin
# Remove=Optim.Remove
# No.Optim=No.Region
# IniState=NULL

      # Load packages
      require(rgdal)
      require(raster)
      require(rgeos)
	    require(ProgGUIinR)
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

      # Filtering calibration region no to optimize
      idx <- 1:length(Parameters)
      idy <- which(No.Optim==rep(region,2))
      Stb <- Parameters[idy]

      # Read and subset input data for the study period
      Data        <- read.table(file.path(Location, 'Inputs', Input), sep='\t', header=T)
      Data$DatesR <- as.POSIXct(Data$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))
      if(is.null(WarmIni)==TRUE){
        Subset    <- seq(which(format(Data$DatesR, format="%m/%Y") == RunIni),
                         which(format(Data$DatesR, format="%m/%Y") == RunEnd))
      } else{
        Subset    <- seq(which(format(Data$DatesR, format="%m/%Y") == WarmIni),
                         which(format(Data$DatesR, format="%m/%Y") == RunEnd))
      }
      Database    <- Data[Subset,]
      time        <- length(Subset)

      # GR2M initial parameters
      nreg        <- length(sort(unique(region)))
      Ini.Param   <- data.frame(Region=sort(unique(region)),
                                X1=Parameters[1:nreg],
                                X2=Parameters[(nreg+1):(2*nreg)],
                                 f=Parameters[((2*nreg)+1):length(Parameters)])

      # Define calibration regions and parameters ranges to optimize
      if (is.null(No.Optim)==TRUE){
        Zone       <- sort(unique(region))
      } else{
        Parameters <- Parameters[!(rep(region,2) %in% No.Optim)]
        Zone       <- sort(unique(region[!(region %in% No.Optim)]))
      }
      Parameters.Min <- rep(Parameters.Min, each=length(Zone))
      Parameters.Max <- rep(Parameters.Max, each=length(Zone))

      # Objective function
      OFUN <- function(Variable){

            # Select model parameters to optimize
            if (is.null(No.Optim)==TRUE){
              Par.Optim <- Variable
            } else{
              dta       <- rbind(cbind(idx2[-idy], Variable), cbind(idy, Stb))
              Par.Optim <- dta[match(sort(dta[,1]), dta[,1]), 2]
            }

            # Auxiliary variables
            qSub       <- matrix(NA, nrow=time, ncol=nsub)
            qOut       <- vector()
            ParamSub   <- list()
            OutModel   <- list()
            States     <- list()
            EndState   <- list()
            Factor     <- list()
            Inputs     <- list()
            FixInputs  <- list()

            # Model parameters to run GR2M model
            Param <- data.frame(Region=sort(unique(region)),
                                    X1=Par.Optim[1:nreg],
                                    X2=Par.Optim[(nreg+1):(2*nreg)],
                                     f=Par.Optim[((2*nreg)+1):length(Par.Optim)])

            # Start loop for each timestep
            for (i in 1:time){
              Date  <- format(Database$DatesR[i], "%m/%Y")
              nDays <- days.in.month(as.numeric(format(Database$DatesR[i],'%Y')),
                                     as.numeric(format(Database$DatesR[i],'%m')))

              foreach (j=1:nsub) %do% {
                  ParamSub[[j]]  <- c(subset(Param$X1, Param$Region==region[j]), subset(Param$X2, Param$Region==region[j]))
                  Factor[[j]]    <- subset(Param$f, Param$Region==region[j])
                  Inputs[[j]]    <- Database[,c(1,j+1,j+1+nsub)]
                  FixInputs[[j]] <- data.frame(DatesR=Inputs[[j]][,1], Factor[[j]]*Inputs[[j]][,c(2,3)])
                  FixInputs[[j]]$DatesR <- as.POSIXct(FixInputs[[j]]$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))
                  if (i==1){
                    OutModel[[j]]  <- GR2MSemiDistr::run_gr2m_step(FixInputs[[j]], ParamSub[[j]], IniState[[j]], Date)
                  }else{
                    States[[j]]    <- OutModel[[j]]$StateEnd
                    OutModel[[j]]  <- GR2MSemiDistr::run_gr2m_step(FixInputs[[j]], ParamSub[[j]], States[[j]], Date)
                  }
                  qSub[i,j]      <- round(OutModel[[j]]$Qsim*area[j]/(86.4*nDays),3)
              }

              # Accumulate streamflow at the basin outlet
              qOut[i]  <- round(sum(qSub[i,]),3)

              # Show message
                cat('\f')
                message('Optimizing with SCE-UA')
                message('======================')
                message('Initial parameters:')
                message(paste0(capture.output(Ini.Param), collapse = "\n"))
                message(' ')
                message('Running Semidistribute GR2M model')
                message(paste0('Time step: ', format(Database$DatesR[i], "%b-%Y")))
                message('Please wait..')
            } #End loop

            # Subset data (without warm-up period)
              Subset2     <- seq(which(format(Database$DatesR, format="%m/%Y") == RunIni),
                                 which(format(Database$DatesR, format="%m/%Y") == RunEnd))
              Database2   <- Database[Subset2,]


            # Evaluation criteria at the outlet
            Qobs <- Database2$Qm3s
            Qsim <- qOut[Subset2]
            if (Remove==TRUE){
              Qsim <- Qsim - qSub[Subset, IdBasin]
            }

            # Evaluation criteria dataframe (only minimizing)
            optim.df <- data.frame(KGE=1-round(KGE(Qsim, Qobs), 3),
                                   NSE=1-round(NSE(Qsim, Qobs), 3),
                                   lnNSE=1-round(NSE(log(Qsim), log(Qobs)), 3),
                                   R=1-round(rPearson(Qsim, Qobs), 3),
                                   RMSE=round(rmse(Qsim, Qobs), 3),
                                   PBIAS=round(pbias(Qsim, Qobs), 3))

          # Return
          OF <- as.numeric(optim.df[colnames(optim.df) %in% Optimization])
          return(OF)

    } # End objective function

    # Optimization with SCE-UA
    Calibration <- sceua(OFUN, pars=Parameters, lower=Parameters.Min, upper=Parameters.Max, maxn=Max.Functions)

    # Extracting results
    if (Optimization == 'PBIAS' | Optimization == 'RMSE'){
    fo <- round(Calibration$value,3)
    } else{
    fo <- round(1-Calibration$value,3)
    }
    Ans <- list(Param=round(Calibration$par,3), Value=fo)

    # Show message
    message("Done!")
    toc()

    # Output
    return(Ans)

} # End (not run)
