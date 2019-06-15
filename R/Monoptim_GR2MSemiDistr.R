#' Optimization of GR2M model parameters with SCE-UA algorithm.
#'
#' @param Parameters       GR2M (X1 and X2) model parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Parameters.Min   Minimum GR2M (X1, X2 and f) model parameters values.
#' @param Parameters.Max   Maximum GR2M (X1, X2 and f) model parameters values.
#' @param Max.Optimization Maximum number of functions used in the optimization loop. 5000 as default.
#' @param Optimization     Mono-objective evaluation criteria for GR2M (NSE, lnNSE, KGE, RMSE, R).
#' @param Region           Calibration region for each subbasin.
#' @param Location     General work directory where data is located.
#' @param Raster       Flow direction raster in GRASS format.
#' @param Shapefile    Subbasins shapefile.
#' @param Input        Model forcing data in airGR format (DatesR,P,T,Qmm). 'Inputs_Basins.txt' as default.
#' @param WarmIni      Initial date 'mm/yyyy' of the warm-up period.
#' @param WarEnd       Final date 'mm/yyyy' of the warm-up period.
#' @param RunIni       Initial date 'mm/yyyy' of the model evaluation period.
#' @param RunEnd       Final date 'mm/yyyy' of the model evaluation period.
#' @param IdBasin      Subbasin ID number to compute outlet model (from shapefile attribute table).
#' @param Remove       Logical value to remove streamflow generated in the IdBasin. FALSE as default.
#' @param No.Optim    Calibration regions not to optimize.
#' @return Best semidistribute GR2M model parameters.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  rtop
#' @import  hydroGOF
#' @import  foreach
#' @import  tictoc
Optim_GR2M_SemiDistr <- function(Parameters, Parameters.Min, Parameters.Max, Max.Optimization=5000,
                                 Optimization='NSE', Region, Location, Shapefile, Input='Inputs_Basins.txt',
								                 WarmIni, WarmEnd, RunIni, RunEnd, IdBasin, Remove=FALSE, No.Optim=NULL){


# Parameters=Model.Param
# Parameters.Min=Model.ParMin
# Parameters.Max=Model.ParMax
# Max.Optimization=Optim.Max
# Optimization=Optim.Eval
# Region=Model.Region
# Location=Location
# Shapefile=File.Shape
# Input='Inputs_Basins.txt'
# WarmIni=WarmUp.Ini
# WarmEnd=WarmUp.End
# RunIni=RunModel.Ini
# RunEnd=RunModel.End
# IdBasin=Optim.Basin
# Remove=Optim.Remove
# No.Optim=No.Region

    # Load packages
      require(rgdal)
      require(raster)
      require(rgeos)
      require(rtop)
      require(hydroGOF)
      require(foreach)
      require(tictoc)
      tic()

    # Filtering calibration region no to optimize
      idx <- 1:length(Parameters)
      idy <- which(No.Optim==rep(Region,2))
      Stb <- Parameters[idy]

    # Load shapefiles
      path.shp   <- file.path(Location,'Inputs', Shapefile)
      area       <- readOGR(path.shp, verbose=F)
      surf       <- area@data$Area
      nsub       <- nrow(area@data)

    # Read input data
      Data        <- read.table(file.path(Location, 'Inputs', Input), sep='\t', header=T)
      Data$DatesR <- as.POSIXct(Data$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))

    # Subset data for the study period
      Subset      <- seq(which(format(Data$DatesR, format="%m/%Y") == WarmIni),
                         which(format(Data$DatesR, format="%m/%Y") == RunEnd))
      Database    <- Data[Subset,]
      time        <- length(Subset)

    # GR2M initial parameters
      nreg      <- length(sort(unique(Region)))
      Ini.Param <- data.frame(Region=sort(unique(Region)),
                              X1=Parameters[1:nreg],
                              X2=Parameters[(nreg+1):(2*nreg)],
                              f=Parameters[((2*nreg)+1):length(Parameters)])

      # Define calibration regions and parameters ranges to optimize
      if (is.null(No.Optim)==TRUE){
        Zone       <- sort(unique(Region))
      } else{
        Parameters <- Parameters[!(rep(Region,2) %in% No.Optim)]
        Zone       <- sort(unique(Region[!(Region %in% No.Optim)]))
      }
      Parameters.Min <- rep(Parameters.Min, each=length(Zone))
      Parameters.Max <- rep(Parameters.Max, each=length(Zone))


    # Objetive function
      OFUN <- function(Variable, Region2, Location2, nsub2, Database2, time2,
                       RunIni2, RunEnd2, WarmIni2, WarmEnd2, IdBasin2,
                       Remove2, Eval, No.Optim2, idx2, idy2, Stb2){

            # Select model parameters to optimize
            if (is.null(No.Optim2)==TRUE){
              Par.Optim <- Variable
            } else{
              dta       <- rbind(cbind(idx2[-idy2], Variable), cbind(idy2, Stb2))
              Par.Optim <- dta[match(sort(dta[,1]), dta[,1]), 2]
            }

            # Auxiliary variables
            qModel     <- matrix(NA, nrow=time2, ncol=nsub2)
            qSub       <- vector()
            ParamSub   <- list()
            OutModel   <- list()
            States     <- list()
            EndState   <- list()
            Factor     <- list()
            Inputs     <- list()
            FixInputs  <- list()

            # Model parameters to run GR2M model
            Zone2 <- sort(unique(Region2))
            nreg  <- length(Zone2)
            Param <- data.frame(Zona=Zone2,
                                X1=Par.Optim[1:nreg],
                                X2=Par.Optim[(nreg+1):(2*nreg)],
                                f=Par.Optim[((2*nreg)+1):length(Par.Optim)])

            # Start loop for each timestep
            for (i in 1:time2){
              Date  <- format(Database2$DatesR[i], "%m/%Y")
              nDays <- days.in.month(as.numeric(format(Database2$DatesR[i],'%Y')),
                                     as.numeric(format(Database2$DatesR[i],'%m')))

                 foreach (j=1:nsub2) %do% {

                    ParamSub[[j]]  <- c(subset(Param$X1, Param$Zona==Region2[j]),
                                        subset(Param$X2, Param$Zona==Region2[j]))
                    Factor[[j]]    <- subset(Param$f, Param$Zona==Region2[j])
                    Inputs[[j]]    <- Database2[,c(1,j+1,j+1+nsub2)]
                    FixInputs[[j]] <- data.frame(DatesR=Inputs[[j]][,1], Factor[[j]]*Inputs[[j]][,c(2,3)])
                    FixInputs[[j]]$DatesR <- as.POSIXct(FixInputs[[j]]$DatesR,"GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))
                    if (i==1){
                      IniState      <- NULL
                      OutModel[[j]] <- GR2MSemiDistr::run_gr2m_step(FixInputs[[j]], ParamSub[[j]], IniState, Date)
                    }else{
                      States[[j]]   <- OutModel[[j]]$StateEnd
                      OutModel[[j]] <- GR2MSemiDistr::run_gr2m_step(FixInputs[[j]], ParamSub[[j]], States[[j]], Date)
                    }
                    qModel[i,j]     <- round(OutModel[[j]]$Qsim*surf[j]/(86.4*nDays),3)
                }
              qSub[i] <- round(sum(qModel[i,], na.rm=T),3)

              # Show message
                cat('\f')
                message('Optimazing with SCE-UA')
                message('======================')
                message('Initial parameters:')
                message(paste0(capture.output(Ini.Param), collapse = "\n"))
                message(' ')
                message('Running Semidistribute GR2M model')
                message(paste0('Time step: ', format(Database$DatesR[i], "%b-%Y")))
                message('Please wait..')
            } #End loop

          # Subset data (without warm-up period)
            Subset2     <- seq(which(format(Database2$DatesR, format="%m/%Y") == RunIni2),
                               which(format(Database2$DatesR, format="%m/%Y") == RunEnd2))
            Database3   <- Database2[Subset2,]

          # Streamflow simulated at the basin outlet and raster streamflows
            Qsim <- qSub[Subset2]
            Qobs <- Database3$Qm3s
            if (Remove2==TRUE){
              Qsub <- qModel
              Qsim <- Qsim - Qsub[,IdBasin2]
            }

          # Evaluation criteria
            if (Eval == 'KGE'){
              H <- 1 - round(KGE(Qsim, Qobs),3)
            }
            if (Eval == 'NSE'){
              H <- 1 - round(NSE(Qsim, Qobs),3)
            }
            if (Eval == 'lnNSE'){
              H <- 1 - round(NSE(ln(Qsim), ln(Qobs)),3)
            }
            if (Eval == 'RMSE'){
              H <- round(rmse(Qsim, Qobs),3)
            }
            if (Eval == 'R'){
              H <- 1 - round(rPearson(Qsim, Qobs),3)
            }
          return(H)
    } # End objectibe function


  # Optimization with SCE-UA
    Ans <- sceua(OFUN, pars=Parameters, lower=Parameters.Min, upper=Parameters.Max, maxn=Max.Optimization,
                 Eval=Optimization, Region2=Region, Location2=Location, nsub2=nsub, Database2=Database, time2=time,
                 RunIni2=RunIni, RunEnd2=RunEnd, WarmIni2=WarmIni, WarmEnd2=WarmEnd, IdBasin2=IdBasin,
                 Remove2=Remove, No.Optim2=No.Optim, idx2=idx, idy2=idy, Stb2=Stb)

  # Show message
    message("Done!")
    toc()

  # Output
    return(Ans)

} # End (not run)
