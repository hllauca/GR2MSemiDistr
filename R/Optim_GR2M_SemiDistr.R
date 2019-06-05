#' Optimization of GR2M model parameters with SCE-UA algorithm.
#'
#' @param Parameters  GR2M (X1 and X2) parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Parameters.Min   Minimum GR2M (X1, X2 and f) model parameters values.
#' @param Parameters.Max   Maximum GR2M (X1, X2 and f) model parameters values.
#' @param Max.Optimization Maximum number of functions used in the optimization loop. 5000 as default.
#' @param Optimization Evaluation criteria for GR2M optimization.
#' @param HRU          Calibration region for each subbasin.
#' @param WorkDir      General work directory where data is located.
#' @param Raster       Flow direction raster in GRASS format.
#' @param Shapefile    Subbasins shapefile.
#' @param Input        Model forcing data in airGR format "[DatesR, P, T, Qobs]".
#' @param WarmIni      Initial date 'mm/yyyy' of the warm-up period.
#' @param WarEnd       Final date 'mm/yyyy' of the warm-up period.
#' @param RunIni       Initial date 'mm/yyyy' of the model evaluation period.
#' @param RunEnd       Final date 'mm/yyyy' of the model evaluation period.
#' @param IdBasin      Subbasin ID number to compute outlet model (from shapefile attribute table).
#' @param Remove       Logical value to remove streamflow generated in the IdBasin. FALSE as default.
#' @param No.Omptim    Calibration regions not to optimize.
#' @return Best semidistribute GR2M model parameters.
#' @export0
Optim_GR2M_SemiDistr <- function(Parameters, Parameters.Min, Parameters.Max, Max.Optimization=5000,
                                  Optimization='NSE', HRU, WorkDir, Shapefile, Input,
								  WarmIni, WarmEnd, RunIni, RunEnd, IdBasin, Remove=FALSE, No.Optim=NULL){

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
      idy <- which(No.Optim==rep(HRU,2))
      Stb <- Parameters[idy]

    # Load shapefiles
      path.shp   <- file.path(WorkDir,'Inputs', Shapefile)
      area       <- readOGR(path.shp, verbose=F)
      nsub       <- nrow(area@data)

    # Read input data
      Data        <- read.table(file.path(WorkDir, 'Inputs', 'Inputs_Basins.txt'), sep='\t', header=T)
      Data$DatesR <- as.POSIXct(Data$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))

    # Subset data for the study period
      Subset      <- seq(which(format(Data$DatesR, format="%m/%Y") == WarmIni),
                         which(format(Data$DatesR, format="%m/%Y") == RunEnd))
      Database    <- Data[Subset,]
      time        <- length(Subset)

    # GR2M initial parameters
      nreg      <- length(sort(unique(HRU)))
      Ini.Param <- data.frame(Region=sort(unique(HRU)),
                          X1=Parameters[1:nreg],
                          X2=Parameters[(nreg+1):(2*nreg)],
                          f=Parameters[((2*nreg)+1):length(Parameters)])

    # Objetive function
      OFUN <- function(Variable, HRU2, WorkDir2, nsub2, Database2, time2,
                        RunIni2, RunEnd2, WarmIni2, WarmEnd2, IdBasin2,
                        Remove2, Eval, No.Optim2, idx2, idy2, Stb2){

            # Select model parameters to optimize
            if (is.null(No.Optim2)==TRUE){
              Optimization <- Variable
            } else{
              dta          <- rbind(cbind(idx2[-idy2], Variable), cbind(idy2, Stb2))
              Optimization <- dta[match(sort(dta[,1]), dta[,1]), 2]
            }

            # Auxiliary variables
            qModel     <- matrix(NA, nrow=time , ncol=nsub2)
            qSub       <- vector()
            ParamSub   <- list()
            OutModel   <- list()
            States     <- list()
            EndState   <- list()
            Factor     <- list()
            Inputs     <- list()
            FixInputs  <- list()

            # Model parameters to run GR2M model
            Zone  <- sort(unique(HRU2))
            nreg  <- length(Zone)
            Param <- data.frame(Zona=Zone,
                                X1=Parameters[1:nreg],
                                X2=Parameters[(nreg+1):(2*nreg)],
                                f=Parameters[((2*nreg)+1):length(Parameters)])

            # Start loop for each timestep
            for (i in 1:time2){
              Date <- format(Database2$DatesR[i], "%m/%Y")

                    foreach (j=1:nsub2) %do% {

                    ParamSub[[j]]  <- c(subset(Param$X1, Param$Zona==HRU[j]),
                                        subset(Param$X2, Param$Zona==HRU[j]))
                    Factor[[j]]    <- subset(Param$f, Param$Zona==HRU[j])
                    Inputs[[j]]    <- Database2[,c(1,j+1,j+1+nsub2)]
                    FixInputs[[j]] <- data.frame(DatesR=Inputs[[j]][,1], Factor[[j]]*Inputs[[j]][,c(2,3)])
                    FixInputs[[j]]$DatesR <- as.POSIXct(FixInputs[[j]]$DatesR,
                                                        "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))

                    if (i==1){
                      IniState      <- NULL
                      OutModel[[j]] <- run_gr2m(FixInputs[[j]],
                                                 ParamSub[[j]],
                                                 IniState[[j]],
                                                 Date)
                    }else{
                      States[[j]]   <- OutModel[[j]]$StateEnd
                      OutModel[[j]] <- run_gr2m(FixInputs[[j]],
                                                 ParamSub[[j]],
                                                 States[[j]],
                                                 Date)
                    }
                    qModel[i,j]     <- round(OutModel[[j]]$Qsim,3)
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
            Subset2     <- seq(which(format(Database$DatesR, format="%m/%Y") == RunIni2),
                               which(format(Database$DatesR, format="%m/%Y") == RunEnd2))
            Database2   <- Database[Subset2,]

          # Streamflow simulated at the basin outlet and raster streamflows
            Qsim <- qSub[Subset2]
            Qobs <- Database2$Qmm
            if (Remove2==TRUE){
              Qsub <- qModel
              Qsim <- Qsim - Qsub[,IdBasin2]
            }

          # Evaluation criteria
            if (Eval == 'KGE'){
              H <- 1 - round(KGE(Qsim, Qobs),3)
              return(H)
            }
            if (Eval == 'NSE'){
              H <- 1 - round(NSE(Qsim, Qobs),3)
              return(H)
            }
            if (Eval == 'lnNSE'){
              H <- 1 - round(NSE(ln(Qsim), ln(Qobs)),3)
              return(H)
            }
            if (Eval == 'RMSE'){
              H <- round(rmse(Qsim, Qobs),3)
              return(H)
            }
            if (Eval == 'R'){
              H <- 1 - round(rPearson(Qsim, Qobs),3)
              return(H)
            }
    }

  # Define calibration regions and parameters to optimize
    if (is.null(No.Optim)==TRUE){
      Zone       <- sort(unique(HRU))
    } else{
      Parameters <- Parameters[!(rep(HRU,2) %in% No.Optim)]
      Zone       <- sort(unique(HRU[!(HRU %in% No.Optim)]))
    }
    Parameters.Min <- rep(Parameters.Min, each=length(Zone))
    Parameters.Max <- rep(Parameters.Max, each=length(Zone))

  # Optimization with SCE-UA
    Ans <- sceua(OFUN, pars=Parameters, lower=Parameters.Min, upper=Parameters.Max, maxn=Max.Optimization,
                 Eval=Optimization, HRU2=HRU, WorkDir2=WorkDir, nsub2=nsub, Database2=Database, time2=time,
                 RunIni2=RunIni, RunEnd2=RunEnd, WarmIni2=WarmIni, WarmEnd2=WarmEnd, IdBasin2=IdBasin,
                 Remove2=Remove, No.Optim2=No.Optim, idx2=idx, idy2=idy, Stb2=Stb)

  # Show message
    message("Done!")
    toc()

  # Output
    return(Ans)

} # End (not run)
