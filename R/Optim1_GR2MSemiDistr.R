#' Optimization of GR2M model parameters with SCE-UA algorithm.
#'
#' @param Parameters      GR2M (X1 and X2) model parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Parameters.Min  Minimum GR2M (X1, X2 and Fpet) model parameters values.
#' @param Parameters.Max  Maximum GR2M (X1, X2 and Fpet) model parameters values.
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
#' @import  parallel
#' @import  tictoc
#' @import  airGR
Optim1_GR2MSemiDistr <- function(Parameters, Parameters.Min, Parameters.Max, Max.Functions=10000,
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
# Max.Functions=1000

      # Load packages
      require(rgdal)
      require(raster)
      require(rgeos)
	    require(ProgGUIinR)
      require(rtop)
      require(hydroGOF)
      require(parallel)
      require(tictoc)
      require(airGR)
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
        Subset    <- seq(which(format(Data$DatesR, format="%m/%Y") == RunIni),
                         which(format(Data$DatesR, format="%m/%Y") == RunEnd))
      } else{
        Subset    <- seq(which(format(Data$DatesR, format="%m/%Y") == WarmIni),
                         which(format(Data$DatesR, format="%m/%Y") == RunEnd))
      }
      Database    <- Data[Subset,]
      time        <- length(Subset)

      # Number of days in a month (to convert mm to m3/s)
      nDays <- c()
      for (j in 1:time){
        nDays[j] <- days.in.month(as.numeric(format(Database$DatesR[j],'%Y')),
                                  as.numeric(format(Database$DatesR[j],'%m')))
      }

      # Filtering calibration region and parameters not to optimize
      Num.Parameters      <- 1:length(Parameters)
      IDStable.Parameters <- which(No.Optim==rep(sort(unique(region)),4))
      Stable.Parameters   <- Parameters[IDStable.Parameters]

      # Define calibration regions and parameters ranges to optimize
      if (is.null(No.Optim)==TRUE){
        Opt.Parameters <- Parameters
        Opt.Region     <- unique(region)
      } else{
        Opt.Parameters <- Parameters[!(rep(sort(unique(region)),4) %in% No.Optim)]
        Opt.Region     <- unique(region[!(region %in% No.Optim)])
      }
      Opt.Parameters.Min <- rep(Parameters.Min, each=length(Opt.Region))
      Opt.Parameters.Min <- rep(Parameters.Max, each=length(Opt.Region))
      Opt.Parameters.Log <- rep(c(TRUE, TRUE, FALSE, FALSE), each=length(Opt.Region))

      # Utils fucntions
      Subset_Param <- function(Param, Region){
        ParamSub  <- c(subset(Param$X1, Param$Region==Region), subset(Param$X2, Param$Region==Region))
        return(ParamSub)
      }

      Forcing_Subbasin <- function(Param, Region, Database, Nsub, ID){
        FactorPP  <- subset(Param$Fpp, Param$Region==Region)
        FactorPET <- subset(Param$Fpet, Param$Region==Region)
        Inputs    <- Database[,c(1,ID+1,ID+1+Nsub)]
        FixInputs <- data.frame(DatesR=Inputs[,1], P=round(FactorPP*Inputs[,2],1), E=round(FactorPET*Inputs[,3],1))
        FixInputs$DatesR <- as.POSIXct(FixInputs$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))
        return(FixInputs)
      }

      # Objective function
      OFUN <- function(Variable){

            # Select model parameters to optimize
            if (is.null(No.Optim)==TRUE){
              All.Parameters <- Variable
            } else{
              New.Parameters <- rbind(cbind(Num.Parameters[-IDStable.Parameters], Variable), cbind(IDStable.Parameters, Stable.Parameters))
              All.Parameters      <- New.Parameters[match(Num.Parameters, New.Parameters[,1]), 2]
            }

            # Model parameters to run GR2M model
            nreg  <- length(sort(unique(region)))
            Param <- data.frame(Region=sort(unique(region)),
                                X1=All.Parameters[1:nreg],
                                X2=All.Parameters[(nreg+1):(2*nreg)],
                                Fpp=All.Parameters[(2*nreg+1):(3*nreg)],
                                Fpet=All.Parameters[(3*nreg+1):length(All.Parameters)])

            cl=makeCluster(detectCores()-1) # Detect and assign a cluster number
            clusterEvalQ(cl,c(library(GR2MSemiDistr))) # Load package to each node
            clusterExport(cl,varlist=c("Param","region","nsub","Database","time",
                                       "IniState","Subset_Param","Forcing_Subbasin"),envir=environment())

            ResModel <- parLapply(cl, 1:nsub, function(i) {

                        # Parameters and factors to run the model
                        ParamSub  <- Subset_Param(Param, region[i])
                        FixInputs <- Forcing_Subbasin(Param, region[i], Database, nsub, i)

                        # Prepare model inputs
                        InputsModel <- CreateInputsModel(FUN_MOD=RunModel_GR2M,
                                                         DatesR=FixInputs$DatesR,
                                                         Precip=FixInputs$P,
                                                         PotEvap=FixInputs$E)

                        # Run GR2M model by an specific initial conditions
                        if(is.null(IniState)==TRUE){

                          # Set-up running options
                          RunOptions <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                                         InputsModel=InputsModel,
                                                         IndPeriod_Run=1:time,
                                                         verbose=FALSE,
                                                         warnings=FALSE)
                        } else{
                          # Set-up running options
                          RunOptions <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                                         InputsModel=InputsModel,
                                                         IniStates=IniState[[i]],
                                                         IndPeriod_Run=1:time,
                                                         verbose=FALSE,
                                                         warnings=FALSE)
                        }

                        # Run GR2M
                        OutputsModel <- RunModel(InputsModel=InputsModel,
                                                 RunOptions=RunOptions,
                                                 Param=ParamSub,
                                                 FUN=RunModel_GR2M)

                        return(OutputsModel)
                      })

            # Close the cluster
            stopCluster(cl)

            # Main model results (Qsim in m3/s and EndState variables)
            if (nsub==1){
              # Streamflow at the basin outlet
              qSub <- (area[1]*ResModel[[1]]$Qsim)/(86.4*nDays)
              qOut <- qSub
            } else{
              # Streamflow at the basin outlet
              Qlist <- list()
              for(w in 1:nsub){Qlist[[w]] <- (area[w]*ResModel[[w]]$Qsim)/(86.4*nDays)}
              qSub <- do.call(cbind, Qlist)
              qOut <- round(apply(qSub, 1, FUN=sum),2)
            }

            # Subset data (without warm-up period)
              Subset2     <- seq(which(format(Database$DatesR, format="%m/%Y") == RunIni),
                                 which(format(Database$DatesR, format="%m/%Y") == RunEnd))
              Database2   <- Database[Subset2,]

            # Evaluation criteria at the outlet
            Qobs <- Database2$Qm3s
            Qsim <- qOut[Subset2]
            if (Remove==TRUE){
              Qsim <- Qsim - qSub[Subset2, IdBasin]
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
    # Show message
    cat('\f')
    message(paste('Optimizing', Optimization, 'with SCE-UA'))
    message('Please wait...')
    Calibration <- sceua(OFUN, pars=Opt.Parameters, lower=Opt.Parameters.Min, upper=Opt.Parameters.Min,
                         plog=Opt.Parameters.Log, maxn=Max.Functions, kstop=3, pcento=0.1)

    # Extracting calibration results
    if (Optimization == 'PBIAS' | Optimization == 'RMSE'){
      fo <- round(Calibration$value,3)
    } else{
      fo <- round(1-Calibration$value,3)
    }
    if (is.null(No.Optim)==TRUE){
      All.Parameters <- round(Calibration$par,3)
    } else{
      New.Parameters <- rbind(cbind(Num.Parameters[-IDStable.Parameters], round(Calibration$par,3)), cbind(IDStable.Parameters, Stable.Parameters))
      All.Parameters <- New.Parameters[match(Num.Parameters, New.Parameters[,1]), 2]
    }
    Ans <- list(Param=All.Parameters, Value=fo)

    # Show message
    message("Done!")
    toc()

    # Output
    return(Ans)

} # End (not run)
