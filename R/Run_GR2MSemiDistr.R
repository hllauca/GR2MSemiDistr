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
#' @import  airGR
#' @import  tictoc
#' @import  ProgGUIinR
#' @import  parallel
Run_GR2MSemiDistr <- function(Parameters, Location, Shapefile, Input='Inputs_Basins.txt',
                              WarmIni=NULL, RunIni, RunEnd, IdBasin, Remove=FALSE,
                              Plot=TRUE, IniState=NULL){

# Parameters=Model.Param
# Location=Location
# Shapefile=File.Shape
# Input='Inputs_Basins.txt'
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
    require(tictoc)
    require(parallel)
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

  # GR2M model parameters
    Zone  <- sort(unique(region))
    nreg  <- length(Zone)
    Param <- data.frame(Region=sort(unique(region)),
                            X1=Parameters[1:nreg],
                            X2=Parameters[(nreg+1):(2*nreg)],
                            Fpp=Parameters[(2*nreg+1):(3*nreg)],
                            Fpet=Parameters[(3*nreg+1):length(Parameters)])

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

  # Show message
    cat('\f')
    message(paste('Running GR2M model', nsub, 'subbasins'))
    message('Please wait...')

  # Run GR2M for each subbasin
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

      # End state variables
      EndState <- list(ResModel[[1]]$StateEnd)

      } else{
      # Streamflow at the basin outlet
        Qlist <- list()
        for(w in 1:nsub){Qlist[[w]] <- (area[w]*ResModel[[w]]$Qsim)/(86.4*nDays)}
        qSub <- do.call(cbind, Qlist)
        qOut <- round(apply(qSub, 1, FUN=sum),2)

      # End state variables
        EndState <- list()
        for(w in 1:nsub){EndState[[w]] <- ResModel[[w]]$StateEnd}
      }

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

  # Forcing data multiplying by a factor Fpp and Fpet
    pp  <- matrix(NA, ncol=nsub, nrow=length(Subset2))
    pet <- matrix(NA, ncol=nsub, nrow=length(Subset2))
    for (w in 1:nsub){
      pp[,w] <- subset(Param$Fpp, Param$Region==region[w])*Database2[,(w+1)]
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
    # return(Ans)
    return(Ans)
}
