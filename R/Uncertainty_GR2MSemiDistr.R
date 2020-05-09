#' Uncertainty analysis of GR2M model parameters with the MCMC algorithm.
#'
#' @param Parameters      GR2M (X1 and X2) model parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Parameters.Min  Minimum GR2M (X1, X2, fprecip and fpet) model parameters values.
#' @param Parameters.Max  Maximum GR2M (X1, X2, fprecip and fpet) model parameters values.
#' @param Niter 	        Number of iterations. 2000 as default.
#' @param Optimization    Mono-objective evaluation criteria for GR2M (NSE, lnNSE, KGE, RMSE, R, PBIAS).
#' @param Location    Directory where 'Inputs' folder is located.
#' @param Shapefile   Subbasin shapefile.
#' @param Input       Forcing data texfile (Dates, Precip, PotEvap, Qobs). 'Inputs_Basins.txt' as default.
#' @param WarmIni     Initial date (in 'mm/yyyy' format) of the warm-up period.
#' @param RunIni      Initial date (in 'mm/yyyy' format) of the model simulation period.
#' @param RunEnd      Final date (in 'mm/yyyy' format) of the model simulation period.
#' @param IdBasin     ID for the outlet subbasin (from shapefile attribute table).
#' @param Remove      Logical value to remove streamflows of the outlet subbasin (IdBasin). FALSE as default.
#' @param No.Optim    Calibration regions not to be optimized.
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @return  Parameter and streamflow uncertanty bounds.
#' @export
#' @import  ProgGUIinR
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  FME
#' @import  hydroGOF
#' @import  parallel
#' @import  tictoc
#' @import  airGR
Uncertainty_GR2MSemiDistr <- function(Parameters, Parameters.Min, Parameters.Max, Niter=2000,
                                 Optimization='NSE', Location, Shapefile, Input='Inputs_Basins.txt',
                                 WarmIni, RunIni, RunEnd, IdBasin, Remove=FALSE, No.Optim=NULL, IniState=NULL){

  # Parameters=Model.Param
  # Parameters.Min=Model.ParMin
  # Parameters.Max=Model.ParMax
  # Niter=5000
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
  require(parallel)
  require(tictoc)
  require(airGR)
  tic()

  # Read sub-basins
  path.shp   <- file.path(Location,'Inputs', Shapefile)
  basin      <- readOGR(path.shp, verbose=F)
  area       <- basin@data$Area
  region     <- basin@data$Region
  nsub       <- nrow(basin@data)

  # Read and subset input data
  Data        <- read.table(file.path(Location, 'Inputs', Input), sep='\t', header=T)
  Data$DatesR <- as.POSIXct(Data$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y","%Y/%m/%d", "%d-%m-%Y"))
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
  Vect.Parameters     <- as.vector(rep(sort(unique(region)),4))
  IDStable.Parameters <- Vect.Parameters %in% No.Optim
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
  Opt.Parameters.Max <- rep(Parameters.Max, each=length(Opt.Region))
  Opt.Parameters.Log <- rep(c(TRUE, TRUE, FALSE, FALSE), each=length(Opt.Region))

  # Auxiliary functions
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

  # Error function
  #================
  EFUN <- function(Variable){

    # Select model parameters to optimize
    if (is.null(No.Optim)==TRUE){
      All.Parameters <- Variable
    } else{
      New.Parameters <- rbind(cbind(Num.Parameters[!IDStable.Parameters], Variable),
                              cbind(Num.Parameters[IDStable.Parameters], Stable.Parameters))
      All.Parameters <- New.Parameters[match(Num.Parameters, New.Parameters[,1]), 2]
    }

    # Model parameters to run GR2M model
    nreg  <- length(sort(unique(region)))
    Param <- data.frame(Region=sort(unique(region)),
                        X1=All.Parameters[1:nreg],
                        X2=All.Parameters[(nreg+1):(2*nreg)],
                        Fpp=All.Parameters[(2*nreg+1):(3*nreg)],
                        Fpet=All.Parameters[(3*nreg+1):length(All.Parameters)])

    cl=makeCluster(detectCores()-1) # Detect and assign a cluster number
    clusterEvalQ(cl,c(library(GR2MSemiDistr),library(airGR))) # Load package to each node
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

    # Calculate residuals
    mRes <- as.vector(na.omit(Qsim-Qobs))
    return(mRes)

  } # End function

  # Show message
  cat('\f')
  message('Parameter uncertainty with MCMC')
  message('Please wait...')
  msr  <- mean((EFUN(Opt.Parameters))^2)
  MCMC <- modMCMC(f=EFUN, p=Opt.Parameters, lower=Opt.Parameters.Min, upper=Opt.Parameters.Max, niter=Niter, var0=msr)


  # Model function
  #===============
  MFUN <- function(Variable){

    # Select model parameters to optimize
    if (is.null(No.Optim)==TRUE){
      All.Parameters <- Variable
    } else{
      New.Parameters <- rbind(cbind(Num.Parameters[!IDStable.Parameters], Variable),
                              cbind(Num.Parameters[IDStable.Parameters], Stable.Parameters))
      All.Parameters <- New.Parameters[match(Num.Parameters, New.Parameters[,1]), 2]
    }

    # Model parameters to run GR2M model
    nreg  <- length(sort(unique(region)))
    Param <- data.frame(Region=sort(unique(region)),
                        X1=All.Parameters[1:nreg],
                        X2=All.Parameters[(nreg+1):(2*nreg)],
                        Fpp=All.Parameters[(2*nreg+1):(3*nreg)],
                        Fpet=All.Parameters[(3*nreg+1):length(All.Parameters)])

    cl=makeCluster(detectCores()-1) # Detect and assign a cluster number
    clusterEvalQ(cl,c(library(GR2MSemiDistr),library(airGR))) # Load package to each node
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
    Qsim <- qOut[Subset2]
    if (Remove==TRUE){
      Qsim <- Qsim - qSub[Subset2, IdBasin]
    }
    return(Qsim)

  } # End model function

  # Show message
  cat('\f')
  message('Generating streamflow uncertainty')
  message('Please wait...')
  pars <- unique(MCMC$pars)
  sens <- list()
  for (w in 1:nrow(pars)){
    sens[[w]] <- MFUN(pars[w,])
  }
  sR   <- do.call(cbind,sens)
  best <- MFUN(Opt.Parameters)
  min  <- apply(sR, 1, min)
  max  <- apply(sR, 1, max)
  mean <- apply(sR, 1, mean)
  std  <- apply(sR, 1, sd)
  q5   <- apply(sR,1, function(x) quantile(x,0.05))
  q90  <- apply(sR,1, function(x) quantile(x,0.9))
  sensStats  <- data.frame(best=best, min=min, max=max, mean=mean, sd=sd, q5, q90)
  sensOutput <- sR
  sen  <- list(stats=sensStats, out=sensOutput)


  # Create output folder ans save simulation
  #=========================================
  dir.create(file.path(Location, 'Outputs'), recursive=T, showWarnings=F)
  Ans <- list(mcmc=MCMC, sens=sen)
  save(Ans, file=file.path(Location,'Outputs','Uncertainty_GR2MSemiDistr.Rda'))

  # Show message
  message("Done!")
  toc()
  return(Ans)

} # End (not run)
