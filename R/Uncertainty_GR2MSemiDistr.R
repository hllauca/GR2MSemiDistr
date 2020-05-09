#' Uncertainty analysis of GR2M model parameters with the MCMC algorithm.
#'
#' @param Parameters      GR2M (X1 and X2) model parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Parameters.Min  Minimum GR2M (X1, X2, fprecip and fpet) model parameters values.
#' @param Parameters.Max  Maximum GR2M (X1, X2, fprecip and fpet) model parameters values.
#' @param Niter 	        Number of iterations. 1000 as default.
#' @param Location    Directory where 'Inputs' folder is located.
#' @param Shapefile   Subbasin shapefile.
#' @param Input       Forcing data texfile (Dates, Precip, PotEvap, Qobs). 'Inputs_Basins.txt' as default.
#' @param WarmIni     Initial date (in 'mm/yyyy' format) of the warm-up period.
#' @param RunIni      Initial date (in 'mm/yyyy' format) of the model simulation period.
#' @param RunEnd      Final date (in 'mm/yyyy' format) of the model simulation period.
#' @param IdBasin     ID for the outlet subbasin (from shapefile attribute table).
#' @param Remove      Logical value to remove streamflows of the outlet subbasin (IdBasin). FALSE as default.
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @return  Parameter and streamflow uncertanty bounds.
#' @export
#' @import  ProgGUIinR
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  FME
#' @import  parallel
#' @import  tictoc
#' @import  airGR
Uncertainty_GR2MSemiDistr <- function(Parameters, Parameters.Min, Parameters.Max, Niter=1000,
                                      Location, Shapefile, Input='Inputs_Basins.txt', WarmIni,
                                      RunIni, RunEnd, IdBasin, Remove=FALSE, IniState=NULL){

  # Parameters=Model.Param
  # Parameters.Min=Model.ParMin
  # Parameters.Max=Model.ParMax
  # Niter=1000
  # Location=Location
  # Shapefile=File.Shape
  # Input='Inputs_Basins.txt'
  # WarmIni=WarmUp.Ini
  # RunIni=RunModel.Ini
  # RunEnd=RunModel.End
  # IdBasin=Optim.Basin
  # Remove=Optim.Remove
  # IniState=NULL


  # Load packages
  require(rgdal)
  require(raster)
  require(rgeos)
  require(ProgGUIinR)
  require(parallel)
  require(tictoc)
  require(airGR)
  require(FME)
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


  # Define calibration regions and parameters ranges to optimize
  Opt.Region     <- unique(region)
  Parameters.Min <- rep(Parameters.Min, each=length(Opt.Region))
  Parameters.Max <- rep(Parameters.Max, each=length(Opt.Region))
  Parameters.Log <- rep(c(TRUE, TRUE, FALSE, FALSE), each=length(Opt.Region))


  # Useful functions
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


  # Residual function
  RFUN <- function(Variable){

    # Model parameters
    nreg  <- length(sort(unique(region)))
    Param <- data.frame(Region=sort(unique(region)),
                        X1=Variable[1:nreg],
                        X2=Variable[(nreg+1):(2*nreg)],
                        Fpp=Variable[(2*nreg+1):(3*nreg)],
                        Fpet=Variable[(3*nreg+1):length(Variable)])

    # Open cluster
    cl=makeCluster(detectCores()-1) # Detect and assign a cluster number
    clusterEvalQ(cl,c(library(GR2MSemiDistr),library(airGR))) # Load package to each node
    clusterExport(cl,varlist=c("Param","region","nsub","Database","time",
                               "IniState","Subset_Param","Forcing_Subbasin"),envir=environment())

    # Run GR2M
    ResModel <- parLapply(cl, 1:nsub, function(i) {

      ParamSub  <- Subset_Param(Param, region[i])
      FixInputs <- Forcing_Subbasin(Param, region[i], Database, nsub, i)

      InputsModel <- CreateInputsModel(FUN_MOD=RunModel_GR2M,
                                       DatesR=FixInputs$DatesR,
                                       Precip=FixInputs$P,
                                       PotEvap=FixInputs$E)

      if(is.null(IniState)==TRUE){
        RunOptions <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                       InputsModel=InputsModel,
                                       IndPeriod_Run=1:time,
                                       verbose=FALSE,
                                       warnings=FALSE)
      } else{
        RunOptions <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                       InputsModel=InputsModel,
                                       IniStates=IniState[[i]],
                                       IndPeriod_Run=1:time,
                                       verbose=FALSE,
                                       warnings=FALSE)
      }

      OutputsModel <- RunModel(InputsModel=InputsModel,
                               RunOptions=RunOptions,
                               Param=ParamSub,
                               FUN=RunModel_GR2M)

      return(OutputsModel)
    })
    stopCluster(cl)

    # Prepare model results
    if (nsub==1){
      qSub <- (area[1]*ResModel[[1]]$Qsim)/(86.4*nDays)
      qOut <- qSub
    } else{
      Qlist <- list()
      for(w in 1:nsub){Qlist[[w]] <- (area[w]*ResModel[[w]]$Qsim)/(86.4*nDays)}
      qSub <- do.call(cbind, Qlist)
      qOut <- round(apply(qSub, 1, FUN=sum),2)
    }

    # Subset model results (exclude warm-up)
    Subset2     <- seq(which(format(Database$DatesR, format="%m/%Y") == RunIni),
                       which(format(Database$DatesR, format="%m/%Y") == RunEnd))
    Database2   <- Database[Subset2,]

    # Calculate model residuals
    Qobs <- Database2$Qm3s
    Qsim <- qOut[Subset2]
    if (Remove==TRUE){
      Qsim <- Qsim - qSub[Subset2, IdBasin]
    }
    mRes <- as.vector(na.omit(Qsim-Qobs))
    return(mRes)

  } # End function


  # Show message
  cat('\f')
  message('Parameter uncertainty with MCMC')
  message('Please wait...')
  msr  <- mean((RFUN(Parameters))^2)
  MCMC <- modMCMC(f=RFUN, p=Parameters, lower=Parameters.Min, upper=Parameters.Max, niter=Niter, var0=msr)


  # Model function
  MFUN <- function(Variable){

    # Model parameters
    nreg  <- length(sort(unique(region)))
    Param <- data.frame(Region=sort(unique(region)),
                        X1=Variable[1:nreg],
                        X2=Variable[(nreg+1):(2*nreg)],
                        Fpp=Variable[(2*nreg+1):(3*nreg)],
                        Fpet=Variable[(3*nreg+1):length(Variable)])

    # Open cluster
    cl=makeCluster(detectCores()-1) # Detect and assign a cluster number
    clusterEvalQ(cl,c(library(GR2MSemiDistr),library(airGR))) # Load package to each node
    clusterExport(cl,varlist=c("Param","region","nsub","Database","time",
                               "IniState","Subset_Param","Forcing_Subbasin"),envir=environment())

    # Run GR2M
    ResModel <- parLapply(cl, 1:nsub, function(i) {

      ParamSub  <- Subset_Param(Param, region[i])
      FixInputs <- Forcing_Subbasin(Param, region[i], Database, nsub, i)

      InputsModel <- CreateInputsModel(FUN_MOD=RunModel_GR2M,
                                       DatesR=FixInputs$DatesR,
                                       Precip=FixInputs$P,
                                       PotEvap=FixInputs$E)

      if(is.null(IniState)==TRUE){
        RunOptions <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                       InputsModel=InputsModel,
                                       IndPeriod_Run=1:time,
                                       verbose=FALSE,
                                       warnings=FALSE)
      } else{
        RunOptions <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                       InputsModel=InputsModel,
                                       IniStates=IniState[[i]],
                                       IndPeriod_Run=1:time,
                                       verbose=FALSE,
                                       warnings=FALSE)
      }

      OutputsModel <- RunModel(InputsModel=InputsModel,
                               RunOptions=RunOptions,
                               Param=ParamSub,
                               FUN=RunModel_GR2M)

      return(OutputsModel)

    })
    stopCluster(cl)

    # Main model results (Qsim in m3/s and EndState variables)
    if (nsub==1){
      qSub <- (area[1]*ResModel[[1]]$Qsim)/(86.4*nDays)
      qOut <- qSub
    } else{
      Qlist <- list()
      for(w in 1:nsub){Qlist[[w]] <- (area[w]*ResModel[[w]]$Qsim)/(86.4*nDays)}
      qSub <- do.call(cbind, Qlist)
      qOut <- round(apply(qSub, 1, FUN=sum),2)
    }

    # Subset model results (exclude warm-up)
    Subset2     <- seq(which(format(Database$DatesR, format="%m/%Y") == RunIni),
                       which(format(Database$DatesR, format="%m/%Y") == RunEnd))
    Database2   <- Database[Subset2,]

    # Streamflow at the basin outlet
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
  best <- MFUN(Parameters)
  min  <- apply(sR, 1, min)
  max  <- apply(sR, 1, max)
  mean <- apply(sR, 1, mean)
  std  <- apply(sR, 1, sd)
  q5   <- apply(sR,1, function(x) quantile(x,0.05))
  q90  <- apply(sR,1, function(x) quantile(x,0.9))
  sensStats  <- data.frame(best=best, min=min, max=max, mean=mean, std=std, q5=q5, q90=q90)
  sensOutput <- sR
  sen  <- list(stats=sensStats, out=sensOutput)


  # Create output folder and save results
  dir.create(file.path(Location, 'Outputs'), recursive=T, showWarnings=F)
  Ans <- list(mcmc=MCMC, sens=sen)
  save(Ans, file=file.path(Location,'Outputs','Uncertainty_GR2MSemiDistr.Rda'))

  # Show message
  message("Done!")
  toc()
  return(Ans)

} # End (not run)
