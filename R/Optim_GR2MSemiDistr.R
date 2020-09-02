#' Model parameter optimization with SCE-UA algorithm.
#'
#' @param Data        File with input data in airGR format (DatesR,P,E,Q)
#' @param Subbasins   Subbasins shapefile.
#' @param RunIni      Initial date of model simulation (in mm/yyyy format).
#' @param RunEnd      Final date of model simulation (in mm/yyyy format).
#' @param WarmUp      Number of months for warm-up.
#' @param Parameters      GR2M model parameters and correction factor of P and E.
#' @param Parameters.Min  Minimum values of GR2M model parameters and correction factor of P and E.
#' @param Parameters.Max  Maximum values of GR2M model parameters and correction factor of P and E.
#' @param Max.Functions 	Maximum number of function evaluation for optimization. 5000 as default.
#' @param Optimization    Objective function (NSE, KGE, RMSE).
#' @param No.Optim    Regions not to be optimized. NULL as default.
#' @return  Optimal GR2M model parameters.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  rtop
#' @import  hydroGOF
#' @import  airGR
#' @import  tictoc
#' @import  parallel
Optim_GR2MSemiDistr <- function(Data,
                                Subbasin,
                                RunIni,
                                RunEnd,
                                WarmUp,
                                Parameters,
                                Parameters.Min,
                                Parameters.Max,
                                Max.Functions=5000,
                                Optimization='NSE',
                                No.Optim=NULL){

  # Data=Data
  # Subbasin=roi
  # RunIni=RunModel.Ini
  # RunEnd=RunModel.End
  # WarmUp=WarmUp
  # Parameters=Model.Param
  # Parameters.Min=Model.Param.Min
  # Parameters.Max=Model.Param.Max
  # Max.Functions=5
  # Optimization='KGE'
  # No.Optim=NULL

  # Load packages
  require(rgdal)
  require(raster)
  require(rgeos)
  require(rtop)
  require(hydroGOF)
  require(airGR)
  require(tictoc)
  require(parallel)
  tic()

  # Load subbasins
  area    <- Subbasin@data$Area
  region  <- Subbasin@data$Region
  nsub    <- nrow(Subbasin@data)

  # Input data
  Data$DatesR <- as.POSIXct(paste0(Data$DatesR,' 00:00:00'),"GMT",
                            tryFormats=c("%Y-%m-%d","%Y/%m/%d","%d-%m-%Y","%d/%m/%Y"))
  Ind_run     <- seq(which(format(Data$DatesR, format="%m/%Y")==RunIni),
                     which(format(Data$DatesR, format="%m/%Y")==RunEnd))
  Database    <- Data[Ind_run,]
  time        <- length(Ind_run)

  # Define calibration regions and parameters ranges to optimize
  if(is.null(No.Optim)==TRUE){
    opt.param  <- Parameters
    opt.region <- unique(region)
  }else{
    n.param    <- 1:length(Parameters)
    v.param    <- as.vector(rep(sort(unique(region)),4))
    id.nopt    <- v.param %in% No.Optim
    nopt.param <- Parameters[id.nopt]
    opt.param  <- Parameters[!id.nopt]
    opt.region <- unique(region[!(region %in% No.Optim)])
  }
  opt.param.min <- rep(Parameters.Min, each=length(opt.region))
  opt.param.max <- rep(Parameters.Max, each=length(opt.region))

  # Useful functions
  Subset_Param <- function(Param, Region){
    ParamSub <- c(subset(Param$X1, Param$Region==Region), subset(Param$X2, Param$Region==Region))
    return(ParamSub)
  }

  Forcing_Subbasin <- function(Param, Region, Database, Nsub, ID){
    FactorPP  <- subset(Param$Fpp, Param$Region==Region)
    FactorPET <- subset(Param$Fpet, Param$Region==Region)
    Inputs    <- Database[,c(1,ID+1,ID+1+Nsub)]
    FixInputs <- data.frame(DatesR=Inputs[,1], P=round(FactorPP*Inputs[,2],1), E=round(FactorPET*Inputs[,3],1))
    FixInputs$DatesR <- as.POSIXct(FixInputs$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%Y/%m/%d", "%d-%m-%Y", "%d/%m/%Y"))
    return(FixInputs)
  }

  numberOfDays <- function(date) {
    m <- format(date, format="%m")
    while (format(date, format="%m")==m) {
      date <- date + 1
    }
    return(as.integer(format(date-1,format="%d")))
  }

  # Number of days in a month (to convert mm to m3/s)
  nDays <- c()
  for(i in 1:time){
    nDays[i] <- numberOfDays(as.Date(Database$DatesR)[i])
  }

  # Objective function
  OFUN <- function(pars){

    # Select model parameters to optimize
    if(is.null(No.Optim)==TRUE){
      all.param <- pars
    }else{
      new.param <- rbind(cbind(n.param[!id.nopt], pars),
                         cbind(n.param[id.nopt], nopt.param))
      all.param <- new.param[match(n.param, new.param[,1]), 2]
    }

    # GR2M model parameters
    Zone  <- sort(unique(region))
    nreg  <- length(Zone)
    Param <- data.frame(Region=sort(unique(region)),
                        X1=Parameters[1:nreg],
                        X2=Parameters[(nreg+1):(2*nreg)],
                        Fpp=Parameters[(2*nreg+1):(3*nreg)],
                        Fpet=Parameters[(3*nreg+1):length(all.param)])

    # Open cluster
    cl=makeCluster(detectCores()-1) # Detect and assign a cluster number
    clusterEvalQ(cl,c(library(GR2MSemiDistr),library(airGR)))
    clusterExport(cl,varlist=c("Param","region","nsub","Database","time",
                               "Subset_Param","Forcing_Subbasin"),envir=environment())

    # Run GR2M
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
      # Set-up running options
      RunOptions <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                     InputsModel=InputsModel,
                                     IndPeriod_Run=1:time,
                                     verbose=FALSE,
                                     warnings=FALSE)

      # Run GR2M
      OutputsModel <- RunModel(InputsModel=InputsModel,
                               RunOptions=RunOptions,
                               Param=ParamSub,
                               FUN_MOD=RunModel_GR2M)

      return(OutputsModel)
    })

    # Close cluster
    stopCluster(cl)

    # Model results
    if(nsub==1){
      qSub <- (area[1]*ResModel[[1]]$Qsim)/(86.4*nDays)
      qOut <- qSub
    }else{
      Qlist <- list()
      for(w in 1:nsub){
        Qlist[[w]] <- (area[w]*ResModel[[w]]$Qsim)/(86.4*nDays)
      }
      qSub <- do.call(cbind, Qlist)
      qOut <- round(apply(qSub, 1, FUN=sum),2)
    }

    # Subset model results (exclude warm-up)
    Qobs        <- Database$Q[-WarmUp:-1]
    Qsim        <- qOut[-WarmUp:-1]

    # Evaluation criteria
    res.df   <- na.omit(data.frame(sim=Qsim, obs=Qobs))
    optim.df <- data.frame(KGE=1-round(KGE(res.df$sim, res.df$obs),3),
                           NSE=1-round(NSE(res.df$sim, res.df$obs),3),
                           RMSE=round(rmse(res.df$sim, res.df$obs), 3))

    # Return
    OF <- as.numeric(optim.df[colnames(optim.df) %in% Optimization])
    return(OF)
  }# End


  # Run optimization
  # Show message
  cat('\f')
  message(paste('Optimizing', Optimization, 'with SCE-UA'))
  message('Please wait...')
  Calibration <- sceua(OFUN,
                       pars=opt.param,
                       lower=opt.param.min,
                       upper=opt.param.max,
                       maxn=Max.Functions)

  # Extracting calibration results
  if(Optimization=='RMSE'){
    fo <- round(Calibration$value,3)
  }else{
    fo <- round(1-Calibration$value,3)
  }
  if(is.null(No.Optim)==TRUE){
    parameter <- round(Calibration$par,3)
  }else{
    order     <- rbind(cbind(n.param[!id.nopt], round(Calibration$par,3)),
                       cbind(n.param[id.nopt], nopt.param))
    parameter <- order[match(n.param, order[,1]), 2]
  }
  Ans <- list(Param=parameter, Value=fo)

  # Create output folder and save simulation
  dir.create(file.path(getwd(), 'Outputs'), recursive=T, showWarnings=F)
  save(Ans, file=file.path(getwd(),'Outputs','Optimization_GR2MSemiDistr.Rda'))

  # Print results
  print("Optimization results:")
  print("======================")
  print(paste0(rep(c('X1=','X2=', 'fprecip=', 'fpet='), each=length(Ans$Param)/4), Ans$Param))
  print(paste0(Optim.Eval,'=', Ans$Value))

  # Show message
  message("Done!")
  toc()

  # Output
  return(Ans)

} # End (not run)
