#' Model parameter optimization with the SCE-UA algorithm.
#' @param Data        Dataframe with model input's data in airGR format from \code{Create_Forcing_Inputs}.
#' (DatesR, P_1, P_2,..,P_n, E_1, E_2, ...E_n, Q). If Q is not available please provide only DatesR, P, and E.
#' @param Subbasins   Subbasins' shapefile. Must contain the following attributes: 'Area' (in km2), 'Region' (in letters), and 'COMID' (identifier number).
#' @param RunIni      Initial date of the model simulation in 'mm/yyyy' format.
#' @param RunEnd      Ending date of the model simulation in 'mm/yyyy' format.
#' @param WarmUp      Number of months for the warm-up period. NULL as default.
#' @param Parameters  Vector of initial model parameters (X1 and X2) and correction factors of P (fp) and E (fpe)
#' in the following order: c(X1, X2, fp, fpe). In the case of existing more than one 'Region'
#' (e.g. regions A and B) please provide model parameters in the following order:
#' c(X1_A, X1_B, X2_A, X2_B, Fp_a, Fp_B, Fpe_A, Fpe_B).
#' @param Parameters.Min  Vector of minimum values of GR2M model parameters and correction factors in the following order: c(X1_min, X2_min, fp_min, fpe_min).
#' @param Parameters.Max  Vector of maximum values of GR2M model parameters and correction factors in the following order: c(X1_max, X2_max, fp_max, fpe_max).
#' @param Max.Functions 	Maximum number of function evaluation for optimization. 5000 as default.
#' @param Optimization    Objective function for optimization (NSE, KGE, or RMSE).
#' @param No.Optim    Regions not to be optimized. NULL as default.
#' @return  List of optimal GR2M model parameters for each 'Region'.
#' @return  Param: Best set of GR2M model parameters (sorted by 'Region').
#' @return  Value: Final value of the objective function.
#' @references Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020). Water. 2021; 13(8):1048. https://doi.org/10.3390/w13081048
#' @export
#' @examples
#' # Optimize GR2M model parameters for a single 'Region' using the KGE metric
#' optim <- Optim_GR2MSemiDistr(Data=data,
#'                              Subbasins=roi,
#'                              RunIni='01/1981',
#'                              RunEnd='12/2002',
#'                              WarmUp=36,
#'                              Parameters=c(1000, 1, 1, 1),
#'                              Parameters.Min=c(1, 0.01, 0.8, 0.8),
#'                              Parameters.Max=c(2000, 2, 1.2, 1.2),
#'                              Max.Functions=1000,
#'                              Optimization='KGE')
#'  best_param <- optim$Param
#' @import  raster
#' @import  rtop
#' @import  hydroGOF
#' @import  airGR
#' @import  tictoc
#' @import  parallel
Optim_GR2MSemiDistr <- function(Data,
                                Subbasins,
                                RunIni,
                                RunEnd,
                                WarmUp=NULL,
                                Parameters,
                                Parameters.Min,
                                Parameters.Max,
                                Max.Functions=5000,
                                Optimization='NSE',
                                No.Optim=NULL){


  # Load packages
  require(raster)
  require(rtop)
  require(hydroGOF)
  require(airGR)
  require(tictoc)
  require(parallel)
  tic()

  # Load subbasins data
  area    <- Subbasins@data$Area
  region  <- Subbasins@data$Region
  nsub    <- nrow(Subbasins@data)

  # Input data
  Data$DatesR <- as.POSIXct(paste0(Data$DatesR,' 00:00:00'),"GMT",
                            tryFormats=c("%Y-%m-%d","%Y/%m/%d",
                                         "%d-%m-%Y","%d/%m/%Y"))
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
    ParamSub <- c(subset(Param$X1, Param$Region==Region),
                  subset(Param$X2, Param$Region==Region))
    return(ParamSub)
  }
  Forcing_Subbasin <- function(Param, Region, Database, Nsub, ID){
    FactorPP  <- subset(Param$Fpp, Param$Region==Region)
    FactorPET <- subset(Param$Fpet, Param$Region==Region)
    Inputs    <- Database[,c(1,ID+1,ID+1+Nsub)]
    FixInputs <- data.frame(DatesR=Inputs[,1],
                            P=round(FactorPP*Inputs[,2],1),
                            E=round(FactorPET*Inputs[,3],1))
    FixInputs$DatesR <- as.POSIXct(FixInputs$DatesR, "GMT",
                                   tryFormats=c("%Y-%m-%d","%Y/%m/%d",
                                                "%d-%m-%Y", "%d/%m/%Y"))
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
  OFUN <- function(parameter_set){

    # Select model parameters to optimize
    if(is.null(No.Optim)==TRUE){
      all.param <- parameter_set
    }else{
      new.param <- rbind(cbind(n.param[!id.nopt], parameter_set),
                         cbind(n.param[id.nopt], nopt.param))
      all.param <- new.param[match(n.param, new.param[,1]), 2]
    }

    # GR2M model parameters
    Zone  <- sort(unique(region))
    nreg  <- length(Zone)
    Param <- data.frame(Region=sort(unique(region)),
                        X1=all.param[1:nreg],
                        X2=all.param[(nreg+1):(2*nreg)],
                        Fpp=all.param[(2*nreg+1):(3*nreg)],
                        Fpet=all.param[(3*nreg+1):length(all.param)])

    # Open cluster
    cl=makeCluster(detectCores()-1) # Detect and assign a cluster number
    clusterEvalQ(cl,c(library(GR2MSemiDistr),library(airGR)))
    clusterExport(cl,varlist=c("Param","region","nsub",
                               "Database","time","Subset_Param",
                               "Forcing_Subbasin"),envir=environment())

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
      qs   <- (area[1]*ResModel[[1]]$Qsim)/(86.4*nDays)
      sink <- qs
    }else{
      QSlist <- list()
      for(w in 1:nsub){
        QSlist[[w]] <- (area[w]*ResModel[[w]]$Qsim)/(86.4*nDays)
      }
      qs   <- do.call(cbind, QSlist)
      sink <- round(apply(qs, 1, FUN=sum),2)
    }

    # Subset model results (exclude warm-up)
    if(is.null(WarmUp)==TRUE){
      Qobs  <- Database$Q
      Qsim  <- sink
    }else{
      Qobs  <- Database$Q[-WarmUp:-1]
      Qsim  <- sink[-WarmUp:-1]
    }


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

  # Print results
  print("Optimization results:")
  print("======================")
  print(paste0(rep(c('X1=','X2=', 'fpr=', 'fpe='), each=length(Ans$Param)/4), Ans$Param))
  print(paste0(Optimization,'=', Ans$Value))

  # Show message
  message("Done!")
  toc()

  # Output
  return(Ans)

} # End (not run)
