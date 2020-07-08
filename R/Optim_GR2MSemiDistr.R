#' Model parameter optimization with SCE-UA algorithm.
#'
#' @param Parameters      GR2M (X1 and X2) model parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Parameters.Min  Minimum GR2M (X1, X2, fprecip and fpet) model parameters values.
#' @param Parameters.Max  Maximum GR2M (X1, X2, fprecip and fpet) model parameters values.
#' @param Max.Functions 	Maximum number of functions used in the optimization loop. 10000 as default.
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
#' @return  Best GR2M model parameters.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  rtop
#' @import  hydroGOF
#' @import  airGR
#' @import  tictoc
#' @import  parallel
#' @import  lubridate
Optim_GR2MSemiDistr <- function(Parameters, Parameters.Min, Parameters.Max, Max.Functions=10000,
                                Optimization='NSE', Location, Shapefile, Input='Inputs_Basins.txt',
									              WarmIni, RunIni, RunEnd, IdBasin, Remove=FALSE, No.Optim=NULL, IniState=NULL){

# Parameters=Model.Param
# Parameters.Min=Model.ParMin
# Parameters.Max=Model.ParMax
# Max.Functions=1500
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
  	require(ProgGUIinR)
  	require(rgdal)
  	require(raster)
  	require(rgeos)
  	require(rtop)
  	require(hydroGOF)
  	require(tictoc)
  	require(parallel)
  	require(lubridate)
    require(airGR)
    tic()

  # Load subbasins
    basin      <- readOGR(file.path(Location,'Inputs', Shapefile), verbose=F)
    area       <- basin@data$Area
    region     <- basin@data$Region
    nsub       <- nrow(basin@data)

  # Input data
    if(is.character(Input)==TRUE){
      Data <- read.table(file.path(Location, 'Inputs', Input), sep='\t', header=T)
    } else{
      Data <- Input
    }
    Data$DatesR <- as.POSIXct(paste0(Data$DatesR,' 00:00:00'), "GMT", tryFormats=c("%Y-%m-%d", "%Y/%m/%d", "%d-%m-%Y", "%d/%m/%Y"))
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

  # Subset calibration regions and parameters not to be optimized
    n.param    <- 1:length(Parameters)
    v.param    <- as.vector(rep(sort(unique(region)),4))
	  id.param   <- v.param %in% No.Optim
    nopt.param <- Parameters[id.param]

  # Define calibration regions and parameters ranges to optimize
    if(is.null(No.Optim)==TRUE){
        opt.param  <- Parameters
        opt.region <- unique(region)
     }else{
        opt.param  <- Parameters[!id.param]
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

  # Objective function
    OFUN <- function(Variable){

		  # Select model parameters to optimize
			if(is.null(No.Optim)==TRUE){
				all.param <- Variable
			}else{
				new.param <- rbind(cbind(n.param[!id.param], Variable),
								     cbind(n.param[id.param], nopt.param))
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
			clusterEvalQ(cl,c(library(GR2MSemiDistr),library(airGR),library(lubridate))) # Load package to each node
			clusterExport(cl,varlist=c("Param","region","nsub","Database","time","IniState","Subset_Param","Forcing_Subbasin"),envir=environment())

		  # Run GR2M
			ResModel <- parLapply(cl, 1:nsub, function(i) {

                # Parameters and factors to run the model
                  ParamSub  <- Subset_Param(Param, region[i])
                  FixInputs <- Forcing_Subbasin(Param, region[i], Database, nsub, i)
                  if(time==1){
                    NewDate   <- as.POSIXct(floor_date(FixInputs$DatesR+months(1),"month"))
                    NewStep   <- data.frame(DatesR=NewDate, P=100, E=100)
                    FixInputs <- rbind(FixInputs,NewStep)
                  }

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
                  }else{
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
    Subset2     <- seq(which(format(Database$DatesR, format="%m/%Y") == RunIni),
                       which(format(Database$DatesR, format="%m/%Y") == RunEnd))
    Database2   <- Database[Subset2,]
    Qobs        <- Database2$Qm3s
    Qsim        <- qOut[Subset2]
    if(Remove==TRUE){
       Qsim <- Qsim - qSub[Subset2, IdBasin]
    }

  # Evaluation criteria
  	res.df   <- na.omit(data.frame(sim=Qsim, obs=Qobs))
    optim.df <- data.frame(KGE=1-round(KGE(res.df$sim, res.df$obs), 3),
                           NSE=1-round(NSE(res.df$sim, res.df$obs), 3),
                           lnNSE=1-round(NSE(log(res.df$sim), log(res.df$obs)), 3),
                           R=1-round(rPearson(res.df$sim, res.df$obs), 3),
                           RMSE=round(rmse(res.df$sim, res.df$obs), 3),
                           PBIAS=round(pbias(res.df$sim, res.df$obs), 3))

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
                         maxn=Max.Functions,
                         kstop=3,
                         pcento=0.1)

  # Extracting calibration results
    if(Optimization == 'PBIAS' | Optimization == 'RMSE'){
      fo <- round(Calibration$value,3)
    }else{
      fo <- round(1-Calibration$value,3)
    }
    if(is.null(No.Optim)==TRUE){
      all.param <- round(Calibration$par,3)
    }else{
      new.param <- rbind(cbind(n.param[!id.param], round(Calibration$par,3)),
                         cbind(n.param[id.param], nopt.param))
      all.param <- new.param[match(n.param, new.param[,1]), 2]
    }
    Ans <- list(Param=all.param, Value=fo)

  # Create output folder and save simulation
    dir.create(file.path(Location, 'Outputs'), recursive=T, showWarnings=F)
    save(Ans, file=file.path(Location,'Outputs','Optimization_GR2MSemiDistr.Rda'))

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
