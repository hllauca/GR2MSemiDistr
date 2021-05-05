#' Run the GR2M model for 'n' subbasins.
#' @param Data        Dataframe with model input's data in airGR format from \code{Create_Forcing_Inputs}.
#' (DatesR, P_1, P_2,..,P_n, E_1, E_2, ...E_n, Q). If Q is not available please provide only DatesR, P, and E.
#' @param Subbasins   Subbasins' shapefile. Must contain the following attributes: 'Area' (in km2), 'Region' (in letters), and 'COMID' (identifier number).
#' @param RunIni      Initial date of the model simulation in 'mm/yyyy' format.
#' @param RunEnd      Ending date of the model simulation in 'mm/yyyy' format.
#' @param WarmUp      Number of months for warm-up. NULL as default.
#' @param Parameters  Vector of model parameters (X1 and X2) and correction factors of P (fp) and E (fpe)
#' in the following order: c(X1, X2, fp, fpe). In the case of existing more than one 'Region'
#' (e.g. regions A and B) please provide model parameters in the following order:
#' c(X1_A, X1_B, X2_A, X2_B, Fp_a, Fp_B, Fpe_A, Fpe_B).
#' @param IniState    Initial states variables. NULL as default.
#' @param Save        Boolean to save results as a text file in the 'Outputs' location. FALSE as default.
#' @param Update      Boolean for the updating mode where only the last month's values will be returned. FALSE as default.
#' @return List of GR2M model outputs.
#' @return PR: Precipitation timeseries for all subbasins in [mm/month].
#' @return AE: Actual evapotranspiration timeseries for all subbasins in [mm/month].
#' @return SM: Soil Moisture timeseries for all subbasins in [mm/month].
#' @return RU: Runoff timeseries for all subbasins in [mm/month].
#' @return QS: Discharge timeseries for all subbasins in [m3/s] (not routed).
#' @return Dates: Vector of dates of the simulation period.
#' @return COMID: Vector of identifier numbers for each subbasin.
#' @return EndState: List of end model states of each subbasin.
#' @return SINK: Basin outlet which contains qsim and qobs data time series in [m3/s].
#' @references Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020). Water. 2021; 13(8):1048. https://doi.org/10.3390/w13081048
#' @export
#' @examples
#' # Run the GR2M model for each subbasin
#' model <- Run_GR2MSemiDistr(Data=data,
#'                            Subbasins=roi,
#'                            RunIni='01/1981',
#'                            RunEnd='12/2016',
#'                            Parameters=c(10.976, 0.665, 1.186, 1.169))
#'
#' # Extract model results
#' View(model$PR) # precipitation [mm/month]
#' View(model$AE) # actual evapotranspiration [mm/month]
#' View(model$SM) # soil moisture [mm/month]
#' View(model$RU) # runoff in [mm/month]
#' print(model$SINK$obs) # observed discharge in [m3/s] at basin outlet
#' print(model$SINK$sim) # simulated discharge in [m3/s] at basin outlet
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  rtop
#' @import  hydroGOF
#' @import  airGR
#' @import  tictoc
#' @import  parallel
#' @import  lubridate
Run_GR2MSemiDistr <- function(Data,
                              Subbasins,
                              RunIni,
                              RunEnd,
                              WarmUp=NULL,
                              Parameters,
                              IniState=NULL,
                              Save=FALSE,
                              Update=FALSE){

  # Load required packages
  require(rgdal)
  require(raster)
  require(rgeos)
  require(rtop)
  require(hydroGOF)
  require(airGR)
  require(tictoc)
  require(parallel)
  require(lubridate)
  tic()

  # Load subbasins data
  area   <- Subbasins@data$Area
  region <- Subbasins@data$Region
  comid  <- as.vector(Subbasins$COMID)
  nsub   <- nrow(Subbasins@data)

  # Subsetting input data
  Data$DatesR <- as.POSIXct(paste0(Data$DatesR,' 00:00:00'),"GMT", tryFormats=c("%Y-%m-%d","%Y/%m/%d","%d-%m-%Y","%d/%m/%Y"))
  Ind_run     <- seq(which(format(Data$DatesR, format="%m/%Y")==RunIni),
                     which(format(Data$DatesR, format="%m/%Y")==RunEnd))
  Database    <- Data[Ind_run,]
  Dates       <- as.Date(Database$DatesR)
  ntime       <- length(Ind_run)

  # Sort model parameters
  Zone  <- sort(unique(region))
  nreg  <- length(Zone)
  Param <- data.frame(Region=sort(unique(region)),
                      X1=Parameters[1:nreg],
                      X2=Parameters[(nreg+1):(2*nreg)],
                      Fpp=Parameters[(2*nreg+1):(3*nreg)],
                      Fpe=Parameters[(3*nreg+1):length(Parameters)])

  # Useful functions
  Subset_Param <- function(Param, Region){
    ParamSub <- c(subset(Param$X1, Param$Region==Region),
                  subset(Param$X2, Param$Region==Region))
    return(ParamSub)
  }
  Forcing_Subbasin <- function(Param, Region, Database, Nsub, ID){
    FactorPP  <- subset(Param$Fpp, Param$Region==Region)
    FactorPET <- subset(Param$Fpe, Param$Region==Region)
    Inputs    <- Database[,c(1,ID+1,ID+1+Nsub)]
    FixInputs <- data.frame(DatesR=Inputs[,1],
                            P=round(FactorPP*Inputs[,2],1),
                            E=round(FactorPET*Inputs[,3],1))
    FixInputs$DatesR <- as.POSIXct(FixInputs$DatesR, "GMT",
                                   tryFormats=c("%Y-%m-%d","%Y/%m/%d",
                                                "%d-%m-%Y","%d/%m/%Y"))
    return(FixInputs)
  }
  numberOfDays <- function(date) {
    m <- format(date, format="%m")
    while (format(date, format="%m")==m) {
      date <- date + 1
    }
    return(as.integer(format(date-1,format="%d")))
  }

  # Number of days in a month (convert mm to m3/s)
  nDays  <- sapply(as.Date(Database$DatesR), numberOfDays)

  # Run hydrological model of each subbasin
    # Show message
    cat('\f')
    message(paste('Running GR2M model for', nsub, 'subbasins'))
    message('Please wait...')

    # Open cluster
    cl=makeCluster(detectCores()-1) # Detect and assign a cluster number
    clusterEvalQ(cl,c(library(GR2MSemiDistr),library(airGR),library(lubridate))) # Load package to each node
    clusterExport(cl,varlist=c("Param","region","nsub","Database",
                               "ntime","IniState","Subset_Param",
                               "Forcing_Subbasin"),envir=environment())

    # Run GR2M
    ResModel <- parLapply(cl, 1:nsub, function(i) {

      # Parameters and factors to run the model
      ParamSub  <- Subset_Param(Param, region[i])
      FixInputs <- Forcing_Subbasin(Param, region[i], Database, nsub, i)
      if(ntime==1){
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
                                       IndPeriod_Run=1:ntime,
                                       verbose=FALSE,
                                       warnings=FALSE)
      } else{
        # Set-up running options
        IniStates <- CreateIniStates(FUN_MOD=RunModel_GR2M,
                                     InputsModel=InputsModel,
                                     ProdStore=IniState[[i]]$Store$Prod,
                                     RoutStore=IniState[[i]]$Store$Rout,
                                     ExpStore=IniState[[i]]$Store$Exp,
                                     UH1=IniState[[i]]$UH$UH1,
                                     UH2=IniState[[i]]$UH$UH2,
                                     GCemaNeigeLayers=NULL,
                                     eTGCemaNeigeLayers=NULL,
                                     GthrCemaNeigeLayers=NULL,
                                     GlocmaxCemaNeigeLayers=NULL)

        RunOptions <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                       InputsModel=InputsModel,
                                       IniStates=IniStates,
                                       IndPeriod_Run=1:ntime,
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


  # Read model results
  if(nsub==1){
    EndState <- list(ResModel[[1]]$StateEnd)
    pr <- round(matrix(ResModel[[1]]$Precip, ncol=1, nrow=ntime),1)
    ae <- round(matrix(ResModel[[1]]$AE, ncol=1, nrow=ntime),1)
    sm <- round(matrix(ResModel[[1]]$Prod, ncol=1, nrow=ntime),1)
    ru <- round(matrix(ResModel[[1]]$Qsim, ncol=1, nrow=ntime),1)
    qs <- round((area[1]*ru)/(86.4*nDays),1)
    qt <- qs
    if(is.null(is.null(Database$Q))==FALSE){
      qo <- Database$Q
    }
    if(is.null(WarmUp)==FALSE){
      pr <- pr[-WarmUp:-1]
      ae <- ae[-WarmUp:-1]
      sm <- sm[-WarmUp:-1]
      ru <- ru[-WarmUp:-1]
      qs <- qs[-WarmUp:-1]
      qt <- qt[-WarmUp:-1]
      if(is.null(is.null(Database$Q))==FALSE){
        qo <- Database$Q[-WarmUp:-1]
      }
      Dates <- Dates[-WarmUp:-1]
    }
  }else{
    EndState <- list()
    PRlist <- list()
    AElist <- list()
    SMlist <- list()
    RUlist <- list()
    QSlist <- list()
    for(w in 1:nsub){
      EndState[[w]] <- ResModel[[w]]$StateEnd
      PRlist[[w]]   <- ResModel[[w]]$Precip
      AElist[[w]]   <- ResModel[[w]]$AE
      SMlist[[w]]   <- ResModel[[w]]$Prod
      RUlist[[w]]   <- ResModel[[w]]$Qsim
      QSlist[[w]]   <- (area[w]*RUlist[[w]])/(86.4*nDays)
    }
    pr <- round(do.call(cbind, PRlist),1)
    ae <- round(do.call(cbind, AElist),1)
    sm <- round(do.call(cbind, SMlist),1)
    ru <- round(do.call(cbind, RUlist),1)
    qs <- round(do.call(cbind, QSlist),1)
    qt <- round(apply(qs,1,FUN=sum),2)
    if(is.null(is.null(Database$Q))==FALSE){
      qo <- Database$Q
    }
    if(is.null(WarmUp)==FALSE){
      pr <- pr[-WarmUp:-1,]
      ae <- ae[-WarmUp:-1,]
      sm <- sm[-WarmUp:-1,]
      ru <- ru[-WarmUp:-1,]
      qs <- qs[-WarmUp:-1,]
      qt <- qt[-WarmUp:-1]
      if(is.null(is.null(Database$Q))==FALSE){
        qo <- Database$Q[-WarmUp:-1]
      }
      Dates <- Dates[-WarmUp:-1]
    }
  }

  if(is.null(Database$Q)==FALSE){
    sink <- data.frame(sim=round(qt,3),obs=round(qo,3))
    rownames(sink) <- Dates
    Ans <- list(SINK=sink,
                QS=qs,
                RU=ru,
                PR=pr,
                AE=ae,
                SM=sm,
                Dates=Dates,
                COMID=comid,
                EndState=EndState)
  }

  if(is.null(Database$Q)==TRUE){
    Ans <- list(QS=qs,
                RU=ru,
                PR=pr,
                AE=ae,
                SM=sm,
                Dates=Dates,
                COMID=comid,
                EndState=EndState)
  }

  # Save model results
  if(Save==TRUE){
    dir.create('./Outputs', recursive=T, showWarnings=F)
    if(Update==TRUE){
      MnYr1    <- format(floor_date(Sys.Date()-months(2), "month"),'%Y%m')
      MnYr2    <- format(floor_date(Sys.Date()-months(1), "month"),'%Y%m')
      PR_name1 <- paste0('PR_GR2MSemiDistr_',MnYr1,'.txt')
      PR_name2 <- paste0('PR_GR2MSemiDistr_',MnYr2,'.txt')
      AE_name1 <- paste0('AE_GR2MSemiDistr_',MnYr1,'.txt')
      AE_name2 <- paste0('AE_GR2MSemiDistr_',MnYr2,'.txt')
      SM_name1 <- paste0('SM_GR2MSemiDistr_',MnYr1,'.txt')
      SM_name2 <- paste0('SM_GR2MSemiDistr_',MnYr2,'.txt')
      RU_name1 <- paste0('RU_GR2MSemiDistr_',MnYr1,'.txt')
      RU_name2 <- paste0('RU_GR2MSemiDistr_',MnYr2,'.txt')
      pr_old   <- read.table(paste0('./Outputs/',PR_name1), header=TRUE, sep='\t')
      ae_old   <- read.table(paste0('./Outputs/',AE_name1), header=TRUE, sep='\t')
      sm_old   <- read.table(paste0('./Outputs/',SM_name1), header=TRUE, sep='\t')
      ru_old   <- read.table(paste0('./Outputs/',RU_name1), header=TRUE, sep='\t')
      pr_new   <- as.data.frame(rbind(as.matrix(pr_old),as.matrix(pr)))
      ae_new   <- as.data.frame(rbind(as.matrix(ae_old),as.matrix(ae)))
      sm_new   <- as.data.frame(rbind(as.matrix(sm_old),as.matrix(sm)))
      ru_new   <- as.data.frame(rbind(as.matrix(ru_old),as.matrix(ru)))
      colnames(pr_new) <- paste0('PR_',comid)
      colnames(ae_new) <- paste0('AE_',comid)
      colnames(sm_new) <- paste0('SM_',comid)
      colnames(ru_new) <- paste0('RU_',comid)
      rownames(pr_new) <- c(as.Date(rownames(pr_old)), Dates)
      rownames(ae_new) <- c(as.Date(rownames(ae_old)), Dates)
      rownames(sm_new) <- c(as.Date(rownames(sm_old)), Dates)
      rownames(ru_new) <- c(as.Date(rownames(ru_old)), Dates)
      write.table(pr_new, file=paste0('./Outputs/',PR_name2), sep='\t')
      write.table(ae_new, file=paste0('./Outputs/',AE_name2), sep='\t')
      write.table(sm_new, file=paste0('./Outputs/',SM_name2), sep='\t')
      write.table(ru_new, file=paste0('./Outputs/',RU_name2), sep='\t')
      file.remove(paste0('./Outputs/',PR_name1))
      file.remove(paste0('./Outputs/',AE_name1))
      file.remove(paste0('./Outputs/',SM_name1))
      file.remove(paste0('./Outputs/',RU_name1))
    } else{
      MnYr    <- format(as.Date(paste0('01/',RunEnd),"%d/%m/%Y"),"%Y%m")
      PR_name <- paste0('PR_GR2MSemiDistr_',MnYr,'.txt')
      AE_name <- paste0('AE_GR2MSemiDistr_',MnYr,'.txt')
      SM_name <- paste0('SM_GR2MSemiDistr_',MnYr,'.txt')
      RU_name <- paste0('RU_GR2MSemiDistr_',MnYr,'.txt')
      pr <- as.data.frame(pr)
      ae <- as.data.frame(ae)
      sm <- as.data.frame(sm)
      ru <- as.data.frame(ru)
      colnames(pr) <- paste0('PR_',comid)
      colnames(ae) <- paste0('AE_',comid)
      colnames(sm) <- paste0('SM_',comid)
      colnames(ru) <- paste0('RU_',comid)
      rownames(pr) <- Dates
      rownames(ae) <- Dates
      rownames(sm) <- Dates
      rownames(ru) <- Dates
      write.table(pr, file=paste0('./Outputs/',PR_name), sep='\t')
      write.table(ae, file=paste0('./Outputs/',AE_name), sep='\t')
      write.table(sm, file=paste0('./Outputs/',SM_name), sep='\t')
      write.table(ru, file=paste0('./Outputs/',RU_name), sep='\t')
    }
  }
  message('Done!')
  toc()
  return(Ans)

} # End (not run)
