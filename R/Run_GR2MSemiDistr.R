#' Run the GR2M model for each subbasins.
#'
#' @param Data        Database in airGR format (DatesR,P,E,Q).
#' @param Subbasins   Subbasins shapefile.
#' @param RunIni      Initial date for model simulation (in mm/yyyy format).
#' @param RunEnd      Final date for model simulation (in mm/yyyy format).
#' @param WarmUp      Number of months for warm-up. NULL as default.
#' @param Parameters  Model parameters and correction factor of P and E.
#' @param IniState    Initial states variables. NULL as default.
#' @param Regional    Boolean to simulate in a regional mode (more than one hydrological station). FALSE as default.
#' @param Save        Boolean to save outputs as text files. FALSE as default.
#' @param Update      Bollean to update previous outputs text files. FALSE as default.
#' @return GR2M model outputs for each subbasin.
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
Run_GR2MSemiDistr <- function(Data,
                              Subbasins,
                              RunIni,
                              RunEnd,
                              WarmUp=NULL,
                              Parameters,
                              IniState=NULL,
                              Regional=FALSE,
                              Save=FALSE,
                              Update=FALSE){

  # Data=Data
  # Subbasins=roi
  # RunIni=RunModel.Ini
  # RunEnd=RunModel.End
  # WarmUp=NULL
  # Parameters=Model.Param
  # IniState=NULL
  # Regional=FALSE
  # Save=TRUE
  # Update=FALSE

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
  nDays <- c()
  for(i in 1:ntime){
    nDays[i] <- numberOfDays(as.Date(Database$DatesR)[i])
  }

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


  # Read model results for each subbasin
  if(nsub==1){
    EndState <- list(ResModel[[1]]$StateEnd)
    pr <- round(matrix(ResModel[[1]]$Precip, ncol=1, nrow=ntime),1)
    ae <- round(matrix(ResModel[[1]]$AE, ncol=1, nrow=ntime),1)
    sm <- round(matrix(ResModel[[1]]$Prod, ncol=1, nrow=ntime),1)
    ru <- round(matrix(ResModel[[1]]$Qsim, ncol=1, nrow=ntime),1)
    qs <- round((area[1]*ru)/(86.4*nDays),1)
    qt <- qs
    if(is.null(WarmUp)==FALSE){
      pr <- pr[-WarmUp:-1]
      ae <- ae[-WarmUp:-1]
      sm <- sm[-WarmUp:-1]
      ru <- ru[-WarmUp:-1]
      qs <- qs[-WarmUp:-1]
      qt <- qt[-WarmUp:-1]
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
    if(is.null(WarmUp)==FALSE){
      pr <- pr[-WarmUp:-1,]
      ae <- ae[-WarmUp:-1,]
      sm <- sm[-WarmUp:-1,]
      ru <- ru[-WarmUp:-1,]
      qs <- qs[-WarmUp:-1,]
      qt <- qt[-WarmUp:-1]
      Dates <- Dates[-WarmUp:-1]
    }
  }

  # Local mode - simulation with a unique hydrological station
  if(Regional==FALSE){
    sink <- data.frame(sim=round(qt,3),obs=round(Database$Q,3))
    if(is.null(WarmUp)==FALSE){
      sink <- sink[-WarmUp:-1,]
    }
    rownames(Qout) <- Dates
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

  # Regional mode - simulation with more than one hydrological station
  if(Regional==TRUE){
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
      colnames(pr_new) <- comid
      colnames(ae_new) <- comid
      colnames(sm_new) <- comid
      colnames(ru_new) <- comid
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
      colnames(pr) <- comid
      colnames(ae) <- comid
      colnames(sm) <- comid
      colnames(ru) <- comid
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
