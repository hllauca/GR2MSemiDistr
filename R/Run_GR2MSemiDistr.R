#' Run the GR2M model for each subbasins.
#'
#' @param Data        File with input data in airGR format (DatesR,P,E,Q)
#' @param Subbasins   Subbasins shapefile.
#' @param RunIni      Initial date of model simulation (in mm/yyyy format).
#' @param RunEnd      Final date of model simulation (in mm/yyyy format).
#' @param Parameters  GR2M model parameters and correction factor of P and E.
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @param Regional    Logical value to simulate in a regional mode (more than one outlet). FALSE as default.
#' @param Save        Logical valute to export simulation results as '.Rda'. TRUE as default.
#' @param Update      Logical value to update a previous outputs text files. FALSE as default.
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
                              Parameters,
                              IniState=NULL,
                              Regional=FALSE,
                              Save=FALSE,
                              Update=FALSE){

  # Data=Data
  # Subbasins=roi
  # RunIni=RunModel.Ini
  # RunEnd=RunModel.End
  # Parameters=Model.Param
  # IniState=NULL
  # Regional=FALSE
  # Save=TRUE
  # Update=FALSE

  # Load packages
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

  # Load Subbasins data
  area    <- Subbasins@data$Area
  region  <- Subbasins@data$Region
  sub.id  <- paste0('GR2M_ID_',as.vector(Subbasins$GR2M_ID))
  nsub    <- nrow(Subbasins@data)

  # Input data
  Data$DatesR <- as.POSIXct(paste0(Data$DatesR,' 00:00:00'),"GMT",
                            tryFormats=c("%Y-%m-%d","%Y/%m/%d",
                                         "%d-%m-%Y","%d/%m/%Y"))
  Ind_run     <- seq(which(format(Data$DatesR, format="%m/%Y")==RunIni),
                     which(format(Data$DatesR, format="%m/%Y")==RunEnd))
  Database    <- Data[Ind_run,]
  ntime       <- length(Ind_run)
  Dates       <- as.Date(Database$DatesR)

  # GR2M model parameters
  Zone  <- sort(unique(region))
  nreg  <- length(Zone)
  Param <- data.frame(Region=sort(unique(region)),
                      X1=Parameters[1:nreg],
                      X2=Parameters[(nreg+1):(2*nreg)],
                      Fpp=Parameters[(2*nreg+1):(3*nreg)],
                      Fpet=Parameters[(3*nreg+1):length(Parameters)])

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

  # Number of days in a month (to convert mm to m3/s)
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
      RunOptions <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                     InputsModel=InputsModel,
                                     IniStates=IniState[[i]],
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

  # Model results
  if(nsub==1){
    prod <- ResModel[[1]]$Prod
    qsub <- (area[1]*ResModel[[1]]$Qsim)/(86.4*nDays)
    qout <- qsub
    EndState <- list(ResModel[[1]]$StateEnd)
  }else{
    Plist <- list()
    Qlist <- list()
    for(w in 1:nsub){
      Plist[[w]] <- ResModel[[w]]$Prod
      Qlist[[w]] <- (area[w]*ResModel[[w]]$Qsim)/(86.4*nDays)
    }
    prod <- do.call(cbind, Plist)
    qsub <- do.call(cbind, Qlist)
    qout <- round(apply(qsub, 1, FUN=sum),2)
    EndState <- list()
    for(w in 1:nsub){
      EndState[[w]] <- ResModel[[w]]$StateEnd
    }
  }

  # Subset model results (exclude warm-up)
  if(nsub==1){
    Qsub <- as.data.frame(matrix(round(qsub,3), ncol=1, nrow=length(qsub)))
    Prod <- as.data.frame(matrix(round(prod,3), ncol=1, nrow=length(prod)))
  }else{
    Qsub <- as.data.frame(round(qsub,3))
    Prod <- as.data.frame(round(prod,3))
  }
  colnames(Qsub) <- sub.id
  rownames(Qsub) <- Dates
  colnames(Prod) <- sub.id
  rownames(Prod) <- Dates

  # Local mode
  if(Regional==FALSE){
    Qout <- data.frame(sim=round(qout,3),
                       obs=round(Database$Q,3))
    rownames(Qout) <- format(Database$DatesR,'%Y-%m-%d')
    Ans <- list(Qout=Qout,
                Qsub=Qsub,
                SM=Prod,
                Dates=Dates,
                GR2M_ID=sub.id,
                EndState=EndState)
  }else{
    Ans <- list(Qsub=Qsub,
                SM=Prod,
                Dates=Dates,
                GR2M_ID=sub.id,
                EndState=EndState)
  }

  # Save results
  if(Save==TRUE){

    # Create output folder
    dir.create(file.path(getwd(),'Outputs'), recursive=T, showWarnings=F)

    # Save dataframes
    if(Update==TRUE){
      MnYr1     <- format(floor_date(Sys.Date()-months(2), "month"),'%b%y')
      MnYr2     <- format(floor_date(Sys.Date()-months(1), "month"),'%b%y')

      ProdName1 <- paste0('SM_GR2MSemiDistr_',MnYr1,'.txt')
      ProdName2 <- paste0('SM_GR2MSemiDistr_',MnYr2,'.txt')
      Prod_Old  <- read.table(file.path(getwd(),'Outputs',ProdName1),
                              header=TRUE, sep='\t')
      Prod_New  <- as.data.frame(rbind(as.matrix(Prod_Old),as.matrix(Prod)))
      colnames(Prod_New) <- sub.id
      rownames(Prod_New) <- c(as.Date(rownames(Prod_Old)), Dates)
      write.table(Prod_New, file=file.path(getwd(),'Outputs',ProdName2),sep='\t')
      file.remove(file.path(getwd(),'Outputs',ProdName1))
    } else{
      MnYr     <- format(as.Date(paste0('01/',RunEnd),"%d/%m/%Y"),"%b%y")
      ProdName <- paste0('SM_GR2MSemiDistr_',MnYr,'.txt')
      write.table(Prod, file=file.path(getwd(),'Outputs',ProdName),sep='\t')
    }
  }
  message('Done!')
  toc()
  return(Ans)

} # End (not run)
