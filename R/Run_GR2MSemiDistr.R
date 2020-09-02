#' Run the GR2M model for each subbasins.
#'
#' @param Data        File with input data in airGR format (DatesR,P,E,Q)
#' @param Subbasins   Subbasins shapefile.
#' @param RunIni      Initial date of model simulation (in mm/yyyy format).
#' @param RunEnd      Final date of model simulation (in mm/yyyy format).
#' @param Parameters  GR2M model parameters and correction factor of P and E.
#' @param Plot        Logical value to plot observed and simulated streamflow timeseries. FALSE as default.
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @param Regional    Logical value to simulate in a regional mode (more than one outlet). FALSE as default.
#' @param Update      Logical value to update a previous production and qsubbasin '.csv' files. FALSE as default.
#' @param Save        Logical valute to export simulation results as '.Rda'. TRUE as default.
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
                              Subbasin,
                              RunIni,
                              RunEnd,
                              Parameters,
                              Plot=TRUE,
                              IniState=NULL,
                              Regional=FALSE,
                              Update=FALSE,
                              Save=TRUE){

  # Data=Data
  # Subbasin=roi
  # RunIni=RunModel.Ini
  # RunEnd=RunModel.End
  # Parameters=Model.Param
  # Plot=TRUE
  # IniState=NULL
  # Regional=FALSE
  # Update=FALSE
  # Save=TRUE

  # Load packages
  require(rgdal)
  require(raster)
  require(rgeos)
  require(rtop)
  require(parallel)
  require(lubridate)
  require(airGR)
  require(hydroGOF)
  require(tictoc)
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

  # Show message
  cat('\f')
  message(paste('Running GR2M model for', nsub, 'subbasins'))
  message('Please wait...')

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

  # Close cluster
  stopCluster(cl)

  # Model results
  if(nsub==1){
    prod <- ResModel[[1]]$Prod
    qSub <- (area[1]*ResModel[[1]]$Qsim)/(86.4*nDays)
    qOut <- qSub
    EndState <- list(ResModel[[1]]$StateEnd)
  }else{
    Plist <- list()
    Qlist <- list()
    for(w in 1:nsub){
      Plist[[w]] <- ResModel[[w]]$Prod
      Qlist[[w]] <- (area[w]*ResModel[[w]]$Qsim)/(86.4*nDays)
    }
    prod <- do.call(cbind, Plist)
    qSub <- do.call(cbind, Qlist)
    qOut <- round(apply(qSub, 1, FUN=sum),2)
    EndState <- list()
    for(w in 1:nsub){
      EndState[[w]] <- ResModel[[w]]$StateEnd
    }
  }

  # Subset model results (exclude warm-up)
  Qsub <- qSub
  Prod <- prod

  # Factors Fpp and Fpet
  pp  <- matrix(NA, ncol=nsub, nrow=nrow(Qsub))
  pet <- matrix(NA, ncol=nsub, nrow=nrow(Qsub))
  for (w in 1:nsub){
    pp[,w]  <- subset(Param$Fpp, Param$Region==region[w])*Database[,(w+1)]
    pet[,w] <- subset(Param$Fpet, Param$Region==region[w])*Database[,(nsub+w)]
  }

  # Local mode
  if(Regional==FALSE){
    Qobs <- Database$Q
    Qsim <- qOut

    if(Plot==TRUE){
      x11()
      ggof(Qsim, Qobs, main='Streamflow at basin outlet', digits=3, gofs=c("NSE","KGE","RMSE"))
    }
    Ans <- list(Qsim=Qsim,
                Qobs=Qobs,
                Qsub=Qsub,
                Precip=pp,
                Evaptr=pet,
                Prod=Prod,
                Dates=format(Database$DatesR,'%Y-%m-%d'),
                EndState=EndState)
  }else{
    Ans <- list(Qsub=Qsub,
                Precip=pp,
                Evaptr=pet,
                Prod=Prod,
                Dates=format(Database$DatesR,'%Y-%m-%d'),
                EndState=EndState)
  }

  # Save results
  if(Save==TRUE){

    # Create output folder
    dir.create(file.path(getwd(),'Outputs'), recursive=T, showWarnings=F)

    # Save simulation
    save(Ans, file='./Outputs/Simulation_GR2MSemiDistr.Rda')

    # Save dataframes
    if(Update==TRUE){
      MnYr1     <- format(floor_date(Sys.Date()-months(2), "month"),'%b%y')
      MnYr2     <- format(floor_date(Sys.Date()-months(1), "month"),'%b%y')

      ProdName1 <- paste0('Production_GR2MSemiDistr_',MnYr1,'.txt')
      ProdName2 <- paste0('Production_GR2MSemiDistr_',MnYr2,'.txt')
      Data      <- read.table(file.path(getwd(),'Outputs',ProdName1), header=T, sep=',')
      Dates     <- as.Date(Data$Dates, tryFormats=c('%Y-%m-%d','%Y/%m/%d','%d/%m/%Y','%d-%m-%Y'))
      Prod_Old  <- Data[,-1]
      Prod_New  <- rbind(as.matrix(Prod_Old),Prod)
      Dates_New <- c(Dates,as.Date(Database$DatesR))
      DataProd  <- data.frame(Dates_New, Prod_New)
      colnames(DataProd) <- c('Dates', paste0('GR2M-ID_',1:nsub))
      write.table(DataProd, file=file.path(getwd,'Outputs',ProdName2), sep='\t', row.names=FALSE)
      file.remove(file.path(getwd(),'Outputs',ProdName1))

      QsubName1 <- paste0('Qsubbasins_GR2MSemiDistr_',MnYr1,'.txt')
      QsubName2 <- paste0('Qsubbasins_GR2MSemiDistr_',MnYr2,'.txt')
      Data      <- read.table(file.path(getwd(),'Outputs',QsubName1), header=T, sep=',')
      Dates     <- as.Date(Data$Dates, tryFormats=c('%Y-%m-%d','%Y/%m/%d','%d/%m/%Y','%d-%m-%Y'))
      Qsub_Old  <- Data[,-1]
      Qsub_New  <- rbind(as.matrix(Qsub_Old), Qsub)
      Dates_New <- c(Dates,as.Date(Database$DatesR))
      DataQsub  <- data.frame(Dates_New, Qsub_New)
      colnames(DataQsub) <- c('Dates', paste0('GR2M-ID_',1:nsub))
      write.table(DataQsub, file=file.path(getwd(),'Outputs',QsubName2), sep='\t', row.names=FALSE)
      file.remove(file.path(getwd(),'Outputs',QsubName1))

    } else{
      MnYr     <- format(as.Date(paste0('01/',RunEnd),"%d/%m/%Y"),"%b%y")
      DataProd <- data.frame(format(Database$DatesR,'%Y-%m-%d'), Prod)
      DataQsub <- data.frame(format(Database$DatesR,'%Y-%m-%d'), Qsub)
      ProdName <- paste0('Production_GR2MSemiDistr_',MnYr,'.txt')
      QsubName <- paste0('Qsubbasins_GR2MSemiDistr_',MnYr,'.txt')
      colnames(DataProd) <- c('Dates', paste0('ID_',1:nsub))
      colnames(DataQsub) <- c('Dates', paste0('ID_',1:nsub))
      write.table(DataProd, file=file.path(getwd(),'Outputs',ProdName), sep='\t', row.names=FALSE)
      write.table(DataQsub, file=file.path(getwd(),'Outputs',QsubName), sep='\t', row.names=FALSE)
    }
  }
  message('Done!')
  toc()
  return(Ans)

} # End (not run)
