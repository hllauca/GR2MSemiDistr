#' Optimization of GR2M model parameters with MOPSOCD algorithm.
#'
#' @param Parameters		   GR2M (X1 and X2) model parameters and a multiplying factor to adjust monthly P and PET values.
#' @param Parameters.Min	 Minimum GR2M (X1, X2 and f) model parameters values.
#' @param Parameters.Max	 Maximum GR2M (X1, X2 and f) model parameters values.
#' @param Optimization		 Multi-objective evaluation criteria (NSE, lnNSE, KGE, RMSE, R, PBIAS)
#' @param Location		 General work directory where data is located.
#' @param Raster			 Flow direction raster in GRASS format.
#' @param Shapefile		 Subbasins shapefile.
#' @param Input				 Model forcing data in airGR format (DatesR,P,T,Qmm). 'Inputs_Basins.txt' as default.
#' @param WarmIni   	 Initial date 'mm/yyyy' of the warm-up period.
#' @param RunIni    	 Initial date 'mm/yyyy' of the model evaluation period.
#' @param RunEnd    	 Final date 'mm/yyyy' of the model evaluation period.
#' @param IdBasin   	 Subbasin ID number to compute outlet model (from shapefile attribute table).
#' @param Remove    	 Logical value to remove streamflow generated in the IdBasin. FALSE as default.
#' @param No.Optim     Calibration regions to exclude
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @return  Best semidistribute GR2M model parameters.
#' @export
#' @import  ProgGUIinR
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  mopsocd
#' @import  hydroGOF
#' @import  foreach
#' @import  tictoc
Optim2_GR2MSemiDistr <- function(Parameters, Parameters.Min, Parameters.Max, Optimization=c('NSE','KGE'),
								                 Location, Shapefile, Input='Inputs_Basins.txt', WarmIni,
								                 RunIni, RunEnd, IdBasin, Remove=FALSE, No.Optim=NULL, IniState=NULL){

# Parameters=Model.Param
# Parameters.Min=Model.ParMin
# Parameters.Max=Model.ParMax
# Optimization=Optimization=c('NSE','R')
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
# Variable=Model.Param

        # Load packages
        require(rgdal)
        require(raster)
        require(rgeos)
        require(ProgGUIinR)
        require(rtop)
        require(hydroGOF)
        require(foreach)
        require(tictoc)
        tic()

        # Load shapefile
        path.shp   <- file.path(Location,'Inputs', Shapefile)
        basin      <- readOGR(path.shp, verbose=F)
        area       <- basin@data$Area
        region     <- basin@data$Region
        nsub       <- nrow(basin@data)

        # Filtering calibration region no to optimize
        idx <- 1:length(Parameters)
        idy <- which(No.Optim==rep(region,2))
        Stb <- Parameters[idy]

        # Read and subset input data for the study period
        Data        <- read.table(file.path(Location, 'Inputs', Input), sep='\t', header=T)
        Data$DatesR <- as.POSIXct(Data$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))
        if(is.null(WarmIni)==TRUE){
          Subset      <- seq(which(format(Data$DatesR, format="%m/%Y") == RunIni),
                             which(format(Data$DatesR, format="%m/%Y") == RunEnd))
        } else{
          Subset      <- seq(which(format(Data$DatesR, format="%m/%Y") == WarmIni),
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
        if (is.null(No.Optim)==TRUE){
          Zone       <- sort(unique(region))
        } else{
          Parameters <- Parameters[!(rep(region,2) %in% No.Optim)]
          Zone       <- sort(unique(region[!(region %in% No.Optim)]))
        }
        Parameters.Min <- rep(Parameters.Min, each=length(Zone))
        Parameters.Max <- rep(Parameters.Max, each=length(Zone))

        # Objective function
        OFUN <- function(Variable){

          # Select model parameters to optimize
          if (is.null(No.Optim)==TRUE){
            Par.Optim <- Variable
          } else{
            dta       <- rbind(cbind(idx2[-idy], Variable), cbind(idy, Stb))
            Par.Optim <- dta[match(sort(dta[,1]), dta[,1]), 2]
          }

          # Model parameters to run GR2M model
          nreg  <- length(sort(unique(region)))
          Param <- data.frame(Region=sort(unique(region)),
                              X1=Par.Optim[1:nreg],
                              X2=Par.Optim[(nreg+1):(2*nreg)],
                              Fpp=Par.Optim[(2*nreg+1):(3*nreg)],
                              Fpet=Par.Optim[(3*nreg+1):length(Par.Optim)])

          # Run GR2M for each subbasin
          cl=makeCluster(detectCores()-1) # Detect and assign a cluster number
          clusterEvalQ(cl,c(library(airGR))) # Load package to each node
          clusterExport(cl,varlist=c("Param","region","nsub","Database","time", "IniState"),envir=environment())

          ResModel <- parLapply(cl, 1:nsub, function(i) {

                      # Parameters and factors to run the model
                      ParamSub  <- c(subset(Param$X1, Param$Region==region[i]), subset(Param$X2, Param$Region==region[i]))
                      FactorPP  <- subset(Param$Fpp, Param$Region==region[i])
                      FactorPET <- subset(Param$Fpet, Param$Region==region[i])
                      Inputs    <- Database[,c(1,i+1,i+1+nsub)]
                      FixInputs <- data.frame(DatesR=Inputs[,1], P=round(FactorPP*Inputs[,2],1), E=round(FactorPET*Inputs[,3],1))
                      FixInputs$DatesR <- as.POSIXct(FixInputs$DatesR, "GMT", tryFormats=c("%Y-%m-%d", "%d/%m/%Y"))

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

          # Evaluation criteria dataframe
          optim.df <- data.frame(KGE=round(KGE(Qsim, Qobs), 3),
                                 NSE=round(NSE(Qsim, Qobs), 3),
                                 lnNSE=round(NSE(log(Qsim), log(Qobs)), 3),
                                 R=round(rPearson(Qsim, Qobs), 3),
                                 RMSE=1-round(rmse(Qsim, Qobs), 3),
                                 PBIAS=1-round(pbias(Qsim, Qobs), 3))

          # Return
          MOF <- as.numeric(optim.df[colnames(optim.df) %in% Optimization])
          return(MOF)

    } # End objective function

  # Optimization with MOPSOCD
  # Show message
    cat('\f')
    message(paste('Optimizing', Optimization, 'with MOPSOCD'))
    message('Please wait...')
    Calibration <- mopsocd(OFUN, varcnt=length(Parameters), fncnt=length(Optimization),
                   lowerbound=Parameters.Min, upperbound=Parameters.Max, opt=1) # Maximizing
    fn   <- data.frame(Calibration$objfnvalues)
    pr   <- data.frame(Calibration$paramvalues)
    colnames(fn) <- Optimization
    fn2  <- fn[colnames(fn) %in% 'NSE']
    mof  <- subset(fn, fn2==max(fn2))
    par  <- pr[which(fn2==max(fn2)),]

  # Show message
    message("Done!")
    toc()

  # Output
    Ans <- list(Param=par, Value=mof, Pareto.Values=fn, Pareto.Param=pr)
    return(Ans)

} # End (not run)
