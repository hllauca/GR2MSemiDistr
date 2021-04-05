#' Uncertainty analysis of GR2M model parameters with the MCMC algorithm.
#'
#' @param Data        File with input data in airGR format (DatesR,P,E,Q).
#' @param Subbasins   Subbasins shapefile.
#' @param Dem         Raster DEM filename.
#' @param RunIni      Initial date of model simulation (in mm/yyyy format).
#' @param RunEnd      Final date of model simulation (in mm/yyyy format).
#' @param WarmUp      Number of months for warm-up. NULL as default.
#' @param Parameters      GR2M model parameters and correction factor of P and E.
#' @param Parameters.Min  Minimum values of GR2M model parameters and correction factor of P and E.
#' @param Parameters.Max  Maximum values of GR2M model parameters and correction factor of P and E.
#' @param Niter 	        Number of iterations. 1000 as default.
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @param MCMC    MCMC data in .Rda format.
#' @return  Lower(Q5) and upper (Q95) streamflows uncertainty bounds.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  FME
#' @import  parallel
#' @import  tictoc
#' @import  airGR
#' @import  abind
Uncertainty_GR2MSemiDistr <- function(Data,
                                      Subbasins,
                                      Dem,
                                      RunIni,
                                      RunEnd,
                                      WarmUp=NULL,
                                      Parameters,
                                      Parameters.Min,
                                      Parameters.Max,
                                      Niter,
                                      IniState=NULL,
                                      MCMC=NULL){

  # Load packages
  require(rgdal)
  require(raster)
  require(rgeos)
  require(parallel)
  require(tictoc)
  require(airGR)
  require(FME)
  require(abind)
  tic()

  # Generate parameters
  if(is.null(MCMC)==TRUE){

    # Show message
    cat('\f')
    message('Calculating parameter uncertainty with MCMC')
    message('Please wait...')

    # Residual function
    RFUN <- function(Variable){

      Ans <- Run_GR2MSemiDistr(Data=Data,
                               Subbasins=Subbasins,
                               RunIni=RunIni,
                               RunEnd=RunEnd,
                               Parameters=Variable,
                               IniState=IniState,
                               Regional=FALSE,
                               Update=FALSE,
                               Save=FALSE)

      # Calculate residuals
      if(is.null(WarmUp)==TRUE){
        Qobs <- Ans$SINK$obs
        Qsim <- Ans$SINK$sim
      }else{
        Qobs <- Ans$SINK$obs[-WarmUp:-1]
        Qsim <- Ans$SINK$sim[-WarmUp:-1]
      }
      mRes <- as.vector(na.omit(Qsim-Qobs))
      return(mRes)
    } # End function

    msr  <- mean((RFUN(Parameters))^2)
    MCMC <- modMCMC(f=RFUN,
                    p=Parameters,
                    lower=Parameters.Min,
                    upper=Parameters.Max,
                    niter=Niter,
                    var0=msr)
    dir.create('./Inputs',recursive=T,showWarnings=F)
    save(MCMC, file='./Inputs/MCMC.Rda')
  }else{
    load('Inputs/MCMC.Rda')
  }


  # Show message
  cat('\f')
  message('Generating streamflow uncertainty bounds')
  message('Please wait...')
  pars <- unique(MCMC$pars)
  Ans  <- list()
  for(w in 1:nrow(pars)){
    # Run model
    Qmod <- Run_GR2MSemiDistr(Data=Data,
                              Subbasins=Subbasins,
                              RunIni=RunIni,
                              RunEnd=RunEnd,
                              Parameters=pars[w,],
                              IniState=IniState,
                              Regional=FALSE,
                              Update=FALSE,
                              Save=FALSE)
    if(is.null(WarmUp)==TRUE){
      Ans[[w]] <- as.matrix(Qmod$QS)
      Dates    <- Qmod$Dates
    }else{
      Ans[[w]] <- as.matrix(Qmod$QS[-WarmUp:-1,])
      Dates    <- Qmod$Dates[-WarmUp:-1]
    }
  }
  qsub <- abind(Ans, along=3)
  q5   <- apply(qsub, c(1,2), function(x) quantile(x,0.05))
  q95  <- apply(qsub, c(1,2), function(x) quantile(x,0.95))

  # Routing streamflow for each subbasin at quantile 5
  M5 <- list(QS=q5,Dates=Dates)
  Q5 <- Routing_GR2MSemiDistr(Model=M5,
                              Subbasins=Subbasins,
                              Dem=Dem,
                              AcumIni=format(as.Date(head(Dates,1)),'%m/%Y'),
                              AcumEnd=format(as.Date(tail(Dates,1)),'%m/%Y'),
                              Save=FALSE,
                              Update=FALSE)

  M95 <- list(QS=q95,Dates=Dates)
  Q95 <- Routing_GR2MSemiDistr(Model=M95,
                               Subbasins=Subbasins,
                               Dem=Dem,
                               AcumIni=format(as.Date(head(Dates,1)),'%m/%Y'),
                               AcumEnd=format(as.Date(tail(Dates,1)),'%m/%Y'),
                               Save=FALSE,
                               Update=FALSE)

  # Prepare data to export
  comid <- paste0('COMID_',as.vector(Subbasins$COMID))
  QR5   <- as.data.frame(round(Q5,3))
  QR95  <- as.data.frame(round(Q95,3))
  colnames(QR5)  <- comid
  colnames(QR95) <- comid
  rownames(QR5)  <- as.Date(Dates)
  rownames(QR95) <- as.Date(Dates)
  Ans <- list(lower=QR5,
              upper=QR95)

  # Show message
  message("Done!")
  toc()
  return(Ans)

} # End (not run)
