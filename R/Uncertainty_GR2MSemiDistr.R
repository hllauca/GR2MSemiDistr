#' Uncertainty analysis of GR2M model parameters with the MCMC algorithm.
#'
#' @param Data        File with input data in airGR format (DatesR,P,E,Q)
#' @param Subbasins   Subbasins shapefile.
#' @param Dem          Raster DEM.
#' @param RunIni      Initial date of model simulation (in mm/yyyy format).
#' @param RunEnd      Final date of model simulation (in mm/yyyy format).
#' @param WarmUp      Number of months for warm-up.
#' @param Parameters      GR2M model parameters and correction factor of P and E.
#' @param Parameters.Min  Minimum values of GR2M model parameters and correction factor of P and E.
#' @param Parameters.Max  Maximum values of GR2M model parameters and correction factor of P and E.
#' @param Niter 	        Number of iterations. 1000 as default.
#' @param IniState    Initial GR2M states variables. NULL as default.
#' @param Positions    Cell numbers to extract data faster for each subbasin. NULL as default.
#' @param MCMC    MCMC data in .Rda format.
#' @return  Lower(Q5) and upper (Q95) streamflows uncertanty bounds.
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
                                      Subbasin,
                                      Dem,
                                      RunIni,
                                      RunEnd,
                                      WarmUp,
                                      Parameters,
                                      Parameters.Min,
                                      Parameters.Max,
                                      Niter,
                                      IniState=NULL,
                                      Positions=NULL,
                                      MCMC=NULL){
  # Data=Ans1
  # Subbasin=roi
  # Dem='Subbasins.tif'
  # RunIni='01/1981'
  # RunEnd='05/1981'
  # WarmUp=2
  # Parameters=BestParam
  # Parameters.Min=c(1, 0.01, 0.8, 0.8)
  # Parameters.Max=c(2000, 2, 1.2, 1.2)
  # Niter=5
  # IniState=NULL
  # Positions=NULL
  # MCMC=NULL


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
                               Subbasin=Subbasin,
                               RunIni=RunIni,
                               RunEnd=RunEnd,
                               Parameters=Variable,
                               IniState=IniState,
                               Regional=FALSE,
                               Update=FALSE,
                               Save=FALSE)

      # Calculate residuals
      Qobs <- Ans$Qobs[-WarmUp:-1]
      Qsim <- Ans$Qsim[-WarmUp:-1]
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
    dir.create(file.path(getwd(),'Inputs'),recursive=T,showWarnings=F)
    save(MCMC, file=file.path(getwd(),'Inputs','MCMC.Rda'))
  }else{
    load(file.path(getwd(),'Inputs','MCMC.Rda'))
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
                              Subbasin=Subbasin,
                              RunIni=RunIni,
                              RunEnd=RunEnd,
                              Parameters=pars[w,],
                              IniState=IniState,
                              Regional=FALSE,
                              Update=FALSE,
                              Save=FALSE)
    Ans[[w]] <- Qmod$Qsub[-WarmUp:-1,]
    Dates    <- Qmod$Dates[-WarmUp:-1]
  }
  qsub <- abind(Ans, along=3)
  q5   <- apply(qsub, c(1,2), function(x) quantile(x,0.05))
  q95  <- apply(qsub, c(1,2), function(x) quantile(x,0.95))

  # Routing streamflow for each subbasin at quantile 5
  M5 <- list(Qsub=q5,Dates=Dates)
  Q5 <- Routing_GR2MSemiDistr(Model=M55,
                              Subbasin=Subbasin,
                              Dem=Dem,
                              AcumIni=format(as.Date(head(Dates,1)),'%m/%Y'),
                              AcumEnd=format(as.Date(tail(Dates,1)),'%m/%Y'),
                              Positions=Positions,
                              Save=FALSE,
                              Update=FALSE,
                              all=FALSE)

  M95 <- list(Qsub=q5,Dates=Dates)
  Q95 <- Routing_GR2MSemiDistr(Model=M95,
                               Subbasin=Subbasin,
                               Dem=Dem,
                               AcumIni=format(as.Date(head(Dates,1)),'%m/%Y'),
                               AcumEnd=format(as.Date(tail(Dates,1)),'%m/%Y'),
                               Positions=Positions,
                               Save=FALSE,
                               Update=FALSE,
                               all=FALSE)

  sens <- list(lower=Q5, upper=Q95)

  # Show message
  message("Done!")
  toc()
  return(sens)

} # End (not run)
