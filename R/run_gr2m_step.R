#' Run lumped GR2M model for each sub-basin and monthly timestep.
#'
#' @param Model.Input Model forcing data in airGR format '(DatesR,P,T)'.
#' @param Model.Parameters GR2M (X1 and X2) model parameters.
#' @param Model.State GR2M state variables for the computing timestep.
#' @param Model.Run Date 'mm/yyyy' to compute timestep.
#' @return GR2M model output for an specific timestep.
#' @export0
run_gr2m_step <- function(Model.Input, Model.Parameter, Model.State, Model.Run){

    # Load packages
    require(airGR)

    # Prepare model inputs
    InputsModel   <- CreateInputsModel(FUN_MOD=RunModel_GR2M,
                                       DatesR=Model.Input$DatesR,
                                       Precip=Model.Input$P,
                                       PotEvap=Model.Input$E)

    # Extract data for an specific timestep
    Ind_Run       <- which(format(Model.Input$DatesR, format="%m/%Y") == Model.Run)

    # Run GR2M model by an specific initial conditions
    if(is.null(Model.State)==TRUE){

      # Set-up running options
      RunOptions    <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                        InputsModel=InputsModel,
                                        IndPeriod_Run=Ind_Run)

      # Run GR2M
      OutputsModel  <- RunModel(InputsModel=InputsModel,
                                RunOptions=RunOptions,
                                Param=Model.Parameter,
                                FUN=RunModel_GR2M)
    } else{

      # Set-up running options
      RunOptions    <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                        InputsModel=InputsModel,
                                        IniStates=Model.State,
                                        IndPeriod_Run=Ind_Run,
                                        verbose=FALSE,
                                        warnings=FALSE)

      # Run GR2M
      OutputsModel  <- RunModel(InputsModel=InputsModel,
                                RunOptions=RunOptions,
                                Param=Model.Parameter,
                                FUN=RunModel_GR2M)
    }

    # Results
    return(OutputsModel)
}
