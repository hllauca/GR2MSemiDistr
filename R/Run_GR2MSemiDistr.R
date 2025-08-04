#' Run the GR2M model for multiple subbasins.
#'
#' This function runs the GR2M monthly water balance model on a semi-distributed basis
#' across multiple subbasins, applying region-specific parameters and correction factors
#' to precipitation and potential evapotranspiration.
#'
#' @param Data Data frame of model inputs in airGR format created with \code{Create_Forcing_Inputs}.
#' Columns should include DatesR, followed by P_1 to P_n, E_1 to E_n, and optionally Q.
#' @param Subbasins A `SpatVector` object containing the subbasins' geometries. Must include the following attributes: `'COMID'` (unique identifier), and `'Region'` (string code for the hydro-climatic region).
#' @param RunIni Initial simulation date in 'mm/yyyy' format.
#' @param RunEnd Final simulation date in 'mm/yyyy' format.
#' @param WarmUp Optional number of months to discard as warm-up. Default is NULL.
#' @param Parameters Numeric vector of model parameters and correction factors.
#' Format: c(X1_region1, ..., X1_regionN, X2_region1, ..., Fpp_regionN, Fpe_region1, ..., Fpe_regionN)
#' @param IniState Optional list of initial state values for each subbasin. Default is NULL.
#' @param Save Logical. If TRUE, outputs will be saved as text files in the 'Outputs/' folder.
#' @param Update Logical. If TRUE and Save = TRUE, overwrites previous month's files before saving new ones.
#'
#' @return A list with the following elements:
#' \item{PR}{Matrix of precipitation [mm/month] per subbasin}
#' \item{AE}{Matrix of actual evapotranspiration [mm/month] per subbasin}
#' \item{SM}{Matrix of soil moisture [mm/month] per subbasin}
#' \item{RU}{Matrix of runoff [mm/month] per subbasin}
#' \item{QS}{Matrix of simulated discharge [m3/s] per subbasin (not routed)}
#' \item{Dates}{Vector of dates used in the simulation}
#' \item{COMID}{Vector of subbasin COMID identifiers}
#' \item{EndState}{List of final model states for each subbasin}
#' \item{SINK}{Optional. Data frame of simulated and observed streamflow at outlet, if Q is provided}
#'
#' @references Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. (2021).
#' PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020). Water, 13(8), 1048.
#' \doi{10.3390/w13081048}
#'
#' @export
#'
#' @examples
#' # Load input forcing data and subbasins (see Create_Forcing_Inputs)
#' model <- Run_GR2MSemiDistr(
#'   Data = data,
#'   Subbasins = roi,
#'   RunIni = '01/1981',
#'   RunEnd = '12/2016',
#'   Parameters = c(10.976, 0.665, 1.186, 1.169)
#' )
#'
#' # View model outputs
#' View(model$PR)  # Precipitation [mm/month]
#' View(model$QS)  # Simulated discharge [m3/s] per subbasin (not routed).
#' print(model$SINK$SIM)  # Simulated outlet discharge
#' print(model$SINK$OBS)  # Observed outlet discharge (if available)
#'
#' @import airGR
#' @import tictoc
#' @import lubridate
#' @importFrom terra expanse

Run_GR2MSemiDistr <- function(Data,
                              Subbasins,
                              RunIni,
                              RunEnd,
                              Parameters,
                              WarmUp = NULL,
                              IniState = NULL,
                              Save = FALSE,
                              Update = FALSE) {

  tic()

  # === Validate Subbasins ===
  if (!inherits(Subbasins, "SpatVector")) {
    stop("Argument 'Subbasins' must be of class 'SpatVector'.")
  }
  if (!"COMID" %in% names(Subbasins)) {
    stop("The 'Subbasins' object must contain a field named 'COMID'.")
  }
  if (!"Region" %in% names(Subbasins)) {
    stop("The 'Subbasins' object must contain a field named 'Region'.")
  }
  comid  <- as.vector(Subbasins$COMID)
  region <- Subbasins$Region
  area   <- terra::expanse(Subbasins, unit = 'km')
  nsub   <- length(comid)

  # === Prepare time indices ===
  Data$DatesR <- as.POSIXct(paste0(Data$DatesR, ' 00:00:00'),
                            tz = "GMT",
                            tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%d-%m-%Y", "%d/%m/%Y"))

  ind_start <- which(format(Data$DatesR, "%m/%Y") == RunIni)
  ind_end   <- which(format(Data$DatesR, "%m/%Y") == RunEnd)

  if (length(ind_start) == 0 || length(ind_end) == 0) {
    stop("RunIni or RunEnd not found in Data$DatesR. Please check date format 'mm/yyyy'.")
  }

  Ind_run  <- ind_start:ind_end
  Database <- Data[Ind_run, ]
  Dates    <- as.Date(Database$DatesR)
  ntime    <- length(Dates)
  nDays    <- lubridate::days_in_month(Dates)

  # === Organize parameters by region ===
  Zone  <- sort(unique(region))
  nreg  <- length(Zone)
  Param <- data.frame(Region = Zone,
                      X1 = Parameters[1:nreg],
                      X2 = Parameters[(nreg + 1):(2 * nreg)],
                      Fpp = Parameters[(2 * nreg + 1):(3 * nreg)],
                      Fpe = Parameters[(3 * nreg + 1):length(Parameters)])

  # === Helper functions ===
  Subset_Param <- function(Param, Region) {
    c(Param$X1[Param$Region == Region],
      Param$X2[Param$Region == Region])
  }

  Forcing_Subbasin <- function(Param, Region, Database, Nsub, ID) {
    factor_pp  <- Param$Fpp[Param$Region == Region]
    factor_pet <- Param$Fpe[Param$Region == Region]
    Inputs     <- Database[, c(1, ID + 1, ID + 1 + Nsub)]
    data.frame(
      DatesR = as.POSIXct(Inputs[, 1], tz = "GMT"),
      P = round(factor_pp * Inputs[, 2], 1),
      E = round(factor_pet * Inputs[, 3], 1)
    )
  }

  # === Run GR2M model ===
  message(sprintf("Running GR2M model for %d subbasins", nsub))

  ResModel <- vector("list", nsub)
  for (i in seq_len(nsub)) {
    ParamSub  <- Subset_Param(Param, region[i])
    FixInputs <- Forcing_Subbasin(Param, region[i], Database, nsub, i)

    if (ntime == 1) {
      next_month <- as.POSIXct(lubridate::floor_date(FixInputs$DatesR + months(1), "month"))
      FixInputs  <- rbind(FixInputs, data.frame(DatesR = next_month, P = 100, E = 100))
    }

    InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                     DatesR = FixInputs$DatesR,
                                     Precip = FixInputs$P,
                                     PotEvap = FixInputs$E)

    if (is.null(IniState)) {
      RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                     InputsModel = InputsModel,
                                     IndPeriod_Run = 1:ntime,
                                     verbose = FALSE,
                                     warnings = FALSE)
    } else {
      IniStates <- CreateIniStates(FUN_MOD = RunModel_GR2M,
                                   InputsModel = InputsModel,
                                   ProdStore = IniState[[i]]$Store$Prod,
                                   RoutStore = IniState[[i]]$Store$Rout,
                                   ExpStore = IniState[[i]]$Store$Exp,
                                   UH1 = IniState[[i]]$UH$UH1,
                                   UH2 = IniState[[i]]$UH$UH2)

      RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                     InputsModel = InputsModel,
                                     IniStates = IniStates,
                                     IndPeriod_Run = 1:ntime,
                                     verbose = FALSE,
                                     warnings = FALSE)
    }

    ResModel[[i]] <- RunModel(InputsModel = InputsModel,
                              RunOptions = RunOptions,
                              Param = ParamSub,
                              FUN = RunModel_GR2M)
  }

  # === Post-processing ===
  extract_outputs <- function(res, area, nDays) {
    list(
      pr = round(res$Precip, 1),
      ae = round(res$AE, 1),
      sm = round(res$Prod, 1),
      ru = round(res$Qsim, 1),
      qs = round((area * res$Qsim) / (86.4 * nDays), 1)
    )
  }

  outputs <- lapply(seq_len(nsub), function(i) extract_outputs(ResModel[[i]], area[i], nDays))
  pr <- do.call(cbind, lapply(outputs, `[[`, "pr"))
  ae <- do.call(cbind, lapply(outputs, `[[`, "ae"))
  sm <- do.call(cbind, lapply(outputs, `[[`, "sm"))
  ru <- do.call(cbind, lapply(outputs, `[[`, "ru"))
  qs <- do.call(cbind, lapply(outputs, `[[`, "qs"))
  qt <- round(rowSums(qs), 2)

  if (!is.null(WarmUp)) {
    indices <- -(seq_len(WarmUp))
    pr <- pr[indices, , drop = FALSE]
    ae <- ae[indices, , drop = FALSE]
    sm <- sm[indices, , drop = FALSE]
    ru <- ru[indices, , drop = FALSE]
    qs <- qs[indices, , drop = FALSE]
    qt <- qt[indices]
    Dates <- Dates[indices]
    if ("Q" %in% names(Database)) {
      qo <- Database$Q[indices]
    }
  } else if ("Q" %in% names(Database)) {
    qo <- Database$Q
  }

  Ans <- list(QS = qs, RU = ru, PR = pr, AE = ae, SM = sm, Dates = Dates, COMID = comid,
              EndState = lapply(ResModel, `[[`, "StateEnd"))

  if (exists("qo")) {
    Ans$SINK <- data.frame(sim = round(qt, 3), obs = round(qo, 3), row.names = Dates)
  }

  if (Save) {
    dir.create("./Outputs", showWarnings = FALSE, recursive = TRUE)

    if (Update) {
      last_month <- format(tail(Dates, 1), "%Y%m")
      previous_month <- format(as.Date(paste0(last_month, "01"), "%Y%m%d") %m-% months(1), "%Y%m")
      vars <- c("PR", "AE", "SM", "RU", "QS")
      for (var in vars) {
        old_file <- paste0("./Outputs/", var, "_GR2MSemiDistr_", previous_month, ".txt")
        if (file.exists(old_file)) file.remove(old_file)
      }
    }

    save_outputs <- function(var, name) {
      df <- as.data.frame(var)
      colnames(df) <- paste0(toupper(name), "_", comid)
      rownames(df) <- Dates
      file <- paste0("./Outputs/", toupper(name), "_GR2MSemiDistr_", format(tail(Dates, 1), "%Y%m"), ".txt")
      write.table(df, file = file)
    }
    save_outputs(pr, "pr")
    save_outputs(ae, "ae")
    save_outputs(sm, "sm")
    save_outputs(ru, "ru")
    save_outputs(ru, "qs")
  }

  message("Done!")
  tictoc::toc()
  return(Ans)
}
