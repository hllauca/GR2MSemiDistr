#' Run the GR2M model for multiple subbasins.
#'
#' This function runs the GR2M monthly water balance model in a semi-distributed manner,
#' applying region-specific parameters and correction factors to precipitation and
#' potential evapotranspiration data for each subbasin.
#'
#' @param Data A data frame containing the model inputs in airGR format. This must include a column DatesR (formatted as date or string), followed by columns P_1 to P_n for precipitation, E_1 to E_n for evapotranspiration, and optionally Q for observed discharge.
#' @param Subbasins A SpatVector object with the geometries of the subbasins. It must contain the attributes COMID (identifier) and Region (hydro-climatic region).
#' @param RunIni Initial date of the simulation period in "mm/yyyy" format.
#' @param RunEnd Final date of the simulation period in "mm/yyyy" format.
#' @param Parameters A data frame with the model parameters and correction factors per region. It must contain the columns: Region, X1, X2, fp, and fe.
#' @param WarmUp Optional. Number of initial months to discard from the outputs as warm-up. Default is NULL.
#' @param IniState Optional list of initial model states per subbasin (used to resume simulations). Default is NULL.
#' @param Save Logical. If TRUE, output files will be saved as .txt files in the ./Outputs folder. Default is FALSE.
#' @param Update Logical. If TRUE and Save = TRUE, removes the output files from the previous month before saving the current ones. Useful for operational updating. Default is FALSE.
#'
#' @return A list containing the following elements:
#' \describe{
#'   \item{PR}{Matrix of precipitation [mm/month] per subbasin.}
#'   \item{AE}{Matrix of actual evapotranspiration [mm/month] per subbasin.}
#'   \item{SM}{Matrix of soil moisture [mm/month] per subbasin.}
#'   \item{RU}{Matrix of runoff [mm/month] per subbasin.}
#'   \item{QS}{Matrix of simulated discharge [m3/s] per subbasin (not routed).}
#'   \item{Dates}{Vector of dates used in the simulation.}
#'   \item{COMID}{Vector of subbasin COMID identifiers.}
#'   \item{EndState}{List of final model states for each subbasin.}
#'   \item{SINK}{Optional. Data frame with simulated and observed streamflow at outlet, if Q is provided.}
#' }
#'
#' @references Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. (2021).
#' PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020). Water, 13(8), 1048.
#' \doi{10.3390/w13081048}
#'
#' @examples
#' model <- Run_GR2MSemiDistr(
#'   Data = data,
#'   Subbasins = roi,
#'   RunIni = '01/1981',
#'   RunEnd = '12/2016',
#'   Parameters = data.frame(
#'     Region = c("R1"),
#'     X1 = 10.5,
#'     X2 = 1.2,
#'     fp = 1.1,
#'     fe = 0.95
#'   )
#' )
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

  tictoc::tic()

  # === Validate Subbasins ===
  if (!inherits(Subbasins, "SpatVector")) {
    stop("Argument 'Subbasins' must be of class 'SpatVector'.")
  }
  if (!all(c("COMID", "Region") %in% names(Subbasins))) {
    stop("The 'Subbasins' object must contain fields 'COMID' and 'Region'.")
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

  # === Organize and validate parameters ===
  Zone  <- sort(unique(region))
  nreg  <- length(Zone)

  if (!is.data.frame(Parameters)) {
    stop("Argument 'Parameters' must be a data frame.")
  }
  required_cols <- c("Region", "X1", "X2", "fp", "fe")
  if (!all(required_cols %in% names(Parameters))) {
    stop("The 'Parameters' data frame must contain the following columns: 'Region', 'X1', 'X2', 'fp', 'fe'.")
  }

  if (nreg != nrow(Parameters)) {
    stop(paste0("The number of regions in 'Subbasins' (", nreg,
                ") does not match the number of rows in 'Parameters' (", nrow(Parameters), ")."))
  }

  # === Helper functions ===
  Subset_Param <- function(Param, Region) {
    c(Param$X1[Param$Region == Region],
      Param$X2[Param$Region == Region])
  }

  Forcing_Subbasin <- function(Param, Region, Database, ID, Nsub) {
    factor_p <- Param$fp[Param$Region == Region]
    factor_e <- Param$fe[Param$Region == Region]
    data.frame(DatesR = as.POSIXct(Database[, 1], tz = "GMT"),
               P = factor_p * Database[, 1 + ID],
               E = factor_e * Database[, 1 + ID + Nsub])
  }

  # === Run GR2M model ===
  message(sprintf("Running GR2M model for %d subbasins", nsub))

  ResModel <- vector("list", nsub)
  for (i in seq_len(nsub)) {
    Param <- Subset_Param(Parameters, region[i])
    Input <- Forcing_Subbasin(Parameters, region[i], Database, i, nsub)

    if (anyNA(Input$P) || anyNA(Input$E)) {
      stop(sprintf("Missing values (NA) detected in P or E for subbasin %s", comid[i]))
    }
    if (any(Input$P < 0)) {
      stop(sprintf("Negative values detected in P for subbasin %s", comid[i]))
    }
    if (any(Input$E < 0)) {
      stop(sprintf("Negative values detected in E for subbasin %s", comid[i]))
    }

    if (ntime == 1) {
      next_month <- lubridate::floor_date(Input$DatesR + months(1), "month")
      Input <- rbind(Input, data.frame(DatesR = next_month, P = 100, E = 100))
    }

    InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                     DatesR = Input$DatesR,
                                     Precip = Input$P,
                                     PotEvap = Input$E)

    RunOptions <- if (is.null(IniState)) {
      CreateRunOptions(FUN_MOD = RunModel_GR2M,
                       InputsModel = InputsModel,
                       IndPeriod_Run = seq_len(ntime),
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

      CreateRunOptions(FUN_MOD = RunModel_GR2M,
                       InputsModel = InputsModel,
                       IniStates = IniStates,
                       IndPeriod_Run = seq_len(ntime),
                       verbose = FALSE,
                       warnings = FALSE)
    }

    ResModel[[i]] <- RunModel(InputsModel = InputsModel,
                              RunOptions = RunOptions,
                              Param = Param,
                              FUN = RunModel_GR2M)
  }

  # === Post-processing ===
  extract_outputs <- function(res, area, nDays) {
    list(
      pr = round(res$Precip, 2),
      ae = round(res$AE, 2),
      sm = round(res$Prod, 2),
      ru = round(res$Qsim, 2),
      qs = round((area * res$Qsim) / (86.4 * nDays), 2)
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
    Ans$SINK <- data.frame(sim = round(qt, 2), obs = round(qo, 2), row.names = Dates)
  }

  if (Save) {
    dir.create("./Outputs", showWarnings = FALSE, recursive = TRUE)

    if (Update) {
      new_month <- format(tail(Dates, 1), "%Y%m")
      old_month <- format(as.Date(paste0(new_month, "01"), "%Y%m%d") %m-% months(1), "%Y%m")
      vars <- c("PR", "AE", "SM", "RU", "QS")
      for (var in vars) {
        old_file <- paste0("./Outputs/", var, "_GR2MSemiDistr_", old_month, ".txt")
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
    save_outputs(qs, "qs")
  }

  message("Processing completed successfully in")
  tictoc::toc()
  return(Ans)
}
