#' Run the GR2M model for multiple subbasins
#'
#' This function executes the GR2M monthly water balance model across multiple subbasins
#' using a semi-distributed approach. It applies region-specific parameters and correction
#' factors for precipitation and evapotranspiration, and optionally saves the results as text files.
#'
#' @param Data A data frame of model input data in the airGR format, as produced by Create_Forcing_Inputs.
#' It must include columns: DatesR, P_1 to P_n, E_1 to E_n, and optionally Q.
#' @param Subbasins A SpatVector object containing the geometries of subbasins. It must include attributes "COMID" (unique subbasin ID) and "Region" (region name/code).
#' @param RunIni Simulation start date in the format "mm/yyyy".
#' @param RunEnd Simulation end date in the format "mm/yyyy".
#' @param Parameters A data frame containing GR2M model parameters and correction factors per region. Must have columns: Region, X1, X2, fp, fe.
#' @param WarmUp Optional number of months to discard from the beginning of the simulation as warm-up. Default is NULL.
#' @param IniState Optional list of initial states for each subbasin. Default is NULL.
#' @param Save Logical. If TRUE, output time series will be saved as text files in the "Outputs/" folder.
#' @param Update Logical. If TRUE and Save = TRUE, it will remove previous month's output files before saving the current ones.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{PR}{Matrix of total precipitation [mm/month] for each subbasin}
#'   \item{AE}{Matrix of actual evapotranspiration [mm/month] for each subbasin}
#'   \item{SM}{Matrix of production store level [mm/month] for each subbasin}
#'   \item{SM}{Matrix of percolation [mm/month] for each subbasin}
#'   \item{RU}{Matrix of runoff [mm/month] for each subbasin}
#'   \item{QS}{Matrix of flow [m3/s] (not routed) for each subbasin}
#'   \item{Dates}{Date vector for the simulation period}
#'   \item{COMID}{Vector of subbasin COMIDs}
#'   \item{EndState}{List of final state variables for each subbasin}
#'   \item{SINK}{Optional. Data frame of simulated and observed outlet discharge, if observed Q is provided}
#' }
#'
#' @references Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. (2021).
#' PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020). Water, 13(8), 1048. \doi{10.3390/w13081048}
#'
#' @examples
#' # Run the GR2M model
#' model <- Run_GR2MSemiDistr(
#'   Data = data,
#'   Subbasins = roi,
#'   RunIni = "01/1981",
#'   RunEnd = "12/2016",
#'   Parameters = data.frame(
#'     Region = c("A"),
#'     X1 = 10.5,
#'     X2 = 0.8,
#'     fp = 1.1,
#'     fe = 1.0
#'   )
#' )
#'
#' @import airGR
#' @import tictoc
#' @import lubridate
#' @importFrom terra expanse
#'
#' @export
#'
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

  Database <- Data[ind_start:ind_end, ]
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
      pc = round(res$Perc, 2),
      ru = round(res$Qsim, 2),
      qs = round((area * res$Qsim) / (86.4 * nDays), 2)
    )
  }

  outputs <- lapply(seq_len(nsub), function(i) extract_outputs(ResModel[[i]], area[i], nDays))
  pr <- do.call(cbind, lapply(outputs, `[[`, "pr"))
  ae <- do.call(cbind, lapply(outputs, `[[`, "ae"))
  sm <- do.call(cbind, lapply(outputs, `[[`, "sm"))
  pc <- do.call(cbind, lapply(outputs, `[[`, "pc"))
  ru <- do.call(cbind, lapply(outputs, `[[`, "ru"))
  qs <- do.call(cbind, lapply(outputs, `[[`, "qs"))
  qt <- round(rowSums(qs), 2)

  if (!is.null(WarmUp)) {
    indices <- -(seq_len(WarmUp))
    pr <- pr[indices, , drop = FALSE]
    ae <- ae[indices, , drop = FALSE]
    sm <- sm[indices, , drop = FALSE]
    pc <- pc[indices, , drop = FALSE]
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

  Ans <- list(Dates = Dates, COMID = comid, QS = qs, RU = ru, PR = pr, AE = ae, SM = sm, PC = pc,
              EndState = lapply(ResModel, `[[`, "StateEnd"))

  if (exists("qo")) {
    Ans$SINK <- data.frame(sim = round(qt, 2), obs = round(qo, 2), row.names = Dates)
  }

  if (Save) {
    dir.create("./Outputs", showWarnings = FALSE, recursive = TRUE)

    if (Update) {
      new_month <- format(tail(Dates, 1), "%Y%m")
      old_month <- format(as.Date(paste0(new_month, "01"), "%Y%m%d") %m-% months(1), "%Y%m")
      vars <- c("PR", "AE", "SM", "PC","RU", "QS")
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
