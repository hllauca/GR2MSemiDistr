#' Run the GR2M model for multiple subbasins
#'
#' This function executes the GR2M monthly water balance model across multiple subbasins
#' using a semi-distributed approach. It applies region-specific parameters and correction
#' factors for precipitation and evapotranspiration, and optionally saves the results as text files.
#'
#' @param Data A data frame of model input data in the airGR format, as produced by Create_Forcing_Inputs.
#' It must include columns: DatesR, P_1 to P_n, E_1 to E_n, and optionally Q.
#' @param Subbasins A SpatVector object containing the geometries of subbasins. It must include attributes 'COMID' (unique subbasin ID) and 'Region' (region name/code).
#' @param RunIni Simulation start date in the format 'mm/yyyy'.
#' @param RunEnd Simulation end date in the format 'mm/yyyy'.
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
#' library(GR2MSemiDistr)
#'
#' # Define region parameters
#' param_df <- data.frame(
#'   Region = unique(cat$Region),
#'   X1 = 2000,
#'   X2 = 1.8,
#'   fp = 0.8,
#'   fe = 0.8
#' )
#'
#' # Run the semi-distributed GR2M model
#' model <- Run_GR2MSemiDistr(
#'   Data = model_inputs,
#'   Subbasins = cat,
#'   RunIni = "01/1981",
#'   RunEnd = "12/2016",
#'   Parameters = param_df,
#' )
#'
#' # Select COMID of the subbasin to plot
#' target_comid <- model$COMID[1]  # Change index as needed
#' idx <- which(model$COMID == target_comid)
#' dates <- as.Date(model$Dates)
#'
#' # Plot model outputs for selected subbasin
#' par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
#'
#' plot(dates, model$PR[, idx], type = "l", col = "blue",
#'      main = paste("Precipitation (PR) - COMID", target_comid),
#'      xlab = "Date", ylab = "mm/month")
#'
#' plot(dates, model$AE[, idx], type = "l", col = "orange",
#'      main = paste("Actual Evapotranspiration (AE) - COMID", target_comid),
#'      xlab = "Date", ylab = "mm/month")
#'
#' plot(dates, model$SM[, idx], type = "l", col = "green4",
#'      main = paste("Production Store (SM) - COMID", target_comid),
#'      xlab = "Date", ylab = "mm/month")
#'
#' plot(dates, model$PC[, idx], type = "l", col = "purple",
#'      main = paste("Percolation (PC) - COMID", target_comid),
#'      xlab = "Date", ylab = "mm/month")
#'
#' plot(dates, model$RU[, idx], type = "l", col = "darkred",
#'      main = paste("Runoff (RU) - COMID", target_comid),
#'      xlab = "Date", ylab = "mm/month")
#'
#' plot(dates, model$QS[, idx], type = "l", col = "blue",
#'      main = paste("Simulated Discharge (QS) - COMID", target_comid),
#'      xlab = "Date", ylab = "m³/s")
#'
#' # Optional: compare simulated vs observed outlet discharge
#' if (!is.null(model$SINK)) {
#'   par(mfrow = c(1, 1))
#'   plot(dates, model$SINK$sim, type = "l", col = "blue",
#'        xlab = "Date", ylab = "Discharge [m³/s]",
#    main = "Simulated vs Observed Discharge (Outlet)")
#'   lines(dates, model$SINK$obs, col = "red", type='o')
#'   legend("topright", legend = c("Simulated", "Observed"),
#'          col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")
#' }
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
  comid  <- as.character(Subbasins$COMID)   # ensure character
  region <- Subbasins$Region
  area   <- terra::expanse(Subbasins, unit = "km")
  nsub   <- length(comid)

  # === Dates: accept Data$Date or Data$DatesR, parse robustly ===
  date_col <- if ("Date" %in% names(Data)) "Date" else if ("DatesR" %in% names(Data)) "DatesR" else NULL
  if (is.null(date_col)) stop("Data must contain a 'Date' or 'DatesR' column.")

  Data[[date_col]] <- as.POSIXct(
    paste0(Data[[date_col]], " 00:00:00"),
    tz = "UTC",
    tryFormats = c("%Y-%m-%d","%Y/%m/%d","%d-%m-%Y","%d/%m/%Y","%Y-%m","%m/%Y")
  )
  if (any(is.na(Data[[date_col]]))) stop("Unrecognized date format in Data date column.")

  # Indices for run window (mm/YYYY)
  ind_start <- which(format(Data[[date_col]], "%m/%Y") == RunIni)
  ind_end   <- which(format(Data[[date_col]], "%m/%Y") == RunEnd)
  if (length(ind_start) == 0 || length(ind_end) == 0) {
    stop("RunIni or RunEnd not found in Data (expected 'mm/YYYY').")
  }
  if (max(ind_end) < min(ind_start)) stop("RunEnd precedes RunIni.")

  Database <- Data[min(ind_start):max(ind_end), , drop = FALSE]
  Dates    <- as.Date(Database[[date_col]])
  ntime    <- length(Dates)
  nDays    <- lubridate::days_in_month(Dates)

  # === Validate Parameters ↔ Regions ===
  if (!is.data.frame(Parameters)) stop("Argument 'Parameters' must be a data frame.")
  required_cols <- c("Region","X1","X2","fp","fe")
  if (!all(required_cols %in% names(Parameters))) {
    stop("Parameters must contain columns: 'Region','X1','X2','fp','fe'.")
  }
  regs_needed <- sort(unique(region))
  regs_have   <- sort(unique(Parameters$Region))
  if (!identical(regs_needed, regs_have)) {
    stop(sprintf("Region mismatch. Needed: %s | Provided: %s",
                 paste(regs_needed, collapse=","), paste(regs_have, collapse=",")))
  }

  # === Check Data columns and map by name (not by position) ===
  p_names <- paste0("P_", comid)
  e_names <- paste0("E_", comid)
  missing_p <- setdiff(p_names, names(Database))
  missing_e <- setdiff(e_names, names(Database))
  if (length(missing_p) || length(missing_e)) {
    stop(sprintf("Data is missing required columns. Missing P: [%s]; Missing E: [%s]",
                 paste(missing_p, collapse=", "), paste(missing_e, collapse=", ")))
  }

  # === Helpers ===
  Subset_Param <- function(Param, RegionID) {
    c(Param$X1[Param$Region == RegionID],
      Param$X2[Param$Region == RegionID])
  }

  Forcing_Subbasin <- function(Param, RegionID, Database, comid_i) {
    factor_p <- Param$fp[Param$Region == RegionID]
    factor_e <- Param$fe[Param$Region == RegionID]
    data.frame(
      DatesR = as.POSIXct(Dates, tz = "UTC"),
      P = factor_p * Database[[paste0("P_", comid_i)]],
      E = factor_e * Database[[paste0("E_", comid_i)]]
    )
  }

  # === Run GR2M model ===
  message(sprintf("Running GR2M model for %d subbasins", nsub))

  ResModel <- vector("list", nsub)
  for (i in seq_len(nsub)) {
    if (i %% max(1, floor(nsub/10)) == 0) {
      message(sprintf("  • Progress: %d/%d subbasins...", i, nsub))
    }
    Param <- Subset_Param(Parameters, region[i])
    Input <- Forcing_Subbasin(Parameters, region[i], Database, comid[i])

    if (anyNA(Input$P) || anyNA(Input$E)) {
      stop(sprintf("NA values detected in P or E for subbasin %s", comid[i]))
    }
    if (any(Input$P < 0)) stop(sprintf("Negative P in subbasin %s", comid[i]))
    if (any(Input$E < 0)) stop(sprintf("Negative E in subbasin %s", comid[i]))

    # Edge case: single month → add a dummy next month (keeps lubridate)
    if (ntime == 1) {
      next_month <- lubridate::floor_date(Input$DatesR[1] + lubridate::months(1), "month")
      Input <- rbind(Input, data.frame(DatesR = next_month, P = 100, E = 100))
    }

    InputsModel <- CreateInputsModel(FUN_MOD  = RunModel_GR2M,
                                     DatesR   = Input$DatesR,
                                     Precip   = Input$P,
                                     PotEvap  = Input$E)

    if (is.null(IniState)) {
      RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                     InputsModel = InputsModel,
                                     IndPeriod_Run = seq_len(ntime),
                                     verbose = FALSE, warnings = FALSE)
    } else {
      IniStates <- CreateIniStates(FUN_MOD   = RunModel_GR2M,
                                   InputsModel = InputsModel,
                                   ProdStore = IniState[[i]]$Store$Prod,
                                   RoutStore = IniState[[i]]$Store$Rout,
                                   ExpStore  = IniState[[i]]$Store$Exp,
                                   UH1 = IniState[[i]]$UH$UH1,
                                   UH2 = IniState[[i]]$UH$UH2)
      RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                     InputsModel = InputsModel,
                                     IniStates = IniStates,
                                     IndPeriod_Run = seq_len(ntime),
                                     verbose = FALSE, warnings = FALSE)
    }

    ResModel[[i]] <- RunModel(InputsModel = InputsModel,
                              RunOptions  = RunOptions,
                              Param       = Subset_Param(Parameters, region[i]),
                              FUN         = RunModel_GR2M)
  }

  # === Post-processing ===
  extract_outputs <- function(res, area_km2, nDays_vec) {
    list(
      pr = round(res$Precip, 2),
      ae = round(res$AE, 2),
      sm = round(res$Prod, 2),
      pc = round(res$Perc, 2),
      ru = round(res$Qsim, 2),
      qs = round((area_km2 * res$Qsim) / (86.4 * nDays_vec), 2)  # m3/s to mm? (convention preserved)
    )
  }

  outputs <- lapply(seq_len(nsub), function(i) extract_outputs(ResModel[[i]], area[i], nDays))
  pr <- do.call(cbind, lapply(outputs, `[[`, "pr"))
  ae <- do.call(cbind, lapply(outputs, `[[`, "ae"))
  sm <- do.call(cbind, lapply(outputs, `[[`, "sm"))
  pc <- do.call(cbind, lapply(outputs, `[[`, "pc"))
  ru <- do.call(cbind, lapply(outputs, `[[`, "ru"))
  qs <- do.call(cbind, lapply(outputs, `[[`, "qs"))
  colnames(pr) <- colnames(ae) <- colnames(sm) <- colnames(pc) <- colnames(ru) <- colnames(qs) <- comid
  qt <- round(rowSums(qs), 2)

  # Warm-up trimming (optional)
  if (!is.null(WarmUp) && WarmUp > 0) {
    idx <- seq_len(WarmUp)
    if (length(idx) >= nrow(qs)) stop("WarmUp is greater than or equal to series length.")
    keep <- setdiff(seq_len(nrow(qs)), idx)
    pr <- pr[keep, , drop = FALSE]
    ae <- ae[keep, , drop = FALSE]
    sm <- sm[keep, , drop = FALSE]
    pc <- pc[keep, , drop = FALSE]
    ru <- ru[keep, , drop = FALSE]
    qs <- qs[keep, , drop = FALSE]
    qt <- qt[keep]
    Dates <- Dates[keep]
    if ("Q" %in% names(Database)) {
      qo <- Database$Q[keep]
    }
  } else if ("Q" %in% names(Database)) {
    qo <- Database$Q
  }

  Ans <- list(Dates = Dates, COMID = comid, QS = qs, RU = ru, PR = pr, AE = ae, SM = sm, PC = pc,
              EndState = lapply(ResModel, `[[`, "StateEnd"))

  if (exists("qo")) {
    Ans$SINK <- data.frame(sim = round(qt, 2), obs = round(qo, 2), row.names = Dates)
  }

  # === Save outputs (tab-separated; keep lubridate for month math) ===
  if (Save) {
    dir.create("./Outputs", showWarnings = FALSE, recursive = TRUE)

    yyyymm_new <- format(tail(Dates, 1), "%Y%m")
    prev_month <- lubridate::`%m-%`(lubridate::floor_date(tail(Dates, 1), "month"),
                                    lubridate::months(1))
    yyyymm_old <- format(prev_month, "%Y%m")

    if (Update) {
      for (tag in c("PR","AE","SM","PC","RU","QS")) {
        old_file <- sprintf("./Outputs/%s_GR2MSemiDistr_%s.txt", tag, yyyymm_old)
        if (file.exists(old_file)) file.remove(old_file)
      }
    }

    save_outputs <- function(var, tag) {
      df <- as.data.frame(var)
      colnames(df) <- paste0(tag, "_", comid)
      rownames(df) <- format(Dates, "%Y-%m-01")
      file <- sprintf("./Outputs/%s_GR2MSemiDistr_%s.txt", tag, yyyymm_new)
      write.table(df, file = file, sep = "\t", quote = FALSE)
    }
    save_outputs(pr, "PR")
    save_outputs(ae, "AE")
    save_outputs(sm, "SM")
    save_outputs(pc, "PC")
    save_outputs(ru, "RU")
    save_outputs(qs, "QS")
  }

  message("Processing completed successfully in...")
  tictoc::toc()
  return(Ans)
}
