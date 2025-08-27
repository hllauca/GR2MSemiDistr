#' Run the GR2M model for multiple subbasins (semi-distributed with routing)
#'
#' This function executes the monthly GR2M water balance model across multiple
#' subbasins using a semi-distributed approach. It applies region-specific
#' parameters and correction factors for precipitation and potential
#' evapotranspiration, routes the simulated flows through the river network via
#' a transfer matrix, and can optionally save results as text files.
#'
#' @param Data A data frame of model input data in the airGR format, as produced by Create_Forcing_Inputs. It must include columns: DatesR, P_1 to P_n, E_1 to E_n, and Q (used for calibration).
#' @param Subbasins A SpatVector object containing the geometries of subbasins. It must include attributes "COMID" (unique subbasin ID) and "Region" (region name/code).
#' @param RunIni Simulation start date in the format "mm/yyyy".
#' @param RunEnd Simulation end date in the format "mm/yyyy".
#' @param Parameters A data frame containing GR2M model parameters and correction factors per region. Must have columns: Region, X1, X2, fp, fe.
#' @param MatrixTransfer Square matrix or data frame (nsub × nsub) of subbasin connectivity. Row and column names must match COMID and are internally reordered to align with P/E inputs and QS.
#' @param Outlet Optional. Outlet subbasin identifier as a COMID code present in Subbasins. If provided, the routed outlet series is extracted from this column.
#' @param WarmUp Optional number of months to discard from the beginning of the simulation as warm-up. Default is NULL.
#' @param IniState Optional list of initial states for each subbasin. Default is NULL.
#' @param Save Logical. If TRUE, output time series are saved as tab-separated text files under "Outputs/".
#' @param Update Logical. If TRUE and Save = TRUE, removes previous month's output files before saving the current ones.
#'
#' @return A list with:
#' \describe{
#'   \item{PR}{Matrix of total precipitation [mm/month] per subbasin.}
#'   \item{AE}{Matrix of actual evapotranspiration [mm/month] per subbasin.}
#'   \item{SM}{Matrix of production store level [mm/month] per subbasin.}
#'   \item{PC}{Matrix of percolation [mm/month] per subbasin.}
#'   \item{RU}{Matrix of runoff [mm/month] per subbasin (model internal output).}
#'   \item{QR}{Matrix of routed discharge [m³/s] per subbasin after network accumulation using MatrixTransfer.}
#'   \item{Dates}{Date vector for the simulation period.}
#'   \item{COMID}{Vector of subbasin COMIDs (column order matches matrices).}
#'   \item{EndState}{List of final state variables for each subbasin.}
#'   \item{SINK}{Only if outlet is provided.
#'               Data frame with outlet series where sim is the routed outlet discharge
#'               (column QR[, Outlet]). If observed outlet discharge QR exists
#'               in Data, obs is included as well. Row names are Dates.}
#' }
#'
#' @details
#' Subbasin-level simulated discharge (in m³/s) is obtained from GR2M runoff
#' (mm/month) using basin area and days per month, and then accumulated/routed
#' through the network using the geometric series of the transfer matrix
#' S = I + T + T^2 + ... + T^(n-1), applied as QR = QS * S'.
#' The outlet series is taken from the routed matrix column corresponding to
#' the provided outlet COMID.
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
#'   MatrixTransfer = matrixT,
#'   Outlet = '8', # basin outlet comid
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
#' plot(dates, model$QR[, idx], type = "l", col = "blue",
#'      main = paste("Simulated Discharge (QR) - COMID", target_comid),
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
                              MatrixTransfer,
                              Outlet = NULL,
                              WarmUp = NULL,
                              IniState = NULL,
                              Save = FALSE,
                              Update = FALSE) {

  tictoc::tic()

  # ==== Basic validation for Subbasins ====
  if (!inherits(Subbasins, "SpatVector")) {
    stop("'Subbasins' must be a 'SpatVector'.")
  }
  if (!all(c("COMID", "Region") %in% names(Subbasins))) {
    stop("'Subbasins' must contain fields 'COMID' and 'Region'.")
  }

  comid  <- as.character(Subbasins$COMID)
  region <- Subbasins$Region
  area   <- terra::expanse(Subbasins, unit = "km")
  nsub   <- length(comid)

  # ==== Dates: detect column, coerce and slice ====
  date_col <- if ("Date" %in% names(Data)) "Date" else if ("DatesR" %in% names(Data)) "DatesR" else NULL
  if (is.null(date_col)) stop("Data must contain 'Date' or 'DatesR'.")

  # Convertir a POSIXct robustamente (si ya es Date/POSIX, se respeta)
  if (inherits(Data[[date_col]], "Date")) {
    Data[[date_col]] <- as.POSIXct(Data[[date_col]], tz = "UTC")
  } else if (!inherits(Data[[date_col]], c("POSIXct", "POSIXt"))) {
    Data[[date_col]] <- as.POSIXct(
      paste0(Data[[date_col]], " 00:00:00"),
      tz = "UTC",
      tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%d-%m-%Y", "%d/%m/%Y", "%Y-%m", "%m/%Y")
    )
  }
  if (any(is.na(Data[[date_col]]))) stop("Unrecognized date format in Data.")

  ind_start <- which(format(Data[[date_col]], "%m/%Y") == RunIni)
  ind_end   <- which(format(Data[[date_col]], "%m/%Y") == RunEnd)
  if (!length(ind_start) || !length(ind_end)) stop("RunIni or RunEnd not found (expected 'mm/YYYY').")
  if (max(ind_end) < min(ind_start))          stop("RunEnd precedes RunIni.")

  Database <- Data[min(ind_start):max(ind_end), , drop = FALSE]
  Dates    <- as.Date(Database[[date_col]])
  ntime    <- length(Dates)
  nDays    <- lubridate::days_in_month(Dates)

  # ==== Parameters: presence, coverage and mapping ====
  if (!is.data.frame(Parameters)) stop("'Parameters' must be a data frame.")
  req_cols <- c("Region", "X1", "X2", "fp", "fe")
  if (!all(req_cols %in% names(Parameters))) {
    stop("Parameters must contain: Region, X1, X2, fp, fe.")
  }

  regs_needed <- sort(unique(region))
  regs_have   <- sort(unique(Parameters$Region))
  if (!identical(regs_needed, regs_have)) {
    stop(sprintf("Region mismatch. Needed: %s | Provided: %s",
                 paste(regs_needed, collapse = ","),
                 paste(regs_have,   collapse = ",")))
  }

  match_idx <- match(region, Parameters$Region)
  X1_vec <- Parameters$X1[match_idx]
  X2_vec <- Parameters$X2[match_idx]
  fp_vec <- Parameters$fp[match_idx]
  fe_vec <- Parameters$fe[match_idx]

  # ==== Forcing columns: check presence and order vs COMID ====
  p_names <- paste0("P_", comid)
  e_names <- paste0("E_", comid)

  miss_p  <- setdiff(p_names, names(Database))
  miss_e  <- setdiff(e_names, names(Database))
  if (length(miss_p) || length(miss_e)) {
    stop(sprintf("Missing columns. P: [%s]; E: [%s]",
                 paste(miss_p, collapse = ", "),
                 paste(miss_e, collapse = ", ")))
  }

  # Verificación estricta de orden (que el orden de P_ y E_ siga a COMID)
  p_idx <- match(p_names, names(Database))
  e_idx <- match(e_names, names(Database))
  if (is.unsorted(p_idx) || is.unsorted(e_idx)) {
    warning("Las columnas de P_ y/o E_ no están en el mismo orden que COMID. Se continuará accediendo por nombre, pero se recomienda mantener el orden para coherencia.")
  }

  # ==== MatrixTransfer: aceptar data.frame/matrix, validar nombres y reordenar ====
  if (is.null(MatrixTransfer)) stop("MatrixTransfer is required.")

  MT <- as.matrix(MatrixTransfer)

  if (any(dim(MT) != c(nsub, nsub))) {
    stop(sprintf("MatrixTransfer must be a square %d x %d object.", nsub, nsub))
  }
  if (is.null(rownames(MT)) || is.null(colnames(MT))) {
    stop("MatrixTransfer must have rownames and colnames corresponding to COMID.")
  }

  rows_mt <- as.character(rownames(MT))
  cols_mt <- as.character(colnames(MT))

  miss_rows <- setdiff(comid, rows_mt)
  miss_cols <- setdiff(comid, cols_mt)
  if (length(miss_rows) || length(miss_cols)) {
    stop(sprintf(
      "MatrixTransfer names do not match COMID. Missing rows: [%s]; Missing cols: [%s]",
      paste(miss_rows, collapse = ", "), paste(miss_cols, collapse = ", ")
    ))
  }
  # Reordenar al orden exacto de COMID (coincidirá con P_, E_, y QS)
  MT <- MT[comid, comid, drop = FALSE]

  # ==== Model run per subbasin ====
  message(sprintf("Running GR2M model for %d subbasins", nsub))
  ResModel <- vector("list", nsub)

  for (i in seq_len(nsub)) {
    if (i %% max(1, floor(nsub / 10)) == 0) {
      message(sprintf("  • Progress: %d/%d...", i, nsub))
    }

    Input <- data.frame(
      DatesR = as.POSIXct(Dates, tz = "UTC"),
      P      = fp_vec[i] * Database[[p_names[i]]],
      E      = fe_vec[i] * Database[[e_names[i]]]
    )

    if (anyNA(Input$P) || anyNA(Input$E)) stop(sprintf("NA in P or E for COMID %s", comid[i]))
    if (any(Input$P < 0))                 stop(sprintf("Negative P in COMID %s", comid[i]))
    if (any(Input$E < 0))                 stop(sprintf("Negative E in COMID %s", comid[i]))

    # Si solo hay 1 paso temporal, rellenar con un mes extra para satisfacer airGR
    if (ntime == 1) {
      d1 <- as.Date(Input$DatesR[1])
      next_month <- seq(d1, by = "month", length.out = 2)[2]
      Input <- rbind(Input, data.frame(DatesR = as.POSIXct(next_month, tz = "UTC"), P = 100, E = 100))
    }

    InputsModel <- CreateInputsModel(
      FUN_MOD = RunModel_GR2M,
      DatesR  = Input$DatesR,
      Precip  = Input$P,
      PotEvap = Input$E
    )

    if (is.null(IniState)) {
      RunOptions <- CreateRunOptions(
        FUN_MOD       = RunModel_GR2M,
        InputsModel   = InputsModel,
        IndPeriod_Run = seq_len(ntime),
        verbose = FALSE, warnings = FALSE
      )
    } else {
      IniStates <- CreateIniStates(
        FUN_MOD     = RunModel_GR2M,
        InputsModel = InputsModel,
        ProdStore   = IniState[[i]]$Store$Prod,
        RoutStore   = IniState[[i]]$Store$Rout,
        ExpStore    = IniState[[i]]$Store$Exp,
        UH1         = IniState[[i]]$UH$UH1,
        UH2         = IniState[[i]]$UH$UH2
      )
      RunOptions <- CreateRunOptions(
        FUN_MOD       = RunModel_GR2M,
        InputsModel   = InputsModel,
        IniStates     = IniStates,
        IndPeriod_Run = seq_len(ntime),
        verbose = FALSE, warnings = FALSE
      )
    }

    ResModel[[i]] <- RunModel(
      InputsModel = InputsModel,
      RunOptions  = RunOptions,
      Param       = c(X1_vec[i], X2_vec[i]),
      FUN         = RunModel_GR2M
    )
  }

  # ==== Helper: bind matrices by key ====
  .to_mat <- function(key) {
    mats <- lapply(ResModel, function(r) round(r[[key]], 2))
    do.call(cbind, mats)
  }

  # ==== Collect per-subbasin outputs (model scale) ====
  PR <- .to_mat("Precip")
  AE <- .to_mat("AE")
  SM <- .to_mat("Prod")
  PC <- .to_mat("Perc")
  RU <- .to_mat("Qsim")

  # ==== Convert to discharge m3/s using basin area and days/month ====
  QS <- do.call(cbind, lapply(seq_len(nsub), function(i) {
    round((area[i] * ResModel[[i]]$Qsim) / (86.4 * nDays), 2)
  }))

  colnames(PR) <- colnames(AE) <- colnames(SM) <- colnames(PC) <- colnames(RU) <- colnames(QS) <- comid

  # ==== Routing with MatrixTransfer (ya reordenado al orden de COMID) ====
  n  <- ncol(QS)
  S  <- diag(n)
  TT <- MT
  for (k in 1:(n - 1)) {
    S  <- S + TT
    TT <- TT %*% MT
  }
  QR <- QS %*% t(S)
  colnames(QR) <- comid
  rownames(QR) <- NULL

  # ==== Observed discharge (optional) ====
  if ("Q" %in% names(Database)) qo <- Database$Q

  # ==== Outlet discharge (optional) ====
  QT <- NULL
  if (!is.null(Outlet)) {
    idx_outlet <- match(Outlet, comid)
    if (is.na(idx_outlet)) stop(sprintf("Outlet COMID '%s' not found in Subbasins.", Outlet))
    QT <- round(QR[, idx_outlet], 2)
  }

  # ==== Warm-up trimming (optional) ====
  if (!is.null(WarmUp) && WarmUp > 0) {
    if (WarmUp >= nrow(QR)) stop("WarmUp is greater than or equal to the series length.")
    keep <- (WarmUp + 1):nrow(QR)
    slice <- function(M) M[keep, , drop = FALSE]
    PR <- slice(PR); AE <- slice(AE); SM <- slice(SM); PC <- slice(PC); RU <- slice(RU); QR <- slice(QR)
    Dates <- Dates[keep]
    if (!is.null(QT)) QT <- QT[keep]
    if (exists("qo")) qo <- qo[keep]
  }

  # ==== Assemble answer ====
  Ans <- list(
    Dates    = Dates,
    COMID    = comid,
    PR       = PR,
    AE       = AE,
    SM       = SM,
    PC       = PC,
    RU       = RU,
    QR       = QR,
    EndState = lapply(ResModel, `[[`, "StateEnd")
  )

  if (!is.null(Outlet)) {
    if (exists("qo")) {
      Ans$SINK <- data.frame(sim = QT, obs = round(qo, 2), row.names = Dates)
    } else {
      Ans$SINK <- data.frame(sim = QT, row.names = Dates)
    }
  }

  # ==== Optional: save to disk ====
  if (Save) {
    if (!dir.exists("./Outputs")) dir.create("./Outputs", recursive = TRUE)

    yyyymm_new <- format(tail(Dates, 1), "%Y%m")
    d0         <- as.Date(format(tail(Dates, 1), "%Y-%m-01"))
    prev_month <- seq(d0, by = "-1 month", length.out = 2)[2]
    yyyymm_old <- format(prev_month, "%Y%m")

    if (Update) {
      for (tag in c("PR", "AE", "SM", "PC", "RU", "QR")) {
        old_file <- sprintf("./Outputs/%s_GR2MSemiDistr_%s.txt", tag, yyyymm_old)
        if (file.exists(old_file)) file.remove(old_file)
      }
    }

    .save_outputs <- function(var, tag) {
      df <- as.data.frame(var)
      colnames(df) <- paste0(tag, "_", comid)
      rownames(df) <- format(Dates, "%Y-%m-01")
      file <- sprintf("./Outputs/%s_GR2MSemiDistr_%s.txt", tag, yyyymm_new)
      write.table(df, file = file, sep = "\t", quote = FALSE)
    }
    .save_outputs(PR, "PR")
    .save_outputs(AE, "AE")
    .save_outputs(SM, "SM")
    .save_outputs(PC, "PC")
    .save_outputs(RU, "RU")
    .save_outputs(QR, "QR")
  }

  message("Processing completed successfully in...")
  tictoc::toc()
  return(Ans)
}
