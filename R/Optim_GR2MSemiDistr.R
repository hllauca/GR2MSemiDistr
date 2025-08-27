#' Parameter optimization for a semi-distributed GR2M model using the SCE-UA algorithm.
#'
#' Optimizes the GR2M monthly water balance model parameters (X1, X2, fp, fe) for multiple calibration regions
#' using the SCE-UA algorithm. The model is run in a semi-distributed setup over multiple subbasins grouped
#' by region. Users can select from various objective functions (e.g., NSE, KGE, RMSE, and linear combinations), including
#' composite metrics. Parameter bounds must be provided, and regions can be excluded from calibration.
#'
#' @param Data A data frame of model input data in the airGR format, as produced by Create_Forcing_Inputs.
#' It must include columns: DatesR, P_1 to P_n, E_1 to E_n, and Q (used for calibration).
#' @param Subbasins A SpatVector object containing the geometries of subbasins. It must include attributes "COMID" (unique subbasin ID) and "Region" (region name/code).
#' @param RunIni Simulation start date in the format "mm/yyyy".
#' @param RunEnd Simulation end date in the format "mm/yyyy".
#' @param WarmUp Optional number of months to discard from the beginning of the simulation as warm-up. Default is NULL.
#' @param MatrixTransfer Square transfer matrix (nsub × nsub) defining upstream-to-downstream connectivity for routing/accumulation of subbasin flows.
#' @param Outlet Optional. Outlet subbasin identifier as a COMID code present in Subbasins. If provided, the routed outlet series is extracted from this column.
#' @param Parameters A data frame containing GR2M model parameters and correction factors per region. Must have columns: Region, X1, X2, fp, fe.
#' @param Parameters.Min Numeric vector of length 4, specifying lower bounds for optimization in the following order:
#' c(X1, X2, fp, fe).
#' @param Parameters.Max Numeric vector of length 4, specifying upper bounds for optimization in the same order as `Parameters.Min`.
#' @param Max.Functions Integer. Maximum number of function evaluations in the optimization. Default is 1500.
#' @param Optimization Character. Objective function to optimize. Available options include:
#' \describe{
#'   \item{"OF1"}{Kling-Gupta Efficiency (KGE).}
#'   \item{"OF2"}{Nash-Sutcliffe Efficiency (NSE).}
#'   \item{"OF3"}{NSE applied to log-transformed flows (log-NSE, Pushpalatha 2012).}
#'   \item{"OF4"}{Root Mean Square Error (RMSE).}
#'   \item{"OF5"}{Percent Bias (PBIAS).}
#'   \item{"OF6"}{FDC Bias (FDC Bias), based on flow duration curves.}
#'   \item{"OF7"}{Pearson correlation coefficient (r).}
#'   \item{"OF8"}{Composite objective: w1 × KGE + w2 × NSE - w3 × RMSE.}
#'   \item{"OF9"}{Composite objective: w1 × KGE + w2 × log-NSE - w3 × RMSE.}
#'   \item{"OF10"}{Composite objective: w1 × KGE.km + w2 × KGE.lf - w3 × RMSE.}
#' }
#' @param No.Optim Character vector (optional). Region names to exclude from the optimization (parameters will be kept fixed).
#' @param w1 Numeric. Weight for the first component in composite objective functions (OF8, OF9, OF10). Default is 0.6.
#' @param w2 Numeric. Weight for the second component in composite objective functions. Default is 0.3.
#' @param w3 Numeric. Weight for the third component in composite objective functions. Default is 0.2.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{Parameters}{Data frame of optimized parameters per region (columns: Region, X1, X2, fp, fe).}
#'   \item{OF}{Final value of the selected objective function.}
#' }
#'
#' @references
#' Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. (2021).
#' PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020).
#' \emph{Water}, 13(8), 1048. \doi{10.3390/w13081048}
#'
#' @examples
#' library(GR2MSemiDistr)
#'
#' # Define initial parameters for each region
#' param_init <- data.frame(
#'   Region = unique(cat$Region),
#'   X1 = rep(500, length(unique(cat$Region))),
#'   X2 = rep(1.5, length(unique(cat$Region))),
#'   fp = rep(1.0, length(unique(cat$Region))),
#'   fe = rep(1.0, length(unique(cat$Region)))
#' )
#'
#' # Optimize parameters using the KGE criterion (OF1)
#' result <- Optim_GR2MSemiDistr(
#'   Data = model_inputs,
#'   Subbasins = cat,
#'   RunIni = "01/1981",
#'   RunEnd = "12/2005",
#'   Parameters = param_init,
#'   Optimization = "OF1"  # KGE
#' )
#'
#' # Extract results
#' best_params <- result$Parameters
#' final_score <- result$OF
#'
#' # Plot observed vs simulated discharge at the outlet
#' model <- Run_GR2MSemiDistr(
#'   Data = model_inputs,
#'   Subbasins = cat,
#'   RunIni = "01/1981",
#'   RunEnd = "12/2005",
#'   Parameters = best_params,
#' )
#'
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
#' @import rtop
#' @import hydroGOF
#' @import terra
#' @import tictoc
#' @import lubridate
#'
#' @export
#'
Optim_GR2MSemiDistr <- function(Data,
                                Subbasins,
                                RunIni,
                                RunEnd,
                                Parameters,
                                MatrixTransfer,
                                Outlet,
                                Parameters.Min = c(100, 0.1, 0.8, 0.8),
                                Parameters.Max = c(2000, 10, 1.2, 1.2),
                                Optimization = 'OF1',
                                Max.Functions = 1000,
                                WarmUp = 36,
                                No.Optim = NULL,
                                w1 = 0.6,
                                w2 = 0.3,
                                w3 = 0.2) {
  tictoc::tic()

  # === Basic checks ===
  if (!inherits(Subbasins, "SpatVector")) stop("Subbasins must be a 'SpatVector'.")
  if (!all(c("COMID", "Region") %in% names(Subbasins))) {
    stop("Subbasins must contain fields 'COMID' and 'Region'.")
  }
  comid  <- as.character(Subbasins$COMID)
  region <- Subbasins$Region
  area   <- terra::expanse(Subbasins, unit = 'km')
  nsub   <- length(comid)

  # Validate MatrixTransfer
  if (!is.matrix(MatrixTransfer) || any(dim(MatrixTransfer) != nsub)) {
    stop("MatrixTransfer must be a square matrix with dimensions equal to the number of subbasins.")
  }
  # Validate Outlet
  idx_outlet <- match(Outlet, comid)
  if (is.na(idx_outlet)) stop(sprintf("Outlet COMID '%s' not found in Subbasins.", Outlet))

  # === Dates: accept Data$Date or Data$DatesR ===
  date_col <- if ("Date" %in% names(Data)) "Date" else if ("DatesR" %in% names(Data)) "DatesR" else NULL
  if (is.null(date_col)) stop("Data must include a 'Date' or 'DatesR' column.")
  Data[[date_col]] <- as.POSIXct(
    paste0(Data[[date_col]], " 00:00:00"),
    tz = "UTC",
    tryFormats = c("%Y-%m-%d","%Y/%m/%d","%d-%m-%Y","%d/%m/%Y","%Y-%m","%m/%Y")
  )
  if (any(is.na(Data[[date_col]]))) stop("Unrecognized date format in Data.")

  # Indices for run window (mm/YYYY)
  ind_start <- which(format(Data[[date_col]], "%m/%Y") == RunIni)
  ind_end   <- which(format(Data[[date_col]], "%m/%Y") == RunEnd)
  if (!length(ind_start) || !length(ind_end)) stop("RunIni or RunEnd not found (mm/YYYY).")
  if (max(ind_end) < min(ind_start)) stop("RunEnd precedes RunIni.")

  Database <- Data[min(ind_start):max(ind_end), , drop = FALSE]
  if (!("Q" %in% names(Database))) stop("Data must include observed flow column 'Q' for calibration.")
  Dates    <- as.Date(Database[[date_col]])
  ntime    <- length(Dates)
  nDays    <- lubridate::days_in_month(Dates)

  # === Regions ↔ Parameters ===
  if (!is.data.frame(Parameters)) stop("Parameters must be a data frame.")
  req_cols <- c("Region","X1","X2","fp","fe")
  if (!all(req_cols %in% names(Parameters))) {
    stop("Parameters must contain columns: Region, X1, X2, fp, fe.")
  }
  all_regions <- sort(unique(region))
  regs_have   <- sort(unique(Parameters$Region))
  if (!identical(all_regions, regs_have)) {
    stop(sprintf("Region mismatch. Needed: %s | Provided: %s",
                 paste(all_regions, collapse=","), paste(regs_have, collapse=",")))
  }

  # Regions to optimize
  if (!is.null(No.Optim) && any(!(No.Optim %in% all_regions))) {
    stop("No.Optim contains regions not present in Subbasins.")
  }
  opt_regions <- if (is.null(No.Optim)) all_regions else setdiff(all_regions, No.Optim)
  if (!length(opt_regions)) stop("No regions left to optimize (check No.Optim).")

  # === Forcing data columns check ===
  p_names <- paste0("P_", comid)
  e_names <- paste0("E_", comid)
  miss_p <- setdiff(p_names, names(Database))
  miss_e <- setdiff(e_names, names(Database))
  if (length(miss_p) || length(miss_e)) {
    stop(sprintf("Missing columns in Data. P: [%s]; E: [%s]",
                 paste(miss_p, collapse=","), paste(miss_e, collapse=",")))
  }

  # === Bounds sanity ===
  if (!(length(Parameters.Min) == 4 && length(Parameters.Max) == 4)) {
    stop("Parameters.Min/Max must be length 4 (X1, X2, fp, fe).")
  }

  # Starting vector and bounds for opt regions (X1,X2,fp,fe) per region
  opt.param     <- unlist(Parameters[Parameters$Region %in% opt_regions, c("X1","X2","fp","fe")], use.names = FALSE)
  opt.param.min <- rep(Parameters.Min, times = length(opt_regions))
  opt.param.max <- rep(Parameters.Max, times = length(opt_regions))
  stopifnot(length(opt.param) == length(opt.param.min),
            length(opt.param) == length(opt.param.max))

  # === Helpers ===
  get_param <- function(param_vec, regions) {
    n <- length(regions)
    data.frame(Region = regions,
               X1 = param_vec[1:n],
               X2 = param_vec[(n + 1):(2 * n)],
               fp = param_vec[(2 * n + 1):(3 * n)],
               fe = param_vec[(3 * n + 1):(4 * n)],
               row.names = NULL)
  }

  forcing_input <- function(Param, RegionID, Database, comid_i) {
    fp <- Param$fp[Param$Region == RegionID]
    fe <- Param$fe[Param$Region == RegionID]
    data.frame(DatesR = as.POSIXct(Dates, tz = "UTC"),
               P = fp * Database[[paste0("P_", comid_i)]],
               E = fe * Database[[paste0("E_", comid_i)]])
  }

  # === Build cumulative routing matrix S = I + T + T^2 + ... + T^(n-1) ===
  build_S <- function(Tmat) {
    n <- ncol(Tmat)
    S <- diag(n)
    TT <- Tmat
    for (k in 1:(n - 1)) {
      S  <- S + TT
      TT <- TT %*% Tmat
    }
    S
  }
  S_route <- build_S(MatrixTransfer)  # precalculated outside OFUN

  # === Objective function using QT (outlet discharge after routing) ===
  OFUN <- function(par) {
    Param <- get_param(par, opt_regions)
    # Merge fixed (No.Optim) if any
    if (!is.null(No.Optim)) {
      fixed <- Parameters[Parameters$Region %in% No.Optim, ]
      Param <- rbind(Param, fixed[, c("Region","X1","X2","fp","fe")])
      Param <- Param[match(all_regions, Param$Region), ]
    } else {
      Param <- Param[match(all_regions, Param$Region), ]
    }

    # Run GR2M for each subbasin (before routing)
    QS <- matrix(NA_real_, nrow = ntime, ncol = nsub)
    for (i in seq_len(nsub)) {
      reg_i   <- region[i]
      param_i <- c(Param$X1[Param$Region == reg_i], Param$X2[Param$Region == reg_i])
      input_i <- forcing_input(Param, reg_i, Database, comid[i])

      model_input <- CreateInputsModel(FUN_MOD = RunModel_GR2M,
                                       DatesR  = input_i$DatesR,
                                       Precip  = input_i$P,
                                       PotEvap = input_i$E)

      run_opt <- CreateRunOptions(FUN_MOD = RunModel_GR2M,
                                  InputsModel   = model_input,
                                  IndPeriod_Run = seq_len(ntime),
                                  verbose = FALSE, warnings = FALSE)

      output <- RunModel(model_input, run_opt, param_i, RunModel_GR2M)

      # Convert to m3/s
      QS[, i] <- (area[i] * output$Qsim) / (86.4 * nDays)
    }

    # Routing: QR = QS %*% t(S), then extract outlet
    QR <- QS %*% t(S_route)
    Qsim <- as.numeric(QR[, idx_outlet])

    # Observed flow
    Qobs <- Database$Q

    # Warm-up trimming
    if (!is.null(WarmUp) && WarmUp > 0) {
      if (WarmUp >= length(Qsim)) return(Inf)
      Qsim <- Qsim[-seq_len(WarmUp)]
      Qobs <- Qobs[-seq_len(WarmUp)]
    }

    # Remove NA pairs
    ok <- is.finite(QT) & is.finite(Qobs)
    if (!any(ok)) return(Inf)
    y <- QT[ok]; x <- Qobs[ok]

    # Metrics
    kge     <- 1 - KGE(y, x)
    nse     <- 1 - NSE(y, x)
    nselog  <- 1 - NSE(y, x, fun = log, epsilon.type = "Pushpalatha2012")
    rmse_v  <- rmse(y, x)
    pbias_v <- abs(pbias(y, x))
    pbiasfdc_v <- abs(pbiasfdc(y, x, plot = FALSE))
    rpear_v <- 1 - rPearson(y, x)
    kgekm_v <- 1 - KGEkm(y, x)
    kgelf_v <- 1 - KGElf(y, x)

    criteria <- c(
      OF1  = kge,
      OF2  = nse,
      OF3  = nselog,
      OF4  = rmse_v,
      OF5  = pbias_v,
      OF6  = pbiasfdc_v,
      OF7  = rpear_v,
      OF8  = w1 * kge + w2 * nse + w3 * rmse_v,
      OF9  = w1 * kge + w2 * nselog + w3 * rmse_v,
      OF10 = w1 * kgekm_v + w2 * kgelf_v + w3 * rmse_v
    )
    criteria[Optimization]
  }

  message(sprintf("Optimizing %s with SCE-UA over %d regions (%d params)...",
                  Optimization, length(opt_regions), length(opt_regions) * 4))

  Calibration <- sceua(OFUN,
                       pars  = opt.param,
                       lower = opt.param.min,
                       upper = opt.param.max,
                       maxn  = Max.Functions)

  # === Organize output ===
  final_params <- get_param(Calibration$par, opt_regions)
  if (!is.null(No.Optim)) {
    fixed <- Parameters[Parameters$Region %in% No.Optim, c("Region","X1","X2","fp","fe")]
    final_params <- rbind(final_params, fixed)
    final_params <- final_params[match(all_regions, final_params$Region), ]
  } else {
    final_params <- final_params[match(all_regions, final_params$Region), ]
  }

  # Report objective value:
  final_obj_raw <- Calibration$value

  message("Optimization complete.")
  cat("\nFinal calibrated parameters per region:\n")
  print(final_params)
  cat(sprintf("Objective Function (%s) = %s\n", Optimization, final_obj_rep))

  tictoc::toc()
  list(Parameters = final_params,
       OF = final_obj_raw)
}
