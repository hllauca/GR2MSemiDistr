#' Parameter optimization for a semi-distributed GR2M model using the SCE-UA algorithm.
#'
#' Optimizes the GR2M monthly water balance model parameters (X1, X2, fp, fe) for multiple calibration regions
#' using the SCE-UA algorithm. The model is run in a semi-distributed setup over multiple subbasins grouped
#' by region. Users can select from various objective functions, including composite metrics.
#' Parameter bounds must be provided, and regions can be excluded from calibration.
#'
#' @param Data data.frame. Model input data in the airGR format, as produced by `Create_Forcing_Inputs`.
#' Must include columns: DatesR, P_1 to P_n, E_1 to E_n, and Q (used for calibration).
#' @param Subbasins SpatVector. Geometries of subbasins. Must include attributes "COMID" (unique subbasin ID) and "Region" (region name/code).
#' @param RunIni character. Simulation start date in the format "mm/yyyy".
#' @param RunEnd character. Simulation end date in the format "mm/yyyy".
#' @param TransferMatrix matrix or dgCMatrix. Subbasin connectivity matrix defining downstream routing.
#' Rows represent receiving (downstream) subbasins, and columns represent donor (upstream) subbasins.
#' Row and column names must match `COMID`. For large river networks (thousands of subbasins),
#' a sparse matrix is strongly recommended for memory efficiency.
#' @param Outlet character. Outlet subbasin identifier as a COMID code present in `Subbasins`. Required for calibration.
#' @param Parameters data.frame. GR2M model parameters and correction factors per region.
#' Must have columns: Region, X1, X2, fp, fe. An optional `k` column (linear-reservoir
#' routing constant per region, in months; see `Run_GR2MSemiDistr()` Details) is used
#' as-is for regions in `No.Optim`; if absent, `k = 0` is assumed (previous,
#' instantaneous-routing behavior).
#' @param Parameters.Min numeric vector of length 5. Lower bounds for optimization
#' in the following order: c(X1, X2, fp, fe, k).
#' @param Parameters.Max numeric vector of length 5. Upper bounds for optimization
#' in the same order as `Parameters.Min`.
#' @param Max.Functions integer. Maximum number of function evaluations in the optimization.
#' Default is 1000.
#' @param Optimization character. Objective function to optimize. Available options include:
#' \describe{
#'   \item{"OF1"}{Kling-Gupta Efficiency (KGE).}
#'   \item{"OF2"}{Nash-Sutcliffe Efficiency (NSE).}
#'   \item{"OF3"}{NSE applied to log-transformed flows (logNSE, Pushpalatha 2012).}
#'   \item{"OF4"}{Root Mean Square Error (RMSE).}
#'   \item{"OF5"}{Percent Bias (PBIAS).}
#'   \item{"OF6"}{Bias in flow duration curves (FDC Bias).}
#'   \item{"OF7"}{Pearson correlation coefficient (r).}
#'   \item{"OF8"}{Composite objective combining KGE, NSE, and RSR (RMSE normalized by
#'   sd(Qobs), Moriasi et al. 2007) with weights (w1, w2, w3).}
#'   \item{"OF9"}{Composite objective combining KGE, logNSE, and RSR with weights.}
#'   \item{"OF10"}{Composite objective combining KGE.km, KGE.lf, and RSR with weights.}
#' }
#' @param WarmUp integer. Number of months discarded as warm-up period when computing the
#' objective function. Default is 36.
#' @param No.Optim character vector (optional). Region names to exclude from the optimization (parameters will be kept fixed).
#' @param w1 numeric. Weight for the first component in composite objective functions (OF8, OF9, OF10). Default is 0.6.
#' @param w2 numeric. Weight for the second component in composite objective functions. Default is 0.3.
#' @param w3 numeric. Weight for the third component in composite objective functions. Default is 0.2.
#' @param Cores integer. Number of parallel workers used to run GR2M across
#' subbasins on every objective-function evaluation. Default is `1`
#' (sequential); values greater than 1 build a `parallel::makeCluster()`
#' PSOCK cluster ONCE before optimization starts and reuse it for every
#' evaluation. Only worth it for large networks (thousands of subbasins).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{Parameters}{data.frame. Optimized parameter set per calibration region.
#'   Always includes all regions present in `Subbasins`; regions excluded from optimization
#'   via `No.Optim` retain their original parameter values. Columns are: Region, X1, X2, fp, fe, k.}
#'   \item{OF}{numeric. Final value of the selected objective function corresponding
#'   to the best parameter set found by the optimization.}
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
#'   X1  = 500,  # Production store capacity [mm]
#'   X2  = 1.5,  # Groundwater exchange coefficient
#'   fp  = 1.0,  # Precipitation correction factor
#'   fe  = 1.0   # Evapotranspiration correction factor
#' )
#'
#' # Calibrate parameters using OF10 as objective function
#' result <- Optim_GR2MSemiDistr(
#'   Data = model_inputs,
#'   Subbasins = cat,
#'   RunIni = "01/1981",
#'   RunEnd = "12/2005",
#'   TransferMatrix = matrixT,
#'   Outlet = '8', # basin outlet comid
#'   Parameters = param_init,
#'   Optimization = "OF10"
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
#'   TransferMatrix = matrixT,
#'   Outlet = '8', # basin outlet comid
#'   Parameters = best_params,
#' )
#'
#'   par(mfrow = c(1, 1))
#'   plot(model$Dates, model$SINK$sim, type = "l", col = "blue", lwd = 1.5,
#'        xlab = "Date", ylab = "Discharge [m³/s]",
#'        main = "Simulated vs Observed Discharge (Outlet)")
#'     lines(model$Dates, model$SINK$obs, col = "red", lty = 2, lwd = 1.5)
#'     legend("topright", legend = c("Simulated", "Observed"),
#'            col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")
#'
#' @importFrom airGR CreateInputsModel CreateRunOptions RunModel RunModel_GR2M
#' @importFrom hydroGOF KGE NSE rPearson KGEkm KGElf rmse pbias pbiasfdc
#'
#' @export
Optim_GR2MSemiDistr <- function(Data,
                                Subbasins,
                                RunIni,
                                RunEnd,
                                Parameters,
                                TransferMatrix,
                                Outlet,
                                Parameters.Min = c(100, 0.1, 0.8, 0.8, 0),
                                Parameters.Max = c(2000, 10, 1.2, 1.2, 6),
                                Optimization = 'OF10',
                                Max.Functions = 1000,
                                WarmUp = 36,
                                No.Optim = NULL,
                                w1 = 0.6,
                                w2 = 0.3,
                                w3 = 0.2,
                                Cores = 1) {
  tictoc::tic()

  # === Validate Subbasins input ===
  if (!inherits(Subbasins, "SpatVector")) stop("Subbasins must be a 'SpatVector'.")
  if (!all(c("COMID", "Region") %in% names(Subbasins))) {
    stop("Subbasins must contain fields 'COMID' and 'Region'.")
  }
  comid  <- as.character(Subbasins$COMID)      # Subbasin IDs as character
  region <- Subbasins$Region                   # Region identifiers for each subbasin
  area   <- terra::expanse(Subbasins, unit = 'km')  # Area of each subbasin in km²
  nsub   <- length(comid)                      # Total number of subbasins

  # === Validate and reorder TransferMatrix ===
  MT <- as.matrix(TransferMatrix)
  if (any(dim(MT) != c(nsub, nsub))) {
    stop(sprintf("TransferMatrix must be a square %d x %d object.", nsub, nsub))
  }
  if (is.null(rownames(MT)) || is.null(colnames(MT))) {
    stop("TransferMatrix must have rownames and colnames corresponding to COMID.")
  }
  miss_rows <- setdiff(comid, rownames(MT))
  miss_cols <- setdiff(comid, colnames(MT))
  if (length(miss_rows) || length(miss_cols)) {
    stop(sprintf(
      "TransferMatrix names do not match COMID. Missing rows: [%s]; Missing cols: [%s]",
      paste(miss_rows, collapse = ", "), paste(miss_cols, collapse = ", ")
    ))
  }
  # Reorder TransferMatrix rows and columns to match COMID order
  MT <- MT[comid, comid, drop = FALSE]

  # === Validate Outlet ===
  idx_outlet <- match(Outlet, comid)
  if (is.na(idx_outlet)) stop(sprintf("Outlet COMID '%s' not found in Subbasins.", Outlet))

  # === Validate Dates column (must be DatesR, as required by airGR) ===
  if (!"DatesR" %in% names(Data)) stop("Data must include a 'DatesR' column.")
  if (inherits(Data$DatesR, "Date")) {
    Data$DatesR <- as.POSIXct(Data$DatesR, tz = "UTC")
  } else if (!inherits(Data$DatesR, c("POSIXct","POSIXt"))) {
    Data$DatesR <- as.POSIXct(
      paste0(Data$DatesR, " 00:00:00"),
      tz = "UTC",
      tryFormats = c("%Y-%m-%d","%Y/%m/%d","%d-%m-%Y","%d/%m/%Y","%Y-%m","%m/%Y")
    )
  }
  if (any(is.na(Data$DatesR))) stop("Unrecognized date format in DatesR.")

  # Identify simulation start and end indices (using mm/YYYY format)
  ind_start <- which(format(Data$DatesR, "%m/%Y") == RunIni)
  ind_end   <- which(format(Data$DatesR, "%m/%Y") == RunEnd)
  if (!length(ind_start) || !length(ind_end)) stop("RunIni or RunEnd not found (mm/YYYY).")
  if (max(ind_end) < min(ind_start)) stop("RunEnd precedes RunIni.")

  # Subset data to simulation period
  Database <- Data[min(ind_start):max(ind_end), , drop = FALSE]
  if (!("Q" %in% names(Database))) stop("Data must include observed flow column 'Q' for calibration.")
  Dates    <- as.Date(Database$DatesR)
  ntime    <- length(Dates)                       # Number of timesteps
  nDays    <- lubridate::days_in_month(Dates)     # Days per month for unit conversion

  # === Validate Parameters data frame ===
  if (!is.data.frame(Parameters)) stop("Parameters must be a data frame.")
  req_cols <- c("Region","X1","X2","fp","fe")
  if (!all(req_cols %in% names(Parameters))) {
    stop("Parameters must contain columns: Region, X1, X2, fp, fe.")
  }
  all_regions <- sort(unique(region))             # Regions present in subbasins
  regs_have   <- sort(unique(Parameters$Region))  # Regions provided in Parameters
  if (!identical(all_regions, regs_have)) {
    stop(sprintf("Region mismatch. Needed: %s | Provided: %s",
                 paste(all_regions, collapse=","), paste(regs_have, collapse=",")))
  }

  # Optional routing constant k (months) per region; default to 0 (instantaneous
  # transfer, i.e. previous behavior) if not supplied, for backward compatibility.
  if (!"k" %in% names(Parameters)) {
    Parameters$k <- 0
    message("Parameters$k not supplied; defaulting to k = 0 (instantaneous routing, previous behavior).")
  }

  # === Determine which regions are optimized ===
  if (!is.null(No.Optim) && any(!(No.Optim %in% all_regions))) {
    stop("No.Optim contains regions not present in Subbasins.")
  }
  opt_regions <- if (is.null(No.Optim)) all_regions else setdiff(all_regions, No.Optim)
  if (!length(opt_regions)) stop("No regions left to optimize (check No.Optim).")

  # === Validate forcing data columns ===
  p_names <- paste0("P_", comid)
  e_names <- paste0("E_", comid)
  miss_p  <- setdiff(p_names, names(Database))
  miss_e  <- setdiff(e_names, names(Database))
  if (length(miss_p) || length(miss_e)) {
    stop(sprintf("Missing columns in Data. P: [%s]; E: [%s]",
                 paste(miss_p, collapse=","), paste(miss_e, collapse=",")))
  }

  # === Validate bounds ===
  if (!(length(Parameters.Min) == 5 && length(Parameters.Max) == 5)) {
    stop("Parameters.Min/Max must be length 5 (X1, X2, fp, fe, k).")
  }

  # Initial parameter vector for optimization (only for opt_regions)
  opt.param     <- unlist(Parameters[Parameters$Region %in% opt_regions, c("X1","X2","fp","fe","k")], use.names = FALSE)
  opt.param.min <- rep(Parameters.Min, times = length(opt_regions))
  opt.param.max <- rep(Parameters.Max, times = length(opt_regions))

  # === Helper function to reconstruct parameter data frame from vector ===
  get_param <- function(param_vec, regions) {
    n <- length(regions)
    data.frame(Region = regions,
               X1 = param_vec[1:n],
               X2 = param_vec[(n + 1):(2 * n)],
               fp = param_vec[(2 * n + 1):(3 * n)],
               fe = param_vec[(3 * n + 1):(4 * n)],
               k  = param_vec[(4 * n + 1):(5 * n)],
               row.names = NULL)
  }

  # === Helper function to build forcing input for a subbasin ===
  forcing_input <- function(Param, RegionID, Database, comid_i) {
    fp <- Param$fp[Param$Region == RegionID]
    fe <- Param$fe[Param$Region == RegionID]
    data.frame(DatesR = as.POSIXct(Dates, tz = "UTC"),
               P = fp * Database[[paste0("P_", comid_i)]],
               E = fe * Database[[paste0("E_", comid_i)]])
  }

  # === Routing topology: built ONCE, outside OFUN ===
  # The network topology never changes during calibration (only the GR2M and
  # routing parameters do), so rebuilding the graph and recomputing the
  # topological order on every objective-function evaluation (potentially
  # Max.Functions times) would be pure waste -- especially for large networks.
  g <- igraph::graph_from_adjacency_matrix(t(MT), mode = "directed")
  if (!igraph::is_dag(g)) {
    stop("TransferMatrix induces a cyclic graph; routing requires a valid dendritic network (DAG).")
  }
  order_sub <- as.integer(igraph::topo_sort(g, mode = "out"))

  # === Parallel cluster: built ONCE, outside OFUN (same reasoning as the
  # routing topology above) -- creating/tearing down a cluster on every
  # objective-function evaluation would cost far more than it saves. Static
  # objects (that don't change across evaluations) are exported once here;
  # only Param (which does change every call) is passed to each parLapply()
  # call directly. No-op when Cores = 1 (plain lapply, no cluster).
  cl <- NULL
  if (Cores > 1) {
    cl <- parallel::makeCluster(Cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(airGR))
    parallel::clusterExport(cl, varlist = c("Database", "region", "area", "nDays",
                                            "comid", "ntime", "forcing_input", "Dates"),
                            envir = environment())
  }

  # === Per-subbasin GR2M run (independent until routing) ===
  run_one_subbasin_calib <- function(i, Param) {
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

    output <- RunModel(InputsModel = model_input,
                       RunOptions = run_opt,
                       Param = param_i,
                       FUN_MOD = RunModel_GR2M)

    # Convert runoff from mm to discharge in m³/s
    (area[i] * output$Qsim) / (86.4 * nDays)
  }

  # === Objective function with routing ===
  OFUN <- function(par) {
    # Build parameter table from current vector
    Param <- get_param(par, opt_regions)

    # Merge fixed parameters (if some regions are excluded from optimization)
    if (!is.null(No.Optim)) {
      fixed <- Parameters[Parameters$Region %in% No.Optim, ]
      Param <- rbind(Param, fixed[, c("Region","X1","X2","fp","fe","k")])
      Param <- Param[match(all_regions, Param$Region), ]
    } else {
      Param <- Param[match(all_regions, Param$Region), ]
    }

    # Run GR2M for each subbasin independently (no routing yet)
    if (Cores > 1) {
      qs_list <- parallel::parLapply(cl, seq_len(nsub), run_one_subbasin_calib, Param = Param)
    } else {
      qs_list <- lapply(seq_len(nsub), run_one_subbasin_calib, Param = Param)
    }
    QS <- do.call(cbind, qs_list)


    # === Routing: propagate flows downstream through a linear reservoir per reach ===
    # (g/order_sub are computed once outside OFUN -- topology doesn't change here.)
    k_vec <- Param$k[match(region, Param$Region)]

    # Initialize routed flows with local runoff
    QR <- QS
    colnames(QR) <- comid

    # Traverse network and propagate flows downstream, transiting each donor's
    # discharge through its own reach (linear reservoir, k in months) before
    # adding it to the receiver -- k = 0 reproduces the previous instantaneous sum.
    for (j in order_sub) {
      # Identify downstream receivers of subbasin j
      rec_ids <- which(MT[, j] != 0)

      if (length(rec_ids) > 0) {
        routed_j <- route_linear_reservoir(QR[, j], k = k_vec[j])
        for (i in rec_ids) {
          QR[, i] <- QR[, i] + routed_j
        }
      }
    }

    # Extract outlet discharge after warm-up period
    Qsim <- as.numeric(QR[-seq_len(WarmUp), idx_outlet])
    # Qsim <- rowSums(QS)[-seq_len(WarmUp)]
    Qobs <- Database$Q[-seq_len(WarmUp)]

    # === Compute performance metrics ===
    kge     <- 1 - KGE(Qsim, Qobs)
    nse     <- 1 - NSE(Qsim, Qobs)
    nselog  <- 1 - NSE(Qsim, Qobs, fun = log, epsilon.type = "Pushpalatha2012")
    rpear_v <- 1 - rPearson(Qsim, Qobs)
    kgekm_v <- 1 - KGEkm(Qsim, Qobs)
    kgelf_v <- 1 - KGElf(Qsim, Qobs)
    rmse_v  <- rmse(Qsim, Qobs)
    pbias_v <- abs(pbias(Qsim, Qobs))
    pbiasfdc_v <- abs(pbiasfdc(Qsim, Qobs, plot = FALSE))

    # RSR (RMSE-observations standard deviation ratio, Moriasi et al. 2007):
    # RMSE normalized by sd(Qobs), used ONLY inside the composite objectives
    # (OF8-OF10) so it is on the same roughly-[0, ~2] scale as the other
    # (1 - efficiency) terms. Using the raw RMSE (in m3/s, unbounded) there
    # would make it dominate the weighted sum regardless of w1/w2/w3 for any
    # basin with non-trivial discharge -- OF4 still reports the raw RMSE.
    rmse_norm <- rmse_v / sd(Qobs)

    # Return selected objective function value
    criteria <- c(
      OF1  = kge,
      OF2  = nse,
      OF3  = nselog,
      OF4  = rmse_v,
      OF5  = pbias_v,
      OF6  = pbiasfdc_v,
      OF7  = rpear_v,
      OF8  = (w1 * kge + w2 * nse + w3 * rmse_norm)/ (w1 + w2 + w3),
      OF9  = (w1 * kge + w2 * nselog + w3 * rmse_norm)/ (w1 + w2 + w3),
      OF10 = (w1 * kgekm_v + w2 * kgelf_v + w3 * rmse_norm)/ (w1 + w2 + w3)
    )
    if (!Optimization %in% names(criteria)) {
      stop(sprintf("Unknown 'Optimization' value '%s'. Must be one of: %s.",
                   Optimization, paste(names(criteria), collapse = ", ")))
    }
    criteria[Optimization]
  }

  # === Run optimization with SCE-UA algorithm ===
  message(sprintf("Optimizing %s with SCE-UA over %d regions (%d params)...",
                  Optimization, length(opt_regions), length(opt_regions) * 4))

  Calibration <- rtop::sceua(OFUN,
                             pars  = opt.param,
                             lower = opt.param.min,
                             upper = opt.param.max,
                             maxn  = Max.Functions)

  # === Organize output ===
  final_params <- get_param(Calibration$par, opt_regions)
  if (!is.null(No.Optim)) {
    fixed <- Parameters[Parameters$Region %in% No.Optim, c("Region","X1","X2","fp","fe","k")]
    final_params <- rbind(final_params, fixed)
    final_params <- final_params[match(all_regions, final_params$Region), ]
  } else {
    final_params <- final_params[match(all_regions, final_params$Region), ]
  }
  final_obj_raw <- Calibration$value

  message("Optimization complete.")
  cat("\nFinal calibrated parameters per region:\n")
  print(final_params)
  cat(sprintf("Objective Function (%s) = %s\n", Optimization, final_obj_raw))
  tictoc::toc()

  # Return
  list(Parameters = final_params,
       OF = final_obj_raw)
}
