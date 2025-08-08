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
#' # Load preprocessed input data
#' data(cat)     # Subbasins (SpatVector with COMID and Region)
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
#'   Data = inputs,
#'   Subbasins = cat,
#'   RunIni = "01/1981",
#'   RunEnd = "12/2005",
#'   Parameters = param_init,
#'   Parameters.Min = c(100, 0.1, 0.8, 0.8),
#'   Parameters.Max = c(2000, 10, 1.2, 1.2),
#'   WarmUp = 12,
#'   Max.Functions = 1000,
#'   Optimization = "OF1"  # KGE
#' )
#'
#' # Extract results
#' best_params <- result$Parameters
#' final_score <- result$OF
#'
#' # Plot observed vs simulated discharge at the outlet
#' model <- Run_GR2MSemiDistr(
#'   Data = inputs,
#'   Subbasins = cat,
#'   RunIni = "01/1981",
#'   RunEnd = "12/2005",
#'   Parameters = best_params,
#'   WarmUp = 12
#' )
#'
#' if (!is.null(model$SINK)) {
#'   dates <- as.Date(model$Dates)
#'   plot(dates, model$SINK$sim, type = "l", col = "blue", lwd = 2,
#'        xlab = "Date", ylab = "Discharge [m³/s]",
#'        main = "Simulated vs Observed Discharge (Outlet)")
#'   lines(dates, model$SINK$obs, col = "darkgreen", lwd = 2, lty = 2)
#'   legend("topright", legend = c("Simulated", "Observed"),
#'          col = c("blue", "darkgreen"), lty = c(1, 2), lwd = 2, bty = "n")
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

  # === Validate inputs ===
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
  Data$DatesR <- as.POSIXct(paste0(Data$DatesR, " 00:00:00"), tz = "GMT",
                            tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%d-%m-%Y", "%d/%m/%Y"))

  ind_start <- which(format(Data$DatesR, "%m/%Y") == RunIni)
  ind_end   <- which(format(Data$DatesR, "%m/%Y") == RunEnd)
  if (length(ind_start) == 0 || length(ind_end) == 0) {
    stop("RunIni or RunEnd not found in 'Data$DatesR'. Use format 'mm/yyyy'.")
  }

  Database <- Data[ind_start:ind_end, ]
  Dates    <- as.Date(Database$DatesR)
  ntime    <- length(Dates)
  nDays    <- lubridate::days_in_month(Dates)

  # === Define calibration regions ===
  all_regions <- sort(unique(region))
  opt_regions <- if (is.null(No.Optim)) all_regions else setdiff(all_regions, No.Optim)
  if (!all(opt_regions %in% Parameters$Region)) stop("Mismatch between optimized regions and parameters.")

  opt.param     <- unlist(Parameters[Parameters$Region %in% opt_regions, c("X1", "X2", "fp", "fe")])
  opt.param.min <- rep(Parameters.Min, each = length(opt_regions))
  opt.param.max <- rep(Parameters.Max, each = length(opt_regions))

  # === Helper functions ===
  get_param <- function(param_vec, regions) {
    n <- length(regions)
    data.frame(Region = regions,
               X1  = param_vec[1:n],
               X2  = param_vec[(n + 1):(2 * n)],
               fp  = param_vec[(2 * n + 1):(3 * n)],
               fe  = param_vec[(3 * n + 1):(4 * n)])
  }

  forcing_input <- function(Param, Region, Database, i, nsub) {
    fp <- Param$fp[Param$Region == Region]
    fe <- Param$fe[Param$Region == Region]
    data.frame(DatesR = Database[, 1],
               P = round(fp * Database[, i + 1], 1),
               E = round(fe * Database[, i + 1 + nsub], 1))
  }

  # === Objective function ===
  OFUN <- function(par) {
    Param <- get_param(par, opt_regions)
    if (!is.null(No.Optim)) {
      fixed <- Parameters[Parameters$Region %in% No.Optim, ]
      Param <- rbind(Param, fixed[match(setdiff(all_regions, opt_regions), fixed$Region), ])
      Param <- Param[match(all_regions, Param$Region), ]
    }

    # Run model
    QS <- matrix(NA, nrow = ntime, ncol = nsub)
    for (i in seq_len(nsub)) {
      reg_i   <- region[i]
      param_i <- c(Param$X1[Param$Region == reg_i], Param$X2[Param$Region == reg_i])
      input_i <- forcing_input(Param, reg_i, Database, i, nsub)

      model_input <- CreateInputsModel(RunModel_GR2M,
                                       DatesR = input_i$DatesR,
                                       Precip = input_i$P,
                                       PotEvap = input_i$E)

      run_opt <- CreateRunOptions(RunModel_GR2M,
                                  InputsModel = model_input,
                                  IndPeriod_Run = seq_len(ntime),
                                  verbose = FALSE,
                                  warnings = FALSE)

      output <- RunModel(model_input, run_opt, param_i, RunModel_GR2M)
      QS[, i] <- (area[i] * output$Qsim) / (86.4 * nDays)
    }

    Qsim <- rowSums(QS, na.rm = TRUE)
    Qobs <- Database$Q
    if (!is.null(WarmUp)) {
      Qsim <- Qsim[-seq_len(WarmUp)]
      Qobs <- Qobs[-seq_len(WarmUp)]
    }

    df <- na.omit(data.frame(y = Qsim, x = Qobs))

    # Calculate all metrics once
    kge   <- 1 - KGE(df$y, df$x)
    nse   <- 1 - NSE(df$y, df$x)
    nselog<- 1 - NSE(df$y, df$x, fun = log, epsilon.type = "Pushpalatha2012")
    rmse  <- rmse(df$y, df$x)
    pbias <- abs(pbias(df$y, df$x))
    pbiasfdc <- abs(pbiasfdc(df$y, df$x, plot = FALSE))
    rpear <- 1 - rPearson(df$y, df$x)
    kgekm <- 1 - KGEkm(df$y, df$x)
    kgelf <- 1 - KGElf(df$y, df$x)

    criteria <- c(
      OF1  = kge,
      OF2  = nse,
      OF3  = nselog,
      OF4  = rmse,
      OF5  = pbias,
      OF6  = pbiasfdc,
      OF7  = rpear,
      OF8  = w1 * kge + w2 * nse + w3 * rmse,
      OF9  = w1 * kge + w2 * nselog + w3 * rmse,
      OF10 = w1 * kgekm + w2 * kgelf + w3 * rmse
    )
    return(criteria[Optimization])
  }

  # === Run optimization ===
  message(paste("Optimizing", Optimization, "using SCE-UA..."))
  Calibration <- sceua(OFUN,
                       pars = opt.param,
                       lower = opt.param.min,
                       upper = opt.param.max,
                       maxn = Max.Functions)

  # === Organize output ===
  final_params <- if (is.null(No.Optim)) {
    get_param(Calibration$par, opt_regions)
  } else {
    fixed <- Parameters[Parameters$Region %in% No.Optim, ]
    all_params <- get_param(Calibration$par, opt_regions)
    full <- rbind(all_params, fixed)
    full[match(all_regions, full$Region), ]
  }

  final_obj <- if (Optimization %in% c("OF4", "OF5", "OF6")) {
    round(Calibration$value, 3)
  } else {
    round(1 - Calibration$value, 3)
  }

  message("Optimization complete.")
  cat("\nFinal calibrated parameters per region:\n")
  print(final_params)
  cat(paste0("Objective Function (", Optimization, ") = ", final_obj, "\n"))

  tictoc::toc()

  return(list(Parameters = final_params, OF = final_obj))
}
