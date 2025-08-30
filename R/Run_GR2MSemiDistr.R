#' Run the GR2M semi-distributed model with river network routing
#'
#' Executes the monthly GR2M water balance model across multiple subbasins in a semi-distributed configuration.
#' Region-specific parameters and correction factors are applied, and discharges are routed through the
#' subbasin connectivity defined in a transfer matrix. The function supports a warm-up period, outlet extraction,
#' saving results to disk, and an update mode for appending new results.
#'
#' @param Data data.frame. Model input data in the airGR format, as produced by `Create_Forcing_Inputs`.
#' Must include columns: DatesR, P_1 to P_n, E_1 to E_n, and optionally Q (observed discharge, used if available).
#' @param Subbasins SpatVector. Geometries of subbasins. Must include attributes "COMID" (unique subbasin ID) and "Region" (region name/code).
#' @param RunIni character. Simulation start date in the format "mm/yyyy".
#' @param RunEnd character. Simulation end date in the format "mm/yyyy".
#' @param Parameters data.frame. GR2M model parameters and correction factors per region. Must include columns: Region, X1, X2, fp, fe.
#' @param TransferMatrix matrix or dgCMatrix. Subbasin connectivity matrix defining downstream routing.
#' Rows represent receiving (downstream) subbasins, and columns represent donor (upstream) subbasins.
#' Row and column names must match `COMID`. The matrix is internally reordered to align with input forcings.
#' For large river networks (thousands of subbasins), a sparse matrix is strongly recommended for memory efficiency.
#' @param Outlet character (optional). COMID of the outlet subbasin to extract routed discharge.
#' If provided, the simulated outlet discharge (`qsim`) is returned, and compared with observed discharge (`Q`) if available in `Data`.
#' @param WarmUp integer (optional). Number of months to discard from the beginning of the simulation as model warm-up.
#' Ignored if `Update = TRUE`. Default is `NULL`.
#' @param IniState list (optional). Initial state variables for GR2M subbasins.
#' If `NULL`, default initial conditions are used.
#' @param Save logical. If `TRUE`, simulation results are saved in tab-separated `.txt` files inside the `./Outputs` directory.
#' Default is `FALSE`.
#' @param Update logical. If `TRUE`, the function runs only for the last month of the input data, appends results to existing
#' output files (if present), and discards the artificial extra month added internally to satisfy airGR's minimum timestep requirement.
#' Default is `FALSE`.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{Dates}{Date vector. Simulation period (monthly).}
#'   \item{COMID}{Character vector. Identifiers of the simulated subbasins.}
#'   \item{PR}{matrix. Precipitation inputs per subbasin (mm). Columns correspond to subbasins.}
#'   \item{AE}{matrix. Actual evapotranspiration per subbasin (mm).}
#'   \item{SM}{matrix. Soil moisture storage per subbasin (mm).}
#'   \item{PC}{matrix. Percolation per subbasin (mm).}
#'   \item{RU}{matrix. Runoff generated per subbasin (mm).}
#'   \item{QS}{matrix. Local discharge per subbasin (m³/s), obtained from runoff conversion by area.}
#'   \item{QR}{matrix. Routed discharge per subbasin (m³/s), after applying the connectivity matrix.}
#'   \item{SINK}{data.frame (optional). Simulated discharge at the outlet subbasin. Contains
#'   column `sim` (simulated) and, if available, column `obs` (observed). Returned only if `Outlet` is provided.}
#'   \item{EndState}{list. Final GR2M states for each subbasin, as returned by `airGR::RunModel_GR2M`.}
#' }
#'
#' @details
#' In `Update = TRUE` mode, only the last month of the input series is processed. Since airGR requires at least two
#' timesteps, the function appends an artificial month with P = 100 mm and E = 100 mm. This artificial record is
#' discarded after the run, ensuring only the actual last month is saved and returned.
#'
#' The transfer matrix describes subbasin connectivity: rows correspond to donor subbasins and columns to receivers.
#' It can be supplied as a `data.frame`, a base R `matrix`, or a sparse `dgCMatrix`. For applications with many
#' subbasins (e.g., >10,000), a sparse representation is strongly advised to reduce memory usage and improve speed.
#'
#' @references
#' Perrin C., Michel C., Andréassian V. (2003). Improvement of a parsimonious model for streamflow simulation.
#' \emph{Journal of Hydrology}, 279(1–4), 275–289. \doi{10.1016/S0022-1694(03)00225-7}
#'
#' Llauca H., Lavado-Casimiro W., Montesinos C., Santini W., Rau P. (2021).
#' PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020).
#' \emph{Water}, 13(8), 1048. \doi{10.3390/w13081048}
#'
#' @examples
#' library(GR2MSemiDistr)
#'
#' # Parameters (not calibrated)
#' parameters <- data.frame(
#'     Region = unique(cat$Region),
#'     X1  = 500,  # Production store capacity [mm]
#'     X2  = 1.5,  # Groundwater exchange coefficient
#'     fp  = 1.0,  # Precipitation correction factor
#'     fe  = 1.0   # Evapotranspiration correction factor
#'     )
#'
#' # Run model with given parameters
#' model <- Run_GR2MSemiDistr2(
#'   Data = model_inputs,
#'   Subbasins = cat,
#'   RunIni = "01/1981",
#'   RunEnd = "12/2016",
#'   Parameters = parameters,
#'   TransferMatrix = matrixT,
#'   Outlet = '8'
#' )
#'
#' # Shown results for a target COMID
#' target_comid <- "10"
#' idx <- which(model$COMID == target_comid)
#'
#' par(mfrow = c(3, 2),
#'     mar = c(4, 4, 3, 1),
#'     cex.main = 0.9,
#'     cex.lab = 0.8,
#'     cex.axis = 0.7)
#'
#' # Precipitation
#' plot(model$Dates, model$PR[, idx], type = "h", lwd = 2, col = "dodgerblue",
#'      main = sprintf("Precipitation (PR) - COMID %s", target_comid),
#'      xlab = "Date", ylab = "mm/month")
#'
#' # Actual evapotranspiration
#' plot(model$Dates, model$AE[, idx], type = "h", lwd = 2, col = "orange",
#'      main = sprintf("Actual Evapotranspiration (AE) - COMID %s", target_comid),
#'      xlab = "Date", ylab = "mm/month")
#'
#' # Soil moisture storage
#' plot(model$Dates, model$SM[, idx], type = "h", lwd = 2, col = "green",
#'      main = sprintf("Soil Moisture (SM) - COMID %s", target_comid),
#'      xlab = "Date", ylab = "mm")
#'
#' # Percolation
#' plot(model$Dates, model$PC[, idx], type = "h", lwd = 2, col = "purple",
#'      main = sprintf("Percolation (PC) - COMID %s", target_comid),
#'      xlab = "Date", ylab = "mm/month")
#'
#' # Runoff
#' plot(model$Dates, model$RU[, idx], type = "h", lwd = 2, col = "red",
#'      main = sprintf("Runoff (RU) - COMID %s", target_comid),
#'      xlab = "Date", ylab = "mm/month")
#'
#' # Routed discharge
#' plot(model$Dates, model$QR[, idx], type = "l", lwd = 1.5, col = "blue",
#'      main = sprintf("Routed Simulated Discharge (QR) - COMID %s", target_comid),
#'      xlab = "Date", ylab = "m³/s")
#'
#' # Compare simulated vs observed discharge at outlet
#' par(mfrow = c(1, 1))
#' plot(model$Dates, model$SINK$sim, type = "l", col = "blue", lwd = 1.5,
#'      xlab = "Date", ylab = "Discharge [m³/s]",
#'      main = "Simulated vs Observed Discharge (Outlet)")
#'   lines(model$Dates, model$SINK$obs, col = "red", lty = 2, lwd = 1.5)
#'   legend("topright", legend = c("Simulated", "Observed"),
#'          col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")
#'
#' @import airGR
#' @import terra
#' @import Matrix
#' @import tictoc
#' @import lubridate
#' @import igraph
#'
#' @export
Run_GR2MSemiDistr <- function(Data,
                              Subbasins,
                              RunIni,
                              RunEnd,
                              Parameters,
                              TransferMatrix,
                              Outlet = NULL,
                              WarmUp = NULL,
                              IniState = NULL,
                              Save = FALSE,
                              Update = FALSE) {

  tictoc::tic()

  # Validate that Subbasins is a SpatVector
  if (!inherits(Subbasins, "SpatVector")) {
    stop("'Subbasins' must be a 'SpatVector'.")
  }

  # Check that Subbasins contains the required attributes
  if (!all(c("COMID", "Region") %in% names(Subbasins))) {
    stop("'Subbasins' must contain attributes 'COMID' and 'Region'.")
  }

  # Extract basic attributes from Subbasins
  comid  <- as.character(Subbasins$COMID)         # Subbasin IDs as character
  region <- Subbasins$Region                      # Region assignment for each subbasin
  area   <- terra::expanse(Subbasins, unit = "km")# Subbasin areas in km²
  nsub   <- length(comid)                         # Number of subbasins

  # Check that input Data contains DatesR column as required by airGR
  if (!"DatesR" %in% names(Data)) {
    stop("Data must contain column 'DatesR' as required by airGR.")
  }

  # Convert DatesR to POSIXct if needed
  if (inherits(Data$DatesR, "Date")) {
    Data$DatesR <- as.POSIXct(Data$DatesR, tz = "UTC")
  } else if (!inherits(Data$DatesR, c("POSIXct", "POSIXt"))) {
    Data$DatesR <- as.POSIXct(
      paste0(Data$DatesR, " 00:00:00"),
      tz = "UTC",
      tryFormats = c("%Y-%m-%d","%Y/%m/%d","%d-%m-%Y","%d/%m/%Y","%Y-%m","%m/%Y")
    )
  }

  # Stop if any date could not be recognized
  if (anyNA(Data$DatesR)) stop("Unrecognized date format in DatesR.")

  # Identify start and end indices for the simulation window
  ind_start <- which(format(Data$DatesR, "%m/%Y") == RunIni)
  ind_end   <- which(format(Data$DatesR, "%m/%Y") == RunEnd)

  if (!length(ind_start) || !length(ind_end)) {
    stop("RunIni or RunEnd not found (expected 'mm/YYYY').")
  }
  if (max(ind_end) < min(ind_start)) {
    stop("RunEnd precedes RunIni.")
  }

  # Subset the database for the selected simulation period
  Database <- Data[min(ind_start):max(ind_end), , drop = FALSE]
  Dates    <- as.Date(Database$DatesR)

  # If Update mode is activated, keep only the last month
  if (Update) {
    Database <- tail(Database, 1)
    Dates    <- tail(Dates, 1)

    # If only one month is available, add an artificial next month
    # with P = 100 and E = 100 to satisfy airGR's minimum timestep requirement
    d1 <- as.Date(Dates[1])
    next_month <- seq(d1, by = "month", length.out = 2)[2]

    Database <- rbind(
      Database,
      cbind(DatesR = as.POSIXct(next_month, tz="UTC"),
            as.data.frame(matrix(100, nrow=1, ncol=ncol(Database)-1,
                                 dimnames=list(NULL, names(Database)[-1]))))
    )

    Dates <- as.Date(Database$DatesR)
  }

  # Number of time steps and days in each month
  ntime    <- length(Dates)
  nDays    <- lubridate::days_in_month(Dates)

  # Validate Parameters structure
  req_cols <- c("Region", "X1", "X2", "fp", "fe")
  if (!all(req_cols %in% names(Parameters))) {
    stop("Parameters must contain: Region, X1, X2, fp, fe.")
  }

  # Ensure that all required regions are covered by Parameters
  regs_needed <- sort(unique(region))
  regs_have   <- sort(unique(Parameters$Region))
  if (!identical(regs_needed, regs_have)) {
    stop("Mismatch between required and provided Regions.")
  }

  # Match parameters to subbasins by Region
  match_idx <- match(region, Parameters$Region)
  X1_vec <- Parameters$X1[match_idx]
  X2_vec <- Parameters$X2[match_idx]
  fp_vec <- Parameters$fp[match_idx]
  fe_vec <- Parameters$fe[match_idx]

  # Validate forcing data presence
  p_names <- paste0("P_", comid)
  e_names <- paste0("E_", comid)
  miss_p  <- setdiff(p_names, names(Database))
  miss_e  <- setdiff(e_names, names(Database))
  if (length(miss_p) || length(miss_e)) {
    stop("Missing forcing columns in Database.")
  }

  # Validate TransferMatrix
  if (is.null(TransferMatrix)) stop("TransferMatrix is required.")

  # Convert data.frame to matrix if needed
  if (is.data.frame(TransferMatrix)) TransferMatrix <- as.matrix(TransferMatrix)

  # Convert to sparse matrix for efficiency
  MT <- as(TransferMatrix, "dgCMatrix")

  # Validate dimensions and row/column names
  if (any(dim(MT) != c(nsub, nsub))) {
    stop(sprintf("TransferMatrix must be %d x %d.", nsub, nsub))
  }
  if (is.null(rownames(MT)) || is.null(colnames(MT))) {
    stop("TransferMatrix must have rownames and colnames corresponding to COMID.")
  }

  rows_mt <- as.character(rownames(MT))
  cols_mt <- as.character(colnames(MT))
  if (!setequal(comid, rows_mt) || !setequal(comid, cols_mt)) {
    stop("TransferMatrix COMID mismatch.")
  }

  # Reorder TransferMatrix to match COMID order
  MT <- MT[comid, comid, drop = FALSE]

  # Run GR2M model for each subbasin
  message(sprintf("Running GR2M for %d subbasins", nsub))
  ResModel <- vector("list", nsub)

  for (i in seq_len(nsub)) {
    # Build input dataframe with precipitation and evapotranspiration
    Input <- data.frame(
      DatesR = Database$DatesR,
      P      = fp_vec[i] * Database[[p_names[i]]],
      E      = fe_vec[i] * Database[[e_names[i]]]
    )

    # Stop if NA values are found in P or E
    if (anyNA(Input$P) || anyNA(Input$E)) {
      stop(sprintf("NA values in P/E for COMID %s", comid[i]))
    }

    # Create input structure for airGR
    InputsModel <- CreateInputsModel(FUN_MOD=RunModel_GR2M,
                                     DatesR=Input$DatesR,
                                     Precip=Input$P,
                                     PotEvap=Input$E)

    # Create run options
    RunOptions  <- CreateRunOptions(FUN_MOD=RunModel_GR2M,
                                    InputsModel=InputsModel,
                                    IndPeriod_Run=seq_len(ntime),
                                    verbose=FALSE,
                                    warnings=FALSE)

    # Run GR2M for the subbasin
    ResModel[[i]] <- RunModel(InputsModel=InputsModel,
                              RunOptions=RunOptions,
                              Param=c(X1_vec[i],X2_vec[i]),
                              FUN=RunModel_GR2M)
  }

  # Helper to extract matrices from ResModel
  mat_from_list <- function(key, round_digits=2) {
    out <- do.call(cbind, lapply(ResModel, function(r) round(r[[key]], round_digits)))
    colnames(out) <- comid
    out
  }

  # Extract outputs: effective precipitation, AE, soil moisture, percolation, runoff
  PR <- mat_from_list("Precip")
  AE <- mat_from_list("AE")
  SM <- mat_from_list("Prod")
  PC <- mat_from_list("Perc")
  RU <- mat_from_list("Qsim")

  # Convert runoff (mm) into discharge (m³/s) using area and number of days
  QS <- matrix(NA, nrow=ntime, ncol=nsub)
  for (i in seq_len(nsub)) {
    QS[,i] <- round((area[i] * ResModel[[i]]$Qsim) / (86.4 * nDays), 2)
  }
  colnames(QS) <- comid

  # === Routing: propagate accumulated flows downstream ===
  message("Performing routing with TransferMatrix...")

  # Build graph with correct orientation:
  # MT[i, j] = 1 means subbasin j drains into i
  g <- igraph::graph_from_adjacency_matrix(t(MT), mode = "directed")

  # Compute topological order (headwaters -> outlet)
  order_sub <- as.integer(igraph::topo_sort(g, mode = "out"))

  # Initialize routed flows with local runoff
  QR <- QS
  colnames(QR) <- comid

  # Traverse network and propagate accumulated flows downstream
  for (j in order_sub) {
    # Identify downstream receivers of subbasin j
    rec_ids <- which(MT[, j] != 0)

    # Pass accumulated discharge from j to each downstream receiver
    if (length(rec_ids) > 0) {
      for (i in rec_ids) {
        QR[, i] <- QR[, i] + QR[, j]
      }
    }
  }

  # If Update mode was used, discard the artificial extra month
  if (Update) {
    PR <- tail(PR, 1)
    AE <- tail(AE, 1)
    SM <- tail(SM, 1)
    PC <- tail(PC, 1)
    RU <- tail(RU, 1)
    QS <- tail(QS, 1)
    QR <- tail(QR, 1)
    Dates <- tail(Dates, 1)
  }

  # Extract outlet discharge if Outlet is provided
  qsim <- NULL
  qobs <- NULL
  if (!is.null(Outlet)) {
    idx_outlet <- match(Outlet, comid)
    if (is.na(idx_outlet)) {
      stop(sprintf("Outlet COMID '%s' not found in Subbasins.", Outlet))
    }
    qsim <- round(QR[, idx_outlet], 2)
    if ("Q" %in% names(Database)) {
      qobs <- Database$Q
    }
  }

  # Save results to disk as tab-separated text files
  if (Save) {
    if (!dir.exists("./Outputs")) {
      dir.create("./Outputs", recursive=TRUE)
    }

    save_one <- function(var, tag) {
      df <- as.data.frame(var)
      rownames(df) <- format(Dates, "%Y-%m-01")
      file <- sprintf("./Outputs/%s_GR2MSemiDistr.txt", tag)

      # In Update mode, append to old file instead of overwriting
      if (Update && file.exists(file)) {
        old_df <- read.table(file, header=TRUE, sep="\t", check.names=FALSE)
        new_df <- rbind(old_df, df)
        write.table(new_df, file=file, sep="\t", quote=FALSE)
      } else {
        write.table(df, file=file, sep="\t", quote=FALSE)
      }
    }

    # Save all variables separately
    save_one(PR,"PR")
    save_one(AE,"AE")
    save_one(SM,"SM")
    save_one(PC,"PC")
    save_one(RU,"RU")
    save_one(QS,"QS")
    save_one(QR,"QR")
  }

  # Report execution time
  message("Processing completed in...")
  tictoc::toc()

  # Assemble results to return
  Ans <- list(
    Dates    = Dates,
    COMID    = comid,
    PR       = PR,
    AE       = AE,
    SM       = SM,
    PC       = PC,
    RU       = RU,
    QS       = QS,
    QR       = QR,
    EndState = lapply(ResModel, `[[`, "StateEnd")
  )

  # Add outlet discharge if available
  if (!is.null(Outlet)) {
    if (!is.null(qobs)) {
      Ans$SINK <- data.frame(sim = qsim, obs = round(qobs,2), row.names = Dates)
    } else {
      Ans$SINK <- data.frame(sim = qsim, row.names = Dates)
    }
  }

  return(Ans)
}
