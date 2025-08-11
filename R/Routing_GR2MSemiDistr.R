#' Routing of simulated discharges for each subbasin using a weighted flow accumulation approach.
#'
#' This function routes monthly simulated streamflows from a semi-distributed GR2M model using a flow accumulation approach based on D8 flow directions.
#' The routing network is derived from a digital elevation model (DEM) via WhiteboxTools, or optionally from a precomputed transfer matrix.
#' Each subbasin contributes its outflow as a weight that is accumulated along the flow direction grid to downstream subbasins.
#'
#' @param Model List containing model results, typically from 'Run_GR2MSemiDistr'. Must include:
#' \describe{
#'   \item{QS}{Matrix of simulated streamflows [m³/s], one column per subbasin.}
#'   \item{Dates}{Vector of dates corresponding to the rows of QS.}
#' }
#' @param Subbasins A 'SpatVector' object with the subbasin boundaries. Must include attributes:
#' \describe{
#'   \item{COMID}{Integer ID of each subbasin.}
#'   \item{Region}{Region code or name.}
#' }
#' @param Dem A 'SpatRaster' representing the digital elevation model. Required unless 'TransferMatrixFile' is used.
#' @param RouIni Character. Start date for routing in "mm/YYYY" format. If NULL, uses the start of the simulation.
#' @param RouEnd Character. End date for routing in "mm/YYYY" format. If NULL, uses the end of the simulation.
#' @param TransferMatrixFile Optional file path (.rds) to a precomputed transfer matrix. If provided, 'Dem' and WhiteboxTools are not used.
#' @param res Numeric. Desired spatial resolution in meters for resampling the DEM (default is 500). Converted to degrees if needed.
#' @param Save Logical. If TRUE, writes the routed discharges (QR) to a tabular .txt file in the './Outputs' directory. Default is FALSE.
#' @param Update Logical. If TRUE, deletes previous month's output and replaces it with the latest. Useful for incremental runs. Default is FALSE.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{QR}{Matrix of routed streamflows [m³/s] per subbasin and time step.}
#'   \item{Dates}{Vector of routing time steps.}
#'   \item{COMID}{Vector of subbasin COMIDs (column names of QR).}
#'   \item{TransferMatrix}{Sparse matrix used for flow routing (rows = outlets, columns = sources).}
#' }
#'
#' @references Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. (2021).
#' PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020). Water, 13(8), 1048. \doi{10.3390/w13081048}
#'
#' @examples
#' library(GR2MSemiDistr)
#'
#' # Run routing using DEM and WhiteboxTools
#' routed <- Routing_GR2MSemiDistr(
#'   Model = model,
#'   Subbasins = cat,
#'   Dem = dem,
#'   RouIni = "01/1981",
#'   RouEnd = "12/2016",
#'  )
#'
#' # Select COMID of the subbasin to plot
#' idx <- which(model$COMID == model$COMID[1]) # Change index as needed
#' dates <- as.Date(model$Dates)
#'
#' plot(dates, routed$QR[, idx], type = "l", col = "darkblue", lwd = 2,
#'      xlab = "Date", ylab = "Discharge [m³/s]",
#'      main = paste("Routed Discharge - Subbasin", routed$COMID[idx]))
#' @import terra
#' @import whitebox
#' @import sf
#' @import Matrix
#' @import lubridate
#' @export
#'
Routing_GR2MSemiDistr <- function(Model,
                                  Subbasins,
                                  Dem,
                                  RouIni = NULL,
                                  RouEnd = NULL,
                                  TransferMatrixFile = NULL,
                                  res = 500,
                                  Save = FALSE,
                                  Update = FALSE) {

  tictoc::tic()

  # --- helpers ---------------------------------------------------------------
  # Build accumulation matrix A = I + T + T^2 + ... (stops within N-1 for a DAG)
  .build_accum_matrix <- function(Tm) {
    n <- nrow(Tm)
    A <- Matrix::Diagonal(n)   # I
    P <- Tm
    it <- 1
    repeat {
      if (Matrix::nnzero(P) == 0) break
      A <- A + P
      P <- P %*% Tm
      it <- it + 1
      if (it >= n) break
    }
    A
  }

  # --- ensure WhiteboxTools ---------------------------------------------------
  if (!file.exists(whitebox::wbt_exe_path())) {
    message("Installing WhiteboxTools..."); whitebox::install_whitebox()
  }
  whitebox::wbt_init()

  # --- validate inputs --------------------------------------------------------
  if (!inherits(Subbasins, "SpatVector")) stop("'Subbasins' must be a SpatVector.")
  if (!all(c("COMID","Region") %in% names(Subbasins)))
    stop("'Subbasins' must contain fields 'COMID' and 'Region'.")
  if (!is.null(Dem)) {
    if (!inherits(Dem, "SpatRaster")) stop("'Dem' must be a SpatRaster.")
    if (!terra::same.crs(Dem, Subbasins)) stop("CRS of 'Dem' and 'Subbasins' must match.")
  }
  if (!is.null(RouIni) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouIni)) stop("RouIni must be 'mm/YYYY'.")
  if (!is.null(RouEnd) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouEnd))   stop("RouEnd must be 'mm/YYYY'.")
  if (!is.null(TransferMatrixFile) && !grepl("\\.rds$", TransferMatrixFile)) stop("TransferMatrixFile must be .rds")

  # --- IDs & QS ---------------------------------------------------------------
  comid <- as.character(Subbasins$COMID)
  nsub  <- length(comid)

  # Subset QS by dates (if requested)
  if ((is.null(RouIni) && is.null(RouEnd)) || nrow(Model$QS) == 1) {
    Dates <- as.Date(Model$Dates)
    QS    <- as.matrix(Model$QS)
  } else {
    all_dates <- as.Date(Model$Dates)
    d0 <- as.Date(paste0("01/", RouIni), "%d/%m/%Y")
    d1 <- as.Date(paste0("01/", RouEnd),  "%d/%m/%Y")
    i0 <- match(d0, all_dates); i1 <- match(d1, all_dates)
    if (any(is.na(c(i0, i1)))) stop("RouIni/RouEnd not found in Model$Dates.")
    Dates <- seq(d0, d1, by = "month")
    QS    <- as.matrix(Model$QS[i0:i1, , drop = FALSE])
  }
  # Align QS columns to COMID order
  if (!is.null(colnames(QS))) {
    miss <- setdiff(comid, colnames(QS))
    if (length(miss)) stop("QS is missing columns for COMID(s): ", paste(miss, collapse=", "))
    QS <- QS[, comid, drop = FALSE]
  } else colnames(QS) <- comid

  # --- build T with a single downstream per subbasin --------------------------
  # If provided, reuse; otherwise derive from DEM + Whitebox D8 + FAC
  if (!is.null(TransferMatrixFile) && file.exists(TransferMatrixFile)) {
    Tm <- readRDS(TransferMatrixFile)
    if (!inherits(Tm, "dgCMatrix")) stop("TransferMatrix must be 'dgCMatrix'.")
    if (!identical(rownames(Tm), comid) || !identical(colnames(Tm), comid))
      stop("TransferMatrix dimnames must match Subbasins COMID order.")
  } else {
    if (is.null(Dem)) stop("Provide 'Dem' when 'TransferMatrixFile' is not provided.")

    # (1) DEM resampling
    if (terra::is.lonlat(Dem)) {
      message("DEM is lon/lat. Converting target 'res' from meters to degrees for resampling...")
      lat_mean <- mean(c(terra::ymin(Dem), terra::ymax(Dem)))
      res_deg  <- res / (111320 * cos(lat_mean * pi / 180))
      target   <- terra::rast(terra::ext(Dem), resolution = res_deg, crs = terra::crs(Dem))
      Dem      <- terra::resample(Dem, target, method = "bilinear")
    } else {
      cur <- terra::res(Dem)
      if (any(abs(cur - c(res,res)) > 1e-6)) {
        target <- terra::rast(terra::ext(Dem), resolution = res, crs = terra::crs(Dem))
        Dem    <- terra::resample(Dem, target, method = "bilinear")
      }
    }

    # (2) Rasterize COMID over DEM grid
    comid_r <- terra::rasterize(Subbasins, Dem, field = "COMID", touches = TRUE)

    # (3) Whitebox preprocessing: Fill → D8 pointer → D8 FAC (cells)
    dem_tmp  <- tempfile(fileext = ".tif")
    dem_fill <- tempfile(fileext = ".tif")
    fdir_tif <- tempfile(fileext = ".tif")
    fac_tif  <- tempfile(fileext = ".tif")
    on.exit(unlink(c(dem_tmp, dem_fill, fdir_tif, fac_tif), force = TRUE), add = TRUE)

    terra::writeRaster(Dem, dem_tmp, overwrite = TRUE)
    whitebox::wbt_fill_depressions(dem_tmp, dem_fill)
    whitebox::wbt_d8_pointer(dem_fill, fdir_tif)
    whitebox::wbt_d8_flow_accumulation(dem_fill, fac_tif, out_type = "cells")

    fdir <- terra::rast(fdir_tif)
    fac  <- terra::rast(fac_tif)

    # (4) D8 code → (dr, dc) offset
    d8_map <- list(`1`=c(0,1), `2`=c(1,1), `4`=c(1,0), `8`=c(1,-1),
                   `16`=c(0,-1), `32`=c(-1,-1), `64`=c(-1,0), `128`=c(-1,1))
    nrows <- terra::nrow(fdir); ncols <- terra::ncol(fdir)

    # (5) Build Tm: one downstream per column (0/1)
    Tm <- Matrix::sparseMatrix(i = integer(0), j = integer(0),
                               dims = c(nsub, nsub),
                               dimnames = list(comid, comid))

    # We will choose, for each COMID j, the boundary cell whose D8 points outside j
    # and has the maximum FAC; the neighbor COMID is the unique downstream i.
    vals <- terra::values(comid_r, mat = FALSE)

    for (k in seq_len(nsub)) {
      id  <- comid[k]
      idx <- which(vals == id)
      if (!length(idx)) next

      best_fac  <- -Inf
      best_dest <- NA_character_

      for (cell in idx) {
        code <- fdir[cell]; if (is.na(code)) next
        key  <- as.character(code); if (!(key %in% names(d8_map))) next
        rc   <- terra::rowColFromCell(fdir, cell)
        off  <- d8_map[[key]]
        r2   <- rc[1] + off[1]; c2 <- rc[2] + off[2]
        if (r2 < 1 || r2 > nrows || c2 < 1 || c2 > ncols) next
        cell2 <- terra::cellFromRowCol(fdir, r2, c2)
        dest  <- comid_r[cell2]

        # Candidate: D8 leaves j and enters another COMID within the domain
        if (is.na(dest) || dest == id) next

        fac_here <- fac[cell]
        if (!is.na(fac_here) && fac_here > best_fac) {
          best_fac  <- fac_here
          best_dest <- as.character(dest)
        }
      }

      # Assign the unique downstream if found in the domain; otherwise j is an outlet
      if (!is.na(best_dest) && (best_dest %in% comid)) {
        Tm[best_dest, id] <- 1
      }
    }

    if (!is.null(TransferMatrixFile)) saveRDS(Tm, TransferMatrixFile)
  }

  # --- Build A and accumulate flows: QR = QS %*% A -----------------------------
  AccumMatrix <- .build_accum_matrix(Tm)  # reachability (includes self)
  QR <- QS %*% AccumMatrix                # accumulated flow at each subbasin

  # --- Optional save -----------------------------------------------------------
  if (Save) {
    dir.create("./Outputs", showWarnings = FALSE, recursive = TRUE)

    yyyymm_new <- format(tail(Dates, 1), "%Y%m")
    if (Update) {
      prev_month <- lubridate::floor_date(tail(Dates, 1), "month") %m-% lubridate::months(1)
      old_file <- sprintf("./Outputs/QR_GR2MSemiDistr_%s.txt", format(prev_month, "%Y%m"))
      if (file.exists(old_file)) file.remove(old_file)
    }

    df <- as.data.frame(round(QR, 1))
    colnames(df) <- paste0("QR_", comid)
    rownames(df) <- format(Dates, "%Y-%m-01")
    out_file <- sprintf("./Outputs/QR_GR2MSemiDistr_%s.txt", yyyymm_new)
    write.table(df, file = out_file, sep = "\t", quote = FALSE)
  }

  message("Processing completed successfully in...")
  tictoc::toc()
  list(QR = QR, Dates = Dates, COMID = comid,
       TransferMatrix = Tm,
       AccumMatrix = AccumMatrix)
}

