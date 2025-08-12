#' Routing of simulated discharges for each subbasin using a weighted flow-accumulation network
#'
#' Routes monthly simulated streamflows from a semi-distributed GR2M model through a
#' single-downstream (D8) river network. The routing network is derived from a Digital
#' Elevation Model (DEM) with WhiteboxTools, or provided directly via a precomputed
#' transfer matrix. Each subbasin's outflow is accumulated downstream using the
#' reachability (accumulation) matrix built from the transfer matrix.
#'
#' @details
#' The transfer matrix T encodes single-downstream connectivity between subbasins:
#' T[i, j] = 1 if subbasin j drains directly to i, and 0 otherwise.
#' The accumulation matrix is A = I + T + T^2 + ... (truncated within N-1
#' steps for a DAG), so that routed flows are computed as QR = QS %*% A, where
#' QS are simulated outflows per subbasin (columns) and time steps (rows).
#'
#' When deriving T from a DEM, D8 flow directions and flow accumulation (FAC) are
#' computed with WhiteboxTools. For each subbasin, the unique downstream neighbor is
#' selected from boundary cells by the candidate with the highest FAC. Subbasin
#' rasterization uses touches = FALSE to avoid spurious edge assignments.
#'
#' @param Model List with GR2M results, typically from Run_GR2MSemiDistr. Must include:
#' \describe{
#'   \item{QS}{Numeric matrix of simulated streamflows [m³/s]; columns are subbasins (COMID).}
#'   \item{Dates}{Vector of dates corresponding to the rows of QS.}
#' }
#' @param Subbasins A SpatVector of subbasin boundaries. Must contain attributes:
#' \describe{
#'   \item{COMID}{Integer ID of each subbasin (unique).}
#'   \item{Region}{Region code or name.}
#' }
#' @param Dem A SpatRaster DEM used to derive the routing network (required unless
#' TransferMatrix is supplied). If the DEM is in lon/lat, res is converted to degrees.
#' @param RouIni Character. Start date for routing in "mm/YYYY" format. If NULL,
#' uses the beginning of Model$Dates.
#' @param RouEnd Character. End date for routing in "mm/YYYY" format. If NULL,
#' uses the end of Model$Dates.
#' @param TransferMatrix Optional numeric square matrix (n x n) encoding direct
#' downstream connectivity (rows = destinations, columns = sources). If provided, Dem
#' and WhiteboxTools are not used. Row/column names (COMID) must match Subbasins$COMID
#' order; if missing, they are set from COMID.
#' @param res Numeric. Target DEM resolution in meters for resampling (default 500).
#' Converted to degrees when Dem is in lon/lat.
#' @param Save Logical. If TRUE, writes routed discharges (QR) as a tabular
#' .txt file into "./Outputs". Default FALSE.
#' @param Update Logical. If TRUE, removes the previous month's output before saving
#' the new one (useful for incremental/operational runs). Default FALSE.
#'
#' @return A list with:
#' \describe{
#'   \item{QR}{Base matrix of routed streamflows [m³/s]; rows = time, columns = COMID.}
#'   \item{Dates}{Vector of routing dates used for QR.}
#'   \item{COMID}{Character vector of subbasin IDs (column order of QR).}
#'   \item{TransferMatrix}{Base matrix used for routing; T[i,j] = 1 if j → i.}
#'   \item{AccumMatrix}{Base matrix A used to accumulate upstream contributions.}
#' }
#'
#' @seealso Run_GR2MSemiDistr
#'
#' @references
#' Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. (2021).
#' PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020).
#' Water, 13(8), 1048. \doi{10.3390/w13081048}
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
                                  TransferMatrix = NULL,  # numeric matrix or Matrix, not a file path
                                  res = 500,
                                  Save = FALSE,
                                  Update = FALSE) {

  tictoc::tic()

  # --- helper -----------------------------------------------------------------
  # Build A = I + T + T^2 + ... (stops within N-1 for a DAG)
  .build_accum_matrix <- function(Tm) {
    n <- nrow(Tm)
    A <- Matrix::Diagonal(n)
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

  if (is.null(TransferMatrix)) {
    if (is.null(Dem)) stop("Provide 'Dem' when 'TransferMatrix' is not provided.")
    if (!inherits(Dem, "SpatRaster")) stop("'Dem' must be a SpatRaster.")
    if (!terra::same.crs(Dem, Subbasins)) stop("CRS of 'Dem' and 'Subbasins' must match.")
  } else {
    if (!is.matrix(TransferMatrix) && !inherits(TransferMatrix, "Matrix"))
      stop("'TransferMatrix' must be a numeric matrix.")
    if (!is.numeric(TransferMatrix)) stop("'TransferMatrix' must be numeric.")
    if (nrow(TransferMatrix) != ncol(TransferMatrix))
      stop("'TransferMatrix' must be square (n x n).")
  }
  if (!is.null(RouIni) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouIni)) stop("RouIni must be 'mm/YYYY'.")
  if (!is.null(RouEnd) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouEnd))   stop("RouEnd must be 'mm/YYYY'.")

  # --- IDs & QS ---------------------------------------------------------------
  comid <- as.character(Subbasins$COMID)
  nsub  <- length(comid)

  # Subset QS by date range (if requested)
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
  storage.mode(QS) <- "double"

  # Align QS columns to COMID order
  if (!is.null(colnames(QS))) {
    miss <- setdiff(comid, colnames(QS))
    if (length(miss)) stop("QS is missing columns for COMID(s): ", paste(miss, collapse=", "))
    QS <- QS[, comid, drop = FALSE]
  } else {
    colnames(QS) <- comid
  }

  # --- transfer matrix (single downstream per subbasin) -----------------------
  if (!is.null(TransferMatrix)) {
    if (!is.null(rownames(TransferMatrix)) && !is.null(colnames(TransferMatrix))) {
      if (!identical(rownames(TransferMatrix), comid) || !identical(colnames(TransferMatrix), comid))
        stop("Dimnames of 'TransferMatrix' must match Subbasins COMID order.")
    } else {
      dimnames(TransferMatrix) <- list(comid, comid)
    }
    Tm_sp <- if (inherits(TransferMatrix, "Matrix")) {
      methods::as(TransferMatrix, "dgCMatrix")
    } else {
      Matrix::Matrix(TransferMatrix, sparse = TRUE)
    }
  } else {
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

    # (2) Rasterize COMID over DEM grid (touches = FALSE for stricter assignment)
    comid_r <- terra::rasterize(Subbasins, Dem, field = "COMID", touches = FALSE)

    # (3) Whitebox preprocessing: Fill -> D8 pointer -> D8 FAC (cells)
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

    # (4) D8 code -> (dr, dc) offset (vectorized lookup)
    codes <- c(1,2,4,8,16,32,64,128)
    dr    <- c(0, 1, 1, 1,  0, -1, -1, -1)
    dc    <- c(1, 1, 0,-1, -1, -1,  0,  1)
    names(dr) <- names(dc) <- as.character(codes)

    nrows <- terra::nrow(fdir); ncols <- terra::ncol(fdir)

    # (5) Pre-extract needed vectors once
    vals_all <- terra::values(comid_r, mat = FALSE)  # COMID per cell (may be NA)
    fac_all  <- terra::values(fac,     mat = FALSE)
    fdir_all <- terra::values(fdir,    mat = FALSE)

    # Build a lookup of cell indices by COMID (single pass)
    idx_valid   <- which(!is.na(vals_all))
    split_by_id <- split(idx_valid, vals_all[idx_valid])  # named by COMID value

    # (6) Prepare triplets to assemble sparse Tm at once
    i_vec <- integer(0)  # row indexes (destination)
    j_vec <- integer(0)  # col indexes (origin)

    # Map COMID to column/row positions once
    comid_pos <- seq_along(comid); names(comid_pos) <- comid

    for (id in comid) {
      cells <- split_by_id[[id]]
      if (is.null(cells) || length(cells) == 0) next

      codes_here <- fdir_all[cells]
      ok         <- !is.na(codes_here) & (codes_here %in% codes)
      if (!any(ok)) next

      cells_ok <- cells[ok]
      codes_ok <- as.character(codes_here[ok])

      # Vectorized (r,c) and neighbor (r2,c2)
      rc  <- terra::rowColFromCell(fdir, cells_ok)
      r2  <- rc[,1] + dr[codes_ok]
      c2  <- rc[,2] + dc[codes_ok]
      inb <- (r2 >= 1 & r2 <= nrows & c2 >= 1 & c2 <= ncols)
      if (!any(inb)) next

      cells_nbr <- terra::cellFromRowCol(fdir, r2[inb], c2[inb])
      dest_id   <- vals_all[cells_nbr]

      # Only transitions that leave this subbasin and land in another valid subbasin
      trans_ok <- !is.na(dest_id) & dest_id != id & (dest_id %in% comid)
      if (!any(trans_ok)) next

      cand_src <- cells_ok[inb][trans_ok]
      cand_dst <- as.character(dest_id[trans_ok])
      cand_fac <- fac_all[cand_src]

      # Pick the single downstream with max FAC among boundary candidates
      if (length(cand_fac)) {
        k_best <- which.max(cand_fac)
        j <- comid_pos[[id]]
        i <- comid_pos[[ cand_dst[k_best] ]]
        i_vec <- c(i_vec, i)
        j_vec <- c(j_vec, j)
      }
    }

    # Assemble sparse transfer matrix in one call
    Tm_sp <- Matrix::sparseMatrix(i = i_vec, j = j_vec,
                                  x = 1, dims = c(nsub, nsub),
                                  dimnames = list(comid, comid))
  }

  # --- route flows ------------------------------------------------------------
  AccumMatrix_sp <- .build_accum_matrix(Tm_sp)
  QR_mat <- as.matrix(QS %*% AccumMatrix_sp)

  # --- optional save ----------------------------------------------------------
  if (Save) {
    dir.create("./Outputs", showWarnings = FALSE, recursive = TRUE)

    yyyymm_new <- format(tail(Dates, 1), "%Y%m")
    if (Update) {
      prev_month <- lubridate::floor_date(tail(Dates, 1), "month") %m-% lubridate::months(1)
      old_file <- sprintf("./Outputs/QR_GR2MSemiDistr_%s.txt", format(prev_month, "%Y%m"))
      if (file.exists(old_file)) file.remove(old_file)
    }

    colnames(QR_mat) <- paste0("QR_", comid)
    rownames(QR_mat) <- format(Dates, "%Y-%m-01")

    QR_export <- as.matrix(round(QR_mat, 1))
    write.table(QR_export,
                file = sprintf("./Outputs/QR_GR2MSemiDistr_%s.txt", yyyymm_new),
                sep = "\t", quote = FALSE, col.names = NA)
  }

  message("Processing completed successfully in...")
  tictoc::toc()

  list(
    QR = QR_mat,
    Dates = Dates,
    COMID = comid,
    TransferMatrix = as.matrix(Tm_sp),
    AccumMatrix = as.matrix(AccumMatrix_sp)
  )
}
