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

  # Ensure WhiteboxTools is available
  if (!file.exists(whitebox::wbt_exe_path())) {
    message("Installing WhiteboxTools...")
    whitebox::install_whitebox()
  }
  whitebox::wbt_init()

  # Basic validations
  if (!inherits(Subbasins, "SpatVector")) stop("Subbasins must be a SpatVector.")
  if (!all(c("COMID", "Region") %in% names(Subbasins))) stop("'Subbasins' needs fields 'COMID' and 'Region'.")
  if (!is.null(Dem)) {
    if (!inherits(Dem, "SpatRaster")) stop("Dem must be a SpatRaster.")
    if (!terra::crs(Dem) == terra::crs(Subbasins)) stop("CRS of 'Dem' and 'Subbasins' must match.")
  }
  if (!is.null(RouIni) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouIni)) stop("RouIni must be 'mm/YYYY'.")
  if (!is.null(RouEnd) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouEnd)) stop("RouEnd must be 'mm/YYYY'.")
  if (!is.null(TransferMatrixFile) && !grepl("\\.rds$", TransferMatrixFile)) stop("TransferMatrixFile must be .rds.")

  # Clean temporary files on exit
  on.exit({
    if (is.null(TransferMatrixFile)) unlink(c("dem_tmp.tif", "dem_filled.tif", "flow_dir.tif"), force = TRUE)
  }, add = TRUE)

  # Clean and prepare Subbasins
  Subbasins <- Subbasins[!is.na(Subbasins$COMID), ]
  if (any(duplicated(Subbasins$COMID))) stop("Duplicated COMID found in 'Subbasins'.")

  # SAFE coercion: preserve real IDs if COMID is factor/character
  if (is.factor(Subbasins$COMID)) {
    Subbasins$COMID <- as.numeric(as.character(Subbasins$COMID))
  } else if (is.character(Subbasins$COMID)) {
    Subbasins$COMID <- as.numeric(Subbasins$COMID)
  } # if already numeric, keep as is

  comid <- Subbasins$COMID

  # Subset QS by time window
  if ((is.null(RouIni) && is.null(RouEnd)) || nrow(Model$QS) == 1) {
    Dates <- Model$Dates
    QS    <- as.matrix(Model$QS)
  } else {
    Dates <- seq(as.Date(paste0('01/', RouIni), format = '%d/%m/%Y'),
                 as.Date(paste0('01/', RouEnd), format = '%d/%m/%Y'), by = 'months')
    Ind   <- seq(which(format(as.Date(Model$Dates), '%d/%m/%Y') == paste0('01/', RouIni)),
                 which(format(as.Date(Model$Dates), '%d/%m/%Y') == paste0('01/', RouEnd)))
    QS <- as.matrix(Model$QS[Ind, ])
  }

  # Branch A: reuse transfer matrix from file
  if (!is.null(TransferMatrixFile) && file.exists(TransferMatrixFile)) {
    transfer_matrix <- readRDS(TransferMatrixFile)
    tm_names <- colnames(transfer_matrix)

    idx_keep <- match(tm_names, as.character(comid))
    if (any(is.na(idx_keep))) stop("Transfer matrix COMID do not fully match Subbasins/Model COMID.")

    QS_aligned <- QS[, idx_keep, drop = FALSE]
    colnames(QS_aligned) <- tm_names

    QR <- QS_aligned %*% transfer_matrix
    colnames(QR) <- rownames(transfer_matrix)
    QR <- round(QR, 1)
    routed_comid <- rownames(transfer_matrix)

  } else {
    if (is.null(Dem)) stop("Provide 'Dem' when 'TransferMatrixFile' is not given.")

    # Resample DEM to requested resolution (meters) or convert meters->degrees if lon/lat
    if (terra::is.lonlat(Dem)) {
      message("DEM is in geographic coordinates (degrees). Converting res from meters to degrees for resampling...")
      lat_mean <- mean(c(terra::ymin(Dem), terra::ymax(Dem)))
      res_deg  <- res / (111320 * cos(lat_mean * pi / 180))
      target_rast <- terra::rast(ext(Dem), resolution = res_deg, crs = crs(Dem))
      Dem <- terra::resample(Dem, target_rast, method = "bilinear")
    } else {
      res_dem <- terra::res(Dem)
      if (any(abs(res_dem - res) > 1e-6)) {
        message("Resampling DEM to ", res, " m resolution...")
        target_rast <- terra::rast(ext(Dem), resolution = res, crs = crs(Dem))
        Dem <- terra::resample(Dem, target_rast, method = "bilinear")
      } else {
        message("DEM already at desired resolution: skipping resampling.")
      }
    }

    # Rasterize COMID; touches=TRUE preserves small polygons
    comid_rast <- terra::rasterize(Subbasins, Dem, field = "COMID", touches = TRUE)

    # Drop categorical encoding and treat 0 as NA (some drivers coerce NA to 0)
    levs <- terra::levels(comid_rast)
    if (!is.null(levs) && length(levs) > 0 && nrow(levs[[1]]) > 0) levels(comid_rast) <- NULL
    comid_rast <- terra::ifel(comid_rast == 0, NA, comid_rast)

    # Flow directions via WhiteboxTools
    terra::writeRaster(Dem, "dem_tmp.tif", overwrite = TRUE)
    whitebox::wbt_fill_depressions("dem_tmp.tif", "dem_filled.tif")
    whitebox::wbt_d8_pointer("dem_filled.tif", "flow_dir.tif")
    flow_dir_rast <- terra::rast("flow_dir.tif")

    # COMID present on the grid
    comid_vals <- sort(unique(na.omit(as.vector(terra::values(comid_rast)))))
    n_sub_grid <- length(comid_vals)
    if (n_sub_grid == 0) stop("No subbasins were rasterized onto the DEM. Check resolution/CRS and polygon sizes.")

    n_sub_shp <- length(unique(Subbasins$COMID))
    if (n_sub_grid != n_sub_shp) {
      message("Notice: ", n_sub_shp, " COMID in 'Subbasins' vs ", n_sub_grid, " COMID rasterized on the grid.")
    }

    # Align QS columns to grid COMID order
    idx_keep <- match(comid_vals, Subbasins$COMID)
    if (any(is.na(idx_keep))) {
      missing_in_model <- comid_vals[is.na(idx_keep)]
      stop("These COMID exist on the grid but not in Model/Subbasins: ", paste(missing_in_model, collapse = ", "))
    }
    QS_aligned <- QS[, idx_keep, drop = FALSE]
    colnames(QS_aligned) <- as.character(comid_vals)

    # Initialize sparse transfer matrix with consistent dims and names
    transfer_matrix <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0),
      dims = c(n_sub_grid, n_sub_grid),
      dimnames = list(as.character(comid_vals), as.character(comid_vals))
    )

    # Fast index maps COMID -> row/col
    row_index <- setNames(seq_len(n_sub_grid), as.character(comid_vals))
    col_index <- row_index

    # D8 direction code to row/col offsets
    d8_map <- list(
      `1` = c(0, 1), `2` = c(1, 1), `4` = c(1, 0), `8` = c(1, -1),
      `16` = c(0, -1), `32` = c(-1, -1), `64` = c(-1, 0), `128` = c(-1, 1)
    )

    nrows <- nrow(flow_dir_rast)
    ncols <- ncol(flow_dir_rast)

    # Populate transfer matrix (rows=receivers, cols=sources)
    for (cell in seq_len(ncell(comid_rast))) {
      current_comid <- comid_rast[cell]
      if (is.na(current_comid)) next

      flow_dir <- flow_dir_rast[cell]
      key <- as.character(flow_dir)
      if (is.na(flow_dir) || !key %in% names(d8_map)) next

      rc <- rowColFromCell(flow_dir_rast, cell)
      off <- d8_map[[key]]
      r2  <- rc[1] + off[1]
      c2  <- rc[2] + off[2]
      if (r2 < 1 || r2 > nrows || c2 < 1 || c2 > ncols) next

      cell2 <- cellFromRowCol(flow_dir_rast, r2, c2)
      next_comid <- comid_rast[cell2]
      if (is.na(next_comid) || next_comid == current_comid) next

      i <- row_index[as.character(next_comid)]
      j <- col_index[as.character(current_comid)]
      if (!is.na(i) && !is.na(j)) transfer_matrix[i, j] <- 1
    }

    # Column-normalize (each source distributes weights to receivers)
    cs <- Matrix::colSums(transfer_matrix)
    cs[cs == 0] <- 1
    transfer_matrix <- transfer_matrix %*% Matrix::Diagonal(x = 1 / cs)
    transfer_matrix@x[is.nan(transfer_matrix@x)] <- 0

    # Routed flows
    QR <- QS_aligned %*% transfer_matrix
    colnames(QR) <- rownames(transfer_matrix)
    QR <- round(QR, 1)
    routed_comid <- rownames(transfer_matrix)
  }

  # Save outputs if requested
  if (Save) {
    dir.create("./Outputs", showWarnings = FALSE, recursive = TRUE)
    last_month <- format(tail(Dates, 1), "%Y%m")
    out_file <- paste0("./Outputs/QR_GR2MSemiDistr_", last_month, ".txt")

    if (Update) {
      last_date <- as.Date(paste0(last_month, "01"), "%Y%m%d")
      prev_month <- format(last_date %m-% months(1), "%Y%m")
      old_file <- paste0("./Outputs/QR_GR2MSemiDistr_", prev_month, ".txt")
      if (file.exists(old_file)) file.remove(old_file)
    }

    df <- as.data.frame(QR)
    colnames(df) <- paste0("QR_", routed_comid)
    rownames(df) <- Dates
    write.table(df, file = out_file)
  }

  message("Processing completed successfully in...")
  tictoc::toc()

  list(QR = QR, Dates = Dates, COMID = colnames(QR), TransferMatrix = transfer_matrix)
}
