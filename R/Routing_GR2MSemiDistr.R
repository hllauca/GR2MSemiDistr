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

  # === Ensure WhiteboxTools is installed and available ===
  if (!file.exists(whitebox::wbt_exe_path())) {
    message("Installing WhiteboxTools...")
    whitebox::install_whitebox()
  }
  # Force initialization after installation to avoid NULL pointer errors
  whitebox::wbt_init()

  # === Validate inputs ===
  if (!inherits(Subbasins, "SpatVector")) {
    stop("Argument 'Subbasins' must be of class 'SpatVector'.")
  }
  if (!all(c("COMID", "Region") %in% names(Subbasins))) {
    stop("The 'Subbasins' object must contain fields 'COMID' and 'Region'.")
  }
  if (!is.null(Dem)) {
    if (!inherits(Dem, "SpatRaster")) stop("Argument 'Dem' must be of class 'SpatRaster'.")
    if (!terra::crs(Dem) == terra::crs(Subbasins)) stop("CRS of 'Dem' and 'Subbasins' must match.")
  }
  if (!is.null(RouIni) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouIni)) stop("RouIni must be in 'mm/YYYY' format.")
  if (!is.null(RouEnd) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouEnd)) stop("RouEnd must be in 'mm/YYYY' format.")
  if (!is.null(TransferMatrixFile) && !grepl("\\.rds$", TransferMatrixFile)) stop("TransferMatrixFile must be an .rds file.")

  # === Remove temporary files on exit or error ===
  on.exit({
    if (is.null(TransferMatrixFile)) {
      unlink(c("dem_tmp.tif", "dem_filled.tif", "comid_tmp.tif", "flow_dir.tif"), force = TRUE)
    }
  }, add = TRUE)

  # === Clean and prepare Subbasins ===
  # Ensure no NA or duplicated COMID to keep a one-to-one mapping
  Subbasins <- Subbasins[!is.na(Subbasins$COMID), ]
  if (any(duplicated(Subbasins$COMID))) {
    stop("Duplicated COMID found in 'Subbasins'. Ensure unique COMID values.")
  }
  comid <- Subbasins$COMID

  # === Subset QS from model time window ===
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

  # === Reuse transfer matrix if available ===
  if (!is.null(TransferMatrixFile) && file.exists(TransferMatrixFile)) {
    transfer_matrix <- readRDS(TransferMatrixFile)

    # Align QS columns to transfer_matrix column names if possible
    tm_names <- colnames(transfer_matrix)
    if (!all(as.character(comid) %in% tm_names)) {
      # Reorder/Subset columns of QS to the transfer matrix naming
      idx_keep <- match(tm_names, as.character(comid))
      if (any(is.na(idx_keep))) {
        stop("Transfer matrix COMID do not fully match Subbasins/Model COMID.")
      }
      QS_aligned <- QS[, idx_keep, drop = FALSE]
    } else {
      QS_aligned <- QS
      colnames(QS_aligned) <- as.character(comid)
      QS_aligned <- QS_aligned[, tm_names, drop = FALSE]
    }

    QR <- QS_aligned %*% transfer_matrix
    QR <- round(QR, 1)
    routed_comid <- rownames(transfer_matrix)

  } else {
    if (is.null(Dem)) stop("You must provide a raster object 'Dem' if 'TransferMatrixFile' is not provided.")

    # === Resample DEM to requested resolution (meters) or convert meters->degrees in lon/lat ===
    if (terra::is.lonlat(Dem)) {
      message("DEM is in geographic coordinates (degrees). Converting res from meters to degrees for resampling...")
      lat_mean <- mean(c(terra::ymin(Dem), terra::ymax(Dem)))  # correct averaging
      res_deg  <- res / (111320 * cos(lat_mean * pi / 180))
      target_rast <- terra::rast(ext(Dem), resolution = res_deg, crs = crs(Dem))
      Dem <- terra::resample(Dem, target_rast, method = "bilinear")
    } else {
      res_dem <- terra::res(Dem)
      # If DEM resolution differs from desired 'res' (tolerate small numeric noise)
      if (any(abs(res_dem - res) > 1e-6)) {
        message("Resampling DEM to ", res, " m resolution...")
        target_rast <- terra::rast(ext(Dem), resolution = res, crs = crs(Dem))
        Dem <- terra::resample(Dem, target_rast, method = "bilinear")
      } else {
        message("DEM already at desired resolution: skipping resampling.")
      }
    }

    # === Rasterize COMID onto DEM grid; touches=TRUE to retain small polygons ===
    comid_rast <- terra::rasterize(Subbasins, Dem, field = "COMID", touches = TRUE)

    # === Write temporary rasters for WhiteboxTools ===
    terra::writeRaster(Dem, "dem_tmp.tif", overwrite = TRUE)
    terra::writeRaster(comid_rast, "comid_tmp.tif", overwrite = TRUE)

    # === Derive flow directions using WhiteboxTools ===
    whitebox::wbt_fill_depressions("dem_tmp.tif", "dem_filled.tif")
    whitebox::wbt_d8_pointer("dem_filled.tif", "flow_dir.tif")

    # === Read back results as SpatRaster ===
    flow_dir_rast <- terra::rast("flow_dir.tif")
    comid_rast    <- terra::rast("comid_tmp.tif")

    # === Determine COMID actually present on the grid ===
    comid_vals <- sort(unique(na.omit(as.vector(terra::values(comid_rast)))))
    n_sub_grid <- length(comid_vals)
    if (n_sub_grid == 0) stop("No subbasins were rasterized onto the DEM. Check resolution/CRS and polygon sizes.")

    # Inform if some subbasins were not rasterized (common with very small polygons)
    n_sub_shp <- length(unique(Subbasins$COMID))
    if (n_sub_grid != n_sub_shp) {
      message("Notice: ", n_sub_shp, " COMID in 'Subbasins' vs ", n_sub_grid, " COMID rasterized on the grid.")
    }

    # === Align QS columns to the order of comid_vals (grid COMID order) ===
    comid_vec_model <- Subbasins$COMID
    idx_keep <- match(comid_vals, comid_vec_model)
    if (any(is.na(idx_keep))) {
      missing_in_model <- comid_vals[is.na(idx_keep)]
      stop("These COMID exist on the grid but not in Model/Subbasins: ",
           paste(missing_in_model, collapse = ", "))
    }
    QS_aligned <- QS[, idx_keep, drop = FALSE]
    colnames(QS_aligned) <- as.character(comid_vals)

    # === Initialize sparse transfer matrix with consistent dims and dimnames ===
    transfer_matrix <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0),
      dims = c(n_sub_grid, n_sub_grid),
      dimnames = list(as.character(comid_vals), as.character(comid_vals))
    )

    # === Build fast index maps COMID -> row/col positions ===
    row_index <- setNames(seq_len(n_sub_grid), as.character(comid_vals))
    col_index <- row_index

    # === D8 direction code to row/col offset mapping ===
    d8_map <- list(
      `1` = c(0, 1), `2` = c(1, 1), `4` = c(1, 0), `8` = c(1, -1),
      `16` = c(0, -1), `32` = c(-1, -1), `64` = c(-1, 0), `128` = c(-1, 1)
    )

    nrows <- nrow(flow_dir_rast)
    ncols <- ncol(flow_dir_rast)

    # === Populate transfer matrix by scanning each cell ===
    # Note: rows = receivers (downstream), cols = sources (upstream)
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
      if (!is.na(i) && !is.na(j)) {
        transfer_matrix[i, j] <- 1
      }
    }

    # === Column-normalize so each source distributes weights to its receivers ===
    cs <- Matrix::colSums(transfer_matrix)
    cs[cs == 0] <- 1
    transfer_matrix <- transfer_matrix %*% Matrix::Diagonal(x = 1 / cs)
    transfer_matrix@x[is.nan(transfer_matrix@x)] <- 0

    # === Optionally save the transfer matrix ===
    if (!is.null(TransferMatrixFile)) {
      saveRDS(transfer_matrix, TransferMatrixFile)
    }

    # === Compute routed flows ===
    QR <- QS_aligned %*% transfer_matrix
    QR <- round(QR, 1)
    routed_comid <- as.character(comid_vals)
  }

  # === Save outputs, if requested ===
  if (Save) {
    dir.create("./Outputs", showWarnings = FALSE, recursive = TRUE)

    # Build filename using last date in 'Dates'
    last_month <- format(tail(Dates, 1), "%Y%m")
    out_file <- paste0("./Outputs/QR_GR2MSemiDistr_", last_month, ".txt")

    # If Update=TRUE, remove previous month's file (if exists)
    if (Update) {
      # Compute previous month safely with lubridate
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

  return(list(
    QR = QR,
    Dates = Dates,
    COMID = colnames(QR),
    TransferMatrix = transfer_matrix
  ))
}
