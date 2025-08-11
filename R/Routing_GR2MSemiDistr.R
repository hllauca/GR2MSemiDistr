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

  # === Ensure WhiteboxTools binary ===
  if (!file.exists(whitebox::wbt_exe_path())) {
    message("Installing WhiteboxTools...")
    whitebox::install_whitebox()
  }
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
    if (!terra::same.crs(Dem, Subbasins)) stop("CRS of 'Dem' and 'Subbasins' must match.")
  }
  if (!is.null(RouIni) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouIni)) stop("RouIni must be 'mm/YYYY'.")
  if (!is.null(RouEnd) && !grepl("^(0[1-9]|1[0-2])/\\d{4}$", RouEnd)) stop("RouEnd must be 'mm/YYYY'.")
  if (!is.null(TransferMatrixFile) && !grepl("\\.rds$", TransferMatrixFile)) stop("TransferMatrixFile must be .rds")

  # === IDs & sizes ===
  comid <- as.character(Subbasins$COMID)
  nsub  <- length(comid)

  # === Subset QS by dates (if requested) ===
  if ((is.null(RouIni) && is.null(RouEnd)) || nrow(Model$QS) == 1) {
    Dates <- as.Date(Model$Dates)
    QS    <- as.matrix(Model$QS)
  } else {
    all_dates <- as.Date(Model$Dates)
    d0 <- as.Date(paste0("01/", RouIni), format = "%d/%m/%Y")
    d1 <- as.Date(paste0("01/", RouEnd),  format = "%d/%m/%Y")
    i0 <- match(d0, all_dates); i1 <- match(d1, all_dates)
    if (any(is.na(c(i0,i1)))) stop("RouIni/RouEnd not found in Model$Dates.")
    Dates <- seq(d0, d1, by = "month")
    QS    <- as.matrix(Model$QS[i0:i1, , drop = FALSE])
  }

  # Ensure QS columns align with COMID order if names exist
  if (!is.null(colnames(QS))) {
    miss <- setdiff(comid, colnames(QS))
    if (length(miss)) stop("QS columns missing COMID(s): ", paste(miss, collapse=", "))
    QS <- QS[, comid, drop = FALSE]
  } else {
    # assume QS already in Subbasins order
    colnames(QS) <- comid
  }

  # === Try to reuse transfer matrix ===
  if (!is.null(TransferMatrixFile) && file.exists(TransferMatrixFile)) {
    transfer_matrix <- readRDS(TransferMatrixFile)
    # basic checks
    if (!inherits(transfer_matrix, "dgCMatrix")) stop("TransferMatrixFile must be a sparse matrix (dgCMatrix).")
    if (!identical(rownames(transfer_matrix), comid) ||
        !identical(colnames(transfer_matrix), comid)) {
      stop("Transfer matrix dimnames must match Subbasins COMID order.")
    }
    QR <- QS %*% transfer_matrix
  } else {
    if (is.null(Dem)) stop("Provide 'Dem' when 'TransferMatrixFile' is not provided.")

    # === DEM resampling ===
    if (terra::is.lonlat(Dem)) {
      message("DEM is lon/lat. Converting target resolution from meters to degrees for resampling...")
      lat_mean <- mean(c(terra::ymin(Dem), terra::ymax(Dem)))
      res_deg  <- res / (111320 * cos(lat_mean * pi / 180))
      target_rast <- terra::rast(terra::ext(Dem), resolution = res_deg, crs = terra::crs(Dem))
      Dem <- terra::resample(Dem, target_rast, method = "bilinear")
    } else {
      res_dem <- terra::res(Dem) # c(xres, yres)
      if (any(abs(res_dem - c(res, res)) > 1e-6)) {
        message("Resampling DEM to ", res, " m...")
        target_rast <- terra::rast(terra::ext(Dem), resolution = res, crs = terra::crs(Dem))
        Dem <- terra::resample(Dem, target_rast, method = "bilinear")
      } else {
        message("DEM already at desired resolution: skipping resampling.")
      }
    }

    # === Rasterize COMID over DEM grid ===
    comid_rast <- terra::rasterize(Subbasins, Dem, field = "COMID", touches = TRUE)

    # === Temp files ===
    dem_tmp  <- tempfile(fileext = ".tif")
    dem_fill <- tempfile(fileext = ".tif")
    com_tmp  <- tempfile(fileext = ".tif")
    fdir_tif <- tempfile(fileext = ".tif")

    # cleanup
    on.exit({
      unlink(c(dem_tmp, dem_fill, com_tmp, fdir_tif), force = TRUE)
    }, add = TRUE)

    terra::writeRaster(Dem, dem_tmp, overwrite = TRUE)
    terra::writeRaster(comid_rast, com_tmp, overwrite = TRUE)

    whitebox::wbt_fill_depressions(dem_tmp, dem_fill)
    whitebox::wbt_d8_pointer(dem_fill, fdir_tif)

    flow_dir_rast <- terra::rast(fdir_tif)
    comid_rast    <- terra::rast(com_tmp)

    # === Build transfer matrix (column-normalized WFA) ===
    transfer_matrix <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0),
      dims = c(nsub, nsub),
      dimnames = list(comid, comid)
    )

    d8_map <- list(
      `1` = c(0, 1), `2` = c(1, 1), `4` = c(1, 0), `8` = c(1, -1),
      `16` = c(0, -1), `32` = c(-1, -1), `64` = c(-1, 0), `128` = c(-1, 1)
    )

    nrows <- terra::nrow(flow_dir_rast)
    ncols <- terra::ncol(flow_dir_rast)
    ncell_tot <- terra::ncell(comid_rast)

    for (cell_id in seq_len(ncell_tot)) {
      current_comid <- comid_rast[cell_id]
      if (is.na(current_comid)) next
      flow_dir <- flow_dir_rast[cell_id]
      key <- as.character(flow_dir)
      if (is.na(flow_dir) || !(key %in% names(d8_map))) next

      rc <- terra::rowColFromCell(flow_dir_rast, cell_id)
      off <- d8_map[[key]]
      r2 <- rc[1] + off[1]; c2 <- rc[2] + off[2]
      if (r2 < 1 || r2 > nrows || c2 < 1 || c2 > ncols) next

      dest_cell  <- terra::cellFromRowCol(flow_dir_rast, r2, c2)
      next_comid <- comid_rast[dest_cell]

      if (!is.na(next_comid) && next_comid != current_comid) {
        # Use Subbasins-defined COMID order
        i_name <- as.character(next_comid)
        j_name <- as.character(current_comid)
        if (i_name %in% comid && j_name %in% comid) {
          transfer_matrix[i_name, j_name] <- 1
        }
      }
    }

    # Normalize columns (split flow equally among multiple downstreams)
    col_sums <- Matrix::colSums(transfer_matrix)
    col_sums[col_sums == 0] <- 1
    transfer_matrix <- transfer_matrix %*% Matrix::Diagonal(x = as.numeric(1 / col_sums))

    if (!is.null(TransferMatrixFile)) {
      saveRDS(transfer_matrix, TransferMatrixFile)
    }

    QR <- QS %*% transfer_matrix
  }

  # === Optional save ===
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
  list(QR = QR, Dates = Dates, COMID = comid, TransferMatrix = transfer_matrix)
}
