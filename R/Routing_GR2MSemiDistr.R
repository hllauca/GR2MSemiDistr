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
#' @param res Numeric. Desired spatial resolution in meters for resampling the DEM (default is 250). Converted to degrees if needed.
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
                                  res = 250,
                                  Save = FALSE,
                                  Update = FALSE) {
  tictoc::tic()

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

  # === Input preparation ===
  comid <- Subbasins$COMID
  nsub  <- length(comid)

  # === Subset QS from model ===
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
    QR <- QS %*% transfer_matrix
    QR <- round(QR, 1)
  } else {
    if (is.null(Dem)) stop("You must provide a raster object 'Dem' if 'TransferMatrixFile' is not provided.")

    # Resample using approximate conversion to degrees if in lon/lat
    if (terra::is.lonlat(Dem)) {
      message("DEM is in geographic coordinates (degrees). Converting res from meters to degrees for resampling...")
      lat_mean <- mean(terra::ymin(Dem), terra::ymax(Dem))
      res_deg <- res / (111320 * cos(lat_mean * pi / 180))
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

    comid_rast <- terra::rasterize(Subbasins, Dem, field = "COMID", touches = FALSE)

    terra::writeRaster(Dem, "dem_tmp.tif", overwrite = TRUE)
    terra::writeRaster(comid_rast, "comid_tmp.tif", overwrite = TRUE)

    whitebox::wbt_fill_depressions("dem_tmp.tif", "dem_filled.tif")
    whitebox::wbt_d8_pointer("dem_filled.tif", "flow_dir.tif")

    flow_dir_rast <- terra::rast("flow_dir.tif")
    comid_rast    <- terra::rast("comid_tmp.tif")
    comid_vals    <- sort(na.omit(unique(terra::values(comid_rast))))

    transfer_matrix <- sparseMatrix(i = integer(0), j = integer(0),
                                    dims = c(nsub, nsub),
                                    dimnames = list(as.character(comid_vals), as.character(comid_vals)))

    d8_map <- list(
      `1` = c(0, 1), `2` = c(1, 1), `4` = c(1, 0), `8` = c(1, -1),
      `16` = c(0, -1), `32` = c(-1, -1), `64` = c(-1, 0), `128` = c(-1, 1)
    )

    nrows <- nrow(flow_dir_rast)
    ncols <- ncol(flow_dir_rast)

    for (i in seq_len(ncell(comid_rast))) {
      current_comid <- comid_rast[i]
      if (is.na(current_comid)) next

      flow_dir <- flow_dir_rast[i]
      if (is.na(flow_dir) || !(as.character(flow_dir) %in% names(d8_map))) next

      rowcol <- rowColFromCell(flow_dir_rast, i)
      offset <- d8_map[[as.character(flow_dir)]]
      dest_row <- rowcol[1] + offset[1]
      dest_col <- rowcol[2] + offset[2]

      if (dest_row < 1 || dest_row > nrows || dest_col < 1 || dest_col > ncols) next

      dest_cell <- cellFromRowCol(flow_dir_rast, dest_row, dest_col)
      next_comid <- comid_rast[dest_cell]

      if (!is.na(next_comid) && next_comid != current_comid) {
        transfer_matrix[as.character(next_comid), as.character(current_comid)] <- 1
      }
    }

    col_sums <- colSums(transfer_matrix)
    col_sums[col_sums == 0] <- 1
    transfer_matrix <- transfer_matrix %*% Diagonal(x = 1 / col_sums)

    if (!is.null(TransferMatrixFile)) {
      saveRDS(transfer_matrix, TransferMatrixFile)
    }

    QR <- QS %*% transfer_matrix
    QR <- round(QR, 1)
  }

  if (Save) {
    dir.create("./Outputs", showWarnings = FALSE, recursive = TRUE)

    if (Update) {
      new_month <- format(tail(Dates, 1), "%Y%m")
      old_month <- format(as.Date(paste0("01", new_month), "%d%Y%m") %m-% months(1), "%Y%m")
      old_file <- paste0("./Outputs/QR_GR2MSemiDistr_", old_month, ".txt")
      if (file.exists(old_file)) file.remove(old_file)
    }

    df <- as.data.frame(QR)
    colnames(df) <- paste0("QR_", comid)
    rownames(df) <- Dates
    out_file <- paste0("./Outputs/QR_GR2MSemiDistr_", format(tail(Dates, 1), "%Y%m"), ".txt")
    write.table(df, file = out_file)
  }

  message("Processing completed successfully in...")
  tictoc::toc()
  return(list(QR = QR, Dates = Dates, COMID = colnames(QR), TransferMatrix = transfer_matrix))
}
