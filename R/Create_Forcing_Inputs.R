#' Extract and prepare model input data in the airGR format (DatesR, P, E, and optionally Q) using gridded monthly precipitation and evapotranspiration.
#'
#' This function extracts spatially averaged monthly precipitation and potential evapotranspiration values for each subbasin,
#' and optionally includes observed streamflow data, to produce a formatted data frame compatible with the airGR hydrological modeling framework.
#'
#' @param Subbasins SpatVector. Subbasin geometries. Must include attributes:
#' 'COMID' (unique subbasin identifier) and 'Region' (character code for the hydro-climatic region).
#' @param Precip SpatRaster. Monthly precipitation fields in [mm/month].
#' Dimensions and time coverage must include at least the requested period (`DateIni`–`DateEnd`).
#' @param PotEvap SpatRaster. Monthly potential evapotranspiration fields in [mm/month].
#' Must align with the precipitation raster in space and time.
#' @param Qobs numeric vector (optional). Observed streamflow series at the basin outlet, in [m³/s].
#' Length must match the number of simulated months if provided. Default is `NULL`.
#' @param DateIni character. Start date of the output time series in `"mm/yyyy"` format.
#' @param DateEnd character. End date of the output time series in `"mm/yyyy"` format.
#' @param IniGrids character. Initial date of the gridded precipitation and evapotranspiration datasets in `"mm/yyyy"` format.
#' Used to align raster time indices with the requested simulation period.
#' @param Save logical. If `TRUE`, saves the resulting data frame to a tab-separated `.txt` file
#' inside the `./Inputs` directory. Default is `FALSE`.
#' @param Update logical. If `TRUE`, returns only the last month of data
#' (useful for operational forecasting/updates). Default is `FALSE`.
#' @param Members integer (optional). Number of ensemble members for streamflow forecasting.
#' If provided, the date vector is replicated accordingly. Default is `NULL`.
#'
#' @return data.frame with model inputs in the airGR format, including the following columns:
#' \describe{
#'   \item{DatesR}{POSIXct. Monthly dates covering the simulation period.}
#'   \item{P_x}{numeric. Mean precipitation for subbasin x [mm/month], where x is the COMID.}
#'   \item{E_x}{numeric. Mean potential evapotranspiration for subbasin x [mm/month], where x is the COMID.}
#'   \item{Q}{numeric (optional). Observed streamflow at the basin outlet [m³/s], included only if `Qobs` was provided.}
#' }
#'
#' @references
#' Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. (2021).
#' PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020).
#' Water, 13(8), 1048. \doi{10.3390/w13081048}
#'
#' @examples
#' library(GR2MSemiDistr)
#'
#' # Define simulation period and grid start date
#' start_date <- "01/1981"
#' end_date   <- "12/2016"
#' ini_grids  <- "01/1981"
#'
#' # Create input data for the model
#' model_inputs <- Create_Forcing_Inputs(
#'   Subbasins = cat,
#'   Precip = grid_pr,
#'   PotEvap = grid_pe,
#'   Qobs = qobs,
#'   DateIni = start_date,
#'   DateEnd = end_date,
#'   IniGrids = ini_grids,
#' )
#'
#' # Select a target COMID to plot
#' target_comid <- "10"
#'
#' # Plot precipitation for the target subbasin
#' P_col <- paste0("P_", target_comid)
#' plot(model_inputs$DatesR, model_inputs[[P_col]], type = "h", lwd = 2, col = "dodgerblue",
#'     ylab = "Precipitation [mm/month]", xlab = "Date",
#'     main = sprintf("Precipitation in Subbasin %s", target_comid))
#'
#' # Plot potential evapotranspiration for the target subbasin
#' E_col <- paste0("E_", target_comid)
#' plot(model_inputs$DatesR, model_inputs[[E_col]], type = "o", lwd = 2, col = "darkgreen",
#'      ylab = "PET [mm/month]", xlab = "Date",
#'      main = sprintf("Potential Evapotranspiration in Subbasin %s", target_comid))
#'
#' # Plot observed streamflow at outlet (if available)
#' if ("Q" %in% names(model_inputs)) {
#'   plot(model_inputs$DatesR, model_inputs$Q, type = "l", col = "darkblue", lwd = 1.5,
#'        ylab = "Q [m³/s]", xlab = "Date",
#'        main = "Observed Streamflow at Basin Outlet")
#' }
#' @import exactextractr
#' @import terra
#' @import sf
#' @import lubridate
#' @import tictoc
#'
#' @export
#'
Create_Forcing_Inputs <- function(Subbasins,
                                  Precip,
                                  PotEvap,
                                  Qobs = NULL,
                                  DateIni,
                                  DateEnd,
                                  IniGrids = "01/1981",
                                  Save = FALSE,
                                  Update = FALSE,
                                  Members = NULL) {

  tictoc::tic()

  # ==== Basic class validation ====
  if (!inherits(Subbasins, "SpatVector")) {
    stop("Argument 'Subbasins' must be of class 'SpatVector'.")
  }
  if (!"COMID" %in% names(Subbasins)) {
    stop("The 'Subbasins' object must contain a field named 'COMID'.")
  }
  if (!is.null(Precip) && !inherits(Precip, "SpatRaster")) {
    stop("Argument 'Precip' must be of class 'SpatRaster'.")
  }
  if (!is.null(PotEvap) && !inherits(PotEvap, "SpatRaster")) {
    stop("Argument 'PotEvap' must be of class 'SpatRaster'.")
  }

  # ==== Date format validation ====
  .check_mmYYYY <- function(x) grepl("^\\d{2}/\\d{4}$", x)
  if (!.check_mmYYYY(IniGrids)) stop("IniGrids must be 'mm/YYYY', e.g. '01/1981'.")
  if (!.check_mmYYYY(DateIni))  stop("DateIni must be 'mm/YYYY', e.g. '01/1990'.")
  if (!.check_mmYYYY(DateEnd))  stop("DateEnd must be 'mm/YYYY', e.g. '12/2020'.")

  # ==== Helper: robust monthly date sequence ====
  .seq_months_from_tag <- function(mmYYYY, n_layers) {
    start <- as.Date(paste0("01/", mmYYYY), "%d/%m/%Y")
    seq(from = start, by = "month", length.out = n_layers)
  }

  # ==== Helper: weighted areal mean extraction ====
  .extract_variable <- function(raster, var_name, IniGrids, DateIni, DateEnd,
                                Subbasins_sf, Update = FALSE) {

    message(sprintf("Computing monthly %s for each subbasin...", var_name))

    nL <- terra::nlyr(raster)
    all_dates <- .seq_months_from_tag(IniGrids, nL)

    ini <- as.Date(paste0("01/", DateIni), "%d/%m/%Y")
    end <- as.Date(paste0("01/", DateEnd), "%d/%m/%Y")

    # Select layer indices within the requested date range
    ind <- which(all_dates >= ini & all_dates <= end)
    if (!length(ind)) stop("Requested range does not intersect raster temporal coverage.")

    # If Update = TRUE, keep only the last date of the range
    if (Update) ind <- tail(ind, 1)

    # Weighted mean using polygon area
    mat <- t(
      exactextractr::exact_extract(raster[[ind]], Subbasins_sf,
                                   fun = "mean", weights = "area", progress = FALSE)
    )
    return(mat)
  }

  # ==== CRS checks and cropping/masking (if rasters are provided) ====
  if (!is.null(Precip)) {
    if (!terra::same.crs(Precip, Subbasins)) {
      stop("CRS mismatch between 'Precip' and 'Subbasins'. Reproject before calling.")
    }
    Precip <- terra::crop(Precip, terra::ext(Subbasins))
    Precip <- terra::mask(Precip, Subbasins)
  }
  if (!is.null(PotEvap)) {
    if (!terra::same.crs(PotEvap, Subbasins)) {
      stop("CRS mismatch between 'PotEvap' and 'Subbasins'. Reproject before calling.")
    }
    PotEvap <- terra::crop(PotEvap, terra::ext(Subbasins))
    PotEvap <- terra::mask(PotEvap, Subbasins)
  }

  # ==== Extract COMID IDs as character ====
  comid <- as.character(Subbasins$COMID)

  # ==== Data extraction ====
  Subbasins_sf <- sf::st_as_sf(Subbasins)

  if (!is.null(Precip)) {
    Prec <- round(.extract_variable(Precip, "precipitation", IniGrids, DateIni, DateEnd,
                              Subbasins_sf, Update),2)
  }

  if (!is.null(PotEvap)) {
    Evap <- round(.extract_variable(PotEvap, "evapotranspiration", IniGrids, DateIni, DateEnd,
                              Subbasins_sf, Update),2)
  }

  # ==== Date vector ====
  if (Update) {
    DatesMonths <- as.Date(paste0("01/", DateEnd), "%d/%m/%Y")
  } else {
    ini <- as.Date(paste0("01/", DateIni),  "%d/%m/%Y")
    end <- as.Date(paste0("01/", DateEnd),  "%d/%m/%Y")
    DatesMonths <- seq(ini, end, by = "month")
    if (!is.null(Members)) {
      DatesMonths <- rep(DatesMonths, times = Members)
      Prec <- Prec[rep(seq_len(nrow(Prec)), times = Members), , drop = FALSE]
      Evap <- Evap[rep(seq_len(nrow(Evap)), times = Members), , drop = FALSE]
    }
    
  }

  # ==== Assemble output data.frame ====
  data_list <- list(DatesR = DatesMonths)

  if (!is.null(Precip)) {
    Prec <- as.data.frame(Prec)
    colnames(Prec) <- paste0("P_", comid)
    data_list <- c(data_list, Prec)
  }

  if (!is.null(PotEvap)) {
    Evap <- as.data.frame(Evap)
    colnames(Evap) <- paste0("E_", comid)
    data_list <- c(data_list, Evap)
  }

  if (!is.null(Qobs)) {
    Flow <- as.data.frame(Qobs)
    if (ncol(Flow) == 1) names(Flow) <- "Q"
    data_list <- c(data_list, Flow)
  }

  Ans <- as.data.frame(data_list)

  # ==== Optional: save to file ====
  if (Save) {
    if (!dir.exists("./Inputs")) dir.create("./Inputs", recursive = TRUE)
    fname <- sprintf("./Inputs/Inputs_model_%s_%s.txt",
                     format(min(Ans$DatesR), "%Y%m"),
                     format(max(Ans$DatesR), "%Y%m"))
    write.table(Ans, file = fname, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  message("Processing completed successfully in...")
  tictoc::toc()
  return(Ans)
}
