#' Extract and prepare model input data in the airGR format (DatesR, P, E, and optionally Q) using gridded monthly precipitation and evapotranspiration.
#'
#' This function extracts spatially averaged monthly precipitation and potential evapotranspiration values for each subbasin,
#' and optionally includes observed streamflow data, to produce a formatted data frame compatible with the airGR hydrological modeling framework.
#'
#' @param Subbasins A SpatVector object containing the subbasins' geometries. Must include the following attributes: 'COMID' (unique identifier), and 'Region' (string code for the hydro-climatic region).
#' @param Precip A SpatRaster object with monthly precipitation data in [mm/month].
#' @param PotEvap A SpatRaster object with monthly potential evapotranspiration data in [mm/month].
#' @param Qobs Optional. A vector with observed streamflow in [m³/s] at the basin outlet.
#' @param DateIni Start date of the output time series in 'mm/yyyy' format.
#' @param DateEnd End date of the output time series in 'mm/yyyy' format.
#' @param IniGrids Initial date of the gridded precipitation and evapotranspiration datasets in mm/yyyy format.
#' @param Save Logical. If TRUE, saves the resulting data frame to a .txt file in the './Inputs' directory. Default is FALSE.
#' @param Update Logical. If TRUE, returns only the last month of data (useful for operational updates). Default is FALSE.
#' @param Members Optional. Integer indicating the number of ensemble members for streamflow forecasting. If provided, the date vector will be repeated accordingly. Default is NULL.
#'
#' @return A data frame containing the model inputs in the airGR format, with columns:
#' \describe{
#'   \item{DatesR}{Monthly dates}
#'   \item{P_x}{Mean precipitation for subbasin x, where x is the COMID}
#'   \item{E_x}{Mean potential evapotranspiration for subbasin x, where x is the COMID}
#'   \item{Q}{Observed streamflow [m³/s], if provided}
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
#' # Create figures to visualize the inputs
#' dates <- as.Date(model_inputs$DatesR)
#'
#' # Plot average precipitation across subbasins
#' P_vars <- grep("^P_", names(model_inputs), value = TRUE)
#' avg_P <- rowMeans(model_inputs[, P_vars], na.rm = TRUE)
#' plot(dates, avg_P, type = "l", col = "blue", lwd = 1.5,
#'      ylab = "Precipitation [mm/month]", xlab = "Date",
#'      main = "Mean Precipitation across Subbasins")
#'
#' # Plot average potential evapotranspiration across subbasins
#' E_vars <- grep("^E_", names(model_inputs), value = TRUE)
#' avg_E <- rowMeans(model_inputs[, E_vars], na.rm = TRUE)
#' plot(dates, avg_E, type = "l", col = "orange", lwd = 1.5,
#'      ylab = "PET [mm/month]", xlab = "Date",
#'      main = "Mean Potential Evapotranspiration across Subbasins")
#'
#' # Plot observed streamflow at outlet (if available)
#' if ("Q" %in% names(model_inputs)) {
#'   plot(dates, model_inputs$Q, type = "l", col = "darkgreen", lwd = 1.5,
#'        ylab = "Observed Streamflow [m³/s]", xlab = "Date",
#'        main = "Observed Streamflow at Outlet")
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
    Prec <- .extract_variable(Precip, "precipitation", IniGrids, DateIni, DateEnd,
                              Subbasins_sf, Update)
  }

  if (!is.null(PotEvap)) {
    Evap <- .extract_variable(PotEvap, "evapotranspiration", IniGrids, DateIni, DateEnd,
                              Subbasins_sf, Update)
  }

  # ==== Date vector ====
  if (Update) {
    DatesMonths <- as.Date(paste0("01/", DateEnd), "%d/%m/%Y")
  } else {
    ini <- as.Date(paste0("01/", DateIni),  "%d/%m/%Y")
    end <- as.Date(paste0("01/", DateEnd),  "%d/%m/%Y")
    DatesMonths <- seq(ini, end, by = "month")
    if (!is.null(Members)) {
      DatesMonths <- rep(DatesMonths, times = Members) # Members handling not modified
    }
  }

  # ==== Assemble output data.frame ====
  data_list <- list(Date = DatesMonths)

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
                     format(min(Ans$Date), "%Y%m"),
                     format(max(Ans$Date), "%Y%m"))
    write.table(Ans, file = fname, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  message("Processing completed successfully in...")
  tictoc::toc()
  return(Ans)
}
