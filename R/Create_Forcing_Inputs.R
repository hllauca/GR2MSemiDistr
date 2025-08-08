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
#' # Load example data
#'
#' data(cat)        # Subbasin geometries (SpatVector)
#' data(pisco_pr)   # Monthly precipitation raster [mm/month]
#' data(pisco_pe)   # Monthly potential evapotranspiration raster [mm/month]
#' data(qobs)       # Observed streamflow [m³/s]
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
#'   Save = FALSE
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
                                  IniGrids = '01/1981',
                                  Save = FALSE,
                                  Update = FALSE,
                                  Members = NULL) {

  tic()

  # ==== Input validation ====
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

  # ==== Format validation ====
  date_format_check <- function(date_string) {
    grepl("^\\d{2}/\\d{4}$", date_string)
  }
  if (!date_format_check(IniGrids)) {
    stop("Argument 'IniGrids' must be in 'mm/yyyy' format (e.g., '01/1981').")
  }
  if (!date_format_check(DateIni)) {
    stop("Argument 'DateIni' must be in 'mm/yyyy' format (e.g., '01/1990').")
  }
  if (!date_format_check(DateEnd)) {
    stop("Argument 'DateEnd' must be in 'mm/yyyy' format (e.g., '12/2020').")
  }

  # ==== Extract subbasin IDs ====
  comid <- as.vector(Subbasins$COMID)
  nsub  <- length(comid)

  # ==== Helper function: mean areal extraction ====
  extract_variable <- function(raster, variable_name) {
    cat('\f')
    cat(paste0("Computing monthly ", variable_name, " for each subbasin\n"))
    cat("Please wait...\n")

    # Convert SpatVector to sf for exactextractr
    Subbasins_sf <- sf::st_as_sf(Subbasins)

    if (is.null(Members)) {
      all_dates <- seq(
        from = as.Date(paste0("01/", IniGrids), "%d/%m/%Y"),
        to   = as.Date(paste0("01/", IniGrids), "%d/%m/%Y") + months(nlyr(raster) - 1),
        by   = "month"
      )
      Ini  <- as.Date(paste0("01/", DateIni), "%d/%m/%Y")
      End  <- as.Date(paste0("01/", DateEnd), "%d/%m/%Y")
      ind  <- which(all_dates >= Ini & all_dates <= End)

      extracted <- t(exact_extract(raster[[ind]], Subbasins_sf, "mean", progress = FALSE))

      if (Update) {
        extracted <- as.data.frame(extracted[nrow(extracted), , drop = FALSE])
      }
    } else {
      extracted <- t(exact_extract(raster, Subbasins_sf, "mean", progress = FALSE))
    }

    return(extracted)
  }

  # ==== Extract data ====
  if (!is.null(Precip)) {
    Prec <- extract_variable(Precip, "precipitation")
  }

  if (!is.null(PotEvap)) {
    Evap <- extract_variable(PotEvap, "evapotranspiration")
  }

  # ==== Build date vector ====
  if (Update) {
    DatesMonths <- as.Date(paste0("01/", DateEnd), "%d/%m/%Y")
  } else {
    Ini <- as.Date(paste0("01/", DateIni), "%d/%m/%Y")
    End <- as.Date(paste0("01/", DateEnd), "%d/%m/%Y")
    DatesMonths <- seq(Ini, End, by = "month")
    if (!is.null(Members)) {
      DatesMonths <- rep(DatesMonths, times = Members)
    }
  }

  # ==== Assemble output dataframe ====
  data_list <- list(DatesR = DatesMonths)

  if (!is.null(Precip)) {
    Prec <- round(Prec, 1)
    colnames(Prec) <- paste0("P_", comid)
    data_list <- c(data_list, as.data.frame(Prec))
  }

  if (!is.null(PotEvap)) {
    Evap <- round(Evap, 1)
    colnames(Evap) <- paste0("E_", comid)
    data_list <- c(data_list, as.data.frame(Evap))
  }

  if (!is.null(Qobs)) {
    Flow <- round(Qobs, 3)
    colnames(Q) <- 'Q'
    data_list <- c(data_list, as.data.frame(Flow))
  }

  Ans <- as.data.frame(data_list)

  # ==== Optional: Save to file ====
  if (Save) {
    if (!dir.exists("./Inputs")) dir.create("./Inputs")
    write.table(Ans, file = "./Inputs/Inputs_model.txt", row.names = FALSE)
  }

  cat("Processing completed successfully in... ")
  toc()
  return(Ans)
}
