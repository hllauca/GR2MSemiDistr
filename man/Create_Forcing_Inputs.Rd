% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Create_Forcing_Inputs.R
\name{Create_Forcing_Inputs}
\alias{Create_Forcing_Inputs}
\title{Extract and prepare model's inputs data in the airGR format (DatesR, P and E) from gridded P and E monthly data.}
\usage{
Create_Forcing_Inputs(
  Subbasins,
  Precip,
  PotEvap,
  Qobs = NULL,
  DateIni,
  DateEnd,
  IniGrids = "01/1981",
  Save = FALSE,
  Update = FALSE,
  Resolution = 0.01,
  Buffer = 1.1,
  Members = NULL
)
}
\arguments{
\item{Subbasins}{Subbasins' shapefile. Must contain the following attributes: 'Area' (in km2), 'Region' (in letters), and 'COMID' (identifier number).}

\item{Precip}{Raster brick of the precipitation data in [mm/month].}

\item{PotEvap}{Raster brick of the evapotranspiration data in [mm/month].}

\item{Qobs}{Observed streamflow data in [m3/s] at the basin outlet. Must have the dates of the output dataset.}

\item{DateIni}{Initial date of the output database in 'mm/yyyy' format.}

\item{DateEnd}{Ending date of the output database in 'mm/yyyy' format.}

\item{IniGrids}{Initial date of the gridded data (P and PE) in 'mm/yyy' format.}

\item{Save}{Boolean to save results as a text file in the 'Outputs' location. FALSE as default.}

\item{Update}{Boolean for the updating mode where only the last month's values will be returned. FALSE as default.}

\item{Resolution}{Resampling resolution for improving subbasins' data extraction. 0.01degrees as default.}

\item{Buffer}{Factor for increase subbasins' limits extents. 1.1 as default.}

\item{Members}{Número de miembros del conjunto modelo. Only for streamflow forecasting purposes. NULL por defecto.}
}
\value{
Return a dataframe of model's inputs data in the airGR format (DatesR, P, E, Q).
}
\description{
Extract and prepare model's inputs data in the airGR format (DatesR, P and E) from gridded P and E monthly data.
}
\examples{
# Load data
require(GR2MSemiDistr)
data(pisco_pr)
data(pisco_pe)
data(qobs)
data(roi)

# Create a database with model's inputs data
data <- Create_Forcing_Inputs(Subbasins=roi,
                              Precip=pisco_pr,
                              PotEvap=pisco_pe,
                              Qobs=qobs,
                              DateIni='01/1981',
                              DateEnd='12/2016',
                              IniGrids='01/1981')
View(data)
}
\references{
Cesar Aybar, Carlos Fernández, Adrian Huerta, Waldo Lavado, Fiorella Vega & Oscar Felipe-Obando (2020) Construction of a high-resolution gridded rainfall dataset for Peru from 1981 to the present day, Hydrological Sciences Journal, 65:5, 770-785, DOI: 10.1080/02626667.2019.1649411

Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020). Water. 2021; 13(8):1048. https://doi.org/10.3390/w13081048
}
