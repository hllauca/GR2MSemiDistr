# GR2MSemiDistr

**GR2MSemiDistr** is an R package for implementing a **semi-distributed monthly water balance model (GR2M)** with spatial discretization into subbasins and flow routing based on a **transfer matrix interpreted as a directed graph**.  

The package is designed for **large-sample and national-scale hydrological applications**, but it can also be applied to **small basins or local studies** with just a few subbasins. It supports both **single-region (shared parameters)** and **multi-region calibration** with correction factors for precipitation and potential evapotranspiration.  

**GR2MSemiDistr** enables **long-term simulations**, **operational updating**, and robust routing schemes that ensure upstream–downstream flow accumulation in river networks ranging from small catchments to large-scale systems (>12,000 subbasins).  

It integrates seamlessly with the `airGR` hydrological modeling framework and the `terra` spatial ecosystem, making it suitable for climate services, water resource assessments, and operational hydrological forecasting.  

The package was developed as part of hydrological research and operational water resource modeling in Peru (see Llauca et al., 2021). Planned extensions include multi-objective calibration methods, automated input preparation from observations, bias correction of gridded forcings, and enhanced routing with linear reservoir schemes.  

> For any issues, bug reports, or suggestions, please contact: **Harold Llauca** (hllauca@senamhi.gob.pe)

---

## Features

- Direct integration with gridded precipitation and evapotranspiration datasets for spatially explicit subbasin inputs.  
- Semi-distributed GR2M monthly water balance modeling over thousands of subbasins.  
- Flexible parameterization: single-region setup (shared parameters across all subbasins) or multi-region calibration with regionalized correction factors for precipitation and potential evapotranspiration.  
- Long-term and operational streamflow simulation for hydrological forecasting and water resources assessment.  
- Routing of subbasin discharges using a transfer matrix interpreted as a directed graph, ensuring correct upstream–downstream accumulation and scalability to large river networks (>12,000 subbasins).  
- Seamless compatibility with the `airGR` hydrological modeling framework and the `terra` spatial ecosystem.  
- Scalable and efficient workflows designed for national-scale hydrological modeling and climate services.  


## Planned developments

- Automated delineation of subbasins and river reaches from DEMs, with connectivity information used to build the transfer matrix.  
- Automated preparation of model input datasets using observed hydro-meteorological records.  
- Bias correction of gridded precipitation and evapotranspiration products based on in-situ observations.  
- Integration of additional calibration methods, including multi-objective optimization algorithms.  
- Enhanced flow routing using a linear reservoir approach with estimable or calibrable routing parameters.  


---

## Installation

```r
# Install devtools if not already installed
install.packages("devtools")

# Install GR2MSemiDistr from GitHub
devtools::install_github("hllauca/GR2MSemiDistr")

# Load the package
library(GR2MSemiDistr)
```

## Example

---

```r
# This example script illustrates the complete workflow for setting up,
# calibrating, and running a semi-distributed GR2M water balance model
# using the GR2MSemiDistr package in R. The workflow includes:
#   - loading example input data (subbasin boundaries, DEM, observed flows, and gridded forcings),
#   - preparing model inputs in the airGR format,
#   - calibrating model parameters with the SCE-UA algorithm,
#   - running the semi-distributed GR2M simulation with routing,
#   - and visualizing hydrological components and outlet performance.
#
# In this simplified case, all subbasins are assumed to belong to a single
# calibration region. Consequently, the same set of optimized parameters
# (X1, X2, fp, fe) is applied uniformly to every subbasin, while precipitation
# and potential evapotranspiration inputs remain spatially distributed.
#
library(GR2MSemiDistr)
library(terra)

# === Load example data ===
data <- GR2MSemiDistr::Load_example_data()
names(data)

cat     <- data$cat       # Subbasin boundaries (SpatVector)
dem     <- data$dem       # DEM (SpatRaster)
qobs    <- data$qobs      # Observed streamflow (data.frame)
grid_pr <- data$grid_pr   # Precipitation (SpatRaster)
grid_pe <- data$grid_pe   # Potential evapotranspiration (SpatRaster)
matrixT <- data$matrixT   # Connectivity matrix (data.frame)

# === Visualize subbasins on DEM ===
plot(mask(dem, cat), alpha = 0.8, main = "Subbasins on DEM")
plot(cat, add = TRUE)
text(cat, labels = cat$COMID, cex = 0.6)

# === Define simulation period ===
start_date <- "01/1981"
end_date   <- "12/2016"
ini_grids  <- "01/1981"

# === Prepare model forcing inputs ===
model_inputs <- Create_Forcing_Inputs(
  Subbasins = cat,
  Precip    = grid_pr,
  PotEvap   = grid_pe,
  Qobs      = qobs,
  DateIni   = start_date,
  DateEnd   = end_date,
  IniGrids  = ini_grids
)

# === Visualize input data ===

# Select a target COMID to plot
target_comid <- "10"

# Plot precipitation for the target subbasin
P_col <- paste0("P_", target_comid)
plot(model_inputs$DatesR, model_inputs[[P_col]], type = "h", lwd = 2, col = "dodgerblue",
     ylab = "Precipitation [mm/month]", xlab = "Date",
     main = sprintf("Precipitation in Subbasin %s", target_comid))

# Plot potential evapotranspiration for the target subbasin
E_col <- paste0("E_", target_comid)
plot(model_inputs$DatesR, model_inputs[[E_col]], type = "o", lwd = 2, col = "darkgreen",
     ylab = "PET [mm/month]", xlab = "Date",
     main = sprintf("Potential Evapotranspiration in Subbasin %s", target_comid))

# Plot observed streamflow at outlet (if available)
if ("Q" %in% names(model_inputs)) {
  plot(model_inputs$DatesR, model_inputs$Q, type = "l", col = "darkblue", lwd = 1.5,
       ylab = "Q [m³/s]", xlab = "Date",
       main = "Observed Streamflow at Basin Outlet")
}

# === Initialize model parameters by region ===
param_init <- data.frame(
  Region = unique(cat$Region),
  X1  = 500,  # Production store capacity [mm]
  X2  = 1.5,  # Groundwater exchange coefficient
  fp  = 1.0,  # Precipitation correction factor
  fe  = 1.0   # Evapotranspiration correction factor
)

# === Calibrate parameters using OF10 as objective function ===
result <- Optim_GR2MSemiDistr2(
  Data          = model_inputs,
  Subbasins     = cat,
  RunIni        = "01/1981",
  RunEnd        = "12/2005",
  TransferMatrix= matrixT,
  Outlet        = '8', # basin outlet COMID
  Parameters    = param_init,
  Optimization  = "OF10"
)

# Extract optimized parameters and performance
best_params <- result$Parameters
final_score <- result$OF

# === Run GR2M model using optimized parameters ===
model <- Run_GR2MSemiDistr2(
  Data          = model_inputs,
  Subbasins     = cat,
  RunIni        = "01/1981",
  RunEnd        = "12/2016",
  TransferMatrix= matrixT,
  Outlet        = '8',
  Parameters    = best_params
)

# === Plot model components for a selected subbasin ===
target_comid <- "10"
idx <- which(model$COMID == target_comid)

par(mfrow = c(3, 2),
    mar = c(4, 4, 3, 1),
    cex.main = 0.9,
    cex.lab = 0.8,
    cex.axis = 0.7)

# Precipitation
plot(model$Dates, model$PR[, idx], type = "h", lwd = 2, col = "dodgerblue",
     main = sprintf("Precipitation (PR) - COMID %s", target_comid),
     xlab = "Date", ylab = "mm/month")

# Actual evapotranspiration
plot(model$Dates, model$AE[, idx], type = "h", lwd = 2, col = "orange",
     main = sprintf("Actual Evapotranspiration (AE) - COMID %s", target_comid),
     xlab = "Date", ylab = "mm/month")

# Soil moisture storage
plot(model$Dates, model$SM[, idx], type = "h", lwd = 2, col = "green",
     main = sprintf("Soil Moisture (SM) - COMID %s", target_comid),
     xlab = "Date", ylab = "mm")

# Percolation
plot(model$Dates, model$PC[, idx], type = "h", lwd = 2, col = "purple",
     main = sprintf("Percolation (PC) - COMID %s", target_comid),
     xlab = "Date", ylab = "mm/month")

# Runoff
plot(model$Dates, model$RU[, idx], type = "h", lwd = 2, col = "red",
     main = sprintf("Runoff (RU) - COMID %s", target_comid),
     xlab = "Date", ylab = "mm/month")

# Routed discharge
plot(model$Dates, model$QR[, idx], type = "l", lwd = 1.5, col = "blue",
     main = sprintf("Routed Simulated Discharge (QR) - COMID %s", target_comid),
     xlab = "Date", ylab = "m³/s")

# === Compare simulated vs observed discharge at outlet ===
par(mfrow = c(1, 1))
plot(model$Dates, model$SINK$sim, type = "l", col = "blue", lwd = 1.5,
     xlab = "Date", ylab = "Discharge [m³/s]",
     main = "Simulated vs Observed Discharge (Outlet)")
lines(model$Dates, model$SINK$obs, col = "red", lty = 2, lwd = 1.5)
legend("topright", legend = c("Simulated", "Observed"),
       col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")
```
