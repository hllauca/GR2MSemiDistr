# GR2MSemiDistr

**GR2MSemiDistr** is an R package for implementing a **semi-distributed monthly water balance model (GR2M)** with spatial subbasin discretization and routing using a **Weighted Flow Accumulation (WFA)** approach. Designed for **large-sample hydrological applications**, it supports spatial parameter regionalization, operational updating, and robust routing based on D8 flow directions.

This package was developed as part of hydrological research and operational water resource modeling in Peru, as described in Llauca et al. (2021).

> For any issues, bug reports, or suggestions, please contact: **Harold Llauca** (hllauca@senamhi.gob.pe)

---

## Features

- Semi-distributed GR2M water balance modeling over thousands of subbasins.
- Regionalized parameterization with correction factors for precipitation and potential evapotranspiration.
- Long-term and operational streamflow simulation for hydrological forecasting and water resource assessments.
- Direct integration with gridded precipitation and evapotranspiration datasets for areal subbasin computations.
- Routing of subbasin discharges via a transfer matrix built from subbasin connectivity information.
- Seamless compatibility with the airGR hydrological modeling framework and the terra spatial ecosystem.
- Scalable workflows designed for national-scale applications and climate services.

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
# This script demonstrates the workflow for setting up, running, and evaluating a hydrological simulation
# using the GR2MSemiDistr package in R. It includes loading input data, preparing model inputs,
# parameter calibration, running the semi-distributed GR2M water balance model, and visualizing outputs.

library(GR2MSemiDistr)

# === Load example data ===
data <- GR2MSemiDistr::Load_example_data()
names(data)
cat     <- data$cat       # Subbasin boundaries (SpatVector)
dem     <- data$dem       # Digital Elevation Model (SpatRaster)
qobs    <- data$qobs      # Observed streamflow at outlet (data.frame)
grid_pr <- data$grid_pr   # Precipitation gridded data (SpatRaster)
grid_pe <- data$grid_pe   # Potential evapotranspiration gridded data (SpatRaster)
matrixT <- data$matrixT   # Connectivity matrix (data.frame)

# === Visualize subbasins on DEM ===
terra::plot(terra::mask(dem, cat), alpha=0.8)
terra::plot(cat, add=TRUE)
terra::text(cat)

# === Define simulation period ===
start_date <- "01/1981"
end_date   <- "12/2016"
ini_grids  <- "01/1981"

# === Prepare model forcing inputs ===
model_inputs <- Create_Forcing_Inputs(
  Subbasins = cat,
  Precip = grid_pr,
  PotEvap = grid_pe,
  Qobs = qobs,
  DateIni = start_date,
  DateEnd = end_date,
  IniGrids = ini_grids,
)

# === Visualize input data ===
dates <- as.Date(model_inputs$DatesR)

# Plot mean precipitation across subbasins
P_vars <- grep("^P_", names(model_inputs), value = TRUE)
avg_P <- rowMeans(model_inputs[, P_vars], na.rm = TRUE)
plot(dates, avg_P, type = "l", col = "blue", lwd = 1.5,
     ylab = "Precipitation [mm/month]", xlab = "Date",
     main = "Mean Precipitation across Subbasins")

# Plot mean PET across subbasins
E_vars <- grep("^E_", names(model_inputs), value = TRUE)
avg_E <- rowMeans(model_inputs[, E_vars], na.rm = TRUE)
plot(dates, avg_E, type = "l", col = "orange", lwd = 1.5,
     ylab = "PET [mm/month]", xlab = "Date",
     main = "Mean Potential Evapotranspiration across Subbasins")

# Plot observed streamflow (if available)
if ("Q" %in% names(model_inputs)) {
  plot(dates, model_inputs$Q, type = "l", col = "darkgreen", lwd = 1.5,
       ylab = "Q [m³/s]", xlab = "Date",
       main = "Observed Streamflow at Basin Outlet")
}

# === Initialize model parameters by region ===
param_init <- data.frame(
  Region = unique(cat$Region),
  X1 = 500,   # Production store capacity [mm]
  X2 = 1.5,   # Groundwater exchange coefficient
  fp = 1.0,   # Precipitation correction factor
  fe = 1.0    # Evapotranspiration correction factor
)

# === Calibrate parameters using OF10 as objective function ===
result <- Optim_GR2MSemiDistr(
  Data = model_inputs,
  Subbasins = cat,
  RunIni = "01/1981",
  RunEnd = "12/2005",
  MatrixTransfer=matrixT,
  Outlet = '8', # basin outlet comid
  Parameters = param_init,
  Optimization = "OF10"
)

# === Extract optimized parameters and performance ===
best_params <- result$Parameters
final_score <- result$OF

# === Run GR2M model using optimized parameters ===
model <- Run_GR2MSemiDistr(
  Data = model_inputs,
  Subbasins = cat,
  RunIni = "01/1981",
  RunEnd = "12/2016",
  MatrixTransfer=matrixT,
  Outlet = '8', # basin outlet comid
  Parameters = best_params
)


# === Plot model components for a selected subbasin ===
target_comid <- model$COMID[8]
idx <- which(model$COMID == target_comid)
dates <- as.Date(model$Dates)

par(mfrow = c(3, 2), mar = c(4, 4, 3, 1))
plot(dates, model$PR[, idx], type = "l", col = "blue",
     main = paste("Precipitation (PR) - COMID", target_comid),
     xlab = "Date", ylab = "mm/month")
plot(dates, model$AE[, idx], type = "l", col = "orange",
     main = paste("Actual Evapotranspiration (AE) - COMID", target_comid),
     xlab = "Date", ylab = "mm/month")
plot(dates, model$SM[, idx], type = "l", col = "green4",
     main = paste("Production Store (SM) - COMID", target_comid),
     xlab = "Date", ylab = "mm")
plot(dates, model$PC[, idx], type = "l", col = "purple",
     main = paste("Percolation (PC) - COMID", target_comid),
     xlab = "Date", ylab = "mm/month")
plot(dates, model$RU[, idx], type = "l", col = "darkred",
     main = paste("Runoff (RU) - COMID", target_comid),
     xlab = "Date", ylab = "mm/month")
plot(dates, model$QR[, idx], type = "l", col = "blue",
     main = paste("Simulated Discharge (QR) - COMID", target_comid),
     xlab = "Date", ylab = "m³/s")

# === Compare simulated vs. observed discharge at outlet ===
if (!is.null(model$SINK)) {
  par(mfrow = c(1, 1))
  plot(dates, model$SINK$sim, type = "l", col = "blue",
       xlab = "Date", ylab = "Discharge [m³/s]",
       main = "Simulated vs Observed Discharge (Outlet)")
  lines(dates, model$SINK$obs, col = "red", type = 'o')
  legend("topright", legend = c("Simulated", "Observed"),
         col = c("blue", "red"), lty = c(1, 2), lwd = 2, bty = "n")
}
```
