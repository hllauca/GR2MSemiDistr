# GR2MSemiDistr

**GR2MSemiDistr** is an R package for implementing a **semi-distributed monthly water balance model (GR2M)** with spatial subbasin discretization and routing using a **Weighted Flow Accumulation (WFA)** approach. Designed for **large-sample hydrological applications**, it supports spatial parameter regionalization, operational updating, and robust routing based on D8 flow directions.

This package was developed as part of hydrological research and operational water resource modeling in Peru, as described in Llauca et al. (2021).

> For any issues, bug reports, or suggestions, please contact: **Harold Llauca** (hllauca@senamhi.gob.pe)

---

## Features

- GR2M water balance modeling in a semi-distributed configuration  
- Regional parameter definition and correction factors  
- Support for long-term and operational streamflow simulation  
- Integration with gridded PISCO precipitation and evapotranspiration datasets  
- Routing of discharges via WFA algorithm using DEM and D8 flow directions  
- Fully compatible with the `airGR` and `terra` ecosystems  

---

## Installation

```r
# Install devtools if not already installed
install.packages("devtools")

# Install GR2MSemiDistr from GitHub
devtools::install_github("hllauca/GR2MSemiDistr")

# Load the package
library(GR2MSemiDistr)
