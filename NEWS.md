# GR2MSemiDistr 5.4

* `Run_GR2MSemiDistr()` gains `K_vec`, an optional named vector (by `COMID`)
  to set the linear-reservoir routing constant `k` per subbasin directly,
  bypassing the per-Region `Parameters$k` lookup -- for setting `k` from an
  externally computed estimate (e.g. reach length and an assumed velocity, or
  an existing routing parameter such as a Muskingum `K` in seconds,
  `K_vec = K_seconds / (30 * 86400)`) without having to calibrate it.

# GR2MSemiDistr 5.3

## New features

* Routing between subbasins is no longer an instantaneous sum: each donor
  subbasin's discharge now transits through a linear reservoir before being
  added to its downstream receiver
  (`Qout(t) = alpha*Qout(t-1) + (1-alpha)*I(t)`), controlled by a new `k`
  (months) column in `Parameters`. `k = 0` (the default when the column is
  absent) reproduces the previous instantaneous behavior exactly, so this is
  fully backward-compatible.
* Routing memory can now carry over across operational `Update = TRUE` runs
  via the new `RouteStatesIni`/`RouteStatesEnd` arguments and return value,
  the same way GR2M's own states already do.
* `k` is now a calibrable parameter in `Optim_GR2MSemiDistr()`
  (`Parameters.Min`/`Parameters.Max` default to length 5: X1, X2, fp, fe, k).
* New `Cores` argument in both `Run_GR2MSemiDistr()` and
  `Optim_GR2MSemiDistr()` to distribute the (mutually independent, per-
  subbasin) GR2M runs over a `parallel::makeCluster()` PSOCK cluster.
  Default `Cores = 1` keeps the previous sequential behavior unchanged.

## Bug fixes

* `Run_GR2MSemiDistr()`: `Database` was not truncated after discarding the
  artificial "fake month" added in single-month `Update = TRUE` runs, which
  broke `SINK` whenever both `Outlet` and observed `Q` were used.
* `Optim_GR2MSemiDistr()`: the routing graph and its topological order were
  rebuilt on every objective-function evaluation instead of once before
  optimization starts.
* Added validation that `TransferMatrix` is acyclic (a valid DAG) before
  routing, in both `Run_GR2MSemiDistr()` and `Optim_GR2MSemiDistr()`.
* `Optim_GR2MSemiDistr()`: the composite objective functions (`OF8`-`OF10`)
  used raw RMSE (in m3/s) inside a weighted sum together with bounded
  `(1 - efficiency)` terms, letting RMSE dominate the objective regardless of
  the weights for basins with non-trivial discharge. They now use RSR (RMSE
  normalized by `sd(Qobs)`, Moriasi et al. 2007) instead; `OF4` still reports
  the raw RMSE, unchanged.
* `Optim_GR2MSemiDistr()`: an invalid `Optimization` value now raises a clear
  error instead of silently propagating `NA` into the SCE-UA optimizer.
* `Create_Forcing_Inputs()`: `Qobs` length is now validated against the
  number of simulated months (a mismatch could previously be silently
  recycled by `data.frame()`, misaligning observed flow against dates).
* `Create_Forcing_Inputs()`: `Members = TRUE` now validates the assumed
  member-blocked raster layer order against `terra::time()` when available,
  and warns explicitly when it cannot be verified.

## Documentation / housekeeping

* `Precip`/`PotEvap` in `Create_Forcing_Inputs()` now have explicit `NULL`
  defaults and are documented as optional, matching how the code already
  treated them.
* Fixed a stale `@return` in `Load_example_data()` (missing `matrixT`) and
  duplicated/contradictory `@param WarmUp` documentation in
  `Optim_GR2MSemiDistr()`.
* Fixed `CAPAS_DIR` in `Build_Transfer_Matrix.R`'s companion build script
  pointing at a nonexistent path.
* `NAMESPACE` now uses precise `importFrom()` per dependency instead of
  whole-package `import()` (except `Matrix`, which needs it for S4 method
  dispatch on base generics such as `t()`), removing several namespace
  name-collision warnings on load.
* Removed committed session/build artifacts (`.RData`, `.RDataTmp`,
  `.Rhistory`, `junk`, `sh.exe.stackdump`, `tar/*.tar.gz`) and excluded them
  going forward via `.gitignore`/`.Rbuildignore`.
