#' Route a monthly inflow series through a single linear reservoir
#'
#' Internal helper used by `Run_GR2MSemiDistr()` and `Optim_GR2MSemiDistr()` to
#' transit the flow contributed by a donor subbasin through the reach that
#' connects it to its downstream receiver, instead of adding it instantly.
#' The reach is modeled as a linear reservoir (`S = k * Qout`), whose mass
#' balance `dS/dt = I(t) - Qout(t)` has the closed-form monthly recursion:
#'
#' `Qout(t) = alpha * Qout(t-1) + (1 - alpha) * I(t)`, with `alpha = exp(-dt / k)`.
#'
#' `k = 0` reproduces the previous instantaneous-transfer behavior exactly
#' (`alpha = exp(-Inf) = 0`, so `Qout(t) = I(t)`), so this is a strict
#' generalization, not a breaking change, when `k` is left at its default of 0.
#'
#' @param I numeric vector. Monthly inflow to the reach (the donor subbasin's
#' already-accumulated upstream discharge).
#' @param k numeric. Transit/storage constant of the reach, in months. `k = 0`
#' means instantaneous transfer (no lag, no attenuation); larger `k` produces
#' more lag and more attenuation (e.g. floodplain/backwater storage).
#' @param dt numeric. Timestep length, in the same units as `k` (months).
#' Fixed at 1 since the model operates at a monthly timestep.
#' @param Qout0 numeric. Initial outflow state (routing memory) before the
#' first value of `I`. Defaults to `I[1]` (steady-state assumption) when not
#' supplied, which avoids an artificial ramp-up at the start of the series.
#'
#' @return numeric vector, same length as `I`: the routed (transited) outflow.
#'
#' @importFrom stats filter
#' @noRd
route_linear_reservoir <- function(I, k, dt = 1, Qout0 = NULL) {
  if (is.na(k) || k < 0) stop("Routing constant 'k' must be a non-negative number.")
  if (is.null(Qout0)) Qout0 <- I[1]

  alpha <- exp(-dt / k)  # k = 0 -> dt/k = Inf -> alpha = 0 (instantaneous transfer)

  as.numeric(stats::filter((1 - alpha) * I, filter = alpha,
                           method = "recursive", init = Qout0))
}
