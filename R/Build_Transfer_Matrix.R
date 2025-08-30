#' Build transfer matrix for subbasins connectivity
#'
#' Create a dense 0/1 transfer matrix (as a data.frame) from a river reach/subbasin
#' layer that contains a unique ID (`COMID`) and its immediate downstream ID
#' (`NextDownID`). **Rows are receivers**, **columns are donors**. The result is a
#' square (nsub × nsub) matrix of subbasin connectivity whose row and column names
#' must match `COMID`.
#'
#' Designed to plug into routing steps used by `Run_GR2MSemiDistr()` and related
#' functions in the package.
#'
#' @param Rivers SpatVector. River reaches or subbasins. Must be a vector object from the `terra` package.
#' Its attribute table must include:
#' \describe{
#'   \item{COMID}{Unique subbasin/reach identifier (character or numeric, internally coerced to character).}
#'   \item{NextDownID}{Identifier of the immediate downstream subbasin/reach. Outlets (no downstream) must be coded with -1.}
#' }
#' Geometry is ignored; only the attribute table is used to build the connectivity.
#'
#' @return dgCMatrix (sparse matrix). Subbasin connectivity matrix where
#' rows represent receiving (downstream) subbasins and columns represent donor (upstream) subbasins.
#' Row and column names are set to `COMID`.
#'
#' @import terra
#' @import Matrix
#'
#' @export
Build_Transfer_Matrix <- function(Rivers) {
  if (!inherits(Rivers, "SpatVector")) {
    stop("'Rivers' must be a 'SpatVector'. Read it first with terra::vect().")
  }

  att <- as.data.frame(Rivers)
  if (!all(c("COMID", "NextDownID") %in% names(att))) {
    stop("Input must contain 'COMID' and 'NextDownID'.")
  }

  df <- att[, c("COMID", "NextDownID")]
  df$COMID      <- as.character(df$COMID)
  df$NextDownID <- as.character(df$NextDownID)

  df <- df[!is.na(df$COMID) & nzchar(df$COMID), , drop = FALSE]
  if (nrow(df) == 0L) stop("No valid COMID rows found after cleaning.")

  # Filter outlets
  is_outlet <- df$NextDownID == "-1"
  donors    <- df$COMID[!is_outlet]
  receivers <- df$NextDownID[!is_outlet]

  ids <- sort(unique(c(df$COMID, receivers)))
  n   <- length(ids)

  if (n == 0L) stop("No node IDs available to build the transfer matrix.")

  row_idx <- match(receivers, ids)
  col_idx <- match(donors,    ids)

  if (anyNA(row_idx) || anyNA(col_idx)) {
    stop("Internal indexing failed. Check 'COMID'/'NextDownID' consistency.")
  }

  # Build sparse connectivity matrix
  Tsp <- Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = 1L,
    dims = c(n, n),
    dimnames = list(ids, ids)
  )

  # Cap duplicates to 1
  if (length(Tsp@x)) Tsp@x[] <- 1L

  return(Tsp)  # dgCMatrix
}
