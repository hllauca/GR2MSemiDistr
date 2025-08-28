#' Build transfer matrix (rows = receivers, cols = donors)
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
#' @param riv A `SpatVector` of river reaches/subbasins. Its attribute table must contain
#'   `COMID` (unique subbasin ID) and `NextDownID` (immediate downstream ID).
#'   Geometry is ignored; only attributes are used. Outlets (no downstream basin/river)
#'   must be encoded with `-1` in `NextDownID`.
#'
#' @return A `data.frame` of size n × n with entries in {0,1}. Row names and
#'   column names are the node IDs (COMID domain). Outlets (no downstream) lead
#'   to rows of zeros.
#'
#' @import terra
#' @import Matrix
#'
#' @export
Build_Transfer_Matrix <- function(riv) {
  # ==== Input validation ====
  if (!inherits(riv, "SpatVector")) {
    stop("'riv' must be a 'SpatVector'. If you have a file path, read it first with terra::vect().")
  }

  # Attributes only
  att <- as.data.frame(riv)
  req_cols <- c("COMID", "NextDownID")
  if (!all(req_cols %in% names(att))) {
    stop("Input must contain columns: 'COMID' and 'NextDownID'.")
  }

  # ==== Basic preparation ====
  df <- att[, req_cols]
  # Coerce to character to avoid numeric/character mismatches
  df$COMID      <- as.character(df$COMID)
  df$NextDownID <- as.character(df$NextDownID)

  # Drop completely missing COMID rows (defensive)
  df <- df[!is.na(df$COMID) & nzchar(df$COMID), , drop = FALSE]
  if (nrow(df) == 0L) stop("No valid COMID rows found after cleaning.")

  # ==== Filter valid donor -> receiver pairs ====
  # Treat -1 as "no downstream" (outlets)
  is_outlet <- df$NextDownID == "-1"
  donors    <- df$COMID[!is_outlet]
  receivers <- df$NextDownID[!is_outlet]

  # If there are no links at all, still return an identity-sized empty matrix over COMIDs
  # (useful for degenerate networks with only sources/outlets)
  if (length(donors) == 0L) {
    ids <- sort(unique(df$COMID))
    n   <- length(ids)
    if (n == 0L) stop("No node IDs available to build the transfer matrix.")
    # dense 0 matrix as data.frame
    out <- as.data.frame(matrix(0L, nrow = n, ncol = n, dimnames = list(ids, ids)))
    return(out)
  }

  # ==== Node set and deterministic order ====
  # Use all donor COMIDs + all valid receivers; sorted for stable reproducibility
  ids <- sort(unique(c(df$COMID, receivers)))
  n   <- length(ids)

  # ==== Map to integer indices (rows = receivers, cols = donors) ====
  row_idx <- match(receivers, ids)  # receivers -> rows
  col_idx <- match(donors,    ids)  # donors    -> cols

  # Sanity check: indices must be within [1, n]
  if (anyNA(row_idx) || anyNA(col_idx)) {
    stop("Internal indexing failed (NA indices). Check 'COMID'/'NextDownID' consistency.")
  }

  # ==== Build sparse matrix; collapse duplicates to 1 ====
  Tsp <- Matrix::sparseMatrix(
    i = row_idx,
    j = col_idx,
    x = 1L,
    dims = c(n, n),
    dimnames = list(ids, ids)
  )
  if (length(Tsp@x)) Tsp@x[] <- 1L  # cap duplicates to 1

  # ==== Convert to dense data.frame (0/1) ====
  Tfull <- as.matrix(Tsp)
  storage.mode(Tfull) <- "integer"

  out <- as.data.frame(Tfull, stringsAsFactors = FALSE)

  # ==== Final structural checks (for downstream compatibility) ====
  # Must be square, named rows/cols, and only contain 0/1
  if (!isTRUE(all.equal(nrow(out), ncol(out)))) {
    stop("Transfer matrix must be square.")
  }
  if (is.null(rownames(out)) || is.null(colnames(out))) {
    stop("Transfer matrix must carry row and column names (COMID domain).")
  }
  return(out)
}
