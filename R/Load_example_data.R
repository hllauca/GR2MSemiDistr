#' Load example hydrological data included in the package
#'
#' Loads spatial and tabular example data (raster and vector) from 'inst/extdata/'.
#'
#' @return A list with objects: `cat`, `dem`, `grid_pr`, `grid_pe`, `qobs`
#' @examples
#' # Load example data
#'
#' data <- GR2MSemiDistr::Load_example_data()
#' names(data)
#' cat     <- data$cat
#' dem     <- data$dem
#' qobs    <- data$qobs
#' grid_pr <- data$grid_pr
#' grid_pe <- data$grid_pe
#'
#' @export
Load_example_data <- function() {
  list(
    cat      = terra::vect(system.file("extdata/cat.gpkg", package = "GR2MSemiDistr")),
    dem      = terra::rast(system.file("extdata/dem.tif", package = "GR2MSemiDistr")),
    grid_pr  = terra::rast(system.file("extdata/grid_pr.tif", package = "GR2MSemiDistr")),
    grid_pe  = terra::rast(system.file("extdata/grid_pe.tif", package = "GR2MSemiDistr")),
    qobs     = read.table(system.file("extdata/qobs.txt", package = "GR2MSemiDistr"), header = FALSE)
  )
}


