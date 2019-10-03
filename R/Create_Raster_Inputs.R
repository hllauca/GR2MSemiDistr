#' Create raster inputs for routing monthly streamflow. Require TauDEM installed.
#'
#' @param Location  Work directory where 'Inputs' folder is located.
#' @param Shapefile Subbasins shapefile.
#' @param Dem Raster DEM.
#' @return Export centroids mask and flow direction rasters.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
Create_Raster_Inputs <- function(Location, Shapefile, Dem){
# Shapefile <- File.Shape
# Dem       <- File.Raster

  # Load packages
    require(rgdal)
    require(raster)
    require(rgeos)

  # Show message
    cat('\f')
    message("Creating Centroids Mask raster")
    message("Please wait...")

  # Load subbasin shapefile and raster DEM
    area  <- readOGR(file.path(Location, 'Inputs', Shapefile), verbose=FALSE)
    dem   <- raster(file.path(Location,'Inputs',Dem))

  # Create a raster mask
    mask  <- dem
    values(mask) <- 1:ncell(dem)

    # Auxiliary function (from https://stackoverflow.com/questions/44327994/calculate-centroid-within-inside-a-spatialpolygon)
    gCentroidWithin <- function(pol) {
      require(rgeos)

      pol$.tmpID <- 1:length(pol)
      # initially create centroid points with gCentroid
      initialCents <- gCentroid(pol, byid = T)

      # add data of the polygons to the centroids
      centsDF <- SpatialPointsDataFrame(initialCents, pol@data)
      centsDF$isCentroid <- TRUE

      # check whether the centroids are actually INSIDE their polygon
      centsInOwnPoly <- sapply(1:length(pol), function(x) {
        gIntersects(pol[x,], centsDF[x, ])
      })

      if(all(centsInOwnPoly) == TRUE){
        return(centsDF)
      } else {
      # substitue outside centroids with points INSIDE the polygon
      newPoints <- SpatialPointsDataFrame(gPointOnSurface(pol[!centsInOwnPoly, ],
                                          byid = T), pol@data[!centsInOwnPoly,])
      newPoints$isCentroid <- FALSE
      centsDF <- rbind(centsDF[centsInOwnPoly,], newPoints)

      # order the points like their polygon counterpart based on `.tmpID`
      centsDF <- centsDF[order(centsDF$.tmpID),]

      # remove `.tmpID` column
      centsDF@data <- centsDF@data[, - which(names(centsDF@data) == ".tmpID")]

      cat(paste(length(pol), "polygons;", sum(centsInOwnPoly), "actual centroids;",
          sum(!centsInOwnPoly), "Points corrected \n"))

      return(centsDF)
     }}

    xycen  <- gCentroidWithin(area)
    val    <- extract(mask, xycen, method='simple')
    mask[mask %in% val] <- 1
    mask[mask!=1] <- 0

  # Save raster
    writeRaster(mask, filename=file.path(getwd(), 'Inputs', 'Centroids_mask.tif'), overwrite=T)


  # Show message
    message("Creating Flow Direction raster")
    message("Please wait...")

  # Create flow direction raster with TauDEM
  # Pitremove
    setwd(file.path(Location,'Inputs'))
    system(paste0("mpiexec -n 8 pitremove -z ",File.Raster," -fel Ras.tif"))

  # D8 flow directions
    system("mpiexec -n 8 D8Flowdir -p Flow_Direction.tif -fel Ras.tif",show.output.on.console=F,invisible=F)
    file.remove('Ras.tif')

  # Show message
  message('Done!')
}
