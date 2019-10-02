#' Create raster inputs (.TIF) to run the Semidistribute GR2M version.
#'
#' @param Location  General work directory where data is located.
#' @param Shapefile Subbasins shapefile.
#' @param Dem Raster DEM for the study area.
#' @return Export centroids mask and flow direction rasters for the study subbasins.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
Create_Raster_Inputs <- function(Location, Shapefile, Dem){

  # Load packages
    require(rgdal)
    require(raster)
    require(rgeos)

  # Show message
    cat('\f')
    message("Creating Centroids Mask and Flow Direction rasters")
    message("Please wait...")

  # Load subbasin shapefile and raster DEM
    area  <- readOGR(file.path(getwd(), 'Inputs', Shapefile), verbose=FALSE)
    nsub  <- nrow(area@data)
    num   <- raster(file.path(getwd(),'Inputs',Dem))

  # Create a raster mask
    qmask         <- num
    values(num)   <- 1:ncell(num)
    num           <- mask(num, area)
    values(qmask) <- 0
    Mask          <- stack(qmask)

  # Extract subbasin rasters centroids
    for (j in 1:nsub){
      xy  <- coordinates(gCentroid(area[j,]))
      val <- extract(num, xy, method='simple')
      if(is.na(val)==TRUE){
        val <- ceiling(median(values(mask(num, area[j,])), na.rm=T))
      }
      qmask[num==val] <- j
      Mask            <- addLayer(Mask, qmask)
      values(qmask)   <- 0
    }
    QMask <- stackApply(Mask, rep(1,nsub+1), fun=sum)

  # Save raster of mask
    writeRaster(QMask, filename=file.path(getwd(), 'Inputs', 'Centroids_mask.tif'), overwrite=T)

  # Pitremove
    setwd(file.path(Location,'Inputs'))
    system(paste0("mpiexec -n 8 pitremove -z ",File.Raster," -fel Ras.tif"))

  # D8 flow directions
    system("mpiexec -n 8 D8Flowdir -p Flow_Direction.tif -fel Ras.tif",show.output.on.console=F,invisible=F)
    file.remove('Ras.tif')
    file.remove('Ši.tif')

    message('Done!')
}
