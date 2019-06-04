#' Create raster inputs (.TIF) to run the Semidistribute GR2M version.
#'
#' @param Shapefile Subbasins shapefile.
#' @param Dem Raster DEM for the study area.
#' @return Export centroids mask and flow direction rasters for the study subbasins.
#' @export0
Create_Raster_Inputs <- function(Shapefile, Dem){

  # Load packages
    require(rgrass7)
    require(rgdal)
    require(raster)
    require(rgeos)

  # Show message
    cat('\f')
    message("Creating Mask and Flow Direction rasters")
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

  # Load GRASS
    loc <- initGRASS('C:/Program Files/GRASS GIS 7.4.4', home=getwd(),
                     gisDbase="GRASS_TEMP", override=TRUE)

  # Importar raster DEM
    execGRASS("r.in.gdal", Sys_show.output.on.console=F, flags=c('o','overwrite'), parameters=list(input=file.path(getwd(),'Inputs',Dem), output="dem"))

  # Setup extention
    execGRASS("g.region", Sys_show.output.on.console=F, parameters=list(raster="dem"))

  # Create flow direction D8
    execGRASS("r.watershed", Sys_show.output.on.console=F, flags=c("overwrite", "s", "a"), parameters=list(elevation="dem", drainage='fdr'))

  # Save raster of flow direction
    execGRASS("r.out.gdal", Sys_show.output.on.console=F, flags='overwrite',
              parameters=list(input='fdr', output=file.path(getwd(),'Inputs','Flow_Direction.tif'), format='GTiff'))

  # Clean GRASS workspace
    unlink(file.path(getwd(), "GRASS_TEMP"), recursive=T)

    message('Done!')
}
