#' Create a TXT file with forcing data inputs (DatesR, Precipitation, Potential Evapotranspiration and Observed streamflow).
#'
#' @param Shapefile Subbasins shapefile.
#' @param Database Path where precipitacion and Pot. evapotranspiration files are located.
#' @param Precip Precipitation filename (in NetCDF format).
#' @param PotEvap Potential evapotranspiration filename (in NetCDF format).
#' @param Qobs Observed streamflow data (in m3/s) for the study period.
#' @param DateIni Initial subset date in 'yyyy/mm/dd' format. '1981/01/01' as default
#' @param DateEnd Final subset date in 'yyyy/mm/dd' format. '2019/02/01' as default
#' @return Export a text file with forcing data inputs to run semidistribute GR2M model.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  tictoc
#' @import  ncdf4
Create_Forcing_Inputs <- function(Shapefile, Database, Precip, PotEvap, Qobs, DateIni='1981/01/01', DateEnd='2019/02/01'){

    # Load packages
      require(rgdal)
      require(raster)
      require(rgeos)
      require(tictoc)
      require(ncdf4)
      tic()

    # Create a vector of dates
      DatesMonths <- seq(as.Date(DateIni), as.Date(DateEnd), by='month')

    # Load subbasins shapefiles
      Basins  <- readOGR(file.path(getwd(),'Inputs', Shapefile))
      nBasins <- nrow(Basins@data)
      proj    <- crs(Basins)

    # Extract monthly precipitation data for each subbasin
    # Because of subbasins with less area than PISCOO pixel, we obtain the time series
    # for the centroid of the subbasin
    # Read precipitation data
      pp      <- brick(file.path(Database, Precip))
      nDataPP <- nlayers(pp)

    # Extract data for each subbasin
      DataPP <- matrix(NA, ncol=nBasins, nrow=length(DatesMonths))
      for (i in 1:nBasins){
        pp.mask <- mask(pp, Basins[i,])
        pp.mean <- as.vector(cellStats(pp.mask, "mean"))
        if (is.na(pp.mean[1]) == TRUE){
          xycen         <- coordinates(gCentroid(Basins[i,]))
          xy.point      <- SpatialPoints(matrix(as.numeric(xycen), ncol=2, nrow=1))
          crs(xy.point) <- crs(Basins)
          pp.out        <- round(as.vector(extract(data_pp, xy.point, method='bilinear')),1)
        }else{
          pp.out        <- round(pp.mean,1)
        }
        cat('\f')
        message('Processing monthly precipitation [mm]')
        message('=====================================')
        message(paste0('Sub-basin ',i,' of ',nBasins))
        message('Please wait...')
        DataPP[,i] <- pp.out
      }
      DataPP[DataPP<0] <- 0

    # Extract monthly potential evapotranspiration for each subbasin
    # In this case we build the climatological series for each subbasin. Because of subbasins with
    # less area than PISCOO pixel, we obtain the time series for the centroid of the subbasin
    # Read potential evapotranspiration data
    pet       <- brick(file.path(Database, PotEvap))
    nDataPET  <- nlayers(pet)
    index     <- rep(1:12, nDataPET/12)
    clim.pet  <- stackApply(pet, index, fun=mean, na.rm=T)
    addYears  <- as.numeric(format(as.Date(DateEnd), '%Y'))-1-2016
    addMonths <- as.numeric(format(as.Date(DateEnd), '%m'))

    # Extract data for each subbasin
      DataPET <- matrix(NA, ncol=nBasins, nrow=length(DatesMonths))
      for (i in 1:nBasins){
        clim.mask <- mask(clim.pet, Basins[i,])
        clim.mean <- as.vector(cellStats(clim.mask, "mean", na.rm=T))
        if (is.na(clim.mean[1]) == TRUE){
          xycen         <- coordinates(gCentroid(Basins[i,]))
          xy.point      <- SpatialPoints(matrix(as.numeric(xycen), ncol=2, nrow=1))
          crs(xy.point) <- proj
          pet.value     <- round(as.vector(extract(clim.pet, xy.point, method='bilinear')),2)
          DataPET[,i]   <- round(c(rep(pet.value, (nDataPET/12+addYears)), pet.value[1:addMonths]), 2)
        }else{
          DataPET[,i]   <- round(c(rep(clim.mean,(nDataPET/12+addYears)), clim.mean[1:addMonths]), 2)
        }
        cat('\f')
        message('Processing monthly potential evap. [mm]')
        message('=======================================')
        message(paste0('Sub-basin ',i,' of ',nBasins))
        message('Please wait...')
      }
      DataPET[DataPET<0] <- 0

    # Save dataframe as .csv (by commas)
      Qobs <- read.table(Qobs, sep='\t', header=F)
      df   <- data.frame(DatesMonths, DataPP, DataPET, Qobs)
      colnames(df) <- c('DatesR', paste0('P',1:nBasins), paste0('E',1:nBasins), 'Qm3s')
      write.table(df, file=file.path(getwd(),'Inputs','Inputs_Basins.txt'), sep='\t', col.names=TRUE, row.names=FALSE)
      toc()
      message('Done!')
}# End (not run)
