#' Create a TXT file with forcing data inputs (DatesR, Precipitation, Potential Evapotranspiration and Observed streamflow).
#'
#' @param Shapefile Subbasins shapefile.
#' @param Database Path where precipitacion and Pot. evapotranspiration files are located.
#' @param Precip Precipitation filename (in NetCDF format).
#' @param PotEvap Potential evapotranspiration filename (in NetCDF format).
#' @param Qobs Observed streamflow data (in m3/s) for the study period.
#' @param Resolution Raster resolution to resample forcing data and extract areal mean values. 0.01 as default.
#' @param DateIni Initial subset date in 'yyyy/mm/dd' format. '1981/01/01' as default
#' @param DateEnd Final subset date in 'yyyy/mm/dd' format. '2016/12/01' as default
#' @return Export a text file with forcing data inputs to run semidistribute GR2M model.
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  tictoc
#' @import  ncdf4
#' @import  parallel
Create_Forcing_Inputs <- function(Shapefile, Database, Precip, PotEvap, Qobs, Resolution=0.01, DateIni='1981/01/01', DateEnd='2016/12/01'){

Shapefile=File.Shape
Database=Database
Precip=File.Precip
PotEvap=File.PotEvap
Qobs=File.Qobs
Resolution=0.01
DateIni='1981/01/01'
DateEnd='2016/12/01'

    # Load packages
      require(rgdal)
      require(raster)
      require(rgeos)
      require(tictoc)
      require(ncdf4)
      require(parallel)
      tic()

    # Useful functions
    # Functions adopted from Cesar Aybar
      generate_mask_geom<-function(cov, geom){
        specialcov=cov
        specialcov[]=1:ncell(cov)
        Position_rowcol<-function(i){
          quad1<-unlist(raster::extract(specialcov,geom[i,],small=T))
        }
        position<-lapply(1:length(geom),Position_rowcol)
        return(position)
      }

      mask_Fast_Extract <- function(cov, positionP, fun=mean, na.rm=T){
        matrix_R <- t(as.matrix(cov))
        sapply(1:length(positionP), function(i){
          Value  <- matrix_R[positionP[[i]]]
          fun(Value, na.rm=na.rm)
        })
      }

    # Create a vector of dates
      DatesMonths <- seq(as.Date(DateIni), as.Date(DateEnd), by='month')

    # Load subbasins shapefiles
      Basins  <- readOGR(file.path(getwd(),'Inputs', Shapefile))
      nBasins <- length(Basins)

    # Extract monthly precipitation data for each subbasin
    # Show message
      cat('\f')
      message('Extracting monthly precipitation [mm]')
      message('Please wait...')

    # Read precipitation data
      pp      <- brick(file.path(Database, Precip))

    # Crop for basin domain
      pp.crop <- crop(pp, extent(Basins)*1.1)

    # Mask for resampling using 'ngb' method
      pp.res      <- raster(extent(pp.crop[[1]]))
      crs(pp.res) <- crs(pp.crop)
      res(pp.res) <- Resolution

    # Cells within each subbasin
      positionPP  <- generate_mask_geom(pp.res[[1]], Basins)

    # Mean areal rainfall for each subbasin
      cl=makeCluster(detectCores()-1) #detectar y asignar numero de cluster
      clusterEvalQ(cl,c(library(raster))) #cargar paquete para aplicar a cada nodo
      clusterExport(cl,varlist=c("pp.res","pp.crop","positionPP","mask_Fast_Extract"),envir=environment())
      # nlayers(pp.crop)
      mean.pp <- parLapply(cl, 1:nlayers(pp.crop), function(z) {
                     res <- resample(pp.crop[[z]], pp.res, method='ngb')
                     ans <- mask_Fast_Extract(cov=res, positionPP, fun=mean, na.rm=T)
                     return(ans)
                 })
      # stopCluster(cl) #Close the cluster
      mean.pp <- do.call(rbind, mean.pp)
      mean.pp <- round(mean.pp[1:length(DatesMonths),],1)


    # Extract monthly potential evapotranspiration for each subbasin
    # Show message
      cat('\f')
      message('Extracting monthly evapotranspiration [mm]')
      message('Please wait...')

    # Read potential evapotranspiration data
      pet      <- brick(file.path(Database, PotEvap))

    # Crop for basin domain
      pet.crop <- crop(pet, extent(Basins)*1.5)

    # Mask for resampling using 'ngb' method
      pet.res      <- raster(extent(pet.crop[[1]]))
      crs(pet.res) <- crs(pet.crop)
      res(pet.res) <- Resolution

    # Cells within each subbasin
      positionPET  <- generate_mask_geom(pet.res[[1]], Basins)

    # Mean areal evapotranspiration for each subbasin
      # cl=makeCluster(detectCores()-1) #detectar y asignar numero de cluster
      clusterEvalQ(cl,c(library(raster))) #cargar paquete para aplicar a cada nodo
      clusterExport(cl,varlist=c("pet.res","pet.crop","positionPET","mask_Fast_Extract"),envir=environment())

      mean.pet <- parLapply(cl, 1:nlayers(pet.crop), function(z) {
        res <- raster::resample(pet.crop[[z]], pet.res, method='ngb')
        ans <- mask_Fast_Extract(cov=res, positionPET, fun=mean, na.rm=T)
        return(ans)
      })
      stopCluster(cl) #Close the cluster
      mean.pet <- do.call(rbind, mean.pet)
      mean.pet <- round(mean.pet[1:length(DatesMonths),],1)

    # Subset streamflow data
      Obs     <- read.table(file.path("Inputs",Qobs), sep='\t', header=F)
      Ind_Obs <- seq(which(as.Date(Obs[,1]) == DateIni),
                     which(as.Date(Obs[,1]) == DateEnd))
      qobs    <- Obs[Ind_Obs,2]
      df      <- data.frame(DatesMonths, mean.pp, mean.pet, qobs)
      colnames(df) <- c('DatesR', paste0('P',1:nBasins), paste0('E',1:nBasins), 'Qm3s')
      write.table(df, file=file.path(getwd(),'Inputs','Inputs_Basins.txt'), sep='\t', col.names=TRUE, row.names=FALSE)
      message('Done!')
      toc()
}# End (not run)
