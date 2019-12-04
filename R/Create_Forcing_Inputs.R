#' Create a text file with data inputs for the GR2M model
#'
#' @param Shapefile Subbasins shapefile.
#' @param Database Directory where precipitation and evapotranspiration data (as netCDF) are located.
#' @param Precip Precipitation filename.
#' @param PotEvap Evapotranspiration filename.
#' @param Qobs Observed streamflow filename (data in m3/s). NULL as default.
#' @param Resolution Raster resolution to resample forcing data and extract areal mean values. 0.01 as default.
#' @param factor Factor from 1 to 1.5 to crop data inputs for an area of interest. 1 as default
#' @param positions Position numbers to mask data by subbasins. NULL as default
#' @param DateIni Initial date 'mm/yyyy' for export data. '01/1981' as default
#' @param DateEnd Final date 'mm/yyyy' for export data. '12/2016' as default
#' @return Export a text file with forcing data inputs (Dates, Precip, Evap, Qobs).
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  tictoc
#' @import  ncdf4
#' @import  parallel
Create_Forcing_Inputs <- function(Shapefile, Database, Precip, PotEvap, Qobs=NULL,
                                  Resolution=0.01, factor=1, positions=NULL, DateIni='01/1981', DateEnd='12/2016'){

# Shapefile=File.Shape
# Database=Database
# Precip=File.Precip
# PotEvap=File.PotEvap
# Qobs=NULL
# Resolution=0.01
# factor=1
# positions=NULL
# DateIni="01/1979"
# DateEnd="11/2019"

    # Load packages
      require(rgdal)
      require(raster)
      require(rgeos)
      require(tictoc)
      require(ncdf4)
      require(parallel)
      tic()


    # Useful functions
    #=================
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


    # Auxiliary variables
    #====================
    # Create a vector of dates
      Ini <- paste0('01/',DateIni)
      End <- paste0('01/',DateEnd)
      DatesMonths <- seq(as.Date(Ini, "%d/%m/%Y"), as.Date(End, "%d/%m/%Y"), by='month')

    # Load subbasins shapefiles
      Basins  <- readOGR(file.path(getwd(),'Inputs', Shapefile))
      nBasins <- length(Basins)


    # Extract monthly precipitation data for each subbasin
    #=====================================================
    # Show message
      cat('\f')
      message('Extracting monthly precipitation [mm]')
      message('Please wait...')

    # Read precipitation data
      pp <- brick(file.path(Database, Precip))

    # Crop for basin domain
      pp.crop <- crop(pp, extent(Basins)*factor)

    # Mask for resampling using 'ngb' method
      pp.res      <- raster(extent(pp.crop[[1]]))
      crs(pp.res) <- crs(pp.crop)
      res(pp.res) <- Resolution

    # Cells within each subbasin
      if(is.null(positions)==TRUE){
          positionPP  <- generate_mask_geom(pp.res[[1]], Basins)
      } else{
          positionPP  <- positions$PP
      }

    # Mean areal rainfall for each subbasin
      cl=makeCluster(detectCores()-1)
      clusterEvalQ(cl,c(library(raster)))
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
    #===============================================================
    # Show message
      cat('\f')
      message('Extracting monthly evapotranspiration [mm]')
      message('Please wait...')

    # Read potential evapotranspiration data
      pet      <- brick(file.path(Database, PotEvap))

    # Crop for basin domain
      pet.crop <- crop(pet, extent(Basins)*factor)

    # Mask for resampling using 'ngb' method
      pet.res      <- raster(extent(pet.crop[[1]]))
      crs(pet.res) <- crs(pet.crop)
      res(pet.res) <- Resolution

    # Cells within each subbasin
      # Cells within each subbasin
      if(is.null(positions)==TRUE){
        positionPET  <- generate_mask_geom(pet.res[[1]], Basins)
      } else{
        positionPET  <- positions$PET
      }

      positionPET  <- generate_mask_geom(pet.res[[1]], Basins)

    # Mean areal evapotranspiration for each subbasin
      clusterEvalQ(cl,c(library(raster)))
      clusterExport(cl,varlist=c("pet.res","pet.crop","positionPET","mask_Fast_Extract"),envir=environment())

      mean.pet <- parLapply(cl, 1:nlayers(pet.crop), function(z) {
        res <- raster::resample(pet.crop[[z]], pet.res, method='ngb')
        ans <- mask_Fast_Extract(cov=res, positionPET, fun=mean, na.rm=T)
        return(ans)
      })
      stopCluster(cl) #Close the cluster
      mean.pet <- do.call(rbind, mean.pet)
      mean.pet <- round(mean.pet[1:length(DatesMonths),],1)


    # Export results in airGR format
    #===============================
    # Data to export
      if(is.null(Qobs)==TRUE){
        df      <- data.frame(DatesMonths, mean.pp, mean.pet)
        colnames(df) <- c('DatesR', paste0('P',1:nBasins), paste0('E',1:nBasins))
      } else{
        # Subsetting streamflow data
        Obs     <- read.table(file.path("Inputs",Qobs), sep='\t', header=F)
        Ind_Obs <- seq(which(format(as.Date(Obs[,1]), "%d/%m/%Y") == Ini),
                       which(format(as.Date(Obs[,1]), "%d/%m/%Y") == End))
        qobs    <- Obs[Ind_Obs,2]
        df      <- data.frame(DatesMonths, mean.pp, mean.pet, qobs)
        colnames(df) <- c('DatesR', paste0('P',1:nBasins), paste0('E',1:nBasins), 'Qm3s')
      }

    # Saving data as text file
      write.table(df, file=file.path(getwd(),'Inputs','Inputs_Basins.txt'), sep='\t', col.names=TRUE, row.names=FALSE)


    # Saveing positions
    if(is.null(positions)==TRUE){
      positions <- list(PP=positionPP, PET=positionPET)
      save(positions, file=file.path(getwd(),'Mask4extract.Rda'))
    }

    # End
    message('Done!')
    toc()

}# End (not run)
