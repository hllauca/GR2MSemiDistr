#' Prepare model data inputs in airGR format.
#'
#' @param Subbasins Subbasins shapefile.
#' @param Precip  Netcdf filename for precipitation dataset.
#' @param PotEvap Netcdf filename for potential evapotranspiration dataset.
#' @param Qobs Observed streamflow (in m3/s). NULL as default.
#' @param DateIni Initial date for subsetting data (in mm/yyyy format).
#' @param DateEnd Final date for subsetting data (in mm/yyyy format).
#' @param Save   Boolean to save database as textfile. FALSE as default.
#' @param Update Boolean to extract the last value for updating model. FALSE as default.
#' @param Positions Cell numbers of subbasins to extract data faster. NULL as default
#' @param Resolution Resolution to resample gridded-datasets. 0.01 as default.
#' @param Buffer Factor to create a buffer of subbasins. 1 as default.
#' @param Members Number of ensemble members for forcasting. NULL as default.
#' @param Horiz Number of months for forcasting. NULL as default.

#' @return Return a dataframe with data inputs in airGR format (DatesR,P,Ep,Q).
#' @export
#' @import  rgdal
#' @import  raster
#' @import  rgeos
#' @import  tictoc
#' @import  parallel
Create_Forcing_Inputs <- function(Subbasins,
                                  Precip,
                                  PotEvap,
                                  Qobs=NULL,
                                  DateIni,
                                  DateEnd,
                                  Save=FALSE,
                                  Update=FALSE,
                                  Positions=NULL,
                                  Resolution=0.01,
                                  Buffer=1,
                                  Members=NULL,
                                  Horiz=NULL){

  # Subbasins=roi
  # Precip=file.path(Inputs,File.Precip)
  # PotEvap=file.path(Inputs,File.PotEvap)
  # Qobs=qdf
  # DateIni=WarmUp.Ini
  # DateEnd=RunModel.End
  # Save=FALSE
  # Update=FALSE
  # Positions=NULL
  # Resolution=0.01
  # Buffer=1
  # Members=NULL
  # Horiz=NULL

  # Load packages
  require(rgdal)
  require(raster)
  require(rgeos)
  require(tictoc)
  require(parallel)
  tic()

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
  if(Update==TRUE){
    DatesMonths <- as.Date(paste0('01/',DateEnd), "%d/%m/%Y")
  }
  if(Update==FALSE){
    Ini <- paste0('01/',DateIni)
    End <- paste0('01/',DateEnd)
    DatesMonths <- seq(as.Date(Ini, "%d/%m/%Y"),
                       as.Date(End, "%d/%m/%Y"),
                       by='month')
    if(is.null(Members)==FALSE){
      DatesMonths <- rep(DatesMonths, length=Horiz*Members)
    }
  }

  # Load subbasin data
  sub.id <- paste0('GR2M_ID_',Subbasin@data@GR2M_ID)
  nsub   <- nrow(Subbasin@data)


  # Extract monthly precipitation data for each subbasin
  # Show message
  cat('\f')
  message('Calculating areal monthly precipitation [mm]')
  message('Please wait...')

  # Read precipitation data
  pp <- brick(Precip)
  if(Update==TRUE){
    pp <- pp[[nlayers(pp)]]
  }

  # Crop for basin domain
  pp.crop <- crop(pp, extent(Subbasins)*Buffer)

  # Mask for resampling using 'ngb' method
  pp.res      <- raster(extent(pp.crop[[1]]))
  crs(pp.res) <- crs(pp.crop)
  res(pp.res) <- Resolution

  # Cells within each subbasin
  if(is.null(Positions)==TRUE){
    positionPP  <- generate_mask_geom(pp.res[[1]], Subbasins)
  } else{
    positionPP  <- Positions$PP
  }

  # Mean areal rainfall for each subbasin
  cl=makeCluster(detectCores()-1)
  clusterEvalQ(cl,c(library(raster)))
  clusterExport(cl,varlist=c("pp.res","pp.crop","positionPP",
                             "mask_Fast_Extract"),envir=environment())
  mean.pp <- parLapply(cl, 1:nlayers(pp.crop), function(z) {
    res <- resample(pp.crop[[z]], pp.res, method='ngb')
    ans <- mask_Fast_Extract(cov=res,
                             positionPP,
                             fun=mean,
                             na.rm=T)
    return(ans)
  })
  # stopCluster(cl) #Close the cluster
  if(Update==TRUE){
    mean.pp <- round(unlist(mean.pp),1)
    mean.pp <- matrix(mean.pp,ncol=length(mean.pp))
  }else{
    mean.pp <- do.call(rbind, mean.pp)
    mean.pp <- round(mean.pp[1:length(DatesMonths),],1)
  }


  # Extract monthly potential evapotranspiration for each subbasin
  # Show message
  cat('\f')
  message('Calcutaling areal monthly pot. evapotranspiration [mm]')
  message('Please wait...')

  # Read potential evapotranspiration data
  pet      <- brick(PotEvap)
  if(Update==TRUE){
    pet <- pet[[nlayers(pet)]]
  }

  # Crop for basin domain
  pet.crop <- crop(pet, extent(Subbasins)*Buffer)

  # Mask for resampling using 'ngb' method
  pet.res      <- raster(extent(pet.crop[[1]]))
  crs(pet.res) <- crs(pet.crop)
  res(pet.res) <- Resolution

  # Cells within each subbasin
  # Cells within each subbasin
  if(is.null(Positions)==TRUE){
    positionPET  <- generate_mask_geom(pet.res[[1]], Subbasins)
  } else{
    positionPET  <- Positions$PET
  }

  # Mean areal evapotranspiration for each subbasin
  clusterEvalQ(cl,c(library(raster)))
  clusterExport(cl,varlist=c("pet.res","pet.crop","positionPET",
                             "mask_Fast_Extract"),envir=environment())

  mean.pet <- parLapply(cl, 1:nlayers(pet.crop), function(z) {
    res <- raster::resample(pet.crop[[z]], pet.res, method='ngb')
    ans <- mask_Fast_Extract(cov=res,
                             positionPET,
                             fun=mean,
                             na.rm=T)
    return(ans)
  })
  stopCluster(cl) #Close the cluster
  if(Update==TRUE){
    mean.pet <- round(unlist(mean.pet),1)
    mean.pet <- matrix(mean.pet,ncol=length(mean.pet))
  }else{
    mean.pet <- do.call(rbind, mean.pet)
    mean.pet <- round(mean.pet[1:length(DatesMonths),],1)
  }

  # Export results in airGR format
  if(is.null(Qobs)==TRUE){
    Ans <- data.frame(DatesMonths, round(mean.pp,1), round(mean.pet,1))
    colnames(Ans) <- c('DatesR', paste0('P_',sub.id), paste0('E_',sub.id))
  } else{
    # Subsetting streamflow data
    Ind_Obs <- seq(which(format(as.Date(Qobs[,1], tryFormats=c('%d/%m/%Y','%Y/%m/%d','%d-%m-%Y','%Y-%m-%d')), "%d/%m/%Y") == Ini),
                   which(format(as.Date(Qobs[,1], tryFormats=c('%d/%m/%Y','%Y/%m/%d','%d-%m-%Y','%Y-%m-%d')), "%d/%m/%Y") == End))
    qobs <- Qobs[Ind_Obs,2]
    Ans  <- data.frame(DatesMonths, round(mean.pp,1), round(mean.pet,1), round(qobs,3))
    colnames(Ans) <- c('DatesR', paste0('P_',sub.id), paste0('E_',sub.id), 'Q')
  }

  # Saving data as text file
  if(Save==TRUE){
    dir.create(file.path(getwd(),'Inputs'))
    write.table(Ans, file=file.path(getwd(),'Inputs','Inputs_model.txt'),
                sep='\t', col.names=TRUE, row.names=FALSE)
  }

  # Saving positions
  if(is.null(Positions)==TRUE){
    Positions_Dat <- list(PP=positionPP, PET=positionPET)
    save(Positions_Dat, file=file.path(getwd(),'Inputs','Positions_Dat.Rda'))
  }

  # End
  message('Done!')
  toc()
  return(Ans)

}# End (not run)
