#' Prepare model inputs at subbasin scales in airGR format from gridded datasets.
#'
#' @param Subbasins Subbasins shapefile.
#' @param Precip  Netcdf file for precipitation (in mm/month).
#' @param PotEvap Netcdf file for potential evapotranspiration (in mm/month).
#' @param Qobs Observed streamflow (in m3/s). NULL as default.
#' @param DateIni Initial date of the data (in mm/yyyy format).
#' @param DateEnd Final date of the data (in mm/yyyy format).
#' @param Save   Boolean to save database as textfile. FALSE as default.
#' @param Update Boolean to extract the last values for model updating. FALSE as default.
#' @param Resolution Resolution to resample gridded data. 0.01 as default.
#' @param Buffer Multiplicative factor to buffer subbasins extents. 1.1 as default.
#' @param Members Number of ensemble members for model forcasting. NULL as default.
#' @param Horiz Number of months for model forcasting. NULL as default.

#' @return Return a dataframe with datavase in airGR format (DatesR,P,E,Q).
#' @export
#' @import rgdal
#' @import raster
#' @import rgeos
#' @import tictoc
#' @import parallel
#' @import exactextractr
#' @import sf
Create_Forcing_Inputs <- function(Subbasins,
                                  Precip,
                                  PotEvap,
                                  Qobs=NULL,
                                  DateIni,
                                  DateEnd,
                                  Save=FALSE,
                                  Update=FALSE,
                                  Resolution=0.01,
                                  Buffer=1.1,
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
  # Resolution=0.01
  # Buffer=1.1
  # Members=NULL
  # Horiz=NULL

  # Load required packages
  require(rgdal)
  require(raster)
  require(rgeos)
  require(parallel)
  require(exactextractr)
  require(sf)
  require(tictoc)
  tic()

  # Load subbasin data
  roi   <- st_as_sf(Subbasins)
  comid <- as.vector(roi$COMID)
  nsub  <- length(comid)

  # Extract monthly precipitation data for each subbasin
  # Show message
  cat('\f')
  message('Calculating monthly mean-areal precipitation [mm]')
  message('Please wait...')

  # Read precipitation data
  pr <- brick(Precip)
  if(Update==TRUE){
    pr <- pr[[nlayers(pr)]]
  }

  # Buffer subbasins and resample raster
  pr_buf      <- crop(pr, extent(roi)*Buffer)
  pr_res      <- raster(extent(pr_buf[[1]]))
  crs(pr_res) <- crs(pr_buf)
  res(pr_res) <- Resolution

  # Mean-areal precipitation for each subbasin
  cl=makeCluster(detectCores()-1)
  clusterEvalQ(cl,c(library(exactextractr), library(raster)))
  clusterExport(cl, varlist=c("pr_buf","pr_res","roi"), envir=environment())
  pr_mean <- parLapply(cl, 1:nlayers(pr_buf), function(z) {
    res <- resample(pr_buf[[z]], pr_res, method='ngb')
    ans <- as.numeric(exact_extract(res, roi, fun='mean'))
    return(ans)
  })
  if(Update==TRUE){
    pr_mean <- round(unlist(pr_mean),1)
    pr_mean <- matrix(pr_mean, ncol=length(pr_mean))
  }else{
    pr_mean <- do.call(rbind, pr_mean)
    pr_mean <- round(pr_mean,1)
  }


  # Extract monthly mean-areal potential evapotranspiration
  # Show message
  cat('\f')
  message('Calcutaling monthly mean-areal pot. evapotranspiration [mm]')
  message('Please wait...')

  # Read potential evapotranspiration data
  pe <- brick(PotEvap)
  if(Update==TRUE){
    pe <- pe[[nlayers(pe)]]
  }

  # Buffer subbasins and resample raster
  pe_buf      <- crop(pe, extent(roi)*Buffer)
  pe_res      <- raster(extent(pe_buf[[1]]))
  crs(pe_res) <- crs(pe_buf)
  res(pe_res) <- Resolution

  # Mean-areal evapotranspiration for each subbasin
  cl=makeCluster(detectCores()-1)
  clusterEvalQ(cl,c(library(exactextractr), library(raster)))
  clusterExport(cl, varlist=c("pe_buf","pe_res","roi"), envir=environment())
  pe_mean <- parLapply(cl, 1:nlayers(pe_buf), function(z) {
    res <- resample(pe_buf[[z]], pe_res, method='ngb')
    ans <- as.numeric(exact_extract(res, roi, fun='mean'))
    return(ans)
  })
  stopCluster(cl) #Close the cluster
  if(Update==TRUE){
    pe_mean <- round(unlist(pe_mean),1)
    pe_mean <- matrix(pe_mean,ncol=length(pe_mean))
  }else{
    pe_mean <- do.call(rbind, pe_mean)
    pe_mean <- round(pe_mean,1)
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

  # Prepare database in airGR format
  if(is.null(Qobs)==TRUE){
    Ans           <- data.frame(DatesMonths, round(pr_mean,1), round(pe_mean,1))
    colnames(Ans) <- c('DatesR', paste0('P_',comid), paste0('E_',comid))
  }else{
    qm_obs        <- Qobs[,2]
    Ans           <- data.frame(DatesMonths, round(pr_mean,1), round(pe_mean,1), round(qm_obs,3))
    colnames(Ans) <- c('DatesR', paste0('P_',comid), paste0('E_',comid), 'Q')
  }

  # Save database as textfile
  if(Save==TRUE){
    dir.create('./Inputs')
    write.table(Ans, file='./Inputs/Inputs_model.txt', sep='\t', col.names=TRUE, row.names=FALSE)
  }

  # Final message
  message('Done!')
  toc()
  return(Ans)

}# End (not run)
