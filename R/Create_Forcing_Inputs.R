#' Extract and prepare model's inputs data in the airGR format (DatesR, P and E) from gridded P and E monthly data.
#' @param Subbasins Subbasins' shapefile. Must contain the following attributes: 'Area' (in km2), 'Region' (in letters), and 'COMID' (identifier number).
#' @param Precip  Raster brick of the precipitation data in [mm/month].
#' @param PotEvap Raster brick of the evapotranspiration data in [mm/month].
#' @param Qobs Observed streamflow data in [m3/s] at the basin outlet. Must have the dates of the output dataset.
#' @param DateIni Initial date of the output database in 'mm/yyyy' format.
#' @param DateEnd Ending date of the output database in 'mm/yyyy' format.
#' @param IniGrids  Initial date of the gridded data (P and PE) in 'mm/yyy' format.
#' @param Save   Boolean to save results as a text file in the 'Outputs' location. FALSE as default.
#' @param Update Boolean for the updating mode where only the last month's values will be returned. FALSE as default.
#' @param Members Número de miembros del conjunto modelo. Only for streamflow forecasting purposes. NULL por defecto.
#' @return Return a dataframe of model's inputs data in the airGR format (DatesR, P, E, Q).
#' @references Cesar Aybar, Carlos Fernández, Adrian Huerta, Waldo Lavado, Fiorella Vega & Oscar Felipe-Obando (2020) Construction of a high-resolution gridded rainfall dataset for Peru from 1981 to the present day, Hydrological Sciences Journal, 65:5, 770-785, DOI: 10.1080/02626667.2019.1649411
#' @references Llauca H, Lavado-Casimiro W, Montesinos C, Santini W, Rau P. PISCO_HyM_GR2M: A Model of Monthly Water Balance in Peru (1981–2020). Water. 2021; 13(8):1048. https://doi.org/10.3390/w13081048
#' @export
#' @examples
#' # Load data
#' library(GR2MSemiDistr)
#' data(pisco_pr)
#' data(pisco_pe)
#' data(qobs)
#' data(roi)
#'
#' # Create a database with model's inputs data
#' data <- Create_Forcing_Inputs(Subbasins=roi,
#'                               Precip=pisco_pr,
#'                               PotEvap=pisco_pe,
#'                               Qobs=qobs,
#'                               DateIni='01/1981',
#'                               DateEnd='12/2016',
#'                               IniGrids='01/1981')
#' View(data)
#' @import exactextractr
#' @import lubridate
#' @import tictoc
Create_Forcing_Inputs <- function(Subbasins,
                                  Precip,
                                  PotEvap,
                                  Qobs=NULL,
                                  DateIni,
                                  DateEnd,
                                  IniGrids='01/1981',
                                  Save=FALSE,
                                  Update=FALSE,
                                  Members=NULL){

  # Load required packages
  library(exactextractr)
  library(lubridate)
  library(tictoc)
  tic()

  # Load subbasin data
  roi   <- Subbasins
  comid <- as.vector(roi$COMID)
  nsub  <- length(comid)

  # Extract monthly precipitation data for each subbasin
  if(is.null(Precip)==FALSE){
    cat('\f')
    message('Calculating monthly mean-areal precipitation [mm/month]')
    message('Please wait...')
    pr <- Precip
    if(Update==TRUE){
      pr_mean <- round(t(exact_extract(pr, roi, 'mean', progress=FALSE)),1)[nlayers(pr),]
    }else{
      if(Members==1 | is.null(Members)==TRUE){
          dates <- seq(from=as.Date(paste0('01/',IniGrids), '%d/%m/%Y'),
                       to=as.Date(paste0('01/',IniGrids), '%d/%m/%Y') + months(nlayers(pr)-1),
                       by='month')
          Ini   <- as.Date(paste0('01/',DateIni),'%d/%m/%Y')
          End   <- as.Date(paste0('01/',DateEnd),'%d/%m/%Y')
          ind   <- which(dates==Ini):which(dates==End)
          dates <- dates[ind]
          pr    <- pr[[ind]]
        }
      pr_mean <- round(t(exact_extract(pr, roi, 'mean', progress=FALSE)),1)
    }
  }


  # Extract monthly mean-areal potential evapotranspiration
  if(is.null(PotEvap)==FALSE){
    cat('\f')
    message('Calcutaling monthly mean-areal evapotranspiration [mm/month]')
    message('Please wait...')
    pe <- PotEvap
    if(Update==TRUE){
      pe_mean <- round(t(exact_extract(pe, roi, 'mean', progress=FALSE)),1)[nlayers(pe),]
    }else{
      if(Members==1 | is.null(Members)==TRUE){
        dates <- seq(from=as.Date(paste0('01/',IniGrids), '%d/%m/%Y'),
                     to=as.Date(paste0('01/',IniGrids), '%d/%m/%Y') + months(nlayers(pe)-1),
                     by='month')
        Ini   <- as.Date(paste0('01/',DateIni),'%d/%m/%Y')
        End   <- as.Date(paste0('01/',DateEnd),'%d/%m/%Y')
        ind   <- which(dates==Ini):which(dates==End)
        dates <- dates[ind]
        pe    <- pe[[ind]]
      }
      pe_mean <- round(t(exact_extract(pe, roi, 'mean', progress=FALSE)),1)
    }
  }


  # Create a vector of dates
  if(Update==TRUE){
    DatesMonths <- as.Date(paste0('01/',DateEnd), "%d/%m/%Y")
  }else{
    Ini <- paste0('01/',DateIni)
    End <- paste0('01/',DateEnd)
    DatesMonths <- seq(as.Date(Ini, "%d/%m/%Y"),
                       as.Date(End, "%d/%m/%Y"),
                       by='month')
    if(is.null(Members)==FALSE){
      DatesMonths <- rep(DatesMonths, times=Members)
    }
  }


  # Prepare database in airGR format
    if(is.null(Qobs)==FALSE & is.null(Precip)==FALSE & is.null(PotEvap)==FALSE){
      qm_obs        <- Qobs[,2]
      Ans           <- data.frame(DatesMonths, round(pr_mean,1), round(pe_mean,1), round(qm_obs,3))
      colnames(Ans) <- c('DatesR', paste0('P_',comid), paste0('E_',comid), 'Q')
    }
    if(is.null(Qobs)==FALSE & is.null(Precip)==FALSE & is.null(PotEvap)==TRUE){
      qm_obs        <- Qobs[,2]
      Ans           <- data.frame(DatesMonths, round(pr_mean,1), round(qm_obs,3))
      colnames(Ans) <- c('DatesR', paste0('P_',comid), 'Q')
    }
    if(is.null(Qobs)==FALSE & is.null(Precip)==TRUE & is.null(PotEvap)==FALSE){
      qm_obs        <- Qobs[,2]
      Ans           <- data.frame(DatesMonths, round(pe_mean,1), round(qm_obs,3))
      colnames(Ans) <- c('DatesR', paste0('E_',comid), 'Q')
    }
    if(is.null(Qobs)==TRUE & is.null(Precip)==FALSE & is.null(PotEvap)==FALSE){
      Ans           <- data.frame(DatesMonths, round(pr_mean,1), round(pe_mean,1))
      colnames(Ans) <- c('DatesR', paste0('P_',comid), paste0('E_',comid))
    }
    if(is.null(Qobs)==TRUE & is.null(Precip)==FALSE & is.null(PotEvap)==TRUE){
      Ans           <- data.frame(DatesMonths, round(pr_mean,1))
      colnames(Ans) <- c('DatesR', paste0('P_',comid))
    }
    if(is.null(Qobs)==TRUE & is.null(Precip)==TRUE & is.null(PotEvap)==FALSE){
      Ans           <- data.frame(DatesMonths, round(pe_mean,1))
      colnames(Ans) <- c('DatesR', paste0('E_',comid))
    }
    if(Save==TRUE){
      dir.create('./Inputs')
      write.table(Ans, file='./Inputs/Inputs_model.txt', sep='\t', col.names=TRUE, row.names=FALSE)
    }

  message('Done!')
  toc()
  return(Ans)

}# End (not run)
