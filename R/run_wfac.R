#' Run a Weighted Flow Accumulation (WFAC) algorithm for a subbasin-streamflow raster
#'
#' @param FlowDir Flow direction raster in GRASS format (directions from 1 to 8).
#' @param Qmodel Subbasin-streamflow raster. Q values only at sub-basin centroids cells.
#' @param Shapefile Subbasins shapefile.
#' @param nSub Numbers of subbasin to be processed.
#' @param t Timestep (for an iteration mode).
#' @return [Qsub] Streamflow values at the oulet of each subbasin. [Qacum] Flow accumulation raster.
#' @export0
run_wfac <- function(FlowDir, Qmodel, Shapefile, nSub, t){

    # Load packages
      require(rgrass7)

      if (t==1){
    # Import raster to GRASS
      writeRAST(as(FlowDir, 'SpatialGridDataFrame'), "fdr", overwrite=T)

    # Set an study extention
      execGRASS("g.region", Sys_show.output.on.console=F, parameters=list(raster="fdr"))
      }

    # Importar Qsim raster
      writeRAST(as(Qmodel, 'SpatialGridDataFrame'), "qweight", overwrite=T)

    # Weighted flow accumulation
      execGRASS("r.accumulate", flags=c("overwrite"),  Sys_show.output.on.console=F,
                parameters=list(direction="fdr",
                                weight='qweight',
                                accumulation="qacum"))

    # Load GRASS raster into R
      Qacum      <- raster(readRAST('qacum'))
      cat('\f')

    # Extract Qsim in each subbasin
      Qsub <- c()
      for(w in 1:nSub){
         Qsub[w] <- maxValue(setMinMax(mask(Qacum, Shapefile[w,])))
      }

    # Results
      ans <- list(Qsub=Qsub, Qacum=Qacum)
      return(ans)
}
