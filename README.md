'GR2MSemiDistr' package
========================
This is a preliminary R package to run a semidistributed version of the GR2M hydrological model applying a Weighted Flow Accumulation algorithm in order to obtain streamflows values at each subbasin outlets.

For any issue or suggestion please write to Harold LLauca (hllauca@gmail.com).

Enjoy it!


Install package
============
In order to use this package please take a look at the following instructions.

install.packages("devtools")

devtools::install_github("hllauca/GR2MSemiDistr")

library(GR2MSemiDistr)


Instructions
============
# Este script ejecuta el modelo GR2M de forma semidistribuida (subcuencas), utilizando los paquetes
# 'airGR' y 'rgrass7'. La topología de la cuenca se identifica automáticamente a partir del raster
# de Flow Direction, utilizando el algorithm de acumulación ponderada del flujo (Weighted Flow
# Accumulation). La calibración de los parámetros del modelo GR2M (X1 y X2) se realiza de forma automática
# mediante el algoritmo de optimización global Shuffle Complex Evolution (SCE-UA) usando el paquete
# 'rtop'. Se requiere instalar previamente GRASS v.7 (https://grass.osgeo.org/download/software/) y
# habilitar la extensión 'r.accumulate' con el comando 'g.extension'.

# ©Harold LLauca
# Email: hllauca@gmail.com

rm(list=ls())    # Remover variables anteriores
options(warn=-1) # Suprimir warnings
cat('\f')        # Limpiar consola


## Configuración del modelo GR2M SemiDistribuido
#===============================================

  ## DATOS DE ENTRADA AL MODELO
  Location     <- 'D:/GR2M_PERU/GR2M_SemiDistr_Test'
  File.Data    <- 'Inputs_Basins.csv'
  File.Shape   <- 'Los_Amigos_prueba.shp'
  File.Raster  <- 'FlowDirection.tif'

  
  ## PERIODO DE EJECUCIÓN DEL MODELO
  WarmUp.Ini   <- '09/2008'
  WarmUp.End   <- '08/2009'
  RunModel.Ini <- '09/2009'
  RunModel.End <- '12/2013'
  
  
  ## PARÁMETROS DEL MODELO
  Model.HRU    <- rep('I',7)         # Región de calibración de cada subcuenca
  Model.Param  <- c(200, 0.2, 0.8)   # Parámetros X1 y X2 para cada región
  No.OptimHRU  <- NULL               # HRUs que no se optimizarán (NULL de no existir)
  WFAC         <- FALSE              # Realizar la acumulación ponderada del flujo
  
  ## OPTIMIZACIÓN AUTOMÁTICA DEL MODELO
  Optim        <- TRUE               # Realizar optimización?
  Optim.Max    <- 1                  # Máx número de iteraciones
  Optim.Eval   <- 'NSE'              # Criterio de desempeño (NSE, lnNSE, R, RMSE, KGE)
  Optim.Basin  <- 7                  # Subcuenca pto. de control
  Optim.Remove <- FALSE              # Elimina Qsim en la subcuenca no deseada (FALSE por defecto)
  Model.ParMin <- c(1, 0.001, 0)     # Mínimos valores de X1 y X2
  Model.ParMax <- c(2000, 2, 1)      # Máximos valores de X1 y X2

  

###################################################################################################
########################################## NO MODIFICAR ###########################################
###################################################################################################
# Directorio de trabajo
  setwd(Location)

  
# Crear carpeta de resultados
  dir.create(file.path(Location, '5_OUTPUT'), recursive=T, showWarnings=F)

  
# Condicional para la optimización
if (Optim == TRUE){
  
# Optimizar parámetros X1 y X2 del modelo GR2M semidistribuido
#=============================================================
  
  # Cargar función
  source(file.path(Location,'1_FUNCIONES','Optim_GR2M_SemiDistr.R'))
  
  # Ejecutar optimización de parámetros del modelo GR2M semidistribuido
  x <- Optim_GR2M_SemiDistr(Parameters=Model.Param,
                            Parameters.Min=Model.ParMin,
                            Parameters.Max=Model.ParMax,
                            Max.Optimization=Optim.Max,
                            Optimization=Optim.Eval,
                            HRU=Model.HRU,
                            WorkDir=Location,
                            Shapefile=File.Shape,
                            Input=File.Data,
                            WarmIni=WarmUp.Ini,
                            WarmEnd=WarmUp.End,
							              RunIni=RunModel.Ini,
                            RunEnd=RunModel.End,
                            IdBasin=Optim.Basin,
                            Remove=Optim.Remove,
                            No.Optim=No.OptimHRU)

	# Mostrar resultados
	print(paste0('Parameters: ',x$par))
	print(paste0('F.O.: ', 1-x$value))
	
	# Guardar resultados
	save(x, file=file.path(Location,'5_OUTPUT','Parameters_Optim.Rda'))
	
	# Reemplazar parámetros por defecto
	Model.Param <- x$par
}

 
# Ejecutar modelo GR2M semidistribuido
#=====================================
  
  # Cargar función
  source(file.path(Location,'1_FUNCIONES','Run_GR2M_SemiDistr.R'))

  # Ejecutar modelo GR2M semidistribuido
  y  <- Run_GR2M_SemiDistr(Parameters=Model.Param,
                           HRU=Model.HRU,
                           WorkDir=Location,
                           Raster=File.Raster,
                           Shapefile=File.Shape,
                           Input=File.Data,
						               WarmIni=WarmUp.Ini,
                           WarmEnd=WarmUp.End,
                           RunIni=RunModel.Ini,
                           RunEnd=RunModel.End,
                           IdBasin=Optim.Basin,
                           Remove=Optim.Remove,
                           Plot=TRUE,
						               IniState=NULL,
						               wfac=WFAC)
  
  # Guardar resultados
  save(y, file=file.path(Location,'5_OUTPUT','Results_GR2M_Semidistr.Rda'))
  
  #Guardar rasters
  dir.create(file.path(Location,'5_OUTPUT','Qsim_rasters'))
  i=1
  rastername <- paste0('Qsim_', y$Dates[i],'.tif')
  writeRaster(y$Qras[[i]], file.path(Location,'5_OUTPUT','Qsim_rasters', rastername), overwrite=TRUE)
  
  
  # Guardar caudales generados en cada subcuenca (en formato .csv)
  Qout           <- data.frame(y$Dates, y$Qsim)
  if (WFAC==TRUE){
    colnames(Qout) <- c('Fecha', paste0('Qsim-', 1:length(Model.HRU)))
  } else{
    colnames(Qout) <- c('Fecha', 'Qoulet')
  }
  write.table(Qout, file=file.path(Location,'5_OUTPUT','Results_GR2M_Semidistr_Qsim.csv'), sep=',', row.names=F)
