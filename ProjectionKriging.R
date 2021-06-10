####--------------- [0] PREAMBLE ----------------------------------------------------
rm(list = ls()) # clearing environment
####--------------- .   PACKAGES ----------------------------------------------------
if("KrigR" %in% rownames(installed.packages()) == FALSE){ # KrigR check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("https://github.com/ErikKusch/KrigR")
}
library(KrigR) 

try(source("PersonalSettings.R")) # I do this here to specify number of cores and API credentials and am thus not sharing this file
source("KrigSleep.R") # version of Kriging that enables sys-sleep at the end of each computation to circumvent overheating of RAM modules

#### CDS API (needed for ERA5-Land downloads)
if(!exists("API_Key") | !exists("API_User")){ # CS API check: if CDS API credentials have not been specified elsewhere
  API_User <- readline(prompt = "Please enter your Climate Data Store API user number and hit ENTER.")
  API_Key <- readline(prompt = "Please enter your Climate Data Store API key number and hit ENTER.")
} # end of CDS API check

#### NUMBER OF CORES
if(!exists("numberOfCores")){ # Core check: if number of cores for parallel processing has not been set yet
  numberOfCores <- as.numeric(readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores())))
} # end of Core check

####--------------- .   DIRECTORIES -------------------------------------------------
mainDir <- getwd() # extract the project folder location
Dir.Projections <- file.path(mainDir, "Projections")
# WORKING DIRECTORY FOR RAW DATA
Dir.Data <- paste(Dir.Projections, "/X - Raw Data", sep="")
if(!dir.exists(Dir.Data)){stop("The raw projection data cannot be found!")}
# WORKING DIRECTORY FOR KRIGING COVARIATES
Dir.COV <- paste(Dir.Projections, "/Z - Covariates", sep="")
if(!dir.exists(Dir.COV)){dir.create(Dir.COV)}
Dir.Mask <- paste(Dir.Projections, "/Y - ShapeFiles", sep="")
if(!dir.exists(Dir.Mask)){dir.create(Dir.Mask)}
Dir.Historical <- file.path(Dir.Projections, "HISTORICAL Aggregate")
if(!dir.exists(Dir.Historical)){dir.create(Dir.Historical)}

####--------------- .   CHECKS & PREPARATIONS ------------------------------------------
#### LAND MASK (for masking species in the sea which are terrestrial and marine)
if(!file.exists(file.path(Dir.Mask, "LandMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_land.zip", destfile = paste(Dir.Mask, "LandMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Mask, "LandMask.zip", sep="/"), exdir = Dir.Mask) # unzip the data
}
Land_shp <- raster::shapefile(file.path(Dir.Mask , "ne_10m_land.shp"))

####--------------- .   CLIMATOLOGY IDENTIFICATION -------------------------------------
Climatologies_fs <- list.files(Dir.Data) # identify file names of climatologies which need downscaling

####--------------- .   RUN SETTINGS ---------------------------------------------------
NMax = 80

####--------------- [1] DOWNSCALING OF INDEPENDENT PROJECTIONS ---------------------------
for(Climatology_Iter in 1:length(Climatologies_fs)){
  ####--------------- .   CLIMATE DATA LOADING -------------------------------------------
  Dir.Clima <- file.path(Dir.Projections, Climatologies_fs[Climatology_Iter])
  dir.create(Dir.Clima)
  FileName <- paste("K", NMax, Climatologies_fs[Climatology_Iter], sep ="_")
  print(paste("########## HANDLING", Climatologies_fs[Climatology_Iter], "with nmax set to", NMax,  "##########################################"))
  if(file.exists(file.path(Dir.Clima, FileName))){
    print(paste(FileName, "already kriged"))
    next()
  }
  Data <- stack(file.path(Dir.Data, Climatologies_fs[Climatology_Iter]))
  extent(Data) <- raster::extent(-180,180,-90,90)
  
  ####--------------- .   COVARIATE DATA DOWNLOAD ----------------------------------------
  print("Loading covariate data. #####################################################")
  if(file.exists(file.path(Dir.COV, "GMTED2010_Target.nc"))){
    print("Covariates already downloaded and processed.")
    Covs_ls <- list(raster::raster(file.path(Dir.COV, "GMTED2010_Train.nc")),
                    raster::raster(file.path(Dir.COV, "GMTED2010_Target.nc")))
  }else{
    print("This process takes about 25min (on the AU server, at least).")
    Covs_ls <- KrigR::download_DEM(Train_ras = raster::raster(file.path(Dir.Data, Climatologies_fs[Climatology_Iter])), # or this want to run off the netcdf intermediate?
                                   Target_res = .01,
                                   Dir = Dir.COV,
                                   Keep_Temporary = TRUE
    )
    Covs_ls[[2]] <- raster::raster(file.path(Dir.COV, "GMTED2010_Target.nc")) # to get around attribute issues
  }
  raster::extent(Covs_ls[[1]]) <- raster::extent(-180,180,-90,90)
  raster::extent(Covs_ls[[2]]) <- raster::extent(-180,180,-90,90)
  
  ####--------------- .   COVARIATE & DATA MASKING ---------------------------------------
  print("Masking data for landmasses. ################################################")
  if(file.exists(file.path(Dir.COV, "GMTED2010_Target_Landmasked.nc"))){
    print("Data already masked.")
    Covs_ls <- list(raster::raster(file.path(Dir.COV, "GMTED2010_Train_Landmasked.nc")),
                    raster::raster(file.path(Dir.COV, "GMTED2010_Target_Landmasked.nc")))
  }else{
    print("This process takes about 20min (on the AU server, at least).")
    Mask_Coarse <- KrigR:::mask_Shape(base.map = Covs_ls[[1]], Shape = Land_shp)
    Covs_ls[[1]] <- mask(Covs_ls[[1]], Mask_Coarse)
    writeRaster(Covs_ls[[1]], file.path(Dir.COV, "GMTED2010_Train_Landmasked.nc"))
    Mask_Fine <- KrigR:::mask_Shape(base.map = Covs_ls[[2]], Shape = Land_shp)
    Covs_ls[[2]] <- mask(Covs_ls[[2]], Mask_Fine)
    writeRaster(Covs_ls[[2]], file.path(Dir.COV, "GMTED2010_Target_Landmasked.nc"))
    rm(list = c("Mask_Fine", "Mask_Coarse"))
  }
  Mask_Data <- KrigR:::mask_Shape(Data[[1]], Land_shp)
  Data <- mask(Data, Mask_Data)
  rm("Mask_Data")
  # raster::extent(Covs_ls[[1]]) <- raster::extent(-180,180,-90,90)
  # raster::extent(Covs_ls[[2]]) <- raster::extent(-180,180,-90,90)
  
  ####--------------- .   KRIGING --------------------------------------------------------
  print("Kriging. ####################################################################")
  Begin <- Sys.time()
  Krigs_ls <- krigRLocalSleep(
    Data = stack(Data),
    Covariates_coarse = Covs_ls[[1]],
    Covariates_fine = Covs_ls[[2]],
    Cores = numberOfCores,
    Dir = Dir.Clima,
    FileName = FileName,
    Keep_Temporary = TRUE,
    nmax = NMax
  )
  End <- Sys.time()
  print(End-Begin)
  gc()
}

####--------------- [2] HISTORICAL ERA5-LAND KRIGING -----------------------------------
####--------------- .   DATA DOWNLOAD AND AGGREGATION ----------------------------------
if(!file.exists(file.path(Dir.Historical, "era5land_tas_1981-2000.nc"))){
  download_ERA(Variable = "2m_temperature",
               DateStart = "1981-01-01",
               DateStop = "1999-12-31",
               TResolution = "month",
               TStep = 1,
               Extent = Land_shp,
               Dir = Dir.Historical,
               API_Key = API_Key,
               API_User = API_User)
  AT_ras <- stack(file.path(Dir.Historical, "2m_temperature_1981-01-01_1999-12-31_month.nc"))
  Index <- rep(1:12, length = nlayers(AT_ras))
  ATClim_ras <- stackApply(AT_ras, indices = Index, fun = mean)
  writeRaster(ATClim_ras, filename = file.path(Dir.Historical, "era5land_tas_1981-2000"), format = "CDF")
}else{
  ATClim_ras <- stack(file.path(Dir.Historical, "era5land_tas_1981-2000.nc"))
}

####--------------- .   COVARIATE DATA DOWNLOAD ----------------------------------------
print("Loading covariate data. #####################################################")
if(file.exists(file.path(Dir.Historical, "GMTED2010_TargetERA.nc"))){
  print("Covariates already downloaded and processed.")
  Covs_ls <- list(raster::raster(file.path(Dir.Historical, "GMTED2010_Train.nc")),
                  raster::raster(file.path(Dir.Historical, "GMTED2010_Target.nc")))
}else{
  print("This process takes about 25min (on the AU server, at least).")
  Covs_ls <- KrigR::download_DEM(Train_ras = ATClim_ras,
                                 Target_res = .01,
                                 Dir = Dir.Historical,
                                 Keep_Temporary = TRUE
  )
  Covs_ls[[2]] <- raster::raster(file.path(Dir.Historical, "GMTED2010_Target.nc")) # to get around attribute issues
}
raster::extent(Covs_ls[[1]]) <- raster::extent(-180,180,-90,90)
raster::extent(Covs_ls[[2]]) <- raster::extent(-180,180,-90,90)

####--------------- .   COVARIATE & DATA MASKING ---------------------------------------
print("Masking data for landmasses. ################################################")
if(file.exists(file.path(Dir.Historical, "GMTED2010_Target_Landmasked.nc"))){
  print("Data already masked.")
  Covs_ls <- list(raster::raster(file.path(Dir.Historical, "GMTED2010_Train_Landmasked.nc")),
                  raster::raster(file.path(Dir.Historical, "GMTED2010_Target_Landmasked.nc")))
}else{
  print("This process takes about 20min (on the AU server, at least).")
  Mask_Coarse <- KrigR:::mask_Shape(base.map = Covs_ls[[1]], Shape = Land_shp)
  Covs_ls[[1]] <- mask(Covs_ls[[1]], Mask_Coarse)
  writeRaster(Covs_ls[[1]], file.path(Dir.Historical, "GMTED2010_Train_Landmasked.nc"))
  Mask_Fine <- KrigR:::mask_Shape(base.map = Covs_ls[[2]], Shape = Land_shp)
  Covs_ls[[2]] <- mask(Covs_ls[[2]], Mask_Fine)
  writeRaster(Covs_ls[[2]], file.path(Dir.Historical, "GMTED2010_Target_Landmasked.nc"))
  rm(list = c("Mask_Fine", "Mask_Coarse"))
}
Mask_Data <- KrigR:::mask_Shape(ATClim_ras[[1]], Land_shp)
Data <- mask(ATClim_ras, Mask_Data)
rm("Mask_Data")
# raster::extent(Covs_ls[[1]]) <- raster::extent(-180,180,-90,90)
# raster::extent(Covs_ls[[2]]) <- raster::extent(-180,180,-90,90)

####--------------- .   KRIGING --------------------------------------------------------
Dir.Clima <- file.path(mainDir, "era5land_tas_1981-2000.nc")
dir.create(Dir.Clima)
FileName <- paste("K", NMax, "era5land_tas_1981-2000.nc", sep ="_")
print("Kriging. ####################################################################")
Begin <- Sys.time()
Krigs_ls <- KrigR::krigR(
  Data = ATClim_ras,
  Covariates_coarse = Covs_ls[[1]],
  Covariates_fine = Covs_ls[[2]],
  Cores = numberOfCores,
  Dir = Dir.Clima,
  FileName = FileName,
  Keep_Temporary = TRUE,
  nmax = NMax
)
End <- Sys.time()
print(End-Begin)