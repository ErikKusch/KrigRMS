####--------------- [0] PREAMBLE ----------------------------------------------------
rm(list=ls())
#### .   PACKAGES ---------------------------------------------------------
print("#### Loading KrigR. ####")
## KrigR Package and Credentials for CDS API
if("KrigR" %in% rownames(installed.packages()) == FALSE){ # KrigR check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("https://github.com/ErikKusch/KrigR")
}
library(KrigR)
try(source("PersonalSettings.R")) # I do this here to specify number of cores and API credentials and am thus not sharing this file
#### CDS API (needed for ERA5-Land downloads)
if(!exists("API_Key") | !exists("API_User")){ # CS API check: if CDS API credentials have not been specified elsewhere
  API_User <- readline(prompt = "Please enter your Climate Data Store API user number and hit ENTER.")
  API_Key <- readline(prompt = "Please enter your Climate Data Store API key number and hit ENTER.")
} # end of CDS API check
#### NUMBER OF CORES
if(!exists("numberOfCores")){ # Core check: if number of cores for parallel processing has not been set yet
  numberOfCores <- as.numeric(readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores())))
} # end of Core check

## Misc Packages
print("#### Loading PACKAGES. ####")
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "httr" # for donwloading zip files of HWSD data base
)
sapply(package_vec, install.load.package)

#### .   DIRECTORIES ---------------------------------------------------------
Dir.Base <- getwd() # read out the project directory
Dir.Figures <- file.path(Dir.Base, "Figures")
if(!dir.exists(Dir.Figures)){dir.create(Dir.Figures)}
Dir.Shapes <- file.path(Dir.Figures, "ShapeFiles")
if(!dir.exists(Dir.Shapes)){dir.create(Dir.Shapes)}
Dir.COV <- file.path(Dir.Figures, "Covariates")
if(!dir.exists(Dir.COV)){dir.create(Dir.COV)}

#### .   FUNCTIONALITY ---------------------------------------------------------
`%nin%` <- Negate(`%in%`)

#### .   COUNTRY MASK (for producing maps with national borders) -----------------------------------------------------------
print("#### Loading COUNTRY MASK. ####")
if(!file.exists(file.path(Dir.Shapes, "CountryMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip", destfile = paste(Dir.Shapes, "CountryMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "CountryMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
CountryMask <- readOGR(Dir.Shapes, "ne_10m_admin_0_countries", verbose = FALSE) # read country mask in
UK_shp <- CountryMask[CountryMask$NAME == "United Kingdom", ] # extracting UK shape

#### .   STATE MASK (for producing maps with state borders) -----------------------------------------------------------
print("#### Loading STATE MASK. ####")
if(!file.exists(file.path(Dir.Shapes, "StateMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip", destfile = paste(Dir.Shapes, "StateMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "StateMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
StateMask <- readOGR(Dir.Shapes, "ne_10m_admin_1_states_provinces", verbose = FALSE) # read state mask in
AK_shp <- StateMask[which(StateMask$name_en == "Alaska"), ] # extracting Alaska shape
AK_ext <- extent(AK_shp) # obtaining extent of Alaska shape
AK_ext[2] <-  -128.95 # crop out some of the westerly Aleutian islands
AK_shp <- crop(AK_shp,  AK_ext) # crop out some Aleutian islands

####--------------- [1] COVARIATES ----------------------------------------------------
#### .   SOIL COVARIATE RETRIEVAL (for kriging with other than DEM) -----------------------------------------------------------
SoilCovs_vec <- c("tkdry", "tksat", "csol", "k_s", "lambda", "psi", "theta_s") # need these names for addressing soil covariates
if(sum(file.exists(c(file.path(Dir.COV, "Covariates_UK.rds"), file.path(Dir.COV, "Covariates_AK.rds")))) != 2){
  print("#### Loading SOIL PROPERTY covariate data. ####") # documentation of these can be found here http://globalchange.bnu.edu.cn/research/soil4.jsp
  # create lists to combine soil data into one
  SoilCovs_ls <- as.list(rep(NA, length(SoilCovs_vec)))
  names(SoilCovs_ls) <- c(SoilCovs_vec)
  ## Downloading, unpacking, and stacking
  for(Soil_Iter in SoilCovs_vec){
    if(!file.exists(file.path(Dir.COV, paste0(Soil_Iter, ".nc")))) { # if not downloaded and processed yet
      print(paste("Handling", Soil_Iter, "data."))
      Dir.Soil <- file.path(Dir.COV, Soil_Iter)
      dir.create(Dir.Soil)
      download.file(paste0("http://globalchange.bnu.edu.cn/download/data/worldptf/", Soil_Iter,".zip"),
                    destfile = file.path(Dir.Soil, paste0(Soil_Iter, ".zip"))
      ) # download data
      unzip(file.path(Dir.Soil, paste0(Soil_Iter, ".zip")), exdir = Dir.Soil) # unzip data
      File <- list.files(Dir.Soil, pattern = ".nc")[1] # only keep first soil layer
      Soil_ras <- raster(file.path(Dir.Soil, File)) # load data
      SoilCovs_ls[[which(names(SoilCovs_ls) == Soil_Iter)]] <- Soil_ras # save to list
      writeRaster(x = Soil_ras, filename = file.path(Dir.COV, Soil_Iter), format = "CDF")
      plot(Soil_ras, main = Soil_Iter, colNA = "black")
      unlink(Dir.Soil, recursive = TRUE)
    }else{
      print(paste(Soil_Iter, "already downloaded and processed."))
      SoilCovs_ls[[which(names(SoilCovs_ls) == Soil_Iter)]] <- raster(file.path(Dir.COV, paste0(Soil_Iter, ".nc")))
    }
  }
  
  if(!file.exists(file.path(Dir.COV, "Soil_AK.nc")) | !file.exists(file.path(Dir.COV, "Soil_UK.nc"))){
    SoilCovs_stack <- stack(SoilCovs_ls)
    SoilCovs_AK <- crop(SoilCovs_stack, AK_shp)
    SoilCovs_AK <- mask(SoilCovs_AK, AK_shp)
    writeRaster(x = SoilCovs_AK, filename = file.path(Dir.COV, "Soil_AK"), format = "CDF")
    SoilCovs_UK <- crop(SoilCovs_stack, UK_shp)
    SoilCovs_UK <- mask(SoilCovs_UK, UK_shp)
    writeRaster(x = SoilCovs_UK, filename = file.path(Dir.COV, "Soil_UK"), format = "CDF")
  }else{
    SoilCovs_AK <- stack( file.path(Dir.COV, "Soil_AK.nc"))
    SoilCovs_UK <- stack( file.path(Dir.COV, "Soil_UK.nc"))
  }
}

#### .   TOPOGRAPHY COVARIATE RETRIEVAL (for kriging with other than DEM) -----------------------------------------------------------
TopoCovs_vec <- c(paste0("SlopesCl", 1:8), paste0("Aspect", c("ClN", "ClE", "ClS", "ClW"))) # need these names for topo covariates
TopoNames_vec <- c("Slopes1", "Slopes2", "Slopes3", "Slopes4", "Slopes5", "Slopes6", "Slopes7", "Slopes8","Slope_aspect_N", "Slope_aspect_E", "Slope_aspect_S", "Slope_aspect_W")
if(sum(file.exists(c(file.path(Dir.COV, "Covariates_UK.rds"), file.path(Dir.COV, "Covariates_AK.rds")))) != 2){
  print("#### Loading TOPOGRAPHY covariate data. ####") # documentation of these can be found here http://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/harmonized-world-soil-database-v12/en/
  # create lists to combine soil data into one
  TopoCovs_ls <- as.list(rep(NA, length(TopoCovs_vec)))
  names(TopoCovs_ls) <- c(TopoNames_vec)
  ## Downloading, unpacking, and stacking
for(Topo_Iter in TopoCovs_vec){
  if(!file.exists(file.path(Dir.COV, paste0(Topo_Iter, ".nc")))) { # if not downloaded and processed yet
    print(paste("Handling", Topo_Iter, "data."))
    Dir.Soil <- file.path(Dir.COV, Topo_Iter)
    dir.create(Dir.Soil)
    ## I used this function to obtain the rar archives from the data base, extracted the data by hand, and recompressed them into .zip for retrieval from a google drive for this script to work without forcing the user to step outside of R
    # httr::GET(paste0("http://www.fao.org/fileadmin/user_upload/soils/HWSD%20Viewer/Glo", Topo_Iter, "_30as.rar?raw=true"),
    #           write_disk(file.path(Dir.Soil, paste0(Topo_Iter, ".rar"))),
    #           progress(), overwrite = TRUE) # download data
    httr::GET("https://www.dropbox.com/sh/dj4avlpkx114pvt/AAAeb5SfGdI2F7lwdrbHS_s6a?raw=true",
              write_disk(file.path(Dir.COV, "HWSD.zip")),
              progress(), overwrite = TRUE) # download data
    unzip(file.path(Dir.COV, "HWSD.zip"), exdir = Dir.Soil, files = paste0(Topo_Iter, ".zip")) # unzip data
    unzip(file.path(Dir.Soil, paste0(Topo_Iter, ".zip")), exdir = Dir.Soil) # unzip data
    File <- list.files(Dir.Soil, pattern = ".asc")[1] # only keep first soil layer
    Soil_ras <- raster(file.path(Dir.Soil, File)) # load data
    TopoCovs_ls[[which(TopoCovs_vec == Topo_Iter)]] <- Soil_ras # save to list
    writeRaster(x = Soil_ras, filename = file.path(Dir.COV, Topo_Iter), format = "CDF")
    plot(Soil_ras, main = Topo_Iter, colNA = "black")
    unlink(Dir.Soil, recursive = TRUE)
  }else{
    print(paste(Topo_Iter, "already downloaded and processed."))
    TopoCovs_ls[[which(TopoCovs_vec == Topo_Iter)]] <- raster(file.path(Dir.COV, paste0(Topo_Iter, ".nc")))
  }
  unlink(file.path(Dir.COV, "HWSD.zip"))
}

if(!file.exists(file.path(Dir.COV, "Topo_AK.nc")) | !file.exists(file.path(Dir.COV, "Topo_UK.nc"))){
  TopoCovs_stack <- stack(TopoCovs_ls)
  TopoCovs_AK <- crop(TopoCovs_stack, AK_shp)
  TopoCovs_AK <- mask(TopoCovs_AK, AK_shp)
  TopoCovs_AK[[13]] <- TopoCovs_AK[[1]]+
    TopoCovs_AK[[2]]*2+TopoCovs_AK[[3]]*3+
    TopoCovs_AK[[4]]*4+TopoCovs_AK[[5]]*5+
    TopoCovs_AK[[6]]*6+TopoCovs_AK[[7]]*7+
    TopoCovs_AK[[8]]*8
  names(TopoCovs_AK)[13] <- "Slopes"
  writeRaster(x = TopoCovs_AK, filename = file.path(Dir.COV, "Topo_AK"), format = "CDF")
  TopoCovs_UK <- crop(TopoCovs_stack, UK_shp)
  TopoCovs_UK <- mask(TopoCovs_UK, UK_shp)
  TopoCovs_UK[[13]] <- TopoCovs_UK[[1]]+
    TopoCovs_UK[[2]]*2+TopoCovs_UK[[3]]*3+
    TopoCovs_UK[[4]]*4+TopoCovs_UK[[5]]*5+
    TopoCovs_UK[[6]]*6+TopoCovs_UK[[7]]*7+
    TopoCovs_UK[[8]]*8
  names(TopoCovs_UK)[13] <- "Slopes"
  writeRaster(x = TopoCovs_UK, filename = file.path(Dir.COV, "Topo_UK"), format = "CDF")
}else{
  TopoCovs_AK <- stack(file.path(Dir.COV, "Topo_AK.nc"))
  TopoCovs_UK <- stack(file.path(Dir.COV, "Topo_UK.nc"))
}
}

#### .   ERA5-LAND REFERENCES (for resampling covariates data) -----------------------------------------------------------
if(sum(file.exists(c(file.path(Dir.COV, "Covariates_UK.rds"), file.path(Dir.COV, "Covariates_AK.rds")))) != 2){
if(!file.exists(file.path(Dir.Figures, "Ref_UK.nc"))){
  REF_UK <- download_ERA(
    Variable = "2m_temperature",
    DateStart = "1984-10-01",
    DateStop = "1984-10-31",
    TResolution = "month",
    TStep = 1,
    Extent = UK_shp,
    Dir = Dir.Figures,
    FileName = "Ref_UK",
    API_Key = API_Key,
    API_User = API_User
  ) 
}else{
  REF_UK <- raster(file.path(Dir.Figures, "Ref_UK.nc"))
}
if(!file.exists(file.path(Dir.Figures, "Ref_AK.nc"))){
  REF_AK <- download_ERA(
    Variable = "2m_temperature",
    DateStart = "1984-10-01",
    DateStop = "1984-10-31",
    TResolution = "month",
    TStep = 1,
    Extent = AK_shp,
    Dir = Dir.Figures,
    FileName = "Ref_AK",
    API_Key = API_Key,
    API_User = API_User
  ) 
}else{
  REF_AK <- raster(file.path(Dir.Figures, "Ref_AK.nc"))
}
}
#### .   GMTED 2010 ELEVATION DATA (for kriging) -----------------------------------------------------------
if(sum(file.exists(c(file.path(Dir.COV, "Covariates_UK.rds"), file.path(Dir.COV, "Covariates_AK.rds")))) != 2){
  GMTED_UK <- download_DEM(
    Train_ras = REF_UK,
    Target_res = TopoCovs_UK[[1]],
    Shape = UK_shp,
    Dir = Dir.COV,
    Keep_Temporary = TRUE
  )
  GMTED_AK <- download_DEM(
    Train_ras = REF_AK,
    Target_res = TopoCovs_AK[[1]],
    Shape = AK_shp,
    Dir = Dir.COV,
    Keep_Temporary = TRUE
  )
}
#### .   COVARIATE RASTERS (combining covariate data for each region) -----------------------------------------------------------
if(sum(file.exists(c(file.path(Dir.COV, "Covariates_UK.rds"), file.path(Dir.COV, "Covariates_AK.rds")))) != 2){
  #### TARGET RESOLUTION (30 arc seconds)
  Covs_UK <- list(stack(SoilCovs_UK, TopoCovs_UK, GMTED_UK[[2]]))
  names(Covs_UK[[1]]) <- c(SoilCovs_vec, TopoNames_vec, "Slopes", "DEM")
  Covs_AK <- list(stack(SoilCovs_AK, TopoCovs_AK, GMTED_AK[[2]]))
  names(Covs_AK[[1]]) <- c(SoilCovs_vec, TopoNames_vec, "Slopes", "DEM")
  
  #### TRAINING RESOLUTION (era5-land)
  Covs_UK[[2]] <- stack(resample(x = stack(SoilCovs_UK, TopoCovs_UK), y = GMTED_UK[[1]]), GMTED_UK[[1]])
  names(Covs_UK[[2]]) <- c(SoilCovs_vec, TopoNames_vec, "Slopes", "DEM")
  Covs_AK[[2]] <- stack(resample(x = stack(SoilCovs_AK, TopoCovs_AK), y = GMTED_AK[[1]]), GMTED_AK[[1]])
  names(Covs_AK[[2]]) <- c(SoilCovs_vec, TopoNames_vec, "Slopes", "DEM")
  
  #### SAVING DATA FOR EASY RETRIEVAL
  saveRDS(object = Covs_UK, file = file.path(Dir.COV, "Covariates_UK.rds"))
  saveRDS(object = Covs_AK, file = file.path(Dir.COV, "Covariates_AK.rds"))
}else{
  Covs_UK <- readRDS(file.path(Dir.COV, "Covariates_UK.rds"))
  Covs_AK <- readRDS(file.path(Dir.COV, "Covariates_AK.rds"))
}

#### [FIGURE 1] (The effect of different Co-Variates on downscaling) -----------------------------------------------------------
Dir.Fig1 <- file.path(Dir.Figures, "Figure1")
if(!dir.exists(Dir.Fig1)){dir.create(Dir.Fig1)}

#### .   TEMPERATURE DATA -----------------------------------------------------------
if(!file.exists(file.path(Dir.Fig1, "ERA5-Land_SAT.nc"))){
  Era5Land_SAT <- download_ERA(
    Variable = "2m_temperature",
    DateStart = "1984-10-01",
    DateStop = "1984-10-31",
    TResolution = "month",
    TStep = 1,
    Extent = UK_shp,
    Dir = Dir.Fig1,
    FileName = "ERA5-Land_SAT",
    API_Key = API_Key,
    API_User = API_User
  )
}else{
  Era5Land_SAT <- raster(file.path(Dir.Fig1, "ERA5-Land_SAT.nc"))
}

if(sum(file.exists(c(file.path(Dir.Fig1, "Fig1A.nc"), file.path(Dir.Fig1, "SE_Fig1A.nc")))) != 2){
  Fig1A_ls <- krigR(
    Data = Era5Land_SAT,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = Covs_UK[[1]],
    KrigingEquation = "ERA ~ DEM",
    Cores = 1,
    Dir = Dir.Fig1,
    FileName = "Fig1A",
    Keep_Temporary = FALSE,
    nmax = 480
  ) 
}

if(sum(file.exists(c(file.path(Dir.Fig1, "Fig1B.nc"), file.path(Dir.Fig1, "SE_Fig1B.nc")))) != 2){
  Fig1B_ls <- krigR(
    Data = Era5Land_SAT,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = Covs_UK[[1]],
    KrigingEquation = "ERA ~ DEM+tkdry+tksat+csol+k_s+lambda+psi+theta_s",
    Cores = 1,
    Dir = Dir.Fig1,
    FileName = "Fig1B",
    Keep_Temporary = FALSE,
    nmax = 480
  ) 
}

if(sum(file.exists(c(file.path(Dir.Fig1, "Fig1C.nc"), file.path(Dir.Fig1, "SE_Fig1C.nc")))) != 2){
  Fig1C_ls <- krigR(
    Data = Era5Land_SAT,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = Covs_UK[[1]],
    KrigingEquation = "ERA ~ DEM+Slopes",
    Cores = 1,
    Dir = Dir.Fig1,
    FileName = "Fig1C",
    Keep_Temporary = FALSE,
    nmax = 480
  ) 
}

if(sum(file.exists(c(file.path(Dir.Fig1, "Fig1D.nc"), file.path(Dir.Fig1, "SE_Fig1D.nc")))) != 2){
  Fig1D_ls <- krigR(
    Data = Era5Land_SAT,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = Covs_UK[[1]],
    KrigingEquation = "ERA ~ DEM+AspectClS+ApsectClN+AspectClW+AspectClE",
    Cores = 1,
    Dir = Dir.Fig1,
    FileName = "Fig1D",
    Keep_Temporary = FALSE,
    nmax = 480
  ) 
}

#### .   SOIL MOISTURE DATA -----------------------------------------------------------
if(!file.exists(file.path(Dir.Fig1, "ERA5-Land_QSOIL.nc"))){
  Era5Land_QSOIL <- download_ERA(
    Variable = "volumetric_soil_water_layer_1",
    DateStart = "1984-10-01",
    DateStop = "1984-10-31",
    TResolution = "month",
    TStep = 1,
    Extent = UK_shp,
    Dir = Dir.Fig1,
    FileName = "ERA5-Land_QSOIL",
    API_Key = API_Key,
    API_User = API_User
  )
}else{
  Era5Land_QSOIL <- raster(file.path(Dir.Fig1, "ERA5-Land_QSOIL.nc"))
}

if(sum(file.exists(c(file.path(Dir.Fig1, "Fig1E.nc"), file.path(Dir.Fig1, "SE_Fig1E.nc")))) != 2){
  Fig1A_ls <- krigR(
    Data = Era5Land_QSOIL,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = Covs_UK[[1]],
    KrigingEquation = "ERA ~ tkdry+tksat+csol+k_s+lambda+psi+theta_s",
    Cores = 1,
    Dir = Dir.Fig1,
    FileName = "Fig1E",
    Keep_Temporary = FALSE,
    nmax = 480
  ) 
}

if(sum(file.exists(c(file.path(Dir.Fig1, "Fig1F.nc"), file.path(Dir.Fig1, "SE_Fig1F.nc")))) != 2){
  Fig1B_ls <- krigR(
    Data = Era5Land_QSOIL,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = Covs_UK[[1]],
    KrigingEquation = "ERA ~ DEM+tkdry+tksat+csol+k_s+lambda+psi+theta_s",
    Cores = 1,
    Dir = Dir.Fig1,
    FileName = "Fig1F",
    Keep_Temporary = FALSE,
    nmax = 480
  ) 
}

if(sum(file.exists(c(file.path(Dir.Fig1, "Fig1G.nc"), file.path(Dir.Fig1, "SE_Fig1G.nc")))) != 2){
  Fig1C_ls <- krigR(
    Data = Era5Land_QSOIL,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = Covs_UK[[1]],
    KrigingEquation = "ERA ~ tkdry+tksat+csol+k_s+lambda+psi+theta_s+Slopes",
    Cores = 1,
    Dir = Dir.Fig1,
    FileName = "Fig1G",
    Keep_Temporary = FALSE,
    nmax = 480
  ) 
}

if(sum(file.exists(c(file.path(Dir.Fig1, "Fig1H.nc"), file.path(Dir.Fig1, "SE_Fig1H.nc")))) != 2){
  Fig1D_ls <- krigR(
    Data = Era5Land_QSOIL,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = Covs_UK[[1]],
    KrigingEquation = "ERA ~ tkdry+tksat+csol+k_s+lambda+psi+theta_s+AspectClS+ApsectClN+AspectClW+AspectClE",
    Cores = 1,
    Dir = Dir.Fig1,
    FileName = "Fig1H",
    Keep_Temporary = FALSE,
    nmax = 480
  ) 
}

#### [FIGURE 2] (Dynamical vs. Statistical Uncertainty) -----------------------------------------------------------
Dir.Fig2 <- file.path(Dir.Figures, "Figure2")
if(!dir.exists(Dir.Fig2)){dir.create(Dir.Fig2)}

#### .   DYNAMICAL -----------------------------------------------------------

# FigA) 9pm 8th September 1984, ensemble members and then sd of all layers
# FigD) 9pm 8th September 1984, ensemble members and then sd of all layers

# FigC+F) 
## Dynamical: 
### Hourly resolution: Download start of data until 31st 12 2000 at hourly resolution, sd across ensemble at each time-step, mean of all time-steps
### Average Monthly ensemble members to to monthly ensemble members, then take sd across ensemble at each time step, mean of all time-steps
### Yearly ensemble mean first per ensemble member, etc. as for monthly


#### .   STATISTICAL -----------------------------------------------------------

# SAT downscaled with GMTED2010
# Qsoil1 downscaled with CHinaData

# FigB) 9pm 8th September 1984, era5-land SAT downscaled to 30arc-sec using elevation and nmax = 480
# FigE) 9pm 8th September 1984, era5-land QSOIL1 downscaled to 30arc-sec using ChinaData and nmax = 480

# FigC+F) 
## Statistical: 
### Hour resolution: 15th January/April/July/October 1200 and 0000 of the year 1984, downscale each of these and average the uncertainty
### Month: Monthly product of January/April/July/October 1984, downscale each and average the uncertainty
### Year: Just 1984 downscale (Just one time-step!)
### Climatology: 1981-2000, aggregate to full climatology with download_ERA() function and then downscale (Just one time-step!)

### Average Monthly ensemble members to to monthly ensemble members, then take sd across ensemble at each time step, mean of all time-steps
### Yearly ensemble mean first per ensemble member, etc. as for monthly

#### [FIGURE 3] (Kriged Products vs. Competitor Climate Products for all months Jan/1981-Dec/2010) -----------------------------------------------------------
Dir.Fig3 <- file.path(Dir.Figures, "Figure3")
if(!dir.exists(Dir.Fig3)){dir.create(Dir.Fig3)}
options(timeout=500)

#### .   ERA5-LAND -----------------------------------------------------------
if(!file.exists(file.path(Dir.Figures, "Era5Land_UK.nc"))){
  Era5Land_UK <- download_ERA(
    Variable = "2m_temperature",
    DateStart = "1981-01-01",
    DateStop = "2000-12-31",
    TResolution = "month",
    TStep = 1,
    Extent = UK_shp,
    Dir = Dir.Figures,
    FileName = "Era5Land_UK",
    API_Key = API_Key,
    API_User = API_User
  ) 
}else{
  Era5Land_UK <- raster(file.path(Dir.Figures, "Era5Land_UK.nc"))
}
if(!file.exists(file.path(Dir.Figures, "Era5Land_AK.nc"))){
  REF_AK <- download_ERA(
    Variable = "2m_temperature",
    DateStart = "1984-10-01",
    DateStop = "1984-10-31",
    TResolution = "month",
    TStep = 1,
    Extent = AK_shp,
    Dir = Dir.Figures,
    FileName = "Era5Land_AK",
    API_Key = API_Key,
    API_User = API_User
  ) 
}else{
  Era5Land_AK <- raster(file.path(Dir.Figures, "Era5Land_AK.nc"))
}

#### .   TERRACLIMATE -----------------------------------------------------------
### TerraClimate 1958 - 2015, monthly, 4x4km ---
## (https://www.northwestknowledge.net/data/5956e20ceb1bc6513f464d11/unzipped/TERRACLIMATE/ & https://climate.northwestknowledge.net/TERRACLIMATE/index_directDownloads.php)
## NOTE: Yes, no average temperature. Weird, I know. Daily records available for Europe only here: https://www.envidat.ch/#/metadata/eur11
Dir.TC <- file.path(Dir.Fig3, "TerraClimate")

if(sum(file.exists(c(file.path(Dir.TC, paste("TC", "tmean", "AK.nc", sep="_")), file.path(Dir.TC, paste("TC", "tmean", "UK.nc", sep="_"))))) != 2){
  dir.create(Dir.TC)
  Vars <- c("tmin", "tmax")
  dates <- 1981:2000
  for(i in 1:length(Vars)){
    if(file.exists(file.path(Dir.TC, paste0("TC_", Vars[i], ".nc")))){next()}
    Dir.Iter <- file.path(Dir.TC, Vars[i])
    dir.create(Dir.Iter)
    ## MONTHLY
    for(k in 1:length(dates)){
      URL <- paste0("https://climate.northwestknowledge.net/TERRACLIMATE-DATA/TerraClimate_", Vars[i], "_", dates[k], ".nc")
      if(!file.exists(file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".nc")))){
        download.file(URL, destfile = file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".nc")), method="libcurl", mode = "wb")
      }
    }
    setwd(Dir.Iter)
    BRICK <- stack(list.files(Dir.Iter, pattern = ".nc"))
    setwd(Dir.TC)
    TC_AK <- crop(BRICK, AK_shp)
    TC_AK <- mask(TC_AK, AK_shp)
    writeRaster(x = TC_AK, filename = file.path(Dir.TC, paste("TC", Vars[i], "AK", sep="_")), format = "CDF")
    TC_UK <- crop(BRICK, UK_shp)
    TC_UK <- mask(TC_UK, UK_shp)
    writeRaster(x = TC_UK, filename = file.path(Dir.TC, paste("TC", Vars[i], "UK", sep="_")), format = "CDF")
    unlink(Dir.Iter, recursive = TRUE)
  }
  ## Creating mean rasters
  AK_fs <- list.files(path = Dir.TC, pattern = "AK")
  AK_max <- stack(file.path(Dir.TC, AK_fs[[1]]))
  AK_min <- stack(file.path(Dir.TC, AK_fs[[2]]))
  AK_mean <- stack(AK_max, AK_min)
  ID <- rep(1:240, 2)
  AK_mean <- stackApply(x = AK_mean, indices = ID, fun = mean)
  writeRaster(x = AK_mean, filename = file.path(Dir.TC, paste("TC", "tmean", "AK", sep="_")), format = "CDF")
  unlink(file.path(Dir.TC, AK_fs))
  UK_fs <- list.files(path = Dir.TC, pattern = "UK")
  UK_max <- stack(file.path(Dir.TC, UK_fs[[1]]))
  UK_min <- stack(file.path(Dir.TC, UK_fs[[2]]))
  UK_mean <- stack(UK_max, UK_min)
  ID <- rep(1:240, 2)
  UK_mean <- stackApply(x = UK_mean, indices = ID, fun = mean)
  writeRaster(x = UK_mean, filename = file.path(Dir.TC, paste("TC", "tmean", "UK", sep="_")), format = "CDF")
  unlink(file.path(Dir.TC, UK_fs))
}
TC_AK <- stack(file.path(Dir.TC, paste("TC", "tmean", "AK.nc", sep="_")))
TC_UK <- stack(file.path(Dir.TC, paste("TC", "tmean", "UK.nc", sep="_")))

if(sum(file.exists(c(file.path(Dir.Fig3, "Fig3A.nc"), file.path(Dir.Fig3, "SE_Fig3A.nc")))) != 2){
  Fig3A_ls <- krigR(
    Data = Era5Land_UK,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = resample(Covs_UK$DEM, TC_UK),
    KrigingEquation = "ERA ~ DEM",
    Cores = numberOfCores,
    Dir = Dir.Fig3,
    FileName = "Fig3A",
    Keep_Temporary = FALSE,
    nmax = 480
  )
}

if(sum(file.exists(c(file.path(Dir.Fig3, "Fig3D.nc"), file.path(Dir.Fig3, "SE_Fig3D.nc")))) != 2){
  Fig3D_ls <- krigR(
    Data = Era5Land_AK,
    Covariates_coarse = Covs_AK[[2]],
    Covariates_fine = resample(Covs_AK$DEM, TC_AK),
    KrigingEquation = "ERA ~ DEM",
    Cores = numberOfCores,
    Dir = Dir.Fig3,
    FileName = "Fig3D",
    Keep_Temporary = FALSE,
    nmax = 480
  )
}

#### .   CHELSA -----------------------------------------------------------
### CHELSA 1979-2013, monthly, 1x1km ---
## (https://chelsa-climate.org/timeseries/)
## NOTE: The technical specifications list soil water, but the download doesn't present it
Dir.CH <- file.path(Dir.Fig3, "CHELSA")
if(!dir.exists(Dir.CH)){dir.create(Dir.CH)}

## MONTHLY
Vars <- c("tmean")
dates <- paste(rep(1981:2000, each = 12), str_pad(string = 1:12, width = 2, "left", 0), sep="_")
for(i in 1:length(Vars)){
  if(sum(file.exists(c(file.path(Dir.CH, paste("CH", Vars[i], "AK.nc", sep="_")), file.path(Dir.CH, paste("CH", Vars[i], "UK.nc", sep="_"))))) == 2){
    CHELSA_UK <- stack(file.path(Dir.CH, paste("CH", Vars[i], "UK.nc", sep="_")))
    CHELSA_AK <- stack(file.path(Dir.CH, paste("CH", Vars[i], "AK.nc", sep="_")))
    next()
  }
  Dir.Iter <- file.path(Dir.CH, Vars[i])
  dir.create(Dir.Iter)
  for(k in 1:length(dates)){
    URL <- paste0("https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/timeseries/", Vars[i],"/CHELSA_", Vars[i],"_", dates[k], "_V1.2.1.tif")
    if(!file.exists(file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".tif")))){
      download.file(URL, destfile = file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".tif")), mode="wb") # download cultural vector
    }
  }
  setwd(Dir.Iter)
  BRICK <- stack(list.files(Dir.Iter, pattern = ".tif"))
  setwd(Dir.CH)
  CHELSA_AK <- crop(BRICK, AK_shp)
  CHELSA_AK <- mask(CHELSA_AK, AK_shp)/10-273.15 # CHELSA reported as 10*C, conver to K
  writeRaster(x = CHELSA_AK, filename = file.path(Dir.CH, paste("CH", Vars[i], "AK", sep="_")), format = "CDF")
  CHELSA_UK <- crop(BRICK, UK_shp)
  CHELSA_UK <- mask(CHELSA_UK, UK_shp)/10-273.15 # CHELSA reported as 10*C, conver to K
  writeRaster(x = CHELSA_UK, filename = file.path(Dir.CH, paste("CH", Vars[i], "UK", sep="_")), format = "CDF")
  unlink(Dir.Iter, recursive = TRUE)
}

if(sum(file.exists(c(file.path(Dir.Fig3, "Fig3B.nc"), file.path(Dir.Fig3, "SE_Fig3B.nc")))) != 2){
  Fig3B_ls <- krigR(
    Data = Era5Land_UK,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = resample(Covs_UK$DEM, CHELSA_UK),
    KrigingEquation = "ERA ~ DEM",
    Cores = numberOfCores,
    Dir = Dir.Fig3,
    FileName = "Fig3B",
    Keep_Temporary = FALSE,
    nmax = 480
  )
}

if(sum(file.exists(c(file.path(Dir.Fig3, "Fig3E.nc"), file.path(Dir.Fig3, "SE_Fig3E.nc")))) != 2){
  Fig3E_ls <- krigR(
    Data = Era5Land_AK,
    Covariates_coarse = Covs_AK[[2]],
    Covariates_fine = resample(Covs_AK$DEM, CHELSA_AK),
    KrigingEquation = "ERA ~ DEM",
    Cores = numberOfCores,
    Dir = Dir.Fig3,
    FileName = "Fig3E",
    Keep_Temporary = FALSE,
    nmax = 480
  )
}

#### .   WORLDCLIM -----------------------------------------------------------
### WORLDCLIM 1960-2018/1970-2000, monthly, 5x5km ---
## (https://www.worldclim.org/data/monthlywth.html)
## It's either 1x1km of climatologies or 2.5minutes (5x5km)
Dir.WC <- file.path(Dir.Fig3, "WorldClim")

if(sum(file.exists(c(file.path(Dir.WC, paste("WC", "tmean", "AK.nc", sep="_")), file.path(Dir.WC, paste("WC", "tmean", "UK.nc", sep="_"))))) != 2){
  dir.create(Dir.WC)
  ## MONTHLY DATA
  Vars <- c("tmin", "tmax")
  dates <- paste(seq(1980, 2000, 10), seq(1989, 2009, 10), sep="-")
  for(i in 1:length(Vars)){
    if(file.exists(file.path(Dir.WC, paste0("WC_", Vars[i], "_UK.nc")))){next()}
    Dir.Iter <- file.path(Dir.WC, Vars[i])
    dir.create(Dir.Iter)
    for(k in 1:length(dates)){
      URL <- paste0("https://data.biogeo.ucdavis.edu/data/worldclim/v2.1/hist/wc2.1_2.5m_", Vars[i], "_", dates[k], ".zip")
      if(!file.exists(file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".zip")))){
        Download <- FALSE
        while(Download == FALSE){
          try(download.file(URL, destfile = file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".zip")))) # download cultural vector
          if(file.exists(file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".zip")))){
            Download <- TRUE  
          }
        }
      }
      if(length(list.files(Dir.Iter, pattern = ".tif")) < 12*9*k){ 
        unzip(file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".zip")), exdir = Dir.Iter) # unzip data
      }
    }
    setwd(Dir.Iter)
    BRICK <- stack(list.files(Dir.Iter, pattern = ".tif"))
    setwd(Dir.WC)
    BRICK <- BRICK[[-1:-12]]
    BRICK <- BRICK[[1:240]]
    WC_AK <- crop(BRICK, AK_shp)
    WC_AK <- mask(WC_AK, AK_shp)
    writeRaster(x = WC_AK, filename = file.path(Dir.WC, paste("WC", Vars[i], "AK", sep="_")), format = "CDF")
    WC_UK <- crop(BRICK, UK_shp)
    WC_UK <- mask(WC_UK, UK_shp)
    writeRaster(x = WC_UK, filename = file.path(Dir.WC, paste("WC", Vars[i], "UK", sep="_")), format = "CDF")
    unlink(Dir.Iter, recursive = TRUE)
  }
  ## Creating mean rasters
  AK_fs <- list.files(path = Dir.WC, pattern = "AK")
  AK_max <- stack(file.path(Dir.WC, AK_fs[[1]]))
  AK_min <- stack(file.path(Dir.WC, AK_fs[[2]]))
  AK_mean <- stack(AK_max, AK_min)
  ID <- rep(1:240, 2)
  AK_mean <- stackApply(x = AK_mean, indices = ID, fun = mean)
  writeRaster(x = AK_mean, filename = file.path(Dir.WC, paste("WC", "tmean", "AK", sep="_")), format = "CDF")
  unlink(file.path(Dir.WC, AK_fs))
  UK_fs <- list.files(path = Dir.WC, pattern = "UK")
  UK_max <- stack(file.path(Dir.WC, UK_fs[[1]]))
  UK_min <- stack(file.path(Dir.WC, UK_fs[[2]]))
  UK_mean <- stack(UK_max, UK_min)
  ID <- rep(1:240, 2)
  UK_mean <- stackApply(x = UK_mean, indices = ID, fun = mean)
  writeRaster(x = UK_mean, filename = file.path(Dir.WC, paste("WC", "tmean", "UK", sep="_")), format = "CDF")
  unlink(file.path(Dir.WC, UK_fs))
}
WC_AK <- stack(file.path(Dir.WC, paste("WC", "tmean", "AK.nc", sep="_")))
WC_UK <- stack(file.path(Dir.WC, paste("WC", "tmean", "UK.nc", sep="_")))

if(sum(file.exists(c(file.path(Dir.Fig3, "Fig3C.nc"), file.path(Dir.Fig3, "SE_Fig3C.nc")))) != 2){
  Fig3C_ls <- krigR(
    Data = Era5Land_UK,
    Covariates_coarse = Covs_UK[[2]],
    Covariates_fine = resample(Covs_UK$DEM, WC_UK),
    KrigingEquation = "ERA ~ DEM",
    Cores = numberOfCores,
    Dir = Dir.Fig3,
    FileName = "Fig3C",
    Keep_Temporary = FALSE,
    nmax = 480
  )
}

if(sum(file.exists(c(file.path(Dir.Fig3, "Fig3F.nc"), file.path(Dir.Fig3, "SE_Fig3F.nc")))) != 2){
  Fig3F_ls <- krigR(
    Data = Era5Land_AK,
    Covariates_coarse = Covs_AK[[2]],
    Covariates_fine = resample(Covs_AK$DEM, WC_AK),
    KrigingEquation = "ERA ~ DEM",
    Cores = numberOfCores,
    Dir = Dir.Fig3,
    FileName = "Fig3F",
    Keep_Temporary = FALSE,
    nmax = 480
  )
}

#### [FIGURE S2] (Localised Kriging and the nmax-trade-off) -----------------------------------------------------------
Dir.FigS2 <- file.path(Dir.Figures, "FigureS2")
if(!dir.exists(Dir.FigS2)){dir.create(Dir.FigS2)}

if(!file.exists(file.path(Dir.FigS2, "ERA5-Land_SAT.nc"))){
  Era5Land_SAT <- download_ERA(
    Variable = "2m_temperature",
    DateStart = "1984-10-01",
    DateStop = "1984-10-31",
    TResolution = "month",
    TStep = 1,
    Extent = UK_shp,
    Dir = Dir.FigS2,
    FileName = "ERA5-Land_SAT",
    API_Key = API_Key,
    API_User = API_User
  )
}else{
  Era5Land_SAT <- raster(file.path(Dir.FigS2, "ERA5-Land_SAT.nc"))
}

if(!file.exists(file.path(Dir.FigS2, "Time_vec.rds"))){
  Nmax_vec <- c(15, 20, 30, 40, 50, 60, 70, 80, 110, 140, 170, 200, 260, 320, 400, 480)
  Time_vec <- rep(NA, length(Nmax_vec))
  counter <- 1
  for(Nmax_Iter in Nmax_vec){
    Begin_time <- Sys.time()
    krigR(
      Data = Era5Land_SAT,
      Covariates_coarse = Covs_UK[[2]],
      Covariates_fine = Covs_UK[[1]],
      KrigingEquation = "ERA ~ DEM",
      Cores = 1,
      Dir = Dir.FigS2,
      FileName = paste0("FigS2_", Nmax_Iter, ".nc"),
      Keep_Temporary = FALSE,
      nmax = Nmax_Iter
    ) 
    End_time <- Sys.time()
    Time_vec[counter] <- End_time-Begin_time
    counter <- counter + 1
}
saveRDS(object = Time_vec, file = file.path(Dir.FigS2, "Time_vec.rds"))
}

#### [FIGURE S3] (Downscaling Uncertainty and Confidence in Kriged Products) -----------------------------------------------------------
Dir.FigS3 <- file.path(Dir.Figures, "FigureS3")
if(!dir.exists(Dir.FigS3)){dir.create(Dir.FigS3)}

# era5-land SAT for every month between 1981 and 2000, aggregate by a factor of and then downscale back down to original era5-land resolution
## SAT
## Qsoil1

###!!! what about other interpolation techniques?

