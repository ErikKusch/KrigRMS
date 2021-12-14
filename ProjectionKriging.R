# PREAMBLE ----------------------------------------------------
rm(list = ls()) # clearing environment
## PACKAGES ---------------------------------------------------
if("KrigR" %in% rownames(installed.packages()) == FALSE){ # KrigR check
  Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
  devtools::install_github("https://github.com/ErikKusch/KrigR")
}
library(KrigR) 

## KrigR API --------------------------------------------------
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

## COUNTRY MASK -----------------------------------------------
print("#### Loading COUNTRY MASK. ####")
if(!file.exists(file.path(Dir.Shapes, "CountryMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip", destfile = paste(Dir.Shapes, "CountryMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "CountryMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
CountryMask <- readOGR(Dir.Shapes, "ne_10m_admin_0_countries", verbose = FALSE) # read country mask in
DE_shp <-  CountryMask[CountryMask$NAME == "Germany", ]


# DATA DOWNLOAD -----------------------------------------------
## HISTORICAL DATA --------------------------------------------
if(!file.exists("historical_ERA_1981-2000.nc")){
  AT_ras <- download_ERA(Variable = "2m_temperature",
               DateStart = "1981-01-01",
               DateStop = "1999-12-31",
               TResolution = "month",
               TStep = 1,
               Extent = DE_shp,
               FileName = "historical_ERA_1981-2000", 
               API_Key = API_Key,
               API_User = API_User)
  Index <- rep(1:12, length = nlayers(AT_ras))
  ATClim_ras <- stackApply(AT_ras, indices = Index, fun = mean)
  writeRaster(ATClim_ras, filename = "historical_ERA_1981-2000.nc", format = "CDF")
}else{
  ATClim_ras <- stack("historical_ERA_1981-2000.nc")
}

## PROJECTION DATA --------------------------------------------
train_ERA <- ATClim_ras[[1]]
### SSP ----
train_SSP <- stack("ssp585_tas_2041-2060.nc")
train_SSP <- train_SSP[[1]]
train_SSP <- crop(train_SSP,extent(train_ERA))

### HISTORICAL ----
train_HIST <- stack("XXX")
train_HIST <- train_HIST[[1]]
train_HIST <- crop(train_HIST,extent(train_ERA))

# KRIGING -----------------------------------------------------

## HISTORICAL DATA --------------------------------------------
GMTED_DE <- download_DEM(
  Train_ras = train_ERA,
  Target_res = 0.008334,
  Keep_Temporary = TRUE
)
Cov_coarse <- GMTED_DE[[1]]
Cov_fine <- GMTED_DE[[2]]

Output_ERA <- krigR(
  Data = train_ERA,
  Covariates_coarse = Cov_coarse, 
  Covariates_fine = Cov_fine,   
  KrigingEquation = "ERA ~ DEM",  
  Cores = 1, 
  Dir = getwd(),  
  FileName = "DE_hist_nmax120", 
  Keep_Temporary = FALSE,
  nmax = 120
)

## PROJECTION DATA --------------------------------------------
### SSP ----
GMTED_DE <- download_DEM(
  Train_ras = train_SSP,
  Target_res = 0.008334,
  Keep_Temporary = TRUE
)
Cov_coarse <- GMTED_DE[[1]]
Cov_fine <- GMTED_DE[[2]]

Output_SSP <- krigR(
  Data = train_SSP,
  Covariates_coarse = Cov_coarse, 
  Covariates_fine = Cov_fine,   
  KrigingEquation = "ERA ~ DEM",  
  Cores = 1, 
  Dir = getwd(),  
  FileName = "DE_SSP585_2041-2060_nmax120", 
  Keep_Temporary = FALSE,
  nmax = 120
)

### HISTORICAL ----
GMTED_DE <- download_DEM(
  Train_ras = train_HIST,
  Target_res = 0.008334,
  Keep_Temporary = TRUE
)
Cov_coarse <- GMTED_DE[[1]]
Cov_fine <- GMTED_DE[[2]]

Output_SSP <- krigR(
  Data = train_HIST,
  Covariates_coarse = Cov_coarse, 
  Covariates_fine = Cov_fine,   
  KrigingEquation = "ERA ~ DEM",  
  Cores = 1, 
  Dir = getwd(),  
  FileName = "DE_CMIP-HIST_nmax120", 
  Keep_Temporary = FALSE,
  nmax = 120
)