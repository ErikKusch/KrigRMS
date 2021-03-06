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
AK_shp <- crop(AK_shp,  c(-180.05, -128.95, 50.95, 72.05)) # crop out some Aleutian islands
plot(AK_shp)

####--------------- [1] COVARIATES ----------------------------------------------------
#### .   SOIL COVARIATE RETRIEVAL (for kriging with other than DEM) -----------------------------------------------------------
SoilCovs_vec <- c("tkdry", "tksat", "csol", "k_s", "lambda", "psi", "theta_s") # need these names for addressing soil covariates
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
SoilCovs_stack <- stack(SoilCovs_ls)

#### .   TOPOGRAPHY COVARIATE RETRIEVAL (for kriging with other than DEM) -----------------------------------------------------------
TopoCovs_vec <- c(paste0("SlopesCl", 1:8), paste0("Aspect", c("ClN", "ClE", "ClS", "ClW"))) # need these names for topo covariates
TopoNames_vec <- c("Slopes1", "Slopes2", "Slopes3", "Slopes4", "Slopes5", "Slopes6", "Slopes7", "Slopes8","Slope_aspect_N", "Slope_aspect_E", "Slope_aspect_S", "Slope_aspect_W")
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
    
    httr::GET(paste0("https://www.dropbox.com/s/aplojivr6rdanf5/", Topo_Iter, ".zip?raw=true"),
              write_disk(file.path(Dir.Soil, paste0(Topo_Iter, ".zip"))),
              progress(), overwrite = TRUE) # download data
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
}
SoilCovs_stack <- stack(TopoCovs_ls)



#### [FIGURE 1] (The effect of different Co-Variates on downscaling) -----------------------------------------------------------
## This is script called 
### SAT + QSOIL1 September 1984, with nmax 480
## --> SAT kriged with Elevation, Elevation+ChinaData, Elevation+Slopes (combination in script for figure 5), slope aspect (ClS, etc.)
## --> QSOIL1 kriged with ChinaData, Elevation+ChinaData, ChinaData+Slopes-Combination, ChinaData+Slopes

#### .   TEMPERATURE DATA -----------------------------------------------------------

#### .   SOIL MOISTURE DATA -----------------------------------------------------------

#### [FIGURE 2] (Dynamical vs. Statistical Uncertainty) -----------------------------------------------------------

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
options(timeout=500)

### downscale each individual month and store those

#### .   CHELSA -----------------------------------------------------------
### CHELSA 1979-2013, monthly, 1x1km ---
## (https://chelsa-climate.org/timeseries/)
## NOTE: The technical specifications list soil water, but the download doesn't present it
Dir.CH <- file.path(Dir.Figures, "CHELSA")
dir.create(Dir.CH)

## MONTHLY
Vars <- c("tmean")
dates <- paste(rep(1981:2000, each = 12), str_pad(string = 1:12, width = 2, "left", 0), sep="_")
for(i in 1:length(Vars)){
  if(file.exists(file.path(Dir.CH, paste0("CH_", Vars[i], sep=".nc")))){next()}
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
    BRICK <- BRICK/10-273.15 # CHELSA reported as 10*C, conver to K
    writeRaster(BRICK, paste("CH", Vars[i], sep="_"), format = "CDF")
    unlink(Dir.Iter, recursive = TRUE)
  }

#### .   TERRACLIMATE -----------------------------------------------------------
### TerraClimate 1958 - 2015, monthly, 4x4km ---
## (https://www.northwestknowledge.net/data/5956e20ceb1bc6513f464d11/unzipped/TERRACLIMATE/ & https://climate.northwestknowledge.net/TERRACLIMATE/index_directDownloads.php)
## NOTE: Yes, no average temperature. Weird, I know. Daily records available for Europe only here: https://www.envidat.ch/#/metadata/eur11
Dir.TC <- file.path(Dir.Figures, "TerraClimate")
dir.create(Dir.TC)

Vars <- c("tmin", "tmax")
dates <- 1981:2000
for(i in 1:length(Vars)){
  if(file.exists(file.path(Dir.TC, paste0("TC_", Vars[i], ".nc")))){next()}
  Dir.Iter <- file.path(Dir.TC, Vars[i])
  dir.create(Dir.Iter)
  ## MONTHLY
  for(k in 1:length(dates)){
    URL <- paste0("http://thredds.northwestknowledge.net:8080/thredds/fileServer/TERRACLIMATE_ALL/data/TerraClimate_", Vars[i], "_", dates[k], ".nc")
    if(!file.exists(file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".nc")))){
      download.file(URL, destfile = file.path(Dir.Iter, paste0(Vars[i],"_", dates[k], ".nc")), method="libcurl", mode = "wb")
    }
  }
  setwd(Dir.Iter)
  BRICK <- stack(list.files(Dir.Iter, pattern = ".nc"))
  setwd(Dir.TC)
  writeRaster(BRICK, paste("TC", Vars[i], sep="_"), format = "CDF")
  unlink(Dir.Iter, recursive = TRUE)
}

#### .   WORLDCLIM -----------------------------------------------------------
### WORLDCLIM 1960-2018/1970-2000, monthly, 5x5km ---
## (https://www.worldclim.org/data/monthlywth.html)
## It's either 1x1km of climatologies or 2.5minutes (5x5km)
Dir.WC <- file.path(Dir.Figures, "WorldClim")
dir.create(Dir.WC)

## MONTHLY DATA
Vars <- c("tmin", "tmax")
dates <- paste(seq(1980, 2000, 10), seq(1989, 2009, 10), sep="-")
for(i in 1:length(Vars)){
  if(file.exists(file.path(Dir.WC, paste0("WC_", Vars[i], ".nc")))){next()}
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
  
  stop("Remove year 1980 and years after the full year 2000")
  
  writeRaster(BRICK, paste("WC", Vars[i], sep="_"), format = "CDF")
  unlink(Dir.Iter, recursive = TRUE)
}

## CLIMATOLOGY
URL <- paste0("https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_tavg.zip")
if(!file.exists(file.path(Dir.WC, "wc2.1_30s_tavg_01.tif"))){
  Download <- FALSE
  while(Download == FALSE){
    try(download.file(URL, destfile = file.path(Dir.WC, "wc2.1_30s_tavg.zip")))
    if(file.exists(file.path(Dir.WC, "wc2.1_30s_tavg.zip"))){
      Download <- TRUE  
    }
  }
  unzip(file.path(Dir.WC, "wc2.1_30s_tavg.zip"), exdir = Dir.WC) # unzip data
  unlink(file.path(Dir.WC, "wc2.1_30s_tavg.zip"), recursive = TRUE)
}

#### .   KRIGR PRODUCT -----------------------------------------------------------


#### [FIGURE S2] (Localised Kriging and the nmax-trade-off) -----------------------------------------------------------

# UK SAT downscaled with elevation and nmax of 15, 20, 30, 40, 50, 60, 70, 80, 100. 110, 140, 170, 200, 260, 320, 400, 480

#### [FIGURE S3] (Downscaling Uncertainty and Confidence in Kriged Products) -----------------------------------------------------------

# era5-land SAT for every month between 1981 and 2000, aggregate by a factor of and then downscale back down to original era5-land resolution
## SAT
## Qsoil1Ñ

