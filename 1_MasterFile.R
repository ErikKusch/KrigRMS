rm(list=ls())
####### PACKAGES ---------------------------------------------------------
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
  numberOfCores <- readline(prompt = paste("How many cores do you want to allocate to these processes? Your machine has", parallel::detectCores()))
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

####### DIRECTORIES ---------------------------------------------------------
Dir.Base <- getwd() # read out the project directory
Dir.Shapes <- file.path(Dir.Base, "ShapeFiles")
if(!dir.exists(Dir.Shapes)){dir.create(Dir.Shapes)}
Dir.COV <- file.path(Dir.Base, "Covariates")
if(!dir.exists(Dir.COV)){dir.create(Dir.COV)}

####### FUNCTIONALITY ---------------------------------------------------------
`%nin%` <- Negate(`%in%`)

#### COUNTRY MASK (for producing maps with national borders) -----------------------------------------------------------
print("#### Loading COUNTRY MASK. ####")
if(!file.exists(file.path(Dir.Shapes, "CountryMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip", destfile = paste(Dir.Shapes, "CountryMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "CountryMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
CountryMask <- readOGR(Dir.Shapes, "ne_10m_admin_0_countries", verbose = FALSE) # read country mask in
UK_shp <- CountryMask[CountryMask$NAME == "United Kingdom", ] # extracting UK shape

#### STATE MASK (for producing maps with state borders) -----------------------------------------------------------
print("#### Loading STATE MASK. ####")
if(!file.exists(file.path(Dir.Shapes, "StateMask.zip"))){ # if land mask has not been downloaded yet
  download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_1_states_provinces.zip", destfile = paste(Dir.Shapes, "StateMask.zip", sep="/")) # download cultural vector
  unzip(paste(Dir.Shapes, "StateMask.zip", sep="/"), exdir = Dir.Shapes) # unzip the data
}
StateMask <- readOGR(Dir.Shapes, "ne_10m_admin_1_states_provinces", verbose = FALSE) # read state mask in
AK_shp <- StateMask[which(StateMask$name_en == "Alaska"), ] # extracting Alaska shape
AK_ext <- extent(AK_shp) # obtaining extent of Alaska shape
AK_ext[2] <- -129.98 # defining eastern-most point of Alaska we want to retain (some Aleutian islands are skipped in our analysis)
AK_shp <- crop(AK_shp,  AK_ext) # crop out some Aleutian islands
plot(AK_shp)

#### SOIL COVARIATE RETRIEVAL (for kriging with other than DEM) -----------------------------------------------------------
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

#### TOPOGRAPHY COVARIATE RETRIEVAL (for kriging with other than DEM) -----------------------------------------------------------
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
    TopoCovs_ls[[which(names(TopoCovs_ls) == Topo_Iter)]] <- Soil_ras # save to list
    writeRaster(x = Soil_ras, filename = file.path(Dir.COV, Topo_Iter), format = "CDF")
    plot(Soil_ras, main = Topo_Iter, colNA = "black")
    unlink(Dir.Soil, recursive = TRUE)
  }else{
    print(paste(Topo_Iter, "already downloaded and processed."))
    TopoCovs_ls[[which(names(TopoCovs_ls) == Topo_Iter)]] <- raster(file.path(Dir.COV, paste0(Topo_Iter, ".nc")))
  }
}
SoilCovs_stack <- stack(TopoCovs_ls)



#### FIGURE 1 (The effect of different Co-Variates on downscaling) -----------------------------------------------------------

#### FIGURE 2 (Dynamical vs. Statistical Uncertainty) -----------------------------------------------------------

#### FIGURE 3 (Kriged Products vs. Competitor Climate Products for all months Jan/1981-Dec/2010) -----------------------------------------------------------




#### FIGURE S2 (Localised Kriging and the nmax-trade-off) -----------------------------------------------------------

#### FIGURE S3 (Downscaling Uncertainty and Confidence in Kriged Products) -----------------------------------------------------------