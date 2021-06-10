krigRLocalSleep <- function(Data = NULL, Covariates_coarse = NULL, Covariates_fine = NULL, KrigingEquation = "ERA ~ DEM", Cores = detectCores(), Dir = getwd(), FileName, Keep_Temporary = TRUE, SingularTry = 10, Variable, PrecipFix = FALSE, Type = "reanalysis", DataSet = "era5-land", DateStart, DateStop, TResolution = "month", TStep = 1, FUN = 'mean', Extent, Buffer = 0.5, ID = "ID", API_Key, API_User, Target_res, Source = "USGS", nmax = Inf,  TryDown = 10, verbose = TRUE, ...){
  ## CALL LIST (for storing how the function as called in the output) ----
  if(is.null(Data)){
    Data_Retrieval <- list(Variable = Variable,
                           Type = Type,
                           PrecipFix = PrecipFix,
                           DataSet = DataSet,
                           DateStart = DateStart,
                           DateStop = DateStop,
                           TResolution = TResolution,
                           TStep = TStep,
                           Extent = Extent)
  }else{
    Data_Retrieval <- "None needed. Data was not queried via krigR function, but supplied by user."
  }
  ## CLIMATE DATA (call to download_ERA function if no Data set is specified) ----
  if(is.null(Data)){ # data check: if no data has been specified
    Data <- download_ERA(Variable = Variable, PrecipFix = PrecipFix, Type = Type, DataSet = DataSet, DateStart = DateStart, DateStop = DateStop, TResolution = TResolution, TStep = TStep, FUN = FUN, Extent = Extent, API_User = API_User, API_Key = API_Key, Dir = Dir, TryDown = TryDown, verbose = verbose)
  } # end of data check
  
  ## COVARIATE DATA (call to download_DEM function when no covariates are specified) ----
  if(is.null(Covariates_coarse) & is.null(Covariates_fine)){ # covariate check: if no covariates have been specified
    if(class(Extent) == "SpatialPolygonsDataFrame" | class(Extent) == "data.frame"){ # Extent check: if Extent input is a shapefile
      Shape <- Extent # save shapefile for use as Shape in masking covariate data
    }else{ # if Extent is not a shape, then extent specification is already baked into Data
      Shape <- NULL # set Shape to NULL so it is ignored in download_DEM function when masking is applied
    } # end of Extent check
    Covs_ls <- download_DEM(Train_ras = Data, Target_res = Target_res, Shape = Shape, Buffer = Buffer, ID = ID, Keep_Temporary = Keep_Temporary, Dir = Dir)
    Covariates_coarse <- Covs_ls[[1]] # extract coarse covariates from download_DEM output
    Covariates_fine <- Covs_ls[[2]] # extract fine covariates from download_DEM output
  } # end of covariate check
  
  ## KRIGING FORMULA (assure that KrigingEquation is a formula object) ----
  KrigingEquation <- as.formula(KrigingEquation)
  
  ## CALL LIST (for storing how the function as called in the output) ----
  Call_ls <- list(Data = KrigR:::SummarizeRaster(Data),
                  Covariates_coarse = KrigR:::SummarizeRaster(Covariates_coarse),
                  Covariates_fine = KrigR:::SummarizeRaster(Covariates_fine),
                  KrigingEquation = KrigingEquation,
                  Cores = Cores,
                  FileName = FileName,
                  Keep_Temporary = Keep_Temporary,
                  nmax = nmax,
                  Data_Retrieval = Data_Retrieval,
                  misc = ...)
  
  ## SANITY CHECKS (step into check_Krig function to catch most common error messages) ----
  Check_Product <- KrigR:::check_Krig(Data = Data, CovariatesCoarse = Covariates_coarse, CovariatesFine = Covariates_fine, KrigingEquation = KrigingEquation)
  KrigingEquation <- Check_Product[[1]] # extract KrigingEquation (this may have changed in check_Krig)
  DataSkips <- Check_Product[[2]] # extract which layers to skip due to missing data (this is unlikely to ever come into action)
  Terms <- unique(unlist(strsplit(labels(terms(KrigingEquation)), split = ":"))) # identify which layers of data are needed
  
  ## DATA REFORMATTING (Kriging requires spatially referenced data frames, reformatting from rasters happens here) ---
  Origin <- raster::as.data.frame(Covariates_coarse, xy = TRUE) # extract covariate layers
  Origin <- Origin[, c(1:2, which(colnames(Origin) %in% Terms))] # retain only columns containing terms
  
  Target <- raster::as.data.frame(Covariates_fine, xy = TRUE) # extract covariate layers
  Target <- Target[, c(1:2, which(colnames(Target) %in% Terms))] # retain only columns containing terms
  Target <- na.omit(Target)
  suppressWarnings(gridded(Target) <- ~x+y) # establish a gridded data product ready for use in kriging
  Target@grid@cellsize[1] <- Target@grid@cellsize[2] # ensure that grid cells are square
  
  ## SET-UP TEMPORARY DIRECTORY (this is where kriged products of each layer will be saved) ----
  Dir.Temp <- file.path(Dir, paste("Kriging", FileName, sep="_"))
  if(!dir.exists(Dir.Temp)){dir.create(Dir.Temp)}
  
  ## KRIGING SPECIFICATION (this will be parsed and evaluated in parallel and non-parallel evaluations further down) ----
  looptext <- "
  OriginK <- cbind(Origin, raster::extract(x = Data[[Iter_Krige]], y = Origin[,1:2], df=TRUE)[, 2]) # combine data of current data layer with training covariate data
  OriginK <- na.omit(OriginK) # get rid of NA cells
  colnames(OriginK)[length(Terms)+3] <- c(terms(KrigingEquation)[[2]]) # assign column names
  suppressWarnings(gridded(OriginK) <-  ~x+y) # generate gridded product
  OriginK@grid@cellsize[1] <- OriginK@grid@cellsize[2] # ensure that grid cells are square

  Iter_Try = 0 # number of tries set to 0
  kriging_result <- NULL
  while(class(kriging_result)[1] != 'autoKrige' & Iter_Try < SingularTry){ # try kriging SingularTry times, this is because of a random process of variogram identification within the automap package that can fail on smaller datasets randomly when it isn't supposed to
    try(invisible(capture.output(kriging_result <- autoKrige(formula = KrigingEquation, input_data = OriginK, new_data = Target, nmax = nmax))), silent = TRUE)
    Iter_Try <- Iter_Try +1
  }
  if(class(kriging_result)[1] != 'autoKrige'){ # give error if kriging fails
    print(paste0('Kriging failed for layer ', Iter_Krige, '. Error message produced by autoKrige function: ', geterrmessage()))
  }

  ## retransform to raster
  try( # try fastest way - this fails with certain edge artefacts in meractor projection and is fixed by using rasterize
    Krig_ras <- raster(x = kriging_result$krige_output, layer = 1), # extract raster from kriging product
    silent = TRUE
  )
  try(
    Var_ras <- raster(x = kriging_result$krige_output, layer = 3), # extract raster from kriging product
    silent = TRUE
  )
  if(!exists('Krig_ras') & !exists('Var_ras')){
    Krig_ras <- rasterize(x = kriging_result$krige_output, y = Covariates_fine[[1]])[[2]] # extract raster from kriging product
    Var_ras <- rasterize(x = kriging_result$krige_output, y = Covariates_fine)[[4]] # extract raster from kriging product
  }
  crs(Krig_ras) <- crs(Data) # setting the crs according to the data
  crs(Var_ras) <- crs(Data) # setting the crs according to the data

  if(Cores == 1){
  Ras_Krig[[Iter_Krige]] <- Krig_ras
  Ras_Var[[Iter_Krige]] <- Var_ras
  } # stack kriged raster into raster list if non-parallel computing
 writeRaster(x = Krig_ras, filename = file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_data.nc')), overwrite = TRUE, format='CDF') # save kriged raster to temporary directory
 writeRaster(x = Var_ras, filename = file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_SE.nc')), overwrite = TRUE, format='CDF') # save kriged raster to temporary directory

  if(Cores == 1){ # core check: if processing non-parallel
    if(Count_Krige == 1){ # count check: if this was the first actual computation
      T_End <- Sys.time() # record time at which kriging was done for current layer
      Duration <- as.numeric(T_End)-as.numeric(T_Begin) # calculate how long it took to krig on layer
      print(paste('Kriging of remaining ', nlayers(Data)-Iter_Krige, ' data layers should finish around: ', as.POSIXlt(T_Begin + Duration*nlayers(Data), tz = Sys.timezone(location=TRUE)), sep='')) # console output with estimate of when the kriging should be done
      ProgBar <- txtProgressBar(min = 0, max = nlayers(Data), style = 3) # create progress bar when non-parallel processing
      Count_Krige <- Count_Krige + 1 # raise count by one so the stimator isn't called again
    } # end of count check
    setTxtProgressBar(ProgBar, Iter_Krige) # update progress bar with number of current layer
  } # end of core check
  Sys.sleep(5400) # sleep for 1.5 hours to allow for RAM to cool down
  "
  
  ## KRIGING PREPARATION (establishing objects which the kriging refers to) ----
  Ras_Krig <- as.list(rep(NA, nlayers(Data))) # establish an empty list which will be filled with kriged layers
  Ras_Var <- as.list(rep(NA, nlayers(Data))) # establish an empty list which will be filled with kriged layers
  
  print("Commencing Kriging")
  ## DATA SKIPS (if certain layers in the data are empty and need to be skipped, this is handled here) ---
  if(!is.null(DataSkips)){ # Skip check: if layers need to be skipped
    for(Iter_Skip in DataSkips){ # Skip loop: loop over all layers that need to be skipped
      Ras_Krig[[Iter_Skip]] <- Data[[Iter_Skip]] # add raw data (which should be empty) to list
      writeRaster(x = Ras_Krig[[Iter_Skip]], filename = file.path(Dir.Temp, str_pad(Iter_Skip,4,'left','0')), overwrite = TRUE, format = 'CDF') # save raw layer to temporary directory, needed for loading back in when parallel processing
    } # end of Skip loop
    Layers_vec <- 1:nlayers(Data) # identify vector of all layers in data
    Compute_Layers <- Layers_vec[which(!Layers_vec %in% DataSkips)] # identify which layers can actually be computed on
  }else{ # if we don't need to skip any layers
    Compute_Layers <- 1:nlayers(Data) # set computing layers to all layers in data
  } # end of Skip check
  
  
  ## ACTUAL KRIGING (carry out kriging according to user specifications either in parallel or on a single core) ----
  if(Cores > 1){ # Cores check: if parallel processing has been specified
    ### PARALLEL KRIGING ---
    ForeachObjects <- c("Dir.Temp", "Cores", "Data", "KrigingEquation", "Origin", "Target", "Covariates_coarse", "Covariates_fine", "Terms", "SingularTry", "nmax") # objects which are needed for each kriging run and are thus handed to each cluster unit
    cl <- makeCluster(Cores) # Assuming Cores node cluster
    registerDoParallel(cl) # registering cores
    foreach(Iter_Krige = Compute_Layers, # kriging loop over all layers in Data, with condition (%:% when(...)) to only run if current layer is not present in Dir.Temp yet
            .packages = c("raster", "stringr", "automap", "ncdf4", "rgdal"), # import packages necessary to each itteration
            .export = ForeachObjects) %:% when(!paste0(str_pad(Iter_Krige,4,"left","0"), '_data.nc') %in% list.files(Dir.Temp)) %dopar% { # parallel kriging loop
              # print("Done")
              Ras_Krig <- eval(parse(text=looptext)) # evaluate the kriging specification per cluster unit per layer
            } # end of parallel kriging loop
    stopCluster(cl) # close down cluster
    Files_krig <- list.files(Dir.Temp)[grep(pattern = "_data.nc", x = list.files(Dir.Temp))]
    Files_var <- list.files(Dir.Temp)[grep(pattern = "_SE.nc", x = list.files(Dir.Temp))]
    for(Iter_Load in 1:length(Files_krig)){ # load loop: load data from temporary files in Dir.Temp
      Ras_Krig[[Iter_Load]] <- raster(file.path(Dir.Temp, Files_krig[Iter_Load])) # load current temporary file and write contents to list of rasters
      Ras_Var[[Iter_Load]] <- raster(file.path(Dir.Temp, Files_var[Iter_Load])) # load current temporary file and write contents to list of rasters
    } # end of load loop
  }else{ # if non-parallel processing has been specified
    ### NON-PARALLEL KRIGING ---
    Count_Krige <- 1 # Establish count variable which is targeted in kriging specification text for producing an estimator
    for(Iter_Krige in Compute_Layers){ # non-parallel kriging loop over all layers in Data
      if(paste0(str_pad(Iter_Krige,4,'left','0'), '_data.nc') %in% list.files(Dir.Temp)){ # file check: if this file has already been produced
        Ras_Krig[[Iter_Krige]] <- raster(file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_data.nc'))) # load already produced kriged file and save it to list of rasters
        Ras_Var[[Iter_Krige]] <- raster(file.path(Dir.Temp, paste0(str_pad(Iter_Krige,4,'left','0'), '_SE.nc')))
        if(!exists("ProgBar")){ProgBar <- txtProgressBar(min = 0, max = nlayers(Data), style = 3)} # create progress bar when non-parallel processing}
        setTxtProgressBar(ProgBar, Iter_Krige) # update progress bar
        next() # jump to next layer
      } # end of file check
      T_Begin <- Sys.time() # record system time when layer kriging starts
      eval(parse(text=looptext)) # evaluate the kriging specification per layer
    } # end of non-parallel kriging loop
  } # end of Cores check
  
  ## SAVING FINAL PRODUCT ----
  if(is.null(DataSkips)){ # Skip check: if no layers needed to be skipped
    Ras_Krig <- brick(Ras_Krig) # convert list of kriged layers in actual rasterbrick of kriged layers
    writeRaster(x = Ras_Krig, filename = file.path(Dir, FileName), overwrite = TRUE, format="CDF") # save final product as raster
    Ras_Var <- brick(Ras_Var) # convert list of kriged layers in actual rasterbrick of kriged layers
    writeRaster(x = Ras_Var, filename = file.path(Dir, paste0("SE_",FileName)), overwrite = TRUE, format="CDF") # save final product as raster
  }else{ # if some layers needed to be skipped
    warning(paste0("Some of the layers in your raster could not be kriged. You will find all the individual layers (kriged and not kriged) in ", Dir, "."))
    Keep_Temporary <- TRUE # keep temporary files so kriged products are not deleted
  } # end of Skip check
  
  ### REMOVE FILES FROM HARD DRIVE ---
  if(Keep_Temporary == FALSE){ # cleanup check
    unlink(Dir.Temp, recursive = TRUE)
  }  # end of cleanup check
  
  Krig_ls <- list(Ras_Krig, Ras_Var, Call_ls)
  names(Krig_ls) <- c("Kriging_Output", "Kriging_SE", "Call")
  return(Krig_ls) # return raster or list of layers
}