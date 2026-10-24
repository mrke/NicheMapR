#' ECMWF SEAS5 seasonal forecast implementation of the microclimate model.
#'
#' Runs the NicheMapR microclimate model from one ensemble member of an ECMWF SEAS5
#' seasonal forecast (Copernicus Climate Data Store dataset 'seasonal-original-single-levels',
#' 1 degree grid, 215 days from the 1st of the start month). Daily forcing is built from:
#' maximum/minimum 2 m air temperature over the last 24 h; 6-hourly 2 m dewpoint (daily mean
#' vapour pressure, giving relative humidity at the maximum and minimum temperatures); 6-hourly 10 m
#' wind (daily min/max speed, adjusted to 'Refhyt'); 6-hourly total cloud cover; accumulated
#' surface solar radiation and total precipitation (daily totals); and orography (grid-cell
#' elevation, for the lapse-rate correction).
#'
#' Data are downloaded with package 'ecmwfr' into 'seas5_path' and reused by later calls, so
#' all members of one forecast come from the same files. The CDS API key is taken from the
#' 'ecmwfr' keyring (under 'cds_user'), the 'CDSAPI_KEY' environment variable or ~/.cdsapirc, and the dataset's
#' Terms of Use must be accepted once at
#' https://cds.climate.copernicus.eu/datasets/seasonal-original-single-levels
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param dstart First day to run, date in format "d/m/Y" e.g. "01/08/2026"
#' @param dfinish Last day to run, date in format "d/m/Y" e.g. "31/12/2026"
#' @param init Forecast start, date in format "d/m/Y" (the 1st of a month); if NA, the 1st of the month of 'dstart'
#' @param member Ensemble member, 0 (control) to 50 (0 to 24 for 1981-2016 hindcasts)
#' @param seas5_path Folder for the downloaded SEAS5 NetCDF files
#' @param cds_time_out Seconds to wait for each CDS request (the CDS queue can take hours)
#' @param cds_user The 'user' label the CDS API key was stored under with \code{ecmwfr::wf_set_key}
#' @param cloud_from_solar Daily mean cloud from SEAS5 solar radiation relative to clear sky (1), keeping the 6-hourly cloud cover range, or 6-hourly cloud cover min/max directly (0)
#' @param dem A digital elevation model produced by microclima function 'get_dem' via R package 'elevatr' (internally generated via same function based on 'loc' if NA)
#' @param dem.res Requested resolution of the DEM from elevatr, m
#' @param pixels Number of pixels along one edge of square requested of DEM requested from elevatr, #
#' @param REFL Soil solar reflectance, decimal \%
#' @param elev Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param minshade Minimum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param maxshade Maximum shade level to us (\%) (can be a single value or a vector of daily values)
#' @param Usrhyt Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
#' @param ... Additional arguments, as for \code{\link{micro_access_s2}}
#' @return As for \code{\link{micro_access_s2}}, plus \code{init}, \code{member} and the daily forcing (\code{TMAXX}, \code{TMINN}, \code{RHMAXX}, \code{RHMINN}, \code{WNMAXX}, \code{WNMINN}, \code{CCMAXX}, \code{CCMINN}, \code{SOLAR} (MJ/m2/d), \code{PRESS} (mean sea level pressure, Pa))
#' @export
#' @examples
#' library(NicheMapR)
#' # forecast started 1 August 2026, ensemble members 0 to 4
#' micros <- lapply(0:4, function(member) micro_seas5(loc = c(144.28, -36.76),
#'   dstart = "01/08/2026", dfinish = "31/12/2026", member = member, seas5_path = "SEAS5/"))
#' soil <- lapply(micros, function(m) m$soil[, "D5cm"])
#' dates <- micros[[1]]$dates
#' plot(dates, soil[[1]], type = 'l', ylim = range(unlist(soil)), ylab = "soil temperature at 5 cm (°C)")
#' for(i in 2:length(soil)) points(dates, soil[[i]], type = 'l', col = i)
micro_seas5 <- function(
  loc = c(144.28, -36.76),
  dstart = "01/08/2026",
  dfinish = "31/12/2026",
  init = NA,
  member = 0,
  seas5_path = "SEAS5/",
  cds_time_out = 6 * 3600,
  cds_user = "ecmwfr",
  cloud_from_solar = 1,
  dem = NA,
  dem.res = 30,
  pixels = 100,
  nyears = as.numeric(substr(dfinish, 7, 10)) - as.numeric(substr(dstart, 7, 10)) + 1,
  REFL = 0.15,
  elev = NA,
  slope = 0,
  aspect = 0,
  lapse_max = 0.0077,
  lapse_min = 0.0039,
  DEP = c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
  minshade = 0,
  maxshade = 90,
  Usrhyt = 0.01,
  Z01 = 0,
  Z02 = 0,
  ZH1 = 0,
  ZH2 = 0,
  runshade = 1,
  clearsky = 0,
  solonly = 0,
  run.gads = 1,
  Soil_Init = NA,
  write_input = 0,
  writecsv = 0,
  windfac = 1,
  adiab_cor = 1,
  warm = 0,
  ERR = 1,
  RUF = 0.004,
  ZH = 0,
  D0 = 0,
  EC = 0.0167238,
  SLE = 0.95,
  Thcond = 2.5,
  Density = 2.56,
  SpecHeat = 870,
  BulkDensity = 1.3,
  PCTWET = 0,
  rainwet = 1.5,
  cap = 1,
  CMH2O = 1,
  hori = rep(0, 24),
  TIMAXS = c(1.0, 1.0, 0.0, 0.0),
  TIMINS = c(0, 0, 1, 1),
  runmoist = 1,
  PE = rep(1.1, 19),
  KS = rep(0.0037, 19),
  BB = rep(4.5, 19),
  BD = rep(BulkDensity, 19),
  DD = rep(Density, 19),
  maxpool = 10000,
  rainmult = 1,
  evenrain = 0,
  SoilMoist_Init = c(0.1, 0.12, 0.15, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
  L = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000,
  R1 = 0.001,
  RW = 2.5e+10,
  RL = 2e+6,
  PC = -1500,
  SP = 10,
  IM = 1e-06,
  MAXCOUNT = 500,
  LAI = 0.1,
  snowmodel = 1,
  snowtemp = 1.5,
  snowdens = 0.375,
  densfun = c(0.5979, 0.2178, 0.001, 0.0038),
  snowmelt = 1,
  undercatch = 1,
  rainmelt = 0.0125,
  shore = 0,
  tides = 0,
  hourly = 0,
  rainhourly = 0,
  rainhour = 0,
  rainoff = 0,
  lamb = 0,
  IUV = 0,
  ndmax = 3,
  soilgrids = 0,
  IR = 0,
  message = 0,
  fail = nyears * 24 * 365,
  runmicro = 1,
  snowcond = 0,
  intercept = max(maxshade) / 100 * 0.3,
  grasshade = 0,
  maxsurf = 85) {

  Refhyt <- 2 # reference height (m) of the forcing
  days <- seq(as.Date(dstart, format = "%d/%m/%Y"), as.Date(dfinish, format = "%d/%m/%Y"), by = "days")
  if(is.na(init)){
    init_date <- as.Date(format(days[1], "%Y-%m-01"))
  }else{
    init_date <- as.Date(init, format = "%d/%m/%Y")
  }
  n_seas5_days <- 215
  last_day <- init_date + n_seas5_days - 1

  # error trapping
  errors <- 0
  if(format(init_date, "%d") != "01"){
    cat("ERROR: SEAS5 forecasts start on the 1st of a month; 'init' must be the 1st.", '\n')
    errors <- 1
  }
  if(days[1] < init_date | days[length(days)] > last_day){
    cat(paste0("ERROR: the SEAS5 forecast started ", init_date, " covers ", init_date, " to ", last_day,
      "; 'dstart' and 'dfinish' must lie within it.", '\n'))
    errors <- 1
  }
  if(member < 0 | member > 50){
    cat("ERROR: 'member' must be between 0 and 50.", '\n')
    errors <- 1
  }
  if(cloud_from_solar %in% c(0, 1) == FALSE){
    cat("ERROR: 'cloud_from_solar' must be 0 or 1.", '\n')
    errors <- 1
  }
  if(DEP[2] - DEP[1] > 3 | DEP[3] - DEP[2] > 3){
    cat("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if(DEP[2] - DEP[1] < 2){
    cat("warning, nodes might be too close near the surface, try a different spacing if the program is crashing \n")
  }
  if(DEP[10] != 200){
    cat("warning, last depth in soil should not be changed from 200 without good reason \n")
  }
  if(loc[1] > 180 | loc[1] < -180 | loc[2] > 90 | loc[2] < -90){
    cat("ERROR: Latitude or longitude (loc) is out of bounds.", '\n')
    errors <- 1
  }
  if(run.gads %in% c(0, 1, 2) == FALSE){
    cat("ERROR: the variable 'run.gads' be either 0, 1 or 2.", '\n')
    errors <- 1
  }
  if(EC < 0.0034 | EC > 0.058){
    cat("ERROR: the eccentricity variable (EC) is out of bounds (0.0034 - 0.058).", '\n')
    errors <- 1
  }
  if(RUF < 0.0001 | RUF > 2){
    cat("ERROR: The roughness height (RUF) must be between 0.0001 and 2 m.", '\n')
    errors <- 1
  }
  if(D0 > 0 & D0 < Usrhyt){
    cat("ERROR: The zero plane displacement height (D0) must be lower than the local height (Usrhyt).", '\n')
    errors <- 1
  }
  if(length(DEP) != 10 | DEP[1] != 0 | any(diff(DEP) <= 0) | DEP[10] > 500){
    cat("ERROR: DEP must be 10 ascending depths from 0 to at most 500 cm.", '\n')
    errors <- 1
  }
  if(min(Thcond) < 0 | min(Density) < 0 | min(SpecHeat) < 0 | min(BulkDensity) < 0){
    cat("ERROR: Thcond, Density, SpecHeat and BulkDensity must be positive.", '\n')
    errors <- 1
  }
  if(REFL < 0 | REFL > 1){
    cat("ERROR: Soil reflectivity value (REFL) must be between 0 and 1.", '\n')
    errors <- 1
  }
  if(slope < 0 | slope > 90 | aspect < 0 | aspect > 365){
    cat("ERROR: slope must be 0 - 90 and aspect 0 - 365 degrees.", '\n')
    errors <- 1
  }
  if(length(hori) != 24 | max(hori) > 90 | min(hori) < 0){
    cat("ERROR: hori must be 24 horizon angles between 0 and 90 degrees.", '\n')
    errors <- 1
  }
  if(SLE < 0.05 | SLE > 1){
    cat("ERROR: Emissivity (SLE) is out of bounds (0.05 - 1.00).", '\n')
    errors <- 1
  }
  if(ERR < 0){
    cat("ERROR: Error bound (ERR) must be > 0.", '\n')
    errors <- 1
  }
  if(Usrhyt < RUF | Usrhyt > Refhyt){
    cat("ERROR: Local height (Usrhyt) must lie between the roughness height (RUF) and the reference height (2 m).", '\n')
    errors <- 1
  }
  if(CMH2O < 0.5 | CMH2O > 2){
    cat("ERROR: Preciptable water in air column (CMH2O) is out of bounds (0.5 - 2 cm).", '\n')
    errors <- 1
  }
  if(max(TIMAXS) > 24 | min(TIMAXS) < 0 | max(TIMINS) > 24 | min(TIMINS) < 0){
    cat("ERROR: TIMAXS and TIMINS must be between 0 and 24.", '\n')
    errors <- 1
  }
  if(max(minshade - maxshade) >= 0 | max(minshade) > 100 | min(minshade) < 0 | max(maxshade) > 100 | min(maxshade) < 0){
    cat("ERROR: shade must be 0 - 100 \\%, with minshade below maxshade.", '\n')
    errors <- 1
  }
  if(errors == 1){
    return(invisible(NULL))
  }

  for(p in c("ecmwfr", "RNetCDF", "terra")){
    if(!requireNamespace(p, quietly = TRUE)){
      stop(paste0("package '", p, "' is needed. Please install it."), call. = FALSE)
    }
  }

  ################## location and terrain #################################
  x <- t(as.matrix(as.numeric(c(loc[1], loc[2]))))
  longlat <- loc
  ALREF <- abs(trunc(x[1])) # local solar noon
  HEMIS <- ifelse(x[2] < 0, 2, 1)
  ALAT <- abs(trunc(x[2]))
  AMINUT <- (abs(x[2]) - ALAT) * 60
  ALONG <- abs(trunc(x[1]))
  ALMINT <- (abs(x[1]) - ALONG) * 60
  azmuth <- aspect

  if(is.na(elev)){
    if(class(dem)[1] == "SpatRaster"){
      cat('using DEM provided to function call \n')
    }else{
      if(!requireNamespace("microclima", quietly = TRUE)){
        stop("package 'microclima' is needed to get the elevation (or supply 'elev'). Please install it.", call. = FALSE)
      }
      cat('downloading DEM via package elevatr \n')
      dem <- microclima::get_dem(lat = loc[2], long = loc[1], resolution = dem.res, xdims = pixels, ydims = pixels)
    }
    xy <- sf::st_set_crs(sf::st_as_sf(data.frame(lon = loc[1], lat = loc[2]), coords = c("lon", "lat")), "EPSG:4326")
    xy <- sf::st_transform(xy, sf::st_crs(dem))
    elev <- as.numeric(terra::extract(dem, terra::vect(xy))[, 2])
  }
  ALTITUDES <- elev

  if(soilgrids == 1){
    cat('extracting soil texture data from SoilGrids \n')
    ov <- jsonlite::fromJSON(paste0('https://rest.isric.org/soilgrids/v2.0/properties/query?lon=', x[1], '&lat=', x[2], '&property=bdod&property=silt&property=clay&property=sand'), flatten = TRUE)
    if(length(ov) > 3){
      soilpro <- cbind(c(0, 5, 15, 30, 60, 100), unlist(ov$properties$layers$depths[[1]]$values.mean) / 100, unlist(ov$properties$layers$depths[[2]]$values.mean) / 10, unlist(ov$properties$layers$depths[[4]]$values.mean) / 10, unlist(ov$properties$layers$depths[[3]]$values.mean) / 10)
      soilpro <- rbind(soilpro, soilpro[6, ])
      soilpro[7, 1] <- 200
      colnames(soilpro) <- c('depth', 'blkdens', 'clay', 'silt', 'sand')
      soil.hydro <- pedotransfer(soilpro = as.data.frame(soilpro), DEP = DEP)
      PE <- soil.hydro$PE
      BB <- soil.hydro$BB
      BD <- soil.hydro$BD
      KS <- soil.hydro$KS
      BulkDensity <- BD[seq(1, 19, 2)]
    }else{
      cat('no SoilGrids data for this site, using user-input soil properties \n')
    }
  }

  VIEWF <- 1 - sum(sin(hori * pi / 180)) / length(hori)

  ################## SEAS5 download #################################
  # whole-degree box: CDS returns the 1-degree cell centres inside it
  area <- c(ceiling(x[2]), floor(x[1]), floor(x[2]), ceiling(x[1])) # N, W, S, E
  if(area[1] == area[3]) area[1] <- area[1] + 1
  if(area[2] == area[4]) area[4] <- area[4] + 1
  daily_leads <- seq(24, 24 * n_seas5_days, 24)
  sixhourly_leads <- seq(6, 24 * n_seas5_days, 6)
  fields <- list(
    tmax = list(cds = "maximum_2m_temperature_in_the_last_24_hours", leads = daily_leads),
    tmin = list(cds = "minimum_2m_temperature_in_the_last_24_hours", leads = daily_leads),
    d2m = list(cds = "2m_dewpoint_temperature", leads = sixhourly_leads),
    u10 = list(cds = "10m_u_component_of_wind", leads = sixhourly_leads),
    v10 = list(cds = "10m_v_component_of_wind", leads = sixhourly_leads),
    tcc = list(cds = "total_cloud_cover", leads = sixhourly_leads),
    msl = list(cds = "mean_sea_level_pressure", leads = sixhourly_leads),
    ssrd = list(cds = "surface_solar_radiation_downwards", leads = daily_leads),
    tp = list(cds = "total_precipitation", leads = daily_leads),
    orog = list(cds = "orography", leads = 0)
  )
  dir.create(seas5_path, showWarnings = FALSE, recursive = TRUE)
  files <- sapply(names(fields), function(f){
    paste0("seas5_", format(init_date, "%Y%m"), "_", f, "_", paste(area, collapse = "_"), ".nc")
  })
  missing <- names(fields)[!file.exists(file.path(seas5_path, files))]
  if(length(missing) > 0){
    if(nchar(Sys.getenv("ecmwfr_PAT")) == 0 & inherits(try(ecmwfr::wf_get_key(user = cds_user), silent = TRUE), "try-error")){
      key <- Sys.getenv("CDSAPI_KEY")
      # R's "~" is Documents on Windows; .cdsapirc usually sits in the user profile
      rc <- c(file.path(Sys.getenv("USERPROFILE"), ".cdsapirc"), path.expand("~/.cdsapirc"))
      rc <- rc[file.exists(rc)]
      if(nchar(key) == 0 & length(rc) > 0){
        line <- grep("^key:", readLines(rc[1], warn = FALSE), value = TRUE)
        if(length(line) > 0) key <- trimws(sub("^key:", "", line[1]))
      }
      if(nchar(key) == 0){
        stop("no CDS API key: set one with ecmwfr::wf_set_key(), ENV CDSAPI_KEY or ~/.cdsapirc", call. = FALSE)
      }
      Sys.setenv(ecmwfr_PAT = key)
    }
    requests <- lapply(missing, function(f){
      list(dataset_short_name = "seasonal-original-single-levels",
        originating_centre = "ecmwf", system = "51",
        variable = fields[[f]]$cds,
        year = format(init_date, "%Y"), month = format(init_date, "%m"), day = "01",
        leadtime_hour = as.character(fields[[f]]$leads),
        area = area, data_format = "netcdf", target = files[[f]])
    })
    cat(paste0("requesting ", length(requests), " SEAS5 field(s) for the forecast started ", init_date,
      " from the CDS (queued requests can take minutes to hours) \n"))
    ecmwfr::wf_request_batch(requests, workers = 4, user = cds_user, path = seas5_path, time_out = cds_time_out)
  }

  # one field for 'member' at the grid cell nearest 'loc', as a vector over lead times
  read_seas5 <- function(f){
    nc <- RNetCDF::open.nc(file.path(seas5_path, files[[f]]))
    on.exit(RNetCDF::close.nc(nc))
    nvars <- RNetCDF::file.inq.nc(nc)$nvars
    dimname <- function(id) RNetCDF::dim.inq.nc(nc, id)$name
    for(v in 0:(nvars - 1)){
      inq <- RNetCDF::var.inq.nc(nc, v)
      if(all(is.na(inq$dimids))) next # scalar, e.g. valid_time with a single lead
      dims <- sapply(inq$dimids, dimname)
      if(any(dims %in% c("longitude", "lon")) & any(dims %in% c("latitude", "lat"))) break
    }
    lon <- RNetCDF::var.get.nc(nc, dims[dims %in% c("longitude", "lon")])
    lat <- RNetCDF::var.get.nc(nc, dims[dims %in% c("latitude", "lat")])
    data <- RNetCDF::var.get.nc(nc, inq$name, unpack = TRUE, collapse = FALSE)
    units <- tryCatch(RNetCDF::att.get.nc(nc, inq$name, "units"), error = function(e) "")
    index <- lapply(dims, function(d){
      if(d %in% c("longitude", "lon")) return(which.min(abs(lon - x[1])))
      if(d %in% c("latitude", "lat")) return(which.min(abs(lat - x[2])))
      if(d %in% c("number", "realization", "member")){
        m <- which(RNetCDF::var.get.nc(nc, d) == member)
        if(length(m) == 0) stop(paste0("SEAS5 ensemble member ", member, " not in ", files[[f]], " (hindcasts before 2017 have members 0-24)"), call. = FALSE)
        return(m)
      }
      seq_len(dim(data)[match(d, dims)])
    })
    values <- as.numeric(do.call(`[`, c(list(data), index)))
    if(length(values) != length(fields[[f]]$leads)){
      stop(paste0("SEAS5 file ", files[[f]], " has ", length(values), " values, expected ", length(fields[[f]]$leads), " lead times"), call. = FALSE)
    }
    attr(values, "units") <- units
    values
  }
  # 6-hourly values to daily (four per UTC day, ending at 24 h leads)
  daily <- function(v, fun) apply(matrix(v, nrow = 4), 2, fun)
  accumulated <- function(v) pmax(diff(c(0, v)), 0) # accumulated from the start; negative steps are packing noise

  Tmax_all <- read_seas5("tmax") - 273.15
  Tmin_all <- read_seas5("tmin") - 273.15
  e_all <- daily(WETAIR(db = read_seas5("d2m") - 273.15, rh = 100)$esat, mean) # vapour pressure, Pa
  wind_all <- sqrt(read_seas5("u10") ^ 2 + read_seas5("v10") ^ 2) * (Refhyt / 10) ^ 0.15 # 10 m to Refhyt
  tcc <- read_seas5("tcc")
  tcc_mean_all <- daily(tcc, mean)
  tcc_min_all <- daily(tcc, min)
  tcc_max_all <- daily(tcc, max)
  press_all <- daily(read_seas5("msl"), mean)
  solar_all <- accumulated(read_seas5("ssrd")) / 1e6 # MJ/m2/d
  rain_all <- accumulated(read_seas5("tp")) * 1000 # mm/d
  orog <- read_seas5("orog")
  if(grepl("s\\*\\*-2|s-2", attr(orog, "units"))) orog <- orog / 9.80665 # geopotential to m

  keep <- as.numeric(days - init_date) + 1
  Tmax <- Tmax_all[keep]
  Tmin <- Tmin_all[keep]
  e <- e_all[keep]
  WNMAXX <- daily(wind_all, max)[keep] * windfac
  WNMINN <- daily(wind_all, min)[keep] * windfac
  SOLAR <- solar_all[keep]
  PRESS <- press_all[keep]
  Rain <- rain_all[keep]
  ndays <- length(days)
  doy <- as.numeric(format(days, "%j"))

  if(is.na(ALTITUDES)){
    ALTITUDES <- as.numeric(orog)
  }
  delta_elev <- as.numeric(orog) - ALTITUDES
  if(adiab_cor == 1){
    TMAXX <- as.matrix(Tmax + delta_elev * lapse_max)
    TMINN <- as.matrix(Tmin + delta_elev * lapse_min)
  }else{
    TMAXX <- as.matrix(Tmax)
    TMINN <- as.matrix(Tmin)
  }
  TMAXX <- TMAXX + warm
  TMINN <- TMINN + warm

  # relative humidity from the daily vapour pressure at the (elevation-corrected) daily extremes
  RHMINN <- pmin(pmax(e / WETAIR(db = TMAXX, rh = 100)$esat * 100, 0.01), 100)
  RHMAXX <- pmin(pmax(e / WETAIR(db = TMINN, rh = 100)$esat * 100, 0.01), 100)

  # cloud (%), from solar relative to clear sky or from cloud cover directly
  if(cloud_from_solar == 1){
    cat("running micro_global to get clear sky solar \n")
    micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, timeinterval = 365, solonly = 1)
    clearsky_daily <- aggregate(micro_clearsky$metout[, 13] / 1e6 * 3600, by = list(micro_clearsky$metout[, 1]), FUN = sum)[, 2]
    cloud <- pmin(pmax(1 - SOLAR / clearsky_daily[pmin(doy, 365)], 0), 1)
    CCMINN <- pmax(cloud - (tcc_mean_all[keep] - tcc_min_all[keep]), 0) * 100
    CCMAXX <- pmin(cloud + (tcc_max_all[keep] - tcc_mean_all[keep]), 1) * 100
  }else{
    CCMINN <- tcc_min_all[keep] * 100
    CCMAXX <- tcc_max_all[keep] * 100
  }
  if(clearsky == 1){
    CCMINN <- CCMINN * 0
    CCMAXX <- CCMAXX * 0
  }

  ################## model setup #################################
  microdaily <- 1
  idayst <- 1
  ida <- ndays
  MINSHADES <- if(length(minshade) == ndays) minshade else rep(minshade[1], ndays)
  MAXSHADES <- if(length(maxshade) == ndays) maxshade else rep(maxshade[1], ndays)

  if(run.gads > 0){
    relhum <- 1
    if(run.gads == 1){
      optdep.summer <- as.data.frame(rungads(longlat[2], longlat[1], relhum, 0))
      optdep.winter <- as.data.frame(rungads(longlat[2], longlat[1], relhum, 1))
    }else{
      optdep.summer <- as.data.frame(gads.r(longlat[2], longlat[1], relhum, 0))
      optdep.winter <- as.data.frame(gads.r(longlat[2], longlat[1], relhum, 1))
    }
    optdep <- as.data.frame(cbind(optdep.winter[, 1], rowMeans(cbind(optdep.summer[, 2], optdep.winter[, 2]))))
    colnames(optdep) <- c("LAMBDA", "OPTDEPTH")
    a <- lm(OPTDEPTH ~ poly(LAMBDA, 6, raw = TRUE), data = optdep)
    LAMBDA <- c(290,295,300,305,310,315,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1380,1400,1420,1440,1460,1480,1500,1540,1580,1600,1620,1640,1660,1700,1720,1780,1800,1860,1900,1950,2000,2020,2050,2100,2120,2150,2200,2260,2300,2320,2350,2380,2400,2420,2450,2490,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000)
    TAI <- predict(a, data.frame(LAMBDA))
  }else{ # as for southern Australia
    TAI <- c(0.0670358341290886,0.0662612704779235,0.065497075238002,0.0647431301168489,0.0639993178022531,0.0632655219571553,0.0625416272145492,0.0611230843885423,0.0597427855962549,0.0583998423063099,0.0570933810229656,0.0558225431259535,0.0545864847111214,0.0533843764318805,0.0522154033414562,0.0499736739981675,0.047855059159556,0.0458535417401334,0.0439633201842001,0.0421788036108921,0.0404946070106968,0.0389055464934382,0.0374066345877315,0.0359930755919066,0.0346602609764008,0.0334037648376212,0.0322193394032758,0.0311029105891739,0.0300505736074963,0.0290585886265337,0.0281233764818952,0.0272415144391857,0.0264097320081524,0.0256249068083005,0.0248840604859789,0.0241843546829336,0.0235230870563317,0.0228976873502544,0.0223057135186581,0.0217448478998064,0.0212128934421699,0.0207077699817964,0.0202275105711489,0.0197702578594144,0.0193342605242809,0.0189178697551836,0.0177713140039894,0.0174187914242432,0.0170790495503944,0.0167509836728154,0.0164335684174899,0.0161258546410128,0.0158269663770596,0.0155360978343254,0.0152525104459325,0.0149755299703076,0.0147045436435285,0.0144389973831391,0.0141783930434343,0.0134220329447663,0.0131772403830191,0.0129356456025128,0.0126970313213065,0.0124612184223418,0.0122280636204822,0.01199745718102,0.0115436048739351,0.0110993711778668,0.0108808815754663,0.0106648652077878,0.0104513876347606,0.0102405315676965,0.00982708969547694,0.00962473896278535,0.00903679230300494,0.00884767454432418,0.0083031278398166,0.00796072474935954,0.00755817587626185,0.00718610751850881,0.00704629977586921,0.00684663903049612,0.00654155580333479,0.00642947339729728,0.00627223096874308,0.00603955966866779,0.00580920937536261,0.00568506186880564,0.00563167068287251,0.00556222005081865,0.00550522989971023,0.00547395763028062,0.0054478983436216,0.00541823364504573,0.00539532163908382,0.00539239864119488,0.00541690124712384,0.00551525885358836,0.00564825853509463,0.00577220185074264,0.00584222986640171,0.00581645238345584,0.00566088137411449,0.00535516862329704,0.00489914757707667,0.00432017939770409,0.0036813032251836,0.00309019064543606,0.00270890436501562,0.00276446109239711,0.00356019862584603)
  }

  RAINFALL <- pmax(Rain + rainoff, 0)
  REFLS <- rep(REFL, ndays)
  soilwet <- ifelse(RAINFALL > rainwet, 90, 0)
  PCTWET <- pmax(soilwet, rep(PCTWET, ndays))

  Numtyps <- 10
  Nodes <- matrix(data = 0, nrow = 10, ncol = ndays)
  Nodes[1:10, ] <- c(1:10)
  ALTT <- as.numeric(ALTITUDES)

  avetemp <- (sum(TMAXX) + sum(TMINN)) / (length(TMAXX) * 2)
  if(is.na(Soil_Init[1])){
    soilinit <- rep(avetemp, 20)
    spinup <- 1
  }else{
    if(snowmodel == 0){
      soilinit <- c(Soil_Init, rep(avetemp, 10))
    }else{
      soilinit <- c(rep(avetemp, 8), Soil_Init[1:10], rep(avetemp, 2))
    }
    spinup <- 0
  }
  tannul <- avetemp
  tannulrun <- rep(avetemp, ndays)

  SLES <- matrix(nrow = ndays, data = SLE)
  moists <- matrix(nrow = 10, ncol = ndays, data = 0)
  if(runmoist == 1){
    moists[1:10, ] <- SoilMoist_Init
  }else{
    moists[1, ndays] <- 0.2
  }
  soilprops <- matrix(data = 0, nrow = 10, ncol = 5)
  soilprops[, 1] <- BulkDensity
  soilprops[, 2] <- pmax(1 - BulkDensity / Density, 0.26) # not used if soil moisture computed
  soilprops[, 3] <- Thcond
  soilprops[, 4] <- SpecHeat
  soilprops[, 5] <- Density
  if(cap == 1){
    soilprops[1:2, 3] <- 0.2
    soilprops[1:2, 4] <- 1920
  }
  if(cap == 2){
    soilprops[1:2, 3] <- 0.1
    soilprops[3:4, 3] <- 0.25
    soilprops[1:4, 4] <- 1920
    soilprops[1:4, 5] <- 1.3
    soilprops[1:4, 1] <- 0.7
  }

  microinput <- c(ndays, RUF, ERR, Usrhyt, Refhyt, Numtyps, Z01, Z02, ZH1, ZH2, idayst, ida, HEMIS, ALAT, AMINUT, ALONG, ALMINT, ALREF, slope, azmuth, ALTT, CMH2O, microdaily, tannul, EC, VIEWF, snowtemp, snowdens, snowmelt, undercatch, rainmult, runshade, runmoist, maxpool, evenrain, snowmodel, rainmelt, writecsv, densfun, hourly, rainhourly, lamb, IUV, RW, PC, RL, SP, R1, IM, MAXCOUNT, IR, message, fail, snowcond, intercept, grasshade, solonly, ZH, D0, TIMAXS, TIMINS, spinup, 0, 360, maxsurf, ndmax)

  TAIRhr <- rep(0, 24 * ndays)
  RHhr <- rep(0, 24 * ndays)
  WNhr <- rep(0, 24 * ndays)
  CLDhr <- rep(0, 24 * ndays)
  SOLRhr <- rep(0, 24 * ndays)
  ZENhr <- rep(-1, 24 * ndays)
  IRDhr <- rep(-1, 24 * ndays)
  RAINhr <- if(rainhourly == 0) rep(0, 24 * ndays) else rainhour
  if(length(LAI) < ndays){
    LAI <- rep(LAI[1], ndays)
  }
  if(shore == 0){
    tides <- matrix(data = 0, nrow = 24 * ndays, ncol = 3)
  }
  hori <- as.numeric(hori)

  micro <- list(tides = tides, microinput = microinput, doy = doy, SLES = SLES, DEP = DEP, Nodes = Nodes, MAXSHADES = MAXSHADES, MINSHADES = MINSHADES, TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN, CCMAXX = CCMAXX, CCMINN = CCMINN, WNMAXX = WNMAXX, WNMINN = WNMINN, TAIRhr = TAIRhr, RHhr = RHhr, WNhr = WNhr, CLDhr = CLDhr, SOLRhr = SOLRhr, RAINhr = RAINhr, ZENhr = ZENhr, IRDhr = IRDhr, REFLS = REFLS, PCTWET = PCTWET, soilinit = soilinit, hori = hori, TAI = TAI, soilprops = soilprops, moists = moists, RAINFALL = RAINFALL, tannulrun = tannulrun, PE = PE, KS = KS, BB = BB, BD = BD, DD = DD, L = L, LAI = LAI)
  if(write_input == 1){
    dir.create("micro csv input", showWarnings = FALSE)
    for(n in setdiff(names(micro), "microinput")){
      write.table(micro[[n]], file = paste0("micro csv input/", n, ".csv"), sep = ",", col.names = NA, qmethod = "double")
    }
    write.table(as.matrix(microinput), file = "micro csv input/microinput.csv", sep = ",", col.names = NA, qmethod = "double")
  }

  # local solar time, labelled UTC so there are no daylight saving jumps
  dates <- as.POSIXct(paste(rep(days, each = 24), sprintf("%02d:00", 0:23)), tz = "UTC")
  dates2 <- as.POSIXct(days, tz = "UTC")
  forcing <- list(RAINFALL = RAINFALL, TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN, WNMAXX = WNMAXX, WNMINN = WNMINN, CCMAXX = CCMAXX, CCMINN = CCMINN, SOLAR = SOLAR, PRESS = PRESS, init = init_date, member = member, dates = dates, dates2 = dates2, elev = ALTT, PE = PE, BD = BD, DD = DD, BB = BB, KS = KS)
  if(runmicro != 1){
    return(forcing)
  }

  cat(paste('running microclimate model for', ndays, 'days from', dstart, 'to', dfinish, 'at site long', loc[1], 'lat', loc[2], 'with SEAS5 member', member, '\n'))
  message('Note: the output column `SOLR` in metout and shadmet is for unshaded horizontal plane solar radiation \n')
  ptm <- proc.time()
  microut <- microclimate(micro)
  print(proc.time() - ptm)
  if(max(microut$metout[, 1] == 0)){
    cat("ERROR: the model crashed - try a different error tolerance (ERR) or a different spacing in DEP", '\n')
  }

  out <- microut[c("metout", "shadmet", "soil", "shadsoil", "tcond", "shadtcond", "specheat", "shadspecheat", "densit", "shaddensit")]
  if(runmoist == 1){
    out <- c(out, microut[c("soilmoist", "shadmoist", "humid", "shadhumid", "soilpot", "shadpot", "plant", "shadplant")])
  }else{
    soil <- microut$soil
    soilpot <- soil; soilpot[, 3:12] <- 0
    soilmoist <- soil; soilmoist[, 3:12] <- 0.5
    humid <- soil; humid[, 3:12] <- 0.99
    plant <- cbind(soil, soil[, 3:4]); plant[, 3:14] <- 0
    out <- c(out, list(soilmoist = soilmoist, shadmoist = soilmoist, humid = humid, shadhumid = humid, soilpot = soilpot, shadpot = soilpot, plant = plant, shadplant = plant))
  }
  if(snowmodel == 1){
    out <- c(out, microut[c("sunsnow", "shdsnow")])
  }
  if(lamb == 1){
    out <- c(out, microut[c("drlam", "drrlam", "srlam")])
  }
  c(out, forcing, list(ndays = ndays, REFL = REFL, longlat = c(x[1], x[2]), nyears = nyears, minshade = MINSHADES, maxshade = MAXSHADES, DEP = DEP, diffuse_frac = NA, dem = dem))
}
