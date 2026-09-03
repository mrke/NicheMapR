#' BARRA implementation of the microclimate model
#'
#' An implementation of the NicheMapR microclimate model driven by BARRA2
#' (Bureau of Meteorology Atmospheric high-resolution Regional Reanalysis for
#' Australia, v2), Australia's ~4 km regional analog to ERA5. Data are read
#' directly from NCI's public THREDDS server via OPeNDAP -- no download or
#' local data folder needed, unlike micro_ncep/micro_era5.
#'
#' Unlike micro_era5 (which relies on package mcera5's ERA5-specific direct/
#' diffuse radiation and net-longwave fields), BARRA only publishes total
#' downwelling shortwave (rsds) and downwelling longwave (rlds) -- the same
#' near-surface subset MicroclimateMapper.jl's BARRA support uses (see
#' src/climate/barra.jl: tas, sfcWind, hurs, psl, rsds, rlds, pr, orog).
#' Cloud cover is therefore derived from rsds via the same clear-sky-ratio
#' method the daily-forcing functions (micro_aust, micro_ncep) use, computed
#' hourly instead of daily, and the Fortran solver's own hourly forcing mode
#' (hourly = 1) is used so BARRA's real sub-daily variation drives the model
#' directly rather than being collapsed to a daily min/max and resynthesised.
#'
#' No terrain/DEM downscaling is attempted by default (flat terrain, BARRA's
#' own `orog` grid supplies point elevation) -- pass `slope`/`aspect`/`hori`
#' directly if you want terrain effects, same as the other micro_*.R
#' functions.
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param dstart First day to run, date in format "d/m/Y" e.g. "01/01/2016"
#' @param dfinish Last day to run, date in format "d/m/Y" e.g. "31/12/2016"
#' @param barra_product "BARRAC2" (~4 km, convection-permitting, default) or "BARRAR2" (~11 km, coarser regional reanalysis)
#' @param barra_domain BARRA domain: "AUST04" (BARRAC2 only, default), "AUS11" or "AUST11" (BARRAR2 only)
#' @param elev Elevation (m). If NA (default), fetched from BARRA's own `orog` grid at `loc`.
#' @param REFL Soil solar reflectance, decimal \%
#' @param slope Slope in degrees (default 0, flat terrain)
#' @param aspect Aspect in degrees, 0 = north (default 0)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param minshade Minimum shade level to use (\%)
#' @param maxshade Maximum shade level to use (\%)
#' @param Usrhyt Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
#' @param ... Additional arguments, see Details
#' @return metout The above ground micrometeorological conditions under the minimum specified shade
#' @return shadmet The above ground micrometeorological conditions under the maximum specified shade
#' @return soil Hourly predictions of the soil temperatures under the minimum specified shade
#' @return shadsoil Hourly predictions of the soil temperatures under the maximum specified shade
#' @return soilmoist Hourly predictions of the soil moisture under the minimum specified shade
#' @return shadmoist Hourly predictions of the soil moisture under the maximum specified shade
#' @usage micro_barra(loc = c(133.8807, -23.6980), dstart = "01/01/2020", dfinish = "31/12/2020",
#' REFL = 0.15, DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minshade = 0, maxshade = 90,
#' Usrhyt = 0.01, ...)
#' @export
#' @details
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#' \code{runshade}{ = 1, populate shadmet/shadsoil/shadmoist etc. under maxshade (1) or just report minshade (0)?}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes (Fortran version), 2=yes (R version), 0=no (use an Australia-typical default)}\cr\cr
#' \code{solonly}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{ndmax}{ = 3, iterations of first day to get a steady periodic}\cr\cr
#' \code{Soil_Init}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{windfac}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#' \code{runmicro}{ = 1, call the microclimate model (1) or not (0), if you just want the downscaled input weather data}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{ERR}{ = 1, Integrator error tolerance for soil temperature calculations}\cr\cr
#' \code{RUF}{ = 0.004, Roughness height (m)}\cr\cr
#' \code{EC}{ = 0.0167238, Eccentricity of the earth's orbit}\cr\cr
#' \code{SLE}{ = 0.95, Substrate longwave IR emissivity (decimal \%)}\cr\cr
#' \code{Thcond}{ = 2.5, Soil minerals thermal conductivity, single value or vector of 10 specific to each depth (W/mK)}\cr\cr
#' \code{Density}{ = 2.56, Soil minerals density, single value or vector of 10 specific to each depth (Mg/m3)}\cr\cr
#' \code{SpecHeat}{ = 870, Soil minerals specific heat, single value or vector of 10 specific to each depth (J/kg-K)}\cr\cr
#' \code{BulkDensity}{ = 1.3, Soil bulk density (Mg/m3), single value or vector of 10 specific to each depth}\cr\cr
#' \code{rainwet}{ = 1.5, mm of rainfall causing the ground to be 90\% wet for the day}\cr\cr
#' \code{cap}{ = 1, organic cap present on soil surface? (lower conductivity, higher specific heat)}\cr\cr
#' \code{CMH2O}{ = 1, Precipitable cm H2O in air column}\cr\cr
#' \code{hori}{ = rep(0,24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#'
#' \strong{ Soil moisture mode parameters:}\cr\cr
#' \code{runmoist}{ = 1, Run soil moisture model? 1=yes, 0=no}\cr\cr
#' \code{PE}{ = rep(1.1,19), Air entry potential (J/kg), 19 values descending through soil}\cr\cr
#' \code{KS}{ = rep(0.0037,19), Saturated conductivity (kg s/m3), 19 values}\cr\cr
#' \code{BB}{ = rep(4.5,19), Campbell's soil 'b' parameter (-), 19 values}\cr\cr
#' \code{BD}{ = rep(1.3,19), Soil bulk density (Mg/m3), 19 values}\cr\cr
#' \code{DD}{ = rep(2.56,19), Soil density (Mg/m3), 19 values}\cr\cr
#' \code{L}{ = c(0,0,8.2,8.0,7.8,7.4,7.1,6.4,5.8,4.8,4.0,1.8,0.9,0.6,0.8,0.4,0.4,0,0)*10000, root density (m/m3), 19 values}\cr\cr
#' \code{maxpool}{ = 10000, Max depth for water pooling on the surface (mm)}\cr\cr
#' \code{rainmult}{ = 1, Rain multiplier for surface soil moisture (-)}\cr\cr
#' \code{evenrain}{ = 0, ignored -- rainfall always comes from BARRA's real hourly distribution (rainhourly is forced to 1)}\cr\cr
#' \code{SoilMoist_Init}{ = c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3), initial soil water content at each soil node, m3/m3}\cr\cr
#' \code{R1}{ = 0.001, root radius, m}\cr\cr
#' \code{RW}{ = 2.5e+10, resistance per unit length of root, m3 kg-1 s-1}\cr\cr
#' \code{RL}{ = 2e+6, resistance per unit length of leaf, m3 kg-1 s-1}\cr\cr
#' \code{PC}{ = -1500, critical leaf water potential for stomatal closure, J kg-1}\cr\cr
#' \code{SP}{ = 10, stability parameter for stomatal closure equation, -}\cr\cr
#' \code{IM}{ = 1e-06, maximum allowable mass balance error, kg}\cr\cr
#' \code{MAXCOUNT}{ = 500, maximum iterations for mass balance, -}\cr\cr
#' \code{LAI}{ = 0.1, leaf area index}\cr\cr
#'
#' \strong{ Snow mode parameters:}\cr\cr
#' \code{snowmodel}{ = 1, run the snow model 1=yes, 0=no}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (Mg/m3), overridden by densfun}\cr\cr
#' \code{densfun}{ = c(0.5979, 0.2178, 0.001, 0.0038), snow density model parameters}\cr\cr
#' \code{snowmelt}{ = 1, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
#' \code{undercatch}{ = 1, undercatch multiplier for converting rainfall to snow}\cr\cr
#' \code{rainmelt}{ = 0.0125, parameter in equation that melts snow with rainfall as a function of air temp}\cr\cr
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)}\cr\cr
#' \code{intercept}{ = max(maxshade) / 100 * 0.3, snow interception fraction for when there's shade (0-1)}\cr\cr
#' \code{grasshade}{ = 0, if 1, shade is removed when snow is present}\cr\cr
#' @examples
#' library(NicheMapR)
#' dstart <- "01/01/2020"
#' dfinish <- "31/01/2020"
#' loc <- c(133.8807, -23.6980) # Alice Springs, NT
#' micro <- micro_barra(loc = loc, dstart = dstart, dfinish = dfinish)
#'
#' metout <- as.data.frame(micro$metout)
#' soil <- as.data.frame(micro$soil)
#' soilmoist <- as.data.frame(micro$soilmoist)
#'
#' dates <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = "UTC"),
#'              as.POSIXct(dfinish, format = "%d/%m/%Y", tz = "UTC") + 23 * 3600, by = "hours")
#' metout <- cbind(dates, metout)
#' soil <- cbind(dates, soil)
#' soilmoist <- cbind(dates, soilmoist)
#'
#' plot(TALOC ~ dates, data = metout, type = "l", ylab = "Air temperature (C)")
#' plot(D0cm ~ dates, data = soil, type = "l", ylab = "Surface soil temperature (C)")
micro_barra <- function(
  loc = c(133.8807, -23.6980),
  dstart = "01/01/2020",
  dfinish = "31/12/2020",
  barra_product = "BARRAC2",
  barra_domain = "AUST04",
  elev = NA,
  REFL = 0.15,
  slope = 0,
  aspect = 0,
  DEP = c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
  minshade = 0,
  maxshade = 90,
  Usrhyt = 0.01,
  Z01 = 0,
  Z02 = 0,
  ZH1 = 0,
  ZH2 = 0,
  runshade = 1,
  run.gads = 1,
  solonly = 0,
  Soil_Init = NA,
  write_input = 0,
  writecsv = 0,
  windfac = 1,
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
  runmoist = 1,
  PE = rep(1.1, 19),
  KS = rep(0.0037, 19),
  BB = rep(4.5, 19),
  BD = rep(BulkDensity, 19),
  DD = rep(Density, 19),
  maxpool = 10000,
  rainmult = 1,
  evenrain = 0,
  SoilMoist_Init = c(0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.3, 0.3, 0.3, 0.3),
  L = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4, 0.4, 0, 0) * 10000,
  R1 = 0.001,
  RW = 2.5e+10,
  RL = 2e+06,
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
  lamb = 0,
  IUV = 0,
  ndmax = 3,
  IR = 0,
  message = 0,
  fail = NA,
  runmicro = 1,
  snowcond = 0,
  intercept = max(maxshade) / 100 * 0.3,
  grasshade = 0,
  maxsurf = 85
) { # end function parameters

  errors <- 0
  Refhyt <- 2 # Reference height (m) at which BARRA's tas/hurs/sfcWind are defined (near-surface = 2 m temp/RH, 10 m wind -- wind is height-corrected below)

  # ── error trapping ──────────────────────────────────────────────────────
  if (DEP[2] - DEP[1] > 3 | DEP[3] - DEP[2] > 3) {
    message("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if (DEP[2] - DEP[1] < 2) {
    message("warning, nodes might be too close near the surface, try a different spacing if the program is crashing \n")
  }
  if (DEP[10] != 200) {
    message("warning, last depth in soil should not be changed from 200 without good reason \n")
  }
  if (loc[1] > 180 | loc[1] < -180 | loc[2] > 90 | loc[2] < -90) {
    message("ERROR: longitude/latitude (loc) is out of bounds. \n")
    errors <- 1
  }
  if (!(barra_product %in% c("BARRAC2", "BARRAR2"))) {
    message("ERROR: barra_product must be 'BARRAC2' or 'BARRAR2'. \n")
    errors <- 1
  }
  if (barra_product == "BARRAC2" & barra_domain != "AUST04") {
    message("ERROR: barra_domain must be 'AUST04' for barra_product = 'BARRAC2'. \n")
    errors <- 1
  }
  if (barra_product == "BARRAR2" & !(barra_domain %in% c("AUS11", "AUST11"))) {
    message("ERROR: barra_domain must be 'AUS11' or 'AUST11' for barra_product = 'BARRAR2'. \n")
    errors <- 1
  }
  if (DEP[1] != 0) {
    message("ERROR: First soil node (DEP[1]) must = 0 cm. \n")
    errors <- 1
  }
  if (length(DEP) != 10) {
    message("ERROR: You must enter 10 different soil depths. \n")
    errors <- 1
  }
  for (i in 1:9) {
    if (DEP[i + 1] <= DEP[i]) {
      message("ERROR: Soil depth (DEP array) is not in ascending size \n")
      errors <- 1
    }
  }
  if (REFL < 0 | REFL > 1) {
    message("ERROR: Soil reflectivity value (REFL) is out of bounds (0-1). \n")
    errors <- 1
  }
  if (Usrhyt < RUF) {
    message("ERROR: Local height (Usrhyt) smaller than roughness height (RUF). \n")
    errors <- 1
  }
  if (Usrhyt > Refhyt) {
    message("ERROR: Reference height (2 m) is less than local height (Usrhyt) \n")
    errors <- 1
  }
  if (max(minshade - maxshade) >= 0) {
    message("ERROR: minshade must be less than maxshade. \n")
    errors <- 1
  }
  # end error trapping

  if (errors == 0) { # continue

    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("package 'ncdf4' is needed. Please install it.", call. = FALSE)
    }

    long <- as.numeric(loc[1])
    lat <- as.numeric(loc[2])
    HEMIS <- ifelse(lat < 0, 2, 1)
    ALAT <- abs(trunc(lat))
    AMINUT <- (abs(lat) - ALAT) * 60
    ALONG <- abs(trunc(long))
    ALMINT <- (abs(long) - ALONG) * 60
    ALREF <- abs(trunc(long))
    azmuth <- aspect

    ystart <- as.numeric(substr(dstart, 7, 10))
    yfinish <- as.numeric(substr(dfinish, 7, 10))
    nyears <- yfinish - ystart + 1
    if (is.na(fail)) fail <- nyears * 24 * 365

    tme <- seq(as.Date(dstart, format = "%d/%m/%Y"), as.Date(dfinish, format = "%d/%m/%Y"), "days")
    doy <- as.numeric(strftime(tme, format = "%j"))
    ndays <- length(doy)
    ida <- ndays
    idayst <- 1
    microdaily <- 1

    MINSHADES <- if (length(minshade) != ndays) rep(minshade[1], ndays) else minshade
    MAXSHADES <- if (runshade == 0) MINSHADES else
      if (length(maxshade) != ndays) rep(maxshade[1], ndays) else maxshade

    # ── BARRA OPeNDAP fetch (hourly, one file per variable per month) ──────
    # https://opus.nci.org.au/pages/viewpage.action?pageId=264241166
    barra_product_name <- switch(barra_product, BARRAC2 = "BARRA-C2", BARRAR2 = "BARRA-R2")
    barra_domain_name  <- switch(barra_domain,  AUST04 = "AUST-04", AUS11 = "AUS-11", AUST11 = "AUST-11")
    barra_base <- "https://thredds.nci.org.au/thredds/dodsC/ob53/output/reanalysis"

    barra_month_url <- function(layer, year, month, freq = "1hr") {
      yyyymm <- sprintf("%d%02d", year, month)
      tag <- paste0(layer, "_", barra_domain_name, "_ERA5_historical_hres_BOM_", barra_product_name, "_v1")
      fname <- paste0(tag, "_", freq, "_", yyyymm, "-", yyyymm, ".nc")
      paste0(barra_base, "/", barra_domain_name, "/BOM/ERA5/historical/hres/", barra_product_name,
             "/v1/", freq, "/", layer, "/latest/", fname)
    }
    barra_static_url <- function(layer) {
      tag <- paste0(layer, "_", barra_domain_name, "_ERA5_historical_hres_BOM_", barra_product_name, "_v1")
      paste0(barra_base, "/", barra_domain_name, "/BOM/ERA5/historical/hres/", barra_product_name,
             "/v1/fx/", layer, "/latest/", tag, ".nc")
    }

    # Reads one variable for one calendar month at the nearest grid cell to
    # (long, lat). Returns a numeric vector, one value per native timestep
    # (hourly for BARRA's "1hr" product).
    barra_read_month <- function(layer, year, month) {
      url <- barra_month_url(layer, year, month)
      nc <- ncdf4::nc_open(url)
      on.exit(ncdf4::nc_close(nc))
      lon_name <- if ("lon" %in% names(nc$dim)) "lon" else "longitude"
      lat_name <- if ("lat" %in% names(nc$dim)) "lat" else "latitude"
      lons <- nc$dim[[lon_name]]$vals
      lats <- nc$dim[[lat_name]]$vals
      loni <- which.min(abs(lons - long))
      lati <- which.min(abs(lats - lat))
      as.numeric(ncdf4::ncvar_get(nc, varid = layer, start = c(loni, lati, 1), count = c(1, 1, -1)))
    }

    barra_read_static <- function(layer) {
      url <- barra_static_url(layer)
      nc <- ncdf4::nc_open(url)
      on.exit(ncdf4::nc_close(nc))
      lon_name <- if ("lon" %in% names(nc$dim)) "lon" else "longitude"
      lat_name <- if ("lat" %in% names(nc$dim)) "lat" else "latitude"
      lons <- nc$dim[[lon_name]]$vals
      lats <- nc$dim[[lat_name]]$vals
      loni <- which.min(abs(lons - long))
      lati <- which.min(abs(lats - lat))
      as.numeric(ncdf4::ncvar_get(nc, varid = layer, start = c(loni, lati), count = c(1, 1)))
    }

    if (is.na(elev)) {
      message("fetching elevation from BARRA's own orog grid via OPeNDAP \n")
      elev <- barra_read_static("orog")
    }
    ALTITUDES <- elev

    year_months <- unique(format(tme, "%Y-%m"))
    message(paste0("extracting BARRA (", barra_product_name, " ", barra_domain_name,
                    ") weather data via OPeNDAP for ", length(year_months), " month(s) \n"))

    barra_vars <- c("tas", "hurs", "sfcWind", "rsds", "rlds", "pr")
    barra_hourly <- setNames(vector("list", length(barra_vars)), barra_vars)
    for (ym in year_months) {
      yr <- as.numeric(substr(ym, 1, 4))
      mo <- as.numeric(substr(ym, 6, 7))
      message(paste0("  ", ym, "\n"))
      for (v in barra_vars) {
        barra_hourly[[v]] <- c(barra_hourly[[v]], barra_read_month(v, yr, mo))
      }
    }
    nhrs <- ndays * 24
    # BARRA's monthly files may run slightly past dfinish (whole months) or
    # the very first file starts at the 1st of that month -- align to the
    # requested [dstart 00:00, dfinish 23:00] window by offsetting from the
    # first day of the first month.
    month_start <- as.Date(paste0(substr(dstart, 4, 10), "-01"), format = "%m/%Y-%d")
    hr_offset <- as.numeric(difftime(as.Date(dstart, format = "%d/%m/%Y"), month_start, units = "days")) * 24
    hr_idx <- (hr_offset + 1):(hr_offset + nhrs)
    for (v in barra_vars) {
      barra_hourly[[v]] <- barra_hourly[[v]][hr_idx]
    }

    TAIRhr <- barra_hourly$tas - 273.15
    RHhr   <- pmin(pmax(barra_hourly$hurs, 0.01), 100)
    if (requireNamespace("microclima", quietly = TRUE)) {
      WNhr <- microclima::windheight(barra_hourly$sfcWind * windfac, 10, 2)
    } else {
      WNhr <- barra_hourly$sfcWind * windfac * (2 / 10) ^ 0.15 # power-law fallback if microclima isn't installed
    }
    SOLRhr <- pmax(barra_hourly$rsds, 0)
    RAINhr <- barra_hourly$pr * 3600 # kg/m^2/s -> mm/hr

    # Cloud cover from the same clear-sky-ratio method micro_aust/micro_ncep
    # use for their daily CCMAXX/CCMINN, computed hourly instead -- BARRA
    # (like AGCD/AWAP) has no direct/diffuse split to derive cloud cover from
    # more directly the way mcera5 does for ERA5's fdir. Clear-sky irradiance
    # uses microclima::solalt (solar altitude) with the standard Meinel &
    # Meinel (1976) clear-sky attenuation -- NicheMapR has no R-level clear-
    # sky function of its own (SOLRAD is Fortran-internal only).
    if (!requireNamespace("microclima", quietly = TRUE)) {
      stop("package 'microclima' is needed for solar geometry (solalt/julday). Please install it.", call. = FALSE)
    }
    tzone <- "UTC"
    hours_seq <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = tzone),
                      as.POSIXct(dfinish, format = "%d/%m/%Y", tz = tzone) + 23 * 3600, by = "hours")
    jd_hr <- microclima::julday(as.numeric(format(hours_seq, "%Y")), as.numeric(format(hours_seq, "%m")), as.numeric(format(hours_seq, "%d")))
    hr_of_day <- as.numeric(format(hours_seq, "%H")) + as.numeric(format(hours_seq, "%M")) / 60
    alt <- microclima::solalt(hr_of_day, lat, long, jd_hr, merid = ALREF)
    sinalt <- sin(alt * pi / 180)
    clearsky_hr <- ifelse(sinalt > 0, 1361 * sinalt * 0.7 ^ ((1 / pmax(sinalt, 0.001)) ^ 0.678), 0)
    CLDhr <- ifelse(clearsky_hr > 1, (1 - SOLRhr / clearsky_hr) * 100, NA)
    CLDhr[!is.na(CLDhr) & CLDhr < 0] <- 0
    CLDhr[!is.na(CLDhr) & CLDhr > 100] <- 100
    CLDhr <- zoo::na.locf(zoo::na.locf(CLDhr, na.rm = FALSE), fromLast = TRUE) # night-time: hold nearest daylight value

    ZENhr <- rep(-1, nhrs) # let the Fortran solver compute zenith angle internally
    IRDhr <- rep(-1, nhrs) # let the Fortran solver compute sky IR internally (from CLDhr + humidity)

    # Daily summaries -- still required by the Fortran solver's microinput/
    # micro list alongside the hourly arrays above (used for e.g. deep-soil
    # boundary conditions), even though hourly = 1 drives the actual physics.
    day_of <- function(h) ((h - 1) %/% 24) + 1
    day_idx <- day_of(1:nhrs)
    TMAXX <- as.matrix(tapply(TAIRhr, day_idx, max))
    TMINN <- as.matrix(tapply(TAIRhr, day_idx, min))
    RHMAXX <- as.numeric(tapply(RHhr, day_idx, max))
    RHMINN <- as.numeric(tapply(RHhr, day_idx, min))
    WNMAXX <- as.numeric(tapply(WNhr, day_idx, max))
    WNMINN <- as.numeric(tapply(WNhr, day_idx, min))
    CCMAXX <- as.numeric(tapply(CLDhr, day_idx, max))
    CCMINN <- as.numeric(tapply(CLDhr, day_idx, min))
    RAINFALL <- as.numeric(tapply(RAINhr, day_idx, sum))
    tannul <- mean(TAIRhr)
    tannulrun <- rep(tannul, ndays)

    hourly <- 1
    rainhourly <- 1

    # ── GADS aerosol optical depth (unchanged from micro_aust.R -- location-
    #    only, not weather-source-dependent) ──────────────────────────────
    if (run.gads > 0) {
      relhum <- 1
      if (run.gads == 1) {
        optdep.summer <- as.data.frame(rungads(lat, long, relhum, 0))
        optdep.winter <- as.data.frame(rungads(lat, long, relhum, 1))
      } else {
        optdep.summer <- as.data.frame(gads.r(lat, long, relhum, 0))
        optdep.winter <- as.data.frame(gads.r(lat, long, relhum, 1))
      }
      optdep <- cbind(optdep.winter[, 1], rowMeans(cbind(optdep.summer[, 2], optdep.winter[, 2])))
      optdep <- as.data.frame(optdep)
      colnames(optdep) <- c("LAMBDA", "OPTDEPTH")
      a <- lm(OPTDEPTH ~ poly(LAMBDA, 6, raw = TRUE), data = optdep)
      LAMBDA <- c(290, 295, 300, 305, 310, 315, 320, 330, 340, 350, 360, 370, 380, 390, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320, 1380, 1400, 1420, 1440, 1460, 1480, 1500, 1540, 1580, 1600, 1620, 1640, 1660, 1700, 1720, 1780, 1800, 1860, 1900, 1950, 2000, 2020, 2050, 2100, 2120, 2150, 2200, 2260, 2300, 2320, 2350, 2380, 2400, 2420, 2450, 2490, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000)
      TAI <- predict(a, data.frame(LAMBDA))
    } else { # Australia-typical default (same fallback micro_aust.R uses)
      TAI <- c(0.0670358341290886, 0.0662612704779235, 0.065497075238002, 0.0647431301168489, 0.0639993178022531, 0.0632655219571553, 0.0625416272145492, 0.0611230843885423, 0.0597427855962549, 0.0583998423063099, 0.0570933810229656, 0.0558225431259535, 0.0545864847111214, 0.0533843764318805, 0.0522154033414562, 0.0499736739981675, 0.047855059159556, 0.0458535417401334, 0.0439633201842001, 0.0421788036108921, 0.0404946070106968, 0.0389055464934382, 0.0374066345877315, 0.0359930755919066, 0.0346602609764008, 0.0334037648376212, 0.0322193394032758, 0.0311029105891739, 0.0300505736074963, 0.0290585886265337, 0.0281233764818952, 0.0272415144391857, 0.0264097320081524, 0.0256249068083005, 0.0248840604859789, 0.0241843546829336, 0.0235230870563317, 0.0228976873502544, 0.0223057135186581, 0.0217448478998064, 0.0212128934421699, 0.0207077699817964, 0.0202275105711489, 0.0197702578594144, 0.0193342605242809, 0.0189178697551836, 0.0177713140039894, 0.0174187914242432, 0.0170790495503944, 0.0167509836728154, 0.0164335684174899, 0.0161258546410128, 0.0158269663770596, 0.0155360978343254, 0.0152525104459325, 0.0149755299703076, 0.0147045436435285, 0.0144389973831391, 0.0141783930434343, 0.0134220329447663, 0.0131772403830191, 0.0129356456025128, 0.0126970313213065, 0.0124612184223418, 0.0122280636204822, 0.01199745718102, 0.0115436048739351, 0.0110993711778668, 0.0108808815754663, 0.0106648652077878, 0.0104513876347606, 0.0102405315676965, 0.00982708969547694, 0.00962473896278535, 0.00903679230300494, 0.00884767454432418, 0.0083031278398166, 0.00796072474935954, 0.00755817587626185, 0.00718610751850881, 0.00704629977586921, 0.00684663903049612, 0.00654155580333479, 0.00642947339729728, 0.00627223096874308, 0.00603955966866779, 0.00580920937536261, 0.00568506186880564, 0.00563167068287251, 0.00556222005081865, 0.00550522989971023, 0.00547395763028062, 0.0054478983436216, 0.00541823364504573, 0.00539532163908382, 0.00539239864119488, 0.00541690124712384, 0.00551525885358836, 0.00564825853509463, 0.00577220185074264, 0.00584222986640171, 0.00581645238345584, 0.00566088137411449, 0.00535516862329704, 0.00489914757707667, 0.00432017939770409, 0.0036813032251836, 0.00309019064543606, 0.00270890436501562, 0.00276446109239711, 0.00356019862584603)
    }

    # ── Soil setup ───────────────────────────────────────────────────────
    Nodes <- matrix(data = 0, nrow = 10, ncol = ndays)
    Nodes[1:10, ] <- 1:10
    REFLS <- rep(REFL, ndays)
    SLES <- matrix(nrow = ndays, data = SLE)

    soilwet <- RAINFALL
    soilwet[soilwet <= rainwet] <- 0
    soilwet[soilwet > 0] <- 90
    PCTWET <- pmax(soilwet, PCTWET)

    soilprops <- matrix(data = 0, nrow = 10, ncol = 5)
    soilprops[, 1] <- BulkDensity
    soilprops[, 2] <- 1 - BulkDensity / Density
    soilprops[soilprops[, 2] < 0.26, 2] <- 0.26
    soilprops[, 3] <- Thcond
    soilprops[, 4] <- SpecHeat
    soilprops[, 5] <- Density
    if (cap == 1) {
      soilprops[1:2, 3] <- 0.2
      soilprops[1:2, 4] <- 1920
    }

    moists <- matrix(nrow = 10, ncol = ndays, data = 0)
    moists[1:10, ] <- SoilMoist_Init

    if (is.na(Soil_Init[1])) {
      soilinit <- rep(mean(TAIRhr), 20)
      spinup <- 1
    } else {
      soilinit <- c(Soil_Init, rep(mean(TAIRhr), 10))
      spinup <- 0
    }

    VIEWF <- 1 - sum(sin(as.data.frame(hori) * pi / 180)) / length(hori)
    TIMAXS <- c(1, 1, 0, 0)
    TIMINS <- c(0, 0, 1, 1)
    tides <- matrix(data = 0, nrow = 24 * ndays, ncol = 3)

    # ── Fortran call ─────────────────────────────────────────────────────
    microinput <- c(
      ndays, RUF, ERR, Usrhyt, Refhyt, 10,
      Z01, Z02, ZH1, ZH2, idayst, ida,
      HEMIS, ALAT, AMINUT, ALONG, ALMINT, ALREF,
      slope, azmuth, ALTITUDES, CMH2O, microdaily, tannul,
      EC, VIEWF, snowtemp, snowdens, snowmelt, undercatch, rainmult,
      1, runmoist, maxpool, evenrain, snowmodel, rainmelt, writecsv,
      densfun, hourly, rainhourly, lamb, IUV,
      RW, PC, RL, SP, R1, IM, MAXCOUNT, IR, message, fail,
      snowcond, intercept, grasshade, solonly, ZH, D0,
      TIMAXS, TIMINS, spinup, 0, 360, maxsurf, ndmax
    )

    micro <- list(
      tides = tides, microinput = microinput, doy = doy, SLES = SLES, DEP = DEP, Nodes = Nodes,
      MAXSHADES = MAXSHADES, MINSHADES = MINSHADES, TMAXX = TMAXX, TMINN = TMINN,
      RHMAXX = RHMAXX, RHMINN = RHMINN, CCMAXX = CCMAXX, CCMINN = CCMINN, WNMAXX = WNMAXX, WNMINN = WNMINN,
      TAIRhr = TAIRhr, RHhr = RHhr, WNhr = WNhr, CLDhr = CLDhr, SOLRhr = SOLRhr, RAINhr = RAINhr,
      ZENhr = ZENhr, IRDhr = IRDhr, REFLS = REFLS, PCTWET = PCTWET, soilinit = soilinit, hori = hori,
      TAI = TAI, soilprops = soilprops, moists = moists, RAINFALL = RAINFALL, tannulrun = tannulrun,
      PE = PE, KS = KS, BB = BB, BD = BD, DD = DD, L = L, LAI = rep(LAI, ndays)
    )

    if (write_input == 1) {
      if (!dir.exists("micro csv input")) dir.create("micro csv input")
      write.table(as.data.frame(microinput), file = "micro csv input/microinput.csv", sep = ",", col.names = NA, qmethod = "double")
    }

    if (runmicro == 1) {
      message("running microclimate model with BARRA hourly forcing \n")
      microut <- microclimate(micro)

      metout <- as.data.frame(microut$metout)
      shadmet <- as.data.frame(microut$shadmet)
      soil <- as.data.frame(microut$soil)
      shadsoil <- as.data.frame(microut$shadsoil)
      soilmoist <- as.data.frame(microut$soilmoist)
      shadmoist <- as.data.frame(microut$shadmoist)
      humid <- as.data.frame(microut$humid)
      shadhumid <- as.data.frame(microut$shadhumid)
      soilpot <- as.data.frame(microut$soilpot)
      shadpot <- as.data.frame(microut$shadpot)
      plant <- as.data.frame(microut$plant)
      shadplant <- as.data.frame(microut$shadplant)
      sunsnow <- as.data.frame(microut$sunsnow)
      shdsnow <- as.data.frame(microut$shdsnow)

      dates <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = "UTC"),
                   as.POSIXct(dfinish, format = "%d/%m/%Y", tz = "UTC") + 23 * 3600, by = "hours")
      dates2 <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = "UTC"),
                    as.POSIXct(dfinish, format = "%d/%m/%Y", tz = "UTC"), by = "days")

      return(list(soil = soil, shadsoil = shadsoil, metout = metout, shadmet = shadmet,
                  soilmoist = soilmoist, shadmoist = shadmoist, humid = humid, shadhumid = shadhumid,
                  soilpot = soilpot, shadpot = shadpot, plant = plant, shadplant = shadplant,
                  sunsnow = sunsnow, shdsnow = shdsnow,
                  RAINFALL = RAINFALL, ndays = ndays, elev = ALTITUDES, REFL = REFL,
                  longlat = c(long, lat), nyears = nyears, minshade = MINSHADES, maxshade = MAXSHADES,
                  DEP = DEP, dates = dates, dates2 = dates2, PE = PE, BD = BD, DD = DD, BB = BB, KS = KS))
    } else {
      return(list(RAINFALL = RAINFALL, TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN,
                  WNMAXX = WNMAXX, WNMINN = WNMINN, CCMAXX = CCMAXX, CCMINN = CCMINN,
                  CLDhr = CLDhr, WNhr = WNhr, TAIRhr = TAIRhr, RHhr = RHhr, RAINhr = RAINhr,
                  SOLRhr = SOLRhr, ZENhr = ZENhr, IRDhr = IRDhr, PE = PE, BD = BD, DD = DD, BB = BB, KS = KS))
    }
  } # end error check
}
