#' ERA5 implementation of the microclimate model, assuming ERA5 grids have been downloaded with package mcera5, and using package microclima downscaling for local topography.
#'
#' An implementation of the NicheMapR microclimate model that integrates the ERA5 hourly weather data and the elevatr package for obtaining DEM using downscaling functions from the microclima package, largely following the methods described in Kearney, M. R., Gillingham, P. K., Bramer, I., Duffy, J. P., & Maclean, I. M. D. (2019). A method for computing hourly, historical, terrain-corrected microclimate anywhere on Earth. Methods in Ecology and Evolution.
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param dstart First day to run, date in format "d/m/Y" e.g. "01/01/2016"
#' @param dfinish Last day to run, date in format "d/m/Y" e.g. "31/12/2016"
#' @param dem A digital elevation model used by microclima for micro-topographic effects, produced by microclima function 'get_dem' via R package 'elevatr' (internally generated via same function based on 'loc' if NA)
#' @param dem2 A digital elevation model used by microclima for meso-climate calculations, produced by microclima function 'get_dem' via R package 'elevatr' (internally generated via same function based on 'loc' if NA)
#' @param dem.res Requested resolution of the DEM from elevatr, m
#' @param zmin minimum elevation of DEM for terrain calculations, m (may need to be made negative if below sea level)
#' @param pixels Number of pixels along one edge of square requested of DEM requested from elevatr, #
#' @param REFL Soil solar reflectance, decimal \%
#' @param slope Slope in degrees (if NA, then derived from DEM with package microclima)
#' @param aspect Aspect in degrees (0 = north) (if NA, then derived from DEM with microclima)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param minshade Minimum shade level to use (can be a single value or a vector of daily values) (\%)
#' @param maxshade Maximum shade level to use (can be a single value or a vector of daily values) (\%)
#' @param Usrhyt Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
#' @param coastal Compute coastal effects with microclima? T (TRUE) or F (FALSE) (can take a while and may have high memory requirements depending on DEM size)
#' @param hourlydata user input of the hourlydata matrix
#' @param dailyprecip user input of daily rainfall
#' @param weather.elev optional value indicating the elevation of values in `hourlydata`. Either a numeric value, corresponding to the elevation in (m) of the location from which `hourlydata` were obtained, or `era5` (default, derived from Copernicus ERA5 climate reanalysis).
#' @param cad.effects optional logical indicating whether to calculate cold air drainage effects (TRUE = Yes, slower. FALSE =  No, quicker)
#' @param ... Additional arguments, see Details
#' @return metout The above ground micrometeorological conditions under the minimum specified shade
#' @return shadmet The above ground micrometeorological conditions under the maximum specified shade
#' @return soil Hourly predictions of the soil temperatures under the minimum specified shade
#' @return shadsoil Hourly predictions of the soil temperatures under the maximum specified shade
#' @return soilmoist Hourly predictions of the soil moisture under the minimum specified shade
#' @return shadmoist Hourly predictions of the soil moisture under the maximum specified shade
#' @return soilpot Hourly predictions of the soil water potential under the minimum specified shade
#' @return shadpot Hourly predictions of the soil water potential under the maximum specified shade
#' @return humid Hourly predictions of the soil humidity under the minimum specified shade
#' @return shadhumid Hourly predictions of the soil humidity under the maximum specified shade
#' @return plant Hourly predictions of plant transpiration, leaf water potential and root water potential under the minimum specified shade
#' @return shadplant Hourly predictions of plant transpiration, leaf water potential and root water potential under the maximum specified shade
#' @return sunsnow Hourly predictions of snow temperature under the minimum specified shade
#' @return shadsnow Hourly predictions snow temperature under the maximum specified shade
#' @usage micro_era5(loc = c(-91.415669, -0.287145), dstart = "01/01/2019", dfinish = "31/07/2019",
#' REFL = 0.15, slope = 0, aspect = 0, DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minshade = 0, maxshade = 90,
#' Usrhyt = 0.01, ...)
#' @export
#' @details
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes (Fortran version), 2=yes (R version), 0=no}\cr\cr
#' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1) or from ERA5 data (2)}\cr\cr
#' \code{solonly}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{ndmax}{ = 3, iterations of first day to get a steady periodic}\cr\cr
#' \code{Soil_Init}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{windfac}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{warm}{ = 0, warming offset vector, °C (negative values mean cooling). Can supply a single value or a vector the length of the number of days to be simulated.}\cr\cr
#' \code{scenario}{ = 0, TerraClimate climate change scenario, either 0, 2 or 4 °C warmer}\cr\cr
#' \code{terra_source}{ = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data", specify location of terraclimate data, goes to the web by default}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#' \code{spatial}{ = 'c:/era5_data/era5', specify folder and file prefix with local ERA5 data extracted via the mcera5 package (no trailing forward slash)}\cr\cr
#' \code{save}{ = 0, don't save forcing data (0), save the forcing data (1) or read previously saved data (2)}\cr\cr
#' \code{runmicro}{ = 1, call the microclimate model (1) or not (0), if you just want the downscaled input weather data}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{ERR}{ = 1, Integrator error tolerance for soil temperature calculations}\cr\cr
#' \code{RUF}{ = 0.004, Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001 (snow) - 0.02 m.}\cr\cr
#' \code{ZH}{ = 0, heat transfer roughness height (m) for Campbell and Norman air temperature/wind speed profile (invoked if greater than 0, 0.02 * canopy height in m if unknown)}\cr\cr
#' \code{D0}{ = 0, zero plane displacement correction factor (m) for Campbell and Norman air temperature/wind speed profile (0.6 * canopy height in m if unknown)}\cr\cr
#' \code{Z01}{ = 0, Top (1st) segment roughness height(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then RUF and Refhyt used)}\cr\cr
#' \code{Z02}{ = 0, 2nd segment roughness height(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then RUF and Refhyt used).}\cr\cr
#' \code{ZH1}{ = 0, Top of (1st) segment, height above surface(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then RUF and Refhyt used).}\cr\cr
#' \code{ZH2}{ = 0, 2nd segment, height above surface(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then RUF and Refhyt used).}\cr\cr
#' \code{EC}{ = 0.0167238, Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)}\cr\cr
#' \code{SLE}{ = 0.95, Substrate longwave IR emissivity (decimal \%), typically close to 1}\cr\cr
#' \code{Thcond}{ = 2.5, Soil minerals thermal conductivity, single value or vector of 10 specific to each depth (W/mK)}\cr\cr
#' \code{Density}{ = 2.56, Soil minerals density, single value or vector of 10 specific to each depth (Mg/m3)}\cr\cr
#' \code{SpecHeat}{ = 870, Soil minerals specific heat, single value or vector of 10 specific to each depth (J/kg-K)}\cr\cr
#' \code{BulkDensity}{ = 1.3, Soil bulk density (Mg/m3), single value or vector of 10 specific to each depth}\cr\cr
#' \code{PCTWET}{ = 0, \% of ground surface area acting as a free water surface (overridden if soil moisture model is running)}\cr\cr
#' \code{rainwet}{ = 1.5, mm of rainfall causing the ground to be 90\% wet for the day}\cr\cr
#' \code{cap}{ = 1, organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)}\cr\cr
#' \code{CMH2O}{ = 1, Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air conditions; 2.0 = humid, tropical conditions (note this is for the whole atmospheric profile, not just near the ground)}\cr\cr
#' \code{hori}{ = rep(NA,24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#'
#' \strong{ Soil moisture mode parameters:}\cr\cr
#'
#' \code{runmoist}{ = 1, Run soil moisture model? 1=yes, 0=no  1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{PE}{ = rep(1.1,19), Air entry potential (J/kg) (19 values descending through soil for specified soil nodes in parameter}
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{KS}{ = rep(0.0037,19), Saturated conductivity, (kg s/m3) (19 values descending through soil for specified soil nodes in parameter}
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{BB}{ = rep(4.5,19), Campbell's soil 'b' parameter (-) (19 values descending through soil for specified soil nodes in parameter}
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{BD}{ = rep(1.3,19), Soil bulk density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter}
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{DD}{ = rep(2.56,19), Soil density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)}\cr\cr
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{maxpool}{ = 10000, Max depth for water pooling on the surface (mm), to account for runoff}\cr\cr
#' \code{rainhourly}{ = 0, Is hourly rain input being supplied (1 = yes, 0 = no)?}\cr\cr
#' \code{rainhour}{ = 0, Vector of hourly rainfall values - overrides daily ERA5 rain if rainhourly = 1}\cr\cr
#' \code{rainmult}{ = 1, Rain multiplier for surface soil moisture (-), used to induce runon}\cr\cr
#' \code{rainoff}{ = 0, Rain offset (mm), used to induce changes in rainfall from ERA5 values. Can be a single value or a vector matching the number of days to simulate. If negative values are used, rainfall will be prevented from becomming negative.}\cr\cr
#' \code{evenrain}{ = 0, Spread daily rainfall evenly across 24hrs (1) or one event at midnight (0)}\cr\cr
#' \code{SoilMoist_Init}{ = c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3), initial soil water content at each soil node, m3/m3}\cr\cr
#' \code{L}{ = c(0,0,8.2,8.0,7.8,7.4,7.1,6.4,5.8,4.8,4.0,1.8,0.9,0.6,0.8,0.4,0.4,0,0)*10000, root density (m/m3), (19 values descending through soil for specified soil nodes in parameter}\cr\cr
#' \code{R1}{ = 0.001, root radius, m}\cr\cr
#' \code{RW}{ = 2.5e+10, resistance per unit length of root, m3 kg-1 s-1}\cr\cr
#' \code{RL}{ = 2e+6, resistance per unit length of leaf, m3 kg-1 s-1}\cr\cr
#' \code{PC}{ = -1500, critical leaf water potential for stomatal closure, J kg-1}\cr\cr
#' \code{SP}{ = 10, stability parameter for stomatal closure equation, -}\cr\cr
#' \code{IM}{ = 1e-06, maximum allowable mass balance error, kg}\cr\cr
#' \code{MAXCOUNT}{ = 500, maximum iterations for mass balance, -}\cr\cr
#' \code{LAI}{ = 0.1, leaf area index (can be a single value or a vector of daily values), used to partition traspiration/evaporation from PET in soil moisture model}\cr\cr
#' \code{microclima.LAI}{ = 0, leaf area index, used by package microclima for radiation calcs}\cr\cr
#' \code{LOR}{ = 1, leaf orientation for package microclima radiation calcs}\cr\cr
#'
#' \strong{ Snow mode parameters:}
#'
#' \code{snowmodel}{ = 1, run the snow model 1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (Mg/m3), overridden by densfun}\cr\cr
#' \code{densfun}{ = c(0.5979, 0.2178, 0.001, 0.0038), slope and intercept of model of snow density as a linear function of snowpack age if first two values are nonzero, and following the exponential function of Sturm et al. 2010 J. of Hydromet. 11:1380-1394 if all values are non-zero; if it is c(0,0,0,0) then fixed density used}\cr\cr
#' \code{snowmelt}{ = 1, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
#' \code{undercatch}{ = 1, undercatch multipier for converting rainfall to snow}\cr\cr
#' \code{rainmelt}{ = 0.0125, paramter in equation that melts snow with rainfall as a function of air temp}\cr\cr
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)}\cr\cr
#' \code{intercept}{ = max(maxshade) / 100 * 0.3, snow interception fraction for when there's shade (0-1)}\cr\cr
#' \code{grasshade}{ = 0, if 1, means shade is removed when snow is present, because shade is cast by grass/low shrubs}\cr\cr
#'
#' \strong{ Intertidal mode parameters:}
#'
#' \code{shore}{ Include tide effects? If 1, the matrix}
#' \code{tides}{ is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0, nrow = length(seq(as.POSIXct(dstart, format = '%d/%m/%Y'), as.POSIXct(dfinish, format = '%d/%m/%Y'), by = 'days')) * 24, ncol = 3), matrix of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#' }
#'
#' \strong{Outputs:}
#'
#' \code{ndays}{ - number of days for which predictions are made}\cr\cr
#' \code{longlat}{ - longitude and latitude for which simulation was run (decimal degrees)}\cr\cr
#' \code{dates}{ - vector of dates (POSIXct, UTC)}\cr\cr
#' \code{nyears}{ - number of years for which predictions are made}\cr\cr
#' \code{RAINFALL}{ - vector of daily rainfall (mm)}\cr\cr
#' \code{elev}{ - elevation at point of simulation (m)}\cr\cr
#' \code{minshade}{ - minimum shade for simulation (\%)}\cr\cr
#' \code{maxshade}{ - maximum shade for simulation (\%)}\cr\cr
#' \code{dem}{ - digital elevation model obtained via 'get_dev' using package 'elevatr' (m)}\cr\cr
#' \code{DEP}{ - vector of depths used (cm)}\cr\cr
#' \code{SLOPE}{ - slope at point of simulation (\%)}\cr\cr
#' \code{ASPECT}{ - aspect at point of simulation (°, 0 is north)}\cr\cr
#' \code{HORIZON}{ - horizon angles at point of simulation (°)}\cr\cr
#'
#' metout/shadmet variables:
#' \itemize{
#' \item 1 DOY - day-of-year
#' \item 2 TIME - time of day (mins)
#' \item 3 TALOC - air temperature (°C) at local height (specified by 'Usrhyt' variable)
#' \item 4 TAREF - air temperature (°C) at reference height (specified by 'Refhyt', 2m default)
#' \item 5 RHLOC - relative humidity (\%) at local height (specified by 'Usrhyt' variable)
#' \item 6 RH  - relative humidity (\%) at reference height (specified by 'Refhyt', 2m default)
#' \item 7 VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)
#' \item 8 VREF - wind speed (m/s) at reference height (specified by 'Refhyt', 2m default)
#' \item 9 SNOWMELT - snowmelt (mm)
#' \item 10 POOLDEP - water pooling on surface (mm)
#' \item 11 PCTWET - soil surface wetness (\%)
#' \item 12 ZEN - zenith angle of sun (degrees - 90 = below the horizon)
#' \item 13 SOLR - solar radiation (W/m2) (unshaded, adjusted for slope, aspect and horizon angle)
#' \item 14 TSKYC - sky radiant temperature (°C)
#' \item 15 DEW - dew fall (mm / h)
#' \item 16 FROST - frost (mm / h)
#' \item 17 SNOWFALL - snow predicted to have fallen (cm)
#' \item 18 SNOWDEP - predicted snow depth (cm)
#' \item 19 SNOWDENS - snow density (g/cm3)
#'}
#' soil and shadsoil variables:
#' \itemize{
#' \item 1 DOY - day-of-year
#' \item 2 TIME - time of day (mins)
#' \item 3-12 D0cm ... - soil temperature (°C) at each of the 10 specified depths
#' }
#'
#' if soil moisture model is run i.e. parameter runmoist = 1\cr
#'
#' soilmoist and shadmoist variables:
#' \itemize{
#' \item 1 DOY - day-of-year
#' \item 2 TIME - time of day (mins)
#' \item 3-12 WC0cm ... - soil moisture (m3/m3) at each of the 10 specified depths
#' }
#' soilpot and shadpot variables:
#' \itemize{
#' \item 1 DOY - day-of-year
#' \item 2 TIME - time of day (mins)
#' \item 3-12 PT0cm ... - soil water potential (J/kg = kPa = bar/100) at each of the 10 specified depths
#' }
#' humid and shadhumid variables:
#' \itemize{
#' \item  1 DOY - day-of-year
#' \item  2 TIME - time of day (mins)
#' \item  3-12 RH0cm ... - soil relative humidity (decimal \%), at each of the 10 specified depths
#' }
#' plant and shadplant variables:
#' \itemize{
#' \item  1 DOY - day-of-year
#' \item  2 TIME - time of day (mins)
#' \item  3 TRANS - plant transpiration rate (g/m2/h)
#' \item  4 LEAFPOT - leaf water potential (J/kg = kPa = bar/100)
#' \item  5-14 RPOT0cm ... - root water potential (J/kg = kPa = bar/100), at each of the 10 specified depths
#' }
#'
#' if snow model is run i.e. parameter snowmodel = 1\cr
#' sunsnow and shdsnow variables:
#' \itemize{
#' \item  1 DOY - day-of-year
#' \item  2 TIME - time of day (mins)
#' \item  3-10 SN1 ... - snow temperature (°C), at each of the potential 8 snow layers (layer 8 is always the bottom - need metout$SNOWDEP to interpret which depth in the snow a given layer represents)
#' }
#'
#' if wavelength-specific solar output is selected i.e. parameter lamb = 1\cr
#' solar output variables
#' drlam (direct solar), drrlam (direct Rayleigh solar) and srlam (scattered solar) variables:
#' \itemize{
#' \item  1 DOY - day-of-year
#' \item  2 TIME - time of day (mins)
#' \item  3-113 290, ..., 4000 - irradiance (W/(m2 nm)) at each of 111 wavelengths from 290 to 4000 nm
#' }
#' @examples
#' library(NicheMapR)
#' library(ecmwfr)
#' library(mcera5)
#' library(lubridate)
#' library(dplyr)
#' library(tidync)
#'
#' # get ERA5 data with package mcera5 (just do once for region and time of interest)
#'
#' # assign your credentials (register here: https://cds.climate.copernicus.eu/user/register)
#' uid <- "$$$$$$"
#' cds_api_key <- "$$$$$$$$-$$$$-$$$$-$$$$-$$$$$$$$$$$$"
#'
#' ecmwfr::wf_set_key(user = uid, key = cds_api_key)
#'
#' # bounding coordinates (in WGS84 / EPSG:4326)
#' xmn <- 130
#' xmx <- 132
#' ymn <- -26
#' ymx <- -24
#'
#' # temporal extent
#' st_time <- lubridate::ymd("2010:07:01")
#' en_time <- lubridate::ymd("2011:12:31")
#'
#' # filename and location for downloaded .nc files
#' file_prefix <- "era5"
#' op <- "C:/Spatial_Data/"
#'
#' # build a request (covering multiple years)
#' req <- build_era5_request(xmin = xmn, xmax = xmx,
#'                           ymin = ymn, ymax = ymx,
#'                           start_time = st_time,
#'                           end_time = en_time,
#'                           outfile_name = file_prefix)
#' str(req)
#' request_era5(request = req, uid = uid, out_path = op)
#'
#' # run micro_era5 for a location (make sure it's within the bounds of your .nc files)
#'
#' dstart <- "01/01/2011"
#' dfinish <- "31/12/2011"
#' loc <- c(131, -25) # somewhere in the middle of Australia
#' micro<-micro_era5(loc = loc, dstart = dstart, dfinish = dfinish, spatial = 'c:/Spatial_Data/era5')
#'
#' metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#' soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
#' soilmoist<-as.data.frame(micro$soilmoist) # soil temperatures, minimum shade
#'
#' # append dates
#' tzone<-paste("Etc/GMT+",0,sep="")
#' dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")
#'
#' metout <- cbind(dates,metout)
#' soil <- cbind(dates,soil)
#' soilmoist <- cbind(dates, soilmoist)
#'
#' # plotting above-ground conditions in minimum shade
#' with(metout,{plot(TALOC ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l",main=paste("air and sky temperature",sep=""), ylim = c(-20, 60))})
#' with(metout,{points(TAREF ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l",lty=2,col='blue')})
#' with(metout,{points(TSKYC ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' ,  type = "l",col='light blue',main=paste("sky temperature",sep=""))})
#' with(metout,{plot(RHLOC ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l",ylim=c(0,100),main=paste("humidity",sep=""))})
#' with(metout,{points(RH ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l",col='blue',lty=2,ylim=c(0,100))})
#' with(metout,{plot(VREF ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#' ,  type = "l",main="wind speed",ylim = c(0, 15))})
#' with(metout,{points(VLOC ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#' ,  type = "l",lty=2,col='blue')})
#' with(metout,{plot(SOLR ~ dates,xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
#' ,  type = "l",main="solar radiation")})
#' with(metout,{plot(SNOWDEP ~ dates,xlab = "Date and Time", ylab = "Snow Depth (cm)"
#' ,  type = "l",main="snow depth")})
#'
#' # plotting soil temperature
#' for(i in 1:10){
#'  if(i==1){
#'    plot(soil[,i+3]~soil[,1],xlab = "Date and Time", ylab = "Soil Temperature (°C)"
#'    ,col=i,type = "l",main=paste("soil temperature",sep=""))
#'  }else{
#'    points(soil[,i+3]~soil[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)",col=i,type = "l")
#'  }
#' }
#'
#' # plotting soil moisture
#' for(i in 1:10){
#'  if(i==1){
#'    plot(soilmoist[,i+3]*100~soilmoist[,1],xlab = "Date and Time", ylab = "Soil Moisture (% volumetric)"
#'    ,col=i,type = "l",main=paste("soil moisture",sep=""))
#'  }else{
#'    points(soilmoist[,i+3]*100~soilmoist[,1],xlab = "Date and Time", ylab = "Soil Moisture
#'     (%)",col=i,type = "l")
#'  }
#' }
micro_era5 <- function(
  loc = c(-5.3, 50.13),
  dstart = "01/01/2017",
  dfinish = "31/12/2017",
  dem = NA,
  dem2 = dem,
  dem.res = 30,
  zmin = 0,
  pixels = 100,
  nyears = as.numeric(substr(dfinish, 7, 10)) - as.numeric(substr(dstart, 7, 10)) + 1,
  REFL = 0.15,
  slope = NA,
  aspect = NA,
  DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200),
  minshade = 0,
  maxshade = 90,
  Usrhyt = 0.01,
  Z01 = 0,
  Z02 = 0,
  ZH1 = 0,
  ZH2 = 0,
  runshade = 1,
  clearsky = 0,
  run.gads = 1,
  solonly = 0,
  Soil_Init = NA,
  write_input = 0,
  writecsv = 0,
  windfac = 1,
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
  hori = rep(NA, 24),
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
  L = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000,
  R1 = 0.001,
  RW = 2.5e+10,
  RL = 2e+06,
  PC = -1500,
  SP = 10,
  IM = 1e-06,
  MAXCOUNT = 500,
  LAI = 0.1,
  microclima.LAI = 0,
  LOR = 1,
  snowmodel = 1,
  snowtemp = 1.5,
  snowdens = 0.375,
  densfun = c(0.5979, 0.2178, 0.001, 0.0038),
  snowmelt = 1,
  undercatch = 1,
  rainmelt = 0.0125,
  shore = 0,
  tides = 0,
  deepsoil = NA,
  rainhour = 0,
  rainhourly = 0,
  rainoff = 0,
  lamb = 0,
  IUV = 0,
  ndmax = 3,
  soilgrids = 0,
  IR = 0,
  message = 0,
  fail = nyears * 24 * 365,
  spatial = 'c:/era5_data/era5',
  save = 0,
  runmicro = 1,
  snowcond = 0,
  intercept = max(maxshade) / 100 * 0.3,
  grasshade = 0,
  scenario = 0,
  terra_source = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data",
  coastal = F,
  hourlydata = NA,
  dailyprecip = NA,
  weather.elev = 'era5',
  cad.effects = TRUE,
  dewrain = 0,
  moiststep = 360,
  maxsurf = 85){ # end function parameters

  # --- 1. Input validation ---
  errors <- 0
  Refhyt <- 2 # Reference height (m): ERA5 data is at 2 m screen height

  # era5-specific checks
  if (is.numeric(loc[1]) && (loc[1] > 180 | loc[2] > 90)) {
    message("ERROR: Latitude or longitude (longlat) is out of bounds. Please enter a correct value.")
    errors <- 1
  }
  if (run.gads == 1) {
    message("If program is crashing, try run.gads = 2.")
  }
  if (!(write_input %in% c(0, 1))) {
    message("ERROR: the variable 'write_input' must be either 0 or 1. Please correct.")
    errors <- 1
  }
  if (!(scenario %in% c(0, 2, 4))) {
    message("ERROR: Scenario can only be 0, 2, or 4 corresponding to the TerraClimate climate change scenarios.")
    errors <- 1
  }

  errors <- errors + validate_micro_inputs(DEP, REFL,
    slope  = ifelse(is.na(slope),  0, slope),
    aspect = ifelse(is.na(aspect), 0, aspect),
    hori   = if (is.na(hori[1])) rep(0, 24) else hori,
    SLE, ERR, RUF, D0, Usrhyt, Refhyt, EC, CMH2O,
    TIMAXS = c(1, 1, 0, 0), TIMINS = c(0, 0, 1, 1),
    minshade, maxshade, Thcond, Density, SpecHeat, BulkDensity,
    run.gads = run.gads)

  if(errors==0){ # continue
    max.date <- as.Date(paste0("01/", format(Sys.time(), "%m/%Y")), format = "%d/%m/%Y", tz = "UTC")
    if(as.Date(dfinish, format = "%d/%m/%Y") > max.date){
      message(paste0("Cannot simulate past ", max.date, "; reducing timespan accordingly \n"))
      dfinish <- as.character(as.Date(max.date, format = "%Y/%m/%d"))
      dfinish <- paste(substr(dfinish, 9, 10), substr(dfinish, 6, 7), substr(dfinish, 1, 4), sep = "/")
    }
    ystart <- as.numeric(substr(dstart, 7, 10))
    yfinish <- as.numeric(substr(dfinish, 7, 10))
    yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)

    #remove trailing forward slash if necessary
    if(is.na(spatial)==FALSE){
      if(substr(x = spatial, start = nchar(spatial), nchar(spatial))=='/'){
        spatial <- substr(spatial, 1, nchar(spatial)-1)
      }
    }
    ################## time related variables #################################

    # for microclima calculations
    tme <- seq(as.Date(dstart, format = "%d/%m/%Y"), as.Date(dfinish, format = "%d/%m/%Y"), "days")

    doy <- as.numeric(strftime(tme, format = "%j"))
    ndays<-length(doy)
    doynum<-ndays
    ida<-ndays
    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial ndmax day burn in
    daystart<-1

    if(length(minshade) != ndays){
      MINSHADES <- rep(0, ndays) + minshade[1] # daily min shade (%)
    }else{
      MINSHADES <- rep(0, ndays) + minshade # daily min shade (%)
    }
    if(length(maxshade) != ndays){
      MAXSHADES <- rep(0, ndays) + maxshade[1] # daily max shade (%)
    }else{
      MAXSHADES <- rep(0, ndays) + maxshade # daily max shade (%)
    }
    idayst <- 1 # start day

    ################## location and terrain #################################
    if (!require("terra", quietly = TRUE)) {
      stop("package 'terra' is needed. Please install it.",
           call. = FALSE)
    }
    if (!require("mcera5", quietly = TRUE)) {
      stop("package 'mcera5' is needed. Please install it.",
           call. = FALSE)
    }
    if (!require("RNetCDF", quietly = TRUE)) {
      stop("package 'RNetCDF' is needed. Please install it.",
           call. = FALSE)
    }
    if (!require("microclima", quietly = TRUE)) {
      stop("package 'microclima' is needed. Please install it via command: devtools::install_github('ilyamaclean/microclima').",
           call. = FALSE)
    }
    require("terra")
    require("mcera5")
    require("RNetCDF")
    require("microclima")
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    # get the local timezone reference longitude
    ALREF <- abs(trunc(x[1]))
    HEMIS <- ifelse(x[2]<0, 2, 1) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2])-ALAT)*60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1])-ALONG)*60
    azmuth<-aspect
    lat <- as.numeric(longlat[2])
    long <- as.numeric(longlat[1])
    loc <- c(long, lat)
    if(class(dem)[1] == "SpatRaster"){
      cat('using DEM provided to function call \n')
    }
    if(save != 2 & class(dem)[1] != "SpatRaster"){
      cat('downloading DEM via package elevatr \n')
      dem <- microclima::get_dem(lat = lat, long = long, resolution = dem.res, xdims = pixels, ydims = pixels) # mercator equal area projection
    }
    if(save == 1){
      save(dem, file = 'dem.Rda')
    }
    if(save == 2){
      load('dem.Rda')
    }
    if(save != 2){
      if(soilgrids == 1){
        sg <- fetch_soilgrids(x, DEP)
        if (!is.null(sg)) { PE <- sg$PE; KS <- sg$KS; BB <- sg$BB; BD <- sg$BD; BulkDensity <- sg$BulkDensity }
      }
    }else{
      if(soilgrids == 1){
        cat("loading SoilGrids data from previous run \n")
        load('PE.Rda')
        load('BB.Rda')
        load('BD.Rda')
        load('KS.Rda')
        load('BulkDensity.Rda')
      }
    }
    if(save == 1 & soilgrids == 1){
      cat("saving SoilGrids data for later \n")
      save(PE, file = 'PE.Rda')
      save(BB, file = 'BB.Rda')
      save(BD, file = 'BD.Rda')
      save(KS, file = 'KS.Rda')
      save(BulkDensity, file = 'BulkDensity.Rda')
    }
    if(is.na(hori[1])){
      hori<-rep(0, 24)
      VIEWF <- 1 # incorporated already by microclima
    }else{
      VIEWF <- 1-sum(sin(as.data.frame(hori)*pi/180))/length(hori) # convert horizon angles to radians and calc view factor(s)
      HORIZON <- spline(x = hori, n = 36, method =  'periodic')$y
      HORIZON[HORIZON < 0] <- 0
      HORIZON[HORIZON > 90] <- 90
    }
    days <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), by = 'days')
    alldays <- seq(as.POSIXct("01/01/1900", format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), Sys.time()-60*60*24, by = 'days')
    startday <- which(as.character(format(alldays, "%d/%m/%Y")) == format(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), "%d/%m/%Y"))
    endday <- which(as.character(format(alldays, "%d/%m/%Y")) == format(as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), "%d/%m/%Y"))
    countday <- endday-startday+1
    dates <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(dfinish, format = "%d/%m/%Y", tz = 'UTC')+23*3600, by = 'hours')
    dates2 <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(dfinish, format = "%d/%m/%Y", tz = 'UTC')+23*3600, by = 'days')
    if(save == 2){
      load('tref.Rda')
      load('SLOPE.Rda')
      load('ASPECT.Rda')
      load('HORIZON.Rda')
      elev <- tref$elev[1] # m
      ALTT <- elev
    }

    if(save != 2){

      # now getting starting point and count for reading netcdf files
      cat(paste0("extracting weather data locally from ", spatial, " \n"))
      dstart.splt <- strsplit(dstart, '/')
      dfinish.splt <- strsplit(dfinish, '/')
      st_time <- lubridate::ymd(paste0(dstart.splt[[1]][3], ":", dstart.splt[[1]][2], ":", dstart.splt[[1]][1]))
      en_time <- lubridate::ymd(paste0(dfinish.splt[[1]][3], ":", dfinish.splt[[1]][2], ":", dfinish.splt[[1]][1]))
      years <- as.numeric(unique(format(tme, "%Y")))
      longitude <- loc[1]
      latitude <- loc[2]
      if(is.logical(hourlydata)){
        if(length(years)==1){
          hourlydata <- mcera5::extract_clim(nc = paste0(spatial, '_', years, '.nc'), long = loc[1], lat = loc[2],
                                             start_time = st_time,  end_time = en_time)
        }else{
          for(i in 1:length(years)){
            if(i==1){
              hourlydata <- mcera5::extract_clim(nc = paste0(spatial, '_', years[i], '.nc'), long = loc[1], lat = loc[2],
                                                 start_time = st_time,
                                                 end_time = as.Date(paste(years[i],'12','31', sep='-')))
            }
            if(i!=1 & i!=length(years)){
              hourlydata.i <- mcera5::extract_clim(nc = paste0(spatial, '_', years[i], '.nc'), long = loc[1], lat = loc[2],
                                                   start_time = as.Date(paste(years[i],'01','01', sep='-')),
                                                   end_time = as.Date(paste(years[i],'12','31', sep='-')))
              hourlydata <- dplyr::bind_rows(hourlydata, hourlydata.i)
            }
            if(i==length(years)){
              hourlydata.i <- mcera5::extract_clim(nc = paste0(spatial, '_', years[i], '.nc'), long = loc[1], lat = loc[2],
                                                   start_time = as.Date(paste(years[i],'01','01', sep='-')),
                                                   end_time = en_time)
              hourlydata <- dplyr::bind_rows(hourlydata, hourlydata.i)
            }
          }
        }
      }
      # gather daily precipitation
      if(is.logical(dailyprecip)){
        if(length(years)==1){
          dailyprecip <- mcera5::extract_precip(nc = paste0(spatial, '_', years, '.nc'), long = loc[1], lat = loc[2],
                                                start_time = st_time,
                                                end_time = en_time)
        }else{
          for(i in 1:length(years)){
            if(i==1){
              dailyprecip <- mcera5::extract_precip(nc = paste0(spatial, '_', years[i], '.nc'), long = loc[1], lat = loc[2],
                                                    start_time = st_time,
                                                    end_time = as.Date(paste(years[i],'12','31', sep='-')))
            }
            if(i!=1 & i!=length(years)){
              dailyprecip.i <- mcera5::extract_precip(nc = paste0(spatial, '_', years[i], '.nc'), long = loc[1], lat = loc[2],
                                                      start_time = as.Date(paste(years[i],'01','01', sep='-')),
                                                      end_time = as.Date(paste(years[i],'12','31', sep='-')))
              dailyprecip <- c(dailyprecip, dailyprecip.i)
            }
            if(i==length(years)){
              dailyprecip.i <- mcera5::extract_precip(nc = paste0(spatial, '_', years[i], '.nc'), long = loc[1], lat = loc[2],
                                                      start_time = as.Date(paste(years[i],'01','01', sep='-')),
                                                      end_time = en_time)
              dailyprecip <- c(dailyprecip, dailyprecip.i)
            }
          }
        }
      }
      cat("computing radiation and elevation effects with package microclima \n")
      microclima.out <- microclima::microclimaforNMR(lat = longlat[2],
                                                     long = longlat[1],
                                                     dstart = dstart,
                                                     dfinish = dfinish,
                                                     l = mean(microclima.LAI),
                                                     x = LOR,
                                                     coastal = coastal,
                                                     hourlydata = as.data.frame(hourlydata),
                                                     dailyprecip = dailyprecip,
                                                     dem = dem,
                                                     demmeso = dem2,
                                                     albr = 0,
                                                     resolution = 30,
                                                     slope = slope,
                                                     aspect = aspect,
                                                     windthresh = 4.5,
                                                     emthresh = 0.78,
                                                     reanalysis2 = TRUE,
                                                     difani = FALSE,
                                                     weather.elev = weather.elev,
                                                     cad.effects = cad.effects,
                                                     zmin = zmin)


      hourlyradwind <- microclima.out$hourlyradwind
      SLOPE <- hourlyradwind$slope[1]
      ASPECT <- hourlyradwind$aspect[1]
      tref <- microclima.out$tref
      ZENhr <- hourlydata$szenith
      ZENhr[ZENhr > 90] <- 90
      HORIZON <- hori
      if(save == 1){
        save(SLOPE, file = 'SLOPE.Rda')
        save(ASPECT, file = 'ASPECT.Rda')
        save(HORIZON, file = 'HORIZON.Rda')
      }
      # NB units for rad = MJ / m^2 / hr (divide by 0.0036 to get to W / m^2)
      # Skyviewfact (time invariant, 1 = complete hemisphere in view)
      # canopyfact (proportion of isotropic radiation blocked out, 1 = no radiation gets in)
      if(save == 1){
        save(tref, file = 'tref.Rda')
      }
      # tref (original hourly ERA5)
      # telev (elevation-corrected temperature)
      # tcad (delta temperature due to cold air drainage)
      elev <- tref$elev[1] # m
      ALTT <- elev
      TAIRhr <- tref$telev + tref$tcad # reference Tair corrected for lapse rate and cold air drainage
      if(scenario != 0){
        TAIRhr_orig <- TAIRhr
        yearstodo <- seq(ystart, yfinish)
        nyears <- yfinish - ystart + 1
        if(yfinish > 2025){
          ystart_terra <- 2025 - nyears + 1
          yfinish_terra <- 2025
          message(paste0("terraclimate climate change scenarios are for 1950 to 2025 - using ", ystart, "-", yfinish), '\n')
        }else{
          ystart_terra <- ystart
          yfinish_terra <- yfinish
        }
        leapyears <- seq(1900, 2060, 4)
        message("generate climate change scenario", '\n')
        # diff spline function
        getdiff <- function(diffs){
          diff1 <- (unlist(diffs[1]) + unlist(diffs[length(diffs)])) / 2
          # generate list of days
          leapcount <- 0
          for(ys in 1:nyears){
            if(ys == 1){
              if(nyears == 1){
                if(yearstodo[ys] %in% leapyears){
                  leapcount <- leapcount + 1
                  day<-c(1, 15.5, 45.5, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5, 366)
                }else{
                  day<-c(1, 15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 365)
                }
              }else{
                if(yearstodo[ys] %in% leapyears){
                  leapcount <- leapcount + 1
                  day<-c(1, 15.5, 45.5, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5)
                }else{
                  day<-c(1, 15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5)
                }
              }
            }else{
              if(ys == nyears){
                if(yearstodo[ys] %in% leapyears){
                  leapcount <- leapcount + 1
                  day<-c(15.5, 45.5, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5, 366)
                }else{
                  day<-c(15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 365)
                }
              }else{
                if(yearstodo[ys] %in% leapyears){
                  leapcount <- leapcount + 1
                  day<-c(15.5, 45.5, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5)
                }else{
                  day<-c(15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5)
                }
              }
            }
            if(ys == 1){
              days2 <- day
              days <- day
            }else{
              if(yearstodo[ys] %in% leapyears){
                days2 <- c(days2, (day + 366 * (ys - 1)) + leapcount)
              }else{
                days2 <- c(days2, (day + 365 * (ys - 1)) + leapcount)
              }
              days <- c(days, day)
            }
          }

          diffs3 <- c(diff1, diffs, diff1)
          days_diffs <- data.frame(matrix(NA, nrow = nyears * 12 + 2, ncol = 3))
          days_diffs[, 1] <- days
          days_diffs[, 3] <- days2
          days_diffs[, 2] <- diffs3
          colnames(days_diffs) <- c("days", "diffs", "new_day")

          # interpolate monthly differences
          f <- approxfun(x = days_diffs$new_day, y = days_diffs$diffs)
          xx <- seq(1, max(days2), 1)
          sp_diff <- f(xx)
          return(sp_diff)
        }

        terra <- as.data.frame(get_terra(x = loc, ystart = ystart_terra, yfinish = yfinish_terra, source = terra_source))
        if(is.infinite(max(terra))){
          message("no TerraClimate data available", '\n')
          stop()
        }
        if(scenario == 4){
          terra_cc <- as.data.frame(get_terra(x = loc, ystart = ystart_terra, yfinish = yfinish_terra, scenario = 4, source = terra_source))
        }
        if(scenario == 2){
          terra_cc <- as.data.frame(get_terra(x = loc, ystart = ystart_terra, yfinish = yfinish_terra, scenario = 2, source = terra_source))
        }
        tmin_diffs <- terra_cc$TMINN - terra$TMINN
        tmax_diffs <- terra_cc$TMAXX - terra$TMAXX
        precip_diffs <- terra_cc$RAINFALL / terra$RAINFALL
        srad_diffs <- terra_cc$SRAD / terra$SRAD
        vpd_diffs <- terra_cc$VPD / terra$VPD
        srad_diffs[is.na(srad_diffs)] <- 1
        vpd_diffs[is.na(vpd_diffs)] <- 1
        precip_diffs[is.na(precip_diffs)] <- 1
        srad_diffs[is.infinite(srad_diffs)] <- 1
        vpd_diffs[is.infinite(vpd_diffs)] <- 1
        precip_diffs[is.infinite(precip_diffs)] <- 1

        TMINN_diff <- getdiff(tmin_diffs)
        TMAXX_diff <- getdiff(tmax_diffs)
        VPD_diff <- getdiff(vpd_diffs)
        SOLAR_diff <- getdiff(srad_diffs)
        RAIN_diff <- getdiff(precip_diffs)

        TMINN_diff <- rep(TMINN_diff, each=24)[1:length(TAIRhr)]
        TMAXX_diff <- rep(TMAXX_diff, each=24)[1:length(TAIRhr)]
        VPD_diff <- rep(VPD_diff, each=24)[1:length(TAIRhr)]
        SOLAR_diff <- rep(SOLAR_diff, each=24)[1:length(TAIRhr)]


        # code to apply changes in min and max air temperature to hourly air temperature data

        # find times of maxima and minima
        mins <- aggregate(TAIRhr, by = list(format(dates, "%Y-%m-%d")), FUN = min)$x
        maxs <- aggregate(TAIRhr, by = list(format(dates, "%Y-%m-%d")), FUN = max)$x
        mins <- rep(mins, each=24)
        maxs <- rep(maxs, each=24)
        mins2 <- TAIRhr / mins[1:length(TAIRhr)]
        maxs2 <- TAIRhr / maxs
        mins2[mins2 != 1] <- 0
        maxs2[maxs2 != 1] <- 0

        # construct weightings so that changes in min and max air temp can be applied
        minweight <- rep(NA, length(mins2))
        maxweight <- minweight
        minweight[mins2 == 1] <- 1
        minweight[maxs2 == 1] <- 0
        maxweight[mins2 == 1] <- 0
        maxweight[maxs2 == 1] <- 1
        minweight <- zoo::na.approx(minweight, na.rm = FALSE)
        maxweight <- zoo::na.approx(maxweight, na.rm = FALSE)
        minweight[is.na(minweight)] <- 0.5
        maxweight[is.na(maxweight)] <- 0.5
        minweight[minweight > 1] <- 1
        minweight[minweight < 0] <- 0
        maxweight[maxweight > 1] <- 1
        maxweight[maxweight < 0] <- 0

        # construct final weighted correction factor and apply
        diff <- TMINN_diff * minweight + TMAXX_diff * maxweight
        TAIRhr <- TAIRhr + diff
      }
      SOLRhr <- hourlyradwind$swrad / 0.0036
      if(scenario != 0){
        SOLRhr <- SOLRhr * SOLAR_diff
      }
      #SOLRhr <- hourlyradwind$swrad / 0.0036
      SOLRhr[SOLRhr < 0] <- 0
      SOLRhr <- zoo::na.approx(SOLRhr)
      cloudhr <- hourlydata$cloudcover
      IRDhr <- hourlydata$downlong / 0.0036
      if(IR != 2){
        IRDhr <- IRDhr * 0 + -1 # make negative so computed internally in the microclimate model
      }
      CLDhr <- hourlydata$cloudcover
      if(scenario != 0){
        CLDhr <- CLDhr * (1 + 1 - SOLAR_diff)
      }
      CLDhr[CLDhr < 0] <- 0
      CLDhr[CLDhr > 100] <- 100
      RHhr <- suppressWarnings(humidityconvert(h = hourlydata$humidity, intype = 'specific', p = hourlydata$pressure, tc = TAIRhr)$relative)
      RHhr[RHhr > 100] <- 100
      RHhr[RHhr < 0] <- 0
      if(scenario != 0){
        VPDhr <- WETAIR(db = TAIRhr_orig, rh = 100)$e - WETAIR(db = TAIRhr_orig, rh = RHhr)$e
        VPDhr <- VPDhr - mean(VPD_diff) # note using mean here because otherwise can lead to unrealistically low RH which then stronly affects TSKY
        e <- WETAIR(db = TAIRhr, rh = 100)$e - VPDhr
        es <- WETAIR(db = TAIRhr, rh = 100)$esat
        RHhr <- (e / es) * 100
        RHhr[RHhr > 100] <- 100
        RHhr[RHhr < 0] <- 0
      }
      WNhr <- hourlyradwind$windspeed
      WNhr[is.na(WNhr)] <- 0.1
      if(rainhourly == 0){
        RAINhr = rep(0, 24 * ndays)
      }else{
        RAINhr = rainhour
      }
      if(length(rainhourly) == length(tme)){ # an hourly rainfall vector has been provided
        aggvec <- rep(1:(length(tme)/24), 24)
        aggvec <- aggvec[order(aggvec)]
        RAINFALL <- aggregate(rainhourly, by = aggvec, FUN = 'sum') # aggregate to daily totals
      }else{
        RAINFALL <- dailyprecip
      }
      PRESShr <- hourlydata$pressure
      if(scenario != 0){
        RAINFALL <- RAINFALL * RAIN_diff[1:length(RAINFALL)] # TODO hack to get around issue with leap years - check this
      }
      RAINFALL[RAINFALL < 0.1] <- 0

      if(rainhourly == 0){ # putting daily rainfall on midnight (account for local solar time)
        midnight <- which(microclima.out$hourlydata$szenith[1:24] == max(microclima.out$hourlydata$szenith[1:24]))
        if(midnight <= 24){
          midnights <- seq(midnight, length(RAINhr / 24) - 1, 24)
        }else{
          midnights <- seq(midnight, length(RAINhr / 24), 24)
        }
        RAINhr[midnights+1] <- RAINFALL/3
        RAINhr[midnights] <- RAINFALL/3
        RAINhr[midnights-1] <- RAINFALL/3
        rainhourly <- 1
      }
      ZENhr2 <- ZENhr
      ZENhr2[ZENhr2!=90] <- 0
      dmaxmin <- function(x, fun) {
        dx <- t(matrix(x, nrow = 24))
        apply(dx, 1, fun)
      }
      TMAXX <- dmaxmin(TAIRhr, max)
      TMINN <- dmaxmin(TAIRhr, min)
      CCMAXX <- dmaxmin(CLDhr, max)
      CCMINN <- dmaxmin(CLDhr, min)
      RHMAXX <- dmaxmin(RHhr, max)
      RHMINN <- dmaxmin(RHhr, min)
      WNMAXX <- dmaxmin(WNhr, max)
      WNMINN <- dmaxmin(WNhr, min)
      PRESS <- dmaxmin(PRESShr, min)

      if(clearsky == 1){
        cat("running micro_global to get clear sky solar \n")
        micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, timeinterval = 365, solonly = 1)
        clearskyrad <- micro_clearsky$metout[, c(1, 13)]
        leapyears <- seq(1900, 2060, 4)
        ystart <- as.numeric(substr(dstart, 7, 10))
        yfinish <- as.numeric(substr(dfinish, 7, 10))
        datestimes <- seq(as.POSIXct(paste0(ystart, '-01-01 00:00:00')), as.POSIXct(paste0(yfinish, '-12-31 23:00:00')), by = 'hours')
        for(year in ystart:yfinish){
          clear <- clearskyrad[, 2]
          if(year %in% leapyears){
            clear <- c(clear[1:1416], clear[1417:(1416 + 24)], clear[(1417 + 1):length(clear)])
          }
          if(year == ystart){
            SOLRhr <- clear
          }else{
            SOLRhr <- c(SOLRhr, clear)
          }
        }
        tzoff <- 12-which(ZENhr[1:24] ==min(ZENhr[1:24]))+1
        if(tzoff != 12){
          if(tzoff > 12){
            SOLRhr <- c(SOLRhr[(length(SOLRhr)-tzoff+1):(length(SOLRhr))], SOLRhr[1:(length(SOLRhr)-tzoff)])
          }else{
            SOLRhr <- c(SOLRhr[(tzoff+1):length(SOLRhr)], SOLRhr[1:tzoff])
          }
        }
      }
      if(save == 1){
        cat("saving met data for later \n")
        save(CCMAXX, file = 'CCMAXX.Rda')
        save(CCMINN, file = 'CCMINN.Rda')
        save(WNMAXX, file = 'WNMAXX.Rda')
        save(WNMINN, file = 'WNMINN.Rda')
        save(TMAXX, file = 'TMAXX.Rda')
        save(TMINN, file = 'TMINN.Rda')
        save(RHMAXX, file = 'RHMAXX.Rda')
        save(RHMINN, file = 'RHMINN.Rda')
        save(RAINFALL, file = 'RAINFALL.Rda')
        save(PRESS, file = 'PRESS.Rda')
        save(CLDhr, file = 'CLDhr.Rda')
        save(WNhr, file = 'WNhr.Rda')
        save(TAIRhr, file = 'TAIRhr.Rda')
        save(RHhr, file = 'RHhr.Rda')
        save(RAINhr, file = 'RAINhr.Rda')
        save(SOLRhr, file = 'SOLRhr.Rda')
        save(ZENhr, file = 'ZENhr.Rda')
        save(IRDhr, file = 'IRDhr.Rda')
        save(microclima.out, file = 'microclima.out.Rda')
      }
    }else{
      cat("loading met data from previous run \n")
      load('CCMAXX.Rda')
      load('CCMINN.Rda')
      load('WNMAXX.Rda')
      load('WNMINN.Rda')
      load('TMAXX.Rda')
      load('TMINN.Rda')
      load('RHMAXX.Rda')
      load('RHMINN.Rda')
      load('RAINFALL.Rda')
      load('PRESS.Rda')
      load('CLDhr.Rda')
      load('WNhr.Rda')
      load('TAIRhr.Rda')
      load('RHhr.Rda')
      load('RAINhr.Rda')
      load('SOLRhr.Rda')
      load('ZENhr.Rda')
      load('IRDhr.Rda')
      load('microclima.out.Rda')
      rainhourly <- 1
    }
    slope <- 0 # already corrected for by microclima
    azmuth <- 0 # already corrected for by microclima

    TAI <- compute_tai(longlat, run.gads, TAI_AUSTRALIA)

    if(max(abs(warm)) != 0){
      if(length(warm) == length(TMAXX)){
       # impose uniform temperature change
       TMAXX <- TMAXX + warm
       TMINN <- TMINN + warm
       warm.hr <- rep(warm, each = 24)
      }else{
       warm.hr <- warm
      }
      TAIRhr <- TAIRhr + warm.hr
      sigma <- 5.67e-8 #Stefan-Boltzman, W/(m.K)
      if(IR == 2){
        IRDhr <- sigma * ((IRDhr / sigma) ^ (1 / 4) + warm.hr) ^ 4 # adjust downward radiation for altered 'sky temperature'
      }
    }
    RAINFALL <- RAINFALL + rainoff
    RAINFALL[RAINFALL < 0] <- 0
    if(length(rainoff) == length(RAINFALL)){
      rainoff.hr <- rep(rainoff, each = 24)
    }else{
      rainoff.hr <- rainoff
    }
    RAINhr <- RAINhr + rainoff.hr
    RAINhr[RAINhr < 0] <- 0
    ALLMINTEMPS <- TMINN
    ALLMAXTEMPS <- TMAXX
    ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS)

    WNMAXX <- WNMAXX * windfac
    WNMINN <- WNMINN * windfac
    WNhr <- WNhr * windfac

    REFLS <- rep(REFL, ndays)
    PCTWET <- rep(PCTWET, ndays)
    soilwet <- RAINFALL
    soilwet[soilwet <= rainwet] <- 0
    soilwet[soilwet > 0] <- 90
    if(ndays < 1){
      PCTWET <- pmax(soilwet, PCTWET)
    }

    Intrvls<-rep(0, ndays)
    Intrvls[1] <- 1 # user-supplied last day-of-year in each time interval sequence
    Numtyps <- 10 # number of substrate types
    Nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
    Nodes[1:10,] <- c(1:10) # deepest nodes for each substrate type
    ALREF <- abs(trunc(x[1]))

    HEMIS <- ifelse(x[2] < 0, 2, 1)
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2]) - ALAT) * 60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1]) - ALONG) * 60

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
    tannul <- mean(unlist(ALLTEMPS))

    if(is.na(deepsoil)){
      if(nyears == 1){
        avetemp <- (sum(TMAXX) + sum(TMINN)) / (length(TMAXX)*2)
        deepsoil <-rep(avetemp, ndays)
      }else{
        avetemp <- rowMeans(cbind(TMAXX, TMINN), na.rm=TRUE)
        if(length(TMAXX) < 365){
          deepsoil <- rep((sum(TMAXX) + sum(TMINN)) / (length(TMAXX) * 2), length(TMAXX))
        }else{
          deepsoil <- terra::roll(avetemp, n = 365, fun = mean, type = 'to')
          yearone <- rep((sum(TMAXX[1:365]) + sum(TMINN[1:365])) / (365 * 2), 365)
          deepsoil[1:365] <- yearone
        }
      }
    }else{
      tannul <- mean(deepsoil)
    }

    SLES <- matrix(nrow = ndays, data = 0)
    SLES <- SLES + SLE

    moists2 <- matrix(nrow = 10, ncol = ndays, data = 0)
    moists2[1, ndays] <- 0.2
    moists <- moists2

    if(runmoist == 1){
      moists2 <- matrix(nrow = 10, ncol = ndays, data = 0) # set up an empty vector for soil moisture values through time
      moists2[1:10, ] <- SoilMoist_Init
      moists <- moists2
    }
    soilprops <- matrix(data = 0, nrow = 10, ncol = 5)
    soilprops[,1] <- BulkDensity
    soilprops[,2] <- 1 - BulkDensity / Density # not used if soil moisture computed
    soilprops[soilprops[,2] < 0.26, 2] <- 0.26
    soilprops[,3] <- Thcond
    soilprops[,4] <- SpecHeat
    soilprops[,5] <- Density
    if(cap == 1){
      soilprops[1:2, 3] <- 0.2
      soilprops[1:2, 4] <- 1920
    }
    if(cap==2){
      soilprops[1:2, 3] <- 0.1
      soilprops[3:4, 3] <- 0.25
      soilprops[1:4, 4] <- 1920
      soilprops[1:4, 5] <- 1.3
      soilprops[1:4, 1] <- 0.7
    }

    ALTT <- as.numeric(ALTT)
    ALREF <- as.numeric(ALREF)
    ALMINT <- as.numeric(ALMINT)
    ALONG <- as.numeric(ALONG)
    AMINUT <- as.numeric(AMINUT)
    ALAT <-as.numeric(ALAT)
    hourly <- 1

    TIMAXS <- c(1, 1, 0, 0)
    TIMINS <- c(0, 0, 1, 1)
    # --- 3. Build microclimate model input list ---
    microinput <- build_microinput(ndays, RUF, ERR, Usrhyt, Refhyt, Numtyps,
      Z01, Z02, ZH1, ZH2, idayst, ida, HEMIS, ALAT, AMINUT, ALONG, ALMINT,
      ALREF, slope, azmuth, ALTT, CMH2O, microdaily, tannul, EC, VIEWF,
      snowtemp, snowdens, snowmelt, undercatch, rainmult, runshade, runmoist,
      maxpool, evenrain, snowmodel, rainmelt, writecsv, densfun, hourly,
      rainhourly, lamb, IUV, RW, PC, RL, SP, R1, IM, MAXCOUNT, IR, message,
      fail, snowcond, intercept, grasshade, solonly, ZH, D0, TIMAXS, TIMINS,
      spinup, dewrain = dewrain, moiststep = moiststep, maxsurf, ndmax)

    if(length(LAI) < ndays){
      LAI<-rep(LAI[1], ndays)
    }
    if(shore == 0){
      tides <- matrix(data = 0, nrow = 24 * ndays, ncol = 3) # make an empty matrix
    }

    micro <- build_micro_list(microinput, doy, SLES, DEP, Nodes, MAXSHADES,
      MINSHADES, TMAXX, TMINN, RHMAXX, RHMINN, CCMAXX, CCMINN, WNMAXX, WNMINN,
      TAIRhr, RHhr, WNhr, CLDhr, SOLRhr, RAINhr, ZENhr, IRDhr, REFLS, PCTWET,
      soilinit, hori, TAI, soilprops, moists, RAINFALL, deepsoil = deepsoil,
      PE, KS, BB, BD, DD, L, LAI, tides)

    # --- 4. Write inputs to CSV (optional) ---
    if (write_input == 1) {
      write_micro_csv(microinput, doy, SLES, DEP, Nodes, MAXSHADES, MINSHADES,
        TIMAXS = c(1, 1, 0, 0), TIMINS = c(0, 0, 1, 1),
        TMAXX, TMINN, RHMAXX, RHMINN, CCMAXX, CCMINN, WNMAXX, WNMINN,
        REFLS, PCTWET, soilinit, hori, TAI, soilprops, moists, RAINFALL,
        deepsoil, PE, BD, DD, BB, KS, L, LAI, tides,
        TAIRhr, RHhr, WNhr, CLDhr, SOLRhr, RAINhr, ZENhr, IRDhr)
    }
    if(is.numeric(loc[1])){
      location <- paste("long", loc[1], "lat", loc[2])
    }else{
      location <- loc
    }
    if(runmicro){
      cat(paste('running microclimate model for ', ndays, ' days from ', dates[1], ' to ', dates[length(dates)], ' at site ', location, '\n'))
    cat('Note: the output column `SOLR` in metout and shadmet is for unshaded solar radiation adjusted for slope, aspect and horizon angle \n')
      ptm <- proc.time() # Start timing
      microut<-microclimate(micro)
      print(proc.time() - ptm) # Stop the clock


    # --- 5. Process and return results ---
    out <- process_micro_output(microut, runmoist, snowmodel, lamb)
    if (max(out$metout[, 1] == 0)) {
      message("ERROR: the model crashed - try a different error tolerance (ERR) or a different spacing in DEP")
    }
    return(build_micro_return(out, RAINFALL, ndays, ALTT, REFL,
      longlat = longlat, nyears, timeinterval = ndays, MINSHADES, MAXSHADES,
      DEP, dates, dates2, PE, BD, DD, BB, KS, dem = dem, diffuse_frac = NA,
      snowmodel, lamb,
      extra = list(SLOPE = SLOPE, ASPECT = ASPECT, HORIZON = HORIZON,
                   microclima.out = microclima.out)))
    }else{
      return(list(RAINFALL = RAINFALL, TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN, WNMAXX = WNMAXX, WNMINN = WNMINN, CCMAXX = CCMAXX, CCMINN = CCMINN, CLDhr = CLDhr, WNhr = WNhr, TAIRhr = TAIRhr, RHhr = RHhr, RAINhr = RAINhr, SOLRhr = SOLRhr, ZENhr = ZENhr, IRDhr = IRDhr, dates = dates, dates2 = dates2, PE=PE,BD=BD,DD=DD,BB=BB,KS=KS))
    }
  } # end error trapping
} # end of micro_era5 function
