#' Global implementation of the microclimate model
#'
#' An implementation of the NicheMapR microclimate model that uses the global climate database
#' derived from "New, M., Lister, D., Hulme, M. and Makin, I., 2002: A high-resolution data
#' set of surface climate over global land areas. Climate Research 21:1-25"
#' It also optionally uses a global monthly soil moisture estimate from NOAA CPC Soil Moisture http://140.172.38.100/psd/thredds/catalog/Datasets/cpcsoil/catalog.html
#' Aerosol attenuation can also be computed based on the Global Aerosol Data Set (GADS)
#' Koepke, P., M. Hess, I. Schult, and E. P. Shettle. 1997. Global Aerosol Data Set. Max-Planck-Institut for Meteorologie, Hamburg
#' by choosing the option 'run.gads<-1'
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param timeinterval The number of time intervals to generate predictions for over a year (must be 12 <= x <=365)
#' @param nyears The number of years to run
#' @param REFL Soil solar reflectance, decimal \%
#' @param elev Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param soiltype Soil type: Rock = 0, sand = 1, loamy sand = 2, sandy loam = 3, loam = 4, silt loam = 5, sandy clay loam = 6, clay loam = 7, silt clay loam = 8, sandy clay = 9, silty clay = 10, clay = 11, user-defined = 12, based on Campbell and Norman 1990 Table 9.1.
#' @param minshade Minimum shade level to use (\%)
#' @param maxshade Maximum shade level to use (\%)
#' @param Usrhyt Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
#' @param ... Additional arguments, see Details
#' @usage micro_global(loc = c(-89.40123, 43.07305), timeinterval = 12, nyears = 1, soiltype = 4,
#' REFL = 0.15, slope = 0, aspect = 0,
#' DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minshade = 0, maxshade = 90,
#' Usrhyt = 0.01, ...)
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
#' @details
#' \itemize{
#' \strong{Parameters controlling how the model runs:}\cr\cr
#'
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes, 0=no}\cr\cr
#' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{solonly}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{ERR}{ = 1.5, Integrator error tolerance for soil temperature calculations}\cr\cr
#' \code{Refhyt}{ = 1.2, Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured}\cr\cr
#' \code{RUF}{ = 0.004, Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001 (snow) - 0.02 m.}\cr\cr
#' \code{ZH}{ = 0, heat transfer roughness height (m) for Campbell and Norman air temperature/wind speed profile (invoked if greater than 1, 0.02 * canopy height in m if unknown)}\cr\cr
#' \code{D0}{ = 0, zero plane displacement correction factor (m) for Campbell and Norman air temperature/wind speed profile (0.6 * canopy height in m if unknown)}\cr\cr
#' \code{Z01}{ = 0, Top (1st) segment roughness height(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then RUF and Refhyt used)}\cr\cr
#' \code{Z02}{ = 0, 2nd segment roughness height(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then RUF and Refhyt used).}\cr\cr
#' \code{ZH1}{ = 0, Top of (1st) segment, height above surface(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then RUF and Refhyt used).}\cr\cr
#' \code{ZH2}{ = 0, 2nd segment, height above surface(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then RUF and Refhyt used).}\cr\cr
#' \code{EC}{ = 0.0167238, Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)}\cr\cr
#' \code{SLE}{ = 0.95, Substrate longwave IR emissivity (decimal \%), typically close to 1}\cr\cr
#' \code{Thcond}{ = 2.5, Soil minerals thermal conductivity (W/mK)}\cr\cr
#' \code{Density}{ = 2.56, Soil minerals density (Mg/m3)}\cr\cr
#' \code{SpecHeat}{ = 870, Soil minerals specific heat (J/kg-K)}\cr\cr
#' \code{BulkDensity}{ = 1.3, Soil bulk density (mg/m3)}\cr\cr
#' \code{PCTWET}{ = 0, \% of ground surface area acting as a free water surface}\cr\cr
#' \code{cap}{ = 1, organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)}\cr\cr
#' \code{CMH2O}{ = 1, Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air conditions; 2.0 = humid, tropical conditions (note this is for the whole atmospheric profile, not just near the ground)}\cr\cr
#' \code{hori}{ = rep(0,24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#' \code{lapse_min}{ = 0.0039 Lapse rate for minimum air temperature (degrees C/m)}
#' \code{lapse_max}{ = 0.0077 Lapse rate for maximum air temperature (degrees C/m)}
#' \code{TIMAXS}{ = c(1, 1, 0, 0), Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise}\cr\cr
#' \code{TIMINS}{ = c(0, 0, 1, 1), Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon}\cr\cr
#' \code{timezone}{ = 0, Use GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)? 1=yes, 0=no}\cr\cr
#' \code{TAI}{ = 0, Vector of 111 values, one per wavelenght bin, for solar attenuation - used to overide GADS}\cr\cr
#'
#' \strong{ Soil moisture mode parameters:}
#'
#' \code{runmoist}{ = 0, Run soil moisture model? 1=yes, 0=no  1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{PE}{ = rep(1.1,19), Air entry potential (J/kg) (19 values descending through soil for specified soil nodes in parameter}
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{KS}{ = rep(0.0037,19), Saturated conductivity, (kg s/m3) (19 values descending through soil for specified soil nodes in parameter}
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{BB}{ = rep(4.5,19), Campbell's soil 'b' parameter (-) (19 values descending through soil for specified soil nodes in parameter}
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{BD}{ = rep(1.3,19), Soil bulk density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)}\cr\cr
#' \code{DD}{ = rep(2.56,19), Soil density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter DEP and points half way between)}\cr\cr
#' \code{maxpool}{ = 10000, Max depth for water pooling on the surface (mm), to account for runoff}\cr\cr
#' \code{rainmult}{ = 1, Rain multiplier for surface soil moisture (-), used to induce runon}\cr\cr
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
#' \code{LAI}{ = 0.1, leaf area index, used to partition traspiration/evaporation from PET}\cr\cr
#'
#' \strong{ Snow mode parameters:}
#'
#' \code{snowmodel}{ = 0, run the snow model 1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (mg/m3), overridden by densfun}\cr\cr
#' \code{densfun}{ = c(0.5979, 0.2178, 0.001, 0.0038), slope and intercept of model of snow density as a linear function of snowpack age if first two values are nonzero, and following the exponential function of Sturm et al. 2010 J. of Hydromet. 11:1380-1394 if all values are non-zero; if it is c(0,0,0,0) then fixed density used}\cr\cr
#' \code{snowmelt}{ = 1, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
#' \code{undercatch}{ = 1, undercatch multipier for converting rainfall to snow}\cr\cr
#' \code{rainmelt}{ = 0.0125, paramter in equation that melts snow with rainfall as a function of air temp}\cr\cr
#' \code{rainfrac}{ = 0.5, fraction of rain that falls on the first day of the month (decimal \% with 0 meaning rain falls evenly) - this parameter allows something other than an even intensity of rainfall when interpolating the montly rainfall data)}\cr\cr
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)}\cr\cr
#' \code{intercept}{ = maxshade / 100 * 0.3, snow interception fraction for when there's shade (0-1)}\cr\cr
#' \code{grasshade}{ = 0, if 1, means shade is removed when snow is present, because shade is cast by grass/low shrubs}\cr\cr
#'
#' \strong{ Intertidal mode parameters:}
#'
#' \code{shore}{ Include tide effects? If 1, the matrix}
#' \code{tides}
#' { is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0., nrow = 24*timeinterval*nyears, ncol = 3), matrix for each how of the simulation of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#' }
#'
#' \strong{Outputs:}
#'
#' \code{ndays}{ - number of days for which predictions are made}\cr\cr
#' \code{longlat}{ - longitude and latitude for which simulation was run (decimal degrees)}\cr\cr
#' \code{nyears}{ - number of years for which predictions are made}\cr\cr
#' \code{RAINFALL}{ - vector of daily rainfall (mm)}\cr\cr
#' \code{elev}{ - elevation at point of simulation (m)}\cr\cr
#' \code{minshade}{ - minimum shade for simulation (\%)}\cr\cr
#' \code{maxshade}{ - maximum shade for simulation (single value - if time varying, in 'MAXSHADES') (\%)}\cr\cr
#' \code{MAXSHADES}{ - vector of maximum shades used (\%)}\cr\cr
#' \code{dem}{ - digital elevation model obtained via 'get_dev' using package 'elevatr' (m)}\cr\cr
#' \code{DEP}{ - vector of depths used (cm)}\cr\cr
#'
#' metout/shadmet variables:
#' \itemize{
#' \item 1 DOY - day-of-year
#' \item 2 TIME - time of day (mins)
#' \item 3 TALOC - air temperature (°C) at local height (specified by 'Usrhyt' variable)
#' \item 4 TAREF - air temperature (°C) at reference height (specified by 'Refhyt', 1.2m default)
#' \item 5 RHLOC - relative humidity (\%) at local height (specified by 'Usrhyt' variable)
#' \item 6 RH  - relative humidity (\%) at reference height (specified by 'Refhyt', 1.2m default)
#' \item 7 VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)
#' \item 8 VREF - wind speed (m/s) at reference height (specified by 'Refhyt', 1.2m default)
#' \item 9 SNOWMELT - snowmelt (mm)
#' \item 10 POOLDEP - water pooling on surface (mm)
#' \item 11 PCTWET - soil surface wetness (\%)
#' \item 12 ZEN - zenith angle of sun (degrees - 90 = below the horizon)
#' \item 13 SOLR - solar radiation (W/m2)
#' \item 14 TSKYC - sky radiant temperature (°C)
#' \item 15 DEW - dew presence (0 or 1)
#' \item 16 FROST - frost presence (0 or 1)
#' \item 17 SNOWFALL - snow predicted to have fallen (cm)
#' \item 18 SNOWDEP - predicted snow depth (cm)
#' \item 19 SNOWDENS - snow density (g/cm3)
#' }
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
#' \item 3-12 WC0cm ... - soil moisuture (m3/m3) at each of the 10 specified depths
#' }
#' soilpot and shadpot variables:
#' \itemize{
#' \item 1 DOY - day-of-year
#' \item 2 TIME - time of day (mins)
#' \item 3-12 PT0cm ... - soil water potential (J/kg = kpa = bar/100) at each of the 10 specified depths
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
#' \item  3 TRANS - plant transpiration rate (kg/m2/s)
#' \item  4 LEAFPOT - leaf water potential (J/kg)
#' \item  5-14 RPOT0cm ... - root water potential (J/kg), at each of the 10 specified depths
#' }
#'
#' if snow model is run i.e. parameter lamb = 1\cr
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
#'micro <- micro_global() # run the model with default location (Madison, Wisconsin) and settings
#'
#'metout <- as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#'shadmet <- as.data.frame(micro$shadmet) # above ground microclimatic conditions, max shade
#'soil <- as.data.frame(micro$soil) # soil temperatures, minimum shade
#'shadsoil <- as.data.frame(micro$shadsoil) # soil temperatures, maximum shade
#'
#'minshade <- micro$minshade
#'maxshade <- micro$maxshade
#'
#'# plotting above-ground conditions in minimum shade
#'with(metout, {plot(TALOC ~ micro$dates, xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l", main = paste("air temperature, ", minshade, "% shade",sep = ""))})
#'with(metout, {points(TAREF ~ micro$dates, xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l",lty = 2, col = 'blue')})
#'with(metout, {plot(RHLOC ~ micro$dates, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",ylim = c(0, 100),main = paste("humidity, ",minshade, "% shade",sep=""))})
#'with(metout, {points(RH ~ micro$dates, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",col = 'blue',lty = 2, ylim = c(0, 100))})
#'with(metout, {plot(TSKYC ~ micro$dates, xlab = "Date and Time", ylab = "Sky Temperature (°C)"
#',  type = "l", main = paste("sky temperature, ", minshade, "% shade", sep=""))})
#'with(metout, {plot(VREF ~ micro$dates, xlab = "Date and Time",  ylab = "Wind Speed (m/s)"
#',  type = "l", main = "wind speed", col = 'blue',ylim = c(0, 15))})
#'with(metout, {points(VLOC ~ micro$dates, xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#',  type = "l", lty = 2)})
#'with(metout, {plot(ZEN ~ micro$dates,xlab = "Date and Time", ylab = "Zenith Angle of Sun (deg)"
#',  type = "l", main = "solar angle, sun")})
#'with(metout, {plot(SOLR ~ micro$dates,xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
#',  type = "l", main = "solar radiation")})
#'
#'# plotting soil temperature for minimum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(soil[,i + 2] ~ micro$dates, xlab = "Date and Time", ylab = "Soil Temperature (°C)"
#'    ,col = i, type = "l", main = paste("soil temperature ", minshade, "% shade", sep=""))
#'  }else{
#'    points(soil[,i + 2] ~ micro$dates, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l")
#'  }
#'}
#'
#'# plotting above-ground conditions in maximum shade
#'with(shadmet,{plot(TALOC ~ micro$dates,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l", main = "air temperature, sun")})
#'with(shadmet,{points(TAREF ~ micro$dates,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l", lty = 2, col = 'blue')})
#'with(shadmet,{plot(RHLOC ~ micro$dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l", ylim = c(0, 100),main = "humidity, shade")})
#'with(shadmet,{points(RH ~ micro$dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l", col = 'blue',lty = 2, ylim = c(0, 100))})
#'with(shadmet,{plot(TSKYC ~ micro$dates,xlab = "Date and Time", ylab = "Sky Temperature (°C)",
#'  type = "l", main = "sky temperature, shade")})
#'
#'# plotting soil temperature for maximum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(shadsoil[,i + 2] ~ micro$dates, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l", main = paste("soil temperature ", maxshade, "% shade", sep=""))
#'  }else{
#'    points(shadsoil[,i + 2] ~ micro$dates, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l")
#'  }
#'}
#' @export
micro_global <- function(
  loc = c(-89.40123, 43.07305),
  timeinterval = 12,
  nyears = 1,
  soiltype = 4,
  REFL = 0.15,
  elev = NA,
  slope = 0,
  aspect = 0,
  lapse_max = 0.0077,
  lapse_min = 0.0039,
  DEP=c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
  minshade = 0,
  maxshade = 90,
  Refhyt = 1.2,
  Usrhyt = 0.01,
  Z01 = 0,
  Z02 = 0,
  ZH1 = 0,
  ZH2 = 0,
  runshade = 1,
  solonly = 0,
  clearsky = 0,
  run.gads = 1,
  write_input = 0,
  writecsv = 0,
  ERR = 1.5,
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
  cap = 1,
  CMH2O = 1,
  hori=rep(0,24),
  TIMAXS = c(1, 1, 0, 0),
  TIMINS = c(0, 0, 1, 1),
  timezone = 0,
  runmoist = 0,
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
  RL = 2e+6,
  PC = -1500,
  SP = 10,
  IM = 1e-06,
  MAXCOUNT = 500,
  LAI=0.1,
  snowmodel = 0,
  snowtemp = 1.5,
  snowdens = 0.375,
  densfun = c(0.5979, 0.2178, 0.001, 0.0038),
  snowmelt = 1,
  undercatch = 1,
  rainmelt = 0.0125,
  rainfrac = 0.5,
  shore = 0,
  tides = 0,
  lamb = 0,
  IUV = 0,
  soilgrids = 0,
  IR = 0,
  message = 0,
  fail = nyears * 24 * 365,
  TAI = 0,
  snowcond = 0,
  intercept = maxshade / 100 * 0.4,
  grasshade = 0
) {

  SoilMoist=SoilMoist_Init
  errors<-0

  # error trapping - originally inside the Fortran code, but now checking before executing Fortran
  if(DEP[2]-DEP[1]>3 | DEP[3]-DEP[2]>3){
    message("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if(DEP[2]-DEP[1]<2){
    cat("warning, nodes might be too close near the surface, try a different spacing if the program is crashing \n")
  }
  if(timeinterval<12 | timeinterval > 365){
    message("ERROR: the variable 'timeinterval' is out of bounds.
        Please enter a correct value (12 - 365).", '\n')
    errors<-1
  }
  if(is.numeric(loc[1])){
    if(loc[1]>180 | loc[2] > 90){
      message("ERROR: Latitude or longitude (longlat) is out of bounds.
        Please enter a correct value.", '\n')
      errors<-1
    }
  }
  if(timezone%in%c(0,1)==FALSE){
    message("ERROR: the variable 'timezone' be either 0 or 1.
      Please correct.", '\n')
    errors<-1
  }
  if(run.gads%in%c(0,1)==FALSE){
    message("ERROR: the variable 'run.gads' be either 0 or 1.
      Please correct.", '\n')
    errors<-1
  }
  if(write_input%in%c(0,1)==FALSE){
    message("ERROR: the variable 'write_input' be either 0 or 1.
      Please correct.", '\n')
    errors<-1
  }
  if(EC<0.0034 | EC > 0.058){
    message("ERROR: the eccentricity variable (EC) is out of bounds.
        Please enter a correct value (0.0034 - 0.058).", '\n')
    errors<-1
  }
  if(RUF<0.0001){
    message("ERROR: The roughness height (RUF) is too small ( < 0.0001).
        Please enter a larger value.", '\n')
    errors<-1
  }
  if(RUF>2){
    message("ERROR: The roughness height (RUF) is too large ( > 2).
        Please enter a smaller value.", '\n')
    errors<-1
  }
  if(DEP[1]!=0){
    message("ERROR: First soil node (DEP[1]) must = 0 cm.
        Please correct", '\n')
    errors<-1
  }
  if(length(DEP)!=10){
    message("ERROR: You must enter 10 different soil depths.", '\n')
    errors<-1
  }
  for(i in 1:9){
    if(DEP[i+1]<=DEP[i]){
      message("ERROR: Soil depth (DEP array) is not in ascending size", '\n')
      errors<-1
    }
  }
  if(DEP[10]>500){
    message("ERROR: Deepest soil depth (DEP array) is too large (<=500 cm)", '\n')
    errors<-1
  }
  if(Thcond<0){
    message("ERROR: Thermal variable conductivity (THCOND) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(Density<0){
    message("ERROR: Density variable (Density) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(SpecHeat<0){
    message("ERROR: Specific heat variable (SpecHeat) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(min(BulkDensity)<0){
    message("ERROR: Bulk density value (BulkDensity) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(REFL<0 | REFL>1){
    message("ERROR: Soil reflectivity value (REFL) is out of bounds.
        Please input a value between 0 and 1.", '\n')
    errors<-1
  }
  if(slope<0 | slope>90){
    message("ERROR: Slope value (slope) is out of bounds.
        Please input a value between 0 and 90.", '\n')
    errors<-1
  }
  if(aspect<0 | aspect>365){
    message("ERROR: Aspect value (aspect) is out of bounds.
        Please input a value between 0 and 365.", '\n')
    errors<-1
  }
  if(max(hori)>90 | min(hori)<0){
    message("ERROR: At least one of your horizon angles (hori) is out of bounds.
        Please input a value between 0 and 90", '\n')
    errors<-1
  }
  if(length(hori)!=24){
    message("ERROR: You must enter 24 horizon angle values.", '\n')
    errors<-1
  }
  if(SLE<0.5 | SLE > 1){
    message("ERROR: Emissivity (SLE) is out of bounds.
        Please enter a correct value (0.05 - 1.00).", '\n')
    errors<-1
  }
  if(ERR<0){
    message("ERROR: Error bound (ERR) is too small.
        Please enter a correct value (> 0.00).", '\n')
    errors<-1
  }
  if(Usrhyt<RUF){
    message("ERROR: Reference height (Usrhyt) smaller than roughness height (RUF).
        Please use a larger height above the surface.", '\n')
    errors<-1
  }
  if(Usrhyt<0.005 | Usrhyt>Refhyt){
    message("ERROR: Reference height (Usrhyt) is out of bounds.
        Please enter a correct value (0.005 - Refhyt).", '\n')
    errors<-1
  }
  if(CMH2O<0.5 | CMH2O>2){
    message("ERROR: Preciptable water in air column (CMH2O) is out of bounds.
        Please enter a correct value (0.1 - 2).", '\n')
    errors<-1
  }
  if(max(TIMAXS)>24 | min(TIMAXS)<0){
    message("ERROR: At least one of your times of weather maxima (TIMAXS) is out of bounds.
        Please input a value between 0 and 24", '\n')
    errors<-1
  }
  if(max(TIMINS)>24 | min(TIMINS)<0){
    message("ERROR: At least one of your times of weather minima (TIMINS) is out of bounds.
        Please input a value between 0 and 24", '\n')
    errors<-1
  }
  if(minshade>maxshade | minshade==maxshade){
    message("ERROR: Your value for minimum shade (minshade) is greater than or equal to the maximum shade (maxshade).
        Please correct this.", '\n')
    errors<-1
  }
  if(minshade>100 | minshade<0){
    message("ERROR: Your value for minimum shade (minshade) is out of bounds.
        Please input a value between 0 and 100.", '\n')
    errors<-1
  }
  if(maxshade>100 | maxshade<0){
    message("ERROR: Your value for maximum shade (maxshade) is out of bounds.
        Please input a value between 0 and 100.", '\n')
    errors<-1
  }
  if(soiltype<0 | soiltype>11){
    message("ERROR: the soil type must range between 1 and 11.
      Please correct.", '\n')
    errors<-1
  }
  # end error trapping

  if(errors==0){ # continue

    ################## time related variables #################################

    doys12<-c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month
    doysn<-doys12 # variable of doys for when doing multiple years
    if(nyears>1 & timeinterval==365){ # create sequence of days for splining across multiple years
      for(i in 1:(nyears-1)){
        doysn<-c(doysn,(doys12+365*i))
      }
    }

    if(timeinterval<365){
      microdaily<-0 # run microclimate model as normal, where each day is iterated 3 times starting with the initial condition of uniform soil temp at mean monthly temperature
    }else{
      microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in
    }

    # now check if doing something other than middle day of each month, and create appropriate vector of Day of Year
    daystart<-as.integer(ceiling(365/timeinterval/2))
    if(timeinterval!=12){
      doys<-seq(daystart,365,as.integer(floor(365/timeinterval)))
    }else{
      doys<-doysn
    }
    doynum <- timeinterval*nyears # total days to do
    doy <- subset(doys, doys!=0) # final vector of Day of Year
    doy<-rep(doy,nyears)
    idayst <- 1 # start day
    ida<-timeinterval*nyears # end day

    ################## location and terrain #################################

    if(is.numeric(loc)==FALSE){ # might not be quite right format, try correcting
      loc=cbind(as.numeric(loc[1]),as.numeric(loc[2]))
    }
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    # get the local timezone reference longitude
    if(timezone==1){ # this now requires registration
      if(!require(geonames, quietly = TRUE)){
        stop('package "geonames" is required to do a specific time zone (timezone=1). Please install it.')
      }
      ALREF<-(geonames::GNtimezone(longlat[2],longlat[1])[4])*-15
    }else{  # just use local solar noon
      ALREF <- abs(trunc(x[1]))
    }
    HEMIS <- ifelse(x[2]<0,2.,1.) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2])-ALAT)*60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1])-ALONG)*60
    azmuth<-aspect

    hori<-as.matrix(hori) #horizon angles
    VIEWF <- 1-sum(sin(as.data.frame(hori)*pi/180))/length(hori) # convert horizon angles to radians and calc view factor(s)
    SLES <- rep(SLE,timeinterval*nyears)
    # creating the shade array
    MAXSHADES <- rep(0,(timeinterval*nyears))+maxshade # daily max shade (%)
    MINSHADES <- rep(0,(timeinterval*nyears))+minshade # daily min shade (%)
    if(soiltype==0){ # simulating rock so turn of soil moisture model and set density equal to bulk density
      BulkDensity<-Density
      cap=0
      runmoist<-0
      PE<-rep(CampNormTbl9_1[1,4],19) #air entry potential J/kg
      KS<-rep(CampNormTbl9_1[1,6],19) #saturated conductivity, kg s/m3
      BB<-rep(CampNormTbl9_1[1,5],19) #soil 'b' parameter
      BD<-rep(BulkDensity,19) # soil bulk density, Mg/m3
      DD<-rep(Density,19) # soil density, Mg/m3
    }else{
      if(soiltype<12){ # use soil properties as specified in Campbell and Norman 1998 Table 9.1
        PE<-rep(CampNormTbl9_1[soiltype,4],19) #air entry potential J/kg
        KS<-rep(CampNormTbl9_1[soiltype,6],19) #saturated conductivity, kg s/m3
        BB<-rep(CampNormTbl9_1[soiltype,5],19) #soil 'b' parameter
        BD<-rep(BulkDensity,19) # soil bulk density, Mg/m3
        DD<-rep(Density,19) # soil density, Mg/m3
      }
    }

    if(soilgrids == 1){
      cat('extracting data from SoilGrids \n')
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("package 'jsonlite' is needed to extract data from SoilGrids, please install it.",
             call. = FALSE)
      }
      require(jsonlite)
      ov <- fromJSON(paste0('https://rest.soilgrids.org/query?lon=',x[1],'&lat=',x[2],',&attributes=BLDFIE,SLTPPT,SNDPPT,CLYPPT'), flatten = TRUE)
      if(length(ov) > 3){
        soilpro <- cbind(c(0,5,15,30,60,100,200), unlist(ov$properties$BLDFIE$M)/1000, unlist(ov$properties$SLTPPT$M), unlist(ov$properties$SNDPPT$M), unlist(ov$properties$CLYPPT$M) )
        colnames(soilpro) <- c('depth', 'blkdens', 'clay', 'silt', 'sand')
        #Now get hydraulic properties for this soil using Cosby et al. 1984 pedotransfer functions.
        soil.hydro<-pedotransfer(soilpro = as.data.frame(soilpro), DEP = DEP)
        PE<-soil.hydro$PE
        BB<-soil.hydro$BB
        BD<-soil.hydro$BD
        KS<-soil.hydro$KS
        BulkDensity <- BD[seq(1,19,2)] #soil bulk density, Mg/m3
      }else{
        cat('no SoilGrids data for this site, using user-input soil properties \n')
      }
    }
    # load global climate files
    gcfolder<-paste(.libPaths()[1],"/gcfolder.rda",sep="")
    if(file.exists(gcfolder)==FALSE){
      folder<-"c:/globalclimate"
      if(file.exists(paste0(folder,"/global_climate.nc"))==FALSE){
        message("You don't appear to have the global climate data set - \n run function get.global.climate(folder = 'folder you want to put it in') .....\n exiting function micro_global")
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))
        stop()
      }
    }else{
      load(gcfolder)
    }
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("package 'raster' is needed. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("package 'ncdf4' is needed. Please install it.",
           call. = FALSE)
    }

    message('extracting climate data \n')
    global_climate<-raster::brick(paste0(folder,"/global_climate.nc"))
    CLIMATE <- raster::extract(global_climate,x)
    ALTT<-as.numeric(CLIMATE[,1])
    delta_elev <- 0
    if(is.na(elev) == FALSE){ # check if user-specified elevation
      delta_elev <- ALTT - elev # get delta for lapse rate correction
      ALTT <- elev # now make final elevation the user-specified one
    }
    adiab_corr_max <- delta_elev * lapse_max
    adiab_corr_min <- delta_elev * lapse_min
    RAINFALL <- CLIMATE[,2:13]
    if(is.na(RAINFALL[1])){
      cat("no climate data for this site, using dummy data so solar is still produced \n")
      CLIMATE <- raster::extract(global_climate,cbind(140,-35))
      ALTT<-as.numeric(CLIMATE[,1])
      delta_elev <- 0
      if(is.na(elev) == FALSE){ # check if user-specified elevation
        delta_elev <- ALTT - elev # get delta for lapse rate correction
        ALTT <- elev # now make final elevation the user-specified one
      }
      adiab_corr_max <- delta_elev * lapse_max
      adiab_corr_min <- delta_elev * lapse_min
      RAINFALL <- CLIMATE[,2:13]*0
      #stop()
    }
    RAINYDAYS <- CLIMATE[,14:25]/10
    WNMAXX <- CLIMATE[,26:37]/10
    WNMINN<-WNMAXX*0.1 # impose diurnal cycle
    TMINN <- CLIMATE[,38:49]/10
    TMAXX <- CLIMATE[,50:61]/10
    TMAXX<-TMAXX+adiab_corr_max
    TMINN<-TMINN+adiab_corr_min
    ALLMINTEMPS<-TMINN
    ALLMAXTEMPS<-TMAXX
    ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)
    RHMINN <- CLIMATE[,62:73]/10
    RHMAXX <- CLIMATE[,74:85]/10
    # correct for potential change in RH with elevation-corrected Tair
    es <- WETAIR(db = TMAXX, rh = 100)$esat
    e <- WETAIR(db = CLIMATE[,50:61]/10, rh = CLIMATE[,62:73]/10)$e
    RHMINN <- (e / es) * 100
    RHMINN[RHMINN>100]<-100
    RHMINN[RHMINN<0]<-0.01
    es <- WETAIR(db = TMINN, rh = 100)$esat
    e <- WETAIR(db = CLIMATE[,38:49]/10, rh = CLIMATE[,74:85]/10)$e
    RHMAXX <- (e / es) * 100
    RHMAXX[RHMAXX>100]<-100
    RHMAXX[RHMAXX<0]<-0.01
    CCMINN <- CLIMATE[,86:97]/10
    if(clearsky==1){
      CCMINN=CCMINN*0
    }
    CCMAXX <- CCMINN
    if(runmoist==0){
      # extract soil moisture
      soilmoisture<-suppressWarnings(raster::brick(paste(folder,"/soilw.mon.ltm.v2.nc",sep="")))
      message("extracting soil moisture data")
      SoilMoist<-raster::extract(soilmoisture,x)/1000 # this is originally in mm/m
    }
    if(is.na(max(SoilMoist, ALTT, CLIMATE))==TRUE){
      message("Sorry, there is no environmental data for this location")
      SoilMoist<-raster::extract(soilmoisture,cbind(140,-35))/1000 # this is originally in mm/m
      #stop()
    }

    # correct for fact that wind is measured at 10 m height
    # wind shear equation v / vo = (h / ho)^a
    #where
    #v = the velocity at height h (m/s)
    #vo = the velocity at height ho (m/s)
    #a = the wind shear exponent
    #Terrain   Wind Shear Exponent
    #- a -
    #  Open water   0.1
    #Smooth, level, grass-covered   0.15 (or more commonly 1/7)
    #Row crops 	0.2
    #Low bushes with a few trees 	0.2
    #Heavy trees 	0.25
    #Several buildings 	0.25
    #Hilly, mountainous terrain 	0.25
    # source http://www.engineeringtoolbox.com/wind-shear-d_1215.html
    WNMINN<-WNMINN*(1.2/10)^0.15
    WNMAXX<-WNMAXX*(1.2/10)^0.15

    if(timeinterval!=12){ # spline from 12 days to chosen time interval
      TMAXX1 <-suppressWarnings(spline(doys12,TMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      TMAXX<-rep(TMAXX1$y,nyears)
      TMINN1 <-suppressWarnings(spline(doys12,TMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      TMINN <- rep(TMINN1$y,nyears)
      RHMAXX1 <-suppressWarnings(spline(doys12,RHMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      RHMAXX <- rep(RHMAXX1$y,nyears)
      RHMINN1 <-suppressWarnings(spline(doys12,RHMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      RHMINN <- rep(RHMINN1$y,nyears)
      CCMAXX1 <-suppressWarnings(spline(doys12,CCMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      CCMAXX <- rep(CCMAXX1$y,nyears)
      CCMINN <- CCMAXX
      WNMAXX1 <-suppressWarnings(spline(doys12,WNMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      WNMAXX<-rep(WNMAXX1$y,nyears)
      WNMINN1 <-suppressWarnings(spline(doys12,WNMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      WNMINN<-rep(WNMINN1$y,nyears)
      if(runmoist==0){
        SoilMoist1 <-suppressWarnings(spline(doys12,SoilMoist,n=timeinterval,xmin=1,xmax=365,method="periodic"))
        SoilMoist<-rep(SoilMoist1$y,nyears)
      }
    }
    if(timeinterval<365){
      TMAXX<-rep(TMAXX,nyears)
      TMINN<-rep(TMINN,nyears)
      RHMAXX<-rep(RHMAXX,nyears)
      RHMINN<-rep(RHMINN,nyears)
      CCMAXX<-rep(CCMAXX,nyears)
      CCMINN<-rep(CCMINN,nyears)
      WNMAXX<-rep(WNMAXX,nyears)
      WNMINN<-rep(WNMINN,nyears)
      if(runmoist==0){
        SoilMoist<-rep(SoilMoist,nyears)
      }
      RAINFALL<-rep(RAINFALL,nyears)
    }
    orig.RAINFALL <- RAINFALL
    # get annual mean temp for creating deep soil (2m) boundary condition
    avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
    soilinit<-rep(avetemp,20)
    tannul<-mean(unlist(ALLTEMPS))
    tannulrun<-rep(tannul,doynum)

    daymon<-c(31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.) # days in each month

    # if doing daily sims, spread rainfall evenly across days based on mean monthly rainfall and the number of rainy days per month
    if(timeinterval==365){
      RAINFALL1<-1:365
      sort<-matrix(data = 0,nrow = 365,ncol = 2)
      daymon<-c(31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.) # days in each month
      m<-1
      b<-0
      for (i in 1:12){ #begin loop throught 12 months of year
        ndays=daymon[i]
        for (k in 1:ndays){
          b<-b+1
          sort[m,1]<-i
          sort[m,2]<-b
          if(k<=RAINYDAYS[i] & rainfrac>0){
            if(k==1){
              RAINFALL1[m]=RAINFALL[i]*rainfrac*rainmult # if first day of month, make user-specified fraction of monthly rainfall fall on first day
            }else{
              RAINFALL1[m]=(RAINFALL[i]*(1-rainfrac)*rainmult)/RAINYDAYS[i] # make remaining rain fall evenly over the remaining number of rainy days for the month, starting at the beginning of the month
            }
          }else{
            if(rainfrac==0){
              RAINFALL1[m]=(RAINFALL[i]*rainmult)/RAINYDAYS[i]
            }else{
              RAINFALL1[m]=0.
            }
          }
          m<-m+1
          if(b>RAINYDAYS[i]){
            b<-0
          }
        }
      }
      RAINFALL2<-as.data.frame(cbind(RAINFALL1,sort))
      #RAINFALL2<-RAINFALL2[order(RAINFALL2$V2,RAINFALL2$V3),] # this line scatters the rainy days evenly across each month - snow predictions better if it is commented out so get rainy days all in a row within the month
      RAINFALL<-rep(as.double(RAINFALL2$RAINFALL1),nyears)
      if(TMINN[1]<snowtemp){
        RAINFALL[1]<-0 # this is needed in some cases to allow the integrator to get started
      }
    }else{
      if(timeinterval!=12){
        RAINFALL<-rep(rep(sum(RAINFALL)/timeinterval,timeinterval),nyears) # just spread evenly across every day
      }else{ # running middle day of each month - divide monthly rain by number of days in month
        RAINFALL<-RAINFALL/rep(daymon,nyears)
      }
    }#end check doing daily sims
    ndays<-length(RAINFALL)
    if(length(TAI) < 111){ # no user supplied values, compute with GADS
      if(run.gads==1){
        ####### get solar attenuation due to aerosols with program GADS #####################
        relhum<-1.
        optdep.summer<-as.data.frame(rungads(longlat[2],longlat[1],relhum,0))
        optdep.winter<-as.data.frame(rungads(longlat[2],longlat[1],relhum,1))
        optdep<-cbind(optdep.winter[,1],rowMeans(cbind(optdep.summer[,2],optdep.winter[,2])))
        optdep<-as.data.frame(optdep)
        colnames(optdep)<-c("LAMBDA","OPTDEPTH")
        a<-lm(OPTDEPTH~poly(LAMBDA, 6, raw=TRUE),data=optdep)
        LAMBDA<-c(290,295,300,305,310,315,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1380,1400,1420,1440,1460,1480,1500,1540,1580,1600,1620,1640,1660,1700,1720,1780,1800,1860,1900,1950,2000,2020,2050,2100,2120,2150,2200,2260,2300,2320,2350,2380,2400,2420,2450,2490,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000)
        TAI<-predict(a,data.frame(LAMBDA))
        ################ end GADS ##################################################
      }else{ # use the original profile from Elterman, L. 1970. Vertical-attenuation model with eight surface meteorological ranges 2 to 13 kilometers. U. S. Airforce Cambridge Research Laboratory, Bedford, Mass.
        TAI<-c(0.42,0.415,0.412,0.408,0.404,0.4,0.395,0.388,0.379,0.379,0.379,0.375,0.365,0.345,0.314,0.3,0.288,0.28,0.273,0.264,0.258,0.253,0.248,0.243,0.236,0.232,0.227,0.223,0.217,0.213,0.21,0.208,0.205,0.202,0.201,0.198,0.195,0.193,0.191,0.19,0.188,0.186,0.184,0.183,0.182,0.181,0.178,0.177,0.176,0.175,0.175,0.174,0.173,0.172,0.171,0.17,0.169,0.168,0.167,0.164,0.163,0.163,0.162,0.161,0.161,0.16,0.159,0.157,0.156,0.156,0.155,0.154,0.153,0.152,0.15,0.149,0.146,0.145,0.142,0.14,0.139,0.137,0.135,0.135,0.133,0.132,0.131,0.13,0.13,0.129,0.129,0.128,0.128,0.128,0.127,0.127,0.126,0.125,0.124,0.123,0.121,0.118,0.117,0.115,0.113,0.11,0.108,0.107,0.105,0.103,0.1)
      } #end check if running gads
    }
    ################ soil properties  ##################################################
    Nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
    if(soilgrids == 1){
      Numtyps <- 10 # number of substrate types
      Nodes[1:10,] <- c(1:10) # deepest nodes for each substrate type
    }else{
      Numtyps <- 2 # number of soil types
      Nodes[1,1:ndays]<-3 # deepest node for first substrate type
      Nodes[2,1:ndays]<-9 # deepest node for second substrate type
    }
    REFLS<-rep(REFL,ndays) # soil reflectances
    PCTWET<-rep(PCTWET,ndays) # soil wetness
    if(runmoist==0){
      moists2<-matrix(nrow= 10, ncol = ndays, data=0) # set up an empty vector for soil moisture values through time
      moists2[1,]<-SoilMoist # fill the first row with monthly soil moisture values
      moists2[2,]<-moists2[1,] # make this row same as first row
      moists<-moists2
    }else{
      moists2<-matrix(nrow=10, ncol = ndays, data=0) # set up an empty vector for soil moisture values through time
      moists2[1:10,]<-SoilMoist_Init
      moists2[moists2>(1-BulkDensity/Density)]<-(1-BulkDensity/Density)
      moists<-moists2
    }

    # now make the soil properties matrix
    # columns are:
    #1) bulk density (Mg/m3)
    #2) volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    #3) thermal conductivity (W/mK)
    #4) specific heat capacity (J/kg-K)
    #5) mineral density (Mg/m3)
    soilprops<-matrix(data = 0, nrow = 10, ncol = 5) # create an empty soil properties matrix
    if(soilgrids == 1){
      soilprops[,1]<-BulkDensity
      soilprops[,2]<-min(0.26, 1 - BulkDensity / Density) # not used if soil moisture computed
      soilprops[,3]<-Thcond
      soilprops[,4]<-SpecHeat
      soilprops[,5]<-Density
      if(cap==1){
        soilprops[1:2,3]<-0.2
        soilprops[1:2,4]<-1920
      }
      if(cap==2){
        soilprops[1:2,3]<-0.1
        soilprops[3:4,3]<-0.25
        soilprops[1:4,4]<-1920
        soilprops[1:4,5]<-1.3
        soilprops[1:4,1]<-0.7
      }
    }else{
      soilprops[1,1]<-BulkDensity # insert soil bulk density to profile 1
      soilprops[2,1]<-BulkDensity # insert soil bulk density to profile 2
      soilprops[1,2]<-min(0.26, 1 - BulkDensity / Density) # insert saturated water content to profile 1
      soilprops[2,2]<-min(0.26, 1 - BulkDensity / Density) # insert saturated water content to profile 2
      if(cap==1){ # insert thermal conductivity to profile 1, and see if 'organic cap' added on top
        soilprops[1,3]<-0.2 # mineral thermal conductivity
      }else{
        soilprops[1,3]<-Thcond # mineral thermal conductivity
      }
      soilprops[2,3]<-Thcond # insert thermal conductivity to profile 2
      if(cap==1){ # insert specific heat to profile 1, and see if 'organic cap' added on top
        soilprops[1,4]<-1920 # mineral heat capacity
      }else{
        soilprops[1,4]<-SpecHeat
      }
      soilprops[2,4]<-SpecHeat # insert specific heat to profile 2
      soilprops[1,5]<-Density # insert mineral density to profile 1
      soilprops[2,5]<-Density # insert mineral density to profile 2
    }
    #########################################################################################

    # Next four parameters are segmented velocity profiles due to bushes, rocks etc. on the surface
    #IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO! (then roughness height is based on the parameter RUF)
    Z01 <- 0 # Top (1st) segment roughness height(m)
    Z02 <- 0 # 2nd segment roughness height(m)
    ZH1 <- 0 # Top of (1st) segment, height above surface(m)
    ZH2 <- 0 # 2nd segment, height above surface(m)

    # hourly option set to 0, so make empty vectors
    hourly<-0
    rainhourly<-0
    TAIRhr<-rep(0,24*ndays)
    RHhr<-rep(0,24*ndays)
    WNhr<-rep(0,24*ndays)
    CLDhr<-rep(0,24*ndays)
    SOLRhr<-rep(0,24*ndays)
    RAINhr<-rep(0,24*ndays)
    ZENhr<-rep(-1,24*ndays)
    IRDhr<-rep(-1,24*ndays)
    # microclimate input parameters list
    microinput<-c(ndays,RUF,ERR,Usrhyt,Refhyt,Numtyps,Z01,Z02,ZH1,ZH2,idayst,ida,HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,CMH2O,microdaily,tannul,EC,VIEWF,snowtemp,snowdens,snowmelt,undercatch,rainmult,runshade,runmoist,maxpool,evenrain,snowmodel,rainmelt,writecsv,densfun,hourly,rainhourly,lamb,IUV,RW,PC,RL,SP,R1,IM,MAXCOUNT,IR,message,fail,snowcond,intercept,grasshade,solonly,ZH,D0)

    doy1<-matrix(data <- 0., nrow = ndays, ncol = 1)
    SLES1<-matrix(data = 0., nrow = ndays, ncol = 1)
    MAXSHADES1<-matrix(data = 0., nrow = ndays, ncol = 1)
    MINSHADES1<-matrix(data = 0., nrow = ndays, ncol = 1)
    TMAXX1<-matrix(data = 0., nrow = ndays, ncol = 1)
    TMINN1<-matrix(data = 0., nrow = ndays, ncol = 1)
    CCMAXX1<-matrix(data = 0., nrow = ndays, ncol = 1)
    CCMINN1<-matrix(data = 0., nrow = ndays, ncol = 1)
    RHMAXX1<-matrix(data = 0., nrow = ndays, ncol = 1)
    RHMINN1<-matrix(data = 0., nrow = ndays, ncol = 1)
    WNMAXX1<-matrix(data = 0., nrow = ndays, ncol = 1)
    WNMINN1<-matrix(data = 0., nrow = ndays, ncol = 1)
    REFLS1<-matrix(data = 0., nrow = ndays, ncol = 1)
    PCTWET1<-matrix(data = 0., nrow = ndays, ncol = 1)
    RAINFALL1<-matrix(data = 0, nrow = ndays, ncol = 1)
    tannul1<-matrix(data = 0, nrow = ndays, ncol = 1)
    moists1<-matrix(data = 0., nrow = 10, ncol = ndays)
    doy1[1:ndays]<-doy
    SLES1[1:ndays]<-SLES
    MAXSHADES1[1:ndays]<-MAXSHADES
    MINSHADES1[1:ndays]<-MINSHADES
    TMAXX1[1:ndays]<-TMAXX
    TMINN1[1:ndays]<-TMINN
    CCMAXX1[1:ndays]<-CCMAXX
    CCMINN1[1:ndays]<-CCMINN
    RHMAXX1[1:ndays]<-RHMAXX
    RHMINN1[1:ndays]<-RHMINN
    WNMAXX1[1:ndays]<-WNMAXX
    WNMINN1[1:ndays]<-WNMINN
    REFLS1[1:ndays]<-REFLS
    PCTWET1[1:ndays]<-PCTWET
    RAINFALL1[1:ndays]<-RAINFALL
    tannul1[1:ndays]<-tannul
    moists1[1:10,1:ndays]<-moists
    if(length(LAI)<ndays){
      LAI<-rep(LAI[1],ndays)
      LAI1 <- LAI
    }
    if(shore==0){
      tides<-matrix(data = 0., nrow = 24*ndays, ncol = 3) # make an empty matrix
    }
    # all microclimate data input list - all these variables are expected by the input argument of the fortran micro2014 subroutine
    micro<-list(tides=tides,microinput=microinput,doy=doy,SLES=SLES1,DEP=DEP,Nodes=Nodes,MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TIMAXS=TIMAXS,TIMINS=TIMINS,TMAXX=TMAXX1,TMINN=TMINN1,RHMAXX=RHMAXX1,RHMINN=RHMINN1,CCMAXX=CCMAXX1,CCMINN=CCMINN1,WNMAXX=WNMAXX1,WNMINN=WNMINN1,TAIRhr=TAIRhr,RHhr=RHhr,WNhr=WNhr,CLDhr=CLDhr,SOLRhr=SOLRhr,RAINhr=RAINhr,ZENhr=ZENhr,IRDhr=IRDhr,REFLS=REFLS1,PCTWET=PCTWET1,soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists1,RAINFALL=RAINFALL1,tannulrun=tannulrun,PE=PE,KS=KS,BB=BB,BD=BD,DD=DD,L=L,LAI=LAI)

    # write all input to csv files in their own folder
    if(write_input==1){
      if(dir.exists("micro csv input")==FALSE){
        dir.create("micro csv input")
      }
      write.table(as.matrix(microinput), file = "micro csv input/microinput.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(doy, file = "micro csv input/doy.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(SLES, file = "micro csv input/SLES.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(DEP, file = "micro csv input/DEP.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(Nodes, file = "micro csv input/Nodes.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(MAXSHADES, file = "micro csv input/Maxshades.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(MINSHADES, file = "micro csv input/Minshades.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TIMAXS, file = "micro csv input/TIMAXS.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TIMINS, file = "micro csv input/TIMINS.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TMAXX, file = "micro csv input/TMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TMINN, file = "micro csv input/TMINN.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(RHMAXX, file = "micro csv input/RHMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(RHMINN, file = "micro csv input/RHMINN.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(CCMAXX, file = "micro csv input/CCMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(CCMINN, file = "micro csv input/CCMINN.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(WNMAXX, file = "micro csv input/WNMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(WNMINN, file = "micro csv input/WNMINN.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(REFLS, file = "micro csv input/REFLS.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(PCTWET, file = "micro csv input/PCTWET.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(soilinit, file = "micro csv input/soilinit.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(hori, file = "micro csv input/hori.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TAI, file = "micro csv input/TAI.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(soilprops, file="micro csv input/soilprop.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(moists,file="micro csv input/moists.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(RAINFALL,file="micro csv input/rain.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(tannulrun,file="micro csv input/tannulrun.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(PE,file="micro csv input/PE.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(BD,file="micro csv input/BD.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(DD,file="micro csv input/DD.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(BB,file="micro csv input/BB.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(KS,file="micro csv input/KS.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(L,file="micro csv input/L.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(LAI,file="micro csv input/LAI.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(tides,file="micro csv input/tides.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TAIRhr,file="micro csv input/TAIRhr.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(RHhr,file="micro csv input/RHhr.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(WNhr,file="micro csv input/WNhr.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(CLDhr,file="micro csv input/CLDhr.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(SOLRhr,file="micro csv input/SOLRhr.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(RAINhr,file="micro csv input/RAINhr.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(ZENhr,file="micro csv input/ZENhr.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(IRDhr,file="micro csv input/IRDhr.csv", sep = ",", col.names = NA, qmethod = "double")
    }

    if(is.numeric(loc[1])){
      location<-paste("long",loc[1],"lat",loc[2])
    }else{
      location<-loc
    }
    message(paste('running microclimate model for',timeinterval,'days by',nyears,'years at site',location,'\n'))
    ptm <- proc.time() # Start timing
    microut<-microclimate(micro)
    message(paste0('runtime ', (proc.time() - ptm)[3], ' seconds')) # Stop the clock

    metout<-microut$metout # retrieve above ground microclimatic conditions, min shade
    shadmet<-microut$shadmet # retrieve above ground microclimatic conditions, max shade
    soil<-microut$soil # retrieve soil temperatures, minimum shade
    shadsoil<-microut$shadsoil # retrieve soil temperatures, maximum shade
    if(runmoist==1){
      soilmoist<-microut$soilmoist # retrieve soil moisture, minimum shade
      shadmoist<-microut$shadmoist # retrieve soil moisture, maximum shade
      humid<-microut$humid # retrieve soil humidity, minimum shade
      shadhumid<-microut$shadhumid # retrieve soil humidity, maximum shade
      soilpot<-microut$soilpot # retrieve soil water potential, minimum shade
      shadpot<-microut$shadpot # retrieve soil water potential, maximum shade
      plant<-microut$plant # retrieve plant output, minimum shade
      shadplant<-microut$shadplant # retrieve plant output, maximum shade
    }else{
      soilpot<-soil
      soilmoist<-soil
      shadpot<-soil
      shadmoist<-soil
      humid<-soil
      shadhumid<-soil
      plant<-cbind(soil,soil[,3:4])
      shadplant<-cbind(soil,soil[,3:4])
      soilpot[,3:12]<-0
      soilmoist[,3:12]<-0.5
      shadpot[,3:12]<-0
      shadmoist[,3:12]<-0.5
      humid[,3:12]<-0.99
      shadhumid[,3:12]<-0.99
      plant[,3:14]<-0
      shadplant[,3:14]<-0
    }
    if(snowmodel == 1){
      sunsnow <- microut$sunsnow
      shdsnow <- microut$shdsnow
    }
    if(timeinterval == 12){
      RAINFALL <- orig.RAINFALL
    }
    if(max(metout[,1] == 0)){
      cat("ERROR: the model crashed - try a different error tolerance (ERR) or a different spacing in DEP", '\n')
    }
    days<-rep(seq(1,timeinterval * nyears), 24)
    days<-days[order(days)]
    dates<-days + metout[, 2] / 60 / 24 - 1 # dates for hourly output
    dates2<-seq(1, timeinterval * nyears) # dates for daily output
    if(lamb == 1){
      drlam<-as.data.frame(microut$drlam) # retrieve direct solar irradiance
      drrlam<-as.data.frame(microut$drrlam) # retrieve direct Rayleigh component solar irradiance
      srlam<-as.data.frame(microut$srlam) # retrieve scattered solar irradiance
      if(snowmodel == 1){
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=minshade,maxshade=maxshade,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam,dates=dates,dates2=dates2))
      }else{
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=minshade,maxshade=maxshade,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam,dates=dates,dates2=dates2))
      }
    }else{
      if(snowmodel == 1){
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=minshade,maxshade=maxshade,DEP=DEP,dates=dates,dates2=dates2))
      }else{
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=minshade,maxshade=maxshade,DEP=DEP,dates=dates,dates2=dates2))
      }
    }
  } # end error trapping
} # end of micro_global function
