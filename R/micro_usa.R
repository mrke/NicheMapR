#' USA implementation of the microclimate model.
#'
#' An implementation of the NicheMapR microclimate model that uses the GRIDMET daily weather database http://www.climatologylab.org/gridmet.html, and specifically uses the following variables: pr, rmax, rmin, srad, tmmn, tmmx, vs. Also uses the following DEM "metdata_elevationdata.nc".
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param dstart First day to run, date in format "d/m/Y" e.g. "01-01-2016"
#' @param dfinish Last day to run, date in format "d/m/Y" e.g. "31-12-2016"
#' @param REFL Soil solar reflectance, decimal \%
#' @param elev Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param minshade Minimum shade level to use (\%)
#' @param maxshade Maximum shade level to us (\%)
#' @param Usrhyt Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
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
#' @usage micro_USA(loc = c(-117.8919161, 36.582016), dstart = "01/01/2016", dfinish = "31/12/2016",
#' REFL = 0.15, slope = 0, aspect = 0, DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minshade = 0, maxshade = 90,
#' Usrhyt = 0.01, ...)
#' @export
#' @details
#' \strong{Parameters controlling how the model runs:}\cr\cr
#'
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes, 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{solonly}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{terrain}{ = 0, Use 250m resolution terrain data? 1=yes, 0=no}\cr\cr
#' \code{dailywind}{ = 1, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{windfac}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{adiab_cor}{ = 1, use adiabatic lapse rate correction? 1=yes, 0=no}\cr\cr
#' \code{warm}{ = 0, uniform warming, °C}\cr\cr
#' \code{spatial}{ = "c:/Australian Environment/", choose location of terrain data}\cr\cr
#' \code{opendap}{ = 1, query met grids via opendap}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#' \code{save}{ = 0, don't save forcing data (0), save the forcing data (1) or read previously saved data (2)}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{ERR}{ = 1.5, Integrator error tolerance for soil temperature calculations}\cr\cr
#' \code{Refhyt}{ = 2, Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured}\cr\cr
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
#' \code{BulkDensity}{ = 1.3, Soil bulk density (Mg/m3)}\cr\cr
#' \code{PCTWET}{ = 0, \% of ground surface area acting as a free water surface}\cr\cr
#' \code{rainwet}{ = 1.5, mm of rainfall causing the ground to be 90\% wet for the day}\cr\cr
#' \code{cap}{ = 1, organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)}\cr\cr
#' \code{CMH2O}{ = 1, Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air conditions; 2.0 = humid, tropical conditions (note this is for the whole atmospheric profile, not just near the ground)}\cr\cr
#' \code{hori}{ = rep(0,24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#' \code{lapse_min}{ = 0.0039 Lapse rate for minimum air temperature (degrees C/m)}
#' \code{lapse_max}{ = 0.0077 Lapse rate for maximum air temperature (degrees C/m)}
#' \code{TIMAXS}{ = c(1.0, 1.0, 0.0, 0.0), Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise}\cr\cr
#' \code{TIMINS}{ = c(0, 0, 1, 1), Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon}\cr\cr
#' \code{timezone}{ = 0, Use GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)? 1=yes, 0=no}\cr\cr
#'
#' \strong{ Soil moisture mode parameters:}
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
#' \code{rainmult}{ = 1, Rain multiplier for surface soil moisture (-), used to induce runon}\cr\cr
#' \code{rainoff}{ = 0, Rain offset (mm), used to induce constant extra input}\cr\cr
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
#' \code{snowmodel}{ = 1, run the snow model 1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (mg/m3), overridden by densfun}\cr\cr
#' \code{densfun}{ = c(0.5979, 0.2178, 0.001, 0.0038), slope and intercept of model of snow density as a linear function of snowpack age if first two values are nonzero, and following the exponential function of Sturm et al. 2010 J. of Hydromet. 11:1380-1394 if all values are non-zero; if it is c(0,0,0,0) then fixed density used}\cr\cr
#' \code{snowmelt}{ = 1, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
#' \code{undercatch}{ = 1, undercatch multipier for converting rainfall to snow}\cr\cr
#' \code{rainmelt}{ = 0.0125, paramter in equation that melts snow with rainfall as a function of air temp}\cr\cr
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)}\cr\cr
#' \code{intercept}{ = maxshade / 100 * 0.3, snow interception fraction for when there's shade (0-1)}\cr\cr
#' \code{grasshade}{ = 0, if 1, means shade is removed when snow is present, because shade is cast by grass/low shrubs}\cr\cr
#'
#' \strong{ Intertidal mode parameters:}
#'
#' \code{shore}{ Include tide effects? If 1, the matrix}
#' \code{tides}
#' { is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0, nrow = length(seq(as.POSIXct(dstart, format = "%d/%m/%Y"), as.POSIXct(dfinish, format = "%d/%m/%Y"), by = "days")) * 24, ncol = 3), matrix for each how of the simulation of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#' }
#'
#' \strong{Outputs:}
#'
#' \code{ndays}{ - number of days for which predictions are made}\cr\cr
#' \code{longlat}{ - longitude and latitude for which simulation was run (decimal degrees)}\cr\cr
#' \code{dates}{ - vector of dates (hourly, POSIXct, timezone = America/Los_Angeles)}\cr\cr
#' \code{dates2}{ - vector of dates (daily, POSIXct, timezone = America/Los_Angeles)}\cr\cr
#' \code{nyears}{ - number of years for which predictions are made}\cr\cr
#' \code{RAINFALL}{ - vector of daily rainfall (mm)}\cr\cr
#' \code{elev}{ - elevation at point of simulation (m)}\cr\cr
#' \code{minshade}{ - minimum shade for simulation (\%)}\cr\cr
#' \code{maxshade}{ - maximum shade for simulation (single value - if time varying, in 'MAXSHADES') (\%)}\cr\cr
#' \code{MAXSHADES}{ - vector of maximum shades used (\%)}\cr\cr
#' \code{DEP}{ - vector of depths used (cm)}\cr\cr
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
#' \item 13 SOLR - solar radiation (W/m2)
#' \item 14 TSKYC - sky radiant temperature (°C)
#' \item 15 DEW - dew presence (0 or 1)
#' \item 16 FROST - frost presence (0 or 1)
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
#' library(NicheMapR)
#' dstart <- "01/01/2016"
#' dfinish <- "31/12/2017"
#' micro<-micro_usa() # run the model at the default location (Madison, Wisconsin) for 2016 to 2017 using opendap
#'
#' metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#' soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
#' soilmoist<-as.data.frame(micro$soilmoist) # soil temperatures, minimum shade
#'
#' # append dates
#' dates <- micro$dates
#'
#' metout <- cbind(dates, metout)
#' soil <- cbind(dates, soil)
#' soilmoist <- cbind(dates, soilmoist)
#' minshade<-micro$minshade
#'
#' # plotting above-ground conditions in minimum shade
#' with(metout, {plot(TALOC ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l", main = paste("air and sky temperature, ", minshade, "% shade", sep = ""), ylim = c(-20, 60))})
#' with(metout, {points(TAREF ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l", lty = 2,col = 'blue')})
#' with(metout, {points(TSKYC ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' ,  type = "l", col = 'light blue', main = paste("sky temperature, ", minshade, "% shade", sep = ""))})
#' with(metout, {plot(RHLOC ~ dates, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l", ylim = c(0, 100), main = paste("humidity, ", minshade, "% shade", sep = ""))})
#' with(metout, {points(RH ~ dates, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l", col = 'blue', lty = 2, ylim = c(0, 100))})
#' with(metout, {plot(VREF ~ dates, xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#' ,  type = "l", main = "wind speed", ylim = c(0, 15))})
#' with(metout, {points(VLOC ~ dates, xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#' ,  type = "l", lty = 2, col = 'blue')})
#' with(metout, {plot(SOLR ~ dates, xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
#' ,  type = "l", main = "solar radiation")})
#' with(metout, {plot(SNOWDEP ~ dates, xlab = "Date and Time", ylab = "Snow Depth (cm)"
#' ,  type = "l", main = "snow depth")})
#'
#' # plotting soil temperature
#' for(i in 1:10){
#'  if(i==1){
#'    plot(soil[,i+3] ~ soil[,1], xlab = "Date and Time", ylab = "Soil Temperature (°C)"
#'    ,col = i,type = "l", main = paste("soil temperature ", minshade, "% shade", sep = ""))
#'  }else{
#'    points(soil[, i + 3] ~ soil[, 1], xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l")
#'  }
#' }
#'
#' # plotting soil moisture
#' for(i in 1:10){
#'  if(i==1){
#'    plot(soilmoist[,i+3] * 100 ~ soilmoist[, 1], xlab = "Date and Time", ylab = "Soil Moisture (% volumetric)"
#'    ,col = i,type = "l", main = paste("soil moisture ", minshade, "% shade", sep = ""))
#'  }else{
#'    points(soilmoist[,i+3] * 100 ~ soilmoist[, 1], xlab = "Date and Time", ylab = "Soil Moisture
#'     (%)", col = i, type = "l")
#'  }
#' }
micro_usa <- function(
  loc = c(-89.40123, 43.07305),
  dstart = "01/01/2016",
  dfinish = "31/12/2017",
  nyears = as.numeric(substr(dfinish, 7, 10)) - as.numeric(substr(dstart, 7, 10)) + 1,
  REFL = 0.15,
  elev = NA,
  slope = 0,
  aspect = 0,
  lapse_max = 0.0077,
  lapse_min = 0.0039,
  DEP=c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
  minshade = 0,
  maxshade = 90,
  Refhyt = 2,
  Usrhyt = 0.01,
  Z01 = 0,
  Z02 = 0,
  ZH1 = 0,
  ZH2 = 0,
  runshade = 1,
  clearsky = 0,
  solonly = 0,
  run.gads = 1,
  write_input = 0,
  writecsv = 0,
  terrain = 0,
  dailywind = 1,
  windfac = 1,
  adiab_cor = 1,
  warm = 0,
  spatial = "P:",
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
  rainwet = 1.5,
  cap = 1,
  CMH2O = 1,
  hori = rep(0,24),
  TIMAXS=c(1.0, 1.0, 0.0, 0.0),
  TIMINS = c(0, 0, 1, 1),
  timezone = 0,
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
  scenario = "",
  year = "",
  hourly = 0,
  rainhourly = 0,
  rainhour = 0,
  rainoff = 0,
  lamb = 0,
  IUV = 0,
  opendap = 1,
  soilgrids = 0,
  IR = 0,
  message = 0,
  fail = nyears * 24 * 365,
  save = 0,
  snowcond = 0,
  intercept = maxshade / 100 * 0.4,
  grasshade = 0) { # end function parameters

  ystart <- as.numeric(substr(dstart, 7, 10))
  yfinish <- as.numeric(substr(dfinish, 7, 10))
  yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
  # error trapping - originally inside the Fortran code, but now checking before executing Fortran
  errors<-0
  if(DEP[2]-DEP[1]>3 | DEP[3]-DEP[2]>3){
    cat("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if(DEP[2]-DEP[1]<2){
    cat("warning, nodes might be too close near the surface, try a different spacing if the program is crashing \n")
  }
  if(is.numeric(loc[1])){
    if(loc[1]>180 | loc[2] > 90){
      cat("ERROR: Latitude or longitude (longlat) is out of bounds.
          Please enter a correct value.", '\n')
      errors<-1
    }
  }
  if(timezone%in%c(0,1)==FALSE){
    cat("ERROR: the variable 'timezone' be either 0 or 1.
        Please correct.", '\n')
    errors<-1
  }
  if(run.gads%in%c(0,1)==FALSE){
    cat("ERROR: the variable 'run.gads' be either 0 or 1.
        Please correct.", '\n')
    errors<-1
  }
  if(write_input%in%c(0,1)==FALSE){
    cat("ERROR: the variable 'write_input' be either 0 or 1.
        Please correct.", '\n')
    errors<-1
  }
  if(EC<0.0034 | EC > 0.058){
    cat("ERROR: the eccentricity variable (EC) is out of bounds.
        Please enter a correct value (0.0034 - 0.058).", '\n')
    errors<-1
  }
  if(RUF<0.0001){
    cat("ERROR: The roughness height (RUF) is too small ( < 0.0001).
        Please enter a larger value.", '\n')
    errors<-1
  }
  if(RUF>2){
    cat("ERROR: The roughness height (RUF) is too large ( > 2).
        Please enter a smaller value.", '\n')
    errors<-1
  }
  if(DEP[1]!=0){
    cat("ERROR: First soil node (DEP[1]) must = 0 cm.
        Please correct", '\n')
    errors<-1
  }
  if(length(DEP)!=10){
    cat("ERROR: You must enter 10 different soil depths.", '\n')
    errors<-1
  }
  for(i in 1:9){
    if(DEP[i+1]<=DEP[i]){
      cat("ERROR: Soil depth (DEP array) is not in ascending size", '\n')
      errors<-1
    }
  }
  if(DEP[10]>500){
    cat("ERROR: Deepest soil depth (DEP array) is too large (<=500 cm)", '\n')
    errors<-1
  }
  if(Thcond<0){
    cat("ERROR: Thermal variable conductivity (THCOND) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(Density<0){
    cat("ERROR: Density variable (Density) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(SpecHeat<0){
    cat("ERROR: Specific heat variable (SpecHeat) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(min(BulkDensity)<0){
    cat("ERROR: Bulk density value (BulkDensity) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(REFL<0 | REFL>1){
    cat("ERROR: Soil reflectivity value (REFL) is out of bounds.
        Please input a value between 0 and 1.", '\n')
    errors<-1
  }
  if(slope<0 | slope>90){
    cat("ERROR: Slope value (slope) is out of bounds.
        Please input a value between 0 and 90.", '\n')
    errors<-1
  }
  if(aspect<0 | aspect>365){
    cat("ERROR: Aspect value (aspect) is out of bounds.
        Please input a value between 0 and 365.", '\n')
    errors<-1
  }
  if(max(hori)>90 | min(hori)<0){
    cat("ERROR: At least one of your horizon angles (hori) is out of bounds.
        Please input a value between 0 and 90", '\n')
    errors<-1
  }
  if(length(hori)!=24){
    cat("ERROR: You must enter 24 horizon angle values.", '\n')
    errors<-1
  }
  if(SLE<0.5 | SLE > 1){
    cat("ERROR: Emissivity (SLE) is out of bounds.
        Please enter a correct value (0.05 - 1.00).", '\n')
    errors<-1
  }
  if(ERR<0){
    cat("ERROR: Error bound (ERR) is too small.
        Please enter a correct value (> 0.00).", '\n')
    errors<-1
  }
  if(Usrhyt<RUF){
    cat("ERROR: Reference height (Usrhyt) smaller than roughness height (RUF).
        Please use a larger height above the surface.", '\n')
    errors<-1
  }
  if(Usrhyt<0.005 | Usrhyt>Refhyt){
    cat("ERROR: Local height (Usrhyt) is out of bounds.
        Please enter a correct value (0.005 - Refhyt).", '\n')
    errors<-1
  }
  if(CMH2O<0.5 | CMH2O>2){
    cat("ERROR: Preciptable water in air column (CMH2O) is out of bounds.
        Please enter a correct value (0.1 - 2cm).", '\n')
    errors<-1
  }
  if(max(TIMAXS)>24 | min(TIMAXS)<0){
    cat("ERROR: At least one of your times of weather maxima (TIMAXS) is out of bounds.
        Please input a value between 0 and 24", '\n')
    errors<-1
  }
  if(max(TIMINS)>24 | min(TIMINS)<0){
    cat("ERROR: At least one of your times of weather minima (TIMINS) is out of bounds.
        Please input a value between 0 and 24", '\n')
    errors<-1
  }
  if(minshade>maxshade | minshade==maxshade){
    cat("ERROR: Your value for minimum shade (minshade) is greater than or equal to the maximum shade (maxshade).
        Please correct this.", '\n')
    errors<-1
  }
  if(minshade>100 | minshade<0){
    cat("ERROR: Your value for minimum shade (minshade) is out of bounds.
        Please input a value between 0 and 100.", '\n')
    errors<-1
  }
  if(maxshade>100 | maxshade<0){
    cat("ERROR: Your value for maximum shade (maxshade) is out of bounds.
        Please input a value between 0 and 100.", '\n')
    errors<-1
  }
  # end error trapping

  if(errors==0){ # continue

    ################## time related variables #################################
    doys12<-c(15.,46.,74.,105.,135.,166.,196.,227.,258.,288.,319.,349.) # middle day of each month

    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in

    daystart<-1
    idayst <- 1 # start day

    ################## location and terrain #################################
    if (!require("raster", quietly = TRUE)) {
      stop("package 'raster' is needed. Please install it.",
        call. = FALSE)
    }
    if (!require("RNetCDF", quietly = TRUE)) {
      stop("package 'RNetCDF' is needed. Please install it.",
        call. = FALSE)
    }
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    require("raster")
    require("RNetCDF")

    # get the local timezone reference longitude
    if(timezone==1){ # this now requires registration
      if(!require(geonames, quietly = TRUE)){
        stop('package "geonames" is required to do a specific time zone (timezone=1). Please install it.')
      }
      ALREF<-(geonames::GNtimezone(longlat[2],longlat[1])[4])*-15
    }else{  # just use local solar noon
      ALREF <- abs(trunc(x[1]))
    }
    HEMIS <- ifelse(x[2]<0, 2, 1) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2])-ALAT)*60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1])-ALONG)*60
    azmuth<-aspect

    if(save != 2){
      #GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.017453292519943295]]
      if(opendap == 1){
        if(!require(futile.logger, quietly = TRUE)){
         stop('package "futile.logger" is required for this opendap extraction process. Please install it.')
        }
        require(utils)
        retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=50, sleep=1) {
          attempts = 0
          retval = try(eval(expr))
          while (isError(retval)) {
            attempts = attempts + 1
            if (attempts >= maxErrors) {
              msg = sprintf("retry: too many retries [[%s]]", capture.output(str(retval)))
              futile.logger::flog.fatal(msg)
              stop(msg)
            } else {
              msg = sprintf("retry: error in attempt %i/%i [[%s]]", attempts, maxErrors,
                capture.output(str(retval)))
              futile.logger::flog.error(msg)
              warning(msg)
            }
            if (sleep > 0) Sys.sleep(sleep)
            retval = try(eval(expr))
          }
          return(retval)
        }

        cat("extracting elevation via opendaps \n")
        baseurl <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/"
        nc <- RNetCDF::open.nc(paste0(baseurl, "/elev/metdata_elevationdata.nc"))
        lon <- RNetCDF::var.get.nc(nc, "lon", unpack = TRUE)
        lat <- RNetCDF::var.get.nc(nc, "lat", unpack = TRUE)
        flat=match(abs(lat-x[2])<1/48,1)
        latindex=which(flat %in% 1)
        flon=match(abs(lon-x[1])<1/48,1)
        lonindex=which(flon %in% 1)
        start <- c(lonindex, latindex, 1)
        count <- c(1, 1, NA)
        USADEM <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "elevation",
          start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
      }else{
        USADEM <- extract(raster::raster(paste0(spatial,"/metdata_elevationdata.nc")), x) # metres above sea level
      }
      if(save == 1){
        cat("saving DEM data for later \n")
        save(USADEM, file = 'USADEM.Rda')
      }
    }
    if(save == 2){
      cat("loading DEM data from previous run \n")
      load('USADEM.Rda')
    }

    ALTITUDES <- NA# extract(raster(paste0(spatial,"/terr50.tif")), x) # to do
    if(is.na(elev) == FALSE){ALTITUDES <- elev} # check if user-specified elevation
    if(is.na(ALTITUDES)==TRUE){ALTITUDES<-USADEM}
    if(save != 2){
      if(soilgrids == 1){
        cat('extracting soil texture data from SoilGrids \n')
        require(jsonlite)
        ov <- jsonlite::fromJSON(paste0('https://rest.soilgrids.org/query?lon=',x[1],'&lat=',x[2],',&attributes=BLDFIE,SLTPPT,SNDPPT,CLYPPT'), flatten = TRUE)
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
    if(terrain==1){ # to do
      cat("extracting terrain data \n")
      # now extract terrain data from elevslpasphori.nc
      # get UTM from dec degrees, NZTM
      HORIZONS <- elevslpasphori[4:27]
      SLOPES <- elevslpasphori[2]
      AZMUTHS <- elevslpasphori[3]
      # the horizons have been arranged so that they go from 0 degrees azimuth (north) clockwise - r.horizon starts
      # in the east and goes counter clockwise!
      HORIZONS <- (ifelse(is.na(HORIZONS),0,HORIZONS))/10 # get rid of na and get back to floating point
      HORIZONS <- data.frame(HORIZONS)
    }else{
      HORIZONS <- hori
      HORIZONS <- data.frame(HORIZONS)
      VIEWF_all <- 1-sum(sin(as.data.frame(hori)*pi/180))/length(hori) # convert horizon angles to radians and calc view factor(s)
      SLOPES<-rep(slope,length(x[,1]))
      AZMUTHS<-rep(aspect,length(x[,1]))
    }
    hori<-HORIZONS
    row.names(hori)<-NULL
    hori<-as.numeric(as.matrix(hori))
    VIEWF<-VIEWF_all

    # setting up for temperature correction using lapse rate given difference between 9sec DEM value and 0.05 deg value
    delta_elev = USADEM - ALTITUDES
    adiab_corr_max <- delta_elev * lapse_max
    adiab_corr_min <- delta_elev * lapse_min
    if(save != 2){
      if(opendap == 1){
        cat("extracting weather data \n")
        dates<-Sys.time()-60*60*24
        curyear<-as.numeric(format(dates,"%Y"))
        days <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900"), as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900"), by = 'days')
        alldays <- seq(as.POSIXct("01/01/1900", format = "%d/%m/%Y", origin = "01/01/1900"), Sys.time()-60*60*24, by = 'days')
        startday <- which(as.character(format(alldays, "%d/%m/%Y")) == format(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900"), "%d/%m/%Y"))
        endday <- which(as.character(format(alldays, "%d/%m/%Y")) == format(as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900"), "%d/%m/%Y"))
        countday <- endday-startday+1
        baseurl <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_met_"
        cat(paste0("reading weather input for ", dstart, " to ", dfinish, " \n"))
        cat(paste0("tmin weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "tmmn_1979_CurrentYear_CONUS.nc"))
        day <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "day", unpack = TRUE)))
        lon <- RNetCDF::var.get.nc(nc, "lon", unpack = TRUE)
        lat <- RNetCDF::var.get.nc(nc, "lat", unpack = TRUE)
        flat <- match(abs(lat-x[2])<1/48,1)
        latindex <- which(flat %in% 1)
        flon <- match(abs(lon-x[1])<1/48,1)
        lonindex <- which(flon %in% 1)
        start <- c(lonindex, latindex, which(day == startday) - 1)
        count <- c(1, 1, countday)
        Tmin <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_minimum_temperature", start = start, count = count, unpack = TRUE))) - 273.15
        RNetCDF::close.nc(nc)
        cat(paste0("tmax weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "tmmx_1979_CurrentYear_CONUS.nc"))
        Tmax <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_maximum_temperature", start = start, count = count, unpack = TRUE))) - 273.15
        RNetCDF::close.nc(nc)
        cat(paste0("rhmin weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "rmin_1979_CurrentYear_CONUS.nc"))
        rhmin <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_minimum_relative_humidity", start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
        cat(paste0("rhmax weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "rmax_1979_CurrentYear_CONUS.nc"))
        rhmax <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_maximum_relative_humidity", start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
        cat(paste0("rain weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "pr_1979_CurrentYear_CONUS.nc"))
        Rain <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "precipitation_amount", start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
        cat(paste0("solar weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "srad_1979_CurrentYear_CONUS.nc"))
        solar <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_mean_shortwave_radiation_at_surface", start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
        cat(paste0("wind weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "vs_1979_CurrentYear_CONUS.nc"))
        Wind <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_mean_wind_speed", start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
      }else{
        cat("extracting weather data \n")
        nc <- RNetCDF::open.nc(paste(spatial, "/tmmx_", yearlist[1], ".nc",
          sep = ""))
        lon <- matrix(RNetCDF::var.get.nc(nc, "lon", unpack = TRUE))
        lat <- matrix(RNetCDF::var.get.nc(nc, "lat", unpack = TRUE))
        lon_1 <- as.numeric(longlat[1])
        lat_1 <- as.numeric(longlat[2])
        dist1 <- abs(lon - lon_1)
        index1 <- which.min(dist1)
        dist2 <- abs(lat - lat_1)
        index2 <- which.min(dist2)
        start <- c(index2, index1, 1)
        count <- c(1, 1, NA)
        for (j in 1:nyears) {
          if (j == 1) {
            cat(paste("reading weather input for ", yearlist[j],
              " \n", sep = ""))
            nc <- RNetCDF::open.nc(paste(spatial, "/tmmn_", yearlist[j],
              ".nc", sep = ""))
            tmin <- as.numeric(RNetCDF::var.get.nc(nc, variable = "air_temperature",
              start = start, count = count, unpack = TRUE))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/tmmx_", yearlist[j],
              ".nc", sep = ""))
            tmax <- as.numeric(RNetCDF::var.get.nc(nc, variable = "air_temperature",
              start = start, count = count, unpack = TRUE))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/rmin_", yearlist[j],
              ".nc", sep = ""))
            rhmin <- as.numeric(RNetCDF::var.get.nc(nc, variable = "relative_humidity",
              start = start, count = count, unpack = TRUE))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/rmax_", yearlist[j],
              ".nc", sep = ""))
            rhmax <- as.numeric(RNetCDF::var.get.nc(nc, variable = "relative_humidity",
              start = start, count = count, unpack = TRUE))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/pr_", yearlist[j],
              ".nc", sep = ""))
            Rain <- as.numeric(RNetCDF::var.get.nc(nc, variable = "precipitation_amount",
              start = start, count = count, unpack = TRUE))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/srad_", yearlist[j],
              ".nc", sep = ""))
            solar <- as.numeric(RNetCDF::var.get.nc(nc, variable = "surface_downwelling_shortwave_flux_in_air",
              start = start, count = count, unpack = TRUE))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/vs_", yearlist[j],
              ".nc", sep = ""))
            Wind <- as.numeric(RNetCDF::var.get.nc(nc, variable = "wind_speed",
              start = start, count = count, unpack = TRUE))
            RNetCDF::close.nc(nc)
            Tmax <- tmax - 273.15
            Tmin <- tmin - 273.15
          } else {
            cat(paste("reading weather input for ", yearlist[j],
              " \n", sep = ""))
            nc <- RNetCDF::open.nc(paste(spatial, "/tmmn_", yearlist[j],
              ".nc", sep = ""))
            tmin <- as.numeric(RNetCDF::var.get.nc(nc, variable = "air_temperature",
              start = start, count = count, unpack = TRUE))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/tmmx_", yearlist[j],
              ".nc", sep = ""))
            tmax <- as.numeric(RNetCDF::var.get.nc(nc, variable = "air_temperature",
              start = start, count = count, unpack = TRUE))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/rmin_", yearlist[j],
              ".nc", sep = ""))
            rhmin <- c(rhmin, as.numeric(RNetCDF::var.get.nc(nc, variable = "relative_humidity",
              start = start, count = count, unpack = TRUE)))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/rmax_", yearlist[j],
              ".nc", sep = ""))
            rhmax <- c(rhmax, as.numeric(RNetCDF::var.get.nc(nc, variable = "relative_humidity",
              start = start, count = count, unpack = TRUE)))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/pr_", yearlist[j],
              ".nc", sep = ""))
            Rain <- c(Rain, as.numeric(RNetCDF::var.get.nc(nc, variable = "precipitation_amount",
              start = start, count = count, unpack = TRUE)))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/srad_", yearlist[j],
              ".nc", sep = ""))
            solar <- c(solar, as.numeric(RNetCDF::var.get.nc(nc, variable = "surface_downwelling_shortwave_flux_in_air",
              start = start, count = count, unpack = TRUE)))
            RNetCDF::close.nc(nc)
            nc <- RNetCDF::open.nc(paste(spatial, "/vs_", yearlist[j],
              ".nc", sep = ""))
            Wind <- c(Wind, as.numeric(RNetCDF::var.get.nc(nc, variable = "wind_speed",
              start = start, count = count, unpack = TRUE)))
            RNetCDF::close.nc(nc)
            Tmax <- c(Tmax, tmax - 273.15)
            Tmin <- c(Tmin, tmin - 273.15)
          }
        }
      }
      # compute clear sky solar for the site of interest, for cloud cover computation below
      cat("running micro_global to get clear sky solar \n")
      micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, timeinterval = 365, solonly = 1)
      clearskyrad <- micro_clearsky$metout[,c(1, 13)]
      clearsky_mean1 <- aggregate(clearskyrad[,2], by = list(clearskyrad[,1]), FUN = mean)[,2]
      leapyears<-seq(1900,2100,4)
      for(j in 1:nyears){
        if(yearlist[j]%in%leapyears){# add day for leap year if needed
          clearsky_mean<-c(clearsky_mean1[1:59],clearsky_mean1[59],clearsky_mean1[60:365])
        }else{
          clearsky_mean <- clearsky_mean1
        }
        if(j == 1){
          allclearsky <- clearsky_mean
        }else{
          allclearsky <- c(allclearsky, clearsky_mean)
        }
      }
      if(opendap == 1 & save != 2){ # truncating if less than whole years requested via opendap
        cut <- as.numeric(days[1] - as.POSIXct(paste0('01/01/', ystart), format = "%d/%m/%Y") + 1)
        allclearsky <- allclearsky[cut:(cut+countday-1)]
      }
      cloud <- (1 - solar / allclearsky) * 100
      cloud[cloud<0]<-0
      cloud[cloud>100]<-100
      if(clearsky == 1){
        cloud <- cloud * 0
      }
      CCMAXX<-as.numeric(cloud)
      CCMINN<-CCMAXX
      CCMINN<-CCMINN*0.5
      CCMAXX<-CCMAXX*2
      CCMINN[CCMINN>100]<-100
      CCMAXX[CCMAXX>100]<-100
      Wind[Wind==0]<-0.1
      if(save == 1){
        cat("saving met data for later \n")
        save(CCMAXX, file = 'CCMAXX.Rda')
        save(CCMINN, file = 'CCMINN.Rda')
        save(Wind, file = 'Wind.Rda')
        save(Tmax, file = 'Tmax.Rda')
        save(Tmin, file = 'Tmin.Rda')
        save(rhmax, file = 'rhmax.Rda')
        save(rhmin, file = 'rhmin.Rda')
        save(Rain, file = 'Rain.Rda')
      }
    }else{
      cat("loading met data from previous run \n")
      load('CCMAXX.Rda')
      load('CCMINN.Rda')
      load('Wind.Rda')
      load('Tmax.Rda')
      load('Tmin.Rda')
      load('rhmax.Rda')
      load('rhmin.Rda')
      load('Rain.Rda')
    }

    ndays<-length(Tmax)
    doynum<-ndays
    leapyears<-seq(1900,2100,4)
    for(k in 1:nyears){
      if(k==1){
        cyear<-ystart
      }else{
        cyear<-cyear+1
      }
      if(cyear %in% leapyears){
        dinyear <- 366
      }else{
        dinyear <- 365
      }
      if(k==1){
        doy <- seq(1,dinyear)
      }else{
        doy <- c(doy, seq(1, dinyear))
      }
    }
    if(opendap == 1 & save != 2){ # could be less than whole years
      doy <- doy[cut:(cut+countday-1)]
    }
    ida<-ndays
    idayst <- 1

    ndays<-length(Tmin)
    maxshades=rep(maxshade,ndays)
    minshades=rep(minshade,ndays)
    shademax<-maxshades

    if(is.na(ALTITUDES)!=TRUE){

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
      }else{ # use a suitable one for Australia (same as around Adelaide/Melbourne)
        TAI<-c(0.0670358341290886,0.0662612704779235,0.065497075238002,0.0647431301168489,0.0639993178022531,0.0632655219571553,0.0625416272145492,0.0611230843885423,0.0597427855962549,0.0583998423063099,0.0570933810229656,0.0558225431259535,0.0545864847111214,0.0533843764318805,0.0522154033414562,0.0499736739981675,0.047855059159556,0.0458535417401334,0.0439633201842001,0.0421788036108921,0.0404946070106968,0.0389055464934382,0.0374066345877315,0.0359930755919066,0.0346602609764008,0.0334037648376212,0.0322193394032758,0.0311029105891739,0.0300505736074963,0.0290585886265337,0.0281233764818952,0.0272415144391857,0.0264097320081524,0.0256249068083005,0.0248840604859789,0.0241843546829336,0.0235230870563317,0.0228976873502544,0.0223057135186581,0.0217448478998064,0.0212128934421699,0.0207077699817964,0.0202275105711489,0.0197702578594144,0.0193342605242809,0.0189178697551836,0.0177713140039894,0.0174187914242432,0.0170790495503944,0.0167509836728154,0.0164335684174899,0.0161258546410128,0.0158269663770596,0.0155360978343254,0.0152525104459325,0.0149755299703076,0.0147045436435285,0.0144389973831391,0.0141783930434343,0.0134220329447663,0.0131772403830191,0.0129356456025128,0.0126970313213065,0.0124612184223418,0.0122280636204822,0.01199745718102,0.0115436048739351,0.0110993711778668,0.0108808815754663,0.0106648652077878,0.0104513876347606,0.0102405315676965,0.00982708969547694,0.00962473896278535,0.00903679230300494,0.00884767454432418,0.0083031278398166,0.00796072474935954,0.00755817587626185,0.00718610751850881,0.00704629977586921,0.00684663903049612,0.00654155580333479,0.00642947339729728,0.00627223096874308,0.00603955966866779,0.00580920937536261,0.00568506186880564,0.00563167068287251,0.00556222005081865,0.00550522989971023,0.00547395763028062,0.0054478983436216,0.00541823364504573,0.00539532163908382,0.00539239864119488,0.00541690124712384,0.00551525885358836,0.00564825853509463,0.00577220185074264,0.00584222986640171,0.00581645238345584,0.00566088137411449,0.00535516862329704,0.00489914757707667,0.00432017939770409,0.0036813032251836,0.00309019064543606,0.00270890436501562,0.00276446109239711,0.00356019862584603)
      } #end check if running gads

      if(adiab_cor==1){
        TMAXX<-as.matrix(Tmax+adiab_corr_max)
        TMINN<-as.matrix(Tmin+adiab_corr_min)
      }else{
        TMAXX<-as.matrix(Tmax)
        TMINN<-as.matrix(Tmin)
      }
      if(warm != 0){
        # impose uniform temperature change
        TMAXX<-TMAXX+seq(0, ndays-1)/(ndays-1)*warm
        TMINN<-TMINN+seq(0, ndays-1)/(ndays-1)*warm
      }
      RAINFALL<-Rain+rainoff

      # correct for potential change in RH with elevation-corrected Tair
      es <- WETAIR(db = TMAXX, rh = 100)$esat
      e <- WETAIR(db = Tmax, rh = rhmin)$e
      RHMINN <- (e / es) * 100
      RHMINN[RHMINN>100]<-100
      RHMINN[RHMINN<0]<-0.01
      es <- WETAIR(db = TMINN, rh = 100)$esat
      e <- WETAIR(db = Tmin, rh = rhmax)$e
      RHMAXX <- (e / es) * 100
      RHMAXX[RHMAXX>100]<-100
      RHMAXX[RHMAXX<0]<-0.01

      ALLMINTEMPS<-TMINN
      ALLMAXTEMPS<-TMAXX
      ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)

      WNMAXX <- Wind * windfac
      WNMINN <- Wind * windfac

      MAXSHADES<-maxshades
      MINSHADES<-minshades

      REFLS <- rep(REFL, ndays)
      PCTWET <- rep(PCTWET, ndays)
      soilwet<-RAINFALL
      soilwet[soilwet<=rainwet] = 0
      soilwet[soilwet>0] = 90
      PCTWET<-pmax(soilwet,PCTWET)

      Intrvls<-rep(0,ndays)
      Intrvls[1] <- 1 # user-supplied last day-of-year in each time interval sequence
      Numtyps <- 10 # number of substrate types
      Numint <- 1  # number of time intervals
      Nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
      Nodes[1:10,] <- c(1:10) # deepest nodes for each substrate type
      ALREF <- abs(trunc(x[1]))

      HEMIS <- ifelse(x[2]<0, 2, 1)
      ALAT <- abs(trunc(x[2]))
      AMINUT <- (abs(x[2])-ALAT)*60
      ALONG <- abs(trunc(x[1]))
      ALMINT <- (abs(x[1])-ALONG)*60
      if(adiab_cor==1){
        ALTT<-ALTITUDES
      }else{
        ALTT<-UKDEM
      }
      SLOPE<-SLOPES
      AZMUTH<-AZMUTHS

      avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
      soilinit<-rep(avetemp,20)
      tannul<-mean(unlist(ALLTEMPS))

      if(nyears==1){
        avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
        tannulrun<-rep(avetemp,ndays)
      }else{
        if(nrow(TMAXX)==1){
          avetemp<-colMeans(cbind(TMAXX, TMINN), na.rm=TRUE)
        }else{
          avetemp<-rowMeans(cbind(TMAXX, TMINN), na.rm=TRUE)
        }
        if(length(TMAXX)<365){
          tannulrun<-rep((sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2),length(TMAXX))
        }else{
          tannulrun<-raster::movingFun(avetemp,n=365,fun=mean,type='to')
          yearone<-rep((sum(TMAXX[1:365])+sum(TMINN[1:365]))/(365*2),365)
          tannulrun[1:365]<-yearone
        }
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
      #Smooth, level, grass-covered   0.15
      #Row crops 	0.2
      #Low bushes with a few trees 	0.2
      #Heavy trees 	0.25
      #Several buildings 	0.25
      #Hilly, mountainous terrain 	0.25
      WNMAXX<-WNMAXX*(2/10)^0.15
      WNMINN<-WNMINN*(2/10)^0.15

      SLES<-matrix(nrow = ndays, data = 0)
      SLES<-SLES+SLE

      moists2<-matrix(nrow=10, ncol = ndays, data=0)
      moists2[1,ndays]<-0.2
      moists<-moists2

      if(runmoist==1){
        moists2<-matrix(nrow=10, ncol = ndays, data=0) # set up an empty vector for soil moisture values through time
        moists2[1:10,]<-SoilMoist_Init
        moists<-moists2
      }
      soilprops<-matrix(data = 0, nrow = 10, ncol = 5)

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

      ALTT<-as.numeric(ALTT)
      ALREF<-as.numeric(ALREF)
      ALMINT<-as.numeric(ALMINT)
      ALONG<-as.numeric(ALONG)
      AMINUT<-as.numeric(AMINUT)
      ALAT<-as.numeric(ALAT)

      # microclimate input parameters list
      microinput<-c(ndays,RUF,ERR,Usrhyt,Refhyt,Numtyps,Z01,Z02,ZH1,ZH2,idayst,ida,HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,CMH2O,microdaily,tannul,EC,VIEWF,snowtemp,snowdens,snowmelt,undercatch,rainmult,runshade,runmoist,maxpool,evenrain,snowmodel,rainmelt,writecsv,densfun,hourly,rainhourly,lamb,IUV,RW,PC,RL,SP,R1,IM,MAXCOUNT,IR,message,fail,snowcond,intercept,grasshade,solonly,ZH,D0)

      # hourly option set to 0, so make empty vectors
      if(hourly==0){
        TAIRhr=rep(0,24*ndays)
        RHhr=rep(0,24*ndays)
        WNhr=rep(0,24*ndays)
        CLDhr=rep(0,24*ndays)
        SOLRhr=rep(0,24*ndays)
        ZENhr=rep(-1,24*ndays)
        IRDhr=rep(-1,24*ndays)
      }
      if(rainhourly==0){
        RAINhr=rep(0,24*ndays)
      }else{
        RAINhr = rainhour
      }

      doy1=matrix(data = 0., nrow = ndays, ncol = 1)
      SLES1=matrix(data = 0., nrow = ndays, ncol = 1)
      MAXSHADES1=matrix(data = 0., nrow = ndays, ncol = 1)
      MINSHADES1=matrix(data = 0., nrow = ndays, ncol = 1)
      TMAXX1=matrix(data = 0., nrow = ndays, ncol = 1)
      TMINN1=matrix(data = 0., nrow = ndays, ncol = 1)
      CCMAXX1=matrix(data = 0., nrow = ndays, ncol = 1)
      CCMINN1=matrix(data = 0., nrow = ndays, ncol = 1)
      RHMAXX1=matrix(data = 0., nrow = ndays, ncol = 1)
      RHMINN1=matrix(data = 0., nrow = ndays, ncol = 1)
      WNMAXX1=matrix(data = 0., nrow = ndays, ncol = 1)
      WNMINN1=matrix(data = 0., nrow = ndays, ncol = 1)
      REFLS1=matrix(data = 0., nrow = ndays, ncol = 1)
      PCTWET1=matrix(data = 0., nrow = ndays, ncol = 1)
      RAINFALL1=matrix(data = 0, nrow = ndays, ncol = 1)
      tannul1=matrix(data = 0, nrow = ndays, ncol = 1)
      moists1=matrix(data = 0., nrow = 10, ncol = ndays)
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
        tides<-matrix(data = 0, nrow = 24*ndays, ncol = 3) # make an empty matrix
      }
      # all microclimate data input list - all these variables are expected by the input argument of the fortran micro2014 subroutine
      micro<-list(tides=tides,microinput=microinput,doy=doy,SLES=SLES1,DEP=DEP,Nodes=Nodes,MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TIMAXS=TIMAXS,TIMINS=TIMINS,TMAXX=TMAXX1,TMINN=TMINN1,RHMAXX=RHMAXX1,RHMINN=RHMINN1,CCMAXX=CCMAXX1,CCMINN=CCMINN1,WNMAXX=WNMAXX1,WNMINN=WNMINN1,TAIRhr=TAIRhr,RHhr=RHhr,WNhr=WNhr,CLDhr=CLDhr,SOLRhr=SOLRhr,RAINhr=RAINhr,ZENhr=ZENhr,IRDhr=IRDhr,REFLS=REFLS1,PCTWET=PCTWET1,soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists1,RAINFALL=RAINFALL1,tannulrun=tannulrun,PE=PE,KS=KS,BB=BB,BD=BD,DD=DD,L=L,LAI=LAI1)
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
      cat(paste('running microclimate model for',ndays,'days from',dstart,' to ', dfinish, ' at site ',location,'\n'))
      ptm <- proc.time() # Start timing
      microut<-microclimate(micro)
      print(proc.time() - ptm) # Stop the clock

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
      if(max(metout[,1] == 0)){
        cat("ERROR: the model crashed - try a different error tolerance (ERR) or a different spacing in DEP", '\n')
      }
      dates <- seq(as.POSIXct(paste0("01/01/",ystart), format = "%d/%m/%Y", tz = 'America/Los_Angeles'), as.POSIXct(paste0("01/01/",yfinish+1), format = "%d/%m/%Y", tz = 'America/Los_Angeles'), by = 'hours')[1:(length(TMAXX) * 24)]
      dates2 <- seq(as.POSIXct(paste0("01/01/",ystart), format = "%d/%m/%Y", tz = 'America/Los_Angeles'), as.POSIXct(paste0("01/01/",yfinish+1), format = "%d/%m/%Y", tz = 'America/Los_Angeles'), by = 'days')[1:(length(TMAXX))]
      if(lamb == 1){
        drlam<-as.data.frame(microut$drlam) # retrieve direct solar irradiance
        drrlam<-as.data.frame(microut$drrlam) # retrieve direct Rayleigh component solar irradiance
        srlam<-as.data.frame(microut$srlam) # retrieve scattered solar irradiance
        if(snowmodel == 1){
          return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,minshade=minshade,maxshade=maxshade,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlamd,dates=dates,dates2=dates2))
        }else{
          return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,minshade=minshade,maxshade=maxshade,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam,dates=dates,dates2=dates2))
        }
      }else{
        if(snowmodel == 1){
          return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,minshade=minshade,maxshade=maxshade,DEP=DEP,dates=dates,dates2=dates2))
        }else{
          return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,minshade=minshade,maxshade=maxshade,DEP=DEP,dates=dates,dates2=dates2))
        }
      }
    } # end of check for na sites
  } # end error trapping
} # end of micro_usa function
