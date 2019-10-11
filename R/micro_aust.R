#' Australian implementation of the microclimate model.
#'
#' An implementation of the NicheMapR microclimate model that uses the AWAP daily weather database
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param ystart First year to run
#' @param yfinish Last year to run
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
#' @usage micro_aust(loc = c(130.5686, -22.6523), ystart = 1990, yfinish = 1990,
#' REFL = 0.15, slope = 0, aspect = 0, DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minshade = 0, maxshade = 90,
#' Usrhyt = 0.01, ...)
#' @export
#' @details
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes, 0=no}\cr\cr
#' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{solonly}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{microclima}{ = 0, Use microclima and elevatr package to adjust solar radiation for terrain? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{manualshade}{ = 1, Use CSIRO Soil and Landscape Grid of Australia? 1=yes, 0=no}\cr\cr
#' \code{soildata}{ = 0, Extract emissivities from gridded data? 1=yes, 0=no}\cr\cr
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
#' \code{BulkDensity}{ = 1.3, Soil bulk density (Mg/m3)}\cr\cr
#' \code{PCTWET}{ = 0, \% of ground surface area acting as a free water surface}\cr\cr
#' \code{rainwet}{ = 1.5, mm of rainfall causing the ground to be 90\% wet for the day}\cr\cr
#' \code{cap}{ = 1, organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)}\cr\cr
#' \code{CMH2O}{ = 1, Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air conditions; 2.0 = humid, tropical conditions (note this is for the whole atmospheric profile, not just near the ground)}\cr\cr
#' \code{hori}{ = rep(0,24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#' \code{lapse_min}{ = 0.0039 Lapse rate for minimum air temperature (degrees C/m)}\cr\cr
#' \code{lapse_max}{ = 0.0077 Lapse rate for maximum air temperature (degrees C/m)}\cr\cr
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
#' \code{microclima.LAI}{ = 0, leaf area index, used by package microclima for radiation calcs}\cr\cr
#' \code{microclima.LOR}{ = 1, leaf orientation for package microclima radiation calcs}\cr\cr
#'
#' \strong{ Snow mode parameters:}
#'
#' \code{snowmodel}{ = 0, run the snow model 1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (mg/m3), overridden by densfun}\cr\cr
#' \code{densfun}{ = c(0,0,0,0), slope and intercept of model of snow density as a linear function of snowpack age if first two values are nonzero, and following the exponential function of Sturm et al. 2010 J. of Hydromet. 11:1380-1394 if all values are non-zero; if it is c(0,0,0,0) then fixed density used}\cr\cr
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
#' \code{tides}{ = matrix(data = 0., nrow = 24*365*nyears, ncol = 3), matrix for each how of the simulation of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#'
#' \strong{Outputs:}
#'
#' \code{ndays}{ - number of days for which predictions are made}\cr\cr
#' \code{longlat}{ - longitude and latitude for which simulation was run (decimal degrees)}\cr\cr
#' \code{dates}{ - vector of dates (hourly, POSIXct, timezone = GMT+10)}\cr\cr
#' \code{dates2}{ - vector of dates (daily, POSIXct, timezone = GMT+10)}\cr\cr
#' \code{nyears}{ - number of years for which predictions are made}\cr\cr
#' \code{RAINFALL}{ - vector of daily rainfall (mm)}\cr\cr
#' \code{elev}{ - elevation at point of simulation (m)}\cr\cr
#' \code{minshade}{ - minimum shade for simulation (\%)}\cr\cr
#' \code{maxshade}{ - maximum shade for simulation (single value - if time varying, in 'MAXSHADES') (\%)}\cr\cr
#' \code{MINSHADES}{ - vector of minimum shades used (\%)}\cr\cr
#' \code{MAXSHADES}{ - vector of maximum shades used (\%)}\cr\cr
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
#'ystart <- 2014
#'yfinish <- 2015
#'nyears <- yfinish - ystart + 1
#'loc <- c(130.5686, -22.6523) # Nyrripi, Northern Territory, Australia
#'micro <- micro_aust(loc = loc, ystart = ystart, yfinish = yfinish, opendap = 1, elev = 0, runshade = 0) # run the model for the middle of the desert in Australia, using opendap
#'
#' metout <- as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#' soil <- as.data.frame(micro$soil) # soil temperatures, minimum shade
#' soilmoist <- as.data.frame(micro$soilmoist) # soil temperatures, minimum shade
#'
#' # append dates
#' dates <- micro$dates
#' metout <- cbind(dates,metout)
#' soil <- cbind(dates,soil)
#' soilmoist <- cbind(dates, soilmoist)
#' minshade<-micro$minshade
#'
#' # plotting above-ground conditions in minimum shade
#' with(metout,{plot(TALOC ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l",main=paste("air and sky temperature, ",minshade,"% shade",sep=""), ylim = c(-20, 60))})
#' with(metout,{points(TAREF ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l",lty=2,col='blue')})
#' with(metout,{points(TSKYC ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' ,  type = "l",col='light blue',main=paste("sky temperature, ",minshade,"% shade",sep=""))})
#' with(metout,{plot(RHLOC ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l",ylim=c(0,100),main=paste("humidity, ",minshade,"% shade",sep=""))})
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
#'    ,col=i,type = "l",main=paste("soil temperature ",minshade,"% shade",sep=""))
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
#'    ,col=i,type = "l",main=paste("soil moisture ",minshade,"% shade",sep=""))
#'  }else{
#'    points(soilmoist[,i+3]*100~soilmoist[,1],xlab = "Date and Time", ylab = "Soil Moisture
#'     (%)",col=i,type = "l")
#'  }
#' }
micro_aust <- function(
  loc = c(130.5686, -22.6523),
  ystart = 1990,
  yfinish = 1990,
  nyears = 1,
  REFL = 0.15,
  elev = NA,
  slope = 0,
  aspect = 0,
  lapse_max = 0.0077,
  lapse_min = 0.0039,
  DEP = c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
  minshade = 0,
  maxshade = 90,
  MINSHADES = NA,
  MAXSHADES = NA,
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
  manualshade = 1,
  soildata = 0,
  terrain = 0,
  dailywind = 1,
  windfac = 1,
  adiab_cor = 1,
  warm = 0,
  spatial = "W:/",
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
  CMH2O=1,
  hori = rep(0, 24),
  TIMAXS = c(1, 1, 0, 0),
  TIMINS=c(0, 0, 1, 1),
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
  SoilMoist_Init=c(0.1, 0.12, 0.15, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
  L = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000,
  R1 = 0.001,
  RW = 2.5e+10,
  RL = 2e+6,
  PC = -1500,
  SP = 10,
  IM = 1e-06,
  MAXCOUNT = 500,
  LAI = 0.1,
  microclima.LAI = 0,
  LOR = 1,
  snowmodel = 0,
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
  barcoo = "",
  quadrangle = 1,
  rainhourly = 0,
  rainhour = 0,
  uid = "",
  pwd = "",
  lamb = 0,
  IUV = 0,
  microclima = 0,
  soilgrids = 0,
  IR = 0,
  opendap = 1,
  message = 0,
  fail = nyears * 24 * 365,
  snowcond = 0,
  intercept = maxshade / 100 * 0.4,
  grasshade = 0
) {

  # function to assist with interpolated data in leap years
  leapfix <- function(indata, yearlist, mult = 1){
    leapyears<-seq(1900,2060,4)
    for(j in 1:length(yearlist)){
      if(yearlist[j] %in% leapyears){# add day for leap year if needed
        if(mult == 1){
          data<-c(indata[1:59], indata[59], indata[60:365])
        }else{
          data<-c(indata[1:(59*mult)], indata[(58*mult+1):(59*mult)], indata[(59*mult+1):(365*mult)])
        }
      }else{
        data <- indata
      }
      if(j==1){
        alldata <- data
      }else{
        alldata <- c(alldata, data)
      }
    }
    return(alldata)
  }

  if(opendap==0){
    require(RODBC)
  }
  if(opendap == 1){
    if(is.na(elev)){
      require(microclima)
      require(raster)
      cat('downloading DEM via package elevatr \n')
      dem <- microclima::get_dem(lat = loc[2], long = loc[1]) # mercator equal area projection
      elev <- extract(dem, loc)[1]
    }
    require(RNetCDF)
    ALTITUDES <- elev
    dbrow <- 1
  }
  errors<-0

  # error trapping - originally inside the Fortran code, but now checking before executing Fortran
  if(opendap == 1 & (ystart < 1990 | ystart > 2015)){
    message("currently no data on the NCI for this time window \n")
    errors<-1
  }
  if(DEP[2]-DEP[1]>3 | DEP[3]-DEP[2]>3){
    message("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if(DEP[2]-DEP[1]<2){
    cat("warning, nodes might be too close near the surface, try a different spacing if the program is crashing \n")
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
  if(is.na(slope) == FALSE & (slope<0 | slope>90)){
    message("ERROR: Slope value (slope) is out of bounds.
        Please input a value between 0 and 90.", '\n')
    errors<-1
  }
  if(is.na(aspect) == FALSE & (aspect<0 | aspect>365)){
    message("ERROR: Aspect value (aspect) is out of bounds.
        Please input a value between 0 and 365.", '\n')
    errors<-1
  }
  if(is.na(hori[1]) == FALSE & (max(hori)>90 | min(hori)<0)){
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
    message("ERROR: Local height (Usrhyt) is out of bounds.
        Please enter a correct value (0.005 - Refhyt).", '\n')
    errors<-1
  }
  if(CMH2O<0.5 | CMH2O>2){
    message("ERROR: Preciptable water in air column (CMH2O) is out of bounds.
        Please enter a correct value (0.1 - 2cm).", '\n')
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
  # end error trapping

  if(errors==0){ # continue

    ################## time related variables #################################
    nyears <- yfinish-ystart+1
    yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
    tzone <- paste("Etc/GMT+",10,sep="")
    dates <- seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days")
    if(yfinish == as.numeric(format(Sys.time(), "%Y"))){ # cut days down if doing current year (incomplete)
      dates <- dates[dates < Sys.time()]
      dailywind <- 0 # no daily wind for current year
    }
    ndays <- length(dates)
    doys12<-c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month
    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in
    daystart<-1
    if(is.na(MAXSHADES)){
      maxshades <- rep(maxshade,ndays)
    }else{
      maxshades <- MAXSHADES
    }
    if(is.na(MINSHADES)){
      minshades <- rep(minshade,ndays)
    }else{
      minshades <- MINSHADES
    }
    leapyears<-seq(1900,2060,4)
    for(mm in 1:nyears){
      if(mm == 1){
        currenty <- ystart
      }else{
        currenty <- ystart + mm - 1
      }
      if(currenty %in% leapyears){
        dayoy <- seq(1,366)
      }else{
        dayoy <- seq(1,365)
      }
      if(mm == 1){
        doy <- dayoy
      }else{
        doy <- c(doy, dayoy)
      }
    }
    idayst <- 1 # start day
    ida<-ndays # end day
    dates<-Sys.time()-60*60*24
    curyear<-as.numeric(format(dates,"%Y"))
    # location and terrain
    f1 <- paste(spatial,"ausclim_rowids.nc",sep="");
    f2 <- paste(spatial,"ausdem_shift1.tif",sep="");
    f3 <- paste(spatial,"agg_9secdem.nc",sep="");
    f4 <- paste(spatial,"Aust9secDEM.tif",sep="");

    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    cat("running micro_global to get clear sky solar \n")
    if(run.gads == 0){
      TAI <- c(0.0670358341290886,0.0662612704779235,0.065497075238002,0.0647431301168489,0.0639993178022531,0.0632655219571553,0.0625416272145492,0.0611230843885423,0.0597427855962549,0.0583998423063099,0.0570933810229656,0.0558225431259535,0.0545864847111214,0.0533843764318805,0.0522154033414562,0.0499736739981675,0.047855059159556,0.0458535417401334,0.0439633201842001,0.0421788036108921,0.0404946070106968,0.0389055464934382,0.0374066345877315,0.0359930755919066,0.0346602609764008,0.0334037648376212,0.0322193394032758,0.0311029105891739,0.0300505736074963,0.0290585886265337,0.0281233764818952,0.0272415144391857,0.0264097320081524,0.0256249068083005,0.0248840604859789,0.0241843546829336,0.0235230870563317,0.0228976873502544,0.0223057135186581,0.0217448478998064,0.0212128934421699,0.0207077699817964,0.0202275105711489,0.0197702578594144,0.0193342605242809,0.0189178697551836,0.0177713140039894,0.0174187914242432,0.0170790495503944,0.0167509836728154,0.0164335684174899,0.0161258546410128,0.0158269663770596,0.0155360978343254,0.0152525104459325,0.0149755299703076,0.0147045436435285,0.0144389973831391,0.0141783930434343,0.0134220329447663,0.0131772403830191,0.0129356456025128,0.0126970313213065,0.0124612184223418,0.0122280636204822,0.01199745718102,0.0115436048739351,0.0110993711778668,0.0108808815754663,0.0106648652077878,0.0104513876347606,0.0102405315676965,0.00982708969547694,0.00962473896278535,0.00903679230300494,0.00884767454432418,0.0083031278398166,0.00796072474935954,0.00755817587626185,0.00718610751850881,0.00704629977586921,0.00684663903049612,0.00654155580333479,0.00642947339729728,0.00627223096874308,0.00603955966866779,0.00580920937536261,0.00568506186880564,0.00563167068287251,0.00556222005081865,0.00550522989971023,0.00547395763028062,0.0054478983436216,0.00541823364504573,0.00539532163908382,0.00539239864119488,0.00541690124712384,0.00551525885358836,0.00564825853509463,0.00577220185074264,0.00584222986640171,0.00581645238345584,0.00566088137411449,0.00535516862329704,0.00489914757707667,0.00432017939770409,0.0036813032251836,0.00309019064543606,0.00270890436501562,0.00276446109239711,0.00356019862584603)
    }else{
      TAI <- 0
    }
    micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, TAI = TAI, timeinterval = 365, solonly = 1)
    clearskyrad <- micro_clearsky$metout[,c(1, 13)]
    clearskysum <- aggregate(clearskyrad[,2], by = list(clearskyrad[,1]), FUN = sum)[,2]

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

    if(soildata==0){
      soilprop<-cbind(0,0)
      # creating the shade array
      if(is.na(MAXSHADES)){
        MAXSHADES <- rep(0,ndays)+maxshade # daily max shade (%)
      }
      if(is.na(MINSHADES)){
        MINSHADES <- rep(0,ndays)+minshade # daily min shade (%)
      }
    }

    if(soildata==1){
      cat('extracting soil emissivities \n')
      static_soil<-paste(spatial,"static_soil.nc",sep="")
      emissivities<-paste(spatial,"aus_emissivities.nc",sep="")
      # read data in from netcdf file
      static_soil_data<-raster::brick(static_soil)
      static_soil_vars <- raster::extract(static_soil_data,x)
      labels<-c('albedo','FAPAR1','FAPAR2','FAPAR3','FAPAR4','FAPAR5','FAPAR6','FAPAR7','FAPAR8','FAPAR9','FAPAR10','FAPAR11','FAPAR12','volwater_Upper','volwater_lower','thick_upper','thick_lower','code')
      colnames(static_soil_vars)<-labels
      emissivities_data<-raster::brick(emissivities)
      SLES2 <- raster::extract(emissivities_data,x)

      # read in other soil related files for working out lumped soil type and properties
      # such as clay % for getting water potential
      filename<-paste(spatial,"ppfInterpAll.txt",sep="")
      ppf<-as.data.frame(read.table(file = filename, sep = ",", header=TRUE))
      filename<-paste(spatial,"Lumped soil types.txt",sep="")
      lumped.soil<-as.data.frame(read.table(file = filename, sep = ","))
      filename<-paste(spatial,"SoilTypeLUT_725_AWAP.csv",sep="")
      soiltype<-as.data.frame(read.table(file = filename, sep = ","))

      soilcode<-subset(soiltype, soiltype[1]==static_soil_vars[18])
      lumped<-subset(lumped.soil, V4==as.character(soilcode[1,2]))
      soiltype<-lumped[1,6]
      soilprop<-subset(ppf, ppf==soilcode[1,2])
    }else{
      SLES2 <- rep(SLE,ndays)
      if(manualshade==0){
        message("extracting shade data \n")
        static_soil<-paste(spatial,"static_soil.nc",sep="")
        emissivities<-paste(spatial,"aus_emissivities.nc",sep="")
        # read data in from netcdf file
        static_soil_data<-raster::brick(static_soil)
        static_soil_vars <- raster::extract(static_soil_data,x)
        labels<-c('albedo','FAPAR1','FAPAR2','FAPAR3','FAPAR4','FAPAR5','FAPAR6','FAPAR7','FAPAR8','FAPAR9','FAPAR10','FAPAR11','FAPAR12','volwater_Upper','volwater_lower','thick_upper','thick_lower','code')
        colnames(static_soil_vars)<-labels
      }
    }
    if(terrain==1){
      message("extracting terrain data \n")
      e<-extent(x[1]-0.05,x[1]+0.05,x[2]-0.05,x[2]+0.05)
      for(i in 1:24){
        horifile<-paste(spatial,'horizon',i,'.tif',sep="")
        horiz<-raster::crop(raster::raster(horifile),e)
        if(i==1){
          horizons_data<-horiz
        }else{
          horizons_data<-raster::stack(horizons_data,horiz)
        }
      }
      HORIZONS <- t(raster::extract(horizons_data,x))
      elev1<-crop(raster::raster(paste(spatial,'elev.tif',sep="")),e)
      slope1<-crop(raster::raster(paste(spatial,'slope.tif',sep="")),e)
      aspect1<-crop(raster::raster(paste(spatial,'aspect.tif',sep="")),e)
      elevslpasp<-raster::stack(elev1,slope1,aspect1)
      ELEVSLPASP <- raster::extract(elevslpasp,x)
      ELEVSLPASP<-as.matrix((ifelse(is.na(ELEVSLPASP),0,ELEVSLPASP)))
      ALTITUDES <- ELEVSLPASP[,1]
      SLOPES <- ELEVSLPASP[,2]
      AZMUTHS <- ELEVSLPASP[,3]
      # the horizons have been arranged so that they go from 0 degrees azimuth (north) clockwise - r.horizon starts
      # in the east and goes counter clockwise!
      HORIZONS <- (ifelse(is.na(HORIZONS),0,HORIZONS))/10 # get rid of na and get back to floating point
      HORIZONS <- data.frame(HORIZONS)
      VIEWF_all <- 1-rowSums(sin(t(HORIZONS)*pi/180))/length(t(HORIZONS)) # convert horizon angles to radians and calc view factor(s)
      r1 <- raster::raster(f1)
      r2 <- raster::raster(f2)
      r3 <- raster::raster(f3)
      dbrow <- raster::extract(r1, x)
      AUSDEM <- raster::extract(r2, x)
      AGG <- raster::extract(r3, x)
    }else{
      if(opendap == 0){
        r1 <- raster::raster(f1)
        r2 <- raster::raster(f2)
        r3 <- raster::raster(f3)
        r4 <- raster::raster(f4)
        dbrow <- raster::extract(r1, x)
        AUSDEM <- raster::extract(r2, x)
        AGG <- raster::extract(r3, x)
        if(is.na(elev) == FALSE){ # check if user-specified elevation
          ALTITUDES <- elev # use user-specified elevation
        }else{
          ALTITUDES <- raster::extract(r4, x) # get elevation from fine res DEM
        }
      }
      HORIZONS <- hori
      HORIZONS <- data.frame(HORIZONS)
      VIEWF_all <- 1-sum(sin(as.data.frame(hori)*pi/180))/length(hori) # convert horizon angles to radians and calc view factor(s)
      SLOPES<-rep(slope,length(x[,1]))
      AZMUTHS<-rep(aspect,length(x[,1]))
    }
    hori<-HORIZONS
    row.names(hori)<-NULL
    hori<-as.numeric(as.matrix(hori))

    if(soildata==1){
      VIEWF<-VIEWF_all
      SLES<-SLES2
    }else{
      VIEWF<-VIEWF_all
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
    # setting up for temperature correction using lapse rate given difference between 9sec DEM value and 0.05 deg value
    if(opendap == 0){
      if(AUSDEM==-9999 | is.na(AUSDEM)=='TRUE'){
        delta_elev = AGG - ALTITUDES
      }else{
        delta_elev = AUSDEM - ALTITUDES
      }
      adiab_corr_max <- delta_elev * lapse_max
      adiab_corr_min <- delta_elev * lapse_min
    }else{
      adiab_corr_max <- 0
      adiab_corr_min <- 0
    }
    if(scenario!=""){
      message("generate climate change scenario", '\n')
      # diff spline function
      getdiff<-function(diffs,grid){
        diff1<-(unlist(diffs[1])+unlist(diffs[12]))/2

        # generate list of days
        for(ys in 1:nyears){
          if(yearstodo[ys]%in%leapyears){
            day<-c(1,15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 366)
          }else{
            day<-c(1,15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 365)
          }
          if(ys==1){
            days2=day
            days=day
          }else{
            if(yearstodo[ys]%in%leapyears){
              days2=c(days2,(day+366*(ys-1)))
            }else{
              days2=c(days2,(day+365*(ys-1)))
            }
            days=c(days,day)
          }
        }

        if(is.na(diffs[1])==TRUE){
          # find the nearest cell with data
          NArem<-grid[[1]]
          NArem<-Which(!is.na(NArem), cells=TRUE)
          dist<-distanceFromPoints(maxTst05[[1]],x)
          distNA<-extract(dist,NArem)
          cellsR<-cbind(distNA,NArem)
          distmin<-which.min(distNA)
          cellrep<-cellsR[distmin,2]
          diffs<-extract(maxTst05,cellrep)
          diff1<-(unlist(diffs[1])+unlist(diffs[12]))/2
        }
        diffs3=rep(c(diff1,diffs,diff1),nyears)
        days_diffs<-data.frame(matrix(NA, nrow = nyears*14, ncol = 3))
        days_diffs[,1]<-days
        days_diffs[,3]<-days2
        days_diffs[,2]<-diffs3
        colnames(days_diffs)<-c("days","diffs","new_day")

        # interpolate monthly differences
        f<-approxfun(x=days_diffs$new_day, y=days_diffs$diffs)
        xx<-seq(1,max(days2),1)
        sp_diff<-f(xx)
        return(sp_diff)
      }
      # Max and Min Air Temps
      ccdir <- "/vlsci/VR0212/mrke"
      load(file = paste0(ccdir, "/Australia Climate Change/",scenario,"/","maxTst05_",scenario,"_",year,".Rda")) #maxTst05
      diffs<-extract(maxTst05,x)
      TMAXX_diff<-getdiff(diffs,maxTst05)
      load(file = paste0(ccdir, "/Australia Climate Change/",scenario,"/","minTst05_",scenario,"_",year,".Rda")) #minTst05
      diffs<-extract(minTst05,x)
      TMINN_diff<-getdiff(diffs,minTst05)
      # RH
      load(file = paste0(ccdir, "/Australia Climate Change/",scenario,"/","RHst05_",scenario,"_",year,".Rda")) #maxTst05
      diffs<-extract(RHst05,x)
      RH_diff<-getdiff(diffs,RHst05)
      # wind
      load(file = paste0(ccdir, "/Australia Climate Change/",scenario,"/","PT_VELst05_",scenario,"_",year,".Rda"))
      diffs<-extract(PT_VELst05,x)
      WIND_diff<-getdiff(diffs,PT_VELst05)
      # SOLAR/CLOUD COVER
      load(file = paste0("c:/Spatial_Data/Australia Climate Change/",scenario,"/","SOLCst05_",scenario,"_",year,".Rda"))
      diffs<-extract(SOLCst05,x)
      SOLAR_diff<-getdiff(diffs,SOLCst05)
    }

    if(opendap == 1){
      message("extracting climate data via opendap - note that there is no wind speed data, so the daily range is assumed to be from 0.5 to 2 m/s \n")
      monstart <- c("0101", "0201", "0301", "0401", "0501", "0601", "0701", "0801", "0901", "1001", "1101", "1201")
      monfinish <- c("0131.nc","0228.nc","0331.nc","0430.nc","0531.nc","0630.nc","0731.nc","0831.nc","0930.nc","1031.nc","1130.nc","1231.nc")
      monfinish2 <- c("0131.nc","0229.nc","0331.nc","0430.nc","0531.nc","0630.nc","0731.nc","0831.nc","0930.nc","1031.nc","1130.nc","1231.nc")
      message("extracting climate data", '\n')
      baseurl<- "http://rs-data1-mel.csiro.au/thredds/dodsC/bawap/"
      get_AWAP <- function(year, var){
        for(i in 1:length(years)){
          year <- years[i]
          cat(paste("reading weather input for variable ", var," for year ", year, " \n", sep = ""))
          for(mm in 1:12){
            mstart <- monstart[mm]
            if(year%in%leapyears){
              mfinish <- monfinish2[mm]
            }else{
              mfinish <- monfinish[mm]
            }
            nc <- RNetCDF::open.nc(paste0(baseurl, var,"/day/",year,"/bom-",var,"_day-",year,mstart,"-",year,mfinish))
            lon <- RNetCDF::var.get.nc(nc, "longitude", unpack = TRUE)
            lat <- RNetCDF::var.get.nc(nc, "latitude", unpack = TRUE)
            flat=match(abs(lat-x[2])<1/40,1)
            latindex=which(flat %in% 1)[1]
            flon=match(abs(lon-x[1])<1/40,1)
            lonindex=which(flon %in% 1)[1]
            start <- c(lonindex, latindex, 1)
            count <- c(1, 1, NA)
            data_1 <- as.numeric(RNetCDF::var.get.nc(nc, variable = paste0(var, "_day"),
                                                     start = start, count, unpack = TRUE))
            if(mm == 1 & i == 1){
              data <- data_1
            }else{
              data <- c(data, data_1)
            }
            RNetCDF::close.nc(nc)
          }
        }
        return(data)
      }
      years <- yearlist
      sol <- zoo::na.approx(get_AWAP(yearlist, "rad"))
      tmax <- get_AWAP(yearlist, "tmax")
      tmin <- get_AWAP(yearlist, "tmin")
      vph09 <- get_AWAP(yearlist, "vph09")
      vph15 <- get_AWAP(yearlist, "vph15")
      Rain <- get_AWAP(yearlist, "rain")

      allclearsky <- leapfix(clearskysum, yearlist)
      allclearsky <- allclearsky[1:ndays]
      # convert from W/d to MJ/d
      allclearsky <- allclearsky * 3600 / 1e6
      cloud <- (1-sol/allclearsky) * 100
      cloud[cloud<0]<-0
      cloud[cloud>100]<-100
      CCMAXX<-as.numeric(cloud)
      CCMINN<-CCMAXX
      CCMINN<-CCMINN*0.5
      CCMAXX<-CCMAXX*2
      CCMINN[CCMINN>100]<-100
      CCMAXX[CCMAXX>100]<-100
      TMAXX <- tmax
      TMINN <- tmin
      VAPRES_max<-apply(cbind(vph09, vph15), FUN = max, MARGIN = 1)*100 # convert from hectopascals to pascals
      VAPRES_min<-apply(cbind(vph09, vph15), FUN = min, MARGIN = 1)*100 # convert from hectopascals to pascals
      es <- WETAIR(db = TMAXX, rh = 100)$esat
      RHMINN <- (VAPRES_min / es) * 100
      RHMINN[RHMINN>100]<-100
      RHMINN[RHMINN<0]<-0.01
      es <- WETAIR(db = TMINN, rh = 100)$esat
      RHMAXX <- (VAPRES_max / es) * 100
      RHMAXX[RHMAXX>100]<-100
      RHMAXX[RHMAXX<0]<-0.01
      WNMAXX <- rep(2, ndays)
      WNMINN <- rep(0.5, ndays)
      RAINFALL <- Rain
    }

    # connect to server
    if(opendap == 0){
      channel2 <- RODBC::odbcConnect("ausclim_predecol", uid = uid, pwd = pwd)
      channel <- RODBC::odbcConnect("AWAPDaily", uid = uid, pwd = pwd)
      # preliminary test for incomplete year, if simulation includes the present year
      for(j in 1:nyears){ # start loop through years
        yeartodo<-yearlist[j]
        lat1<-x[2]-0.024
        lat2<-x[2]+0.025
        lon1<-x[1]-0.024
        lon2<-x[1]+0.025
        query<-paste("SELECT a.latitude, a.longitude, b.* FROM [AWAPDaily].[dbo].[latlon] as a
                   , [AWAPDaily].[dbo].[",yeartodo,"] as b where (a.id = b.id) and
                   (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
                   order by b.day",sep="")
        output<- sqlQuery(channel,query)
        if(yearlist[j]<1971){
          output$vpr<-output$tmin/output$tmin-1
        }
        if(yearlist[j]>1989){
          output$sol<-as.numeric(as.character(output$sol))
        }else{
          output$sol<-output$tmin/output$tmin-1
        }
        if(j==1){
          results<-output
        }else{
          results<-rbind(results,output)
        }
      }
    }
    if(is.na(MAXSHADES[1])){
      maxshades=rep(maxshade,ndays)
    }
    if(is.na(MINSHADES[1])){
      minshades=rep(minshade,ndays)
    }
    doys<-seq(daystart,ndays,1)
    leapyears<-seq(1900,2060,4)
    for(mm in 1:nyears){
      if(mm == 1){
        currenty <- ystart
      }else{
        currenty <- ystart + mm
      }
      if(currenty %in% leapyears){
        dayoy <- seq(1,366)
      }else{
        dayoy <- seq(1,365)
      }
      if(mm == 1){
        doy <- dayoy
      }else{
        doy <- c(doy, dayoy)
      }
    }
    doy<-leapfix(doy, yearlist)
    doy<-doy[1:ndays]
    ida<-ndays
    idayst <- 1 # start month
    # end preliminary test for incomplete year, if simulation includes the present year


    if((soildata==1 & nrow(soilprop)>0) | soildata == 0){
      if(soildata==1){
        # get static soil data into arrays
        REFL <- static_soil_vars[,1]  # albedo/reflectances
        maxshades <- static_soil_vars[,2:13] # assuming FAPAR represents shade
        shademax<-maxshades
      }else{
        if(manualshade==0){
          maxshades <- static_soil_vars[,2:13] # assuming FAPAR represents shade
        }
        shademax<-maxshades
      }
      if((is.na(dbrow)!=TRUE & is.na(ALTITUDES)!=TRUE) | opendap == 1){
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

        if(opendap == 0){
          for(j in 1:nyears){ # start loop through years
            yeartodo<-yearlist[j]
            lat1<-x[2]-0.024
            lat2<-x[2]+0.025
            lon1<-x[1]-0.024
            lon2<-x[1]+0.025
            query<-paste("SELECT a.latitude, a.longitude, b.*
                       FROM [AWAPDaily].[dbo].[latlon] as a
                       , [AWAPDaily].[dbo].[",yeartodo,"] as b
                       where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
                       order by b.day",sep="")
            output<- sqlQuery(channel,query)
            if(yearlist[j]<1971){
              output$vpr<-output$tmin/output$tmin-1
            }
            if(yearlist[j]>1989){
              output$sol<-as.numeric(as.character(output$sol))
            }else{
              output$sol<-output$tmin/output$tmin-1
            }
            if(j==1){
              results<-output
            }else{
              results<-rbind(results,output)
            }
          }
          if(dailywind==1){
            channel <- RODBC::odbcConnect("dailywind", uid = uid, pwd = pwd)
            if(min(yearlist)<1975){
              # get mean of 1975-1984
              for(j in 1:10){ # start loop through years
                yeartodo<-1974+j
                lat1<-x[2]-0.024
                lat2<-x[2]+0.025
                lon1<-x[1]-0.024
                lon2<-x[1]+0.025
                query<-paste("SELECT a.latitude, a.longitude, b.*
                         FROM [dailywind].[dbo].[latlon] as a
                         , [dailywind].[dbo].[",yeartodo,"] as b
                         where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
                         order by b.day",sep="")
                output<- sqlQuery(channel,query)
                if(j==1){
                  dwindmean<-output
                }else{
                  dwindmean<-cbind(dwindmean,output[,5])
                }
              }
              dwindmean<-cbind(dwindmean[,1:4],rowMeans(dwindmean[,5:14]))
              colnames(dwindmean)[5]<-'wind'
            }
            for(j in 1:nyears){ # start loop through years
              yeartodo<-yearlist[j]
              if(yeartodo<1975){
                output<-dwindmean
              }else{
                lat1<-x[2]-0.024
                lat2<-x[2]+0.025
                lon1<-x[1]-0.024
                lon2<-x[1]+0.025
                query<-paste("SELECT a.latitude, a.longitude, b.*
                         FROM [dailywind].[dbo].[latlon] as a
                         , [dailywind].[dbo].[",yeartodo,"] as b
                         where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
                         order by b.day",sep="")
                output<- sqlQuery(channel,query)
              }
              if(j==1){
                dwind<-output
              }else{
                dwind<-rbind(dwind,output)
              }
            }
            dwind<-dwind$wind/15.875
          }
        }
        if(opendap == 0){
          if(adiab_cor==1){
            TMAXX.orig <- results$tmax
            TMINN.orig <- results$tmin
            TMAXX<-as.matrix(results$tmax+adiab_corr_max)
            TMINN<-as.matrix(results$tmin+adiab_corr_min)

          }else{
            TMAXX<-as.matrix(results$tmax)
            TMINN<-as.matrix(results$tmin)
          }
          if(scenario!=""){
            TMAXX=TMAXX+TMAXX_diff
            TMINN=TMINN+TMINN_diff
          }
          RAINFALL<-results$rr
          output_AWAPDaily<-results
        }else{
          if(adiab_cor==1){
            TMAXX.orig <- TMAXX
            TMINN.orig <- TMINN
            TMAXX<-as.matrix(TMAXX+adiab_corr_max)
            TMINN<-as.matrix(TMINN+adiab_corr_min)
          }
        }
        if(scenario!=""){
          # first work out for each site the new predicted rainfall amount for each month - use this to adjust for fact that will underestimate chcange
          # using proportion becasue 0 x % is still 0
          # add columns with days, months and years
          RAIN_current<-as.data.frame(RAINFALL)
          dates<-seq(ISOdate(ystart,1,1,tz=paste("Etc/GMT-",10,sep=""))-3600*12, ISOdate((ystart+nyears),1,1,tz=paste("Etc/GMT-",10,sep=""))-3600*13, by="days")
          dates<-subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
          RAINFALL_sum<-aggregate(RAIN_current, by=list(format(dates,"%m-%Y")),FUN=sum)
          dates2<-RAINFALL_sum$Group.1
          RAINFALL_sum<-RAINFALL_sum[order(as.Date(paste("01-",RAINFALL_sum$Group.1,sep=""),"%m-%Y")),2]

          load(file = paste("c:/Spatial_Data/Australia Climate Change/",scenario,"/","RAINst05_",scenario,"_",year,".Rda",sep=""))
          diffs<-rep(extract(RAINst05,x),nyears)

          if(is.na(diffs[1])==TRUE){
            print("no data")
            # find the nearest cell with data
            NArem<-RAINst05[[1]] # don't need to re-do this within bioregion loop
            NArem<-Which(!is.na(NArem), cells=TRUE)
            dist<-distanceFromPoints(RAINst05[[1]],y)
            distNA<-extract(dist,NArem)
            cellsR<-cbind(distNA,NArem)
            distmin<-which.min(distNA)
            cellrep<-cellsR[distmin,2]
            diffs<-rep(extract(RAINst05,cellrep),nyears)
          }

          rainfall_new<-(RAINFALL_sum*diffs)

          rainfall_new[rainfall_new < 0 ]= 0 # get rid of any negative rainfall values

          ## Now extract predicted change in mm
          load(file = paste("c:/Spatial_Data/Australia Climate Change/",scenario,"/","RAINst05_mm_",scenario,"_",year,".Rda",sep=""))
          rainfall_change_mm<-rep(extract(RAINst05_mm,x),nyears)

          #########Now get predicted change in rainfall (could also get this from OzClim or ClimSim layer)#############
          Diff_prop<-rainfall_new/RAINFALL_sum # proportion change
          Diff_prop[Diff_prop=='NaN']= 0
          Diff_prop[Diff_prop=='Inf']= 0 ## If was previously no rainfall and now is rainfall need to alter code so this is simply added

          newRAINFALL=rep(NA,length(RAINFALL))
          for (k in 1:length(RAINFALL)){ # loop through each sites applying monthly % changes
            month<-which(dates2==format(dates[k],"%m-%Y"))
            # Test if predicted rainfall matches up - use rainfall_change_mm layer
            Rain_adj<-rainfall_change_mm[month]

            # test for if proportional change is 0 (because current rainfall is 0)
            #but rainfall predicted to increase
            if(Diff_prop[month]==0 & Rain_adj > 1){ # couldn't get proportion as no current rain days but need to add rain
              print('new rain days needed')
              # previously no rain, randomly select a day and put all rain on it
              listD<-seq(1,length(newRAINFALL),1)
              altD<-sample(listD,1)
              newRAINFALL[altD]<-Rain_adj/30
            }else{
              newRAINFALL[k]<-RAINFALL[k]*Diff_prop[month]
            }
          } # end of loop through each day
          newRAINFALL[newRAINFALL<0.1]<-0
          newRAINFALL[is.na(newRAINFALL)]<-0
          RAINFALL=newRAINFALL
        }
        # cloud cover
        if(opendap==0){
          if(ystart>1989 & sum(results[,9],na.rm=TRUE)>0){ # solar radiation data available
            allclearsky <- leapfix(clearskysum, yearlist)
            allclearsky <- allclearsky[1:ndays]
            # convert from W/d to MJ/d
            allclearsky <- allclearsky * 3600 / 1e6
            if(is.na(output_AWAPDaily[1,9])==TRUE){
              output_AWAPDaily[1,9]=mean(output_AWAPDaily[,9],na.rm=TRUE)
            }
            if(is.na(output_AWAPDaily[ndays,9])==TRUE){
              output_AWAPDaily[nrow(output_AWAPDaily),9]=mean(output_AWAPDaily[,9],na.rm=TRUE)
            }
            solar<-zoo::na.approx(output_AWAPDaily[,9])
            if(scenario!=""){
              solar=solar*SOLAR_diff
            }
            cloud <- (1-solar/allclearsky) * 100
            cloud[cloud<0]<-0
            cloud[cloud>100]<-100
            CCMAXX<-as.numeric(cloud)
            CCMINN<-CCMAXX
          }else{
            datestart1<-"01/01/1990" # day, month, year
            datefinish1<-"31/12/2014" # day, month, year
            datestart1<-strptime(datestart1, "%d/%m/%Y") # convert to date format
            datefinish1<-strptime(datefinish1, "%d/%m/%Y") # convert to date format
            yearstart1<-as.numeric(format(datestart1, "%Y")) # yet year start
            yearfinish1<-as.numeric(format(datefinish1, "%Y")) # yet year finish
            years1<-seq(yearstart1,yearfinish1,1) # get sequence of years to d0
            doystart<-datestart1$yday+1 # get day-of-year at start
            doyfinish<-datefinish1$yday+1 # get day-of-year at finish
            years1<-seq(yearstart1,yearfinish1,1) # get sequence of years to do

            for(i in 1:length(years1)){ # start loop through years
              # syntax for query
              if(length(years1)==1){ # doing a period within a year
                query<-paste0("SELECT a.latitude, a.longitude, b.* FROM [AWAPDaily].[dbo].[latlon] as a
                    , [AWAPDaily].[dbo].[",years1[i],"] as b where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day between ",doystart," and ",doyfinish,")
                  order by b.day")
              }else{
                if(i==1){ # doing first year, start at day requested
                  query<-paste0("SELECT a.latitude, a.longitude, b.* FROM [AWAPDaily].[dbo].[latlon] as a
                     , [AWAPDaily].[dbo].[",years1[i],"] as b where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day >= ",doystart,")
                     order by b.day")
                }else{
                  if(i==length(years1)){ # doing last year, only go up to last day requested
                    query<-paste0("SELECT a.latitude, a.longitude, b.* FROM [AWAPDaily].[dbo].[latlon] as a
                       , [AWAPDaily].[dbo].[",years1[i],"] as b where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day <= ",doyfinish,")
                       order by b.day")
                  }else{ # doing in between years, so get all data for this year
                    query<-paste0("SELECT a.latitude, a.longitude, b.* FROM [AWAPDaily].[dbo].[latlon] as a
                       , [AWAPDaily].[dbo].[",years1[i],"] as b where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
                      order by b.day")
                  }}}
              # exectue query and concatenate if necessary
              if(i==1){
                output1<- sqlQuery(channel,query)
              }else{
                output1<-rbind(output1,sqlQuery(channel,query))
              }
            } # end loop through years
            output1$clearsky<-leapfix(clearskysum, seq(1990, 2014)) * 3600 / 1e6
            glm_sol<-coefficients(with(output1,glm(sol~rr+tmax+tmin+day+clearsky)))
            output_AWAPDaily$clearsky<-leapfix(clearskysum, yearlist) * 3600 / 1e6
            output_AWAPDaily[,9]<-glm_sol[1]+glm_sol[2]*output_AWAPDaily$rr+glm_sol[3]*output_AWAPDaily$tmax+glm_sol[4]*output_AWAPDaily$tmin+glm_sol[5]*output_AWAPDaily$day+glm_sol[6]*output_AWAPDaily$clearsky
            if(scenario!=""){
              output_AWAPDaily[,9]=output_AWAPDaily[,9]*SOLAR_diff
            }
            cloud<-(1-as.data.frame(output_AWAPDaily$sol)/as.data.frame(output_AWAPDaily$clearsky))*100
            cloud[cloud<0]<-0
            cloud[cloud>100]<-100
            cloud<-as.matrix(cbind(output_AWAPDaily[,4],cloud))
            CCMAXX<-cloud[,2]
            CCMINN<-CCMAXX
          }# end check for year 1990 or later
          if(ystart>1970){ #vapour pressure data available
            if(is.na(output_AWAPDaily[1,8])==TRUE){
              output_AWAPDaily[1,8]=mean(output_AWAPDaily[,8],na.rm=TRUE)
            }
            VAPRES<-zoo::na.approx(output_AWAPDaily[,8])
            VAPRES<-VAPRES*100 # convert from hectopascals to pascals
            es <- WETAIR(db = TMAXX, rh = 100)$esat
            RHMINN <- (VAPRES / es) * 100
            RHMINN[RHMINN>100]<-100
            RHMINN[RHMINN<0]<-0.01
            es <- WETAIR(db = TMINN, rh = 100)$esat
            RHMAXX <- (VAPRES / es) * 100
            RHMAXX[RHMAXX>100]<-100
            RHMAXX[RHMAXX<0]<-0.01
            if(scenario!=""){
              RHMINN<-RHMINN+RH_diff
            }
            if(scenario!=""){
              RHMAXX<-RHMAXX+RH_diff
            }
          }else{
            if(exists("output1")==FALSE){
              datestart1<-"01/01/1990" # day, month, year
              datefinish1<-"31/12/2014" # day, month, year
              datestart1<-strptime(datestart1, "%d/%m/%Y") # convert to date format
              datefinish1<-strptime(datefinish1, "%d/%m/%Y") # convert to date format
              yearstart1<-as.numeric(format(datestart1, "%Y")) # yet year start
              yearfinish1<-as.numeric(format(datefinish1, "%Y")) # yet year finish
              years1<-seq(yearstart1,yearfinish1,1) # get sequence of years to d0
              doystart<-datestart1$yday+1 # get day-of-year at start
              doyfinish<-datefinish1$yday+1 # get day-of-year at finish
              years1<-seq(yearstart1,yearfinish1,1) # get sequence of years to do
              for(i in 1:length(years1)){ # start loop through years
                # syntax for query
                if(length(years1)==1){ # doing a period within a year
                  query<-paste0("SELECT a.latitude, a.longitude, b.* FROM [AWAPDaily].[dbo].[latlon] as a
                    , [AWAPDaily].[dbo].[",years1[i],"] as b where (a.id = b.id) and (a.latitude between ",
                                lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day between
                    ",doystart," and ",doyfinish,") order by b.day")
                }else{
                  if(i==1){ # doing first year, start at day requested
                    query<-paste0("SELECT a.latitude, a.longitude, b.* FROM [AWAPDaily].[dbo].[latlon] as a
                       , [AWAPDaily].[dbo].[",years1[i],"] as b where (a.id = b.id) and (a.latitude between ",
                                  lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day >= ",
                                  doystart,") order by b.day")
                  }else{
                    if(i==length(years1)){ # doing last year, only go up to last day requested
                      query<-paste0("SELECT a.latitude, a.longitude, b.* FROM [AWAPDaily].[dbo].[latlon] as a
                        , [AWAPDaily].[dbo].[",years1[i],"] as b where (a.id = b.id) and (a.latitude between "
                                    ,lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day <= ",
                                    doyfinish,") order by b.day")
                    }else{ # doing in between years, so get all data for this year
                      query<-paste0("SELECT a.latitude, a.longitude, b.* FROM [AWAPDaily].[dbo].[latlon] as a
                        , [AWAPDaily].[dbo].[",years1[i],"] as b where (a.id = b.id) and (a.latitude between
                        ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") order by b.day")
                    }}}
                # exectue query and concatenate if necessary
                if(i==1){
                  output1<- sqlQuery(channel,query)
                }else{
                  output1<-rbind(output1,sqlQuery(channel,query))
                }
              } # end loop through years
            }
            glm_vpr<-coefficients(with(output1,glm(vpr~rr+tmax+tmin+day)))
            output_AWAPDaily[,8]<-glm_vpr[1]+glm_vpr[2]*output_AWAPDaily$rr+glm_vpr[3]*output_AWAPDaily$tmax+glm_vpr[4]*output_AWAPDaily$tmin+glm_vpr[5]*output_AWAPDaily$day
            VAPRES<-zoo::na.approx(output_AWAPDaily[,8])
            VAPRES<-VAPRES*100 # convert from hectopascals to pascals
            # correct for potential change in RH with elevation-corrected Tair
            es <- WETAIR(db = TMAXX, rh = 100)$esat
            RHMINN <- (VAPRES / es) * 100
            RHMINN[RHMINN>100]<-100
            RHMINN[RHMINN<0]<-0.01
            es <- WETAIR(db = TMINN, rh = 100)$esat
            RHMAXX <- (VAPRES / es) * 100
            RHMAXX[RHMAXX>100]<-100
            RHMAXX[RHMAXX<0]<-0.01
            if(scenario!=""){
              RHMINN<-RHMINN+RH_diff
            }
            if(scenario!=""){
              RHMAXX<-RHMAXX+RH_diff
            }
          }#end check for year is 1971 or later
        }
        if(adiab_cor==1){
          RHMAXX.orig <- RHMAXX
          RHMINN.orig <- RHMINN
          # correct for potential change in RH with elevation-corrected Tair
          es <- WETAIR(db = TMAXX, rh = 100)$esat
          e <- WETAIR(db = TMAXX.orig, rh = RHMINN.orig)$e
          RHMINN <- (e / es) * 100
          RHMINN[RHMINN>100]<-100
          RHMINN[RHMINN<0]<-0.01
          es <- WETAIR(db = TMINN, rh = 100)$esat
          e <- WETAIR(db = TMINN.orig, rh = RHMAXX.orig)$e
          RHMAXX <- (e / es) * 100
          RHMAXX[RHMAXX>100]<-100
          RHMAXX[RHMAXX<0]<-0.01
        }
        # AUSCLIM query statements
        clouds<-paste("select cloud1,cloud2,cloud3,cloud4,cloud5,cloud6,cloud7,cloud8,cloud9,cloud10,cloud11,cloud12 FROM cloudcover WHERE i = ",dbrow,sep="")
        maxwinds<-paste("select maxwind1,maxwind2,maxwind3,maxwind4,maxwind5,maxwind6,maxwind7,maxwind8,maxwind9,maxwind10,maxwind11,maxwind12 FROM maxwind WHERE i = ",dbrow,sep="")
        minwinds<-paste("select minwind1,minwind2,minwind3,minwind4,minwind5,minwind6,minwind7,minwind8,minwind9,minwind10,minwind11,minwind12 FROM minwind WHERE i = ",dbrow,sep="")
        maxhumidities<-paste("select maxhum1,maxhum2,maxhum3,maxhum4,maxhum5,maxhum6,maxhum7,maxhum8,maxhum9,maxhum10,maxhum11,maxhum12 FROM maxhum WHERE i = ",dbrow,sep="")
        minhumidities<-paste("select minhum1,minhum2,minhum3,minhum4,minhum5,minhum6,minhum7,minhum8,minhum9,minhum10,minhum11,minhum12 FROM minhum WHERE i = ",dbrow,sep="")
        rainfall<-paste("select rainfall1,rainfall2,rainfall3,rainfall4,rainfall5,rainfall6,rainfall7,rainfall8,rainfall9,rainfall10,rainfall11,rainfall12 FROM rainfall WHERE i = ",dbrow,sep="")
        rainydays<-paste("select rainy1,rainy2,rainy3,rainy4,rainy5,rainy6,rainy7,rainy8,rainy9,rainy10,rainy11,rainy12 FROM rainydays WHERE i = ",dbrow,sep="")
        ALLMINTEMPS<-TMINN
        ALLMAXTEMPS<-TMAXX
        ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)
        if(opendap == 0){
          WNMAXX <- sqlQuery(channel2,maxwinds)
          WNMINN <- sqlQuery(channel2,minwinds)
          if(dailywind!=1 ){
            WNMAXX1 <-suppressWarnings(spline(doys12,WNMAXX,n=365,xmin=1,xmax=365,method="periodic"))
            WNMAXX<-leapfix(WNMAXX1$y,yearlist)
            WNMINN1 <-suppressWarnings(spline(doys12,WNMINN,n=365,xmin=1,xmax=365,method="periodic"))
            WNMINN<-leapfix(WNMINN1$y,yearlist)
            if(scenario!=""){
              WNMAXX=WNMAXX*WIND_diff
              WNMINN=WNMINN*WIND_diff
            }
          }
        }
        if(soildata==1){
          SLES<-suppressWarnings(spline(doys12,SLES,n=365,xmin=1,xmax=365,method="periodic"))$y
          SLES <- leapfix(SLES, yearlist)
          maxshades1 <-suppressWarnings(spline(doys12,shademax,n=365,xmin=1,xmax=365,method="periodic"))
          MAXSHADES<-leapfix(maxshades1$y*100,yearlist)
          MAXSHADES<-MAXSHADES[1:ndays]
          if(manualshade==1){
            if(is.na(MAXSHADES)){
              maxshades <- rep(maxshade,ndays)
              MAXSHADES<-maxshades
            }
            if(is.na(MINSHADES)){
              minshades <- rep(minshade,ndays)
              MINSHADES<-minshades
            }
          }
        }else{
          if(manualshade==0){
            maxshades1 <-suppressWarnings(spline(doys12,shademax,n=365,xmin=1,xmax=365,method="periodic"))
            MAXSHADES<-leapfix(maxshades1$y*100,yearlist)
            minshades <- rep(minshade,ndays)
            MINSHADES<-minshades
          }else{
            MAXSHADES<-maxshades
            MINSHADES<-minshades
          }
        }
        REFLS <- rep(REFL, ndays)
        if((soildata==1)&(length(RAINFALL)>0)){
          soilwet<-RAINFALL
          soilwet[soilwet<=rainwet] = 0
          soilwet[soilwet>0] = 90
          PCTWET<-pmax(soilwet,PCTWET)
        }else{
          REFLS <- rep(REFL, ndays)
          PCTWET <- rep(PCTWET, ndays)
          soilwet<-RAINFALL
          soilwet[soilwet<=rainwet] = 0
          soilwet[soilwet>0] = 90
          PCTWET<-pmax(soilwet,PCTWET)
        }
        Numtyps <- 10 # number of substrate types
        Nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
        Nodes[1:10,] <- c(1:10) # deepest nodes for each substrate type
        if(timezone==1){
          if(!require(geonames)){
            stop('package "geonames" is required.')
          }
          ALREF<-(GNtimezone(longlat[2],longlat[1])[4])*-15
        }else{
          ALREF <- abs(trunc(x[1]))
        }
        HEMIS <- ifelse(x[2]<0,2.,1.)
        ALAT <- abs(trunc(x[2]))
        AMINUT <- (abs(x[2])-ALAT)*60
        ALONG <- abs(trunc(x[1]))
        ALMINT <- (abs(x[1])-ALONG)*60
        ALTT<-ALTITUDES
        SLOPE<-SLOPES
        AZMUTH<-AZMUTHS

        avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
        soilinit<-rep(avetemp,20)
        tannul<-mean(unlist(ALLTEMPS))

        if(nyears==1){
          avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
          if(ystart %in% leapyears){
            tannulrun<-rep(avetemp,366)
          }else{
            tannulrun<-rep(avetemp,365)
          }
        }else{
          if(nrow(TMAXX)==1){
            avetemp<-rowMeans(t(rbind(TMAXX, TMINN)), na.rm=TRUE)
          }else{
            avetemp<-rowMeans(cbind(TMAXX, TMINN), na.rm=TRUE)
          }
          if(length(TMAXX)<365){
            tannulrun<-rep((sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2),length(TMAXX))
          }else{
            tannulrun<-raster::movingFun(avetemp,n=365,fun=mean,type='to')
            yearone<-rep((sum(TMAXX[1:365])+sum(TMINN[1:365]))/(365*2),365)
            tannulrun[1:365]<-yearone
            # SST
          }
        }
        hourly <- 0
        if(microclima == 1){
          hourly <- 2
          cat('using microclima and elevatr to adjust solar for topographic and vegetation effects \n')
          if (!require("microclima", quietly = TRUE)) {
            stop("package 'microclima' is needed. Please install it.",
                 call. = FALSE)
          }
          if (!require("zoo", quietly = TRUE)) {
            stop("package 'zoo' is needed. Please install it.",
                 call. = FALSE)
          }
          cat("Downloading digital elevation data \n")
          lat <- x[2]
          long <- x[1]
          tt <- seq(as.POSIXct(paste0('01/01/',ystart), format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(paste0('31/12/',yfinish), format = "%d/%m/%Y", tz = 'UTC')+23*3600, by = 'hours')
          timediff <- x[1]/15
          hour.microclima <- as.numeric(format(tt, "%H")) + timediff-floor(timediff)
          jd <- julday(as.numeric(format(tt, "%Y")), as.numeric(format(tt, "%m")), as.numeric(format(tt, "%d")))
          dem <- microclima::get_dem(r = NA, lat = lat, long = long, resolution = 100, zmin = -20)
          xy <- data.frame(x = long, y = lat)
          coordinates(xy) = ~x + y
          proj4string(xy) = "+init=epsg:4326"
          xy <- as.data.frame(spTransform(xy, crs(dem)))
          if (class(slope) == "logical") {
            slope <- terrain(dem, unit = "degrees")
            slope <- extract(slope, xy)
          }
          if (class(aspect) == "logical") {
            aspect <- terrain(dem, opt = "aspect", unit = "degrees")
            aspect <- extract(aspect, xy)
          }
          ha <- 0
          if(is.na(hori[1]) == "TRUE"){
            ha36 <- 0
            for (i in 0:35) {
              har <- horizonangle(dem, i * 10, res(dem)[1])
              ha36[i + 1] <- atan(extract(har, xy)) * (180/pi)
            }
          }else{
            ha36 <- spline(x = hori, n = 36, method =  'periodic')$y
            ha36[ha36 < 0] <- 0
            ha36[ha36 > 90] <- 90
          }
          for (i in 1:length(hour.microclima)) {
            saz <- solazi(hour.microclima[i], lat, long, jd[i], merid = long)
            saz <- round(saz/10, 0) + 1
            saz <- ifelse(saz > 36, 1, saz)
            ha[i] <- ha36[saz]
          }
          #demmeso <- dem
          #info <- .eleveffects(hourlydata, demmeso, lat, long, windthresh = 4.5, emthresh = 0.78)
          #elev <- info$tout
          cloudhr <- cbind(rep(seq(1,length(cloud)),24), rep(cloud, 24))
          cloudhr <- cloudhr[order(cloudhr[,1]),]
          cloudhr <- cloudhr[,2]
          cloudhr <- leapfix(cloudhr, yearlist, 24)
          dsw2 <- leapfix(clearskyrad[,2], yearlist, 24) *(0.36+0.64*(1-cloudhr/100)) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
          # partition total solar into diffuse and direct using code from microclima::hourlyNCEP
          si <- microclima::siflat(hour.microclima, lat, long, jd, merid = long)
          am <- microclima::airmasscoef(hour.microclima, lat, long, jd, merid = long)
          dp <- vector(length = length(jd))
          for (i in 1:length(jd)) {
            dp[i] <- microclima:::difprop(dsw2[i], jd[i], hour.microclima[i], lat, long, watts = TRUE, hourly = TRUE, merid = long)
          }
          dp[dsw2 == 0] <- NA
          dnir <- (dsw2 * (1 - dp))/si
          dnir[si == 0] <- NA
          difr <- (dsw2 * dp)
          edni <- dnir/((4.87/0.0036) * (1 - dp))
          edif <- difr/((4.87/0.0036) * dp)
          bound <- function(x, mn = 0, mx = 1) {
            x[x > mx] <- mx
            x[x < mn] <- mn
            x
          }
          odni <- bound((log(edni)/-am), mn = 0.001, mx = 1.7)
          odif <- bound((log(edif)/-am), mn = 0.001, mx = 1.7)
          nd <- length(odni)
          sel <- which(is.na(am * dp * odni * odif) == F)
          dp[1] <- dp[min(sel)]
          odni[1] <- odni[min(sel)]
          odif[1] <- odif[min(sel)]
          dp[nd] <- dp[max(sel)]
          odni[nd] <- odni[max(sel)]
          odif[nd] <- odif[max(sel)]
          dp[nd] <- dp[max(sel)]
          odni[nd] <- odni[max(sel)]
          odif[nd] <- odif[max(sel)]
          if (!require("raster", quietly = TRUE)) {
            stop("package 'raster' is needed. Please install it.",
                 call. = FALSE)
          }
          dp <- na.approx(dp, na.rm = F)
          odni <- na.approx(odni, na.rm = F)
          odif <- na.approx(odif, na.rm = F)
          h_dp <- bound(dp)
          h_oi <- bound(odni, mn = 0.24, mx = 1.7)
          h_od <- bound(odif, mn = 0.24, mx = 1.7)
          afi <- exp(-am * h_oi)
          afd <- exp(-am * h_od)
          h_dni <- (1 - h_dp) * afi * 4.87/0.0036
          h_dif <- h_dp * afd * 4.87/0.0036
          h_dni[si == 0] <- 0
          h_dif[is.na(h_dif)] <- 0
          radwind2 <- .shortwave.ts(h_dni * 0.0036, h_dif * 0.0036, jd, hour.microclima, lat, long, slope, aspect, ha = ha, svv = 1, x = LOR, l = mean(microclima.LAI), albr = 0, merid = long, dst = 0, difani = FALSE)
          #microclima.out$hourlyradwind <- radwind2
          SOLRhr <- radwind2$swrad / 0.0036
          VIEWF <- 1 # accounted for already in microclima cals
          hori <- rep(0, 24) # accounted for already in microclima calcs
        }

        if(opendap == 0){
          # correct for fact that wind is measured at 2 m height
          # wind shear equation v / vo = (h / ho)^a
          #where
          #v = the velocity at height h (m/s)
          #vo = the velocity at height ho (m/s)
          #a = the wind shear exponent
          #Terrain   Wind Shear Exponent
          #- a -
          #  Open water   0.1
          #Smooth, level, grass-covered   0.15
          #Row crops   0.2
          #Low bushes with a few trees 	0.2
          #Heavy trees 	0.25
          #Several buildings 	0.25
          #Hilly, mountainous terrain 	0.25
          if(dailywind!=1){
            WNMINN<-WNMINN*(1.2/10)^0.15*.1 # reduce min wind further because have only 9am/3pm values to get max/min
            WNMAXX<-WNMAXX*(1.2/10)^0.15
            message('min wind * 0.1 \n')
          }else{
            WNMAXX<-dwind*(1.2/2)^0.15
            WNMINN<-WNMAXX
            WNMAXX<-WNMAXX*2
            WNMINN<-WNMINN*0.5
            WNMINN[WNMINN<0.1]<-0.1
            message('min wind * 0.5 \n')
            message('max wind * 2 \n')
          }
          CCMINN<-CCMINN*0.5
          CCMAXX<-CCMAXX*2
          CCMINN[CCMINN>100]<-100
          CCMAXX[CCMAXX>100]<-100
          if(clearsky==1){
            CCMINN=CCMINN*0
            CCMAXX=CCMAXX*0
            message('running for clear sky conditions')
          }else{
            message('min cloud * 0.5 \n')
            message('max cloud * 2 \n')
          }
        }else{
          if(clearsky==1){
            CCMINN=CCMINN*0
            CCMAXX=CCMAXX*0
            message('running for clear sky conditions')
          }
        }
        # impose uniform warming
        TMAXX<-TMAXX+warm
        TMINN<-TMINN+warm
        # impose wind multiplication factor
        WNMAXX <- WNMAXX * windfac
        WNMINN <- WNMINN * windfac
        if(soildata!=1){
          SLES<-matrix(nrow=ndays,data=0)
          SLES<-SLES+SLE
        }
        #quick fix to make it so that MINSHADES is at the user-specified value and MAXSHADES is from the FAPAR database
        if(soildata==1 & manualshade==0){
          MINSHADES<-MAXSHADES
          MINSHADES[1:length(MINSHADES)]<-minshade
          MINSHADES<-MAXSHADES # this is to make one shade level effectively, as dictated by FAPAR
          MAXSHADES<-MINSHADES+0.1
        }

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

        # microclimate input parameters listALTT,ALREF,ALMINT,ALONG,AMINUT,ALAT
        ALTT<-as.numeric(ALTT)
        ALREF<-as.numeric(ALREF)
        ALMINT<-as.numeric(ALMINT)
        ALONG<-as.numeric(ALONG)
        AMINUT<-as.numeric(AMINUT)
        ALAT<-as.numeric(ALAT)

        # microclimate input parameters list
        microinput<-c(ndays,RUF,ERR,Usrhyt,Refhyt,Numtyps,Z01,Z02,ZH1,ZH2,idayst,ida,HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,CMH2O,microdaily,tannul,EC,VIEWF,snowtemp,snowdens,snowmelt,undercatch,rainmult,runshade,runmoist,maxpool,evenrain,snowmodel,rainmelt,writecsv,densfun,hourly,rainhourly,lamb,IUV,RW,PC,RL,SP,R1,IM,MAXCOUNT,IR,message,fail,snowcond,intercept,grasshade,solonly,ZH,D0)

        # hourly option set to 0, so make empty vectors
        if(hourly!=1){
          TAIRhr=rep(0,24*ndays)
          RHhr=rep(0,24*ndays)
          WNhr=rep(0,24*ndays)
          CLDhr=rep(0,24*ndays)
          ZENhr=rep(-1,24*ndays)
          IRDhr=rep(-1,24*ndays)
        }
        if(hourly==0){
          SOLRhr=rep(0,24*ndays)
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
        SLES1[1:ndays]<-SLES[1:ndays]
        MAXSHADES1[1:ndays]<-MAXSHADES
        MINSHADES1[1:ndays]<-MINSHADES
        TMAXX1[1:ndays]<-TMAXX
        TMINN1[1:ndays]<-TMINN
        CCMAXX1[1:ndays]<-CCMAXX
        CCMINN1[1:ndays]<-CCMINN
        RHMAXX1[1:ndays]<-RHMAXX
        RHMINN1[1:ndays]<-RHMINN
        WNMAXX1[1:ndays]<-WNMAXX[1:ndays]
        WNMINN1[1:ndays]<-WNMINN[1:ndays]
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
        message(paste('running microclimate model for',ndays,'days from ',ystart, ' to ', yfinish, ' at site',location,'\n'))
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
          plant<-microut$plant # retrieve plant , minimum shade
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
        dates <- seq(as.POSIXct(paste0("01/01/",ystart), format = "%d/%m/%Y", tz = 'Etc/GMT+10'), as.POSIXct(paste0("01/01/",yfinish+1), format = "%d/%m/%Y ", tz = 'Etc/GMT+10'), by = 'hours')[1:(length(TMAXX1)*24)]
        dates2 <- seq(as.POSIXct(paste0("01/01/",ystart), format = "%d/%m/%Y", tz = 'Etc/GMT+10'), as.POSIXct(paste0("01/01/",yfinish+1), format = "%d/%m/%Y", tz = 'Etc/GMT+10'), by = 'days')[1:length(TMAXX1)]
        if(lamb == 1){
          drlam<-as.data.frame(microut$drlam) # retrieve direct solar irradiance
          drrlam<-as.data.frame(microut$drrlam) # retrieve direct Rayleigh component solar irradiance
          srlam<-as.data.frame(microut$srlam) # retrieve scattered solar irradiance
          if(snowmodel == 1){
            return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,minshade=minshade,maxshade=maxshade,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam,dates=dates,dates2=dates2))
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
    } # end of check if soil data is being used but no soil data returned
  } # end error trapping
} # end of micro_aust function
