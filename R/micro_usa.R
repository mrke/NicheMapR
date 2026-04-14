#' USA implementation of the microclimate model.
#'
#' An implementation of the NicheMapR microclimate model that uses the GRIDMET daily weather database http://www.climatologylab.org/gridmet.html, and specifically uses the following variables: pr, rmax, rmin, srad, tmmn, tmmx, vs. Also uses the following DEM "metdata_elevationdata.nc".
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param dstart First day to run, date in format "d/m/Y" e.g. "01/01/2016"
#' @param dfinish Last day to run, date in format "d/m/Y" e.g. "31/12/2016"
#' @param REFL Soil solar reflectance, decimal \%
#' @param elev Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param minshade Minimum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param maxshade Maximum shade level to us (\%) (can be a single value or a vector of daily values)
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
#' @return tcond Hourly predictions of the soil thermal conductivity under the minimum specified shade
#' @return shadtcond Hourly predictions of the soil thermal conductivity under the maximum specified shade
#' @return specheat Hourly predictions of the soil specific heat capacity under the minimum specified shade
#' @return shadspecheat Hourly predictions of soil specific heat capacity under the maximum specified shade
#' @return densit Hourly predictions of the soil density under the minimum specified shade
#' @return shaddensit Hourly predictions of the soil density under the maximum specified shade
#' @usage micro_USA(loc = c(-117.8919161, 36.582016), dstart = "01/01/2016", dfinish = "31/12/2016",
#' REFL = 0.15, slope = 0, aspect = 0, DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minshade = 0, maxshade = 90,
#' Usrhyt = 0.01, ...)
#' @export
#' @details
#' \strong{Parameters controlling how the model runs:}\cr\cr
#'
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes (Fortran version), 2=yes (R version), 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{solonly}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{ndmax}{ = 3, iterations of first day to get a steady periodic}\cr\cr
#' \code{Soil_Init}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{microclima}{ = 0, Use microclima and elevatr package to adjust solar radiation for terrain? 1 = yes, 0 = no}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{windfac}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{adiab_cor}{ = 1, use adiabatic lapse rate correction? 1=yes, 0=no}\cr\cr
#' \code{warm}{ = 0, warming offset vector, °C (negative values mean cooling). Can supply a single value or a vector the length of the number of days to be simulated.}\cr\cr
#' \code{spatial}{ = "c:/Australian Environment/", choose location of terrain data}\cr\cr
#' \code{opendap}{ = 1, query met grids via opendap}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
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
#' \code{rainoff}{ = 0, Rain offset (mm), used to induce changes in rainfall from GRIDMET values. Can be a single value or a vector matching the number of days to simulate. If negative values are used, rainfall will be prevented from becomming negative.}\cr\cr
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
#' \code{LAI}{ = 0.1, leaf area index (can be a single value or a vector of daily values), used to partition traspiration/evaporation from PET}\cr\cr
#' \code{microclima.LAI}{ = 0, leaf area index, used by package microclima for radiation calcs}\cr\cr
#' \code{microclima.LOR}{ = 1, leaf orientation for package microclima radiation calcs}\cr\cr
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
#' \code{tides}
#' { is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0, nrow = length(seq(as.POSIXct(dstart, format = '%d/%m/%Y'), as.POSIXct(dfinish, format = '%d/%m/%Y'), by = 'days')) * 24, ncol = 3), matrix of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
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
#' \code{minshade}{ - minimum shade for each day of simulation (\%)}\cr\cr
#' \code{maxshade}{ - maximum shade for each day of simulation (\%)}\cr\cr
#' \code{DEP}{ - vector of depths used (cm)}\cr\cr
#' \code{diffuse_frac}{ - vector of hourly values of the fraction of total solar radiation that is diffuse (-), computed by microclima if microclima > 0}\cr\cr
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
#' \item 13 SOLR - solar radiation (W/m2) (unshaded, horizontal plane)
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
#' tcond and shadtcond variables:
#' \itemize{
#' \item  1 DOY - day-of-year
#' \item  2 TIME - time of day (mins)
#' \item  3-12 TC0cm ... - soil thermal conductivity (W/m-K), at each of the 10 specified depths
#' }
#' specheat and shadspecheat variables:
#' \itemize{
#' \item  1 DOY - day-of-year
#' \item  2 TIME - time of day (mins)
#' \item  3-12 SP0cm ... - soil specific heat capacity (J/kg-K), at each of the 10 specified depths
#' }
#' densit and shaddensit variables:
#' \itemize{
#' \item  1 DOY - day-of-year
#' \item  2 TIME - time of day (mins)
#' \item  3-12 DE0cm ... - soil density (Mg/m3), at each of the 10 specified depths
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
#' dstart <- "01/01/2016"
#' dfinish <- "31/12/2017"
#' micro<-micro_usa(dstart = dstart, dfinish = dfinish) # run the model at the default location (Madison, Wisconsin) for 2016 to 2017 using opendap
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
  #terrain = 0,
  windfac = 1,
  adiab_cor = 1,
  warm = 0,
  spatial = "P:",
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
  microclima.LAI = 0,
  microclima.LOR = 1,
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
  ndmax = 3,
  microclima = 0,
  microclima.dem.res = 100,
  microclima.zmin = -20,
  opendap = 1,
  soilgrids = 0,
  IR = 0,
  message = 0,
  fail = nyears * 24 * 365,
  save = 0,
  runmicro = 1,
  snowcond = 0,
  intercept = max(maxshade) / 100 * 0.3,
  grasshade = 0,
  maxsurf = 85) { # end function parameters

  ystart <- as.numeric(substr(dstart, 7, 10))
  yfinish <- as.numeric(substr(dfinish, 7, 10))
  yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
  # error trapping - originally inside the Fortran code, but now checking before executing Fortran
  # --- 1. Input validation ---
  Refhyt <- 2  # GRIDMET met data reference height (m) — 2 m standard screen height
  errors <- 0
  # usa-specific checks
  if(is.numeric(loc[1])){
    if(loc[1] > 180 | loc[2] > 90){
      message("ERROR: Latitude or longitude (longlat) is out of bounds. Please enter a correct value.")
      errors <- 1
    }
  }
  if(timezone %in% c(0, 1) == FALSE){
    message("ERROR: the variable 'timezone' must be either 0 or 1. Please correct.")
    errors <- 1
  }
  if(run.gads == 1){
    message("If program is crashing, try run.gads = 2.")
  }
  if(write_input %in% c(0, 1) == FALSE){
    message("ERROR: the variable 'write_input' must be either 0 or 1. Please correct.")
    errors <- 1
  }
  # common parameter checks
  errors <- errors + validate_micro_inputs(DEP, REFL, slope, aspect, hori, SLE, ERR,
    RUF, D0, Usrhyt, Refhyt, EC, CMH2O, TIMAXS, TIMINS, minshade, maxshade,
    Thcond, Density, SpecHeat, BulkDensity, run.gads = run.gads)

  if(errors==0){ # continue

    ################## time related variables #################################
    doys12<-c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month

    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial ndmax day burn in

    daystart<-1
    idayst <- 1 # start day

    ################## location and terrain #################################
    if (!require("terra", quietly = TRUE)) {
      stop("package 'terra' is needed. Please install it.",
        call. = FALSE)
    }
    if (!require("RNetCDF", quietly = TRUE)) {
      stop("package 'RNetCDF' is needed. Please install it.",
        call. = FALSE)
    }
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    require("terra")
    require("RNetCDF")

    # reference meridian for solar noon calculation
    ALREF <- get_timezone_alref(lon = x[1], lat = x[2], timezone = timezone)
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
        USADEM <- as.numeric(terra::extract(terra::rast(paste0(spatial,"/metdata_elevationdata.nc")), x)[, 2]) # metres above sea level
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

    ALTITUDES <- NA# as.numeric(terra::extract(rast(paste0(spatial,"/terr50.tif")), x)) # to do
    if(is.na(elev) == FALSE){ALTITUDES <- elev} # check if user-specified elevation
    if(is.na(ALTITUDES)==TRUE){ALTITUDES<-USADEM}
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
    #if(terrain==1){ # to do
      # cat("extracting terrain data \n")
      # # now extract terrain data from elevslpasphori.nc
      # # get UTM from dec degrees, NZTM
      # HORIZONS <- elevslpasphori[4:27]
      # SLOPES <- elevslpasphori[2]
      # AZMUTHS <- elevslpasphori[3]
      # # the horizons have been arranged so that they go from 0 degrees azimuth (north) clockwise - r.horizon starts
      # # in the east and goes counter clockwise!
      # HORIZONS <- (ifelse(is.na(HORIZONS),0,HORIZONS))/10 # get rid of na and get back to floating point
      # HORIZONS <- data.frame(HORIZONS)
    #}else{
      HORIZONS <- hori
      HORIZONS <- data.frame(HORIZONS)
      VIEWF_all <- 1-sum(sin(as.data.frame(hori)*pi/180))/length(hori) # convert horizon angles to radians and calc view factor(s)
      SLOPES<-rep(slope,length(x[,1]))
      AZMUTHS<-rep(aspect,length(x[,1]))
    #}
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
        nc <- RNetCDF::open.nc(paste0(baseurl, "tmmn_1979_CurrentYear_CONUS.nc#fillmismatch"))
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
        nc <- RNetCDF::open.nc(paste0(baseurl, "tmmx_1979_CurrentYear_CONUS.nc#fillmismatch"))
        Tmax <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_maximum_temperature", start = start, count = count, unpack = TRUE))) - 273.15
        RNetCDF::close.nc(nc)
        cat(paste0("rhmin weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "rmin_1979_CurrentYear_CONUS.nc#fillmismatch"))
        rhmin <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_minimum_relative_humidity", start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
        cat(paste0("rhmax weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "rmax_1979_CurrentYear_CONUS.nc#fillmismatch"))
        rhmax <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_maximum_relative_humidity", start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
        cat(paste0("rain weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "pr_1979_CurrentYear_CONUS.nc#fillmismatch"))
        Rain <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "precipitation_amount", start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
        cat(paste0("solar weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "srad_1979_CurrentYear_CONUS.nc#fillmismatch"))
        solar <- retry(as.numeric(RNetCDF::var.get.nc(nc, variable = "daily_mean_shortwave_radiation_at_surface", start = start, count = count, unpack = TRUE)))
        RNetCDF::close.nc(nc)
        cat(paste0("wind weather input \n"))
        nc <- RNetCDF::open.nc(paste0(baseurl, "vs_1979_CurrentYear_CONUS.nc#fillmismatch"))
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
    ida <- ndays
    idayst <- 1

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

    if(is.na(ALTITUDES)!=TRUE){

      TAI <- compute_tai(longlat, run.gads, TAI_AUSTRALIA)

      if(adiab_cor==1){
        TMAXX<-as.matrix(Tmax+adiab_corr_max)
        TMINN<-as.matrix(Tmin+adiab_corr_min)
      }else{
        TMAXX<-as.matrix(Tmax)
        TMINN<-as.matrix(Tmin)
      }
      if(warm != 0){
        # impose uniform temperature change
        TMAXX<-TMAXX+warm
        TMINN<-TMINN+warm
      }
      RAINFALL<-Rain+rainoff
      RAINFALL[RAINFALL < 0] <- 0

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
      ALREF <- get_timezone_alref(lon = x[1], lat = x[2], timezone = timezone)

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
      tannul<-mean(unlist(ALLTEMPS))

      if(nyears==1){
        avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
        deepsoil<-rep(avetemp,ndays)
      }else{
        if(nrow(TMAXX)==1){
          avetemp<-colMeans(cbind(TMAXX, TMINN), na.rm=TRUE)
        }else{
          avetemp<-rowMeans(cbind(TMAXX, TMINN), na.rm=TRUE)
        }
        if(length(TMAXX)<365){
          deepsoil<-rep((sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2),length(TMAXX))
        }else{
          deepsoil<-terra::roll(avetemp,n=365,fun=mean,type='to')
          yearone<-rep((sum(TMAXX[1:365])+sum(TMINN[1:365]))/(365*2),365)
          deepsoil[1:365]<-yearone
        }
      }

      if(microclima == 1){
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
        timediff <- x[1] / 15
        hour.microclima <- as.numeric(format(tt, "%H")) + timediff-floor(timediff)
        jd <- julday(as.numeric(format(tt, "%Y")), as.numeric(format(tt, "%m")), as.numeric(format(tt, "%d")))
        dem_result <- fetch_dem(loc[1], loc[2],
                                dem.res = microclima.dem.res, zmin = microclima.zmin,
                                pixels = NULL,
                                terrain = 1, horizon.step = 10, nangles = 36,
                                user_hori = if (!is.na(hori[1])) hori else NULL,
                                slope = slope, aspect = aspect)
        dem <- dem_result$dem; slope <- dem_result$slope; aspect <- dem_result$aspect
        ha36 <- dem_result$hori
        for (i in 1:length(hour.microclima)) {
          saz <- solazi(hour.microclima[i], lat, long, jd[i], merid = long)
          saz <- round(saz/10, 0) + 1
          saz <- ifelse(saz > 36, 1, saz)
          ha[i] <- ha36[saz]
        }
        #demmeso <- dem
        #info <- .eleveffects(hourlydata, demmeso, lat, long, windthresh = 4.5, emthresh = 0.78)
        #elev <- info$tout
        cloudhr <- cbind(rep(seq(1, length(cloud)),24), rep(cloud, 24))
        cloudhr <- cloudhr[order(cloudhr[,1]),]
        cloudhr <- cloudhr[,2]
        cloudhr <- leapfix(cloudhr, yearlist, 24)
        dsw2 <- leapfix(clearskyrad[,2], yearlist, 24) *(0.36+0.64*(1-cloudhr/100)) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
        # partition total solar into diffuse and direct using code from microclima::hourlyNCEP
        sol <- compute_solar_partition(dsw2, jd, hour.microclima, lat, long,
                                       slope, aspect, ha, microclima.LOR, microclima.LAI)
        SOLRhr <- sol$SOLRhr_all
        diffuse_frac <- sol$diffuse_frac_all
        VIEWF <- 1 # accounted for already in microclima cals
        hori <- rep(0, 24) # accounted for already in microclima calcs
      }else{
        diffuse_frac <- NA
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
      soilprops[,2] <- 1 - BulkDensity / Density # not used if soil moisture computed
      soilprops[soilprops[,2] < 0.26, 2] <- 0.26
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

      # --- 3. Prepare model inputs ---
      microinput <- build_microinput(ndays, RUF, ERR, Usrhyt, Refhyt, Numtyps,
        Z01, Z02, ZH1, ZH2, idayst, ida, HEMIS, ALAT, AMINUT, ALONG, ALMINT,
        ALREF, slope, azmuth, ALTT, CMH2O, microdaily, tannul, EC, VIEWF,
        snowtemp, snowdens, snowmelt, undercatch, rainmult, runshade, runmoist,
        maxpool, evenrain, snowmodel, rainmelt, writecsv, densfun, hourly,
        rainhourly, lamb, IUV, RW, PC, RL, SP, R1, IM, MAXCOUNT, IR, message,
        fail, snowcond, intercept, grasshade, solonly, ZH, D0, TIMAXS, TIMINS,
        spinup, maxsurf, ndmax)

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

      if(length(LAI)<ndays){
        LAI<-rep(LAI[1],ndays)
      }
      if(shore==0){
        tides<-matrix(data = 0, nrow = 24*ndays, ncol = 3) # make an empty matrix
      }
      # all microclimate data input list - all these variables are expected by the input argument of the fortran micro2014 subroutine
      micro <- build_micro_list(microinput, doy, SLES, DEP, Nodes,
        MAXSHADES, MINSHADES, TMAXX, TMINN, RHMAXX, RHMINN, CCMAXX, CCMINN,
        WNMAXX, WNMINN, TAIRhr, RHhr, WNhr, CLDhr, SOLRhr, RAINhr, ZENhr,
        IRDhr, REFLS, PCTWET, soilinit, hori, TAI, soilprops, moists, RAINFALL,
        deepsoil, PE, KS, BB, BD, DD, L, LAI, tides)
      # write all input to csv files in their own folder
      if(write_input == 1){
        write_micro_csv(microinput, doy, SLES, DEP, Nodes, MAXSHADES, MINSHADES,
          TIMAXS, TIMINS, TMAXX, TMINN, RHMAXX, RHMINN, CCMAXX, CCMINN,
          WNMAXX, WNMINN, REFLS, PCTWET, soilinit, hori, TAI, soilprops, moists,
          RAINFALL, deepsoil, PE, BD, DD, BB, KS, L, LAI, tides,
          TAIRhr, RHhr, WNhr, CLDhr, SOLRhr, RAINhr, ZENhr, IRDhr)
      }
      if(is.numeric(loc[1])){
        location<-paste("long",loc[1],"lat",loc[2])
      }else{
        location<-loc
      }
      if(runmicro){
        cat(paste('running microclimate model for',ndays,'days from',dstart,' to ', dfinish, ' at site ',location,'\n'))
      message('Note: the output column `SOLR` in metout and shadmet is for unshaded horizontal plane solar radiation \n')
      ptm <- proc.time() # Start timing
      microut<-microclimate(micro)
      print(proc.time() - ptm) # Stop the clock

      # --- 5. Process and return results ---
      out <- process_micro_output(microut, runmoist, snowmodel, lamb)
      if(max(out$metout[,1] == 0)){
        message("ERROR: the model crashed - try a different error tolerance (ERR) or a different spacing in DEP")
      }
      dates <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = 'America/Los_Angeles'),
                   as.POSIXct(dfinish, format = "%d/%m/%Y", tz = 'America/Los_Angeles') + 23 * 3600,
                   by = 'hours')[1:(length(TMAXX) * 24)]
      dates2 <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = 'America/Los_Angeles'),
                    as.POSIXct(dfinish, format = "%d/%m/%Y", tz = 'America/Los_Angeles'),
                    by = 'days')[1:(length(TMAXX))]
      return(build_micro_return(out, RAINFALL, ndays, ALTT, REFL,
        longlat = c(x[1], x[2]), nyears, timeinterval = ndays, MINSHADES, MAXSHADES,
        DEP, dates, dates2, PE, BD, DD, BB, KS, dem = NA, diffuse_frac,
        snowmodel, lamb))
    }else{
      return(list(RAINFALL = RAINFALL, TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN, WNMAXX = WNMAXX, WNMINN = WNMINN, CCMAXX = CCMAXX, CCMINN = CCMINN, CLDhr = CLDhr, WNhr = WNhr, TAIRhr = TAIRhr, RHhr = RHhr, RAINhr = RAINhr, SOLRhr = SOLRhr, ZENhr = ZENhr, IRDhr = IRDhr, dates = dates, dates2 = dates2, PE=PE,BD=BD,DD=DD,BB=BB,KS=KS))
    }
    } # end of check for na sites
  } # end error trapping
} # end of micro_usa function
