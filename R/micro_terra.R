#' Global historical monthly implementation of the microclimate model
#'
#' An implementation of the NicheMapR microclimate model that uses the global monthly climate database
#' described in "Abatzoglou, J. T., Dobrowski, S. Z., Parks, S. A., & Hegewisch, K. C. (2018). TerraClimate,
#' a high-resolution global dataset of monthly climate and climatic water balance from 1958–2015.
#' Scientific Data, 5(1), 170191. https://doi.org/10.1038/sdata.2017.191"
#' This dataset includes climate change scenarios for +2 and +4 deg C warming (via parameter 'scenario')
#' Aerosol attenuation can also be computed based on the Global Aerosol Data Set (GADS)
#' Koepke, P., M. Hess, I. Schult, and E. P. Shettle. 1997. Global Aerosol Data Set. Max-Planck-Institut for Meteorologie, Hamburg
#' by choosing the option 'run.gads<-1' (Fortran version, quicker but may crash on some systems) or 'run.gads<-2' (R version)
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param timeinterval The number of time intervals to generate predictions for over a year (must be 12 <= x <=365)
#' @param ystart First year to run
#' @param yfinish Last year to run
#' @param REFL Soil solar reflectance, decimal \%
#' @param elev Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param minshade Minimum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param maxshade Maximum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param Usrhyt Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
#' @param ... Additional arguments, see Details
#' @usage micro_terra(loc = c(-89.40123, 43.07305), timeinterval = 12, ystart = 2000, yfinish = 2015,
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
#' @return tcond Hourly predictions of the soil thermal conductivity under the minimum specified shade
#' @return shadtcond Hourly predictions of the soil thermal conductivity under the maximum specified shade
#' @return specheat Hourly predictions of the soil specific heat capacity under the minimum specified shade
#' @return shadspecheat Hourly predictions of soil specific heat capacity under the maximum specified shade
#' @return densit Hourly predictions of the soil density under the minimum specified shade
#' @return shaddensit Hourly predictions of the soil density under the maximum specified shade
#' @export
#' @details
#' \itemize{
#' \strong{Parameters controlling how the model runs:}\cr\cr
#'
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes (Fortran version), 2=yes (R version), 0=no}\cr\cr
#' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{solonly}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{Soil_Init}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{elevatr}{ = 0, Use elevatr package to get high resolution elevation for location? 1 = yes, 0 = no}\cr\cr
#' \code{terrain}{ = 0, Use elevatr package to adjust horizon angles, slope and aspect? 1 = yes, 0 = no}\cr\cr
#' \code{microclima}{ = 0, Use microclima and elevatr package to compute diffuse fraction of solar radiation (1) and adjust solar radiation for terrain (2)? 0 = no}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{ERR}{ = 1.5, Integrator error tolerance for soil temperature calculations}\cr\cr
#' \code{Refhyt}{ = 2, Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured}\cr\cr
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
#' \code{cap}{ = 1, organic cap present on soil surface? (cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)}\cr\cr
#' \code{CMH2O}{ = 1, Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air conditions; 2.0 = humid, tropical conditions (note this is for the whole atmospheric profile, not just near the ground)}\cr\cr
#' \code{hori}{ = rep(0,24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#' \code{lapse_min}{ = 0.0039 Lapse rate for minimum air temperature (degrees C/m)}\cr\cr
#' \code{lapse_max}{ = 0.0077 Lapse rate for maximum air temperature (degrees C/m)}\cr\cr
#' \code{TIMAXS}{ = c(1, 1, 0, 0), Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise}\cr\cr
#' \code{TIMINS}{ = c(0, 0, 1, 1), Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon}\cr\cr
#' \code{timezone}{ = 0, Use GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)? 1=yes, 0=no}\cr\cr
#' \code{TAI}{ = 0, Vector of 111 values, one per wavelenght bin, for solar attenuation - used to overide GADS}\cr\cr
#' \code{windfac}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{warm}{ = 0, warming offset vector, °C (negative values mean cooling). Can supply a single value or a vector the length of the number of days to be simulated.}\cr\cr
#' \code{terra_source}{ = NA, specify location of terraclimate data, goes to the web by default}\cr\cr
#' \code{scenario}{ = 0, TerraClimate climate change scenario, either 0, 2 or 4 °C warmer}\cr\cr
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
#' \code{LAI}{ = 0.1, leaf area index (can be a single value or a vector of daily values), used to partition traspiration/evaporation from PET}\cr\cr
#' \code{microclima.LAI}{ = 0, leaf area index, used by package microclima for radiation calcs}\cr\cr
#' \code{microclima.LOR}{ = 1, leaf orientation for package microclima radiation calcs}\cr\cr
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
#' \code{rainfrac}{ = 0.5, fraction of rain that falls on the first day of the month (decimal \% with 0 meaning rain falls evenly) - this parameter allows something other than an even intensity of rainfall when interpolating the montly rainfall data)}\cr\cr
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)}\cr\cr
#' \code{intercept}{ = max(maxshade) / 100 * 0.3, snow interception fraction for when there's shade (0-1)}\cr\cr
#' \code{grasshade}{ = 0, if 1, means shade is removed when snow is present, because shade is cast by grass/low shrubs}\cr\cr
#'
#' \strong{ Intertidal mode parameters:}
#'
#' \code{shore}{ Include tide effects? If 1, the matrix}
#' \code{tides}
#' { is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0., nrow = 24*365*nyears, ncol = 3), matrix of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#' }
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
#' \code{minshade}{ - minimum shade for each day of simulation (\%)}\cr\cr
#' \code{maxshade}{ - maximum shade for each day of simulation (\%)}\cr\cr
#' \code{dem}{ - digital elevation model obtained via 'get_dem' using package 'elevatr' (m)}\cr\cr
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
#'micro <- micro_terra() # run the model with default location (Madison, Wisconsin) from 2000 to 2015 and settings
#'
#'metout <- as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#'shadmet <- as.data.frame(micro$shadmet) # above ground microclimatic conditions, max shade
#'soil <- as.data.frame(micro$soil) # soil temperatures, minimum shade
#'shadsoil <- as.data.frame(micro$shadsoil) # soil temperatures, maximum shade
#'
#'minshade <- micro$minshade[1]
#'maxshade <- micro$maxshade[1]
#'
#'# plotting above-ground conditions in minimum shade
#'with(metout, {plot(TALOC ~ micro$dates3, xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l", main = paste("air temperature, ", minshade, "% shade",sep = ""))})
#'with(metout, {points(TAREF ~ micro$dates3, xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l",lty = 2, col = 'blue')})
#'with(metout, {plot(RHLOC ~ micro$dates3, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",ylim = c(0, 100),main = paste("humidity, ",minshade, "% shade",sep=""))})
#'with(metout, {points(RH ~ micro$dates3, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",col = 'blue',lty = 2, ylim = c(0, 100))})
#'with(metout, {plot(TSKYC ~ micro$date3s, xlab = "Date and Time", ylab = "Sky Temperature (°C)"
#',  type = "l", main = paste("sky temperature, ", minshade, "% shade", sep=""))})
#'with(metout, {plot(VREF ~ micro$dates3, xlab = "Date and Time",  ylab = "Wind Speed (m/s)"
#',  type = "l", main = "wind speed", col = 'blue',ylim = c(0, 15))})
#'with(metout, {points(VLOC ~ micro$dates3, xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#',  type = "l", lty = 2)})
#'with(metout, {plot(ZEN ~ micro$dates3,xlab = "Date and Time", ylab = "Zenith Angle of Sun (deg)"
#',  type = "l", main = "solar angle, sun")})
#'with(metout, {plot(SOLR ~ micro$dates3,xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
#',  type = "l", main = "solar radiation")})
#'
#'# plotting soil temperature for minimum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(soil[,i + 2] ~ micro$dates3, xlab = "Date and Time", ylab = "Soil Temperature (°C)"
#'    ,col = i, type = "l", main = paste("soil temperature ", minshade, "% shade", sep=""))
#'  }else{
#'    points(soil[,i + 2] ~ micro$dates3, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l")
#'  }
#'}
#'points(metout$SNOWDEP ~ micro$dates, type = 'h', col = 'light blue')
#
#'# plotting above-ground conditions in maximum shade
#'with(shadmet,{plot(TALOC ~ micro$dates3,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l", main = "air temperature, sun")})
#'with(shadmet,{points(TAREF ~ micro$dates3,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l", lty = 2, col = 'blue')})
#'with(shadmet,{plot(RHLOC ~ micro$dates3,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l", ylim = c(0, 100),main = "humidity, shade")})
#'with(shadmet,{points(RH ~ micro$dates3,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l", col = 'blue',lty = 2, ylim = c(0, 100))})
#'with(shadmet,{plot(TSKYC ~ micro$dates3,xlab = "Date and Time", ylab = "Sky Temperature (°C)",
#'  type = "l", main = "sky temperature, shade")})
#'
#'# plotting soil temperature for maximum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(shadsoil[,i + 2] ~ micro$dates3, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l", main = paste("soil temperature ", maxshade, "% shade", sep=""))
#'  }else{
#'    points(shadsoil[,i + 2] ~ micro$dates3, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l")
#'  }
#'}
#'points(shadmet$SNOWDEP ~ micro$dates, type = 'h', col = 'light blue')
micro_terra <- function(
  loc = c(-89.4557, 43.1379),
  ystart = 2000,
  yfinish = 2015,
  timeinterval = 12,
  nyears = yfinish - ystart + 1,
  REFL = 0.15,
  elev = NA,
  slope = 0,
  aspect = 0,
  lapse_max = 0.0077,
  lapse_min = 0.0039,
  DEP=c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
  minshade = 0,
  maxshade = 90,
  dem = NA,
  Refhyt = 2,
  Usrhyt = 0.01,
  Z01 = 0,
  Z02 = 0,
  ZH1 = 0,
  ZH2 = 0,
  runshade = 1,
  solonly = 0,
  clearsky = 0,
  run.gads = 1,
  Soil_Init = NA,
  write_input = 0,
  writecsv = 0,
  elevatr = 0,
  terrain = 0,
  microclima = 0,
  adiab_cor = 1,
  warm = 0,
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
  microclima.LAI = 0,
  microclima.LOR = 1,
  snowmodel = 1,
  snowtemp = 1.5,
  snowdens = 0.375,
  densfun = c(0.5979, 0.2178, 0.001, 0.0038),
  snowmelt = 1,
  undercatch = 1,
  rainmelt = 0.0125,
  rainfrac = 0.5,
  shore = 0,
  tides = 0,
  deepsoil = NA,
  lamb = 0,
  IUV = 0,
  soilgrids = 0,
  IR = 0,
  message = 0,
  fail = nyears * 24 * 365,
  TAI = 0,
  windfac = 1,
  snowcond = 0,
  intercept = max(maxshade) / 100 * 0.3,
  grasshade = 0,
  terra_source = NA,
  scenario = 0
) {

  # function to assist with interpolated data in leap years
  leapfix <- function(indata, yearlist, mult = 1){
    leapyears <- seq(1900, 2060, 4)
    for(j in 1:length(yearlist)){
      if(yearlist[j] %in% leapyears){# add day for leap year if needed
        if(mult == 1){
          data <- c(indata[1:59], indata[59], indata[60:365])
        }else{
          data <- c(indata[1:(59 * mult)], indata[(58*mult+1):(59 * mult)], indata[(59 * mult + 1):(365 * mult)])
        }
      }else{
        data <- indata
      }
      if(j == 1){
        alldata <- data
      }else{
        alldata <- c(alldata, data)
      }
    }
    return(alldata)
  }

  errors <- 0

  # error trapping - originally inside the Fortran code, but now checking before executing Fortran
  if(ystart < 1958){
    cat("ERROR: TerraClimate climate data is not available prior to 1958", '\n')
    errors<-1
  }
  curdate <- Sys.time() - 60 * 60 * 24
  curyear <- as.numeric(format(curdate, "%Y"))
  if(ystart > curyear - 1){
    cat(paste0("warning: TerraClimate climate data is only available until ", curyear - 1, ")", '\n'))
  }
  if(DEP[2]-DEP[1]>3 | DEP[3]-DEP[2]>3){
    message("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if(DEP[2]-DEP[1]<2){
    cat("warning, nodes might be too close near the surface, try a different spacing if the program is crashing \n")
  }
  if(DEP[10] != 200){
    cat("warning, last depth in soil should not be changed from 200 without good reason \n")
  }
  if(!(timeinterval %in% c(12, 365))){
    message("ERROR: the variable 'timeinterval' is out of bounds.
        Please enter a correct value (12 or 365).", '\n')
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
  if(run.gads == 1){
    message("If program is crashing, try run.gads = 2.", '\n')
  }
  if(run.gads%in%c(0,1,2)==FALSE){
    message("ERROR: the variable 'run.gads' be either 0, 1 or 2.
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
  if(min(Thcond)<0){
    cat("ERROR: Thermal variable conductivity (THCOND) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(min(Density)<0){
    cat("ERROR: Density variable (Density) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(min(SpecHeat)<0){
    cat("ERROR: Specific heat variable (SpecHeat) is negative.
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
  if(SLE<0.05 | SLE > 1){
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
  if(max(minshade-maxshade) >= 0){
    cat("ERROR: Your value(s) for minimum shade (minshade) is greater than or equal to the maximum shade (maxshade).
        Please correct this.", '\n')
    errors<-1
  }
  if(max(minshade)>100 | min(minshade)<0){
    cat("ERROR: Your value(s) for minimum shade (minshade) is out of bounds.
        Please input a value between 0 and 100.", '\n')
    errors<-1
  }
  if(max(maxshade)>100 | min(maxshade)<0){
    cat("ERROR: Your value(s) for maximum shade (maxshade) is out of bounds.
        Please input a value between 0 and 100.", '\n')
    errors<-1
  }
  if(!(scenario %in% c(0, 2, 4))){
    cat("ERROR: Scenario can only be 0, 2, or 4 corresponding to the TerraClimate climate change scenarios", '\n')
    errors<-1
  }
  if(scenario > 0 & (yfinish > 2015 | ystart < 1985)){
    cat("ERROR: TerraClimate climate change scenarios are only for years 1985 to 2015", '\n')
    errors<-1
  }
  # end error trapping

  if(errors==0){ # continue

    ################## time related variables #################################
    tzone <- paste("Etc/GMT+", 10, sep = "")
    dates <- seq(ISOdate(ystart, 1, 1, tz = tzone) - 3600 * 12, ISOdate((ystart + nyears), 1, 1, tz = tzone) - 3600 * 13, by = "days")
    if(timeinterval == 365){
      ndays <- length(dates)
    }else{
      ndays <- 12 * nyears
    }
    doys12 <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month
    microdaily <- 1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in
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
    if(timeinterval == 365){
      leapyears <- seq(1900, 2060, 4)
      for(mm in 1:nyears){
        if(mm == 1){
          currenty <- ystart
        }else{
          currenty <- ystart + mm - 1
        }
        if(currenty %in% leapyears){
          dayoy <- seq(1, 366)
        }else{
          dayoy <- seq(1, 365)
        }
        if(mm == 1){
          doy <- dayoy
        }else{
          doy <- c(doy, dayoy)
        }
      }
    }else{
      doy <- rep(doys12, nyears)
    }

    if(timeinterval<365){
      microdaily<-0 # run microclimate model as normal, where each day is iterated 3 times starting with the initial condition of uniform soil temp at mean monthly temperature
    }else{
      microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in
    }

    # now check if doing something other than middle day of each month, and create appropriate vector of Day of Year
    ndays<-length(doy)
    ida<-ndays

    idayst <- 1 # start day
    yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)

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

    SLES <- rep(SLE,ndays)
    SoilMoist <- SoilMoist_Init

    if(soilgrids == 1){
      cat('extracting data from SoilGrids \n')
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("package 'json' is needed to extract data from SoilGrids, please install it.",
             call. = FALSE)
      }
      require(jsonlite)
      #ov <- fromJSON(paste0('https://rest.isric.org/query?lon=',x[1],'&lat=',x[2],',&attributes=BLDFIE,SLTPPT,SNDPPT,CLYPPT'), flatten = TRUE)
      ov <- jsonlite::fromJSON(paste0('https://rest.isric.org/soilgrids/v2.0/properties/query?lon=',x[1],'&lat=',x[2],'&property=bdod&property=silt&property=clay&property=sand'), flatten = TRUE)
      if(length(ov) > 3){
        soilpro <- cbind(c(0, 5, 15, 30, 60, 100), unlist(ov$properties$layers$depths[[1]]$values.mean) / 100, unlist(ov$properties$layers$depths[[2]]$values.mean) / 10, unlist(ov$properties$layers$depths[[4]]$values.mean) / 10, unlist(ov$properties$layers$depths[[3]]$values.mean) / 10)
        soilpro <- rbind(soilpro, soilpro[6, ])
        soilpro[7, 1] <- 200
        #soilpro <- cbind(c(0,5,15,30,60,100,200), unlist(ov$properties$BLDFIE$M)/1000, unlist(ov$properties$CLYPPT$M), unlist(ov$properties$SLTPPT$M), unlist(ov$properties$SNDPPT$M) )
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

    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("package 'raster' is needed. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("package 'ncdf4' is needed. Please install it.",
           call. = FALSE)
    }
    library(ncdf4)
    require(utils)
    retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=100, sleep=30) {
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
    if(is.na(terra_source)){
     var <- "tmax"
     baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc#fillmismatch")

    nc <- retry(nc_open(baseurlagg))
    lon <- retry(ncvar_get(nc, "lon"))
    lat <- retry(ncvar_get(nc, "lat"))
    flat <- match(abs(lat - x[2]) < 1/48, 1)
    latindex <- which(flat %in% 1)
    if(length(latindex) == 0){
      flat <- match(abs(lat - x[2]) < 1/47.9, 1)
      latindex <- which(flat %in% 1)[1]
    }
    flon <- match(abs(lon - x[1]) < 1/48, 1)
    lonindex <- which(flon %in% 1)
    if(length(lonindex) == 0){
      flon <- match(abs(lon - x[1]) < 1/47.9, 1)
      lonindex <- which(flon %in% 1)[1]
    }
    time <- retry(ncvar_get(nc, "time"))
    alldays <- head(seq(as.Date('1900-01-01'), as.Date(paste0(yfinish + 1,'-01-01')), 'days'), -1)
    days1900 <- seq(1:length(alldays))
    allmonths <- head(seq(as.Date('1900-02-01'), as.Date(paste0(yfinish + 1,'-02-01')), '1 month'), -1) - 1
    allmonth.days <- head(days1900[which(alldays %in% allmonths & alldays >= as.Date('1957-12-01'))], -1)
    month.dates.to.do <- days1900[which(alldays %in% allmonths & alldays > as.Date(paste0(ystart - 1, '-12-01')) & alldays < as.Date(paste0(yfinish, '-12-31')))]
    timeindex <- which(time %in% month.dates.to.do)
    start <- c(lonindex, latindex, timeindex[1])
    count <- c(1, 1, length(timeindex))
    if(scenario == 0){
      message('extracting maximum air temperature data from TerraClimate\n')
      TMAXX <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
      var <- 'tmin'
      message('extracting minimum air temperature data from TerraClimate \n')
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      TMINN <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
      message('extracting precipitation data from TerraClimate \n')
      var <- 'ppt'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      RAINFALL <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
      message('extracting wind speed data from TerraClimate \n')
      var <- 'ws'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      WIND <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
      message('extracting vapour pressure deficit data from TerraClimate \n')
      var <- 'vpd'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      VPD <- retry(as.numeric(ncvar_get(nc, varid = var,start = start, count)))
      nc_close(nc)
      message('extracting solar radiation data from TerraClimate \n')
      var <- 'srad'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      SRAD <- retry(as.numeric(ncvar_get(nc, varid = var,start = start, count)))
      if(runmoist == 0){
        # extract soil moisture
        #soilmoisture <- suppressWarnings(raster::brick(paste(folder, "/soilw.mon.ltm.v2.nc", sep = "")))
        message("extracting soil moisture data from TerraClimate")
        #SoilMoist <- raster::extract(soilmoisture, x) / 1000 # this is originally in mm/m
        var <- 'soil'
        baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        SoilMoist <- retry(as.numeric(ncvar_get(nc, varid = var,start = start, count))) / 1000 # this is originally in mm/m
        nc_close(nc)
      }
    }else{
      if(scenario == 2){
        base <- '/thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data_plus2C/TerraClimate_2c'
      }else{
        base <- '/thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data_plus4C/TerraClimate_4c'
      }
      start <- c(lonindex, latindex, 1)
      count <- c(1, 1, -1)
      for(i in 1:length(yearlist)){
        var <- "tmax"
        baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' maximum air temperature data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          TMAXX <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          TMAXX <- c(TMAXX, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        var <- "tmin"
        baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' minimum air temperature data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          TMINN <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          TMINN <- c(TMINN, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        var <- "ppt"
        baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' precipitation data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          RAINFALL <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          RAINFALL <- c(RAINFALL, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        var <- "vpd"
        baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' vapour pressure deficit data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          VPD <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          VPD <- c(VPD, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        var <- "srad"
        baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data_plus2C/TerraClimate_2c_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' solar radiation data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          SRAD <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          SRAD <- c(SRAD, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        if(runmoist == 0){
          var <- "soil"
          baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
          nc <- retry(nc_open(baseurlagg))
          message(paste0('extracting plus ', scenario,' soil moisture data from TerraClimate for ', yearlist[i], '\n'))
          if(i == 1){
            SoilMoist <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))) / 1000 # this is originally in mm/m
          }else{
            SoilMoist <- c(SoilMoist, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))) / 1000)# this is originally in mm/m
          }
          nc_close(nc)
        }
      }
      timeindex <- which(time %in% month.dates.to.do)
      start <- c(lonindex, latindex, timeindex[1])
      count <- c(1, 1, length(timeindex))
      message('extracting wind speed data from TerraClimate \n')
      var <- 'ws'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      WIND <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
    }
    }else{
      alldays <- head(seq(as.Date('1900-01-01'), as.Date(paste0(yfinish + 1,'-01-01')), 'days'), -1)
      days1900 <- seq(1:length(alldays))
      allmonths <- head(seq(as.Date('1900-02-01'), as.Date(paste0(yfinish + 1,'-02-01')), '1 month'), -1) - 1
      allmonth.days <- head(days1900[which(alldays %in% allmonths & alldays >= as.Date('1957-12-01'))], -1)
      month.dates.to.do <- days1900[which(alldays %in% allmonths & alldays > as.Date(paste0(ystart - 1, '-12-01')) & alldays < as.Date(paste0(yfinish, '-12-31')))]
      terra <- as.data.frame(get_terra(x = loc, ystart = ystart, yfinish = yfinish, scenario = 0, source = terra_source))
      TMINN <- terra$TMINN
      TMAXX <- terra$TMAXX
      RAINFALL <- terra$RAINFALL
      VPD <- terra$VPD
      SRAD <- terra$SRAD
      SoilMoist <- terra$SoilMoist
      WIND <- terra$WIND
      if(scenario == 4){
        terra_cc <- as.data.frame(get_terra(x = loc, ystart = ystart, yfinish = yfinish, scenario = 4, source = terra_source))
        TMINN <- terra_cc$TMINN
        TMAXX <- terra_cc$TMAXX
        RAINFALL <- terra_cc$RAINFALL
        VPD <- terra_cc$VPD
        SRAD <- terra_cc$SRAD
        SoilMoist <- terra_cc$SoilMoist
      }
      if(scenario == 2){
        terra_cc <- as.data.frame(get_terra(x = loc, ystart = ystart, yfinish = yfinish, scenario = 2, source = terra_source))
        TMINN <- terra_cc$TMINN
        TMAXX <- terra_cc$TMAXX
        RAINFALL <- terra_cc$RAINFALL
        VPD <- terra_cc$VPD
        SRAD <- terra_cc$SRAD
        SoilMoist <- terra_cc$SoilMoist
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


    global_climate <- raster::brick(paste0(folder, "/global_climate.nc"))
    CLIMATE <- raster::extract(global_climate, x)
    #ALTT <- as.numeric(CLIMATE[, 1])
    library(dismo)
    elevation <- getData("worldclim", var = "alt", res = 2.5)
    ALTT <- raster::extract(elevation, x)
    if(is.na(ALTT)){
      ALTT <- 0
    }
    if(terrain == 1){
      elevatr <- 1
    }
    if(is.na(elev) & elevatr == 1){
      require(microclima)
      require(raster)
      cat('downloading DEM via package elevatr \n')
      dem <- microclima::get_dem(lat = loc[2], long = loc[1]) # mercator equal area projection
      xy <- data.frame(x = loc[1], y = loc[2])
      coordinates(xy) = ~x + y
      proj4string(xy) = "+init=epsg:4326"
      xy <- as.data.frame(spTransform(xy, crs(dem)))
      elev <- raster::extract(dem, xy)[1]
      if(terrain == 1){
        cat('computing slope, aspect and horizon angles \n')
        slope <- terrain(dem, unit = "degrees")
        slope <- raster::extract(slope, xy)
        aspect <- terrain(dem, opt = "aspect", unit = "degrees")
        aspect <- raster::extract(aspect, xy)
        ha24 <- 0
        for (i in 0:23) {
          har <- horizonangle(dem, i * 10, res(dem)[1])
          ha24[i + 1] <- atan(raster::extract(har, xy)) * (180/pi)
        }
        hori <- ha24
      }
    }
    hori<-as.matrix(hori) #horizon angles
    VIEWF <- 1-sum(sin(as.data.frame(hori) * pi / 180)) / length(hori) # convert horizon angles to radians and calc view factor(s)

    require(RNetCDF)
    ALTITUDES <- elev
    dbrow <- 1
    if(is.na(max(SoilMoist, ALTT, TMAXX)) == TRUE){
      message("Sorry, there is no environmental data for this location")
      SoilMoist <- raster::extract(soilmoisture, cbind(140, -35)) / 1000 # this is originally in mm/m
    }
    delta_elev <- 0
    if(is.na(elev) == FALSE){ # check if user-specified elevation
      delta_elev <- ALTT - elev # get delta for lapse rate correction
      ALTT <- elev # now make final elevation the user-specified one
    }
    adiab_corr_max <- delta_elev * lapse_max
    adiab_corr_min <- delta_elev * lapse_min

    if(is.na(RAINFALL[1])){
      cat("no climate data for this site, using dummy data so solar is still produced \n")
      CLIMATE <- raster::extract(global_climate, cbind(140, -35))
      CLIMATE[2:97] <- 0
      ALTT<-as.numeric(CLIMATE[, 1])
      delta_elev <- 0
      if(is.na(elev) == FALSE){ # check if user-specified elevation
        delta_elev <- ALTT - elev # get delta for lapse rate correction
        ALTT <- elev # now make final elevation the user-specified one
      }
      adiab_corr_max <- delta_elev * lapse_max
      adiab_corr_min <- delta_elev * lapse_min
      RAINFALL <- CLIMATE[, 2:13] * 0
      #stop()
    }
    RAINYDAYS <- CLIMATE[, 14:25] / 10
    # WNMAXX <- CLIMATE[, 26:37] / 10 * windfac
    # WNMINN<-WNMAXX * 0.1 # impose diurnal cycle
    # TMINN <- CLIMATE[, 38:49] / 10
    # TMAXX <- CLIMATE[, 50:61] / 10
    WNMAXX <- WIND
    WNMINN<-WNMAXX * 0.1 # impose diurnal cycle
    TMAXX <- TMAXX + adiab_corr_max
    TMINN <- TMINN + adiab_corr_min
    ALLMINTEMPS <- TMINN
    ALLMAXTEMPS <- TMAXX
    ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)
    MEANTEMPS <- (TMAXX + TMINN) / 2
    e <- WETAIR(db = MEANTEMPS, rh = 100)$e - VPD * 1000
    es <- WETAIR(db = TMAXX, rh = 100)$esat
    RHMINN <- (e / es) * 100
    es <- WETAIR(db = TMINN, rh = 100)$esat
    RHMAXX <- (e / es) * 100
    RHMINN[RHMINN>100]<-100
    RHMINN[RHMINN<0]<-0.01
    RHMAXX[RHMAXX>100]<-100
    RHMAXX[RHMAXX<0]<-0.01

    cat("running micro_global to get clear sky solar \n")
    if(run.gads == 0){
      TAI <- c(0.0670358341290886, 0.0662612704779235, 0.065497075238002, 0.0647431301168489, 0.0639993178022531, 0.0632655219571553, 0.0625416272145492, 0.0611230843885423, 0.0597427855962549, 0.0583998423063099, 0.0570933810229656, 0.0558225431259535, 0.0545864847111214, 0.0533843764318805, 0.0522154033414562, 0.0499736739981675, 0.047855059159556, 0.0458535417401334, 0.0439633201842001, 0.0421788036108921, 0.0404946070106968, 0.0389055464934382, 0.0374066345877315, 0.0359930755919066, 0.0346602609764008, 0.0334037648376212, 0.0322193394032758, 0.0311029105891739, 0.0300505736074963, 0.0290585886265337, 0.0281233764818952, 0.0272415144391857, 0.0264097320081524, 0.0256249068083005, 0.0248840604859789, 0.0241843546829336, 0.0235230870563317, 0.0228976873502544, 0.0223057135186581, 0.0217448478998064, 0.0212128934421699, 0.0207077699817964, 0.0202275105711489, 0.0197702578594144, 0.0193342605242809, 0.0189178697551836, 0.0177713140039894, 0.0174187914242432, 0.0170790495503944, 0.0167509836728154, 0.0164335684174899, 0.0161258546410128, 0.0158269663770596, 0.0155360978343254, 0.0152525104459325, 0.0149755299703076, 0.0147045436435285, 0.0144389973831391, 0.0141783930434343, 0.0134220329447663, 0.0131772403830191, 0.0129356456025128, 0.0126970313213065, 0.0124612184223418, 0.0122280636204822, 0.01199745718102, 0.0115436048739351, 0.0110993711778668, 0.0108808815754663, 0.0106648652077878, 0.0104513876347606, 0.0102405315676965, 0.00982708969547694, 0.00962473896278535, 0.00903679230300494, 0.00884767454432418, 0.0083031278398166, 0.00796072474935954, 0.00755817587626185, 0.00718610751850881, 0.00704629977586921, 0.00684663903049612, 0.00654155580333479, 0.00642947339729728, 0.00627223096874308, 0.00603955966866779, 0.00580920937536261, 0.00568506186880564, 0.00563167068287251, 0.00556222005081865, 0.00550522989971023, 0.00547395763028062, 0.0054478983436216, 0.00541823364504573, 0.00539532163908382, 0.00539239864119488, 0.00541690124712384, 0.00551525885358836, 0.00564825853509463, 0.00577220185074264, 0.00584222986640171, 0.00581645238345584, 0.00566088137411449, 0.00535516862329704, 0.00489914757707667, 0.00432017939770409, 0.0036813032251836, 0.00309019064543606, 0.00270890436501562, 0.00276446109239711, 0.00356019862584603)
    }else{
      TAI <- 0
    }
    micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, TAI = TAI, timeinterval = 365, solonly = 1)
    clearskyrad <- micro_clearsky$metout[, c(1, 13)][, 2]
    clearskymean <- leapfix(clearskyrad, yearlist, 24)

    dates2 <- head(seq(as.Date(paste0(ystart, '-01-01')), as.Date(paste0(yfinish + 1, '-01-01')), by = "days"), -1)
    dates <- head(seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'Etc/GMT+10'), as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y ", tz = 'Etc/GMT+10'), by = 'hours'), -1)
    ndays <- length(dates2)

    clearskymean <- aggregate(clearskymean, by = list(format(dates, '%Y-%m')), FUN = mean)[, 2]

    #allclearsky <- allclearsky[1:ndays]
    # convert from W/d to MJ/d
    #allclearsky <- allclearsky * 3600 / 1e6
    cloud <- (1 - SRAD / clearskymean) * 100
    cloud[cloud < 0] <- 0
    cloud[cloud > 100] <- 100
    CCMINN <- cloud * 0.5
    CCMAXX <- cloud * 2
    CCMINN[CCMINN > 100] <- 100
    CCMAXX[CCMAXX > 100] <- 100
    CCMINN[CCMINN < 0] <- 0
    CCMAXX[CCMAXX < 0] <- 0
    #CCMINN <- CLIMATE[, 86:97] / 10
    if(clearsky == 1){
      CCMINN <- CCMINN * 0
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
    WNMINN <- WNMINN * (2 / 10) ^ 0.15
    WNMAXX <- WNMAXX * (2 / 10) ^ 0.15
    # impose uniform warming
    TMAXX <- TMAXX + warm
    TMINN <- TMINN + warm
    methspline <- 'periodic'
    if(timeinterval != 12){ # spline from 12 days to chosen time interval
      for(i in 1:nyears){
        if(yearlist[i] %in% seq(1900, 2060, 4)){
          xmax <- 366
        }else{
          xmax <- 365
        }
        if(i == 1){
          start <- 1
          end <- 12
          TMAXX1 <-suppressWarnings(spline(doys12,TMAXX[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          TMAXX2<-TMAXX1$y
          TMINN1 <-suppressWarnings(spline(doys12,TMINN[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          TMINN2 <- TMINN1$y
          RHMAXX1 <-suppressWarnings(spline(doys12,RHMAXX[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          RHMAXX2 <- RHMAXX1$y
          RHMINN1 <-suppressWarnings(spline(doys12,RHMINN[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          RHMINN2 <- RHMINN1$y
          CCMAXX1 <-suppressWarnings(spline(doys12,CCMAXX[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          CCMAXX2 <- CCMAXX1$y
          CCMINN1 <-suppressWarnings(spline(doys12,CCMINN[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          CCMINN2 <- CCMINN1$y
          WNMAXX1 <-suppressWarnings(spline(doys12,WNMAXX[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          WNMAXX2 <- WNMAXX1$y
          WNMINN1 <-suppressWarnings(spline(doys12,WNMINN[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          WNMINN2 <- WNMINN1$y
          if(runmoist==0){
            SoilMoist1 <- suppressWarnings(spline(doys12,SoilMoist[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
            SoilMoist2 <- SoilMoist1$y
          }
        }else{
          start <- end + 1
          end <- end + 12
          TMAXX1 <-suppressWarnings(spline(c(0, doys12), c(tail(TMAXX2, 1), TMAXX[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          TMAXX2<-c(TMAXX2, TMAXX1$y)
          TMINN1 <-suppressWarnings(spline(c(0, doys12), c(tail(TMINN2, 1), TMINN[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          TMINN2<-c(TMINN2, TMINN1$y)
          RHMAXX1 <-suppressWarnings(spline(c(0, doys12), c(tail(RHMAXX2, 1), RHMAXX[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          RHMAXX2<-c(RHMAXX2, RHMAXX1$y)
          RHMINN1 <-suppressWarnings(spline(c(0, doys12), c(tail(RHMINN2, 1), RHMINN[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          RHMINN2<-c(RHMINN2, RHMINN1$y)
          CCMAXX1 <-suppressWarnings(spline(c(0, doys12), c(tail(CCMAXX2, 1), CCMAXX[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          CCMAXX2<-c(CCMAXX2, CCMAXX1$y)
          CCMINN1 <-suppressWarnings(spline(c(0, doys12), c(tail(CCMINN2, 1), CCMINN[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          CCMINN2<-c(CCMINN2, CCMINN1$y)
          WNMAXX1 <-suppressWarnings(spline(c(0, doys12), c(tail(WNMAXX2, 1), WNMAXX[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          WNMAXX2<-c(WNMAXX2, WNMAXX1$y)
          WNMINN1 <-suppressWarnings(spline(c(0, doys12), c(tail(WNMINN2, 1), WNMINN[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          WNMINN2<-c(WNMINN2, WNMINN1$y)
          if(runmoist==0){
            SoilMoist1 <-suppressWarnings(spline(c(0, doys12), c(tail(SoilMoist2, 1), SoilMoist[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
            SoilMoist2<-c(SoilMoist2, SoilMoist1$y)
          }
        }
      }
      TMAXX <- TMAXX2
      TMINN <- TMINN2
      RHMAXX <- RHMAXX2
      RHMINN <- RHMINN2
      WNMAXX <- WNMAXX2
      WNMINN <- WNMINN2
      CCMAXX <- CCMAXX2
      CCMINN <- CCMINN2
      if(runmoist==0){
        SoilMoist <- SoilMoist2
      }
    }
    CCMINN[CCMINN > 100] <- 100
    CCMAXX[CCMAXX > 100] <- 100
    CCMINN[CCMINN < 0] <- 0
    CCMAXX[CCMAXX < 0] <- 0

    orig.RAINFALL <- RAINFALL

    # get annual mean temp for creating deep soil (2m) boundary condition
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

    if(timeinterval == 12){
      if(is.na(deepsoil)){
        if(nyears == 1){
          avetemp <- (sum(TMAXX) + sum(TMINN)) / (length(TMAXX)*2)
          deepsoil <-rep(avetemp, length(TMAXX))
        }else{
          avetemp <- rowMeans(cbind(TMAXX, TMINN), na.rm=TRUE)
          deepsoil <- raster::movingFun(avetemp, n = 12, fun = mean, type = 'to')
          yearone <- rep((sum(TMAXX[1:12]) + sum(TMINN[1:12])) / (12 * 2), 12)
          deepsoil[1:12] <- yearone
        }
      }else{
        tannul <- mean(deepsoil)
      }
    }else{
      if(is.na(deepsoil)){
        if(nyears == 1){
          avetemp <- (sum(TMAXX) + sum(TMINN)) / (length(TMAXX)*2)
          deepsoil <-rep(avetemp, ndays)
        }else{
          avetemp <- rowMeans(cbind(TMAXX, TMINN), na.rm=TRUE)
          if(length(TMAXX) < 365){
            deepsoil <- rep((sum(TMAXX) + sum(TMINN)) / (length(TMAXX) * 2), length(TMAXX))
          }else{
            deepsoil <- raster::movingFun(avetemp, n = 365, fun = mean, type = 'to')
            yearone <- rep((sum(TMAXX[1:365]) + sum(TMINN[1:365])) / (365 * 2), 365)
            deepsoil[1:365] <- yearone
          }
        }
      }else{
        tannul <- mean(deepsoil)
      }
    }
    tannulrun <- deepsoil
    daymon<-c(31,28,31,30,31,30,31,31,30,31,30,31) # days in each month
    bb <- 0
    # if doing daily sims, spread rainfall evenly across days based on mean monthly rainfall and the number of rainy days per month
    if(timeinterval==365){
      for(j in 1:length(yearlist)){
        if(yearlist[j] %in% seq(1900, 2060, 4)){
          dys <- 366
          daymon2 <- c(31,29,31,30,31,30,31,31,30,31,30,31) # days in each month
        }else{
          dys <- 365
          daymon2 <- daymon
        }
        RAINFALL1<-1:dys
        sort<-matrix(data = 0,nrow = dys,ncol = 2)
        m<-1
        b<-0
        for (i in 1:12){ #begin loop through 12 months of year
          ndays=daymon2[i]
          for (k in 1:ndays){
            b<-b+1
            sort[m,1]<-i
            sort[m,2]<-b
            if(k<=RAINYDAYS[i] & rainfrac>0){
              if(k==1){
                RAINFALL1[m]<-RAINFALL[i + (j - 1) * 12]*rainfrac*rainmult # if first day of month, make user-specified fraction of monthly rainfall fall on first day
              }else{
                RAINFALL1[m]<-(RAINFALL[i + (j - 1) * 12]*(1-rainfrac)*rainmult)/RAINYDAYS[i] # make remaining rain fall evenly over the remaining number of rainy days for the month, starting at the beginning of the month
              }
            }else{
              if(rainfrac==0){
                RAINFALL1[m]<-(RAINFALL[i + (j - 1) * 12]*rainmult)/RAINYDAYS[i]
              }else{
                RAINFALL1[m]<-0
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
        if(j == 1){
          RAINFALL3<-RAINFALL2$RAINFALL1
          if(TMINN[1]<snowtemp){
            RAINFALL3[1]<-0 # this is needed in some cases to allow the integrator to get started
          }
        }else{
          RAINFALL3 <- c(RAINFALL3, RAINFALL2$RAINFALL1)
        }
      }
    }else{
      if(timeinterval!=12){
        RAINFALL3<-rep(rep(sum(RAINFALL)/timeinterval,timeinterval),nyears) # just spread evenly across every day
      }else{ # running middle day of each month - divide monthly rain by number of days in month
        RAINFALL3<-RAINFALL/rep(daymon,nyears)
      }
    }#end check doing daily sims
    RAINFALL <- RAINFALL3

    SOLRhr<-rep(0,24*ndays)

    hourly <- 0
    if(microclima > 0){
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
      if(!is.raster(dem)){
        dem <- microclima::get_dem(r = NA, lat = lat, long = long, resolution = 100, zmin = -20)
      }
      xy <- data.frame(x = long, y = lat)
      coordinates(xy) = ~x + y
      proj4string(xy) = "+init=epsg:4326"
      xy <- as.data.frame(spTransform(xy, crs(dem)))
      if (class(slope) == "logical") {
        slope <- terrain(dem, unit = "degrees")
        slope <- raster::extract(slope, xy)
      }
      if (class(aspect) == "logical") {
        aspect <- terrain(dem, opt = "aspect", unit = "degrees")
        aspect <- raster::extract(aspect, xy)
      }
      ha <- 0
      ha36 <- 0
      for (i in 0:35) {
        har <- horizonangle(dem, i * 10, res(dem)[1])
        ha36[i + 1] <- atan(raster::extract(har, xy)) * (180/pi)
      }
      for (i in 1:length(hour.microclima)) {
        saz <- solazi(hour.microclima[i], lat, long, jd[i], merid = long)
        saz <- round(saz/10, 0) + 1
        saz <- ifelse(saz > 36, 1, saz)
        ha[i] <- ha36[saz]
      }

      methspline <- 'periodic'
      for(i in 1:nyears){
        if(yearlist[i] %in% seq(1900, 2060, 4)){
          xmax <- 366
        }else{
          xmax <- 365
        }
        if(i == 1){
          start <- 1
          end <- 12
          cloud1 <-suppressWarnings(spline(doys12,cloud[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          cloud2<-cloud1$y
        }else{
          start <- end + 1
          end <- end + 12
          cloud1 <-suppressWarnings(spline(c(0, doys12), c(tail(cloud2, 1), cloud[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          cloud2<-c(cloud2, cloud1$y)
        }
      }
      cloudhr <- cbind(rep(seq(1, length(cloud2)),24), rep(cloud2, 24))
      cloudhr <- cloudhr[order(cloudhr[,1]),]
      cloudhr <- cloudhr[,2]
      #cloudhr <- leapfix(cloudhr, yearlist, 24)
      dsw2 <- leapfix(clearskyrad, yearlist, 24) *(0.36+0.64*(1-cloudhr/100)) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
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
      diffuse_frac_all <- h_dif / (h_dni + h_dif) # calculated diffuse fraction
      diffuse_frac_all[is.na(diffuse_frac_all)] <- 1
      radwind2 <- .shortwave.ts(h_dni * 0.0036, h_dif * 0.0036, jd, hour.microclima, lat, long, slope, aspect, ha = ha, svv = 1, x = microclima.LOR, l = mean(microclima.LAI), albr = 0, merid = long, dst = 0, difani = FALSE)
      #microclima.out$hourlyradwind <- radwind2
      SOLRhr_all <- radwind2$swrad / 0.0036
      if(microclima == 2){ # use hourly solar from microclima
        hourly <- 2
        VIEWF <- 1 # accounted for already in microclima cals
        hori <- rep(0, 24) # accounted for already in microclima calcs
      }
    }else{
      diffuse_frac <- NA
    }
    if(timeinterval == 12 & microclima > 0){
      dates_all <- head(seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y ", tz = 'UTC'), by = 'hours'), -1)
      month.dates.to.do2 <- days1900[which(alldays %in% allmonths & alldays > as.Date(paste0(ystart, '-01-01')) & alldays < as.Date(paste0(yfinish + 1, '-01-01')))]
      mon <- alldays[month.dates.to.do2]
      dates3 <- dates_all[which(format(dates_all, "%Y-%m-%d") %in% as.character(mon))]
      SOLRhr <- SOLRhr_all[which(dates_all %in% dates3)]
      diffuse_frac <- diffuse_frac_all[which(dates_all %in% dates3)]
    }

    ndays<-length(RAINFALL)
    if(length(TAI) < 111){ # no user supplied values, compute with GADS
      if(run.gads > 0){
        ####### get solar attenuation due to aerosols with program GADS #####################
        relhum <- 1
        if(run.gads == 1){ # fortran version
          optdep.summer <- as.data.frame(rungads(longlat[2], longlat[1], relhum, 0))
          optdep.winter <- as.data.frame(rungads(longlat[2], longlat[1], relhum, 1))
        }else{ # r version
          optdep.summer <- as.data.frame(gads.r(longlat[2], longlat[1], relhum, 0))
          optdep.winter <- as.data.frame(gads.r(longlat[2], longlat[1], relhum, 1))
        }
        optdep <-cbind(optdep.winter[,1],rowMeans(cbind(optdep.summer[,2],optdep.winter[,2])))
        optdep <-as.data.frame(optdep)
        colnames(optdep) <- c("LAMBDA", "OPTDEPTH")
        a <-lm(OPTDEPTH ~ poly(LAMBDA, 6, raw = TRUE), data = optdep)
        LAMBDA <- c(290, 295, 300, 305, 310, 315, 320, 330, 340, 350, 360, 370, 380, 390, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1080, 1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280, 1300, 1320, 1380, 1400, 1420, 1440, 1460, 1480, 1500, 1540, 1580, 1600, 1620, 1640, 1660, 1700, 1720, 1780, 1800, 1860, 1900, 1950, 2000, 2020, 2050, 2100, 2120, 2150, 2200, 2260, 2300, 2320, 2350, 2380, 2400, 2420, 2450, 2490, 2500, 2600, 2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000)
        TAI <- predict(a, data.frame(LAMBDA))
        ################ end GADS ##################################################
      }else{ # use the original profile from Elterman, L. 1970. Vertical-attenuation model with eight surface meteorological ranges 2 to 13 kilometers. U. S. Airforce Cambridge Research Laboratory, Bedford, Mass.
        TAI <-c(0.42, 0.415, 0.412, 0.408, 0.404, 0.4, 0.395, 0.388, 0.379, 0.379, 0.379, 0.375, 0.365, 0.345, 0.314, 0.3, 0.288, 0.28, 0.273, 0.264, 0.258, 0.253, 0.248, 0.243, 0.236, 0.232, 0.227, 0.223, 0.217, 0.213, 0.21, 0.208, 0.205, 0.202, 0.201, 0.198, 0.195, 0.193, 0.191, 0.19, 0.188, 0.186, 0.184, 0.183, 0.182, 0.181, 0.178, 0.177, 0.176, 0.175, 0.175, 0.174, 0.173, 0.172, 0.171, 0.17, 0.169, 0.168, 0.167, 0.164, 0.163, 0.163, 0.162, 0.161, 0.161, 0.16, 0.159, 0.157, 0.156, 0.156, 0.155, 0.154, 0.153, 0.152, 0.15, 0.149, 0.146, 0.145, 0.142, 0.14, 0.139, 0.137, 0.135, 0.135, 0.133, 0.132, 0.131, 0.13, 0.13, 0.129, 0.129, 0.128, 0.128, 0.128, 0.127, 0.127, 0.126, 0.125, 0.124, 0.123, 0.121, 0.118, 0.117, 0.115, 0.113, 0.11, 0.108, 0.107, 0.105, 0.103, 0.1)
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
    #hourly<-0
    rainhourly<-0
    TAIRhr<-rep(0,24*ndays)
    RHhr<-rep(0,24*ndays)
    WNhr<-rep(0,24*ndays)
    CLDhr<-rep(0,24*ndays)
    RAINhr<-rep(0,24*ndays)
    ZENhr<-rep(-1,24*ndays)
    IRDhr<-rep(-1,24*ndays)
    # microclimate input parameters list
    microinput<-c(ndays,RUF,ERR,Usrhyt,Refhyt,Numtyps,Z01,Z02,ZH1,ZH2,idayst,ida,HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,CMH2O,microdaily,tannul,EC,VIEWF,snowtemp,snowdens,snowmelt,undercatch,rainmult,runshade,runmoist,maxpool,evenrain,snowmodel,rainmelt,writecsv,densfun,hourly,rainhourly,lamb,IUV,RW,PC,RL,SP,R1,IM,MAXCOUNT,IR,message,fail,snowcond,intercept,grasshade,solonly,ZH,D0,TIMAXS,TIMINS,spinup)

    if(length(LAI)<ndays){
      LAI<-rep(LAI[1],ndays)
    }
    if(shore==0){
      tides<-matrix(data = 0., nrow = 24*ndays, ncol = 3) # make an empty matrix
    }
    # all microclimate data input list - all these variables are expected by the input argument of the fortran micro2014 subroutine
    micro<-list(tides=tides,microinput=microinput,doy=doy,SLES=SLES,DEP=DEP,Nodes=Nodes,MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TMAXX=TMAXX,TMINN=TMINN,RHMAXX=RHMAXX,RHMINN=RHMINN,CCMAXX=CCMAXX,CCMINN=CCMINN,WNMAXX=WNMAXX,WNMINN=WNMINN,TAIRhr=TAIRhr,RHhr=RHhr,WNhr=WNhr,CLDhr=CLDhr,SOLRhr=SOLRhr,RAINhr=RAINhr,ZENhr=ZENhr,IRDhr=IRDhr,REFLS=REFLS,PCTWET=PCTWET,soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists,RAINFALL=RAINFALL,tannulrun=tannulrun,PE=PE,KS=KS,BB=BB,BD=BD,DD=DD,L=L,LAI=LAI)

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
    if(timeinterval == 12){
      message(paste('running microclimate model for middle day of each month for',ystart,'to ', yfinish, 'at site',location,'\n'))
    }else{
      message(paste('running microclimate model, splining monthly data to daily, for',ystart,'to ', yfinish, 'at site',location,'\n'))
    }
    message('Note: the output column `SOLR` in metout and shadmet is for unshaded horizontal plane solar radiation \n')
    ptm <- proc.time() # Start timing
    microut<-microclimate(micro)
    message(paste0('runtime ', (proc.time() - ptm)[3], ' seconds')) # Stop the clock

    metout<-microut$metout # retrieve above ground microclimatic conditions, min shade
    shadmet<-microut$shadmet # retrieve above ground microclimatic conditions, max shade
    soil<-microut$soil # retrieve soil temperatures, minimum shade
    shadsoil<-microut$shadsoil # retrieve soil temperatures, maximum shade
    tcond <- microut$tcond
    shadtcond <- microut$shadtcond
    specheat <- microut$specheat
    shadspecheat <- microut$shadspecheat
    densit <- microut$densit
    shaddensit <- microut$shaddensit
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
    days <- rep(seq(1, ndays), 24)
    days <- days[order(days)]
    dates <- days + rep(seq(0, 1380, 60), ndays)  / 60 / 24 - 1 # dates for hourly output
    #dates2<-seq(1, timeinterval * nyears) # dates for daily output
    dates2 <- month.dates.to.do
    tzone <- paste("Etc/GMT+", 0, sep = "")
    dates3 <- head(seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y ", tz = 'UTC'), by = 'hours'), -1)
    if(timeinterval == 12){
      dates_all <- head(seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y ", tz = 'UTC'), by = 'hours'), -1)
      month.dates.to.do2 <- days1900[which(alldays %in% allmonths & alldays > as.Date(paste0(ystart, '-01-01')) & alldays < as.Date(paste0(yfinish + 1, '-01-01')))]
      mon <- alldays[month.dates.to.do2]
      dates3 <- dates_all[which(format(dates_all, "%Y-%m-%d") %in% as.character(mon))]
    }else{
      dates <- dates3
    }

    if(lamb == 1){
      drlam<-as.data.frame(microut$drlam) # retrieve direct solar irradiance
      drrlam<-as.data.frame(microut$drrlam) # retrieve direct Rayleigh component solar irradiance
      srlam<-as.data.frame(microut$srlam) # retrieve scattered solar irradiance
      if(snowmodel == 1){
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,tcond=tcond,shadtcond=shadtcond,specheat=specheat,shadspecheat=shadspecheat,densit=densit,shaddensit=shaddensit,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=MINSHADES,maxshade=MAXSHADES,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam,dates=dates,dates2=dates2,dates3=dates3,PE=PE,BD=BD,DD=DD,BB=BB,KS=KS,slope=slope,aspect=aspect,hori=hori,dem=dem, diffuse_frac = diffuse_frac))
      }else{
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,tcond=tcond,shadtcond=shadtcond,specheat=specheat,shadspecheat=shadspecheat,densit=densit,shaddensit=shaddensit,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=MINSHADES,maxshade=MAXSHADES,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam,dates=dates,dates2=dates2,dates3=dates3,PE=PE,BD=BD,DD=DD,BB=BB,KS=KS,slope=slope,aspect=aspect,hori=hori,dem=dem, diffuse_frac = diffuse_frac))
      }
    }else{
      if(snowmodel == 1){
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,tcond=tcond,shadtcond=shadtcond,specheat=specheat,shadspecheat=shadspecheat,densit=densit,shaddensit=shaddensit,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=MINSHADES,maxshade=MAXSHADES,DEP=DEP,dates=dates,dates2=dates2,dates3=dates3,PE=PE,BD=BD,DD=DD,BB=BB,KS=KS,slope=slope,aspect=aspect,hori=hori,dem=dem, diffuse_frac = diffuse_frac))
      }else{
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,tcond=tcond,shadtcond=shadtcond,specheat=specheat,shadspecheat=shadspecheat,densit=densit,shaddensit=shaddensit,RAINFALL=RAINFALL,ndays=ndays,elev=ALTT,REFL=REFL[1],longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=MINSHADES,maxshade=MAXSHADES,DEP=DEP,dates=dates,dates2=dates2,dates3=dates3,PE=PE,BD=BD,DD=DD,BB=BB,KS=KS,slope=slope,aspect=aspect,hori=hori,dem=dem, diffuse_frac = diffuse_frac))
      }
    }
  } # end error trapping
} # end of micro_terra function
