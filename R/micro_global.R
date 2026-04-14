#' Global implementation of the microclimate model
#'
#' An implementation of the NicheMapR microclimate model that uses the global climate database
#' derived from "New, M., Lister, D., Hulme, M. and Makin, I., 2002: A high-resolution data
#' set of surface climate over global land areas. Climate Research 21:1-25"
#' It also optionally uses a global monthly soil moisture estimate from NOAA CPC Soil Moisture http://140.172.38.100/psd/thredds/catalog/Datasets/cpcsoil/catalog.html
#' Aerosol attenuation can also be computed based on the Global Aerosol Data Set (GADS)
#' Koepke, P., M. Hess, I. Schult, and E. P. Shettle. 1997. Global Aerosol Data Set. Max-Planck-Institut for Meteorologie, Hamburg
#' by choosing the option 'run.gads <- 1' (Fortran version, quicker but may crash on some systems) or 'run.gads <- 2' (R version)
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param timeinterval The number of time intervals to generate predictions for over a year (must be 12 <= x <=365)
#' @param nyears The number of years to run
#' @param dem A digital elevation model used produced by microclima function 'get_dem' via R package 'elevatr' (internally generated via same function based on 'loc' if dem = NA and terrain/elevatr/microclima != 0)
#' @param REFL Soil solar reflectance, decimal \%
#' @param elev Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param soiltype Soil type: Rock = 0, sand = 1, loamy sand = 2, sandy loam = 3, loam = 4, silt loam = 5, sandy clay loam = 6, clay loam = 7, silt clay loam = 8, sandy clay = 9, silty clay = 10, clay = 11, user-defined = 12, based on Campbell and Norman 1990 Table 9.1.
#' @param minshade Minimum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param maxshade Maximum shade level to use (\%) (can be a single value or a vector of daily values)
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
#' @return tcond Hourly predictions of the soil thermal conductivity under the minimum specified shade
#' @return shadtcond Hourly predictions of the soil thermal conductivity under the maximum specified shade
#' @return specheat Hourly predictions of the soil specific heat capacity under the minimum specified shade
#' @return shadspecheat Hourly predictions of soil specific heat capacity under the maximum specified shade
#' @return densit Hourly predictions of the soil density under the minimum specified shade
#' @return shaddensit Hourly predictions of the soil density under the maximum specified shade
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
#' \code{ndmax}{ = 3, iterations of each day to get a steady periodic}\cr\cr
#' \code{Soil_Init}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{elevatr}{ = 0, Use elevatr package to get high resolution elevation for location? 1 = yes, 0 = no}\cr\cr
#' \code{terrain}{ = 0, Use elevatr package to adjust horizon angles, slope and aspect? 1 = yes, 0 = no}\cr\cr
#' \code{dem.res}{ = 30, Requested resolution of the DEM from elevatr, m}\cr\cr
#' \code{zmin minimum}{ = -20, elevation of DEM for terrain calculations, m (may need to be made negative if below sea level)}\cr\cr
#' \code{pixels}{ = 100, Number of pixels along one edge of square requested of DEM requested from elevatr}\cr\cr
#' \code{microclima}{ = 0, Use microclima and elevatr package to compute diffuse fraction of solar radiation (1) and adjust solar radiation for terrain (2)? 0 = no}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
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
#'
#' \strong{ Soil moisture mode parameters:}
#'
#' \code{runmoist}{ = 0, Run soil moisture model? 1=yes, 0=no  1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{PE}{ = rep(1.1,19), Air entry potential (J/kg) (19 values descending through soil for specified soil nodes in parameter}
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{KS}{ = rep(0.0037,19), Saturated conductivity, (kg s/m^3) (19 values descending through soil for specified soil nodes in parameter}
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
#' \code{L}{ = c(0,0,8.2,8.0,7.8,7.4,7.1,6.4,5.8,4.8,4.0,1.8,0.9,0.6,0.8,0.4,0.4,0,0)*10000, root density (m/m^3), (19 values descending through soil for specified soil nodes in parameter}\cr\cr
#' \code{R1}{ = 0.001, root radius, m}\cr\cr
#' \code{RW}{ = 2.5e+10, resistance per unit length of root, m^3 kg-1 s-1}\cr\cr
#' \code{RL}{ = 2e+6, resistance per unit length of leaf, m^3 kg-1 s-1}\cr\cr
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
#' \code{snowmodel}{ = 0, run the snow model 1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (Mg/m3), overridden by densfun}\cr\cr
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
#' \item 4 TAREF - air temperature (°C) at reference height (specified by 'Refhyt', 1.2m default)
#' \item 5 RHLOC - relative humidity (\%) at local height (specified by 'Usrhyt' variable)
#' \item 6 RH  - relative humidity (\%) at reference height (specified by 'Refhyt', 1.2m default)
#' \item 7 VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)
#' \item 8 VREF - wind speed (m/s) at reference height (specified by 'Refhyt', 1.2m default)
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
#'micro <- micro_global() # run the model with default location (Madison, Wisconsin) and settings
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
#'  if(i == 1){
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
#'  if(i == 1){
#'    plot(shadsoil[,i + 2] ~ micro$dates, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l", main = paste("soil temperature ", maxshade, "% shade", sep=""))
#'  }else{
#'    points(shadsoil[,i + 2] ~ micro$dates, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l")
#'  }
#'}
#' @export
micro_global <- function(
  loc = c(-89.4557, 43.1379),
  timeinterval = 12,
  nyears = 1,
  soiltype = 4,
  REFL = 0.15,
  elev = NA,
  slope = 0,
  aspect = 0,
  lapse_max = 0.0077,
  lapse_min = 0.0039,
  DEP = c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
  minshade = 0,
  maxshade = 90,
  dem = NA,
  dem.res = 30,
  zmin = -20,
  pixels = 100,
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
  microclima.dem.res = 100,
  microclima.zmin = -20,
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
  cap = 1,
  CMH2O = 1,
  hori = rep(0,24),
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
  LAI = 0.1,
  microclima.LAI = 0,
  microclima.LOR = 1,
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
  ndmax = 3,
  soilgrids = 0,
  IR = 0,
  message = 0,
  fail = nyears * 24 * 365,
  runmicro = 1,
  TAI = 0,
  warm = 0,
  windfac = 1,
  snowcond = 0,
  intercept = max(maxshade) / 100 * 0.3,
  grasshade = 0,
  maxsurf = 85
) {

  SoilMoist <- SoilMoist_Init
  # Reference height (m) at which input air temperature, wind speed and
  # relative humidity are measured
  Refhyt <- 1.2

  # --- 1. Input validation --------------------------------------------------
  # File-specific checks first (timeinterval, soiltype, loc bounds, run options)
  errors <- 0
  if (timeinterval < 12 | timeinterval > 365) {
    message("ERROR: 'timeinterval' is out of bounds. Please enter a value between 12 and 365.")
    errors <- 1
  }
  if (is.numeric(loc[1])) {
    if (loc[1] > 180 | loc[2] > 90) {
      message("ERROR: Latitude or longitude (loc) is out of bounds. Please enter correct values.")
      errors <- 1
    }
  }
  if (timezone %in% c(0, 1) == FALSE) {
    message("ERROR: 'timezone' must be either 0 or 1. Please correct.")
    errors <- 1
  }
  if (run.gads == 1) {
    message("If program is crashing, try run.gads = 2.")
  }
  if (write_input %in% c(0, 1) == FALSE) {
    message("ERROR: 'write_input' must be 0 or 1. Please correct.")
    errors <- 1
  }
  if (soiltype < 0 | soiltype > 11) {
    message("ERROR: soil type must range between 1 and 11. Please correct.")
    errors <- 1
  }
  # Shared checks common to all micro_*() functions
  errors <- errors + validate_micro_inputs(
    DEP = DEP, REFL = REFL, slope = slope, aspect = aspect, hori = hori,
    SLE = SLE, ERR = ERR, RUF = RUF, D0 = D0, Usrhyt = Usrhyt,
    Refhyt = Refhyt, EC = EC, CMH2O = CMH2O,
    TIMAXS = TIMAXS, TIMINS = TIMINS,
    minshade = minshade, maxshade = maxshade,
    Thcond = Thcond, Density = Density, SpecHeat = SpecHeat,
    BulkDensity = BulkDensity,
    run.gads = run.gads
  )

  if(errors == 0){ # continue

    ################## time related variables #################################

    doys12 <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month
    doysn <- doys12 # variable of doys for when doing multiple years
    if(nyears>1 & timeinterval == 365){ # create sequence of days for splining across multiple years
      for(i in 1:(nyears-1)){
        doysn <- c(doysn,(doys12+365*i))
      }
    }

    if(timeinterval < 365){
      microdaily <- 0 # run microclimate model as normal, where each day is iterated ndmax times starting with the initial condition of uniform soil temp at mean monthly temperature
    }else{
      microdaily <- 1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in
    }

    # now check if doing something other than middle day of each month, and create appropriate vector of Day of Year
    daystart <- as.integer(ceiling(365 / timeinterval / 2))
    if(timeinterval!=12){
      doys <- seq(daystart, 365, as.integer(floor(365 / timeinterval)))
    }else{
      doys <- doysn
    }
    doynum <- timeinterval * nyears # total days to do
    doy <- subset(doys, doys != 0) # final vector of Day of Year
    doy <- rep(doy, nyears)
    idayst  <-  1 # start day
    ida <- timeinterval * nyears # end day

    ################## location and terrain #################################

    if(is.numeric(loc) == FALSE){ # might not be quite right format, try correcting
      loc=cbind(as.numeric(loc[1]),as.numeric(loc[2]))
    }
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    # Get reference longitude for solar noon correction (ALREF).
    # timezone = 0: use local apparent solar noon (simple truncation).
    # timezone = 1: look up civil UTC offset via geonames web service.
    ALREF <- get_timezone_alref(lon = x[1], lat = x[2], timezone = timezone)
    HEMIS <- ifelse(x[2]<0,2.,1.) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2])-ALAT)*60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1])-ALONG)*60
    azmuth <- aspect

    if(terrain == 1){
      elevatr <- 1
    }
    if(is.na(elev) & elevatr == 1){
      dem_result <- fetch_dem(loc[1], loc[2], dem.res = dem.res, zmin = zmin,
                              pixels = pixels, terrain = terrain, horizon.step = 15)
      dem <- dem_result$dem
      elev <- dem_result$elev
      if (terrain == 1) { slope <- dem_result$slope; aspect <- dem_result$aspect; hori <- dem_result$hori }
    }

    hori <- as.matrix(hori) #horizon angles
    VIEWF <- 1-sum(sin(as.data.frame(hori) * pi / 180)) / length(hori) # convert horizon angles to radians and calc view factor(s)
    SLES <- rep(SLE, timeinterval * nyears)
    # Expand scalar or short shade inputs to one value per simulation day
    shades <- setup_shade_vectors(minshade, maxshade, ndays = timeinterval * nyears)
    MINSHADES <- shades$MINSHADES
    MAXSHADES <- shades$MAXSHADES
    if(soiltype == 0){ # simulating rock so turn of soil moisture model and set density equal to bulk density
      BulkDensity <- Density
      cap=0
      runmoist <- 0
      PE <- rep(CampNormTbl9_1[1,4],19) #air entry potential J/kg
      KS <- rep(CampNormTbl9_1[1,6],19) #saturated conductivity, kg s/m3
      BB <- rep(CampNormTbl9_1[1,5],19) #soil 'b' parameter
      BD <- rep(BulkDensity,19) # soil bulk density, Mg/m3
      DD <- rep(Density,19) # soil density, Mg/m3
    }else{
      if(soiltype<12){ # use soil properties as specified in Campbell and Norman 1998 Table 9.1
        PE <- rep(CampNormTbl9_1[soiltype,4],19) #air entry potential J/kg
        KS <- rep(CampNormTbl9_1[soiltype,6],19) #saturated conductivity, kg s/m3
        BB <- rep(CampNormTbl9_1[soiltype,5],19) #soil 'b' parameter
        BD <- rep(BulkDensity,19) # soil bulk density, Mg/m3
        DD <- rep(Density,19) # soil density, Mg/m3
      }
    }

    if(soilgrids == 1){
      sg <- fetch_soilgrids(x, DEP)
      if (!is.null(sg)) { PE <- sg$PE; KS <- sg$KS; BB <- sg$BB; BD <- sg$BD; BulkDensity <- sg$BulkDensity }
    }
    # load global climate files
    gcfolder <- paste(.libPaths()[1],"/gcfolder.rda",sep="")
    if(file.exists(gcfolder) == FALSE){
      folder <- "c:/globalclimate"
      if(file.exists(paste0(folder,"/global_climate.nc")) == FALSE){
        message("You don't appear to have the global climate data set - \n run function get.global.climate(folder = 'folder you want to put it in') .....\n exiting function micro_global")
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))
        stop()
      }
    }else{
      load(gcfolder)
    }
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("package 'terra' is needed. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("package 'ncdf4' is needed. Please install it.",
           call. = FALSE)
    }

    message('extracting climate data \n')
    global_climate <- terra::rast(paste0(folder, "/global_climate.nc"))
    CLIMATE <- t(as.numeric(terra::extract(global_climate, x)))
    ALTT <- as.numeric(CLIMATE[1])

    delta_elev <- 0
    if(is.na(elev) == FALSE){ # check if user-specified elevation
      delta_elev <- ALTT - elev # get delta for lapse rate correction
      ALTT <- elev # now make final elevation the user-specified one
    }
    adiab_corr_max <- delta_elev * lapse_max
    adiab_corr_min <- delta_elev * lapse_min
    RAINFALL <- CLIMATE[2:13]
    if(is.na(RAINFALL[1])){
      cat("no climate data for this site, using dummy data so solar is still produced \n")
      CLIMATE <- t(as.numeric(terra::extract(global_climate, cbind(140, -35))))
      CLIMATE[2:97] <- 0
      ALTT<-as.numeric(CLIMATE[1])
      delta_elev <- 0
      if(is.na(elev) == FALSE){ # check if user-specified elevation
        delta_elev <- ALTT - elev # get delta for lapse rate correction
        ALTT <- elev # now make final elevation the user-specified one
      }
      adiab_corr_max <- delta_elev * lapse_max
      adiab_corr_min <- delta_elev * lapse_min
      RAINFALL <- CLIMATE[2:13] * 0
      #stop()
    }
    RAINYDAYS <- CLIMATE[14:25] / 10
    WNMAXX <- CLIMATE[26:37] / 10 * windfac
    WNMINN <- WNMAXX * 0.1 # impose diurnal cycle
    TMINN <- CLIMATE[38:49] / 10
    TMAXX <- CLIMATE[50:61] / 10
    TMAXX <- TMAXX + adiab_corr_max
    TMINN <- TMINN + adiab_corr_min
    ALLMINTEMPS <- TMINN
    ALLMAXTEMPS <- TMAXX
    ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)
    RHMINN <- CLIMATE[62:73] / 10
    RHMAXX <- CLIMATE[74:85] / 10
    # correct for potential change in RH with elevation-corrected Tair
    es <- WETAIR(db = TMAXX, rh = 100)$esat
    e <- WETAIR(db = CLIMATE[50:61] / 10, rh = CLIMATE[62:73] / 10)$e
    RHMINN <- (e / es) * 100
    RHMINN[RHMINN>100] <- 100
    RHMINN[RHMINN<0] <- 0.01
    es <- WETAIR(db = TMINN, rh = 100)$esat
    e <- WETAIR(db = CLIMATE[38:49] / 10, rh = CLIMATE[74:85] / 10)$e
    RHMAXX <- (e / es) * 100
    RHMAXX[RHMAXX>100] <- 100
    RHMAXX[RHMAXX<0] <- 0.01
    CCMINN <- CLIMATE[86:97] / 10
    if(clearsky == 1){
      CCMINN <- CCMINN * 0
    }
    CCMAXX <- CCMINN
    if(runmoist == 0){
      # extract soil moisture
      soilmoisture <- suppressWarnings(terra::rast(paste(folder, "/soilw.mon.ltm.v2.nc", sep = "")))
      message("extracting soil moisture data")
      SoilMoist <- t(as.numeric(terra::extract(soilmoisture, x))) / 1000 # this is originally in mm/m
    }
    if(is.na(max(SoilMoist, ALTT, CLIMATE)) == TRUE){
      message("Sorry, there is no environmental data for this location")
      SoilMoist <- t(as.numeric(terra::extract(soilmoisture, cbind(140, -35)))) / 1000 * (1 - BulkDensity / Density) # this is originally in mm/m
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
    WNMINN <- WNMINN * (1.2 / 10) ^ 0.15
    WNMAXX <- WNMAXX * (1.2 / 10) ^ 0.15
    # impose uniform warming
    TMAXX <- TMAXX + warm
    TMINN <- TMINN + warm
    if(timeinterval != 12){ # spline from 12 days to chosen time interval
      TMAXX1 <- suppressWarnings(spline(doys12,TMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      TMAXX <- rep(TMAXX1$y,nyears)
      TMINN1 <- suppressWarnings(spline(doys12,TMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      TMINN <- rep(TMINN1$y,nyears)
      RHMAXX1 <-suppressWarnings(spline(doys12,RHMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      RHMAXX <- rep(RHMAXX1$y,nyears)
      RHMINN1 <-suppressWarnings(spline(doys12,RHMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      RHMINN <- rep(RHMINN1$y,nyears)
      CCMAXX1 <-suppressWarnings(spline(doys12,CCMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      CCMAXX <- rep(CCMAXX1$y,nyears)
      CCMINN <- CCMAXX
      WNMAXX1 <- suppressWarnings(spline(doys12,WNMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      WNMAXX <- rep(WNMAXX1$y,nyears)
      WNMINN1 <- suppressWarnings(spline(doys12,WNMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      WNMINN <- rep(WNMINN1$y,nyears)
      if(runmoist == 0){
        SoilMoist1 <- suppressWarnings(spline(doys12,SoilMoist,n=timeinterval,xmin=1,xmax=365,method="periodic"))
        SoilMoist <- rep(SoilMoist1$y,nyears)
      }
    }
    if(timeinterval<365){
      TMAXX <- rep(TMAXX,nyears)
      TMINN <- rep(TMINN,nyears)
      RHMAXX <- rep(RHMAXX,nyears)
      RHMINN <- rep(RHMINN,nyears)
      CCMAXX <- rep(CCMAXX,nyears)
      CCMINN <- rep(CCMINN,nyears)
      WNMAXX <- rep(WNMAXX,nyears)
      WNMINN <- rep(WNMINN,nyears)
      if(runmoist == 0){
        SoilMoist <- rep(SoilMoist,nyears)
      }
      RAINFALL <- rep(RAINFALL,nyears)
    }
    orig.RAINFALL <- RAINFALL
    # get annual mean temp for creating deep soil (2m) boundary condition
    avetemp <- (sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
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
    deepsoil <- rep(tannul, doynum)  # deep soil boundary temperature (annual mean)

    daymon <- c(31,28,31,30,31,30,31,31,30,31,30,31) # days in each month

    # if doing daily sims, spread rainfall evenly across days based on mean monthly rainfall and the number of rainy days per month
    if(timeinterval == 365){
      RAINFALL1 <- 1:365
      sort <- matrix(data = 0,nrow = 365,ncol = 2)
      m <- 1
      b <- 0
      for (i in 1:12){ #begin loop through 12 months of year
        ndays=daymon[i]
        for (k in 1:ndays){
          b <- b+1
          sort[m,1] <- i
          sort[m,2] <- b
          if(k<=RAINYDAYS[i] & rainfrac>0){
            if(k == 1){
              RAINFALL1[m] <- RAINFALL[i]*rainfrac*rainmult # if first day of month, make user-specified fraction of monthly rainfall fall on first day
            }else{
              RAINFALL1[m] <- (RAINFALL[i]*(1-rainfrac)*rainmult)/RAINYDAYS[i] # make remaining rain fall evenly over the remaining number of rainy days for the month, starting at the beginning of the month
            }
          }else{
            if(rainfrac == 0){
              RAINFALL1[m] <- (RAINFALL[i]*rainmult)/RAINYDAYS[i]
            }else{
              RAINFALL1[m] <- 0
            }
          }
          m <- m+1
          if(b>RAINYDAYS[i]){
            b <- 0
          }
        }
      }
      RAINFALL2 <- as.data.frame(cbind(RAINFALL1,sort))
      #RAINFALL2 <- RAINFALL2[order(RAINFALL2$V2,RAINFALL2$V3),] # this line scatters the rainy days evenly across each month - snow predictions better if it is commented out so get rainy days all in a row within the month
      RAINFALL <- rep(as.double(RAINFALL2$RAINFALL1),nyears)
      RAINFALL[!is.finite(RAINFALL)] <- 0
      if(TMINN[1]<snowtemp){
        RAINFALL[1] <- 0 # this is needed in some cases to allow the integrator to get started
      }
    }else{
      if(timeinterval!=12){
        RAINFALL <- rep(rep(sum(RAINFALL)/timeinterval,timeinterval),nyears) # just spread evenly across every day
      }else{ # running middle day of each month - divide monthly rain by number of days in month
        RAINFALL <- RAINFALL/rep(daymon,nyears)
      }
    }#end check doing daily sims

    ndays <- length(RAINFALL)
    SOLRhr <- rep(0,24*ndays)

    hourly <- 0
    if(microclima > 0 & timeinterval %in% c(12, 365)){

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
      yearlist <- seq(1960, 1960 + nyears - 1)
      tt <- seq(as.POSIXct(paste0('01/01/', yearlist[1]), format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(paste0('31/12/', yearlist[nyears]), format = "%d/%m/%Y", tz = 'UTC')+23*3600, by = 'hours')
      timediff <- x[1] / 15
      hour.microclima <- as.numeric(format(tt, "%H")) + timediff-floor(timediff)
      jd <- julday(as.numeric(format(tt, "%Y")), as.numeric(format(tt, "%m")), as.numeric(format(tt, "%d")))
      dem_result <- fetch_dem(loc[1], loc[2],
                              dem.res = microclima.dem.res, zmin = microclima.zmin,
                              pixels = NULL,
                              existing_dem = if (inherits(dem, c("SpatRaster", "RasterLayer"))) dem else NULL,
                              terrain = 1, horizon.step = 10, nangles = 36,
                              slope = slope, aspect = aspect)
      dem <- dem_result$dem; slope <- dem_result$slope; aspect <- dem_result$aspect
      ha36 <- dem_result$hori
      for (i in 1:length(hour.microclima)) {
        saz <- solazi(hour.microclima[i], lat, long, jd[i], merid = long)
        saz <- round(saz/10, 0) + 1
        saz <- ifelse(saz > 36, 1, saz)
        ha[i] <- ha36[saz]
      }
      cloud <- rep(CCMAXX / 2, nyears)
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
          cloud1 <- suppressWarnings(spline(doys12,cloud[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          cloud2 <- cloud1$y
        }else{
          start <- end + 1
          end <- end + 12
          cloud1 <- suppressWarnings(spline(c(0, doys12), c(tail(cloud2, 1), cloud[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          cloud2 <- c(cloud2, cloud1$y)
        }
      }
      cloudhr <- cbind(rep(seq(1, length(cloud2)),24), rep(cloud2, 24))
      cloudhr <- cloudhr[order(cloudhr[,1]),]
      cloudhr <- cloudhr[,2]
      #cloudhr <- leapfix(cloudhr, yearlist, 24)
      micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, TAI = TAI, timeinterval = 365, solonly = 1)
      clearskyrad <- micro_clearsky$metout[, c(1, 13)][, 2]
      dsw2 <- leapfix(clearskyrad, yearlist, 24) *(0.36+0.64*(1-cloudhr/100)) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
      # partition total solar into diffuse and direct using code from microclima::hourlyNCEP
      sol <- compute_solar_partition(dsw2, jd, hour.microclima, lat, long,
                                     slope, aspect, ha, microclima.LOR, microclima.LAI)
      SOLRhr_all <- sol$SOLRhr_all
      diffuse_frac_all <- sol$diffuse_frac_all
      diffuse_frac <- diffuse_frac_all
      if(microclima == 2){ # use hourly solar from microclima
        hourly <- 2
        VIEWF <- 1 # accounted for already in microclima cals
        hori <- rep(0, 24) # accounted for already in microclima calcs
      }
    }else{
      diffuse_frac <- NA
    }
    if(timeinterval == 12 & microclima > 0){
      dates_all <- head(seq(as.POSIXct(paste0("01/01/", yearlist[1]), format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(paste0("01/01/", yearlist[nyears] + 1), format = "%d/%m/%Y ", tz = 'UTC'), by = 'hours'), -1)
      dates_15th <- which(format(dates_all, "%d") == "15")
      diffuse_frac <- diffuse_frac_all[dates_15th]
    }

    if(length(TAI) < 111){ # no user supplied values, compute with GADS
      TAI <- compute_tai(longlat, run.gads, TAI_ELTERMAN)
    }
    ################ soil properties  ##################################################
    Nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
    if(soilgrids == 1){
      Numtyps <- 10 # number of substrate types
      Nodes[1:10,] <- c(1:10) # deepest nodes for each substrate type
    }else{
      Numtyps <- 2 # number of soil types
      Nodes[1,1:ndays] <- 3 # deepest node for first substrate type
      Nodes[2,1:ndays] <- 9 # deepest node for second substrate type
    }
    REFLS <- rep(REFL,ndays) # soil reflectances
    PCTWET <- rep(PCTWET,ndays) # soil wetness
    if(runmoist == 0){
      moists2 <- matrix(nrow= 10, ncol = ndays, data=0) # set up an empty vector for soil moisture values through time
      moists2[1,] <- SoilMoist # fill the first row with monthly soil moisture values
      moists2[2,] <- moists2[1,] # make this row same as first row
      moists <- moists2
    }else{
      moists2 <- matrix(nrow=10, ncol = ndays, data=0) # set up an empty vector for soil moisture values through time
      moists2[1:10,] <- SoilMoist_Init
      moists2[moists2>(1-BulkDensity/Density)] <- (1-BulkDensity/Density)
      moists <- moists2
    }

    # now make the soil properties matrix
    # columns are:
    #1) bulk density (Mg/m3)
    #2) volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    #3) thermal conductivity (W/mK)
    #4) specific heat capacity (J/kg-K)
    #5) mineral density (Mg/m3)
    soilprops <- matrix(data = 0, nrow = 10, ncol = 5) # create an empty soil properties matrix
    if(soilgrids == 1){
      soilprops[,1] <- BulkDensity
      soilprops[,2] <- 1 - BulkDensity / Density # not used if soil moisture computed
      soilprops[soilprops[,2] < 0.26, 2] <- 0.26
      soilprops[,3] <- Thcond
      soilprops[,4] <- SpecHeat
      soilprops[,5] <- Density
      if(cap == 1){
        soilprops[1:2,3] <- 0.2
        soilprops[1:2,4] <- 1920
      }
      if(cap == 2){
        soilprops[1:2,3] <- 0.1
        soilprops[3:4,3] <- 0.25
        soilprops[1:4,4] <- 1920
        soilprops[1:4,5] <- 1.3
        soilprops[1:4,1] <- 0.7
      }
    }else{
      soilprops[1,1] <- BulkDensity # insert soil bulk density to profile 1
      soilprops[2,1] <- BulkDensity # insert soil bulk density to profile 2
      soilprops[1,2] <- min(0.26, 1 - BulkDensity / Density) # insert saturated water content to profile 1
      soilprops[2,2] <- min(0.26, 1 - BulkDensity / Density) # insert saturated water content to profile 2
      if(cap == 1){ # insert thermal conductivity to profile 1, and see if 'organic cap' added on top
        soilprops[1,3] <- 0.2 # mineral thermal conductivity
      }else{
        soilprops[1,3] <- Thcond # mineral thermal conductivity
      }
      soilprops[2,3] <- Thcond # insert thermal conductivity to profile 2
      if(cap == 1){ # insert specific heat to profile 1, and see if 'organic cap' added on top
        soilprops[1,4] <- 1920 # mineral heat capacity
      }else{
        soilprops[1,4] <- SpecHeat
      }
      soilprops[2,4] <- SpecHeat # insert specific heat to profile 2
      soilprops[1,5] <- Density # insert mineral density to profile 1
      soilprops[2,5] <- Density # insert mineral density to profile 2
    }
    #########################################################################################

    # --- 3. Prepare model inputs ----------------------------------------------

    # Hourly input override vectors — unused in micro_global (daily dataset),
    # so fill with zeros / sentinel values that tell the Fortran to ignore them
    hourly     <- 0
    rainhourly <- 0
    TAIRhr <- rep(0,  24 * ndays)
    RHhr   <- rep(0,  24 * ndays)
    WNhr   <- rep(0,  24 * ndays)
    CLDhr  <- rep(0,  24 * ndays)
    SOLRhr <- rep(0,  24 * ndays)
    RAINhr <- rep(0,  24 * ndays)
    ZENhr  <- rep(-1, 24 * ndays)   # -1 = not supplied
    IRDhr  <- rep(-1, 24 * ndays)   # -1 = not supplied

    # Build the flat numeric vector consumed by the Fortran microclimate solver
    microinput <- build_microinput(
      ndays = ndays, RUF = RUF, ERR = ERR, Usrhyt = Usrhyt, Refhyt = Refhyt,
      Numtyps = Numtyps, Z01 = Z01, Z02 = Z02, ZH1 = ZH1, ZH2 = ZH2,
      idayst = idayst, ida = ida, HEMIS = HEMIS, ALAT = ALAT, AMINUT = AMINUT,
      ALONG = ALONG, ALMINT = ALMINT, ALREF = ALREF, slope = slope,
      azmuth = azmuth, ALTT = ALTT, CMH2O = CMH2O, microdaily = microdaily,
      tannul = tannul, EC = EC, VIEWF = VIEWF, snowtemp = snowtemp,
      snowdens = snowdens, snowmelt = snowmelt, undercatch = undercatch,
      rainmult = rainmult, runshade = runshade, runmoist = runmoist,
      maxpool = maxpool, evenrain = evenrain, snowmodel = snowmodel,
      rainmelt = rainmelt, writecsv = writecsv, densfun = densfun,
      hourly = hourly, rainhourly = rainhourly, lamb = lamb, IUV = IUV,
      RW = RW, PC = PC, RL = RL, SP = SP, R1 = R1, IM = IM,
      MAXCOUNT = MAXCOUNT, IR = IR, message = message, fail = fail,
      snowcond = snowcond, intercept = intercept, grasshade = grasshade,
      solonly = solonly, ZH = ZH, D0 = D0, TIMAXS = TIMAXS, TIMINS = TIMINS,
      spinup = spinup, maxsurf = maxsurf, ndmax = ndmax
      # dewrain / moiststep left at defaults (0 / 360) for this daily dataset
    )

    if (length(LAI) < ndays) {
      LAI <- rep(LAI[1], ndays)
    }
    if (shore == 0) {
      tides <- matrix(data = 0, nrow = 24 * ndays, ncol = 3)
    }

    # Assemble the complete input list for microclimate()
    micro <- build_micro_list(
      microinput = microinput, doy = doy, SLES = SLES, DEP = DEP,
      Nodes = Nodes, MAXSHADES = MAXSHADES, MINSHADES = MINSHADES,
      TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN,
      CCMAXX = CCMAXX, CCMINN = CCMINN, WNMAXX = WNMAXX, WNMINN = WNMINN,
      TAIRhr = TAIRhr, RHhr = RHhr, WNhr = WNhr, CLDhr = CLDhr,
      SOLRhr = SOLRhr, RAINhr = RAINhr, ZENhr = ZENhr, IRDhr = IRDhr,
      REFLS = REFLS, PCTWET = PCTWET, soilinit = soilinit, hori = hori,
      TAI = TAI, soilprops = soilprops, moists = moists, RAINFALL = RAINFALL,
      deepsoil = deepsoil, PE = PE, KS = KS, BB = BB, BD = BD, DD = DD,
      L = L, LAI = LAI, tides = tides
    )

    # Optionally dump all inputs to CSV files for debugging
    if (write_input == 1) {
      write_micro_csv(
        microinput = microinput, doy = doy, SLES = SLES, DEP = DEP,
        Nodes = Nodes, MAXSHADES = MAXSHADES, MINSHADES = MINSHADES,
        TIMAXS = TIMAXS, TIMINS = TIMINS,
        TMAXX = TMAXX, TMINN = TMINN, RHMAXX = RHMAXX, RHMINN = RHMINN,
        CCMAXX = CCMAXX, CCMINN = CCMINN, WNMAXX = WNMAXX, WNMINN = WNMINN,
        REFLS = REFLS, PCTWET = PCTWET, soilinit = soilinit, hori = hori,
        TAI = TAI, soilprops = soilprops, moists = moists, RAINFALL = RAINFALL,
        deepsoil = deepsoil, PE = PE, BD = BD, DD = DD, BB = BB, KS = KS,
        L = L, LAI = LAI, tides = tides,
        TAIRhr = TAIRhr, RHhr = RHhr, WNhr = WNhr, CLDhr = CLDhr,
        SOLRhr = SOLRhr, RAINhr = RAINhr, ZENhr = ZENhr, IRDhr = IRDhr
      )
    }

    # --- 4. Run the microclimate model ----------------------------------------
    if (is.numeric(loc[1])) {
      location <- paste("long", loc[1], "lat", loc[2])
    } else {
      location <- loc
    }

    if (runmicro) {
      message(paste('running microclimate model for', timeinterval, 'days by',
                    nyears, 'years at site', location))
      message('Note: the output column `SOLR` in metout and shadmet is for unshaded horizontal plane solar radiation')
      ptm <- proc.time()
      microut <- microclimate(micro)
      message(paste0('runtime ', (proc.time() - ptm)[3], ' seconds'))

      if (timeinterval == 12) {
        RAINFALL <- orig.RAINFALL  # restore monthly rainfall for the return value
      }
      # Compute hourly and daily date sequences for the output
      day_seq <- rep(seq(1, timeinterval * nyears), 24)
      day_seq <- day_seq[order(day_seq)]
      metout  <- microut$metout

      if (max(metout[, 1] == 0)) {
        message("ERROR: the model crashed - try a different error tolerance (ERR) or a different spacing in DEP")
      }
      dates   <- day_seq + metout[, 2] / 60 / 24 - 1  # fractional day index
      dates2  <- seq(1, timeinterval * nyears)          # integer day index

      # --- 5. Process and return results ------------------------------------
      out <- process_micro_output(microut,
                                  runmoist  = runmoist,
                                  snowmodel = snowmodel,
                                  lamb      = lamb)

      return(build_micro_return(
        out          = out,
        RAINFALL     = RAINFALL,
        ndays        = ndays,
        ALTT         = ALTT,
        REFL         = REFL,
        longlat      = c(x[1], x[2]),
        nyears       = nyears,
        timeinterval = timeinterval,
        MINSHADES    = MINSHADES,
        MAXSHADES    = MAXSHADES,
        DEP          = DEP,
        dates        = dates,
        dates2       = dates2,
        PE           = PE,
        BD           = BD,
        DD           = DD,
        BB           = BB,
        KS           = KS,
        dem          = dem,
        diffuse_frac = diffuse_frac,
        snowmodel    = snowmodel,
        lamb         = lamb
      ))

    } else {
      # runmicro = FALSE: return only the prepared weather inputs, not model output
      days_hr <- rep(seq(1, timeinterval * nyears), 24)
      days_hr <- days_hr[order(days_hr)]
      dates   <- days_hr
      dates2  <- seq(1, timeinterval * nyears)
      return(list(
        RAINFALL = RAINFALL, TMAXX = TMAXX, TMINN = TMINN,
        RHMAXX = RHMAXX, RHMINN = RHMINN, WNMAXX = WNMAXX, WNMINN = WNMINN,
        CCMAXX = CCMAXX, CCMINN = CCMINN, CLDhr = CLDhr, WNhr = WNhr,
        TAIRhr = TAIRhr, RHhr = RHhr, RAINhr = RAINhr, SOLRhr = SOLRhr,
        ZENhr = ZENhr, IRDhr = IRDhr, dates = dates, dates2 = dates2,
        PE = PE, BD = BD, DD = DD, BB = BB, KS = KS
      ))
    }
  } # end if(errors == 0)
} # end of micro_global function
