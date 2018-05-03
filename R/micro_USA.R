#' USA implementation of the microclimate model.
#'
#' An implementation of the Niche Mapper microclimate model that uses the GRIDMET daily weather database http://www.climatologylab.org/gridmet.html, and specifically uses the following variables: pr, rmax, rmin, srad, tmmn, tmmx, vs. Also uses the following DEM "metdata_elevationdata.nc".
#' @param loc Either a longitude and latitude (decimal degrees) or a place name to search for on Google Earth
#' @param timeinterval The number of time intervals to generate predictions for over a year (must be 12 <= x <=365)
#' @param ystart First year to run
#' @param yfinish Last year to run
#' @param REFL Soil solar reflectance, decimal \%
#' @param elev Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param soiltype Soil type: Rock = 0, sand = 1, loamy sand = 2, sandy loam = 3, loam = 4, silt loam = 5, sandy clay loam = 6, clay loam = 7, silt clay loam = 8, sandy clay = 9, silty clay = 10, clay = 11, user-defined = 12, based on Campbell and Norman 1990 Table 9.1.
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
#' @usage micro_USA(loc = "Madison Wisconsin, USA", timeinterval = 365, ystart = 1990, yfinish = 1990, soiltype = 4,
#' REFL = 0.15, slope = 0, aspect = 0, DEP = c(0., 2.5,  5.,  10.,  15,  20,  30,  50,  100,  200), minshade = 0, maxshade = 90,
#' Usrhyt = 0.01, ...)
#' @export
#'
#' @details
#' \itemize{
#' \strong{ Parameters controling how the model runs:}\cr\cr
#'
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{rungads}{ = 1, Use the Global Aerosol Database? 1=yes, 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{terrain}{ = 0, Use 250m resolution terrain data? 1=yes, 0=no}\cr\cr
#' \code{dailywind}{ = 1, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{windfac}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{adiab_cor}{ = 1, use adiabatic lapse rate correction? 1=yes, 0=no}\cr\cr
#' \code{warm}{ = 0, uniform warming, deg C}\cr\cr
#' \code{spatial}{ = "c:/Australian Environment/", choose location of terrain data}\cr\cr
#' \code{opendap}{ = 1, query met grids via opendap (does not work on PC unless you compile ncdf4 - see https://github.com/pmjherman/r-ncdf4-build-opendap-windows)}\cr\cr
#' \code{soilgrids}{ = 1, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#' \code{save}{ = 0, don't save forcing data (0), save the forcing data (1) or read previously saved data (2)}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{ERR}{ = 2.5, Integrator error tolerance for soil temperature calculations}\cr\cr
#' \code{Refhyt}{ = 2, Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured}\cr\cr
#' \code{RUF}{ = 0.004, Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001 (snow) - 0.02 m.}\cr\cr
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
#' \code{DEP}
#' { and points half way between)}\cr\cr
#' \code{LAI}{ = 0.1, leaf area index, used to partition traspiration/evaporation from PET}\cr\cr
#'
#' \strong{ Snow mode parameters:}
#'
#' \code{snowmodel}{ = 1, run the snow model 1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (deg C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (mg/m3), overridden by densfun}\cr\cr
#' \code{densfun}{ = c(0,0), slope and intercept of linear model of snow density as a function of day-of-year - if it is c(0,0) then fixed density used}\cr\cr
#' \code{snowmelt}{ = 0.9, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
#' \code{undercatch}{ = 1, undercatch multipier for converting rainfall to snow}\cr\cr
#' \code{rainmelt}{ = 0.0125, paramter in equation that melts snow with rainfall as a function of air temp}\cr\cr
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)}\cr\cr
#' \code{intercept}{ = maxshade / 100 * 0.3, snow interception fraction for when there's shade (0-1)}\cr\cr
#'
#' \strong{ Intertidal mode parameters:}
#'
#' \code{shore}{ Include tide effects? If 1, the matrix}
#' \code{tides}
#' { is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0., nrow = 24*timeinterval*nyears, ncol = 3), matrix for each how of the simulation of 1. tide state (0=out, 1=in), 2. Water temperature (deg C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#' }
#'
#' \strong{Outputs:}
#' metout/shadmet variables:
#' \itemize{
#' \item 1 DOY - day-of-year
#' \item 2 TIME - time of day (mins)
#' \item 3 TALOC - air temperature (deg C) at local height (specified by 'Usrhyt' variable)
#' \item 4 TAREF - air temperature (deg C) at reference height (specified by 'Refhyt', 2m default)
#' \item 5 RHLOC - relative humidity (\%) at local height (specified by 'Usrhyt' variable)
#' \item 6 RH  - relative humidity (\%) at reference height (specified by 'Refhyt', 2m default)
#' \item 7 VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)
#' \item 8 VREF - wind speed (m/s) at reference height (specified by 'Refhyt', 2m default)
#' \item 9 SNOWMELT - snowmelt (mm)
#' \item 10 POOLDEP - water pooling on surface (mm)
#' \item 11 PCTWET - soil surface wetness (\%)
#' \item 12 ZEN - zenith angle of sun (degrees - 90 = below the horizon)
#' \item 13 SOLR - solar radiation (W/m2)
#' \item 14 TSKYC - sky radiant temperature (deg C)
#' \item 15 DEW - dew presence (0 or 1)
#' \item 16 FROST - frost presence (0 or 1)
#' \item 17 SNOWFALL - snow predicted to have fallen (cm)
#' \item 18 SNOWDEP - predicted snow depth (cm)
#'}
#' soil and shadsoil variables:
#' \itemize{
#' \item 1 DOY - day-of-year
#' \item 2 TIME - time of day (mins)
#' \item 3-12 D0cm ... - soil temperature (deg C) at each of the 10 specified depths
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
#' \item  4 LEAFPOT - leaf water potentail (J/kg)
#' \item  5-14 RPOT0cm ... - root water potentail (J/kg), at each of the 10 specified depths
#' }
#'
#' if snow model is run i.e. parameter lamb = 1\cr
#' sunsnow and shdsnow variables:
#' \itemize{
#' \item  1 DOY - day-of-year
#' \item  2 TIME - time of day (mins)
#' \item  3-10 SN1 ... - snow temperature (deg C), at each of the potential 8 snow layers (layer 8 is always the bottom - need metout$SNOWDEP to interpret which depth in the snow a given layer represents)
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
#' micro<-micro_usa(loc = 'Madison Wisconsin, USA', runshade = 0, ystart = 2014, yfinish = 2016) # run the model using SoilGrids data at Madison for 2014 to 2016
#'
#' metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#' soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
#' soilmoist<-as.data.frame(micro$soilmoist) # soil temperatures, minimum shade
#'
#' # append dates
#' ystart <- 2014
#' yfinish <- 2016
#' nyears <- yfinish-ystart+1
#' tzone<-paste("Etc/GMT+",0,sep="")
#' dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
#'
#' metout <- cbind(dates,metout)
#' soil <- cbind(dates,soil)
#' soilmoist <- cbind(dates, soilmoist)
#' minshade<-micro$minshade
#'
#' # plotting above-ground conditions in minimum shade
#' with(metout,{plot(TALOC ~ dates,xlab = "Date and Time", ylab = "Temperature (deg C)"
#' , type = "l",main=paste("air and sky temperature, ",minshade,"% shade",sep=""), ylim = c(-20, 60))})
#' with(metout,{points(TAREF ~ dates,xlab = "Date and Time", ylab = "Temperature (deg C)"
#' , type = "l",lty=2,col='blue')})
#' with(metout,{points(TSKYC ~ dates,xlab = "Date and Time", ylab = "Temperature (deg C)"
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
#'    plot(soil[,i+3]~soil[,1],xlab = "Date and Time", ylab = "Soil Temperature (deg C)"
#'    ,col=i,type = "l",main=paste("soil temperature ",minshade,"% shade",sep=""))
#'  }else{
#'    points(soil[,i+3]~soil[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (deg C)",col=i,type = "l")
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
micro_usa <- function(loc = "Madison, Wisconsin", timeinterval = 365, ystart = 2016,
  yfinish = 2016, nyears = 1, soiltype = 4, REFL = 0.15, elev = NA, slope = 0,
  aspect = 0, lapse_max = 0.0077, lapse_min = 0.0039,
  DEP=c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200), minshade = 0, maxshade = 90,
  Refhyt = 2, Usrhyt = 0.01, Z01 = 0, Z02 = 0, ZH1 = 0, ZH2 = 0, runshade = 1,
  clearsky = 0, rungads = 1, write_input = 0, writecsv = 0,
  terrain = 0, dailywind = 1, windfac = 1, adiab_cor = 1, warm = 0,
  spatial = "N:/USA", ERR = 2.5, RUF = 0.004, EC = 0.0167238, SLE = 0.95,
  Thcond = 2.5, Density = 2.56, SpecHeat = 870, BulkDensity = 1.3, PCTWET = 0,
  rainwet = 1.5, cap = 1, CMH2O = 1, hori = rep(0,24), TIMAXS=c(1.0, 1.0, 0.0, 0.0),
  TIMINS = c(0, 0, 1, 1), timezone = 0, runmoist = 1, PE = rep(1.1, 19),
  KS = rep(0.0037, 19), BB = rep(4.5, 19), BD = rep(BulkDensity, 19),
  DD = rep(Density, 19), maxpool = 10000, rainmult = 1, evenrain = 0,
  SoilMoist_Init = c(0.1, 0.12, 0.15, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
  L = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000,
  R1 = 0.001, RW = 2.5e+10, RL = 2e+6, PC = -1500, SP = 10, IM = 1e-06, MAXCOUNT = 500,
  LAI = 0.1, snowmodel = 1, snowtemp = 1.5, snowdens = 0.375, densfun = c(0, 0),
  snowmelt = 0.9, undercatch = 1, rainmelt = 0.0125, shore = 0,
  tides = matrix(data = 0, nrow = 24 * timeinterval * nyears, ncol = 3),
  scenario = "", year = "", hourly = 0, rainhourly = 0, rainhour = 0,
  rainoff = 0, lamb = 0, IUV = 0, opendap = 1, soilgrids = 1, IR = 0, message = 0,
  fail = nyears * 24 * 365, save = 0, snowcond = 0, intercept = maxshade / 100 * 0.3) {

  # loc="Madison, Wisconsin"
  # timeinterval=365
  # ystart=2014
  # yfinish=2016
  # nyears=1
  # soiltype=4
  # REFL=0.15
  # slope=0
  # aspect=0
  # DEP=c(0., 2.5,  5.,  10.,  15.,  20.,  30.,  50.,  100.,  200.)
  # minshade=0
  # maxshade=90
  # Refhyt=2
  # Usrhyt=.01
  # Z01=0
  # Z02=0
  # ZH1=0
  # ZH2=0
  # runshade=0
  # clearsky=0
  # rungads=1
  # write_input=0
  # writecsv=0
  # terrain=0
  # dailywind=1
  # adiab_cor=1
  # warm=0
  # spatial="C:/USA"
  # ERR=1.5
  # RUF=0.004
  # EC=0.0167238
  # SLE=0.95
  # Thcond=2.5
  # Density=2.56
  # SpecHeat=870
  # BulkDensity=1.3
  # PCTWET=0
  # rainwet=1.5
  # cap=1
  # CMH2O=1
  # hori=rep(0,24)
  # TIMAXS=c(1.0, 1.0, 0.0, 0.0)
  # TIMINS=c(0, 0, 1, 1)
  # timezone=0
  # runmoist=1
  # PE=rep(1.1,19)
  # KS=rep(0.0037,19)
  # BB=rep(4.5,19)
  # BD=rep(1.3,19)
  # DD=rep(2.56,19)
  # maxpool=10000
  # rainmult=1
  # evenrain=0
  # SoilMoist_Init=c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3)
  # L = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000
  # LAI=0.1
  # snowmodel=1
  # snowtemp=1.5
  # snowdens=0.375
  # densfun=c(0,0)
  # snowmelt=0.9
  # undercatch=1
  # rainmelt=0.0125
  # shore=0
  # tides=matrix(data = 0., nrow = 24*timeinterval*nyears, ncol = 3)
  # scenario=""
  # hourly=0
  # rainhour = 0
  # rainoff=0
  # lamb = 0
  # IUV = 0
  # R1 = 0.001
  # RW = 2.5e+10
  # RL = 2e+6
  # PC = -1500
  # SP = 10
  # IM = 1e-06
  # MAXCOUNT = 500
  # windfac=1
  # rainhourly = 0
  # opendap = 1
  # soilgrids = 1
  # IR = 0
  # message = 0
  # fail = nyears * 24 * 365
  # elev = NA
  # lapse_max = 0.0077
  # lapse_min = 0.0039
  # snowcond = 0
  # intercept = maxshade / 100 * 0.3

  # error trapping - originally inside the Fortran code, but now checking before executing Fortran
  errors<-0
  if(DEP[2]-DEP[1]>3 | DEP[3]-DEP[2]>3){
    cat("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if(DEP[2]-DEP[1]<2){
    cat("warning, nodes might be too close near the surface, try a different spacing if the program is crashing \n")
  }
  if(timeinterval<12 | timeinterval > 365){
    cat("ERROR: the variable 'timeinterval' is out of bounds.
        Please enter a correct value (12 - 365).", '\n')
    errors<-1
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
  if(rungads%in%c(0,1)==FALSE){
    cat("ERROR: the variable 'rungads' be either 0 or 1.
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
  if(soiltype<0 | soiltype>11){
    cat("ERROR: the soil type must range between 1 and 11.
        Please correct.", '\n')
    errors<-1
  }
  # end error trapping

  if(errors==0){ # continue

    ################## time related variables #################################
    nyears<-yfinish-ystart+1

    doys12<-c(15.,46.,74.,105.,135.,166.,196.,227.,258.,288.,319.,349.) # middle day of each month

    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in

    daystart<-as.integer(ceiling(365/timeinterval/2))
    idayst <- 1 # start day
    dates<-Sys.time()-60*60*24
    curyear<-as.numeric(format(dates,"%Y"))

    ################## location and terrain #################################
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("package 'raster' is needed. Please install it.",
        call. = FALSE)
    }
    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("package 'ncdf4' is needed. Please install it.",
        call. = FALSE)
    }
    if(is.numeric(loc)==FALSE){ # use geocode to get location from site name via googlemaps
      if (!requireNamespace("dismo", quietly = TRUE)) {
        stop("dismo needed for the place name geocode function to work. Please install it.",
          call. = FALSE)
      }
      if (!requireNamespace("XML", quietly = TRUE)) {
        stop("XML needed for the place name geocode function to work. Please install it.",
          call. = FALSE)
      }
      if (!requireNamespace("proj4", quietly = TRUE)) {
        stop("package 'proj4' is needed. Please install it.",
          call. = FALSE)
      }
      longlat <- dismo::geocode(loc)[3:4] # assumes first geocode match is correct
      if(nrow(longlat>1)){longlat<-longlat[1,]}
      x <- t(as.matrix(as.numeric(c(longlat[1,1],longlat[1,2]))))
    }else{
      longlat <- loc
      x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))
    }
    library(raster)
    library(proj4)
    library(ncdf4)
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
        require(utils)
        require(futile.logger)
        retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=20, sleep=0) {
          attempts = 0
          retval = try(eval(expr))
          while (isError(retval)) {
            attempts = attempts + 1
            if (attempts >= maxErrors) {
              msg = sprintf("retry: too many retries [[%s]]", capture.output(str(retval)))
              flog.fatal(msg)
              stop(msg)
            } else {
              msg = sprintf("retry: error in attempt %i/%i [[%s]]", attempts, maxErrors,
                capture.output(str(retval)))
              flog.error(msg)
              warning(msg)
            }
            if (sleep > 0) Sys.sleep(sleep)
            retval = try(eval(expr))
          }
          return(retval)
        }

        cat("extracting elevation via opendaps \n")
        baseurl <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/"
        nc <- nc_open(paste0(baseurl, "/elev/metdata_elevationdata.nc"))
        lon <- ncvar_get(nc, "lon")
        lat <- ncvar_get(nc, "lat")
        flat=match(abs(lat-x[2])<1/48,1)
        latindex=which(flat %in% 1)
        flon=match(abs(lon-x[1])<1/48,1)
        lonindex=which(flon %in% 1)
        start <- c(lonindex,latindex,1)
        count <- c(1, 1, 1)
        USADEM <- retry(as.numeric(ncvar_get(nc, varid = "elevation",
          start = start, count)))
        nc_close(nc)
      }else{
        USADEM <- extract(raster(paste0(spatial,"/metdata_elevationdata.nc")), x) # metres above sea level
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
        cat('extracting data from SoilGrids \n')
        require(rjson)
        require(sp)
        require(GSIF)
        pnts <- data.frame(lon=x[1], lat=x[2], id=c("p1"))
        coordinates(pnts) <- ~lon+lat
        proj4string(pnts) <- CRS("+proj=longlat +datum=WGS84")
        soilgrids.r <- REST.SoilGrids(c("BLDFIE", "SLTPPT","SNDPPT", "CLYPPT"))
        ov <- over(soilgrids.r, pnts)
        if(length(ov) > 3){
          soilpro <- cbind(c(0,5,15,30,60,100,200), t(ov[3:9])/1000, t(ov[11:17]), t(ov[19:25]), t(ov[27:33]) )
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
      cat("loading SoilGrids data from previous run \n")
      load('PE.Rda')
      load('BB.Rda')
      load('BD.Rda')
      load('KS.Rda')
      load('BulkDensity.Rda')
    }
    if(save == 1){
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
        baseurl <- "http://thredds.northwestknowledge.net:8080/thredds/dodsC/MET/"
        yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
        for (j in 1:nyears) {
          if (j == 1) {
            cat(paste("reading weather input for ", yearlist[j], " \n", sep = ""))
            cat(paste("tmin weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "/tmmn/tmmn_", yearlist[j],
              ".nc"))
            lon <- ncvar_get(nc, "lon")
            lat <- ncvar_get(nc, "lat")
            flat=match(abs(lat-x[2])<1/48,1)
            latindex=which(flat %in% 1)
            flon=match(abs(lon-x[1])<1/48,1)
            lonindex=which(flon %in% 1)
            start <- c(latindex,lonindex,1)
            count <- c(1, 1, -1)
            tmin <- retry(as.numeric(ncvar_get(nc, varid = "air_temperature",
              start = start, count)))
            nc_close(nc)
            cat("\014")
            cat(paste("tmax weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "tmmx/tmmx_", yearlist[j],
              ".nc"))
            tmax <- retry(as.numeric(ncvar_get(nc, varid = "air_temperature",
              start = start, count)))
            nc_close(nc)
            cat("\014")
            cat(paste("rhmin weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "rmin/rmin_", yearlist[j],
              ".nc"))
            rhmin <- retry(as.numeric(ncvar_get(nc, varid = "relative_humidity",
              start = start, count)))
            nc_close(nc)
            cat("\014")
            cat(paste("rhmax weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "rmax/rmax_", yearlist[j],
              ".nc"))
            rhmax <- retry(as.numeric(ncvar_get(nc, varid = "relative_humidity",
              start = start, count)))
            nc_close(nc)
            cat("\014")
            cat(paste("Rain weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "pr/pr_", yearlist[j],
              ".nc"))
            Rain <- retry(as.numeric(ncvar_get(nc, varid = "precipitation_amount",
              start = start, count)))
            nc_close(nc)
            cat("\014")
            cat(paste("solar weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "srad/srad_", yearlist[j],
              ".nc"))
            solar <- retry(as.numeric(ncvar_get(nc, varid = "surface_downwelling_shortwave_flux_in_air",
              start = start, count)))
            nc_close(nc)
            cat("\014")
            cat(paste("Wind weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "vs/vs_", yearlist[j],
              ".nc"))
            cat("\014")
            Wind <- retry(as.numeric(ncvar_get(nc, varid = "wind_speed",
              start = start, count)))
            nc_close(nc)
            cat("\014")
            Tmax <- tmax - 273.15
            Tmin <- tmin - 273.15
          }else{
            cat(paste("reading weather input for ", yearlist[j], " \n", sep = ""))
            cat(paste("tmin weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "/tmmn/tmmn_", yearlist[j],
              ".nc"))
            tmin <- retry(as.numeric(ncvar_get(nc, varid = "air_temperature",
              start = start, count)))
            nc_close(nc)
            cat("\014")
            cat(paste("tmax weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "tmmx/tmmx_", yearlist[j],
              ".nc"))
            tmax <- retry(as.numeric(ncvar_get(nc, varid = "air_temperature",
              start = start, count)))
            nc_close(nc)
            cat("\014")
            cat(paste("rmin weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "rmin/rmin_", yearlist[j],
              ".nc"))
            rhmin <- c(rhmin, retry(as.numeric(ncvar_get(nc, varid = "relative_humidity",
              start = start, count))))
            nc_close(nc)
            cat("\014")
            cat(paste("rmax weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "rmax/rmax_", yearlist[j],
              ".nc"))
            rhmax <- c(rhmax, retry(as.numeric(ncvar_get(nc, varid = "relative_humidity",
              start = start, count))))
            nc_close(nc)
            cat("\014")
            cat(paste("Rain weather input for ", yearlist[j], " \n", sep = ""))
            nc <- nc_open(paste0(baseurl, "pr/pr_", yearlist[j],
              ".nc"))
            Rain <- c(Rain, retry(as.numeric(ncvar_get(nc, varid = "precipitation_amount",
              start = start, count))))
            nc_close(nc)
            nc <- nc_open(paste0(baseurl, "srad/srad_", yearlist[j],
              ".nc"))
            cat("\014")
            cat(paste("solar weather input for ", yearlist[j], " \n", sep = ""))
            solar <- c(solar, retry(as.numeric(ncvar_get(nc, varid = "surface_downwelling_shortwave_flux_in_air",
              start = start, count))))
            nc_close(nc)
            nc <- nc_open(paste0(baseurl, "vs/vs_", yearlist[j],
              ".nc"))
            cat("\014")
            cat(paste("wind weather input for ", yearlist[j], " \n", sep = ""))
            Wind <- c(Wind, retry(as.numeric(ncvar_get(nc, varid = "wind_speed",
              start = start, count))))
            nc_close(nc)
            cat("\014")
            Tmax <- c(Tmax, tmax - 273.15)
            Tmin <- c(Tmin, tmin - 273.15)
          }
        }
      }else{
        cat("extracting weather data \n")
        yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
        nc <- nc_open(paste(spatial, "/tmmx_", yearlist[1], ".nc",
          sep = ""))
        lon <- matrix(ncvar_get(nc, "lon"))
        lat <- matrix(ncvar_get(nc, "lat"))
        lon_1 <- as.numeric(longlat[1])
        lat_1 <- as.numeric(longlat[2])
        dist1 <- abs(lon - lon_1)
        index1 <- which.min(dist1)
        dist2 <- abs(lat - lat_1)
        index2 <- which.min(dist2)
        start <- c(index2, index1, 1)
        count <- c(1, 1, -1)
        for (j in 1:nyears) {
          if (j == 1) {
            cat(paste("reading weather input for ", yearlist[j],
              " \n", sep = ""))
            nc <- nc_open(paste(spatial, "/tmmn_", yearlist[j],
              ".nc", sep = ""))
            tmin <- as.numeric(ncvar_get(nc, varid = "air_temperature",
              start = start, count))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/tmmx_", yearlist[j],
              ".nc", sep = ""))
            tmax <- as.numeric(ncvar_get(nc, varid = "air_temperature",
              start = start, count))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/rmin_", yearlist[j],
              ".nc", sep = ""))
            rhmin <- as.numeric(ncvar_get(nc, varid = "relative_humidity",
              start = start, count))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/rmax_", yearlist[j],
              ".nc", sep = ""))
            rhmax <- as.numeric(ncvar_get(nc, varid = "relative_humidity",
              start = start, count))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/pr_", yearlist[j],
              ".nc", sep = ""))
            Rain <- as.numeric(ncvar_get(nc, varid = "precipitation_amount",
              start = start, count))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/srad_", yearlist[j],
              ".nc", sep = ""))
            solar <- as.numeric(ncvar_get(nc, varid = "surface_downwelling_shortwave_flux_in_air",
              start = start, count))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/vs_", yearlist[j],
              ".nc", sep = ""))
            Wind <- as.numeric(ncvar_get(nc, varid = "wind_speed",
              start = start, count))
            nc_close(nc)
            Tmax <- tmax - 273.15
            Tmin <- tmin - 273.15
          } else {
            cat(paste("reading weather input for ", yearlist[j],
              " \n", sep = ""))
            nc <- nc_open(paste(spatial, "/tmmn_", yearlist[j],
              ".nc", sep = ""))
            tmin <- as.numeric(ncvar_get(nc, varid = "air_temperature",
              start = start, count))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/tmmx_", yearlist[j],
              ".nc", sep = ""))
            tmax <- as.numeric(ncvar_get(nc, varid = "air_temperature",
              start = start, count))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/rmin_", yearlist[j],
              ".nc", sep = ""))
            rhmin <- c(rhmin, as.numeric(ncvar_get(nc, varid = "relative_humidity",
              start = start, count)))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/rmax_", yearlist[j],
              ".nc", sep = ""))
            rhmax <- c(rhmax, as.numeric(ncvar_get(nc, varid = "relative_humidity",
              start = start, count)))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/pr_", yearlist[j],
              ".nc", sep = ""))
            Rain <- c(Rain, as.numeric(ncvar_get(nc, varid = "precipitation_amount",
              start = start, count)))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/srad_", yearlist[j],
              ".nc", sep = ""))
            solar <- c(solar, as.numeric(ncvar_get(nc, varid = "surface_downwelling_shortwave_flux_in_air",
              start = start, count)))
            nc_close(nc)
            nc <- nc_open(paste(spatial, "/vs_", yearlist[j],
              ".nc", sep = ""))
            Wind <- c(Wind, as.numeric(ncvar_get(nc, varid = "wind_speed",
              start = start, count)))
            nc_close(nc)
            Tmax <- c(Tmax, tmax - 273.15)
            Tmin <- c(Tmin, tmin - 273.15)
          }
        }
      }
      # compute clear sky solar for the site of interest, for cloud cover computation below
      cat("running micro_global to get clear sky solar \n")
      micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, timeinterval = 365)
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
      cloud <- (1 - solar / allclearsky) * 100
      cloud[cloud<0]<-0
      cloud[cloud>100]<-100
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
    ida<-ndays
    idayst <- 1

    tzone<-paste("Etc/GMT-12",sep="") # doing it this way ignores daylight savings!
    dim<-length(seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days"))
    maxshades=rep(maxshade,dim)
    minshades=rep(minshade,dim)
    shademax<-maxshades

    if(is.na(ALTITUDES)!=TRUE){

      if(rungads==1){
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
        TMAXX<-as.matrix(Tmaxx)
        TMINN<-as.matrix(Tminn)
      }
      if(scenario!=""){ # to do - climate change scenarios
        TMAXX=TMAXX+TMAXX_diff
        TMINN=TMINN+TMINN_diff
        VP<-VP*VP_diff2 # modify the predicted VP by this factor
      }
      RAINFALL<-Rain+rainoff

      # correct for potential change in RH with elevation-corrected Tair
      es <- WETAIR(db = TMAXX, rh = 100)$esat
      e <- WETAIR(db = Tmax, rh = rhmin)$e
      RHMINN <- (e / es) * 100
      RHMINN[RHMINN>100]<-100
      RHMINN[RHMINN<0]<-0.01
      es <- WETAIR(db = TMINN, rh = 100)$esat
      e <- WETAIR(db = Tmin, rh = rhmin)$e
      RHMAXX <- (e / es) * 100
      RHMAXX[RHMAXX>100]<-100
      RHMAXX[RHMAXX<0]<-0.01
      TMAXX<-TMAXX+adiab_corr_max
      TMINN<-TMINN+adiab_corr_min
      ALLMINTEMPS<-TMINN
      ALLMAXTEMPS<-TMAXX
      ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)

      WNMAXX <- Wind * windfac
      WNMINN <- Wind * windfac

      MAXSHADES<-maxshades
      MINSHADES<-minshades


      REFLS <- rep(REFL, dim)
      PCTWET <- rep(PCTWET, dim)
      soilwet<-RAINFALL
      soilwet[soilwet<=rainwet] = 0
      soilwet[soilwet>0] = 90
      PCTWET<-pmax(soilwet,PCTWET)

      Intrvls<-rep(0,dim)
      Intrvls[1] <- 1 # user-supplied last day-of-year in each time interval sequence
      Numtyps <- 1 # number of substrate types
      Numint <- 1  # number of time intervals
      Nodes <- matrix(data = 0, nrow = 10, ncol = dim) # deepest nodes for each substrate type
      Nodes[1,1] <- 10. # deepest nodes for each substrate type
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
          tannulrun<-movingFun(avetemp,n=365,fun=mean,type='to')
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


      # impose uniform warming
      TMAXX<-TMAXX+warm
      TMINN<-TMINN+warm


      SLES<-matrix(nrow = dim, data = 0)
      SLES<-SLES+SLE

      moists2<-matrix(nrow=10, ncol = ndays, data=0)
      moists2[1,ndays]<-0.2
      moists<-moists2

      if(runmoist==1){
        if(timeinterval==365){
          moists2<-matrix(nrow=10, ncol = dim, data=0) # set up an empty vector for soil moisture values through time
        }else{
          moists2<-matrix(nrow=10, ncol = timeinterval, data=0) # set up an empty vector for soil moisture values through time
        }
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
      microinput<-c(dim,RUF,ERR,Usrhyt,Refhyt,Numtyps,Z01,Z02,ZH1,ZH2,idayst,ida,HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,CMH2O,microdaily,tannul,EC,VIEWF,snowtemp,snowdens,snowmelt,undercatch,rainmult,runshade,runmoist,maxpool,evenrain,snowmodel,rainmelt,writecsv,densfun,hourly,rainhourly,lamb,IUV,RW,PC,RL,SP,R1,IM,MAXCOUNT,IR,message,fail,snowcond,intercept)

      # hourly option set to 0, so make empty vectors
      if(hourly==0){
        TAIRhr=rep(0,24*dim)
        RHhr=rep(0,24*dim)
        WNhr=rep(0,24*dim)
        CLDhr=rep(0,24*dim)
        SOLRhr=rep(0,24*dim)
        ZENhr=rep(-1,24*dim)
      }
      if(rainhourly==0){
        RAINhr=rep(0,24*dim)
      }else{
        RAINhr = rainhour
      }

      doy1=matrix(data = 0., nrow = dim, ncol = 1)
      SLES1=matrix(data = 0., nrow = dim, ncol = 1)
      MAXSHADES1=matrix(data = 0., nrow = dim, ncol = 1)
      MINSHADES1=matrix(data = 0., nrow = dim, ncol = 1)
      TMAXX1=matrix(data = 0., nrow = dim, ncol = 1)
      TMINN1=matrix(data = 0., nrow = dim, ncol = 1)
      CCMAXX1=matrix(data = 0., nrow = dim, ncol = 1)
      CCMINN1=matrix(data = 0., nrow = dim, ncol = 1)
      RHMAXX1=matrix(data = 0., nrow = dim, ncol = 1)
      RHMINN1=matrix(data = 0., nrow = dim, ncol = 1)
      WNMAXX1=matrix(data = 0., nrow = dim, ncol = 1)
      WNMINN1=matrix(data = 0., nrow = dim, ncol = 1)
      REFLS1=matrix(data = 0., nrow = dim, ncol = 1)
      PCTWET1=matrix(data = 0., nrow = dim, ncol = 1)
      RAINFALL1=matrix(data = 0, nrow = dim, ncol = 1)
      tannul1=matrix(data = 0, nrow = dim, ncol = 1)
      moists1=matrix(data = 0., nrow = 10, ncol = dim)
      doy1[1:dim]<-doy
      SLES1[1:dim]<-SLES
      MAXSHADES1[1:dim]<-MAXSHADES
      MINSHADES1[1:dim]<-MINSHADES
      TMAXX1[1:dim]<-TMAXX
      TMINN1[1:dim]<-TMINN
      CCMAXX1[1:dim]<-CCMAXX
      CCMINN1[1:dim]<-CCMINN
      RHMAXX1[1:dim]<-RHMAXX
      RHMINN1[1:dim]<-RHMINN
      WNMAXX1[1:dim]<-WNMAXX
      WNMINN1[1:dim]<-WNMINN
      REFLS1[1:dim]<-REFLS
      PCTWET1[1:dim]<-PCTWET
      RAINFALL1[1:dim]<-RAINFALL
      tannul1[1:dim]<-tannul
      moists1[1:10,1:dim]<-moists
      if(length(LAI)<dim){
        LAI<-rep(LAI[1],dim)
        LAI1 <- LAI
      }
      if(shore==0){
        tides<-matrix(data = 0, nrow = 24*dim, ncol = 3) # make an empty matrix
      }
      # all microclimate data input list - all these variables are expected by the input argument of the fortran micro2014 subroutine
      micro<-list(tides=tides,microinput=microinput,doy=doy,SLES=SLES1,DEP=DEP,Nodes=Nodes,MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TIMAXS=TIMAXS,TIMINS=TIMINS,TMAXX=TMAXX1,TMINN=TMINN1,RHMAXX=RHMAXX1,RHMINN=RHMINN1,CCMAXX=CCMAXX1,CCMINN=CCMINN1,WNMAXX=WNMAXX1,WNMINN=WNMINN1,TAIRhr=TAIRhr,RHhr=RHhr,WNhr=WNhr,CLDhr=CLDhr,SOLRhr=SOLRhr,RAINhr=RAINhr,ZENhr=ZENhr,REFLS=REFLS1,PCTWET=PCTWET1,soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists1,RAINFALL=RAINFALL1,tannulrun=tannulrun,PE=PE,KS=KS,BB=BB,BD=BD,DD=DD,L=L,LAI=LAI1)
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
      }
      if(is.numeric(loc[1])){
        location<-paste("long",loc[1],"lat",loc[2])
      }else{
        location<-loc
      }
      cat(paste('running microclimate model for',timeinterval,'days by',nyears,'years at site',location,'\n'))
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
      if(lamb == 1){
        drlam<-as.data.frame(microut$drlam) # retrieve direct solar irradiance
        drrlam<-as.data.frame(microut$drrlam) # retrieve direct Rayleigh component solar irradiance
        srlam<-as.data.frame(microut$srlam) # retrieve scattered solar irradiance
        if(snowmodel == 1){
          return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,dim=dim,ALTT=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=minshade,maxshade=maxshade,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam))
        }else{
          return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,dim=dim,ALTT=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=minshade,maxshade=maxshade,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam))
        }
      }else{
        if(snowmodel == 1){
          return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,dim=dim,ALTT=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=minshade,maxshade=maxshade,DEP=DEP))
        }else{
          return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,dim=dim,ALTT=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=c(x[1],x[2]),nyears=nyears,timeinterval=timeinterval,minshade=minshade,maxshade=maxshade,DEP=DEP))
        }
      }
    } # end of check for na sites
  } # end error trapping
} # end of micro_usa function
