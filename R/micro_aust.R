#' Australian implementation of the microclimate model.
#'
#' An implementation of the Niche Mapper microclimate model that uses the AWAP daily weather database
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
#' @usage micro_aust(loc = "Melbourne, Australia", timeinterval = 365, ystart = 1990, yfinish = 1990, soiltype = 4,
#' REFL = 0.15, slope = 0, aspect = 0, DEP = c(0., 2.5,  5.,  10.,  15,  20,  30,  50,  100,  200), minshade = 0, maxshade = 90,
#' Usrhyt = 0.01, ...)
#' @export
#' @details
#' \strong{ Parameters controling how the model runs:}\cr\cr
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes, 0=no}\cr\cr
#' \code{IR}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{manualshade}{ = 1, Use CSIRO Soil and Landscape Grid of Australia? 1=yes, 0=no}\cr\cr
#' \code{soildata}{ = 1, Use CSIRO Soil and Landscape Grid of Australia? 1=yes, 0=no}\cr\cr
#' \code{terrain}{ = 0, Use 250m resolution terrain data? 1=yes, 0=no}\cr\cr
#' \code{dailywind}{ = 1, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{windfac}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{adiab_cor}{ = 1, use adiabatic lapse rate correction? 1=yes, 0=no}\cr\cr
#' \code{warm}{ = 0, uniform warming, deg C}\cr\cr
#' \code{spatial}{ = "c:/Australian Environment/", choose location of terrain data}\cr\cr
#' \code{vlsci}{ = 0, running on the VLSCI system? 1=yes, 0=no}\cr\cr
#' \code{opendap}{ = 1, query met grids via opendap (does not work on PC unless you compile ncdf4 - see https://github.com/pmjherman/r-ncdf4-build-opendap-windows)}\cr\cr
#' \code{loop}{ = 0, if doing multiple years, this shifts the starting year by the integer value}\cr\cr

#' \code{soilgrids}{ = 1, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr

#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{ERR}{ = 1.5, Integrator error tolerance for soil temperature calculations}\cr\cr
#' \code{Refhyt}{ = 1.2, Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured}\cr\cr
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
#' \code{lapse_min}{ = 0.0039 Lapse rate for minimum air temperature (degrees C/m)}\cr\cr
#' \code{lapse_max}{ = 0.0077 Lapse rate for maximum air temperature (degrees C/m)}\cr\cr
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
#' \code{snowdens}{ = 0.375, snow density (mg/m3), overridden by }
#' \code{densfun}\cr\cr
#' \code{densfun}{ = c(0,0), slope and intercept of linear model of snow density as a function of day-of-year - if it is c(0,0) then fixed density used}\cr\cr
#' \code{snowmelt}{ = 0.9, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
#' \code{undercatch}{ = 1, undercatch multipier for converting rainfall to snow}\cr\cr
#' \code{rainmelt}{ = 0.0125, paramter in equation that melts snow with rainfall as a function of air temp}\cr\cr
#'
#' \strong{ Intertidal mode parameters:}
#'
#' \code{shore}{ Include tide effects? If 1, the matrix}
#' \code{tides}
#' { is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0., nrow = 24*timeinterval*nyears, ncol = 3), matrix for each how of the simulation of 1. tide state (0=out, 1=in), 2. Water temperature (deg C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#'
#' \strong{Outputs:}
#' metout/shadmet variables:
#' \itemize{
#' \item 1 DOY - day-of-year
#' \item 2 TIME - time of day (mins)
#' \item 3 TALOC - air temperature (deg C) at local height (specified by 'Usrhyt' variable)
#' \item 4 TAREF - air temperature (deg C) at reference height (specified by 'Refhyt', 1.2m default)
#' \item 5 RHLOC - relative humidity (\%) at local height (specified by 'Usrhyt' variable)
#' \item 6 RH  - relative humidity (\%) at reference height (specified by 'Refhyt', 1.2m default)
#' \item 7 VLOC - wind speed (m/s) at local height (specified by 'Usrhyt' variable)
#' \item 8 VREF - wind speed (m/s) at reference height (specified by 'Refhyt', 1.2m default)
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
#' \item  3-10 SN0cm ... - snow temperature (deg C), at the soil surface and each of the potential 8 layers
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
#'micro<-micro_aust() # run the model with default location and settings
#'
#'metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#'shadmet<-as.data.frame(micro$shadmet) # above ground microclimatic conditions, max shade
#'soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
#'shadsoil<-as.data.frame(micro$shadsoil) # soil temperatures, maximum shade
#'
#'# append dates
#'days<-rep(seq(1,12),24)
#'days<-days[order(days)]
#'dates<-days+metout$TIME/60/24-1 # dates for hourly output
#'dates2<-seq(1,12,1) # dates for daily output
#'
#'plotmetout<-cbind(dates,metout)
#'plotsoil<-cbind(dates,soil)
#'plotshadmet<-cbind(dates,shadmet)
#'plotshadsoil<-cbind(dates,shadsoil)
#'
#'minshade<-micro$minshade
#'maxshade<-micro$maxshade
#'
#'# plotting above-ground conditions in minimum shade
#'with(plotmetout,{plot(TALOC ~ dates,xlab = "Date and Time", ylab = "Air Temperature (deg C)"
#', type = "l",main=paste("air temperature, ",minshade,"% shade",sep=""))})
#'with(plotmetout,{points(TAREF ~ dates,xlab = "Date and Time", ylab = "Air Temperature (deg C)"
#', type = "l",lty=2,col='blue')})
#'with(plotmetout,{plot(RHLOC ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",ylim=c(0,100),main=paste("humidity, ",minshade,"% shade",sep=""))})
#'with(plotmetout,{points(RH ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",col='blue',lty=2,ylim=c(0,100))})
#'with(plotmetout,{plot(TSKYC ~ dates,xlab = "Date and Time", ylab = "Sky Temperature (deg C)"
#',  type = "l",main=paste("sky temperature, ",minshade,"% shade",sep=""))})
#'with(plotmetout,{plot(VREF ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#',  type = "l",main="wind speed",ylim = c(0, 15))})
#'with(plotmetout,{points(VLOC ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#',  type = "l",lty=2,col='blue')})
#'with(plotmetout,{plot(ZEN ~ dates,xlab = "Date and Time", ylab = "Zenith Angle of Sun (deg)"
#',  type = "l",main="solar angle, sun")})
#'with(plotmetout,{plot(SOLR ~ dates,xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
#',  type = "l",main="solar radiation")})
#'
#'# plotting soil temperature for minimum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(plotsoil[,i+3]~plotsoil[,1],xlab = "Date and Time", ylab = "Soil Temperature (deg C)"
#'    ,col=i,type = "l",main=paste("soil temperature ",minshade,"% shade",sep=""))
#'  }else{
#'    points(plotsoil[,i+3]~plotsoil[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (deg C)",col=i,type = "l")
#'  }
#'}
#'
#'# plotting above-ground conditions in maximum shade
#'with(plotshadmet,{plot(TALOC ~ dates,xlab = "Date and Time", ylab = "Air Temperature (deg C)"
#', type = "l",main="air temperature, sun")})
#'with(plotshadmet,{points(TAREF ~ dates,xlab = "Date and Time", ylab = "Air Temperature (deg C)"
#', type = "l",lty=2,col='blue')})
#'with(plotshadmet,{plot(RHLOC ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",ylim=c(0,100),main="humidity, shade")})
#'with(plotshadmet,{points(RH ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",col='blue',lty=2,ylim=c(0,100))})
#'with(plotshadmet,{plot(TSKYC ~ dates,xlab = "Date and Time", ylab = "Sky Temperature (deg C)",
#'  type = "l",main="sky temperature, shade")})
#'
#'# plotting soil temperature for maximum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(plotshadsoil[,i+3]~plotshadsoil[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (deg C)",col=i,type = "l",main=paste("soil temperature ",maxshade,"% shade",sep=""))
#'  }else{
#'    points(plotshadsoil[,i+3]~plotshadsoil[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (deg C)",col=i,type = "l")
#'  }
#'}
micro_aust <- function(loc="Nyrripi, Northern Territory",timeinterval=365,ystart=1990,yfinish=1990,
  nyears=1,soiltype=4,REFL=0.15, elev = NA, slope=0,aspect=0, lapse_max = 0.0077, lapse_min = 0.0039,
  DEP=c(0., 2.5,  5.,  10.,  15,  20,  30,  50,  100,  200),
  minshade=0,maxshade=90,Refhyt=1.2,Usrhyt=0.01,Z01=0,Z02=0,ZH1=0,ZH2=0,
  runshade=1,clearsky=0,run.gads=1,write_input=0,writecsv=0,manualshade=1,
  soildata=1,terrain=0,dailywind=1,windfac=1,adiab_cor=1,warm=0,spatial="c:/Australian Environment/",vlsci=0,
  ERR=1.5,RUF=0.004,EC=0.0167238,SLE=0.95,Thcond=2.5,Density=2.56,SpecHeat=870,BulkDensity=1.3,
  PCTWET=0,rainwet=1.5,cap=1,CMH2O=1.,hori=rep(0,24),
  TIMAXS=c(1.0, 1.0, 0.0, 0.0),TIMINS=c(0, 0, 1, 1),timezone=0,
  runmoist=1,PE=rep(1.1,19),KS=rep(0.0037,19),BB=rep(4.5,19),BD=rep(BulkDensity,19),DD=rep(Density,19),
  maxpool=10000,rainmult=1,evenrain=0,
  SoilMoist_Init=c(0.1,0.12,0.15,0.3,0.4,0.4,0.4,0.4,0.4,0.4),
  L=c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,
    4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,
    0,0)*10000, R1 = 0.001, RW = 2.5e+10, RL = 2e+6, PC = -1500, SP = 10, IM = 1e-06, MAXCOUNT = 500,
  LAI=0.1,
  snowmodel=0,snowtemp=1.5,snowdens=0.375,densfun=c(0,0),snowmelt=0.9,undercatch=1,rainmelt=0.0125,
  shore=0,tides=matrix(data = 0., nrow = 24*timeinterval*nyears, ncol = 3),loop=0, scenario="",year="",
  barcoo="",quadrangle=1,hourly=0,rainhourly=0,rainhour=0, uid = "", pwd = "",
  lamb = 0, IUV = 0, soilgrids = 1, IR = 0, opendap = 1, message = 0, fail = nyears * 24 * 365) {

  #   loc="Nyrripi, Northern Territory"
  #   timeinterval=365
  #   ystart=2000
  #   yfinish=2001
  #   nyears=yfinish-ystart+1
  #   soiltype=4
  #   REFL=0.15
  #   slope=0
  #   aspect=0
  #   DEP=c(0., 2.5,  5.,  10.,  15.,  20.,  30.,  50.,  100.,  200.)
  #   minshade=0
  #   maxshade=90
  #   Usrhyt=.01
  #   Z01=0
  #   Z02=0
  #   ZH1=0
  #   ZH2=0
  #   runshade=1
  #   clearsky=0
  #   run.gads=1
  #   write_input=0
  #   writecsv=0
  #   manualshade=1
  #   soildata=0
  #   terrain=0
  #   dailywind=1
  #   adiab_cor=1
  #   warm=0
  #   spatial="w:/"
  #   vlsci=0
  #   loop=0
  #   ERR=1.5
  #   RUF=0.004
  #   EC=0.0167238
  #   SLE=0.95
  #   Thcond=2.5
  #   Density=2.560
  #   SpecHeat=870
  #   BulkDensity=1.300
  #   PCTWET=0
  #   rainwet=1.5
  #   cap=1
  #   CMH2O=1
  #   hori=rep(0,24)
  #   TIMAXS=c(1.0, 1.0, 0.0, 0.0)
  #   TIMINS=c(0, 0, 1, 1)
  #   timezone=0
  #   runmoist=1
  #   PE=rep(1.1,19)
  #   KS=rep(0.0037,19)
  #   BB=rep(4.5,19)
  #   BD=rep(1.3,19)
  #   DD=rep(2.56,19)
  #   maxpool=10000
  #   rainmult=1
  #   evenrain=0
  #   SoilMoist_Init=c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3)
  #   L=c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,
  #       4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,
  #       0,0)*10000
  #   LAI=0.1
  #   snowmodel=0
  #   snowtemp=1.5
  #   snowdens=0.375
  #   densfun=c(0,0)
  #   snowmelt=0.9
  #   undercatch=1
  #   rainmelt=0.0125
  #   shore=0
  #   tides=matrix(data = 0., nrow = 24*timeinterval*nyears, ncol = 3)
  #   scenario=""
  #   year=2070
  #   barcoo=""
  #   quadrangle=1
  #   hourly=0
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
  # Refhyt <- 1.2

  if(vlsci==0 | opendap==0){
    library(RODBC)
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
  if(soiltype<0 | soiltype>11){
    message("ERROR: the soil type must range between 1 and 11.
      Please correct.", '\n')
    errors<-1
  }
  # end error trapping

  if(errors==0){ # continue

    ################## time related variables #################################
    nyears<-yfinish-ystart+1
    if(vlsci==1){
      ndays<-365*nyears
    }

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

    # now check if doing something other than middle day of each month, and create appropriate vector of day-of-year
    daystart<-as.integer(ceiling(365/timeinterval/2))
    if(timeinterval!=12){
      doys<-seq(daystart,365,as.integer(floor(365/timeinterval)))
    }else{
      doys<-doysn
    }
    doynum <- timeinterval*nyears # total days to do
    doy <- subset(doys, doys!=0) # final vector of day-of-year
    doy<-rep(doy,nyears)
    dim<-doynum
    maxshades=rep(maxshade,dim)
    minshades=rep(minshade,dim)
    doys<-seq(daystart,dim,1)
    if(opendap == 1){
    leapyears<-seq(1972,2060,4)
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
    dim<-length(doy)
    }
    idayst <- 1 # start day
    ida<-timeinterval*nyears # end day
    dates<-Sys.time()-60*60*24
    curyear<-as.numeric(format(dates,"%Y"))

    ################## location and terrain #################################
    if(vlsci==0){
      f1 <- paste(spatial,"ausclim_rowids.nc",sep="");
      f2 <- paste(spatial,"ausdem_shift1.tif",sep="");
      f3 <- paste(spatial,"agg_9secdem.nc",sep="");
      f4 <- paste(spatial,"Aust9secDEM.tif",sep="");
    }else{
      f1 <- "/vlsci/VR0212/shared/Spatial_Data/ausclim_rowids.nc";
      f2 <- "/vlsci/VR0212/shared/Spatial_Data/ausdem_shift1.tif";
      f3 <- "/vlsci/VR0212/shared/Spatial_Data/agg_9secdem.nc";
      f4 <- "/vlsci/VR0212/shared/Spatial_Data/Aust9secDEM.grd";
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
      longlat <- dismo::geocode(loc)[3:4] # assumes first geocode match is correct
      if(nrow(longlat>1)){longlat<-longlat[1,]}
      x <- t(as.matrix(as.numeric(c(longlat[1,1],longlat[1,2]))))
    }else{
      longlat <- loc
      x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))
    }

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
      MAXSHADES <- rep(0,dim)+maxshade # daily max shade (%)
      MINSHADES <- rep(0,dim)+minshade # daily min shade (%)
    }

    if(soildata==1){
      message("extracting soil data \n")
      if(vlsci==0){
        static_soil<-paste(spatial,"static_soil.nc",sep="")
        emissivities<-paste(spatial,"aus_emissivities.nc",sep="")
      }else{
        static_soil<-'/vlsci/VR0212/shared/Spatial_Data/static_soil.nc'
        emissivities<-'/vlsci/VR0212/shared/Spatial_Data/aus_emissivities.nc'
      }
      # read data in from netcdf file
      static_soil_data<-raster::brick(static_soil)
      static_soil_vars <- raster::extract(static_soil_data,x)
      labels<-c('albedo','FAPAR1','FAPAR2','FAPAR3','FAPAR4','FAPAR5','FAPAR6','FAPAR7','FAPAR8','FAPAR9','FAPAR10','FAPAR11','FAPAR12','volwater_Upper','volwater_lower','thick_upper','thick_lower','code')
      colnames(static_soil_vars)<-labels
      emissivities_data<-raster::brick(emissivities)
      SLES2 <- raster::extract(emissivities_data,x)

      # read in other soil related files for working out lumped soil type and properties
      # such as clay % for getting water potential
      if(vlsci==0){
        filename<-paste(spatial,"ppfInterpAll.txt",sep="")
        ppf<-as.data.frame(read.table(file = filename, sep = ",", header=TRUE))
        filename<-paste(spatial,"Lumped soil types.txt",sep="")
        lumped.soil<-as.data.frame(read.table(file = filename, sep = ","))
        filename<-paste(spatial,"SoilTypeLUT_725_AWAP.csv",sep="")
        soiltype<-as.data.frame(read.table(file = filename, sep = ","))
      }else{
        filename<-'/vlsci/VR0212/shared/Spatial_Data/ppfInterpAll.txt'
        ppf<-as.data.frame(read.table(file = filename, sep = ",", header=TRUE))
        filename<-'/vlsci/VR0212/shared/Spatial_Data/Lumped soil types.txt'
        lumped.soil<-as.data.frame(read.table(file = filename, sep = ","))
        filename<-'/vlsci/VR0212/shared/Spatial_Data/SoilTypeLUT_725_AWAP.csv'
        soiltype<-as.data.frame(read.table(file = filename, sep = ","))
      }
      soilcode<-subset(soiltype, soiltype[1]==static_soil_vars[18])
      lumped<-subset(lumped.soil, V4==as.character(soilcode[1,2]))
      soiltype<-lumped[1,6]
      soilprop<-subset(ppf, ppf==soilcode[1,2])
    }else{
      SLES2 <- rep(SLE,timeinterval*nyears)
      if(manualshade==0){
        message("extracting shade data \n")
        if(vlsci==0){
          static_soil<-paste(spatial,"static_soil.nc",sep="")
          emissivities<-paste(spatial,"aus_emissivities.nc",sep="")
        }else{
          static_soil<-'/vlsci/VR0212/shared/Spatial_Data/static_soil.nc'
          emissivities<-'/vlsci/VR0212/shared/Spatial_Data/aus_emissivities.nc'
        }
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
        if(vlsci==0){
          horifile<-paste(spatial,'horizon',i,'.tif',sep="")
        }else{
          horifile<-paste('/vlsci/VR0212/shared/Spatial_Data/','horizon',i,'.tif',sep="")
        }
        horiz<-raster::crop(raster::raster(horifile),e)
        if(i==1){
          horizons_data<-horiz
        }else{
          horizons_data<-raster::stack(horizons_data,horiz)
        }
      }
      HORIZONS <- t(raster::extract(horizons_data,x))
      if(vlsci==0){
        elev1<-crop(raster::raster(paste(spatial,'elev.tif',sep="")),e)
        slope1<-crop(raster::raster(paste(spatial,'slope.tif',sep="")),e)
        aspect1<-crop(raster::raster(paste(spatial,'aspect.tif',sep="")),e)
        elevslpasp<-raster::stack(elev1,slope1,aspect1)
      }else{
        elev1<-raster::crop(raster(paste('/vlsci/VR0212/shared/Spatial_Data/','elev.tif',sep="")),e)
        slope1<-raster::crop(raster(paste('/vlsci/VR0212/shared/Spatial_Data/','slope.tif',sep="")),e)
        aspect1<-raster::crop(raster(paste('/vlsci/VR0212/shared/Spatial_Data/','aspect.tif',sep="")),e)
        elevslpasp<-raster::stack(elev1,slope1,aspect1)
      }
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
      #ALTITUDES <- AUSDEM
      #message("using 0.05 res DEM!")
      #}
      HORIZONS <- hori
      HORIZONS <- data.frame(HORIZONS)
      VIEWF_all <- rep(1,length(x[,1]))
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
    # setting up for temperature correction using lapse rate given difference between 9sec DEM value and 0.05 deg value
    if(AUSDEM==-9999 | is.na(AUSDEM)=='TRUE'){
      delta_elev = AGG - ALTITUDES
    }else{
      delta_elev = AUSDEM - ALTITUDES
    }
    adiab_corr_max <- delta_elev * lapse_max
    adiab_corr_min <- delta_elev * lapse_min

    if(scenario!=""){
      message("generate climate change scenario", '\n')
      # diff spline function
      getdiff<-function(diffs,grid){
        diff1<-(unlist(diffs[1])+unlist(diffs[12]))/2

        # generate list of days
        for(ys in 1:nyears){
          day<-c(1,15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 365)
          day.leap<-c(1,15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 366)
          if(ys==1){
            days2=day
            days=day
          }else{
            days2=c(days2,(day+365*(ys-1)))
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

      ########### Max and Min Air Temps ################

      load(file = paste("c:/Spatial_Data/Australia Climate Change/",scenario,"/","maxTst05_",scenario,"_",year,".Rda",sep="")) #maxTst05

      diffs<-extract(maxTst05,x)
      TMAXX_diff<-getdiff(diffs,maxTst05)

      load(file = paste("c:/Spatial_Data/Australia Climate Change/",scenario,"/","minTst05_",scenario,"_",year,".Rda",sep="")) #minTst05

      diffs<-extract(minTst05,x)
      TMINN_diff<-getdiff(diffs,minTst05)

      ################ RH ############################

      load(file = paste("c:/Spatial_Data/Australia Climate Change/",scenario,"/","RHst05_",scenario,"_",year,".Rda",sep="")) #maxTst05

      diffs<-extract(RHst05,x)
      RH_diff<-getdiff(diffs,RHst05)

      ################ wind ############################

      load(file = paste("c:/Spatial_Data/Australia Climate Change/",scenario,"/","PT_VELst05_",scenario,"_",year,".Rda",sep=""))

      diffs<-extract(PT_VELst05,x)
      WIND_diff<-getdiff(diffs,PT_VELst05)

      ############# SOLAR/CLOUD COVER ##################

      load(file = paste("c:/Spatial_Data/Australia Climate Change/",scenario,"/","SOLCst05_",scenario,"_",year,".Rda",sep=""))

      diffs<-extract(SOLCst05,x)
      SOLAR_diff<-getdiff(diffs,SOLCst05)
    }

    if(opendap == 1){
      message("extracting climate data", '\n')

      baseurl<-"http://dapds00.nci.org.au/thredds/dodsC/ub8/au/climate/"

      yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
      for (j in 1:nyears) {
        if (j == 1) {
          cat(paste("reading weather input for ", yearlist[j]," \n", sep = ""))
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.rad.", yearlist[j],
            ".nc"))
          lon <- ncvar_get(nc, "longitude")
          lat <- ncvar_get(nc, "latitude")
          flat=match(abs(lat-x[2])<1/44,1)
          latindex=which(flat %in% 1)
          flon=match(abs(lon-x[1])<1/44,1)
          lonindex=which(flon %in% 1)
          start <- c(latindex,lonindex,1)
          count <- c(1, 1, -1)

          sol <- as.numeric(ncvar_get(nc, varid = "rad",
            start = start, count))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.tmax.", yearlist[j],
            ".nc"))
          tmax <- as.numeric(ncvar_get(nc, varid = "tmax",
            start = start, count))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.tmin.", yearlist[j],
            ".nc"))
          tmin <- as.numeric(ncvar_get(nc, varid = "tmin",
            start = start, count))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.vph09.", yearlist[j],
            ".nc"))
          vph09 <- as.numeric(ncvar_get(nc, varid = "vph09",
            start = start, count))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.vph15.", yearlist[j],
            ".nc"))
          vph15 <- as.numeric(ncvar_get(nc, varid = "vph15",
            start = start, count))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.rain.", yearlist[j],
            ".nc"))
          Rain <- as.numeric(ncvar_get(nc, varid = "rain",
            start = start, count))
          nc_close(nc)
        }else{

          cat(paste("reading weather input for ", yearlist[j],
            " \n", sep = ""))
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.rad.", yearlist[j],
            ".nc"))
          sol <- c(sol, as.numeric(ncvar_get(nc, varid = "rad",
            start = start, count)))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.tmax.", yearlist[j],
            ".nc"))
          tmax <- c(tmax, as.numeric(ncvar_get(nc, varid = "tmax",
            start = start, count)))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.tmin.", yearlist[j],
            ".nc"))
          tmin <- c(tmin, as.numeric(ncvar_get(nc, varid = "tmin",
            start = start, count)))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.vph09.", yearlist[j],
            ".nc"))
          vph09 <- c(vph09, as.numeric(ncvar_get(nc, varid = "vph09",
            start = start, count)))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.vph15.", yearlist[j],
            ".nc"))
          vph15 <- c(vph15, as.numeric(ncvar_get(nc, varid = "vph15",
            start = start, count)))
          nc_close(nc)
          nc <- nc_open(paste0(baseurl, "AGCD.BoM.daily.rain.", yearlist[j],
            ".nc"))
          Rain <- c(Rain, as.numeric(ncvar_get(nc, varid = "rain",
            start = start, count)))
          nc_close(nc)
        }
      }
      # compute clear sky solar for the site of interest, for cloud cover computation below
      cat("running micro_global to get clear sky solar \n")
      TAI<-c(0.0670358341290886,0.0662612704779235,0.065497075238002,0.0647431301168489,0.0639993178022531,0.0632655219571553,0.0625416272145492,0.0611230843885423,0.0597427855962549,0.0583998423063099,0.0570933810229656,0.0558225431259535,0.0545864847111214,0.0533843764318805,0.0522154033414562,0.0499736739981675,0.047855059159556,0.0458535417401334,0.0439633201842001,0.0421788036108921,0.0404946070106968,0.0389055464934382,0.0374066345877315,0.0359930755919066,0.0346602609764008,0.0334037648376212,0.0322193394032758,0.0311029105891739,0.0300505736074963,0.0290585886265337,0.0281233764818952,0.0272415144391857,0.0264097320081524,0.0256249068083005,0.0248840604859789,0.0241843546829336,0.0235230870563317,0.0228976873502544,0.0223057135186581,0.0217448478998064,0.0212128934421699,0.0207077699817964,0.0202275105711489,0.0197702578594144,0.0193342605242809,0.0189178697551836,0.0177713140039894,0.0174187914242432,0.0170790495503944,0.0167509836728154,0.0164335684174899,0.0161258546410128,0.0158269663770596,0.0155360978343254,0.0152525104459325,0.0149755299703076,0.0147045436435285,0.0144389973831391,0.0141783930434343,0.0134220329447663,0.0131772403830191,0.0129356456025128,0.0126970313213065,0.0124612184223418,0.0122280636204822,0.01199745718102,0.0115436048739351,0.0110993711778668,0.0108808815754663,0.0106648652077878,0.0104513876347606,0.0102405315676965,0.00982708969547694,0.00962473896278535,0.00903679230300494,0.00884767454432418,0.0083031278398166,0.00796072474935954,0.00755817587626185,0.00718610751850881,0.00704629977586921,0.00684663903049612,0.00654155580333479,0.00642947339729728,0.00627223096874308,0.00603955966866779,0.00580920937536261,0.00568506186880564,0.00563167068287251,0.00556222005081865,0.00550522989971023,0.00547395763028062,0.0054478983436216,0.00541823364504573,0.00539532163908382,0.00539239864119488,0.00541690124712384,0.00551525885358836,0.00564825853509463,0.00577220185074264,0.00584222986640171,0.00581645238345584,0.00566088137411449,0.00535516862329704,0.00489914757707667,0.00432017939770409,0.0036813032251836,0.00309019064543606,0.00270890436501562,0.00276446109239711,0.00356019862584603)
      micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, TAI = TAI, timeinterval = 365)
      clearskyrad <- micro_clearsky$metout[,c(1, 13)]
      clearsky_mean1 <- aggregate(clearskyrad[,2], by = list(clearskyrad[,1]), FUN = sum)[,2]
      leapyears<-seq(1972,2060,4)
      for(j in 1:nyears){
        if(yearlist[j]%in%leapyears){# add day for leap year if needed
          clearsky_mean<-c(clearsky_mean1[1:59],clearsky_mean1[59],clearsky_mean1[60:365])
        }else{
          clearsky_mean <- clearsky_mean1
        }
        if(j==1){
          allclearsky <- clearsky_mean
        }else{
          allclearsky <- c(allclearsky, clearsky_mean)
        }
      }
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
      ndays<-length(TMAXX)
      dim<-ndays
      VAPRES_max<-apply(cbind(vph09, vph15), FUN = max, MARGIN = 1)*100 # convert from hectopascals to pascals
      VAPRES_min<-apply(cbind(vph09, vph15), FUN = min, MARGIN = 1)*100 # convert from hectopascals to pascals
      TMAXK<-TMAXX+273.15
      loge<-TMAXK
      loge2<-loge
      loge[loge2>273.16]<- -7.90298*(373.16/TMAXK[loge2>273.16]-1.)+5.02808*log10(373.16/TMAXK[loge2>273.16])-1.3816E-07*(10.^(11.344*(1.-TMAXK[loge2>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/TMAXK[loge2>273.16]-1.))-1.)+log10(1013.246)
      loge[loge2<=273.16]<- -9.09718*(273.16/TMAXK[loge2<=273.16]-1.)-3.56654*log10(273.16/TMAXK[loge2<=273.16])+.876793*(1.-TMAXK[loge2<=273.16]/273.16)+log10(6.1071)
      estar<-(10.^loge)*100.
      RHMINN<-(VAPRES_min/estar)*100
      RHMINN[RHMINN>100]<-100
      RHMINN[RHMINN<0]<-0.01
      #RHMINN
      TMINK<-TMINN+273.15
      loge<-TMINK
      loge2<-loge
      loge[loge2>273.16]<- -7.90298*(373.16/TMAXK[loge2>273.16]-1.)+5.02808*log10(373.16/TMAXK[loge2>273.16])-1.3816E-07*(10.^(11.344*(1.-TMAXK[loge2>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/TMAXK[loge2>273.16]-1.))-1.)+log10(1013.246)
      loge[loge2<=273.16]<- -9.09718*(273.16/TMAXK[loge2<=273.16]-1.)-3.56654*log10(273.16/TMAXK[loge2<=273.16])+.876793*(1.-TMAXK[loge2<=273.16]/273.16)+log10(6.1071)
      estar<-(10.^loge)*100.
      RHMAXX<-(VAPRES_max/estar)*100
      RHMAXX[RHMAXX>100]<-100
      RHMAXX[RHMAXX<0]<-0.01
      WNMAXX <- rep(2, length(ndays))
      WNMINN <- rep(0.5, length(ndays))
      RAINFALL <- Rain
    }

    # connect to server
    if(vlsci==0 & opendap == 0){
      channel2 <- RODBC::odbcConnect("ausclim_predecol", uid = uid, pwd = pwd)
      channel <- RODBC::odbcConnect("AWAPDaily", uid = uid, pwd = pwd)


      # preliminary test for incomplete year, if simulation includes the present year
      yearlist<-seq(ystart,(ystart+(nyears-1)),1)
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
        if(nrow(output)>365){
          # fix leap years
          output<-output[-60,]
        }
        if(j==1){
          results<-output
        }else{
          results<-rbind(results,output)
        }
      }

      nyears2<-nrow(results)/365
      ndays<-nrow(results)
      doysn2<-doysn[doysn<=ndays]
      doysn2<-doysn[1:round(nyears2*12)]
    } #end vlsci check
    dim<-ndays
    maxshades=rep(maxshade,dim)
    minshades=rep(minshade,dim)
    doys<-seq(daystart,dim,1)
    leapyears<-seq(1972,2060,4)
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
    #doy<-rep(seq(1,365),nyears)
    doy<-doy[1:dim]
    ida<-ndays
    idayst <- 1 # start month
    # end preliminary test for incomplete year, if simulation includes the present year

    if((soildata==1 & nrow(soilprop)>0)|soildata==0){

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

      if(is.na(dbrow)!=TRUE & is.na(ALTITUDES)!=TRUE){

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


        if(vlsci==0 & opendap==0){
          yearlist<-seq(ystart,(ystart+(nyears-1)),1)
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

            if(nrow(output)>365){
              # fix leap years
              output<-output[-60,]
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

                if(nrow(output)>365){
                  # fix leap years
                  output<-output[-60,]
                }
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

                if(nrow(output)>365){
                  # fix leap years
                  output<-output[-60,]
                }
              }
              if(j==1){
                dwind<-output
              }else{
                dwind<-rbind(dwind,output)
              }
            }
            dwind<-dwind$wind/15.875
          }
        }else{ #vlsci==1
          if(vlsci == 1){
            load(paste(barcoo,'TMAXX.bin',sep=''))
            TMAXX <- as.matrix(data[quadrangle,2:7301])
            load(paste(barcoo,'TMINN.bin',sep=''))
            TMINN <- as.matrix(data[quadrangle,2:7301])
            load(paste(barcoo,'RHMAXX.bin',sep=''))
            RHMAXX <- as.numeric(as.matrix(data[quadrangle,2:7301]))
            load(paste(barcoo,'RHMINN.bin',sep=''))
            RHMINN <- as.numeric(as.matrix(data[quadrangle,2:7301]))
            load(paste(barcoo,'CCMAXX.bin',sep=''))
            CCMAXX <- as.numeric(as.matrix(data[quadrangle,2:7301]))
            load(paste(barcoo,'CCMINN.bin',sep=''))
            CCMINN <- as.numeric(as.matrix(data[quadrangle,2:7301]))

            load(paste(barcoo,'RAINFALL.bin',sep=''))
            RAINFALL <- as.matrix(data[quadrangle,2:7301])

            RHMAXX<-t(RHMAXX)
            RHMINN<-t(RHMINN)
            CCMAXX<-t(CCMAXX)
            CCMINN<-t(CCMINN)
            TMAXX<-t(as.data.frame(TMAXX[((ystart-1989)*365-364):((yfinish-1989)*365)]))
            TMINN<-t(as.data.frame(TMINN[((ystart-1989)*365-364):((yfinish-1989)*365)]))
            RHMAXX<-t(as.data.frame(RHMAXX[((ystart-1989)*365-364):((yfinish-1989)*365)]))
            RHMINN<-t(as.data.frame(RHMINN[((ystart-1989)*365-364):((yfinish-1989)*365)]))
            CCMAXX<-t(as.data.frame(CCMAXX[((ystart-1989)*365-364):((yfinish-1989)*365)]))
            CCMINN<-t(as.data.frame(CCMINN[((ystart-1989)*365-364):((yfinish-1989)*365)]))
            RAINFALL<-t(as.data.frame(RAINFALL[((ystart-1989)*365-364):((yfinish-1989)*365)]))

            if(dailywind==1){

              load(paste(barcoo,'dwind.bin',sep=''))
              dwind<- as.matrix(data[quadrangle,2:7301])
              dwind<-t(as.data.frame(dwind[((ystart-1989)*365-364):((yfinish-1989)*365)]))
              #WNMINN<-WNMAXX
            }
            if(adiab_cor==1){
              TMAXX<-TMAXX+adiab_corr_max
              TMINN<-TMINN+adiab_corr_min
            }
          }
        } #end vlsci check
        if(vlsci==0){
          if(opendap == 0){
            if(adiab_cor==1){
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
            dim<-length(RAINFALL)
            output_AWAPDaily<-results
          }else{
            if(adiab_cor==1){
              TMAXX<-as.matrix(TMAXX+adiab_corr_max)
              TMINN<-as.matrix(TMINN+adiab_corr_min)
            }
            dim <- length(RAINFALL)
          }
          if(scenario!=""){
            # first work out for each site the new predicted rainfall amount for each month - use this to adjust for fact that will underestimate chcange
            # using proportion becasue 0 x % is still 0
            # add columns with days, months and years
            RAIN_current<-as.data.frame(RAINFALL)
            dates<-seq(ISOdate(ystart,1,1,tz=paste("Etc/GMT-",10,sep=""))-3600*12, ISOdate((ystart+nyears),1,1,tz=paste("Etc/GMT-",10,sep=""))-3600*13, by="days")
            dates=subset(dates, format(dates, "%m/%d")!= "02/29") # remove leap years
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

        }

        # cloud cover
        if(opendap==0){
          if(vlsci==0){
            if(ystart>1989 & sum(results[,9],na.rm=TRUE)>0){ # solar radiation data available
              query<-paste("SELECT a.*
                       FROM [ausclim].[dbo].[clearskysol] as a
                       where (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
                       ",sep="")
              output_ausclim<- sqlQuery(channel,query)

              if(nrow(output_ausclim)==0){ #no satellite coverage, get data from ausclim
                clouds<-paste("select cloud1,cloud2,cloud3,cloud4,cloud5,cloud6,cloud7,cloud8,cloud9,cloud10,cloud11,cloud12 FROM cloudcover WHERE i = ",dbrow,sep="")
                #CCMAXX <- dbGetQuery(con,statement=clouds)*100
                if(vlsci==0){
                  CCMAXX<- sqlQuery(channel2,clouds)*100
                }
                CCMINN <- CCMAXX
                CCMAXX1 <-suppressWarnings(spline(doys12,CCMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
                CCMAXX <- rep(CCMAXX1$y,nyears)
                CCMINN <- CCMAXX
              }else{
                weekly_sol<-cbind(1:52,t(output_ausclim[3:54]))
                daily_sol <-suppressWarnings(spline(seq(3,361,7),weekly_sol[,2],n=365,xmin=1,xmax=365,method="periodic"))
                daily_sol<-as.numeric(daily_sol$y)
                if(yfinish==curyear){ # if it's the current year, then make sure daily_sol goes however far into the year we are
                  daily_sol<-c(rep(daily_sol,nyears-1),daily_sol[1:(ndays-(nyears-1)*365)])
                }else{
                  daily_sol<-rep(daily_sol,nyears)
                }
                if(is.na(output_AWAPDaily[1,9])==TRUE){
                  output_AWAPDaily[1,9]=mean(output_AWAPDaily[,9],na.rm=TRUE)
                }
                if(is.na(output_AWAPDaily[dim,9])==TRUE){
                  output_AWAPDaily[nrow(output_AWAPDaily),9]=mean(output_AWAPDaily[,9],na.rm=TRUE)
                }
                solar<-zoo::na.approx(output_AWAPDaily[,9])
                if(scenario!=""){
                  solar=solar*SOLAR_diff
                }
                cloud<-(1-as.data.frame(solar)/as.data.frame(daily_sol))*100
                cloud[cloud<0]<-0
                cloud[cloud>100]<-100
                cloud<-as.matrix(cbind(output_AWAPDaily[,4],cloud))
                CCMAXX<-cloud[,2]
                CCMINN<-CCMAXX
              }
            }else{
              #           clouds<-paste("select cloud1,cloud2,cloud3,cloud4,cloud5,cloud6,cloud7,cloud8,cloud9,cloud10,cloud11,cloud12 FROM cloudcover WHERE i = ",dbrow,sep="")
              #           #CCMAXX <- dbGetQuery(con,statement=clouds)*100
              #           CCMAXX<- sqlQuery(channel2,clouds)*100
              #           CCMINN <- CCMAXX
              #           CCMAXX1 <-spline(doys12,CCMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic")
              #           CCMAXX <- rep(CCMAXX1$y,nyears)
              #           CCMINN <- CCMAXX
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
                  query<-paste("SELECT a.latitude, a.longitude, b.*
    FROM [AWAPDaily].[dbo].[latlon] as a
  , [AWAPDaily].[dbo].[",years1[i],"] as b
  where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day between ",doystart," and ",doyfinish,")
  order by b.day",sep="")
                }else{
                  if(i==1){ # doing first year, start at day requested
                    query<-paste("SELECT a.latitude, a.longitude, b.*
  FROM [AWAPDaily].[dbo].[latlon] as a
  , [AWAPDaily].[dbo].[",years1[i],"] as b
  where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day >= ",doystart,")
  order by b.day",sep="")
                  }else{
                    if(i==length(years1)){ # doing last year, only go up to last day requested
                      query<-paste("SELECT a.latitude, a.longitude, b.*
  FROM [AWAPDaily].[dbo].[latlon] as a
  , [AWAPDaily].[dbo].[",years1[i],"] as b
  where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day <= ",doyfinish,")
  order by b.day",sep="")
                    }else{ # doing in between years, so get all data for this year
                      query<-paste("SELECT a.latitude, a.longitude, b.*
  FROM [AWAPDaily].[dbo].[latlon] as a
  , [AWAPDaily].[dbo].[",years1[i],"] as b
  where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
  order by b.day",sep="")
                    }}}


                # exectue query and concatenate if necessary
                if(i==1){
                  output1<- sqlQuery(channel,query)
                }else{
                  output1<-rbind(output1,sqlQuery(channel,query))
                }
              } # end loop through years
              query<-paste("SELECT a.*
                       FROM [ausclim].[dbo].[clearskysol] as a
                       where (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
                       ",sep="")
              output_ausclim<- sqlQuery(channel,query)

              weekly_sol<-cbind(1:52,t(output_ausclim[3:54]))
              daily_sol <-suppressWarnings(spline(seq(3,361,7),weekly_sol[,2],n=365,xmin=1,xmax=365,method="periodic"))
              daily_sol<-as.numeric(daily_sol$y)
              dates11<-seq(datestart1,datefinish1,"days")
              output1$dates<-dates11
              for(yr in 1:length(years1)){
                ndays<-nrow(subset(output1,format(output1$dates,"%Y")==as.character(yearstart1-1+yr)))
                clear<-daily_sol
                if(ndays==366){# add day for leap year if needed
                  clear<-c(daily_sol[1:59],daily_sol[59],daily_sol[60:365])
                }
                if(yr==1){
                  all_clear<-clear
                }else{
                  all_clear<-c(all_clear,clear)
                }
              }
              output1$clearsky<-all_clear
              glm_sol<-coefficients(with(output1,glm(sol~rr+tmax+tmin+day+clearsky)))

              dates12<-seq(as.Date(paste(ystart,"-","01-01",sep="")),as.Date(paste(yfinish,"-","12-31",sep="")),"days")
              dates12<-subset(dates12, format(dates12, "%m/%d")!= "02/29") # remove leap years

              results$dates<-dates12
              for(yr in 1:length(nyears)){
                ndays<-nrow(subset(results,format(results$dates,"%Y")==as.character(ystart-1+yr)))
                clear<-daily_sol
                if(ndays==366){# add day for leap year if needed
                  clear<-c(daily_sol[1:59],daily_sol[59],daily_sol[60:365])
                }
                if(yr==1){
                  all_clear<-clear
                }else{
                  all_clear<-c(all_clear,clear)
                }
              }
              output_AWAPDaily$clearsky<-all_clear

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
              TMAXK<-TMAXX+273.15
              loge<-TMAXK
              loge[loge>273.16]<- -7.90298*(373.16/TMAXK[loge>273.16]-1.)+5.02808*log10(373.16/TMAXK[loge>273.16])-1.3816E-07*(10.^(11.344*(1.-TMAXK[loge>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/TMAXK[loge>273.16]-1.))-1.)+log10(1013.246)
              loge[loge<=273.16]<- -9.09718*(273.16/TMAXK[loge<=273.16]-1.)-3.56654*log10(273.16/TMAXK[loge<=273.16])+.876793*(1.-TMAXK[loge<=273.16]/273.16)+log10(6.1071)
              estar<-(10.^loge)*100.
              RHMINN<-(VAPRES/estar)*100
              if(scenario!=""){
                RHMINN<-RHMINN+RH_diff
              }
              RHMINN[RHMINN>100]<-100
              RHMINN[RHMINN<0]<-0.01
              #RHMINN
              TMINK<-TMINN+273.15
              loge<-TMINK
              loge[loge>273.16]<- -7.90298*(373.16/TMINK[loge>273.16]-1.)+5.02808*log10(373.16/TMINK[loge>273.16])-1.3816E-07*(10.^(11.344*(1.-TMINK[loge>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/TMINK[loge>273.16]-1.))-1.)+log10(1013.246)
              loge[loge<=273.16]<- -9.09718*(273.16/TMINK[loge<=273.16]-1.)-3.56654*log10(273.16/TMINK[loge<=273.16])+.876793*(1.-TMINK[loge<=273.16]/273.16)+log10(6.1071)
              estar<-(10.^loge)*100.
              RHMAXX<-(VAPRES/estar)*100
              if(scenario!=""){
                RHMAXX<-RHMAXX+RH_diff
              }
              RHMAXX[RHMAXX>100]<-100
              RHMAXX[RHMAXX<0]<-0.01
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
                    query<-paste("SELECT a.latitude, a.longitude, b.*
    FROM [AWAPDaily].[dbo].[latlon] as a
  , [AWAPDaily].[dbo].[",years1[i],"] as b
  where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day between ",doystart," and ",doyfinish,")
  order by b.day",sep="")
                  }else{
                    if(i==1){ # doing first year, start at day requested
                      query<-paste("SELECT a.latitude, a.longitude, b.*
  FROM [AWAPDaily].[dbo].[latlon] as a
  , [AWAPDaily].[dbo].[",years1[i],"] as b
  where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day >= ",doystart,")
  order by b.day",sep="")
                    }else{
                      if(i==length(years1)){ # doing last year, only go up to last day requested
                        query<-paste("SELECT a.latitude, a.longitude, b.*
  FROM [AWAPDaily].[dbo].[latlon] as a
  , [AWAPDaily].[dbo].[",years1[i],"] as b
  where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,") and (b.day <= ",doyfinish,")
  order by b.day",sep="")
                      }else{ # doing in between years, so get all data for this year
                        query<-paste("SELECT a.latitude, a.longitude, b.*
  FROM [AWAPDaily].[dbo].[latlon] as a
  , [AWAPDaily].[dbo].[",years1[i],"] as b
  where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
  order by b.day",sep="")
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
              TMAXK<-TMAXX+273.15
              loge<-TMAXK
              loge2<-loge
              loge[loge2>273.16]<- -7.90298*(373.16/TMAXK[loge2>273.16]-1.)+5.02808*log10(373.16/TMAXK[loge2>273.16])-1.3816E-07*(10.^(11.344*(1.-TMAXK[loge2>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/TMAXK[loge2>273.16]-1.))-1.)+log10(1013.246)
              loge[loge2<=273.16]<- -9.09718*(273.16/TMAXK[loge2<=273.16]-1.)-3.56654*log10(273.16/TMAXK[loge2<=273.16])+.876793*(1.-TMAXK[loge2<=273.16]/273.16)+log10(6.1071)
              estar<-(10.^loge)*100.
              RHMINN<-(VAPRES/estar)*100
              if(scenario!=""){
                RHMINN<-RHMINN+RH_diff
              }
              RHMINN[RHMINN>100]<-100
              RHMINN[RHMINN<0]<-0.01
              #RHMINN
              TMINK<-TMINN+273.15
              loge<-TMINK
              loge2<-loge
              loge[loge2>273.16]<- -7.90298*(373.16/TMINK[loge2>273.16]-1.)+5.02808*log10(373.16/TMINK[loge2>273.16])-1.3816E-07*(10.^(11.344*(1.-TMINK[loge2>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/TMINK[loge2>273.16]-1.))-1.)+log10(1013.246)
              loge[loge2<=273.16]<- -9.09718*(273.16/TMINK[loge2<=273.16]-1.)-3.56654*log10(273.16/TMINK[loge2<=273.16])+.876793*(1.-TMINK[loge2<=273.16]/273.16)+log10(6.1071)
              estar<-(10.^loge)*100.
              RHMAXX<-(VAPRES/estar)*100
              if(scenario!=""){
                RHMAXX<-RHMAXX+RH_diff
              }
              RHMAXX[RHMAXX>100]<-100
              RHMAXX[RHMAXX<0]<-0.01


            }#end check for year is 1971 or later
          } #end vlsci check
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
          if(vlsci==0){
            WNMAXX <- sqlQuery(channel2,maxwinds)
            WNMINN <- sqlQuery(channel2,minwinds)
          }else{
            if(dailywind!=1){

              load(paste(barcoo,'WNMAXX.bin',sep=''))
              WNMAXX <- as.numeric(as.matrix(data[quadrangle,2:7301]))
              load(paste(barcoo,'WNMINN.bin',sep=''))
              WNMINN<-as.numeric(as.matrix(data[quadrangle,2:7301]))

              WNMAXX<-t(WNMAXX)
              WNMINN<-t(WNMINN)
              WNMAXX<-t(as.data.frame(WNMAXX[((ystart-1989)*365-364):((yfinish-1989)*365)]))
              WNMINN<-t(as.data.frame(WNMINN[((ystart-1989)*365-364):((yfinish-1989)*365)]))
            }
          }
          if(dailywind!=1 ){
            WNMAXX1 <-suppressWarnings(spline(doys12,WNMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
            WNMAXX<-rep(WNMAXX1$y,nyears)
            WNMINN1 <-suppressWarnings(spline(doys12,WNMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
            WNMINN<-rep(WNMINN1$y,nyears)
            if(scenario!=""){
              WNMAXX=WNMAXX*WIND_diff
              WNMINN=WNMINN*WIND_diff
            }
          }
        }
        if(soildata==1){
          SLES1<-suppressWarnings(spline(doys12,SLES,n=timeinterval,xmin=1,xmax=365,method="periodic"))
          SLES<-rep(SLES1$y,nyears)
          SLES<-SLES[1:ndays]
          maxshades1 <-suppressWarnings(spline(doys12,shademax,n=timeinterval,xmin=1,xmax=365,method="periodic"))
          MAXSHADES<-rep(maxshades1$y*100,nyears)
          MAXSHADES<-MAXSHADES[1:ndays]
          if(manualshade==1){
            maxshades <- rep(maxshade,dim)
            MAXSHADES<-maxshades
            minshades <- rep(minshade,dim)
            MINSHADES<-minshades
          }
        }else{
          if(manualshade==0){
            cat("need to sort out leap years here \n")
            maxshades1 <-suppressWarnings(spline(doys12,shademax,n=timeinterval,xmin=1,xmax=365,method="periodic"))
            MAXSHADES<-rep(maxshades1$y*100,nyears)
            minshades <- rep(minshade,dim)
            MINSHADES<-minshades
          }else{
            MAXSHADES<-maxshades
            MINSHADES<-minshades
          }
        }

        REFLS <- rep(REFL, dim)
        if((soildata==1)&(length(RAINFALL)>0)){
          soilwet<-RAINFALL
          soilwet[soilwet<=rainwet] = 0
          soilwet[soilwet>0] = 90
          PCTWET<-pmax(soilwet,PCTWET)
        }else{
          REFLS <- rep(REFL, dim)
          PCTWET <- rep(PCTWET, dim)
          soilwet<-RAINFALL
          soilwet[soilwet<=rainwet] = 0
          soilwet[soilwet>0] = 90
          PCTWET<-pmax(soilwet,PCTWET)
        }

        Numtyps <- 10 # number of substrate types
        Nodes <- matrix(data = 0, nrow = 10, ncol = dim) # deepest nodes for each substrate type
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
          tannulrun<-rep(avetemp,365)
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

        if(opendap == 0){
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
          #Row crops   0.2
          #Low bushes with a few trees 	0.2
          #Heavy trees 	0.25
          #Several buildings 	0.25
          #Hilly, mountainous terrain 	0.25
          if(dailywind!=1){
            WNMINN<-WNMINN*(1.2/10)^0.15*.1 # reduce min wind further because have only 9am/3pm values to get max/min
            WNMAXX<-WNMAXX*(1.2/10)^0.15
            WNMINN<-WNMINN#*3.25 # for snow
            WNMAXX<-WNMAXX#*3.25 # for snow
            message('min wind * 0.1 \n')
            #message('max wind * 2.0 for snow ')
          }else{
            if(snowmodel==0){
              WNMAXX<-dwind*(1.2/2)^0.15
              WNMINN<-WNMAXX
              WNMAXX<-WNMAXX*2#*3.5#*5
              WNMINN<-WNMINN*0.5#1.5#*3.5#*2
              message('min wind * 0.5 \n')
              message('max wind * 2 \n')
            }else{
              WNMAXX<-dwind*(1.2/2)^0.15#*3.25
              WNMINN<-WNMAXX
              WNMAXX<-WNMAXX*2#*2.5#*3.5#*5
              WNMINN<-WNMINN*0.5#*3#1.5#*3.5#*2
              WNMINN[WNMINN<0.1]<-0.1
              #message('min wind * 0.5 * 1.5 for snow \n')
              #message('max wind * 2 * 2.5 for snow \n')
            }
          }
          if(vlsci==0){
            CCMINN<-CCMINN*0.5
            CCMAXX<-CCMAXX*2
            CCMINN[CCMINN>100]<-100
            CCMAXX[CCMAXX>100]<-100
          }
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
          SLES<-matrix(nrow=dim,data=0)
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
          if(timeinterval==365){
            moists2<-matrix(nrow=10, ncol = dim, data=0) # set up an empty vector for soil moisture values through time
          }else{
            moists2<-matrix(nrow=10, ncol = timeinterval, data=0) # set up an empty vector for soil moisture values through time
          }
          moists2[1:10,]<-SoilMoist_Init
          moists<-moists2
        }
        soilprops<-matrix(data = 0, nrow = 10, ncol = 6)

        soilprops[,1]<-BulkDensity
        soilprops[,2]<-min(0.26, 1 - BulkDensity / Density) # not used if soil moisture computed
        soilprops[,3]<-20 # not used
        soilprops[,4]<-Thcond
        soilprops[,5]<-SpecHeat
        soilprops[,6]<-Density
        #         }
        if(cap==1){
          soilprops[1:2,4]<-0.2
          soilprops[1:2,5]<-1920
        }
        if(cap==2){
          soilprops[1:2,4]<-0.1
          soilprops[3:4,4]<-0.25
          soilprops[1:4,5]<-1920
          soilprops[1:4,6]<-1.3
          soilprops[1:4,1]<-0.7
        }

        if(loop>0){
          TMAXX<-c(TMAXX[((loop)*365+1):(nyears*365)],TMAXX[1:((loop)*365)])
          TMINN<-c(TMINN[((loop)*365+1):(nyears*365)],TMINN[1:((loop)*365)])
          RHMAXX<-c(RHMAXX[((loop)*365+1):(nyears*365)],RHMAXX[1:((loop)*365)])
          RHMINN<-c(RHMINN[((loop)*365+1):(nyears*365)],RHMINN[1:((loop)*365)])
          CCMAXX<-c(CCMAXX[((loop)*365+1):(nyears*365)],CCMAXX[1:((loop)*365)])
          CCMINN<-c(CCMINN[((loop)*365+1):(nyears*365)],CCMINN[1:((loop)*365)])
          WNMAXX<-c(WNMAXX[((loop)*365+1):(nyears*365)],WNMAXX[1:((loop)*365)])
          WNMINN<-c(WNMINN[((loop)*365+1):(nyears*365)],WNMINN[1:((loop)*365)])
          PCTWET<-c(PCTWET[((loop)*365+1):(nyears*365)],PCTWET[1:((loop)*365)])
          moists<-cbind(moists[,((loop)*365+1):(nyears*365)],moists[,1:((loop)*365)])
          RAINFALL<-c(RAINFALL[((loop)*365+1):(nyears*365)],RAINFALL[1:((loop)*365)])

        }

        # microclimate input parameters listALTT,ALREF,ALMINT,ALONG,AMINUT,ALAT
        ALTT<-as.numeric(ALTT)
        ALREF<-as.numeric(ALREF)
        ALMINT<-as.numeric(ALMINT)
        ALONG<-as.numeric(ALONG)
        AMINUT<-as.numeric(AMINUT)
        ALAT<-as.numeric(ALAT)

        # microclimate input parameters list
        microinput<-c(dim,RUF,ERR,Usrhyt,Refhyt,Numtyps,Z01,Z02,ZH1,ZH2,idayst,ida,HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,CMH2O,microdaily,tannul,EC,VIEWF,snowtemp,snowdens,snowmelt,undercatch,rainmult,runshade,runmoist,maxpool,evenrain,snowmodel,rainmelt,writecsv,densfun,hourly,rainhourly,lamb,IUV,RW,PC,RL,SP,R1,IM,MAXCOUNT,IR,message,fail)

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
        }
        if(shore==0){
          tides<-matrix(data = 0., nrow = 24*dim, ncol = 3) # make an empty matrix
        }
        # all microclimate data input list - all these variables are expected by the input argument of the fortran micro2014 subroutine
        micro<-list(tides=tides,microinput=microinput,doy=doy,SLES=SLES1,DEP=DEP,Nodes=Nodes,MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TIMAXS=TIMAXS,TIMINS=TIMINS,TMAXX=TMAXX1,TMINN=TMINN1,RHMAXX=RHMAXX1,RHMINN=RHMINN1,CCMAXX=CCMAXX1,CCMINN=CCMINN1,WNMAXX=WNMAXX1,WNMINN=WNMINN1,TAIRhr=TAIRhr,RHhr=RHhr,WNhr=WNhr,CLDhr=CLDhr,SOLRhr=SOLRhr,RAINhr=RAINhr,ZENhr=ZENhr,REFLS=REFLS1,PCTWET=PCTWET1,soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists1,RAINFALL=RAINFALL1,tannulrun=tannulrun,PE=PE,KS=KS,BB=BB,BD=BD,DD=DD,L=L,LAI=LAI)
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
    } # end of check if soil data is being used but no soil data returned
  } # end error trapping
} # end of micro_aust function
