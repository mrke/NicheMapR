#' NCEP implementation of the microclimate model, with package microclima computing hourly forcing.
#'
#' An implementation of the NicheMapR microclimate model that uses the RNCEP package daily weather database https://sites.google.com/site/michaelukemp/rncep, and specifically uses the following variables: air.2m, prate.sfc, shum.2m, pres.sfc, tcdc.eatm, uwnd.10m, vwnd.10m, dswrf.sfc. At the moment uses the same DEM from the CRU global climate data set.
#' @encoding UTF-8
#' @param loc Either a longitude and latitude (decimal degrees) or a place name to search for on Google Earth
#' @param dstart First day to run, date in format "d-m-Y" e.g. "01-01-2016"
#' @param dfinish Last day to run, date in format "d-m-Y" e.g. "31-12-2016"
#' @param dem a digital elevation model produce by microclima function 'get_dem' via R package 'elevatr' (internally generated via same function based on 'loc' if NA)
#' @param REFL Soil solar reflectance, decimal \%
#' @param slope Slope in degrees (if negative, then derived from DEM with package microclima)
#' @param aspect Aspect in degrees (0 = north) (if negative, then derived from DEM with microclima)
#' @param DEP Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
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
#' @usage micro_clima(loc = 'Galapagos', dstart = "01-01-2017", dfinish = "31-12-2017",
#' REFL = 0.15, slope = 0, aspect = 0, DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200),
#' Usrhyt = 0.01, ...)
#' @export
#' @details
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#' \code{run.gads}{ = 1, Use the Global Aerosol Database? 1=yes, 0=no}\cr\cr
#' \code{lamb}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{IUV}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{writecsv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{reanalysis}{ = TRUE, Use reanalysis2 NCEP data? TRUE/FALSE}\cr\cr
#' \code{windfac}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{warm}{ = 0, uniform warming, °C}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#' \code{spatial}{ = NA, specify folder with local NCEP data (no trailing forward slash), goes to the web if NA }\cr\cr
#' \code{save}{ = 0, don't save forcing data (0), save the forcing data (1) or read previously saved data (2)}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{ERR}{ = 1.5, Integrator error tolerance for soil temperature calculations}\cr\cr
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
#' \code{LAI}{ = 0.1, leaf area index, used to partition traspiration/evaporation from PET in soil moisture model and for package microclima radiation calcs}\cr\cr
#' \code{LOR}{ = 1, leaf orientation for package microclima radiation calcs}\cr\cr
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
#' \code{dim}{ - number of days for which predictions are made}\cr\cr
#' \code{longlat}{ - longitude and latitude for which simulation was run (decimal degrees)}\cr\cr
#' \code{dates}{ - vector of dates (POSIXct, UTC)}\cr\cr
#' \code{nyears}{ - number of years for which predictions are made}\cr\cr
#' \code{RAINFALL}{ - vector of daily rainfall (mm)}\cr\cr
#' \code{elev}{ - elevation at point of simulation (m)}\cr\cr
#' \code{minshade}{ - minimum shade for simulation (\%)}\cr\cr
#' \code{maxshade}{ - maximum shade for simulation (single value - if time varying, in 'MAXSHADES') (\%)}\cr\cr
#' \code{MAXSHADES}{ - vector of maximum shades used (\%)}\cr\cr
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
#' \item  4 LEAFPOT - leaf water potentail (J/kg)
#' \item  5-14 RPOT0cm ... - root water potentail (J/kg), at each of the 10 specified depths
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
#' dstart <- "01/01/1976"
#' dfinish <- "31/12/1978"
#' micro<-micro_clima(loc = 'Galapagos', soilgrids = 0, dstart = dstart, dfinish = dfinish, reanalysis = FALSE)
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
micro_ncep <- function(
  loc = c(-5.3, 50.13),
  dstart = "01/01/2017",
  dfinish = "31/12/2017",
  dem = NA,
  nyears=as.numeric(substr(dfinish, 7, 10)) - as.numeric(substr(dstart, 7, 10)) + 1,
  REFL=0.15,
  slope=0,
  aspect=0,
  DEP=c(0., 2.5,  5.,  10.,  15.,  20.,  30.,  50.,  100.,  200.),
  Refhyt=1,
  Usrhyt=.01,
  Z01=0,
  Z02=0,
  ZH1=0,
  ZH2=0,
  run.gads=1,
  write_input=0,
  writecsv=0,
  reanalysis=TRUE,
  windfac = 1,
  warm=0,
  ERR=1.5,
  RUF=0.004,
  EC=0.0167238,
  SLE=0.95,
  Thcond=2.5,
  Density=2.56,
  SpecHeat=870,
  BulkDensity=1.3,
  PCTWET=0,
  rainwet=1.5,
  cap=1,
  CMH2O=1,
  hori=rep(0,24),
  runmoist=1,
  PE=rep(1.1,19),
  KS=rep(0.0037,19),
  BB=rep(4.5,19),
  BD=rep(1.3,19),
  DD=rep(2.56,19),
  maxpool=10000,
  rainmult=1,
  evenrain=0,
  SoilMoist_Init=c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3),
  L = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000,
  R1 = 0.001,
  RW = 2.5e+10,
  RL = 2e+6,
  PC = -1500,
  SP = 10,
  IM = 1e-06,
  MAXCOUNT = 500,
  LAI=0.1,
  LOR=1,
  snowmodel=1,
  snowtemp=1.5,
  snowdens=0.375,
  densfun=c(0.5979, 0.2178, 0.001, 0.0038),
  snowmelt=1,
  undercatch=1,
  rainmelt=0.0125,
  shore=0,
  tides = 0,
  hourly=1,
  rainhourly = 0,
  rainhour = 0,
  rainoff=0,
  lamb = 0,
  IUV = 0,
  soilgrids = 0,
  message = 0,
  fail = nyears * 24 * 365,
  spatial = NA,
  save = 0,
  snowcond = 0,
  intercept = 0 / 100 * 0.3,
  grasshade = 0){ # end function parameters

  #
  #   loc = c(-5.3, 50.13)
  #   dstart = "01/01/2010"
  #   dfinish = "31/12/2010"
  #   dem = NA
  #   nyears=as.numeric(substr(dfinish, 7, 10)) - as.numeric(substr(dstart, 7, 10)) + 1
  #   REFL=0.15
  #   slope=0
  #   aspect=0
  #   DEP=c(0., 2.5,  5.,  10.,  15.,  20.,  30.,  50.,  100.,  200.)
  #   Refhyt=2
  #   Usrhyt=.01
  #   Z01=0
  #   Z02=0
  #   ZH1=0
  #   ZH2=0
  #   run.gads=1
  #   write_input=0
  #   writecsv=0
  #   reanalysis=TRUE
  #   windfac = 1
  #   warm=0
  #   ERR=1.5
  #   RUF=0.004
  #   EC=0.0167238
  #   SLE=0.95
  #   Thcond=2.5
  #   Density=2.56
  #   SpecHeat=870
  #   BulkDensity=1.3
  #   PCTWET=0
  #   rainwet=1.5
  #   cap=1
  #   CMH2O=1
  #   hori=rep(0,24)
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
  #   L = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000
  #   R1 = 0.001
  #   RW = 2.5e+10
  #   RL = 2e+6
  #   PC = -1500
  #   SP = 10
  #   IM = 1e-06
  #   MAXCOUNT = 500
  #   LAI=0.1
  #   LOR=1
  #   snowmodel=1
  #   snowtemp=1.5
  #   snowdens=0.375
  #   densfun=c(0.5979, 0.2178, 0.001, 0.0038)
  #   snowmelt=1
  #   undercatch=1
  #   rainmelt=0.0125
  #   shore=0
  #   tides = 0
  #   hourly=1
  #   rainhourly = 0
  #   rainhour = 0
  #   rainoff=0
  #   lamb = 0
  #   IUV = 0
  #   soilgrids = 0
  #   message = 0
  #   fail = nyears * 24 * 365
  #   spatial = "c:/Spatial_Data/ncep/"
  #   save = 0
  #   snowcond = 0
  #   intercept = 0 / 100 * 0.3
  #   grasshade = 0
  library(microclima)
  library(elevatr)
  library(RNCEP)
  library(zoo)
  library(raster)
  library(proj4)
  library(ncdf4)
  # Inputs:
  # (1) r - a raster object with projection system specified
  # (2) lat and long in decimal degrees (Only used if r is NA)
  # (3) resolution (resolution required)
  #  (4) returned dataset includes bathymetry. zmin sets values below this
  # in returned in zmin
  # https://mapzen.com/documentation/terrain-tiles/data-sources/
  get_dem <- function(r = NA, lat, long, resolution = 30, zmin = 0) {
    if (resolution < 30) {
      warning("Higher resolution data only available for some locations. DEM
            may be derived by interpolation")
    }
    if (class(r) != "RasterLayer") {
      xy <- data.frame(x = long, y = lat)
      coordinates(xy) = ~x + y
      proj4string(xy) = "+init=epsg:4326"
      if (lat >= -80 & lat <= 84)
        xy <- as.data.frame(spTransform(xy, CRS("+init=epsg:3395")))
      if (lat > 84)
        xy <- as.data.frame(spTransform(xy, CRS("+init=epsg:3413")))
      if (lat < -80)
        xy <- as.data.frame(spTransform(xy, CRS("+init=epsg:3976")))
      e <- extent(c(xy$x - 100 * resolution, xy$x + 100 * resolution,
        xy$y - 100 * resolution, xy$y + 100 * resolution))
      r <- raster(e)
      res(r) <- resolution
      if (lat >= -80 & lat <= 84)
        crs(r) <- "+init=epsg:3395"
      if (lat > 84)
        crs(r) <- "+init=epsg:3413"
      if (lat < -80)
        crs(r) <- "+init=epsg:3976"
    } else {
      lat <- latlongfromraster(r)$lat
      long <- latlongfromraster(r)$long
      res(r) <- resolution
    }
    z = ceiling(log((cos(lat * pi/180) * 2 * pi * 6378137) / (256 * resolution), 2))
    z <- ifelse(z > 14, 14, z)
    r2 <-get_elev_raster(r, z = z, src = "aws")
    r2<- resample(r2, r)
    m2 <- getValues(r2, format = "matrix")
    m2[m2 < zmin] <- zmin
    m2[is.na(m2)] <- zmin
    r2 <- if_raster(m2, r2)
    return(r2)
  }

  # get nearest land pixel from ncep
  landlatlong <- function(lat, long, msk) {
    xy <- data.frame(xyFromCell(msk, 1:ncell(msk)))
    long <- ifelse(long < 0, long  + 360, long)
    xydif <- data.frame(x = xy$x - long, y = xy$y - lat)
    xydif$dist <- sqrt(xydif$x^2 + xydif$y^2)
    sel <- which(xydif$dist == min(xydif$dist))
    return(xy[sel[1],])
  }
  get_inputs <- function(lat, long, tme, reanalysis2 = TRUE, spatial = NA) {
    sorttimes <- function(tme) {
      tme2 <- as.POSIXct(tme)
      tme2 <- c((tme2 - 24 * 3600), tme2, (tme2 + 24 * 3600), (tme2 + 48 * 3600))
      tme2 <- as.POSIXlt(tme2)
      yrs <- unique(tme2$year + 1900)
      mths <- unique(tme2$mon + 1)
      yrs <-yrs[order(yrs)]
      mths <-mths[order(mths)]
      tma <- 0
      dm <- c(31,28,31,30,31,30,31,31,30,31,30,31) * 4 - 1
      for (yr in 1:length(yrs)) {
        dm[2] <- ifelse(yrs[yr]%%4 == 0, 115, 111)
        for (mth in 1:length(mths)) {
          tmym <- as.POSIXct(c(0:dm[mth]) * 3600 * 6,
            origin = paste0(yrs[yr],"-",mths[mth],"-01 00:00"),
            tz = 'GMT')
          tma <- c(tma, tmym)
        }
      }
      tma <- as.POSIXlt(tma[-1], origin = '1970-01-01', tz = 'GMT')
      sel <- which(tma >= min(tme2) & tma <= max(tme2))
      sel <- sel[-length(sel)]
      return(list(tme = tma, sel = sel))
    }
    ncepget1 <- function(climvar, tme2, ll, sel) {
      yrs <- unique(tme2$year + 1900)
      mths <- unique(tme2$mon + 1)
      v <- NCEP.gather(climvar, level = 'gaussian',
        years.minmax = c(min(yrs),max(yrs)),
        months.minma = c(min(mths):max(mths)),
        lat.southnorth = c(ll$y,ll$y), lon.westeast = c(ll$x,ll$x),
        reanalysis2 = reanalysis2, return.units = FALSE, status.bar = FALSE)
      v <- apply(v, 3, mean, na.rm = T)
      v[sel]
    }
    ncquery <- function(filename, var, start, count, year){
      nc <- nc_open(paste(spatial, "/",filename,year,".nc", sep = ""))
      out <- as.numeric(ncvar_get(nc, varid = var, start = start, count = count))
      nc_close(nc)
      out
    }
    ll <- data.frame(x = long, y = lat)
    tme2 <- sorttimes(tme)$tme
    sel <- sorttimes(tme)$sel
    if(is.na(spatial) == TRUE){
      # These variables are forecasts valid 6 hours after the reference time.
      Tk <- ncepget1('air.2m', tme2, ll, sel)
      Tkmin <- ncepget1('tmin.2m', tme2, ll, sel)
      Tkmax <- ncepget1('tmax.2m', tme2, ll, sel)
      sh <- ncepget1('shum.2m', tme2, ll, sel)
      pr <- ncepget1('pres.sfc', tme2, ll, sel)
      wu <- ncepget1('uwnd.10m', tme2, ll, sel)
      wv <- ncepget1('vwnd.10m', tme2, ll, sel)
      # These variables are 6 hour averages starting at the reference time.
      dlw <- ncepget1('dlwrf.sfc', tme2, ll, sel)
      ulw <- ncepget1('ulwrf.sfc', tme2, ll, sel)
      dsw <- ncepget1('dswrf.sfc', tme2, ll, sel)
      tcdc <- ncepget1('tcdc.eatm', tme2, ll, sel)
      prate <- ncepget1('prate.sfc', tme2, ll, sel)
    }else{
      cat(paste0("extracting weather data locally from ", spatial, " \n"))
      years <- as.numeric(unique(format(tme, "%Y")))
      nyears <- length(years)
      nc <- nc_open(paste(spatial, "/tmax.2m.gauss.", years[1], ".nc",
        sep = ""))
      lon2 <- matrix(ncvar_get(nc, "lon"))
      #lon2[lon2 > 180] <- - 180 + (lon2[lon2 > 180] - 180)
      lat2 <- matrix(ncvar_get(nc, "lat"))
      lon_1 <- long
      if(lon_1 < 0){lon_1 <- 180 - (long*-1) + 180}
      lat_1 <- lat
      dist1 <- abs(lon2 - lon_1)
      index1 <- which.min(dist1)
      dist2 <- abs(lat2 - lat_1)
      index2 <- which.min(dist2)
      start <- c(index1, index2, 1, 1)
      count <- c(1, 1, 1, -1)
      start2 <- c(index1, index2, 1, 1460-4)
      count2 <- c(1, 1, 1, 4)
      nc_close(nc)

      for (j in 1:(nyears+2)) {
        if (j == 1) {
          Tkmin <- ncquery("tmin.2m.gauss.", "tmin", start2, count2, years[j]-1)
          Tkmax <- ncquery("tmax.2m.gauss.", "tmax", start2, count2, years[j]-1)
          Tk <- ncquery("air.2m.gauss.", "air", start2, count2, years[j]-1)
          sh <- ncquery("shum.2m.gauss.", "shum", start2, count2, years[j]-1)
          pr <- ncquery("pres.sfc.gauss.", "pres", start2[c(1,2,4)], count2[c(1,2,4)], years[j]-1)
          tcdc <- ncquery("tcdc.eatm.gauss.", "tcdc", start2[c(1,2,4)], count2[c(1,2,4)], years[j]-1)
          dsw <- ncquery("dswrf.sfc.gauss.", "dswrf", start2[c(1,2,4)], count2[c(1,2,4)], years[j]-1)
          dlw <- ncquery("dlwrf.sfc.gauss.", "dlwrf", start2[c(1,2,4)], count2[c(1,2,4)], years[j]-1)
          ulw <- ncquery("ulwrf.sfc.gauss.", "ulwrf", start2[c(1,2,4)], count2[c(1,2,4)], years[j]-1)
          wu <- ncquery("uwnd.10m.gauss.", "uwnd", start2, count2, years[j]-1)
          wv <- ncquery("vwnd.10m.gauss.", "vwnd", start2, count2, years[j]-1)
          prate <- ncquery("prate.sfc.gauss.", "prate", start2[c(1,2,4)], count2[c(1,2,4)], years[j]-1)
        }else{
          if(j <= nyears+1){
            cat(paste("reading weather input for ", years[j-1]," \n", sep = ""))
            Tkmin <- c(Tkmin, ncquery("tmin.2m.gauss.", "tmin", start, count, years[j-1]))
            Tkmax <- c(Tkmax, ncquery("tmax.2m.gauss.", "tmax", start, count, years[j-1]))
            Tk <- c(Tk, ncquery("air.2m.gauss.", "air", start, count, years[j-1]))
            sh <- c(sh, ncquery("shum.2m.gauss.", "shum", start, count, years[j-1]))
            pr <- c(pr, ncquery("pres.sfc.gauss.", "pres", start[c(1,2,4)], count[c(1,2,4)], years[j-1]))
            tcdc <- c(tcdc, ncquery("tcdc.eatm.gauss.", "tcdc", start[c(1,2,4)], count[c(1,2,4)], years[j-1]))
            dsw <- c(dsw, ncquery("dswrf.sfc.gauss.", "dswrf", start[c(1,2,4)], count[c(1,2,4)], years[j-1]))
            dlw <- c(dlw, ncquery("dlwrf.sfc.gauss.", "dlwrf", start[c(1,2,4)], count[c(1,2,4)], years[j-1]))
            ulw <- c(ulw, ncquery("ulwrf.sfc.gauss.", "ulwrf", start[c(1,2,4)], count[c(1,2,4)], years[j-1]))
            wu <- c(wu, ncquery("uwnd.10m.gauss.", "uwnd", start, count, years[j-1]))
            wv <- c(wv, ncquery("vwnd.10m.gauss.", "vwnd", start, count, years[j-1]))
            prate <- c(prate, ncquery("prate.sfc.gauss.", "prate", start[c(1,2,4)], count[c(1,2,4)], years[j-1]))
          } else {
            Tkmin <- c(Tkmin, ncquery("tmin.2m.gauss.", "tmin", start, count2, years[j-2]+1))
            Tkmax <- c(Tkmax, ncquery("tmax.2m.gauss.", "tmax", start, count2, years[j-2]+1))
            Tk <- c(Tk, ncquery("air.2m.gauss.", "air", start, count2, years[j-2]+1))
            sh <- c(sh, ncquery("shum.2m.gauss.", "shum", start, count2, years[j-2]+1))
            pr <- c(pr, ncquery("pres.sfc.gauss.", "pres", start[c(1,2,4)], count2[c(1,2,4)], years[j-2]+1))
            tcdc <- c(tcdc, ncquery("tcdc.eatm.gauss.", "tcdc", start[c(1,2,4)], count2[c(1,2,4)], years[j-2]+1))
            dsw <- c(dsw, ncquery("dswrf.sfc.gauss.", "dswrf", start[c(1,2,4)], count2[c(1,2,4)], years[j-2]+1))
            dlw <- c(dlw, ncquery("dlwrf.sfc.gauss.", "dlwrf", start[c(1,2,4)], count2[c(1,2,4)], years[j-2]+1))
            ulw <- c(ulw, ncquery("ulwrf.sfc.gauss.", "ulwrf", start[c(1,2,4)], count2[c(1,2,4)], years[j-2]+1))
            wu <- c(wu, ncquery("uwnd.10m.gauss.", "uwnd", start, count2, years[j-2]+1))
            wv <- c(wv, ncquery("vwnd.10m.gauss.", "vwnd", start, count2, years[j-2]+1))
            prate <- c(prate, ncquery("prate.sfc.gauss.", "prate", start[c(1,2,4)], count2[c(1,2,4)], years[j-2]+1))
          }
        }
      }
      dsw[dsw<0] <- 0
      prate[prate<0] <- 0
    }
    dfout <- data.frame(obs_time = tme2[sel], Tk, Tkmin, Tkmax, sh, pr, wu, wv, dlw, ulw, dsw, tcdc, prate)
    rownames(dfout) <- NULL
    return(dfout)
  }
  # needs a time offset
  get_rain <- function(lat, long, tme) {
    yrs <- unique(tme$year + 1900)
    mths <- unique(tme$mon + 1)
    ll <- data.frame(x = long, y = lat)
    #ll <- landlatlong(lat, long, msk)
    precip <- NCEP.gather('prate.sfc', level = 'gaussian',
      years.minmax = c(min(yrs),max(yrs)),
      months.minma = c(min(mths):max(mths)),
      lat.southnorth = c(ll$y,ll$y), lon.westeast = c(ll$x,ll$x),
      status.bar = FALSE, reanalysis2 = reanalysis)
    precip <- apply(precip, 3, mean, na.rm = T)
    precip <- precip * 6 * 3600
    mpre <- t(matrix(as.numeric(precip), nrow = 4))
    dpre <- apply(mpre, 1, sum)
    dpre
  }

  siflat <- function(localtime, lat, long, julian){
    saltitude <- solalt(localtime, lat, long, julian)
    alt <- saltitude * (pi/180)
    index <- cos(pi/2 - alt)
    index[index < 0] <- 0
    index
  }
  difprop2 <- function (rad, julian, localtime, lat, long, hourly = FALSE,
    watts = TRUE, merid = 0, dst = 0)
  {
    if (watts)
      rad <- rad * 0.0036
    lt <- c(localtime:(localtime + 5)) + 0.5
    jd <- rep(julian, 6)
    sel <- which(lt > 24)
    lt[sel] <- lt[sel] -24
    jd[sel] <- jd[sel] + 1
    sa <- solalt(lt, lat, long, jd, merid, dst)
    alt <- sa * (pi/180)
    k1 <- 0.83 - 0.56 * exp(-0.06 * sa)
    si <- cos(pi/2 - alt)
    k <- rad/(4.87 * si)
    k[is.na(k)] <- 0
    k <- ifelse(k > k1, k1, k)
    k[k < 0] <- 0
    rho <- k/k1
    sigma3a <- 0.021 + 0.397 * rho - 0.231 * rho^2 - 0.13 *
      exp(-1 * (((rho - 0.931)/0.134)^2)^0.834)
    sigma3b <- 0.12 + 0.65 * (rho - 1.04)
    sigma3 <- ifelse(rho <= 1.04, sigma3a, sigma3b)

    k2 <- 0.95 * k1
    d1 <- ifelse(sa > 1.4, 0.07 + 0.046 * (90 - sa)/(sa + 3),
      1)
    d1 <- mean(d1)
    K <- 0.5 * (1 + sin(pi * (k - 0.22)/(k1 - 0.22) - pi/2))
    d2 <- 1 - ((1 - d1) * (0.11 * sqrt(K) + 0.15 * K + 0.74 *
        K^2))
    d3 <- (d2 * k2) * (1 - k)/(k * (1 - k2))
    alpha <- (1/sin(alt))^0.6
    kbmax <- 0.81^alpha
    kmax <- (kbmax + d2 * k2/(1 - k2))/(1 + d2 * k2/(1 - k2))
    dmax <- (d2 * k2) * (1 - kmax)/(kmax * (1 - k2))
    d4 <- 1 - kmax * (1 - dmax)/k
    d <- ifelse(k <= kmax, d3, d4)
    d <- ifelse(k <= k2, d2, d)
    d <- ifelse(k <= 0.22, 1, d)
    kX <- 0.56 - 0.32 * exp(-0.06 * sa)
    kL <- (k - 0.14)/(kX - 0.14)
    kR <- (k - kX)/0.71
    delta <- ifelse(k >= 0.14 & k < kX, -3 * kL^2 * (1 - kL) *
        sigma3^1.3, 0)
    delta <- ifelse(k >= kX & k < (kX + 0.71), 3 * kR * (1 -
        kR)^2 * sigma3^0.6, delta)
    d[sigma3 > 0.01] <- d[sigma3 > 0.01] + delta[sigma3 > 0.01]
    d[si <= 0] <- NA
    d[d > 1] <- 1
    mean(d, na.rm = T)
  }

  hourlytemp2 <- function (julian, dni, dif, em, mintemp, maxtemp,
    lat, long, merid = 0, tz = 0, dst = 0)
  {
    solarday <- function(julian, localtime, lat, long, tz = 0,
      dst = 0) {
      src <- suntimes(julian, lat, long, tz, dst)$sunrise
      ssc <- suntimes(julian, lat, long, tz, dst)$sunset
      ssp <- suntimes(julian - 1, lat, long, tz, dst)$sunset
      srn <- suntimes(julian + 1, lat, long, tz, dst)$sunrise
      st <- ifelse(localtime >= src & localtime <= ssc, ((localtime -
          src)/(ssc - src)) * 12 + 6, localtime)
      st <- ifelse(localtime > ssc, ((localtime - ssc)/(srn +
          24 - ssc)) * 12 + 18, st)
      st <- ifelse(localtime < src, ((localtime + 24 - ssp)/(src +
          24 - ssp)) * 12 - 6, st)
      st <- (st/24)%%1
      st
    }
    propmaxrad <- function(dif, dct, am) {
      theta <- 1.1 * (0.7^(am^0.678))
      mxr <- 4.87 * theta
      pr <- (dif + dct)/mxr
      pr <- ifelse(pr > 1, 1, pr)
      pr
    }
    l1 <- unlist(lapply(list(em, dni, dif), length))
    l2 <- unlist(lapply(list(julian, mintemp, maxtemp), length))
    if (length(unique(l1)) != 1) {
      warning("Number of hourly values not identical")
    }
    if (length(unique(l2)) != 1) {
      warning("Number of daily values not identical")
    }
    if (length(em)%%24 != 0)
      warning("Hourly values must be multiple of 24")
    if (l1[1] != l2[1] * 24) {
      warning("Number of hourly values not 24 times number of daily values")
    }
    o <- order(rep(c(1:length(julian)), 24))
    jd <- rep(julian, 24)[o]
    localtime <- rep(c(0:23), length(mintemp))
    solfrac <- solarday(julian, localtime, lat, long, tz, dst)
    am <- airmasscoef(localtime, lat, long, julian, merid, dst)
    dct <- 0
    for (i in 1:length(jd)) {
      dct[i] <- dni[i] * solarindex(0, 0, localtime[i], lat,
        long, jd[i], shadow = FALSE)
    }
    pr <- propmaxrad(dif, dct, am)
    tfrac <- 110.42201 * sin(solfrac) - 38.64802 * cos(solfrac) -
      79.82963 * sin(2 * solfrac) + 65.39122 * cos(2 * solfrac) +
      15.54387 * sin(3 * solfrac) - 26.30047 * cos(3 * solfrac)
    mns <- rep(mintemp, 24)[o]
    mxs <- rep(maxtemp, 24)[o]
    tc <- (mxs - mns) * tfrac + mns
    day <- which(is.na(pr) == F)
    ngt <- which(is.na(pr))
    tfrac[day] <- -0.5154639 + 1.1060481 * tfrac[day] + 6.1147121 *
      pr[day] - 7.5826432 * tfrac[day] * pr[day]
    tfrac[ngt] <- -4.614759 + 5.909548 * tfrac[ngt] + 4.308512 *
      em[ngt] - 5.020017 * tfrac[ngt] * em[ngt]
    tfrac <- 1/(1 + exp(-1 * tfrac))
    td <- array(tfrac, dim = c(24, length(tfrac)/24))
    tmns <- apply(td, 2, min)
    tmxs <- apply(td, 2, max)
    for (i in 1:dim(td)[2]) {
      td[, i] <- (td[, i] - tmns[i])/(tmxs[i] - tmns[i])
    }
    tfrac <- as.vector(td)
    tc <- (mxs - mns) * tfrac + mns
    tc
  }
  hourlyNCEP <- function(tme, lat, long, reanalysis2 = TRUE, spatial = NA) {
    int <- as.numeric(tme[2]) - as.numeric(tme[1])
    lgth <- (length(tme) * int) / (24 * 3600)
    tme2 <- as.POSIXlt(c(0:(lgth - 1)) * 3600 * 24, origin = min(tme), tz = 'GMT')
    ncepdata <- get_inputs(lat, long, tme2, reanalysis2, spatial)
    prate <- ncepdata$prate
    tme6 <- as.POSIXlt(ncepdata$obs_time)
    n <- (length(tme6) - 1) * 6 + 1
    h_pr <- spline(tme6, ncepdata$pr, n = n)$y
    h_sh <- spline(tme6, ncepdata$sh, n = n)$y
    h_tcdc <- spline(tme6, ncepdata$tcdc, n = n)$y
    tmo  <- spline(tme6, ncepdata$sh, n = n)$x
    tmo <- as.POSIXlt(tmo, origin = "1970-01-01 00:00",
      tz = "GMT")
    # tmoend <- as.POSIXlt(c(1:6) * 3600, origin = max(tmo),
    #   tz = "GMT")
    # tmo2 <- c(tmo, tmoend)
    tmo2 <- as.POSIXlt(seq(0,(lgth + 2), by = 1/24) * 3600 * 24, origin = min(tme)-3600*24, tz = 'GMT')
    jd <- julday(tmo2$year + 1900,  tmo2$mon + 1, tmo2$mday)
    si <- siflat(tmo2$hour, lat, long, jd)
    amc <- airmasscoef(tmo2$hour, lat, long, jd)
    sa2<-solalt(tmo2$hour, lat, long, jd, 0, 0)
    z <- 90 - sa2
    amc <- 1/(cos(z * pi/180) + 0.50572 * (96.07995 - z)^(-1.6364))
    amc[z >= 88] <- 0
    si_m <- 0
    am_m <- 0
    for (i in 1:(length(tme6))) {
      st <- (i - 1) * 6 + 1
      ed <- st + 5
      si_m[i] <- mean(si[st:ed], na.rm = T)
      am_m[i] <- mean(amc[st:ed], na.rm = T)
    }
    af <- 1.1 * 0.7 ^ (am_m^0.678)
    jd <- julday(tme6$year + 1900, tme6$mon + 1, tme6$mday)
    dp <- 0
    for (i in 1:length(jd)) {
      dp[i] <- difprop2(ncepdata$dsw[i], jd[i], tme6$hour[i], lat, long)
    }
    dp[ncepdata$dsw == 0] <- NA
    dirr <- (ncepdata$dsw * (1 - dp)) / si_m
    difr <- (ncepdata$dsw * dp) / af
    dirr[dirr > 1367] <- 0 # ensuring low solar angles don't cause values > solar constant
    difr[difr > 1367] <- 0 # ensuring low solar angles don't cause values > solar constant
    globr <- dirr + difr
    #globr[globr > (4.87 / 0.0036)] <- (4.87 / 0.0036)
    nd <- length(globr)
    globr[1] <- mean(globr[1:4], na.rm = T)
    globr[nd] <- mean(globr[(nd - 4):nd], na.rm = T)
    dp[1] <- mean(dp[1:4], na.rm = T)
    dp[nd] <- mean(dp[(nd - 4):nd], na.rm = T)
    globr <- na.approx(globr, na.rm = F)
    dp <- na.approx(dp, na.rm = F)
    h_dp <- spline(tme6, dp, n = n)$y
    h_gr <- spline(tme6, globr, n = n)$y
    tmorad <- as.POSIXlt(tmo + 3600 * 3)
    jd <- julday(tmorad$year + 1900, tmorad$mon + 1, tmorad$mday)
    si <- siflat(tmorad$hour, lat, long, jd)
    szenith <- 90 - solalt(tmorad$hour, lat, long, jd)
    am <- airmasscoef(tmorad$hour, lat, long, jd)
    af <- 1.1 * 0.7 ^ (am^0.678)
    h_dni <- h_gr * (1 - h_dp)
    h_dni[si == 0] <- 0
    h_dif <- h_gr * h_dp * af
    h_dif[is.na(h_dif)] <- 0
    em <- ncepdata$dlw / ncepdata$ulw
    em_h <- spline(tme6, em, n = n)$y
    t6sel <- c(4:(length(tme6)-5))
    thsel <-c(22:(length(tmorad)-22))
    thsel2 <-c(19:(length(tmo) - 25))
    tcmn <- ncepdata$Tkmin[t6sel] - 273.15
    tcmx <- ncepdata$Tkmax[t6sel] - 273.15
    tcmn <-t(matrix(tcmn, nrow = 4))
    tcmx <-t(matrix(tcmx, nrow = 4))
    tmin <- apply(tcmn, 1, min)
    tmax <- apply(tcmx, 1, max)
    jd <- julday(tme2$year + 1900, tme2$mon + 1, tme2$mday)
    h_tc<-hourlytemp2(julian = jd, em = em_h[thsel], dni = h_dni[thsel],
      dif = h_dif[thsel], mintemp = tmin, maxtemp = tmax,
      lat = lat, long = long)
    hlwu <- 2.043e-10 * (h_tc + 273.15)^4
    hlwd <- em_h[thsel] *  hlwu
    h_nlw <- hlwu - hlwd
    h_uw <- spline(tme6, ncepdata$wu, n = n)$y
    h_vw <- spline(tme6, ncepdata$wv, n = n)$y
    h_ws <- sqrt(h_uw^2 + h_vw^2)
    h_ws <- windheight(h_ws, 10, 1)
    h_wd <- atan2(h_uw, h_vw) * 180/pi + 180
    h_wd <- h_wd%%360
    hourlyout <- data.frame(obs_time = tmorad[thsel], temperature = h_tc, humidity = h_sh[thsel2],
      pressure = h_pr[thsel2], windspeed = h_ws[thsel2], winddir = h_wd[thsel2],
      emissivity = em_h[thsel2], netlong = h_nlw, uplong = hlwu, downlong = hlwd,
      rad_dni = h_dni[thsel] * 0.0036, rad_dif = h_dif[thsel] * 0.0036,
      szenith = szenith[thsel], h_tcdc = h_tcdc[thsel])
    return(list(hourlyout = hourlyout, prate = prate))
  }
  # allow for horizon angle and slope/aspect override, as well as horizon angle, slope and aspect output, substituted 'lor' for 'x'
  pointradwind <- function(hourlydata, dem, lat, long, l, x, albr = 0.15, zmin = 0, slope = NA, aspect = NA) {
    m <- is_raster(dem)
    m[is.na(m)] <- zmin
    m[m < zmin] <- zmin
    dem <- if_raster(m, dem)
    xy <- data.frame(x = long, y = lat)
    coordinates(xy) = ~x + y
    proj4string(xy) = "+init=epsg:4326"
    xy <- as.data.frame(spTransform(xy, crs(dem)))
    reso <- res(dem)[1]
    wsc36 <- 0
    wsc36atground <- 0
    for (i in 0:35) {
      wscr <- windcoef(dem, i*10, hgt = 1, res = reso)
      wscr2 <- windcoef(dem, i*10, hgt = 0, res = reso)
      wsc36[i + 1] <- extract(wscr, xy)
      wsc36atground[i + 1] <- extract(wscr2, xy)
    }
    wshelt <- 0
    wsheltatground <- 0
    for (i in 1:length(hourlydata$winddir)) {
      daz <- round(hourlydata$winddir[i] / 10, 0) + 1
      wshelt[i] <- wsc36[daz]
      wsheltatground[i] <- wsc36atground[daz]
    }
    if (class(slope) == "logical") {
      slope <- terrain(dem, unit = 'degrees')
      slope <- extract(slope, xy)
    }
    if (class(aspect) == "logical") {
      aspect <- terrain(dem, opt = 'aspect', unit = 'degrees')
      aspect <- extract(aspect, xy)
    }
    svf <- skyviewveg(dem, array(l, dim = dim(dem)[1:2]),
      array(x, dim = dim(dem)[1:2]), res = reso)
    fr <- canopy(array(l, dim = dim(dem)[1:2]), array(x, dim = dim(dem)[1:2]))
    svf <- extract(svf, xy)
    fr <- mean(fr)
    mslope <- mean_slope(dem, res = reso)
    mha <- extract(mslope, xy)
    ha36 <- 0
    for (i in 0:35) {
      har <- horizonangle(dem, i*10, reso)
      ha36[i + 1] <- atan(extract(har, xy)) * (180/pi)
    }
    tme <- as.POSIXlt(hourlydata$obs_time)
    ha <- 0
    jd <- julday(tme$year + 1900, tme$mon + 1, tme$mday)
    for (i in 1:length(tme)) {
      saz <- solazi(tme$hour[i], lat, long, jd[i])
      saz <- round(saz / 10, 0) + 1
      saz <- ifelse(saz > 36, 1, saz)
      ha[i] <- ha36[saz]
    }
    si <- siflat(tme$hour, lat, long, jd)
    sa <- solalt(tme$hour, lat, long, jd)
    si[sa < ha] <- 0
    dirr <- si * hourlydata$rad_dni
    a <- slope * (pi/180)
    k <- hourlydata$rad_dni / 4.87
    k <- ifelse(k > 1, 1, k)
    isor <- 0.5 * hourlydata$rad_dif * (1 + cos(a)) * (1 - k)
    cisr <- k * hourlydata$rad_dif * si
    sdi <- (slope + mha) * (pi/180)
    refr <- 0.5 * albr * (1 - cos(sdi)) * hourlydata$rad_dif
    fd <- dirr + cisr
    fdf <- isor + refr
    kk <- ((x^2 + 1/(tan(sa * (pi/180))^2))^0.5)/(x + 1.774 * (x + 1.182)^(-0.733))
    trd <- exp(-kk * l)
    trf <- (1 - fr)
    fgd <- fd * trd
    fged <- fdf * trf * svf
    swrad <- fgd + fged
    hourlyrad <- data.frame(swrad = swrad, skyviewfact = svf, canopyfact = fr,
      whselt = wsheltatground, windspeed = wshelt * hourlydata$windspeed,
      slope = slope, aspect = aspect)
    return(list(hourlyr = hourlyrad, slp = slope, asp = aspect, hori = ha36))
  }

  cadconditions2 <- function (em, wind, startjul, lat, long,
    starttime = 0, hourint = 1, windthresh = 4.5, emthresh = 0.725,
    tz = 0, dst = 0, con = TRUE)
  {
    jd <- floor(c(1:length(em)) * hourint/24 - hourint/24 + startjul +
        starttime/24)

    st <- suntimes(jd, lat, long, tz, dst)
    hrs <- (c(1:length(em)) * hourint - hourint + starttime)%%24
    dn <- ifelse(hrs > (st$sunrise + 3) & hrs < st$sunset, 1,
      0)
    cad <- ifelse(dn < 1 & wind < windthresh & em < emthresh,
      1, 0)
    if (hourint == 1 & con) {
      cad[2] <- ifelse(cad[1] & cad[2] == 1, 1, 0)
      cad <- c(cad[1:2], cad[3:length(cad)] * cad[2:(length(cad) -
          1)] * cad[1:(length(cad) - 2)])
    }
    cad
  }
  flowdir <- function (dem)
  {
    md <- is_raster(dem)
    fd <- md * 0
    md2 <- array(NA, dim = c(dim(md)[1] + 2, dim(md)[2] + 2))
    md2[2:(dim(md)[1] + 1), 2:(dim(md)[2] + 1)] <- md
    v <- c(1:length(md))
    v <- v[is.na(md) == F]
    x <- arrayInd(v, dim(md))[, 1]
    y <- arrayInd(v, dim(md))[, 2]
    for (i in 1:length(x)) {
      md9 <- md2[x[i]:(x[i] + 2), y[i]:(y[i] + 2)]
      fd[x[i], y[i]] <- round(mean(which(md9 == min(md9, na.rm = T))),
        0)
    }
    if_raster(fd, dem)
  }

  flowacc <- function (dem)
  {
    dm <- is_raster(dem)
    fd <- flowdir(dm)
    fa <- fd * 0 + 1
    o <- order(dm, decreasing = T, na.last = NA)
    for (i in 1:length(o)) {
      x <- arrayInd(o[i], dim(dm))[1]
      y <- arrayInd(o[i], dim(dm))[2]
      f <- fd[x, y]
      x2 <- x + (f - 1)%%3 - 1
      y2 <- y + (f - 1)%/%3 - 1
      if (x2 > 0 & x2 < dim(dm)[1] & y2 > 0 & y2 < dim(dm)[2])
        fa[x2, y2] <- fa[x, y] + 1
    }
    if_raster(fa, dem)
  }
  integertobinary8 <- function(i) {
    a <- 2 ^ (0:9)
    b <- 2 * a
    binc <- format(sapply(i, function(x) sum(10 ^ (0:9)[(x %% b) >= a])),
      scientific = FALSE)
    if (nchar(binc) > 8) warning("Integer > 8 bit binary")
    binc <- str_pad(binc, 8, pad = "0")
    binc
  }
  eleveffects <- function(hourlydata, dem, lat, long, windthresh = 4.5,
    emthresh = 0.78) {
    xy <- data.frame(x = long, y = lat)
    elevncep <- extract(demworld, xy)
    coordinates(xy) = ~x + y
    proj4string(xy) = "+init=epsg:4326"
    xy <- as.data.frame(spTransform(xy, crs(dem)))
    elev <- extract(dem, xy)
    # elevation effect
    lr <- lapserate(hourlydata$temperature, hourlydata$humidity,
      hourlydata$pressure)
    elevt <- lr * (elev - elevncep) + hourlydata$temperature
    tme <- as.POSIXlt(hourlydata$obs_time)
    jds <- julday(tme$year[1] + 1900, tme$mday[1] + 1, tme$mday[1])
    cad <- cadconditions2(hourlydata$emissivity, hourlydata$windspeed,
      jds, lat, long, starttime = tme$hour[1], hourint = 1,
      windthresh = windthresh, emthresh = emthresh)
    pxls <- dim(dem)[1] * dim(dem)[2]
    if (pxls > 300 * 300) {
      basins <- basindelin_big(dem)
    } else {
      basins <- basindelin(dem)
    }
    basins <- basinmerge(dem, basins, 2)
    basins <- basinsort(dem, basins)
    fa <- flowacc(dem)
    pfa <- is_raster(fa) * 0
    td <- is_raster(fa) * 0
    bm <- is_raster(basins)
    dm <- is_raster(dem)
    for (b in 1:max(bm, na.rm = TRUE)) {
      sel <- which(bm == b)
      fao <- log(is_raster(fa)[sel])
      pfa[sel] <- fao/max(fao, na.rm = TRUE)
      ed <- max(dm[sel], na.rm = TRUE) - dm[sel]
      td[sel] <- ed
    }
    cdif <- pfa * td
    cdif <- if_raster(cdif, dem)
    cdif <- extract(cdif, xy)
    cadt <- cdif * lr * cad
    tout <- data.frame(tref = hourlydata$temperature,
      elev = elev, elevncep = elevncep,
      telev = lr * elev,
      tcad = cadt)
    tout
  }

  ystart <- as.numeric(substr(dstart, 7, 10))
  yfinish <- as.numeric(substr(dfinish, 7, 10))
  yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
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

  # error trapping - originally inside the Fortran code, but now checking before executing Fortran
  errors<-0
  if(as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900") > Sys.time() | as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900") > Sys.time()){
    cat("sorry, no NCEP data for these times - please choose a different date range \n")
  }
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
  if(slope>90){
    cat("ERROR: Slope value (slope) is out of bounds.
        Please input a value between 0 and 90.", '\n')
    errors<-1
  }
  if(aspect>365){
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
  # end error trapping

  if(errors==0){ # continue

    moist.lai <- 0.1 # to do - option to decouple shade from LAI used in soil moisture
    # for microclima calculations
    tme <- as.POSIXlt(seq(ISOdate(ystart, 1, 1, tz = "UTC") - 3600 * 12, ISOdate((ystart + nyears), 1, 1, tz = "UTC") - 3600 * 13, by = "days"))



    ################## time related variables #################################
    doys12<-c(15.,46.,74.,105.,135.,166.,196.,227.,258.,288.,319.,349.) # middle day of each month

    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in

    daystart<-1
    idayst <- 1 # start day

    ################## location and terrain #################################
    if (!requireNamespace("raster", quietly = TRUE)) {
      stop("package 'raster' is needed. Please install it.",
        call. = FALSE)
    }
    if (!requireNamespace("RNCEP", quietly = TRUE)) {
      stop("package 'RNCEP' is needed. Please install it.",
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

    # get the local timezone reference longitude
    ALREF <- abs(trunc(x[1]))
    HEMIS <- ifelse(x[2]<0, 2, 1) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2])-ALAT)*60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1])-ALONG)*60
    azmuth<-aspect
    azmuth <- 0
    slope <- 0
    lat <- as.numeric(longlat[2])
    long <- as.numeric(longlat[1])
    loc <- c(long, lat)
    if(class(dem)[1] == "RasterLayer"){
      cat('using DEM provided to function call \n')
    }
    if(save != 2 & class(dem)[1] != "RasterLayer"){
      cat('downloading DEM via package elevatr \n')
      dem <- get_dem(lat = lat, long = long) # mercator equal area projection
    }
    if(save == 1){
      save(dem, file = 'dem.Rda')
    }
    if(save == 2){
      load('dem.Rda')
    }
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
    hori<-rep(0, 24)
    slope <- 0 # incorporated already by microclima
    azmuth <- 0 # incorporated already by microclima
    VIEWF <- 1 # incorporated already by microclima

    # setting up for temperature correction using lapse rate given difference between 9sec DEM value and 0.05 deg value
    days <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900"), as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900"), by = 'days')
    alldays <- seq(as.POSIXct("01/01/1900", format = "%d/%m/%Y", origin = "01/01/1900"), Sys.time()-60*60*24, by = 'days')
    startday <- which(as.character(format(alldays, "%d/%m/%Y")) == format(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900"), "%d/%m/%Y"))
    endday <- which(as.character(format(alldays, "%d/%m/%Y")) == format(as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900"), "%d/%m/%Y"))
    countday <- endday-startday+1
    #tt <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900"), as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900")+23*3600, by = 'hours')
    tt <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(dfinish, format = "%d/%m/%Y", tz = 'UTC')+23*3600, by = 'hours')
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
      cat("extracting weather data with RNCEP \n")
      #ncepdata <- get_inputs(lat, long, tme)
      #save(ncepdata, file = 'ncepdata.Rda')
      #hourlydata <- hourlyNCEP(ncepdata, lat, long)
      hourlyd <- hourlyNCEP(tme, longlat[2], longlat[1], reanalysis2, spatial)
      hourlydata <- hourlyd$hourlyout
      precip <- hourlyd$prate[-c(1:4)]
      precip <- precip[1:(length(precip)-4)]
      precip <- precip * 6 * 3600
      precip <- aggregate(precip, by = list(format(tt[seq(1, length(tt), 6)], "%Y-%m-%d")), sum)$x
      dailyrain <- precip
      ZENhr <- hourlydata$szenith
      ZENhr[ZENhr > 90] <- 90
      #dailyrain <- get_rain(lat, long, tme)
      #save(dailyrain, file = 'dailyrain.Rda')
      l <- mean(LAI) # leaf area index
      lor <- LOR # leaf orientation
      albr <- REFL # mean albedo of surfaces from which radiation is reflected

      cat("computing radiation and elevation effects with package microclima \n")
      hourlyr <- pointradwind(hourlydata, dem, lat, long, l, lor, albr)
      hourlyradwind <- hourlyr$hourlyr
      SLOPE <- hourlyr$slp
      ASPECT <- hourlyr$asp
      HORIZON <- hourlyr$hori
      if(save == 1){
        save(SLOPE, file = 'SLOPE.Rda')
        save(ASPECT, file = 'ASPECT.Rda')
        save(HORIZON, file = 'HORIZON.Rda')
      }
      # NB units for rad = MJ / m^2 / hr (divide by 0.0036 to get to W / m^2)
      # Skyviewfact (time invariant, 1 = complete hemisphere in view)
      # canopyfact (proportion of isotropic radiation blocked out, 1 = no radiation gets in)
      tref <- eleveffects(hourlydata, dem, lat, long)
      if(save == 1){
        save(tref, file = 'tref.Rda')
      }
      # tref (Original hourly NCEP)
      # telev (delta temperature due to elev)
      # tcad (delta temperature due to cold air drainage)
      elev <- tref$elev[1] # m
      ALTT <- elev
      TAIRhr <- tref$tref + tref$telev + tref$tcad # reference Tair plus offsets due to lapse rate and cold air drainage
      SOLRhr <- hourlyradwind$swrad / 0.0036
      SOLRhr[SOLRhr < 0] <- 0
      CLDhr <- hourlydata$h_tcdc
      CLDhr[CLDhr < 0] <- 0
      CLDhr[CLDhr > 100] <- 100
      #if(clearsky == 1){
      #  CLDhr <- CLDhr * 0
      #}
      IRDhr <- hourlydata$downlong / .0036
      RHhr <- suppressWarnings(humidityconvert(h = hourlydata$humidity, intype = 'specific', p = hourlydata$pressure, tc = TAIRhr)$relative)
      RHhr[RHhr > 100] <- 100
      RHhr[RHhr < 0] <- 0
      WNhr <- hourlyradwind$windspeed#*(1/10)^0.15 # adjust from 10m to 2m
      WNhr[is.na(WNhr)] <- 0.1
      RAINhr <- WNhr * 0 # using daily rain for now
      PRESShr <- hourlydata$pressure
      RAINFALL <- dailyrain
      RAINFALL[RAINFALL < 0.1] <- 0
      ZENhr2 <- ZENhr
      ZENhr2[ZENhr2!=90] <- 0
      # rleb <- rle(x = ZENhr2)$lengths[1:3]
      # length.orig <- length(PRESShr)
      # crop <- rleb[1] + (24-rleb[3])/2
      # crop2 <- length.orig-(25-crop)
      # CLDhr<-CLDhr[crop:crop2]
      # WNhr<-WNhr[crop:crop2]
      # TAIRhr<-TAIRhr[crop:crop2]
      # RHhr<-RHhr[crop:crop2]
      # RAINhr<-RAINhr[crop:crop2]
      # SOLRhr<-SOLRhr[crop:crop2]
      # ZENhr<-ZENhr[crop:crop2]
      # PRESShr<-PRESShr[crop:crop2]
      # tt<-tt[crop:crop2]
      # tt<-tt+(25-crop)*3600
      TMAXX1<-aggregate(TAIRhr,by=list(format(tt, "%d/%m/%Y")),max) # maximum air temperatures (°C)
      TMAXX<-TMAXX1$x[order(as.POSIXct(TMAXX1$Group.1, format = "%d/%m/%Y"))]
      TMINN1<-aggregate(TAIRhr,by=list(format(tt, "%d/%m/%Y")),min) # minimum air temperatures (°C)
      TMINN<-TMINN1$x[order(as.POSIXct(TMINN1$Group.1, format = "%d/%m/%Y"))]
      #RAINFALL1<-aggregate(data.frame(RAINhr2 * 6 * 3600),by=list(format(tointerp, "%d/%m/%Y")),sum) # monthly mean rainfall (mm)
      #RAINFALL<-RAINFALL1[order(as.POSIXct(RAINFALL1$Group.1, format = "%d/%m/%Y")),2]
      CCMAXX1<-aggregate(CLDhr,by=list(format(tt, "%d/%m/%Y")),max) # max cloud cover (%)
      CCMAXX<-CCMAXX1$x[order(as.POSIXct(CCMAXX1$Group.1, format = "%d/%m/%Y"))]
      CCMINN1<-aggregate(CLDhr,by=list(format(tt, "%d/%m/%Y")),min) # min cloud cover (%)
      CCMINN<-CCMINN1$x[order(as.POSIXct(CCMINN1$Group.1, format = "%d/%m/%Y"))]
      RHMAXX1<-aggregate(RHhr,by=list(format(tt, "%d/%m/%Y")),max) # max relative humidity (%)
      RHMAXX<-RHMAXX1$x[order(as.POSIXct(RHMAXX1$Group.1, format = "%d/%m/%Y"))]
      RHMINN1<-aggregate(RHhr,by=list(format(tt, "%d/%m/%Y")),min) # min relative humidity (%)
      RHMINN<-RHMINN1$x[order(as.POSIXct(RHMINN1$Group.1, format = "%d/%m/%Y"))]
      WNMAXX1<-aggregate(WNhr,by=list(format(tt, "%d/%m/%Y")),max) # max wind speed (m/s)
      WNMAXX<-WNMAXX1$x[order(as.POSIXct(WNMAXX1$Group.1, format = "%d/%m/%Y"))]
      WNMINN1<-aggregate(WNhr,by=list(format(tt, "%d/%m/%Y")),min) # min wind speed (m/s)
      WNMINN<-WNMINN1$x[order(as.POSIXct(WNMINN1$Group.1, format = "%d/%m/%Y"))]
      PRESS1<-aggregate(PRESShr,by=list(format(tt, "%d/%m/%Y")),min) # min wind speed (m/s)
      PRESS<-PRESS1$x[order(as.POSIXct(PRESS1$Group.1, format = "%d/%m/%Y"))]

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
    }

    ndays<-length(TMAXX)
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
    #cut <- as.numeric(days[1] - as.POSIXct(paste0('01/01/', ystart), format = "%d/%m/%Y") + 1)
    #doy <- doy[(cut+1):(cut+countday-1)]
    #RAINFALL<-RAINFALL[(cut+1):(cut+countday-1)]
    ida<-ndays
    idayst <- 1

    dim<-length(TMAXX)
    maxshades<-rep(0.1,dim)
    minshades<-rep(0,dim)
    shademax<-maxshades
    maxshade<-.1
    minshade<-0

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

    if(warm != 0){
      # impose uniform temperature change
      TMAXX<-TMAXX+seq(0, dim-1)/(dim-1)*warm
      TMINN<-TMINN+seq(0, dim-1)/(dim-1)*warm
      TAIRhr<-TAIRhr+seq(0, dim-1)/(dim-1)*warm
    }
    RAINFALL<-RAINFALL+rainoff
    RAINhr<-RAINhr+rainoff
    ALLMINTEMPS<-TMINN
    ALLMAXTEMPS<-TMAXX
    ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)

    WNMAXX <- WNMAXX * windfac
    WNMINN <- WNMINN * windfac
    WNhr <- WNhr * windfac

    MAXSHADES<-maxshades
    MINSHADES<-minshades

    REFLS <- rep(REFL, dim)
    PCTWET <- rep(PCTWET, dim)
    #RAINFALL <- RAINFALL[(1+cut):(dim+cut)]
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
    slope <- 0 # microclima deals with this
    azmuth <- 0 # microclima deals with this

    avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
    soilinit<-rep(avetemp,20)
    tannul<-mean(unlist(ALLTEMPS))

    if(nyears==1){
      avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
      tannulrun<-rep(avetemp,ndays)
    }else{
      avetemp<-rowMeans(cbind(TMAXX, TMINN), na.rm=TRUE)
      if(length(TMAXX)<365){
        tannulrun<-rep((sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2),length(TMAXX))
      }else{
        tannulrun<-movingFun(avetemp,n=365,fun=mean,type='to')
        yearone<-rep((sum(TMAXX[1:365])+sum(TMINN[1:365]))/(365*2),365)
        tannulrun[1:365]<-yearone
      }
    }

    SLES<-matrix(nrow = dim, data = 0)
    SLES<-SLES+SLE

    moists2<-matrix(nrow=10, ncol = ndays, data=0)
    moists2[1,ndays]<-0.2
    moists<-moists2

    if(runmoist==1){
      moists2<-matrix(nrow=10, ncol = dim, data=0) # set up an empty vector for soil moisture values through time
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
    runshade <- 0
    IR <- 0
    # microclimate input parameters list
    microinput<-c(dim,RUF,ERR,Usrhyt,Refhyt,Numtyps,Z01,Z02,ZH1,ZH2,idayst,ida,HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,CMH2O,microdaily,tannul,EC,VIEWF,snowtemp,snowdens,snowmelt,undercatch,rainmult,runshade,runmoist,maxpool,evenrain,snowmodel,rainmelt,writecsv,densfun,hourly,rainhourly,lamb,IUV,RW,PC,RL,SP,R1,IM,MAXCOUNT,IR,message,fail,snowcond,intercept,grasshade)

    # hourly option set to 0, so make empty vectors
    if(hourly==0){
      TAIRhr=rep(0,24*dim)
      RHhr=rep(0,24*dim)
      WNhr=rep(0,24*dim)
      CLDhr=rep(0,24*dim)
      SOLRhr=rep(0,24*dim)
      ZENhr=rep(-1,24*dim)
      IRDhr=rep(-1,24*dim)
    }else{
      CLDhr=rep(0,24*dim)
    }
    if(rainhourly==0){
      RAINhr=rep(0,24*dim)
    }else{
      RAINhr = RAINhr
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
    LAI <- moist.lai
    if(length(LAI)<dim){
      LAI<-rep(LAI[1],dim)
      LAI1 <- LAI
    }
    if(shore==0){
      tides<-matrix(data = 0, nrow = 24*dim, ncol = 3) # make an empty matrix
    }
    TIMAXS <- c(1.0, 1.0, 0.0, 0.0)
    TIMINS <- c(0, 0, 1, 1)
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
    cat(paste('running microclimate model for ',dim,' days from ',tt[1],' to ', tt[length(tt)], ' at site ',location,'\n'))
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
    if(lamb == 1){
      drlam<-as.data.frame(microut$drlam) # retrieve direct solar irradiance
      drrlam<-as.data.frame(microut$drrlam) # retrieve direct Rayleigh component solar irradiance
      srlam<-as.data.frame(microut$srlam) # retrieve scattered solar irradiance
      if(snowmodel == 1){
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,dim=dim,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=longlat,nyears=nyears,minshade=minshade,maxshade=maxshade,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam, SLOPE = SLOPE, ASPECT = ASPECT, HORIZON = HORIZON, dates=tt, dem = dem,dates2=dates2))
      }else{
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,dim=dim,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=longlat,nyears=nyears,minshade=minshade,maxshade=maxshade,DEP=DEP,drlam=drlam,drrlam=drrlam,srlam=srlam, SLOPE = SLOPE, ASPECT = ASPECT, HORIZON = HORIZON, dates=tt, dem = dem,dates2=dates2))
      }
    }else{
      if(snowmodel == 1){
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,sunsnow=sunsnow,shdsnow=shdsnow,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,dim=dim,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=longlat,nyears=nyears,minshade=minshade,maxshade=maxshade,DEP=DEP, SLOPE = SLOPE, ASPECT = ASPECT, HORIZON = HORIZON, dates=tt, dem = dem,dates2=dates2))
      }else{
        return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,humid=humid,shadhumid=shadhumid,soilpot=soilpot,shadpot=shadpot,plant=plant,shadplant=shadplant,RAINFALL=RAINFALL,dim=dim,elev=ALTT,REFL=REFL[1],MAXSHADES=MAXSHADES,longlat=longlat,nyears=nyears,minshade=minshade,maxshade=maxshade,DEP=DEP, SLOPE = SLOPE, ASPECT = ASPECT, HORIZON = HORIZON, dates=tt, dem = dem,dates2=dates2))
      }
    }
  } # end error trapping
} # end of micro_clima function
