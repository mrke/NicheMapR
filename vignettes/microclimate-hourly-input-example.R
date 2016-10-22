## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ------------------------------------------------------------------------
library(NicheMapR)
library(zoo)

## ------------------------------------------------------------------------
head(SCANsites)

## ------------------------------------------------------------------------
sitenum='2184' # Ford Dry Lake
site=subset(SCANsites,id==sitenum) # subset the SCANsites dataset for Ford Dry Lake
name=site$name # the name of the sites
Latitude<-site$lat # the latitude in decimal degrees
Longitude<-site$lon # the longitude in decimal degrees
Elevation<-site$elev/3.28084 # elevation, converted from feet to metres
TZoffset<-site$`GMT offset` # the offset from Greenwich Mean Time, in hours
ystart=2015 # start yeaar
yfinish=2015 # end year
nyears=yfinish-ystart+1 # number of years to run


## ------------------------------------------------------------------------
weather<-SCAN_FordDryLake_2015 # make SCAN_FordDrylake_2015 supplied package data the weather input variable

## ------------------------------------------------------------------------
writecsv<-0 # make Fortran code write output as csv files
runshade<-1 # run the model twice, once for each shade level (1) or just for the first shade level (0)?
runmoist<-1 # run soil moisture model (0=no, 1=yes)?
snowmodel<-1 # run the snow model (0=no, 1=yes)? - note that this runs slower
hourly<-1 # run the model with hourly input data
microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day

## ------------------------------------------------------------------------
longlat<-c(Longitude,Latitude) # decimal degrees longitude and latitude from the SCAN site data table
julnum<-floor(nrow(weather)/24) # number of days to run, determined by counting the number of rows in the weather dataset and dividing by 24 to get days, but keeping it as a whole number
idayst <- 1 # start month
ida<-julnum # end month
HEMIS <- ifelse(longlat[2]<0,2.,1.) # chose hemisphere based on latitude
ALAT <- abs(trunc(longlat[2])) # degrees latitude
AMINUT <- (abs(longlat[2])-ALAT)*60 # minutes latitude
ALONG <- abs(trunc(longlat[1])) # degrees longitude
ALMINT <- (abs(longlat[1])-ALONG)*60 # minutes latitude
ALREF <- ALONG # reference longitude for time zone

## ------------------------------------------------------------------------
EC <- 0.0167238 # Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)
RUF <- 0.004 # Roughness height (m), , e.g. sand is 0.05, grass may be 2.0, current allowed range: 0.001 (snow) - 2.0 cm.
# Next four parameters are segmented velocity profiles due to bushes, rocks etc. on the surface
#IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO! (then roughness height is based on the parameter RUF)
Z01 <- 0. # Top (1st) segment roughness height(m)
Z02 <- 0. # 2nd segment roughness height(m)
ZH1 <- 0. # Top of (1st) segment, height above surface(m)
ZH2 <- 0. # 2nd segment, height above surface(m)
SLE <- 0.96 # Substrate longwave IR emissivity (decimal %), typically close to 1
ERR <- 1.5 # Integrator error for soil temperature calculations
Refhyt <- 2 # Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured
DEP <- c(0., 2.5,  5., 10., 15., 20.,  30.,  50.,  100.,  200.) # Soil nodes (cm) - keep spacing close near the surface, last value is where it is assumed that the soil temperature is at the annual mean air temperature
Thcond <- 2.5 # soil minerals thermal conductivity (W/mC)
Density <- 2560. # soil minerals density (kg/m3)
SpecHeat <- 870. # soil minerals specific heat (J/kg-K)
BulkDensity <- 1300 # soil bulk density (kg/m3)
SatWater <- 0.26 # volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
Clay <- 20 # clay content for matric potential calculations (%)
REFL<-0.20 # soil reflectance (decimal %)
ALTT<-Elevation # altitude (m)
slope<-0. # slope (degrees, range 0-90)
azmuth<-180. # aspect (degrees, 0 = North, range 0-360)
hori<-rep(0,24) # enter the horizon angles (degrees) so that they go from 0 degrees azimuth (north) clockwise in 15 degree intervals
VIEWF <- 1-sum(sin(hori*pi/180))/length(hori) # convert horizon angles to radians and calc view factor(s)
PCTWET<-0 # percentage of surface area acting as a free water surface (%)
CMH2O <- 1. # precipitable cm H2O in air column, 0.1 = VERY DRY; 1.0 = MOIST AIR CONDITIONS; 2.0 = HUMID, TROPICAL CONDITIONS (note this is for the whole atmospheric profile, not just near the ground)
TIMAXS <- c(1.0, 1.0, 0.0, 0.0)   # Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise
TIMINS <- c(0.0, 0.0, 1.0, 1.0)   # Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon
minshade<-0. # minimum available shade (%)
maxshade<-90. # maximum available shade (%)
Usrhyt <- 0.01# local height (m) at which air temperature, relative humidity and wind speed calculatinos will be made
# Aerosol profile using GADS
relhum<-1.
optdep.summer<-as.data.frame(rungads(longlat[2],longlat[1],relhum,0))
optdep.winter<-as.data.frame(rungads(longlat[2],longlat[1],relhum,1))
optdep<-cbind(optdep.winter[,1],rowMeans(cbind(optdep.summer[,2],optdep.winter[,2])))
optdep<-as.data.frame(optdep)
colnames(optdep)<-c("LAMBDA","OPTDEPTH")
a<-lm(OPTDEPTH~poly(LAMBDA, 6, raw=TRUE),data=optdep)
LAMBDA<-c(290,295,300,305,310,315,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1380,1400,1420,1440,1460,1480,1500,1540,1580,1600,1620,1640,1660,1700,1720,1780,1800,1860,1900,1950,2000,2020,2050,2100,2120,2150,2200,2260,2300,2320,2350,2380,2400,2420,2450,2490,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000)
TAI<-predict(a,data.frame(LAMBDA))

## ------------------------------------------------------------------------
# check if first element is NA and, if so, use next non-NA value for na.approx function
if(is.na(weather$TAVG.H[1])==TRUE){ # mean hourly air temperature
  weather$TAVG.H[1]<-weather$TAVG.H[which(!is.na(weather$TAVG.H))[1]]
}

if(is.na(weather$WSPDV.H[1])==TRUE){ # mean hourly wind speed
  weather$WSPDV.H[1]<-weather$WSPDV.H[which(!is.na(weather$WSPDV.H))[1]]
}
if(is.na(weather$RHUM[1])==TRUE){ # mean hourly relative humidity
  weather$RHUM[1]<-weather$RHUM[which(!is.na(weather$RHUM))[1]]
}
if(is.na(weather$SRADV.H[1])==TRUE){ # mean hourly solar radiation
  weather$SRADV.H[1]<-weather$SRADV.H[which(!is.na(weather$SRADV.H))[1]]
}
if(is.na(weather$PRCP.H[1])==TRUE){ # hourly precipitation
  weather$PRCP.H[1]<-weather$PRCP.H[which(!is.na(weather$PRCP.H))[1]]
}
# use na.approx function from zoo package to fill in missing data
TAIRhr<-weather$TAVG.H<-na.approx(weather$TAVG.H)
RHhr<-weather$RHUM.I<-na.approx(weather$RHUM.I)
SOLRhr<-weather$SRADV.H<-na.approx(weather$SRADV.H)
RAINhr<-weather$PRCP.H<-na.approx(weather$PRCP.H*25.4) # convert rainfall from inches to mm
WNhr<-weather$WSPDV.H<-na.approx(weather$WSPDV.H*0.44704) # convert wind speed from miles/hour to m/s

## ------------------------------------------------------------------------
# run global microclimate model in clear sky mode to get clear sky radiation
micro<-micro_global(loc=c(Longitude,Latitude),timeinterval = 365, clearsky = 1)
# append dates
tzone=paste("Etc/GMT",TZoffset,sep="") # doing it this way ignores daylight savings!
dates=seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
clear<-as.data.frame(cbind(dates,as.data.frame(rep(micro$metout[1:(365*24),13],nyears))),stringsAsFactors = FALSE)
julday<-rep(seq(1,365),nyears)[1:floor(nrow(weather)/24)] # julian days to run
clear=as.data.frame(clear,stringsAsFactors = FALSE)
colnames(clear)=c("datetime","sol")

# find the maximum observed solar and adjust the clear sky prediction to this value 
maxsol=max(SOLRhr)
clear2<-clear[,2]*(maxsol/max(clear[,2])) # get ratio of max observed to predicted max clear sky solar

# compute cloud cover from ratio of max to observed solar
sol=SOLRhr
sol[sol<5]<-0 # remove very low values
clr<-clear2
clr[clr<5]<-0 # remove very low values
a=((clr-sol)/clr)*100 # get ratio of observed to predicted solar, convert to %
a[a>100]<-100 # cap max 100%
a[a<0]<-0 # cap min at 0%
a[is.na(a)]=0 # replace NA with zero
a[is.infinite(a)]=0 # replace infinity with zero
a[a==0]=NA # change all zeros to NA for na.approx
a=na.approx(a,na.rm = FALSE) # apply na.approx, but leave any trailing NAs
a[is.na(a)]=0 # make trailing NAs zero
CLDhr<-a # now we have hourly cloud cover

## ------------------------------------------------------------------------
# aggregate hourly data to daily min/max
CCMAXX<-aggregate(CLDhr,by=list(weather$Date), FUN=max)[,2]#c(100,100) # max cloud cover (%)
CCMINN<-aggregate(CLDhr,by=list(weather$Date), FUN=min)[,2]#c(0,15.62) # min cloud cover (%)
TMAXX<-aggregate(TAIRhr,by=list(weather$Date), FUN=max)[,2]#c(40.1,31.6) # maximum air temperatures (deg C)
TMINN<-aggregate(TAIRhr,by=list(weather$Date), FUN=min)[,2]#c(19.01,19.57) # minimum air temperatures (deg C)
RAINFALL<-aggregate(RAINhr,by=list(weather$Date), FUN=sum)[,2]#c(19.01,19.57) # minimum air temperatures (deg C)
RHMAXX<-aggregate(RHhr,by=list(weather$Date), FUN=max)[,2]#c(90.16,80.92) # max relative humidity (%)
RHMINN<-aggregate(RHhr,by=list(weather$Date), FUN=min)[,2]#c(11.05,27.9) # min relative humidity (%)
WNMAXX<-aggregate(WNhr,by=list(weather$Date), FUN=max)[,2]#c(1.35,2.0) # max wind speed (m/s)
WNMINN<-aggregate(WNhr,by=list(weather$Date), FUN=min)[,2]#c(0.485,0.610) # min wind speed (m/s)

## ------------------------------------------------------------------------
tannul<-mean(c(TMAXX,TMINN)) # annual mean temperature for getting monthly deep soil temperature (deg C)
tannulrun<-rep(tannul,julnum) # monthly deep soil temperature (2m) (deg C)
# creating the arrays of environmental variables that are assumed not to change with month for this simulation
MAXSHADES <- rep(maxshade,julnum) # daily max shade (%)
MINSHADES <- rep(minshade,julnum) # daily min shade (%)
SLES <- rep(SLE,julnum) # set up vector of ground emissivities for each day
REFLS<-rep(REFL,julnum) # set up vector of soil reflectances for each day
PCTWET<-rep(PCTWET,julnum) # set up vector of soil wetness for each day

## ------------------------------------------------------------------------
# set up a profile of soil properites with depth for each day to be run
Numtyps <- 1 # number of soil types
Nodes <- matrix(data = 0, nrow = 10, ncol = julnum) # array of all possible soil nodes for max time span of 20 years
Nodes[1,1:julnum]<-10 # deepest node for first substrate type
Density<-Density/1000 # density of minerals - convert to Mg/m3
BulkDensity<-BulkDensity/1000 # density of minerals - convert to Mg/m3

# now make the soil properties matrix
# columns are:
#1) bulk density (Mg/m3)
#2) volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
#3) clay content (%)
#4) thermal conductivity (W/mK)
#5) specific heat capacity (J/kg-K)
#6) mineral density (Mg/m3)
soilprops<-matrix(data = 0, nrow = 10, ncol = 6) # create an empty soil properties matrix
soilprops[1,1]<-BulkDensity # insert soil bulk density to profile 1
soilprops[1,2]<-SatWater # insert saturated water content to profile 1
soilprops[1,3]<-Clay     # insert percent clay to profile 1
soilprops[1,4]<-Thcond # insert thermal conductivity to profile 1
soilprops[1,5]<-SpecHeat # insert specific heat to profile 1
soilprops[1,6]<-Density # insert mineral density to profile 1
soilinit<-rep(tannul,20) # make iniital soil temps equal to mean annual

## ------------------------------------------------------------------------
#use Campbell and Norman Table 9.1 soil moisture properties
soiltype=3 # 3 = sandy loam
PE<-rep(CampNormTbl9_1[soiltype,4],19) #air entry potential J/kg
KS<-rep(CampNormTbl9_1[soiltype,6],19) #saturated conductivity, kg s/m3
BB<-rep(CampNormTbl9_1[soiltype,5],19) #soil 'b' parameter
BD<-rep(BulkDensity,19) # soil bulk density, Mg/m3
soiltype=5 # change deeper nodes to 5 = a silt loam
PE[10:19]<-CampNormTbl9_1[soiltype,4] #air entry potential J/kg
KS[10:19]<-CampNormTbl9_1[soiltype,6] #saturated conductivity, kg s/m3
BB[10:19]<-CampNormTbl9_1[soiltype,5] #soil 'b' parameter

L<-c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,0,0)*10000 # root density at each node, mm/m3 (from Campell 1985 Soil Physics with Basic, p. 131)
LAI<-0.1 # leaf area index, used to partition traspiration/evaporation from PET
rainmult<-1 # rainfall multiplier to impose catchment
maxpool<-10 # max depth for water pooling on the surface, mm (to account for runoff)
evenrain<-0 # spread daily rainfall evenly across 24hrs (1) or one event at midnight (2)
SoilMoist_Init<-rep(0.2,10) # initial soil water content for each node, m3/m3
moists<-matrix(nrow=10, ncol = julnum, data=0) # set up an empty vector for soil moisture values through time
moists[1:10,]<-SoilMoist_Init # insert inital soil moisture

## ------------------------------------------------------------------------
snowtemp<-1.5 # temperature at which precipitation falls as snow (used for snow model)
snowdens<-0.375 # snow density (mg/m3)
densfun<-c(0,0) # slope and intercept of linear model of snow density as a function of day of year - if it is c(0,0) then fixed density used
snowmelt<-0.9 # proportion of calculated snowmelt that doesn't refreeze
undercatch<-1. # undercatch multipier for converting rainfall to snow
rainmelt<-0.0125 # parameter in equation from Anderson's SNOW-17 model that melts snow with rainfall as a function of air temp

## ------------------------------------------------------------------------
# intertidal simulation input vector (col 1 = tide in(1)/out(0), col 2 = sea water temperature in deg C, col 3 = % wet from wave splash)
tides<-matrix(data = 0., nrow = 24*julnum, ncol = 3) # matrix for tides

