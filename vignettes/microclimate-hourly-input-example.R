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

