## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ---- message=FALSE, warnings=FALSE--------------------------------------
library(NicheMapR) 
longlat <- c(-89.40123, 43.07305) # Madison, Wisconsin, USA
micro <- micro_global(loc = longlat)

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE----------
knitr::kable(head(micro$metout[, 1:9], 2), digits = 2)
knitr::kable(head(micro$metout[, 10:18], 2), digits = 2)

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE----------
knitr::kable(head(micro$soil[, 1:12], 2), digits = 2)

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
require(lattice, quietly = TRUE)
soil <- as.data.frame(micro$soil) # get the soil data
minshade <- micro$minshade # get the value for minimum shade
with(subset(soil, DOY==196 | DOY==349), {xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(DOY), ylim = c(-20, 70), xlab = "Time of Day (min)", ylab = "Soil Temperature (deg C)", auto.key = list(columns = 5), as.table = TRUE, type = "b", main = paste(minshade,"% shade"))})

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
micro <- micro_global(loc = longlat, runshade = 0, minshade = 50)
soil <- as.data.frame(micro$soil) # get the soil data
minshade<-micro$minshade # get the value for minimum shade
with(subset(soil, DOY==196 | DOY==349),{xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(DOY), ylim = c(-20, 70), xlab = "Time of Day (min)", ylab = "Soil Temperature (deg C)", auto.key=list(columns = 5), as.table = TRUE, type = "b", main = paste(minshade, "% shade"))})

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
micro <- micro_global(loc = longlat, runshade = 0, minshade = 0, slope = 45, aspect = 180)
soil <- as.data.frame(micro$soil) # get the soil data
minshade <- micro$minshade # get the value for minimum shade
with(subset(soil, DOY==196 | DOY==349), {xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(DOY), ylim = c(-20, 70), xlab = "Time of Day (min)", ylab = "Soil Temperature (deg C)", auto.key=list(columns = 5), as.table = TRUE, type = "b", main = paste(minshade, "% shade, 45 degree slope, 180 degrees aspect"))})

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
micro <- micro_global(loc = longlat, runshade = 0, minshade = 0, hori = c(0, 0, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 0, 0, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65))
soil <- as.data.frame(micro$soil) # get the soil data
minshade <- micro$minshade # get the value for minimum shade
with(subset(soil,DOY==196 | DOY==349),{xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(DOY), ylim = c(-20, 70), xlab = "Time of Day (min)", ylab = "Soil Temperature (deg C)", auto.key = list(columns = 5), as.table = TRUE, type = "b", main=paste(minshade,"% shade, north-south gully"))})

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
micro <- micro_global(runshade = 0, minshade = 0, soiltype = 0)
soil <- as.data.frame(micro$soil) # get the soil data
minshade <- micro$minshade # get the value for minimum shade
with(subset(soil, DOY==196 | DOY==349), {xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(DOY), ylim = c(-20, 70), xlab = "Time of Day (min)", ylab = "Soil Temperature (deg C)",  auto.key = list(columns = 5), as.table = TRUE, type = "b", main = paste(minshade, "% shade, rock substrate"))})

## ---- message=FALSE, warnings=FALSE--------------------------------------
micro <- micro_global(loc = longlat, runmoist = 1)

## ---- message=FALSE, warnings=FALSE, echo=FALSE, results='asis'----------
knitr::kable(head(micro$soilmoist[, 1:9], 2))
knitr::kable(head(micro$shadmoist[, 10:12], 2))

## ---- message=FALSE, warnings=FALSE, echo=FALSE, results='asis'----------
knitr::kable(head(micro$soilpot[, 1:9], 2))
knitr::kable(head(micro$shadpot[, 10:12], 2))

## ---- message=FALSE, warnings=FALSE, echo=FALSE, results='asis'----------
knitr::kable(head(micro$humid[, 1:9], 2))
knitr::kable(head(micro$shadhumid[, 10:12], 2))

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
soilmoist <- as.data.frame(micro$soilmoist) # get the minimum shade soil moisture output
minshade <- micro$minshade # get the value for minimum shade
# append dates
days <- rep(seq(1, 12), 24)
days <- days[order(days)]
dates <- days + soilmoist$TIME / 60 / 24 - 1 # dates for hourly output
soilmoist<-cbind(dates, soilmoist)
for(i in 1:10){
  if(i == 1){
    plot(soilmoist[, i + 3] ~ soilmoist[, 1], ylim = c(0, 0.5), xlab = "Month",  ylab = "Soil Moisture (% vol)", col = i, type = "l", main = paste("soil moisture, ", minshade, "% shade"))
  }else{
    points(soilmoist[, i + 3] ~ soilmoist[, 1], col = i, type = "l")
  }
}

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
micro<-micro_global(loc = longlat, runmoist = 1, soiltype = 11)
soilmoist<-as.data.frame(micro$soilmoist) # get the minimum shade soil moisture output
minshade<-micro$minshade # get the value for minimum shade
# append dates
days<-rep(seq(1,12),24)
days<-days[order(days)]
dates<-days+soilmoist$TIME/60/24-1 # dates for hourly output
soilmoist<-cbind(dates,soilmoist)
for(i in 1:10){
  if(i==1){
    plot(soilmoist[,i+3]~soilmoist[,1], ylim=c(0,0.5),xlab = "Month", ylab = "Soil Moisture (%
      vol)",col=i,type = "l",main=paste("soil moisture, ",minshade,"% shade"))
  }else{
    points(soilmoist[,i+3]~soilmoist[,1] ,col=i,type = "l")
  }
}

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
nyears<-5
micro<-micro_global(loc = longlat, runmoist = 1, soiltype = 11, nyears = nyears)
soilmoist<-as.data.frame(micro$soilmoist) # get the minimum shade soil moisture output
minshade<-micro$minshade # get the value for minimum shade
# append dates
days<-rep(seq(1,12*nyears)/12,24)
days<-days[order(days)]
dates<-days+soilmoist$TIME/60/24/(12) # dates for hourly output
soilmoist<-cbind(dates,soilmoist)
for(i in 1:10){
  if(i==1){
    plot(soilmoist[,i+3]~soilmoist[,1], ylim=c(0,0.5),xlab = "Year", ylab = "Soil Moisture (%
      vol)",col=i,type = "l",main=paste("soil moisture, ",minshade,"% shade"))
  }else{
    points(soilmoist[,i+3]~soilmoist[,1],col=i,type = "l")
  }
}

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
timeinterval <- 365
micro <- micro_global(loc = longlat, runmoist = 1, soiltype = 11, timeinterval = timeinterval)
soilmoist <- as.data.frame(micro$soilmoist)  # get the minimum shade soil moisture output
minshade <- micro$minshade  # get the value for minimum shade
# append dates
days <- rep(seq(1, timeinterval), 24)
days <- days[order(days)]
dates <- days + soilmoist$TIME/60/24 - 1  # dates for hourly output
soilmoist <- cbind(dates, soilmoist)
for (i in 1:10) {
    if (i == 1) {
        plot(soilmoist[, i + 3] ~ soilmoist[, 1], ylim = c(0, 0.5), xlab = "Day of Year", ylab = "Soil Moisture (% vol)", col = i, type = "l", main = paste("soil moisture, ", minshade, "% shade"))
    } else {
        points(soilmoist[, i + 3] ~ soilmoist[, 1], col = i, type = "l")
    }
}

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
soil <- as.data.frame(micro$soil)  # get the minimum shade soil temperature output
minshade <- micro$minshade  # get the value for minimum shade
# append dates
days <- rep(seq(1, timeinterval), 24)
days <- days[order(days)]
dates <- days + soil$TIME/60/24 - 1  # dates for hourly output
soil <- cbind(dates, soil)
for (i in 1:10) {
    if (i == 1) {
        plot(soil[, i + 3] ~ soil[, 1], ylim = c(-20, 70), xlab = "Day of Year", ylab = "Soil Temperature (deg C)", col = i, type = "l", main = paste("soil temperature, ", minshade, "% shade"))
    } else {
        points(soil[, i + 3] ~ soil[, 1], col = i, type = "l")
    }
}

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
timeinterval <- 365
nyears <- 2  # running two years, first one acting as a 'burn in' year and discarded
micro <- micro_global(loc = longlat, runmoist = 1, snowmodel = 1, timeinterval = timeinterval, nyears = 2)
soil <- as.data.frame(micro$soil)[(365 * 24 + 1):(365 * 24 * nyears), ]  # get the minimum shade soil temperature output, discarding the first year
metout <- as.data.frame(micro$metout)[(365 * 24 + 1):(365 * 24 * nyears), ]  # get the minimum shade above ground conditions, discarding the first year
minshade <- micro$minshade  # get the value for minimum shade
# append dates
days <- rep(seq(366, timeinterval * nyears), 24)
days <- days[order(days)]
dates <- days + soil$TIME / 60 / 24 - 1  # dates for hourly output
soil <- cbind(dates, soil)
metout <- cbind(dates, metout)
for (i in 1:10) {
    if (i == 1) {
        plot(soil[, i + 3] ~ soil[, 1], ylim = c(-20, 70), xlab = "Day of Year", ylab = "Soil Temperature (deg C)", col = i, type = "l", main = paste("soil  temperature, ", minshade, "% shade"))
    } else {
        points(soil[, i + 3] ~ soil[, 1], col = i, type = "l")
    }
}
plot(metout$SNOWDEP ~ metout$dates, xlab = "Time of Day (min)", ylab = "snow depth, cm / snow fall, mm", type = "h", main = paste("snow depth (cm) and snow fall (mm) ", minshade, "% shade", sep = ""), col = "light blue")
points(metout$SNOWFALL ~ metout$dates, xlab = "Time of Day (min)", type = "h", col = "blue")

## ---- message=FALSE, warnings=FALSE, fig.width=7, fig.height=6-----------
nyears<-1
timeinterval<-365
# if you have data on tides for a shoreline, use the code below to specify the hours the tide is in, the water temperature and wave splash
  tides<-matrix(data = 0., nrow = 24*timeinterval*nyears, ncol = 3) # make an empty matrix
  tides[,2]<-5 # set a constant sea water temperature
  tides[,3]<-0 # zero the wave splash column
  tides[12,3]<-90 # make the 12th value (hour 12 on 1st day) get a spash, resulting in the rock evaporating as if it is 90% wetted
  mocktides<-rep(c(rep(0,12),rep(1,13)),timeinterval*nyears) # made a sequence of tides, where 1 means tide is in and 2 means out, and offset it to 24 hour cycle
  mocktides<-mocktides[1:8760] # subset it to a year long of 24 hour tides
  tides[1:(timeinterval*nyears*24),1]<-mocktides # put the mock tides in the tides vector, column 1
micro<-micro_global(loc=longlat, shore = 1, soiltype = 0, timeinterval = timeinterval, nyears = nyears, tides = tides)
soil<-as.data.frame(micro$soil) # get the minimum shade soil temperature output first year
minshade<-micro$minshade # get the value for minimum shade
# append dates
days<-rep(seq(1,timeinterval),24)
days<-days[order(days)]
dates<-days+soil$TIME/60/24-1 # dates for hourly output
soil<-cbind(dates,soil)
for(i in 1:10){
  if(i==1){
    plot(soil[,i+3]~soil[,1], ylim=c(-20,70),xlab = "Day of Year", ylab = "Soil Temperature (deg
      C)",col=i,type = "l",main=paste("soil temperature, ",minshade,"% shade"))
  }else{
    points(soil[,i+3]~soil[,1], col=i,type = "l")
  }
}  

