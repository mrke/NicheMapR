## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ---- echo=FALSE---------------------------------------------------------
#dyn.load(paste(lib.loc = .libPaths()[1],'/NicheMapR/libs/x64/microclimate.dll',sep=""))
#library(GADS) # for some reason it won't automatically load GADS via the require statement

## ------------------------------------------------------------------------
library(NicheMapR)

micro<-micro_global(loc="Madison, Wisconsin, USA")

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(micro$metout[,1:9], 2))
knitr::kable(head(micro$metout[,10:18], 2))

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(micro$soil[,1:9], 2))
knitr::kable(head(micro$soil[,10:12], 2))

## ---- fig.width=7, fig.height=6------------------------------------------
require(lattice,quietly = TRUE)
soil<-as.data.frame(micro$soil) # get the soil data
minshade<-micro$minshade # get the value for minimum shade
with(subset(soil,JULDAY==196 | JULDAY==349),{xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm
  + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(JULDAY),ylim=c(-20,70),xlab = "Time of Day
  (min)", ylab = "Soil Temperature (deg C)", auto.key=list(columns = 5), as.table = TRUE, type =
    "b", main=paste(minshade,"% shade"))})

## ---- fig.width=7, fig.height=6------------------------------------------
micro<-micro_global(loc = "Madison, Wisconsin, USA", runshade = 0, minshade = 50)

soil<-as.data.frame(micro$soil) # get the soil data
minshade<-micro$minshade # get the value for minimum shade
with(subset(soil,JULDAY==196 | JULDAY==349),{xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm
  + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(JULDAY),ylim=c(-20,70),xlab = "Time of Day
  (min)", ylab = "Soil Temperature (deg C)", auto.key=list(columns = 5), as.table = TRUE, type =
    "b", main=paste(minshade,"% shade"))})

## ---- fig.width=7, fig.height=6------------------------------------------
micro<-micro_global(loc = "Madison, Wisconsin, USA", runshade = 0, minshade = 0, slope = 45, aspect
  = 180)

soil<-as.data.frame(micro$soil) # get the soil data
minshade<-micro$minshade # get the value for minimum shade
with(subset(soil,JULDAY==196 | JULDAY==349),{xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm
  + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(JULDAY),ylim=c(-20,70),xlab = "Time of Day
  (min)", ylab = "Soil Temperature (deg C)", auto.key=list(columns = 5), as.table = TRUE, type =
    "b", main=paste(minshade,"% shade, 45 degree slope, 180 degrees aspect"))})

## ---- fig.width=7, fig.height=6------------------------------------------
micro<-micro_global(loc="Madison, Wisconsin, USA", runshade = 0, minshade = 0, hori = c(0, 0, 65,
  65, 65, 65, 65, 65, 65, 65, 65, 65, 0, 0, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65))

soil<-as.data.frame(micro$soil) # get the soil data
minshade<-micro$minshade # get the value for minimum shade
with(subset(soil,JULDAY==196 | JULDAY==349),{xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm
  + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(JULDAY),ylim=c(-20,70),xlab = "Time of Day
  (min)", ylab = "Soil Temperature (deg C)", auto.key=list(columns = 5), as.table = TRUE, type =
    "b", main=paste(minshade,"% shade, north-south gully"))})

## ---- fig.width=7, fig.height=6------------------------------------------
micro<-micro_global(loc="Madison, Wisconsin, USA", runshade = 0, minshade = 0, soiltype = 0)

soil<-as.data.frame(micro$soil) # get the soil data
minshade<-micro$minshade # get the value for minimum shade
with(subset(soil,JULDAY==196 | JULDAY==349),{xyplot(D0cm + D2.5cm + D5cm + D10cm + D15cm + D20cm
  + D30cm + D50cm + D100cm + D200cm ~ TIME | as.factor(JULDAY),ylim=c(-20,70),xlab = "Time of Day
  (min)", ylab = "Soil Temperature (deg C)", auto.key=list(columns = 5), as.table = TRUE, type =
    "b", main=paste(minshade,"% shade, rock substrate"))})

