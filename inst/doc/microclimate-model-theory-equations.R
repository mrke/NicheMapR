## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ----fig.width=8, fig.height=3,echo=FALSE, results='asis', fig.align='center', fig.cap='Stucture of the NicheMapR microclimate model Fortran program. Solid boxes are subroutines and dashed boxes are functions.', fig.pos='!h', message=FALSE, warnings=FALSE----
library(png)
library(grid)
img <- readPNG("Fig1.png")
 grid.raster(img)

## ----echo=FALSE, message=FALSE, warnings=FALSE, results='asis', fig.width=7, fig.height=6, fig.align='center', fig.cap='Solar spectral irradiance as a function of wavelength'----
optical<-read.csv("OPTICAL.csv")
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
par(mar = c(5,5,4,2) + 0.1)
ylab=expression(paste("solar spectral irradiance,  W m"^{-2}))
plot(optical[,2]~optical[,1],type='l',ylab=expression(paste("solar spectral irradiance,  W m"^{-2})),xlab="wavelength, nm")

## ----echo=FALSE, message=FALSE, warnings=FALSE, results='asis', fig.width=7, fig.height=6, fig.align='center', fig.cap='Wavelength-specific optical thickness for different atmospheric components'----
par(mfrow = c(2,2)) # set up for 2 plots in 2 columns
plot(optical[,3]~optical[,1],type='l',ylab="optical thickness, cm",xlab="wavelength, nm",main="a) molecular")
plot(optical[,5]~optical[,1],type='l',ylab="optical thickness, cm",xlab="wavelength, nm",main="b) ozone")
plot(optical[,6]~optical[,1],type='l',ylab="optical thickness, cm",xlab="wavelength, nm",main="c) water")
plot(optical[,4]~optical[,1],type='l',ylab="optical thickness, cm",xlab="wavelength, nm",main="d) aerosol")

## ----echo=FALSE, message=FALSE, warnings=FALSE, results='asis', fig.width=7, fig.height=10, fig.align='center', fig.cap='Seasonal, latitudinal and altitudinal correction factors for atmospheric ozone', fig.subcap=c('seasonal and latitudinal correction factors','altitudinal correction factors')----
altfct<-read.csv("ALTFCT.csv")
oz<-read.csv("OZ.csv")
par(mfrow = c(2,1)) # set up for 2 plots in 1 columns
for(i in 2:13){
  if(i==2){
    plot(oz[,i]~oz[,1],type='b',ylab="correction factor",xlab="latitude",main="a) latitude and season", col=i-1,pch=i-1,ylim=c(0.2,0.5))
  }else{
    points(oz[,i]~oz[,1],type='b',ylab="correction factor",xlab="latitude",pch=i-1,col=i-1)
  }
}
legend(-95,.5, c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'),pch=seq(1,12), col=seq(1,12),bty="n",horiz=TRUE,cex=1)

plot(altfct[,2]~altfct[,1],type='b',ylab="correction factor",xlab="altitude, km", main="b) altitude",col=1,pch=1,ylim=c(0,1))
points(altfct[,3]~altfct[,1],type='b',col=1,pch=2,ylim=c(0,1))
points(altfct[,4]~altfct[,1],type='b',col=1,pch=3,ylim=c(0,1))
legend(16,1.05, c('molecular', 'aerosol', 'ozone'),pch=seq(1,3), col=1,bty="n",horiz=FALSE,cex=1)

