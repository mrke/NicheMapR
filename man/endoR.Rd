% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/endoR_solvendo.R
\encoding{UTF-8}
\name{endoR}
\alias{endoR}
\title{endoR - the endotherm model of NicheMapR}
\usage{
endoR(AMASS = 1, NGEOM = 4, GMREF = 3, FURTHRMK = 0, ZFURD = 2E-03, ZFURV = 2E-03, TC = 37, TCMAX = 45, TA = 20, TGRD = TA, TSKY = TA, VEL = 0.1, RH = 5, QSOLR = 0, Z = 20, SHADE = 0, NITESHAD = 0,...)
}
\arguments{
\item{TA}{= 20, air temperature at local height (°C)}

\item{TGRD}{= TA, ground temperature (°C)}

\item{TSKY}{= TA, sky temperature (°C)}

\item{VEL}{= 0.1, wind speed (m/s)}

\item{RH}{= 5, relative humidity (\%)}

\item{QSOLR}{= 0, solar radiation, horizontal plane (W/m2)}

\item{Z}{= 20, zenith angle of sun (degrees from overhead)}

\item{SHADE}{= 0, shade level (\%)}

\item{AMASS}{= 1, # kg}

\item{NGEOM}{= 4, # cylinder (ngeom = 1), sphere (ngeom = 2) and ellipsoid (ngeom = 4). If a truncated cone (5) or ellipsoidal cylinder (3), we will use the cylinder equations (ngeom=1).}

\item{GMREF}{= 3, # initial ratio between long and short axis (-)}

\item{FURTHRMK}{= 0, # user-specified fur thermal conductivity (W/mK), not used if 0}

\item{ZFURD}{= 2E-03, # fur depth, dorsal (m)}

\item{ZFURV}{= 2E-03, # fur depth, ventral (m)}

\item{TC}{= 37, # core temperature (°C)}

\item{TCMAX}{= 45, # maximum core temperature (°C)}
}
\description{
This model implements postural and physiological thermoregulation under a given
environmental scenario for an organism of a specified shape and no extra body
parts. In this function the sequence of thermoregulatory events in the face of
heat stress is to first uncurl, second change flesh conductivity, third raise
core temprature, fourth pant and fifth sweat.
}
\examples{
library(NicheMapR)
# environment
TAs <- seq(0, 50, 2) # air temperatures (Â°C)
VEL <- 0.002 # wind speed (m/s)
RH <- 10 # relative humidity (\%)
SOL <- 100 # solar radiation (W/m2)

# core temperature
TC <- 38 # core temperature (deg C)
TCMAX <- 43 # maximum core temperature (Â°C)
RAISETC <- 0.25 # increment by which TC is elevated (Â°C)

# size and shape
AMASS <- 0.0337 # mass (kg)
GMREF <- 1.1 # start off near to a sphere (-)
GMULTMAX <- 5 # maximum ratio of length to width/depth

# fur/feather properties
DHAIRD = 30E-06 # hair diameter, dorsal (m)
DHAIRV = 30E-06 # hair diameter, ventral (m)
LHAIRD = 23.1E-03 # hair length, dorsal (m)
LHAIRV = 22.7E-03 # hair length, ventral (m)
ZFURD = 5.8E-03 # fur depth, dorsal (m)
ZFURV = 5.6E-03 # fur depth, ventral (m)
RHOD = 8000E+04 # hair density, dorsal (1/m2)
RHOV = 8000E+04 # hair density, ventral (1/m2)
REFLD = 0.248  # fur reflectivity dorsal (fractional, 0-1)
REFLV = 0.351  # fur reflectivity ventral (fractional, 0-1)

# physiological responses
SKINW <- 0.1 # base skin wetness (\%)
MXWET <- 20 # maximum skin wetness (\%)
SWEAT <- 0.25 # intervals by which skin wetness is increased (\%)
Q10 <- 2 # A10 effect of body temperature on metabolic rate
QBASAL <- 10 ^ (-1.461 + 0.669 * log10(AMASS * 1000)) # basal heat generation (W) (bird formula from McKechnie and Wolf 2004 Phys. & Biochem. Zool. 77:502-521)
DELTAR <- 5 # offset between air temeprature and breath (Â°C)
EXTREF <- 15 # O2 extraction efficiency (\%)
PANTING <- 0.1 # turns on panting, the value being the increment by which the panting multiplier is increased up to the maximum value, PANTMAX
PANTMAX <- 3# maximum panting rate - multiplier on air flow through the lungs above that determined by metabolic rate

ptm <- proc.time() # start timing
endo.out <- lapply(1:length(TAs), function(x){endoR(TA = TAs[x], SOL = SOL, VEL = VEL, TC = TC, TCMAX = TCMAX, RH = RH, AMASS = AMASS, GMREF = GMREF, GMULTMAX = GMULTMAX, SKINW = SKINW, SWEAT = SWEAT, MXWET = MXWET, Q10 = Q10, QBASAL = QBASAL, DELTAR = DELTAR, DHAIRD = DHAIRD, DHAIRV = DHAIRV, LHAIRD = LHAIRD, LHAIRV = LHAIRV, ZFURD = ZFURD, ZFURV = ZFURV, RHOD = RHOD, RHOV = RHOV, REFLD = REFLD, RAISETC = RAISETC, PANTING = PANTING, PANTMAX = PANTMAX, EXTREF = EXTREF)}) # run endoR across environments
proc.time() - ptm # stop timing

endo.out <- do.call("rbind", lapply(endo.out, data.frame)) # turn results into data frame

QGEN <- endo.out$RESPGEN # metabolic rate (W)
H2O <- endo.out$EVAP.G.H # g/h water evaporated
TFA_D <- endo.out$TFA_D # dorsal fur surface temperature
TFA_V <- endo.out$TFA_V # ventral fur surface temperature
TskinD <- endo.out$TSKIN_D # dorsal skin temperature
TskinV <- endo.out$TSKIN_V # ventral skin temperature
TCs <- endo.out$TC # core temperature

par(mfrow = c(2, 2))
par(oma = c(2, 1, 2, 2) + 0.1)
par(mar = c(3, 3, 1.5, 1) + 0.1)
par(mgp = c(2, 1, 0))
plot(QGEN ~ TAs, type = 'l', ylab = 'metabolic rate, W', xlab = 'air temperature, deg C', ylim = c(0.2, 1.2))
plot(H2O ~ TAs, type = 'l', ylab = 'water loss, g/h', xlab = 'air temperature, deg C', ylim = c(0, 1.5))
plot(TFA_D ~ TAs, type = 'l', col = 'grey', ylab = 'fur, skin and core temperature, deg C', xlab = 'air temperature, deg C', ylim = c(10, 50))
points(TFA_V ~ TAs, type = 'l', col = 'grey', lty = 2)
points(TskinD ~ TAs, type = 'l', col = 'orange')
points(TskinV ~ TAs, type = 'l', col = 'orange', lty = 2)
points(TCs ~ TAs, type = 'l', col = 'red')
plot(endo.out$AIRVOL * 1e6 / 60 ~ TAs, ylim=c(0,250),  lty = 1, xlim=c(-5,50), main = "minute volume", ylab = "ml / min", xlab=paste("air temperature (deg C)"), type = 'l')
}