## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE, tidy.opts=list(width.cutoff=60), tidy=TRUE  
)

## ------------------------------------------------------------------------
library(NicheMapR)
library(knitr)

endo.out <- endoR(TA = 0)

treg <- as.data.frame(endo.out$treg)
morph <- as.data.frame(endo.out$morph)
enbal <- as.data.frame(endo.out$enbal)
masbal <- as.data.frame(endo.out$masbal)

kable(treg[, 1:10])
kable(treg[, 11:15])
kable(morph[, 1:10])
kable(morph[, 11:20])
kable(enbal)
kable(masbal)

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
library(NicheMapR)
# environment
TAs <- seq(0, 50, 1) # air temperature (deg C)
VEL <- 0.002 # wind speed (m/s)
vd <- WETAIR(rh = 30, db = 40)$vd # Weather and Schoenbaechler had 16.7 mm Hg above 40 deg C = 30% RH at 40 deg C
vd_sat <- WETAIR(rh = 100, db = TAs)$vd # Weather and Schoenbaechler had 16.7 mm Hg above 40 deg C = 30% RH at 40 deg C
exp_rh <- vd / vd_sat * 100
exp_rh[exp_rh > 100] <- 100
exp_rh[TAs < 30] <- 15
hum <- exp_rh#rep(humidity,96)

# core temperature
TC <- 38 # core temperature (deg C)
TCMAX <- 43 # maximum core temperature (deg C)
RAISETC <- 0.25 # increment by which TC is elevated (deg C)

# size and shape
AMASS <- 0.0337 # mass (kg)
SHAPE_B_REF <- 1.1 # start off near to a sphere (-)
SHAPE_B_MAX <- 5 # maximum ratio of length to width/depth
UNCURL <- 0.1 # allows the animal to uncurl to SHAPE_B_MAX, the value being the increment SHAPE_B is increased per iteration
SHAPE <- 4 # use ellipsoid geometry
SAMODE <- 1 # use bird surface area relations (2 is mammal, 0 is based on shape specified in GEOM)

# feather properties
DHAIRD <- 30E-06 # hair diameter, dorsal (m)
DHAIRV <- 30E-06 # hair diameter, ventral (m)
LHAIRD <- 23.1E-03 # hair length, dorsal (m)
LHAIRV <- 22.7E-03 # hair length, ventral (m)
ZFURD <- 5.8E-03 # fur depth, dorsal (m)
ZFURV <- 5.6E-03 # fur depth, ventral (m)
RHOD <- 8000E+04 # hair density, dorsal (1/m2)
RHOV <- 8000E+04 # hair density, ventral (1/m2)
REFLD <- 0.248  # fur reflectivity dorsal (fractional, 0-1)
REFLV <- 0.351  # fur reflectivity ventral (fractional, 0-1)

# physiological responses
SKINW <- 0.5 # base skin wetness (%)
MXWET <- 5 # maximum skin wetness (%)
SWEAT <- 0.25 # intervals by which skin wetness is increased (%)
Q10s <- rep(1, length(TAs))
Q10s[TAs >= TCMAX] <- 2 # assuming Q10 effect kicks in only after air temp rises above TCMAX
QBASAL <- 10 ^ (-1.461 + 0.669 * log10(AMASS * 1000)) # basal heat generation (W)
DELTAR <- 5 # offset between air temeprature and breath (°C)
EXTREF <- 25 # O2 extraction efficiency (%)
PANTING <- 0.1 # turns on panting, the value being the increment by which the panting multiplier is increased up to the maximum value, PANTMAX
PANTMAX <- 15# maximum panting rate - multiplier on air flow through the lungs above that determined by metabolic rate

ptm <- proc.time() # start timing
endo.out <- lapply(1:length(TAs), function(x){endoR(TA = TAs[x], VEL = VEL, TC = TC, TCMAX = TCMAX, RH = hum[x], AMASS = AMASS, SHAPE = SHAPE, SHAPE_B_REF = SHAPE_B_REF, SHAPE_B_MAX = SHAPE_B_MAX, SKINW = SKINW, SWEAT = SWEAT, MXWET = MXWET, Q10 = Q10s[x], QBASAL = QBASAL, DELTAR = DELTAR, DHAIRD = DHAIRD, DHAIRV = DHAIRV, LHAIRD = LHAIRD, LHAIRV = LHAIRV, ZFURD = ZFURD, ZFURV = ZFURV, RHOD = RHOD, RHOV = RHOV, REFLD = REFLD, RAISETC = RAISETC, PANTING = PANTING, PANTMAX = PANTMAX, EXTREF = EXTREF, UNCURL = UNCURL, SAMODE = SAMODE)}) # run endoR across environments
proc.time() - ptm # stop timing

endo.out1 <- do.call("rbind", lapply(endo.out, data.frame)) # turn results into data frame
treg <- endo.out1[, grep(pattern = "treg", colnames(endo.out1))]
colnames(treg) <- gsub(colnames(treg), pattern = "treg.", replacement = "")
morph <- endo.out1[, grep(pattern = "morph", colnames(endo.out1))]
colnames(morph) <- gsub(colnames(morph), pattern = "morph.", replacement = "")
enbal <- endo.out1[, grep(pattern = "enbal", colnames(endo.out1))]
colnames(enbal) <- gsub(colnames(enbal), pattern = "enbal.", replacement = "")
masbal <- endo.out1[, grep(pattern = "masbal", colnames(endo.out1))]
colnames(masbal) <- gsub(colnames(masbal), pattern = "masbal.", replacement = "")

QGEN <- enbal$QMET # metabolic rate (W)
H2O <- masbal$H2OResp_g + masbal$H2OCut_g # g/h water evaporated
TFA_D <- treg$TFA_D # dorsal fur surface temperature
TFA_V <- treg$TFA_V # ventral fur surface temperature
TskinD <- treg$TSKIN_D # dorsal skin temperature
TskinV <- treg$TSKIN_V # ventral skin temperature
TCs <- treg$TC # core temperature

par(mfrow = c(2, 2))
par(oma = c(2, 1, 2, 2) + 0.1)
par(mar = c(3, 3, 1.5, 1) + 0.1)
par(mgp = c(2, 1, 0))
plot(QGEN ~ TAs, type = 'l', ylab = 'metabolic rate, W', xlab = 'air temperature, deg C', ylim = c(0.2, 1.2))
points(Weathers1976Fig1$Tair, Weathers1976Fig1$mlO2gh * 20.1 / 3600 * (AMASS * 1000), pch = 16, col = 2)
plot(H2O ~ TAs, type = 'l', ylab = 'water loss, g/h', xlab = 'air temperature, deg C', ylim = c(0, 1.5))
points(Weathers1976Fig3$Tair, Weathers1976Fig3$mgH2Ogh * AMASS, pch = 16, col = 2)
plot(TFA_D ~ TAs, type = 'l', col = 'grey', ylab = 'fur, skin and core temperature, deg C', xlab = 'air temperature, deg C', ylim = c(10, 50))
points(TFA_V ~ TAs, type = 'l', col = 'grey', lty = 2)
points(TskinD ~ TAs, type = 'l', col = 'orange')
points(TskinV ~ TAs, type = 'l', col = 'orange', lty = 2)
points(TCs ~ TAs, type = 'l', col = 'red')
points(Weathers1976Fig2$Tair, Weathers1976Fig2$Tb, pch = 16, col = 2)
plot(masbal$AIR_L * 1000 / 60 ~ TAs, ylim=c(0,250),  lty = 1, xlim=c(-5,50), main = "minute volume", ylab = "ml / min", xlab=paste("air temperature (deg C)"), type = 'l') 
points(Weathers1976Fig5$breaths_min * (13.2 * AMASS ^ 1.08) * ((Weathers1976Fig5$Tair + 273.15) / 273.15) ~ Weathers1976Fig5$Tair, col='red', pch = 16) # tidal volume allometry from Lasiewski, R. C., and W. A. Calder. 1971, correcting volume according to PV = nRT equation, where V_2 = T_2 * V_1 / T_2, and T_1 is at STP, so 0 deg C

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
library(NicheMapR)
loc <- c(131.05, -22.75)
Usrhyt <- 0.07
maxshade <- 100
micro <- micro_global(loc = loc, Usrhyt = Usrhyt, maxshade = maxshade)

metout <- as.data.frame(micro$metout)
soil <- as.data.frame(micro$soil)
shadmet <- as.data.frame(micro$shadmet)
shadsoil <- as.data.frame(micro$shadsoil)
dates <- micro$dates

TAs <- metout$TALOC
TAREFs <- metout$TAREF
TSKYs <- metout$TSKYC
TGRDs <- soil$D0cm
VELs <- metout$VLOC
RHs <- metout$RHLOC
QSOLRs <- metout$SOLR 
Zs <- metout$ZEN
ELEV <- micro$elev
ABSSB <- 1-micro$REFL

# core temperature
TC <- 38 # core temperature (deg C)
TCMAX <- 43 # maximum core temperature (deg C)
RAISETC <- 0.25 # increment by which TC is elevated (deg C)

# size and shape
AMASS <- 0.0337 # mass (kg)
SHAPE <- 4 # use ellipsoid geometry
SHAPE_B_REF <- 1.1 # start off near to a sphere (-)
SHAPE_B_MAX <- 5 # maximum ratio of length to width/depth
UNCURL <- 0.1 # allows the animal to uncurl to SHAPE_B_MAX, the value being the increment SHAPE_B is increased per iteration

# feather properties
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
SKINW <- 0.1 # base skin wetness (%)
MXWET <- 0.1 # maximum skin wetness (%)
SWEAT <- 0.25 # intervals by which skin wetness is increased (%)
Q10 <- 1 # Q10 effect of body temperature on metabolic rate (-)
QBASAL <- 10 ^ (-1.461 + 0.669 * log10(AMASS * 1000)) # basal heat generation (W)
DELTAR <- 5 # offset between air temeprature and breath (°C)
EXTREF <- 25 # O2 extraction efficiency (%)
PANTING <- 0.1 # turns on panting, the value being the increment by which the panting multiplier is increased up to the maximum value, PANTMAX
PANTMAX <- 15# maximum panting rate - multiplier on air flow through the lungs above that determined by metabolic rate

ptm <- proc.time() # start timing
endo.out <- lapply(1:length(TAs), function(x){endoR(TA = TAs[x], TAREF = TAREFs[x], TSKY = TSKYs[x], TGRD = TGRDs[x], VEL = VELs[x], RH = RHs[x], QSOLR = QSOLRs[x], Z = Zs[x], ELEV = ELEV, ABSSB = ABSSB, TC = TC, TCMAX = TCMAX, AMASS = AMASS, SHAPE = SHAPE, SHAPE_B_REF = SHAPE_B_REF, SHAPE_B_MAX = SHAPE_B_MAX, SKINW = SKINW, SWEAT = SWEAT, Q10 = Q10, QBASAL = QBASAL, DELTAR = DELTAR, DHAIRD = DHAIRD, DHAIRV = DHAIRV, LHAIRD = LHAIRD, LHAIRV = LHAIRV, ZFURD = ZFURD, ZFURV = ZFURV, RHOD = RHOD, RHOV = RHOV, REFLD = REFLD, RAISETC = RAISETC, PANTING = PANTING, PANTMAX = PANTMAX, EXTREF = EXTREF, UNCURL = UNCURL, SAMODE = SAMODE, SHADE = 0)})
proc.time() - ptm
endo.out1 <- do.call("rbind", lapply(endo.out, data.frame))

treg <- endo.out1[, grep(pattern = "treg", colnames(endo.out1))]
colnames(treg) <- gsub(colnames(treg), pattern = "treg.", replacement = "")
morph <- endo.out1[, grep(pattern = "morph", colnames(endo.out1))]
colnames(morph) <- gsub(colnames(morph), pattern = "morph.", replacement = "")
enbal <- endo.out1[, grep(pattern = "enbal", colnames(endo.out1))]
colnames(enbal) <- gsub(colnames(enbal), pattern = "enbal.", replacement = "")
masbal <- endo.out1[, grep(pattern = "masbal", colnames(endo.out1))]
colnames(masbal) <- gsub(colnames(masbal), pattern = "masbal.", replacement = "")

QGEN <- enbal$QMET # metabolic rate (W)
H2O <- masbal$H2OResp_g + masbal$H2OCut_g # g/h water evaporated
TFA_D <- treg$TFA_D # dorsal fur surface temperature
TFA_V <- treg$TFA_V # ventral fur surface temperature
TskinD <- treg$TSKIN_D # dorsal skin temperature
TskinV <- treg$TSKIN_V # ventral skin temperature
TCs <- treg$TC # core temperature
SkinW <- treg$SKINWET # skin wetness (%)
Pant <- treg$PANT # panting multiplier (-)

par(mfrow = c(2, 2))
par(oma = c(2, 1, 2, 2) + 0.1)
par(mar = c(3, 3, 1.5, 1) + 0.1)
par(mgp = c(2, 1, 0))
plot(QGEN ~ dates, type = 'l', ylab = 'metabolic rate, W', xlab = 'time', ylim = c(0, QBASAL * 5))
plot(H2O ~ dates, type = 'l', ylab = 'water loss, g / h', xlab = 'time', ylim = c(0, 2))
plot(TFA_D ~ dates, type = 'l', col = 'grey', ylab = 'fur, skin and core temperature, deg C', xlab = 'time', ylim = c(0, 60))
points(TFA_V ~ dates, type = 'l', col = 'grey', lty = 2)
points(TskinD ~ dates, type = 'l', col = 'orange')
points(TskinV ~ dates, type = 'l', col = 'orange', lty = 2)
points(TCs ~ dates, type = 'l', col = 'red')
plot((H2O / (AMASS * 1000)) * 1000 ~ dates, type = 'l', ylab = 'H2O loss, % body mass / h', xlab = 'time', ylim = c(0, 50))

