## ----echo = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE, tidy.opts=list(width.cutoff=60), tidy=TRUE  
)

## ----message=FALSE, warnings=FALSE--------------------------------------------
library(NicheMapR)
library(knitr)

endo.out <- endoR(TA = 0)

treg <- as.data.frame(endo.out$treg)
morph <- as.data.frame(endo.out$morph)
enbal <- as.data.frame(endo.out$enbal)
masbal <- as.data.frame(endo.out$masbal)

## ----echo=FALSE, message=FALSE, warnings=FALSE--------------------------------
treg.lab <- c("core temperature (°C)", "lung temperature (°C)", "dorsal skin temperature (°C)", "ventral skin temperature (°C)", "dorsal fur/feather-air interface temperature (°C)", "ventral fur-air interface temperature (°C)", "current ratio between long and short axis due to postural change (-)", "breathing rate multiplier (-)", "part of the skin surface that is wet (%)", "thermal conductivity of flesh (W/mC)", "effective thermal conductivity of fur/feathers (W/mC)", "thermal conductivity of dorsal fur/feathers (W/mC)", "thermal conductivity of ventral fur/feathers (W/mC)", "thermal conductivity of compressed fur/feathers (W/mC)", "Q10 multiplier on metabolic rate (-)")

morph.lab <- c("total outer surface area (m2)", "total volume (m3)","characteristic dimension for convection (m)", "fat mass (kg)", "thickness of fat layer (m)", "flesh volume (m3)", "length (m)", "width (m)", "height (m)", "diameter, core to skin (m)", "diameter, core to fur/feathers (m)", "silhouette area (m2)", "silhouette area normal to sun's rays (m2)", "silhouette area parallel to sun's rays (m2)", "total skin area (m2)", "skin area available for evaporation (m2)", "area for convection (m2)", "area for conduction (m2)", "configuration factor to sky (-)", "configuration factor to ground (-)")

enbal.lab <- c("solar radiation absorbed (W)", "longwave (infra-red) radiation absorbed (W)", "metabolic heat generation (W)", "evaporation (W)", "longwave (infra-red) radiation lost (W)", "convection (W)", "conduction (W)", "energy balance (W)", "iterations required for a solution (-)", "was a solution found (0=no, 1=yes)")

masbal.lab <- c("breating rate (L/h)", "oxgyen consumption rate (L/h)", "respiratory water loss (g/h)", "cutaneous water loss (g/h)", "oxygen inhaled (mol/h)", "oxygen expelled (mol/h)", "nitrogen inhaled (mol/h)", "nitrogen expelled (mol/h)", "air inhaled (mol/h)", "air expelled (mol/h)")

kable(cbind(t(round(treg, 4)), treg.lab), col.names = c("value", "description and units"), caption = 'treg')
kable(cbind(t(round(morph, 4)), morph.lab), col.names = c("value", "description and units"), caption = 'morph')
kable(cbind(t(round(enbal, 4)), enbal.lab), col.names = c("value", "description and units"), caption = 'enbal')
kable(cbind(t(round(masbal, 4)), masbal.lab), col.names = c("value", "description and units"), caption = 'masbal')

## ----fig.width=10, fig.height=7, fig.show = "hold", message=FALSE, warnings=FALSE----
library(NicheMapR)
# environment
TAs <- seq(0, 50, 1) # air temperature (°C)
VEL <- 0.1 # wind speed (m/s)
vd <- WETAIR(rh = 30, db = 40)$vd # Weathers and Schoenbaechler had 16.7 mm Hg 
# above 40 °C = 30% RH at 40 °C
vd_sat <- WETAIR(rh = 100, db = TAs)$vd # Weathers and Schoenbaechler had 16.7 mm Hg 
# above 40 °C = 30% RH at 40 °C
exp_rh <- vd / vd_sat * 100
exp_rh[exp_rh > 100] <- 100
exp_rh[TAs < 30] <- 15
hum <- exp_rh

# core temperature
TC <- 38 # core temperature (°C)
TC_MAX <- 43 # maximum core temperature (°C)
TC_INC <- 0.25 # increment by which TC is elevated (°C)

# size and shape
AMASS <- 0.0337 # mass (kg)
SHAPE_B <- 1.1 # start off near to a sphere (-)
SHAPE_B_MAX <- 5 # maximum ratio of length to width/depth
UNCURL <- 0.1 # allows the animal to uncurl to SHAPE_B_MAX, the value being the increment 
# SHAPE_B is increased per iteration
SHAPE <- 4 # use ellipsoid geometry
SAMODE <- 1 # use bird surface area relations (2 is mammal, 0 is based on shape specified 
# in GEOM)

# feather properties
DHAIRD <- 30E-06 # hair diameter, dorsal (m)
DHAIRV <- 30E-06 # hair diameter, ventral (m)
LHAIRD <- 23.1E-03 # hair length, dorsal (m)
LHAIRV <- 22.7E-03 # hair length, ventral (m)
ZFURD <- 5.9E-03 # fur depth, dorsal (m)
ZFURV <- 5.7E-03 # fur depth, ventral (m)
RHOD <- 5000E+04 # hair density, dorsal (1/m2)
RHOV <- 5000E+04 # hair density, ventral (1/m2)
REFLD <- 0.248  # fur reflectivity dorsal (fractional, 0-1)
REFLV <- 0.351  # fur reflectivity ventral (fractional, 0-1)

# physiological responses
PCTWET <- 0.5 # base skin wetness (%)
PCTWET_MAX <- 5 # maximum skin wetness (%)
PCTWET_INC <- 0.25 # intervals by which skin wetness is increased (%)
Q10s <- rep(1, length(TAs))
Q10s[TAs >= TC_MAX] <- 2 # assuming Q10 effect kicks in only after air temp rises above TC_MAX
QBASAL <- 10 ^ (-1.461 + 0.669 * log10(AMASS * 1000)) # basal heat generation (W)
DELTAR <- 5 # offset between air temperature and breath (°C)
EXTREF <- 25 # O2 extraction efficiency (%)
PANT_INC <- 0.1 # turns on panting, the value being the increment by which the panting multiplier
# is increased up to the maximum value, PANT_MAX
PANT_MAX <- 15 # maximum panting rate - multiplier on air flow through the lungs above 
# that determined by metabolic rate
PANT_MULT <- 1 # multiplier on basal metabolic rate at maximum panting level

ptm <- proc.time() # start timing
endo.out <- lapply(1:length(TAs), function(x){endoR(TA = TAs[x], VEL = VEL, TC = TC, TC_MAX = TC_MAX, RH = hum[x], AMASS = AMASS, SHAPE = SHAPE, SHAPE_B = SHAPE_B, SHAPE_B_MAX = SHAPE_B_MAX, PCTWET = PCTWET, PCTWET_INC = PCTWET_INC, PCTWET_MAX = PCTWET_MAX, Q10 = Q10s[x], QBASAL = QBASAL, DELTAR = DELTAR, DHAIRD = DHAIRD, DHAIRV = DHAIRV, LHAIRD = LHAIRD, LHAIRV = LHAIRV, ZFURD = ZFURD, ZFURV = ZFURV, RHOD = RHOD, RHOV = RHOV, REFLD = REFLD, TC_INC = TC_INC, PANT_INC = PANT_INC, PANT_MAX = PANT_MAX, EXTREF = EXTREF, UNCURL = UNCURL, SAMODE = SAMODE, PANT_MULT = PANT_MULT)}) # run endoR across environments
proc.time() - ptm # stop timing

# extract the output
endo.out1 <- do.call("rbind", lapply(endo.out, data.frame))

# thermoregulation output
treg <- endo.out1[, grep(pattern = "treg", colnames(endo.out1))]
colnames(treg) <- gsub(colnames(treg), pattern = "treg.", replacement = "")

# morphometric output
morph <- endo.out1[, grep(pattern = "morph", colnames(endo.out1))]
colnames(morph) <- gsub(colnames(morph), pattern = "morph.", replacement = "")

# heat balance
enbal <- endo.out1[, grep(pattern = "enbal", colnames(endo.out1))]
colnames(enbal) <- gsub(colnames(enbal), pattern = "enbal.", replacement = "")

# mass aspects
masbal <- endo.out1[, grep(pattern = "masbal", colnames(endo.out1))]
colnames(masbal) <- gsub(colnames(masbal), pattern = "masbal.", replacement = "")

QGEN <- enbal$QGEN # metabolic rate (W)
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
plot(QGEN ~ TAs, type = 'l', ylab = 'metabolic rate, W', xlab = expression("air temperature, "*~degree*C*""), ylim = c(0.2, 1.2))
points(Weathers1976Fig1$Tair, Weathers1976Fig1$mlO2gh * 20.1 / 3600 * (AMASS * 1000), pch = 16, col = 2)
legend(x = 30, y = 1.2, legend = c("observed", "predicted"), col = c("red", "black"), lty = c(NA, 1), bty = "n", pch = c(16, NA))
plot(H2O ~ TAs, type = 'l', ylab = 'water loss, g/h', xlab = expression("air temperature, "*~degree*C*""), ylim = c(0, 1.5))
points(masbal$H2OResp_g ~ TAs, type = 'l', lty = 2)
points(masbal$H2OCut_g ~ TAs, type = 'l', lty = 2, col = 'blue')
legend(x = 3, y = 1.5, legend = c("total (obs)", "total (pred)", "respiratory", "cutaneous"), col = c("red", "black", "black", "blue"), lty = c(NA, 1, 2, 2), bty = "n", pch = c(16, NA, NA, NA))
points(Weathers1976Fig3$Tair, Weathers1976Fig3$mgH2Ogh * AMASS, pch = 16, col = 2)
plot(TFA_D ~ TAs, type = 'l', col = 'grey', ylab = expression("feather, skin and core temperature, "*~degree*C*""), xlab = expression("air temperature, "*~degree*C*""), ylim = c(10, 50))
points(TFA_V ~ TAs, type = 'l', col = 'grey', lty = 2)
points(TskinD ~ TAs, type = 'l', col = 'orange')
points(TskinV ~ TAs, type = 'l', col = 'orange', lty = 2)
points(TCs ~ TAs, type = 'l', col = 'red')
legend(x = 30, y = 33, legend = c("core (obs)", "core (pred)", "skin dorsal", "skin ventral", "feathers dorsal", "feathers ventral"), col = c("red", "red", "orange", "orange", "grey", "grey"), lty = c(NA, 1, 1, 2, 1, 2), bty = "n", pch = c(16, NA, NA, NA, NA, NA))
points(Weathers1976Fig2$Tair, Weathers1976Fig2$Tb, pch = 16, col = 2)
plot(masbal$AIR_L * 1000 / 60 ~ TAs, ylim = c(0, 250),  lty = 1, xlim = c(-5, 50), ylab = "ml air / min", xlab=expression("air temperature, "*~degree*C*""), type = 'l') 
legend(x = 0, y = 250, legend = c("observed", "predicted"), col = c("red", "black"), lty = c(NA, 1), bty = "n", pch = c(16, NA))
points(Weathers1976Fig5$breaths_min * (13.2 * AMASS ^ 1.08) * ((Weathers1976Fig5$Tair + 273.15) / 273.15) ~ Weathers1976Fig5$Tair, col = 2, pch = 16) # tidal volume allometry from Lasiewski, R. C., and W. A. Calder. 1971, correcting volume according to PV = nRT equation, where V_2 = T_2 * V_1 / T_2, and T_1 is at STP, so 0 °C

## ----fig.width=10, fig.height=7, fig.show = "hold", message=FALSE, warnings=FALSE----
library(NicheMapR)
loc <- c(131.05, -22.75)
Usrhyt <- 0.07
maxshade <- 100
micro <- micro_global(loc = loc, Usrhyt = Usrhyt, maxshade = maxshade)

metout <- as.data.frame(micro$metout) # unshaded above-ground conditions
soil <- as.data.frame(micro$soil) # unshaded below-ground soil temperatures
shadmet <- as.data.frame(micro$shadmet) # shaded above-ground conditions
shadsoil <- as.data.frame(micro$shadsoil) # shaded below-ground soil temperatures
dates <- micro$dates

TAs <- metout$TALOC # air temperatures at height of animal (°C)
TAREFs <- metout$TAREF # air temperatures at reference height (°C)
TSKYs <- metout$TSKYC # sky temperatures (°C)
TGRDs <- soil$D0cm # surface temperatures (°C)
VELs <- metout$VLOC # wind speeds at animal height (m/s)
RHs <- metout$RHLOC # relative humidity at animal height (%)
QSOLRs <- metout$SOLR # solar radiation (W/m2)
Zs <- metout$ZEN # zenith angle of the sun (degrees)
ELEV <- micro$elev # elevation (m)
ABSSB <- 1 - micro$REFL # substrate solar absorptivity (%)

# core temperature
TC <- 38 # core temperature (°C)
TC_MAX <- 43 # maximum core temperature (°C)
TC_INC <- 0.25 # increment by which TC is elevated (°C)

# size and shape
AMASS <- 0.0337 # mass (kg)
SHAPE <- 4 # use ellipsoid geometry
SHAPE_B <- 1.1 # start off near to a sphere (-)
SHAPE_B_MAX <- 5 # maximum ratio of length to width/depth
UNCURL <- 0.1 # allows the animal to uncurl to SHAPE_B_MAX, the value being the increment 
# SHAPE_B is increased per iteration

# feather properties
DHAIRD <- 30E-06 # hair diameter, dorsal (m)
DHAIRV <- 30E-06 # hair diameter, ventral (m)
LHAIRD <- 23.1E-03 # hair length, dorsal (m)
LHAIRV <- 22.7E-03 # hair length, ventral (m)
ZFURD <- 5.9E-03 # fur depth, dorsal (m)
ZFURV <- 5.7E-03 # fur depth, ventral (m)
RHOD <- 5000E+04 # hair density, dorsal (1/m2)
RHOV <- 5000E+04 # hair density, ventral (1/m2)
REFLD <- 0.248  # fur reflectivity dorsal (fractional, 0-1)
REFLV <- 0.351  # fur reflectivity ventral (fractional, 0-1)

# physiological responses
PCTWET <- 0.1 # base skin wetness (%)
PCTWET_MAX <- 0.1 # maximum skin wetness (%)
PCTWET_INC <- 0.25 # intervals by which skin wetness is increased (%)
Q10 <- 1 # Q10 effect of body temperature on metabolic rate (-)
QBASAL <- 10 ^ (-1.461 + 0.669 * log10(AMASS * 1000)) # basal heat generation (W)
DELTAR <- 5 # offset between air temperature and breath (°C)
EXTREF <- 25 # O2 extraction efficiency (%)
PANT_INC <- 0.1 # turns on panting, the value being the increment by which the panting 
# multiplier is increased up to the maximum value, PANT_MAX
PANT_MAX <- 15 # maximum panting rate - multiplier on air flow through the lungs above 
# that determined by metabolic rate
PANT_MULT <- 1 # multiplier on basal metabolic rate at maximum panting level

# run the model
ptm <- proc.time() # start timing
endo.out <- lapply(1:length(TAs), function(x){endoR(TA = TAs[x], TAREF = TAREFs[x], TSKY = TSKYs[x], TGRD = TGRDs[x], VEL = VELs[x], RH = RHs[x], QSOLR = QSOLRs[x], Z = Zs[x], ELEV = ELEV, ABSSB = ABSSB, TC = TC, TC_MAX = TC_MAX, AMASS = AMASS, SHAPE = SHAPE, SHAPE_B = SHAPE_B, SHAPE_B_MAX = SHAPE_B_MAX, PCTWET = PCTWET, PCTWET_INC = PCTWET_INC, Q10 = Q10, QBASAL = QBASAL, DELTAR = DELTAR, DHAIRD = DHAIRD, DHAIRV = DHAIRV, LHAIRD = LHAIRD, LHAIRV = LHAIRV, ZFURD = ZFURD, ZFURV = ZFURV, RHOD = RHOD, RHOV = RHOV, REFLD = REFLD, TC_INC = TC_INC, PANT_INC = PANT_INC, PANT_MAX = PANT_MAX, EXTREF = EXTREF, UNCURL = UNCURL, SAMODE = SAMODE, SHADE = 0, PANT_MULT = PANT_MULT)})
proc.time() - ptm # end timing

# extract the output
endo.out1 <- do.call("rbind", lapply(endo.out, data.frame))

# thermoregulation output
treg <- endo.out1[, grep(pattern = "treg", colnames(endo.out1))]
colnames(treg) <- gsub(colnames(treg), pattern = "treg.", replacement = "")

# morphometric output
morph <- endo.out1[, grep(pattern = "morph", colnames(endo.out1))]
colnames(morph) <- gsub(colnames(morph), pattern = "morph.", replacement = "")

# heat balance
enbal <- endo.out1[, grep(pattern = "enbal", colnames(endo.out1))]
colnames(enbal) <- gsub(colnames(enbal), pattern = "enbal.", replacement = "")

# mass aspects
masbal <- endo.out1[, grep(pattern = "masbal", colnames(endo.out1))]
colnames(masbal) <- gsub(colnames(masbal), pattern = "masbal.", replacement = "")

# extract variables for plotting
QGEN <- enbal$QGEN # metabolic rate (W)
H2O <- masbal$H2OResp_g + masbal$H2OCut_g # g/h water evaporated
TFA_D <- treg$TFA_D # dorsal fur surface temperature
TFA_V <- treg$TFA_V # ventral fur surface temperature
TskinD <- treg$TSKIN_D # dorsal skin temperature
TskinV <- treg$TSKIN_V # ventral skin temperature
TCs <- treg$TC # core temperature
SkinW <- treg$SKINWET # skin wetness (%)
Pant <- treg$PANT # panting multiplier (-)

# plot the results
par(mfrow = c(2, 2))
par(oma = c(2, 1, 2, 2) + 0.1)
par(mar = c(3, 3, 1.5, 1) + 0.1)
par(mgp = c(2, 1, 0))
plot(QGEN ~ dates, type = 'l', ylab = 'metabolic rate, W', xlab = 'time', ylim = c(0, QBASAL * 5))
plot(H2O ~ dates, type = 'l', ylab = 'water loss, g / h', xlab = 'time', ylim = c(0, 2))
points(masbal$H2OResp_g ~ dates, type = 'l', lty = 2)
points(masbal$H2OCut_g ~ dates, type = 'l', lty = 2, col = 'blue')
legend(x = 4, y = 2, legend = c("total", "respiratory", "cutaneous"), col = c("black", "black", "blue"), lty = c(1, 2, 2), bty = "n")  
plot(TFA_D ~ dates, type = 'l', col = 'grey', ylab = 'temperature, °C', xlab = 'time', ylim = c(0, 60))
points(TFA_V ~ dates, type = 'l', col = 'grey', lty = 2)
points(TskinD ~ dates, type = 'l', col = 'orange')
points(TskinV ~ dates, type = 'l', col = 'orange', lty = 2)
points(TCs ~ dates, type = 'l', col = 'red')
legend(x = 1, y = 65, legend = c("core", "skin dorsal", "skin ventral", "feathers dorsal", "feathers ventral"), col = c("red", "orange", "orange", "grey", "grey"), lty = c(1, 1, 2, 1, 2), bty = "n", ncol = 2)  
plot((H2O / (AMASS * 1000)) * 100 ~ dates, type = 'l', ylab = 'H2O loss, % body mass / h', xlab = 'time', ylim = c(0, 5))

