## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE, tidy.opts=list(width.cutoff=60), tidy=TRUE  
)

## ------------------------------------------------------------------------
library(NicheMapR)
library(knitr)
endo.out <- endoR(TA = 0)
kable(endo.out[, 1:9])
kable(endo.out[, 10:19])
kable(endo.out[, 20:29])
kable(endo.out[, 30:39])
kable(endo.out[, 40:49])
kable(endo.out[, 50:57])

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
TAs <- seq(0,50) # air temperature (deg C)
TC <- 37 # core temperature (deg C)
AMASS <- 1 # mass (kg)
ZFURD <- 2E-03 # fur depth, dorsal (m)
ZFURV <- 2E-03 # fur depth, dorsal (m)
GMREF <- 4
SKINW <- 0.05
SWEAT <- 0.25
QBASAL <- (70 * AMASS ^ 0.75) * (4.185 / (24 * 3.6)) # basal heat generation (W)

ptm <- proc.time()
endo.out <- lapply(1:length(TAs), function(x){endoR(TA = TAs[x], AMASS = AMASS, ZFURD = ZFURD, ZFURV = ZFURV, GMREF = GMREF, SKINW = SKINW, SWEAT = SWEAT, QBASAL = QBASAL)})
proc.time() - ptm
endo.out <- do.call("rbind", lapply(endo.out, data.frame))

QGEN <- endo.out$RESPGEN # metabolic rate (W)
H2O <- endo.out$EVAP.G.H # g/h water evaporated
TFA_D <- endo.out$TFA_D # dorsal fur surface temperature
TFA_V <- endo.out$TFA_V # ventral fur surface temperature
TskinD <- endo.out$TSKIN_D # dorsal skin temperature
TskinV <- endo.out$TSKIN_V # ventral skin temperature
TC <- endo.out$TC # core temperature
plot(QGEN ~ TAs, type = 'l', ylab = 'metabolic rate, W', xlab = 'air temperature, deg C')
plot(H2O ~ TAs, type = 'l', ylab = 'water loss, g/h', xlab = 'air temperature, deg C')
plot(TFA_D ~ TAs, type = 'l', col = 'grey', ylab = 'fur, skin and core temperature, deg C', xlab = 'air temperature, deg C')
points(TFA_V ~ TAs, type = 'l', col = 'grey', lty = 2)
points(TskinD ~ TAs, type = 'l', col = 'orange')
points(TskinV ~ TAs, type = 'l', col = 'orange', lty = 2)
points(TC ~ TAs, type = 'l', col = 'red')

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
micro <- micro_global(loc = c(145, -37), runshade = 0, Usrhyt = 0.01)

metout <- as.data.frame(micro$metout)
soil <- as.data.frame(micro$soil)
days<-rep(seq(1,12),24)
days<-days[order(days)]
dates<-days+metout$TIME/60/24-1 # dates for hourly output

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


endo.out <- lapply(1:length(TAs), function(x){endoR(TA = TAs[x], TAREF = TAREFs[x], TSKY = TSKYs[x], 
  TGRD = TGRDs[x], VEL = VELs[x], RH = RHs[x], QSOLR = QSOLRs[x], Z = Zs[x], ELEV = ELEV, ABSSB = ABSSB, AMASS = AMASS, ZFURD = ZFURD, ZFURV = ZFURV, GMREF = GMREF, SKINW = SKINW, SWEAT = SWEAT)})
endo.out <- do.call("rbind", lapply(endo.out, data.frame))

QGEN <- endo.out$RESPGEN # metabolic rate (W)
H2O <- endo.out$GEVAP * 3600 # g/h water evaporated
TFA_D <- endo.out$TFA_D # dorsal fur surface temperature
TFA_V <- endo.out$TFA_V # ventral fur surface temperature
TskinD <- endo.out$TSKIN_D # dorsal skin temperature
TskinV <- endo.out$TSKIN_V # ventral skin temperature
TC <- endo.out$TC # core temperature

plot(QGEN ~ dates, type = 'l', ylab = 'metabolic rate, W', xlab = 'time')
plot(H2O ~ dates, type = 'l', ylab = 'water loss, g/h', xlab = 'time')
plot(TFA_D ~ dates, type = 'l', col = 'grey', ylab = 'fur, skin and core temperature, deg C', xlab = 'time')
points(TFA_V ~ dates, type = 'l', col = 'grey', lty = 2)
points(TskinD ~ dates, type = 'l', col = 'orange')
points(TskinV ~ dates, type = 'l', col = 'orange', lty = 2)
points(TC ~ dates, type = 'l', col = 'red')

