## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE, tidy.opts=list(width.cutoff=60), tidy=TRUE  
)

## ---- warning=FALSE-----------------------------------------------------------
library(NicheMapR)
library(knitr)
library(plotrix)

## -----------------------------------------------------------------------------
# environmental input
TA <- 20 # air temperature, for calculation of conductivity of air (°C)

# shape input
PVEN <- 0.3 # maximum fraction of surface area that is ventral (fractional, 0-1)

# fur properties
DHAIRD <- 30E-06 # hair diameter, dorsal (m)
DHAIRV <- 30E-06 # hair diameter, ventral (m)
LHAIRD <- 23.9E-03 # hair length, dorsal (m)  
LHAIRV <- 23.9E-03 # hair length, ventral (m)  
ZFURD <- 9E-03 # fur depth, dorsal (m)
ZFURV <- 9E-03 # fur depth, ventral (m)
ZFURCOMP <- ZFURV # depth of compressed fur (for conduction) (m)
RHOD <- 3968E+04 # hair density, dorsal (1/m2) 
RHOV <- 2781E+04 # hair density, ventral (1/m2)
REFLD <- 0.301  # fur reflectivity dorsal (fractional, 0-1) 
REFLV <- 0.301  # fur reflectivity ventral (fractional, 0-1)
KHAIR <- 0.209  # hair thermal conductivity (W/m°C)

# call the subroutine
IRPROP.out <- IRPROP(TA, DHAIRD, DHAIRV, LHAIRD, LHAIRV, ZFURD, ZFURV, RHOD, RHOV, REFLD, REFLV, ZFURCOMP, PVEN, KHAIR)

# output
KEFARA <- IRPROP.out[1:3] # effective thermal conductivity of fur array, mean, dorsal, ventral (W/mK)
BETARA <- IRPROP.out[4:6] # term involved in computing optical thickness (1/m)
B1ARA <- IRPROP.out[7:9] # optical thickness array, mean, dorsal, ventral (-)
DHAR <- IRPROP.out[10:12] # fur diameter array, mean, dorsal, ventral (m)
LHAR <- IRPROP.out[13:15] # fur length array, mean, dorsal, ventral (m)
RHOAR <- IRPROP.out[16:18] # fur density array, mean, dorsal, ventral (fibers/m2)  
ZZFUR <- IRPROP.out[19:21] # fur depth array, mean, dorsal, ventral (m)  
REFLFR <- IRPROP.out[22:24] # fur reflectivity array, mean, dorsal, ventral (fractional, 0-1)
FURTST <- IRPROP.out[25] # test of presence of fur (length x diameter x density x depth) (-)
KFURCMPRS <- IRPROP.out[26] # effective thermal conductivity of compressed ventral fur (W/mK)

IRPROP.lab <- c("KEFARA mean", "KEFARA dorsal", "KEFARA ventral", "BETARA mean", "BETARA dorsal", "BETARA ventral", "B1ARA mean", "B1ARA dorsal", "B1ARA ventral", "DHAR mean", "DHAR dorsal", "DHAR ventral", "LHAR mean", "LHAR dorsal", "LHAR ventral", "RHOAR mean", "RHOAR dorsal", "RHOAR ventral", "ZZFUR mean", "ZZFUR dorsal", "ZZFUR ventral", "REFLFR mean", "REFLFR dorsal", "REFLFR ventral", "FURTST", "KFURCMPRS")
kable(cbind(IRPROP.lab, IRPROP.out[1:26]))

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
DHAIRs <- seq(0, 150, 2) # hair diameters (micrometers)
KEFARAs <- NULL
for(i in 1:length(DHAIRs)){
  KEFARAs[i] <- IRPROP(TA, DHAIRs[i] * 1E-06, DHAIRs[i] * 1E-06, LHAIRD, LHAIRV, ZFURD, ZFURV, RHOD, RHOV, REFLD, REFLV, ZFURCOMP, PVEN, KHAIR)[1]
}
plot(KEFARAs ~ DHAIRs, type = 'p', pch = 16, ylab = 'effective fur conductivity, W K-1 m-1', xlab = 'hair diameter, um')

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
RHOs <- seq(0, 50000, 500) # hair densities (1/cm2)
KEFARAs <- NULL
for(i in 1:length(RHOs)){
  KEFARAs[i] <- IRPROP(TA, DHAIRD, DHAIRV, LHAIRD, LHAIRV, ZFURD, ZFURV, RHOs[i] * 1E+04, RHOs[i] * 1E+04, REFLD, REFLV, ZFURCOMP, PVEN, KHAIR)[1]
}
plot(KEFARAs ~ RHOs, type = 'p', pch = 16, ylab = 'effective fur conductivity, W K-1 m-1', xlab = 'hair density, 1/cm2')

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
ZFURs <- seq(0, 50, 1) # hair depths (mm)
KEFARAs <- NULL
for(i in 1:length(ZFURs)){
  KEFARAs[i] <- IRPROP(TA, DHAIRD, DHAIRV, LHAIRD, LHAIRV, ZFURs[i] * 1E-03, ZFURs[i] * 1E-03, RHOD, RHOV, REFLD, REFLV, ZFURCOMP, PVEN, KHAIR)[1]
}
plot(KEFARAs ~ ZFURs, type = 'p', pch = 16, ylab = 'effective fur conductivity, W K-1 m-1', xlab = 'fur depth, mm')

## -----------------------------------------------------------------------------
# input
AMASS <- 10 # animal mass (kg)
ANDENS <- 1000 # animal density (kg/m3)
SUBQFAT <- 1 # is subcutaneous fat present? (0 is no, 1 is yes)
FATPCT <- 20 # body fat (%)
SHAPE <- 4 # shape, 1 is cylinder, 2 is sphere, 3 is plate, 4 is ellipsoid
SHAPE_B <- 2.7 # current ratio between long and short axis (-)
SHAPE_C <- SHAPE_B # current ratio of length:height (for plate)
DHARA <- DHAR[1] # fur diameter, mean (m) (from IRPROP)
RHOARA <- RHOAR[1] # hair density, mean (1/m2) (from IRPROP)
ZFUR <- ZZFUR[1] # fur depth, mean (m) (from IRPROP)
PCOND <- 0 # fraction of body in contact with substrate (fractional, 0-1)
SAMODE <- 0 # if 1, uses bird skin surface area scaling from Walsberg
# and King 1978. JEB Biology 76:185–189, if 2, uses mammal surface area
# scaling from Stahl (1967) J. of App. Physiology, 453–460.
ORIENT <- 0 # if 1, largest surface area normal to sun's ray's, if 2, largest surface parallel to sun's rays, if 0, average of normal/parallel posture
Z <- 20 # zenith angle of the sun
ZEN <- pi/180*Z # convert degrees to radians

# call the subroutine
GEOM.out <- GEOM_ENDO(AMASS, ANDENS, FATPCT, SHAPE, ZFUR, SUBQFAT, SHAPE_B, SHAPE_C, DHARA, RHOARA, PCOND, SAMODE, ORIENT, ZEN)

# output
VOL <- GEOM.out[1] # volume (m3)
D <- GEOM.out[2] # characteristic dimension for convection (m)
MASFAT <- GEOM.out[3] # mass body fat (kg)
VOLFAT <- GEOM.out[4] # volume body fat (m3)
ALENTH <- GEOM.out[5] # length (m)
AWIDTH <- GEOM.out[6] # width (m)
AHEIT <- GEOM.out[7] # height (m)
ATOT <- GEOM.out[8] # total area at fur/feathers-air interface (m2)
ASIL <- GEOM.out[9] # silhouette area to use in solar calcs (m2) may be normal, parallel
# or average set via ORIENT
ASILN <- GEOM.out[10] # silhouette area normal to sun (m2)
ASILP <- GEOM.out[11] # silhouette area parallel to sun (m2)
GMASS <- GEOM.out[12] # mass (g)
AREASKIN <- GEOM.out[13] # area of skin (m2)
FLSHVL <- GEOM.out[14] # flesh volume (m3)
FATTHK <- GEOM.out[15] # fat layer thickness (m)
ASEMAJ <- GEOM.out[16] # semimajor axis length (m)
BSEMIN <- GEOM.out[17] # b semiminor axis length (m)
CSEMIN <- GEOM.out[18] # c semiminor axis length (m)
#(currently only prolate spheroid)
CONVSK <- GEOM.out[19] # area of skin for evaporation (total skin area - hair area) (m2)
CONVAR <- GEOM.out[20] # area for convection (total area minus ventral area in contact 
#with the substrate, as determined by PCOND), m2
R1 <- GEOM.out[21] # shape-specific core-skin radius in shortest dimension (m)
R2 <- GEOM.out[22] # shape-specific core-fur radius in shortest dimension (m)

GEOM.lab <- c("VOL", "D", "MASFAT", "VOLFAT", "ALENTH", "AWIDTH", "AHEIT", "ATOT", "ASIL", "ASILN", "ASILP", "GMASS", "AREASKIN", "FLSHVL", "FATTHK", "ASEMAJ", "BSEMIN", "CSEMIN", "CONVSK", "CONVAR", "R1", "R2")
kable(cbind(GEOM.lab, t(GEOM.out)))

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
if(SHAPE == 4){ #ellipsoid
par(mfrow=c(1,2))
plot(c(0,ASEMAJ*2+ZFUR*2), c(0,ASEMAJ*2+ZFUR*2), type="n", main="ellipsoid, sagittal section", ylab = 'minor axis, m', xlab = 'major axis, m', asp=1)
draw.ellipse(ASEMAJ+ZFUR, ASEMAJ+ZFUR, col="black", border = "black", a = ASEMAJ+ZFUR, b = BSEMIN+ZFUR)
draw.ellipse(ASEMAJ+ZFUR, ASEMAJ+ZFUR, col="pink", border = "pink", a = ASEMAJ, b = BSEMIN)
draw.ellipse(ASEMAJ+ZFUR, ASEMAJ+ZFUR, col="red", border = "red", a = ASEMAJ-FATTHK, b = BSEMIN-FATTHK)

plot(c(0,ASEMAJ*2+ZFUR*2), c(0,ASEMAJ*2+ZFUR*2), type="n", main="ellipsoid, transverse section", ylab = 'minor axis, m', xlab = 'minor axis, m', asp=1)
draw.ellipse(ASEMAJ+ZFUR, ASEMAJ+ZFUR, col="black", border = "black", a = BSEMIN+ZFUR, b = CSEMIN+ZFUR)
draw.ellipse(ASEMAJ+ZFUR, ASEMAJ+ZFUR, col="pink", border = "pink", a = BSEMIN, b = CSEMIN)
draw.ellipse(ASEMAJ+ZFUR, ASEMAJ+ZFUR, col="red", border = "red", a = BSEMIN-FATTHK, b = CSEMIN-FATTHK)
par(mfrow=c(1,1))
}

## -----------------------------------------------------------------------------
# environmental variables
QSOLR <- 1000 # solar radiation, horizontal plane (W/m2)
SHADE <- 0 # shade (%, 0-100)
Z <- 20 # zenith angle of sun (degrees from overhead)
ABSSB <- 0.8 # solar absorptivity of substrate (fractional, 0-1)

# traits
AREA <- ATOT # surface area for solar exchange, m2 (from GEOM)
ABSAND <- 0.8 # solar absorptivity of dorsal fur (fractional, 0-1)
ABSANV <- 0.8 # solar absorptivity of ventral fur (fractional, 0-1)
ASIL <- ASIL # silhouette area to sun, m2 (from GEOM)
PDIF <- 0.15 # proportion of solar radiation that is diffuse (fractional, 0-1)
FASKY <- 0.4 # configuration factor to sky (-)
FAVEG <- 0 # configuration factor to vegetation (-)

# solar radiation normal to sun's rays
ZEN <- pi / 180 * Z # convert degrees to radians
if(Z < 90){ # compute solar radiation on a surface normal to the direct rays of the sun
  CZ <- cos(ZEN)
  QNORM <- QSOLR / CZ
}else{ # diffuse skylight only
  QNORM <- QSOLR
}

SOLAR.out <- SOLAR_ENDO(AREA, ABSAND, ABSANV, ABSSB, ASIL, PDIF, QNORM, SHADE, QSOLR, FASKY, FAVEG)

QSOLAR <- SOLAR.out[1]  # total (global) solar radiation (W)
QSDIR <- SOLAR.out[2]  # direct solar radiation (W)
QSSKY <- SOLAR.out[3]  # diffuse solar radiation from sky (W)
QSRSB <- SOLAR.out[4]  # diffuse solar radiation reflected from substrate (W)
QSDIFF <- SOLAR.out[5]  # total diffuse solar radiation (W)
QDORSL <- SOLAR.out[6]  # dorsal direct solar radiation (W)
QVENTR <- SOLAR.out[7]  # ventral diffuse solar radiation (W)

SOLAR.lab <- c("QSOLAR", "QSDIR", "QSSKY", "QSRSB", "QSDIFF", "QDORSL", "QVENTR")
kable(cbind(SOLAR.lab, t(SOLAR.out)))

## -----------------------------------------------------------------------------
# input
TS <- 33 # skin temperature (°C)
TENV <- 20 # air temperature (°C)
TFA <- 10 # fur/air interface temperature (°C)
SHAPE <- 4 # shape, 1 is cylinder, 2 is sphere, 3 is plate, 4 is ellipsoid
SURFAR <- CONVAR # surface area for convection, m2 (from GEOM)
FLTYPE <- 0 # FLUID TYPE: 0 = AIR; 1 = FRESH WATER; 2 = SALT WATER
FURTST <- FURTST # test of fur presence (-) from IRPROP 
VEL <- 1 # wind speed (m/s)
ELEV <- 0 # elevation (m)
BP <- -1 # barometric pressure (Pa), negative means altitude is used

# run subroutine
CONV.out <- CONV_ENDO(TS, TENV, SHAPE, SURFAR, FLTYPE, FURTST, D, TFA, VEL, ZFUR, BP, ELEV)

QCONV <- CONV.out[1] # convective heat loss (W)
HC <- CONV.out[2] # combined convection coefficient (W/m2 C)
HCFREE <- CONV.out[3] # free convection coefficient (W/m2 C)
HCFOR <- CONV.out[4] # forced convection coefficient (W/m2 C)
HD <- CONV.out[5] # mass transfer coefficient (kg/m2 s)
HDFREE <- CONV.out[6] # free mass transfer coefficient (kg/m2 s)
HDFORC <- CONV.out[7] # forced mass transfer coefficient (kg/m2 s)
ANU <- CONV.out[8] # Nusselt number (-)
RE <- CONV.out[9] # Reynold's number (-)
GR <- CONV.out[10] # Grasshof number (-)
PR <- CONV.out[11] # Prandtl number (-)
RA <- CONV.out[12] # Rayleigh number (-)
SC <- CONV.out[13] # Schmidt number (-)
BP <- CONV.out[14] # barometric pressure (Pa)

CONV.lab <- c("QCONV", "HC", "HCFREE", "HCFOR", "HD", "HDFREE", "HDFORC", "ANU", "RE", "GR", "PR", "RA", "SC", "BP")
kable(cbind(CONV.lab, t(CONV.out)))

## -----------------------------------------------------------------------------
BP <- BP # barometric pressure (Pa) (from CONV)
TA <- 20 # air temperature (°C)
RELHUM <- 20 # relative humidity (%)
TC <- 37 # core temperature (°C)
TSKIN <- 33 # skin temperature (°C)
PCTWET <- 11 # part of the skin surface that is wet (%)
FLYHR <- 0 # is flight occurring this hour? (imposes forced evaporative loss)
BAREVAP <- 0 # is evaporation partly from bare skin? (0 = no, 1 = yes, % defined with PCTBAREVAP)
PCTBAREVAP <- 0 # surface area for evaporation that is skin, e.g. licking paws (%)
PCTEYES <- 0.03 # surface area made up by the eye (%) - make zero if sleeping
ZFUR <- ZFUR # fur depth (m)
FURWET <- 0 # part of the fur surface that is wet (%)

SEVAP.out <- SEVAP_ENDO(BP, TA, RELHUM, VEL, TC, TSKIN, ELEV, PCTWET, FLYHR, CONVSK, HD, HDFREE, PCTBAREVAP, PCTEYES, ZFUR, FURWET, TFA, CONVAR)

QSEVAP <- SEVAP.out[1] # skin evaporative heat loss (W)
WEYES <- SEVAP.out[2] # ocular evaporation (kg/s)
WCUTHF <- SEVAP.out[3] # forced cutaneous evaporation (kg/s)
WCUTF <- SEVAP.out[4] # free cutaneous evaporation (kg/s)
WCUT <- SEVAP.out[5] # total cutaneous evaporation (kg/s)
WTFUR <- SEVAP.out[6] # total fur evaporation (kg/s)
QFSEVAP <- SEVAP.out[7] # fur evaporative heat loss (W))

SEVAP.lab <- c("QSEVAP", "WEYES", "WCUTHF", "WCUTF", "WCUT", "WTFUR", "QFSEVAP")
kable(cbind(SEVAP.lab, t(SEVAP.out)))

## -----------------------------------------------------------------------------
# environment
FLTYPE <- 0 # fluid type: 0 = air; 1 = fresh water; 2 = salt water
TAREF <- TA # reference air temperature (e.g. at 1.2 or 2m) (°C)
TLOWER <- TA # ground temperature (°C)
TSKY <- TA # sky temperature (°C)
TCONDSB <- TA # surface temperature for conduction (°C)
TBUSH <- TA # bush temperature (°C)
TVEG <- TAREF # assume vegetation casting shade is at reference (e.g. 1.2m or 2m) air temperature (°C)
RH <- 20 # relative humidity (%)
VEL <- 1 # wind speed (m/s)
BP <- BP # barometric pressure (Pa), negative means altitude is used (from CONV)
ELEV <- 0 # elevation (m)
KSUB <- 2.79 # substrate thermal conductivity (W/m°C)

# physiology and morphology 
PCTWET <- 0 # part of the skin surface that is wet (%)
AK1 <- 0.9 # initial thermal conductivity of flesh (0.412 - 2.8 W/mK)
AK2 <- 0.230 # conductivity of fat (W/mK)
XR <- 1 # fractional depth of fur at which longwave radiation is exchanged (0-1)
EMISAN <- 0.99 # animal emissivity (-)
BAREVAP <- 0 # is evaporation partly from bare skin? (0 = no, 1 = yes, % defined with PCTBAREVAP)
PCTBAREVAP <- 0 # surface area for evaporation that is skin, e.g. licking paws (%)
PCTEYES <- 0.03 # surface area made up by the eye (%) - make zero if sleeping
# behaviour
FLYHR <- 0 # is flight occurring this hour? (imposes forced evaporative loss)

# configuration factors
FAGRD <- 0.5 # configuration factor to ground (-)
FASKY <- 0.5 # configuration factor to sky (-)
FAVEG <- 0 # this is for overhead veg (at TAREF) (-)
FABUSH <- 0 # this is for veg below/around animal (at TALOC) (-)

# reference configuration factors
FABUSHREF <- FABUSH # nearby bush (-)
FASKYREF <- FASKY # sky (-)
FAGRDREF <- FAGRD # ground (-)
FAVEGREF <- FAVEG # vegetation (-)

FURTHRMK <- 0 # user-specified fur thermal conductivity (W/m C), not used if 0
SHAPE <- 4 # shape, 1 is cylinder, 2 is sphere, 3 is plate, 4 is ellipsoid

# Initial values
TS <- TC # current guess of skin temperature (deg C)
TFA <- TA # current guess of fur/air interface temperature (deg C)

DIFTOL <- 0.001 # tolerance for SIMULSOL

# vector to hold the SIMULSOL results for dorsal and ventral side
SIMULSOL.out <- matrix(data = 0, nrow = 2, ncol = 16) 

# repeat for each side, dorsal and ventral, of the animal
for(S in 1:2){ 
  
 # Calculating solar intensity entering fur. This will depend 
 # on whether we are calculating the fur temperature for the 
 # dorsal side or the ventral side. The dorsal side will have
 # solar inputs from the direct beam hitting the silhouette 
 # area as well as diffuse solar scattered from the sky. The
 # ventral side will have diffuse solar scattered off the 
 # substrate.
  
# Setting config factors and solar depending on whether the dorsal side (S=1) or ventral side (S=2) is being estimated.

  if(QSOLAR > 0.0){
    if(S == 1){
      # proportion of upward view that is sky
      FASKY <- FASKYREF /(FASKYREF + FAVEGREF)
      # proportion of upward view that is vegetation (shade)
      FAVEG <- FAVEGREF / (FASKYREF + FAVEGREF) 
      FAGRD <- 0.0
      FABUSH <- 0.0
      QSLR <- 2 * QSDIR + ((QSSKY / FASKYREF) * FASKY) # direct x 2
      # because assuming sun in both directions, and un-adjusting QSSKY
      # for config factor imposed in SOLAR_ENDO and back to new larger
      # one in both directions
    }else{  # doing ventral side
      FASKY <- 0.0
      FAVEG <- 0.0
      FAGRD <- FAGRDREF / (FAGRDREF + FABUSHREF)
      FABUSH <- FABUSHREF / (FAGRDREF + FABUSHREF)
      QSLR <- (QVENTR / (1 - FASKYREF - FAVEGREF)) * (1 - (2 * PCOND)) # un-adjust by 
      # config factor imposed in SOLAR_ENDO to have it coming in both 
      # directions, but also cutting down according to fractional area
      # conducting to ground (in both directions)
    }
  }else{
    QSLR <- 0.0
    if(S==1){
      FASKY <- FASKYREF / (FASKYREF + FAVEGREF)
      FAVEG <- FAVEGREF / (FASKYREF + FAVEGREF)
      FAGRD <- 0.0
      FABUSH <- 0.0
    }else{
      FASKY <- 0.0
      FAVEG <- 0.0
      FAGRD <- FAGRDREF / (FAGRDREF + FABUSHREF)
      FABUSH <- FABUSHREF / (FAGRDREF + FABUSHREF)
    }
  }
  
  # set fur depth and conductivity
  # index for KEFARA, the conductivity, is the average (1), front/dorsal (2), 
  # back/ventral(3) of the body part
  if(QSOLR > 0 | ZZFUR[2] != ZZFUR[3]){
    if(S == 1){
      ZL <- ZZFUR[2]
      KEFF <- KEFARA[2]
    }else{
      ZL <- ZZFUR[3]
      KEFF <- KEFARA[3]
    }
  }else{
    ZL <- ZZFUR[1]
    KEFF <- KEFARA[1]
  }
  
  RSKIN <- R1 # body radius (including fat), m
  RFLESH <- R1 - FATTHK # body radius flesh only (no fat), m
  RFUR <- R1 + ZL # body radius including fur, m
  D <- 2 * RFUR # diameter, m
  RRAD <- RSKIN + (XR * ZL) # effective radiation radius, m
  LEN <- ALENTH # length, m
  
  if(SHAPE != 4){ #! For cylinder and sphere geometries
   RFURCMP <- RSKIN + ZFURCOMP
  }else{
   RFURCMP <- RFUR #! Note that this value is never used if conduction not being modeled,
   # but need to have a value for the calculations
   }

   if(SHAPE == 4){  #! For ellipsoid geometry
    BLCMP <- BSEMIN + ZFURCOMP
   }else{
    BLCMP <- RFUR #! Note that this value is never used if conduction not being modeled, 
    # but need to have a value for the calculations
   }
  
  # Correcting volume to account for subcutaneous fat
  if(SUBQFAT == 1 & FATTHK > 0.0){
    VOL <- FLSHVL
  }
  
  # Calculating the "Cd" variable: Qcond = Cd(Tskin-Tsub), where Cd = Conduction area*ksub/subdepth
  if(S == 2){
   AREACND <- ATOT * (PCOND * 2)
   CD <- (AREACND * KSUB) / 0.025 # assume conduction happens from 2.5 cm depth
  }else{ #doing dorsal side, no conduction. No need to adjust areas used for convection.
   AREACND <- 0
   CD <- 0
  }
  
  
  # package up inputs
  FURVARS <- c(LEN,ZFUR,FURTHRMK,KEFF,BETARA,FURTST,ZL,LHAR[S+1],DHAR[S+1],RHOAR[S+1],REFLFR[S+1],KHAIR,S)
  GEOMVARS <- c(SHAPE,SUBQFAT,CONVAR,VOL,D,CONVAR,CONVSK,RFUR,RFLESH,RSKIN,XR,RRAD,ASEMAJ,BSEMIN,CSEMIN,CD,PCOND,RFURCMP,BLCMP,KFURCMPRS)
  ENVVARS <- c(FLTYPE,TA,TS,TBUSH,TVEG,TLOWER,TSKY,TCONDSB,RH,VEL,BP,ELEV,FASKY,FABUSH,FAVEG,FAGRD,QSLR)
  TRAITS <- c(TC,AK1,AK2,EMISAN,FATTHK,FLYHR,FURWET,PCTBAREVAP,PCTEYES)
  
  # set IPT, the geometry assumed in SIMULSOL: 1 = cylinder, 2 = sphere, 3 = ellipsoid
  if(SHAPE %in% c(1, 3)){
    IPT <- 1
  }
  if(SHAPE == 2){
    IPT <- 2
  }
  if(SHAPE == 4){
    IPT <- 3
  }
  
  # call SIMULSOL
  SIMULSOL.out[S,] <- SIMULSOL(DIFTOL, IPT, FURVARS, GEOMVARS, ENVVARS, TRAITS, TFA, PCTWET, TS)
}

SIMULSOL.out <- cbind(c(1,2), SIMULSOL.out)
colnames(SIMULSOL.out) <- c("SIDE", "TFA", "TSKIN", "QCONV", "QCOND", "QGENNET", "QSEVAP", "QRAD", "QSLR", "QRSKY", "QRBSH", "QRVEG", "QRGRD", "QFSEVAP", "NTRY", "SUCCESS", "KFUR")
tSIMULSOL.out <- t(SIMULSOL.out)
colnames(tSIMULSOL.out) <- c("DORSAL", "VENTRAL")
kable(tSIMULSOL.out)

## -----------------------------------------------------------------------------
# define basal metabolic rate
QBASAL <- (70 * AMASS ^ 0.75) * (4.185 / (24 * 3.6)) # heat generation (W), from Kleiber 1947
DELTAR <- 0 # offset between air temperature and breath (°C)
O2GAS <- 20.95 # oxygen concentration of air (%)
N2GAS <- 79.02 # nitrogen concentration of air (%)
CO2GAS <- 0.0412 # carbon dioxide concentration of air (%)
R_PCO2 <- CO2GAS / 100 # reference atmospheric dioxide concentration of air (proportion),
# to allow for anthropogenic change (%)
RQ <- 0.80 # respiratory quotient (fractional, typically 0.7 (fats)
# to 1 (carbs), with 0.8 typical for protein)
EXTREF <- 20 # O2 extraction efficiency (%)
RELXIT <- 100 # relative humidity of exhaled air (%)
PANT <- 1 # multiplier on breathing rate (-)

# Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(FASKY+FAVEG))+(ventral value*FAGRD)
GEND <- SIMULSOL.out[1, 6]
GENV <- SIMULSOL.out[2, 6]
DMULT <- FASKYREF + FAVEGREF
VMULT <- 1 - DMULT # Assume that reflectivity of veg below equals reflectivity of soil 
# so VMULT left as 1 - DMULT
X <- GEND * DMULT + GENV * VMULT # weighted estimate of metabolic heat generation

# lung temperature and temperature of exhaled air
TS <- (SIMULSOL.out[1, 2] + SIMULSOL.out[2, 2]) * 0.5
TLUNG <- (TC + TS) * 0.5 # average of skin and core
TAEXIT <- min(TA + DELTAR, TLUNG) # temperature of exhaled air, °C
QMIN <- QBASAL
QM1 <- X - (5 * QMIN) 
QM2 <- X + (10 * QMIN)
QSUM <- X
TOL <- AMASS * 0.01
ZBRENT.in <- c(TA, O2GAS, N2GAS, CO2GAS, BP, QMIN, RQ, TLUNG, GMASS, EXTREF, RH, RELXIT, 1, TAEXIT, QSUM, PANT, R_PCO2)

# call ZBRENT subroutine which calls RESPFUN
ZBRENT.out <- ZBRENT_ENDO(QM1, QM2, TOL, ZBRENT.in)
colnames(ZBRENT.out) <- c("RESPFN","QRESP","GEVAP", "PCTO2", "PCTN2", "PCTCO2", "RESPGEN", "O2STP", "O2MOL1", "N2MOL1", "AIRML1", "O2MOL2", "N2MOL2", "AIRML2", "AIRVOL")
tZBRENT.out <- t(ZBRENT.out)
colnames(tZBRENT.out) <- c("OUTPUT")
kable(tZBRENT.out)

