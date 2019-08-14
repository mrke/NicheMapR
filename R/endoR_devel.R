#' endoR_devel - the development version of the endotherm model of NicheMapR
#'
#' This model uses R code to implement postural and physiological
#' thermoregulation under a given environmental scenario for an organism of
#' a specified shape and no extra body parts. In this function the sequence of
#' thermoregulatory events in the face of heat stress is to first uncurl,
#' second change flesh conductivity, third raise core temprature, fourth
#' pant and fifth sweat. This can be modified to be more specific to the
#' species of interest, e.g. by changing the sequence of responses or having
#' some happen in parallel. However it is very slow. endoR_solvendo implements
#' the endoR sequence of thermoregulation within FORTRAN and is 100x faster.
#' Thus it is best to use endoR as a basis for prototyping and refining and
#' then to adjust the FORTRAN code of endoR_solvendo (SOLVENDO.f) accordingly.
#' @encoding UTF-8
#' @param AMASS = 1, # kg
#' @param NGEOM = 4, # cylinder (ngeom = 1), sphere (ngeom = 2), plate (ngeom = 3) and ellipsoid (ngeom = 4)
#' @param GMREF = 3, # initial ratio between long and short axis (-)
#' @param FURTHRMK = 0, # user-specified fur thermal conductivity (W/mK), not used if 0
#' @param ZFURD = 2E-03, # fur depth, dorsal (m)
#' @param ZFURV = 2E-03, # fur depth, ventral (m)
#' @param TC = 37, # core temperature (°C)
#' @param TCMAX = 45, # maximum core temperature (°C)
#' @param TA = 20, air temperature at local height (°C)
#' @param TGRD = TA, ground temperature (°C)
#' @param TSKY = TA, sky temperature (°C)
#' @param VEL = 0.1, wind speed (m/s)
#' @param RH = 5, relative humidity (\%)
#' @param QSOLR = 0, solar radiation, horizontal plane (W/m2)
#' @param Z = 20, zenith angle of sun (degrees from overhead)
#' @param SHADE = 0, shade level (\%)
#' @usage endoR(AMASS = 1, NGEOM = 4, GMREF = 3, FURTHRMK = 0, ZFURD = 2E-03, ZFURV = 2E-03, TC = 37, TCMAX = 45, TA = 20, TGRD = TA, TSKY = TA, VEL = 0.1, RH = 5, QSOLR = 0, Z = 20, SHADE = 0, NITESHAD = 0,...)
#' @details
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#' \code{DIFTOL}{ = 0.001, error tolerance for SIMULSOL (°C)}\cr\cr
#'
#' \strong{ Environment:}\cr\cr
#' \code{TAREF}{ = TA, air temperature at reference height (°C)}\cr\cr
#' \code{ELEV}{ = 0, elevation (m)}\cr\cr
#' \code{ABSSB}{ = 0.8, solar absorptivity of substrate (fractional, 0-1)}\cr\cr
#' \code{FLTYPE}{ = 0, FLUID TYPE: 0 = AIR; 1 = FRESH WATER; 2 = SALT WATER - needs to be looked at - only invoked in main program when the dive table is set up}\cr\cr
#' \code{TCONDSB}{ = TGRD, surface temperature for conduction (°C)}\cr\cr
#' \code{TBUSH}{ = TA, bush temperature (°C)}\cr\cr
#' \code{BP}{ = -1, Pa, negatve means elevation is used}\cr\cr
#' \code{O2GAS}{ = 20.95, oxygen concentration of air (\%)}\cr\cr
#' \code{N2GAS}{ = 79.02, nitrogen concetration of air (\%)}\cr\cr
#' \code{CO2GAS}{ = 0.03, carbon dioxide concentration of air (\%)}\cr\cr
#' \code{PCTDIF}{ = 0.15, proportion of solar radiation that is diffuse (fractional, 0-1)}\cr\cr
#'
#' \strong{ Behaviour:}\cr\cr
#' \code{SHADE}{ = 0, shade level (\%)}\cr\cr
#' \code{NITESHAD}{ = 0, flag for if animal is behaviourally seeking shade for warmth at night - remove?}\cr\cr
#' \code{FLYHR}{ = 0, is flight occuring this hour? (imposes forced evaporative loss)}\cr\cr
#' \code{UNCURL}{ = 1, allows the animal to uncurl to GMULTMAX, the value being the increment GMULT is increased per iteration}\cr\cr
#' \code{RAISETC}{ = 1, turns on core temperature elevation, the value being the increment by which TC is increased per iteration}\cr\cr
#' \code{SWEAT}{ = 0.25, turns on sweating, the value being the increment by which SKINW is increased per iteration (\%)}\cr\cr
#' \code{MXWET}{ = 100, maximum surface area that can be wet (\%)}\cr\cr
#' \code{AK1inc}{ = 0.5, turns on thermal conductivity increase (W/mK), the value being the increment by which AK1 is increased per iteration (W/mC)}\cr\cr
#' \code{AKMAX}{ = 2.8, maximum flesh conductivity (W/mK)}\cr\cr
#' \code{PANT}{ = 1, multiplier on breathing rate to simulate panting (-)}\cr\cr
#' \code{PANTING}{ = 0.1, increment for multiplier on breathing rate to simulate panting (-)}\cr\cr
#'
#' \strong{ General morphology:}\cr\cr
#' \code{ANDENS}{ = 1000, body density (kg/m3)}\cr\cr
#' \code{SUBQFAT}{ = 0, is subcutaneous fat present? (0 is no, 1 is yes)}\cr\cr
#' \code{FATPCT}{ = 20, \% body fat}\cr\cr
#' \code{GMULT}{ = GMREF, current ratio between long and short axis (-)}\cr\cr
#' \code{GMULTMAX}{ = GMREF, max possible ratio between long and short axis (-)}\cr\cr
#' \code{MAXPTVEN}{ = 0.5, maxium fraction of surface area that is ventral (fractional, 0-1)}\cr\cr
#' \code{PTCOND}{ = 0, \% of body area touching the substrate}\cr\cr
#' \code{BIRD}{ = 0, if 1, uses bird skin surface area scaling from Walsberg, G. E., and J. E. King. 1978. The Relationship of the External Surface Area of Birds to Skin Surface Area and Body Mass. Journal of Experimental Biology 76:185–189}\cr\cr
#' \code{MAMMAL}{ = 0, if 1, uses mammal surface area scaling from Stahl W. R. (1967) Scaling of respiratory variables in mammals. Journal of Applied Physiology 22 , 453–460.}\cr\cr
#' \code{ORIENT}{ = 0, if 0, long axis parallel to ground, if 1, long axis is perpendicular to the ground}\cr\cr
#'
#' \strong{ Fur properties:}\cr\cr
#' \code{DHAIRD}{ = 30E-06, hair diameter, dorsal (m)}\cr\cr
#' \code{DHAIRV}{ = 30E-06, hair diameter, ventral (m)}\cr\cr
#' \code{LHAIRD}{ = 23.9E-03, hair length, dorsal (m)}\cr\cr
#' \code{LHAIRV}{ = 23.9E-03, hair length, dorsal (m)}\cr\cr
#' \code{RHOD}{ = 3000E+04, hair density, dorsal (1/m2)}\cr\cr
#' \code{RHOV}{ = 3000E+04, hair density, ventral (1/m2)}\cr\cr
#' \code{REFLD}{ = 0.2, fur reflectivity dorsal (fractional, 0-1)}\cr\cr
#' \code{REFLV}{ = 0.2, fur reflectivity ventral (fractional, 0-1)}\cr\cr
#'
#' \strong{ Radiation exchange:}\cr\cr
#' \code{EMISAN}{ = 0.99, animal emissivity (-)}\cr\cr
#' \code{FATOBJ}{ = 0, configuration factor to nearby object}\cr\cr
#' \code{FABUSH}{ = 0, this is for veg below/around animal (at TALOC)}\cr\cr
#' \code{FGDREF}{ = 0.5, reference configuration factor to ground}\cr\cr
#' \code{FSKREF}{ = 0.5, configuration factor to sky}\cr\cr
#'
#' \strong{ Nest properties:}\cr\cr
#' \code{NESTYP}{ = 0, # for nest calculations, to do)}\cr\cr
#' \code{RoNEST}{ = 0, # for nest calculations, to do}\cr\cr
#'
#' \strong{ Physiology:}\cr\cr
#' \code{AK1}{ = 0.9, # initial thermal conductivity of flesh (0.412 - 2.8 W/mK)}\cr\cr
#' \code{AK2}{ = 0.230, # conductivity of fat (W/mK)}\cr\cr
#' \code{QBASAL}{ = (70 \* AMASS ^ 0.75) \* (4.185 / (24 \* 3.6)), # basal heat generation (W)}\cr\cr
#' \code{SKINW}{ = 0, # part of the skin surface that is wet (\%)}\cr\cr
#' \code{BAREVAP}{ = 0, # is evaporation partly from bare skin? (0 = no, 1 = yes, \% defined with PCTSKINEVAP)}\cr\cr
#' \code{PCTBAREVAP}{ = 2.8, maximum flesh conductivity (W/mK)}\cr\cr
#' \code{PCTEYES}{ = 0, # surface area made up by the eye (\%) - make zero if sleeping}\cr\cr
#' \code{DELTAR}{ = 0, # offset between air temperature and breath (°C)}\cr\cr
#' \code{RELXIT}{ = 100, # relative humidity of exhaled air, \%}\cr\cr
#' \code{TIMACT}{ = 1, # multiplier on metabolic rate for activity costs}\cr\cr
#' \code{RQ}{ = 0.80, # respiratory quotient (fractional, 0-1)}\cr\cr
#' \code{EXTREF}{ = 20, # O2 extraction efficiency (\%)}\cr\cr
#' \code{PANTMAX}{ = 10, # maximum breathing rate multiplier to simulate panting (-)}\cr\cr
#' \code{Q10}{ = 1, # Q10 factor for adjusting BMR for TC}\cr\cr
#'
#' \strong{ Initial conditions:}\cr\cr
#' \code{TS}{ = TC - 3, # initial skin temperature (°C)}\cr\cr
#' \code{TFA}{ = TA, # initial fur/air interface temperature (°C)}\cr\cr
#'
#' \strong{Outputs:}
#' \code{TC}{core temperature (°C)}\cr\cr
#' \code{TFA_D}{dorsal fur-air interface temperature (°C)}\cr\cr
#' \code{TFA_V}{ventral fur-air interface temperature (°C)}\cr\cr
#' \code{TSKIN_D}{dorsal skin temperature (°C)}\cr\cr
#' \code{TSKIN_V}{ventral skin temperature (°C)}\cr\cr
#' \code{QCONV_D}{dorsal convection heat exchange (W)}\cr\cr
#' \code{QCONV_V}{ventral convection heat exchange (W)}\cr\cr
#' \code{QCOND_D}{dorsal conduction heat exchange (W)}\cr\cr
#' \code{QCOND_V}{ventral conduction heat exchange (W)}\cr\cr
#' \code{QGENNET_D}{dorsal net heat generation (W)}\cr\cr
#' \code{QGENNET_V}{ventral net heat generation (W)}\cr\cr
#' \code{QSEVAP_D}{dorsal evaporative heat exchange (W)}\cr\cr
#' \code{QSEVAP_V}{ventral evaporative heat exchange (W)}\cr\cr
#' \code{QRAD_D}{dorsal radiant heat loss (W)}\cr\cr
#' \code{QRAD_V}{ventral radiant heat loss (W)}\cr\cr
#' \code{QSLR_D}{dorsal solar heat gain (W)}\cr\cr
#' \code{QSLR_V}{ventral solar heat gain (W)}\cr\cr
#' \code{QRSKY_D}{dorsal radiant heat incomming from sky (W)}\cr\cr
#' \code{QRSKY_V}{ventral radiant heat incomming from sky (W)}\cr\cr
#' \code{QRBSH_D}{dorsal radiant heat incomming from nearby bush (W)}\cr\cr
#' \code{QRBSH_V}{ventral radiant heat incomming from nearby bush (W)}\cr\cr
#' \code{QRVEG_D}{dorsal radiant heat incomming from vegetation (W)}\cr\cr
#' \code{QRVEG_V}{ventral radiant heat incomming from vegetation (W)}\cr\cr
#' \code{QRGRD_D}{dorsal radiant heat incomming from ground (W)}\cr\cr
#' \code{QRGRD_V}{ventral radiant heat incomming from ground (W)}\cr\cr
#' \code{NTRY_D}{number of iterations need for convergence of dorsal heat budget}\cr\cr
#' \code{NTRY_V}{number of iterations need for convergence of ventral heat budget}\cr\cr
#' \code{SUCCESS_D}{test of success convergence for dorsal heat budget}\cr\cr
#' \code{SUCCESS_V}{test of success convergence for ventral heat budget}\cr\cr
#' \code{RESPFN}{energy balance test after call to RESPFUN (W)}\cr\cr
#' \code{QRESP}{respiratory heat exchange (W)}\cr\cr
#' \code{GEVAP}{respiratory water loss (g/s)}\cr\cr
#' \code{PCTO2}{ambient oxygen gas concentration (\%)}\cr\cr
#' \code{PCTN2}{ambient nitgrogen gas concentration (\%)}\cr\cr
#' \code{PCTCO2}{ambient carbon dioxide gas concentration (\%)}\cr\cr
#' \code{RESPGEN}{total metabolic rate (W)}\cr\cr
#' \code{O2STP}{oxygen consumption at standard temperature and pressure (L/s)}\cr\cr
#' \code{O2MOL1}{oxygen entering lungs (moles/s)}\cr\cr
#' \code{N2MOL1}{nitrogen entering lungs (moles/s)}\cr\cr
#' \code{AIRML1}{air entering lungs (moles/s)}\cr\cr
#' \code{O2MOL2}{oxygen leaving lungs (moles/s)}\cr\cr
#' \code{N2MOL2}{nitrogen leaving lungs (moles/s)}\cr\cr
#' \code{AIRML2}{air leaving lungs (moles/s)}\cr\cr
#' \code{AIRVOL}{air entering lungs (L/s)}\cr\cr
#' \code{GMULT}{shape multiplier for postural change (-)}\cr\cr
#' \code{SKINW}{skin area that is wet (\%)}\cr\cr
#' \code{SWEAT.G.H}{sweating rate (g/h)}\cr\cr
#' \code{EVAP.G.H}{evaporation rate (g/h)}\cr\cr
#' \code{EXTREF}{oxygen extraction efficiency (\%)}\cr\cr
#' \code{AK}{skin thermal conductivity (W/m°C)}\cr\cr
#' \code{TA}{air temperature (°C)}\cr\cr
#' \code{TGRD}{ground temperature, driving longwave heat gain (°C)}\cr\cr
#' \code{TCONDSB}{substrate temperature, driving conductive heat exchange (°C)}\cr\cr
#' \code{TSKY}{sky temperature (°C)}\cr\cr
#' \code{VEL}{wind speed (m/s)}\cr\cr
#' \code{RH}{relative humidity (\%)}\cr\cr
#' \code{QSOLR}{solar radiation (W/m2)}\cr\cr
#' @examples
#' library(NicheMapR)
#' # environment (central Australia)
#' micro <- micro_global(loc = c(131.05, -22.75), runshade = 0, Usrhyt = 0.01)
#'
#' metout <- as.data.frame(micro$metout)
#' soil <- as.data.frame(micro$soil)
#' days<-rep(seq(1,12),24)
#' days<-days[order(days)]
#' dates<-days+metout$TIME/60/24-1 # dates for hourly output
#'
#' TAs <- metout$TALOC
#' TAREFs <- metout$TAREF
#' TSKYs <- metout$TSKYC
#' TGRDs <- soil$D0cm
#' VELs <- metout$VLOC
#' RHs <- metout$RHLOC
#' QSOLRs <- metout$SOLR
#' Zs <- metout$ZEN
#' ELEV <- micro$elev
#' ABSSB <- 1-micro$REFL
#'
#' # core temperature
#' TC <- 38 # core temperature (deg C)
#' TCMAX <- 43 # maximum core temperature (°C)
#' RAISETC <- 0.25 # increment by which TC is elevated (°C)
#'
#' # size and shape
#' AMASS <- 0.0337 # mass (kg)
#' GMREF <- 1.1 # start off near to a sphere (-)
#' GMULTMAX <- 5 # maximum ratio of length to width/depth
#'
#' # fur/feather properties
#' DHAIRD = 30E-06 # hair diameter, dorsal (m)
#' DHAIRV = 30E-06 # hair diameter, ventral (m)
#' LHAIRD = 23.1E-03 # hair length, dorsal (m)
#' LHAIRV = 22.7E-03 # hair length, ventral (m)
#' ZFURD = 5.8E-03 # fur depth, dorsal (m)
#' ZFURV = 5.6E-03 # fur depth, ventral (m)
#' RHOD = 8000E+04 # hair density, dorsal (1/m2)
#' RHOV = 8000E+04 # hair density, ventral (1/m2)
#' REFLD = 0.248  # fur reflectivity dorsal (fractional, 0-1)
#' REFLV = 0.351  # fur reflectivity ventral (fractional, 0-1)
#'
#' # physiological responses
#' SKINW <- 0.1 # base skin wetness (%)
#' MXWET <- 20 # maximum skin wetness (%)
#' SWEAT <- 0.25 # intervals by which skin wetness is increased (%)
#' Q10 <- 2 # A10 effect of body temperature on metabolic rate
#' QBASAL <- 10 ^ (-1.461 + 0.669 * log10(AMASS * 1000)) # basal heat generation (W) (bird formula from McKechnie and Wolf 2004 Phys. & Biochem. Zool. 77:502-521)
#' DELTAR <- 5 # offset between air temeprature and breath (°C)
#' EXTREF <- 15 # O2 extraction efficiency (%)
#' PANTING <- 0.1 # turns on panting, the value being the increment by which the panting multiplier is increased up to the maximum value, PANTMAX
#' PANTMAX <- 3# maximum panting rate - multiplier on air flow through the lungs above that determined by metabolic rate
#'
#' ptm <- proc.time() # start timing
#' endo.out <- lapply(1:length(TAs), function(x){endoR(TA = TAs[x], TAREF = TAREFs[x], TSKY = TSKYs[x],
#'                                                     TGRD = TGRDs[x], VEL = VELs[x], RH = RHs[x], QSOLR = QSOLRs[x], Z = Zs[x], ELEV = ELEV, ABSSB = ABSSB, TC = TC, TCMAX = TCMAX, AMASS = AMASS, GMREF = GMREF, GMULTMAX = GMULTMAX, SKINW = SKINW, SWEAT = SWEAT, Q10 = Q10, QBASAL = QBASAL, DELTAR = DELTAR, DHAIRD = DHAIRD, DHAIRV = DHAIRV, LHAIRD = LHAIRD, LHAIRV = LHAIRV, ZFURD = ZFURD, ZFURV = ZFURV, RHOD = RHOD, RHOV = RHOV, REFLD = REFLD, RAISETC = RAISETC, PANTING = PANTING, PANTMAX = PANTMAX, EXTREF = EXTREF)})
#' proc.time() - ptm
#' endo.out <- do.call("rbind", lapply(endo.out, data.frame))
#'
#' QGEN <- endo.out$RESPGEN # metabolic rate (W)
#' H2O <- endo.out$GEVAP * 3600 # g/h water evaporated
#' TFA_D <- endo.out$TFA_D # dorsal fur surface temperature
#' TFA_V <- endo.out$TFA_V # ventral fur surface temperature
#' TskinD <- endo.out$TSKIN_D # dorsal skin temperature
#' TskinV <- endo.out$TSKIN_V # ventral skin temperature
#' TCs <- endo.out$TC # core temperature
#' SkinW <- endo.out$SKINW # skin wetness (%)
#' Pant <- endo.out$PANT # panting multiplier (-)
#'
#' par(mfrow = c(2, 2))
#' par(oma = c(2, 1, 2, 2) + 0.1)
#' par(mar = c(3, 3, 1.5, 1) + 0.1)
#' par(mgp = c(2, 1, 0))
#' plot(QGEN ~ dates, type = 'l', ylab = 'metabolic rate, W', xlab = 'time')
#' plot(H2O ~ dates, type = 'l', ylab = 'water loss, g/h', xlab = 'time')
#' plot(TFA_D ~ dates, type = 'l', col = 'grey', ylab = 'fur, skin and core temperature, deg C', xlab = 'time', ylim = c(0, 60))
#' points(TFA_V ~ dates, type = 'l', col = 'grey', lty = 2)
#' points(TskinD ~ dates, type = 'l', col = 'orange')
#' points(TskinV ~ dates, type = 'l', col = 'orange', lty = 2)
#' points(TCs ~ dates, type = 'l', col = 'red')
#' plot(SkinW ~ dates, type = 'l', col = 'black', ylab = 'skin wetness (%)/panting rate (-)', xlab = 'time', ylim = c(0, 20))
#' points(Pant ~ dates, type = 'l', col = 'grey', lty = 2)
#' @export
endoR_devel <- function(
  TA = 20, # air temperature at local height (°C)
  TAREF = TA, # air temeprature at reference height (°C)
  TGRD = TA, # ground temperature (°C)
  TSKY = TA, # sky temperature (°C)
  VEL = 0.1, # wind speed (m/s)
  RH = 5, # relative humidity (%)
  QSOLR = 0, # solar radiation, horizontal plane (W/m2)
  Z = 20, # zenith angle of sun (degrees from overhead)
  ELEV = 0, # elevation (m)
  ABSSB = 0.8, # solar absorptivity of substrate (fractional, 0-1)

  # other environmental variables
  FLTYPE = 0, # FLUID TYPE: 0 = AIR; 1 = FRESH WATER; 2 = SALT WATER - need's to be looked at - only invoked in main program when the dive table is set up
  TCONDSB = TGRD, # surface temperature for conduction (°C)
  TBUSH = TA, # bush temperature (°C)
  BP = -1, # Pa, negatve means elevation is used
  O2GAS = 20.95, # oxygen concentration of air (%)
  N2GAS = 79.02, # nitrogen concetration of air (%)
  CO2GAS = 0.03, # carbon dioxide concentration of air (%)
  PCTDIF = 0.15, # proportion of solar radiation that is diffuse (fractional, 0-1)

  # BEHAVIOUR

  SHADE = 0, # shade level (%)
  NITESHAD = 0, # flag for if animal is behaviourally seeking shade for warmth at night - remove?
  FLYHR = 0, # is flight occuring this hour? (imposes forced evaporative loss)
  UNCURL = 1, # allows the animal to uncurl to GMULTMAX, the value being the increment GMULT is increased per iteration
  RAISETC = 1, # turns on core temperature elevation, the value being the increment by which TC is increased per iteration
  SWEAT = 0.25, # turns on sweating, the value being the increment by which SKINW is increased per iteration
  MXWET = 100, # maximum surface area that can be wet (%)
  AK1inc = 0.5, # turns on thermal conductivity increase (W/mK), the value being the increment by which AK1 is increased per iteration
  AKMAX = 2.8, # maximum flesh conductivity (W/mK)
  PANT = 1, # multiplier on breathing rate to simulate panting (-)
  PANTING = 0.1, # increment for multiplier on breathing rate to simulate panting (-)

  # MORPHOLOGY

  # geometry
  AMASS = 1, # kg
  ANDENS = 1000, # kg/m3
  SUBQFAT = 0, # is subcutaneous fat present? (0 is no, 1 is yes)
  FATPCT = 20, # % body fat
  NGEOM = 4, # cylinder (ngeom = 1), sphere (ngeom = 2) and ellipsoid (ngeom = 4). If a truncated cone (5) or ellipsoidal cylinder (3), we will use the cylinder equations (ngeom=1).
  GMREF = 3, # initial ratio between long and short axis (-)
  GMULT = GMREF, # current ratio between long and short axis (-)
  GMULTMAX = GMREF, # max possible ratio between long and short axis (-)
  MAXPTVEN = 0.5, # maxium fraction of surface area that is ventral (fractional, 0-1)
  PTCOND = 0, # % of body area touching the substrate
  BIRD = 0, # if 1, uses bird skin surface area allometry from Walsberg, G. E., and J. E. King. 1978. The Relationship of the External Surface Area of Birds to Skin Surface Area and Body Mass. Journal of Experimental Biology 76:185–189.
  MAMMAL = 0, # if 1, uses mammal surface area from Stahl W. R. (1967) Scaling of respiratory variables in mammals. Journal of Applied Physiology 22 , 453–460.
  ORIENT = 0, # if 0, long axis parallel to ground, if 1, long axis is perpendicular to the ground

  # fur properties
  FURTHRMK = 0, # user-specified fur thermal conductivity (W/mK), not used if 0
  DHAIRD = 30E-06, # hair diameter, dorsal (m)
  DHAIRV = 30E-06, # hair diameter, ventral (m)
  LHAIRD = 23.9E-03, # hair length, dorsal (m)
  LHAIRV = 23.9E-03, # hair length, ventral (m)
  ZFURD = 2E-03, # fur depth, dorsal (m)
  ZFURV = 2E-03, # fur depth, ventral (m)
  RHOD = 3000E+04, # hair density, dorsal (1/m2)
  RHOV = 3000E+04, # hair density, ventral (1/m2)
  REFLD = 0.2,  # fur reflectivity dorsal (fractional, 0-1)
  REFLV = 0.2,  # fur reflectivity ventral (fractional, 0-1)

  # radiation exchange
  EMISAN = 0.99, # animal emissivity (-)
  FATOBJ = 0, # configuration factor to nearby object
  FABUSH = 0, # this is for veg below/around animal (at TALOC)
  FGDREF = 0.5, # reference configuration factor to ground
  FSKREF = 0.5, # configuration factor to sky

  # nest properties
  NESTYP = 0, # for nest calculations, to do
  RoNEST = 0, # for nest calculations, to do

  # PHYSIOLOGY

  # thermal
  TC = 37, # core temperature (°C)
  TCMAX = 45, # maximum core temperature (°C)
  AK1 = 0.9, # initial thermal conductivity of flesh (0.412 - 2.8 W/mC)
  AK2 = 0.230, # conductivity of fat (W/mK)

  # evaporation
  SKINW = 0, # part of the skin surface that is wet (%)
  BAREVAP = 0, # is evaporation partly from bare skin? (0 = no, 1 = yes, % defined with PCTSKINEVAP)
  PCTBAREVAP = 0, # surface area for evaporation that is skin, e.g. licking paws (%)
  PCTEYES = 0, # surface area made up by the eye (%) - make zero if sleeping
  DELTAR = 0, # offset between air temeprature and breath (°C)
  RELXIT = 100, # relative humidity of exhaled air, %

  # metabolism/respiration
  QBASAL = (70 * AMASS ^ 0.75) * (4.185 / (24 * 3.6)), # basal heat generation (W)
  TIMACT = 1, # multiplier on metabolic rate for activity costs
  RQ = 0.80, # respiratory quotient (fractional, 0-1)
  EXTREF = 20, # O2 extraction efficiency (%)
  PANTMAX = 10, # maximum breathing rate multiplier to simulate panting (-)
  Q10 = 1, # Q10 factor for adjusting BMR for TC

  # initial conditions
  TS = TC - 3, # skin temperature (°C)
  TFA = TA, # fur/air interface temperature (°C)

  # other model settings
  DIFTOL = 0.001 # tolerance for SIMULSOL
){
  if(PANTING == 0){
    PANTMAX <- PANT # can't pant, so panting level set to current value
  }
  if(SWEAT == 0){
    MXWET <- SKINW # can't sweat, so max maximum skin wetness equal to current value
  }
  if(RAISETC == 0){
    TCMAX <- TC # can't raise Tc, so max value set to current value
  }
  if(AK1inc == 0){
    AKMAX <- AK1 # can't change thermal conductivity, so max value set to current value
  }
  if(UNCURL == 0){
    GMULTMAX <- GMULT # can't change posture, so max multiplier of dimension set to current value
  }
  PANTSTEP <- 0
  QGEN <- 0
  TCREF <- TC
  QBASREF <- QBASAL
  # check if heat stressed already
  if(TA >= TC){
    # set core temperature, flesh thermal conductivity and shape to
    # extreme heat loss values, adjusting basal metabolic
    # rate for temperature increase
    if(TA > TCMAX){
      TC <- TCMAX
      Q10mult <- Q10^((TC - TCREF)/10)
      QBASAL = QBASREF * Q10mult
    }
    AK1 <- AKMAX
    GMULT <- GMULTMAX
  }

  while(QGEN < QBASAL){

    ### IRPROP, infrared radiation properties of fur

    # call the IR properties subroutine
    IRPROP.out <- IRPROP(TA, GMULTMAX, GMREF, GMULT, DHAIRD, DHAIRV, LHAIRD, LHAIRV, ZFURD, ZFURV, RHOD, RHOV, REFLD, REFLV, MAXPTVEN)

    # output
    KEFARA <- IRPROP.out[2:4] # effective thermal conductivity of fur array, mean, dorsal, ventral (W/mK)
    BETARA <- IRPROP.out[5:7] # term involved in computing optical thickess (1/mK2)
    B1ARA <- IRPROP.out[8:10] # optical thickness array, mean, dorsal, ventral (m)
    DHAR <- IRPROP.out[11:13] # fur diameter array, mean, dorsal, ventral (m)
    LHAR <- IRPROP.out[14:16] # fur length array, mean, dorsal, ventral (m)
    RHOAR <- IRPROP.out[17:19] # fur density array, mean, dorsal, ventral (1/m2)
    ZZFUR <- IRPROP.out[20:22] # fur depth array, mean, dorsal, ventral (m)
    REFLFR <- IRPROP.out[23:25] # fur reflectivity array, mean, dorsal, ventral (fractional, 0-1)
    FURTST <- IRPROP.out[26] # test of presence of fur (length x diamater x density x depth) (-)

    ### GEOM, geometry

    # input
    DHARA <- DHAR[1] # fur diameter, mean (m) (from IRPROP)
    RHOARA <- RHOAR[1] # hair density, mean (1/m2) (from IRPROP)
    ZFUR <- ZZFUR[1] # fur depth, mean (m) (from IRPROP)

    # call the subroutine
    GEOM.out <- GEOM(AMASS, ANDENS, FATPCT, NGEOM, ZFUR, SUBQFAT, GMULT, GMREF, DHARA, RHOARA, PTCOND, BIRD, MAMMAL, ORIENT)

    # output
    R <- GEOM.out[1] # radius as determined assumming the volume as a sphere, m
    VOL <- GEOM.out[2] # volume, m3
    D <- GEOM.out[3] # diameter as determined assumming the volume as a sphere, m
    MASFAT <- GEOM.out[4] # mass body fat, kg
    VOLFAT <- GEOM.out[5] # volume body fat, m3
    ALENTH <- GEOM.out[6] # length, m
    AWIDTH <- GEOM.out[7] # width, m
    AHEIT <- GEOM.out[8] # height, m
    ATOT <- GEOM.out[9] # total area, m2
    ASILN <- GEOM.out[10] # silhouette area normal to sun, m2
    ASILP <- GEOM.out[11] # silhouette area parallel to sun, m2
    AL <- GEOM.out[12] # effective lenght for convection, m
    GMASS <- GEOM.out[13] # mass, g
    AREASKIN <- GEOM.out[14] # area of skin, m2
    AREA <- GEOM.out[15] # total area at fur/feathers-air interface, m2
    FLSHVL <- GEOM.out[16] # flesh volume, m3
    FATTHK <- GEOM.out[17] # fat layer thickness, m
    ASEMAJ <- GEOM.out[18] # semimajor axis length, m
    BSEMIN <- GEOM.out[19] # b semiminor axis length, m
    CSEMIN <- GEOM.out[20] # c semiminor axis length, m (currently only prolate spheroid)
    CONVSK <- GEOM.out[21] # area of skin for evaporation (total skin area - hair area), m2
    CONVAR <- GEOM.out[22] # area for convection (total area minus ventral area, as determined by PTCOND), m2
    R1 <- GEOM.out[23] # shape-specific core-skin radius in shortest dimension, m

    ### F_FACTOR, radiation configuration factors
    # at this stage make sure NESTYP = 0 to get correct configuration factors
    NESTYP <- 0
    F_FACTOR.out <- F_FACTOR(SHADE, NITESHAD, QSOLR, FATOBJ, NESTYP, RoNEST, R1, FGDREF, FSKREF, AREASKIN, EMISAN)

    FAVEG <- F_FACTOR.out[1] # configuration factor to vegetation
    FASKY <- F_FACTOR.out[2] # configuration factor to sky
    FAGRD <- F_FACTOR.out[3] # configuration factor to ground
    FANEST <- F_FACTOR.out[4] # configuration factor to nest wall
    # constants for infra-red exchange calculatiosn AREASKIN*CONFIG*EMISAN*SIG
    C3 <- F_FACTOR.out[5] # sky
    C4 <- F_FACTOR.out[6] # ground
    C5 <- F_FACTOR.out[7] # object
    C6 <- F_FACTOR.out[8] # vegetation (shade)
    C7 <- F_FACTOR.out[9] # nest

    ### SOLAR, solar radiation

    # solar radiation normal to sun's rays
    ZEN <- pi/180*Z # convert degrees to radians
    if(Z < 90){ # compute solar radiation on a surface normal to the direct rays of the sun
      CZ = cos(ZEN)
      QNORM = QSOLR/CZ
    }else{ # diffuse skylight only
      QNORM = QSOLR
    }

    ABSAND <- 1 - REFLFR[2] # solar absorptivity of dorsal fur (fractional, 0-1)
    ABSANV <- 1 - REFLFR[3] # solar absorptivity of ventral fur (fractional, 0-1)

    SOLAR.out <- SOLAR(AREA, ABSAND, ABSANV, ABSSB, ASILN, PCTDIF, QNORM, SHADE,
      QSOLR, FASKY, FATOBJ, FAVEG)

    QSOLAR <- SOLAR.out[1] # total (global) solar radiation (W) QSOLAR,QSDIR,QSOBJ,QSSKY,QSRSB,QSDIFF,QDORSL,QVENTR
    QSDIR <- SOLAR.out[2] # direct solar radiaton (W)
    QSOBJ <- SOLAR.out[3] # lateral diffuse solar radiation (W)
    QSSKY <- SOLAR.out[4] # diffuse solar radiation from sky (W)
    QSRSB <- SOLAR.out[5] # diffuse solar radiation reflected from substrate (W)
    QSDIFF <- SOLAR.out[6] # total diffuse solar radiation (W)
    QDORSL <- SOLAR.out[7] # dorsal direct solar radiation (W)
    QVENTR <- SOLAR.out[8] # ventral diffuse solar radiaton (W)

    ### CONV, convection

    # input
    SURFAR <- CONVAR # surface area for convection, m2 (from GEOM)
    TENV <- TA # fluid temperature (°C)

    # run subroutine
    CONV.out <- CONV(TS, TENV, NGEOM, SURFAR, FLTYPE, FURTST, D, TFA, VEL, ZFUR, BP, ELEV)

    QCONV <- CONV.out[1] # convective heat loss (W)
    HC <- CONV.out[2] # combined convection coefficient
    HCFREE <- CONV.out[3] # free convection coefficient
    HCFOR <- CONV.out[4] # forced convection coefficient
    HD <- CONV.out[5] # mass transfer coefficient
    HDFREE <- CONV.out[6] # free mass transfer coefficient
    HDFORC <- CONV.out[7] # forced mass transfer coefficient
    ANU <- CONV.out[8] # Nusselt number (-)
    RE <- CONV.out[9] # Reynold's number (-)
    GR <- CONV.out[10] # Grasshof number (-)
    PR <- CONV.out[11] # Prandlt number (-)
    RA <- CONV.out[12] # Rayleigh number (-)
    SC <- CONV.out[13] # Schmidt number (-)
    BP <- CONV.out[14] # barometric pressure (Pa)

    ### SIMULSOL, simultaneous solution of heat balance
    SIMULSOL.out <- matrix(data = 0, nrow = 2, ncol = 14) # vector to hold the SIMULSOL results for dorsal and ventral side

    # reference configuration factors
    FABUSHREF <- FABUSH # nearby bush
    FATOBJREF <- FATOBJ # nearby object
    FASKYREF <- FASKY # sky
    FAGRDREF <- FAGRD # ground
    FAVEGREF <- FAVEG # vegetation
    # repeat for each side, dorsal and ventral, of the animal

    for(S in 1:2){

      # set infrared environment
      TVEG <- TAREF # assume vegetation casting shade is at 1.2 m (reference) air temperature (°C)
      SKYIR <- C3 * (TSKY + 273.15) ^ 4 # sky infrared incoming (W)
      VEGIR <- C6 * (TVEG + 273.15) ^ 4 # vegetation infrared incomming (W)
      SKYRAD <- SKYIR + VEGIR
      #TLOCUP <- (((SKYIN) / (C3 + C6)) ^ 0.25) - 273.15
      SKYIN <- SKYRAD
      GRDIN <- C4 * (TGRD + 273.15) ^ 4 # note, MK put C4 here wherease before it was just SIG
      TLOWER <- TGRD

      # Calculating solar intensity entering fur. This will depend on whether we are calculating the fur temperature for the dorsal side or the ventral side. The dorsal side will have solar inputs from the direct beam hitting the silhouette area as well as diffuse solar scattered from the sky and objects. The ventral side will have diffuse solar scattered off the substrate.

      # Resetting config factors and solar depending on whether the dorsal side (S=1) or ventral side (S=2) is being estimated.
      if(QSOLAR > 0.0){
        if(S==1){
          FASKY <- FASKYREF/(FASKYREF+FATOBJREF+FAVEGREF)
          FATOBJ <- FATOBJREF/(FASKYREF+FATOBJREF+FAVEGREF)
          FAVEG <- FAVEGREF/(FASKYREF+FATOBJREF+FAVEGREF)
          FAGRD <- 0.0
          FABUSH <- 0.0
          if(FATOBJ == 0.0){
            QSLR <- 2*QSDIR+((QSSKY/FASKYREF)*FASKY)
          }else{
            QSLR <- 2*QSDIR+((QSSKY/FASKYREF)*FASKY)+((QSOBJ/FATOBJREF)*FATOBJ)
          }
        }else{
          FASKY <- 0.0
          FATOBJ <- 0.0
          FAVEG <- 0.0
          FAGRD <- FAGRDREF/(1 - FAGRDREF - FATOBJREF - FABUSHREF)
          FABUSH <- FABUSHREF/(1 - FAGRDREF - FATOBJREF - FABUSHREF)
          QSLR <- QVENTR/(1 - FASKYREF - FATOBJREF - FAVEGREF)
        }
      }else{
        QSLR <- 0.0
        if(S==1){
          FASKY <- FASKYREF/(FASKYREF+FATOBJREF+FAVEGREF)
          FATOBJ <- FATOBJREF/(FASKYREF+FATOBJREF+FAVEGREF)
          FAVEG <- FAVEGREF/(FASKYREF+FATOBJREF+FAVEGREF)
          FAGRD <- 0.0
          FABUSH <- 0.0
        }else{
          FASKY <- 0.0
          FATOBJ <- 0.0
          FAVEG <- 0.0
          FAGRD <- FAGRDREF/(1 - FAGRDREF - FATOBJREF - FAVEGREF)
          FABUSH <- FABUSHREF/(1 - FAGRDREF - FATOBJREF - FAVEGREF)
        }
      }

      # set fur depth and conductivity
      # index for KEFARA, the conductivity, is the average (1), front/dorsal (2), back/ventral(3) of the body part
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

      RDXDEP <- 1 # not used yet - relates to radiation through fur
      XR <- RDXDEP # not used yet - relates to radiation through fur
      X <- RDXDEP # not used yet - relates to radiation through fur
      RSKIN <- R1 # body radius (including fat), m
      RFLESH <- R1 - FATTHK # body radius flesh only (no fat), m
      RFUR <- R1 + ZL # body radius including fur, m
      D <- 2 * RFUR # diameter, m
      RRAD <- RSKIN + (XR * ZL) # effective radiation radius, m
      LEN <- ALENTH # length, m

      # Correcting volume to account for subcutaneous fat
      if(SUBQFAT == 1 & FATTHK > 0.0){
        VOL <- FLSHVL
      }

      # Getting compressed fur thermal conductivity (outputs a variable called KFURCMPRS)
      AREACND <- ATOT * PTCOND
      # CALL COMPRSKEFF # to do
      KFURCMPRS <- 1
      ZFURCOMP <- 1
      # Calculating the "Cd" variable: Qcond = Cd(Tskin-Tsub), where Cd = Conduction area*((kfur/zfur)+(ksub/subdepth))
      CD <- AREACND * ((KFURCMPRS / ZFURCOMP))

      # package up inputs
      FURVARS <- c(LEN,ZFUR,FURTHRMK,KEFF,BETARA,FURTST,ZL)
      GEOMVARS <- c(NGEOM,SUBQFAT,CONVAR,VOL,D,CONVAR,CONVSK,RFUR,RFLESH,RSKIN,XR,RRAD,ASEMAJ,BSEMIN,CSEMIN,CD)
      ENVVARS <- c(FLTYPE,TA,TS,TBUSH,TVEG,TLOWER,TSKY,TCONDSB,RH,VEL,BP,ELEV,FASKY,FABUSH,FAVEG,FAGRD,QSLR)
      TRAITS <- c(TC,AK1,AK2,EMISAN,FATTHK,FLYHR,BAREVAP,PCTBAREVAP,PCTEYES)

      # set IPT, the geometry assumed in SIMULSOL: 1 = cylinder, 2 = sphere, 3 = ellipsoid
      if(NGEOM %in% c(1,3,5)){
        IPT <- 1
      }
      if(NGEOM == 2){
        IPT <- 2
      }
      if(NGEOM == 4){
        IPT <- 3
      }

      # call SIMULSOL
      SIMULSOL.out[S,] <- SIMULSOL(DIFTOL, IPT, FURVARS, GEOMVARS, ENVVARS, TRAITS, TFA, SKINW, TS)
    }
    TSKINMAX <- max(SIMULSOL.out[1,2], SIMULSOL.out[2,2])
    ### ZBRENT and RESPFUN

    # Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(FASKY+FAVEG+FATOBJ))+(ventral value*FAGRD)
    GEND <- SIMULSOL.out[1, 5]
    GENV <- SIMULSOL.out[2, 5]
    DMULT <- FASKYREF + FAVEGREF + FATOBJ
    VMULT <- 1 - DMULT # Assume that reflectivity of veg below = ref of soil so VMULT left as 1 - DMULT
    X <- GEND * DMULT + GENV * VMULT # weighted estimate of metabolic heat generation

    # reset configuration factors
    FABUSH <- FABUSHREF # nearby bush
    FATOBJ <- FATOBJREF # nearby object
    FASKY <- FASKYREF # sky
    FAGRD <- FAGRDREF # ground
    FAVEG <- FAVEGREF # vegetation

    # lung temperature and temperature of exhaled air
    TLUNG <- (TC + (SIMULSOL.out[1, 2] + SIMULSOL.out[1, 2]) * 0.5) * 0.5 # average of skin and core
    TAEXIT <- min(TA + DELTAR, TLUNG) # temperature of exhaled air, °C

    # now guess for metabolic rate that balances the heat budget while allowing metabolic rate
    # to remain at or above QBASAL, via 'shooting method' ZBRENT
    QMIN <- QBASAL
    if(TA < TC & TSKINMAX < TC){
     QM1 <- QBASAL * 2 * -1
     QM2 <- QBASAL * 50
    }else{
     QM1 <- QBASAL * 50* -1
     QM2 <- QBASAL * 2
    }
    QSUM <- X
    TOL <- AMASS * 0.01

    ZBRENT.in <- c(TA, O2GAS, N2GAS, CO2GAS, BP, QMIN, RQ, TLUNG, GMASS, EXTREF, RH,
      RELXIT, TIMACT, TAEXIT, QSUM, PANT)
    # call ZBRENT subroutine which calls RESPFUN
    ZBRENT.out <- ZBRENT(QM1, QM2, TOL, ZBRENT.in)
    colnames(ZBRENT.out) <- c("RESPFN","QRESP","GEVAP", "PCTO2", "PCTN2", "PCTCO2", "RESPGEN", "O2STP", "O2MOL1", "N2MOL1", "AIRML1", "O2MOL2", "N2MOL2", "AIRML2", "AIRVOL")

    QGEN <- ZBRENT.out[7]
    GMULTLAST <- GMULT
    AK1LAST <- AK1
    TCLAST <- TC
    PANTLAST <- PANT
    SKINWLAST <- SKINW

    if(GMULT < GMULTMAX){
      GMULT <- GMULT + UNCURL
    }else{
      GMULT <- GMULTMAX
      if(AK1 < AKMAX){
        AK1 <- AK1 + AK1inc
      }else{
        AK1 <- AKMAX
        if(TC < TCMAX){
          TC <- TC + RAISETC
          Q10mult <- Q10^((TC - TCREF)/10)
          QBASAL = QBASREF * Q10mult
        }else{
          TC <- TCMAX
          Q10mult <- Q10^((TC - TCREF)/10)
          QBASAL = QBASREF * Q10mult
          if(PANT < PANTMAX){
            PANT <- PANT + PANTING
            PANTSTEP <- PANTSTEP + 1
            #PANT <- round(PANTMAX - (PANTMAX - 1) * exp(-0.02 / (PANTMAX / 10) * PANTSTEP), 1)
          }else{
            PANT <- PANTMAX
            SKINW <- SKINW + SWEAT
            if(SKINW > MXWET | SWEAT == 0){
              SKINW <- MXWET
                break
            }
          }
        }
      }
    }
  }

  HTOVPR <- 2.5012E+06 - 2.3787E+03 * TA
  SWEAT.G.H <- (SIMULSOL.out[1,6] + SIMULSOL.out[2,6]) * 0.5 / HTOVPR * 1000 * 3600
  EVAP.G.H <- ZBRENT.out[3] * 3600 + SWEAT.G.H
  endo.out <- as.matrix(cbind(TCLAST, TLUNG, SIMULSOL.out[1,1], SIMULSOL.out[2,1], SIMULSOL.out[1,2], SIMULSOL.out[2,2], SIMULSOL.out[1,3], SIMULSOL.out[2,3], SIMULSOL.out[1,4], SIMULSOL.out[2,4], SIMULSOL.out[1,5], SIMULSOL.out[2,5], SIMULSOL.out[1,6], SIMULSOL.out[2,6], SIMULSOL.out[1,7], SIMULSOL.out[2,7], SIMULSOL.out[1,8], SIMULSOL.out[2,8], SIMULSOL.out[1,9], SIMULSOL.out[2,9], SIMULSOL.out[1,10], SIMULSOL.out[2,10], SIMULSOL.out[1,11], SIMULSOL.out[2,11], SIMULSOL.out[1,12], SIMULSOL.out[2,12], SIMULSOL.out[1,13], SIMULSOL.out[2,13], SIMULSOL.out[1,14], SIMULSOL.out[2,14], ZBRENT.out, GMULTLAST, PANTLAST, SKINWLAST, SWEAT.G.H, EVAP.G.H, AK1LAST, TA, TGRD, TCONDSB, TSKY, VEL, RH, QSOLR))
  colnames(endo.out) <- c("TC", "TLUNG", "TFA_D", "TFA_V", "TSKIN_D", "TSKIN_V", "QCONV_D", "QCONV_V", "QCOND_D", "QCOND_V", "QGENNET_D", "QGENNET_V", "QSEVAP_D", "QSEVAP_V", "QRAD_D", "QRAD_V", "QSLR_D", "QSLR_V", "QRSKY_D", "QRSKY_V", "QRBSH_D", "QRBSH_V", "QRVEG_D", "QRVEG_V", "QRGRD_D", "QRGRD_V", "NTRY_D", "NTRY_V", "SUCCESS_D", "SUCCESS_V", "RESPFN","QRESP","GEVAP", "PCTO2", "PCTN2", "PCTCO2", "RESPGEN", "O2STP", "O2MOL1", "N2MOL1", "AIRML1", "O2MOL2", "N2MOL2", "AIRML2", "AIRVOL", "GMULT", "PANT", "SKINW", "SWEAT.G.H", "EVAP.G.H", "AK", "TA", "TGRD", "TCONDSB", "TSKY", "VEL", "RH", "QSOLR")

  return(endo.out)
}
