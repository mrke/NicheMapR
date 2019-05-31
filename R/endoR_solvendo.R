#' endoR
#'
#' Endotherm model of NicheMapR
#' @encoding UTF-8
#' @param AMASS = 1, # kg
#' @param NGEOM = 4, # cylinder (ngeom = 1), sphere (ngeom = 2) and ellipsoid (ngeom = 4). If a truncated cone (5) or ellipsoidal cylinder (3), we will use the cylinder equations (ngeom=1).
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
#' @export
#' @details
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#' \code{DIFTOL}{ = 0.001, error tolerance for SIMULSOL (°C)}\cr\cr
#'
#' \strong{ Environment:}\cr\cr
#' \code{TAREF}{ = TA, air temperature at reference height (°C)}\cr\cr
#' \code{ELEV}{ = 0, elevation (m)}\cr\cr
#' \code{ABSSB}{ = 0.8, solar absorptivity of substrate (fractional, 0-1)}\cr\cr
#' \code{FLTYPE}{ = 0, FLUID TYPE: 0 = AIR; 1 = FRESH WATER; 2 = SALT WATER - need's to be looked at - only invoked in main program when the dive table is set up}\cr\cr
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
#' \code{MAXPTVEN}{ = 0, maxium fraction of surface area that is ventral (fractional, 0-1)}\cr\cr
#' \code{AWING}{ = 0, area of wing, to do}\cr\cr
#' \code{PTCOND}{ = 0, \% of body area touching the substrate}\cr\cr
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
#' \code{TC}{ }\cr\cr
#' \code{TFA_D}{ }\cr\cr
#' \code{TFA_V}{ }\cr\cr
#' \code{TSKIN_D}{ }\cr\cr
#' \code{TSKIN_V}{ }\cr\cr
#' \code{QCONV_D}{ }\cr\cr
#' \code{QCONV_V}{ }\cr\cr
#' \code{QCOND_D}{ }\cr\cr
#' \code{QCOND_V}{ }\cr\cr
#' \code{QGENNET_D}{ }\cr\cr
#' \code{QGENNET_V}{ }\cr\cr
#' \code{QSEVAP_D}{ }\cr\cr
#' \code{QSEVAP_V}{ }\cr\cr
#' \code{QRAD_D}{ }\cr\cr
#' \code{QRAD_V}{ }\cr\cr
#' \code{QSLR_D}{ }\cr\cr
#' \code{QSLR_V}{ }\cr\cr
#' \code{QRSKY_D}{ }\cr\cr
#' \code{QRSKY_V}{ }\cr\cr
#' \code{QRBSH_D}{ }\cr\cr
#' \code{QRBSH_V}{ }\cr\cr
#' \code{QRVEG_D}{ }\cr\cr
#' \code{QRVEG_V}{ }\cr\cr
#' \code{QRGRD_D}{ }\cr\cr
#' \code{QRGRD_V}{ }\cr\cr
#' \code{NTRY_D}{ }\cr\cr
#' \code{NTRY_V}{ }\cr\cr
#' \code{SUCCESS_D}{ }\cr\cr
#' \code{SUCCESS_V}{ }\cr\cr
#' \code{RESPFN}{ }\cr\cr
#' \code{QRESP}{ }\cr\cr
#' \code{GEVAP}{ }\cr\cr
#' \code{PCTO2}{ }\cr\cr
#' \code{PCTN2}{ }\cr\cr
#' \code{PCTCO2}{ }\cr\cr
#' \code{RESPGEN}{ }\cr\cr
#' \code{O2STP}{ }\cr\cr
#' \code{O2MOL1}{ }\cr\cr
#' \code{N2MOL1}{ }\cr\cr
#' \code{AIRML1}{ }\cr\cr
#' \code{O2MOL2}{ }\cr\cr
#' \code{N2MOL2}{ }\cr\cr
#' \code{AIRML2}{ }\cr\cr
#' \code{AIRVOL}{ }\cr\cr
#' \code{GMULT}{ }\cr\cr
#' \code{SKINW}{ }\cr\cr
#' \code{SWEAT.G.H}{ }\cr\cr
#' \code{EVAP.G.H}{ }\cr\cr
#' \code{EXTREF}{ }\cr\cr
#' \code{AK}{ }\cr\cr
#' \code{TA}{ }\cr\cr
#' \code{TGRD}{ }\cr\cr
#' \code{TCONDSB}{ }\cr\cr
#' \code{TSKY}{ }\cr\cr
#' \code{VEL}{ }\cr\cr
#' \code{RH}{ }\cr\cr
#' \code{QSOLR}{ }\cr\cr
#' @examples
#' library(NicheMapR)
#' dstart <- "02/01/2017"
#' dfinish <- "30/12/2017"
endoR_solvendo <- function(
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
  MAXPTVEN = 0, # maxium fraction of surface area that is ventral (fractional, 0-1)
  AWING = 0, # area of wing, to do
  PTCOND = 0, # % of body area touching the substrate

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
#
#   TA = 50 # air temperature at local height (°C)
#   TAREF = TA # air temeprature at reference height (°C)
#   TGRD = TA # ground temperature (°C)
#   TSKY = TA # sky temperature (°C)
#   VEL = 0.1 # wind speed (m/s)
#   RH = 5 # relative humidity (%)
#   QSOLR = 0 # solar radiation, horizontal plane (W/m2)
#   Z = 20 # zenith angle of sun (degrees from overhead)
#   ELEV = 0 # elevation (m)
#   ABSSB = 0.8 # solar absorptivity of substrate (fractional, 0-1)
#
#   # other environmental variables
#   FLTYPE = 0 # FLUID TYPE: 0 = AIR; 1 = FRESH WATER; 2 = SALT WATER - need's to be looked at - only invoked in main program when the dive table is set up
#   TCONDSB = TGRD # surface temperature for conduction (°C)
#   TBUSH = TA # bush temperature (°C)
#   BP = -1 # Pa, negatve means elevation is used
#   O2GAS = 20.95 # oxygen concentration of air (%)
#   N2GAS = 79.02 # nitrogen concetration of air (%)
#   CO2GAS = 0.03 # carbon dioxide concentration of air (%)
#   PCTDIF = 0.15 # proportion of solar radiation that is diffuse (fractional, 0-1)
#
#   # BEHAVIOUR
#
#   SHADE = 0 # shade level (%)
#   NITESHAD = 0 # flag for if animal is behaviourally seeking shade for warmth at night - remove?
#   FLYHR = 0 # is flight occuring this hour? (imposes forced evaporative loss)
#   UNCURL = 1 # allows the animal to uncurl to GMULTMAX, the value being the increment GMULT is increased per iteration
#   RAISETC = 1 # turns on core temperature elevation, the value being the increment by which TC is increased per iteration
#   SWEAT = 0.25 # turns on sweating, the value being the increment by which SKINW is increased per iteration
#   MXWET = 100 # maximum surface area that can be wet (%)
#   AK1inc = 0.5 # turns on thermal conductivity increase (W/mK), the value being the increment by which AK1 is increased per iteration
#   AKMAX = 2.8 # maximum flesh conductivity (W/mK)
#   PANT = 1 # multiplier on breathing rate to simulate panting (-)
#   PANTING = 0.1 # increment for multiplier on breathing rate to simulate panting (-)
#
#   # MORPHOLOGY
#
#   # geometry
#   AMASS = 1 # kg
#   ANDENS = 1000 # kg/m3
#   SUBQFAT = 0 # is subcutaneous fat present? (0 is no, 1 is yes)
#   FATPCT = 20 # % body fat
#   NGEOM = 4 # cylinder (ngeom = 1), sphere (ngeom = 2) and ellipsoid (ngeom = 4). If a truncated cone (5) or ellipsoidal cylinder (3), we will use the cylinder equations (ngeom=1).
#   GMREF = 3 # initial ratio between long and short axis (-)
#   GMULT = GMREF # current ratio between long and short axis (-)
#   GMULTMAX = GMREF # max possible ratio between long and short axis (-)
#   MAXPTVEN = 0 # maxium fraction of surface area that is ventral (fractional, 0-1)
#   AWING = 0 # area of wing, to do
#   PTCOND = 0 # % of body area touching the substrate
#
#   # fur properties
#   FURTHRMK = 0 # user-specified fur thermal conductivity (W/mK), not used if 0
#   DHAIRD = 30E-06 # hair diameter, dorsal (m)
#   DHAIRV = 30E-06 # hair diameter, ventral (m)
#   LHAIRD = 23.9E-03 # hair length, dorsal (m)
#   LHAIRV = 23.9E-03 # hair length, ventral (m)
#   ZFURD = 2E-03 # fur depth, dorsal (m)
#   ZFURV = 2E-03 # fur depth, ventral (m)
#   RHOD = 3000E+04 # hair density, dorsal (1/m2)
#   RHOV = 3000E+04 # hair density, ventral (1/m2)
#   REFLD = 0.2  # fur reflectivity dorsal (fractional, 0-1)
#   REFLV = 0.2  # fur reflectivity ventral (fractional, 0-1)
#
#   # radiation exchange
#   EMISAN = 0.99 # animal emissivity (-)
#   FATOBJ = 0 # configuration factor to nearby object
#   FABUSH = 0 # this is for veg below/around animal (at TALOC)
#   FGDREF = 0.5 # reference configuration factor to ground
#   FSKREF = 0.5 # configuration factor to sky
#
#   # nest properties
#   NESTYP = 0 # for nest calculations, to do
#   RoNEST = 0 # for nest calculations, to do
#
#   # PHYSIOLOGY
#
#   # thermal
#   TC = 37 # core temperature (°C)
#   TCMAX = 45 # maximum core temperature (°C)
#   AK1 = 0.9 # initial thermal conductivity of flesh (0.412 - 2.8 W/mC)
#   AK2 = 0.230# conductivity of fat (W/mK)
#
#   # evaporation
#   SKINW = 0 # part of the skin surface that is wet (%)
#   BAREVAP = 0 # is evaporation partly from bare skin? (0 = no, 1 = yes, % defined with PCTSKINEVAP)
#   PCTBAREVAP = 0 # surface area for evaporation that is skin, e.g. licking paws (%)
#   PCTEYES = 0 # surface area made up by the eye (%) - make zero if sleeping
#   DELTAR = 0 # offset between air temeprature and breath (°C)
#   RELXIT = 100 # relative humidity of exhaled air, %
#
#   # metabolism/respiration
#   QBASAL = (70 * AMASS ^ 0.75) * (4.185 / (24 * 3.6)) # basal heat generation (W)
#   TIMACT = 1 # multiplier on metabolic rate for activity costs
#   RQ = 0.80 # respiratory quotient (fractional, 0-1)
#   EXTREF = 20 # O2 extraction efficiency (%)
#   PANTMAX = 2 # maximum breathing rate multiplier to simulate panting (-)
#   Q10 = 1 # Q10 factor for adjusting BMR for TC
#
#   # initial conditions
#   TS = TC - 3 # skin temperature (°C)
#   TFA = TA # fur/air interface temperature (°C)
#
#   # other model settings
#   DIFTOL = 0.001 # tolerance for SIMULSOL

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
  TVEG <- TA
  SOLVENDO.input <- c(QGEN, QBASAL,TA, GMULTMAX, GMREF, GMULT, DHAIRD, DHAIRV, LHAIRD, LHAIRV, ZFURD, ZFURV, RHOD, RHOV, REFLD, REFLV, MAXPTVEN, NGEOM,EMISAN,FATOBJ,FSKREF,FGDREF,NESTYP,PCTDIF,ABSSB,AWING,FLTYPE,ELEV,BP,NITESHAD,SHADE,QSOLR,RoNEST,Z,VEL, TS,TFA,FABUSH,FURTHRMK,RH,TCONDSB,TBUSH,TC,PCTBAREVAP,FLYHR,BAREVAP,AK1,AK2,PCTEYES,DIFTOL,SKINW,TSKY,TVEG,TAREF,DELTAR,RQ, TIMACT, O2GAS, N2GAS, CO2GAS,RELXIT,PANT,EXTREF,UNCURL,AKMAX,AK1inc,TCMAX,RAISETC,TCREF,Q10,QBASREF,PANTMAX,PANTSTEP,MXWET, SWEAT,TGRD,AMASS,ANDENS,SUBQFAT,FATPCT,PTCOND,PANTING)
  endo.out <- SOLVENDO(SOLVENDO.input)
  colnames(endo.out) <- c("TC", "TLUNG", "TFA_D", "TFA_V", "TSKIN_D", "TSKIN_V", "QCONV_D", "QCONV_V", "QCOND_D", "QCOND_V", "QGENNET_D", "QGENNET_V", "QSEVAP_D", "QSEVAP_V", "QRAD_D", "QRAD_V", "QSLR_D", "QSLR_V", "QRSKY_D", "QRSKY_V", "QRBSH_D", "QRBSH_V", "QRVEG_D", "QRVEG_V", "QRGRD_D", "QRGRD_V", "NTRY_D", "NTRY_V", "SUCCESS_D", "SUCCESS_V", "RESPFN","QRESP","GEVAP", "PCTO2", "PCTN2", "PCTCO2", "RESPGEN", "O2STP", "O2MOL1", "N2MOL1", "AIRML1", "O2MOL2", "N2MOL2", "AIRML2", "AIRVOL", "GMULT", "PANT", "SKINW", "SWEAT.G.H", "EVAP.G.H", "EXTREF", "AK", "TA", "TGRD", "TCONDSB", "TSKY", "VEL", "RH", "QSOLR")
  endo.out
  dyn.unload('C:/Users/mrke/OneDrive - The University of Melbourne/Documents/R/win-library/3.5/NicheMapR/libs/win/x64/SOLVENDO.dll')
  return(endo.out)
}
