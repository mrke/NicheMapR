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
#' @param SHAPE = 4, # shape, 1 is cylinder, 2 is sphere, 3 is plate, 4 is ellipsoid
#' @param SHAPE_B_REF = 3, # initial ratio between long and short axis (-)
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
#' @usage endoR(AMASS = 1, SHAPE = 4, SHAPE_B_REF = 3, FURTHRMK = 0, ZFURD = 2E-03, ZFURV = 2E-03, TC = 37, TCMAX = 45, TA = 20, TGRD = TA, TSKY = TA, VEL = 0.1, RH = 5, QSOLR = 0, Z = 20, SHADE = 0, NITESHAD = 0,...)
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
#' \code{UNCURL}{ = 1, allows the animal to uncurl to SHAPE_B_MAX, the value being the increment SHAPE_B is increased per iteration}\cr\cr
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
#' \code{SHAPE_B}{ = SHAPE_B_REF, current ratio between long and short axis (-)}\cr\cr
#' \code{SHAPE_B_MAX}{ = SHAPE_B_REF, max possible ratio between long and short axis (-)}\cr\cr
#' \code{SHAPE_C}{ = SHAPE_B, current ratio of length:height (plate)}\cr\cr
#' \code{MAXPTVEN}{ = 0.5, maxium fraction of surface area that is ventral (fractional, 0-1)}\cr\cr
#' \code{PCOND}{ = 0, fraction of surface area that is touching the substrate (fractional, 0-1)}\cr\cr
#' \code{MAXPCOND}{ = 0, maximum fraction of surface area that is touching the substrate (fractional, 0-1)}\cr\cr
#' \code{SAMODE}{ = 0, if 0, uses surface area for SHAPE geometry, if 1, uses bird skin surface area allometry from Walsberg & King. 1978. JEB 76:185–189, if 2 uses mammal surface area from Stahl 1967.J. App. Physiol. 22, 453–460.}\cr\cr
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
#' \code{ZFURCOMP}{ = ZFURV, depth of compressed fur (for conduction) (m)}\cr\cr
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
#' \code{SKINW}{ = 0.5, # part of the skin surface that is wet (\%)}\cr\cr
#' \code{FURWET}{ = 0, # Area of fur/feathers that is wet after rain (\%)}\cr\cr
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
#'
#' treg variables (thermoregulatory response):
#' \itemize{
#' \item 1 TC - core temperature (°C)
#' \item 2 TLUNG - lung temperature (°C)
#' \item 3 TSKIN_D  - dorsal skin temperature (°C)
#' \item 4 TSKIN_V - ventral skin temperature (°C)
#' \item 5 TFA_D - dorsal fur-air interface temperature (°C)
#' \item 6 TFA_V - ventral fur-air interface temperature (°C)
#' \item 7 SHAPE_B - current ratio between long and short axis due to postural change (-)
#' \item 8 PANT - breathing rate multiplier (-)
#' \item 9 SKINWET - part of the skin surface that is wet (\%)
#' \item 10 K_FLESH - thermal conductivity of flesh (W/mC)
#' \item 11 K_FUR - thermal conductivity of flesh (W/mC)
#' \item 12 K_FUR_D - thermal conductivity of dorsal fur (W/mC)
#' \item 13 K_FUR_V - thermal conductivity of ventral fur (W/mC)
#' \item 14 K_COMPFUR - thermal conductivity of compressed fur (W/mC)
#' \item 15 Q10 - Q10 multiplier on metabolic rate (-)
#' }
#' morph variables (morphological traits):
#' \itemize{
#' \item 1 AREA - total outer surface area (m2)
#' \item 2 VOLUME - total volume (m3)
#' \item 3 CHAR_DIM  - characteristic dimension for convection (m)
#' \item 4 MASS_FAT - fat mass (kg)
#' \item 5 FAT_THICK - thickness of fat layer (m)
#' \item 6 FLESH_VOL - flesh volume (m3)
#' \item 7 LENGTH - length (m)
#' \item 8 WIDTH - width (m)
#' \item 9 HEIGHT - height (m)
#' \item 10 DIAM_FLESH - diameter, core to skin (m)
#' \item 11 DIAM_FUR - diameter, core to fur (m)
#' \item 12 AREA_SIL - silhouette area (m2)
#' \item 13 AREA_SILN - silhouette area normal to sun's rays (m2)
#' \item 14 AREA_ASILP - silhouette area parallel to sun's rays (m2)
#' \item 15 AREA_SKIN - total skin area (m2)
#' \item 16 AREA_SKIN_EVAP - skin area available for evaporation (m2)
#' \item 17 AREA_CONV - area for convection (m2)
#' \item 18 AREA_COND - area for conduction (m2)
#' \item 19 F_SKY - configuration factor to sky (-)
#' \item 20 F_GROUND - configuration factor to ground (-)
#' }
#' enbal variables (energy balance):
#' \itemize{
#' \item 1 QSOL - solar radiation absorbed (W)
#' \item 2 QIRIN - longwave (infra-red) radiation absorbed (W)
#' \item 3 QMET  - characteristic dimension for convection (W)
#' \item 4 QEVAP - evaporation (W)
#' \item 5 QIROUT - longwave (infra-red) radiation lost (W)
#' \item 6 QCONV - convection (W)
#' \item 7 QCOND - conduction (W)
#' \item 8 ENB - energy balance (W)
#' \item 9 NTRY - iterations required for a solution (-)
#' \item 10 SUCCESS - was a solution found (0=no, 1=yes)
#' }
#' masbal variables (mass exchanges):
#' \itemize{
#' \item 1 AIR_L - breating rate (L/h)
#' \item 2 O2_L - oxgyen consumption rate (L/h)
#' \item 3 H2OResp_g - respiratory water loss (g/h)
#' \item 4 H2OCut_g - cutaneous water loss (g/h)
#' \item 5 O2_mol_in - oxygen inhaled (mol/h)
#' \item 6 O2_mol_out - oxygen expelled (mol/h)
#' \item 7 N2_mol_in - nitrogen inhaled (mol/h)
#' \item 8 N2_mol_out - nitrogen expelled (mol/h)
#' \item 9 AIR_mol_in - air inhaled (mol/h)
#' \item 10 AIR_mol_out - air expelled (mol/h)
#' }
#' @examples
#' library(NicheMapR)
#' # environment
#' TAs <- seq(0, 50, 2) # air temperatures (°C)
#' VEL <- 0.002 # wind speed (m/s)
#' RH <- 10 # relative humidity (\%)
#' QSOLR <- 100 # solar radiation (W/m2)
#'
#' # core temperature
#' TC <- 38 # core temperature (deg C)
#' TCMAX <- 43 # maximum core temperature (°C)
#' RAISETC <- 0.25 # increment by which TC is elevated (°C)
#'
#' # size and shape
#' AMASS <- 0.0337 # mass (kg)
#' SHAPE_B_REF <- 1.1 # start off near to a sphere (-)
#' SHAPE_B_MAX <- 5 # maximum ratio of length to width/depth
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
#' SKINW <- 0.1 # base skin wetness (\%)
#' MXWET <- 20 # maximum skin wetness (\%)
#' SWEAT <- 0.25 # intervals by which skin wetness is increased (\%)
#' Q10 <- 2 # Q10 effect of body temperature on metabolic rate
#' QBASAL <- 10 ^ (-1.461 + 0.669 * log10(AMASS * 1000)) # basal heat generation (W) (bird formula from McKechnie and Wolf 2004 Phys. & Biochem. Zool. 77:502-521)
#' DELTAR <- 5 # offset between air temeprature and breath (°C)
#' EXTREF <- 15 # O2 extraction efficiency (\%)
#' PANTING <- 0.1 # turns on panting, the value being the increment by which the panting multiplier is increased up to the maximum value, PANTMAX
#' PANTMAX <- 3# maximum panting rate - multiplier on air flow through the lungs above that determined by metabolic rate
#'
#' ptm <- proc.time() # start timing
#' endo.out <- lapply(1:length(TAs), function(x){endoR_devel(TA = TAs[x], QSOLR = QSOLR, VEL = VEL, TC = TC, TCMAX = TCMAX, RH = RH, AMASS = AMASS, SHAPE_B_REF = SHAPE_B_REF, SHAPE_B_MAX = SHAPE_B_MAX, SKINW = SKINW, SWEAT = SWEAT, MXWET = MXWET, Q10 = Q10, QBASAL = QBASAL, DELTAR = DELTAR, DHAIRD = DHAIRD, DHAIRV = DHAIRV, LHAIRD = LHAIRD, LHAIRV = LHAIRV, ZFURD = ZFURD, ZFURV = ZFURV, RHOD = RHOD, RHOV = RHOV, REFLD = REFLD, RAISETC = RAISETC, PANTING = PANTING, PANTMAX = PANTMAX, EXTREF = EXTREF)}) # run endoR across environments
#' proc.time() - ptm # stop timing
#'
#' endo.out1 <- do.call("rbind", lapply(endo.out, data.frame)) # turn results into data frame
#' treg <- endo.out1[, grep(pattern = "treg", colnames(endo.out1))]
#' colnames(treg) <- gsub(colnames(treg), pattern = "treg.", replacement = "")
#' morph <- endo.out1[, grep(pattern = "morph", colnames(endo.out1))]
#' colnames(morph) <- gsub(colnames(morph), pattern = "morph.", replacement = "")
#' enbal <- endo.out1[, grep(pattern = "enbal", colnames(endo.out1))]
#' colnames(enbal) <- gsub(colnames(enbal), pattern = "enbal.", replacement = "")
#' masbal <- endo.out1[, grep(pattern = "masbal", colnames(endo.out1))]
#' colnames(masbal) <- gsub(colnames(masbal), pattern = "masbal.", replacement = "")
#'
#' QGEN <- enbal$QMET # metabolic rate (W)
#' H2O <- masbal$H2OResp_g + masbal$H2OCut_g # g/h water evaporated
#' TFA_D <- treg$TFA_D # dorsal fur surface temperature
#' TFA_V <- treg$TFA_V # ventral fur surface temperature
#' TskinD <- treg$TSKIN_D # dorsal skin temperature
#' TskinV <- treg$TSKIN_V # ventral skin temperature
#' TCs <- treg$TC # core temperature
#'
#' par(mfrow = c(2, 2))
#' par(oma = c(2, 1, 2, 2) + 0.1)
#' par(mar = c(3, 3, 1.5, 1) + 0.1)
#' par(mgp = c(2, 1, 0))
#' plot(QGEN ~ TAs, type = 'l', ylab = 'metabolic rate, W', xlab = 'air temperature, deg C', ylim = c(0.2, 1.2))
#' plot(H2O ~ TAs, type = 'l', ylab = 'water loss, g/h', xlab = 'air temperature, deg C', ylim = c(0, 1.5))
#' points(masbal$H2OResp_g ~ TAs, type = 'l', lty = 2)
#' points(masbal$H2OCut_g ~ TAs, type = 'l', lty = 2, col = 'blue')
#' legend(x = 3, y = 1.5, legend = c("total", "respiratory", "cutaneous"), col = c("black", "black", "blue"), lty = c(1, 2, 2), bty = "n")
#' plot(TFA_D ~ TAs, type = 'l', col = 'grey', ylab = 'temperature, deg C', xlab = 'air temperature, deg C', ylim = c(10, 50))
#' points(TFA_V ~ TAs, type = 'l', col = 'grey', lty = 2)
#' points(TskinD ~ TAs, type = 'l', col = 'orange')
#' points(TskinV ~ TAs, type = 'l', col = 'orange', lty = 2)
#' points(TCs ~ TAs, type = 'l', col = 'red')
#' legend(x = 30, y = 33, legend = c("core", "skin dorsal", "skin ventral", "feathers dorsal", "feathers ventral"), col = c("red", "orange", "orange", "grey", "grey"), lty = c(1, 1, 2, 1, 2), bty = "n")
#' plot(masbal$AIR_L * 1000 / 60 ~ TAs, ylim=c(0,250),  lty = 1, xlim=c(-5,50), ylab = "ml / min", xlab=paste("air temperature (deg C)"), type = 'l')
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
  UNCURL = 1, # allows the animal to uncurl to SHAPE_B_MAX, the value being the increment SHAPE_B is increased per iteration
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
  SHAPE = 4, # shape, 1 is cylinder, 2 is sphere, 3 is plate, 4 is ellipsoid
  SHAPE_B_REF = 3, # initial ratio between long and short axis (-)
  SHAPE_B = SHAPE_B_REF, # current ratio between long and short axis (-)
  SHAPE_B_MAX = SHAPE_B_REF, # max possible ratio between long and short axis (-)
  SHAPE_C = SHAPE_B, # current ratio of length:height (plate)
  MAXPTVEN = 0.5, # maxium fraction of surface area that is ventral (fractional, 0-1)
  PCOND = 0, # fraction of surface area that is touching the substrate (fractional, 0-1)
  MAXPCOND = 0, # maximum fraction of surface area that is touching the substrate (fractional, 0-1)
  SAMODE = 0, # if 0, uses surface area for SHAPE parameter geometry, if 1, uses bird skin surface area allometry from Walsberg & King. 1978. JEB 76:185–189, if 2 uses mammal surface area from Stahl 1967.J. App. Physiol. 22, 453–460.
  ORIENT = 0, # if 1 = normal to sun's rays (heat maximising), if 2 = parallel to sun's rays (heat minimising), or 0 = average

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
  ZFURCOMP = ZFURV, # depth of compressed fur (for conduction) (m)

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
  SKINW = 0.5, # part of the skin surface that is wet (%)
  FURWET = 0, # part of the fur/feathers that is wet after rain (%)
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
    SHAPE_B_MAX <- SHAPE_B # can't change posture, so max multiplier of dimension set to current value
  }
  Q10mult <- 1
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
    SHAPE_B <- SHAPE_B_MAX
  }

  while(QGEN < QBASAL){

    ### IRPROP, infrared radiation properties of fur

    # call the IR properties subroutine
    IRPROP.out <- IRPROP(TA, SHAPE_B_MAX, SHAPE_B_REF, SHAPE_B, DHAIRD, DHAIRV, LHAIRD, LHAIRV, ZFURD, ZFURV, RHOD, RHOV, REFLD, REFLV, MAXPTVEN, ZFURCOMP)

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
    KFURCMPRS <- IRPROP.out[27] # effictive thermal conductivity of compressed ventral fur (W/mK)

    ### GEOM, geometry

    # input
    DHARA <- DHAR[1] # fur diameter, mean (m) (from IRPROP)
    RHOARA <- RHOAR[1] # hair density, mean (1/m2) (from IRPROP)
    ZFUR <- ZZFUR[1] # fur depth, mean (m) (from IRPROP)

    # call the subroutine
    GEOM.out <- GEOM(AMASS, ANDENS, FATPCT, SHAPE, ZFUR, SUBQFAT, SHAPE_B, SHAPE_B_REF, SHAPE_C, DHARA, RHOARA, PCOND, SAMODE, ORIENT)

    # output
    VOL <- GEOM.out[1] # volume, m3
    D <- GEOM.out[2] # characteristic dimension for convection, m
    MASFAT <- GEOM.out[3] # mass body fat, kg
    VOLFAT <- GEOM.out[4] # volume body fat, m3
    ALENTH <- GEOM.out[5] # length, m
    AWIDTH <- GEOM.out[6] # width, m
    AHEIT <- GEOM.out[7] # height, m
    ATOT <- GEOM.out[8] # total area at fur/feathers-air interface, m2
    ASIL <- GEOM.out[9] # silhouette area to use in solar calcs, m2 may be normal, parallel or average set via ORIENT
    ASILN <- GEOM.out[10] # silhouette area normal to sun, m2
    ASILP <- GEOM.out[11] # silhouette area parallel to sun, m2
    GMASS <- GEOM.out[12] # mass, g
    AREASKIN <- GEOM.out[13] # area of skin, m2
    FLSHVL <- GEOM.out[14] # flesh volume, m3
    FATTHK <- GEOM.out[15] # fat layer thickness, m
    ASEMAJ <- GEOM.out[16] # semimajor axis length, m
    BSEMIN <- GEOM.out[17] # b semiminor axis length, m
    CSEMIN <- GEOM.out[18] # c semiminor axis length, m (currently only prolate spheroid)
    CONVSK <- GEOM.out[19] # area of skin for evaporation (total skin area - hair area), m2
    CONVAR <- GEOM.out[20] # area for convection (total area minus ventral area, as determined by PCOND), m2
    R1 <- GEOM.out[21] # shape-specific core-skin radius in shortest dimension, m
    R2 <- GEOM.out[22] # shape-specific core-fur radius in shortest dimension, m

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

    ABSAND <- 1 - REFLD # solar absorptivity of dorsal fur (fractional, 0-1)
    ABSANV <- 1 - REFLV # solar absorptivity of ventral fur (fractional, 0-1)

    SOLAR.out <- SOLAR(ATOT, ABSAND, ABSANV, ABSSB, ASILN, PCTDIF, QNORM, SHADE,
      QSOLR, FASKY, FATOBJ, FAVEG)

    QSOLAR <- SOLAR.out[1] # total (global) solar radiation (W) QSOLAR,QSDIR,QSOBJ,QSSKY,QSRSB,QSDIFF,QDORSL,QVENTR
    QSDIR <- SOLAR.out[2] # direct solar radiaton (W)
    QSOBJ <- SOLAR.out[3] # lateral diffuse solar radiation (W)
    QSSKY <- SOLAR.out[4] # diffuse solar radiation from sky (W)
    QSRSB <- SOLAR.out[5] # diffuse solar radiation reflected from substrate (W)
    QSDIFF <- SOLAR.out[6] # total diffuse solar radiation (W)
    QDORSL <- SOLAR.out[7] # total dorsal solar radiation (W)
    QVENTR <- SOLAR.out[8] # total ventral solar radiaton (W)

    ### CONV, convection

    # input
    SURFAR <- CONVAR # surface area for convection, m2 (from GEOM)
    TENV <- TA # fluid temperature (°C)

    # run subroutine
    CONV.out <- CONV(TS, TENV, SHAPE, SURFAR, FLTYPE, FURTST, D, TFA, VEL, ZFUR, BP, ELEV)

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
    SIMULSOL.out <- matrix(data = 0, nrow = 2, ncol = 15) # vector to hold the SIMULSOL results for dorsal and ventral side

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
        }else{  # doing ventral side. NB edit - adjust QSLR for PCOND here.
          FASKY <- 0.0
          FATOBJ <- 0.0
          FAVEG <- 0.0
          FAGRD <- FAGRDREF/(1 - FAGRDREF - FATOBJREF - FABUSHREF)
          FABUSH <- FABUSHREF/(1 - FAGRDREF - FATOBJREF - FABUSHREF)
          QSLR <- (QVENTR/(1 - FASKYREF - FATOBJREF -
                             FAVEGREF))*(1-(2*PCOND))
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

      # Calculating the "Cd" variable: Qcond = Cd(Tskin-Tsub), where Cd = Conduction area*((kfur/zfur)+(ksub/subdepth))
      if(S == 2){
        AREACND <- ATOT * (PCOND *2)
        CD <- AREACND * ((KFURCMPRS/ZFURCOMP))
        CONVAR<-CONVAR - AREACND #NB edit - Adjust area used for convection to account for PCOND. This is sent in to simulsol & then conv (unpacked as SURFAR)
      } else{ #doing dorsal side, no conduction. No need to adjust areas used for convection.
        AREACND = 0
        CD <- AREACND * ((KFURCMPRS/ZFURCOMP))
      }


      # package up inputs
      FURVARS <- c(LEN,ZFUR,FURTHRMK,KEFF,BETARA,FURTST,ZL)
      GEOMVARS <- c(SHAPE,SUBQFAT,CONVAR,VOL,D,CONVAR,CONVSK,RFUR,RFLESH,RSKIN,XR,RRAD,ASEMAJ,BSEMIN,CSEMIN,CD)
      ENVVARS <- c(FLTYPE,TA,TS,TBUSH,TVEG,TLOWER,TSKY,TCONDSB,RH,VEL,BP,ELEV,FASKY,FABUSH,FAVEG,FAGRD,QSLR)
      TRAITS <- c(TC,AK1,AK2,EMISAN,FATTHK,FLYHR,FURWET,PCTBAREVAP,PCTEYES)

      # set IPT, the geometry assumed in SIMULSOL: 1 = cylinder, 2 = sphere, 3 = ellipsoid
      if(SHAPE %in% c(1,3,5)){
        IPT <- 1
      }
      if(SHAPE == 2){
        IPT <- 2
      }
      if(SHAPE == 4){
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
    SHAPE_B_LAST <- SHAPE_B
    AK1LAST <- AK1
    TCLAST <- TC
    PANTLAST <- PANT
    SKINWLAST <- SKINW

    if(SHAPE_B < SHAPE_B_MAX){
      SHAPE_B <- SHAPE_B + UNCURL
    }else{
      SHAPE_B <- SHAPE_B_MAX
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

  # SIMULSOL output, dorsal
  TFA.D <- SIMULSOL.out[1, 1] # temperature of feathers/fur-air interface, deg C
  TSKCALCAV.D <- SIMULSOL.out[1, 2] # averagek skin temperature, deg C
  QCONV.D <- SIMULSOL.out[1, 3] # convection, W
  QCOND.D <- SIMULSOL.out[1, 4] # conduction, W
  QGENNET.D <- SIMULSOL.out[1, 5] # heat generation from flesh, W
  QSEVAP.D <- SIMULSOL.out[1, 6] # cutaneous evaporative heat loss, W
  QRAD.D <- SIMULSOL.out[1, 7] # radiation lost at fur/feathers/bare skin, W
  QSLR.D <- SIMULSOL.out[1, 8] # solar radiation, W
  QRSKY.D <- SIMULSOL.out[1, 9] # sky radiation, W
  QRBSH.D <- SIMULSOL.out[1, 10] # bush/object radiation, W
  QRVEG.D <- SIMULSOL.out[1, 11] # overhead vegetation radiation (shade), W
  QRGRD.D <- SIMULSOL.out[1, 12] # ground radiation, W
  QFSEVAP.D <- SIMULSOL.out[1, 13] # fur evaporative heat loss, W
  NTRY.D <- SIMULSOL.out[1, 14] # solution attempts, #
  SUCCESS.D <- SIMULSOL.out[1, 15] # successful solution found? (0 no, 1 yes)

  # SIMULSOL output, ventral
  TFA.V <- SIMULSOL.out[2, 1] # temperature of feathers/fur-air interface, deg C
  TSKCALCAV.V <- SIMULSOL.out[2, 2] # averagek skin temperature, deg C
  QCONV.V <- SIMULSOL.out[2, 3] # convection, W
  QCOND.V <- SIMULSOL.out[2, 4] # conduction, W
  QGENNET.V <- SIMULSOL.out[2, 5] # heat generation from flesh, W
  QSEVAP.V <- SIMULSOL.out[2, 6] # cutaneous evaporative heat loss, W
  QRAD.V <- SIMULSOL.out[2, 7] # radiation lost at fur/feathers/bare skin, W
  QSLR.V <- SIMULSOL.out[2, 8] # solar radiation, W
  QRSKY.V <- SIMULSOL.out[2, 9] # sky radiation, W
  QRBSH.V <- SIMULSOL.out[2, 10] # bush/object radiation, W
  QRVEG.V <- SIMULSOL.out[2, 11] # overhead vegetation radiation (shade), W
  QRGRD.V <- SIMULSOL.out[2, 12] # ground radiation, W
  QFSEVAP.V <- SIMULSOL.out[2, 13] # fur evaporative heat loss, W
  NTRY.V <- SIMULSOL.out[2, 14] # solution attempts, #
  SUCCESS.V <- SIMULSOL.out[2, 15] # successful solution found? (0 no, 1 yes)

  RESPFN <- ZBRENT.out[1] # heat sum (should be near zero), W
  QRESP <- ZBRENT.out[2] # respiratory heat loss, W
  GEVAP <- ZBRENT.out[3] # respiratory evaporation (g/s)
  PCTO2 <- ZBRENT.out[4] # O2 concentration (%)
  PCTN2 <- ZBRENT.out[5] # N2 concentration (%)
  PCTCO2 <- ZBRENT.out[6] # CO2 concentration (%)
  RESPGEN <- ZBRENT.out[7] # metabolic heat (W)
  O2STP <- ZBRENT.out[8] # O2 in rate at STP (L/s)
  O2MOL1 <- ZBRENT.out[9] # O2 in (mol/s)
  N2MOL1 <- ZBRENT.out[10] # N2 in (mol/s)
  AIRML1 <- ZBRENT.out[11] # air in (mol/s)
  O2MOL2 <- ZBRENT.out[12] # O2 out (mol/s)
  N2MOL2 <- ZBRENT.out[13] # N2 out (mol/s)
  AIRML2 <- ZBRENT.out[14] # air out (mol/s)
  AIRVOL <- ZBRENT.out[15] # air out at STP (L/s)

  HTOVPR <- 2.5012E+06 - 2.3787E+03 * TA # latent heat of vapourisation, W/kg/C
  SWEAT.G.S <- (QSEVAP.D + QSEVAP.V) * 0.5 / HTOVPR * 1000 # water lost from skin, g/s
  EVAP.G.S <- GEVAP + SWEAT.G.S # total evaporative water loss, g/s
  sigma <- 5.6697E-8
  QIROUT.D <- sigma * EMISAN * AREASKIN * (TSKCALCAV.D + 273.15) ^ 4
  QIRIN.D <- QRAD.D * -1 + QIROUT.D
  QIROUT.V <- sigma * EMISAN * AREASKIN * (TSKCALCAV.D + 273.15) ^ 4
  QIRIN.V <- QRAD.V * -1 + QIROUT.V

  QSOL <- QSLR.D * DMULT + QSLR.V * VMULT # solar, W
  QIRIN <- QIRIN.D * DMULT + QIRIN.V * VMULT # infrared in, W
  QMET <- RESPGEN # metabolism, W
  QEVAP <- QSEVAP.D * DMULT + QSEVAP.V * VMULT + QFSEVAP.D * DMULT + QFSEVAP.V * VMULT + QRESP # evaporation, W
  QIROUT <- QIROUT.D * DMULT + QIROUT.V * VMULT # infrared out, W
  QCONV <- QCONV.D * DMULT + QCONV.V * VMULT # convection, W
  QCOND <- QCOND.D * DMULT + QCOND.V * VMULT # conduction, W

  treg <- c(TC, TLUNG, TSKCALCAV.D, TSKCALCAV.V, TFA.D, TFA.V, SHAPE_B, PANT, SKINW, AK1, KEFARA[1], KEFARA[2], KEFARA[3], KFURCMPRS, Q10mult)
  names(treg) <- c("TC", "TLUNG", "TSKIN_D", "TSKIN_V", "TFA_D", "TFA_V", "SHAPE_B", "PANT", "SKINWET", "K_FLESH", "K_FUR", "K_FUR_D", "K_FUR_V", "K_COMPFUR", "Q10")

  morph <- c(ATOT, VOL, D, MASFAT, FATTHK, FLSHVL, ALENTH, AWIDTH, AHEIT, R1, R2, ASIL, ASILN, ASILP, AREASKIN, CONVSK, CONVAR, AREACND, FASKY, FAGRD)
  names(morph) <- c("AREA", "VOLUME", "CHAR_DIM", "MASS_FAT", "FAT_THICK", "FLESH_VOL", "LENGTH", "WIDTH", "HEIGHT", "DIAM_FLESH", "DIAM_FUR", "AREA_SIL", "AREA_SILN", "AREA_ASILP", "AREA_SKIN", "AREA_SKIN_EVAP", "AREA_CONV", "AREA_COND", "F_SKY", "F_GROUND")

  enbal <- c(QSOL, QIRIN, QMET, QEVAP, QIROUT, QCONV, QCOND, RESPFN, max(NTRY.D, NTRY.V), min(SUCCESS.D, SUCCESS.V))
  names(enbal) <- c("QSOL", "QIRIN", "QMET", "QEVAP", "QIROUT", "QCONV", "QCOND", "ENB", "NTRY", "SUCCESS")

  masbal <- c(AIRVOL, O2STP, GEVAP, SWEAT.G.S, O2MOL1, O2MOL2, N2MOL1, N2MOL2, AIRML1, AIRML2) * 3600
  names(masbal) <- c("AIR_L", "O2_L", "H2OResp_g", "H2OCut_g", "O2_mol_in", "O2_mol_out", "N2_mol_in", "N2_mol_out", "AIR_mol_in", "AIR_mol_out")

  return(list(treg = treg, morph = morph, enbal = enbal, masbal = masbal))
}
