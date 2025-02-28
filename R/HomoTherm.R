#' HomoTherm - the human model of NicheMapR
#'
#' This is a multi-part application of the endoR model to a human. It simulates a single
#' environmental scenario (use HomoTherm_var for a sequence of environments).
#' @encoding UTF-8
#' @param MASS = 70, mass of person (kg)
#' @param QMETAB_REST = 105, resting metabolic rate (W)
#' @param ACTIVE = FALSE, activity state (-)
#' @param MET = 1, MET units of activity (-)
#' @param INSDEPDs = c(1e-02, rep(6e-03, 3)), clothing depth, dorsal (m)
#' @param INSDEPVs = c(1e-09, rep(6e-03, 3)), clothing depth, ventral (m)
#' @param TA = 21, air temperature at local height (°C)
#' @param TGRD = TA, ground temperature (°C)
#' @param TSKY = TA, sky temperature (°C)
#' @param VEL = 0.1, wind speed (m/s)
#' @param RH = 50, relative humidity (\%)
#' @param QSOLR = 0, solar radiation, horizontal plane (W/m2)
#' @param Z = 20, zenith angle of sun (degrees from overhead)
#' @usage HomoTherm(MASS = 70, QMETAB_REST = 105, ACTIVE = FALSE, MET = 1, INSDEPDs = c(1e-02, rep(6e-03, 3)), INSDEPVs = c(1e-09, rep(6e-03, 3)), TA = 21, TGRD = TA, TSKY = TA, VEL = 0.1, RH = 50, QSOLR = 0, Z = 20,...)
#' @export
#' @details
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#' \code{EXCEED.TCMAX}{ = TRUE, allow the mode to continue increasing core temperature? (-)}\cr\cr
#' \code{MAXITER }{ = 500, maximum iterations beyond TC_MAX allowed when EXCEED.TMAX = TRUE}\cr\cr
#'
#' \strong{ Environment:}\cr\cr
#' \code{TAREF}{ = TA, air temperature at reference height (°C)}\cr\cr
#' \code{SHADE}{ = 0, shade on person (radiates at reference height temperature) (\%)}\cr\cr
#' \code{ELEV}{ = 0, elevation (m)}\cr\cr
#' \code{ABSSB}{ = 0.85, solar absorptivity of substrate (fractional, 0-1)}\cr\cr
#' \code{BP}{ = -1, Pa, negative means elevation is used}\cr\cr
#' \code{O2GAS}{ = 20.95, oxygen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)}\cr\cr
#' \code{N2GAS}{ = 79.02, nitrogen concetration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)}\cr\cr
#' \code{CO2GAS}{ = 0.0412, carbon dioxide concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)}\cr\cr
#' \code{PDIF}{ = 0.15, proportion of solar radiation that is diffuse (fractional, 0-1)}\cr\cr
#' \code{GRAV}{ = 9.80665, acceleration due to gravity, (m/s^2)}\cr\cr
#' \code{CONV_ENHANCE}{ = 1, convective enhancement factor (> 1 for turbulent outdoor conditions) (-)}\cr\cr
#'
#' \strong{ Whole body parameters:}\cr\cr
#' \code{MAXSWEAT}{ = 0.75, maximum sweating rate (L/h/m^2)}\cr\cr
#' \code{Q10}{ = 2, Q10 factor for adjusting BMR for TC}\cr\cr
#' \code{RQ}{ = 0.80, respiratory quotient (fractional, 0-1)}\cr\cr
#' \code{EXTREF}{ = 20, O2 extraction efficiency (\%)}\cr\cr
#'
#' \strong{ Part-specific morphological parameters (head, torso, arms, legs):}\cr\cr
#' \code{DENSITYs}{ = rep(1050, 4), body density (kg/m^3)}\cr\cr
#' \code{MASSFRACs}{ = c(0.0761, 0.501, 0.049, 0.162), fraction of total mass (-)}\cr\cr
#' \code{AREAFRACs}{ = c(0.0829, 0.327, 0.110, 0.185), fraction of total surface area (-)}\cr\cr
#' \code{PJOINs}{ = c(0.0275, 0.0824, 0.02174, 0.0333), fraction of part joined with rest of body (-)}\cr\cr
#' \code{SUBQFATs}{ = rep(1, 4), is subcutaneous fat present? (0 is no, 1 is yes)}\cr\cr
#' \code{FATPCT}{ = c(5, 36, 10, 23), \% body fat}\cr\cr
#' \code{SHAPE_Bs}{ = c(1.6, 1.9, 11, 7.0), ratio between long and short axis (-)}\cr\cr
#' \code{FSKREFs}{ = c(0.50, 0.42, 0.35, 0.35), configuration factor to sky}\cr\cr
#' \code{FGDREFs}{ = c(0.38, 0.42, 0.35, 0.35), reference configuration factor to ground}\cr\cr
#' \code{EMISANs}{ = rep(0.98, 4), emissivity each body part (-)}\cr\cr
#' \code{REFLD}{ = rep(0.3, 4), solar reflectivity dorsal (fractional, 0-1)}\cr\cr
#' \code{REFLV}{ = rep(0.3, 4), solar reflectivity ventral (fractional, 0-1)}\cr\cr
#'
#' \strong{ Part-specific physiological parameters (head, torso, arms, legs):}\cr\cr
#' \code{TC_RESTs}{ = rep(36.8, 4), resting core temperature (°C)}\cr\cr
#' \code{TC_ACTIVEs}{ = rep(37.5, 4), active core temperature (°C)}\cr\cr
#' \code{TC_INCs}{ = rep(0.04, 4), core temperature increment (°C)}\cr\cr
#' \code{TC_MAXs}{ = rep(38, 4), maximum tolerated core temperature (°C)}\cr\cr
#' \code{PCTWETs}{ = rep(4, 4), skin wettedness (\%)}\cr\cr
#' \code{PCTWET_INCs}{ = rep(0.5, 4), skin wettedness increment (\%)}\cr\cr
#' \code{PCTWET_MAXs}{ = rep(100, 4), maximum skin surface area that can be wet (\%)}\cr\cr
#' \code{CLOWETs}{ = rep(0, 4), insulation wettedness (\%)}\cr\cr
#' \code{PCTBAREVAPs}{ = c(60, 0, 0, 0), bare area where free and forced evaporation can occur (\%)}\cr\cr
#' \code{KFLESHs}{ = c(0.9, 0.9, 0.5, 0.5), flesh thermal conductivity (W/m°C)}\cr\cr
#' \code{KFLESH_INCs}{ = rep(0.05, 4), surface thermal conductivity increment (W/m°C)}\cr\cr
#' \code{KFLESH_MAXs}{ = rep(5, 4), maximum flesh conductivity (W/m°C)}\cr\cr
#' \code{KFATs}{ = rep(0.23, 4), fat conductivity (W/m°C)}\cr\cr
#'
#' \strong{ Insulation properties:}\cr\cr
#' \code{KCLOs}{ = rep(0, 4), insulation thermal conductivity manual override values (computed internally if zero) (W/mC)}\cr\cr
#' \code{DHAIRDs}{ = c(7.5e-5, rep(1E-06, 3)), fibre diameter, dorsal (m)}\cr\cr
#' \code{DHAIRVs}{ = c(7.5e-5, rep(1E-06, 3)), fibre diameter, ventral (m)}\cr\cr
#' \code{LHAIRDs}{ = c(50e-3, 50e-3, 50e-3, 50e-3), fibre length, dorsal (m)}\cr\cr
#' \code{LHAIRVs}{ = c(1e-9, 50e-3, 50e-3, 50e-3), fibre length, ventral (m)}\cr\cr
#' \code{INSDENDs}{ = rep(3e+08, 4), fibre density, dorsal (1/m2)}\cr\cr
#' \code{INSDENVs}{ = c(3e+05, rep(3e+08, 3)), fibre density, ventral (1/m2)}\cr\cr
#'
#' \strong{Outputs:}
#'
#' balance variables (general, whole-body output):
#' \itemize{
#' \item 1 T_CORE - core temperature (°C)
#' \item 2 T_LUNG - lung temperature (°C)
#' \item 3 T_SKIN - skin temperature (°C)
#' \item 4 T_CLO - insulation temperature (°C)
#' \item 5 PCTWET - skin wettedness (\%)
#' \item 6 K_FLESH - thermal conductivity of flesh (W/m°C)
#' \item 7 EVAP_CUT_L - cutaneous water loss (L/h)
#' \item 8 EVAP_RESP_L - respiratory water loss (L/h)
#' \item 9 SWEAT_L - water lost as sweat (may be higher than EVAP_CUT_L due to dripping) (L/h)
#' \item 10 K_FLESH - thermal conductivity of flesh (W/m°C)
#' \item 11 QMETAB  - metabolic heat production (W)
#' \item 12 QSLR - solar radiation absorbed (W)
#' \item 13 QIRIN - longwave (infra-red) radiation absorbed (W)
#' \item 14 QIROUT - longwave (infra-red) radiation lost (W)
#' \item 15 QCONV_RESP - respiratory sensible heat (W)
#' \item 16 QEVAP_RESP - respiratory evaporative heat (W)
#' \item 17 QEVAP_CUT - cutaneous evaporation (W)
#' \item 18 QCONV - convection (W)
#' \item 19 AREA - total surface area (m^2)
#' \item 20 AREA_RAD - total area for radiation exchange (m^2)
#' }
#' respire variables (respiratory response):
#' \itemize{
#' \item 1 AIR_L - air flowing through the lungs (L/h)
#' \item 2 O2_L - O2 consumed (L/h)
#' \item 3 O2_mol_in - inspired O2 (mol/h)
#' \item 4 O2_mol_out - expired O2 (mol/h)
#' \item 5 AIR_mol_in - inspired air (mol/h)
#' \item 6 AIR_mol_out - expired air (mol/h)
#' }
#' treg (thermoregulatory response variables, one table per body part):
#' \itemize{
#' \item 1 T_CORE - core temperature (°C)
#' \item 2 TSKIN_D - dorsal skin temperature (°C)
#' \item 3 TSKIN_V - ventral skin temperature (°C)
#' \item 4 TCLO_D - dorsal fur-air interface temperature (°C)
#' \item 5 TCLO_V - ventral fur-air interface temperature (°C)
#' \item 6 PCTWET - part of the skin surface that is wet (\%)
#' \item 7 K_FLESH - thermal conductivity of flesh (W/m°C)
#' \item 8 K_CLO_D - thermal conductivity of dorsal fur (W/m°C)
#' \item 9 Q10 - Q10 multiplier on metabolic rate (-)
#' }
#' morph variables (morphological traits, one table per body part):
#' \itemize{
#' \item 1 MASS - mass (kg)
#' \item 2 AREA - total outer surface area (m2)
#' \item 3 VOLUME - total volume (m3)
#' \item 4 CHAR_DIMENSION  - characteristic dimension for convection (m)
#' \item 5 MASS_FAT - fat mass (kg)
#' \item 6 FAT_THICK - thickness of fat layer (m)
#' \item 7 FLESH_VOL - flesh volume (m3)
#' \item 8 LENGTH - length (without fur) (m)
#' \item 9 WIDTH - width (without fur) (m)
#' \item 10 HEIGHT - height (without fur) (m)
#' \item 11 R_SKIN - radius, core to skin (m)
#' \item 12 R_FUR - radius, core to fur (m)
#' \item 13 AREA_SILHOUETTE - silhouette area (m2)
#' \item 14 AREA_SKIN - total skin area (m2)
#' \item 15 AREA_SKIN_EVAP - skin area available for evaporation (m2)
#' \item 16 AREA_CONV - area for convection (m2)
#' \item 17 AREA_JOIN - area for conduction (m2)
#' \item 18 F_SKY - configuration factor to sky (-)
#' \item 19 F_GROUND - configuration factor to ground (-)
#' }
#' enbal variables (energy balance, one table per body part):
#' \itemize{
#' \item 1 QSOL - solar radiation absorbed (W)
#' \item 2 QIRIN - longwave (infra-red) radiation absorbed (W)
#' \item 3 QGEN  - metabolic heat production (W)
#' \item 4 QEVAP - evaporation (W)
#' \item 5 QIROUT - longwave (infra-red) radiation lost (W)
#' \item 6 QCONV - convection (W)
#' \item 7 ENB - energy balance (W)
#' \item 8 NTRY - iterations required for a solution (-)
#' \item 9 SUCCESS - was a solution found (0=no, 1=yes)
#' }
#' @examples
#' library(NicheMapR)
#' # environment
#' TA <- 0 # air temperature, °C
#' VEL <- 0.1 # wind speeds, m/s
#' RH <- 50 # humidity, %
#' # set insulation depth, flesh conductivity and fat
#' INSDEPDs <- c(1e-02, rep(6.15e-03, 3)) # 'dorsal' clothing depth, m
#' INSDEPVs <- c(1e-09, rep(6.15e-03, 3)) # 'ventral' clothing depth, m
#' KCLOs <- rep(0.04, 4) # clothing thermal conductivity, W/m·K
#' FATPCTs <- c(5, 36, 10, 23) # body fat %
# simulate
#' HomoTherm.out <- HomoTherm(INSDEPDs = INSDEPDs * 0,
#'                                INSDEPVs = INSDEPVs * 0,
#'                                KCLOs = KCLOs,
#'                                FATPCTs = FATPCTs,
#'                                TA = TA,
#'                                VEL = VEL,
#'                                RH = RH,
#'                                EXCEED.TCMAX = TRUE)
#'balance <- HomoTherm.out$balance
#'balance # report output
HomoTherm <- function(MASS = 70,
                      QMETAB_REST = 105,
                      ACTIVE = FALSE,
                      MET = 1,
                      Q10 = 2,
                      RQ = 0.8,
                      EXTREF = 25,
                      TC_RESTs = rep(36.8, 4),
                      TC_ACTIVEs = rep(37.5, 4),
                      TC_MAXs = rep(38, 4),
                      EXCEED.TCMAX = TRUE,
                      TC_INCs = rep(0.04, 4),
                      DENSITYs = rep(1050, 4),
                      MASSFRACs = c(0.07609801, 0.50069348, 0.04932963, 0.16227462),
                      AREAFRACs = c(0.08291887, 0.32698460, 0.11025155, 0.18479669),
                      SHAPE_Bs = c(1.6, 1.9, 11, 7.0),
                      PJOINs = c(0.02753623, 0.08239728, 0.02173913, 0.03333333),
                      SUBQFATs =  rep(1, 4),
                      FATPCTs =  c(5, 36, 10, 23),
                      KFLESHs = c(0.9, 0.9, 0.5, 0.5),
                      KFLESH_MAXs = rep(5, 4),
                      KFLESH_INCs = rep(0.05, 4), #rep(0.2, 4),
                      KFATs = rep(0.23, 4),
                      KCLOs = rep(0, 4),
                      DHAIRDs = c(7.5e-5, rep(1E-06, 3)),
                      DHAIRVs = c(7.5e-5, rep(1E-06, 3)),
                      LHAIRDs = c(50e-3, 50e-3, 50e-3, 50e-3),
                      LHAIRVs = c(1e-9, 50e-3, 50e-3, 50e-3),
                      INSDEPDs = c(1e-02, rep(6e-03, 3)),
                      INSDEPVs = c(1e-09, rep(6e-03, 3)),
                      INSDENDs = rep(3e+08, 4),
                      INSDENVs = c(3e+05, rep(3e+08, 3)),
                      REFLDs = rep(0.3, 4),
                      REFLVs = rep(0.3, 4),
                      EMISANs = rep(0.98, 4),
                      FSKREFs = c(0.50, 0.42, 0.35, 0.35),
                      FGDREFs = c(0.38, 0.42, 0.35, 0.35),
                      PCTWETs = rep(4, 4),
                      PCTWET_INCs = rep(0.5, 4),
                      PCTWET_MAXs = rep(100, 4),
                      PCTBAREVAPs = c(60, 0, 0, 0),
                      MAXSWEAT = 0.75,
                      CLOWETs = rep(0, 4),
                      CONV_ENHANCE = 1,
                      TA = 21,
                      TAREF = TA[1],
                      TSKY = TA[1],
                      TGRD = TA[1],
                      VEL = 1,
                      RH = 50,
                      QSOLR = 0,
                      SHADE = 0,
                      PDIF = 0.15,
                      Z = 20,
                      ELEV = 0,
                      BP = 101325,
                      GRAV = 9.80665,
                      ABSSB = 0.85,
                      O2GAS = 20.95,
                      N2GAS = 79.02,
                      CO2GAS = 0.0422,
                      MAXITER = 500){
  # check if only one air temp / wind speed / relative humidity specified
  # needs one value for each body part so gradients with height can be
  # accounted for
  SHAPEs <- c(4, 1, 1, 1) # ellipsoid (head) and cylinders (trunk, arms, legs)
  NPARTs <- c(1, 1, 2, 2) # one head, one trunk two arms and two legs
  PVENs <-  rep(0.5, 4) # make fraction ventral even
  ORIENTs = rep(0, 4) # don't orient perp or parallel to sun
  AK1_SUBQFAT <- 0.51 # flesh conductivity value at which subcutaneous fat influence declines
  SUBQFAT_REDUCE <- 0.9 # factor by which subcutaneous fat influence is reduced
  if(length(TA == 1)){
    TA <- rep(TA, 4)
  }
  if(length(TA == 1)){
    VEL <- rep(VEL, 4)
  }
  if(length(TA == 1)){
    RH <- rep(RH, 4)
  }
  get.vec <- function(sublist, subitem, x){
    unlist(lapply(lapply(x, `[[`, sublist), `[[`, subitem))
  }
  TC_REF <- mean(TC_RESTs)
  if(ACTIVE){
    QMETAB_MIN <- QMETAB_REST * MET
    TCs <- TC_ACTIVEs
  }else{
    QMETAB_MIN <- QMETAB_REST
    TCs <- TC_RESTs
  }
  QMET_REF <- QMETAB_MIN
  QGEN <- 0
  parts2do <- length(SHAPEs)
  parts <- vector(mode = "list", length = parts2do)
  counter <- 0
  counter2 <- 0
  SWEAT <- 0
  HEIGHT_GUESS <- 100 * (MASS / 23.5) ^ 0.5 # BMI equation, assuming 23.5 BMI (this gets update after first iteration)
  AREA <- 0.00718 * MASS ^ 0.425 * HEIGHT_GUESS ^ 0.725 # DuBois area, m2

  while(QGEN < QMETAB_MIN){

    for(i in 1:parts2do){  # cycle through each body part
      ANDENS <- DENSITYs[i]
      MASSFRAC <- MASSFRACs[i]
      SHAPE <- SHAPEs[i]
      SHAPE_B <- SHAPE_Bs[i]
      SHAPE_B_MAX <- SHAPE_B
      PVEN <- PVENs[i]
      PCOND <- PJOINs[i]
      SUBQFAT <- SUBQFATs[i]
      FATPCT <- FATPCTs[i]
      TC <- TCs[i]
      TC_INC <- TC_INCs[i]
      TC_MAX <- TC_MAXs[i]
      AK1 <- KFLESHs[i]
      #AK1_MAX <- KFLESH_MAXs[i]
      #AK1_INC <- KFLESH_INCs[i]
      AK2 <- KFATs[i]
      FURTHRMK <- KCLOs[i]
      DHAIRD <- DHAIRDs[i]
      DHAIRV <- DHAIRVs[i]
      LHAIRD <- LHAIRDs[i]
      LHAIRV <- LHAIRVs[i]
      ZFURD <- INSDEPDs[i]
      ZFURV <- INSDEPVs[i]
      RHOD <- INSDENDs[i]
      RHOV <- INSDENVs[i]
      REFLD <- REFLDs[i]
      REFLV <- REFLVs[i]
      EMISAN <- EMISANs[i]
      FGDREF <- FGDREFs[i]
      FSKREF <- FSKREFs[i]
      ORIENT <- ORIENTs[i]
      PCTWET <- PCTWETs[i]
      PCTWET_INC <- PCTWET_INCs[i]
      PCTWET_MAX <- PCTWET_MAXs[i]
      FURWET <- CLOWETs[i]
      PCTBAREVAP <- PCTBAREVAPs[i]
      AMASS <- MASSFRAC * MASS
      QBASPART <- QMETAB_MIN * MASSFRAC
      TREGMODE <- 2 # irrelevant because THERMOREG = 0
      THERMOREG <- 0 # turn off thermoregulation for body part
      TCONDSB <- TC # conduction is to other limbs via joins
      parts[[i]] <- endoR(ANDENS = ANDENS,
                          FATDEN = ANDENS, # doing this to prevent changes in area as fat is cut down during thermoreg response to simulate vasodilation
                          EMISAN = EMISAN,
                          THERMOREG = THERMOREG,
                          TREGMODE = TREGMODE,
                          RESPIRE = 0,
                          PZFUR = 0,
                          TC_MIN = 35,
                          ZFURD_MAX = ZFURD,
                          ZFURV_MAX = ZFURV,
                          PCOND = PCOND,
                          FATPCT = FATPCT,
                          SUBQFAT = SUBQFAT,
                          AK1 = AK1,
                          PCTBAREVAP = PCTBAREVAP,
                          TA = TA[i], # may vary per body part
                          TAREF = TAREF,
                          TSKY = TSKY,
                          TGRD = TGRD,
                          VEL = VEL[i], # may vary per body part
                          RH = RH[i], # may vary per body part
                          QSOLR = QSOLR,
                          PDIF = PDIF,
                          Z = Z,
                          ELEV = ELEV,
                          ABSSB = ABSSB,
                          AMASS = AMASS,
                          SHAPE = SHAPE,
                          PVEN = PVEN,
                          SHAPE_B = SHAPE_B,
                          SHAPE_B_MAX = SHAPE_B_MAX,
                          DHAIRD = DHAIRD,
                          DHAIRV = DHAIRV,
                          LHAIRD = LHAIRD,
                          LHAIRV = LHAIRV,
                          ZFURD = ZFURD,
                          ZFURV = ZFURV,
                          RHOD = RHOD,
                          RHOV = RHOV,
                          REFLD = REFLD,
                          REFLV = REFLV,
                          PCTWET = PCTWET,
                          PCTWET_INC = PCTWET_INC,
                          FGDREF = FGDREF,
                          FSKREF = FSKREF,
                          TC = TC,
                          TCONDSB = TCONDSB,
                          PANT_INC = 0.0001,
                          PANT_MAX = 1.0001,
                          TC_INC = TC_INC,
                          TC_MAX = TC_MAX,
                          BP = BP,
                          GRAV = GRAV,
                          QBASAL = QBASPART,
                          FURTHRMK = FURTHRMK,
                          SHADE = SHADE,
                          FURWET = FURWET,
                          ZFURCOMP = 0.5,
                          ORIENT = ORIENT,
                          CONV_ENHANCE = CONV_ENHANCE,
                          Q10 = Q10)
    }

    # if first iteration, get surface area
    if(counter2 == 0){
      mass.parts <- MASSFRACs * MASS
      head.morph <- c(mass.parts[1], parts[[1]][2]$morph[c(1:12, 15:20)])
      trunk.morph <- c(mass.parts[2], parts[[2]][2]$morph[c(1:12, 15:20)])
      arm.morph <- c(mass.parts[3], parts[[3]][2]$morph[c(1:12, 15:20)])
      leg.morph <- c(mass.parts[4], parts[[4]][2]$morph[c(1:12, 15:20)])
      names(head.morph) <- c("MASS", "AREA", "VOLUME", "CHAR_DIMENSION", "MASS_FAT", "FAT_THICK", "FLESH_VOL", "LENGTH", "WIDTH", "HEIGHT", "R_SKIN", "R_INS", "AREA_SILHOUETTE", "AREA_SKIN", "AREA_SKIN_EVAP", "AREA_CONV", "AREA_JOIN", "F_SKY", "F_GROUND")
      names(trunk.morph) <- names(head.morph)
      names(arm.morph) <- names(head.morph)
      names(leg.morph) <- names(head.morph)
      AREA <- head.morph[2] - head.morph[17] + trunk.morph[2] - trunk.morph[17] + (arm.morph[2] - arm.morph[17]) * 2 + (leg.morph[2] - leg.morph[17]) * 2
      AREA_RAD <- (head.morph[2] - head.morph[17]) * (FSKREFs[1] + FGDREFs[1]) + (trunk.morph[2] - trunk.morph[17]) * (FSKREFs[2] + FGDREFs[2]) + (arm.morph[2] - arm.morph[17]) * (FSKREFs[3] + FGDREFs[3]) * 2 + (leg.morph[2] - leg.morph[17]) * (FSKREFs[4] + FGDREFs[4]) * 2
    }

    # weighted (by surface area) mean skin and fur temperature
    # get.vec(a, b, c) gets element b of sub-list a from list c
    TSKIN <- sum((get.vec(1, 3, parts) + get.vec(1, 4, parts)) / 2 * AREAFRACs  * NPARTs)
    TC <- sum(unlist(lapply(lapply(parts, `[[`, 1), `[[`, 1)) * AREAFRACs  * NPARTs)
    TSKINDs <- c(parts[[1]]$treg[3], parts[[2]]$treg[3], parts[[3]]$treg[3], parts[[4]]$treg[3])
    TSKINVs <- c(parts[[1]]$treg[4], parts[[2]]$treg[4], parts[[3]]$treg[4], parts[[4]]$treg[4])
    QGEN <- sum(get.vec(3, 3, parts) * NPARTs)
    TCORE <- parts[[2]]$treg[1]

    # compute lung temp as temperature half way between core to fat
    trunk.morph <- c(0, parts[[2]][2]$morph[c(1:12, 15:20)])
    trunk.enbal <- parts[[2]][3]$enbal[c(1:7, 9:10)][c(3, 1, 2, 5, 4, 6:9)]
    Q_gen <- trunk.enbal[1]
    V_gen <- trunk.morph[7]
    R_skin <- trunk.morph[11]
    R_gen <- R_skin - trunk.morph[6]
    R_lung <- R_gen / 2
    T_fat <- (TSKINDs[2] + TSKINVs[2]) / 2 + ((Q_gen * R_gen ^ 2 ) / (2 * KFATs[2] * V_gen)) * log((R_skin / (R_gen)))
    #TLUNG <- (TC + TSKIN) / 2 # average of skin and core, but no greater than core, °C
    #TLUNG <- (TCORE + T_fat) / 2 # average of fat and core, °C
    TLUNG <- T_fat + (Q_gen * (R_gen ^ 2 - R_lung ^ 2)) / (4 * KFLESHs[2] * V_gen)
    #TAEXIT <- min(mean(TA[1], TLUNG), TLUNG) # temperature of exhaled air, no greater than lung, °C
    #TAEXIT <- min(0.25 * TA[1] + 26, TLUNG) # from Ferrus et al. 1980. Respiratory water loss. Respiration Physiology 39:367–381.
    TAEXIT <- min(0.0837 * TA[1] + 32, TLUNG) # from Ferrus et al. 1980. Respiratory water loss. Respiration Physiology 39:367–381.
    TFA <- sum((unlist(lapply(lapply(parts, `[[`, 1), `[[`, 5)) + unlist(lapply(lapply(parts, `[[`, 1), `[[`, 6))) / 2 * AREAFRACs  * NPARTs)
    PCTWET <- min(100, sum(unlist(lapply(lapply(parts, `[[`, 1), `[[`, 9)) * AREAFRACs  * NPARTs))
    AK1 <- sum(unlist(lapply(lapply(parts, `[[`, 1), `[[`, 10)) * AREAFRACs  * NPARTs)

    if(mean(TA[1:2]) < TC){
      QM1 <- QMETAB_MIN * 2 * -1
      QM2 <- QMETAB_MIN * 50
    }else{
      QM1 <- QMETAB_MIN * 10 * -1
      QM2 <- QMETAB_MIN * 2
    }
    QSUM <- QGEN
    TOL <- MASS * 0.001
    PANT <- 1
    TIMACT <- 1
    RELXIT <- 92 # assumed, based on Berry, E. (1914). Relative Humidity of Expired Air. American Physical Education Review, 19(6), 452–454. https://doi.org/10.1080/23267224.1914.10651422
    ZBRENT.in <- c(TA[1], O2GAS, N2GAS, CO2GAS, BP, QMETAB_MIN, RQ, TLUNG, GMASS = MASS * 1000, EXTREF, RH[1],
                   RELXIT, TIMACT, TAEXIT, QSUM, PANT)

    # call ZBRENT subroutine which calls RESPFUN
    ZBRENT.out <- ZBRENT_ENDO(QM1, QM2, TOL, ZBRENT.in)
    colnames(ZBRENT.out) <- c("RESPFN", "QRESP", "GEVAP", "PCTO2", "PCTN2", "PCTCO2", "RESPGEN", "O2STP", "O2MOL1", "N2MOL1", "AIRML1", "O2MOL2", "N2MOL2", "AIRML2", "AIRVOL")

    #RESPFN <- ZBRENT.out[1] # heat sum (should be near zero), W
    QRESP <- ZBRENT.out[2] # respiratory heat loss, W
    GEVAP <- ZBRENT.out[3] # respiratory evaporation (g/s)
    #PCTO2 <- ZBRENT.out[4] # O2 concentration (%)
    #PCTN2 <- ZBRENT.out[5] # N2 concentration (%)
    #PCTCO2 <- ZBRENT.out[6] # CO2 concentration (%)
    QGEN <- ZBRENT.out[7] # metabolic heat (W)
    O2STP <- ZBRENT.out[8] # O2 in rate at STP (L/s)
    O2MOL1 <- ZBRENT.out[9] # O2 in (mol/s)
    N2MOL1 <- ZBRENT.out[10] # N2 in (mol/s)
    AIRML1 <- ZBRENT.out[11] # air in (mol/s)
    O2MOL2 <- ZBRENT.out[12] # O2 out (mol/s)
    N2MOL2 <- ZBRENT.out[13] # N2 out (mol/s)
    AIRML2 <- ZBRENT.out[14] # air out (mol/s)
    AIRVOL <- ZBRENT.out[15] # air out at STP (L/s)

    # compute required evaporation rate and sweating efficiency, and thus required sweat rate
    HTOVPR <- 2.5012E+06 - 2.3787E+03 * mean(TA)
    EVAP_RESP_L <- GEVAP * 3600 / 1000 # L/h
    QEVAP <- sum(unlist(lapply(lapply(parts, `[[`, 3), `[[`, 4)) * NPARTs)
    EVAP_CUT_L <- QEVAP / HTOVPR * 3600 # L/h
    r <- 1 - (PCTWET / 100)^2 / 2 # Parsons p. 43, sweating efficiency
    SWEAT <- EVAP_CUT_L / r # L/h

    # cap sweating at maximum rate
    if(SWEAT > MAXSWEAT * AREA){
      PCTWET_MAXs <- PCTWETs
    }
    if(QGEN < QMETAB_MIN){
      for(i in 1:4){
      if(KFLESHs[i] < KFLESH_MAXs[i]){ # vasodilation
        KFLESHs[i] <- KFLESHs[i] + KFLESH_INCs[i]
        if(KFLESHs[i] > AK1_SUBQFAT){
          FATPCTs[i] <- FATPCTs[i] * SUBQFAT_REDUCE # blood flow bypasses fat layer
        }
      }else{
        if(PCTWETs[i] < PCTWET_MAXs[i]){ # sweating
          PCTWETs[i] <- PCTWETs[i] + PCTWET_INCs[i]
        }else{
          TCs[i] <- TCs[i] + TC_INCs[i] # let core rise (adjust metabolic rate accordingly)
          Q10mult <- Q10 ^ ((mean(TCs) - TC_REF) / 10)
          QMETAB_MIN <- QMET_REF * Q10mult
        }
      }
      }
    }
    if(max(TSKINDs, TSKINVs) > 33){ # augment current vasoconstriction with sweating and core temperature rise
      PCTWETs <- PCTWETs + PCTWET_INC / 6
      for(i in 1:length(TCs)){
        if(TCs[i] < TC_MAX){
          TCs[i] <- TCs[i] + TC_INC / 1.5
          Q10mult <- Q10 ^ ((mean(TCs) - TC_REF) / 10)
          QMETAB_MIN <- QMET_REF * Q10mult
        }
      }
    }
    if(max(TSKINDs, TSKINVs) > TC - 2){ # augment sweating rate further
      PCTWETs <- PCTWETs + PCTWET_INC * 2
    }
    for(i in 1:4){
      if(PCTWETs[i] > PCTWET_MAXs[i]){
        PCTWETs[i] <- PCTWET_MAXs[i]
      }
      if(KFLESHs[i] > KFLESH_MAXs[i]){
        KFLESHs[i] <- KFLESH_MAXs[i]
      }
    }
    if(max(TCs) >= TC_MAX + TC_INC | TC_INC == 0){
      counter <- counter + 1
      if(!EXCEED.TCMAX | counter > MAXITER){
       break # dead
      }
    }
    counter2 <- counter2 + 1
  }

  ZBRENT.out <- as.data.frame(ZBRENT.out)
  respire1 <- c(AIRVOL, O2STP, O2MOL1, O2MOL2, AIRML1, AIRML2) * 3600
  respire <- matrix(data = respire1, nrow = 1, ncol = 6)
  respire.names <- c("AIR_L", "O2_L", "O2_mol_in", "O2_mol_out", "AIR_mol_in", "AIR_mol_out")
  colnames(respire) <- respire.names

  # sum heat exchange terms
  QSLR <- sum(unlist(lapply(lapply(parts, `[[`, 3), `[[`, 1)) * NPARTs)
  QRAD_IN <- sum(unlist(lapply(lapply(parts, `[[`, 3), `[[`, 2)) * NPARTs)
  QEVAP <- sum(unlist(lapply(lapply(parts, `[[`, 3), `[[`, 4)) * NPARTs)
  QRAD_OUT <- sum(unlist(lapply(lapply(parts, `[[`, 3), `[[`, 5)) * NPARTs)
  QCONV <- sum(unlist(lapply(lapply(parts, `[[`, 3), `[[`, 6)) * NPARTs)
  QCOND <- sum(unlist(lapply(lapply(parts, `[[`, 3), `[[`, 7)) * NPARTs)

  mass.parts <- MASSFRACs * MASS
  head.treg <- parts[[1]][1]$treg[c(1, 3:6, 9, 10, 12:13, 17)]
  head.morph <- c(mass.parts[1], parts[[1]][2]$morph[c(1:12, 15:20)])
  head.enbal <- parts[[1]][3]$enbal[c(1:7, 9:10)][c(3, 1, 2, 5, 4, 6:9)]
  trunk.treg <- parts[[2]][1]$treg[c(1, 3:6, 9, 10, 12:13, 17)]
  trunk.morph <- c(mass.parts[2], parts[[2]][2]$morph[c(1:12, 15:20)])
  trunk.enbal <- parts[[2]][3]$enbal[c(1:7, 9:10)][c(3, 1, 2, 5, 4, 6:9)]
  arm.treg <- parts[[3]][1]$treg[c(1, 3:6, 9, 10, 12:13, 17)]
  arm.morph <- c(mass.parts[3], parts[[3]][2]$morph[c(1:12, 15:20)])
  arm.enbal <- parts[[3]][3]$enbal[c(1:7, 9:10)][c(3, 1, 2, 5, 4, 6:9)]
  leg.treg <- parts[[4]][1]$treg[c(1, 3:6, 9, 10, 12:13, 17)]
  leg.morph <- c(mass.parts[4], parts[[4]][2]$morph[c(1:12, 15:20)])
  leg.enbal <- parts[[4]][3]$enbal[c(1:7, 9:10)][c(3, 1, 2, 5, 4, 6:9)]
  names(head.treg) <- c("T_CORE","TSKIN_D", "TSKIN_V", "TCLO_D", "TCLO_V", "PCTWET", "K_FLESH", "K_CLO_D", "K_CLO_V", "Q10")
  names(trunk.treg) <- names(head.treg)
  names(arm.treg) <- names(head.treg)
  names(leg.treg) <- names(head.treg)
  #names(head.enbal) <- c("QSOL", "QIRIN", "QGEN", "QEVAP", "QIROUT", "QCONV", "ENB", "NTRY", "SUCCESS")
  names(head.enbal) <- c("QMETAB", "QSLR", "QRAD_IN", "QRAD_OUT", "QEVAP", "QCONV", "ENB", "NTRY", "SUCCESS")
  names(trunk.enbal) <- names(head.enbal)
  names(arm.enbal) <- names(head.enbal)
  names(leg.enbal) <- names(head.enbal)
  names(head.morph) <- c("MASS", "AREA", "VOLUME", "CHAR_DIMENSION", "MASS_FAT", "FAT_THICK", "FLESH_VOL", "LENGTH", "WIDTH", "HEIGHT", "R_SKIN", "R_INS", "AREA_SILHOUETTE", "AREA_SKIN", "AREA_SKIN_EVAP", "AREA_CONV", "AREA_JOIN", "F_SKY", "F_GROUND")
  names(trunk.morph) <- names(head.morph)
  names(arm.morph) <- names(head.morph)
  names(leg.morph) <- names(head.morph)
  AREA <- head.morph[2] - head.morph[17] + trunk.morph[2] - trunk.morph[17] + (arm.morph[2] - arm.morph[17]) * 2 + (leg.morph[2] - leg.morph[17]) * 2
  AREA_RAD <- (head.morph[2] - head.morph[17]) * (FSKREFs[1] + FGDREFs[1]) + (trunk.morph[2] - trunk.morph[17]) * (FSKREFs[2] + FGDREFs[2]) + (arm.morph[2] - arm.morph[17]) * (FSKREFs[3] + FGDREFs[3]) * 2 + (leg.morph[2] - leg.morph[17]) * (FSKREFs[4] + FGDREFs[4]) * 2
  CP <- WETAIR(db = TAEXIT, rh = RH[1], bp = BP)$cp # heat capacity of air (J/KG/K)
  QCONV_RESP <- CP * AIRML1 * 0.0289647 * (TA[1] - TLUNG) # HEAT CAPCITY (J/KG/K) * MOLES AIR (MOL / S) * MOLAR MASS OF AIR (KG/MOL) * DELTA TEMPERATURE
  QEVAP_RESP <- QRESP + QCONV_RESP
  balance <- c(TCORE, TLUNG, TSKIN, TFA, PCTWET, AK1, EVAP_CUT_L, EVAP_RESP_L, SWEAT, QGEN, QSLR, QRAD_IN, QRAD_OUT, QCONV_RESP, QEVAP_RESP * -1, QEVAP * -1, QCONV * -1, AREA, AREA_RAD)
  names(balance) <- c("T_CORE", "T_LUNG", "T_SKIN", "T_CLO", "PCTWET", "K_FLESH", "EVAP_CUT_L", "EVAP_RESP_L", "SWEAT_L", "QMETAB", "QSLR", "QRAD_IN", "QRAD_OUT", "QCONV_RESP", "QEVAP_RESP", "QEVAP_CUT", "QCONV", "AREA", "AREA_RAD")
  return(list(balance = balance,
              respire = respire,
              head.treg = head.treg,
              head.morph = head.morph,
              head.enbal = head.enbal,
              trunk.treg = trunk.treg,
              trunk.morph = trunk.morph,
              trunk.enbal = trunk.enbal,
              arm.treg = arm.treg,
              arm.morph = arm.morph,
              arm.enbal = arm.enbal,
              leg.treg = leg.treg,
              leg.morph = leg.morph,
              leg.enbal = leg.enbal))
}
