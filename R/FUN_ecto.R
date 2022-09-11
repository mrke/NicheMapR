#' FUN_ecto
#'
#' R version of Fortran FUN.f (ectotherm model) for guessing a core body temperature that balances the heat budget.
#'
#' @encoding UTF-8
#' @param AMASS body mass, kg
#' @param GEOMETRY organism shape, 0-5, 0=plate, 1=cyl, 2=ellips, 3=lizard (desert iguana), 4=frog (leopard frog), 5=custom
#' @param ATOT total body surface area (m2)
#' @param AV ventral surface area (m2)
#' @param AT body surface area contacting another organism of same temperature (m2)
#' @param AL animal characteristic dimension (length) (m)
#' @param VOL body volume (m3)
#' @param R total body radius (m)
#' @param R1 flesh radius (i.e. radius to start of insulation if present) (m)
#' @param RINSUL depth of insulation (m) (not yet used)
#' @param ASEMAJR length of semi-major radius for ellipsoid (m)
#' @param BSEMINR length of semi-minor radius 1 for ellipsoid (m)
#' @param CSEMINR length of semi-minor radius 1 for ellipsoid (m)
#' @param M_1 metabolic rate parameter 1 V_O2=M_1*M^M_2*10^(M_3*Tb), in ml O2 / h, default parameters for lizards based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231
#' @param M_2 metabolic rate parameter 2
#' @param M_3 metabolic rate parameter 3
#' @param EXTREF oxygen extraction efficiency (\%)
#' @param PANT multiplier on breathing rate, for respiratory water loss via panting
#' @param RQ respiratory quotient (fractional, 0-1)
#' @param FLSHCOND thermal conductivity of flesh (W/mK)
#' @param PSI_BODY water potential of body (J/kg) - affects skin humidity for water vapour exchange
#' @param SKINW fraction of total surface area acting as a free water surface (fractional, 0-1)
#' @param AEFF effective area acting as a free-water exchanger, drives cutaneous evaporation (m2)
#' @param PEYES proportion of total surface area that is 'wet' eye (fractional, 0-1)
#' @param FATOSK configuration factor to sky (-)
#' @param FATOSB configuration factor to substrate (-)
#' @param FATOBJ configuration factor to nearby object (-) (not functional at the moment)
#' @param EMISAN emissivity of animal (fractional, 0-1)
#' @param EMISSB emissivity of substrate (fractional, 0-1)
#' @param EMISSK emissivity of sky (fractional, 0-1)
#' @param FLTYPE fluid type, air (0) or water (1)
#' @param TA air temperature (°C)
#' @param TSKY sky temperature (°C)
#' @param TSUBST substrate temperature (driving conduction) (°C)
#' @param TGRD ground temperature (driving radiation gain) (°C)
#' @param VEL wind speed (m/s)
#' @param QSOLAR solar radiation in, computed by SOLAR_ecto (W)
#' @param QIRIN longwave radiation in, computed by RADIN_ecto (W)
#' @param RELHUM relative humidity (\%)
#' @param BP air pressure (Pa)
#' @param ALT elevation (m)
#' @param SUBTK substrate thermal conductivity (W/mK)
#' @param O2GAS oxygen concentration of atmosphere (\%)
#' @param CO2GAS carbon dioxide concentration of atmosphere (\%)
#' @param X current guess of core body temperature (°C)
#' @export
FUN_ecto <- function(AMASS = AMASS,
                     GEOMETRY = GEOMETRY,
                     ATOT = AREA,
                     AV = AV,
                     AT = AT,
                     AL = AL,
                     VOL = VOL,
                     R = R,
                     R1 = R1,
                     RINSUL = RINSUL,
                     ASEMAJR = ASEMAJR,
                     BSEMINR = BSEMINR,
                     CSEMINR = CSEMINR,
                     M_1 = M_1,
                     M_2 = M_2,
                     M_3 = M_3,
                     EXTREF = EXTREF,
                     PANT = PANT,
                     RQ = RQ,
                     FLSHCOND = FLSHCOND,
                     PSI_BODY = PSI_BODY,
                     SKINW = SKINW,
                     AEFF = AEFF,
                     PEYES = PEYES,
                     FATOSK = FATOSK,
                     FATOSB = FATOSB,
                     FATOBJ = FATOBJ,
                     EMISAN = EMISAN,
                     EMISSB = EMISSB,
                     EMISSK = EMISSK,
                     FLTYPE = FLTYPE,
                     TA = TA,
                     TSKY = TSKY,
                     TSUBST = TSUBST,
                     TGRD = TGRD,
                     VEL = VEL,
                     QSOLAR = QSOLAR,
                     QIRIN = QIRIN,
                     RELHUM = RELHUM,
                     BP = BP,
                     ALT = ALT,
                     SUBTK = SUBTK,
                     O2GAS = O2GAS,
                     CO2GAS = CO2GAS,
                     N2GAS = N2GAS,
                     X){

  # C     NICHEMAPR: SOFTWARE FOR BIOPHYSICAL MECHANISTIC NICHE MODELLING
  #
  # C     COPYRIGHT (C) 2018 MICHAEL R. KEARNEY AND WARREN P. PORTER
  #
  # C     THIS PROGRAM IS FREE SOFTWARE: YOU CAN REDISTRIBUTE IT AND/OR MODIFY
  # C     IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
  # C     THE FREE SOFTWARE FOUNDATION, EITHER VERSION 3 OF THE LICENSE, OR (AT
  # C     YOUR OPTION) ANY LATER VERSION.
  #
  # C     THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT
  # C     WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
  # C     MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. SEE THE GNU
  # C     GENERAL PUBLIC LICENSE FOR MORE DETAILS.
  #
  # C     YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
  # C     ALONG WITH THIS PROGRAM. IF NOT, SEE HTTP://WWW.GNU.ORG/LICENSES/.
  #
  # C     EQUATIONS FOR STEADY STATE HEAT BUDGET, USED TO FIND TB VIA ROOT
  # C     FINDING ALGORITHM



  PI <- 3.14159265

  # C     THE GUESSED VARIABLE, X, IS CORE TEMPERATURE (C)
  # C     THIS ASSUMES UNIFORM BODY TEMPERATURE.

  #     CONTROL OF BODY TEMPERATURE GUESSES FOR STABILITY PURPOSES
  if(X > 100){
    X <- 100
  }

  TC <- X
  XTRY <- X

  #C     GET THE METABOLIC RATE
  #C     CHECKING FOR INANIMATE OBJECT
  #C      ALIVE, BUT IS IT TOO COLD?
  if(TC >= 0){
    MET.out <- MET_ecto(AMASS = AMASS,
                        XTRY = XTRY,
                        M_1 = M_1,
                        M_2 = M_2,
                        M_3 = M_3)
    QMETAB <- MET.out
  }else{
    #C       TOO COLD, SUPER LOW METABOLISM
    QMETAB <- 0.0001
    TC <- X
  }

  #C     GET THE RESPIRATORY WATER LOSS
  #C     CHECKING FOR FLUID TYPE
  if(FLTYPE == 0){
    #C      AIR
    #C      CALL FOR RESPIRATORY WATER & ENERGY LOSS
    if(QMETAB >= 0){
      RESP.out <- RESP_ecto(XTRY = XTRY,
                            AMASS = AMASS,
                            TC = TC,
                            QMETAB = QMETAB,
                            EXTREF = EXTREF,
                            PANT = PANT,
                            RQ = RQ,
                            TA = TA,
                            RELHUM = RELHUM,
                            BP = BP,
                            O2GAS = O2GAS,
                            CO2GAS = CO2GAS,
                            N2GAS = N2GAS)
      QRESP <- RESP.out$QRESP
    }else{
      #C       NEGATIVE METABOLIC RATE. NO PHYSIOLOGICAL MEANING - DEAD.
      QRESP <- 0
      QMETAB <- 0
    }
  }

  #C     NET INTERNAL HEAT GENERATION
  QGENET <- QMETAB - QRESP
  #C     NET INTERNAL HEAT GENERATION/UNIT VOLUME. USE FOR ESTIMATING SKIN TEMP.
  GN <- QGENET / VOL

  #C     COMPUTING SURFACE TEMPERATURE AS DICTATED BY GEOMETRY
  if(GEOMETRY == 0){
    #C      FLAT PLATE
    TSKIN <- TC - GN * R ^ 2 / (2 * FLSHCOND)
    RFLESH <- R
    TLUNG <- TC
  }

  #C     FIRST SET AVERAGE BODY TEMPERATURE FOR ESTIMATION OF AVEARAGE LUNG TEMPERATURE
  if(GEOMETRY == 1){
    #C      CYLINDER: FROM P. 270 BIRD, STEWART & LIGHTFOOT. 1960. TRANSPORT PHENOMENA.
    #C      TAVE = (GR ^ 2/(8K)) + TSKIN, WHERE TSKIN = TCORE - GR ^ 2/(4K)
    #C      NOTE:  THESE SHOULD ALL BE SOLVED SIMULTANEOUSLY.  THIS IS AN APPROXIMATION
    #C      USING CYLINDER GEOMETRY. SUBCUTANEOUS FAT IS ALLOWED IN CYLINDER & SPHERE
    #C      CALCULATIONS.
    RFLESH <- R1 - RINSUL
    TSKIN <- TC - GN * RFLESH ^ 2 / (4 * FLSHCOND)
    #C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN
    TLUNG <- (GN * RFLESH ^ 2) / (8 * FLSHCOND) + TSKIN
  }

  if(GEOMETRY == 2){
    #C      ELLIPSOID: DERIVED 24 OCTOBER, 1993  W. PORTER
    A <- ASEMAJR
    B <- BSEMINR
    C <- CSEMINR
    ASQ <- A ^ 2
    BSQ <- B ^ 2
    CSQ <- C ^ 2
    TSKIN <- TC - (GN / (2 * FLSHCOND)) * ((ASQ * BSQ * CSQ) / (ASQ * BSQ + ASQ * CSQ + BSQ * CSQ))
    #C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN
    TLUNG <- (GN / (4 * FLSHCOND)) * ((ASQ * BSQ * CSQ) / (ASQ * BSQ + ASQ * CSQ + BSQ * CSQ)) + TSKIN
  }

  if(GEOMETRY == 4){
    #C      SPHERE:
    RFLESH <- R1 - RINSUL
    RSKIN <- R1
    #C      FAT LAYER, IF ANY
    S1 <- (QGENET / (4 * PI * FLSHCOND)) * ((RFLESH - RSKIN) / (RFLESH * RSKIN))
    TSKIN <- TC - (GN * RFLESH ^ 2) / (6 * FLSHCOND) + S1
    #C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN (12 BECAUSE TLUNG IS 1/2 THE TC-TSKIN DIFFERENCE, 6*AK1)
    TLUNG <- (GN * RFLESH ^ 2) / (12 * FLSHCOND) + TSKIN
  }

  if((GEOMETRY == 3) | (GEOMETRY == 5)){
    # C      MODEL LIZARD/CUSTOM SHAPE AS CYLINDER
    # C      CYLINDER: FROM P. 270 BIRD, STEWART & LIGHTFOOT. 1960. TRANSPORT PHENOMENA.
    # C      TAVE = (GR ^ 2/(8K)) + TSKIN, WHERE TSKIN = TCORE - GR ^ 2/(4K)
    # C      NOTE:  THESE SHOULD ALL BE SOLVED SIMULTANEOUSLY.  THIS IS AN APPROXIMATION
    # C      USING CYLINDER GEOMETRY. SUBCUTANEOUS FAT IS ALLOWED IN CYLINDER & SPHERE
    # C      CALCULATIONS.
    RFLESH <- R1 - RINSUL
    TSKIN <- TC - GN * RFLESH ^ 2 / (4 * FLSHCOND)
    #C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN
    TLUNG <- (GN * RFLESH ^ 2) / (8 * FLSHCOND) + TSKIN
  }

  #C     LIMITING LUNG TEMPERATURE EXTREMES
  if(TLUNG > TC){
    TLUNG <- TC
  }

  CONV.out <- CONV_ecto(GEOMETRY = GEOMETRY,
                        ATOT = ATOT,
                        AV = AV,
                        AL = AL,
                        AT = AT,
                        BP = BP,
                        ALT = ALT,
                        TA = TA,
                        VEL = VEL,
                        FLTYPE = FLTYPE,
                        TSKIN = TSKIN)
  QCONV <- CONV.out$QCONV
  HD <- CONV.out$HD
  RESP.out <- RESP_ecto(XTRY = XTRY,
                        AMASS = AMASS,
                        TC = TC,
                        QMETAB = QMETAB,
                        EXTREF = EXTREF,
                        PANT = PANT,
                        RQ = RQ,
                        TA = TA,
                        RELHUM = RELHUM,
                        BP = BP,
                        O2GAS = O2GAS,
                        CO2GAS = CO2GAS,
                        N2GAS = N2GAS)
  QRESP <- RESP.out$QRESP
  GEVAP <- RESP.out$GEVAP
  SEVAP.out <- SEVAP_ecto(TC = TC,
                          TSKIN = TSKIN,
                          GEVAP = GEVAP,
                          PSI_BODY = PSI_BODY,
                          SKINW = SKINW,
                          AEFF = AEFF,
                          ATOT = ATOT,
                          HD = HD,
                          PEYES = PEYES,
                          TA = TA,
                          RELHUM = RELHUM,
                          VEL = VEL,
                          BP = BP)
  QSEVAP <- SEVAP.out$QSEVAP

  QIROUT <- RADOUT_ecto(TSKIN = TSKIN,
                        ATOT = ATOT,
                        AV = AV,
                        AT = AT,
                        FATOSK = FATOSK,
                        FATOSB = FATOSB,
                        EMISAN = EMISAN)

  QCOND <- COND_ecto(AV = AV,
                     TSKIN = TSKIN,
                     TSUBST = TSUBST,
                     SUBTK = SUBTK)

  if(FLTYPE == 1){
    #C      WATER ENVIRONMENT
    QSEVAP <- 0
    WEVAP <- 0
    QIRIN <- 0
    QIROUT <- 0
    QCOND <- 0
  }

  QIN <- QSOLAR + QIRIN + QMETAB
  QOUT <- QRESP + QSEVAP + QIROUT + QCONV + QCOND
  #C     FINDING THE DEVIATION FROM ZERO IN GUESSING THE SOLUTION
  ENB <- QIN - QOUT
  FUN <- ENB
  #message(paste(QSOLAR, QIRIN, QMETAB, QRESP, QSEVAP, QIROUT, QCONV, QCOND, TSKIN), '\n')

  return(FUN)
}
