#' Ectotherm heat budget model
#'
#' A modularised implementation of the heat budget of the Niche Mapper ectotherm model that computes body temperature,
#' and water loss for a single environment with no behaviour. Like endoR_devel, it can be used to develop customised
#' versions of the NicheMapR ectotherm model.
#'
#' @encoding UTF-8
#' @param Ww_g = 40, Wet weight of animal (g), note this model is 'steady state' so no lags in heating/cooling due to mass
#' @param shape = 3, Organism shape, 0-5, Determines whether standard or custom shapes/surface area/volume relationships are used: 0=plate, 1=cyl, 2=ellips, 3=lizard (desert iguana), 4=frog (leopard frog), 5=custom (see details)
#' @param alpha = 0.85, Solar absorptivity, 0-1
#' @param M_1 = 0.013, Metabolic rate parameter 1 V_O2=M_1*M^M_2*10^(M_3*Tb), in ml O2 / h, default parameters for lizards based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231
#' @param M_2 = 0.800, Metabolic rate parameter 2
#' @param M_3 = 0.038, Metabolic rate parameter 3
#' @param pct_wet = 0.1, \% of surface area acting as a free-water exchanger, for computing cutaneous water loss
#' @param pct_eyes = 0, \% of surface area taken up by open eyes, for computing ocular water loss (only when active)
#' @param pct_mouth = 0, \% of surface area taken up by open mouth, for computing panting water loss
#' @param pantmax = 1, maximum multiplier on breathing rate, for respiratory water loss via panting (>1 invokes panting)
#' @param F_O2 = 20, \% oxygen extraction efficiency, for respiratory water loss
#' @param O2gas = 20.95, \% O2 in air
#' @param CO2gas = 0.03, \% CO2 in air
#' @param N2gas = 79.02, \% nitrogen in air
#' @param psi_body = -7.07 * 100, water potential of body (J/kg) - affects skin humidity for water vapour exchange
#' @param delta_air = 0.1, temperature difference (°C) between expired and inspired air, for computing respiratory water loss
#' @param TA = 20, air temperature at local height (°C)
#' @param TGRD = 40, ground temperature (°C)
#' @param TSKY = -5, sky temperature (°C)
#' @param VEL = 1, wind speed (m/s)
#' @param RH = 30, relative humidity (\%)
#' @param QSOLR = 800, solar radiation, horizontal plane (W/m2)
#' @param Z = 20, zenith angle of sun (degrees from overhead)
#' @param SHADE = 0, shade level (\%)
#' @usage ectotherm(Ww_g = 40, shape = 3, alpha = 0.85, postur = 0, TA = 20, TGRD = 40, TSKY = -5, VEL = 1, RH = 30, QSOLR = 800, Z = 20 ...)
#' @details
#' \strong{ Environmental inputs:}
#'
#' \itemize{
#' \item{\code{elev}{ = 0, elevation (m)}\cr}
#' \item{\code{alpha_sub}{ = 0.8, solar absorptivity of substrate (fractional, 0-1)}\cr}
#' \item{\code{fluid}{ = 0, fluid type: 0 = air; 1 = water}\cr}
#' \item{\code{TSUBST}{ = TGRD, surface temperature for conduction (°C)}\cr}
#' \item{\code{K_sub}{ = 2.79, substrate thermal conductivity (W/m°C)}\cr}
#' \item{\code{pres}{ = 101325, atmospheric pressure (Pa)}\cr}
#' \item{\code{O2gas}{ = 20.95, oxygen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)}\cr}
#' \item{\code{N2gas}{ = 79.02, nitrogen concetration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)}\cr}
#' \item{\code{CO2gas}{ = 0.0412, carbon dioxide concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)}\cr}
#' \item{\code{PDIF}{ = 0.15, proportion of solar radiation that is diffuse (fractional, 0-1)}\cr}
#'}
#' \strong{ Morphological parameters:}
#'
#' \itemize{
#' \item{\code{custom_shape}{ = c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743), Custom shape coefficients. Operates if shape=5, and consists of 4 pairs of values representing the parameters a and b of a relationship AREA=a*mass^b, where AREA is in cm2 and mass is in g. The first pair are a and b for total surface area, then a and b for ventral area, then for sillhouette area normal to the sun, then sillhouette area perpendicular to the sun}\cr}
#' \item{\code{shape_a}{ = 1, Proportionality factor (-) for going from volume to area, keep this 1 (redundant parameter that should be removed)}\cr}
#' \item{\code{shape_b}{ = 3, Proportionality factor (-) for going from volume to area, represents ratio of width:height for a plate, length:diameter for cylinder, b axis:a axis for ellipsoid }\cr}
#' \item{\code{shape_c}{ = 0.6666666667, Proportionality factor (-) for going from volume to area, represents ratio of length:height for a plate, c axis:a axis for ellipsoid}\cr}
#' \item{\code{fatosk}{ = 0.4, Configuration factor to sky (-) for infrared calculations}\cr}
#' \item{\code{fatosb}{ = 0.4, Configuration factor to subsrate for infrared calculations}\cr}
#' \item{\code{pct_cond}{ = 10, Percentage of animal surface contacting the substrate (\%)}\cr}
#' \item{\code{pct_touch}{ = 0, Percentage of animal surface area contacting another animal of same temperature (\%)}\cr}
#' \item{\code{k_flesh}{ = 0.5, Thermal conductivity of flesh (W/mC, range: 0.412-2.8)}\cr}
#' \item{\code{rho_body}{ = 1000, Density of flesh (kg/m3)}\cr}
#' \item{\code{epsilon}{ = 0.95, Emissivity of animal (0-1)}\cr}
#' \item{\code{postur}{ = 1, postural orientation to sun, 1 = perpendicular, 2 = parallel, 0 = half way between (foraging)}\cr}
#'}
#' \strong{Outputs:}
#'
#' temperature variables:
#'
#' TC - Body temperature (°C)
#'
#' TSKIN - Skin temperature (°C)
#'
#' TLUNG - Lung temperature (°C)
#'
#' enbal variables:
#' \itemize{
#' \item 1 QSOL - Solar radiation absorbed (W)
#' \item 2 QIRIN - Infrared radiation absorbed (W)
#' \item 3 QMET - Metabolic heat production (W)
#' \item 4 QEVAP - Evaporative heat loss (W)
#' \item 5 QIROUT - Infrared radiation lost (W)
#' \item 6 QCONV - Heat lost by convection (W)
#' \item 7 QCOND - Heat lost by conduction (W)
#' \item 8 ENB - Energy balance (W)
#' \item 9 NTRY - Iterations that were required for solution to heat balance equation
#'}
#' masbal variables:
#' \itemize{
#' \item 1 O2_ml - Oxygen consumption rate (ml/h)
#' \item 2 H2OResp_g - Respiratory water loss (g/h)
#' \item 3 H2OCut_g - Cutaneous water loss (g/h)
#' \item 4 H2OEye_g - Ocular water loss (g/h)
#'}
#' @examples
#' library(NicheMapR)
#'
#' # parameters
#' Ww_g <- 40
#' shape <- 3
#' alpha <- 0.85
#' ectoR.out <- ectoR_devel(Ww_g = Ww_g, # wet weight, g
#'                          shape = shape, # using lizard geometry
#'                          alpha = alpha, # solar absorptivity
#'                          postur = 0, # average posture, half way between normal and parallel to sun
#'                          TA = 20, # air temperature at lizard height, deg C
#'                          TGRD = 40, # ground temperature, deg C
#'                          TSKY = -5, # sky temperature, deg C
#'                          VEL = 1, # wind speed, m/s
#'                          RH = 30, # relative humidity, %
#'                          QSOLR = 800, # total horizontal plane solar radiation, W/m2
#'                          Z = 20 # solar zenith angle, degrees
#' )
#' # return body temperature
#' ectoR.out$TC
#' # return skin temperature
#' ectoR.out$TS
#' # return heat budget
#' ectoR.out$enbal
#' # return O2 consumption rate and water loss
#' ectoR.out$masbal
#'
#' # run microclimate model in monthly mode at default site (Madison, Wisconsin, USA)
#' micro <- micro_global()
#'
#' # extract full sun conditions
#' metout <- as.data.frame(micro$metout)
#' soil <- as.data.frame(micro$soil)
#'
#' # get required inputs
#' TAs <- metout$TALOC
#' TGRDs <- soil$D0cm
#' TSKYs <- metout$TSKYC
#' VELs <- metout$VLOC
#' RHs <- metout$RHLOC
#' QSOLRs <- metout$SOLR
#' Zs <- metout$ZEN
#'
#' # use ectoR_devel to compute body temperature in open without respiratory heat loss,
#' # conduction, or metabolic heat gain
#' TC <- unlist(lapply(1:length(TAs), function(x){ectoR_devel(
#'   Ww_g = Ww_g, # wet weight, g
#'   shape = shape, # using lizard geometry
#'   alpha = alpha, # solar absorptivity
#'   M_1 = 0, # turn of metabolic heat
#'   postur = 0, # average posture, half way between normal and parallel to sun
#'   pantmax = 0, # turn off respiratory heat exchange
#'   pct_cond = 0, # negligible conduction to substrate
#'   alpha_sub = (1 - micro$REF), # substrate absorptivity used in microclimate model
#'   elev = micro$elev, # elevation from microclimate model
#'   TA = TAs[x], # air temperature at lizard height from microclimate model, deg C
#'   TGRD = TGRDs[x], # ground temperature from microclimate model, deg C
#'   TSKY = TSKYs[x], # sky temperature from microclimate model, deg C
#'   VEL = VELs[x], # wind speed from microclimate model, m/s
#'   RH = RHs[x], # relative humidity from microclimate model, %
#'   QSOLR = QSOLRs[x], # total horizontal plane solar radiation from microclimate model, W/m2
#'   Z = Zs[x] # solar zenith angle from microclimate model, degrees
#' )$TC})) # run ectoR_devel across environments
#'
#' # run ectotherm model for a non-behaving animal without respiratory heat loss,
#' # conduction or metabolic heat gain
#' ecto <- ectotherm(
#'   Ww_g = Ww_g, # wet weight, g
#'   shape = shape, # using lizard geometry
#'   alpha_min = alpha, # minimum solar absorptivity
#'   alpha_max = alpha, # maximum solar absorptivity
#'   M_1 = 0, # turn of metabolic heat
#'   postur = 0, # average posture, half way between normal and parallel to sun
#'   pantmax = 0, # turn off respiratory heat exchange
#'   pct_cond = 0, # negligible conduction to substrate
#'   live = 0
#' )
#'
#' # extract results
#' environ <- as.data.frame(ecto$environ)
#' enbal <- as.data.frame(ecto$enbal)
#' masbal <- as.data.frame(ecto$masbal)
#' TC_ectotherm <- environ$TC
#'
#' # compare
#' time <- micro$dates
#' plot(time, TC, type = 'l', ylim = c(-25, 55), ylab = 'body temperature, deg C', xlab = 'month of year')
#' points(time, TC_ectotherm, type = 'l', col = 2)
#'
#' @export
ectoR_devel <- function(
    Ww_g = 40,
    alpha = 0.85,
    epsilon = 0.95,
    rho_body = 1000,
    fatosk = 0.4,
    fatosb = 0.4,
    shape = 3,
    shape_a = 1,
    shape_b = 3,
    shape_c = 2 / 3,
    custom_shape = c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743),
    pct_cond = 10,
    pct_touch = 0,
    postur = 1,
    k_flesh = 0.5,
    M_1 = 0.013,
    M_2 = 0.8,
    M_3 = 0.038,
    pct_wet = 0.1,
    pct_eyes = 0,
    pct_mouth = 0,
    psi_body = -700,
    pantmax = 1,
    F_O2 = 20,
    RQ = 0.8,
    delta_air = 0.1,
    elev = 0,
    alpha_sub = 0.2,
    epsilon_sub = 1,
    epsilon_sky = 1,
    pres = 101325,
    fluid = 0,
    O2gas = 20.95,
    CO2gas = 0.03,
    N2gas = 79.02,
    K_sub = 0.5,
    PDIF = 0.1,
    SHADE = 0,
    QSOLR = 1000,
    Z = 20,
    TA = 20,
    TGRD = 30,
    TSUBST = 30,
    TSKY = -5,
    VEL = 1,
    RH = 5
    ){ # end function parameters

  errors <- 0

  # error trapping
  if(shape < 0 | shape > 5  | shape%%1 != 0){
    message("error: shape can only be an integer from 0 to 5 \n")
    errors<-1
  }
  if(alpha < 0 | alpha > 1){
    message("error: alpha can only be from 0 to 1 \n")
    errors<-1
  }
  if(pct_wet < 0 | pct_wet > 100){
    message("error: pct_wet can only be from 0 to 100 \n")
    errors<-1
  }
  if(pct_eyes < 0 | pct_eyes > 100){
    message("error: pct_eyes can only be from 0 to 100 \n")
    errors<-1
  }
  if((pantmax < 0) | (pantmax > 0 & pantmax < 1)){
    message("error: pantmax should be greater than or equal to 1, or zero if you want to simluate the effect of no respiratory water loss\n")
    errors<-1
  }
  if(F_O2 < 0 | F_O2 > 100){
    message("error: F_O2 can only be from 0 to 100 \n")
    errors<-1
  }
  if(!fluid %in% c(0, 1)){
    message("error: fluid must be 0 or 1 \n")
    errors<-1
  }
  if(alpha_sub < 0 | alpha_sub > 1){
    message("error: alpha_sub can only be from 0 to 1 \n")
    errors<-1
  }
  if(shape_a < 0){
    message("error: shape_a can't be negative \n")
    errors<-1
  }
  if(shape_b < 0){
    message("error: shape_b can't be negative \n")
    errors<-1
  }
  if(shape_c < 0){
    message("error: shape_c can't be negative \n")
    errors<-1
  }
  if(fatosk < 0 | fatosk > 1){
    message("error: fatosk can only be from 0 to 1 \n")
    errors<-1
  }
  if(fatosb < 0 | fatosb > 1){
    message("error: fatosb can only be from 0 to 1 \n")
    errors<-1
  }
  if(fatosk + fatosb > 1){
    message("error: the sum of fatosb and fatosk can't exceed 1 \n")
    errors<-1
  }
  if(pct_cond < 0 | pct_cond > 100){
    message("error: pct_cond can only be from 0 to 100 \n")
    errors<-1
  }
  if(k_flesh < 0){
    message("error: k_flesh can't be negative \n")
    errors<-1
  }
  if(rho_body < 0){
    message("error: rho_body can't be negative \n")
    errors<-1
  }
  if(epsilon < 0 | epsilon > 1){
    message("error: epsilon can only be from 0 to 1 \n")
    errors<-1
  }
  if(epsilon < 0.9){
    message("warning: epsilon is rarely below 0.9 for living things \n")
    errors<-0
  }
  if(errors == 0){

    # constants
    PI <- 3.14159265 # defining to be same precision as Fortran code

    # parameter name translations from input arguments, and value conversions

    ZEN <- Z / 180 * pi
    AMASS <- Ww_g / 1000 # animal wet weight (kg)
    SKINW <- pct_wet / 100 # fractional skin wetness
    PEYES <- pct_eyes / 100 # fractional of surface area that is wet eyes
    PMOUTH <- pct_mouth / 100 # fraction of surface area that is wet mouth
    SKINT <- pct_touch / 100 # fraction of surface area that is touching another individual at the same temperature
    PTCOND <- pct_cond / 100 # fraction of surface area conducting to the ground
    ALT <- elev # 'altitude' (m) (technically correct term is elevation)
    EMISSK <- epsilon_sky # emissivity of the sky (0-1)
    EMISSB <- epsilon_sub # emissivity of the substrate (0-1)
    ABSSB <- alpha_sub # solar absorbtivity of the substrate (0-1)
    BP <- pres # barometric pressure
    RELHUM <- RH # relative humidity
    PSI_BODY <- psi_body #
    ABSAN <- alpha # animal solar absorbtivity
    EXTREF <- F_O2 # oxygen extraction efficiency (0-1)
    GEOMETRY <- shape # animal shape
    FATOSK <- fatosk # radiation configuration factor to sky
    FATOSB <- fatosb # radiation configuration factor to ground
    O2GAS <- O2gas # O2 gas concentration (%)
    CO2GAS <- CO2gas # CO2 gas concentration (%)
    N2GAS <- N2gas # N2 gas concentration (%)
    FLSHCOND <- k_flesh # flesh thermal conductivity (W/mK)
    ANDENS <- rho_body # body density (kg/m3)
    PANT <- pantmax # panting modifier, >=1 (or 0 if cutting out respiration)
    FLTYPE <- fluid # air or water? (0 or 1)
    EMISAN <- epsilon # emissivity of animal skin (0-1)
    EMISSB <- epsilon_sub # emissivity of substrate (0-1)
    EMISSK <- epsilon_sky # emissivity of sky (0-1)
    SUBTK <- K_sub # substrate thermal conductivity (W/mK)
    DELTAR <- delta_air # temperature difference between inspired and expired air
    CUSTOMGEOM <- custom_shape # parameters for customised geometry

    # unused parameters
    FATOBJ <- 0 # configuration factor to nearby object of different temp to sky and ground (e.g. warm rock, fire), not yet used
    RINSUL <- 0 # radius of insulation, not used yet

    # shape setup for custom animals
    if(shape == 3){ # lizard proportions
      shape_a <- 1
      shape_b <- 1
      shape_c <- 4
    }
    if(shape == 4){ # frog proportions
      shape_a <- 1
      shape_b <- 1
      shape_c <- 0.5
    }
    SHP <- c(shape_a, shape_b, shape_c)
    if(SHP[1] == SHP[2] & SHP[2] == SHP[3]){
      SHP[3]=SHP[3]-0.0000001
    }

    # call GEOM_ecto to get lengths, areas and volume
    GEOM.out <- GEOM_ecto(AMASS = AMASS,
                          GEOMETRY = GEOMETRY,
                          SHP = SHP,
                          CUSTOMGEOM = CUSTOMGEOM,
                          ANDENS = ANDENS,
                          SKINW = SKINW,
                          SKINT = SKINT,
                          RINSUL = RINSUL,
                          PTCOND = PTCOND,
                          PMOUTH = PMOUTH,
                          PANT = PANT)
    VOL <- GEOM.out$VOL
    AREA <- GEOM.out$AREA
    ATOT <- AREA
    AV <- GEOM.out$AV
    AT <- GEOM.out$AT
    AL <- GEOM.out$AL
    ASILN <- GEOM.out$ASILN
    ASILP <- GEOM.out$ASILP
    AEFF <- GEOM.out$AEFF
    R1 <- GEOM.out$R1
    R <- GEOM.out$R
    ASEMAJR <- GEOM.out$ASEMAJR
    BSEMINR <- GEOM.out$BSEMINR
    CSEMINR <- GEOM.out$CSEMINR

    # set silhouette area
    if(postur == 1){
      ASIL <- ASILN
    }
    if(postur == 2){
      ASIL <- ASILP
    }
    if(postur == 0){
      ASIL <- (ASILN + ASILP) / 2
    }

    # compute solar load
    SOLAR.out <- SOLAR_ecto(ATOT = AREA,
                            ASIL = ASIL,
                            AV = AV,
                            AT = AT,
                            ABSAN = ABSAN,
                            ABSSB = ABSSB,
                            FATOSK = FATOSK,
                            FATOSB = FATOSB,
                            FATOBJ = FATOBJ,
                            ZEN = ZEN,
                            QSOLR = QSOLR,
                            PDIF = PDIF,
                            SHADE = SHADE,
                            postur = postur)
    QSOLAR <- SOLAR.out$QSOLAR

    # compute infrared radiation in
    QIRIN <- RADIN_ecto(ATOT = ATOT,
                        AV = AV,
                        AT = AT,
                        FATOSK = FATOSK,
                        FATOSB = FATOSB,
                        FATOBJ = FATOBJ,
                        EMISAN = EMISAN,
                        EMISSB = EMISSB,
                        EMISSK = EMISSK,
                        TSKY = TSKY,
                        TGRD = TGRD)

    # call FUN_ecto to find a solution for core body temperature
    TC <- uniroot(function(x) FUN_ecto(AMASS = AMASS,
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
                                       x), lower = -50, upper = 100)$root

    # recalculate everything with known TC, to get all outputs (run FUN_ecto again)

    # metabolic rate and heat production
    MET.out <- MET_ecto(AMASS = AMASS,
                        XTRY = TC,
                        M_1 = M_1,
                        M_2 = M_2,
                        M_3 = M_3)
    QMETAB <-  max(0.0001, MET.out)

    # respiration
    RESP.out <- RESP_ecto(XTRY = TC,
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

    # compute skin temperature given metabolic rate

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

    # compute convection now that TSKIN is known
    CONV.out <- CONV_ecto(GEOMETRY = GEOMETRY,
                          ATOT = AREA,
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
    GEVAP <- RESP.out$GEVAP

    # compute cutaneous evaporation now that mass transfer coefficient, HD, and TSKIN are known
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

    # compute radiation out know that TSKIN is known
    QIROUT <- RADOUT_ecto(TSKIN = TSKIN,
                          ATOT = ATOT,
                          AV = AV,
                          AT = AT,
                          FATOSK = FATOSK,
                          FATOSB = FATOSB,
                          EMISAN = EMISAN)

    # compute conduction know that TSKIN is known
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

    # balance the budget
    QIN <- QSOLAR + QIRIN + QMETAB
    QOUT <- QRESP + QSEVAP + QIROUT + QCONV + QCOND
    #C     FINDING THE DEVIATION FROM ZERO RESULTING FROM GUESSING THE SOLUTION
    ENB <- QIN - QOUT

    # outputs
    enbal <- t(matrix(c(QSOLAR, QIRIN, QMETAB, QRESP, QSEVAP, QIROUT, QCONV, QCOND, ENB)))
    colnames(enbal) <- c("QSOL", "QIRIN", "QMET", "QRESP", "QEVAP", "QIROUT", "QCONV", "QCOND", "ENB")
    O2_ml <- 10 ^ (M_3 * TC) * M_1 * Ww_g ^ M_2 # ml/h
    H2OResp_g <- SEVAP.out$WRESP * 3600 # g/h
    H2OCut_g <- SEVAP.out$WCUT * 3600 # g/h
    H2OEyes_g <- SEVAP.out$WEYES * 3600 # g/h
    masbal <- t(matrix(c(O2_ml, H2OResp_g, H2OCut_g, H2OEyes_g)))
    colnames(masbal) <- c("O2_ml", "H2OResp_g", "H2OCut_g", "H2OEyes_g")

    return(list(TC = TC, TSKIN = TSKIN, TLUNG = TLUNG, enbal = enbal, masbal = masbal, GEOM.out = GEOM.out, SOLAR.out = SOLAR.out, RESP.out = RESP.out, CONV.out = CONV.out, SEVAP.out = SEVAP.out))
  } # end error check
}
