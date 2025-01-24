#' HomoTherm_var - the human model of NicheMapR
#'
#' This is a multi-part application of the endoR model to a human. It simulates runs
#' the HomoTherm function across a sequence of environments.
#' @encoding UTF-8
#' @param MASS = 70, mass of person (kg)
#' @param QMETAB_REST = 105, resting metabolic rate (W)
#' @param ACTIVE = FALSE, activity state (-)
#' @param MET = 1, MET units of activity (-)
#' @param INSDEPDs = c(1e-02, rep(6e-03, 3)), clothing depth, dorsal (m)
#' @param INSDEPVs = c(1e-09, rep(6e-03, 3)), clothing depth, ventral (m)
#' @param TAs = 21, air temperature at local height (°C)
#' @param TGRDs = TAs, ground temperature (°C)
#' @param TSKYs = TAs, sky temperature (°C)
#' @param VELs = 0.1, wind speed (m/s)
#' @param RHs = 50, relative humidity (\%)
#' @param QSOLRs = 0, solar radiation, horizontal plane (W/m2)
#' @param Zs = 20, zenith angle of sun (degrees from overhead)
#' @usage HomoTherm_var(MASS = 70, QMETAB_REST = 105, ACTIVE = FALSE, MET = 1, INSDEPDs = c(1e-02, rep(6e-03, 3)), INSDEPVs = c(1e-09, rep(6e-03, 3)), TAs = 21, TGRDs = TAs, TSKYs = TAs, VELs = 0.1, RHs = 50, QSOLRs = 0, Zs = 20,...)
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
#' \code{Refhyt}{ = 2, height of weather observations (m)}\cr\cr
#' \code{RUF}{ = 0.004, roughness height for calculating air/wind/humidity profile with get_profile (m)}\cr\cr
#' \code{ZH}{ = 1, heat transfer roughness height for calculating air/wind/humidity profile with get_profile (m) (-)}\cr\cr
#' \code{D0}{ = 1, zero plane displacement correction factor for calculating air/wind/humidity profile with get_profile (m)}\cr\cr
#' \code{CONV_ENHANCE}{ = 1, convective enhancement factor (> 1 for turbulent outdoor conditions) (-)}\cr\cr
#'
#' \strong{ Whole body parameters:}\cr\cr
#' \code{MAXSWEAT}{ = 1500, maximum sweating rate (g/h)}\cr\cr
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
#' \code{FATPCT}{ = c(5, 36, 10, 23) * 0.5, \% body fat}\cr\cr
#' \code{SHAPE_Bs}{ = c(1.6, 1.9, 11, 7.0), ratio between long and short axis (-)}\cr\cr
#' \code{FSKREFs}{ = c(0.50, 0.42, 0.35, 0.35), configuration factor to sky}\cr\cr
#' \code{FGDREFs}{ = c(0.38, 0.42, 0.35, 0.35), reference configuration factor to ground}\cr\cr
#' \code{EMISANs}{ = rep(0.95, 4), emissivity each body part (-)}\cr\cr
#' \code{REFLD}{ = rep(0.3, 4), solar reflectivity dorsal (fractional, 0-1)}\cr\cr
#' \code{REFLV}{ = rep(0.3, 4), solar reflectivity ventral (fractional, 0-1)}\cr\cr
#' \code{heights}{ = rep(NA, 4), height of mid-point of each body part (m), can be calculated with the 'get_heights' function}\cr\cr
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
#' TAs <- seq(-5, 50) # sequence of air temperatures, °C
#' VELs <- rep(0.1, length(TAs)) # keep wind speeds constant, m/s
#' RHs <- rep(50, length(TAs)) # keep humidity constant, %
#' # set insulation depth, flesh conductivity and fat
#' INSDEPDs <- c(1e-02, rep(6.15e-03, 3)) # 'dorsal' clothing depth, m
#' INSDEPVs <- c(1e-09, rep(6.15e-03, 3)) # 'ventral' clothing depth, m
#' KCLOs <- rep(0.04, 4) # clothing thermal conductivity, W/m·K
#' FATPCTs <- c(5, 36, 10, 23) # body fat %
# simulate across the air temperatures
#' HomoTherm.out <- HomoTherm_var(INSDEPDs = INSDEPDs * 0,
#'                                INSDEPVs = INSDEPVs * 0,
#'                                KCLOs = KCLOs,
#'                                FATPCTs = FATPCTs,
#'                                TAs = TAs,
#'                                VELs = VELs,
#'                                RHs = RHs,
#'                                EXCEED.TCMAX = TRUE)
#'balance <- HomoTherm.out$balance
#'LCT <- TAs[which(balance$K_FLESH > balance$K_FLESH[1])[1]-1]
#'UCT <- TAs[balance$T_SKIN > 34][1]
#'plot(TAs, balance$QMETAB, type = 'l', col = 'red', lwd = 1.5, ylim = c(0, 300), ylab = 'watts', xlab = 'air temperature, °C', xaxs = 'i', yaxs = 'i')
#'points(TAs, (balance$QEVAP_RESP + balance$QEVAP_CUT) * -1, type = 'l', col = 'blue', lwd = 1.5)
#'points(TAs, balance$PCTWET, type = 'l', col = 'lightblue', lwd = 1.5)
#'points(TAs, balance$K_FLESH * 10, type = 'l', lty = 2, col = 'red')
#'legend(x = -5, y = 300, legend = c("Q_metab", "Q_evap", "% wet", "k_flesh×10"), col = c("red", "blue", "lightblue", "red"), lty = c(1, 1, 1, 2), bty = "n", cex = 0.75)
#'abline(v = LCT, lty = 2)
#'abline(v = UCT, lty = 2)
HomoTherm_var <- function(MASS = 70,
                          QMETAB_REST = 105,
                          ACTIVE = 0,
                          MET = 1,
                          Q10 = 2,
                          RQ = 0.8,
                          EXTREF = 25,
                          TC_RESTs = rep(36.8, 4),
                          TC_ACTIVEs = rep(37.5, 4),
                          TC_MAXs = rep(38, 4),
                          TC_INCs = rep(0.04, 4),
                          DENSITYs = rep(1050, 4),
                          MASSFRACs = c(0.07609801, 0.50069348, 0.04932963, 0.16227462),
                          AREAFRACs = c(0.08291887, 0.32698460, 0.11025155, 0.18479669),
                          SHAPE_Bs = c(1.6, 1.9, 11, 7.0),
                          heights = rep(NA, 4),
                          PJOINs = c(0.02753623, 0.08239728, 0.02173913, 0.03333333),
                          SUBQFATs = rep(1, 4),
                          FATPCTs = c(5, 36, 10, 23),
                          KFLESHs = c(0.9, 0.9, 0.5, 0.5),
                          KFLESH_MAXs = rep(5, 4),
                          KFLESH_INCs = rep(0.05, 4), # rep(0.2, 4),
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
                          MAXSWEATs = rep(1500, 4),
                          CLOWETs = rep(0, 4),
                          EXCEED.TCMAX = TRUE,
                          CONV_ENHANCE = 1,
                          Refhyt = 2,
                          RUF = 0.004,
                          ZH = 0,
                          D0 = 0,
                          TAs = 21,
                          TAREFs = TAs,
                          TSKYs = TAs,
                          TGRDs = TAs,
                          VELs = rep(1, length(TAs)),
                          VREFs = VELs,
                          SPEED = 0,
                          RHs = rep(50, length(TAs)),
                          RHREFs = RHs,
                          QSOLRs = rep(0, length(TAs)),
                          SHADEs = rep(0, length(TAs)),
                          PDIFs = rep(0.15, length(TAs)),
                          Zs = rep(20, length(TAs)),
                          ELEV = 0,
                          BPs = rep(101325, length(TAs)),
                          GRAV = 9.80665,
                          ABSSB = 0.85,
                          O2GAS = 20.95,
                          N2GAS = 79.02,
                          CO2GAS = 0.0422,
                          SIL_adjust = TRUE){

  NPARTs <- c(1, 1, 2, 2) # one head, one trunk two arms and two legs

  # avoid body parts with mixtures of fur and no fur
  INSDEPDs[INSDEPDs == 0 & INSDEPVs > 0] <- 1e-9
  INSDEPVs[INSDEPVs == 0 & INSDEPDs > 0] <- 1e-9

  check_same_length <- function(...) {
    vectors <- list(...)
    lengths <- sapply(vectors, length)  # Get lengths of all vectors
    all(lengths == lengths[1])          # Check if all lengths are equal
  }
  if(!check_same_length(TAs, TAREFs, TSKYs, TGRDs, RHs, RHREFs, VELs, VREFs, QSOLRs, PDIFs, Zs, BPs, SHADEs)){
    message('warning, not all environmental vectors are the same length - will now check if just TAs is a vector and make all others the same length')
    # check if user has entered a vector of air temps and just a single
    # value for other variables - fix
    if(length(TAs) > 1 & length(TAREFs == 1)){
      TAREFs <- rep(TAREFs, length(TAs))
    }
    if(length(TAs) > 1 & length(TSKYs == 1)){
      TSKYs <- rep(TSKYs, length(TAs))
    }
    if(length(TAs) > 1 & length(TSKYs == 1)){
      TSKYs <- rep(TSKYs, length(TAs))
    }
    if(length(TAs) > 1 & length(TGRDs == 1)){
      TGRDs <- rep(TGRDs, length(TAs))
    }
    if(length(TAs) > 1 & length(RHs == 1)){
      RHs <- rep(RHs, length(TAs))
    }
    if(length(TAs) > 1 & length(RHREFs == 1)){
      RHREFs <- rep(RHREFs, length(TAs))
    }
    if(length(TAs) > 1 & length(VELs == 1)){
      VELs <- rep(VELs, length(TAs))
    }
    if(length(TAs) > 1 & length(VREFs == 1)){
      VREFs <- rep(VREFs, length(TAs))
    }
    if(length(TAs) > 1 & length(QSOLRs == 1)){
      QSOLRs <- rep(QSOLRs, length(TAs))
    }
    if(length(TAs) > 1 & length(PDIFs == 1)){
      PDIFs <- rep(PDIFs, length(TAs))
    }
    if(length(TAs) > 1 & length(Zs == 1)){
      Zs <- rep(Zs, length(TAs))
    }
    if(length(TAs) > 1 & length(BPs == 1)){
      BPs <- rep(BPs, length(TAs))
    }
    if(length(TAs) > 1 & length(SHADEs == 1)){
      SHADEs <- rep(SHADEs, length(TAs))
    }
  }

  if(SIL_adjust){
    # work out silhouette area relation and correction factor
    get_silhouette <- function(ZENs = seq(0, 90, 1), # degrees, zenith angles
                               AZIMUTHs = rep(0, length(ZENs)),
                               MASSs = c(5.32, 35.07, 3.43, 11.34), # kg, masses per part
                               DENSITYs = rep(1050, 4), # kg/m3, densities per part
                               INSDEPDs = c(1e-02, rep(6e-03, 3)),
                               SHAPE_Bs = c(1.75, 1.87, 6.65, 6.70),
                               SUBQFATs = rep(1, 4),
                               FATPCTs = c(5, 36, 10, 23),
                               PJOINs = c(0.0348, 0.0965, 0.0345, 0.0347),
                               plot.sil = TRUE){

      SHAPEs <- c(4, 1, 1, 1)
      ORIENTs <- rep(3, 4)
      for(j in 1:length(ZENs)){
        GEOM.head <- c(ZENs[j], GEOM_ENDO(MASSs[1], DENSITYs[1], FATPCTs[1], SHAPEs[1], INSDEPDs[1], SUBQFATs[1], SHAPE_Bs[1], SHAPE_Bs[1], 0, 0, PJOINs[1], 0, ORIENTs[1], ZEN = ZENs[j]))
        GEOM.trunk <- c(ZENs[j], GEOM_ENDO(MASSs[2], DENSITYs[2], FATPCTs[2], SHAPEs[2], INSDEPDs[2], SUBQFATs[2], SHAPE_Bs[2], SHAPE_Bs[2], 0, 0, PJOINs[2], 0, ORIENTs[2], ZEN = ZENs[j]))
        GEOM.arm <- c(ZENs[j], GEOM_ENDO(MASSs[3], DENSITYs[3], FATPCTs[3], SHAPEs[3], INSDEPDs[3], SUBQFATs[3], SHAPE_Bs[3], SHAPE_Bs[3], 0, 0, PJOINs[3], 0, ORIENTs[3], ZEN = ZENs[j]))
        GEOM.leg <- c(ZENs[j], GEOM_ENDO(MASSs[4], DENSITYs[4], FATPCTs[4], SHAPEs[4], INSDEPDs[4], SUBQFATs[4], SHAPE_Bs[4], SHAPE_Bs[4], 0, 0, PJOINs[4], 0, ORIENTs[4], ZEN = ZENs[j]))
        if(j == 1){
          GEOM.heads <- GEOM.head
          GEOM.trunks <- GEOM.trunk
          GEOM.arms <- GEOM.arm
          GEOM.legs <- GEOM.leg
        }else{
          GEOM.heads <- rbind(GEOM.heads, GEOM.head)
          GEOM.trunks <- rbind(GEOM.trunks, GEOM.trunk)
          GEOM.arms <- rbind(GEOM.arms, GEOM.arm)
          GEOM.legs <- rbind(GEOM.legs, GEOM.leg)
        }
      }
      GEOM.heads <- as.data.frame(GEOM.heads)
      GEOM.trunks <- as.data.frame(GEOM.trunks)
      GEOM.arms <- as.data.frame(GEOM.arms)
      GEOM.legs <- as.data.frame(GEOM.legs)

      GEOM.lab <- c("ZEN", "VOL", "D", "MASFAT", "VOLFAT", "ALENTH", "AWIDTH", "AHEIT", "ATOT", "ASIL", "ASILN", "ASILP", "GMASS", "AREASKIN", "FLSHVL", "FATTHK", "ASEMAJ", "BSEMIN", "CSEMIN", "CONVSK", "CONVAR", "R1", "R2")

      colnames(GEOM.heads) <- GEOM.lab
      colnames(GEOM.trunks) <- GEOM.lab
      colnames(GEOM.arms) <- GEOM.lab
      colnames(GEOM.legs) <- GEOM.lab
      AREA <- max(GEOM.heads$ATOT + GEOM.trunks$ATOT + GEOM.arms$ATOT * 2 + GEOM.legs$ATOT * 2)

      # compare with Underwood and Ward calculation
      thetas <- (90 - ZENs) * pi / 180
      psis <- AZIMUTHs * pi / 180
      sil.underwood <- (0.043 * sin(thetas) + 2.997 * cos(thetas) * (0.02133 * cos(psis) ^ 2 + 0.0091 * sin(psis)^2) ^ 0.5) / 1.81 * AREA
      sil.endoR <- GEOM.heads$ASIL + GEOM.trunks$ASIL + GEOM.arms$ASIL * 2 + GEOM.legs$ASIL * 2
      if(plot.sil){
        par(mfrow = c(2, 1))
        par(oma = c(2, 1, 2, 2) + 0.1)
        par(mar = c(3, 3, 1.5, 1) + 0.1)
        par(mgp = c(2, 1, 0))
        plot(ZENs, sil.endoR, type = 'l', ylim = c(0, max(sil.endoR)), ylab = 'silhouette area, m2', xlab = 'zenith angle, degrees')
        points(ZENs, sil.underwood, type = 'l', lty = 2)
        legend(c(0, max(sil.endoR)), cex = 0.75, legend = c('endoR', 'Underwood & Wang'), lty = c(1, 2), bty = 'n')
      }
      # model difference as a function of zenith and adjust diffuse solar fraction
      delta.calc <- sil.underwood / sil.endoR
      fit <- lm(delta.calc ~ poly(ZENs, 8))
      if(plot.sil){
        plot(ZENs, delta.calc, type = 'l', ylab = 'Underwood/endoR silhouette area', xlab = 'zenith angle, degrees')
        points(predict(fit), type = 'l', lty = 2)
        legend(max(ZENs) * 0.5, .7, cex = 0.75, legend = c('observed', 'polynomial fit'), lty = c(1, 2), bty = 'n')
        par(mfrow = c(1, 1))
      }
      return(list(fit = fit,
                  sil.underwood = sil.underwood,
                  sil.HomoTherm = sil.endoR,
                  AREA = AREA))
    }
    # adjust solar radiation for silhouette area

    # run get_silhouette function to obtain correction factor
    # to match empirical relationship in Underwood and Ward (1966)
    MASSs <- MASS * MASSFRACs
    sil.out <- get_silhouette(ZENs = seq(0, 90),
                              MASSs = MASSs,
                              DENSITYs = DENSITYs,
                              SHAPE_Bs = SHAPE_Bs,
                              FATPCTs = FATPCTs,
                              PJOINs = PJOINs,
                              plot.sil = FALSE)
    fit <- sil.out$fit # multiple regression pars for adjusting silhouette area
    AREA <- sil.out$AREA # area of human based on current pars
    adjust <- predict(fit, newdata = data.frame(ZENs = Zs))

    thetas <- (90 - Zs) * pi / 180
    psis <- 0 * pi / 180
    # Underwood and Wang model prediction
    sil.underwood <- (0.043 * sin(thetas) + 2.997 * cos(thetas) * (0.02133 * cos(psis) ^ 2 + 0.0091 * sin(psis)^2) ^ 0.5) / 1.81 * AREA
    Qnorm <- QSOLRs/cos(Zs * pi / 180) # solar normal to sun
    Qnorm[Qnorm > 1367] <- 0

    # correct calculation of direct and diffuse solar
    F_SKY <- FSKREFs[1] * AREAFRACs[1] + FSKREFs[2] * AREAFRACs[2] + FSKREFs[3] * AREAFRACs[3] * 2 + FSKREFs[4] * AREAFRACs[4] * 2
    F_GRD <- FGDREFs[1] * AREAFRACs[1] + FGDREFs[2] * AREAFRACs[2] + FGDREFs[3] * FGDREFs[3] * 2 + FGDREFs[4] * AREAFRACs[4] * 2
    Qdir <- sil.underwood * Qnorm # correct
    Qdir1 <- (sil.underwood / adjust) * Qnorm # biased
    Qdif <- F_SKY * AREA * QSOLRs * PDIFs + F_GRD * AREA * QSOLRs * (1 - ABSSB) * PDIFs
    Qsol <- Qdir + Qdif
    Qsol1 <- Qdir1 + Qdif
    SOLR_adj <- Qsol/Qsol1
    SOLR_adj[is.na(SOLR_adj)] <- 1
    QSOLRs <- QSOLRs * SOLR_adj
  }

  balance <- matrix(data = 0, nrow = length(TAs), ncol = 19)
  respire <- matrix(data = 0, nrow = length(TAs), ncol = 6)
  head.treg <- matrix(data = 0, nrow = length(TAs), ncol = 10)
  head.morph <- matrix(data = 0, nrow = length(TAs), ncol = 19)
  head.enbal <- matrix(data = 0, nrow = length(TAs), ncol = 9)
  colnames(head.treg) <- c("T_CORE", "TSKIN_D", "TSKIN_V", "TCLO_D", "TCLO_V", "PCTWET", "K_FLESH", "K_CLO_D", "K_CLO_V", "Q10")
  colnames(head.morph) <- c("MASS", "AREA", "VOLUME", "CHAR_DIMENSION", "MASS_FAT", "FAT_THICK", "FLESH_VOL", "LENGTH", "WIDTH", "HEIGHT", "R_SKIN", "R_FUR", "AREA_SILHOUETTE", "AREA_SKIN", "AREA_SKIN_EVAP", "AREA_CONV", "AREA_JOIN", "F_SKY", "F_GROUND")
  colnames(head.enbal) <- c("QSOL", "QIRIN", "QGEN", "QEVAP", "QIROUT", "QCONV", "ENB", "NTRY", "SUCCESS")

  trunk.treg <- head.treg
  trunk.morph <- head.morph
  trunk.enbal <- head.enbal
  arm.treg <- head.treg
  arm.morph <- head.morph
  arm.enbal <- head.enbal
  leg.treg <- head.treg
  leg.morph <- head.morph
  leg.enbal <- head.enbal

  for(j in 1:length(TAs)){

    TAREF <- TAREFs[j]
    if(is.na(heights[1])){
      TA <- rep(TAs[j], length(NPARTs))
      VEL <- rep(VELs[j], length(NPARTs))
      RH <- rep(RHs[j], length(NPARTs))
    }else{ # adjust for height of each body part
      profile.out <- get_profile(Refhyt = Refhyt,
                                 RUF = RUF,
                                 heights = heights,
                                 TAREF = TAREFs[j],
                                 VREF = VREFs[j],
                                 RH = RHREFs[j],
                                 D0cm = TGRDs[j],
                                 ZEN = Zs[j],
                                 warn = FALSE)
      TA <- profile.out$TAs[2:(length(NPARTs) + 1)] # first output is at ground level
      VEL <- profile.out$VELs[2:(length(NPARTs) + 1)] # first output is at ground level
      RH <- profile.out$RHs[2:(length(NPARTs) + 1)] # first output is at ground level
    }
    TSKY <- TSKYs[j]
    TGRD <- TGRDs[j]
    TBUSH <- TAREFs[j]
    QSOLR <- QSOLRs[j]
    Z <- Zs[j]
    BP <- BPs[j]
    PDIF <- PDIFs[j]
    SHADE <- SHADEs[j]
    VEL[VEL < SPEED] <- SPEED
    out <- HomoTherm(MASS = MASS,
                     QMETAB_REST = QMETAB_REST,
                     ACTIVE = ACTIVE,
                     MET = MET,
                     Q10 = Q10,
                     RQ = RQ,
                     EXTREF = EXTREF,
                     TC_RESTs = TC_RESTs,
                     TC_ACTIVEs = TC_ACTIVEs,
                     TC_MAXs = TC_MAXs,
                     EXCEED.TCMAX = EXCEED.TCMAX,
                     TC_INCs = TC_INCs,
                     DENSITYs = DENSITYs,
                     MASSFRACs = MASSFRACs,
                     AREAFRACs = AREAFRACs,
                     SHAPE_Bs = SHAPE_Bs,
                     PJOINs = PJOINs,
                     SUBQFATs = SUBQFATs,
                     FATPCTs = FATPCTs,
                     KFLESHs = KFLESHs,
                     KFLESH_MAXs = KFLESH_MAXs,
                     KFLESH_INCs = KFLESH_INCs,
                     KFATs = KFATs,
                     KCLOs = KCLOs,
                     DHAIRDs = DHAIRDs,
                     DHAIRVs = DHAIRVs,
                     LHAIRDs = LHAIRDs,
                     LHAIRVs = LHAIRVs,
                     INSDEPDs = INSDEPDs,
                     INSDEPVs = INSDEPVs,
                     INSDENDs = INSDENDs,
                     INSDENVs = INSDENVs,
                     REFLDs = REFLDs,
                     REFLVs = REFLVs,
                     EMISANs = EMISANs,
                     FGDREFs = FGDREFs,
                     FSKREFs = FSKREFs,
                     PCTWETs = PCTWETs,
                     PCTWET_INCs = PCTWET_INCs,
                     PCTWET_MAXs = PCTWET_MAXs,
                     PCTBAREVAPs = PCTBAREVAPs,
                     MAXSWEATs = MAXSWEATs,
                     CLOWETs = rep(0, length(PCTWETs)),
                     CONV_ENHANCE = CONV_ENHANCE,
                     TA = TA,
                     TAREF = TAREF,
                     TSKY = TSKY,
                     TGRD = TGRD,
                     VEL = VEL,
                     RH = RH,
                     QSOLR = QSOLR,
                     SHADE = SHADE,
                     PDIF = PDIF,
                     Z = Z,
                     ELEV = ELEV,
                     BP = BP,
                     ABSSB = ABSSB,
                     O2GAS = O2GAS,
                     N2GAS = N2GAS,
                     CO2GAS = CO2GAS)
    balance[j, ] <- out$balance
    respire[j, ] <- out$respire
    head.treg[j, ] <- out$head.treg
    head.morph[j, ] <- out$head.morph
    head.enbal[j, ] <- out$head.enbal
    trunk.treg[j, ] <- out$trunk.treg
    trunk.morph[j, ] <- out$trunk.morph
    trunk.enbal[j, ] <- out$trunk.enbal
    arm.treg[j, ] <- out$arm.treg
    arm.morph[j, ] <- out$arm.morph
    arm.enbal[j, ] <- out$arm.enbal
    leg.treg[j, ] <- out$leg.treg
    leg.morph[j, ] <- out$leg.morph
    leg.enbal[j, ] <- out$leg.enbal
  }

  balance <- as.data.frame(balance)
  colnames(balance) <- c("T_CORE", "T_LUNG", "T_SKIN", "T_CLO", "PCTWET", "K_FLESH", "EVAP_CUT_L", "EVAP_RESP_L", "SWEAT_L", "QMETAB", "QSLR", "QRAD_IN", "QRAD_OUT", "QCONV_RESP", "QEVAP_RESP", "QEVAP_CUT", "QCONV", "AREA", "AREA_RAD")
  colnames(respire) <- c("AIR_L", "O2_L", "O2_mol_in", "O2_mol_out", "AIR_mol_in", "AIR_mol_out")
  return(list(balance = balance,
              respire = as.data.frame(respire),
              head.treg = as.data.frame(head.treg),
              head.enbal = as.data.frame(head.enbal),
              head.morph = as.data.frame(head.morph),
              trunk.treg = as.data.frame(trunk.treg),
              trunk.enbal = as.data.frame(trunk.enbal),
              trunk.morph = as.data.frame(trunk.morph),
              arm.treg = as.data.frame(arm.treg),
              arm.enbal = as.data.frame(arm.enbal),
              arm.morph = as.data.frame(arm.morph),
              leg.treg = as.data.frame(leg.treg),
              leg.enbal = as.data.frame(leg.enbal),
              leg.morph = as.data.frame(leg.morph)))
}
