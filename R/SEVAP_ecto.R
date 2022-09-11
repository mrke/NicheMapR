#' SEVAP_ecto
#'
#' R version of Fortran SEVAP.f (ectotherm model)
#' @encoding UTF-8
#' @param TC core body temperature (°C)
#' @param TSKIN skin temperature (°C)
#' @param GEVAP respiratory water loss rate (just used for getting total water loss) (kg/s)
#' @param PSI_BODY water potential of body (J/kg) - affects skin humidity for water vapour exchange
#' @param SKINW fraction of total surface area acting as a free water surface (fractional, 0-1)
#' @param AEFF effective area acting as a free-water exchanger, drives cutaneous evaporation (m2)
#' @param ATOT total body surface area (m2)
#' @param HD mass transfer coefficient calculated with CONV_ecto (m/s)
#' @param PEYES proportion of total surface area that is 'wet' eye (fractional, 0-1)
#' @param TA air temperature (°C)
#' @param RELHUM relative humidity (\%)
#' @param VEL wind speed (m/s)
#' @param BP air pressure (Pa)
#' @export
SEVAP_ecto <- function(
    TC = 25,
    TSKIN = 25.1,
    GEVAP = 1.177235e-09,
    PSI_BODY = -7.07 * 100,
    SKINW = 0.1,
    AEFF = 1.192505e-05,
    ATOT = 0.01325006,
    HD = 0.02522706,
    PEYES = 0.03 / 100,
    TA = 20,
    RELHUM = 50,
    VEL = 0.1,
    BP = 101325){

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
  # C     THIS SUBROUTINE COMPUTES SKIN EVAPORATION BASED ON THE MASS TRANSFER
  # C     COEFFICIENT, % OF SURFACE OF THE SKIN ACTING AS A FREE WATER SURFACE
  # C     AND EXPOSED TO THE AIR, AND THE VAPOR DENSITY GRADIENT BETWEEN THE
  # C     SURFACE AND THE AIR, EACH AT THEIR OWN TEMPERATURE.

  TAIR <- TA
  V <- VEL
  XTRY <- TC
  MW <- 0.018 #! molar mass of water, kg/mol
  RG <- 8.314 #! gas constant, J/mol/K
  #C     CALCULATING SKIN SURFACE SATURATION VAPOR DENSITY
  #C      RH = 100.
  RH <- exp(PSI_BODY / (RG / MW * (TSKIN + 273.15))) * 100 #
  DB <- TSKIN

  #C     SETTING 3 PARAMETERS FOR WETAIR, SINCE RH IS KNOWN (SEE WETAIR LISTING)
  WB <- 0
  DP <- 999
  #C     BP CALCULATED FROM ALTITUDE USING THE STANDARD ATMOSPHERE
  #C     EQUATIONS FROM SUBROUTINE DRYAIR    (TRACY ET AL,1972)
  PSTD <- 101325
  #C     PATMOS=PSTD*((1.-(.0065*ALT/288.))**(1./.190284))
  PATMOS <- BP
  BP <- PATMOS

  WETAIR.out <- WETAIR(DB, WB, RH, DP, BP)
  VDSURF <- WETAIR.out$vd

  #C     AIR VAPOR DENSITY
  RH <- RELHUM
  DB <- TAIR

  WETAIR.out <- WETAIR(DB, WB, RH, DP, BP)
  VDAIR <- WETAIR.out$vd

  WEYES <- HD * PEYES * ATOT * (VDSURF - VDAIR)

  WRESP <- GEVAP / 1000

  if(WEYES > 0){
    WCUT <- (AEFF - PEYES * ATOT * SKINW) * HD * (VDSURF - VDAIR)
  }else{
    WCUT <- AEFF * HD * (VDSURF - VDAIR)
  }
  WATER <- WEYES + WRESP + WCUT
  #C     END OF COMPUTING AEFF FOR SURFACE OR NOT

  #10 CONTINUE

  #C     FROM DRYAIR: LATENT HEAT OF VAPORIZATION
  HTOVPR <- 2.5012E+06 - 2.3787E+03 * TAIR
  QSEVAP <- (WEYES + WCUT) * HTOVPR
  #C     KG/S TO G/S
  WEYES <- WEYES * 1000
  WRESP <- WRESP * 1000
  WCUT <- WCUT * 1000
  WEVAP <- WATER * 1000

  return(list(QSEVAP = QSEVAP, WEVAP = WEVAP, WRESP = WRESP, WCUT = WCUT, WEYES = WEYES))
}
