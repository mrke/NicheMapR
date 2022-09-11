#' RADIN_ecto
#'
#' R version of Fortran RADIN.f (ectotherm model) for calculating absorbed long-wave radiation.
#'
#' @encoding UTF-8
#' @param ATOT total body surface area (m2)
#' @param AV ventral surface area (m2)
#' @param AT body surface area contacting another organism of same temperature (m2)
#' @param FATOSK configuration factor to sky (-)
#' @param FATOSB configuration factor to substrate (-)
#' @param FATOBJ configuration factor to nearby object (-) (not functional at the moment)
#' @param EMISAN emissivity of animal (fractional, 0-1)
#' @param EMISSB emissivity of substrate (fractional, 0-1)
#' @param EMISSK emissivity of sky (fractional, 0-1)
#' @param TSKY sky temperature (°C)
#' @param TGRD ground temperature (°C)
#' @export
RADIN_ecto <- function(
    ATOT = 0.01325006,
    AV = 0.001325006,
    AT = 0,
    FATOSK = 0.4,
    FATOSB = 0.4,
    FATOBJ = 0,
    EMISAN = 0.95,
    EMISSB = 0.95,
    EMISSK = 0.8,
    TSKY = 10,
    TGRD = 30){

  # C     NICHEMAPR: SOFTWARE FOR BIOPHYSICAL MECHANISTIC NICHE MODELLING
  #
  # C     COPYRIGHT (C) 2018 MICHAEL R. KEARNEY AND WARREN P. PORTER
  #
  # C     THIS PROGRAM IS FREE SOFTWARE: YOU CAN REDISTRIBUTE IT AND/OR MODIFY
  # C     IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
  # C     THE FREE SOFTWARE FOUNDATION, EITHER VERSION 3 OF THE LICENSE, OR (AT
  #C      YOUR OPTION) ANY LATER VERSION.
  #
  # C     THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT
  # C     WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
  # C     MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. SEE THE GNU
  # C     GENERAL PUBLIC LICENSE FOR MORE DETAILS.
  #
  # C     YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
  # C     ALONG WITH THIS PROGRAM. IF NOT, SEE HTTP://WWW.GNU.ORG/LICENSES/.
  #
  # C     COMPUTES LONGWAVE RADIATION ABSORBED
  SIG <- 5.669E-08

  TKSKY <- TSKY + 273.15
  TKSUB <- TGRD + 273.15
  TKOBJ <- TGRD + 273.15
  FSKY <- FATOSK - FATOBJ
  if(FSKY < 0){
    FSKY <- 0
  }
  QIROBJ <- EMISAN * FATOBJ * (ATOT-AT / 2) * EMISSB * SIG * TKOBJ ^ 4
  QIRSKY <- EMISAN * FSKY * (ATOT-AT / 2) * EMISSK * SIG * TKSKY ^ 4
  QIRSUB <- EMISAN * FATOSB * (ATOT-AV-AT / 2) * EMISSB * SIG * TKSUB ^ 4
  QIRIN <- QIRSKY + QIRSUB + QIROBJ
  return(QIRIN)
}

