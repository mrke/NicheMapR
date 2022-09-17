#' SOLAR_ecto
#'
#' R version of Fortran SOLAR.f (ectotherm model) for calculating absorbed solar radiation.
#'
#' @encoding UTF-8
#' @param ATOT surface area for solar exchange, m2
#' @param ASIL silouette area for solar exchange, m2
#' @param AV ventral area for solar exchange with substrate, m2
#' @param AT total surface area touching another animal, m2
#' @param ABSAN solar absorptivity of animal (fractional, 0-1)
#' @param ABSSB solar absorptivity of substrate (fractional, 0-1)
#' @param FATOSK radiation configuration factor to sky (fractional, 0-1)
#' @param FATOBJ radiation configuration factor to ground (fractional, 0-1)
#' @param ZEN solar zenith angle (radians)
#' @param PDIF proportion of solar radiation that is diffuse (fractional, 0-1)
#' @param SHADE shade (\%)
#' @param postur postural orientation to sun, 1 = perpendicular, 2 = parallel, 0 = half way between
#' @export
SOLAR_ecto <- function(
    ATOT = 0.01325006,
    ASIL= 0.004718043,
    AV = 0.001325006,
    AT = 0,
    ABSAN = 0.85,
    ABSSB = 0.8,
    FATOSK = 0.4,
    FATOSB = 0.4,
    FATOBJ = 0,
    ZEN = 20 * pi / 180,
    QSOLR = 1000,
    PDIF = 0.1,
    SHADE = 0,
    postur = 1){
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
  # C     COMPUTES SOLAR RADIATION ABSORBED

  PI <- 3.14159265

  #C     CHECKING FOR SCATTERED SKYLIGHT ONLY WHEN SUN BELOW HORIZON
  ZENITH <- ZEN * 180 / PI
  #C     DIRECT BEAM COMPONENT
  if(ZENITH < 90){
    #C      DIRECT BEAM (NORMAL TO THE DIRECT BEAM)
    QNORM <- (QSOLR / cos(ZEN))
    if(postur != 1){
      QNORM <- QSOLR
    }
    if(QNORM > 1367){
      #C       MAKING SURE THAT LOW SUN ANGLES DON'T LEAD TO SOLAR VALUES
      #C       GREATER THAN THE SOLAR CONSTANT
      QNORM <- 1367
    }
    if(ZENITH >= 90){
      QNORM <- 0
    }
    QSDIR <- ABSAN * ASIL * (1 - PDIF) * QNORM * (1 - (SHADE / 100))
  }else{
    QSDIR <- 0.00
    QNORM <- 0.00
  }

  #C      DIFFUSE COMPONENTS (SKY AND SUBSTRATE)
  QSOBJ <- ABSAN * FATOBJ * (ATOT - AT / 2) * PDIF * QNORM
  QSSKY <- ABSAN * FATOSK * (ATOT - AT / 2) * PDIF * QNORM * (1 - (SHADE / 100))
  QSRSB <- ABSAN * FATOSB * (ATOT - AV - AT / 2) * (1 - ABSSB) * QSOLR * (1 - (SHADE / 100))
  QSDIFF <- QSSKY + QSRSB + QSOBJ
  QSOLAR <- QSDIR + QSDIFF

  return(list(QSOLAR = QSOLAR, QSDIFF = QSDIFF, QSRSB = QSRSB, QSSKY = QSSKY, QSOBJ = QSOBJ))
}
