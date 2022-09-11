#' MET_ecto
#'
#' R version of Fortran MET.f (ectotherm model) for computing metabolic heat production using an allometric function.
#'
#' @encoding UTF-8
#' @param AMASS body mass, kg
#' @param XTRY current guess of core body temperature (°C)
#' @param M_1 metabolic rate parameter 1 V_O2=M_1*M^M_2*10^(M_3*Tb), in ml O2 / h, default parameters for lizards
#' @param M_2 metabolic rate parameter 2
#' @param M_3 metabolic rate parameter 3
#' @export
MET_ecto <- function(
    AMASS = 0.04,
    XTRY = 25,
    M_1 = 0.013,
    M_2 = 0.8,
    M_3 = 0.038){

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
  # C      THIS SUBROUTINE COMPUTES METABOLIC HEAT

  GMASS <- AMASS * 1000
  TC <- XTRY

  #C      USE REGRESSION
  #C      LIZARD REGRESSION
  if(TC > 50){
    #C       CAP METABOLIC RATE EQUATION WITH MAX OF TC = 50
    QMETAB <- 0.0056 * 10 ^ (M_3 * 50)  * M_1 * GMASS ^ M_2
  }

  if(TC >= 1){
    #C        ACCEPTABLE TEMPERATURE RANGE
    QMETAB <- 0.0056 * 10 ^ (M_3 * TC) * M_1 * GMASS ^ M_2
  }else{
    QMETAB <- 0.01
  }
  return(QMETAB)
}
