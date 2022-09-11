#' RADOUT_ecto
#'
#' R version of Fortran RADOUT.f (ectotherm model) for calculating outgoing long-wave radiation.
#'
#' @encoding UTF-8
#' @param TSKIN skin temperature (°C)
#' @param ATOT total body surface area (m2)
#' @param AV ventral surface area (m2)
#' @param AT body surface area contacting another organism of same temperature (m2)
#' @param FATOSK configuration factor to sky (-)
#' @param FATOSB configuration factor to substrate (-)
#' @param EMISAN emissivity of animal (fractional, 0-1)
#' @export
RADOUT_ecto <- function(
    TSKIN = 25.1,
    ATOT = 0.01325006,
    AV = 0.001325006,
    AT = 0,
    FATOSK = 0.4,
    FATOSB = 0.4,
    EMISAN = 0.95){

# C     NICHEMAPR: SOFTWARE FOR BIOPHYSICAL MECHANISTIC NICHE MODELLING
#
# C     COPYRIGHT (C) 2018 MICHAEL R. KEARNEY AND WARREN P. PORTER
#
# C     THIS PROGRAM IS FREE SOFTWARE: YOU CAN REDISTRIBUTE IT AND/OR MODIFY
# C     IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
# C     THE FREE SOFTWARE FOUNDATION, EITHER VERSION 3 OF THE LICENSE, OR (AT
#                                                                          C      YOUR OPTION) ANY LATER VERSION.
#
# C     THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT
# C     WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
# C     MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. SEE THE GNU
# C     GENERAL PUBLIC LICENSE FOR MORE DETAILS.
#
# C     YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
# C     ALONG WITH THIS PROGRAM. IF NOT, SEE HTTP://WWW.GNU.ORG/LICENSES/.
#
# C     COMPUTES LONGWAVE RADIATION LOST
SIG <- 5.669E-08
X <- TSKIN
XK <- X + 273.15
QIR2SK <- (ATOT - AT / 2) * FATOSK * EMISAN * SIG * XK ^ 4
QIR2SB <- (ATOT - AV - AT / 2) * FATOSB * EMISAN * SIG * XK ^ 4
QIROUT <- QIR2SK + QIR2SB

return(QIROUT)
}
