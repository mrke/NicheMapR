#' IRPROP
#'
#' R wrapper for Fortran binary of IRPROP (endotherm model)
#' @encoding UTF-8
#' @param TA air temperature at local height (°C)
#' @param DHAIRD hair diameter, dorsal (m)
#' @param DHAIRV hair diameter, ventral (m)
#' @param LHAIRD hair length, dorsal (m)
#' @param LHAIRV hair length, ventral (m)
#' @param ZFURD fur depth, dorsal (m)
#' @param ZFURV fur depth, ventral (m)
#' @param RHOD hair density, dorsal (1/m2)
#' @param RHOV hair density, ventral (1/m2)
#' @param REFLD fur reflectivity dorsal (fractional, 0-1)
#' @param REFLV fur reflectivity ventral (fractional, 0-1)
#' @param ZFURCOMP depth of compressed fur (for conduction) (m)
#' @param PVEN fraction of surface area that is ventral fur (fractional, 0-1)
#' @param KHAIR hair thermal conductivity (W/m°C)
#' @export
IRPROP <- function(TA, DHAIRD, DHAIRV, LHAIRD, LHAIRV, ZFURD, ZFURV, RHOD, RHOV, REFLD, REFLV, ZFURCOMP, PVEN, KHAIR){
  a <- .Fortran("IRPROP",
                as.double(TA),
                as.double(DHAIRD),
                as.double(DHAIRV),
                as.double(LHAIRD),
                as.double(LHAIRV),
                as.double(ZFURD),
                as.double(ZFURV),
                as.double(RHOD),
                as.double(RHOV),
                as.double(REFLD),
                as.double(REFLV),
                as.double(ZFURCOMP),
                as.double(PVEN),
                as.double(KHAIR),
                results=matrix(data = 0., nrow = 1, ncol = 26),
                PACKAGE = "NicheMapR")

  results <- matrix(data = 0., nrow = 1, ncol = 26)

  storage.mode(results)<-"double"
  results <- a$results

  return (results)
}
