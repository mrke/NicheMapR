#' CONV_ENDO
#'
#' R wrapper for Fortran binary of CONV_ENDO (endotherm model)
#' @encoding UTF-8
#' @param TS skin temperature (°C)
#' @param TENV fluid temperature at local height (°C)
#' @param SHAPE shape, 1 is cylinder, 2 is sphere, 3 is plate, 4 is ellipsoid
#' @param SURFAR surface area for convection, m2
#' @param FLTYPE FLUID TYPE: 0 = AIR; 1 = FRESH WATER; 2 = SALT WATER
#' @param FURTST test of fur presence (-) from IRPROP
#' @param D characteristic dimension for convection, m
#' @param TFA fur/air interface temperature (°C)
#' @param CONV_ENHANCE convective enhancement factor (-)
#' @param VEL wind speed (m/s)
#' @param ZFUR fur depth, mean (m) (from IRPROP)
#' @param BP barometric pressure (Pa), negative means altitude is used
#' @param ELEV elevation (m)
#' @param GRAV acceleration due to gravity, (m/s^2)
#' @export
CONV_ENDO <- function(TS, TENV, SHAPE, SURFAR, FLTYPE, FURTST, D, TFA, VEL, ZFUR, BP, ELEV, CONV_ENHANCE, GRAV, O2GAS, N2GAS, CO2GAS){
  a <- .Fortran("CONV_ENDO",
    as.double(TS),
    as.double(TENV),
    as.double(SHAPE),
    as.double(SURFAR),
    as.double(FLTYPE),
    as.double(FURTST),
    as.double(D),
    as.double(TFA),
    as.double(VEL),
    as.double(ZFUR),
    as.double(BP),
    as.double(ELEV),
    as.double(CONV_ENHANCE),
    as.double(GRAV),
    as.double(O2GAS),
    as.double(N2GAS),
    as.double(CO2GAS),
    results=matrix(data = 0, nrow = 1, ncol = 14),
    PACKAGE = "NicheMapR")

  results <- matrix(data = 0, nrow = 1, ncol = 14)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
