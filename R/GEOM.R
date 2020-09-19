#' GEOM_ENDO
#'
#' R wrapper for Fortran binary of GEOM_ENDO (endotherm model)
#' @param AMASS mass, kg
#' @param ANDENS volume, kg/m3
#' @param FATPCT fat percentage, \%
#' @param POSTUR shape, -
#' @param ZFUR fur depth, m
#' @param SUBQFAT subcutaneous fat presence, -
#' @param SHAPE_B multiplier on long axis for shape change, -
#' @param SHAPE_B_REF reference multiplier on long axis for shape change, -
#' @param SHAPE_C multiplier on width relative to height for plate shape change, -
#' @param DHARA hair diameter, m
#' @param RHOARA hair density, 1/m2
#' @param PTCOND percentage of surface area conducting to substrate, \%
#' @param SAMODE surface area mode, -
#' @param ORIENT orientation to solar, -
#' @export
GEOM_ENDO <- function(AMASS, ANDENS, FATPCT, POSTUR, ZFUR, SUBQFAT, SHAPE_B, SHAPE_B_REF, SHAPE_C, DHARA, RHOARA, PTCOND, SAMODE, ORIENT){
  a <- .Fortran("GEOM_ENDO",
    as.double(AMASS),
    as.double(ANDENS),
    as.double(FATPCT),
    as.double(POSTUR),
    as.double(ZFUR),
    as.double(SUBQFAT),
    as.double(SHAPE_B),
    as.double(SHAPE_B_REF),
    as.double(SHAPE_C),
    as.double(DHARA),
    as.double(RHOARA),
    as.double(PTCOND),
    as.double(SAMODE),
    as.double(ORIENT),
    results=matrix(data = 0, nrow = 1, ncol = 22),
    PACKAGE = "NicheMapR")

  results <- matrix(data = 0, nrow = 1, ncol = 22)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
