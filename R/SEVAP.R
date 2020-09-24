#' SEVAP_ENDO
#'
#' R wrapper for Fortran binary of SEVAP_ENDO (endotherm model)
#' @encoding UTF-8
#' @param BP barometric pressure (Pa), negative means altitude is used
#' @param TA air temperature (°C)
#' @param RELHUM relative humidity (\%)
#' @param VEL wind speed (m/s)
#' @param TC core temperature (°C)
#' @param TSURF air temperature (°C)
#' @param ELEV elevation (m)
#' @param SKINW part of the skin surface that is wet (\%)
#' @param FLYHR is flight occurring this hour? (imposes forced evaporative loss)
#' @param CONVSK area of skin for evaporation (total skin area - hair area), m2
#' @param HD mass transfer coefficient
#' @param HDFREE free mass transfer coefficient
#' @param PCTBAREVAP surface area for evaporation that is skin, e.g. licking paws (\%)
#' @param PCTEYES surface area made up by the eye (\%) - make zero if sleeping
#' @param ZFUR fur depth (m)
#' @param FURWET part of the fur surface that is wet (\%)
#' @param TFA fur/air interface temperature (°C)
#' @param CONVAR area for convection (total area minus ventral area, as determined by PCOND), m2
#' @export
SEVAP_ENDO <- function(BP, TA, RELHUM, VEL, TC, TSURF, ELEV, SKINW, FLYHR,
  CONVSK, HD, HDFREE, PCTBAREVAP, PCTEYES, ZFUR, FURWET, TFA, CONVAR){
  a <- .Fortran("SEVAP_ENDO",
    as.double(BP),
    as.double(TA),
    as.double(RELHUM),
    as.double(VEL),
    as.double(TC),
    as.double(TSURF),
    as.double(ELEV),
    as.double(SKINW),
    as.double(FLYHR),
    as.double(CONVSK),
    as.double(HD),
    as.double(HDFREE),
    as.double(PCTBAREVAP),
    as.double(PCTEYES),
    as.double(ZFUR),
    as.double(FURWET),
    as.double(TFA),
    as.double(CONVAR),
    results=matrix(data = 0., nrow = 1, ncol = 7),
    PACKAGE = "NicheMapR")

  results <- matrix(data = 0., nrow = 1, ncol = 7)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
