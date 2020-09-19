#' SEVAP_ENDO
#'
#' R wrapper for Fortran binary of SEVAP_ENDO (endotherm model)
#' @param BP A
#' @param TA A
#' @param RELHUM A
#' @param VEL A
#' @param TC A
#' @param TSURF A
#' @param ALT A
#' @param SKINW A
#' @param FLYHR A
#' @param CONVSK A
#' @param HD A
#' @param HDFREE A
#' @param PCTBAREVAP A
#' @param PCTEYES A
#' @param ZFUR A
#' @param FURWET A
#' @param TFA A
#' @param CONVAR A
#' @export
SEVAP_ENDO <- function(BP, TA, RELHUM, VEL, TC, TSURF, ALT, SKINW, FLYHR,
  CONVSK, HD, HDFREE, PCTBAREVAP, PCTEYES, ZFUR, FURWET, TFA, CONVAR){
  a <- .Fortran("SEVAP_ENDO",
    as.double(BP),
    as.double(TA),
    as.double(RELHUM),
    as.double(VEL),
    as.double(TC),
    as.double(TSURF),
    as.double(ALT),
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
