#' CONV_ENDO
#'
#' R wrapper for Fortran binary of CONV_ENDO (endotherm model)
#' @param TS A
#' @param TENV A
#' @param NGEOM A
#' @param SURFAR A
#' @param FLTYPE A
#' @param FURTST A
#' @param D A
#' @param TFA A
#' @param VEL A
#' @param ZFUR A
#' @param BP A
#' @param ELEV A
#' @export
CONV_ENDO <- function(TS, TENV, NGEOM, SURFAR, FLTYPE, FURTST, D, TFA, VEL, ZFUR, BP, ELEV){
  a <- .Fortran("CONV_ENDO",
    as.double(TS),
    as.double(TENV),
    as.double(NGEOM),
    as.double(SURFAR),
    as.double(FLTYPE),
    as.double(FURTST),
    as.double(D),
    as.double(TFA),
    as.double(VEL),
    as.double(ZFUR),
    as.double(BP),
    as.double(ELEV),
    results=matrix(data = 0, nrow = 1, ncol = 14),
    PACKAGE = "NicheMapR")

  results <- matrix(data = 0, nrow = 1, ncol = 14)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
