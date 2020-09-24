#' ZBRENT_ENDO
#'
#' R wrapper for Fortran binary of ZBRENT_ENDO (endotherm model)
#' @encoding UTF-8
#' @param QM1 bracketing guess 1 of metabolic heat generation, W
#' @param QM2 bracketing guess 2 of metabolic heat generation, W
#' @param TOL tolerance of solution, W
#' @param ZBRENT.in vector of input parameters for RESPFUN
#' @export
ZBRENT_ENDO <- function(QM1, QM2, TOL, ZBRENT.in, results){

  a <- .Fortran("ZBRENT_ENDO",
    as.double(QM1),
    as.double(QM2),
    as.double(TOL),
    as.double(ZBRENT.in),
    results=matrix(data = 0., nrow = 1, ncol = 15),
    PACKAGE = "NicheMapR")
  results <- matrix(data = 0., nrow = 1, ncol = 15)
  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
