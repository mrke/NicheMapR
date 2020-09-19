#' ZBRENT_ENDO
#'
#' R wrapper for Fortran binary of ZBRENT_ENDO (endotherm model)
#' @param QM1 A
#' @param QM2 A
#' @param TOL A
#' @param ZBRENT.in A
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
