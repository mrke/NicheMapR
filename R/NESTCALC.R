#' NESTCALC
#'
#' R wrapper for Fortran binary of NESTCALC (endotherm model)
#' @encoding UTF-8
#' @param NESTCALC.input description
NESTCALC <- function(NESTCALC.input){
  a <- .Fortran("NESTCALC",
    as.double(NESTCALC.input),
    results=matrix(data = 0, nrow = 1, ncol = 2),
    PACKAGE = "NicheMapR")
  results <- matrix(data = 0, nrow = 1, ncol = 2)
  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
