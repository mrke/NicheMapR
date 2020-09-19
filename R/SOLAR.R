#' SOLAR_ENDO
#'
#' R wrapper for Fortran binary of SOLAR_ENDO (endotherm model)
#' @param AREATOTL A
#' @param ABSAND A
#' @param ABSANV A
#' @param ABSSB A
#' @param ASILN A
#' @param PCTDIF A
#' @param QNORM A
#' @param SHADE A
#' @param QSOLR A
#' @param FASKY A
#' @param FATOBJ A
#' @param FAVEG A
#' @export
SOLAR_ENDO <- function(AREATOTL, ABSAND, ABSANV, ABSSB, ASILN, PCTDIF, QNORM, SHADE,
  QSOLR, FASKY, FAVEG){
  a <- .Fortran("SOLAR_ENDO",
    as.double(AREATOTL),
    as.double(ABSAND),
    as.double(ABSANV),
    as.double(ABSSB),
    as.double(ASILN),
    as.double(PCTDIF),
    as.double(QNORM),
    as.double(SHADE),
    as.double(QSOLR),
    as.double(FASKY),
    as.double(FAVEG),
    results=matrix(data = 0., nrow = 1, ncol = 7),
    PACKAGE = "NicheMapR")

  results <- matrix(data = 0., nrow = 1, ncol = 7)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
