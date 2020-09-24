#' SOLAR_ENDO
#'
#' R wrapper for Fortran binary of SOLAR_ENDO (endotherm model)
#' @encoding UTF-8
#' @param AREATOTL surface area for solar exchange, m2
#' @param ABSAND solar absorptivity of dorsal fur (fractional, 0-1)
#' @param ABSANV solar absorptivity of ventral fur (fractional, 0-1)
#' @param ABSSB  solar absorptivity of substrate (fractional, 0-1)
#' @param ASIL silhouette area to sun, m2
#' @param PCTDIF proportion of solar radiation that is diffuse (fractional, 0-1)
#' @param QNORM compute solar radiation on a surface normal to the direct rays of the sun (W/m2)
#' @param SHADE shade (fractional, 0-1)
#' @param QSOLR solar radiation, horizontal plane (W/m2)
#' @param FASKY configuration factor to sky (-)
#' @param FAVEG configuration factor to vegetation (-)
#' @export
SOLAR_ENDO <- function(AREATOTL, ABSAND, ABSANV, ABSSB, ASIL, PCTDIF, QNORM, SHADE,
  QSOLR, FASKY, FAVEG){
  a <- .Fortran("SOLAR_ENDO",
    as.double(AREATOTL),
    as.double(ABSAND),
    as.double(ABSANV),
    as.double(ABSSB),
    as.double(ASIL),
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
