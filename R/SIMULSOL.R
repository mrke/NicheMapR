#' SIMULSOL
#'
#' R wrapper for Fortran binary of SIMULSOL (endotherm model)
#' @encoding UTF-8
#' @param DIFTOL error tolerance for SIMULSOL
#' @param IPT geometry for SIMULSOL, 1 = cylinder, 2 = sphere, 3 = ellipsoid
#' @param FURVARS fur input variable vector
#' @param GEOMVARS shape and size input variable vector
#' @param ENVVARS environmental input vector
#' @param TRAITS other trait inputs vector
#' @param TFA current guess of fur/air-interface temperature (°C)
#' @param SKINW part of the skin surface that is wet (\%)
#' @param TSKIN current guess of skin temperature (°C)
#' @export
SIMULSOL <- function(DIFTOL, IPT, FURVARS, GEOMVARS, ENVVARS, TRAITS, TFA,
                     SKINW, TSKIN, results){
  a <- .Fortran("SIMULSOL",
    as.double(DIFTOL),
    as.double(IPT),
    as.double(FURVARS),
    as.double(GEOMVARS),
    as.double(ENVVARS),
    as.double(TRAITS),
    as.double(TFA),
    as.double(SKINW),
    as.double(TSKIN),
    results=matrix(data = 0., nrow = 1, ncol = 16),
    PACKAGE = "NicheMapR")

  results <- matrix(data = 0., nrow = 1, ncol = 16)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
