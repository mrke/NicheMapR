#' SIMULSOL
#'
#' R wrapper for Fortran binary of SIMULSOL (endotherm model)
#' @encoding UTF-8
#' @param DIFTOL A
#' @param IPT A
#' @param FURVARS A
#' @param GEOMVARS A
#' @param ENVVARS A
#' @param TRAITS A
#' @param TFA A
#' @param SKINW A
#' @param TSKIN A
#' @export
SIMULSOL <- function(DIFTOL, IPT, FURVARS, GEOMVARS, ENVVARS, TRAITS, TFA,
                     SKINW, TSKIN, results){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/SIMULSOL.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/SIMULSOL.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/SIMULSOL.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/SIMULSOL.so'
  }
  if (!is.loaded('SIMULSOL')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
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
    results=matrix(data = 0., nrow = 1, ncol = 15),
    PACKAGE = "SIMULSOL")
  #dyn.unload("SIMULSOL.dll")

  results <- matrix(data = 0., nrow = 1, ncol = 15)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
