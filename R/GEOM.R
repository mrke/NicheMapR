#' GEOM
#'
#' R wrapper for Fortran binary of GEOM (endotherm model)
#' @param AMASS mass, kg
#' @param ANDENS volume, kg/m3
#' @param FATPCT fat percentage, \%
#' @param POSTUR shape, -
#' @param ZFUR fur depth, m
#' @param SUBQFAT subcutaneous fat presence, -
#' @param GMULT multiplier on long axis for shape change, -
#' @param GMREF reference multiplier on long axis for shape change, -
#' @param DHARA hair diameter, m
#' @param RHOARA hair density, 1/m2
#' @param PTCOND percentage of surface area conducting to substrate, \%
#' @param SAMODE surface area mode, -
#' @param ORIENT orientation to solar, -
#' @export
GEOM <- function(AMASS, ANDENS, FATPCT, POSTUR, ZFUR, SUBQFAT, GMULT, GMREF, DHARA, RHOARA, PTCOND, SAMODE, ORIENT){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/GEOM.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/GEOM.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/GEOM.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/GEOM.so'
  }
  if (!is.loaded('GEOM')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("GEOM",
    as.double(AMASS),
    as.double(ANDENS),
    as.double(FATPCT),
    as.double(POSTUR),
    as.double(ZFUR),
    as.double(SUBQFAT),
    as.double(GMULT),
    as.double(GMREF),
    as.double(DHARA),
    as.double(RHOARA),
    as.double(PTCOND),
    as.double(SAMODE),
    as.double(ORIENT),
    results=matrix(data = 0, nrow = 1, ncol = 22),
    PACKAGE = "GEOM")
  #dyn.unload("GEOM.dll")

  results <- matrix(data = 0, nrow = 1, ncol = 22)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
