#' GEOM
#'
#' R wrapper for Fortran binary of GEOM (endotherm model)
#' @param AMASS A
#' @param ANDENS A
#' @param FATPCT A
#' @param POSTUR A
#' @param ZFUR A
#' @param SUBQFAT A
#' @param GMULT A
#' @param GMREF A
#' @param DHARA A
#' @param RHOARA A
#' @param PTCOND A
#' @param BIRD A
#' @param MAMMAL A
#' @param ORIENT A
#' @export
GEOM <- function(AMASS, ANDENS, FATPCT, POSTUR, ZFUR, SUBQFAT, GMULT,
  GMREF, DHARA, RHOARA, PTCOND, BIRD, MAMMAL, ORIENT){
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
    as.double(BIRD),
    as.double(MAMMAL),
    as.double(ORIENT),
    results=matrix(data = 0, nrow = 1, ncol = 25),
    PACKAGE = "GEOM")
  #dyn.unload("GEOM.dll")

  results <- matrix(data = 0, nrow = 1, ncol = 25)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
