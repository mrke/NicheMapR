#' F_FACTOR
#'
#' R wrapper for Fortran binary of F_FACTOR (endotherm model)
#' @param ASHADE A
#' @param NITESHAD A
#' @param QSOLR A
#' @param FATOB A
#' @param NESTYP A
#' @param RoNEST A
#' @param R1 A
#' @param FGDREF A
#' @param FSKREF A
#' @param AREASKIN A
#' @param EMISAN A
#' @export
F_FACTOR <- function(ASHADE, NITESHAD, QSOLR, FATOB, NESTYP, RoNEST, R1, FGDREF,
  FSKREF, AREASKIN, EMISAN){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/F_FACTOR.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/F_FACTOR.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/F_FACTOR.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/F_FACTOR.so'
  }
  if (!is.loaded('F_FACTOR')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("F_FACTOR",
    as.double(ASHADE),
    as.double(NITESHAD),
    as.double(QSOLR),
    as.double(FATOB),
    as.double(NESTYP),
    as.double(RoNEST),
    as.double(R1),
    as.double(FGDREF),
    as.double(FSKREF),
    as.double(AREASKIN),
    as.double(EMISAN),
    results=matrix(data = 0., nrow = 1, ncol = 9),
    PACKAGE = "F_FACTOR")
  #dyn.unload("F_FACTOR.dll")

  results <- matrix(data = 0., nrow = 1, ncol = 9)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
