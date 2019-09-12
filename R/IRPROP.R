#' IRPROP
#'
#' R wrapper for Fortran binary of IRPROP (endotherm model)
#' @param TA A
#' @param GMULTMAX A
#' @param GMREF A
#' @param GMULT A
#' @param DHAIRD A
#' @param DHAIRV A
#' @param LHAIRD A
#' @param LHAIRV A
#' @param ZFURD A
#' @param ZFURV A
#' @param RHOD A
#' @param RHOV A
#' @param REFLD A
#' @param REFLV A
#' @param MAXPTVEN A
#' @param ZFURCOMP A
#' @export
IRPROP <- function(TA, GMULTMAX, GMREF, GMULT, DHAIRD, DHAIRV, LHAIRD, LHAIRV,
  ZFURD, ZFURV, RHOD, RHOV, REFLD, REFLV, MAXPTVEN, ZFURCOMP){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/IRPROP.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/IRPROP.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/IRPROP.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/IRPROP.so'
  }
  if (!is.loaded('IRPROP')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("IRPROP",
                as.double(TA),
                as.double(GMULTMAX),
                as.double(GMREF),
                as.double(GMULT),
                as.double(DHAIRD),
                as.double(DHAIRV),
                as.double(LHAIRD),
                as.double(LHAIRV),
                as.double(ZFURD),
                as.double(ZFURV),
                as.double(RHOD),
                as.double(RHOV),
                as.double(REFLD),
                as.double(REFLV),
                as.double(MAXPTVEN),
                as.double(ZFURCOMP),
                results=matrix(data = 0., nrow = 1, ncol = 27),
                PACKAGE = "IRPROP")
  #dyn.unload("IRPROP.dll")

  results <- matrix(data = 0., nrow = 1, ncol = 27)

  storage.mode(results)<-"double"
  results <- a$results

  return (results)
}
