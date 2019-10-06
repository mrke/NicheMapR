#' CONV
#'
#' R wrapper for Fortran binary of CONV (endotherm model)
#' @param TS A
#' @param TENV A
#' @param NGEOM A
#' @param SURFAR A
#' @param FLTYPE A
#' @param FURTST A
#' @param D A
#' @param TFA A
#' @param VEL A
#' @param ZFUR A
#' @param BP A
#' @param ELEV A
#' @export
CONV <- function(TS, TENV, NGEOM, SURFAR, FLTYPE, FURTST, D, TFA, VEL, ZFUR, BP, ELEV){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/CONV.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/CONV.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/CONV.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/CONV.so'
  }
  if (!is.loaded('CONV')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("CONV",
    as.double(TS),
    as.double(TENV),
    as.double(NGEOM),
    as.double(SURFAR),
    as.double(FLTYPE),
    as.double(FURTST),
    as.double(D),
    as.double(TFA),
    as.double(VEL),
    as.double(ZFUR),
    as.double(BP),
    as.double(ELEV),
    results=matrix(data = 0, nrow = 1, ncol = 14),
    PACKAGE = "CONV")
  #dyn.unload("CONV.dll")

  results <- matrix(data = 0, nrow = 1, ncol = 14)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
