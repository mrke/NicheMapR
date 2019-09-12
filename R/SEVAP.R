#' SEVAP
#'
#' R wrapper for Fortran binary of SEVAP (endotherm model)
#' @param BP A
#' @param TA A
#' @param RELHUM A
#' @param VEL A
#' @param TC A
#' @param TSURF A
#' @param ALT A
#' @param SKINW A
#' @param FLYHR A
#' @param CONVSK A
#' @param HD A
#' @param HDFREE A
#' @param PCTBAREVAP A
#' @param PCTEYES A
#' @param ZFUR A
#' @param FURWET A
#' @param TFA A
#' @param CONVAR A
#' @export
SEVAP <- function(BP, TA, RELHUM, VEL, TC, TSURF, ALT, SKINW, FLYHR,
  CONVSK, HD, HDFREE, PCTBAREVAP, PCTEYES, ZFUR, FURWET, TFA, CONVAR){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/SEVAP.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/SEVAP.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/SEVAP.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/SEVAP.so'
  }
  if (!is.loaded('SEVAP')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("SEVAP",
    as.double(BP),
    as.double(TA),
    as.double(RELHUM),
    as.double(VEL),
    as.double(TC),
    as.double(TSURF),
    as.double(ALT),
    as.double(SKINW),
    as.double(FLYHR),
    as.double(CONVSK),
    as.double(HD),
    as.double(HDFREE),
    as.double(PCTBAREVAP),
    as.double(PCTEYES),
    as.double(ZFUR),
    as.double(FURWET),
    as.double(TFA),
    as.double(CONVAR),
    results=matrix(data = 0., nrow = 1, ncol = 7),
    PACKAGE = "SEVAP")
  #dyn.unload("SEVAP.dll")

  results <- matrix(data = 0., nrow = 1, ncol = 7)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
