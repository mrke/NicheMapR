#' ZBRENT
#'
#' R wrapper for Fortran binary of ZBRENT (endotherm model)
#' @param QM1 A
#' @param QM2 A
#' @param TOL A
#' @param ZBRENT.in A
#' @export
ZBRENT <- function(QM1, QM2, TOL, ZBRENT.in, results){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/ZBRENT.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/ZBRENT.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/ZBRENT.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/ZBRENT.so'
  }
  if (!is.loaded('ZBRENT')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("ZBRENT",
    as.double(QM1),
    as.double(QM2),
    as.double(TOL),
    as.double(ZBRENT.in),
    results=matrix(data = 0., nrow = 1, ncol = 15),
    PACKAGE = "ZBRENT")
  #dyn.unload("ZBRENT.dll")
  results <- matrix(data = 0., nrow = 1, ncol = 15)
  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
