#' ZBRAC_SHELTER
#'
#' R wrapper for Fortran binary of ZBRAC_SHELTER (endotherm model)
#' @param X1 A
#' @param X2 A
#' @param TOL A
#' @param ZBRENT.in A
#' @export
ZBRAC_SHELTER <- function(X1, X2, ZBRENT.in, results){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/ZBRAC_SHELTER.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/ZBRAC_SHELTER.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/ZBRAC_SHELTER.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/ZBRAC_SHELTER.so'
  }
  if (!is.loaded('ZBRAC_SHELTER')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("ZBRAC_SHELTER",
    as.double(X1),
    as.double(X2),
    as.double(ZBRENT.in),
    results=matrix(data = 0, nrow = 1, ncol = 2),
    PACKAGE = "ZBRAC_SHELTER")
  #dyn.unload("ZBRAC_SHELTER.dll")
  results <- matrix(data = 0, nrow = 1, ncol = 2)
  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
