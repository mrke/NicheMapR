#' ZBRENT_SHELTER
#'
#' R wrapper for Fortran binary of ZBRENT_SHELTER (endotherm model)
#' @param X1 A
#' @param X2 A
#' @param TOL A
#' @param ZBRENT.in A
#' @export
ZBRENT_SHELTER <- function(X1, X2, TOL, ZBRENT.in){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/ZBRENT_SHELTER.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/ZBRENT_SHELTER.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/ZBRENT_SHELTER.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/ZBRENT_SHELTER.so'
  }
  if (!is.loaded('ZBRENT_SHELTER')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("ZBRENT_SHELTER",
                as.double(X1),
                as.double(X2),
                as.double(TOL),
                as.double(ZBRENT.in),
                results=matrix(data = 0, nrow = 1, ncol = 3),
                PACKAGE = "ZBRENT_SHELTER")
  #dyn.unload("ZBRENT_SHELTER.dll")
  results <- a[[1]][1]
}
