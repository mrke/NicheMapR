#' SOLVENDO
#'
#' R wrapper for Fortran binary of SOLVENDO (endotherm model)
#' @param AMASS A
#' @export
SOLVENDO <- function(SOLVENDO.input){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/SOLVENDO.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/SOLVENDO.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/SOLVENDO.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/SOLVENDO.so'
  }
  if (!is.loaded('SOLVENDO')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("SOLVENDO",
                as.double(SOLVENDO.input),
                 results=matrix(data = 0., nrow = 1, ncol = 58),
                PACKAGE = "SOLVENDO")
  #dyn.unload("SOLVENDO.dll")

  results <- matrix(data = 0., nrow = 1, ncol = 58)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
