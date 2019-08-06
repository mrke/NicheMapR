#' SOLAR
#'
#' R wrapper for Fortran binary of SOLAR (endotherm model)
#' @param AREATOTL A
#' @param ABSAND A
#' @param ABSANV A
#' @param ABSSB A
#' @param ASILN A
#' @param PCTDIF A
#' @param QNORM A
#' @param SHADE A
#' @param QSOLR A
#' @param FASKY A
#' @param FATOBJ A
#' @param FAVEG A
#' @export
SOLAR <- function(AREATOTL, ABSAND, ABSANV, ABSSB, ASILN, PCTDIF, QNORM, SHADE,
  QSOLR, FASKY, FATOBJ, FAVEG){
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/win/x64/SOLAR.dll'
    } else {
      libpath='/NicheMapR/libs/win/i386/SOLAR.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/linux/SOLAR.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/mac/SOLAR.so'
  }
  if (!is.loaded('SOLAR')) {
    dyn.load(paste0(lib.loc = .libPaths()[1],libpath))
  }
  a <- .Fortran("SOLAR",
    as.double(AREATOTL),
    as.double(ABSAND),
    as.double(ABSANV),
    as.double(ABSSB),
    as.double(ASILN),
    as.double(PCTDIF),
    as.double(QNORM),
    as.double(SHADE),
    as.double(QSOLR),
    as.double(FASKY),
    as.double(FATOBJ),
    as.double(FAVEG),
    results=matrix(data = 0., nrow = 1, ncol = 8),
    PACKAGE = "SOLAR")
  #dyn.unload("SOLAR.dll")

  results <- matrix(data = 0., nrow = 1, ncol = 8)

  storage.mode(results)<-"double"
  results <- a$results
  return (results)
}
