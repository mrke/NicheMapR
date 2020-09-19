#' IRPROP
#'
#' R wrapper for Fortran binary of IRPROP (endotherm model)
#' @param TA A
#' @param SHAPE_B_MAX A
#' @param SHAP_B_REF A
#' @param SHAPE_B A
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
#' @param KHAIR A
#' @export
IRPROP <- function(TA, SHAPE_B_MAX, SHAP_B_REF, SHAPE_B, DHAIRD, DHAIRV, LHAIRD, LHAIRV,
  ZFURD, ZFURV, RHOD, RHOV, REFLD, REFLV, MAXPTVEN, ZFURCOMP, KHAIR){
  a <- .Fortran("IRPROP",
                as.double(TA),
                as.double(SHAPE_B_MAX),
                as.double(SHAP_B_REF),
                as.double(SHAPE_B),
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
                as.double(KHAIR),
                results=matrix(data = 0., nrow = 1, ncol = 26),
                PACKAGE = "NicheMapR")

  results <- matrix(data = 0., nrow = 1, ncol = 26)

  storage.mode(results)<-"double"
  results <- a$results

  return (results)
}
