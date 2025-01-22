#' get_clo - auxiliary function for HomoTherm
#'
#' This function takes the output from the HomoTherm model and computes the
#' surface area-weighted insulation, reported back in 'clo' units.
#' @encoding UTF-8
#' @param indata the HomoTherm output list
#' @param INSDEPDs = c(1e-02, rep(6e-03, 3)), clothing depth, dorsal (m)
#' @param INSDEPVs = c(1e-09, rep(6e-03, 3)), clothing depth, ventral (m)
#' @param AREAFRACs = c(0.08291887, 0.32698460, 0.11025155, 0.18479669), fractional areas (-)
#' @usage get_clo(HomoTherm.out, INSDEPDs = c(1e-02, rep(4.1e-03, 3)), INSDEPVs = c(1e-09, rep(4.1e-03, 3)))
#' @export
get_clo <- function(indata, # the output from HomoTherm
                    INSDEPDs = c(1e-02, rep(4.1e-03, 3)),
                    INSDEPVs = c(1e-09, rep(4.1e-03, 3)),
                    AREAFRACs = c(0.08291887, 0.32698460, 0.11025155, 0.18479669)){
  # surface-area weighted mean clothing thermal conductivity
  k_clo <- (indata$head.treg[8] + indata$head.treg[9]) / 2 * AREAFRACs[1] +
    (indata$trunk.treg[8] + indata$trunk.treg[9]) / 2 * AREAFRACs[2] +
    (indata$arm.treg[8] + indata$arm.treg[9]) / 2 * AREAFRACs[3] * 2 +
    (indata$leg.treg[8] + indata$leg.treg[9]) / 2 * AREAFRACs[4] * 2 # W/m/K

  # surface-area weighted mean clothing depth
  INSDEP <- mean(c(INSDEPDs[1], INSDEPVs[1]) * AREAFRACs[1] +
                   mean(c(INSDEPDs[2], INSDEPVs[2])) * AREAFRACs[2] +
                   mean(c(INSDEPDs[3], INSDEPVs[3])) * AREAFRACs[3] * 2 +
                   mean(c(INSDEPDs[4], INSDEPVs[4])) * AREAFRACs[4] * 2) # m

  # calculate implied CLO
  clo <- INSDEP / (k_clo * 0.155)
  names(clo) <- 'clo'
  clo
}
