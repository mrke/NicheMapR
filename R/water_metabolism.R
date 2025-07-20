#' Estimate Metabolic Energy and Water Flux from Diet Composition
#'
#' Calculates energy availability and water gain/loss from digestion and excretion
#' based on dietary composition and metabolic rate.
#'
#' @param p_prot Proportion of protein in dry food (default = 0.44).
#' @param p_fat Proportion of fat in dry food (default = 0.48).
#' @param p_carb Proportion of carbohydrate in dry food (default = 0.08 * 0.8).
#' @param p_dry Proportion of dry matter in wet food (default = 0.3).
#' @param kap_X Digestive efficiency (default = 0.8).
#' @param mrate Metabolic rate in J/time.
#' @param p_X Feeding rate in J/time (default = `mrate * 2`).
#' @param p_urea Proportion of urea in nitrogenous waste (default = 0.1).
#' @param p_fecalH2O Proportion of water in faeces (default = 0.8).
#'
#' @return A list containing:
#' \describe{
#'   \item{Q_avail}{Net energy available after metabolism (J/time).}
#'   \item{H2O_avail}{Net water available after accounting for excretion (g/time).}
#'   \item{H2OMet_g}{Metabolic water generated (g/time).}
#'   \item{H2OFaeces_g}{Water lost via faeces (g/time).}
#'   \item{WetFaeces_G}{Total wet mass of faeces (g/time).}
#'   \item{H2OUrine_g}{Water lost via urine (g/time).}
#'   \item{H2OFree_g}{Free water absorbed from wet food (g/time).}
#'   \item{WetFood_g}{Total wet food consumed (g/time).}
#'   \item{DryFood_g}{Total dry food consumed (g/time).}
#' }
#'
#' @details
#' Macronutrient metabolism produces metabolic water and contributes to energy
#' gain. This function estimates energy available after subtracting metabolic
#' requirements and calculates water gained from food and lost through excretion.
#'
#' @examples
#' water_metabolism(mrate = 5000)
#'
#' @export
water_metabolism <- function(p_prot = 0.44,
                             p_fat = 0.48,
                             p_carb = 0.08 * 0.8,
                             p_dry = 0.3,
                             kap_X = 0.8,
                             mrate = mrate,
                             p_X = mrate * 2,
                             p_urea = 0.1,
                             p_fecalH2O = 0.8) {

  gfatpg <- p_fat
  gprotpg <- p_prot
  gcarbpg <- p_carb
  gsum <- gfatpg + gprotpg + gcarbpg
  gundig <- 1.00 - gsum
  totcarb <- gcarbpg + gundig

  fatjpg <- gfatpg * (9400 * 4.185)
  protjpg <- gprotpg * (4199 * 4.185)
  carbjpg <- totcarb * (4200 * 4.185)

  jabspgr <- fatjpg + protjpg + kap_X * carbjpg
  dryabs <- p_X / jabspgr
  DryFood_g <- dryabs / (1 - gundig)
  WetFood_g <- DryFood_g / p_dry
  Q_avail <- jabspgr * dryabs - mrate

  gprot <- dryabs * gprotpg
  gfat  <- dryabs * gfatpg
  gcarb <- dryabs * gcarbpg
  gtot <- gprot + gfat + gcarb

  H2OFree_g <- gtot / p_dry - gtot
  H2OMet_g <- gprot * 0.40 + gfat * 1.07 + gcarb * 0.56

  gurea <- gprot * 0.343
  gurine <- gurea / p_urea
  H2OUrine_g <- (gurine - gurea) / 1.0474

  DryFaeces_g <- DryFood_g * (1 - kap_X)
  WetFaeces_G <- DryFaeces_g / (1 - p_fecalH2O)
  H2OFaeces_g <- WetFaeces_G * p_fecalH2O - (WetFood_g - gtot / p_dry) * (1 - p_dry)

  H2O_avail <- H2OFree_g + H2OMet_g - H2OUrine_g - H2OFaeces_g

  return(list(Q_avail = Q_avail,
              H2O_avail = H2O_avail,
              H2OMet_g = H2OMet_g,
              H2OFaeces_g = H2OFaeces_g,
              WetFaeces_G = WetFaeces_G,
              H2OUrine_g = H2OUrine_g,
              H2OFree_g = H2OFree_g,
              WetFood_g = WetFood_g,
              DryFood_g = DryFood_g))
}
