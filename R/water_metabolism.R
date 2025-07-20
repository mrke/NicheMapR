#' Estimate Metabolic Energy and Water Flux from Diet Composition
#'
#' Calculates energy availability and water gain/loss from digestion and excretion
#' based on dietary composition and metabolic rate.
#'
#' @param mrate Metabolic rate in J/time.
#' @param p_X Feeding rate in J/time (default = `mrate * 2`).
#' @param pct_prot Percentage of protein in dry food (default = 44).
#' @param pct_fat Percentage of fat in dry food (default = 48).
#' @param pct_carb Percentage of carbohydrate in dry food (default = 0.08 * 80).
#' @param kap_X Digestive efficiency (default = 0.8).
#' @param pct_H_X Percentage of water in food (default = 82).
#' @param pct_H_P Percentage of water in faeces (default = 73).
#' @param pct_H_N Percentage of water excreted in nitrogenous waste (default = 0, i.e. all uric acid).
#' @param pct_urea Percentage of urea in nitrogenous waste (default = 10).
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
#' Assumptions of energy content and metabolic water based on the following
#' (taken from original ectotherm model code)
#'
#' PROTEINS
#' Approximate gram molecular weight of amino acids is 137 g/mole
#' Based on info from Handbook of Chemistry & Physics
#' Alanine = 89.1 g/mole (with water - 18 g/mole)
#' Arginine = 174.2 "          "
#' Aspargine =132.1 "          "
#' Aspartic acid = 133.1       "
#' Cystine = 121.2            etc.
#' Glutamic acid = 147.1
#' Glutamine = 146.1
#' Glycine = 75.1
#' Histidine = 155.2
#' Isoleucine = 131.2
#' Leucine = 131.2
#' Lysine = 146.2
#' Methionine = 149.2
#' Phenylalanine = 165.2
#' Proline = 115.1
#' Serine = 105.1
#' Threonine = 119.1
#' Tryptophan = 204.2
#' Tyrosine = 181.2
#' Valine = 117.1        Average = 137 g/mole
#' 4300 calories/g, calories/ml O2 = 4.5 Hainsworth, F.R. 1981
#' 0.40 g WATER/g protein oxidized: Hainsworth, F.R. 1981.
#' Animal Physiology.Addison-Wesley Publ. Co., Reading, MA. 669 p.
#'
#' LIPIDS
#' Data from Guyton, A.C. 1991.  Textbook of Medical Physiology.
#' 8th ed.  W.B. Saunders, Philadelphia.  1014 p.
#' Triglycerides used for energy.  They are stearic acid (880 g/mol),
#' oleic acid (879 g/mol) and palmitic acid (754 g/mol).
#' Assume as an average 850 g/mol.
#' calories/ml O2 = 4.7; 9400 calories/g fat  (Kleiber, 1961)
#' 1.07 g WATER/g lipid oxidized: Hainsworth, F.R. 1981.
#'
#' CARBOHYDRATES
#' Glucose (180 g/mol); 5.0 cal/ml O2; 4200 calories/g  (Kleiber, 1961)
#' 0.56 g WATER/g carbohydrate oxidized: Hainsworth, F.R. 1981.
#' @examples
#' water_metabolism(mrate = 5000)
#'
#' @export
water_metabolism <- function(mrate = mrate,
                             p_X = mrate * 2,
                             pct_prot = 0.44,
                             pct_fat = 0.48,
                             pct_carb = 0.08 * 0.8,
                             kap_X = 0.8,
                             pct_H_X = 82,
                             pct_H_P = 73,
                             pct_H_N = 0,
                             pct_urea = 10) {

  gfatpg <- pct_fat / 100
  gprotpg <- pct_prot / 100
  gcarbpg <- pct_carb / 100
  p_dry <- 1 - (pct_H_X / 100)
  p_fecalH2O <- pct_H_P / 100
  p_urea <- pct_urea / 100
  p_H_N <- pct_H_N / 100

  gsum <- gfatpg + gprotpg + gcarbpg
  gundig <- 1.00 - gsum
  totcarb <- gcarbpg + gundig

  fatjpg <- gfatpg * (9400 * 4.185)
  protjpg <- gprotpg * (4300 * 4.185)
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
  if(p_urea > 0){
    gurine <- gurea / p_urea
    H2OUrine_g <- (gurine - gurea) / 1.0474
  }else{
    H2OUrine_g <- 0
  }

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
