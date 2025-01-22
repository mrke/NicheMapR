#' get_heights - auxiliary function for HomoTherm
#'
#' This function gets the heights at the mid-point of each body part of a human.
#' @encoding UTF-8
#' @param MASSs = c(5.32, 35.07, 3.43, 11.34), mass of each part (kg)
#' @param HEIGHT = 170, total height of person (m)
#' @param DENSITYs = rep(1050, 4), density of each part (kg/m3)
#' @param SHAPE_Bs = c(1.75, 1.87, 6.65, 6.70), ratios of long to short axes of each part (-)
#' @param INSDEPDs = c(1e-02, rep(6e-03, 3)), clothing depth, dorsal (m)
#' @param INSDEPVs = c(1e-09, rep(6e-03, 3)), clothing depth, ventral (m)
#' @usage get_heights(MASSs = c(5.32, 35.07, 3.43, 11.34), HEIGHT = 170, DENSITYs = rep(1050, 4), SHAPE_Bs = c(1.75, 1.87, 6.65, 6.70))
#' @export
get_heights <- function(MASSs = c(5.32, 35.07, 3.43, 11.34),
                        HEIGHT = 170,
                        DENSITYs = rep(1050, 4),
                        SHAPE_Bs = c(1.75, 1.87, 6.65, 6.70)
){
  SAMODE <- 0
  SHAPEs <- c(4, 1, 1, 1)
  FATPCTs <- c(5, 36, 10, 23) * 0.5 # 22.6%, typical human
  GEOM.head <- GEOM_ENDO(MASSs[1], DENSITYs[1], FATPCTs[1], SHAPEs[1], 0, 0, SHAPE_Bs[1], SHAPE_Bs[1], 0, 0, 0, 0, 0, 0)
  GEOM.trunk <- GEOM_ENDO(MASSs[2], DENSITYs[2], FATPCTs[2], SHAPEs[2], 0, 0, SHAPE_Bs[2], SHAPE_Bs[2], 0, 0, 0, 0, 0, 0)
  GEOM.arm <- GEOM_ENDO(MASSs[3], DENSITYs[3], FATPCTs[3], SHAPEs[3], 0, 0, SHAPE_Bs[3], SHAPE_Bs[3], 0, 0, 0, 0, 0, 0)
  GEOM.leg <- GEOM_ENDO(MASSs[4], DENSITYs[4], FATPCTs[4], SHAPEs[4], 0, 0, SHAPE_Bs[4], SHAPE_Bs[4], 0, 0, 0, 0, 0, 0)

  GEOM.heads <- as.data.frame(GEOM.head)
  GEOM.trunks <- as.data.frame(GEOM.trunk)
  GEOM.arms <- as.data.frame(GEOM.arm)
  GEOM.legs <- as.data.frame(GEOM.leg)
  GEOM.lab <- c("VOL", "D", "MASFAT", "VOLFAT", "ALENTH", "AWIDTH",
                "AHEIT", "ATOT", "ASIL", "ASILN", "ASILP", "GMASS", "AREASKIN",
                "FLSHVL", "FATTHK", "ASEMAJ", "BSEMIN", "CSEMIN", "CONVSK",
                "CONVAR", "R1", "R2")
  colnames(GEOM.heads) <- GEOM.lab
  colnames(GEOM.trunks) <- GEOM.lab
  colnames(GEOM.arms) <- GEOM.lab
  colnames(GEOM.legs) <- GEOM.lab
  head.height <- HEIGHT / 100 - GEOM.heads$ALENTH / 2
  trunk.height <- HEIGHT / 100 - GEOM.heads$ALENTH - GEOM.trunks$ALENTH / 2
  arm.height <- HEIGHT / 100 - GEOM.heads$ALENTH - GEOM.arms$ALENTH / 2
  leg.height <- GEOM.legs$ALENTH / 2
  heights <- c(head.height, trunk.height, arm.height, leg.height)
  return(heights)
}
