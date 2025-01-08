#' GET_SHAPES
#'
#' R wrapper for Fortran binary of GET_SHAPES (HomoTherm model)
#' @encoding UTF-8
#' @param MASSs mass per body part (head, trunk, arm, leg), kg
#' @param HEIGHT target height, m
#' @param AREA body surface area, m2
#' @param SHAPE_Bs shape factors per body part, -
#' @param SHAPE_Bs.min min allowable SHAPE_Bs, -
#' @param SHAPE_Bs.max max allowable SHAPE_Bs, -
#' @param SUBQFAT subcutaneous fat presence, -
#' @param AREAFRACs area fraction per body part, -
#' @param DENSITYs density per body part, kg/m3
#' @param FATPCTs fat percentage per body part, /%
#' @param SHAPEs shape type per body part, -
#' @param SUBQFATs presence of subcutaneous fat per body part, -
#' @param DHAIRDs dorsal hair diameters, m
#' @param DHAIRVs ventral hair diameters, m
#' @param INSDENDs dorsal insulation depth, m
#' @param INSDENvs ventral insulation depth, m
#' @param tol tolerance, -
#' @param maxiter maximum iterations,
#' @export
GET_SHAPES <- function(MASSs = c(5.0987, 33.5670,  3.2830, 10.8540),
                       HEIGHT = 170,
                       AREA = 1.947321,
                       SHAPE_Bs = c(1.75, 1.87, 6.65, 6.7),
                       SHAPE_Bs.min = c(1.3, 2.3, 7.5, 6.5),
                       SHAPE_Bs.max = c(1.7, 4, 21, 20),
                       AREAFRACs = c(0.08291887, 0.32698460, 0.11025155, 0.18479669),
                       DENSITYs = c(1050, 1050, 1050, 1050),
                       FATPCTs = c(5, 36, 10, 23),
                       SHAPEs = c(4, 1, 1, 1),
                       SUBQFATs = c(1, 1, 1, 1),
                       DHAIRDs = c(3e-05, 3e-05, 3e-05, 3e-05),
                       DHAIRVs = c(3e-05, 3e-05, 3e-05, 3e-05),
                       INSDENDs = c(3e+08, 3e+08, 3e+08, 3e+08),
                       INSDENVs = c(3e+08, 3e+08, 3e+08, 3e+08),
                       tol = 1e-6,
                       maxiter = 10000){
  a <- .Fortran("GET_SHAPES",
                as.double(MASSs),
                as.double(HEIGHT),
                as.double(AREA),
                as.double(SHAPE_Bs),
                as.double(SHAPE_Bs.min),
                as.double(SHAPE_Bs.max),
                as.double(AREAFRACs),
                as.double(DENSITYs),
                as.double(FATPCTs),
                as.double(SHAPEs),
                as.double(SUBQFATs),
                as.double(DHAIRDs),
                as.double(DHAIRVs),
                as.double(INSDENDs),
                as.double(INSDENVs),
                as.double(tol),
                as.double(maxiter),
                results=matrix(data = 0, nrow = 1, ncol = 10),
                PACKAGE = "NicheMapR")

  results <- matrix(data = 0, nrow = 1, ncol = 10)

  storage.mode(results)<-"double"
  results <- a$results
  SHAPE_Bs <- results[1:4]
  PJOINs <- results[5:8]
  HEIGHT_out <- results[9]
  AREA_out <- results[10]
  return (list(SHAPE_Bs=SHAPE_Bs, PJOINs=PJOINs, HEIGHT_out = HEIGHT_out, AREA_out = AREA_out))
}
