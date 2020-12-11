#' Pedotransfer functions for soil hydraulic properties
#'
#' A function to compute soil hydrological properties from information on bulk
#' density and soil composition (clay/silt/sand composition) at particular depths,
#' with capacity to spline results to other depths.
#' Calculations of Campbell's b, air-entry potential and hydraulic conductivity
#' are based on either on equations in Campbell, G. S. 1985. Soil Physics
#' with Basic: Transport Models for Soil-Plant Systems. Elsevier, Amsterdam,
#' or Cosby, B. J., G. M. Hornberger, R. B. Clapp, and T. R. Ginn. 1984.
#' A Statistical Exploration of the Relationships of Soil Moisture Characteristics
#' to the Physical Properties of Soils. Water Resources Research 20:682-690.
#' Field capacity and permanent wilting point are based on Rab, M. A., S. Chandra, P.
#' D. Fisher, N. J. Robinson, M. Kitching, C. D. Aumann, and M. Imhof. 2011. Modelling
#' and prediction of soil water contents at field capacity and permanent wilting
#' point of dryland cropping soils. Soil Research 49:389-407.
#' @param soilpro Matrix of n x 5 matrix of soil composition with the following columns 1. depth (cm), 2. bulk density (Mg/m3), 3. clay (\%), 4. silt (\%), 5. sand (\%)
#' @param model Choice of equation to compute soil hydraulic parameters (see details)
#' @param DEP sequence of depths at which results are required, within the range provided by column 1 of input table 'soilpro'
#' @return PE air entry water potential (J/kg), Campbell (1985) eq. 5.12, p. 46 or Cosby et al. (1984) Table 5
#' @return BB Campbell's b parameter, Campbell (1985) eq. 5.11, p. 45 or Cosby et al. (1984) Table 5
#' @return BD bulk density, Mg/m3
#' @return KS saturated hydraulic conductivity (kg s / m3), Campbell (1985) eq. 6.12, p. 54 or Cosby et al. (1984) Table 5
#' @return FC Field capacity (m3/m3, \%) Based on model 6 in Table 6 of Rab, M. A., S. Chandra, P. D. Fisher, N. J. Robinson, M. Kitching, C. D. Aumann, and M. Imhof. 2011. Modelling and prediction of soil water contents at field capacity and permanent wilting point of dryland cropping soils. Soil Research 49:389-407.
#' @return PWP Permanent Wilting Point (m3/m3, \%) Based on model 2 in Table 7 of Rab et al. 2011 (cited above)
#' @usage pedotransfer(soilpro, 1, DEP)
#' @details
#' \code{model}{ = 1, choose the pedotransfer function to use; 0 for Cosby et al. (1984) univariate regression (their Table 5), for Cosby et al. (1984) multivariate regression (their Table 4), or 2 for Campbell (1985)}\cr\cr
#' @export
pedotransfer <- function(soilpro = as.data.frame(soilpro), model = 1, DEP = soilpro[,1]){

  soil_depths <- soilpro[, 1]
  if(length(DEP) != nrow(soilpro)){
   nodeout <- length(DEP) * 2 - 2
  }else{
   nodeout <- length(DEP)
  }
  # spline to new depths if needed
  if(nodeout != length(DEP) & nodeout != 0){

    # find half-way points between given depths
    DEP2 <- rep(0, nodeout)
    j <- 1
    for(i in 1:length(DEP2)){ # loop through all depths
      if(i%%2 == 0){ # alternates between nodes to get half way points
        DEP2[i] <- DEP2[i - 1] + (DEP[j] - DEP2[i - 1]) / 2
      }else{
        DEP2[i] <- DEP[j]
        j <- j + 1
      }
    }
    DEP2 <- as.data.frame(floor(DEP2))
    colnames(DEP2) <- "DEPTH"

    BD_spline <- spline(soil_depths,soilpro[, 2], n = 201, xmin = 0, xmax = 200, method = 'natural')
    BD_spline <- as.data.frame(cbind(BD_spline$x,BD_spline$y))
    colnames(BD_spline) <- c('DEPTH', 'VALUE')
    BD <- merge(DEP2, BD_spline)
    BD <- c(BD[1, 2], BD[, 2])

    sand_spline <-spline(soil_depths,soilpro[, 5], n = 201, xmin = 0, xmax = 200, method = 'natural')
    sand_spline <- as.data.frame(cbind(sand_spline$x,sand_spline$y))
    colnames(sand_spline) <- c('DEPTH', 'VALUE')
    sand <- merge(DEP2, sand_spline)
    sand <- c(sand[1,2], sand[,2])
    sand[sand < 0] <- 0
    sand[sand > 100] <- 100

    silt_spline <-spline(soil_depths,soilpro[, 4], n = 201, xmin = 0, xmax = 200, method = 'natural')
    silt_spline <- as.data.frame(cbind(silt_spline$x, silt_spline$y))
    colnames(silt_spline) <- c('DEPTH', 'VALUE')
    silt <- merge(DEP2, silt_spline)
    silt <- c(silt[1, 2], silt[, 2])
    silt[silt < 0] <- 0
    silt[silt > 100] <- 100

    clay_spline <-spline(soil_depths,soilpro[, 3], n = 201, xmin = 0, xmax = 200, method = 'natural')
    clay_spline <- as.data.frame(cbind(clay_spline$x, clay_spline$y))
    colnames(clay_spline) <- c('DEPTH', 'VALUE')
    clay <- merge(DEP2, clay_spline)
    clay <- c(clay[1, 2], clay[, 2])
    clay[clay < 0] <- 0
    clay[clay > 100] <- 100
  }else{
    BD <- soilpro[, 2]
    sand <- soilpro[, 5]
    silt <- soilpro[, 4]
    clay <- soilpro[, 3]
  }# end check if splining

  # field capacity, model 6 from Table 6 of Rab et al. (2011).
  FC <- (7.561 + 1.176 * clay - 0.009843 * clay ^ 2 + 0.2132 * silt) / 100

  # permanent wilting point, model 2, Table 7 of Rab et al. (2011)
  PWP <- (-1.304 + 1.117 * clay - 0.009309 * clay ^ 2) / 100

  if(model == 0){ # use Cosby et al. (1984) equations, from Table 5

    # Campbell's b parameter
    BB <- clay * Cosby1984Tbl5$slope[1] + Cosby1984Tbl5$intercept[1]
    BB <- log10(10^(BB) / 10.2) # un-log, convert to J/kg, relog


    # air entry water potential (J/kg)
    PE <- sand * Cosby1984Tbl5$slope[2] + Cosby1984Tbl5$intercept[2]
    PE <- 10 ^ (PE) / 10.2 * -1 # un-log, convert to J/kg

    # saturated hydraulic conductivity (kg s / m3), Campbell (1985) eq. 6.12, p. 54
    KS <- sand * Cosby1984Tbl5$slope[3] + Cosby1984Tbl5$intercept[3]
    KS <- 10 ^ (KS) * 0.0007196666 # un-log, convert to kg s m-3
  }

  if(model == 1){ # use Cosby et al. (1984) equations, from Table 4

    # Campbell's b parameter
    BB <- clay * Cosby1984Tbl4$slope[1] + sand * Cosby1984Tbl4$slope[2] + Cosby1984Tbl4$intercept[1]
    BB <- log10(10 ^ (BB) / 10.2) # un-log, convert to J/kg, relog

    # air entry water potential (J/kg)
    PE <- sand * Cosby1984Tbl4$slope[3] + silt * Cosby1984Tbl4$slope[4] + Cosby1984Tbl4$intercept[3]
    PE <- 10 ^ (PE) / 10.2 * -1 # un-log, convert to J/kg

    # saturated hydraulic conductivity (kg s / m3), Campbell (1985) eq. 6.12, p. 54
    KS <- sand * Cosby1984Tbl4$slope[5] + clay * Cosby1984Tbl4$slope[6] + Cosby1984Tbl4$intercept[5]
    KS <- 10 ^ (KS) * 0.0007196666 # un-log, convert to kg s m-3
  }
  if(model == 2){ # use Campbell (1985) equations

    # particle diameters from Campbell (1985) p. 10
    dclay <- 0.001 #mm
    dsilt <- 0.026 #mm
    dsand <- 1.05 #mm

    # Campbell (1985) eq. 2.17, p. 9
    a <- (clay / 100) * log(dclay) + (sand / 100) * log(dsand) + (silt / 100) * log(dsilt)

    # Campbell (1985) eq. 2.18, p. 9
    b.1 <- (((clay / 100) * log(dclay) ^ 2 + (sand / 100) * log(dsand) ^ 2 + (silt / 100) * log(dsilt) ^ 2) - a ^ 2) ^ (1 / 2)

    # Campbell (1985) eq. 2.15, p. 9
    dg <- exp(a) # geometric mean of particle diameter

    # Campbell (1985) eq. 2.16, p. 9
    sigma_g <- exp(b.1) # geometric standard deviation of particle diameter

    # Reference air entry water potential Campbell (1985) eq. 5.10, p. 45
    PES <- (0.5 * dg ^ (-1 / 2)) * -1 # air entry water potential reference value at bulk density of 1.3 Mg/m3

    # Campbell's b parameter, Campbell (1985) eq. 5.11, p. 45
    BB <- -2 * PES + 0.2 * sigma_g # slope of ln of water potential against ln of volumetric water content

    # air entry water potential (J/kg), Campbell (1985) eq. 5.12, p. 46
    PE <- PES * (BD / 1.3) ^ (0.67 * BB) #

    # saturated hydraulic conductivity (kg s / m3), Campbell (1985) eq. 6.12, p. 54
    KS <- 0.004 * (1.3 / BD) ^ (1.3 * BB) * exp(-6.9 * clay / 100 - 3.7 * silt / 100)
  }
  return(list(BD = BD, BB = BB, KS = KS, PE = PE, FC = FC, PWP = PWP))
}
