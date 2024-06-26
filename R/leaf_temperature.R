#' Leaf temperature calculation
#'
#' Campbell & Norman's (1998) leaf temperature calculation based on their equation 14.1
#' solved using R's root-finding function 'uniroot'. Uses the functions SOLAR_ecto and
#' RADIN_ecto to calculate incoming radiation based on TGRD, TSKY, QSOLR, Z and PDIF.
#' The ectotherm function (and ectoR_devel) can also be used to calculate leaf temperature
#' and is faster and more general (e.g. can simulate different shapes) than the Campbell
#' & Norman function. Both methods are demonstrated in the example for this function.
#' Michael Kearney developed this R function and its sub-function leaf_bal and example in July 2023.
#' Campbell, G. S., & Norman, J. M. (1998). Environmental Biophysics. Springer.
#' @param w = 0.1, leaf width, m
#' @param A = 0.01, leaf surface area, m^2
#' @param A_sil = 0.005, leaf silhouette area, m^2
#' @param alpha_L = 0.5, leaf solar absorptivity, -
#' @param alpha_S = 0.8, ground solar absorptivity, -
#' @param epsilon_L = 0.97, leaf emissivity, -
#' @param epsilon_sub = 0.97, substrate emissivity, -
#' @param epsilon_sky = 0.78, sky emissivity, -
#' @param g_vs_ab = 0.5, vapour conductance, abaxial (top of leaf), mol/m2/s
#' @param g_vs_ad = 0, vapour conductance, adaxial (bottom of leaf), mol/m2/s
#' @param TA = 30, air temperature at leaf height from microclimate model, deg C
#' @param TGRD = 60, ground temperature from microclimate model, deg C
#' @param TSKY = 10, sky temperature from microclimate model, deg C
#' @param VEL = 2, wind speed from microclimate model, m/s
#' @param RH = 20, relative humidity from microclimate model, \%
#' @param QSOLR = 1000, total horizontal plane solar radiation from microclimate model, W/m2
#' @param Z = 0, solar zenith angle, degrees
#' @param PRESS = 101325, atmospheric pressure, Pa
#' @param PDIF = 0.15, proportion of solar radiation that is diffuse, -
#' @param conv_enhance = 1.4, outdoor convective enhancement factor, Mitchell 1976
#' @param heliotropic = 0, orient towards sun (1) or not (0)
#' @return leaf temperature (°C)
#' @usage leaf_temperature(w = 0.1, A = 0.01, A_sil = 0.005, alpha_L = 0.5, alpha_S = 0.85, g_vs_ab = 0.15, g_vs_ad = 0.15, TA = 30, TGRD = 60, TSKY = 10, VEL = 2, RH = 20, QSOLR = 1000, Z = 0, PDIF = 15)
#' @examples
#' library(NicheMapR)
#'
#' # simulate a microclimate
#' loc <- c(141.86, -29.05) # outback Australia
#' micro <<- micro_global(loc = loc, microclima = 1) # setting microclima = 1 to get partitioned diffuse and direct solar - needs internet connection for the DEM download
#' dates <- micro$dates # mock dates
#' metout <- as.data.frame(micro$metout) # microclimate aboveground conditions
#' soil <- as.data.frame(micro$soil) # soil temperature
#'
#' # obtain relevant microclimate conditions
#' P_a <- get_pressure(micro$elev/288) # hourlydata$pressure
#' PRESSs <- rep(P_a, nrow(metout)) # atmospheric pressure, Pa
#' TAs <- metout$TALOC # air temperature at leaf height from microclimate model, deg C
#' TGRDs <- soil$D0cm # ground temperature from microclimate model, deg C
#' TSKYs <- metout$TSKY # sky temperature from microclimate model, deg C
#' VELs <- metout$VLOC # wind speed from microclimate model, m/s
#' RHs <- metout$RHLOC # relative humidity from microclimate model, %
#' QSOLRs <- metout$SOLR # total horizontal plane solar radiation from microclimate model, W/m2
#' Zs <- metout$ZEN # solar zenith angle, degrees
#' PDIFs <- micro$diffuse_frac # use variable diffuse fraction
#' epsilon_sky <- 1 # emissivity has already been incorporated in calculation of TSKYs
#' alpha_S <- 1 - micro$REF # substrate solar absorptivity
#'
#' # leaf functional traits required for the heat and water budget
#' Ww_g <- 1 # wet weight, g
#' shape <- 2 # 0=plate, 1=cylinder, 2=ellipsoid
#' g_vs_ab <- 0.3 # leaf vapour conductance, abaxial (bottom of leaf), mol/m2/s
#' g_vs_ad <- 0.0 # leaf vapour conductance, adaxial (top of leaf), mol/m2/s
#' g_vs_base <- 0.01 # base leaf vapour conductance when stomata closed, equivalent to 0.1 (µmol H2O) / (m^2 s Pa) from tealeaves
#' shape_b <- 0.0025 # ratio of b axis:a axis for ellipsoid
#' shape_c <- 0.1176 # ratio of c axis:a axis for ellipsoid
#' epsilon_sub <- 0.95 # emissivity of the substrate (-)
#' epsilon_L <- 0.97 # emissivity of the leaf (-)
#' alpha_L <- 0.5 # solar absorptivity of the leaf (-)
#' fatosk <- 0.5 # radiation configuration factor to sky (-)
#' fatosb <- 0.5 # radiation configuration factor to substrate (-)
#' conv_enhance <- 1.4 # convective enhancement factor (-)
#' pct_cond <- 0 # percent of leaf conducting to the ground (%)
#'
#' # set up vapour conductance vectors and simulate stomatal closure at night
#' g_vs_abs <- rep(g_vs_ab, nrow(metout))
#' g_vs_ads <- rep(g_vs_ad, nrow(metout))
#' g_vs_abs[metout$ZEN == 90] <- 0 # close stomata when the sun is set
#' g_vs_ads[metout$ZEN == 90] <- 0 # close stomata when the sun is set
#' g_vs_abs <- g_vs_abs + g_vs_base / 2
#' g_vs_ads <- g_vs_ads + g_vs_base / 2
#'
#' # get characteristic dimension and areas using ecto_devel function GEOM_ecto
#' GEOM.out <- GEOM_ecto(AMASS = Ww_g / 1000, GEOMETRY = shape, SHP = c(1, shape_b, shape_c), PTCOND = pct_cond / 100, PMOUTH = 0, SKINW = 0 / 100)
#' w <- GEOM.out$AL / 0.7 # leaf width, m
#' A <- GEOM.out$AREA # total leaf surface area, m2
#' A_sil <- (GEOM.out$ASILN + GEOM.out$ASILP) / 2 # leaf silhoette area to direct solar radiation, m2
#'
#' # Campbell and Norman Calculation
#' T_leaf_CN <- unlist(lapply(1:length(TAs), function(x){leaf_temperature(
#'   w = w, # leaf width, m
#'   A = A, # leaf surface area, m^2
#'   A_sil = A_sil, # leaf silhouette area, m^2
#'   alpha_L = alpha_L, # leaf solar absorptivity, -
#'   alpha_S = alpha_S, # ground solar absorptivity, -
#'   epsilon_L = epsilon_L, # leaf emissivity, -
#'   epsilon_sub = epsilon_sub, # ground emissivity, -
#'   epsilon_sky = epsilon_sky, # sky emissivity, -
#'   g_vs_ab = g_vs_abs[x], # vapour conductance, abaxial (bottom of leaf), mol/m2/s
#'   g_vs_ad = g_vs_ads[x], # vapour conductance, adaxial (top of leaf), mol/m2/s
#'   TA = TAs[x], # air temperature at lizard height from microclimate model, deg C
#'   TGRD = TGRDs[x], # ground temperature from microclimate model, deg C
#'   TSKY = TSKYs[x], # sky temperature from microclimate model, deg C
#'   VEL = VELs[x], # wind speed from microclimate model, m/s
#'   RH = RHs[x], # relative humidity from microclimate model, %
#'   QSOLR = QSOLRs[x], # total horizontal plane solar radiation from microclimate model, W/m2
#'   Z = Zs[x], # solar zenith angle from microclimate model, degrees
#'   PRESS = PRESSs[x],
#'   PDIF = PDIFs[x], # proportion solar radiation that is diffuse, -
#'   conv_enhance = conv_enhance # convective enhancement factor
#' )})) # run leaf_temperature across environments
#'
#' month <- 1 # choose a month to plot
#' sub <- which(floor(dates) + 1 == month)
#' T_air_2m <- metout$TAREF
#' T_air_1cm <- metout$TALOC
#' plot(T_air_1cm[sub], type = 'l', col = 'blue', ylab = 'temperature, deg C', xlab = 'hour of day', ylim = c(15, 50))
#' points(T_air_2m[sub], type = 'l', col = 'blue', lty = 2)
#' points(T_leaf_CN[sub], type = 'l')
#'
#' # compare to ectotherm model calculation
#'
#' # model settings
#' live <- 0 # don't simulate behaviour or respiration
#' leaf <- 1 # use the stomatal conductance formulation for evaporation
#' postur <- 3 # no heliotropy
#'
#' ecto <- ectotherm(leaf = leaf,
#'                   pct_wet = 0,
#'                   g_vs_ab = g_vs_abs,
#'                   g_vs_ad = g_vs_ads,
#'                   preshr = PRESSs,
#'                   PDIF = PDIFs,
#'                   conv_enhance = conv_enhance,
#'                   Ww_g = Ww_g,
#'                   alpha_max = alpha_L,
#'                   alpha_min = alpha_L,
#'                   fatosk = fatosk,
#'                   fatosb = fatosb,
#'                   epsilon_sub = epsilon_sub,
#'                   epsilon = epsilon_L,
#'                   shape_b = shape_b,
#'                   shape_c = shape_c,
#'                   live = live,
#'                   pct_cond = pct_cond,
#'                   shape = shape,
#'                   postur = postur)
#' environ <- as.data.frame(ecto$environ)
#' T_leaf_NMR <- environ$TC
#'
#' points(T_leaf_NMR[sub], type = 'l', lty = 2)
#' legend(x = 3, y = 50, legend = c("T_air 1.2m", "T_air 1cm", "T_leaf_CN", "T_leaf_NMR"), col = c('blue', 'blue', 'black', 'black'), lty = c(2, 1, 1, 2), bty = "n")
#' @export
leaf_temperature <- function(w = 0.1, # leaf width, m
                   A = 0.01, # leaf surface area, m^2
                   A_sil = 0.005, # leaf silhouette area, m^2
                   alpha_L = 0.5, # leaf solar absorptivity, -
                   alpha_S = 0.85, # ground solar absorptivity, -
                   epsilon_L = 0.97, # leaf emissivity, -
                   epsilon_sub = 0.97, # substrate emissivity, -
                   epsilon_sky = 1, # sky emissivity, -
                   g_vs_ab = 0.15, # vapour conductance, abaxial (top of leaf), mol/m2/s
                   g_vs_ad = 0.15, # vapour conductance, adaxial (bottom of leaf), mol/m2/s
                   TA = 30, # air temperature at leaf height from microclimate model, deg C
                   TGRD = 60, # ground temperature from microclimate model, deg C
                   TSKY = 10, # sky temperature from microclimate model, deg C
                   VEL = 2, # wind speed from microclimate model, m/s
                   RH = 20, # relative humidity from microclimate model, %
                   QSOLR = 1000, # total horizontal plane solar radiation from microclimate model, W/m2
                   Z = 0, # solar zenith angle, degrees
                   PRESS = 101325, # atmospheric pressure, Pa
                   PDIF = 0.15, # proportion of solar radiation that is diffuse, -
                   conv_enhance = 1.4, # outdoor convective enhancement factor, Mitchell 1976
                   heliotropic = 0 # orient towards the sun?
){

  # constants
  sigma <- 5.67e-8 # Stefan-Boltzmann constant, W/m2/K
  ZEN <- Z * pi / 180
  Q_sol <- SOLAR_ecto(ATOT = A, ASIL = A_sil, AV = 0, AT = 0, ABSAN = alpha_L, ABSSB = alpha_S, FATOSK = 0.5, FATOSB = 0.5, FATOBJ = 0, ZEN = ZEN, QSOLR = QSOLR, PDIF = PDIF, SHADE = 0, LIVE = heliotropic)$QSOLAR
  #Q_sol <- 0.5 * A * QSOLR * alpha_L * (1 + (1 - alpha_S))
  Q_IR <- RADIN_ecto(TSKY = TSKY, TGRD = TGRD, ATOT = A, AV = 0, AT = 0, FATOSK = 0.5, FATOSB = 0.5, FATOBJ = 0, EMISAN = epsilon_L, EMISSB = epsilon_sub, EMISSK = epsilon_sky)
  #Q_IR <- 0.5 * A * sigma * epsilon_sub * epsilon_L * (TGRD + 273.15) ^ 4 + 0.5 * A * sigma * epsilon_L * epsilon_sky * (TSKY + 273.15) ^ 4
  R_abs <- (Q_sol + Q_IR) / A # total radiation absorbed W per m2
  T_L <- uniroot(function(x) leaf_bal(w = w, # leaf width, m
                                      epsilon_L = epsilon_L, # leaf emissivity, -
                                      g_vs_ab = g_vs_ab, # vapour conductance, abaxial (top of leaf), mol/m2/s
                                      g_vs_ad = g_vs_ad, # vapour conductance, adaxial (bottom of leaf), mol/m2/s
                                      T_a = TA, # air temperature, deg C
                                      P_a = PRESS, # atmospheric pressure, Pa
                                      RH = RH, # relative humidity, %
                                      u = VEL, # wind speed, m/s
                                      R_abs = R_abs, # total radiation absorbed, W/m2
                                      conv_enhance = conv_enhance, # outdoor convective enhancement factor, Mitchell 1976
                                      x), lower = -100, upper = 100)$root
return(T_L)
}

leaf_bal <- function(w = 0.1, # leaf width, m
                     epsilon_L = 0.97, # leaf emissivity, -
                     g_vs_ab = 0.5, # vapour conductance, abaxial (top of leaf), mol/m2/s
                     g_vs_ad = 0, # vapour conductance, adaxial (bottom of leaf), mol/m2/s
                     T_a = 38, # air temperature, deg C
                     P_a = 100000, # atmospheric pressure, Pa
                     RH = 16.6, # relative humidity, %
                     u = 1.5, # wind speed, m/s
                     R_abs = 750, # total radiation absorbed, W/m2
                     conv_enhance = 1.4, # outdoor convective enhancement factor, Mitchell 1976
                     x){
  # constants
  M_air <- 28.96 # molar mass of air, J/g
  sigma <- 5.67e-8 # Stefan-Boltzmann constant, W/m2/K
  lambda <- 40650 # latent heat of vaporisation of water, J/mol

  # derived traits - characteristic dimension and vapour conductance
  d <- 0.7 * w # leaf characteristic dimension, m
  g_va <- conv_enhance * 0.147 * (u / d) ^ (1 / 2) # boundary conductance, mol/m2/s, Campbell & Norman EQ 7.33
  g_v <- (0.5 * g_vs_ab * g_va) / (g_vs_ab + g_va) + (0.5 * g_vs_ad * g_va) / (g_vs_ad + g_va) # vapour conductance, mol/m2/s

  # air properties
  dryair.out <- DRYAIR(T_a, P_a, 0) # get properties of dry air using NicheMapR DRYAIR function
  rho_air <- dryair.out$densty
  mu <- dryair.out$visdyn
  RE <- rho_air * u * d / mu
  wetair.out <- WETAIR(db = T_a, rh = RH, bp = P_a) # get properties of wet air using NicheMapR WETAIR function
  e_a <- wetair.out$e
  c_p <- wetair.out$cp * (M_air / 1000) # specific heat of air, J/kg/K to J/mol/C

  # sensible heat exchange parameters
  g_Ha_forced <- conv_enhance * 0.135 * (u / d) ^ (1 / 2) # forced heat transfer conductance, mol/m2/s, Campbell & Norman Table 7.6
  g_Ha_free <- 0.05 * ((abs(x - T_a)) / d) ^ (1 / 4) # free heat transfer conductance, mol/m2/s, Campbell & Norman Table 7.6
  if(RE > 1000){
    g_Ha <- g_Ha_forced
  }else{
    if(RE < 10){
      g_Ha <- g_Ha_free
    }else{
      g_Ha <- (g_Ha_free ^ 3 + g_Ha_forced ^ 3) ^ (1 / 3)
    }
  }
  g_Ha <- (g_Ha_free ^ 3 + g_Ha_forced ^ 3) ^ (1 / 3)

  # heat balance to be solved for x, based on equation 14.1 in Campbell and Norman (1990)
  R_abs - epsilon_L * sigma * (x + 273.15) ^ 4 - c_p * g_Ha * (x - T_a) - lambda * g_v * (WETAIR(db = x)$esat - e_a) / P_a
}
