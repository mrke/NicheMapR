#' One-lump Transient Heat Budget for Constant Environment
#'
#' Transient, 'one-lump', heat budget for computing rate of change of temperature
#' under a constant environment
#' Michael Kearney, Raymond Huey and Warren Porter developed this R function and example in September 2017.
#' @param t = seq(1,3600,60), time intervals (s) at which output is required
#' @param Tc_init = 5, initial temperature (°C)
#' @param Ww_g = 500, animal weight (g)
#' @param rho_body = 932, animal density (kg/m3)
#' @param q = 0, metabolic heat production rate W/m3
#' @param c_body = 3073, Specific heat of flesh J/(kg-K)
#' @param k_flesh = 0.5, Thermal conductivity of flesh (W/mK, range: 0.412-2.8)
#' @param emis = 0.95, Emissivity of animal (0-1)
#' @param abs = 0.85, solar absorptivity, decimal percent
#' @param geom = 2, Organism shape, 0-5, Determines whether standard or custom shapes/surface area/volume relationships are used: 0=plate, 1=cyl, 2=ellips, 3=lizard (desert iguana), 4=frog (leopard frog), 5=custom (see parameter 'shape_coefs')
#' @param shape_b = 1/5, Proportionality factor (-) for going from volume to area, represents ratio of width:height for a plate, length:diameter for cylinder, b axis:a axis for ellipsoid
#' @param shape_c = 1/5, Proportionality factor (-) for going from volume to area, represents ratio of length:height for a plate, c axis:a axis for ellipsoid
#' @param shape_coefs = c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743), Custom shape coefficients. Operates if geom=5, and consists of 4 pairs of values representing the parameters a and b of a relationship AREA=a*Ww_g^b, where AREA is in cm2 and Ww_g is in g. The first pair are a and b for total surface area, then a and b for ventral area, then for sillhouette area normal to the sun, then sillhouette area perpendicular to the sun
#' @param posture = 'n' pointing normal 'n', parallel 'p' to the sun's rays, or 'a' in between?
#' @param orient = 1, does the object orient toward the sun? (0,1)
#' @param fatosk = 0.4, Configuration factor to sky (-) for infrared calculations
#' @param fatosb = 0.4, Configuration factor to subsrate for infrared calculations
#' @param alpha_sub = 0.2, substrate solar reflectivity, decimal percent
#' @param pdif = 0.1, proportion of solar energy that is diffuse (rather than direct beam)
#' @param Tair = 30, air temperature (°C)
#' @param Trad = 30, radiant temperature (°C), averaging ground and sky
#' @param vel = 0.1, wind speed (m/s)
#' @param Qsol = 500, solar radiation (W/m2)
#' @param Zen = 20, zenith angle of sun (90 is below horizon), degrees
#' @param press = 101325, air pressure (Pa)
#' @return Tc Core temperature (°C)
#' @return Tcf Final (steady state) temperature (°C), if conditions remained constant indefinately
#' @return tau Time constant (s)
#' @return dTc Rate of change of core temperature (°C/s)
#' @usage onelump(t, Tc_init, Ww_g, geom, Tair, Trad, vel, Qsol, Zen, ...)
#' @examples
#' library(NicheMapR)
#'
#' # compare heating rates of a variety of different-sized objects
#'
#' t=seq(1,3600*2,60) # times (in seconds) to report back values - every minute for 2 hrs
#' tmins <- t/60
#'
#' par(mfrow = c(1,2))
#' Ww_g <- 5 # body weight, g
#' Tc_init <- 20 # initial body temperature, °C
#' geom <- 2 # shape (2 = ellipsoid)
#' Tair <- 20 # air temperature, °C
#' Trad <- Tair # radiant temperature, °C
#' vel <- 1 # wind speed, m/s
#' Qsol <- 500 # horizontal plane solar radiation, W/m2
#' Zen <- 20 # zenith angle of sun, degrees
#' alpha <- 0.85 # solar absorptivity, -
#'
#' Tbs<-onelump(t=t, alpha = alpha, Tc_init = Tc_init, Ww_g = Ww_g,
#'   geom = geom, Tair = Tair, Trad = Trad, vel = vel, Qsol = Qsol, Zen = Zen)
#' plot(Tbs$Tc ~ tmins, type= 'l' ,col = 1, ylim = c(20, 30), ylab = 'Temperature, °C',xlab='time, min', las = 1)
#' text(80, 27, "    500 g")
#' text(80, 24.5, "5 g")
#' text(90, 20.5, "Tair for both sizes", col = "blue")
#' text(90, 26, "   vel = 1.0 m/s")
#' text(90, 23.5, "     vel = 1.0 m/s")"     vel = 1.0 m/s")
#'
#' Ww_g <- 500 # body weight, g
#' Tbs<-onelump(t=t, alpha = alpha, Tc_init = Tc_init, Ww_g = Ww_g,
#'   geom = geom, Tair = Tair, Trad = Trad, vel = vel, Qsol = Qsol, Zen = Zen)
#' points(Tbs$Tc~tmins,type='l',lty = 2, col=1)
#' abline(Tair,0, lty = 1, col = 'blue')
#' abline(h = Tair + .1, lty = 2, col = 'blue')
#'
#' Ww_g <- 5 # body weight, g
#' Tair <- 25 # air temperature, °C
#' vel <-0.5 # wind speed, m/s
#'
#' Tbs<-onelump(t=t, alpha = alpha, Tc_init = Tc_init, Ww_g = Ww_g,
#'   geom = geom, Tair = Tair, Trad = Trad, vel = vel, Qsol = Qsol, Zen = Zen)
#' plot(Tbs$Tc~tmins,type='l',col=1,ylim=c(20,30),ylab='Temperature, °C',xlab='time, min', las = 1)
#' abline(h = Tair, lty = 1, col = 'blue')
#'
#' Ww_g <- 500 # body weight, g
#' Tair <- 20 # air temperature, °C
#' vel <-1 # wind speed, m/s
#'
#' Tbs<-onelump(t=t, alpha = alpha, Tc_init = Tc_init, Ww_g = Ww_g,
#'   geom = geom, Tair = Tair, Trad = Trad, vel = vel, Qsol = Qsol, Zen = Zen)
#' points(Tbs$Tc~tmins,type='l',lty = 2, col=1)
#' abline(h = Tair, lty = 2, col = 'blue')
#'
#' text(65, 29.65, "5 g")
#' text(85, 27, "500 g")
#' text(90, 20.5, "Tair for large animal", col = "blue")
#' text(90, 25.4, "Tair for small animal", col = "blue")
#' text(80, 28.65, "vel = 0.5 m/s")
#' text(93, 26, "vel = 1.0 m/s")
#' @export
onelump<-function(t = seq(1, 3600, 60), Tc_init = 5, Ww_g = 500,
  geom = 2, Tair = 30, Trad = 30, vel = 0.1, Qsol = 500, Zen = 20, k_flesh = 0.5,
  q = 0, c_body = 3073, emis = 0.95, rho_body = 932, alpha = 0.85,
  shape_coefs = c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743),
  shape_b = 1/5, shape_c = 1/5, posture = 'n', orient = 1, fatosk = 0.4, fatosb = 0.4,
  alpha_sub = 0.8, pdif = 0.1, press = 101325){

  sigma <- 5.67e-8 #Stefan-Boltzman, W/(m.K)
  Zenith <- Zen * pi / 180 # zenith angle in radians
  Tc <- Tc_init
  Tskin <- Tc + 0.1
  vel[vel < 0.1] <- 0.1 # don't let wind speed go too low
  S2 <- 0.0001 # shape factor, arbitrary initialization, because not defined except for ellipsoid
  DENSTY <- press / (287.04 * (Tair + 273.15)) # air density, kg/m3
  THCOND <- 0.02425 + (7.038 * 10 ^ -5 * Tair) # air thermal conductivity, W/(m.K)
  VISDYN <- (1.8325 * 10 ^ -5 * ((296.16 + 120) / ((Tair + 273.15) + 120))) * (((Tair + 273.15) / 296.16) ^ 1.5) # dynamic viscosity of air, kg/(m.s)

  # geometry section ############################################################
  m <- Ww_g / 1000 # convert weight to kg
  C <- m * c_body # thermal capacitance, J/K
  V <- m / rho_body # volume, m3
  Q_gen <- q * V # total metabolic heat, J
  L <- V ^ (1 / 3) # characteristic dimension, m

  # FLAT PLATE geometry
  if (geom == 0) {
    AHEIT <- (V / (shape_b * shape_c)) ^ (1 / 3) # length, m
    AWIDTH <- AHEIT * shape_b # width, m
    ALENTH <- AHEIT * shape_c # height, m
    ATOT <- ALENTH * AWIDTH * 2 + ALENTH * AHEIT * 2 + AWIDTH * AHEIT * 2 # total area, m2
    ASILN <- ALENTH * AWIDTH # max silhouette area, m2
    ASILP <- AWIDTH * AHEIT # min silhouette area, m2
    L <- AHEIT # characteristic dimension, m
    if (AWIDTH <= ALENTH) {
      L <- AWIDTH
    }else{
      L <- ALENTH
    }
    R <- ALENTH / 2 # 'radius', m
  }

  # CYLINDER geometry
  if (geom == 1) {
    R1 <- (V / (pi * shape_b * 2)) ^ (1 / 3) # radius, m
    ALENTH <- 2 * R1 * shape_b # length, m
    ATOT <- 2 * pi * R1 ^ 2 + 2 * pi * R1 * ALENTH # total surface area, m2
    AWIDTH <- 2 * R1 # width, m
    ASILN <- AWIDTH * ALENTH # max silhouette area, m2
    ASILP <- pi * R1 ^ 2 # min silhouette area, m2
    L <- ALENTH # characteristic dimension, m
    R2 <- L / 2
    if (R1 > R2) {
      # choose shortest dimension as R
      R <- R2
    } else{
      R <- R1
    }
  }

  # Ellipsoid geometry
  if (geom == 2) {
    A1 <- ((3 / 4) * V / (pi * shape_b * shape_c)) ^ 0.333 # axis A, m
    B1 <- A1 * shape_b # axis B, m
    C1 <- A1 * shape_c # axis C, m
    P1 <- 1.6075 # a constant
    ATOT <- (4 * pi * (((A1 ^ P1 * B1 ^ P1 + A1 ^ P1 * C1 ^ P1 + B1 ^ P1 * C1 ^ P1)        ) / 3) ^ (1 / P1)) # total surface area, m2
    ASILN <- max(pi * A1 * C1, pi * B1 * C1) # max silhouette area, m2
    ASILP <- min(pi * A1 * C1, pi * B1 * C1) # min silhouette area, m2
    S2 <- (A1 ^ 2 * B1 ^ 2 * C1 ^ 2) / (A1 ^ 2 * B1 ^ 2 + A1 ^ 2 * C1 ^ 2 + B1 ^ 2 * C1 ^ 2) # fraction of semi-major and minor axes, see Porter and Kearney 2009 supp1
    k_flesh <- 0.5# + 6.14 * B1 + 0.439 # thermal conductivity of flesh as a function of radius, see Porter and Kearney 2009
  }

  # Lizard geometry - DESERT IGUANA (PORTER ET AL. 1973 OECOLOGIA)
  if (geom == 3) {
    ATOT <- (10.4713 * Ww_g ^ .688) / 10000. # total surface area, m2
    AV <- (0.425 * Ww_g ^ .85) / 10000. # ventral surface area, m2
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984
    ASILN <- (3.798 * Ww_g ^ .683) / 10000. # Max. silhouette area (normal to the sun), m2
    ASILP <- (0.694 * Ww_g ^ .743) / 10000. # Min. silhouette area (pointing toward the sun), m2
    R <- L
  }

  # Frog geometry - LEOPARD FROG (C.R. TRACY 1976 ECOL. MONOG.)
  if (geom == 4) {
    ATOT <- (12.79 * Ww_g ^ 0.606) / 10000. # total surface area, m2
    AV <- (0.425 * Ww_g ^ 0.85) / 10000. # ventral surface area, m2
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: EQ'N 11 TRACY 1976
    ZEN <- 0
    PCTN <- 1.38171E-06 * ZEN ^ 4 - 1.93335E-04 * ZEN ^ 3 + 4.75761E-03 * ZEN ^ 2 - 0.167912 * ZEN + 45.8228
    ASILN <- PCTN * ATOT / 100 # Max. silhouette area (normal to the sun), m2
    ZEN <- 90
    PCTP <- 1.38171E-06 * ZEN ^ 4 - 1.93335E-04 * ZEN ^ 3 + 4.75761E-03 * ZEN ^ 2 - 0.167912 * ZEN + 45.8228
    ASILP <- PCTP * ATOT / 100 # Min. silhouette area (pointing toward the sun), m2
    R <- L
  }

  # user defined geometry
  if (geom == 5) {
    ATOT <- (shape_coefs[1] * Ww_g ^ shape_coefs[2]) / 10000. # total surface area, m2
    AV <- (shape_coefs[3] * Ww_g ^ shape_coefs[4]) / 10000 # ventral surface area, m2
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984
    # User must define Max. silhouette area (normal to the sun)
    ASILN <- (shape_coefs[5] * Ww_g ^ shape_coefs[6]) / 10000 # Max. silhouette area (normal to the sun), m2
    # User must define Min. silhouette area (pointing toward the sun)
    ASILP <- (shape_coefs[7] * Ww_g ^ shape_coefs[8]) / 10000 # Min. silhouette area (pointing toward the sun), m2
    R <- L
  }
  # end geometry section ############################################################

  if (max(Zen) >= 89) {
    Q_norm <- 0
  } else{
    if(orient == 1){
      Q_norm <- (Qsol / cos(Zenith))
    }else{
      Q_norm <- Qsol
    }
  }
  if (Q_norm > 1367) {
    Q_norm <- 1367 #making sure that low sun angles don't lead to solar values greater than the solar constant
  }
  if (posture == 'p') {
    Q_abs <- (Q_norm * (1 - pdif) * ASILP + Qsol * pdif * fatosk * ATOT + Qsol * (1 - alpha_sub) * fatosb * ATOT) * alpha
  }
  if (posture == 'n') {
    Q_abs <- (Q_norm * (1 - pdif) * ASILN + Qsol * pdif * fatosk * ATOT + Qsol * (1 - alpha_sub) * fatosb * ATOT) * alpha
  }
  if (posture == 'a') {
    Q_abs <- (Q_norm * (1 - pdif) * (ASILN + ASILP) / 2 + Qsol * pdif * fatosk * ATOT + Qsol * (1 - alpha_sub) * fatosb * ATOT) * alpha
  }

  Re <- DENSTY * vel * L / VISDYN # Reynolds number
  PR <- 1005.7 * VISDYN / THCOND # Prandlt number

  if (geom == 0) {
    NUfor <- 0.102 * Re ^ 0.675 * PR ^ (1. / 3.)
  }
  if (geom == 3 | geom == 5) {
    NUfor <- 0.35 * Re ^ 0.6
  }
  if (geom == 1) {
    #       FORCED CONVECTION OF A CYLINDER
    #       ADJUSTING NU - RE CORRELATION FOR RE NUMBER (P. 260 MCADAMS,1954)
    if (Re < 4) {
      NUfor = 0.891 * Re ^ 0.33
    } else{
      if (Re < 40) {
        NUfor = 0.821 * Re ^ 0.385
      } else{
        if (Re < 4000) {
          NUfor = 0.615 * Re ^ 0.466
        } else{
          if (Re < 40000) {
            NUfor = 0.174 * Re ^ 0.618
          } else{
            if (Re < 400000) {
              NUfor = 0.0239 * Re ^ 0.805
            } else{
              NUfor = 0.0239 * Re ^ 0.805
            }
          }
        }
      }
    }
  }
  if (geom == 2 | geom == 4) {
    NUfor <- 0.35 * Re ^ (0.6) # Nusselt number, forced convection
  }
  h_conv_forced <- NUfor * THCOND / L # convection coefficent, forced

  GR <- abs(DENSTY ^ 2 * (1 / (Tair + 273.15)) * 9.80665 * L ^ 3 * (Tskin - Tair) / VISDYN ^ 2) # Grashof number
  Raylei <- GR * PR # Rayleigh number

  # get Nusselt for Free Convection
  if (geom == 0) {
    NUfre = 0.55 * Raylei ^ 0.25
  }
  if (geom == 1 | geom == 3 | geom == 5) {
    if (Raylei < 1.0e-05) {
      NUfre = 0.4
    } else{
      if (Raylei < 0.1) {
        NUfre = 0.976 * Raylei ^ 0.0784
      } else{
        if (Raylei < 100) {
          NUfre = 1.1173 * Raylei ^ 0.1344
        } else{
          if (Raylei < 10000.) {
            NUfre = 0.7455 * Raylei ^ 0.2167
          } else{
            if (Raylei < 1.0e+09) {
              NUfre = 0.5168 * Raylei ^ 0.2501
            } else{
              if (Raylei < 1.0e+12) {
                NUfre = 0.5168 * Raylei ^ 0.2501
              } else{
                NUfre = 0.5168 * Raylei ^ 0.2501
              }
            }
          }
        }
      }
    }
  }

  if (geom == 2 | geom == 4) {
    Raylei = (GR ^ 0.25) * (PR ^ 0.333)
    NUfre = 2 + 0.60 * Raylei
  }
  h_conv_free <- NUfre * THCOND / L # convection coefficent, forced
  h_conv <- h_conv_free + h_conv_forced # combined convection coefficient
  Nu <- h_conv * L / THCOND # Nu combined
  #R_conv <- 1 / (h_conv * ATOT) # convective resistance, eq. 3 of Transient Equations Derivation vignette
  h_rad <- 4 * emis * sigma * ((Tc + Trad) / 2 + 273.15) ^ 3 # radiation heat transfer coefficient, eq. 46 of Transient Equations Derivation vignette

  if(geom == 2){ # ellipsoid
    j <- (Q_abs + Q_gen + h_conv * ATOT * ((q * S2) / (2 * k_flesh) + Tair) + h_rad * ATOT * ((q * S2) / (2 * k_flesh) + Trad)) / C #based on eq. 48 of Transient Equations Derivation vignette
  }else{ # assume cylinder
    j <- (Q_abs + Q_gen + h_conv * ATOT * ((q * R ^ 2) / (4 * k_flesh) + Tair) + h_rad * ATOT * ((q * R ^ 2) / (2 * k_flesh) + Trad)) / C #based on eq. 48 of Transient Equations Derivation vignette
  }
  theta_Tc <- ATOT * (Tc * h_conv + Tc * h_rad) / C #based on eq. 48 of Transient Equations Derivation vignette
  theta <- ATOT * (h_conv + h_rad) / C #based on eq. 48 of Transient Equations Derivation vignette
  Tcf <- j / theta # final Tc = j/theta, based on eq. 23 of Transient Equations Derivation vignette
  Tci <- Tc # initial temperature
  Tc <- (Tci - Tcf) * exp(-1 * theta * t) + Tcf # Tc at time t, based on eq. 33 of Transient Equations Derivation vignette
  tau <- 1 / theta # time constant
  dTc <- j - theta_Tc # rate of temperature change (°C/sec)
  return(list(Tc = Tc, Tcf = Tcf, tau = tau, dTc = dTc))
}
