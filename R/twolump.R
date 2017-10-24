#' two-lump Transient Heat Budget (for use with deSolve package)
#'
#' Transient, 'one-lump', heat budget for computing rate of change of temperature
#' under a constant environment
#' Michael Kearney, Raymond Huey and Warren Porter developed this R function and example in September 2017.
#' @param t = seq(1,3600,60), time intervals (s) at which output is required
#' @param Tc_init = 5, initial core temperature (°C) 
#' @param Tsh_init = 5.1, initial shell temperature (°C) 
#' @param Ts_init = 5.2, initial surface temperature (°C) 
#' @param amass = 500, animal mass (g)
#' @param x_shell = 0.001, shell thickness, m
#' @param geom = 2, Organism shape, 0-5, Determines whether standard or custom shapes/surface area/volume relationships are used: 0=plate, 1=cyl, 2=ellips, 3=lizard (desert iguana), 4=frog (leopard frog), 5=custom (see parameter 'customallom')
#' @param k_inner = 0.5, Thermal conductivity of inner shell (W/mK, range: 0.412-2.8)
#' @param k_outer = 0.5, Thermal conductivity of outer shell (W/mK, range: 0.412-2.8)
#' @param q = 0, metabolic heat production rate W/m3
#' @param cp_inner = 3073, Specific heat of flesh J/(kg-K)
#' @param cp_outer = 3073, Specific heat of outer shell J/(kg-K)
#' @param rho = 932, animal density (kg/m3)
#' @param emis = 0.95, Emissivity of animal (0-1)
#' @param abs = 0.85, solar absorptivity, decimal percent
#' @param shape_coefs = c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743), Custom allometry coefficients. Operates if lometry=5, and consists of 4 pairs of values representing the parameters a and b of a relationship AREA=a*mass^b, where AREA is in cm2 and mass is in g. The first pair are a and b for total surface area, then a and b for ventral area, then for sillhouette area normal to the sun, then sillhouette area perpendicular to the sun
#' @param shape_b = 3, Proportionality factor (-) for going from volume to area, represents ratio of width:height for a plate, length:diameter for cylinder, b axis:a axis for ellipsoid
#' @param shape_c = 0.6666666667, Proportionality factor (-) for going from volume to area, represents ratio of length:height for a plate, c axis:a axis for ellipsoid
#' @param posture = 'n' pointing normal 'n', parallel 'p' to the sun's rays, or 'b' in between?
#' @param orient does the object orient toward the sun? (0,1)
#' @param fatosk = 0.4, Configuration factor to sky (-) for infrared calculations
#' @param fatosb = 0.4, Configuration factor to subsrate for infrared calculations
#' @param abs_sub = 0.2, substrate solar reflectivity, decimal percent
#' @param pctdif = 0.1, proportion of solar energy that is diffuse (rather than direct beam)
#' @param Tair = 30, air temperature (°C)
#' @param Trad = 30, radiant temperature (°C), averaging ground and sky
#' @param vel = 0.1, wind speed (m/s)
#' @param Qsol = 500, solar radiation (W/m2)
#' @param Zen = 20, zenith angle of sun (90 is below horizon), degrees
#' @param press = 101325, air pressure (Pa)
#' @return Tc Core temperature (deg C)
#' @return Tcf Final (steady state) temperature (deg C), if conditions remained constant indefinately
#' @return tau Time constant (s)
#' @return dTc Rate of change of core temperature (deg C/s)
#' @usage twolump(t, Tc_init, Tsh_init, Ts_init, mass, geom, Tair, Trad, vel, Qsol, Zen, ...)
#' @examples
#' library(NicheMapR)
#'
#' # compare heating rates of a variety of different-sized objects
#'
#' t=seq(1,3600*2,60) # times (in seconds) to report back values - every minute for 2 hrs
#' tmins <- t/60
#'
#' par(mfrow = c(1,2))
#' mass <- 5
#' x_shell <- 0.001
#' Tc_init <- 20
#' Tsh_init <- 20.1
#' Ts_init <- 20.2
#' geom <- 3
#' Tair <- 20
#' Trad <- Tair
#' vel <- 1
#' Qsol <- 500
#' Zen <- 20
#' abs <- 0.85
#'
#' Tbs<-twolump(t=t, abs = abs, Tc_init = Tc_init, Tsh_init = Tsh_init, Ts_init = Ts_init, mass = mass, x_shell = x_shell,
#'   geom = geom, Tair = Tair, Trad = Trad, vel = vel, Qsol = Qsol, Zen = Zen)
#' plot(Tbs$Tc ~ tmins, type= 'l' ,col = 1, ylim = c(20, 32), ylab = 'Temperature, °C',xlab='time, s', las = 1)
#' text(80, 29, "    5000 g")
#' text(80, 25, "5 g")
#' text (40, 20.5, "Tair for both sizes", col = "blue")
#' text(90, 28, "   vel = 1.0 m/s")
#' text(90, 24, "     vel = 1.0 m/s")
#'
#' mass <- 5000
#' Tbs<-twolump(t=t, abs = abs, Tc_init = Tc_init, Tsh_init = Tsh_init, Ts_init = Ts_init, mass = mass, x_shell = x_shell,
#'   geom = geom, Tair = Tair, Trad = Trad, vel = vel, Qsol = Qsol, Zen = Zen)
#' points(Tbs$Tc~tmins,type='l',lty = 2, col=1)
#' abline(Tair,0, lty = 1, col = 'blue')
#' abline(h = Tair + .1, lty = 2, col = 'blue')
#'
#' mass <- 5
#' Tair <- 25
#' Trad <- Tair
#' vel <-0.5
#'
#' Tbs<-twolump(t=t, abs = abs, Tc_init = Tc_init, Tsh_init = Tsh_init, Ts_init = Ts_init, mass = mass, x_shell = x_shell,
#'   geom = geom, Tair = Tair, Trad = Trad, vel = vel, Qsol = Qsol, Zen = Zen)
#' plot(Tbs$Tc~tmins,type='l',col=1,ylim=c(20,32),ylab='Temperature, °C',xlab='time, s', las = 1)
#' abline(h = Tair, lty = 1, col = 'blue')
#'
#' mass <- 5000
#' Tair <- 20
#' Trad <- Tair
#' vel <-1
#'
#' Tbs<-twolump(t=t, abs = abs, Tc_init = Tc_init, Tsh_init = Tsh_init, Ts_init = Ts_init, mass = mass, x_shell = x_shell,
#'   geom = geom, Tair = Tair, Trad = Trad, vel = vel, Qsol = Qsol, Zen = Zen)
#' points(Tbs$Tc~tmins,type='l',lty = 2, col=1)
#' abline(h = Tair, lty = 2, col = 'blue')
#'
#' text(65, 31, "5 g")
#' text(80, 29, "5000 g")
#' text(60, 20.5, "Tair for large animal", col = "blue")
#' text(63, 25.4, "Tair for small animal", col = "blue")
#' text(80, 30.1, "vel = 0.5 m/s")
#' text(93, 28, "vel = 1.0 m/s")
#' @export
twolump<-function(t = seq(1, 3600, 60), Tc_init = 5, Tsh_init = 5.1, Ts_init = 5.2, mass = 500, x_shell = 0.001,
  geom = 2, Tair = 30, Trad=30, vel = 0.1, Qsol = 500, Zen = 20, k_inner = 0.5, k_outer = 0.5,
  q = 0, cp_inner = 3073, cp_outer = 3073, emis = 0.95, rho = 932, abs = 0.85,
  customallom = c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743),
  shape_b = 0.5, shape_c = 0.5, posture = 'n', orient = 1, fatosk = 0.4, fatosb = 0.4,
  abs_sub = 0.8, pctdif = 0.1, press = 101325){
  
  sigma <- 5.67e-8 #Stefan-Boltzman, W/(m.K)
  Zenith <- Zen * pi / 180 # zenith angle in radians
  Tc <- Tc_init # core temperature, deg C
  Tsh <- Tsh_init # skin temperature
  Ts <- Ts_init # surface temperature
  vel[vel < 0.01] <- 0.01 # don't let wind speed go too low - always some free convection
  S2 <- 0.0001 # shape factor, arbitrary initialization, because not defined except for ellipsoid
  DENSTY <- press / (287.04 * (Tair + 273)) # air density, kg/m3
  THCOND <- 0.02425 + (7.038 * 10 ^ -5 * Tair) # air thermal conductivity, W/(m.K)
  VISDYN <- (1.8325 * 10 ^ -5 * ((296.16 + 120) / ((Tair + 273.15) + 120))) * (((Tair + 273.15) / 296.16) ^ 1.5) # dynamic viscosity of air, kg/(m.s)
  
  # geometry section ############################################################
  m <- mass / 1000 # convert mass to kg  C<-m*cp # thermal capacitance, J/K
  C <- m * cp_inner # thermal capacitance, J/K
  V <- m / rho # volume, m3
  Qgen <- q * V # total metabolic heat, J
  L <- V ^ (1 / 3) # characteristic dimension, m
  V_inner <- (L - shell) * 3
  V_shell <- V - V_inner
  Cs <- V_shell * rho * cp_shell
  Cc <- V_inner * rho * cp_inner  
  
  # FLAT PLATE geometry
  if (geom == 0) {
    AHEIT <- (V / (shape_b * shape_c)) ^ (1 / 3) # length, m
    AWIDTH <- AHEIT * shape_b # width, m
    ALENTH <- AHEIT * shape_c # height, m
    V_inner <- (AWIDTH - x_shell) * (ALENTH - x_shell) * (AHEIT - x_shell)
    V_shell <- V - V_inner
    Cs <- V_shell * rho * cp
    Cc <- V_inner * rho * cp
    ATOT <- ALENTH * AWIDTH * 2 + ALENTH * AHEIT * 2 + AWIDTH * AHEIT * 2 # total area, m2
    ASILN <- ALENTH * AWIDTH # max silhouette area, m2
    ASILP <- AWIDTH * AHEIT # min silhouette area, m2
    L <- AHEIT # characteristic dimension, m
    if (AWIDTH <= ALENTH) {
      L <- AWIDTH
    } else{
      L <- ALENTH
    }
    R <- ALENTH / 2 # 'radius', m
  }
  
  # CYLINDER geometry      
  if(geom == 1){
    R1 <- (V / (pi * shape_b * 2)) ^ (1 / 3) # radius, m
    ALENTH <- 2 * R1 * shape_b # length, m
    V_inner <- pi * (R1 - x_shell) ^ 2 * (ALENTH - x_shell)
    V_shell <- V - V_inner
    Cs <- V_shell * rho * cp
    Cc <- V_inner * rho * cp  
    ATOT<- 2 * pi * R1 ^ 2 + 2 * pi * R1 * ALENTH # total surface area, m2
    AWIDTH <- 2 * R1 # width, m
    ASILN <- AWIDTH * ALENTH # max silhouette area, m2
    ASILP <- pi * R1 ^ 2 # min silhouette area, m2
    L <- ALENTH # characteristic dimension, m
    R2 <- L / 2
    if(R1 > R2){ # choose shortest dimension as R
      R <- R2
    }else{
      R <- R1
    }
  }
  
  # Ellipsoid geometry
  if(geom == 2){
    A1<-((3 /4 ) * V / (pi * shape_b * shape_c)) ^ (1 / 3) # axis A, m   
    B1 <- A1 * shape_b # axis B, m
    C1 <- A1 * shape_c # axis C, m
    V_inner <- (4 / 3) * pi * (A1 - x_shell) * (B1 - x_shell) * (C1 - x_shell)
    V_shell <- V - V_inner
    Cs <- V_shell * rho * cp
    Cc <- V_inner * rho * cp 
    P1 <- 1.6075 # a constant
    ATOT <- (4 * pi * (((A1 ^ P1 * B1 ^ P1 + A1 ^ P1 * C1 ^ P1 + B1 ^ P1 * C1 ^ P1)) / 3) ^ (1 / P1)) # total surface area, m2
    ASILN <- max(pi * A1 * C1, pi * B1 * C1) # max silhouette area, m2
    ASILP <- min(pi * A1 * C1, pi * B1 * C1) # min silhouette area, m2
    S2 <- (A1 ^ 2 * B1 ^ 2 * C1 ^ 2) / (A1 ^ 2 * B1 ^ 2 + A1 ^ 2 * C1 ^ 2 + B1 ^ 2 * C1 ^ 2) # fraction of semi-major and minor axes, see Porter and Kearney 2009 supp1
    kflesh <- 0.5 # + 6.14 * B1 + 0.439 # thermal conductivity of flesh as a function of radius, see Porter and Kearney 2009
  }              
  
  # Lizard geometry - DESERT IGUANA (PORTER ET AL. 1973 OECOLOGIA)
  if (geom == 3) {
    ATOT <- (10.4713 * mass ^ .688) / 10000. # total surface area, m2
    AV <- (0.425 * mass ^ .85) / 10000. # ventral surface area, m2
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984
    ASILN <- (3.798 * mass ^ .683) / 10000. # Max. silhouette area (normal to the sun), m2
    ASILP <- (0.694 * mass ^ .743) / 10000. # Min. silhouette area (pointing toward the sun), m2
    R <- L
  }
  
  # Frog geometry - LEOPARD FROG (C.R. TRACY 1976 ECOL. MONOG.)
  if (geom == 4) {
    ATOT <- (12.79 * mass ^ 0.606) / 10000. # total surface area, m2
    AV <- (0.425 * mass ^ 0.85) / 10000. # ventral surface area, m2
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
    ATOT <- (shape_coefs[1] * mass ^ shape_coefs[2]) / 10000. # total surface area, m2
    AV <- (shape_coefs[3] * mass ^ shape_coefs[4]) / 10000 # ventral surface area, m2
    # NORMAL AND POINTING @ SUN SILHOUETTE AREA: PORTER & TRACY 1984
    # User must define Max. silhouette area (normal to the sun)
    ASILN <- (shape_coefs[5] * mass ^ shape_coefs[6]) / 10000 # Max. silhouette area (normal to the sun), m2
    # User must define Min. silhouette area (pointing toward the sun)
    ASILP <- (shape_coefs[7] * mass ^ shape_coefs[8]) / 10000 # Min. silhouette area (pointing toward the sun), m2
    R <- L
  }
  # end geometry section ############################################################
  
  if (max(Zen) >= 89) {
    Qnorm <- 0
  } else{
    if(orient == 1){
      Qnorm <- (Qsol / cos(Zenith))
    }else{
      Qnorm <- Qsol
    }
  }
  if (Qnorm > 1367) {
    Qnorm <- 1367 #making sure that low sun angles don't lead to solar values greater than the solar constant
  }
  if (posture == 'p') {
    Qabs <- (Qnorm * (1 - pctdif) * ASILP + Qsol * pctdif * fatosk * ATOT + Qsol * (1 - abs_sub) * fatosb * ATOT) * abs
  }
  if (posture == 'n') {
    Qabs <- (Qnorm * (1 - pctdif) * ASILN + Qsol * pctdif * fatosk * ATOT + Qsol * (1 - abs_sub) * fatosb * ATOT) * abs
  }
  if (posture == 'b') {
    Qabs <- (Qnorm * (1 - pctdif) * (ASILN + ASILP) / 2 + Qsol * pctdif * fatosk * ATOT + Qsol * (1 - abs_sub) * fatosb * ATOT) * abs
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
  hc_forced <- NUfor * THCOND / L # convection coefficent, forced
  
  GR <- abs(DENSTY ^ 2 * (1 / (Tair + 273.15)) * 9.80665 * L ^ 3 * (Ts - Tair) / VISDYN ^ 2) # Grashof number
  Raylei <- GR * PR # Rayleigh number
  
  # get Nusselt for Free Convect
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
  hc_free <- NUfre * THCOND / L # convection coefficent, forced
  hc <- hc_free + hc_forced # combined convection coefficient
  Nu <- hc * L / THCOND # Nu combined
  Rconv <- 1 / (hc * ATOT) # convective resistance, eq. 5 of Kearney, Huey and Porter 2017 Appendix 1
  
  hr <- 4 * emis * sigma * ((Tc + Trad) / 2 + 273.15) ^ 3 # radiation resistance, eq. 49 of Kearney, Huey and Porter 2017 Appendix 1
  
  Rrad <- 1 / (hr * ATOT)
  Rconv <- 1 / (hc * ATOT)
  
  Rs <- x_shell / (k_shell * ATOT) # resistance of skin
  Qgen <- 0
  Qresp <- 0
  m_bl <- 0.012 / 1000 / 60 / 1e-6 # blood flow rate kg/s/m3
  V_bl <- ATOT * 0.1 # blood volume
  Rb <- 1 / (m_bl * cp * V_bl) # blood resistance
  F <- 1 / (Cc * Rb) # from eq 110 and 111
  H <- (Qgen - Qresp) / Cc # from eq 110 and 111
  Ts <- (-1 * H + F * Tc) / F # from eq 112
  Tsh <- (Qabs + Trad / Rrad + Tair / Rconv + 2 * Ts / Rs) / (1 / Rrad + 1 / Rconv + 2 / Rs) # eq 106
  A <- (1 / Cs) * (1 / Rb + 2 / Rs - (2 / Rs) / (Rs / (2 * Rrad) + Rs / (2 * Rconv) + 1)) #eq 109 part 2
  B <- 1 / (Rb * Cs) # eq 109 part 3
  P1 <- (-1 * (F + A) + ((F + A) ^ 2 - 4 * F * A + 4 * F * B) ^ (1 / 2)) / 2 # eq 123
  P2 <- (-1 * (F + A) - ((F + A) ^ 2 - 4 * F * A + 4 * F * B) ^ (1 / 2)) / 2 # eq 123
  D <- (Qabs + Trad / Rrad + Tair / Rconv) / ((Rs / (2 * Rrad) + Rs / (2 * Rconv) + 1) * Cs) #eq 109 part 4
  alpha <- (Tco + (A + D) / (B - A)) # from eq 132
  beta <- F * Tso + (H * B - F * D) / (B - A) # from eq 133
  C1 <- (alpha * (P2 + F) - beta) / (P2 - P1) # from eq 134
  C2 <- (beta - alpha * (P1 + F)) / (P2 - P1) # from eq 135
  Tcf <- -1 * (F * D + A * H) / (F * (B - A)) # from eq 128
  Tc <- C1 * exp(P1 * t) + C2 * exp(P2 * t) + Tcf # from eq 128
  dTc <- (Qgen - Qresp - (Tc - Ts) / Rb) / Cc # from eq 103
  dTs <- ((Tc - Ts) / Rb - (Ts - Tsh) / (Rs / 2)) / Cs # from eq 104
  
  return(list(Tc=Tc,Tcf=Tcf,Tsh=Tsh,dTc=dTc,dTs=dTs))
}
