#' micro subroutine from microclimate model to compute wind speed, air temperature and relative humidity profiles
#'
#' Function for computing wind speed, air temperature and humidity profiles at a range of heights
#' beyond what comes out of the microclimate model, using the same functions that
#' the microclimate model uses in the micro.f subroutine. It is useful for situations where
#' the organism grows through time and thus the local height needs to change, or when simulating
#' a multi-part organism where the parts are at different heights (e.g., leg vs. head of a human)
#' @param Refhyt = 1.2, Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured (must match the micro_* function you are using, e.g. 1.2 for micro_global, 2 for micro_era5)
#' @param RUF = 0.004, Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001 (snow) - 0.02 m. (match to the value used in the original microclimate simulation)
#' @param ZH = 0, heat transfer roughness height (m) for Campbell and Norman air temperature/wind speed profile (invoked if greater than 0, 0.02 * canopy height in m if unknown
#' @param D0 = 0, zero plane displacement correction factor (m) for Campbell and Norman air temperature/wind speed profile (0.6 * canopy height in m if unknown)
#' @param TAREF = 27.8, air temperature (deg C) at reference height
#' @param VREF = 2.75, wind speed (m/s) at reference height
#' @param RH = 49.0415, relative humidity (%) at reference height
#' @param D0cm = 48.6, soil surface temperature (deg C)
#' @param maxsurf = 95, maximum allowed soil surface temperature - this is the default value in all micro_* functions
#' @param ZEN = 21.5, zenith angle (degrees) of sun - used in determining if free convection conditions or not
#' @param heights = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1), vector of heights (m) for which the profile is desired (make them between zero and the reference height, don't include the reference height, and make the minimum greater than the roughness height
#' @return heights Heights (m) at which results are reported, with 0 and Refhyt added on
#' @return VELs Wind speeds (m/s) with increasing height
#' @return TAs Air temperatures (deg C) with increasing height
#' @return RHs Relative humidity (%) with increasing height
#' @return QCONV Convective heat loss (W/m2) from surface
#' @return USTAR Friction velocity (m/s)
#' @usage get_profile(ZH = 0.004, D0 = 0.012)
#' @examples
#' library(NicheMapR)
#' RUF <- 0.004 # choose a roughness height
#' micro <- micro_global(RUF = RUF) # run with defaults other than roughness height (Madison, Wisconsin)
#' dates <- micro$dates # extract mock dates (units of months)
#' metout <- as.data.frame(micro$metout) # above ground min shade conditions
#' soil <- as.data.frame(micro$soil) # below ground min shade conditions
#' newheights <- c(0.1, 0.4) # m, new height needed (can be a single value or a vector of heights)
#' profile.out <- lapply(1:length(metout$TALOC),
#'                      function(x){get_profile(Refhyt = 1.2, # needs to be what micro_global uses as Refhyt
#'                                              RUF = RUF,
#'                                              heights = newheights,
#'                                              TAREF = metout$TAREF[x],
#'                                              VREF = metout$VREF[x],
#'                                              RH = metout$RH[x],
#'                                              D0cm = soil$D0cm[x],
#'                                              ZEN = metout$ZEN[x])}) # run get_profile across all times at new height
#' profile.out1 <- do.call("rbind", lapply(profile.out, data.frame)) # turn results into data frame
#' newheight.out <- subset(profile.out1, heights == newheights[2])
#' plot(dates, metout$TALOC, ylab = 'air temperature, deg C', type = 'l')
#' points(dates, profile.out1$TAs[profile.out1$heights == newheights[1]], type = 'l', lty = 2)
#' points(dates, profile.out1$TAs[profile.out1$heights == newheights[2]], type = 'l', lty = 3)
#' plot(dates, metout$VLOC, ylab = 'wind speed, m/s', type = 'l', ylim = c(0, 3))
#' points(dates, profile.out1$VELs[profile.out1$heights == newheights[1]], type = 'l', lty = 2)
#' points(dates, profile.out1$VELs[profile.out1$heights == newheights[2]], type = 'l', lty = 3)
#' plot(dates, metout$RHLOC, ylab = 'relative humidity, %', type = 'l', ylim = c(0, 100))
#' points(dates, profile.out1$RHs[profile.out1$heights == newheights[1]], type = 'l', lty = 2)
#' points(dates, profile.out1$RHs[profile.out1$heights == newheights[2]], type = 'l', lty = 3)
#'@export
get_profile <- function(Refhyt = 1.2,
                        RUF = 0.004,
                        ZH = 0, #0.002
                        D0 = 0, #0.06
                        TAREF = 27.77818,
                        VREF = 2.749575,
                        RH = 49.0415,
                        D0cm = 48.58942,
                        maxsurf = 95,
                        ZEN = 21.50564,
                        heights = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1)){
  errors <- 0
  if(min(heights) < RUF){
    message("ERROR: the minimum height is not greater than the roughness height, 'RUF'.
        Please enter a correct value.", '\n')
    errors <- 1
  }
  if(max(heights) >= Refhyt){
    message("warning: there are heights greater than or equal to the reference height, 'Refhyt'.
        Removing these values.", '\n')
    heights <- heights[heights < Refhyt]
  }

  # 1 SEGMENT VELOCITY PROFILE - W. PORTER
  # VELOCITY PROFILE - Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971). Flux-Profile Relationships in the Atmospheric Surface Layer. Journal of the Atmospheric Sciences, 28(2), 181–189. doi:10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2
  # SUBLAYER MODEL - Garratt, J. R., & Hicks, B. B. (1973). Momentum, heat and water vapour transfer to and from natural and artificial surfaces. Quarterly Journal of the Royal Meteorological Society, 99(422), 680–687. doi:10.1002/qj.49709942209
  # Z=REFERENCE HEIGHT
  # Z0=ROUGHNESS HEIGHT
  # T1=TEMPERATURE AT REFERENCE HEIGHT
  # T3=CURRENT ESTIMATE OF GROUND SURFACE TEMP.
  # V=VELOCITY AT REF. HEIGHT
  # QC=COMPUTED (HERE) CONVECTIVE HEAT TRANSFER AT THE SURFACE
  # AMOL=MONIN-OBUKHOV LENGTH
  # NAIR=NO. OF HEIGHTS FOR AIR TEMP'S AND VELOCITIES
  # ZZ=ARRAY OF HEIGHT VALUES
  # VV=ARRAY OF COMPUTED (HERE) VELOCITIES FOR EACH HEIGHT
  # T=ARRAY OF AIR TEMP'S COMPUTED HERE ASSUMING A LOG PROFILE
  #
  # THIS SUBROUTINE IS MODIFIED (FEB. 1979)FOR SHEAR OCCURRING ABOVE
  # THE SURFACE DUE TO VEGETATION SPACED OVER THE SURFACE.
  # TEMP. PROFILE REMAINS LOGARITHMIC. VEL. PROFILE LOGARITHMIC IN SEGMENTS

  get_Obukhov <- function(TA, TS, VEL, z, z0){
    AMOL <- -30 #! initial Monin-Obukhov length cm
    GAM <- 16
    kappa <- 0.4 # Kármán constant
    RCPTKG <- 6.003e-8 # RHO*CP*T/(K*G) = 6.003D-8 IN CAL-MIN-CM-C UNITS
    Z <- z * 100
    Z0 <- z0 * 100
    ZRATIO <- Z / Z0 + 1 # ratio of reference to roughness height
    DUM <- log(ZRATIO)
    PHI <- function(Z, GAM, AMOL){
      (1-min(1, GAM * Z / AMOL)) ^ 0.25
    }
    PSI1 <- function(X){
      2 * log((1 + X) / 2) + log((1 + X ^ 2) / 2) - 2 * atan(X) + pi / 2
    }
    PSI2 <- function(X){
      2 * log((1 + X ^ 2) / 2)
    }
    DIFFT <- TA - TS # temp at reference height minus ground temp
    TAVE <- (TA + TS + 546.3) / 2 # ave temp in Kelvin
    RCP <- 0.08472 / TAVE
    count <- 0
    DEL <- 1
    while(DEL > 1e-2 & count < 500 & !is.na(AMOL)){
      count <- count + 1
      # ITERATING TO FIND THE MONIN-OBUKHOV LENGTH (AMOL)
      X <- PHI(Z, GAM, AMOL)
      Y <- PSI1(X)
      YY <- PSI2(X)
      USTAR <- kappa * (VEL * 100 * 60) / (log(Z / Z0) - Y)
      if(AMOL > 0){
        STS <- 0.62 / (Z0 * USTAR / 12) ^ 0.45 # SUBLAYER STANTON NO.
        STB <- 0.64 / DUM # BULK STANTON NO.
        QC <- RCP * DIFFT * USTAR * STB / (1 + STB / STS) # convective heat transfer at the surface
      }else{
        STS <- 0.62 / (Z0 * USTAR / 12) ^ 0.45
        STB <- (0.64 / DUM) * (1- 0.1 * Z / AMOL) # BULK STANTON NO.
        STO <- STB / (1 + STB / STS)
        QC <- RCP * DIFFT * USTAR  *STO
      }
      AMOLN <- RCPTKG * USTAR ^ 3/ QC
      DEL <- abs((AMOLN - AMOL) / AMOL)
      AMOL <- AMOLN
    }
    return(list(AMOL = AMOL/100, STS = STS, STO = STO, STB = STB, USTAR = USTAR, QC = QC)) # convert back to m
  }

  T1 <- TAREF
  T3 <- D0cm
  # unit conversions
  Z <- Refhyt * 100 # m to cm
  Z0 <- RUF * 100 # m to cm
  ZH_cm <- ZH * 100 # m to cm
  D0_cm <- D0 * 100 # m to cm
  V <- VREF * 100 * 60 # m/s to cm/min
  # DEFINE AIR HEIGHTS (cm)
  AIRDP <- c(Z, rev(heights) * 100)
  ZZ <- AIRDP
  NAIR <- length(AIRDP)
  VV <- rep(0, NAIR) # output wind speeds
  T <- VV # output temperatures
  RHs <- VV # output relative humidities
  VV[1] <- V
  T[1] <- T1

  # some necessary functions
  RHOCP <- function(TAVE){0.08472 / TAVE}
  PHI <- function(Z){(1 - GAM * Z / AMOL) ^ 0.25}
  PSI1 <- function(X){2 * log((1 + X) / 2) + log((1 + X * X) / 2) - 2 * atan(X) + pi / 2}
  PSI2 <- function(X){2 * log((1 +X * X) / 2)}

  GAM <- 16
  RCPTKG <- 6.003e-8 # RHO * CP * T / (K * G) = 6.003D-8 IN CAL-MIN-CM-C UNITS
  kappa <- 0.4 # Kármán constant

  # COMPUTING VEL. PROFILE PARAMETERS FROM 200 CM REFERENCE VELOCITY
  ZRATIO <- Z / Z0 + 1 # ratio of reference to roughness height
  DUM <- log(ZRATIO)
  USTAR <- kappa * V / DUM # friction velocity
  DIFFT <- T1 - T3 # temp at reference height minus ground temp
  TAVE <- (T3 + T1 + 546) / 2 # ave temp in Kelvin
  RCP <- RHOCP(TAVE)
  AMOL <- -30 # initial Monin-Obukhov length (unstable conditions because negative)
  ITER <- 0 # initialise counter

  # alternative Campbell and Norman 1998 vertical air temperature profile calculation option
  if(ZH > 0){
    STS <- 0.62 / (Z0 * USTAR / 12) ^ 0.45 #SUBLAYER STANTON NO.
    STB <- 0.64 / DUM # BULK STANTON NO.
    QC <- RCP * DIFFT * USTAR * STB / (1 + STB / STS) # convective heat transfer at the surface
    # Use vertical temperature profile from Campbell and Norman 1998
    if(NAIR <= 0){break}
    for(i in 2:NAIR) {
      # FILL OUT VEL. AND TEMP. PROFILES
      if(T1 >= T3 | T3 <= maxsurf | ZEN >= 90){
        VV[i] <- (USTAR / kappa) * log(ZZ[i] / Z0 + 1)
      }else{
        X1 <- PHI(ZZ[i])
        Y1 <- PSI1(X1)
        ADUM <- ZZ[i] / Z0 - Y1
        VV[i] <- (USTAR / kappa) * log(ADUM)
      }
      A <- (T1 - T3)/(1 - log((Z - D0_cm) / ZH_cm))
      T0 <- T1 + A * log((Z - D0_cm) / ZH_cm)
      T[i] <- T0 - A * log((ZZ[i] - D0_cm) / ZH_cm)
    }
  }else{
    if(T1 >= T3 | T3 <= maxsurf | ZEN >= 90){
      STS <- 0.62 / (Z0 * USTAR / 12) ^ 0.45 #SUBLAYER STANTON NO.
      STB <- 0.64 / DUM # BULK STANTON NO.
      QC <- RCP * DIFFT * USTAR * STB / (1 + STB / STS) # convective heat transfer at the surface
      if(NAIR <= 0){break}
      for(i in 2:NAIR) {
        # FILL OUT VEL. AND TEMP. PROFILES
        VV[i] <- (USTAR / kappa) * log(ZZ[i] / Z0 + 1)
        # COMPUTING FICTITIOUS TEMP. AT TOP OF SUBLAYER
        TZO <- (T1 * STB + T3 * STS) / (STB + STS)
        T[i] <- TZO + (T1 - TZO) * log(ZZ[i] / Z0 + 1) / DUM
      }
      # CHECK FOR FREE CONVECTION (LAPSE) CONDITIONS
    }else{
      if(NAIR <= 0){break}
      for(i in 2:NAIR) {
        X1 <- PHI(ZZ[i])
        Y1 <- PSI1(X1)
        YY2 <- PSI2(X1)
        X <- PHI(Z)
        Y <- PSI1(X)
        YY <- PSI2(X)
        # FILL OUT VELOCITY AND TEMP. PROFILES
        ADUM <- ZZ[i]/Z0-Y1
        VV[i] <- (USTAR / kappa)*log(ADUM)
        Obukhov.out <- get_Obukhov(T1, T3, V / 100 / 60, ZZ[i] / 100, Z0 / 100)
        TZO <- (T1 * Obukhov.out$STB + T3 * Obukhov.out$STS) / (Obukhov.out$STB + Obukhov.out$STS)
        T[i] <- TZO + (T1 - TZO) * log(ZZ[i] / Z0 - YY2) / log(Z / Z0-YY)
      }
    }
  }
  # add zero, ref height and reorder
  heights <- c(0, heights, Refhyt)
  VV <- c(0, rev(VV))
  T <- c(T3, rev(T))
  e <- WETAIR(db = T1, rh = RH)$e
  es <- WETAIR(db = T)$esat
  RHs <- (e / es) * 100
  RHs[RHs > 100] <- 100
  RHs[RHs < 0] <- 0
  return(list(heights = heights,
              VELs = VV / 6000, # m/s
              TAs = T, # deg C
              RHs = RHs, # humidity, %
              QCONV = QC * 4.185 * 10000 / 60, # W
              USTAR = USTAR / 6000 # m/s
  )
  )
}
