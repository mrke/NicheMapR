#' micro subroutine from microclimate model to compute wind speed, air temperature and relative humidity profiles
#'
#' Function for computing wind speed, air temperature and humidity profiles at a range of heights
#' beyond what comes out of the microclimate model, using the same functions that
#' the microclimate model uses in the micro.f subroutine. It is useful for situations where
#' the organism grows through time and thus the local height needs to change, or when simulating
#' a multi-part organism where the parts are at different heights (e.g., leg vs. head of a human)
#' @param reference_height = 1.2, Reference height (m), reference height at which air temperature, wind speed and relative humidity input data are measured (must match the micro_* function you are using, e.g. 1.2 for micro_global, 2 for micro_era5)
#' @param roughness_height = 0.004, Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001 (snow) - 0.02 m. (match to the value used in the original microclimate simulation)
#' @param canopy_roughness_height = 0, heat transfer roughness height (m) for Campbell and Norman air temperature/wind speed profile (invoked if greater than 0, 0.02 * canopy height in m if unknown
#' @param zero_plane_displacement = 0, zero plane displacement correction factor (m) for Campbell and Norman air temperature/wind speed profile (0.6 * canopy height in m if unknown)
#' @param air_temperature_reference = 27.8, air temperature (deg C) at reference height
#' @param wind_speed_reference = 2.75, wind speed (m/s) at reference height
#' @param relative_humidity_reference = 49.0415, relative humidity (pct) at reference height
#' @param surface_temperature = 48.6, soil surface temperature (deg C)
#' @param maximum_surface_temperature = 95, maximum allowed soil surface temperature - this is the default value in all micro_* functions
#' @param zenith_angle = 21.5, zenith angle (degrees) of sun - used in determining if free convection conditions or not
#' @param wind_shear_exponent = 0.15, wind shear exponent for extending above reference height (open water 0.1, Smooth, level, grass-covered 0.15 (or more commonly 1/7), row crops 0.2, low bushes with a few trees 0.2, heavy trees or several buildings or mountainous terrain 0.25, (source http://www.engineeringtoolbox.com/wind-shear-d_1215.html)
#' @param heights = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1), vector of heights (m) for which the profile is desired (make them between zero and the reference height, don't include the reference height, and make the minimum greater than the roughness height
#' @param warn = TRUE, show warning messages?
#' @return heights Heights (m) at which results are reported, with 0 and reference_height added on
#' @return VELs Wind speeds (m/s) with increasing height
#' @return TAs Air temperatures (deg C) with increasing height
#' @return RHs Relative humidity (pct) with increasing height
#' @return QCONV Convective heat loss (W/m2) from surface
#' @return USTAR Friction velocity (m/s)
#' @usage get_profile(canopy_roughness_height = 0.004, zero_plane_displacement = 0.012)
#' @examples
#' library(NicheMapR)
#' roughness_height <- 0.004 # choose a roughness height
#' micro <- micro_global(roughness_height = roughness_height) # run with defaults other than roughness height (Madison, Wisconsin)
#' dates <- micro$dates # extract mock dates (units of months)
#' micromet_lowshade <- as.data.frame(micro$micromet_lowshade) # above ground min shade conditions
#' soil_temperature_lowshade <- as.data.frame(micro$soil_temperature_lowshade) # below ground min shade conditions
#' newheights <- c(0.1, 0.4) # m, new height needed (can be a single value or a vector of heights)
#' profile.out <- lapply(1:length(micromet_lowshade$air_temperature_local),
#'                      function(x){get_profile(reference_height = 1.2, # needs to be what micro_global uses as reference_height
#'                                              roughness_height = roughness_height,
#'                                              heights = newheights,
#'                                              air_temperature_reference = micromet_lowshade$air_temperature_reference[x],
#'                                              wind_speed_reference = micromet_lowshade$wind_speed_reference[x],
#'                                              relative_humidity_reference = micromet_lowshade$relative_humidity_reference[x],
#'                                              surface_temperature = soil_temperature_lowshade$depth_0cm[x],
#'                                              zenith_angle = micromet_lowshade$zenith_angle[x])}) # run get_profile across all times at new height
#' profile.out1 <- do.call("rbind", lapply(profile.out, data.frame)) # turn results into data frame
#' newheight.out <- subset(profile.out1, heights == newheights[2])
#' plot(dates, micromet_lowshade$air_temperature_local, ylab = 'air temperature, deg C', type = 'l')
#' points(dates, profile.out1$TAs[profile.out1$heights == newheights[1]], type = 'l', lty = 2)
#' points(dates, profile.out1$TAs[profile.out1$heights == newheights[2]], type = 'l', lty = 3)
#' plot(dates, micromet_lowshade$wind_speed_local, ylab = 'wind speed, m/s', type = 'l', ylim = c(0, 3))
#' points(dates, profile.out1$VELs[profile.out1$heights == newheights[1]], type = 'l', lty = 2)
#' points(dates, profile.out1$VELs[profile.out1$heights == newheights[2]], type = 'l', lty = 3)
#' plot(dates, micromet_lowshade$relative_humidity_local, ylab = 'relative humidity, pct', type = 'l', ylim = c(0, 100))
#' points(dates, profile.out1$RHs[profile.out1$heights == newheights[1]], type = 'l', lty = 2)
#' points(dates, profile.out1$RHs[profile.out1$heights == newheights[2]], type = 'l', lty = 3)
#'@export
get_profile <- function(reference_height = 1.2,
                        roughness_height = 0.004,
                        canopy_roughness_height = 0,
                        zero_plane_displacement = 0,
                        air_temperature_reference = 27.77818,
                        wind_speed_reference = 2.749575,
                        relative_humidity_reference = 49.0415,
                        surface_temperature = 48.58942,
                        maximum_surface_temperature = 95,
                        zenith_angle = 21.50564,
                        wind_shear_exponent = 0.15,
                        heights = c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1),
                        warn = FALSE) {
  errors <- 0
  addheight <- FALSE
  if (min(heights) < roughness_height) {
    message(
      "ERROR: the minimum height is not greater than the roughness height, 'roughness_height'.
        Please enter a correct value.",
      '\n'
    )
    errors <- 1
  }
  if (min(heights) > reference_height) {
    addheight <- TRUE
    heights <- c(0.01, heights)
  }
  if (max(heights) >= reference_height) {
    if (warn) {
      message(
        "warning: there are heights greater than or equal to the reference height, 'reference_height'.
        Assuming constant air temperature above the reference height and adjusting wind speed according to shear parameter.",
        '\n'
      )
    }
    heights_orig <- heights
    heights <- heights_orig[heights_orig < reference_height]
    heights_extra <- heights_orig[heights_orig > reference_height]
  }
  if (errors == 0) {
    # continue
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

    get_Obukhov <- function(TA, TS, VEL, z, z0) {
      AMOL <- -30 #! initial Monin-Obukhov length cm
      GAM <- 16
      kappa <- 0.4 # Kármán constant
      RCPTKG <- 6.003e-8 # RHO*CP*T/(K*G) = 6.003D-8 IN CAL-MIN-CM-C UNITS
      Z <- z * 100
      Z0 <- z0 * 100
      ZRATIO <- Z / Z0 + 1 # ratio of reference to roughness height
      DUM <- log(ZRATIO)
      PHI <- function(Z, GAM, AMOL) {
        (1 - min(1, GAM * Z / AMOL))^0.25
      }
      PSI1 <- function(X) {
        2 * log((1 + X) / 2) + log((1 + X^2) / 2) - 2 * atan(X) + pi / 2
      }
      PSI2 <- function(X) {
        2 * log((1 + X^2) / 2)
      }
      DIFFT <- TA - TS # temp at reference height minus ground temp
      TAVE <- (TA + TS + 546.3) / 2 # ave temp in Kelvin
      RCP <- 0.08472 / TAVE
      count <- 0
      DEL <- 1
      while (DEL > 1e-2 & count < 500 & !is.na(AMOL)) {
        count <- count + 1
        # ITERATING TO FIND THE MONIN-OBUKHOV LENGTH (AMOL)
        X <- PHI(Z, GAM, AMOL)
        Y <- PSI1(X)
        YY <- PSI2(X)
        USTAR <- kappa * (VEL * 100 * 60) / (log(Z / Z0) - Y)
        if (AMOL > 0) {
          STS <- 0.62 / (Z0 * USTAR / 12)^0.45 # SUBLAYER STANTON NO.
          STB <- 0.64 / DUM # BULK STANTON NO.
          QC <- RCP * DIFFT * USTAR * STB / (1 + STB / STS) # convective heat transfer at the surface
        } else{
          STS <- 0.62 / (Z0 * USTAR / 12)^0.45
          STB <- (0.64 / DUM) * (1 - 0.1 * Z / AMOL) # BULK STANTON NO.
          STO <- STB / (1 + STB / STS)
          QC <- RCP * DIFFT * USTAR  * STO
        }
        AMOLN <- RCPTKG * USTAR^3 / QC
        DEL <- abs((AMOLN - AMOL) / AMOL)
        AMOL <- AMOLN
      }
      return(list(
        AMOL = AMOL / 100,
        STS = STS,
        STO = STO,
        STB = STB,
        USTAR = USTAR,
        QC = QC
      )) # convert back to m
    }

    T1 <- air_temperature_reference
    T3 <- surface_temperature
    # unit conversions
    Z <- reference_height * 100 # m to cm
    Z0 <- roughness_height * 100 # m to cm
    ZH_cm <- canopy_roughness_height * 100 # m to cm
    D0_cm <- zero_plane_displacement * 100 # m to cm
    V <- wind_speed_reference * 100 * 60 # m/s to cm/min
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
    RHOCP <- function(TAVE) {
      0.08472 / TAVE
    }
    PHI <- function(Z) {
      (1 - GAM * Z / AMOL)^0.25
    }
    PSI1 <- function(X) {
      2 * log((1 + X) / 2) + log((1 + X * X) / 2) - 2 * atan(X) + pi / 2
    }
    PSI2 <- function(X) {
      2 * log((1 + X * X) / 2)
    }

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

    STS <- 0.62 / (Z0 * USTAR / 12)^0.45 #SUBLAYER STANTON NO.
    STB <- 0.64 / DUM # BULK STANTON NO.
    QC <- RCP * DIFFT * USTAR * STB / (1 + STB / STS) # convective heat transfer at the surface

    # alternative Campbell and Norman 1998 vertical air temperature profile calculation option
    if (canopy_roughness_height > 0) {
      for (i in 2:NAIR) {
        A <- (T1 - T3) / (1 - log((Z - D0_cm) / ZH_cm))
        T0 <- T1 + A * log((Z - D0_cm) / ZH_cm)
        T[i] <- T0 - A * log((ZZ[i] - D0_cm) / ZH_cm)
      }
    }
    if (T1 >= T3 | zenith_angle >= 90) {
      STS <- 0.62 / (Z0 * USTAR / 12)^0.45 #SUBLAYER STANTON NO.
      STB <- 0.64 / DUM # BULK STANTON NO.
      QC <- RCP * DIFFT * USTAR * STB / (1 + STB / STS) # convective heat transfer at the surface
      for (i in 2:NAIR) {
        # FILL OUT VEL. AND TEMP. PROFILES
        VV[i] <- (USTAR / kappa) * log(ZZ[i] / Z0 + 1)
        # COMPUTING FICTITIOUS TEMP. AT TOP OF SUBLAYER
        if (canopy_roughness_height == 0) {
          TZO <- (T1 * STB + T3 * STS) / (STB + STS)
          T[i] <- TZO + (T1 - TZO) * log(ZZ[i] / Z0 + 1) / DUM
        }
      }
      # CHECK FOR FREE CONVECTION (LAPSE) CONDITIONS
    } else{
      for (i in 2:NAIR) {
        X1 <- PHI(ZZ[i])
        Y1 <- PSI1(X1)
        YY2 <- PSI2(X1)
        X <- PHI(Z)
        Y <- PSI1(X)
        YY <- PSI2(X)
        # FILL OUT VELOCITY AND TEMP. PROFILES
        ADUM <- ZZ[i] / Z0 - Y1
        VV[i] <- (USTAR / kappa) * log(ADUM)
        Obukhov.out <- get_Obukhov(T1, T3, V / 100 / 60, ZZ[i] / 100, Z0 / 100)
        if (canopy_roughness_height == 0) {
          TZO <- (T1 * Obukhov.out$STB + T3 * Obukhov.out$STS) / (Obukhov.out$STB + Obukhov.out$STS)
          T[i] <- TZO + (T1 - TZO) * log(ZZ[i] / Z0 - YY2) / log(Z / Z0 -
                                                                   YY)
        }
      }
    }

    # add zero, ref height and reorder
    heights <- c(0, heights, reference_height)
    VV <- c(0, rev(VV))
    T <- c(T3, rev(T))
    e <- WETAIR(db = T1, rh = relative_humidity_reference)$e
    es <- WETAIR(db = T)$esat
    RHs <- (e / es) * 100
    RHs[RHs > 100] <- 100
    RHs[RHs < 0] <- 0
    # add heights above reference height if any
    if (exists("heights_extra")) {
      VV_extra <- V * (heights_extra / reference_height)^wind_shear_exponent
      T_extra <- VV_extra * 0 + air_temperature_reference
      RH_extra <- VV_extra * 0 + relative_humidity_reference
      heights <- c(heights, heights_extra)
      VV <- c(VV, VV_extra)
      T <- c(T, T_extra)
      RHs <- c(RHs, RH_extra)
    }
    if (addheight) {
      VV <- VV[-2]
      T <- T[-2]
      RHs <- RHs[-2]
      heights <- heights[-2]
    }
    return(
      list(
        heights = heights,
        VELs = VV / 6000,
        # m/s
        TAs = T,
        # deg C
        RHs = RHs,
        # humidity, pct
        QCONV = QC * 4.185 * 10000 / 60,
        # W
        USTAR = USTAR / 6000 # m/s
      )
    )
  }
}
