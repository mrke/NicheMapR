# Shared helper functions used by all micro_*.R microclimate model drivers.
#
# These internal helpers encapsulate boilerplate that is common across the 11
# climate-dataset-specific wrappers (micro_global, micro_era5, micro_ncep,
# micro_openmeteo, micro_aust, micro_nz, micro_silo, micro_uk, micro_usa,
# micro_terra, micro_access_s2).  They are NOT exported; users interact only
# with the per-dataset functions.

# ---------------------------------------------------------------------------
# 1.  Input validation
# ---------------------------------------------------------------------------

#' Validate common microclimate model inputs
#'
#' Checks the parameters that are shared by every \code{micro_*()} function
#' and returns an error flag.  File-specific checks (e.g. \code{timeinterval},
#' \code{soiltype}) are handled inside each calling function.
#'
#' @param DEP         Numeric vector (length 10) of soil node depths (cm).
#' @param REFL        Soil solar reflectance (0–1).
#' @param slope       Slope in degrees (0–90).
#' @param aspect      Aspect in degrees (0–365).
#' @param hori        Numeric vector (length 24) of horizon angles (degrees).
#' @param SLE         Substrate longwave emissivity (0.05–1).
#' @param ERR         Integrator error tolerance (must be > 0).
#' @param RUF         Roughness height, m (0.0001–2).
#' @param D0          Zero-plane displacement height, m.
#' @param Usrhyt      Local height for met output, m (must exceed RUF and be
#'                    no greater than Refhyt).
#' @param Refhyt      Reference height of input met data, m.
#' @param EC          Orbital eccentricity (0.0034–0.058).
#' @param CMH2O       Precipitable water in air column, cm (0.1–2).
#' @param TIMAXS      Numeric vector (length 4): time offsets for weather maxima.
#' @param TIMINS      Numeric vector (length 4): time offsets for weather minima.
#' @param minshade    Minimum shade level(s), \% (0–100).
#' @param maxshade    Maximum shade level(s), \% (0–100, must exceed minshade).
#' @param Thcond      Soil mineral thermal conductivity, W/mK (must be >= 0).
#' @param Density     Soil mineral density, Mg/m3 (must be >= 0).
#' @param SpecHeat    Soil mineral specific heat, J/kg-K (must be >= 0).
#' @param BulkDensity Soil bulk density, Mg/m3 (must be >= 0).
#'
#' @return Integer: \code{0} if all checks pass, \code{1} if any check fails.
#'   Error messages are emitted via \code{message()}.
#'
#' @keywords internal
validate_micro_inputs <- function(DEP, REFL, slope, aspect, hori, SLE, ERR,
                                  RUF, D0, Usrhyt, Refhyt, EC, CMH2O,
                                  TIMAXS, TIMINS, minshade, maxshade,
                                  Thcond, Density, SpecHeat, BulkDensity) {
  errors <- 0

  # --- Soil depth node spacing (warnings only, do not set errors) -----------
  if (DEP[2] - DEP[1] > 3 | DEP[3] - DEP[2] > 3) {
    message("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing")
  }
  if (DEP[2] - DEP[1] < 2) {
    message("warning, nodes might be too close near the surface, try a different spacing if the program is crashing")
  }
  if (DEP[10] != 200) {
    message("warning, last depth in soil should not be changed from 200 without good reason")
  }

  # --- DEP vector structure -------------------------------------------------
  if (DEP[1] != 0) {
    message("ERROR: First soil node (DEP[1]) must = 0 cm. Please correct.")
    errors <- 1
  }
  if (length(DEP) != 10) {
    message("ERROR: You must enter 10 different soil depths.")
    errors <- 1
  }
  for (i in 1:9) {
    if (DEP[i + 1] <= DEP[i]) {
      message("ERROR: Soil depth (DEP array) is not in ascending order.")
      errors <- 1
    }
  }
  if (DEP[10] > 500) {
    message("ERROR: Deepest soil depth (DEP array) is too large (must be <= 500 cm).")
    errors <- 1
  }

  # --- Roughness and displacement heights -----------------------------------
  if (RUF < 0.0001) {
    message("ERROR: The roughness height (RUF) is too small (< 0.0001). Please enter a larger value.")
    errors <- 1
  }
  if (RUF > 2) {
    message("ERROR: The roughness height (RUF) is too large (> 2). Please enter a smaller value.")
    errors <- 1
  }
  if (D0 > 0 & D0 < Usrhyt) {
    message("ERROR: The zero plane displacement height (D0) must be lower than the local height (Usrhyt). Please enter a smaller value.")
    errors <- 1
  }
  if (Usrhyt < RUF) {
    message("ERROR: Local height (Usrhyt) smaller than roughness height (RUF). Please use a larger height above the surface.")
    errors <- 1
  }
  if (Usrhyt > Refhyt) {
    message("ERROR: Reference height is less than local height (Usrhyt).")
    errors <- 1
  }

  # --- Soil physical properties ---------------------------------------------
  if (min(Thcond) < 0) {
    message("ERROR: Thermal conductivity (Thcond) is negative. Please input a positive value.")
    errors <- 1
  }
  if (min(Density) < 0) {
    message("ERROR: Soil mineral density (Density) is negative. Please input a positive value.")
    errors <- 1
  }
  if (min(SpecHeat) < 0) {
    message("ERROR: Specific heat (SpecHeat) is negative. Please input a positive value.")
    errors <- 1
  }
  if (min(BulkDensity) < 0) {
    message("ERROR: Bulk density (BulkDensity) is negative. Please input a positive value.")
    errors <- 1
  }

  # --- Radiation and site geometry ------------------------------------------
  if (REFL < 0 | REFL > 1) {
    message("ERROR: Soil reflectivity (REFL) is out of bounds. Please input a value between 0 and 1.")
    errors <- 1
  }
  if (slope < 0 | slope > 90) {
    message("ERROR: Slope value is out of bounds. Please input a value between 0 and 90.")
    errors <- 1
  }
  if (aspect < 0 | aspect > 365) {
    message("ERROR: Aspect value is out of bounds. Please input a value between 0 and 365.")
    errors <- 1
  }
  if (max(hori) > 90 | min(hori) < 0) {
    message("ERROR: At least one horizon angle (hori) is out of bounds. Please input values between 0 and 90.")
    errors <- 1
  }
  if (length(hori) != 24) {
    message("ERROR: You must enter 24 horizon angle values.")
    errors <- 1
  }
  if (SLE < 0.05 | SLE > 1) {
    message("ERROR: Emissivity (SLE) is out of bounds. Please enter a value between 0.05 and 1.00.")
    errors <- 1
  }

  # --- Atmospheric and orbital parameters -----------------------------------
  if (ERR < 0) {
    message("ERROR: Error bound (ERR) is too small. Please enter a value > 0.")
    errors <- 1
  }
  if (EC < 0.0034 | EC > 0.058) {
    message("ERROR: Eccentricity (EC) is out of bounds. Please enter a value between 0.0034 and 0.058.")
    errors <- 1
  }
  if (CMH2O < 0.1 | CMH2O > 2) {
    message("ERROR: Precipitable water (CMH2O) is out of bounds. Please enter a value between 0.1 and 2.")
    errors <- 1
  }

  # --- Weather timing offsets -----------------------------------------------
  if (max(TIMAXS) > 24 | min(TIMAXS) < 0) {
    message("ERROR: At least one time of weather maximum (TIMAXS) is out of bounds. Please input values between 0 and 24.")
    errors <- 1
  }
  if (max(TIMINS) > 24 | min(TIMINS) < 0) {
    message("ERROR: At least one time of weather minimum (TIMINS) is out of bounds. Please input values between 0 and 24.")
    errors <- 1
  }

  # --- Shade levels ---------------------------------------------------------
  if (max(minshade - maxshade) >= 0) {
    message("ERROR: Minimum shade (minshade) must be less than maximum shade (maxshade). Please correct.")
    errors <- 1
  }
  if (max(minshade) > 100 | min(minshade) < 0) {
    message("ERROR: Minimum shade (minshade) is out of bounds. Please input a value between 0 and 100.")
    errors <- 1
  }
  if (max(maxshade) > 100 | min(maxshade) < 0) {
    message("ERROR: Maximum shade (maxshade) is out of bounds. Please input a value between 0 and 100.")
    errors <- 1
  }

  return(errors)
}


# ---------------------------------------------------------------------------
# 2.  Shade vector construction
# ---------------------------------------------------------------------------

#' Create daily shade fraction vectors
#'
#' Expands scalar or short \code{minshade}/\code{maxshade} inputs to a full
#' vector of length \code{ndays}, as required by the Fortran model.
#'
#' @param minshade Scalar or vector of minimum shade values (\%).
#' @param maxshade Scalar or vector of maximum shade values (\%).
#' @param ndays    Total number of simulation days.
#'
#' @return Named list with elements \code{MINSHADES} and \code{MAXSHADES},
#'   each a numeric vector of length \code{ndays}.
#'
#' @keywords internal
setup_shade_vectors <- function(minshade, maxshade, ndays) {
  if (length(minshade) != ndays) {
    MINSHADES <- rep(minshade[1], ndays)
  } else {
    MINSHADES <- minshade
  }
  if (length(maxshade) != ndays) {
    MAXSHADES <- rep(maxshade[1], ndays)
  } else {
    MAXSHADES <- maxshade
  }
  list(MINSHADES = MINSHADES, MAXSHADES = MAXSHADES)
}


# ---------------------------------------------------------------------------
# 3.  Timezone reference longitude (ALREF)
# ---------------------------------------------------------------------------

#' Compute the reference longitude for solar noon correction
#'
#' Returns \code{ALREF}, the longitude (degrees) used by the Fortran model to
#' correct solar time.  Two modes are supported:
#' \itemize{
#'   \item \code{timezone = 0} — simple truncation of the site longitude to
#'     the nearest whole degree (uses local apparent solar noon).
#'   \item \code{timezone = 1} — queries the GeoNames web service via the
#'     \pkg{geonames} package to obtain the civil UTC offset, then converts it
#'     to a reference longitude.  Requires a registered GeoNames account (see
#'     \code{?geonames::GNtimezone}).
#' }
#'
#' @param lon      Site longitude (decimal degrees, negative = west).
#' @param lat      Site latitude (decimal degrees, negative = south).
#' @param timezone Integer: \code{0} = local solar noon, \code{1} = civil
#'   timezone lookup via GeoNames.
#'
#' @return Numeric scalar: reference longitude for solar time correction (degrees).
#'
#' @keywords internal
get_timezone_alref <- function(lon, lat, timezone) {
  if (timezone == 1) {
    if (!requireNamespace("geonames", quietly = TRUE)) {
      stop('package "geonames" is required for timezone = 1. Please install it.')
    }
    # GNtimezone returns a data frame; column 4 is the UTC offset in hours.
    # Multiply by -15 to convert hours to the equivalent reference longitude.
    ALREF <- as.numeric(geonames::GNtimezone(lat, lon)[4]) * -15
  } else {
    # Use the nearest whole-degree longitude as the solar noon reference
    ALREF <- abs(trunc(lon))
  }
  return(ALREF)
}


# ---------------------------------------------------------------------------
# 4.  Microinput vector construction
# ---------------------------------------------------------------------------

#' Build the microinput parameter vector for the Fortran microclimate model
#'
#' Assembles the numeric vector consumed by the Fortran \code{microclimate()}
#' routine.  All \code{micro_*()} drivers use the same vector layout; the only
#' difference between the baseline dataset drivers and the ERA5/NCEP/OpenMeteo
#' drivers is the \code{dewrain} and \code{moiststep} parameters (defaulting to
#' \code{0} and \code{360} respectively for the baseline group).
#'
#' @param ndays     Total simulation days.
#' @param RUF       Roughness height (m).
#' @param ERR       Integrator error tolerance.
#' @param Usrhyt    Local height for met output (m).
#' @param Refhyt    Reference height of input met data (m).
#' @param Numtyps   Number of soil types.
#' @param Z01       Top segment roughness height (m).
#' @param Z02       Second segment roughness height (m).
#' @param ZH1       Top of first segment height (m).
#' @param ZH2       Second segment height (m).
#' @param idayst    Start day index.
#' @param ida       End day index.
#' @param HEMIS     Hemisphere flag (1 = north, 2 = south).
#' @param ALAT      Latitude degrees (absolute value, integer part).
#' @param AMINUT    Latitude minutes (decimal).
#' @param ALONG     Longitude degrees (absolute value, integer part).
#' @param ALMINT    Longitude minutes (decimal).
#' @param ALREF     Reference longitude for solar noon correction (degrees).
#' @param slope     Surface slope (degrees).
#' @param azmuth    Surface aspect/azimuth (degrees).
#' @param ALTT      Site elevation (m).
#' @param CMH2O     Precipitable water in air column (cm H2O).
#' @param microdaily  Flag: 1 = daily iteration mode (365-day runs), 0 = normal.
#' @param tannul    Mean annual air temperature used to initialise deep soil (°C).
#' @param EC        Orbital eccentricity.
#' @param VIEWF     Sky view factor (0–1).
#' @param snowtemp  Temperature threshold for snowfall (°C).
#' @param snowdens  Snow density (Mg/m3).
#' @param snowmelt  Proportion of calculated melt that does not refreeze.
#' @param undercatch Undercatch multiplier for rain-to-snow conversion.
#' @param rainmult  Rain multiplier for surface soil moisture.
#' @param runshade  Flag: 1 = run both shade levels, 0 = minimum shade only.
#' @param runmoist  Flag: 1 = run soil moisture model.
#' @param maxpool   Maximum water pooling depth (mm).
#' @param evenrain  Flag: 1 = spread rainfall evenly over 24 h, 0 = midnight event.
#' @param snowmodel Flag: 1 = run snow model.
#' @param rainmelt  Rain melt parameter.
#' @param writecsv  Flag: 1 = write Fortran output as CSV files.
#' @param densfun   Snow density function coefficients (vector of 4).
#' @param hourly    Flag: 1 = hourly input data provided.
#' @param rainhourly Flag: 1 = hourly rainfall provided.
#' @param lamb      Flag: 1 = return wavelength-specific solar output.
#' @param IUV       Flag: 1 = use gamma function for scattered solar radiation.
#' @param RW        Resistance per unit root length (m3 kg-1 s-1).
#' @param PC        Critical leaf water potential for stomatal closure (J kg-1).
#' @param RL        Resistance per unit leaf length (m3 kg-1 s-1).
#' @param SP        Stability parameter for stomatal closure.
#' @param R1        Root radius (m).
#' @param IM        Maximum allowable mass balance error (kg).
#' @param MAXCOUNT  Maximum iterations for mass balance.
#' @param IR        Longwave radiation formula: 0 = Campbell & Norman, 1 = Swinbank.
#' @param message   Flag: 1 = allow Fortran integrator warnings.
#' @param fail      Maximum integrator restarts before quitting.
#' @param snowcond  Effective snow thermal conductivity W/mC (0 = use density function).
#' @param intercept Snow interception fraction under shade (0–1).
#' @param grasshade Flag: 1 = remove shade when snow present (shade from grass/low shrubs).
#' @param solonly   Flag: 1 = only run SOLRAD solar radiation sub-model.
#' @param ZH        Heat transfer roughness height (m).
#' @param D0        Zero-plane displacement height (m).
#' @param TIMAXS    Vector of 4: time offsets for weather maxima (h).
#' @param TIMINS    Vector of 4: time offsets for weather minima (h).
#' @param spinup    Number of spin-up days.
#' @param maxsurf   Maximum surface temperature (°C).
#' @param ndmax     Number of day iterations for steady periodic solution.
#' @param dewrain   Dew/rain parameter (ERA5/NCEP/OpenMeteo only; default 0).
#' @param moiststep Moisture time step (ERA5/NCEP/OpenMeteo only; default 360).
#'
#' @return Numeric vector of microclimate model input parameters.
#'
#' @keywords internal
build_microinput <- function(ndays, RUF, ERR, Usrhyt, Refhyt, Numtyps,
                             Z01, Z02, ZH1, ZH2, idayst, ida,
                             HEMIS, ALAT, AMINUT, ALONG, ALMINT, ALREF,
                             slope, azmuth, ALTT, CMH2O, microdaily, tannul,
                             EC, VIEWF, snowtemp, snowdens, snowmelt,
                             undercatch, rainmult, runshade, runmoist,
                             maxpool, evenrain, snowmodel, rainmelt, writecsv,
                             densfun, hourly, rainhourly, lamb, IUV,
                             RW, PC, RL, SP, R1, IM, MAXCOUNT, IR,
                             message, fail, snowcond, intercept, grasshade,
                             solonly, ZH, D0, TIMAXS, TIMINS, spinup,
                             maxsurf, ndmax,
                             dewrain = 0, moiststep = 360) {
  c(ndays, RUF, ERR, Usrhyt, Refhyt, Numtyps, Z01, Z02, ZH1, ZH2,
    idayst, ida, HEMIS, ALAT, AMINUT, ALONG, ALMINT, ALREF,
    slope, azmuth, ALTT, CMH2O, microdaily, tannul,
    EC, VIEWF, snowtemp, snowdens, snowmelt, undercatch, rainmult,
    runshade, runmoist, maxpool, evenrain, snowmodel, rainmelt, writecsv,
    densfun, hourly, rainhourly, lamb, IUV, RW, PC, RL, SP, R1, IM,
    MAXCOUNT, IR, message, fail, snowcond, intercept, grasshade, solonly,
    ZH, D0, TIMAXS, TIMINS, spinup, dewrain, moiststep, maxsurf, ndmax)
}


# ---------------------------------------------------------------------------
# 5.  Micro input list construction
# ---------------------------------------------------------------------------

#' Assemble the micro input list for the microclimate() Fortran wrapper
#'
#' Builds the named list that is passed directly to \code{microclimate()}.
#' Every element corresponds to an array expected by the Fortran subroutine.
#'
#' @param microinput Numeric vector from \code{\link{build_microinput}}.
#' @param doy        Integer vector of day-of-year values for each simulation day.
#' @param SLES       Numeric vector of surface emissivity values (one per day).
#' @param DEP        Numeric vector of 10 soil node depths (cm).
#' @param Nodes      Matrix of node weights/depths for the Fortran solver.
#' @param MAXSHADES  Daily maximum shade fractions (\%).
#' @param MINSHADES  Daily minimum shade fractions (\%).
#' @param TMAXX      Daily maximum air temperatures (°C).
#' @param TMINN      Daily minimum air temperatures (°C).
#' @param RHMAXX     Daily maximum relative humidity (\%).
#' @param RHMINN     Daily minimum relative humidity (\%).
#' @param CCMAXX     Daily maximum cloud cover (oktas or fraction).
#' @param CCMINN     Daily minimum cloud cover.
#' @param WNMAXX     Daily maximum wind speed (m/s).
#' @param WNMINN     Daily minimum wind speed (m/s).
#' @param TAIRhr     Hourly air temperature override (°C); zeros if not used.
#' @param RHhr       Hourly relative humidity override (\%); zeros if not used.
#' @param WNhr       Hourly wind speed override (m/s); zeros if not used.
#' @param CLDhr      Hourly cloud cover override; zeros if not used.
#' @param SOLRhr     Hourly solar radiation override (W/m2); zeros if not used.
#' @param RAINhr     Hourly rainfall override (mm); zeros if not used.
#' @param ZENhr      Hourly zenith angle override (degrees); -1 if not used.
#' @param IRDhr      Hourly longwave radiation override (W/m2); -1 if not used.
#' @param REFLS      Daily soil reflectance values.
#' @param PCTWET     Daily fraction of surface acting as free water.
#' @param soilinit   Matrix of initial soil temperatures at each node.
#' @param hori       Horizon angles (24 values, degrees).
#' @param TAI        Solar attenuation vector (111 wavelength bins); 0 = use GADS.
#' @param soilprops  Matrix of soil thermal properties at each node.
#' @param moists     Matrix of initial soil moisture at each node.
#' @param RAINFALL   Daily rainfall (mm).
#' @param deepsoil   Annual mean temperature for deep soil boundary condition (°C).
#' @param PE         Air entry potential (J/kg) at each of 19 soil levels.
#' @param KS         Saturated conductivity (kg s/m3) at each of 19 soil levels.
#' @param BB         Campbell's soil b parameter at each of 19 soil levels.
#' @param BD         Bulk density (Mg/m3) at each of 19 soil levels.
#' @param DD         Mineral density (Mg/m3) at each of 19 soil levels.
#' @param L          Root density (m/m3) at each of 19 soil levels.
#' @param LAI        Leaf area index (one value per simulation day).
#' @param tides      Matrix (ndays*24 x 3) of tidal state, temperature and splash.
#'
#' @return Named list ready to pass to \code{microclimate()}.
#'
#' @keywords internal
build_micro_list <- function(microinput, doy, SLES, DEP, Nodes,
                             MAXSHADES, MINSHADES,
                             TMAXX, TMINN, RHMAXX, RHMINN,
                             CCMAXX, CCMINN, WNMAXX, WNMINN,
                             TAIRhr, RHhr, WNhr, CLDhr, SOLRhr,
                             RAINhr, ZENhr, IRDhr,
                             REFLS, PCTWET, soilinit, hori, TAI,
                             soilprops, moists, RAINFALL, deepsoil,
                             PE, KS, BB, BD, DD, L, LAI, tides) {
  list(
    tides      = tides,
    microinput = microinput,
    doy        = doy,
    SLES       = SLES,
    DEP        = DEP,
    Nodes      = Nodes,
    MAXSHADES  = MAXSHADES,
    MINSHADES  = MINSHADES,
    TMAXX      = TMAXX,
    TMINN      = TMINN,
    RHMAXX     = RHMAXX,
    RHMINN     = RHMINN,
    CCMAXX     = CCMAXX,
    CCMINN     = CCMINN,
    WNMAXX     = WNMAXX,
    WNMINN     = WNMINN,
    TAIRhr     = TAIRhr,
    RHhr       = RHhr,
    WNhr       = WNhr,
    CLDhr      = CLDhr,
    SOLRhr     = SOLRhr,
    RAINhr     = RAINhr,
    ZENhr      = ZENhr,
    IRDhr      = IRDhr,
    REFLS      = REFLS,
    PCTWET     = PCTWET,
    soilinit   = soilinit,
    hori       = hori,
    TAI        = TAI,
    soilprops  = soilprops,
    moists     = moists,
    RAINFALL   = RAINFALL,
    tannulrun  = deepsoil,   # Fortran expects 'tannulrun' as the list element name
    PE         = PE,
    KS         = KS,
    BB         = BB,
    BD         = BD,
    DD         = DD,
    L          = L,
    LAI        = LAI
  )
}


# ---------------------------------------------------------------------------
# 6.  Write inputs to CSV
# ---------------------------------------------------------------------------

#' Write all microclimate model inputs to CSV files
#'
#' Creates a folder \code{"micro csv input"} in the working directory and
#' writes every model input array to a separate CSV file.  Useful for
#' debugging or passing inputs to other tools.
#'
#' Invoked only when the calling function's \code{write_input} parameter is 1.
#'
#' @param microinput Numeric microinput vector (from \code{\link{build_microinput}}).
#' @param doy        Day-of-year vector.
#' @param SLES       Surface emissivity vector.
#' @param DEP        Soil depths vector (cm).
#' @param Nodes      Soil node weight matrix.
#' @param MAXSHADES  Daily maximum shade vector.
#' @param MINSHADES  Daily minimum shade vector.
#' @param TIMAXS     Time-of-maximum offsets vector.
#' @param TIMINS     Time-of-minimum offsets vector.
#' @param TMAXX      Daily maximum air temperature vector.
#' @param TMINN      Daily minimum air temperature vector.
#' @param RHMAXX     Daily maximum relative humidity vector.
#' @param RHMINN     Daily minimum relative humidity vector.
#' @param CCMAXX     Daily maximum cloud cover vector.
#' @param CCMINN     Daily minimum cloud cover vector.
#' @param WNMAXX     Daily maximum wind speed vector.
#' @param WNMINN     Daily minimum wind speed vector.
#' @param REFLS      Daily soil reflectance vector.
#' @param PCTWET     Daily surface wetness fraction vector.
#' @param soilinit   Initial soil temperature matrix.
#' @param hori       Horizon angles vector (24 values).
#' @param TAI        Solar attenuation vector.
#' @param soilprops  Soil thermal properties matrix.
#' @param moists     Initial soil moisture matrix.
#' @param RAINFALL   Daily rainfall vector (mm).
#' @param deepsoil   Deep soil boundary temperature vector/scalar.
#' @param PE         Air entry potential vector (19 levels).
#' @param BD         Bulk density vector (19 levels).
#' @param DD         Mineral density vector (19 levels).
#' @param BB         Soil b parameter vector (19 levels).
#' @param KS         Saturated conductivity vector (19 levels).
#' @param L          Root density vector (19 levels).
#' @param LAI        Leaf area index vector.
#' @param tides      Tidal state matrix.
#' @param TAIRhr     Hourly air temperature vector.
#' @param RHhr       Hourly relative humidity vector.
#' @param WNhr       Hourly wind speed vector.
#' @param CLDhr      Hourly cloud cover vector.
#' @param SOLRhr     Hourly solar radiation vector.
#' @param RAINhr     Hourly rainfall vector.
#' @param ZENhr      Hourly zenith angle vector.
#' @param IRDhr      Hourly longwave radiation vector.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
write_micro_csv <- function(microinput, doy, SLES, DEP, Nodes,
                            MAXSHADES, MINSHADES, TIMAXS, TIMINS,
                            TMAXX, TMINN, RHMAXX, RHMINN,
                            CCMAXX, CCMINN, WNMAXX, WNMINN,
                            REFLS, PCTWET, soilinit, hori, TAI,
                            soilprops, moists, RAINFALL, deepsoil,
                            PE, BD, DD, BB, KS, L, LAI, tides,
                            TAIRhr, RHhr, WNhr, CLDhr, SOLRhr,
                            RAINhr, ZENhr, IRDhr) {
  if (!dir.exists("micro csv input")) {
    dir.create("micro csv input")
  }
  wt <- function(x, fname) {
    write.table(x, file = file.path("micro csv input", fname),
                sep = ",", col.names = NA, qmethod = "double")
  }
  wt(as.matrix(microinput), "microinput.csv")
  wt(doy,        "doy.csv")
  wt(SLES,       "SLES.csv")
  wt(DEP,        "DEP.csv")
  wt(Nodes,      "Nodes.csv")
  wt(MAXSHADES,  "Maxshades.csv")
  wt(MINSHADES,  "Minshades.csv")
  wt(TIMAXS,     "TIMAXS.csv")
  wt(TIMINS,     "TIMINS.csv")
  wt(TMAXX,      "TMAXX.csv")
  wt(TMINN,      "TMINN.csv")
  wt(RHMAXX,     "RHMAXX.csv")
  wt(RHMINN,     "RHMINN.csv")
  wt(CCMAXX,     "CCMAXX.csv")
  wt(CCMINN,     "CCMINN.csv")
  wt(WNMAXX,     "WNMAXX.csv")
  wt(WNMINN,     "WNMINN.csv")
  wt(REFLS,      "REFLS.csv")
  wt(PCTWET,     "PCTWET.csv")
  wt(soilinit,   "soilinit.csv")
  wt(hori,       "hori.csv")
  wt(TAI,        "TAI.csv")
  wt(soilprops,  "soilprop.csv")
  wt(moists,     "moists.csv")
  wt(RAINFALL,   "rain.csv")
  wt(deepsoil,   "tannulrun.csv")
  wt(PE,         "PE.csv")
  wt(BD,         "BD.csv")
  wt(DD,         "DD.csv")
  wt(BB,         "BB.csv")
  wt(KS,         "KS.csv")
  wt(L,          "L.csv")
  wt(LAI,        "LAI.csv")
  wt(tides,      "tides.csv")
  wt(TAIRhr,     "TAIRhr.csv")
  wt(RHhr,       "RHhr.csv")
  wt(WNhr,       "WNhr.csv")
  wt(CLDhr,      "CLDhr.csv")
  wt(SOLRhr,     "SOLRhr.csv")
  wt(RAINhr,     "RAINhr.csv")
  wt(ZENhr,      "ZENhr.csv")
  wt(IRDhr,      "IRDhr.csv")
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# 7.  Process Fortran output
# ---------------------------------------------------------------------------

#' Unpack and post-process microclimate model output
#'
#' Extracts all output arrays from the list returned by \code{microclimate()},
#' fills in placeholder values when the soil moisture model is not run
#' (\code{runmoist = 0}), and optionally extracts snow and spectral outputs.
#'
#' @param microut   List returned by \code{microclimate()}.
#' @param runmoist  Flag: 1 = soil moisture model was run, 0 = it was not.
#' @param snowmodel Flag: 1 = snow model was run.
#' @param lamb      Flag: 1 = wavelength-specific solar output was requested.
#'
#' @return Named list containing all post-processed output arrays:
#'   \code{metout}, \code{shadmet}, \code{soil}, \code{shadsoil},
#'   \code{soilmoist}, \code{shadmoist}, \code{humid}, \code{shadhumid},
#'   \code{soilpot}, \code{shadpot}, \code{plant}, \code{shadplant},
#'   \code{tcond}, \code{shadtcond}, \code{specheat}, \code{shadspecheat},
#'   \code{densit}, \code{shaddensit}, and (conditionally) \code{sunsnow},
#'   \code{shdsnow}, \code{drlam}, \code{drrlam}, \code{srlam}.
#'
#' @keywords internal
process_micro_output <- function(microut, runmoist, snowmodel, lamb) {
  # --- Above-ground micromet conditions -------------------------------------
  metout    <- microut$metout    # min-shade above-ground conditions
  shadmet   <- microut$shadmet   # max-shade above-ground conditions

  # --- Soil temperatures ----------------------------------------------------
  soil      <- microut$soil      # min-shade soil temperatures
  shadsoil  <- microut$shadsoil  # max-shade soil temperatures

  # --- Soil thermal properties (always available) ---------------------------
  tcond        <- microut$tcond
  shadtcond    <- microut$shadtcond
  specheat     <- microut$specheat
  shadspecheat <- microut$shadspecheat
  densit       <- microut$densit
  shaddensit   <- microut$shaddensit

  # --- Soil moisture, humidity and plant outputs ----------------------------
  if (runmoist == 1) {
    soilmoist  <- microut$soilmoist
    shadmoist  <- microut$shadmoist
    humid      <- microut$humid
    shadhumid  <- microut$shadhumid
    soilpot    <- microut$soilpot
    shadpot    <- microut$shadpot
    plant      <- microut$plant
    shadplant  <- microut$shadplant
  } else {
    # Fill placeholders so downstream code always has these arrays available
    soilpot   <- soil;    soilpot[, 3:12]   <- 0     # water potential = 0
    shadpot   <- soil;    shadpot[, 3:12]   <- 0
    soilmoist <- soil;    soilmoist[, 3:12] <- 0.5   # arbitrary mid-range moisture
    shadmoist <- soil;    shadmoist[, 3:12] <- 0.5
    humid     <- soil;    humid[, 3:12]     <- 0.99  # near-saturated humidity
    shadhumid <- soil;    shadhumid[, 3:12] <- 0.99
    plant     <- cbind(soil, soil[, 3:4]);  plant[, 3:14]     <- 0
    shadplant <- cbind(soil, soil[, 3:4]);  shadplant[, 3:14] <- 0
  }

  # --- Snow outputs (only when snowmodel = 1) --------------------------------
  sunsnow <- NULL
  shdsnow <- NULL
  if (snowmodel == 1) {
    sunsnow <- microut$sunsnow
    shdsnow <- microut$shdsnow
  }

  # --- Wavelength-specific solar outputs (only when lamb = 1) ---------------
  drlam  <- NULL
  drrlam <- NULL
  srlam  <- NULL
  if (lamb == 1) {
    drlam  <- as.data.frame(microut$drlam)   # direct solar irradiance
    drrlam <- as.data.frame(microut$drrlam)  # direct Rayleigh component
    srlam  <- as.data.frame(microut$srlam)   # scattered solar irradiance
  }

  list(
    metout       = metout,       shadmet       = shadmet,
    soil         = soil,         shadsoil      = shadsoil,
    soilmoist    = soilmoist,    shadmoist     = shadmoist,
    humid        = humid,        shadhumid     = shadhumid,
    soilpot      = soilpot,      shadpot       = shadpot,
    plant        = plant,        shadplant     = shadplant,
    tcond        = tcond,        shadtcond     = shadtcond,
    specheat     = specheat,     shadspecheat  = shadspecheat,
    densit       = densit,       shaddensit    = shaddensit,
    sunsnow      = sunsnow,      shdsnow       = shdsnow,
    drlam        = drlam,        drrlam        = drrlam,
    srlam        = srlam
  )
}


# ---------------------------------------------------------------------------
# 8.  Build the final return list
# ---------------------------------------------------------------------------

#' Assemble the final return list for a micro_*() function
#'
#' Constructs the named list returned to the user when the full microclimate
#' model has been run (\code{runmicro = TRUE}).  Handles the four combinations
#' of \code{snowmodel} and \code{lamb} flags.  Dataset-specific extra outputs
#' (e.g., \code{SILO.data}, \code{SLOPE}, \code{microclima.out}) are appended
#' via the \code{extra} argument.
#'
#' @param out         Named list from \code{\link{process_micro_output}}.
#' @param RAINFALL    Daily rainfall vector (mm).
#' @param ndays       Total simulation days.
#' @param ALTT        Site elevation (m).
#' @param REFL        Soil reflectance used (scalar).
#' @param longlat     Numeric vector \code{c(longitude, latitude)}.
#' @param nyears      Number of simulated years.
#' @param timeinterval Number of time steps per year.
#' @param MINSHADES   Daily minimum shade vector (\%).
#' @param MAXSHADES   Daily maximum shade vector (\%).
#' @param DEP         Soil depth nodes vector (cm).
#' @param dates       POSIXct or numeric vector of hourly date/time values.
#' @param dates2      POSIXct or numeric vector of daily date values.
#' @param PE          Air entry potential vector (19 levels).
#' @param BD          Bulk density vector (19 levels).
#' @param DD          Mineral density vector (19 levels).
#' @param BB          Soil b parameter vector (19 levels).
#' @param KS          Saturated conductivity vector (19 levels).
#' @param dem         Digital elevation model raster (or \code{NA}).
#' @param diffuse_frac Hourly diffuse fraction of solar radiation (or \code{NA}).
#' @param snowmodel   Flag: 1 = snow model was run.
#' @param lamb        Flag: 1 = wavelength-specific solar output available.
#' @param extra       Optional named list of additional dataset-specific outputs
#'   to append to the return list (e.g., \code{list(SILO.data = ...)}).
#'
#' @return Named list of microclimate model outputs.
#'
#' @keywords internal
build_micro_return <- function(out, RAINFALL, ndays, ALTT, REFL, longlat,
                               nyears, timeinterval, MINSHADES, MAXSHADES,
                               DEP, dates, dates2, PE, BD, DD, BB, KS,
                               dem, diffuse_frac, snowmodel, lamb,
                               extra = list()) {
  # Core outputs common to all variants
  base <- list(
    soil         = out$soil,
    shadsoil     = out$shadsoil,
    metout       = out$metout,
    shadmet      = out$shadmet,
    soilmoist    = out$soilmoist,
    shadmoist    = out$shadmoist,
    humid        = out$humid,
    shadhumid    = out$shadhumid,
    soilpot      = out$soilpot,
    shadpot      = out$shadpot,
    plant        = out$plant,
    shadplant    = out$shadplant,
    tcond        = out$tcond,
    shadtcond    = out$shadtcond,
    specheat     = out$specheat,
    shadspecheat = out$shadspecheat,
    densit       = out$densit,
    shaddensit   = out$shaddensit,
    RAINFALL     = RAINFALL,
    ndays        = ndays,
    elev         = ALTT,
    REFL         = REFL[1],
    longlat      = c(longlat[1], longlat[2]),
    nyears       = nyears,
    timeinterval = timeinterval,
    minshade     = MINSHADES,
    maxshade     = MAXSHADES,
    DEP          = DEP,
    dates        = dates,
    dates2       = dates2,
    PE           = PE,
    BD           = BD,
    DD           = DD,
    BB           = BB,
    KS           = KS,
    dem          = dem,
    diffuse_frac = diffuse_frac
  )

  # Add snow outputs if the snow model was run
  if (snowmodel == 1) {
    base$sunsnow <- out$sunsnow
    base$shdsnow <- out$shdsnow
  }

  # Add wavelength-specific solar outputs if requested
  if (lamb == 1) {
    base$drlam  <- out$drlam
    base$drrlam <- out$drrlam
    base$srlam  <- out$srlam
  }

  # Append any dataset-specific extras (e.g. SILO.data, SLOPE, microclima.out)
  if (length(extra) > 0) {
    base <- c(base, extra)
  }

  return(base)
}
