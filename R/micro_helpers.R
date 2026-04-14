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
#' @param run.gads    Integer: aerosol model selector. Must be 0 (use default
#'   profile), 1 (run Fortran GADS), or 2 (run R GADS).
#'
#' @return Integer: \code{0} if all checks pass, \code{1} if any check fails.
#'   Error messages are emitted via \code{message()}.
#'
#' @keywords internal
validate_micro_inputs <- function(DEP, REFL, slope, aspect, hori, SLE, ERR,
                                  RUF, D0, Usrhyt, Refhyt, EC, CMH2O,
                                  TIMAXS, TIMINS, minshade, maxshade,
                                  Thcond, Density, SpecHeat, BulkDensity,
                                  run.gads = NULL) {
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

  # --- GADS aerosol model selector ------------------------------------------
  if (!is.null(run.gads) && !(run.gads %in% c(0, 1, 2))) {
    message("ERROR: run.gads must be 0 (default profile), 1 (Fortran GADS), or 2 (R GADS).")
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


# ---------------------------------------------------------------------------
# 9.  GADS aerosol optical depth / TAI computation
# ---------------------------------------------------------------------------

#' Default TAI profile: Elterman (1970)
#'
#' Aerosol transmittance profile from Elterman, L. 1970. Vertical-attenuation
#' model with eight surface meteorological ranges 2 to 13 kilometers. U.S.
#' Airforce Cambridge Research Laboratory, Bedford, Mass.  Used as the fallback
#' in \code{micro_global} and \code{micro_terra} when \code{run.gads = 0}.
#'
#' @keywords internal
TAI_ELTERMAN <- c(0.42, 0.415, 0.412, 0.408, 0.404, 0.4, 0.395, 0.388, 0.379, 0.379, 0.379, 0.375, 0.365, 0.345, 0.314, 0.3, 0.288, 0.28, 0.273, 0.264, 0.258, 0.253, 0.248, 0.243, 0.236, 0.232, 0.227, 0.223, 0.217, 0.213, 0.21, 0.208, 0.205, 0.202, 0.201, 0.198, 0.195, 0.193, 0.191, 0.19, 0.188, 0.186, 0.184, 0.183, 0.182, 0.181, 0.178, 0.177, 0.176, 0.175, 0.175, 0.174, 0.173, 0.172, 0.171, 0.17, 0.169, 0.168, 0.167, 0.164, 0.163, 0.163, 0.162, 0.161, 0.161, 0.16, 0.159, 0.157, 0.156, 0.156, 0.155, 0.154, 0.153, 0.152, 0.15, 0.149, 0.146, 0.145, 0.142, 0.14, 0.139, 0.137, 0.135, 0.135, 0.133, 0.132, 0.131, 0.13, 0.13, 0.129, 0.129, 0.128, 0.128, 0.128, 0.127, 0.127, 0.126, 0.125, 0.124, 0.123, 0.121, 0.118, 0.117, 0.115, 0.113, 0.11, 0.108, 0.107, 0.105, 0.103, 0.1)

#' Default TAI profile: Australia (Adelaide/Melbourne representative)
#'
#' Aerosol transmittance profile representative of southern Australia,
#' generated by GADS at Adelaide/Melbourne conditions.  Used as the fallback
#' in all dataset-specific \code{micro_*()} functions (other than
#' \code{micro_global} and \code{micro_terra}) when \code{run.gads = 0}.
#'
#' @keywords internal
TAI_AUSTRALIA <- c(0.0670358341290886, 0.0662612704779235, 0.065497075238002, 0.0647431301168489, 0.0639993178022531, 0.0632655219571553, 0.0625416272145492, 0.0611230843885423, 0.0597427855962549, 0.0583998423063099, 0.0570933810229656, 0.0558225431259535, 0.0545864847111214, 0.0533843764318805, 0.0522154033414562, 0.0499736739981675, 0.047855059159556, 0.0458535417401334, 0.0439633201842001, 0.0421788036108921, 0.0404946070106968, 0.0389055464934382, 0.0374066345877315, 0.0359930755919066, 0.0346602609764008, 0.0334037648376212, 0.0322193394032758, 0.0311029105891739, 0.0300505736074963, 0.0290585886265337, 0.0281233764818952, 0.0272415144391857, 0.0264097320081524, 0.0256249068083005, 0.0248840604859789, 0.0241843546829336, 0.0235230870563317, 0.0228976873502544, 0.0223057135186581, 0.0217448478998064, 0.0212128934421699, 0.0207077699817964, 0.0202275105711489, 0.0197702578594144, 0.0193342605242809, 0.0189178697551836, 0.0177713140039894, 0.0174187914242432, 0.0170790495503944, 0.0167509836728154, 0.0164335684174899, 0.0161258546410128, 0.0158269663770596, 0.0155360978343254, 0.0152525104459325, 0.0149755299703076, 0.0147045436435285, 0.0144389973831391, 0.0141783930434343, 0.0134220329447663, 0.0131772403830191, 0.0129356456025128, 0.0126970313213065, 0.0124612184223418, 0.0122280636204822, 0.01199745718102, 0.0115436048739351, 0.0110993711778668, 0.0108808815754663, 0.0106648652077878, 0.0104513876347606, 0.0102405315676965, 0.00982708969547694, 0.00962473896278535, 0.00903679230300494, 0.00884767454432418, 0.0083031278398166, 0.00796072474935954, 0.00755817587626185, 0.00718610751850881, 0.00704629977586921, 0.00684663903049612, 0.00654155580333479, 0.00642947339729728, 0.00627223096874308, 0.00603955966866779, 0.00580920937536261, 0.00568506186880564, 0.00563167068287251, 0.00556222005081865, 0.00550522989971023, 0.00547395763028062, 0.0054478983436216, 0.00541823364504573, 0.00539532163908382, 0.00539239864119488, 0.00541690124712384, 0.00551525885358836, 0.00564825853509463, 0.00577220185074264, 0.00584222986640171, 0.00581645238345584, 0.00566088137411449, 0.00535516862329704, 0.00489914757707667, 0.00432017939770409, 0.0036813032251836, 0.00309019064543606, 0.00270890436501562, 0.00276446109239711, 0.00356019862584603)

#' Compute aerosol transmittance (TAI) via GADS or a default profile
#'
#' Runs the GADS aerosol model (Fortran or R version) to compute
#' wavelength-specific aerosol optical depth, fits a 6th-order polynomial, and
#' predicts TAI at the 111 standard wavelength bands used by the microclimate
#' model.  When \code{run.gads = 0} the supplied \code{default_tai} vector is
#' returned unchanged.
#'
#' @param longlat     Numeric vector \code{c(longitude, latitude)}.
#' @param run.gads    Integer: \code{0} = use default profile, \code{1} = run
#'   Fortran GADS, \code{2} = run R GADS.
#' @param default_tai Numeric vector of length 111 returned when
#'   \code{run.gads = 0}.  Use \code{TAI_ELTERMAN} or \code{TAI_AUSTRALIA}.
#'
#' @return Numeric vector of length 111: aerosol transmittance at each
#'   standard wavelength band.
#'
#' @keywords internal
# ---------------------------------------------------------------------------
# 10.  SoilGrids fetch
# ---------------------------------------------------------------------------

#' Fetch soil hydraulic properties from SoilGrids v2
#'
#' Downloads bulk density, silt, clay, and sand fractions from the ISRIC
#' SoilGrids v2 REST API for a given location, then applies Cosby et al. 1984
#' pedotransfer functions to derive hydraulic parameters for the microclimate
#' model.
#'
#' @param x   Numeric vector \code{c(longitude, latitude)}.
#' @param DEP Soil depth node vector (cm), length 10.
#'
#' @return Named list with elements \code{PE}, \code{KS}, \code{BB},
#'   \code{BD}, \code{BulkDensity} when data are available, or \code{NULL}
#'   when SoilGrids returns no data for the location.
#'
#' @keywords internal
fetch_soilgrids <- function(x, DEP) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("package 'jsonlite' is needed to extract data from SoilGrids, please install it.",
         call. = FALSE)
  }
  message("extracting soil texture data from SoilGrids")
  ov <- jsonlite::fromJSON(
    paste0('https://rest.isric.org/soilgrids/v2.0/properties/query?lon=', x[1],
           '&lat=', x[2],
           '&property=bdod&property=silt&property=clay&property=sand'),
    flatten = TRUE)
  if (length(ov) > 3) {
    soilpro <- cbind(
      c(0, 5, 15, 30, 60, 100),
      unlist(ov$properties$layers$depths[[1]]$values.mean) / 100,
      unlist(ov$properties$layers$depths[[2]]$values.mean) / 10,
      unlist(ov$properties$layers$depths[[4]]$values.mean) / 10,
      unlist(ov$properties$layers$depths[[3]]$values.mean) / 10
    )
    soilpro <- rbind(soilpro, soilpro[6, ])
    soilpro[7, 1] <- 200
    colnames(soilpro) <- c('depth', 'blkdens', 'clay', 'silt', 'sand')
    soil.hydro <- pedotransfer(soilpro = as.data.frame(soilpro), DEP = DEP)
    return(list(
      PE          = soil.hydro$PE,
      KS          = soil.hydro$KS,
      BB          = soil.hydro$BB,
      BD          = soil.hydro$BD,
      BulkDensity = soil.hydro$BD[seq(1, 19, 2)]
    ))
  } else {
    message("no SoilGrids data for this site, using user-input soil properties")
    return(NULL)
  }
}


compute_tai <- function(longlat, run.gads, default_tai) {
  if (run.gads > 0) {
    relhum <- 1
    if (run.gads == 1) {
      optdep.summer <- as.data.frame(rungads(longlat[2], longlat[1], relhum, 0))
      optdep.winter <- as.data.frame(rungads(longlat[2], longlat[1], relhum, 1))
    } else {
      optdep.summer <- as.data.frame(gads.r(longlat[2], longlat[1], relhum, 0))
      optdep.winter <- as.data.frame(gads.r(longlat[2], longlat[1], relhum, 1))
    }
    optdep <- cbind(optdep.winter[, 1],
                    rowMeans(cbind(optdep.summer[, 2], optdep.winter[, 2])))
    optdep <- as.data.frame(optdep)
    colnames(optdep) <- c("LAMBDA", "OPTDEPTH")
    a <- lm(OPTDEPTH ~ poly(LAMBDA, 6, raw = TRUE), data = optdep)
    LAMBDA <- c(290, 295, 300, 305, 310, 315, 320, 330, 340, 350, 360, 370,
                380, 390, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580,
                600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820,
                840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1080,
                1100, 1120, 1140, 1160, 1180, 1200, 1220, 1240, 1260, 1280,
                1300, 1320, 1380, 1400, 1420, 1440, 1460, 1480, 1500, 1540,
                1580, 1600, 1620, 1640, 1660, 1700, 1720, 1780, 1800, 1860,
                1900, 1950, 2000, 2020, 2050, 2100, 2120, 2150, 2200, 2260,
                2300, 2320, 2350, 2380, 2400, 2420, 2450, 2490, 2500, 2600,
                2700, 2800, 2900, 3000, 3100, 3200, 3300, 3400, 3500, 3600,
                3700, 3800, 3900, 4000)
    return(predict(a, data.frame(LAMBDA)))
  } else {
    return(default_tai)
  }
}


#' Download DEM and optionally compute terrain attributes
#'
#' Downloads a digital elevation model via \code{microclima::get_dem()},
#' extracts site elevation, and optionally computes slope, aspect, and horizon
#' angles. Accepts an existing raster to skip the download.
#'
#' @param lon Longitude (decimal degrees).
#' @param lat Latitude (decimal degrees).
#' @param dem.res DEM resolution (m) passed to \code{microclima::get_dem()}.
#'   \code{NULL} to omit (uses the function's own default).
#' @param zmin Minimum elevation filter (m) passed to
#'   \code{microclima::get_dem()}. \code{NULL} to omit.
#' @param pixels Number of pixels in x and y dimensions passed as
#'   \code{xdims}/\code{ydims}. \code{NULL} to omit.
#' @param existing_dem Optional existing \code{RasterLayer} or
#'   \code{SpatRaster}. When non-\code{NULL}, the download is skipped and this
#'   raster is used directly (converted to SpatRaster if needed).
#' @param terrain Integer; \code{1} to compute slope, aspect, and horizon
#'   angles.
#' @param horizon.step Angle step (degrees) between successive horizon angle
#'   readings. \code{nangles} readings are returned at
#'   \code{0, horizon.step, 2*horizon.step, ...} degrees.
#' @param nangles Number of horizon angle readings to return. Use \code{24}
#'   (default) for the elevatr/Fortran-model path and \code{36} for the
#'   microclima solar-adjustment path.
#' @param user_hori Optional user-supplied horizon angles. When non-\code{NULL}
#'   and non-\code{NA}, they are spline-interpolated to \code{nangles} rather
#'   than computed from the DEM.
#' @param slope Current slope value. A \code{logical} (default \code{FALSE})
#'   triggers computation from the DEM when \code{terrain == 1}; a numeric
#'   value is passed through unchanged.
#' @param aspect Current aspect value; same logic as \code{slope}.
#' @return Named list: \code{dem} (SpatRaster), \code{elev} (numeric),
#'   \code{slope}, \code{aspect}, \code{hori} (\code{nangles}-element vector or
#'   \code{NA}).
#' @keywords internal
fetch_dem <- function(lon, lat, dem.res = 10, zmin = NULL, pixels = 100,
                      existing_dem = NULL, terrain = 0, horizon.step = 15,
                      nangles = 24, user_hori = NULL, slope = FALSE,
                      aspect = FALSE) {
  if (is.null(existing_dem)) {
    if (!requireNamespace("microclima", quietly = TRUE)) {
      stop("package 'microclima' is needed to download a DEM via elevatr, please install it.",
           call. = FALSE)
    }
    message("downloading DEM via package elevatr")
    dem_args <- list(lat = lat, long = lon)
    if (!is.null(dem.res)) dem_args$resolution <- dem.res
    if (!is.null(zmin)) dem_args$zmin <- zmin
    if (!is.null(pixels)) { dem_args$xdims <- pixels; dem_args$ydims <- pixels }
    existing_dem <- do.call(microclima::get_dem, dem_args)
  }
  dem <- terra::rast(existing_dem)  # normalise to SpatRaster (no-op if already SpatRaster)
  xy <- data.frame(lon = lon, lat = lat) |>
    sf::st_as_sf(coords = c("lon", "lat"))
  xy <- sf::st_set_crs(xy, "EPSG:4326")
  xy <- sf::st_transform(xy, sf::st_crs(dem))
  elev <- as.numeric(terra::extract(dem, xy)[, 2])
  hori_out <- NA_real_
  if (terrain == 1) {
    message("computing slope, aspect and horizon angles")
    if (is.logical(slope)) {
      slope_r <- terra::terrain(dem, v = "slope", unit = "degrees")
      slope <- as.numeric(terra::extract(slope_r, xy)[, 2])
    }
    if (is.logical(aspect)) {
      aspect_r <- terra::terrain(dem, v = "aspect", unit = "degrees")
      aspect <- as.numeric(terra::extract(aspect_r, xy)[, 2])
    }
    if (!is.null(user_hori) && !is.na(user_hori[1])) {
      hori_out <- spline(x = user_hori, n = nangles, method = 'periodic')$y
      hori_out[hori_out < 0] <- 0
      hori_out[hori_out > 90] <- 90
    } else {
      ha <- numeric(nangles)
      for (i in seq_len(nangles)) {
        har <- microclima::horizonangle(dem, (i - 1L) * horizon.step, terra::res(dem)[1])
        ha[i] <- atan(as.numeric(terra::extract(har, xy)[, 2])) * (180 / pi)
      }
      hori_out <- ha
    }
  }
  list(dem = dem, elev = elev, slope = slope, aspect = aspect, hori = hori_out)
}


#' Compute terrain-adjusted solar radiation and diffuse fraction
#'
#' Partitions hourly total solar radiation into direct and diffuse components
#' and applies topographic correction via \code{microclima::.shortwave.ts()}.
#' This block is identical across all \code{micro_*()} functions that use the
#' microclima package for solar adjustment.
#'
#' @param dsw2 Numeric vector; hourly total downwelling shortwave radiation
#'   (W m\eqn{^{-2}}).
#' @param jd Integer vector; Julian day for each hourly time step.
#' @param hour.microclima Numeric vector; local decimal hour for each time step.
#' @param lat Latitude (decimal degrees).
#' @param long Longitude (decimal degrees).
#' @param slope Site slope (degrees).
#' @param aspect Site aspect (degrees, clockwise from north).
#' @param ha Numeric vector; horizon angle for each hourly time step (degrees),
#'   derived from the 36-angle lookup table via solar azimuth.
#' @param LOR Leaf optical reflectance passed to \code{.shortwave.ts()} as
#'   the \code{x} argument (\code{microclima.LOR} user parameter).
#' @param LAI Numeric vector; leaf area index. \code{mean(LAI)} is passed to
#'   \code{.shortwave.ts()} as \code{l}.
#' @return Named list: \code{SOLRhr_all} (terrain-adjusted hourly solar
#'   radiation, W m\eqn{^{-2}}) and \code{diffuse_frac_all} (hourly diffuse
#'   fraction, 0–1).
#' @keywords internal
compute_solar_partition <- function(dsw2, jd, hour.microclima, lat, long,
                                    slope, aspect, ha, LOR, LAI) {
  bound <- function(x, mn = 0, mx = 1) {
    x[x > mx] <- mx
    x[x < mn] <- mn
    x
  }
  si  <- microclima::siflat(hour.microclima, lat, long, jd, merid = long)
  am  <- microclima::airmasscoef(hour.microclima, lat, long, jd, merid = long)
  dp  <- vector(length = length(jd))
  for (i in seq_along(jd)) {
    dp[i] <- microclima:::difprop(dsw2[i], jd[i], hour.microclima[i], lat, long,
                                  watts = TRUE, hourly = TRUE, merid = long)
  }
  dp[dsw2 == 0] <- NA
  dnir <- (dsw2 * (1 - dp)) / si;  dnir[si == 0] <- NA
  difr  <- dsw2 * dp
  edni  <- dnir / ((4.87 / 0.0036) * (1 - dp))
  edif  <- difr / ((4.87 / 0.0036) * dp)
  odni  <- bound((log(edni) / -am), mn = 0.001, mx = 1.7)
  odif  <- bound((log(edif) / -am), mn = 0.001, mx = 1.7)
  nd    <- length(odni)
  sel   <- which(!is.na(am * dp * odni * odif))
  dp[1]   <- dp[min(sel)];   odni[1]   <- odni[min(sel)];   odif[1]   <- odif[min(sel)]
  dp[nd]  <- dp[max(sel)];   odni[nd]  <- odni[max(sel)];   odif[nd]  <- odif[max(sel)]
  dp    <- zoo::na.approx(dp,   na.rm = FALSE)
  odni  <- zoo::na.approx(odni, na.rm = FALSE)
  odif  <- zoo::na.approx(odif, na.rm = FALSE)
  h_dp  <- bound(dp)
  h_oi  <- bound(odni, mn = 0.24, mx = 1.7)
  h_od  <- bound(odif, mn = 0.24, mx = 1.7)
  afi   <- exp(-am * h_oi);  afd <- exp(-am * h_od)
  h_dni <- (1 - h_dp) * afi * 4.87 / 0.0036
  h_dif <- h_dp * afd * 4.87 / 0.0036
  h_dni[si == 0]     <- 0
  h_dif[is.na(h_dif)] <- 0
  diffuse_frac_all <- h_dif / (h_dni + h_dif)
  diffuse_frac_all[is.na(diffuse_frac_all)] <- 1
  radwind2 <- .shortwave.ts(h_dni * 0.0036, h_dif * 0.0036, jd,
                            hour.microclima, lat, long, slope, aspect,
                            ha = ha, svv = 1, x = LOR, l = mean(LAI),
                            albr = 0, merid = long, dst = 0, difani = FALSE)
  SOLRhr_all <- radwind2$swrad / 0.0036
  list(SOLRhr_all = SOLRhr_all, diffuse_frac_all = diffuse_frac_all)
}
