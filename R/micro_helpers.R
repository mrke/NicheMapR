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
#' @param depths         Numeric vector (length 10) of soil node depths (cm).
#' @param surface_reflectivity        Soil solar reflectance (0–1).
#' @param slope       Slope in degrees (0–90).
#' @param aspect      Aspect in degrees (0–365).
#' @param horizon_angles        Numeric vector (length 24) of horizon angles (degrees).
#' @param surface_emissivity         Substrate longwave emissivity (0.05–1).
#' @param tolerance         Integrator error tolerance (must be > 0).
#' @param roughness_height         Roughness height, m (0.0001–2).
#' @param zero_plane_displacement          Zero-plane displacement height, m.
#' @param local_height      Local height for met output, m (must exceed roughness_height and be
#'                    no greater than reference_height).
#' @param reference_height      Reference height of input met data, m.
#' @param orbital_eccentricity          Orbital eccentricity (0.0034–0.058).
#' @param precipitable_water       Precipitable water in air column, cm (0.1–2).
#' @param maxima_times      Numeric vector (length 4): time offsets for weather maxima.
#' @param minima_times      Numeric vector (length 4): time offsets for weather minima.
#' @param minimum_shade    Minimum shade level(s), \% (0–100).
#' @param maximum_shade    Maximum shade level(s), \% (0–100, must exceed minimum_shade).
#' @param mineral_conductivity      Soil mineral thermal conductivity, W/mK (must be >= 0).
#' @param mineral_density     Soil mineral density, Mg/m3 (must be >= 0).
#' @param mineral_heat_capacity    Soil mineral specific heat, J/kg-K (must be >= 0).
#' @param bulk_density Soil bulk density, Mg/m3 (must be >= 0).
#' @param global_aerosol_database    Integer: aerosol model selector. Must be 0 (use default
#'   profile), 1 (run Fortran GADS), or 2 (run R GADS).
#'
#' @return Integer: \code{0} if all checks pass, \code{1} if any check fails.
#'   Error messages are emitted via \code{message()}.
#'
#' @keywords internal
validate_micro_inputs <- function(depths, surface_reflectivity, slope, aspect, horizon_angles, surface_emissivity, tolerance,
                                  roughness_height, zero_plane_displacement, local_height, reference_height, orbital_eccentricity, precipitable_water,
                                  maxima_times, minima_times, minimum_shade, maximum_shade,
                                  mineral_conductivity, mineral_density, mineral_heat_capacity, bulk_density,
                                  global_aerosol_database = NULL) {
  errors <- 0

  # --- Soil depth node spacing (warnings only, do not set errors) -----------
  if (depths[2] - depths[1] > 3 | depths[3] - depths[2] > 3) {
    message("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing")
  }
  if (depths[2] - depths[1] < 2) {
    message("warning, nodes might be too close near the surface, try a different spacing if the program is crashing")
  }
  if (depths[10] != 200) {
    message("warning, last depth in soil_temperature_lowshade should not be changed from 200 without good reason")
  }

  # --- depths vector structure -------------------------------------------------
  if (depths[1] != 0) {
    message("ERROR: First soil_temperature_lowshade node (depths[1]) must = 0 cm. Please correct.")
    errors <- 1
  }
  if (length(depths) != 10) {
    message("ERROR: You must enter 10 different soil_temperature_lowshade depths.")
    errors <- 1
  }
  for (i in 1:9) {
    if (depths[i + 1] <= depths[i]) {
      message("ERROR: Soil depth (depths array) is not in ascending order.")
      errors <- 1
    }
  }
  if (depths[10] > 500) {
    message("ERROR: Deepest soil_temperature_lowshade depth (depths array) is too large (must be <= 500 cm).")
    errors <- 1
  }

  # --- Roughness and displacement heights -----------------------------------
  if (roughness_height < 0.0001) {
    message("ERROR: The roughness height (roughness_height) is too small (< 0.0001). Please enter a larger value.")
    errors <- 1
  }
  if (roughness_height > 2) {
    message("ERROR: The roughness height (roughness_height) is too large (> 2). Please enter a smaller value.")
    errors <- 1
  }
  if (zero_plane_displacement > 0 & zero_plane_displacement < local_height) {
    message("ERROR: The zero plane displacement height (zero_plane_displacement) must be lower than the local height (local_height). Please enter a smaller value.")
    errors <- 1
  }
  if (local_height < roughness_height) {
    message("ERROR: Local height (local_height) smaller than roughness height (roughness_height). Please use a larger height above the surface.")
    errors <- 1
  }
  if (local_height > reference_height) {
    message("ERROR: Reference height is less than local height (local_height).")
    errors <- 1
  }

  # --- Soil physical properties ---------------------------------------------
  if (min(mineral_conductivity) < 0) {
    message("ERROR: Thermal conductivity (mineral_conductivity) is negative. Please input a positive value.")
    errors <- 1
  }
  if (min(mineral_density) < 0) {
    message("ERROR: Soil mineral density (mineral_density) is negative. Please input a positive value.")
    errors <- 1
  }
  if (min(mineral_heat_capacity) < 0) {
    message("ERROR: Specific heat (mineral_heat_capacity) is negative. Please input a positive value.")
    errors <- 1
  }
  if (min(bulk_density) < 0) {
    message("ERROR: Bulk density (bulk_density) is negative. Please input a positive value.")
    errors <- 1
  }

  # --- Radiation and site geometry ------------------------------------------
  if (surface_reflectivity < 0 | surface_reflectivity > 1) {
    message("ERROR: Soil reflectivity (surface_reflectivity) is out of bounds. Please input a value between 0 and 1.")
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
  if (max(horizon_angles) > 90 | min(horizon_angles) < 0) {
    message("ERROR: At least one horizon angle (horizon_angles) is out of bounds. Please input values between 0 and 90.")
    errors <- 1
  }
  if (length(horizon_angles) != 24) {
    message("ERROR: You must enter 24 horizon angle values.")
    errors <- 1
  }
  if (surface_emissivity < 0.05 | surface_emissivity > 1) {
    message("ERROR: Emissivity (surface_emissivity) is out of bounds. Please enter a value between 0.05 and 1.00.")
    errors <- 1
  }

  # --- Atmospheric and orbital parameters -----------------------------------
  if (tolerance < 0) {
    message("ERROR: Error bound (tolerance) is too small. Please enter a value > 0.")
    errors <- 1
  }
  if (orbital_eccentricity < 0.0034 | orbital_eccentricity > 0.058) {
    message("ERROR: Eccentricity (orbital_eccentricity) is out of bounds. Please enter a value between 0.0034 and 0.058.")
    errors <- 1
  }
  if (precipitable_water < 0.1 | precipitable_water > 2) {
    message("ERROR: Precipitable water (precipitable_water) is out of bounds. Please enter a value between 0.1 and 2.")
    errors <- 1
  }

  # --- Weather timing offsets -----------------------------------------------
  if (max(maxima_times) > 24 | min(maxima_times) < 0) {
    message("ERROR: At least one time of weather maximum (maxima_times) is out of bounds. Please input values between 0 and 24.")
    errors <- 1
  }
  if (max(minima_times) > 24 | min(minima_times) < 0) {
    message("ERROR: At least one time of weather minimum (minima_times) is out of bounds. Please input values between 0 and 24.")
    errors <- 1
  }

  # --- GADS aerosol model selector ------------------------------------------
  if (!is.null(global_aerosol_database) && !(global_aerosol_database %in% c(0, 1, 2))) {
    message("ERROR: global_aerosol_database must be 0 (default profile), 1 (Fortran GADS), or 2 (R GADS).")
    errors <- 1
  }

  # --- Shade levels ---------------------------------------------------------
  if (max(minimum_shade - maximum_shade) >= 0) {
    message("ERROR: Minimum shade (minimum_shade) must be less than maximum shade (maximum_shade). Please correct.")
    errors <- 1
  }
  if (max(minimum_shade) > 100 | min(minimum_shade) < 0) {
    message("ERROR: Minimum shade (minimum_shade) is out of bounds. Please input a value between 0 and 100.")
    errors <- 1
  }
  if (max(maximum_shade) > 100 | min(maximum_shade) < 0) {
    message("ERROR: Maximum shade (maximum_shade) is out of bounds. Please input a value between 0 and 100.")
    errors <- 1
  }

  return(errors)
}


# ---------------------------------------------------------------------------
# 2.  Shade vector construction
# ---------------------------------------------------------------------------

#' Create daily shade fraction vectors
#'
#' Expands scalar or short \code{minimum_shade}/\code{maximum_shade} inputs to a full
#' vector of length \code{ndays}, as required by the Fortran model.
#'
#' @param minimum_shade Scalar or vector of minimum shade values (\%).
#' @param maximum_shade Scalar or vector of maximum shade values (\%).
#' @param ndays    Total number of simulation days.
#'
#' @return Named list with elements \code{minimum_shade_daily} and \code{maximum_shade_daily},
#'   each a numeric vector of length \code{ndays}.
#'
#' @keywords internal
setup_shade_vectors <- function(minimum_shade, maximum_shade, ndays) {
  if (length(minimum_shade) != ndays) {
    minimum_shade_daily <- rep(minimum_shade[1], ndays)
  } else {
    minimum_shade_daily <- minimum_shade
  }
  if (length(maximum_shade) != ndays) {
    maximum_shade_daily <- rep(maximum_shade[1], ndays)
  } else {
    maximum_shade_daily <- maximum_shade
  }
  list(minimum_shade_daily = minimum_shade_daily, maximum_shade_daily = maximum_shade_daily)
}


# ---------------------------------------------------------------------------
# 3.  Timezone reference longitude (solar_noon_longitude)
# ---------------------------------------------------------------------------

#' Compute the reference longitude for solar noon correction
#'
#' Returns \code{solar_noon_longitude}, the longitude (degrees) used by the Fortran model to
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
    solar_noon_longitude <- as.numeric(geonames::GNtimezone(lat, lon)[4]) * -15
  } else {
    # Use the nearest whole-degree longitude as the solar noon reference
    solar_noon_longitude <- abs(trunc(lon))
  }
  return(solar_noon_longitude)
}


# ---------------------------------------------------------------------------
# 4.  Microinput vector construction
# ---------------------------------------------------------------------------

#' Build the micro_input parameter vector for the Fortran microclimate model
#'
#' Assembles the numeric vector consumed by the Fortran \code{microclimate()}
#' routine.  All \code{micro_*()} drivers use the same vector layout; the only
#' difference between the baseline dataset drivers and the ERA5/NCEP/OpenMeteo
#' drivers is the \code{dewrain} and \code{moiststep} parameters (defaulting to
#' \code{0} and \code{360} respectively for the baseline group).
#'
#' @param ndays     Total simulation days.
#' @param roughness_height       Roughness height (m).
#' @param tolerance       Integrator error tolerance.
#' @param local_height    Local height for met output (m).
#' @param reference_height    Reference height of input met data (m).
#' @param Numtyps   Number of soil types.
#' @param roughness_height_1       Top segment roughness height (m).
#' @param roughness_height_2       Second segment roughness height (m).
#' @param wind_profile_height_1       Top of first segment height (m).
#' @param wind_profile_height_2       Second segment height (m).
#' @param start_day    Start day index.
#' @param end_day       End day index.
#' @param hemisphere     Hemisphere flag (1 = north, 2 = south).
#' @param latitude_degrees      Latitude degrees (absolute value, integer part).
#' @param latitude_minutes    Latitude minutes (decimal).
#' @param longitude_degrees     Longitude degrees (absolute value, integer part).
#' @param longitude_minutes    Longitude minutes (decimal).
#' @param solar_noon_longitude     Reference longitude for solar noon correction (degrees).
#' @param slope     Surface slope (degrees).
#' @param azmuth    Surface aspect/azimuth (degrees).
#' @param elevation      Site elevation (m).
#' @param precipitable_water     Precipitable water in air column (cm H2O).
#' @param microdaily  Flag: 1 = daily iteration mode (365-day runs), 0 = normal.
#' @param mean_annual_temperature    Mean annual air temperature used to initialise deep soil (°C).
#' @param orbital_eccentricity        Orbital eccentricity.
#' @param VIEWF     Sky view factor (0–1).
#' @param snowtemp  Temperature threshold for snowfall (°C).
#' @param snowdens  Snow density (Mg/m3).
#' @param snowmelt  Proportion of calculated melt that does not refreeze.
#' @param undercatch Undercatch multiplier for rain-to-snow conversion.
#' @param rain_multiplier  Rain multiplier for surface soil moisture.
#' @param runshade  Flag: 1 = run both shade levels, 0 = minimum shade only.
#' @param soil_moisture_model  Flag: 1 = run soil moisture model.
#' @param maximum_pooling_depth   Maximum water pooling depth (mm).
#' @param evenrain  Flag: 1 = spread rainfall evenly over 24 h, 0 = midnight event.
#' @param snow_model Flag: 1 = run snow model.
#' @param rainmelt  Rain melt parameter.
#' @param output_to_csv  Flag: 1 = write Fortran output as CSV files.
#' @param densfun   Snow density function coefficients (vector of 4).
#' @param hourly    Flag: 1 = hourly input data provided.
#' @param rainhourly Flag: 1 = hourly rainfall provided.
#' @param radiation_per_wavelength      Flag: 1 = return wavelength-specific solar output.
#' @param scattered_uv       Flag: 1 = use gamma function for scattered solar radiation.
#' @param root_resistance        Resistance per unit root length (m3 kg-1 s-1).
#' @param stomatal_closure_potential        Critical leaf water potential for stomatal closure (J kg-1).
#' @param leaf_resistance        Resistance per unit leaf length (m3 kg-1 s-1).
#' @param stomatal_stability_parameter        Stability parameter for stomatal closure.
#' @param root_radius        Root radius (m).
#' @param moist_error        Maximum allowable mass balance error (kg).
#' @param moist_count  Maximum iterations for mass balance.
#' @param longwave_radiation_model        Longwave radiation formula: 0 = Campbell & Norman, 1 = Swinbank.
#' @param message   Flag: 1 = allow Fortran integrator warnings.
#' @param fail      Maximum integrator restarts before quitting.
#' @param snowcond  Effective snow thermal conductivity W/mC (0 = use density function).
#' @param intercept Snow interception fraction under shade (0–1).
#' @param grasshade Flag: 1 = remove shade when snow present (shade from grass/low shrubs).
#' @param solar_model_only   Flag: 1 = only run SOLRAD solar radiation sub-model.
#' @param canopy_roughness_height        Heat transfer roughness height (m).
#' @param zero_plane_displacement        Zero-plane displacement height (m).
#' @param maxima_times    Vector of 4: time offsets for weather maxima (h).
#' @param minima_times    Vector of 4: time offsets for weather minima (h).
#' @param spinup    Number of spin-up days.
#' @param maxsurf   Maximum surface temperature (°C).
#' @param max_iterations_per_day     Number of day iterations for steady periodic solution.
#' @param dewrain   Dew/rain parameter (ERA5/NCEP/OpenMeteo only; default 0).
#' @param moiststep Moisture time step (ERA5/NCEP/OpenMeteo only; default 360).
#'
#' @return Numeric vector of microclimate model input parameters.
#'
#' @keywords internal
build_microinput <- function(ndays, roughness_height, tolerance, local_height, reference_height, Numtyps,
                             roughness_height_1, roughness_height_2, wind_profile_height_1, wind_profile_height_2, start_day, end_day,
                             hemisphere, latitude_degrees, latitude_minutes, longitude_degrees, longitude_minutes, solar_noon_longitude,
                             slope, azmuth, elevation, precipitable_water, microdaily, mean_annual_temperature,
                             orbital_eccentricity, VIEWF, snowtemp, snowdens, snowmelt,
                             undercatch, rain_multiplier, runshade, soil_moisture_model,
                             maximum_pooling_depth, evenrain, snow_model, rainmelt, output_to_csv,
                             densfun, hourly, rainhourly, radiation_per_wavelength, scattered_uv,
                             root_resistance, stomatal_closure_potential, leaf_resistance, stomatal_stability_parameter, root_radius, moist_error, moist_count, longwave_radiation_model,
                             message, fail, snowcond, intercept, grasshade,
                             solar_model_only, canopy_roughness_height, zero_plane_displacement, maxima_times, minima_times, spinup,
                             maxsurf, max_iterations_per_day,
                             dewrain = 0, moiststep = 360) {
  c(ndays, roughness_height, tolerance, local_height, reference_height, Numtyps, roughness_height_1, roughness_height_2, wind_profile_height_1, wind_profile_height_2,
    start_day, end_day, hemisphere, latitude_degrees, latitude_minutes, longitude_degrees, longitude_minutes, solar_noon_longitude,
    slope, azmuth, elevation, precipitable_water, microdaily, mean_annual_temperature,
    orbital_eccentricity, VIEWF, snowtemp, snowdens, snowmelt, undercatch, rain_multiplier,
    runshade, soil_moisture_model, maximum_pooling_depth, evenrain, snow_model, rainmelt, output_to_csv,
    densfun, hourly, rainhourly, radiation_per_wavelength, scattered_uv, root_resistance, stomatal_closure_potential, leaf_resistance, stomatal_stability_parameter, root_radius, moist_error,
    moist_count, longwave_radiation_model, message, fail, snowcond, intercept, grasshade, solar_model_only,
    canopy_roughness_height, zero_plane_displacement, maxima_times, minima_times, spinup, dewrain, moiststep, maxsurf, max_iterations_per_day)
}


# ---------------------------------------------------------------------------
# 5.  Micro input list construction
# ---------------------------------------------------------------------------

#' Assemble the micro input list for the microclimate() Fortran wrapper
#'
#' Builds the named list that is passed directly to \code{microclimate()}.
#' Every element corresponds to an array expected by the Fortran subroutine.
#'
#' @param micro_input Numeric vector from \code{\link{build_microinput}}.
#' @param day_of_year        Integer vector of day-of-year values for each simulation day.
#' @param surface_emissivity       Numeric vector of surface emissivity values (one per day).
#' @param depths        Numeric vector of 10 soil node depths (cm).
#' @param soil_nodes      Matrix of node weights/depths for the Fortran solver.
#' @param maximum_shade_daily  Daily maximum shade fractions (\%).
#' @param minimum_shade_daily  Daily minimum shade fractions (\%).
#' @param reference_temperature_max      Daily maximum air temperatures (°C).
#' @param reference_temperature_min      Daily minimum air temperatures (°C).
#' @param reference_humidity_max     Daily maximum relative humidity (\%).
#' @param reference_humidity_min     Daily minimum relative humidity (\%).
#' @param cloud_max     Daily maximum cloud cover (oktas or fraction).
#' @param cloud_min     Daily minimum cloud cover.
#' @param reference_wind_max     Daily maximum wind speed (m/s).
#' @param reference_wind_min     Daily minimum wind speed (m/s).
#' @param reference_temperature     Hourly air temperature override (°C); zeros if not used.
#' @param reference_humidity       Hourly relative humidity override (\%); zeros if not used.
#' @param reference_wind_speed       Hourly wind speed override (m/s); zeros if not used.
#' @param cloud_cover      Hourly cloud cover override; zeros if not used.
#' @param global_radiation     Hourly solar radiation override (W/m2); zeros if not used.
#' @param rainfall_hourly     Hourly rainfall override (mm); zeros if not used.
#' @param zenith_angle_hourly      Hourly zenith angle override (degrees); -1 if not used.
#' @param longwave_radiation      Hourly longwave radiation override (W/m2); -1 if not used.
#' @param albedo      Daily soil reflectance values.
#' @param soil_wetness     Daily fraction of surface acting as free water.
#' @param soilinit   Matrix of initial soil temperatures at each node.
#' @param horizon_angles       Horizon angles (24 values, degrees).
#' @param aerosol_optical_depth        Solar attenuation vector (111 wavelength bins); 0 = use GADS.
#' @param soil_properties  Matrix of soil thermal properties at each node.
#' @param soil_moisture_profile     Matrix of initial soil moisture at each node.
#' @param rainfall   Daily rainfall (mm).
#' @param deep_soil_temperature   Annual mean temperature for deep soil boundary condition (°C).
#' @param air_entry_water_potential         Air entry potential (J/kg) at each of 19 soil levels.
#' @param saturated_hydraulic_conductivity         Saturated conductivity (kg s/m3) at each of 19 soil levels.
#' @param campbell_b_parameter         Campbell's soil b parameter at each of 19 soil levels.
#' @param soil_bulk_density         Bulk density (Mg/m3) at each of 19 soil levels.
#' @param soil_mineral_density         Mineral density (Mg/m3) at each of 19 soil levels.
#' @param root_density          Root density (m/m3) at each of 19 soil levels.
#' @param leaf_area_index        Leaf area index (one value per simulation day).
#' @param tides      Matrix (ndays*24 x 3) of tidal state, temperature and splash.
#'
#' @return Named list ready to pass to \code{microclimate()}.
#'
#' @keywords internal
build_micro_list <- function(micro_input, day_of_year, surface_emissivity, depths, soil_nodes,
                             maximum_shade_daily, minimum_shade_daily,
                             reference_temperature_max, reference_temperature_min, reference_humidity_max, reference_humidity_min,
                             cloud_max, cloud_min, reference_wind_max, reference_wind_min,
                             reference_temperature, reference_humidity, reference_wind_speed, cloud_cover, global_radiation,
                             rainfall_hourly, zenith_angle_hourly, longwave_radiation,
                             albedo, soil_wetness, soilinit, horizon_angles, aerosol_optical_depth,
                             soil_properties, soil_moisture_profile, rainfall, deep_soil_temperature,
                             air_entry_water_potential, saturated_hydraulic_conductivity, campbell_b_parameter, soil_bulk_density, soil_mineral_density, root_density, leaf_area_index, tides) {
  list(
    tides      = tides,
    micro_input = micro_input,
    day_of_year        = day_of_year,
    surface_emissivity       = surface_emissivity,
    depths        = depths,
    soil_nodes      = soil_nodes,
    maximum_shade_daily  = maximum_shade_daily,
    minimum_shade_daily  = minimum_shade_daily,
    reference_temperature_max      = reference_temperature_max,
    reference_temperature_min      = reference_temperature_min,
    reference_humidity_max     = reference_humidity_max,
    reference_humidity_min     = reference_humidity_min,
    cloud_max     = cloud_max,
    cloud_min     = cloud_min,
    reference_wind_max     = reference_wind_max,
    reference_wind_min     = reference_wind_min,
    reference_temperature     = reference_temperature,
    reference_humidity       = reference_humidity,
    reference_wind_speed       = reference_wind_speed,
    cloud_cover      = cloud_cover,
    global_radiation     = global_radiation,
    rainfall_hourly     = rainfall_hourly,
    zenith_angle_hourly      = zenith_angle_hourly,
    longwave_radiation      = longwave_radiation,
    albedo      = albedo,
    soil_wetness     = soil_wetness,
    soilinit   = soilinit,
    horizon_angles       = horizon_angles,
    aerosol_optical_depth        = aerosol_optical_depth,
    soil_properties  = soil_properties,
    soil_moisture_profile     = soil_moisture_profile,
    rainfall   = rainfall,
    deep_soil_temperature  = deep_soil_temperature,
    air_entry_water_potential         = air_entry_water_potential,
    saturated_hydraulic_conductivity         = saturated_hydraulic_conductivity,
    campbell_b_parameter         = campbell_b_parameter,
    soil_bulk_density         = soil_bulk_density,
    soil_mineral_density         = soil_mineral_density,
    root_density          = root_density,
    leaf_area_index        = leaf_area_index
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
#' @param micro_input Numeric micro_input vector (from \code{\link{build_microinput}}).
#' @param day_of_year        Day-of-year vector.
#' @param surface_emissivity       Surface emissivity vector.
#' @param depths        Soil depths vector (cm).
#' @param soil_nodes      Soil node weight matrix.
#' @param maximum_shade_daily  Daily maximum shade vector.
#' @param minimum_shade_daily  Daily minimum shade vector.
#' @param maxima_times     Time-of-maximum offsets vector.
#' @param minima_times     Time-of-minimum offsets vector.
#' @param reference_temperature_max      Daily maximum air temperature vector.
#' @param reference_temperature_min      Daily minimum air temperature vector.
#' @param reference_humidity_max     Daily maximum relative humidity vector.
#' @param reference_humidity_min     Daily minimum relative humidity vector.
#' @param cloud_max     Daily maximum cloud cover vector.
#' @param cloud_min     Daily minimum cloud cover vector.
#' @param reference_wind_max     Daily maximum wind speed vector.
#' @param reference_wind_min     Daily minimum wind speed vector.
#' @param albedo      Daily soil reflectance vector.
#' @param soil_wetness     Daily surface wetness fraction vector.
#' @param soilinit   Initial soil temperature matrix.
#' @param horizon_angles       Horizon angles vector (24 values).
#' @param aerosol_optical_depth        Solar attenuation vector.
#' @param soil_properties  Soil thermal properties matrix.
#' @param soil_moisture_profile     Initial soil moisture matrix.
#' @param rainfall   Daily rainfall vector (mm).
#' @param deep_soil_temperature   Deep soil boundary temperature vector/scalar.
#' @param air_entry_water_potential         Air entry potential vector (19 levels).
#' @param soil_bulk_density         Bulk density vector (19 levels).
#' @param soil_mineral_density         Mineral density vector (19 levels).
#' @param campbell_b_parameter         Soil b parameter vector (19 levels).
#' @param saturated_hydraulic_conductivity         Saturated conductivity vector (19 levels).
#' @param root_density          Root density vector (19 levels).
#' @param leaf_area_index        Leaf area index vector.
#' @param tides      Tidal state matrix.
#' @param reference_temperature     Hourly air temperature vector.
#' @param reference_humidity       Hourly relative humidity vector.
#' @param reference_wind_speed       Hourly wind speed vector.
#' @param cloud_cover      Hourly cloud cover vector.
#' @param global_radiation     Hourly solar radiation vector.
#' @param rainfall_hourly     Hourly rainfall vector.
#' @param zenith_angle_hourly      Hourly zenith angle vector.
#' @param longwave_radiation      Hourly longwave radiation vector.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @keywords internal
write_micro_csv <- function(micro_input, day_of_year, surface_emissivity, depths, soil_nodes,
                            maximum_shade_daily, minimum_shade_daily, maxima_times, minima_times,
                            reference_temperature_max, reference_temperature_min, reference_humidity_max, reference_humidity_min,
                            cloud_max, cloud_min, reference_wind_max, reference_wind_min,
                            albedo, soil_wetness, soilinit, horizon_angles, aerosol_optical_depth,
                            soil_properties, soil_moisture_profile, rainfall, deep_soil_temperature,
                            air_entry_water_potential, soil_bulk_density, soil_mineral_density, campbell_b_parameter, saturated_hydraulic_conductivity, root_density, leaf_area_index, tides,
                            reference_temperature, reference_humidity, reference_wind_speed, cloud_cover, global_radiation,
                            rainfall_hourly, zenith_angle_hourly, longwave_radiation) {
  if (!dir.exists("micro csv input")) {
    dir.create("micro csv input")
  }
  wt <- function(x, fname) {
    write.table(x, file = file.path("micro csv input", fname),
                sep = ",", col.names = NA, qmethod = "double")
  }
  wt(as.matrix(micro_input), "micro_input.csv")
  wt(day_of_year,        "day_of_year.csv")
  wt(surface_emissivity,       "surface_emissivity.csv")
  wt(depths,        "depths.csv")
  wt(soil_nodes,      "soil_nodes.csv")
  wt(maximum_shade_daily,  "Maxshades.csv")
  wt(minimum_shade_daily,  "Minshades.csv")
  wt(maxima_times,     "maxima_times.csv")
  wt(minima_times,     "minima_times.csv")
  wt(reference_temperature_max,      "reference_temperature_max.csv")
  wt(reference_temperature_min,      "reference_temperature_min.csv")
  wt(reference_humidity_max,     "reference_humidity_max.csv")
  wt(reference_humidity_min,     "reference_humidity_min.csv")
  wt(cloud_max,     "cloud_max.csv")
  wt(cloud_min,     "cloud_min.csv")
  wt(reference_wind_max,     "reference_wind_max.csv")
  wt(reference_wind_min,     "reference_wind_min.csv")
  wt(albedo,      "albedo.csv")
  wt(soil_wetness,     "soil_wetness.csv")
  wt(soilinit,   "soilinit.csv")
  wt(horizon_angles,       "horizon_angles.csv")
  wt(aerosol_optical_depth,        "aerosol_optical_depth.csv")
  wt(soil_properties,  "soilprop.csv")
  wt(soil_moisture_profile,     "soil_moisture_profile.csv")
  wt(rainfall,   "rain.csv")
  wt(deep_soil_temperature,   "deep_soil_temperature.csv")
  wt(air_entry_water_potential,         "air_entry_water_potential.csv")
  wt(soil_bulk_density,         "soil_bulk_density.csv")
  wt(soil_mineral_density,         "soil_mineral_density.csv")
  wt(campbell_b_parameter,         "campbell_b_parameter.csv")
  wt(saturated_hydraulic_conductivity,         "saturated_hydraulic_conductivity.csv")
  wt(root_density,          "root_density.csv")
  wt(leaf_area_index,        "leaf_area_index.csv")
  wt(tides,      "tides.csv")
  wt(reference_temperature,     "reference_temperature.csv")
  wt(reference_humidity,       "reference_humidity.csv")
  wt(reference_wind_speed,       "reference_wind_speed.csv")
  wt(cloud_cover,      "cloud_cover.csv")
  wt(global_radiation,     "global_radiation.csv")
  wt(rainfall_hourly,     "rainfall_hourly.csv")
  wt(zenith_angle_hourly,      "zenith_angle_hourly.csv")
  wt(longwave_radiation,      "longwave_radiation.csv")
  invisible(NULL)
}


# ---------------------------------------------------------------------------
# 7.  Process Fortran output
# ---------------------------------------------------------------------------

#' Unpack and post-process microclimate model output
#'
#' Extracts all output arrays from the list returned by \code{microclimate()},
#' fills in placeholder values when the soil moisture model is not run
#' (\code{soil_moisture_model = 0}), and optionally extracts snow and spectral outputs.
#'
#' @param microut   List returned by \code{microclimate()}.
#' @param soil_moisture_model  Flag: 1 = soil moisture model was run, 0 = it was not.
#' @param snow_model Flag: 1 = snow model was run.
#' @param radiation_per_wavelength      Flag: 1 = wavelength-specific solar output was requested.
#'
#' @return Named list containing all post-processed output arrays:
#'   \code{micromet_lowshade}, \code{micromet_highshade}, \code{soil}, \code{soil_temperature_highshade},
#'   \code{soil_moisture_lowshade}, \code{soil_moisture_highshade}, \code{humid}, \code{soil_humidity_highshade},
#'   \code{soil_water_potential_lowshade}, \code{soil_water_potential_highshade}, \code{plant}, \code{plant_output_highshade},
#'   \code{soil_conductivity_lowshade}, \code{soil_conductivity_highshade}, \code{soil_specific_heat_lowshade}, \code{soil_specific_heat_highshade},
#'   \code{soil_density_lowshade}, \code{soil_density_highshade}, and (conditionally) \code{snow_output_lowshade},
#'   \code{snow_output_highshade}, \code{direct_solar_spectrum}, \code{rayleigh_solar_spectrum}, \code{diffuse_solar_spectrum}.
#'
#' @keywords internal
process_micro_output <- function(microut, soil_moisture_model, snow_model, radiation_per_wavelength) {
  # --- Above-ground micromet conditions -------------------------------------
  micromet_lowshade    <- microut$micromet_lowshade    # min-shade above-ground conditions
  micromet_highshade   <- microut$micromet_highshade   # max-shade above-ground conditions

  # --- Soil temperatures ----------------------------------------------------
  soil_temperature_lowshade      <- microut$soil_temperature_lowshade      # min-shade soil_temperature_lowshade temperatures
  soil_temperature_highshade  <- microut$soil_temperature_highshade  # max-shade soil_temperature_lowshade temperatures

  # --- Soil thermal properties (always available) ---------------------------
  soil_conductivity_lowshade        <- microut$soil_conductivity_lowshade
  soil_conductivity_highshade    <- microut$soil_conductivity_highshade
  soil_specific_heat_lowshade     <- microut$soil_specific_heat_lowshade
  soil_specific_heat_highshade <- microut$soil_specific_heat_highshade
  soil_density_lowshade       <- microut$soil_density_lowshade
  soil_density_highshade   <- microut$soil_density_highshade

  # --- Soil moisture, humidity and plant outputs ----------------------------
  if (soil_moisture_model == 1) {
    soil_moisture_lowshade  <- microut$soil_moisture_lowshade
    soil_moisture_highshade  <- microut$soil_moisture_highshade
    soil_humidity_lowshade      <- microut$soil_humidity_lowshade
    soil_humidity_highshade  <- microut$soil_humidity_highshade
    soil_water_potential_lowshade    <- microut$soil_water_potential_lowshade
    soil_water_potential_highshade    <- microut$soil_water_potential_highshade
    plant_output_lowshade      <- microut$plant_output_lowshade
    plant_output_highshade  <- microut$plant_output_highshade
  } else {
    # Fill placeholders so downstream code always has these arrays available
    soil_water_potential_lowshade   <- soil_temperature_lowshade;    soil_water_potential_lowshade[, 3:12]   <- 0     # water potential = 0
    soil_water_potential_highshade   <- soil_temperature_lowshade;    soil_water_potential_highshade[, 3:12]   <- 0
    soil_moisture_lowshade <- soil_temperature_lowshade;    soil_moisture_lowshade[, 3:12] <- 0.5   # arbitrary mid-range moisture
    soil_moisture_highshade <- soil_temperature_lowshade;    soil_moisture_highshade[, 3:12] <- 0.5
    soil_humidity_lowshade     <- soil_temperature_lowshade;    soil_humidity_lowshade[, 3:12]     <- 0.99  # near-saturated humidity
    soil_humidity_highshade <- soil_temperature_lowshade;    soil_humidity_highshade[, 3:12] <- 0.99
    plant_output_lowshade     <- cbind(soil_temperature_lowshade, soil_temperature_lowshade[, 3:4]);  plant_output_lowshade[, 3:14]     <- 0
    plant_output_highshade <- cbind(soil_temperature_lowshade, soil_temperature_lowshade[, 3:4]);  plant_output_highshade[, 3:14] <- 0
  }

  # --- Snow outputs (only when snow_model = 1) --------------------------------
  snow_output_lowshade <- NULL
  snow_output_highshade <- NULL
  if (snow_model == 1) {
    snow_output_lowshade <- microut$snow_output_lowshade
    snow_output_highshade <- microut$snow_output_highshade
  }

  # --- Wavelength-specific solar outputs (only when radiation_per_wavelength = 1) ---------------
  direct_solar_spectrum  <- NULL
  rayleigh_solar_spectrum <- NULL
  diffuse_solar_spectrum  <- NULL
  if (radiation_per_wavelength == 1) {
    direct_solar_spectrum  <- as.data.frame(microut$direct_solar_spectrum)   # direct solar irradiance
    rayleigh_solar_spectrum <- as.data.frame(microut$rayleigh_solar_spectrum)  # direct Rayleigh component
    diffuse_solar_spectrum  <- as.data.frame(microut$diffuse_solar_spectrum)   # scattered solar irradiance
  }

  list(
    micromet_lowshade       = micromet_lowshade,       micromet_highshade       = micromet_highshade,
    soil_temperature_lowshade         = soil_temperature_lowshade,         soil_temperature_highshade      = soil_temperature_highshade,
    soil_moisture_lowshade    = soil_moisture_lowshade,    soil_moisture_highshade     = soil_moisture_highshade,
    soil_humidity_lowshade        = soil_humidity_lowshade,        soil_humidity_highshade     = soil_humidity_highshade,
    soil_water_potential_lowshade      = soil_water_potential_lowshade,      soil_water_potential_highshade       = soil_water_potential_highshade,
    plant_output_lowshade        = plant_output_lowshade,        plant_output_highshade     = plant_output_highshade,
    soil_conductivity_lowshade        = soil_conductivity_lowshade,        soil_conductivity_highshade     = soil_conductivity_highshade,
    soil_specific_heat_lowshade     = soil_specific_heat_lowshade,     soil_specific_heat_highshade  = soil_specific_heat_highshade,
    soil_density_lowshade       = soil_density_lowshade,       soil_density_highshade    = soil_density_highshade,
    snow_output_lowshade      = snow_output_lowshade,      snow_output_highshade       = snow_output_highshade,
    direct_solar_spectrum        = direct_solar_spectrum,        rayleigh_solar_spectrum        = rayleigh_solar_spectrum,
    diffuse_solar_spectrum        = diffuse_solar_spectrum
  )
}


# ---------------------------------------------------------------------------
# 8.  Build the final return list
# ---------------------------------------------------------------------------

#' Assemble the final return list for a micro_*() function
#'
#' Constructs the named list returned to the user when the full microclimate
#' model has been run (\code{runmicro = TRUE}).  Handles the four combinations
#' of \code{snow_model} and \code{radiation_per_wavelength} flags.  Dataset-specific extra outputs
#' (e.g., \code{SILO.data}, \code{SLOPE}, \code{microclima.out}) are appended
#' via the \code{extra} argument.
#'
#' @param out         Named list from \code{\link{process_micro_output}}.
#' @param rainfall    Daily rainfall vector (mm).
#' @param ndays       Total simulation days.
#' @param elevation        Site elevation (m).
#' @param surface_reflectivity        Soil reflectance used (scalar).
#' @param longlat     Numeric vector \code{c(longitude, latitude)}.
#' @param nyears      Number of simulated years.
#' @param timeinterval Number of time steps per year.
#' @param minimum_shade_daily   Daily minimum shade vector (\%).
#' @param maximum_shade_daily   Daily maximum shade vector (\%).
#' @param depths         Soil depth nodes vector (cm).
#' @param dates       POSIXct or numeric vector of hourly date/time values.
#' @param dates2      POSIXct or numeric vector of daily date values.
#' @param air_entry_water_potential          Air entry potential vector (19 levels).
#' @param soil_bulk_density          Bulk density vector (19 levels).
#' @param soil_mineral_density          Mineral density vector (19 levels).
#' @param campbell_b_parameter          Soil b parameter vector (19 levels).
#' @param saturated_hydraulic_conductivity          Saturated conductivity vector (19 levels).
#' @param dem         Digital elevation model raster (or \code{NA}).
#' @param diffuse_frac Hourly diffuse fraction of solar radiation (or \code{NA}).
#' @param snow_model   Flag: 1 = snow model was run.
#' @param radiation_per_wavelength        Flag: 1 = wavelength-specific solar output available.
#' @param extra       Optional named list of additional dataset-specific outputs
#'   to append to the return list (e.g., \code{list(SILO.data = ...)}).
#'
#' @return Named list of microclimate model outputs.
#'
#' @keywords internal
build_micro_return <- function(out, rainfall, ndays, elevation, surface_reflectivity, longlat,
                               nyears, timeinterval, minimum_shade_daily, maximum_shade_daily,
                               depths, dates, dates2, air_entry_water_potential, soil_bulk_density, soil_mineral_density, campbell_b_parameter, saturated_hydraulic_conductivity,
                               dem, diffuse_frac, snow_model, radiation_per_wavelength,
                               extra = list()) {
  # Core outputs common to all variants
  base <- list(
    soil_temperature_lowshade         = out$soil_temperature_lowshade,
    soil_temperature_highshade     = out$soil_temperature_highshade,
    micromet_lowshade       = out$micromet_lowshade,
    micromet_highshade      = out$micromet_highshade,
    soil_moisture_lowshade    = out$soil_moisture_lowshade,
    soil_moisture_highshade    = out$soil_moisture_highshade,
    soil_humidity_lowshade        = out$soil_humidity_lowshade,
    soil_humidity_highshade    = out$soil_humidity_highshade,
    soil_water_potential_lowshade      = out$soil_water_potential_lowshade,
    soil_water_potential_highshade      = out$soil_water_potential_highshade,
    plant_output_lowshade        = out$plant_output_lowshade,
    plant_output_highshade    = out$plant_output_highshade,
    soil_conductivity_lowshade        = out$soil_conductivity_lowshade,
    soil_conductivity_highshade    = out$soil_conductivity_highshade,
    soil_specific_heat_lowshade     = out$soil_specific_heat_lowshade,
    soil_specific_heat_highshade = out$soil_specific_heat_highshade,
    soil_density_lowshade       = out$soil_density_lowshade,
    soil_density_highshade   = out$soil_density_highshade,
    rainfall     = rainfall,
    ndays        = ndays,
    elevation         = elevation,
    surface_reflectivity         = surface_reflectivity[1],
    longlat      = c(longlat[1], longlat[2]),
    nyears       = nyears,
    timeinterval = timeinterval,
    minimum_shade     = minimum_shade_daily,
    maximum_shade     = maximum_shade_daily,
    depths          = depths,
    dates        = dates,
    dates2       = dates2,
    air_entry_water_potential           = air_entry_water_potential,
    soil_bulk_density           = soil_bulk_density,
    soil_mineral_density           = soil_mineral_density,
    campbell_b_parameter           = campbell_b_parameter,
    saturated_hydraulic_conductivity           = saturated_hydraulic_conductivity,
    dem          = dem,
    diffuse_frac = diffuse_frac
  )

  # Add snow outputs if the snow model was run
  if (snow_model == 1) {
    base$snow_output_lowshade <- out$snow_output_lowshade
    base$snow_output_highshade <- out$snow_output_highshade
  }

  # Add wavelength-specific solar outputs if requested
  if (radiation_per_wavelength == 1) {
    base$direct_solar_spectrum  <- out$direct_solar_spectrum
    base$rayleigh_solar_spectrum <- out$rayleigh_solar_spectrum
    base$diffuse_solar_spectrum  <- out$diffuse_solar_spectrum
  }

  # Append any dataset-specific extras (e.g. SILO.data, SLOPE, microclima.out)
  if (length(extra) > 0) {
    base <- c(base, extra)
  }

  return(base)
}


# ---------------------------------------------------------------------------
# 9.  GADS aerosol optical depth / aerosol_optical_depth computation
# ---------------------------------------------------------------------------

#' Default aerosol_optical_depth profile: Elterman (1970)
#'
#' Aerosol transmittance profile from Elterman, L. 1970. Vertical-attenuation
#' model with eight surface meteorological ranges 2 to 13 kilometers. U.S.
#' Airforce Cambridge Research Laboratory, Bedford, Mass.  Used as the fallback
#' in \code{micro_global} and \code{micro_terra} when \code{global_aerosol_database = 0}.
#'
#' @keywords internal
TAI_ELTERMAN <- c(0.42, 0.415, 0.412, 0.408, 0.404, 0.4, 0.395, 0.388, 0.379, 0.379, 0.379, 0.375, 0.365, 0.345, 0.314, 0.3, 0.288, 0.28, 0.273, 0.264, 0.258, 0.253, 0.248, 0.243, 0.236, 0.232, 0.227, 0.223, 0.217, 0.213, 0.21, 0.208, 0.205, 0.202, 0.201, 0.198, 0.195, 0.193, 0.191, 0.19, 0.188, 0.186, 0.184, 0.183, 0.182, 0.181, 0.178, 0.177, 0.176, 0.175, 0.175, 0.174, 0.173, 0.172, 0.171, 0.17, 0.169, 0.168, 0.167, 0.164, 0.163, 0.163, 0.162, 0.161, 0.161, 0.16, 0.159, 0.157, 0.156, 0.156, 0.155, 0.154, 0.153, 0.152, 0.15, 0.149, 0.146, 0.145, 0.142, 0.14, 0.139, 0.137, 0.135, 0.135, 0.133, 0.132, 0.131, 0.13, 0.13, 0.129, 0.129, 0.128, 0.128, 0.128, 0.127, 0.127, 0.126, 0.125, 0.124, 0.123, 0.121, 0.118, 0.117, 0.115, 0.113, 0.11, 0.108, 0.107, 0.105, 0.103, 0.1)

#' Default aerosol_optical_depth profile: Australia (Adelaide/Melbourne representative)
#'
#' Aerosol transmittance profile representative of southern Australia,
#' generated by GADS at Adelaide/Melbourne conditions.  Used as the fallback
#' in all dataset-specific \code{micro_*()} functions (other than
#' \code{micro_global} and \code{micro_terra}) when \code{global_aerosol_database = 0}.
#'
#' @keywords internal
TAI_AUSTRALIA <- c(0.0670358341290886, 0.0662612704779235, 0.065497075238002, 0.0647431301168489, 0.0639993178022531, 0.0632655219571553, 0.0625416272145492, 0.0611230843885423, 0.0597427855962549, 0.0583998423063099, 0.0570933810229656, 0.0558225431259535, 0.0545864847111214, 0.0533843764318805, 0.0522154033414562, 0.0499736739981675, 0.047855059159556, 0.0458535417401334, 0.0439633201842001, 0.0421788036108921, 0.0404946070106968, 0.0389055464934382, 0.0374066345877315, 0.0359930755919066, 0.0346602609764008, 0.0334037648376212, 0.0322193394032758, 0.0311029105891739, 0.0300505736074963, 0.0290585886265337, 0.0281233764818952, 0.0272415144391857, 0.0264097320081524, 0.0256249068083005, 0.0248840604859789, 0.0241843546829336, 0.0235230870563317, 0.0228976873502544, 0.0223057135186581, 0.0217448478998064, 0.0212128934421699, 0.0207077699817964, 0.0202275105711489, 0.0197702578594144, 0.0193342605242809, 0.0189178697551836, 0.0177713140039894, 0.0174187914242432, 0.0170790495503944, 0.0167509836728154, 0.0164335684174899, 0.0161258546410128, 0.0158269663770596, 0.0155360978343254, 0.0152525104459325, 0.0149755299703076, 0.0147045436435285, 0.0144389973831391, 0.0141783930434343, 0.0134220329447663, 0.0131772403830191, 0.0129356456025128, 0.0126970313213065, 0.0124612184223418, 0.0122280636204822, 0.01199745718102, 0.0115436048739351, 0.0110993711778668, 0.0108808815754663, 0.0106648652077878, 0.0104513876347606, 0.0102405315676965, 0.00982708969547694, 0.00962473896278535, 0.00903679230300494, 0.00884767454432418, 0.0083031278398166, 0.00796072474935954, 0.00755817587626185, 0.00718610751850881, 0.00704629977586921, 0.00684663903049612, 0.00654155580333479, 0.00642947339729728, 0.00627223096874308, 0.00603955966866779, 0.00580920937536261, 0.00568506186880564, 0.00563167068287251, 0.00556222005081865, 0.00550522989971023, 0.00547395763028062, 0.0054478983436216, 0.00541823364504573, 0.00539532163908382, 0.00539239864119488, 0.00541690124712384, 0.00551525885358836, 0.00564825853509463, 0.00577220185074264, 0.00584222986640171, 0.00581645238345584, 0.00566088137411449, 0.00535516862329704, 0.00489914757707667, 0.00432017939770409, 0.0036813032251836, 0.00309019064543606, 0.00270890436501562, 0.00276446109239711, 0.00356019862584603)

#' Compute aerosol transmittance (aerosol_optical_depth) via GADS or a default profile
#'
#' Runs the GADS aerosol model (Fortran or R version) to compute
#' wavelength-specific aerosol optical depth, fits a 6th-order polynomial, and
#' predicts aerosol_optical_depth at the 111 standard wavelength bands used by the microclimate
#' model.  When \code{global_aerosol_database = 0} the supplied \code{default_tai} vector is
#' returned unchanged.
#'
#' @param longlat     Numeric vector \code{c(longitude, latitude)}.
#' @param global_aerosol_database    Integer: \code{0} = use default profile, \code{1} = run
#'   Fortran GADS, \code{2} = run R GADS.
#' @param default_tai Numeric vector of length 111 returned when
#'   \code{global_aerosol_database = 0}.  Use \code{TAI_ELTERMAN} or \code{TAI_AUSTRALIA}.
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
#' @param depths Soil depth node vector (cm), length 10.
#'
#' @return Named list with elements \code{air_entry_water_potential}, \code{saturated_hydraulic_conductivity}, \code{campbell_b_parameter},
#'   \code{soil_bulk_density}, \code{bulk_density} when data are available, or \code{NULL}
#'   when SoilGrids returns no data for the location.
#'
#' @keywords internal
fetch_soilgrids <- function(x, depths) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("package 'jsonlite' is needed to extract data from SoilGrids, please install it.",
         call. = FALSE)
  }
  message("extracting soil_temperature_lowshade texture data from SoilGrids")
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
    soil.hydro <- pedotransfer(soilpro = as.data.frame(soilpro), depths = depths)
    return(list(
      air_entry_water_potential          = soil.hydro$air_entry_water_potential,
      saturated_hydraulic_conductivity          = soil.hydro$saturated_hydraulic_conductivity,
      campbell_b_parameter          = soil.hydro$campbell_b_parameter,
      soil_bulk_density          = soil.hydro$soil_bulk_density,
      bulk_density = soil.hydro$soil_bulk_density[seq(1, 19, 2)]
    ))
  } else {
    message("no SoilGrids data for this site, using user-input soil_temperature_lowshade properties")
    return(NULL)
  }
}


compute_tai <- function(longlat, global_aerosol_database, default_tai) {
  if (global_aerosol_database > 0) {
    relhum <- 1
    if (global_aerosol_database == 1) {
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
#' @param dem_resolution DEM resolution (m) passed to \code{microclima::get_dem()}.
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
#' @return Named list: \code{dem} (SpatRaster), \code{elevation} (numeric),
#'   \code{slope}, \code{aspect}, \code{horizon_angles} (\code{nangles}-element vector or
#'   \code{NA}).
#' @keywords internal
fetch_dem <- function(lon, lat, dem_resolution = 10, zmin = NULL, pixels = 100,
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
    if (!is.null(dem_resolution)) dem_args$resolution <- dem_resolution
    if (!is.null(zmin)) dem_args$zmin <- zmin
    if (!is.null(pixels)) { dem_args$xdims <- pixels; dem_args$ydims <- pixels }
    existing_dem <- do.call(microclima::get_dem, dem_args)
  }
  dem <- terra::rast(existing_dem)  # normalise to SpatRaster (no-op if already SpatRaster)
  xy <- data.frame(lon = lon, lat = lat) |>
    sf::st_as_sf(coords = c("lon", "lat"))
  xy <- sf::st_set_crs(xy, "EPSG:4326")
  xy <- sf::st_transform(xy, sf::st_crs(dem))
  elevation <- as.numeric(terra::extract(dem, xy)[, 2])
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
  list(dem = dem, elevation = elevation, slope = slope, aspect = aspect, horizon_angles = hori_out)
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
#' @param leaf_area_index Numeric vector; leaf area index. \code{mean(leaf_area_index)} is passed to
#'   \code{.shortwave.ts()} as \code{l}.
#' @return Named list: \code{SOLRhr_all} (terrain-adjusted hourly solar
#'   radiation, W m\eqn{^{-2}}) and \code{diffuse_frac_all} (hourly diffuse
#'   fraction, 0–1).
#' @keywords internal
compute_solar_partition <- function(dsw2, jd, hour.microclima, lat, long,
                                    slope, aspect, ha, LOR, leaf_area_index) {
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
                            ha = ha, svv = 1, x = LOR, l = mean(leaf_area_index),
                            albr = 0, merid = long, dst = 0, difani = FALSE)
  SOLRhr_all <- radwind2$swrad / 0.0036
  list(SOLRhr_all = SOLRhr_all, diffuse_frac_all = diffuse_frac_all)
}
