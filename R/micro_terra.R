#' Global historical monthly implementation of the microclimate model
#'
#' An implementation of the NicheMapR microclimate model that uses the global monthly climate database
#' described in "Abatzoglou, J. T., Dobrowski, S. Z., Parks, S. A., & Hegewisch, K. C. (2018). TerraClimate,
#' a high-resolution global dataset of monthly climate and climatic water balance from 1958–2023.
#' Scientific Data, 5(1), 170191. https://doi.org/10.1038/sdata.2017.191"
#' This dataset includes climate change scenarios for +2 and +4 deg C warming for 1950-2025 (via parameter 'scenario')
#' Aerosol attenuation can also be computed based on the Global Aerosol Data Set (GADS)
#' Koepke, P., M. Hess, I. Schult, and E. P. Shettle. 1997. Global Aerosol Data Set. Max-Planck-Institut for Meteorologie, Hamburg
#' by choosing the option 'global_aerosol_database <- 1' (Fortran version, quicker but may crash on some systems) or 'global_aerosol_database <- 2' (R version)
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param timeinterval The number of time intervals to generate predictions for over a year (must be 12 <= x <=365)
#' @param ystart First year to run
#' @param yfinish Last year to run
#' @param surface_reflectivity Soil solar reflectance, decimal \%
#' @param elevation Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param depths Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param minimum_shade Minimum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param maximum_shade Maximum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param local_height Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
#' @param ... Additional arguments, see Details
#' @usage micro_terra(loc = c(-89.40123, 43.07305), timeinterval = 12, ystart = 2000, yfinish = 2015,
#' surface_reflectivity = 0.15, slope = 0, aspect = 0,
#' depths = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minimum_shade = 0, maximum_shade = 90,
#' local_height = 0.01, ...)
#' @return micromet_lowshade The above ground micrometeorological conditions under the minimum specified shade
#' @return micromet_highshade The above ground micrometeorological conditions under the maximum specified shade
#' @return soil_temperature_lowshade Hourly predictions of the soil temperatures under the minimum specified shade
#' @return soil_temperature_highshade Hourly predictions of the soil temperatures under the maximum specified shade
#' @return soil_moisture_lowshade Hourly predictions of the soil moisture under the minimum specified shade
#' @return soil_moisture_highshade Hourly predictions of the soil moisture under the maximum specified shade
#' @return soil_water_potential_lowshade Hourly predictions of the soil water potential under the minimum specified shade
#' @return soil_water_potential_highshade Hourly predictions of the soil water potential under the maximum specified shade
#' @return soil_humidity_lowshade Hourly predictions of the soil humidity under the minimum specified shade
#' @return soil_humidity_highshade Hourly predictions of the soil humidity under the maximum specified shade
#' @return plant_output_lowshade Hourly predictions of plant transpiration, leaf water potential and root water potential under the minimum specified shade
#' @return plant_output_highshade Hourly predictions of plant transpiration, leaf water potential and root water potential under the maximum specified shade
#' @return snow_output_lowshade Hourly predictions of snow temperature under the minimum specified shade
#' @return snow_output_highshade Hourly predictions snow temperature under the maximum specified shade
#' @return soil_conductivity_lowshade Hourly predictions of the soil thermal conductivity under the minimum specified shade
#' @return soil_conductivity_highshade Hourly predictions of the soil thermal conductivity under the maximum specified shade
#' @return soil_specific_heat_lowshade Hourly predictions of the soil specific heat capacity under the minimum specified shade
#' @return soil_specific_heat_highshade Hourly predictions of soil specific heat capacity under the maximum specified shade
#' @return soil_density_lowshade Hourly predictions of the soil density under the minimum specified shade
#' @return soil_density_highshade Hourly predictions of the soil density under the maximum specified shade
#' @export
#' @details
#' \itemize{
#' \strong{Parameters controlling how the model runs:}\cr\cr
#'
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{global_aerosol_database}{ = 1, Use the Global Aerosol Database? 1=yes (Fortran version), 2=yes (R version), 0=no}\cr\cr
#' \code{longwave_radiation_model}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{solar_model_only}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{radiation_per_wavelength}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{scattered_uv}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{max_iterations_per_day}{ = 3, iterations per day to get a steady periodic}\cr\cr
#' \code{initial_soil_temperature}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{output_to_csv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{elevatr}{ = 0, Use elevatr package to get high resolution elevation for location? 1 = yes, 0 = no}\cr\cr
#' \code{terrain}{ = 0, Use elevatr package to adjust horizon angles, slope and aspect? 1 = yes, 0 = no}\cr\cr
#' \code{microclima}{ = 0, Use microclima and elevatr package to compute diffuse fraction of solar radiation (1) and adjust solar radiation for terrain (2)? 0 = no}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#' \code{runmicro}{ = 1, call the microclimate model (1) or not (0), if you just want the downscaled input weather data}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#' \code{tolerance}{ = 1, Integrator error tolerance for soil temperature calculations}\cr\cr
#' \code{roughness_height}{ = 0.004, Roughness height (m), e.g. smooth desert is 0.0003, closely mowed grass may be 0.001, bare tilled soil 0.002-0.006, current allowed range: 0.00001 (snow) - 0.02 m.}\cr\cr
#' \code{canopy_roughness_height}{ = 0, heat transfer roughness height (m) for Campbell and Norman air temperature/wind speed profile (invoked if greater than 0, 0.02 * canopy height in m if unknown)}\cr\cr
#' \code{zero_plane_displacement}{ = 0, zero plane displacement correction factor (m) for Campbell and Norman air temperature/wind speed profile (0.6 * canopy height in m if unknown)}\cr\cr
#' \code{roughness_height_1}{ = 0, Top (1st) segment roughness height(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then roughness_height and reference_height used)}\cr\cr
#' \code{roughness_height_2}{ = 0, 2nd segment roughness height(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then roughness_height and reference_height used).}\cr\cr
#' \code{wind_profile_height_1}{ = 0, Top of (1st) segment, height above surface(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then roughness_height and reference_height used).}\cr\cr
#' \code{wind_profile_height_2}{ = 0, 2nd segment, height above surface(m) - IF NO EXPERIMENTAL WIND PROFILE DATA SET THIS TO ZERO! (then roughness_height and reference_height used).}\cr\cr
#' \code{orbital_eccentricity}{ = 0.0167238, Eccenricity of the earth's orbit (current value 0.0167238, ranges between 0.0034 to 0.058)}\cr\cr
#' \code{surface_emissivity}{ = 0.95, Substrate longwave IR emissivity (decimal \%), typically close to 1}\cr\cr
#' \code{mineral_conductivity}{ = 2.5, Soil minerals thermal conductivity, single value or vector of 10 specific to each depth (W/mK)}\cr\cr
#' \code{mineral_density}{ = 2.56, Soil minerals density, single value or vector of 10 specific to each depth (Mg/m3)}\cr\cr
#' \code{mineral_heat_capacity}{ = 870, Soil minerals specific heat, single value or vector of 10 specific to each depth (J/kg-K)}\cr\cr
#' \code{bulk_density}{ = 1.3, Soil bulk density (Mg/m3), single value or vector of 10 specific to each depth}\cr\cr
#' \code{soil_wetness}{ = 0, \% of ground surface area acting as a free water surface (overridden if soil moisture model is running)}\cr\cr
#' \code{organic_soil_cap}{ = 1, organic cap present on soil surface? (organic_soil_cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)}\cr\cr
#' \code{precipitable_water}{ = 1, Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air conditions; 2.0 = humid, tropical conditions (note this is for the whole atmospheric profile, not just near the ground)}\cr\cr
#' \code{horizon_angles}{ = rep(0,24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#' \code{minimum_temperature_lapse_rate}{ = 0.0039 Lapse rate for minimum air temperature (degrees C/m)}\cr\cr
#' \code{maximum_temperature_lapse_rate}{ = 0.0077 Lapse rate for maximum air temperature (degrees C/m)}\cr\cr
#' \code{maxima_times}{ = c(1, 1, 0, 0), Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise}\cr\cr
#' \code{minima_times}{ = c(0, 0, 1, 1), Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon}\cr\cr
#' \code{timezone}{ = 0, Use GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)? 1=yes, 0=no}\cr\cr
#' \code{aerosol_optical_depth}{ = 0, Vector of 111 values, one per wavelenght bin, for solar attenuation - used to overide GADS}\cr\cr
#' \code{wind_multiplier}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{air_temperature_offset}{ = 0, warming offset vector, °C (negative values mean cooling). Can supply a single value or a vector the length of the number of days to be simulated.}\cr\cr
#' \code{terra_source}{ = NA, specify location of terraclimate data, goes to the web by default}\cr\cr
#' \code{scenario}{ = 0, TerraClimate climate change scenario, either 0, 2 or 4 °C warmer}\cr\cr
#'
#' \strong{ Soil moisture mode parameters:}
#'
#' \code{soil_moisture_model}{ = 0, Run soil moisture model? 1=yes, 0=no  1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{air_entry_water_potential}{ = rep(1.1,19), Air entry potential (J/kg) (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{saturated_hydraulic_conductivity}{ = rep(0.0037,19), Saturated conductivity, (kg s/m3) (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{campbell_b_parameter}{ = rep(4.5,19), Campbell's soil 'b' parameter (-) (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{soil_bulk_density}{ = rep(1.3,19), Soil bulk density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter depths and points half way between)}\cr\cr
#' \code{soil_mineral_density}{ = rep(2.56,19), Soil density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter depths and points half way between)}\cr\cr
#' \code{maximum_pooling_depth}{ = 10000, Max depth for water pooling on the surface (mm), to account for runoff}\cr\cr
#' \code{rain_multiplier}{ = 1, Rain multiplier for surface soil moisture (-), used to induce runon}\cr\cr
#' \code{evenrain}{ = 0, Spread daily rainfall evenly across 24hrs (1) or one event at midnight (0)}\cr\cr
#' \code{initial_soil_moisture}{ = c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3), initial soil water content at each soil node, m3/m3}\cr\cr
#' \code{root_density}{ = c(0,0,8.2,8.0,7.8,7.4,7.1,6.4,5.8,4.8,4.0,1.8,0.9,0.6,0.8,0.4,0.4,0,0)*10000, root density (m/m3), (19 values descending through soil for specified soil nodes in parameter}\cr\cr
#' \code{root_radius}{ = 0.001, root radius, m}\cr\cr
#' \code{root_resistance}{ = 2.5e+10, resistance per unit length of root, m3 kg-1 s-1}\cr\cr
#' \code{leaf_resistance}{ = 2e+6, resistance per unit length of leaf, m3 kg-1 s-1}\cr\cr
#' \code{stomatal_closure_potential}{ = -1500, critical leaf water potential for stomatal closure, J kg-1}\cr\cr
#' \code{stomatal_stability_parameter}{ = 10, stability parameter for stomatal closure equation, -}\cr\cr
#' \code{moist_error}{ = 1e-06, maximum allowable mass balance error, kg}\cr\cr
#' \code{moist_count}{ = 500, maximum iterations for mass balance, -}\cr\cr
#' \code{leaf_area_index}{ = 0.1, leaf area index (can be a single value or a vector of daily values), used to partition traspiration/evaporation from PET}\cr\cr
#' \code{microclima.leaf_area_index}{ = 0, leaf area index, used by package microclima for radiation calcs}\cr\cr
#' \code{microclima.LOR}{ = 1, leaf orientation for package microclima radiation calcs}\cr\cr
#'
#' \strong{ Snow mode parameters:}
#'
#' \code{snow_model}{ = 1, run the snow model 1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (Mg/m3), overridden by densfun}\cr\cr
#' \code{densfun}{ = c(0.5979, 0.2178, 0.001, 0.0038), slope and intercept of model of snow density as a linear function of snowpack age if first two values are nonzero, and following the exponential function of Sturm et al. 2010 J. of Hydromet. 11:1380-1394 if all values are non-zero; if it is c(0,0,0,0) then fixed density used}\cr\cr
#' \code{snowmelt}{ = 1, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
#' \code{undercatch}{ = 1, undercatch multipier for converting rainfall to snow}\cr\cr
#' \code{rainmelt}{ = 0.0125, paramter in equation that melts snow with rainfall as a function of air temp}\cr\cr
#' \code{rainfrac}{ = 0.5, fraction of rain that falls on the first day of the month (decimal \% with 0 meaning rain falls evenly) - this parameter allows something other than an even intensity of rainfall when interpolating the montly rainfall data)}\cr\cr
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)}\cr\cr
#' \code{intercept}{ = max(maximum_shade) / 100 * 0.3, snow interception fraction for when there's shade (0-1)}\cr\cr
#' \code{grasshade}{ = 0, if 1, means shade is removed when snow is present, because shade is cast by grass/low shrubs}\cr\cr
#'
#' \strong{ Intertidal mode parameters:}
#'
#' \code{shore}{ Include tide effects? If 1, the matrix}
#' \code{tides}
#' { is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0., nrow = 24*365*nyears, ncol = 3), matrix of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#' }
#'
#' \strong{Outputs:}
#'
#' \code{ndays}{ - number of days for which predictions are made}\cr\cr
#' \code{longlat}{ - longitude and latitude for which simulation was run (decimal degrees)}\cr\cr
#' \code{dates}{ - vector of dates (hourly, POSIXct, timezone = GMT+10)}\cr\cr
#' \code{dates2}{ - vector of dates (daily, POSIXct, timezone = GMT+10)}\cr\cr
#' \code{nyears}{ - number of years for which predictions are made}\cr\cr
#' \code{rainfall}{ - vector of daily rainfall (mm)}\cr\cr
#' \code{elevation}{ - elevation at point of simulation (m)}\cr\cr
#' \code{minimum_shade}{ - minimum shade for each day of simulation (\%)}\cr\cr
#' \code{maximum_shade}{ - maximum shade for each day of simulation (\%)}\cr\cr
#' \code{dem}{ - digital elevation model obtained via 'get_dem' using package 'elevatr' (m)}\cr\cr
#' \code{depths}{ - vector of depths used (cm)}\cr\cr
#' \code{diffuse_frac}{ - vector of hourly values of the fraction of total solar radiation that is diffuse (-), computed by microclima if microclima > 0}\cr\cr
#'
#' micromet_lowshade/micromet_highshade variables:
#' \itemize{
#' \item 1 day_of_year - day-of-year
#' \item 2 time - time of day (mins)
#' \item 3 air_temperature_local - air temperature (°C) at local height (specified by 'local_height' variable)
#' \item 4 air_temperature_reference - air temperature (°C) at reference height (specified by 'reference_height', 2m default)
#' \item 5 relative_humidity_local - relative humidity (\%) at local height (specified by 'local_height' variable)
#' \item 6 relative_humidity_reference  - relative humidity (\%) at reference height (specified by 'reference_height', 2m default)
#' \item 7 wind_speed_local - wind speed (m/s) at local height (specified by 'local_height' variable)
#' \item 8 wind_speed_reference - wind speed (m/s) at reference height (specified by 'reference_height', 2m default)
#' \item 9 snow_melt - snowmelt (mm)
#' \item 10 pool_depth - water pooling on surface (mm)
#' \item 11 soil_wetness - soil surface wetness (\%)
#' \item 12 zenith_angle - zenith angle of sun (degrees - 90 = below the horizon)
#' \item 13 solar_radiation - solar radiation (W/m2) (unshaded, horizontal plane)
#' \item 14 sky_temperature - sky radiant temperature (°C)
#' \item 15 dew - dew fall (mm / h)
#' \item 16 frost - frost (mm / h)
#' \item 17 snow_fall - snow predicted to have fallen (cm)
#' \item 18 snow_depth - predicted snow depth (cm)
#' \item 19 snow_density - snow density (g/cm3)
#' }
#' soil_temperature_lowshade and soil_temperature_highshade variables:
#' \itemize{
#' \item 1 day_of_year - day-of-year
#' \item 2 time - time of day (mins)
#' \item 3-12 depth_0cm ... - soil temperature (°C) at each of the 10 specified depths
#' }
#'
#' if soil moisture model is run i.e. parameter soil_moisture_model = 1\cr
#'
#' soil_moisture_lowshade and soil_moisture_highshade variables:
#' \itemize{
#' \item 1 day_of_year - day-of-year
#' \item 2 time - time of day (mins)
#' \item 3-12 depth_0cm ... - soil moisture (m3/m3) at each of the 10 specified depths
#' }
#' soil_water_potential_lowshade and soil_water_potential_highshade variables:
#' \itemize{
#' \item 1 day_of_year - day-of-year
#' \item 2 time - time of day (mins)
#' \item 3-12 depth_0cm ... - soil water potential (J/kg = kPa = bar/100) at each of the 10 specified depths
#' }
#' soil_humidity_lowshade and soil_humidity_highshade variables:
#' \itemize{
#' \item  1 day_of_year - day-of-year
#' \item  2 time - time of day (mins)
#' \item  3-12 depth_0cm ... - soil relative humidity (decimal \%), at each of the 10 specified depths
#' }
#' plant_output_lowshade and plant_output_highshade variables:
#' \itemize{
#' \item  1 day_of_year - day-of-year
#' \item  2 time - time of day (mins)
#' \item  3 transpiration - plant transpiration rate (g/m2/h)
#' \item  4 leaf_water_potential - leaf water potential (J/kg = kPa = bar/100)
#' \item  5-14 depth_0cm ... - root water potential (J/kg = kPa = bar/100), at each of the 10 specified depths
#' }
#' soil_conductivity_lowshade and soil_conductivity_highshade variables:
#' \itemize{
#' \item  1 day_of_year - day-of-year
#' \item  2 time - time of day (mins)
#' \item  3-12 depth_0cm ... - soil thermal conductivity (W/m-K), at each of the 10 specified depths
#' }
#' soil_specific_heat_lowshade and soil_specific_heat_highshade variables:
#' \itemize{
#' \item  1 day_of_year - day-of-year
#' \item  2 time - time of day (mins)
#' \item  3-12 SP_depth_0cm ... - soil specific heat capacity (J/kg-K), at each of the 10 specified depths
#' }
#' soil_density_lowshade and soil_density_highshade variables:
#' \itemize{
#' \item  1 day_of_year - day-of-year
#' \item  2 time - time of day (mins)
#' \item  3-12 depth_0cm ... - soil density (Mg/m3), at each of the 10 specified depths
#' }
#'
#' if snow model is run i.e. parameter snow_model = 1\cr
#' snow_output_lowshade and snow_output_highshade variables:
#' \itemize{
#' \item  1 day_of_year - day-of-year
#' \item  2 time - time of day (mins)
#' \item  3-10 SN1 ... - snow temperature (°C), at each of the potential 8 snow layers (layer 8 is always the bottom - need micromet_lowshade$snow_depth to interpret which depth in the snow a given layer represents)
#' }
#'
#' if wavelength-specific solar output is selected i.e. parameter radiation_per_wavelength = 1\cr
#' solar output variables
#' direct_solar_spectrum (direct solar), rayleigh_solar_spectrum (direct Rayleigh solar) and diffuse_solar_spectrum (scattered solar) variables:
#' \itemize{
#' \item  1 day_of_year - day-of-year
#' \item  2 time - time of day (mins)
#' \item  3-113 290, ..., 4000 - irradiance (W/(m2 nm)) at each of 111 wavelengths from 290 to 4000 nm
#' }
#' @examples
#'library(geodata)
#'elevation <- geodata::worldclim_global(var = 'elevation', res = 2.5, path=tempdir()) # download DEM in advance (needed for lapse rate correction - automatically downloaded otherwise)
#'micro <- micro_terra(elevation=elevation) # run the model with default location (Madison, Wisconsin) from 2000 to 2015 and settings
#'
#'micromet_lowshade <- as.data.frame(micro$micromet_lowshade) # above ground microclimatic conditions, min shade
#'micromet_highshade <- as.data.frame(micro$micromet_highshade) # above ground microclimatic conditions, max shade
#'soil_temperature_lowshade <- as.data.frame(micro$soil_temperature_lowshade) # soil temperatures, minimum shade
#'soil_temperature_highshade <- as.data.frame(micro$soil_temperature_highshade) # soil temperatures, maximum shade
#'
#'minimum_shade <- micro$minimum_shade[1]
#'maximum_shade <- micro$maximum_shade[1]
#'
#'# plotting above-ground conditions in minimum shade
#'with(micromet_lowshade, {plot(air_temperature_local ~ micro$dates3, xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l", main = paste("air temperature, ", minimum_shade, "% shade",sep = ""))})
#'with(micromet_lowshade, {points(air_temperature_reference ~ micro$dates3, xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l",lty = 2, col = 'blue')})
#'with(micromet_lowshade, {plot(relative_humidity_local ~ micro$dates3, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",ylim = c(0, 100),main = paste("humidity, ",minimum_shade, "% shade",sep=""))})
#'with(micromet_lowshade, {points(relative_humidity_reference ~ micro$dates3, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",col = 'blue',lty = 2, ylim = c(0, 100))})
#'with(micromet_lowshade, {plot(sky_temperature ~ micro$date3s, xlab = "Date and Time", ylab = "Sky Temperature (°C)"
#',  type = "l", main = paste("sky temperature, ", minimum_shade, "% shade", sep=""))})
#'with(micromet_lowshade, {plot(wind_speed_reference ~ micro$dates3, xlab = "Date and Time",  ylab = "Wind Speed (m/s)"
#',  type = "l", main = "wind speed", col = 'blue',ylim = c(0, 15))})
#'with(micromet_lowshade, {points(wind_speed_local ~ micro$dates3, xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#',  type = "l", lty = 2)})
#'with(micromet_lowshade, {plot(zenith_angle ~ micro$dates3,xlab = "Date and Time", ylab = "Zenith Angle of Sun (deg)"
#',  type = "l", main = "solar angle, sun")})
#'with(micromet_lowshade, {plot(solar_radiation ~ micro$dates3,xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
#',  type = "l", main = "solar radiation")})
#'
#'# plotting soil temperature for minimum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(soil_temperature_lowshade[,i + 2] ~ micro$dates3, xlab = "Date and Time", ylab = "Soil Temperature (°C)"
#'    ,col = i, type = "l", main = paste("soil temperature ", minimum_shade, "% shade", sep=""))
#'  }else{
#'    points(soil_temperature_lowshade[,i + 2] ~ micro$dates3, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l")
#'  }
#'}
#'points(micromet_lowshade$snow_depth ~ micro$dates, type = 'h', col = 'light blue')
#
#'# plotting above-ground conditions in maximum shade
#'with(micromet_highshade,{plot(air_temperature_local ~ micro$dates3,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l", main = "air temperature, sun")})
#'with(micromet_highshade,{points(air_temperature_reference ~ micro$dates3,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l", lty = 2, col = 'blue')})
#'with(micromet_highshade,{plot(relative_humidity_local ~ micro$dates3,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l", ylim = c(0, 100),main = "humidity, shade")})
#'with(micromet_highshade,{points(relative_humidity_reference ~ micro$dates3,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l", col = 'blue',lty = 2, ylim = c(0, 100))})
#'with(micromet_highshade,{plot(sky_temperature ~ micro$dates3,xlab = "Date and Time", ylab = "Sky Temperature (°C)",
#'  type = "l", main = "sky temperature, shade")})
#'
#'# plotting soil temperature for maximum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(soil_temperature_highshade[,i + 2] ~ micro$dates3, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l", main = paste("soil temperature ", maximum_shade, "% shade", sep=""))
#'  }else{
#'    points(soil_temperature_highshade[,i + 2] ~ micro$dates3, xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l")
#'  }
#'}
#'points(micromet_highshade$snow_depth ~ micro$dates, type = 'h', col = 'light blue')
micro_terra <- function(
  loc = c(-89.4557, 43.1379),
  ystart = 2000,
  yfinish = 2015,
  timeinterval = 12,
  elevation = NA,
  nyears = yfinish - ystart + 1,
  surface_reflectivity = 0.15,
  slope = 0,
  aspect = 0,
  maximum_temperature_lapse_rate = 0.0077,
  minimum_temperature_lapse_rate = 0.0039,
  depths = c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
  minimum_shade = 0,
  maximum_shade = 90,
  dem = NA,
  local_height = 0.01,
  roughness_height_1 = 0,
  roughness_height_2 = 0,
  wind_profile_height_1 = 0,
  wind_profile_height_2 = 0,
  runshade = 1,
  solar_model_only = 0,
  clearsky = 0,
  global_aerosol_database = 1,
  initial_soil_temperature = NA,
  write_input = 0,
  output_to_csv = 0,
  elevatr = 0,
  terrain = 0,
  microclima = 0,
  microclima.dem_resolution = 100,
  microclima.zmin = -20,
  adiab_cor = 1,
  air_temperature_offset = 0,
  tolerance = 1,
  roughness_height = 0.004,
  canopy_roughness_height = 0,
  zero_plane_displacement = 0,
  orbital_eccentricity = 0.0167238,
  surface_emissivity = 0.95,
  mineral_conductivity = 2.5,
  mineral_density = 2.56,
  mineral_heat_capacity = 870,
  bulk_density = 1.3,
  soil_wetness = 0,
  organic_soil_cap = 1,
  precipitable_water = 1,
  horizon_angles = rep(0,24),
  maxima_times = c(1, 1, 0, 0),
  minima_times = c(0, 0, 1, 1),
  timezone = 0,
  soil_moisture_model = 0,
  air_entry_water_potential = rep(1.1, 19),
  saturated_hydraulic_conductivity = rep(0.0037, 19),
  campbell_b_parameter = rep(4.5, 19),
  soil_bulk_density = rep(bulk_density, 19),
  soil_mineral_density = rep(mineral_density, 19),
  maximum_pooling_depth = 10000,
  rain_multiplier = 1,
  evenrain = 0,
  initial_soil_moisture = c(0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.3, 0.3, 0.3, 0.3),
  root_density = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4 ,0.4, 0, 0) * 10000,
  root_radius = 0.001,
  root_resistance = 2.5e+10,
  leaf_resistance = 2e+6,
  stomatal_closure_potential = -1500,
  stomatal_stability_parameter = 10,
  moist_error = 1e-06,
  moist_count = 500,
  leaf_area_index = 0.1,
  microclima.leaf_area_index = 0,
  microclima.LOR = 1,
  snow_model = 1,
  snowtemp = 1.5,
  snowdens = 0.375,
  densfun = c(0.5979, 0.2178, 0.001, 0.0038),
  snowmelt = 1,
  undercatch = 1,
  rainmelt = 0.0125,
  rainfrac = 0.5,
  shore = 0,
  tides = 0,
  deep_soil_temperature = NA,
  radiation_per_wavelength = 0,
  scattered_uv = 0,
  max_iterations_per_day = 3,
  soilgrids = 0,
  longwave_radiation_model = 0,
  message = 0,
  fail = nyears * 24 * 365,
  runmicro = 1,
  aerosol_optical_depth = 0,
  wind_multiplier = 1,
  snowcond = 0,
  intercept = max(maximum_shade) / 100 * 0.3,
  grasshade = 0,
  terra_source = NA,
  scenario = 0,
  maxsurf = 85
) {

  # Reference height (m) at which input air temperature, wind speed and
  # relative humidity are measured (TerraClimate data are at 2 m height)
  reference_height <- 2

  # --- 1. Input validation --------------------------------------------------
  # Dataset-specific checks first
  errors <- 0
  if (ystart < 1950) {
    message("ERROR: TerraClimate climate data is not available prior to 1950.")
    errors <- 1
  }
  curdate <- Sys.time() - 60 * 60 * 24
  curyear <- as.numeric(format(curdate, "%Y"))
  if (ystart > curyear - 1) {
    message(paste0("warning: TerraClimate climate data is only available until ", curyear - 1, "."))
  }
  if (!(timeinterval %in% c(12, 365))) {
    message("ERROR: 'timeinterval' must be 12 or 365 for TerraClimate.")
    errors <- 1
  }
  if (is.numeric(loc[1])) {
    if (loc[1] > 180 | loc[2] > 90) {
      message("ERROR: Latitude or longitude (loc) is out of bounds. Please enter correct values.")
      errors <- 1
    }
  }
  if (timezone %in% c(0, 1) == FALSE) {
    message("ERROR: 'timezone' must be either 0 or 1. Please correct.")
    errors <- 1
  }
  if (global_aerosol_database == 1) {
    message("If program is crashing, try global_aerosol_database = 2.")
  }
  if (write_input %in% c(0, 1) == FALSE) {
    message("ERROR: 'write_input' must be 0 or 1. Please correct.")
    errors <- 1
  }
  if (!(scenario %in% c(0, 2, 4))) {
    message("ERROR: 'scenario' can only be 0, 2, or 4 (TerraClimate climate change scenarios).")
    errors <- 1
  }
  if (scenario > 0 & (yfinish > 2025 | ystart < 1950)) {
    message("ERROR: TerraClimate climate change scenarios are only for years 1950 to 2025.")
    errors <- 1
  }
  # Shared checks common to all micro_*() functions
  errors <- errors + validate_micro_inputs(
    depths = depths, surface_reflectivity = surface_reflectivity, slope = slope, aspect = aspect, horizon_angles = horizon_angles,
    surface_emissivity = surface_emissivity, tolerance = tolerance, roughness_height = roughness_height, zero_plane_displacement = zero_plane_displacement, local_height = local_height,
    reference_height = reference_height, orbital_eccentricity = orbital_eccentricity, precipitable_water = precipitable_water,
    maxima_times = maxima_times, minima_times = minima_times,
    minimum_shade = minimum_shade, maximum_shade = maximum_shade,
    mineral_conductivity = mineral_conductivity, mineral_density = mineral_density, mineral_heat_capacity = mineral_heat_capacity,
    bulk_density = bulk_density,
    global_aerosol_database = global_aerosol_database
  )

  if(errors == 0){ # continue

    ################## time related variables #################################
    tzone <- paste("Etc/GMT+", 10, sep = "")
    dates <- seq(ISOdate(ystart, 1, 1, tz = tzone) - 3600 * 12, ISOdate((ystart + nyears), 1, 1, tz = tzone) - 3600 * 13, by = "days")
    if(timeinterval == 365){
      ndays <- length(dates)
    }else{
      ndays <- 12 * nyears
    }
    doys12 <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month
    # Expand scalar or short shade inputs to one value per simulation day
    shades <- setup_shade_vectors(minimum_shade, maximum_shade, ndays)
    minimum_shade_daily <- shades$minimum_shade_daily
    maximum_shade_daily <- shades$maximum_shade_daily
    if(timeinterval == 365){
      leapyears <- seq(1900, 2060, 4)
      for(mm in 1:nyears){
        if(mm == 1){
          currenty <- ystart
        }else{
          currenty <- ystart + mm - 1
        }
        if(currenty %in% leapyears){
          dayoy <- seq(1, 366)
        }else{
          dayoy <- seq(1, 365)
        }
        if(mm == 1){
          day_of_year <- dayoy
        }else{
          day_of_year <- c(day_of_year, dayoy)
        }
      }
    }else{
      day_of_year <- rep(doys12, nyears)
    }

    if(timeinterval<365){
      microdaily <- 0 # run microclimate model as normal, where each day is iterated max_iterations_per_day times starting with the initial condition of uniform soil temp at mean monthly temperature
    }else{
      microdaily <- 1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in
    }

    # now check if doing something other than middle day of each month, and create appropriate vector of Day of Year
    ndays <- length(day_of_year)
    end_day <- ndays

    start_day <- 1 # start day
    yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)

    ################## location and terrain #################################

    if(is.numeric(loc)==FALSE){ # might not be quite right format, try correcting
      loc <- cbind(as.numeric(loc[1]), as.numeric(loc[2]))
    }
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    # Get reference longitude for solar noon correction
    solar_noon_longitude <- get_timezone_alref(lon = x[1], lat = x[2], timezone = timezone)
    hemisphere <- ifelse(x[2] < 0, 2, 1) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    latitude_degrees <- abs(trunc(x[2]))
    latitude_minutes <- (abs(x[2])-latitude_degrees)*60
    longitude_degrees <- abs(trunc(x[1]))
    longitude_minutes <- (abs(x[1])-longitude_degrees)*60
    azmuth <- aspect

    surface_emissivity <- rep(surface_emissivity,ndays)
    SoilMoist <- initial_soil_moisture

    if(soilgrids == 1){
      sg <- fetch_soilgrids(x, depths)
      if (!is.null(sg)) { air_entry_water_potential <- sg$air_entry_water_potential; saturated_hydraulic_conductivity <- sg$saturated_hydraulic_conductivity; campbell_b_parameter <- sg$campbell_b_parameter; soil_bulk_density <- sg$soil_bulk_density; bulk_density <- sg$bulk_density }
    }

    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("package 'terra' is needed. Please install it.",
           call. = FALSE)
    }
    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("package 'ncdf4' is needed. Please install it.",
           call. = FALSE)
    }
    library(ncdf4)
    require(utils)
    retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=100, sleep=30) {
      attempts = 0
      retval = try(eval(expr))
      while (isError(retval)) {
        attempts = attempts + 1
        if (attempts >= maxErrors) {
          msg = sprintf("retry: too many retries [[%s]]", capture.output(str(retval)))
          futile.logger::flog.fatal(msg)
          stop(msg)
        } else {
          msg = sprintf("retry: error in attempt %i/%i [[%s]]", attempts, maxErrors,
                        capture.output(str(retval)))
          futile.logger::flog.error(msg)
          warning(msg)
        }
        if (sleep > 0) Sys.sleep(sleep)
        retval = try(eval(expr))
      }
      return(retval)
    }
    if(is.na(terra_source)){
     var <- "tmax"
     baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1950_CurrentYear_GLOBE.nc#fillmismatch")

    nc <- retry(nc_open(baseurlagg))
    lon <- retry(ncvar_get(nc, "lon"))
    lat <- retry(ncvar_get(nc, "lat"))
    flat <- match(abs(lat - x[2]) < 1/48, 1)
    latindex <- which(flat %in% 1)
    if(length(latindex) == 0){
      flat <- match(abs(lat - x[2]) < 1/47.9, 1)
      latindex <- which(flat %in% 1)[1]
    }
    flon <- match(abs(lon - x[1]) < 1/48, 1)
    lonindex <- which(flon %in% 1)
    if(length(lonindex) == 0){
      flon <- match(abs(lon - x[1]) < 1/47.9, 1)
      lonindex <- which(flon %in% 1)[1]
    }
    time <- retry(ncvar_get(nc, "time"))
    alldays <- head(seq(as.Date('1900-01-01'), as.Date(paste0(yfinish + 1,'-01-01')), 'days'), -1)
    days1900 <- seq(1:length(alldays))
    allmonths <- head(seq(as.Date('1900-02-01'), as.Date(paste0(yfinish + 1,'-02-01')), '1 month'), -1) - 1
    allmonth.days <- head(days1900[which(alldays %in% allmonths & alldays >= as.Date('1957-12-01'))], -1)
    month.dates.to.do <- days1900[which(alldays %in% allmonths & alldays > as.Date(paste0(ystart - 1, '-12-01')) & alldays < as.Date(paste0(yfinish, '-12-31')))]
    timeindex <- which(time %in% month.dates.to.do)
    start <- c(lonindex, latindex, timeindex[1])
    count <- c(1, 1, length(timeindex))
    if(scenario == 0){
      message('extracting maximum air temperature data from TerraClimate\n')
      reference_temperature_max <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
      var <- 'tmin'
      message('extracting minimum air temperature data from TerraClimate \n')
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1950_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      reference_temperature_min <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
      message('extracting precipitation data from TerraClimate \n')
      var <- 'ppt'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1950_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      rainfall <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
      message('extracting wind speed data from TerraClimate \n')
      var <- 'ws'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1950_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      WIND <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
      message('extracting vapour pressure deficit data from TerraClimate \n')
      var <- 'vpd'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1950_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      VPD <- retry(as.numeric(ncvar_get(nc, varid = var,start = start, count)))
      nc_close(nc)
      message('extracting solar radiation data from TerraClimate \n')
      var <- 'srad'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1950_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      SRAD <- retry(as.numeric(ncvar_get(nc, varid = var,start = start, count)))
      if(soil_moisture_model == 0){
        # extract soil moisture
        #soilmoisture <- suppressWarnings(terra::rast(paste(folder, "/soilw.mon.ltm.v2.nc", sep = "")))
        message("extracting soil moisture data from TerraClimate")
        #SoilMoist <- as.numeric(terra::extract(soilmoisture, x)) / 1000 # this is originally in mm/m
        var <- 'soil'
        baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1950_CurrentYear_GLOBE.nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        SoilMoist <- retry(as.numeric(ncvar_get(nc, varid = var,start = start, count))) / 1000 * (1 - bulk_density / mineral_density) # this is originally in mm/m
        nc_close(nc)
      }
    }else{
      if(scenario == 2){
        base <- '/thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data_plus2C/TerraClimate_plus2C'
      }else{
        base <- '/thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data_plus4C/TerraClimate_plus4C'
      }
      start <- c(lonindex, latindex, 1)
      count <- c(1, 1, -1)
      for(i in 1:length(yearlist)){
        var <- "tmax"
        baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' maximum air temperature data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          reference_temperature_max <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          reference_temperature_max <- c(reference_temperature_max, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        var <- "tmin"
        baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' minimum air temperature data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          reference_temperature_min <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          reference_temperature_min <- c(reference_temperature_min, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        var <- "ppt"
        baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' precipitation data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          rainfall <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          rainfall <- c(rainfall, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        var <- "vpd"
        baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' vapour pressure deficit data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          VPD <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          VPD <- c(VPD, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        var <- "srad"
        baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
        nc <- retry(nc_open(baseurlagg))
        message(paste0('extracting plus ', scenario,' solar radiation data from TerraClimate for ', yearlist[i], '\n'))
        if(i == 1){
          SRAD <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
        }else{
          SRAD <- c(SRAD, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))))
        }
        nc_close(nc)
        if(soil_moisture_model == 0){
          var <- "soil"
          baseurlagg <- paste0(paste0("http:/", base, "_", var),"_", yearlist[i], ".nc#fillmismatch")
          nc <- retry(nc_open(baseurlagg))
          message(paste0('extracting plus ', scenario,' soil moisture data from TerraClimate for ', yearlist[i], '\n'))
          if(i == 1){
            SoilMoist <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))) / 1000 # this is originally in mm/m
          }else{
            SoilMoist <- c(SoilMoist, retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))) / 1000)# this is originally in mm/m
          }
          nc_close(nc)
        }
      }
      timeindex <- which(time %in% month.dates.to.do)
      start <- c(lonindex, latindex, timeindex[1])
      count <- c(1, 1, length(timeindex))
      message('extracting wind speed data from TerraClimate \n')
      var <- 'ws'
      baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1950_CurrentYear_GLOBE.nc#fillmismatch")
      nc <- retry(nc_open(baseurlagg))
      WIND <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      nc_close(nc)
    }
    }else{
      alldays <- head(seq(as.Date('1900-01-01'), as.Date(paste0(yfinish + 1,'-01-01')), 'days'), -1)
      days1900 <- seq(1:length(alldays))
      allmonths <- head(seq(as.Date('1900-02-01'), as.Date(paste0(yfinish + 1,'-02-01')), '1 month'), -1) - 1
      allmonth.days <- head(days1900[which(alldays %in% allmonths & alldays >= as.Date('1957-12-01'))], -1)
      month.dates.to.do <- days1900[which(alldays %in% allmonths & alldays > as.Date(paste0(ystart - 1, '-12-01')) & alldays < as.Date(paste0(yfinish, '-12-31')))]
      terra <- as.data.frame(get_terra(x = loc, ystart = ystart, yfinish = yfinish, scenario = 0, source = terra_source))
      reference_temperature_min <- terra$reference_temperature_min
      reference_temperature_max <- terra$reference_temperature_max
      rainfall <- terra$rainfall
      VPD <- terra$VPD
      SRAD <- terra$SRAD
      SoilMoist <- terra$SoilMoist
      WIND <- terra$WIND
      if(scenario == 4){
        terra_cc <- as.data.frame(get_terra(x = loc, ystart = ystart, yfinish = yfinish, scenario = 4, source = terra_source))
        reference_temperature_min <- terra_cc$reference_temperature_min
        reference_temperature_max <- terra_cc$reference_temperature_max
        rainfall <- terra_cc$rainfall
        VPD <- terra_cc$VPD
        SRAD <- terra_cc$SRAD
        SoilMoist <- terra_cc$SoilMoist
      }
      if(scenario == 2){
        terra_cc <- as.data.frame(get_terra(x = loc, ystart = ystart, yfinish = yfinish, scenario = 2, source = terra_source))
        reference_temperature_min <- terra_cc$reference_temperature_min
        reference_temperature_max <- terra_cc$reference_temperature_max
        rainfall <- terra_cc$rainfall
        VPD <- terra_cc$VPD
        SRAD <- terra_cc$SRAD
        SoilMoist <- terra_cc$SoilMoist
      }
    }
    # load global climate files
    gcfolder<-paste(.libPaths()[1],"/gcfolder.rda",sep="")
    if(file.exists(gcfolder)==FALSE){
      folder<-"c:/globalclimate"
      if(file.exists(paste0(folder,"/global_climate.nc"))==FALSE){
        message("You don't appear to have the global climate data set - \n run function get.global.climate(folder = 'folder you want to put it in') .....\n exiting function micro_global")
        opt <- options(show.error.messages=FALSE)
        on.exit(options(opt))
        stop()
      }
    }else{
      load(gcfolder)
    }


    global_climate <- terra::rast(paste0(folder, "/global_climate.nc"))
    CLIMATE <- t(as.numeric(terra::extract(global_climate, x)))
    #elevation <- as.numeric(CLIMATE[, 1])
    if(class(elevation)[1] == "SpatRaster"){
      cat('using elevation DEM provided to function call \n')
    }else{
      if (!requireNamespace("geodata", quietly = TRUE)) {
        stop("package 'geodata' is needed. Please install it.",
             call. = FALSE)
      }
      library(geodata)
      elevation <- geodata::worldclim_global(var = 'elevation', res = 2.5, path=tempdir())
    }
    elevation <- as.numeric(terra::extract(elevation, x))
    if(is.na(elevation)){
      elevation <- 0
    }
    if(terrain == 1){
      elevatr <- 1
    }
    if(is.na(elevation) & elevatr == 1){
      dem_result <- fetch_dem(loc[1], loc[2], terrain = terrain, horizon.step = 10)
      dem <- dem_result$dem
      elevation <- dem_result$elevation
      if (terrain == 1) { slope <- dem_result$slope; aspect <- dem_result$aspect; horizon_angles <- dem_result$horizon_angles }
    }
    horizon_angles<-as.matrix(horizon_angles) #horizon angles
    VIEWF <- 1-sum(sin(as.data.frame(horizon_angles) * pi / 180)) / length(horizon_angles) # convert horizon angles to radians and calc view factor(s)

    require(RNetCDF)
    ALTITUDES <- elevation
    dbrow <- 1
    if(is.na(max(SoilMoist, elevation, reference_temperature_max)) == TRUE){
      message("Sorry, there is no environmental data for this location")
      SoilMoist <- as.numeric(terra::extract(soilmoisture, cbind(140, -35))) / 1000 # this is originally in mm/m
    }
    delta_elev <- 0
    if(is.na(elevation) == FALSE){ # check if user-specified elevation
      delta_elev <- elevation - elevation # get delta for lapse rate correction
      elevation <- elevation # now make final elevation the user-specified one
    }
    adiab_corr_max <- delta_elev * maximum_temperature_lapse_rate
    adiab_corr_min <- delta_elev * minimum_temperature_lapse_rate

    if(is.na(rainfall[1])){
      cat("no climate data for this site, using dummy data so solar is still produced \n")
      CLIMATE <- as.numeric(terra::extract(global_climate, cbind(140, -35)))
      CLIMATE[2:97] <- 0
      elevation<-as.numeric(CLIMATE[, 1])
      delta_elev <- 0
      if(is.na(elevation) == FALSE){ # check if user-specified elevation
        delta_elev <- elevation - elevation # get delta for lapse rate correction
        elevation <- elevation # now make final elevation the user-specified one
      }
      adiab_corr_max <- delta_elev * maximum_temperature_lapse_rate
      adiab_corr_min <- delta_elev * minimum_temperature_lapse_rate
      rainfall <- CLIMATE[, 2:13] * 0
      #stop()
    }
    RAINYDAYS <- CLIMATE[, 14:25] / 10
    # reference_wind_max <- CLIMATE[, 26:37] / 10 * wind_multiplier
    # reference_wind_min<-reference_wind_max * 0.1 # impose diurnal cycle
    # reference_temperature_min <- CLIMATE[, 38:49] / 10
    # reference_temperature_max <- CLIMATE[, 50:61] / 10
    reference_wind_max <- WIND
    reference_wind_min<-reference_wind_max * 0.1 # impose diurnal cycle
    reference_temperature_max <- reference_temperature_max + adiab_corr_max
    reference_temperature_min <- reference_temperature_min + adiab_corr_min
    ALLMINTEMPS <- reference_temperature_min
    ALLMAXTEMPS <- reference_temperature_max
    ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)
    MEANTEMPS <- (reference_temperature_max + reference_temperature_min) / 2
    e <- WETAIR(db = MEANTEMPS, rh = 100)$e - VPD * 1000
    es <- WETAIR(db = reference_temperature_max, rh = 100)$esat
    reference_humidity_min <- (e / es) * 100
    es <- WETAIR(db = reference_temperature_min, rh = 100)$esat
    reference_humidity_max <- (e / es) * 100
    reference_humidity_min[reference_humidity_min>100]<-100
    reference_humidity_min[reference_humidity_min<0]<-0.01
    reference_humidity_max[reference_humidity_max>100]<-100
    reference_humidity_max[reference_humidity_max<0]<-0.01

    cat("running micro_global to get clear sky solar \n")
    if(global_aerosol_database == 0){
      aerosol_optical_depth <- c(0.0670358341290886, 0.0662612704779235, 0.065497075238002, 0.0647431301168489, 0.0639993178022531, 0.0632655219571553, 0.0625416272145492, 0.0611230843885423, 0.0597427855962549, 0.0583998423063099, 0.0570933810229656, 0.0558225431259535, 0.0545864847111214, 0.0533843764318805, 0.0522154033414562, 0.0499736739981675, 0.047855059159556, 0.0458535417401334, 0.0439633201842001, 0.0421788036108921, 0.0404946070106968, 0.0389055464934382, 0.0374066345877315, 0.0359930755919066, 0.0346602609764008, 0.0334037648376212, 0.0322193394032758, 0.0311029105891739, 0.0300505736074963, 0.0290585886265337, 0.0281233764818952, 0.0272415144391857, 0.0264097320081524, 0.0256249068083005, 0.0248840604859789, 0.0241843546829336, 0.0235230870563317, 0.0228976873502544, 0.0223057135186581, 0.0217448478998064, 0.0212128934421699, 0.0207077699817964, 0.0202275105711489, 0.0197702578594144, 0.0193342605242809, 0.0189178697551836, 0.0177713140039894, 0.0174187914242432, 0.0170790495503944, 0.0167509836728154, 0.0164335684174899, 0.0161258546410128, 0.0158269663770596, 0.0155360978343254, 0.0152525104459325, 0.0149755299703076, 0.0147045436435285, 0.0144389973831391, 0.0141783930434343, 0.0134220329447663, 0.0131772403830191, 0.0129356456025128, 0.0126970313213065, 0.0124612184223418, 0.0122280636204822, 0.01199745718102, 0.0115436048739351, 0.0110993711778668, 0.0108808815754663, 0.0106648652077878, 0.0104513876347606, 0.0102405315676965, 0.00982708969547694, 0.00962473896278535, 0.00903679230300494, 0.00884767454432418, 0.0083031278398166, 0.00796072474935954, 0.00755817587626185, 0.00718610751850881, 0.00704629977586921, 0.00684663903049612, 0.00654155580333479, 0.00642947339729728, 0.00627223096874308, 0.00603955966866779, 0.00580920937536261, 0.00568506186880564, 0.00563167068287251, 0.00556222005081865, 0.00550522989971023, 0.00547395763028062, 0.0054478983436216, 0.00541823364504573, 0.00539532163908382, 0.00539239864119488, 0.00541690124712384, 0.00551525885358836, 0.00564825853509463, 0.00577220185074264, 0.00584222986640171, 0.00581645238345584, 0.00566088137411449, 0.00535516862329704, 0.00489914757707667, 0.00432017939770409, 0.0036813032251836, 0.00309019064543606, 0.00270890436501562, 0.00276446109239711, 0.00356019862584603)
    }else{
      aerosol_optical_depth <- 0
    }
    micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, aerosol_optical_depth = aerosol_optical_depth, timeinterval = 365, solar_model_only = 1)
    clearskyrad <- micro_clearsky$micromet_lowshade[, c(1, 13)][, 2]
    clearskymean <- leapfix(clearskyrad, yearlist, 24)

    dates2 <- head(seq(as.Date(paste0(ystart, '-01-01')), as.Date(paste0(yfinish + 1, '-01-01')), by = "days"), -1)
    dates <- head(seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'Etc/GMT+10'), as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y ", tz = 'Etc/GMT+10'), by = 'hours'), -1)
    #ndays <- length(dates2)

    clearskymean <- aggregate(clearskymean, by = list(format(dates, '%Y-%m')), FUN = mean)[, 2]

    #allclearsky <- allclearsky[1:ndays]
    # convert from W/d to MJ/d
    #allclearsky <- allclearsky * 3600 / 1e6
    cloud <- (1 - SRAD / clearskymean) * 100
    cloud[cloud < 0] <- 0
    cloud[cloud > 100] <- 100
    cloud_min <- cloud * 0.5
    cloud_max <- cloud * 2
    cloud_min[cloud_min > 100] <- 100
    cloud_max[cloud_max > 100] <- 100
    cloud_min[cloud_min < 0] <- 0
    cloud_max[cloud_max < 0] <- 0
    #cloud_min <- CLIMATE[, 86:97] / 10
    if(clearsky == 1){
      cloud_min <- cloud_min * 0
      cloud_max <- cloud_min
    }

    # correct for fact that wind is measured at 10 m height
    # wind shear equation v / vo = (h / ho)^a
    #where
    #v = the velocity at height h (m/s)
    #vo = the velocity at height ho (m/s)
    #a = the wind shear exponent
    #Terrain   Wind Shear Exponent
    #- a -
    #  Open water   0.1
    #Smooth, level, grass-covered   0.15 (or more commonly 1/7)
    #Row crops 	0.2
    #Low bushes with a few trees 	0.2
    #Heavy trees 	0.25
    #Several buildings 	0.25
    #Hilly, mountainous terrain 	0.25
    # source http://www.engineeringtoolbox.com/wind-shear-d_1215.html
    reference_wind_min <- reference_wind_min * (2 / 10) ^ 0.15
    reference_wind_max <- reference_wind_max * (2 / 10) ^ 0.15
    # impose uniform warming
    reference_temperature_max <- reference_temperature_max + air_temperature_offset
    reference_temperature_min <- reference_temperature_min + air_temperature_offset
    methspline <- 'periodic'
    if(timeinterval != 12){ # spline from 12 days to chosen time interval
      for(i in 1:nyears){
        if(yearlist[i] %in% seq(1900, 2060, 4)){
          xmax <- 366
        }else{
          xmax <- 365
        }
        if(i == 1){
          start <- 1
          end <- 12
          TMAXX1 <-suppressWarnings(spline(doys12,reference_temperature_max[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          TMAXX2<-TMAXX1$y
          TMINN1 <-suppressWarnings(spline(doys12,reference_temperature_min[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          TMINN2 <- TMINN1$y
          RHMAXX1 <-suppressWarnings(spline(doys12,reference_humidity_max[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          RHMAXX2 <- RHMAXX1$y
          RHMINN1 <-suppressWarnings(spline(doys12,reference_humidity_min[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          RHMINN2 <- RHMINN1$y
          CCMAXX1 <-suppressWarnings(spline(doys12,cloud_max[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          CCMAXX2 <- CCMAXX1$y
          CCMINN1 <-suppressWarnings(spline(doys12,cloud_min[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          CCMINN2 <- CCMINN1$y
          WNMAXX1 <-suppressWarnings(spline(doys12,reference_wind_max[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          WNMAXX2 <- WNMAXX1$y
          WNMINN1 <-suppressWarnings(spline(doys12,reference_wind_min[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          WNMINN2 <- WNMINN1$y
          if(soil_moisture_model==0){
            SoilMoist1 <- suppressWarnings(spline(doys12,SoilMoist[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
            SoilMoist2 <- SoilMoist1$y
          }
        }else{
          start <- end + 1
          end <- end + 12
          TMAXX1 <-suppressWarnings(spline(c(0, doys12), c(tail(TMAXX2, 1), reference_temperature_max[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          TMAXX2<-c(TMAXX2, TMAXX1$y)
          TMINN1 <-suppressWarnings(spline(c(0, doys12), c(tail(TMINN2, 1), reference_temperature_min[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          TMINN2<-c(TMINN2, TMINN1$y)
          RHMAXX1 <-suppressWarnings(spline(c(0, doys12), c(tail(RHMAXX2, 1), reference_humidity_max[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          RHMAXX2<-c(RHMAXX2, RHMAXX1$y)
          RHMINN1 <-suppressWarnings(spline(c(0, doys12), c(tail(RHMINN2, 1), reference_humidity_min[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          RHMINN2<-c(RHMINN2, RHMINN1$y)
          CCMAXX1 <-suppressWarnings(spline(c(0, doys12), c(tail(CCMAXX2, 1), cloud_max[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          CCMAXX2<-c(CCMAXX2, CCMAXX1$y)
          CCMINN1 <-suppressWarnings(spline(c(0, doys12), c(tail(CCMINN2, 1), cloud_min[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          CCMINN2<-c(CCMINN2, CCMINN1$y)
          WNMAXX1 <-suppressWarnings(spline(c(0, doys12), c(tail(WNMAXX2, 1), reference_wind_max[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          WNMAXX2<-c(WNMAXX2, WNMAXX1$y)
          WNMINN1 <-suppressWarnings(spline(c(0, doys12), c(tail(WNMINN2, 1), reference_wind_min[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          WNMINN2<-c(WNMINN2, WNMINN1$y)
          if(soil_moisture_model==0){
            SoilMoist1 <-suppressWarnings(spline(c(0, doys12), c(tail(SoilMoist2, 1), SoilMoist[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
            SoilMoist2<-c(SoilMoist2, SoilMoist1$y)
          }
        }
      }
      reference_temperature_max <- TMAXX2
      reference_temperature_min <- TMINN2
      reference_humidity_max <- RHMAXX2
      reference_humidity_min <- RHMINN2
      reference_wind_max <- WNMAXX2
      reference_wind_min <- WNMINN2
      cloud_max <- CCMAXX2
      cloud_min <- CCMINN2
      if(soil_moisture_model==0){
        SoilMoist <- SoilMoist2
      }
    }
    cloud_min[cloud_min > 100] <- 100
    cloud_max[cloud_max > 100] <- 100
    cloud_min[cloud_min < 0] <- 0
    cloud_max[cloud_max < 0] <- 0

    orig.rainfall <- rainfall

    # get annual mean temp for creating deep soil (2m) boundary condition
    avetemp<-(sum(reference_temperature_max)+sum(reference_temperature_min))/(length(reference_temperature_max)*2)
    if(is.na(initial_soil_temperature[1])){
      soilinit <- rep(avetemp, 20)
      spinup <- 1
    }else{
      if(snow_model == 0){
        soilinit <- c(initial_soil_temperature, rep(avetemp, 10))
      }else{
        soilinit <- c(rep(avetemp, 8), initial_soil_temperature[1:10], rep(avetemp, 2))
      }
      spinup <- 0
    }

    mean_annual_temperature<-mean(unlist(ALLTEMPS))

    if(timeinterval == 12){
      if(is.na(deep_soil_temperature)){
        if(nyears == 1){
          avetemp <- (sum(reference_temperature_max) + sum(reference_temperature_min)) / (length(reference_temperature_max)*2)
          deep_soil_temperature <-rep(avetemp, length(reference_temperature_max))
        }else{
          avetemp <- rowMeans(cbind(reference_temperature_max, reference_temperature_min), na.rm=TRUE)
          deep_soil_temperature <- terra::roll(avetemp, n = 12, fun = mean, type = 'to')
          yearone <- rep((sum(reference_temperature_max[1:12]) + sum(reference_temperature_min[1:12])) / (12 * 2), 12)
          deep_soil_temperature[1:12] <- yearone
        }
      }else{
        mean_annual_temperature <- mean(deep_soil_temperature)
      }
    }else{
      if(is.na(deep_soil_temperature)){
        if(nyears == 1){
          avetemp <- (sum(reference_temperature_max) + sum(reference_temperature_min)) / (length(reference_temperature_max)*2)
          deep_soil_temperature <- rep(avetemp, ndays)
        }else{
          avetemp <- rowMeans(cbind(reference_temperature_max, reference_temperature_min), na.rm=TRUE)
          if(length(reference_temperature_max) < 365){
            deep_soil_temperature <- rep((sum(reference_temperature_max) + sum(reference_temperature_min)) / (length(reference_temperature_max) * 2), length(reference_temperature_max))
          }else{
            deep_soil_temperature <- terra::roll(avetemp, n = 365, fun = mean, type = 'to')
            yearone <- rep((sum(reference_temperature_max[1:365]) + sum(reference_temperature_min[1:365])) / (365 * 2), 365)
            deep_soil_temperature[1:365] <- yearone
          }
        }
      }else{
        mean_annual_temperature <- mean(deep_soil_temperature)
      }
    }
    deep_soil_temperature <- deep_soil_temperature  # deep_soil_temperature already computed above; alias kept for clarity
    daymon<-c(31,28,31,30,31,30,31,31,30,31,30,31) # days in each month
    bb <- 0
    # if doing daily sims, spread rainfall evenly across days based on mean monthly rainfall and the number of rainy days per month
    if(timeinterval==365){
      for(j in 1:length(yearlist)){
        if(yearlist[j] %in% seq(1900, 2060, 4)){
          dys <- 366
          daymon2 <- c(31,29,31,30,31,30,31,31,30,31,30,31) # days in each month
        }else{
          dys <- 365
          daymon2 <- daymon
        }
        RAINFALL1<-1:dys
        sort<-matrix(data = 0,nrow = dys,ncol = 2)
        m<-1
        b<-0
        for (i in 1:12){ #begin loop through 12 months of year
          ndays2=daymon2[i]
          for (k in 1:ndays2){
            b<-b+1
            sort[m,1]<-i
            sort[m,2]<-b
            if(k<=RAINYDAYS[i] & rainfrac>0){
              if(k==1){
                RAINFALL1[m]<-rainfall[i + (j - 1) * 12]*rainfrac*rain_multiplier # if first day of month, make user-specified fraction of monthly rainfall fall on first day
              }else{
                RAINFALL1[m]<-(rainfall[i + (j - 1) * 12]*(1-rainfrac)*rain_multiplier)/RAINYDAYS[i] # make remaining rain fall evenly over the remaining number of rainy days for the month, starting at the beginning of the month
              }
            }else{
              if(rainfrac==0){
                RAINFALL1[m]<-(rainfall[i + (j - 1) * 12]*rain_multiplier)/RAINYDAYS[i]
              }else{
                RAINFALL1[m]<-0
              }
            }
            m<-m+1
            if(b>RAINYDAYS[i]){
              b<-0
            }
          }
        }
        RAINFALL2<-as.data.frame(cbind(RAINFALL1,sort))
        #RAINFALL2<-RAINFALL2[order(RAINFALL2$V2,RAINFALL2$V3),] # this line scatters the rainy days evenly across each month - snow predictions better if it is commented out so get rainy days all in a row within the month
        if(j == 1){
          RAINFALL3<-RAINFALL2$RAINFALL1
          if(reference_temperature_min[1]<snowtemp){
            RAINFALL3[1]<-0 # this is needed in some cases to allow the integrator to get started
          }
        }else{
          RAINFALL3 <- c(RAINFALL3, RAINFALL2$RAINFALL1)
        }
      }
    }else{
      if(timeinterval!=12){
        RAINFALL3<-rep(rep(sum(rainfall)/timeinterval,timeinterval),nyears) # just spread evenly across every day
      }else{ # running middle day of each month - divide monthly rain by number of days in month
        RAINFALL3<-rainfall/rep(daymon,nyears)
      }
    }#end check doing daily sims
    rainfall <- RAINFALL3
    rainfall[!is.finite(rainfall)] <- 0

    global_radiation<-rep(0,24*ndays)

    hourly <- 0
    if(microclima > 0){
      cat('using microclima and elevatr to adjust solar for topographic and vegetation effects \n')
      if (!require("microclima", quietly = TRUE)) {
        stop("package 'microclima' is needed. Please install it.",
             call. = FALSE)
      }
      if (!require("zoo", quietly = TRUE)) {
        stop("package 'zoo' is needed. Please install it.",
             call. = FALSE)
      }
      cat("Downloading digital elevation data \n")
      lat <- x[2]
      long <- x[1]
      tt <- seq(as.POSIXct(paste0('01/01/',ystart), format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(paste0('31/12/',yfinish), format = "%d/%m/%Y", tz = 'UTC')+23*3600, by = 'hours')
      timediff <- x[1] / 15
      hour.microclima <- as.numeric(format(tt, "%H")) + timediff-floor(timediff)
      jd <- julday(as.numeric(format(tt, "%Y")), as.numeric(format(tt, "%m")), as.numeric(format(tt, "%d")))
      dem_result <- fetch_dem(loc[1], loc[2],
                              dem_resolution = microclima.dem_resolution, zmin = microclima.zmin,
                              pixels = NULL,
                              existing_dem = if (inherits(dem, c("SpatRaster", "RasterLayer"))) dem else NULL,
                              terrain = 1, horizon.step = 10, nangles = 36,
                              slope = slope, aspect = aspect)
      dem <- dem_result$dem; slope <- dem_result$slope; aspect <- dem_result$aspect
      ha36 <- dem_result$horizon_angles
      for (i in 1:length(hour.microclima)) {
        saz <- solazi(hour.microclima[i], lat, long, jd[i], merid = long)
        saz <- round(saz/10, 0) + 1
        saz <- ifelse(saz > 36, 1, saz)
        ha[i] <- ha36[saz]
      }

      methspline <- 'periodic'
      for(i in 1:nyears){
        if(yearlist[i] %in% seq(1900, 2060, 4)){
          xmax <- 366
        }else{
          xmax <- 365
        }
        if(i == 1){
          start <- 1
          end <- 12
          cloud1 <-suppressWarnings(spline(doys12,cloud[start:end],n=xmax,xmin=1,xmax=xmax,method=methspline))
          cloud2<-cloud1$y
        }else{
          start <- end + 1
          end <- end + 12
          cloud1 <-suppressWarnings(spline(c(0, doys12), c(tail(cloud2, 1), cloud[start:end]), n = xmax, xmin = 1, xmax = xmax, method = methspline))
          cloud2<-c(cloud2, cloud1$y)
        }
      }
      cloudhr <- cbind(rep(seq(1, length(cloud2)),24), rep(cloud2, 24))
      cloudhr <- cloudhr[order(cloudhr[,1]),]
      cloudhr <- cloudhr[,2]
      #cloudhr <- leapfix(cloudhr, yearlist, 24)
      dsw2 <- leapfix(clearskyrad, yearlist, 24) *(0.36+0.64*(1-cloudhr/100)) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
      # partition total solar into diffuse and direct using code from microclima::hourlyNCEP
      sol <- compute_solar_partition(dsw2, jd, hour.microclima, lat, long,
                                     slope, aspect, ha, microclima.LOR, microclima.leaf_area_index)
      SOLRhr_all <- sol$SOLRhr_all
      diffuse_frac_all <- sol$diffuse_frac_all
      if(microclima == 2){ # use hourly solar from microclima
        hourly <- 2
        VIEWF <- 1 # accounted for already in microclima cals
        horizon_angles <- rep(0, 24) # accounted for already in microclima calcs
      }
    }else{
      diffuse_frac_all <- NA
    }
    if(timeinterval == 12 & microclima > 0){
      dates_all <- head(seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y ", tz = 'UTC'), by = 'hours'), -1)
      month.dates.to.do2 <- days1900[which(alldays %in% allmonths & alldays > as.Date(paste0(ystart, '-01-01')) & alldays < as.Date(paste0(yfinish + 1, '-01-01')))]
      mon <- alldays[month.dates.to.do2]
      dates3 <- dates_all[which(format(dates_all, "%Y-%m-%d") %in% as.character(mon))]
      global_radiation <- SOLRhr_all[which(dates_all %in% dates3)]
      diffuse_frac <- diffuse_frac_all[which(dates_all %in% dates3)]
    }else{
      diffuse_frac <- diffuse_frac_all
    }

    ndays<-length(rainfall)
    if(length(aerosol_optical_depth) < 111){ # no user supplied values, compute with GADS
      aerosol_optical_depth <- compute_tai(longlat, global_aerosol_database, TAI_ELTERMAN)
    }
    ################ soil properties  ##################################################
    soil_nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
    if(soilgrids == 1){
      Numtyps <- 10 # number of substrate types
      soil_nodes[1:10,] <- c(1:10) # deepest nodes for each substrate type
    }else{
      Numtyps <- 2 # number of soil types
      soil_nodes[1,1:ndays]<-3 # deepest node for first substrate type
      soil_nodes[2,1:ndays]<-9 # deepest node for second substrate type
    }
    albedo<-rep(surface_reflectivity,ndays) # soil reflectances
    soil_wetness<-rep(soil_wetness,ndays) # soil wetness
    if(soil_moisture_model==0){
      moists2<-matrix(nrow= 10, ncol = ndays, data=0) # set up an empty vector for soil moisture values through time
      moists2[1,]<-SoilMoist # fill the first row with monthly soil moisture values
      moists2[2,]<-moists2[1,] # make this row same as first row
      soil_moisture_profile<-moists2
    }else{
      moists2<-matrix(nrow=10, ncol = ndays, data=0) # set up an empty vector for soil moisture values through time
      moists2[1:10,]<-initial_soil_moisture
      moists2[moists2>(1-bulk_density/mineral_density)]<-(1-bulk_density/mineral_density)
      soil_moisture_profile<-moists2
    }

    # now make the soil properties matrix
    # columns are:
    #1) bulk density (Mg/m3)
    #2) volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    #3) thermal conductivity (W/mK)
    #4) specific heat capacity (J/kg-K)
    #5) mineral density (Mg/m3)
    soil_properties<-matrix(data = 0, nrow = 10, ncol = 5) # create an empty soil properties matrix
    if(soilgrids == 1){
      soil_properties[,1]<-bulk_density
      soil_properties[,2] <- 1 - bulk_density / mineral_density # not used if soil moisture computed
      soil_properties[soil_properties[,2] < 0.26, 2] <- 0.26
      soil_properties[,3]<-mineral_conductivity
      soil_properties[,4]<-mineral_heat_capacity
      soil_properties[,5]<-mineral_density
      if(organic_soil_cap==1){
        soil_properties[1:2,3]<-0.2
        soil_properties[1:2,4]<-1920
      }
      if(organic_soil_cap==2){
        soil_properties[1:2,3]<-0.1
        soil_properties[3:4,3]<-0.25
        soil_properties[1:4,4]<-1920
        soil_properties[1:4,5]<-1.3
        soil_properties[1:4,1]<-0.7
      }
    }else{
      soil_properties[1,1]<-bulk_density # insert soil bulk density to profile 1
      soil_properties[2,1]<-bulk_density # insert soil bulk density to profile 2
      soil_properties[1,2]<-min(0.26, 1 - bulk_density / mineral_density) # insert saturated water content to profile 1
      soil_properties[2,2]<-min(0.26, 1 - bulk_density / mineral_density) # insert saturated water content to profile 2
      if(organic_soil_cap==1){ # insert thermal conductivity to profile 1, and see if 'organic cap' added on top
        soil_properties[1,3]<-0.2 # mineral thermal conductivity
      }else{
        soil_properties[1,3]<-mineral_conductivity # mineral thermal conductivity
      }
      soil_properties[2,3]<-mineral_conductivity # insert thermal conductivity to profile 2
      if(organic_soil_cap==1){ # insert specific heat to profile 1, and see if 'organic cap' added on top
        soil_properties[1,4]<-1920 # mineral heat capacity
      }else{
        soil_properties[1,4]<-mineral_heat_capacity
      }
      soil_properties[2,4]<-mineral_heat_capacity # insert specific heat to profile 2
      soil_properties[1,5]<-mineral_density # insert mineral density to profile 1
      soil_properties[2,5]<-mineral_density # insert mineral density to profile 2
    }
    #########################################################################################

    # Next four parameters are segmented velocity profiles due to bushes, rocks etc. on the surface
    #IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO! (then roughness height is based on the parameter roughness_height)
    roughness_height_1 <- 0 # Top (1st) segment roughness height(m)
    roughness_height_2 <- 0 # 2nd segment roughness height(m)
    wind_profile_height_1 <- 0 # Top of (1st) segment, height above surface(m)
    wind_profile_height_2 <- 0 # 2nd segment, height above surface(m)

    # hourly option set to 0, so make empty vectors
    #hourly <- 0
    rainhourly <- 0
    reference_temperature <- rep(0,24*ndays)
    reference_humidity <- rep(0,24*ndays)
    reference_wind_speed <- rep(0,24*ndays)
    cloud_cover <- rep(0,24*ndays)
    rainfall_hourly <- rep(0,24*ndays)
    zenith_angle_hourly <- rep(-1,24*ndays)
    longwave_radiation <- rep(-1,24*ndays)
    # --- 3. Prepare model inputs ---
    micro_input <- build_microinput(ndays, roughness_height, tolerance, local_height, reference_height, Numtyps,
      roughness_height_1, roughness_height_2, wind_profile_height_1, wind_profile_height_2, start_day, end_day, hemisphere, latitude_degrees, latitude_minutes, longitude_degrees, longitude_minutes,
      solar_noon_longitude, slope, azmuth, elevation, precipitable_water, microdaily, mean_annual_temperature, orbital_eccentricity, VIEWF,
      snowtemp, snowdens, snowmelt, undercatch, rain_multiplier, runshade, soil_moisture_model,
      maximum_pooling_depth, evenrain, snow_model, rainmelt, output_to_csv, densfun, hourly,
      rainhourly, radiation_per_wavelength, scattered_uv, root_resistance, stomatal_closure_potential, leaf_resistance, stomatal_stability_parameter, root_radius, moist_error, moist_count, longwave_radiation_model, message,
      fail, snowcond, intercept, grasshade, solar_model_only, canopy_roughness_height, zero_plane_displacement, maxima_times, minima_times,
      spinup, maxsurf, max_iterations_per_day)

    if(length(leaf_area_index)<ndays){
      leaf_area_index <- rep(leaf_area_index[1],ndays)
    }
    if(shore==0){
      tides <- matrix(data = 0, nrow = 24 * ndays, ncol = 3) # make an empty matrix
    }
    micro <- build_micro_list(micro_input, day_of_year, surface_emissivity, depths, soil_nodes,
      maximum_shade_daily, minimum_shade_daily, reference_temperature_max, reference_temperature_min, reference_humidity_max, reference_humidity_min, cloud_max, cloud_min,
      reference_wind_max, reference_wind_min, reference_temperature, reference_humidity, reference_wind_speed, cloud_cover, global_radiation, rainfall_hourly, zenith_angle_hourly,
      longwave_radiation, albedo, soil_wetness, soilinit, horizon_angles, aerosol_optical_depth, soil_properties, soil_moisture_profile, rainfall,
      deep_soil_temperature, air_entry_water_potential, saturated_hydraulic_conductivity, campbell_b_parameter, soil_bulk_density, soil_mineral_density, root_density, leaf_area_index, tides)

    # write all input to csv files in their own folder
    if(write_input == 1){
      write_micro_csv(micro_input, day_of_year, surface_emissivity, depths, soil_nodes, maximum_shade_daily, minimum_shade_daily,
        maxima_times, minima_times, reference_temperature_max, reference_temperature_min, reference_humidity_max, reference_humidity_min, cloud_max, cloud_min,
        reference_wind_max, reference_wind_min, albedo, soil_wetness, soilinit, horizon_angles, aerosol_optical_depth, soil_properties, soil_moisture_profile,
        rainfall, deep_soil_temperature, air_entry_water_potential, soil_bulk_density, soil_mineral_density, campbell_b_parameter, saturated_hydraulic_conductivity, root_density, leaf_area_index, tides,
        reference_temperature, reference_humidity, reference_wind_speed, cloud_cover, global_radiation, rainfall_hourly, zenith_angle_hourly, longwave_radiation)
    }

    if(is.numeric(loc[1])){
      location <- paste("long",loc[1],"lat",loc[2])
    }else{
      location <- loc
    }
    if(timeinterval == 12){
      message(paste('running microclimate model for middle day of each month for',ystart,'to ', yfinish, 'at site',location,'\n'))
    }else{
      message(paste('running microclimate model, splining monthly data to daily, for',ystart,'to ', yfinish, 'at site',location,'\n'))
    }
    if(runmicro){
      message('Note: the output column `SOLR` in micromet_lowshade and micromet_highshade is for unshaded horizontal plane solar radiation \n')
    ptm <- proc.time() # Start timing
    microut <- microclimate(micro)
    message(paste0('runtime ', (proc.time() - ptm)[3], ' seconds')) # Stop the clock

    # --- 5. Process and return results ---
    out <- process_micro_output(microut, soil_moisture_model, snow_model, radiation_per_wavelength)

    # TerraClimate: restore original monthly rainfall for return list
    if(timeinterval == 12){
      rainfall <- orig.rainfall
    }
    if(max(out$micromet_lowshade[,1] == 0)){
      message("ERROR: the model crashed - try a different error tolerance (tolerance) or a different spacing in depths")
    }

    # Build date sequences for output
    days <- rep(seq(1, ndays), 24)
    days <- days[order(days)]
    dates <- days + rep(seq(0, 1380, 60), ndays) / 60 / 24 - 1  # numeric day-fraction for hourly output
    dates2 <- month.dates.to.do                                   # daily output dates
    tzone <- paste("Etc/GMT+", 0, sep = "")
    dates3 <- head(seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'UTC'),
                       as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y ", tz = 'UTC'),
                       by = 'hours'), -1)
    if(timeinterval == 12){
      dates_all <- head(seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'UTC'),
                            as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y ", tz = 'UTC'),
                            by = 'hours'), -1)
      month.dates.to.do2 <- days1900[which(alldays %in% allmonths &
                                           alldays > as.Date(paste0(ystart, '-01-01')) &
                                           alldays < as.Date(paste0(yfinish + 1, '-01-01')))]
      mon <- alldays[month.dates.to.do2]
      dates3 <- dates_all[which(format(dates_all, "%Y-%m-%d") %in% as.character(mon))]
    }else{
      dates <- dates3  # for daily data, dates is the full POSIXct hourly sequence
    }

    return(build_micro_return(out, rainfall, ndays, elevation, surface_reflectivity,
      longlat = c(x[1], x[2]), nyears, timeinterval, minimum_shade_daily, maximum_shade_daily,
      depths, dates, dates2, air_entry_water_potential, soil_bulk_density, soil_mineral_density, campbell_b_parameter, saturated_hydraulic_conductivity, dem, diffuse_frac,
      snow_model, radiation_per_wavelength,
      extra = list(dates3 = dates3, slope = slope, aspect = aspect, horizon_angles = horizon_angles)))
  }else{
    # Weather-only return (runmicro = FALSE): climate inputs without model run
    return(list(rainfall = rainfall, reference_temperature_max = reference_temperature_max, reference_temperature_min = reference_temperature_min,
                reference_humidity_max = reference_humidity_max, reference_humidity_min = reference_humidity_min, reference_wind_max = reference_wind_max,
                reference_wind_min = reference_wind_min, cloud_max = cloud_max, cloud_min = cloud_min,
                cloud_cover = cloud_cover, reference_wind_speed = reference_wind_speed, reference_temperature = reference_temperature, reference_humidity = reference_humidity,
                rainfall_hourly = rainfall_hourly, global_radiation = global_radiation, zenith_angle_hourly = zenith_angle_hourly, longwave_radiation = longwave_radiation,
                dates = dates, dates2 = dates2, dates3 = dates3,
                air_entry_water_potential = air_entry_water_potential, soil_bulk_density = soil_bulk_density, soil_mineral_density = soil_mineral_density, campbell_b_parameter = campbell_b_parameter, saturated_hydraulic_conductivity = saturated_hydraulic_conductivity))
  }
  } # end error trapping
} # end of micro_terra function
