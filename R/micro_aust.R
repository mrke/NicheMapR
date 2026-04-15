#' Australian implementation of the microclimate model.
#'
#' An implementation of the NicheMapR microclimate model that uses the AWAP daily weather database
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param ystart First year to run
#' @param yfinish Last year to run
#' @param surface_reflectivity Soil solar reflectance, decimal \%
#' @param elevation Elevation, if to be user specified (m)
#' @param slope Slope in degrees
#' @param aspect Aspect in degrees (0 = north)
#' @param depths Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param minimum_shade Minimum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param maximum_shade Maximum shade level to us (\%) (can be a single value or a vector of daily values)
#' @param local_height Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
#' @param ... Additional arguments, see Details
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
#' @usage micro_aust(loc = c(130.5686, -22.6523), ystart = 1990, yfinish = 1990,
#' surface_reflectivity = 0.15, slope = 0, aspect = 0, depths = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minimum_shade = 0, maximum_shade = 90,
#' local_height = 0.01, ...)
#' @export
#' @details
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{global_aerosol_database}{ = 1, Use the Global Aerosol Database? 1=yes (Fortran version), 2=yes (R version), 0=no}\cr\cr
#' \code{longwave_radiation_model}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{solar_model_only}{ = 0, Only run SOLRAD to get solar radiation? 1 = yes, 0 = no}\cr\cr
#' \code{radiation_per_wavelength}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{scattered_uv}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{max_iterations_per_day}{ = 3, iterations of first day to get a steady periodic}\cr\cr
#' \code{initial_soil_temperature}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{microclima}{ = 0, Use microclima and elevatr package to adjust solar radiation for terrain? 1 = yes, 0 = no}\cr\cr
#' \code{output_to_csv}{ = 0, Make Fortran code write output as csv files? 1 = yes, 0 = no}\cr\cr
#' \code{manualshade}{ = 1, Use CSIRO Soil and Landscape Grid of Australia? 1 = yes, 0 = no}\cr\cr
#' \code{soildata}{ = 0, Extract emissivities from gridded data? 1 = yes, 0 = no}\cr\cr
#' \code{terrain}{ = 0, Use 250m resolution terrain data? 1 = yes, 0 = no}\cr\cr
#' \code{dailywind}{ = 1, use McVicar grids for wind speed (does not work for opendap = 1)? 1 = yes, 0 = no}\cr\cr
#' \code{wind_multiplier}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{adiab_cor}{ = 1, use adiabatic lapse rate correction? 1 = yes, 0 = no}\cr\cr
#' \code{air_temperature_offset}{ = 0, warming offset vector, °C (negative values mean cooling). Can supply a single value or a vector the length of the number of days to be simulated.}\cr\cr
#' \code{spatial}{ = "c:/Australian Environment/", choose location of terrain data}\cr\cr
#' \code{opendap}{ = 1, query met grids via opendap}\cr\cr
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
#' \code{rainwet}{ = 1.5, mm of rainfall causing the ground to be 90\% wet for the day}\cr\cr
#' \code{organic_soil_cap}{ = 1, organic cap present on soil surface? (organic_soil_cap has lower conductivity - 0.2 W/mC - and higher specific heat 1920 J/kg-K)}\cr\cr
#' \code{precipitable_water}{ = 1, Precipitable cm H2O in air column, 0.1 = very dry; 1.0 = moist air conditions; 2.0 = humid, tropical conditions (note this is for the whole atmospheric profile, not just near the ground)}\cr\cr
#' \code{horizon_angles}{ = rep(0, 24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#' \code{minimum_temperature_lapse_rate}{ = 0.0039 Lapse rate for minimum air temperature (degrees C/m)}\cr\cr
#' \code{maximum_temperature_lapse_rate}{ = 0.0077 Lapse rate for maximum air temperature (degrees C/m)}\cr\cr
#' \code{maxima_times}{ = c(1.0, 1.0, 0.0, 0.0), Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise}\cr\cr
#' \code{minima_times}{ = c(0, 0, 1, 1), Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon}\cr\cr
#' \code{timezone}{ = 0, Use GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)? 1 = yes, 0 = no}\cr\cr
#'
#' \strong{ Soil moisture mode parameters:}
#'
#' \code{soil_moisture_model}{ = 1, Run soil moisture model? 1 = yes, 0 = no  1 = yes, 0 = no (note that this may cause slower runs)}\cr\cr
#' \code{air_entry_water_potential}{ = rep(1.1, 19), Air entry potential (J/kg) (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{saturated_hydraulic_conductivity}{ = rep(0.0037, 19), Saturated conductivity, (kg s/m3) (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{campbell_b_parameter}{ = rep(4.5, 19), Campbell's soil 'b' parameter (-) (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{soil_bulk_density}{ = rep(1.3, 19), Soil bulk density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{soil_mineral_density}{ = rep(2.56, 19), Soil density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter depths and points half way between)}\cr\cr
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{maximum_pooling_depth}{ = 10000, Max depth for water pooling on the surface (mm), to account for runoff}\cr\cr
#' \code{rain_multiplier}{ = 1, Rain multiplier for surface soil moisture (-), used to induce runon}\cr\cr
#' \code{evenrain}{ = 0, Spread daily rainfall evenly across 24hrs (1) or one event at midnight (0)}\cr\cr
#' \code{initial_soil_moisture}{ = c(0.1, 0.12, 0.15, 0.2, 0.25, 0.3, 0.3, 0.3, 0.3, 0.3), initial soil water content at each soil node, m3/m3}\cr\cr
#' \code{root_density}{ = c(0, 0, 8.2, 8.0, 7.8, 7.4, 7.1, 6.4, 5.8, 4.8, 4.0, 1.8, 0.9, 0.6, 0.8, 0.4, 0.4, 0, 0)*10000, root density (m/m3), (19 values descending through soil for specified soil nodes in parameter}\cr\cr
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
#' \code{snow_model}{ = 0, run the snow model 1 = yes, 0 = no (note that this may cause slower runs)}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (Mg/m3), overridden by densfun}\cr\cr
#' \code{densfun}{ = c(0, 0, 0, 0), slope and intercept of model of snow density as a linear function of snowpack age if first two values are nonzero, and following the exponential function of Sturm et al. 2010 J. of Hydromet. 11:1380-1394 if all values are non-zero; if it is c(0, 0, 0, 0) then fixed density used}\cr\cr
#' \code{snowmelt}{ = 1, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
#' \code{undercatch}{ = 1, undercatch multipier for converting rainfall to snow}\cr\cr
#' \code{rainmelt}{ = 0.0125, paramter in equation that melts snow with rainfall as a function of air temp}\cr\cr
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)}\cr\cr
#' \code{intercept}{ = max(maximum_shade) / 100 * 0.3, snow interception fraction for when there's shade (0-1)}\cr\cr
#' \code{grasshade}{ = 0, if 1, means shade is removed when snow is present, because shade is cast by grass/low shrubs}\cr\cr
#'
#' \strong{ Intertidal mode parameters:}
#'
#' \code{shore}{ Include tide effects? If 1, the matrix}
#' \code{tides}
#' { is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0, nrow = 24 \* 365 \* nyears, ncol = 3), matrix of 1. tide state (0 = out, 1 = in), 2. Water temperature (°C) and 3. Wave splash (0 = yes, 1 = no)}\cr\cr
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
#' \code{depths}{ - vector of depths used (cm)}\cr\cr
#' \code{diffuse_frac}{ - vector of hourly values of the fraction of total solar radiation that is diffuse (-), computed by microclima if microclima > 0}\cr\cr
#'
#' micromet_lowshade/micromet_highshade variables:
#' \itemize{
#' \item 1 day_of_year - day-of-year
#' \item 2 time - time of day (mins)
#' \item 3 air_temperature_local - air temperature (°C) at local height (specified by 'local_height' variable)
#' \item 4 air_temperature_reference - air temperature (°C) at reference height (specified by 'reference_height', 1.2m default)
#' \item 5 relative_humidity_local - relative humidity (\%) at local height (specified by 'local_height' variable)
#' \item 6 relative_humidity_reference  - relative humidity (\%) at reference height (specified by 'reference_height', 1.2m default)
#' \item 7 wind_speed_local - wind speed (m/s) at local height (specified by 'local_height' variable)
#' \item 8 wind_speed_reference - wind speed (m/s) at reference height (specified by 'reference_height', 1.2m default)
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
#'}
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
#'ystart <- 2014
#'yfinish <- 2015
#'nyears <- yfinish - ystart + 1
#'loc <- c(130.5686, -22.6523) # Nyrripi, Northern Territory, Australia
#'micro <- micro_aust(loc = loc, ystart = ystart, yfinish = yfinish, opendap = 1, elevation = 0, runshade = 0) # run the model for the middle of the desert in Australia, using opendap
#'
#' micromet_lowshade <- as.data.frame(micro$micromet_lowshade) # above ground microclimatic conditions, min shade
#' soil_temperature_lowshade <- as.data.frame(micro$soil_temperature_lowshade) # soil temperatures, minimum shade
#' soil_moisture_lowshade <- as.data.frame(micro$soil_moisture_lowshade) # soil temperatures, minimum shade
#'
#' # append dates
#' dates <- micro$dates
#' micromet_lowshade <- cbind(dates, micromet_lowshade)
#' soil_temperature_lowshade <- cbind(dates, soil_temperature_lowshade)
#' soil_moisture_lowshade <- cbind(dates, soil_moisture_lowshade)
#' minimum_shade <- micro$minimum_shade[1]
#'
#' # plotting above-ground conditions in minimum shade
#' with(micromet_lowshade, {plot(air_temperature_local ~ dates, xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l" ,main = paste("air and sky temperature, ", minimum_shade, "% shade", sep = ""), ylim = c(-20, 60))})
#' with(micromet_lowshade, {points(air_temperature_reference ~ dates, xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l", lty = 2, col = 'blue')})
#' with(micromet_lowshade,{points(sky_temperature ~ dates, xlab = "Date and Time", ylab = "Temperature (°C)"
#' ,  type = "l", col = 'light blue', main = paste("sky temperature, ", minimum_shade, "% shade", sep = ""))})
#' with(micromet_lowshade, {plot(relative_humidity_local ~ dates, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l", ylim = c(0, 100), main = paste("humidity, ", minimum_shade, "% shade", sep = ""))})
#' with(micromet_lowshade, {points(relative_humidity_reference ~ dates, xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l", col = 'blue', lty = 2, ylim = c(0, 100))})
#' with(micromet_lowshade, {plot(wind_speed_reference ~ dates, xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#' ,  type = "l", main = "wind speed", ylim = c(0, 15))})
#' with(micromet_lowshade, {points(wind_speed_local ~ dates, xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#' ,  type = "l", lty = 2, col = 'blue')})
#' with(micromet_lowshade, {plot(solar_radiation ~ dates, xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
#' ,  type = "l", main="solar radiation")})
#' with(micromet_lowshade, {plot(snow_depth ~ dates, xlab = "Date and Time", ylab = "Snow Depth (cm)"
#' ,  type = "l", main = "snow depth")})
#'
#' # plotting soil temperature
#' for(i in 1:10){
#'  if(i == 1){
#'    plot(soil_temperature_lowshade[, i + 3] ~ soil_temperature_lowshade[, 1], xlab = "Date and Time", ylab = "Soil Temperature (°C)"
#'    , col = i, type = "l", main = paste("soil temperature ", minimum_shade, "% shade",sep=""))
#'  }else{
#'    points(soil_temperature_lowshade[, i + 3] ~ soil_temperature_lowshade[, 1], xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)", col = i, type = "l")
#'  }
#' }
#'
#' # plotting soil moisture
#' for(i in 1:10){
#'  if(i == 1){
#'    plot(soil_moisture_lowshade[, i + 3] * 100 ~ soil_moisture_lowshade[, 1], xlab = "Date and Time", ylab = "Soil Moisture (% volumetric)"
#'    ,col = i, type = "l", main = paste("soil moisture ", minimum_shade, "% shade", sep = ""))
#'  }else{
#'    points(soil_moisture_lowshade[, i + 3] * 100 ~ soil_moisture_lowshade[, 1], xlab = "Date and Time", ylab = "Soil Moisture
#'     (%)", col = i, type = "l")
#'  }
#' }
micro_aust <- function(
  loc = c(130.5686, -22.6523),
  ystart = 1990,
  yfinish = 1990,
  nyears = 1,
  surface_reflectivity = 0.15,
  elevation = NA,
  slope = 0,
  aspect = 0,
  maximum_temperature_lapse_rate = 0.0077,
  minimum_temperature_lapse_rate = 0.0039,
  depths = c(0, 2.5, 5, 10, 15, 20, 30, 50, 100, 200),
  minimum_shade = 0,
  maximum_shade = 90,
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
  manualshade = 1,
  soildata = 0,
  terrain = 0,
  dailywind = 1,
  wind_multiplier = 1,
  adiab_cor = 1,
  air_temperature_offset = 0,
  spatial = "root_density:/",
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
  rainwet = 1.5,
  organic_soil_cap = 1,
  precipitable_water = 1,
  horizon_angles = rep(0, 24),
  maxima_times = c(1, 1, 0, 0),
  minima_times = c(0, 0, 1, 1),
  timezone = 0,
  soil_moisture_model = 1,
  air_entry_water_potential = rep(1.1, 19),
  saturated_hydraulic_conductivity = rep(0.0037, 19),
  campbell_b_parameter = rep(4.5, 19),
  soil_bulk_density = rep(bulk_density, 19),
  soil_mineral_density = rep(mineral_density, 19),
  maximum_pooling_depth = 10000,
  rain_multiplier = 1,
  evenrain = 0,
  initial_soil_moisture = c(0.1, 0.12, 0.15, 0.3, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4),
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
  snow_model = 0,
  snowtemp = 1.5,
  snowdens = 0.375,
  densfun = c(0.5979, 0.2178, 0.001, 0.0038),
  snowmelt = 1,
  undercatch = 1,
  rainmelt = 0.0125,
  shore = 0,
  tides = 0,
  scenario = "",
  year = "",
  rainhourly = 0,
  rainhour = 0,
  uid = "",
  pwd = "",
  host = "",
  radiation_per_wavelength = 0,
  scattered_uv = 0,
  max_iterations_per_day = 3,
  microclima = 0,
  microclima.dem_resolution = 100,
  microclima.zmin = -20,
  soilgrids = 0,
  longwave_radiation_model = 0,
  opendap = 1,
  message = 0,
  fail = nyears * 24 * 365,
  runmicro = 1,
  snowcond = 0,
  intercept = max(maximum_shade) / 100 * 0.3,
  grasshade = 0,
  maxsurf = 85
) {

  if(opendap == 0){
    require(RMySQL)
  }
  if(opendap == 1){
    if(is.na(elevation)){
      require(microclima)
      require(terra)
      cat('downloading DEM via package elevatr \n')
      dem <- microclima::get_dem(lat = loc[2], long = loc[1]) # mercator equal area projection
      xy = data.frame(lon = loc[1], lat = loc[2]) |>
        sf::st_as_sf(coords = c("lon", "lat"))
      xy <- sf::st_set_crs(xy, "EPSG:4326")
      xy <- sf::st_transform(xy, sf::st_crs(dem))
      elevation <- as.numeric(terra::extract(dem, xy))
      #xy <- data.frame(x = loc[1], y = loc[2])
      #coordinates(xy) = ~x + y
      #proj4string(xy) = "+init=epsg:4326"
      #xy <- as.data.frame(spTransform(xy, crs(dem)))
    }
    require(RNetCDF)
    ALTITUDES <- elevation
    dbrow <- 1
  }
  # --- 1. Input validation ---
  reference_height <- 1.2  # met data reference height (m) for AWAP
  errors <- 0
  # aust-specific checks
  if(opendap == 1 & (ystart < 1990 | ystart > 2017)){
    message("might not be data on the NCI for this time window")
  }
  if(is.numeric(loc[1])){
    if(loc[1] > 180 | loc[2] > 90){
      message("ERROR: Latitude or longitude (longlat) is out of bounds. Please enter a correct value.")
      errors <- 1
    }
  }
  if(timezone %in% c(0, 1) == FALSE){
    message("ERROR: the variable 'timezone' must be either 0 or 1. Please correct.")
    errors <- 1
  }
  if(global_aerosol_database == 1){
    message("If program is crashing, try global_aerosol_database = 2.")
  }
  if(write_input %in% c(0, 1) == FALSE){
    message("ERROR: the variable 'write_input' must be either 0 or 1. Please correct.")
    errors <- 1
  }
  # common parameter checks
  errors <- errors + validate_micro_inputs(depths, surface_reflectivity, slope, aspect, horizon_angles, surface_emissivity, tolerance,
    roughness_height, zero_plane_displacement, local_height, reference_height, orbital_eccentricity, precipitable_water, maxima_times, minima_times, minimum_shade, maximum_shade,
    mineral_conductivity, mineral_density, mineral_heat_capacity, bulk_density, global_aerosol_database = global_aerosol_database)

  if(errors == 0){ # continue

    ################## time related variables #################################
    nyears <- yfinish - ystart + 1
    yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
    tzone <- paste("Etc/GMT+", 10, sep = "")
    dates <- seq(ISOdate(ystart, 1, 1, tz = tzone) - 3600 * 12, ISOdate((ystart + nyears), 1, 1, tz = tzone) - 3600 * 13, by = "days")
    if(yfinish == as.numeric(format(Sys.time(), "%Y"))){ # cut days down if doing current year (incomplete)
      dates <- dates[dates < Sys.time()]
      dailywind <- 0 # no daily wind for current year
    }
    ndays <- length(dates)
    doys12 <- c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month
    microdaily <- 1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial max_iterations_per_day day burn in
    daystart <- 1
    if(length(minimum_shade) != ndays){
      minimum_shade_daily <- rep(0, ndays) + minimum_shade[1] # daily min shade (%)
    }else{
      minimum_shade_daily <- rep(0, ndays) + minimum_shade # daily min shade (%)
    }
    if(length(maximum_shade) != ndays){
      maximum_shade_daily <- rep(0, ndays) + maximum_shade[1] # daily max shade (%)
    }else{
      maximum_shade_daily <- rep(0, ndays) + maximum_shade # daily max shade (%)
    }
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
    start_day <- 1 # start day
    end_day <- ndays # end day
    curdate <- Sys.time() - 60 * 60 * 24
    curyear <- as.numeric(format(curdate, "%Y"))
    # location and terrain
    f1 <- paste0(spatial, "ausclim_rowids.nc")
    f2 <- paste0(spatial, "ausdem_shift1.tif")
    f3 <- paste0(spatial, "agg_9secdem.nc")
    f4 <- paste0(spatial, "Aust9secDEM.tif")

    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1], loc[2]))))

    cat("running micro_global to get clear sky solar \n")
    if(global_aerosol_database == 0){
      aerosol_optical_depth <- c(0.0670358341290886, 0.0662612704779235, 0.065497075238002, 0.0647431301168489, 0.0639993178022531, 0.0632655219571553, 0.0625416272145492, 0.0611230843885423, 0.0597427855962549, 0.0583998423063099, 0.0570933810229656, 0.0558225431259535, 0.0545864847111214, 0.0533843764318805, 0.0522154033414562, 0.0499736739981675, 0.047855059159556, 0.0458535417401334, 0.0439633201842001, 0.0421788036108921, 0.0404946070106968, 0.0389055464934382, 0.0374066345877315, 0.0359930755919066, 0.0346602609764008, 0.0334037648376212, 0.0322193394032758, 0.0311029105891739, 0.0300505736074963, 0.0290585886265337, 0.0281233764818952, 0.0272415144391857, 0.0264097320081524, 0.0256249068083005, 0.0248840604859789, 0.0241843546829336, 0.0235230870563317, 0.0228976873502544, 0.0223057135186581, 0.0217448478998064, 0.0212128934421699, 0.0207077699817964, 0.0202275105711489, 0.0197702578594144, 0.0193342605242809, 0.0189178697551836, 0.0177713140039894, 0.0174187914242432, 0.0170790495503944, 0.0167509836728154, 0.0164335684174899, 0.0161258546410128, 0.0158269663770596, 0.0155360978343254, 0.0152525104459325, 0.0149755299703076, 0.0147045436435285, 0.0144389973831391, 0.0141783930434343, 0.0134220329447663, 0.0131772403830191, 0.0129356456025128, 0.0126970313213065, 0.0124612184223418, 0.0122280636204822, 0.01199745718102, 0.0115436048739351, 0.0110993711778668, 0.0108808815754663, 0.0106648652077878, 0.0104513876347606, 0.0102405315676965, 0.00982708969547694, 0.00962473896278535, 0.00903679230300494, 0.00884767454432418, 0.0083031278398166, 0.00796072474935954, 0.00755817587626185, 0.00718610751850881, 0.00704629977586921, 0.00684663903049612, 0.00654155580333479, 0.00642947339729728, 0.00627223096874308, 0.00603955966866779, 0.00580920937536261, 0.00568506186880564, 0.00563167068287251, 0.00556222005081865, 0.00550522989971023, 0.00547395763028062, 0.0054478983436216, 0.00541823364504573, 0.00539532163908382, 0.00539239864119488, 0.00541690124712384, 0.00551525885358836, 0.00564825853509463, 0.00577220185074264, 0.00584222986640171, 0.00581645238345584, 0.00566088137411449, 0.00535516862329704, 0.00489914757707667, 0.00432017939770409, 0.0036813032251836, 0.00309019064543606, 0.00270890436501562, 0.00276446109239711, 0.00356019862584603)
    }else{
      aerosol_optical_depth <- 0
    }
    micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, aerosol_optical_depth = aerosol_optical_depth, timeinterval = 365, solar_model_only = 1)
    clearskyrad <- micro_clearsky$micromet_lowshade[, c(1, 13)]
    clearskysum <- aggregate(clearskyrad[, 2], by = list(clearskyrad[, 1]), FUN = sum)[, 2]

    # reference meridian for solar noon calculation
    solar_noon_longitude <- get_timezone_alref(lon = x[1], lat = x[2], timezone = timezone)
    hemisphere <- ifelse(x[2] < 0, 2, 1) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    latitude_degrees <- abs(trunc(x[2]))
    latitude_minutes <- (abs(x[2]) - latitude_degrees) * 60
    longitude_degrees <- abs(trunc(x[1]))
    longitude_minutes <- (abs(x[1]) - longitude_degrees) * 60
    azmuth <- aspect

    if(soildata == 0){
      soilprop <- cbind(0, 0)
    }

    if(soildata == 1){
      cat('extracting soil emissivities \n')
      static_soil <- paste0(spatial, "static_soil.nc")
      emissivities <- paste0(spatial, "aus_emissivities.nc")
      # read data in from netcdf file
      static_soil_data <- terra::rast(static_soil)
      static_soil_vars <- as.numeric(terra::extract(static_soil_data, x))
      labels <- c('albedo', 'FAPAR1', 'FAPAR2', 'FAPAR3', 'FAPAR4', 'FAPAR5', 'FAPAR6', 'FAPAR7', 'FAPAR8', 'FAPAR9', 'FAPAR10', 'FAPAR11', 'FAPAR12', 'volwater_Upper', 'volwater_lower', 'thick_upper', 'thick_lower', 'code')
      colnames(static_soil_vars) <- labels
      emissivities_data <- terra::rast(emissivities)
      SLES2 <- as.numeric(terra::extract(emissivities_data, x))

      # read in other soil related files for working out lumped soil type and properties
      # such as clay % for getting water potential
      filename <- paste0(spatial, "ppfInterpAll.txt")
      ppf <- as.data.frame(read.table(file = filename, sep = ",", header = TRUE))
      filename <- paste0(spatial,"Lumped soil types.txt")
      lumped.soil <- as.data.frame(read.table(file = filename, sep = ","))
      filename <- paste0(spatial, "SoilTypeLUT_725_AWAP.csv")
      soiltype <- as.data.frame(read.table(file = filename, sep = ","))
      soilcode <- subset(soiltype, soiltype[1] == static_soil_vars[18])
      lumped <- subset(lumped.soil, V4 == as.character(soilcode[1, 2]))
      soiltype <- lumped[1, 6]
      soilprop <- subset(ppf, ppf == soilcode[1, 2])
    }else{
      SLES2 <- rep(surface_emissivity, ndays)
      if(manualshade == 0){
        message("extracting shade data \n")
        static_soil <- paste0(spatial, "static_soil.nc")
        emissivities <- paste0(spatial, "aus_emissivities.nc")
        # read data in from netcdf file
        static_soil_data <- terra::rast(static_soil)
        static_soil_vars <- as.numeric(terra::extract(static_soil_data, x))
        labels <- c('albedo', 'FAPAR1', 'FAPAR2', 'FAPAR3', 'FAPAR4', 'FAPAR5', 'FAPAR6', 'FAPAR7', 'FAPAR8', 'FAPAR9', 'FAPAR10', 'FAPAR11', 'FAPAR12', 'volwater_Upper', 'volwater_lower', 'thick_upper', 'thick_lower', 'code')
        colnames(static_soil_vars) <- labels
      }
    }
    if(terrain == 1){
      message("extracting terrain data \n")
      e <- extent(x[1] - 0.05, x[1] + 0.05, x[2] - 0.05, x[2] + 0.05)
      for(i in 1:24){
        horifile <- paste0(spatial,'horizon', i, '.tif')
        horiz <- terra::crop(terra::rast(horifile), e)
        if(i == 1){
          horizons_data <- horiz
        }else{
          horizons_data <- terra::rast(horizons_data, horiz)
        }
      }
      HORIZONS <- as.numeric(t(terra::extract(horizons_data, x)))
      elev1 <- crop(terra::rast(paste0(spatial,'elevation.tif')), e)
      slope1 <- crop(terra::rast(paste0(spatial,'slope.tif')), e)
      aspect1 <- crop(terra::rast(paste0(spatial,'aspect.tif')), e)
      elevslpasp <- terra::rast(elev1, slope1, aspect1)
      ELEVSLPASP <- as.numeric(terra::extract(elevslpasp, x))
      ELEVSLPASP <- as.matrix((ifelse(is.na(ELEVSLPASP), 0, ELEVSLPASP)))
      ALTITUDES <- ELEVSLPASP[, 1]
      SLOPES <- ELEVSLPASP[, 2]
      AZMUTHS <- ELEVSLPASP[, 3]
      # the horizons have been arranged so that they go from 0 degrees azimuth (north) clockwise - r.horizon starts
      # in the east and goes counter clockwise!
      HORIZONS <- (ifelse(is.na(HORIZONS), 0, HORIZONS)) / 10 # get rid of na and get back to floating point
      HORIZONS <- data.frame(HORIZONS)
      VIEWF_all <- 1 - rowSums(sin(t(HORIZONS) * pi / 180)) / length(t(HORIZONS)) # convert horizon angles to radians and calc view factor(s)
      r1 <- terra::rast(f1)
      r2 <- terra::rast(f2)
      r3 <- terra::rast(f3)
      dbrow <- as.numeric(terra::extract(r1, x))
      AUSDEM <- as.numeric(terra::extract(r2, x))
      AGG <- as.numeric(terra::extract(r3, x))
    }else{
      if(opendap == 0){
        r1 <- terra::rast(f1)
        r2 <- terra::rast(f2)
        r3 <- terra::rast(f3)
        r4 <- terra::rast(f4)
        dbrow <- as.numeric(terra::extract(r1, x))
        AUSDEM <- as.numeric(terra::extract(r2, x))
        AGG <- as.numeric(terra::extract(r3, x))
        if(is.na(elevation) == FALSE){ # check if user-specified elevation
          ALTITUDES <- elevation # use user-specified elevation
        }else{
          ALTITUDES <- as.numeric(terra::extract(r4, x)) # get elevation from fine res DEM
        }
      }
      HORIZONS <- horizon_angles
      HORIZONS <- data.frame(HORIZONS)
      VIEWF_all <- 1 - sum(sin(as.data.frame(horizon_angles) * pi / 180)) / length(horizon_angles) # convert horizon angles to radians and calc view factor(s)
      SLOPES <- rep(slope, length(x[, 1]))
      AZMUTHS <- rep(aspect, length(x[, 1]))
    }
    horizon_angles <- HORIZONS
    row.names(horizon_angles) <- NULL
    horizon_angles <- as.numeric(as.matrix(horizon_angles))

    if(soildata == 1){
      VIEWF <- VIEWF_all
      surface_emissivity <- SLES2
    }else{
      VIEWF <- VIEWF_all
    }

    if(soilgrids == 1){
      sg <- fetch_soilgrids(x, depths)
      if (!is.null(sg)) { air_entry_water_potential <- sg$air_entry_water_potential; saturated_hydraulic_conductivity <- sg$saturated_hydraulic_conductivity; campbell_b_parameter <- sg$campbell_b_parameter; soil_bulk_density <- sg$soil_bulk_density; bulk_density <- sg$bulk_density }
    }
    # setting up for temperature correction using lapse rate given difference between 9sec DEM value and 0.05 deg value
    if(opendap == 0){
      if(AUSDEM == -9999 | is.na(AUSDEM) == 'TRUE'){
        delta_elev = AGG - ALTITUDES
      }else{
        delta_elev = AUSDEM - ALTITUDES
      }
      adiab_corr_max <- delta_elev * maximum_temperature_lapse_rate
      adiab_corr_min <- delta_elev * minimum_temperature_lapse_rate
    }else{
      adiab_corr_max <- 0
      adiab_corr_min <- 0
    }
    if(scenario!=""){
      message("generate climate change scenario", '\n')
      # diff spline function
      getdiff <- function(diffs, grid){
        diff1 <- (unlist(diffs[1]) + unlist(diffs[12])) / 2

        # generate list of days
        for(ys in 1:nyears){
          if(yearstodo[ys] %in% leapyears){
            day<-c(1, 15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 366)
          }else{
            day<-c(1, 15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 365)
          }
          if(ys == 1){
            days2 <- day
            days <- day
          }else{
            if(yearstodo[ys] %in% leapyears){
              days2 <- c(days2, (day + 366 * (ys - 1)))
            }else{
              days2 <- c(days2, (day +365 * (ys - 1)))
            }
            days <- c(days, day)
          }
        }

        if(is.na(diffs[1]) == TRUE){
          # find the nearest cell with data
          NArem <- grid[[1]]
          NArem <- Which(!is.na(NArem), cells = TRUE)
          dist <- distanceFromPoints(maxTst05[[1]], x)
          distNA <- as.numeric(terra::extract(dist, NArem))
          cellsR <- cbind(distNA, NArem)
          distmin <- which.min(distNA)
          cellrep <- cellsR[distmin, 2]
          diffs <- as.numeric(terra::extract(maxTst05, cellrep))
          diff1 <- (unlist(diffs[1]) + unlist(diffs[12])) / 2
        }
        diffs3 <- rep(c(diff1, diffs, diff1), nyears)
        days_diffs <- data.frame(matrix(NA, nrow = nyears * 14, ncol = 3))
        days_diffs[, 1] <- days
        days_diffs[, 3] <- days2
        days_diffs[, 2] <- diffs3
        colnames(days_diffs) <- c("days", "diffs", "new_day")

        # interpolate monthly differences
        f <- approxfun(x = days_diffs$new_day, y = days_diffs$diffs)
        xx <- seq(1, max(days2), 1)
        sp_diff <- f(xx)
        return(sp_diff)
      }
      # Max and Min Air Temps
      ccdir <- "/vlsci/VR0212/mrke"
      load(file = paste0(ccdir, "/Australia Climate Change/", scenario, "/", "maxTst05_", scenario, "_", year, ".Rda")) #maxTst05
      diffs <- as.numeric(terra::extract(maxTst05,  x))
      TMAXX_diff <- getdiff(diffs, maxTst05)
      load(file = paste0(ccdir, "/Australia Climate Change/", scenario, "/", "minTst05_", scenario, "_", year, ".Rda")) #minTst05
      diffs <- as.numeric(terra::extract(minTst05, x))
      TMINN_diff <- getdiff(diffs, minTst05)
      # RH
      load(file = paste0(ccdir,  "/Australia Climate Change/", scenario, "/", "RHst05_", scenario, "_", year, ".Rda")) #maxTst05
      diffs <- as.numeric(terra::extract(RHst05, x)[, 2])
      RH_diff <- getdiff(diffs, RHst05)
      # wind
      load(file = paste0(ccdir, "/Australia Climate Change/", scenario, "/", "PT_VELst05_", scenario, "_", year, ".Rda"))
      diffs <- as.numeric(terra::extract(PT_VELst05, x))
      WIND_diff <- getdiff(diffs, PT_VELst05)
      # SOLAR/CLOUD COVER
      load(file = paste0("c:/Spatial_Data/Australia Climate Change/", scenario, "/", "SOLCst05_", scenario, "_", year, ".Rda"))
      diffs <- as.numeric(terra::extract(SOLCst05, x))
      SOLAR_diff <- getdiff(diffs, SOLCst05)
    }

    if(opendap == 1){

      message("extracting climate data via opendap - note that there is no wind speed data, so the daily range is assumed to be from 0.5 to 2 m/s \n")
      baseurl<- "http://opendap.bom.gov.au:8080/thredds/dodsC/agcd/"
      get_AGCD <- function(dstart, dfinish, var, x, baseurl = "http://opendap.bom.gov.au:8080/thredds/dodsC/agcd/") {
        # Convert date strings to Date objects
        tstart <- as.Date(dstart, format = "%d/%m/%Y")
        tfinish <- as.Date(dfinish, format = "%d/%m/%Y")

        # Generate sequence of dates
        dates <- seq(tstart, tfinish, by = "1 day")

        # Precompute file paths
        if(var == 'solar'){
          filepaths <- paste0(
            baseurl, var, "_2018/total/r005/01day/",
            format(dates, "%Y"), "/",
            sprintf("%s_dailyExposure.nc", format(dates, "%Y%m%d"))
          )
          var2 <- "dailyExposure"
        }else{
          if(var == 'precip'){
            filepaths <- paste0(
              baseurl, var, "/total/r005/01day/",
              format(dates, "%Y"), "/",
              sprintf("%s_total_r005_%s_%s.nc", var,
                      format(dates, "%Y%m%d"), format(dates, "%Y%m%d"))
            )
            var2 <- "precip"
          }else{
            filepaths <- paste0(
              baseurl, var, "/mean/r005/01day/",
              format(dates, "%Y"), "/",
              sprintf("%s_mean_r005_%s_%s.nc", var,
                      format(dates, "%Y%m%d"), format(dates, "%Y%m%d"))
            )
            var2 <- var
          }
        }
        if(var == 'vapourpres_h09' | var == 'vapourpres_h15'){
          var2 <- "vapourpres"
        }
        cat(paste0("Reading weather input for variable: ", var, "\n"))

        # Helper function to extract the value from one file
        extract_value <- function(ncfile) {
          nc <- RNetCDF::open.nc(ncfile)
          if(length(grep("dailyExposure", ncfile)) != 0){
            lon <- RNetCDF::var.get.nc(nc, "longitude", unpack = TRUE)
            lat <- RNetCDF::var.get.nc(nc, "latitude", unpack = TRUE)
          }else{
            lon <- RNetCDF::var.get.nc(nc, "lon", unpack = TRUE)
            lat <- RNetCDF::var.get.nc(nc, "lat", unpack = TRUE)
          }
          # Find nearest index (efficient)
          latindex <- which.min(abs(lat - x[2]))
          lonindex <- which.min(abs(lon - x[1]))

          start <- c(lonindex, latindex, 1)
          count <- c(1, 1, NA)
          value <- as.numeric(RNetCDF::var.get.nc(nc, variable = var2,
                                                  start = start, count = count,
                                                  unpack = TRUE))
          RNetCDF::close.nc(nc)
          return(value)
        }

        # Read all data efficiently (can also try parallel later)
        data <- vapply(filepaths, extract_value, numeric(1))

        return(data)
      }
      dstart <- paste0('01/01/', ystart)
      dfinish <- paste0('31/12/', yfinish)
      sol <- zoo::na.approx(get_AGCD(dstart, dfinish, "solar", x))
      tmax <- get_AGCD(dstart, dfinish, "tmax", x)
      tmin <- get_AGCD(dstart, dfinish, "tmin", x)
      vph09 <- get_AGCD(dstart, dfinish, "vapourpres_h09", x)
      vph15 <- get_AGCD(dstart, dfinish, "vapourpres_h15", x)
      Rain <- get_AGCD(dstart, dfinish, "precip", x)

      allclearsky <- leapfix(clearskysum, yearlist)
      allclearsky <- allclearsky[1:ndays]
      # convert from W/d to MJ/d
      allclearsky <- allclearsky * 3600 / 1e6
      cloud <- (1 - sol / allclearsky) * 100
      cloud[cloud < 0] <- 0
      cloud[cloud > 100] <- 100
      cloud_max <- as.numeric(cloud)
      cloud_min <- cloud_max
      cloud_min <- cloud_min * 0.5
      cloud_max <- cloud_max * 2
      cloud_min[cloud_min > 100] <- 100
      cloud_max[cloud_max > 100] <- 100
      reference_temperature_max <- tmax
      reference_temperature_min <- tmin
      VAPRES_max <- apply(cbind(vph09, vph15), FUN = max, MARGIN = 1) * 100 # convert from hectopascals to pascals
      VAPRES_min <- apply(cbind(vph09, vph15), FUN = min, MARGIN = 1) * 100 # convert from hectopascals to pascals
      es <- WETAIR(db = reference_temperature_max, rh = 100)$esat
      reference_humidity_min <- (VAPRES_min / es) * 100
      reference_humidity_min[reference_humidity_min > 100] <- 100
      reference_humidity_min[reference_humidity_min < 0] <- 0.01
      es <- WETAIR(db = reference_temperature_min, rh = 100)$esat
      reference_humidity_max <- (VAPRES_max / es) * 100
      reference_humidity_max[reference_humidity_max > 100] <- 100
      reference_humidity_max[reference_humidity_max < 0] <- 0.01
      reference_wind_max <- rep(2, ndays)
      reference_wind_min <- rep(0.5, ndays)
      rainfall <- Rain
    }

    # connect to server
    if(opendap == 0){
      #channel2 <- RODBC::odbcConnect("ausclim_predecol", uid = uid, pwd = pwd)
      #channel <- RODBC::odbcConnect("AWAPDaily", uid = uid, pwd = pwd)
      lat1 <- x[2] - 0.024
      lat2 <- x[2] + 0.025
      lon1 <- x[1] - 0.024
      lon2 <- x[1] + 0.025
      channel <- RMySQL::dbConnect(MySQL(), user = uid, password = pwd, host = host, dbname = "AWAPDaily", port = 3306)
      # preliminary test for incomplete year, if simulation includes the present year
      for(j in 1:nyears){ # start loop through years
        yeartodo <- yearlist[j]
        query <- paste("SELECT a.latitude, a.longitude, b.* FROM AWAPDaily.latlon as a
                       , AWAPDaily.", yeartodo," as b where (a.id = b.id) and
                       (a.latitude between ", lat1, " and ",lat2,") and (a.longitude between ", lon1, " and ", lon2, ")
                       order by b.day", sep = "")
        output<- dbGetQuery(channel, query)
        if(yearlist[j] < 1971){
          output$vpr <- output$tmin/output$tmin-1
        }
        if(yearlist[j] > 1989){
          output$sol <- as.numeric(as.character(output$sol))
        }else{
          output$sol <- output$tmin / output$tmin - 1
        }
        if(j == 1){
          results <- output
        }else{
          results <- rbind(results, output)
        }
      }
      dbDisconnect(channel)
    }
    doys<-seq(daystart, ndays, 1)
    leapyears <- seq(1900, 2060, 4)
    for(mm in 1:nyears){
      if(mm == 1){
        currenty <- ystart
      }else{
        currenty <- ystart + mm
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
    day_of_year <- leapfix(day_of_year, yearlist)
    day_of_year <- day_of_year[1:ndays]
    end_day <- ndays
    start_day <- 1 # start month
    # end preliminary test for incomplete year, if simulation includes the present year

    if((soildata == 1 & nrow(soilprop) > 0) | soildata == 0){
      if(soildata == 1){
        # get static soil data into arrays
        surface_reflectivity <- static_soil_vars[, 1]  # albedo/reflectances
        maxshades <- static_soil_vars[, 2:13] # assuming FAPAR represents shade
        shademax <- maxshades
      }else{
        if(manualshade == 0){
          maxshades <- static_soil_vars[, 2:13] # assuming FAPAR represents shade
        }
        shademax <- maximum_shade
      }
      if((is.na(dbrow) != TRUE & is.na(ALTITUDES) != TRUE) | opendap == 1){
        aerosol_optical_depth <- compute_tai(longlat, global_aerosol_database, TAI_AUSTRALIA)

        if(opendap == 0){
          channel <- RMySQL::dbConnect(MySQL(), user = uid, password = pwd, host = host, dbname = "AWAPDaily", port = 3306)
          for(j in 1:nyears){ # start loop through years
            yeartodo<-yearlist[j]
            lat1 <- x[2] - 0.024
            lat2 <- x[2] + 0.025
            lon1 <- x[1] - 0.024
            lon2 <- x[1] + 0.025
            query<-paste("SELECT a.latitude, a.longitude, b.*
                       FROM AWAPDaily.latlon as a
                       , AWAPDaily.", yeartodo, " as b
                       where (a.id = b.id) and (a.latitude between ", lat1, " and ", lat2, ") and (a.longitude between ",lon1," and ",lon2,")
                       order by b.day", sep = "")
            output<- dbGetQuery(channel, query)
            if(yearlist[j] < 1971){
              output$vpr <- output$tmin / output$tmin - 1
            }
            if(yearlist[j] > 1989){
              output$sol <- as.numeric(as.character(output$sol))
            }else{
              output$sol <- output$tmin / output$tmin - 1
            }
            if(j==1){
              results <- output
            }else{
              results <- rbind(results, output)
            }
          }
          dbDisconnect(channel)
          if(dailywind == 1){
            channel3 <- RMySQL::dbConnect(MySQL(), user = uid, password = pwd, host = host, dbname = "dailywind", port = 3306)
            if(min(yearlist) < 1975){
              # get mean of 1975-1984
              for(j in 1:10){ # start loop through years
                yeartodo <- 1974 + j
                lat1 <- x[2] - 0.024
                lat2 <- x[2] + 0.025
                lon1 <- x[1] - 0.024
                lon2 <- x[1] + 0.025
                query<-paste("SELECT a.latitude, a.longitude, b.*
                         FROM dailywind.latlon as a
                         , dailywind.", yeartodo, " as b
                         where (a.id = b.id) and (a.latitude between ",lat1," and ",lat2,") and (a.longitude between ",lon1," and ",lon2,")
                         order by b.day", sep = "")
                output <- dbGetQuery(channel3, query)
                if(j == 1){
                  dwindmean <- output
                }else{
                  dwindmean <- cbind(dwindmean, output[, 5])
                }
              }
              dwindmean<-cbind(dwindmean[, 1:4], rowMeans(dwindmean[, 5:14]))
              colnames(dwindmean)[5] <- 'wind'
            }
            for(j in 1:nyears){ # start loop through years
              yeartodo <- yearlist[j]
              if(yeartodo < 1975){
                output <- dwindmean
              }else{
                lat1 <- x[2] - 0.024
                lat2 <- x[2] + 0.025
                lon1 <- x[1] - 0.024
                lon2 <- x[1] + 0.025
                query <- paste("SELECT a.latitude, a.longitude, b.*
                         FROM dailywind.latlon as a
                             , dailywind.", yeartodo, " as b
                             where (a.id = b.id) and (a.latitude between ", lat1, " and ", lat2, ") and (a.longitude between ", lon1, " and ", lon2, ")
                             order by b.day", sep = "")
                output <- dbGetQuery(channel3, query)
                }
              if(j == 1){
                dwind <- output
              }else{
                dwind <- rbind(dwind, output)
              }
            }
            dwind <- dwind$wind / 15.875 # conversion byte (i.e., an 8-bit unsigned integer ranging in value from 0 to 255) to m/s
            dbDisconnect(channel3)
          }
        }
        if(opendap == 0){
          if(adiab_cor==1){
            reference_temperature_max.orig <- results$tmax
            reference_temperature_min.orig <- results$tmin
            reference_temperature_max <- as.matrix(results$tmax + as.numeric(adiab_corr_max))
            reference_temperature_min <- as.matrix(results$tmin + as.numeric(adiab_corr_min))

          }else{
            reference_temperature_max <- as.matrix(results$tmax)
            reference_temperature_min <- as.matrix(results$tmin)
          }
          if(scenario!=""){
            reference_temperature_max <- reference_temperature_max + TMAXX_diff
            reference_temperature_min <- reference_temperature_min + TMINN_diff
          }
          rainfall <- results$rr
          output_AWAPDaily <- results
        }else{
          if(adiab_cor == 1){
            reference_temperature_max.orig <- reference_temperature_max
            reference_temperature_min.orig <- reference_temperature_min
            reference_temperature_max<-as.matrix(reference_temperature_max + adiab_corr_max)
            reference_temperature_min<-as.matrix(reference_temperature_min + adiab_corr_min)
          }
        }
        if(scenario != ""){
          # first work out for each site the new predicted rainfall amount for each month - use this to adjust for fact that will underestimate chcange
          # using proportion because 0 x % is still 0
          # add columns with days, months and years
          RAIN_current <- as.data.frame(rainfall)
          dates <- seq(ISOdate(ystart, 1, 1, tz = paste("Etc/GMT+", 10, sep="")) - 3600 * 12, ISOdate((ystart+nyears),1, 1, tz = paste("Etc/GMT+",10, sep=""))-3600*13, by="days")
          dates <- subset(dates, format(dates, "%m/%d") != "02/29") # remove leap years
          RAINFALL_sum <- aggregate(RAIN_current, by = list(format(dates, "%m-%Y")), FUN = sum)
          dates2 <- RAINFALL_sum$Group.1
          RAINFALL_sum <- RAINFALL_sum[order(as.Date(paste0("01-",RAINFALL_sum$Group.1), "%m-%Y")), 2]

          load(file = paste0("c:/Spatial_Data/Australia Climate Change/", scenario, "/", "RAINst05_", scenario, "_", year, ".Rda"))
          diffs <- rep(as.numeric(terra::extract(RAINst05, x)), nyears)

          if(is.na(diffs[1]) == TRUE){
            print("no data")
            # find the nearest cell with data
            NArem <-RAINst05[[1]] # don't need to re-do this within bioregion loop
            NArem <- Which(!is.na(NArem), cells = TRUE)
            dist <- distanceFromPoints(RAINst05[[1]], y)
            distNA <- as.numeric(terra::extract(dist, NArem))
            cellsR <- cbind(distNA, NArem)
            distmin <- which.min(distNA)
            cellrep <- cellsR[distmin, 2]
            diffs <- rep(as.numeric(terra::extract(RAINst05, cellrep)), nyears)
          }

          rainfall_new <- (RAINFALL_sum * diffs)

          rainfall_new[rainfall_new < 0 ] <- 0 # get rid of any negative rainfall values

          ## Now extract predicted change in mm
          load(file = paste0("c:/Spatial_Data/Australia Climate Change/", scenario, "/" ,"RAINst05_mm_", scenario, "_", year, ".Rda"))
          rainfall_change_mm <- rep(as.numeric(terra::extract(RAINst05_mm, x)), nyears)

          #########Now get predicted change in rainfall (could also get this from OzClim or ClimSim layer)#############
          Diff_prop <- rainfall_new/RAINFALL_sum # proportion change
          Diff_prop[Diff_prop == 'NaN'] <- 0
          Diff_prop[Diff_prop == 'Inf'] <- 0 ## If was previously no rainfall and now is rainfall need to alter code so this is simply added

          newRAINFALL <- rep(NA, length(rainfall))
          for (k in 1:length(rainfall)){ # loop through each sites applying monthly % changes
            month <- which(dates2 == format(dates[k], "%m-%Y"))
            # Test if predicted rainfall matches up - use rainfall_change_mm layer
            Rain_adj <- rainfall_change_mm[month]

            # test for if proportional change is 0 (because current rainfall is 0)
            #but rainfall predicted to increase
            if(Diff_prop[month] == 0 & Rain_adj > 1){ # couldn't get proportion as no current rain days but need to add rain
              print('new rain days needed')
              # previously no rain, randomly select a day and put all rain on it
              listD<-seq(1, length(newRAINFALL), 1)
              altD<-sample(listD, 1)
              newRAINFALL[altD] <- Rain_adj / 30
            }else{
              newRAINFALL[k] <- rainfall[k] * Diff_prop[month]
            }
          } # end of loop through each day
          newRAINFALL[newRAINFALL < 0.1] <- 0
          newRAINFALL[is.na(newRAINFALL)] <- 0
          rainfall <- newRAINFALL
        }
        # cloud cover
        if(opendap == 0){
          if(ystart > 1989 & sum(results[, 9], na.rm = TRUE) > 0){ # solar radiation data available
            allclearsky <- leapfix(clearskysum, yearlist)
            allclearsky <- allclearsky[1:ndays]
            # convert from W/d to MJ/d
            allclearsky <- allclearsky * 3600 / 1e6
            if(is.na(output_AWAPDaily[1, 9]) == TRUE){
              output_AWAPDaily[1, 9] = mean(output_AWAPDaily[, 9], na.rm = TRUE)
            }
            if(is.na(output_AWAPDaily[ndays, 9]) == TRUE){
              output_AWAPDaily[nrow(output_AWAPDaily), 9]=mean(output_AWAPDaily[, 9], na.rm = TRUE)
            }
            solar <- zoo::na.approx(output_AWAPDaily[, 9])
            if(scenario != ""){
              solar <- solar*SOLAR_diff
            }
            cloud <- (1 - solar / allclearsky) * 100
            cloud[cloud < 0] <- 0
            cloud[cloud > 100] <- 100
            cloud_max <- as.numeric(cloud)
            cloud_min <- cloud_max
          }else{
            datestart1 <- "01/01/1990" # day, month, year
            datefinish1 <- "31/12/2014" # day, month, year
            datestart1 <- strptime(datestart1, "%d/%m/%Y") # convert to date format
            datefinish1 <- strptime(datefinish1, "%d/%m/%Y") # convert to date format
            yearstart1 <- as.numeric(format(datestart1, "%Y")) # yet year start
            yearfinish1 <- as.numeric(format(datefinish1, "%Y")) # yet year finish
            years1 <- seq(yearstart1, yearfinish1, 1) # get sequence of years to d0
            doystart <- datestart1$yday + 1 # get day-of-year at start
            doyfinish <- datefinish1$yday + 1 # get day-of-year at finish
            years1 <- seq(yearstart1, yearfinish1, 1) # get sequence of years to do
            channel <- RMySQL::dbConnect(MySQL(), user = uid, password = pwd, host = host, dbname = "AWAPDaily", port = 3306)
            for(i in 1:length(years1)){ # start loop through years
              # syntax for query
              if(length(years1) == 1){ # doing a period within a year
                query<-paste0("SELECT a.latitude, a.longitude, b.* FROM AWAPDaily.latlon as a
                    , AWAPDaily.", years1[i], " as b where (a.id = b.id) and (a.latitude between ", lat1, " and ",lat2, ") and (a.longitude between ", lon1, " and ",lon2, ") and (b.day between ",doystart, " and ",doyfinish, ")
                              order by b.day")
               }else{
                if(i==1){ # doing first year, start at day requested
                  query<-paste0("SELECT a.latitude, a.longitude, b.* FROM AWAPDaily.latlon as a
                     , AWAPDaily.", years1[i], " as b where (a.id = b.id) and (a.latitude between ", lat1, " and ",lat2, ") and (a.longitude between ", lon1, " and ", lon2, ") and (b.day >= ", doystart,")
                     order by b.day")
                 }else{
                  if(i==length(years1)){ # doing last year, only go up to last day requested
                    query<-paste0("SELECT a.latitude, a.longitude, b.* FROM AWAPDaily.latlon as a
                       , AWAPDaily.", years1[i], " as b where (a.id = b.id) and (a.latitude between ", lat1, " and ",lat2, ") and (a.longitude between ", lon1," and ", lon2, ") and (b.day <= ", doyfinish,")
                       order by b.day")
                  }else{ # doing in between years, so get all data for this year
                    query<-paste0("SELECT a.latitude, a.longitude, b.* FROM AWAPDaily.latlon as a
                       , AWAPDaily.", years1[i], " as b where (a.id = b.id) and (a.latitude between ", lat1, " and ", lat2, ") and (a.longitude between ", lon1," and ", lon2, ")
                      order by b.day")
                  }}}
              if(i==1){
                output1 <- dbGetQuery(channel, query)
              }else{
                output1 <- rbind(output1, dbGetQuery(channel, query))
              }
            } # end loop through years
            dbDisconnect(channel)
            output1$sol <- as.numeric(output1$sol)
            output1$clearsky <- leapfix(clearskysum, seq(1990, 2014)) * 3600 / 1e6
            glm_sol <- coefficients(with(output1, glm(sol ~ rr + tmax + tmin + day + clearsky)))
            output_AWAPDaily$clearsky <- leapfix(clearskysum, yearlist) * 3600 / 1e6
            output_AWAPDaily[, 9] <- glm_sol[1] + glm_sol[2] * output_AWAPDaily$rr + glm_sol[3] * output_AWAPDaily$tmax + glm_sol[4] * output_AWAPDaily$tmin + glm_sol[5] * output_AWAPDaily$day + glm_sol[6] * output_AWAPDaily$clearsky
            if(scenario != ""){
              output_AWAPDaily[, 9] <- output_AWAPDaily[, 9] * SOLAR_diff
            }
            cloud <- (1 - as.data.frame(output_AWAPDaily$sol) / as.data.frame(output_AWAPDaily$clearsky)) * 100
            cloud[cloud < 0] <- 0
            cloud[cloud > 100] <- 100
            cloud <- as.matrix(cbind(output_AWAPDaily[, 4], cloud))
            cloud_max <- cloud[, 2]
            cloud_min <- cloud_max
          }# end check for year 1990 or later
          if(ystart > 1970){ #vapour pressure data available
            if(is.na(output_AWAPDaily[1, 8]) == TRUE){
              output_AWAPDaily[1, 8] = mean(output_AWAPDaily[, 8], na.rm = TRUE)
            }
            VAPRES <- zoo::na.approx(output_AWAPDaily[, 8])
            VAPRES <- VAPRES * 100 # convert from hectopascals to pascals
            es <- WETAIR(db = reference_temperature_max, rh = 100)$esat
            reference_humidity_min <- (VAPRES / es) * 100
            reference_humidity_min[reference_humidity_min > 100] <- 100
            reference_humidity_min[reference_humidity_min < 0] <- 0.01
            es <- WETAIR(db = reference_temperature_min, rh = 100)$esat
            reference_humidity_max <- (VAPRES / es) * 100
            reference_humidity_max[reference_humidity_max > 100] <- 100
            reference_humidity_max[reference_humidity_max < 0] <- 0.01
            if(scenario != "" ){
              reference_humidity_min <- reference_humidity_min + RH_diff
            }
            if(scenario != ""){
              reference_humidity_max <- reference_humidity_max + RH_diff
            }
          }else{
            if(exists("output1") == FALSE){
              datestart1 <- "01/01/1990" # day, month, year
              datefinish1 <- "31/12/2014" # day, month, year
              datestart1 <- strptime(datestart1, "%d/%m/%Y") # convert to date format
              datefinish1 <- strptime(datefinish1, "%d/%m/%Y") # convert to date format
              yearstart1 <- as.numeric(format(datestart1, "%Y")) # yet year start
              yearfinish1 <- as.numeric(format(datefinish1, "%Y")) # yet year finish
              years1 <- seq(yearstart1, yearfinish1, 1) # get sequence of years to d0
              doystart <- datestart1$yday + 1 # get day-of-year at start
              doyfinish <- datefinish1$yday + 1 # get day-of-year at finish
              years1 <- seq(yearstart1, yearfinish1, 1) # get sequence of years to do
              channel <- RMySQL::dbConnect(MySQL(), user = uid, password = pwd, host = host, dbname = "AWAPDaily", port = 3306)
              for(i in 1:length(years1)){ # start loop through years
                 if(length(years1) == 1){ # doing a period within a year
                  query <- paste0("SELECT a.latitude, a.longitude, b.* FROM AWAPDaily.latlon as a
                    , AWAPDaily.", years1[i], " as b where (a.id = b.id) and (a.latitude between ",
                                lat1, " and ", lat2, ") and (a.longitude between ", lon1, " and ", lon2, ") and (b.day between
                    ", doystart, " and ", doyfinish, ") order by b.day")
                }else{
                  if(i == 1){ # doing first year, start at day requested
                    query <- paste0("SELECT a.latitude, a.longitude, b.* FROM AWAPDaily.latlon as a
                       , AWAPDaily.", years1[i]," as b where (a.id = b.id) and (a.latitude between ",
                                  lat1, " and ", lat2, ") and (a.longitude between ", lon1, " and ", lon2, ") and (b.day >= ",
                                  doystart, ") order by b.day")
                  }else{
                    if(i == length(years1)){ # doing last year, only go up to last day requested
                      query <- paste0("SELECT a.latitude, a.longitude, b.* FROM AWAPDaily.latlon as a
                        , AWAPDaily.", years1[i]," as b where (a.id = b.id) and (a.latitude between "
                                    , lat1, " and ", lat2, ") and (a.longitude between ", lon1, " and ", lon2, ") and (b.day <= ",
                                    doyfinish, ") order by b.day")
                    }else{ # doing in between years, so get all data for this year
                      query <- paste0("SELECT a.latitude, a.longitude, b.* FROM AWAPDaily.latlon as a
                        , AWAPDaily.", years1[i], " as b where (a.id = b.id) and (a.latitude between
                        ", lat1, " and ", lat2, ") and (a.longitude between ", lon1, " and ", lon2, ") order by b.day")
                    }}}
                if(i==1){
                  output1 <- dbGetQuery(channel, query)
                }else{
                  output1 <- rbind(output1, dbGetQuery(channel, query))
                }
              } # end loop through years
              dbDisconnect(channel)
            }
            glm_vpr <- coefficients(with(output1, glm(vpr ~ rr + tmax + tmin + day)))
            output_AWAPDaily[, 8]<- glm_vpr[1] + glm_vpr[2] * output_AWAPDaily$rr + glm_vpr[3] * output_AWAPDaily$tmax + glm_vpr[4] * output_AWAPDaily$tmin + glm_vpr[5] * output_AWAPDaily$day
            VAPRES <- zoo::na.approx(output_AWAPDaily[, 8])
            VAPRES<-VAPRES * 100 # convert from hectopascals to pascals
            # correct for potential change in RH with elevation-corrected Tair
            es <- WETAIR(db = reference_temperature_max, rh = 100)$esat
            reference_humidity_min <- (VAPRES / es) * 100
            reference_humidity_min[reference_humidity_min > 100] <- 100
            reference_humidity_min[reference_humidity_min < 0] <- 0.01
            es <- WETAIR(db = reference_temperature_min, rh = 100)$esat
            reference_humidity_max <- (VAPRES / es) * 100
            reference_humidity_max[reference_humidity_max > 100] <- 100
            reference_humidity_max[reference_humidity_max < 0] <- 0.01
            if(scenario != ""){
              reference_humidity_min <- reference_humidity_min + RH_diff
            }
            if(scenario != ""){
              reference_humidity_max <- reference_humidity_max + RH_diff
            }
          }#end check for year is 1971 or later
        }
        if(adiab_cor == 1){
          reference_humidity_max.orig <- reference_humidity_max
          reference_humidity_min.orig <- reference_humidity_min
          # correct for potential change in RH with elevation-corrected Tair
          es <- WETAIR(db = reference_temperature_max, rh = 100)$esat
          e <- WETAIR(db = reference_temperature_max.orig, rh = reference_humidity_min.orig)$e
          reference_humidity_min <- (e / es) * 100
          reference_humidity_min[reference_humidity_min > 100] <- 100
          reference_humidity_min[reference_humidity_min < 0] <- 0.01
          es <- WETAIR(db = reference_temperature_min, rh = 100)$esat
          e <- WETAIR(db = reference_temperature_min.orig, rh = reference_humidity_max.orig)$e
          reference_humidity_max <- (e / es) * 100
          reference_humidity_max[reference_humidity_max > 100] <- 100
          reference_humidity_max[reference_humidity_max < 0] <- 0.01
        }
        # AUSCLIM query statements
        clouds <- paste("select cloud1,cloud2,cloud3,cloud4,cloud5,cloud6,cloud7,cloud8,cloud9,cloud10,cloud11,cloud12 FROM cloudcover WHERE i = ",dbrow,sep="")
        maxwinds <- paste("select maxwind1,maxwind2,maxwind3,maxwind4,maxwind5,maxwind6,maxwind7,maxwind8,maxwind9,maxwind10,maxwind11,maxwind12 FROM maxwind WHERE i = ",dbrow,sep="")
        minwinds <- paste("select minwind1,minwind2,minwind3,minwind4,minwind5,minwind6,minwind7,minwind8,minwind9,minwind10,minwind11,minwind12 FROM minwind WHERE i = ",dbrow,sep="")
        maxhumidities <- paste("select maxhum1,maxhum2,maxhum3,maxhum4,maxhum5,maxhum6,maxhum7,maxhum8,maxhum9,maxhum10,maxhum11,maxhum12 FROM maxhum WHERE i = ",dbrow,sep="")
        minhumidities <- paste("select minhum1,minhum2,minhum3,minhum4,minhum5,minhum6,minhum7,minhum8,minhum9,minhum10,minhum11,minhum12 FROM minhum WHERE i = ",dbrow,sep="")
        rainfall <- paste("select rainfall1,rainfall2,rainfall3,rainfall4,rainfall5,rainfall6,rainfall7,rainfall8,rainfall9,rainfall10,rainfall11,rainfall12 FROM rainfall WHERE i = ",dbrow,sep="")
        rainydays <- paste("select rainy1,rainy2,rainy3,rainy4,rainy5,rainy6,rainy7,rainy8,rainy9,rainy10,rainy11,rainy12 FROM rainydays WHERE i = ",dbrow,sep="")
        ALLMINTEMPS <- reference_temperature_min
        ALLMAXTEMPS <- reference_temperature_max
        ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS)
        if(opendap == 0){
          channel2 <- RMySQL::dbConnect(MySQL(), user = uid, password = pwd, host = host, dbname = "ausclim", port = 3306)
          reference_wind_max <- dbGetQuery(channel2, maxwinds)
          reference_wind_min <- dbGetQuery(channel2, minwinds)
          dbDisconnect(channel2)
          if(dailywind != 1 ){
            WNMAXX1 <- suppressWarnings(spline(doys12, reference_wind_max, n = 365, xmin = 1, xmax = 365, method = "periodic"))
            reference_wind_max <- leapfix(WNMAXX1$y, yearlist)
            WNMINN1 <- suppressWarnings(spline(doys12, reference_wind_min, n = 365, xmin = 1, xmax = 365, method = "periodic"))
            reference_wind_min <- leapfix(WNMINN1$y, yearlist)
            if(scenario != ""){
              reference_wind_max <- reference_wind_max * WIND_diff
              reference_wind_min <- reference_wind_min * WIND_diff
            }
          }
        }
        if(soildata == 1){
          surface_emissivity <- suppressWarnings(spline(doys12, surface_emissivity, n = 365, xmin = 1, xmax = 365, method = "periodic"))$y
          surface_emissivity <- leapfix(surface_emissivity, yearlist)
          maxshades1 <- suppressWarnings(spline(doys12, shademax, n = 365, xmin = 1, xmax = 365, method = "periodic"))
          maximum_shade_daily <- leapfix(maxshades1$y * 100, yearlist)
          maximum_shade_daily <- maximum_shade_daily[1:ndays]
        }else{
          if(manualshade == 0){
            maxshades1 <-suppressWarnings(spline(doys12, shademax, n = 365, xmin = 1, xmax = 365, method = "periodic"))
            maximum_shade_daily <- leapfix(maxshades1$y * 100, yearlist)
          }
        }
        albedo <- rep(surface_reflectivity, ndays)
        if((soildata == 1) & (length(rainfall) > 0)){
          soilwet <- rainfall
          soilwet[soilwet <= rainwet] <- 0
          soilwet[soilwet > 0] <- 90
          soil_wetness <- pmax(soilwet, soil_wetness)
        }else{
          albedo <- rep(surface_reflectivity, ndays)
          soil_wetness <- rep(soil_wetness, ndays)
          soilwet <- rainfall
          soilwet[soilwet <= rainwet] <- 0
          soilwet[soilwet > 0] <- 90
          soil_wetness <- pmax(soilwet, soil_wetness)
        }
        Numtyps <- 10 # number of substrate types
        soil_nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
        soil_nodes[1:10, ] <- c(1:10) # deepest nodes for each substrate type
        solar_noon_longitude <- get_timezone_alref(lon = x[1], lat = x[2], timezone = timezone)
        hemisphere <- ifelse(x[2]<0, 2, 1)
        latitude_degrees <- abs(trunc(x[2]))
        latitude_minutes <- (abs(x[2]) - latitude_degrees) * 60
        longitude_degrees <- abs(trunc(x[1]))
        longitude_minutes <- (abs(x[1]) - longitude_degrees) * 60
        elevation <- ALTITUDES
        SLOPE <- SLOPES
        AZMUTH <- AZMUTHS

        avetemp <- (sum(reference_temperature_max) + sum(reference_temperature_min)) / (length(reference_temperature_max) * 2)
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
        mean_annual_temperature <- mean(unlist(ALLTEMPS))

        if(nyears == 1){
          avetemp<-(sum(reference_temperature_max)+sum(reference_temperature_min)) / (length(reference_temperature_max) * 2)
          if(ystart %in% leapyears){
            deep_soil_temperature <- rep(avetemp, 366)
          }else{
            deep_soil_temperature <- rep(avetemp, 365)
          }
        }else{
          if(nrow(reference_temperature_max) == 1){
            avetemp <- rowMeans(t(rbind(reference_temperature_max, reference_temperature_min)), na.rm = TRUE)
          }else{
            avetemp <- rowMeans(cbind(reference_temperature_max, reference_temperature_min), na.rm = TRUE)
          }
          if(length(reference_temperature_max) < 365){
            deep_soil_temperature <- rep((sum(reference_temperature_max) + sum(reference_temperature_min)) / (length(reference_temperature_max) * 2), length(reference_temperature_max))
          }else{
            deep_soil_temperature <- terra::roll(avetemp, n = 365, fun = mean, type = 'to')
            yearone <- rep((sum(reference_temperature_max[1:365]) + sum(reference_temperature_min[1:365])) / (365 * 2), 365)
            deep_soil_temperature[1:365] <- yearone
            # SST
          }
        }
        hourly <- 0
        if(microclima == 1){
          hourly <- 2
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
                                  terrain = 1, horizon.step = 10, nangles = 36,
                                  user_hori = if (!is.na(horizon_angles[1])) horizon_angles else NULL,
                                  slope = slope, aspect = aspect)
          dem <- dem_result$dem; slope <- dem_result$slope; aspect <- dem_result$aspect
          ha36 <- dem_result$horizon_angles
          for (i in 1:length(hour.microclima)) {
            saz <- solazi(hour.microclima[i], lat, long, jd[i], merid = long)
            saz <- round(saz/10, 0) + 1
            saz <- ifelse(saz > 36, 1, saz)
            ha[i] <- ha36[saz]
          }
          #demmeso <- dem
          #info <- .eleveffects(hourlydata, demmeso, lat, long, windthresh = 4.5, emthresh = 0.78)
          #elevation <- info$tout
          cloudhr <- cbind(rep(seq(1, length(cloud)),24), rep(cloud, 24))
          cloudhr <- cloudhr[order(cloudhr[,1]),]
          cloudhr <- cloudhr[,2]
          cloudhr <- leapfix(cloudhr, yearlist, 24)
          dsw2 <- leapfix(clearskyrad[,2], yearlist, 24) *(0.36+0.64*(1-cloudhr/100)) # Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
          # partition total solar into diffuse and direct using code from microclima::hourlyNCEP
          sol <- compute_solar_partition(dsw2, jd, hour.microclima, lat, long,
                                         slope, aspect, ha, microclima.LOR, microclima.leaf_area_index)
          global_radiation <- sol$SOLRhr_all
          diffuse_frac <- sol$diffuse_frac_all
          VIEWF <- 1 # accounted for already in microclima cals
          horizon_angles <- rep(0, 24) # accounted for already in microclima calcs
        }else{
          diffuse_frac <- NA
        }

        if(opendap == 0){
          # correct for fact that wind is measured at 2 m height
          # wind shear equation v / vo = (h / ho)^a
          #where
          #v = the velocity at height h (m/s)
          #vo = the velocity at height ho (m/s)
          #a = the wind shear exponent
          #Terrain   Wind Shear Exponent
          #- a -
          #  Open water   0.1
          #Smooth, level, grass-covered   0.15
          #Row crops   0.2
          #Low bushes with a few trees 	0.2
          #Heavy trees 	0.25
          #Several buildings 	0.25
          #Hilly, mountainous terrain 	0.25
          if(dailywind != 1){
            reference_wind_min <- reference_wind_min * (1.2 / 10) ^ 0.15 * 0.1 # reduce min wind further because have only 9am/3pm values to get max/min
            reference_wind_max <- reference_wind_max * (1.2 / 10) ^ 0.15
            message('min wind * 0.1 \n')
          }else{
            reference_wind_max<-dwind * (1.2 / 2) ^ 0.15
            reference_wind_min<-reference_wind_max
            reference_wind_max<-reference_wind_max * 2
            reference_wind_min<-reference_wind_min * 0.5
            reference_wind_min[reference_wind_min < 0.1] <- 0.1
            message('min wind * 0.5 \n')
            message('max wind * 2 \n')
          }
          cloud_min <- cloud_min * 0.5
          cloud_max <- cloud_max * 2
          cloud_min[cloud_min > 100] <- 100
          cloud_max[cloud_max > 100] <- 100
          if(clearsky == 1){
            cloud_min <- cloud_min * 0
            cloud_max <- cloud_max * 0
            message('running for clear sky conditions')
          }else{
            message('min cloud * 0.5 \n')
            message('max cloud * 2 \n')
          }
        }else{
          if(clearsky == 1){
            cloud_min <- cloud_min * 0
            cloud_max <- cloud_max * 0
            message('running for clear sky conditions')
          }
        }
        # impose uniform warming
        reference_temperature_max <- reference_temperature_max + air_temperature_offset
        reference_temperature_min <- reference_temperature_min + air_temperature_offset
        # impose wind multiplication factor
        reference_wind_max <- reference_wind_max * wind_multiplier
        reference_wind_min <- reference_wind_min * wind_multiplier
        if(soildata != 1){
          surface_emissivity <- matrix(nrow = ndays, data = 0)
          surface_emissivity <- surface_emissivity + surface_emissivity
        }
        #quick fix to make it so that minimum_shade_daily is at the user-specified value and maximum_shade_daily is from the FAPAR database
        if(soildata == 1 & manualshade == 0){
          minimum_shade_daily <- maximum_shade_daily
          minimum_shade_daily[1:length(minimum_shade_daily)] <- minimum_shade
          minimum_shade_daily <- maximum_shade_daily # this is to make one shade level effectively, as dictated by FAPAR
          maximum_shade_daily <- minimum_shade_daily + 0.1
        }

        moists2 <- matrix(nrow = 10, ncol = ndays, data = 0)
        moists2[1, ndays] <- 0.2
        soil_moisture_profile <- moists2

        if(soil_moisture_model == 1){
          moists2 <- matrix(nrow = 10, ncol = ndays, data = 0) # set up an empty vector for soil moisture values through time
          moists2[1:10,] <- initial_soil_moisture
          soil_moisture_profile <- moists2
        }
        soil_properties<-matrix(data = 0, nrow = 10, ncol = 5)

        soil_properties[, 1] <- bulk_density
        soil_properties[,2] <- 1 - bulk_density / mineral_density # not used if soil moisture computed
        soil_properties[soil_properties[,2] < 0.26, 2] <- 0.26
        soil_properties[, 3] <- mineral_conductivity
        soil_properties[, 4] <- mineral_heat_capacity
        soil_properties[, 5] <- mineral_density

        if(organic_soil_cap == 1){
          soil_properties[1:2, 3] <- 0.2
          soil_properties[1:2, 4] <- 1920
        }

        if(organic_soil_cap == 2){
          soil_properties[1:2, 3] <- 0.1
          soil_properties[3:4, 3] <- 0.25
          soil_properties[1:4, 4] <- 1920
          soil_properties[1:4, 5] <- 1.3
          soil_properties[1:4, 1] <- 0.7
        }

        # microclimate input parameters list elevation, solar_noon_longitude, longitude_minutes, longitude_degrees, latitude_minutes, latitude_degrees
        elevation <- as.numeric(elevation)
        solar_noon_longitude <- as.numeric(solar_noon_longitude)
        longitude_minutes <- as.numeric(longitude_minutes)
        longitude_degrees <- as.numeric(longitude_degrees)
        latitude_minutes <- as.numeric(latitude_minutes)
        latitude_degrees <- as.numeric(latitude_degrees)

        # --- 3. Prepare model inputs ---
        micro_input <- build_microinput(ndays, roughness_height, tolerance, local_height, reference_height, Numtyps,
          roughness_height_1, roughness_height_2, wind_profile_height_1, wind_profile_height_2, start_day, end_day, hemisphere, latitude_degrees, latitude_minutes, longitude_degrees, longitude_minutes,
          solar_noon_longitude, slope, azmuth, elevation, precipitable_water, microdaily, mean_annual_temperature, orbital_eccentricity, VIEWF,
          snowtemp, snowdens, snowmelt, undercatch, rain_multiplier, runshade, soil_moisture_model,
          maximum_pooling_depth, evenrain, snow_model, rainmelt, output_to_csv, densfun, hourly,
          rainhourly, radiation_per_wavelength, scattered_uv, root_resistance, stomatal_closure_potential, leaf_resistance, stomatal_stability_parameter, root_radius, moist_error, moist_count, longwave_radiation_model, message,
          fail, snowcond, intercept, grasshade, solar_model_only, canopy_roughness_height, zero_plane_displacement, maxima_times, minima_times,
          spinup, maxsurf, max_iterations_per_day)

        # hourly option set to 0, so make empty vectors
        if(hourly != 1){
          reference_temperature <- rep(0, 24 * ndays)
          reference_humidity <- rep(0, 24 * ndays)
          reference_wind_speed <- rep(0, 24 * ndays)
          cloud_cover <- rep(0, 24 * ndays)
          zenith_angle_hourly <- rep(-1, 24 * ndays)
          longwave_radiation <- rep(-1, 24 * ndays)
        }
        if(hourly == 0){
          global_radiation <- rep(0, 24 * ndays)
        }
        if(rainhourly == 0){
          rainfall_hourly <- rep(0, 24 * ndays)
        }else{
          rainfall_hourly <- rainhour
        }

        if(length(leaf_area_index) <ndays){
          leaf_area_index <- rep(leaf_area_index[1], ndays)
        }
        if(shore == 0){
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
          location <- paste("long", loc[1], "lat", loc[2])
        }else{
          location <- loc
        }
        if(opendap == 0){
         all_cons <- dbListConnections(MySQL())
         for(con in all_cons) +  dbDisconnect(con)
        }
        if(runmicro){
        message(paste('running microclimate model for', ndays, 'days from ', ystart, ' to ', yfinish, ' at site', location, '\n'))
        message('Note: the output column `SOLR` in micromet_lowshade and micromet_highshade is for unshaded horizontal plane solar radiation \n')
        ptm <- proc.time() # Start timing
        microut <- microclimate(micro)
        message(paste0('runtime ', (proc.time() - ptm)[3], ' seconds')) # Stop the clock

        # --- 5. Process and return results ---
        out <- process_micro_output(microut, soil_moisture_model, snow_model, radiation_per_wavelength)
        if(max(out$micromet_lowshade[,1] == 0)){
          message("ERROR: the model crashed - try a different error tolerance (tolerance) or a different spacing in depths")
        }
        dates <- seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'Etc/GMT+10'),
                     as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y ", tz = 'Etc/GMT+10'),
                     by = 'hours')[1:(length(reference_temperature_max) * 24)]
        dates2 <- seq(as.POSIXct(paste0("01/01/", ystart), format = "%d/%m/%Y", tz = 'Etc/GMT+10'),
                      as.POSIXct(paste0("01/01/", yfinish + 1), format = "%d/%m/%Y", tz = 'Etc/GMT+10'),
                      by = 'days')[1:length(reference_temperature_max)]
        return(build_micro_return(out, rainfall, ndays, elevation, surface_reflectivity,
          longlat = c(x[1], x[2]), nyears, timeinterval = ndays, minimum_shade_daily, maximum_shade_daily,
          depths, dates, dates2, air_entry_water_potential, soil_bulk_density, soil_mineral_density, campbell_b_parameter, saturated_hydraulic_conductivity, dem = NA, diffuse_frac,
          snow_model, radiation_per_wavelength))
      }else{
        # Weather-only return (runmicro = FALSE): climate inputs without model run
        return(list(rainfall = rainfall, reference_temperature_max = reference_temperature_max, reference_temperature_min = reference_temperature_min,
                    reference_humidity_max = reference_humidity_max, reference_humidity_min = reference_humidity_min, reference_wind_max = reference_wind_max,
                    reference_wind_min = reference_wind_min, cloud_max = cloud_max, cloud_min = cloud_min,
                    cloud_cover = cloud_cover, reference_wind_speed = reference_wind_speed, reference_temperature = reference_temperature, reference_humidity = reference_humidity,
                    rainfall_hourly = rainfall_hourly, global_radiation = global_radiation, zenith_angle_hourly = zenith_angle_hourly,
                    longwave_radiation = longwave_radiation, dates = dates, dates2 = dates2,
                    air_entry_water_potential = air_entry_water_potential, soil_bulk_density = soil_bulk_density, soil_mineral_density = soil_mineral_density, campbell_b_parameter = campbell_b_parameter, saturated_hydraulic_conductivity = saturated_hydraulic_conductivity))
      }
      } # end of check for na sites
    } # end of check if soil data is being used but no soil data returned
  } # end error trapping
} # end of micro_aust function
