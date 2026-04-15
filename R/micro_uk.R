#' United Kingdom implementation of the microclimate model.
#'
#' An implementation of the NicheMapR microclimate model that uses the CHESS daily weather database. Specifically uses the following files, and needs all 12 for a given year: dtr, tas, huss, precip, rsds, sfcWind, psurf. Also uses the OS Terrain 50 DEM (assumes file name is 'terr50.tif'), which was aggregated to 1km to match met data for lapse-rate correction (assumes file name is 'terr1000.tif'). Citation for CHESS data: Robinson, E.L. ; Blyth, E. ; Clark, D.B. ; Comyn-Platt, E. ; Finch, J. ; Rudd, A.C. (2017). Climate hydrology and ecology research support system meteorology dataset for Great Britain (1961-2015) [CHESS-met] v1.2. NERC Environmental Information Data Centre. https://catalogue.ceh.ac.uk/documents/2ab15bf0-ad08-415c-ba64-831168be7293
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
#' @usage micro_UK(loc = "London, UK", timeinterval = 365, ystart = 2015, yfinish = 2015,
#' surface_reflectivity = 0.15, slope = 0, aspect = 0, depths = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minimum_shade = 0, maximum_shade = 90,
#' local_height = 0.01, ...)
#' @export
#'
#' @details
#' \itemize{
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#'
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{global_aerosol_database}{ = 1, Use the Global Aerosol Database? 1=yes (Fortran version), 2=yes (R version), 0=no}\cr\cr
#' \code{longwave_radiation_model}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1)}\cr\cr
#' \code{solar_model_only}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{radiation_per_wavelength}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{scattered_uv}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{max_iterations_per_day}{ = 3, iterations of first day to get a steady periodic}\cr\cr
#' \code{initial_soil_temperature}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{output_to_csv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{manualshade}{ = 1, Use CSIRO Soil and Landscape Grid of Australia? 1=yes, 0=no}\cr\cr
#' \code{terrain}{ = 0, Use 250m resolution terrain data? 1=yes, 0=no}\cr\cr
#' \code{dailywind}{ = 1, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{wind_multiplier}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{adiab_cor}{ = 1, use adiabatic lapse rate correction? 1=yes, 0=no}\cr\cr
#' \code{air_temperature_offset}{ = 0, warming offset vector, °C (negative values mean cooling). Can supply a single value or a vector the length of the number of days to be simulated.}\cr\cr
#' \code{spatial}{ = "c:/Australian Environment/", choose location of terrain data}\cr\cr
#' \code{opendap}{ = 1, query met grids via opendap (does not work on stomatal_closure_potential unless you compile ncdf4 - see https://github.com/pmjherman/r-ncdf4-build-opendap-windows)}\cr\cr
#' \code{soilgrids}{ = 1, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#' \code{save}{ = 0, don't save forcing data (0), save the forcing data (1) or read previously saved data (2)}\cr\cr
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
#' \code{horizon_angles}{ = rep(0,24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#' \code{minimum_temperature_lapse_rate}{ = 0.0039 Lapse rate for minimum air temperature (degrees C/m)}\cr\cr
#' \code{maximum_temperature_lapse_rate}{ = 0.0077 Lapse rate for maximum air temperature (degrees C/m)}\cr\cr
#' \code{maxima_times}{ = c(1.0, 1.0, 0.0, 0.0), Time of Maximums for Air Wind RelHum Cloud (h), air & Wind max's relative to solar noon, humidity and cloud cover max's relative to sunrise}\cr\cr
#' \code{minima_times}{ = c(0, 0, 1, 1), Time of Minimums for Air Wind RelHum Cloud (h), air & Wind min's relative to sunrise, humidity and cloud cover min's relative to solar noon}\cr\cr
#' \code{timezone}{ = 0, Use GNtimezone function in package geonames to correct to local time zone (excluding daylight saving correction)? 1=yes, 0=no}\cr\cr
#'
#' \strong{ Soil moisture mode parameters:}
#'
#' \code{soil_moisture_model}{ = 1, Run soil moisture model? 1=yes, 0=no  1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{air_entry_water_potential}{ = rep(1.1,19), Air entry potential (J/kg) (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{saturated_hydraulic_conductivity}{ = rep(0.0037,19), Saturated conductivity, (kg s/m3) (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{campbell_b_parameter}{ = rep(4.5,19), Campbell's soil 'b' parameter (-) (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{soil_bulk_density}{ = rep(1.3,19), Soil bulk density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter}
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{soil_mineral_density}{ = rep(2.56,19), Soil density (Mg/m3)  (19 values descending through soil for specified soil nodes in parameter depths and points half way between)}\cr\cr
#' \code{depths}
#' { and points half way between)}\cr\cr
#' \code{maximum_pooling_depth}{ = 10000, Max depth for water pooling on the surface (mm), to account for runoff}\cr\cr
#' \code{rain_multiplier}{ = 1, Rain multiplier for surface soil moisture (-), used to induce runon}\cr\cr
#' \code{rainoff}{ = 0, Rain offset (mm), used to induce changes in rainfall from CHESS values. Can be a single value or a vector matching the number of days to simulate. If negative values are used, rainfall will be prevented from becomming negative.}\cr\cr
#' \code{evenrain}{ = 0, Spread daily rainfall evenly across 24hrs (1) or one event at midnight (0)}\cr\cr
#' \code{initial_soil_moisture}{ = c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3), initial soil water content at each soil node, m3/m3}\cr\cr
#' \code{root_density}{ = c(0,0,8.2,8.0,7.8,7.4,7.1,6.4,5.8,4.8,4.0,1.8,0.9,0.6,0.8,0.4,0.4,0,0)*10000, root density (m/m3), (19 values descending through soil for specified soil nodes in parameter}\cr\cr
#' \code{root_radius}{ = 0.001, root radius, m}\cr\cr
#' \code{root_resistance}{ = 2.5e+10, resistance per unit length of root, m3 kg-1 s-1}\cr\cr
#' \code{leaf_resistance}{ = 2e+6, resistance per unit length of leaf, m3 kg-1 s-1}\cr\cr
#' \code{stomatal_closure_potential}{ = -1500, critical leaf water potential for stomatal closure, J kg-1}\cr\cr
#' \code{stomatal_stability_parameter}{ = 10, stability parameter for stomatal closure equation, -}\cr\cr
#' \code{moist_error}{ = 1e-06, maximum allowable mass balance error, kg}\cr\cr
#' { and points half way between)}\cr\cr
#' \code{leaf_area_index}{ = 0.1, leaf area index (can be a single value or a vector of daily values), used to partition traspiration/evaporation from PET}\cr\cr
#'
#' \strong{ Snow mode parameters:}
#'
#' \code{snow_model}{ = 1, run the snow model 1=yes, 0=no (note that this may cause slower runs)}\cr\cr
#' \code{snowtemp}{ = 1.5, Temperature (°C) at which precipitation falls as snow}\cr\cr
#' \code{snowdens}{ = 0.375, snow density (Mg/m3), overridden by densfun}\cr\cr
#' \code{densfun}{ = c(0.5979, 0.2178, 0.001, 0.0038), slope and intercept of model of snow density as a linear function of snowpack age if first two values are nonzero, and following the exponential function of Sturm et al. 2010 J. of Hydromet. 11:1380-1394 if all values are non-zero; if it is c(0,0,0,0) then fixed density used}\cr\cr
#' \code{snowmelt}{ = 0.9, proportion of calculated snowmelt that doesn't refreeze}\cr\cr
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
#' \code{tides}{ = matrix(data = 0., nrow = 24*timeinterval*nyears, ncol = 3), matrix of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#' }
#'
#' \strong{Outputs:}
#'
#' \code{ndays}{ - number of days for which predictions are made}\cr\cr
#' \code{longlat}{ - longitude and latitude for which simulation was run (decimal degrees)}\cr\cr
#' \code{dates}{ - vector of dates (hourly, POSIXct, timezone = GMT)}\cr\cr
#' \code{dates2}{ - vector of dates (daily, POSIXct, timezone = GMT)}\cr\cr
#' \code{nyears}{ - number of years for which predictions are made}\cr\cr
#' \code{rainfall}{ - vector of daily rainfall (mm)}\cr\cr
#' \code{elevation}{ - elevation at point of simulation (m)}\cr\cr
#' \code{minimum_shade}{ - minimum shade for each day of simulation (\%)}\cr\cr
#' \code{maximum_shade}{ - maximum shade for each day of simulation (\%)}\cr\cr
#' \code{depths}{ - vector of depths used (cm)}\cr\cr
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
#' library(NicheMapR)
#' micro<-micro_uk(runshade = 0) # run the model with default location (London) and settings
#'
#' micromet_lowshade<-as.data.frame(micro$micromet_lowshade) # above ground microclimatic conditions, min shade
#' soil_temperature_lowshade<-as.data.frame(micro$soil_temperature_lowshade) # soil temperatures, minimum shade
#' soil_moisture_lowshade<-as.data.frame(micro$soil_moisture_lowshade) # soil temperatures, minimum shade
#'
#' # append dates
#' ystart <- 2015
#' yfinish <- 2015
#' nyears <- yfinish-ystart+1
#' tzone<-paste("Etc/GMT+",0,sep="")
#' dates<-seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="hours")
#'
#' micromet_lowshade <- cbind(dates,micromet_lowshade)
#' soil_temperature_lowshade <- cbind(dates,soil_temperature_lowshade)
#' soil_moisture_lowshade <- cbind(dates, soil_moisture_lowshade)
#' minimum_shade<-micro$minimum_shade
#'
#' # plotting above-ground conditions in minimum shade
#' with(micromet_lowshade,{plot(air_temperature_local ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l",main=paste("air and sky temperature, ",minimum_shade,"% shade",sep=""), ylim = c(-20, 60))})
#' with(micromet_lowshade,{points(air_temperature_reference ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l",lty=2,col='blue')})
#' with(micromet_lowshade,{points(sky_temperature ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' ,  type = "l",col='light blue',main=paste("sky temperature, ",minimum_shade,"% shade",sep=""))})
#' with(micromet_lowshade,{plot(relative_humidity_local ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l",ylim=c(0,100),main=paste("humidity, ",minimum_shade,"% shade",sep=""))})
#' with(micromet_lowshade,{points(relative_humidity_reference ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l",col='blue',lty=2,ylim=c(0,100))})
#' with(micromet_lowshade,{plot(wind_speed_reference ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#' ,  type = "l",main="wind speed",ylim = c(0, 15))})
#' with(micromet_lowshade,{points(wind_speed_local ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#' ,  type = "l",lty=2,col='blue')})
#' with(micromet_lowshade,{plot(solar_radiation ~ dates,xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
#' ,  type = "l",main="solar radiation")})
#' with(micromet_lowshade,{plot(snow_depth ~ dates,xlab = "Date and Time", ylab = "Snow Depth (cm)"
#' ,  type = "l",main="snow depth")})
#'
#' # plotting soil temperature
#' for(i in 1:10){
#'  if(i==1){
#'    plot(soil_temperature_lowshade[,i+3]~soil_temperature_lowshade[,1],xlab = "Date and Time", ylab = "Soil Temperature (°C)"
#'    ,col=i,type = "l",main=paste("soil temperature ",minimum_shade,"% shade",sep=""))
#'  }else{
#'    points(soil_temperature_lowshade[,i+3]~soil_temperature_lowshade[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)",col=i,type = "l")
#'  }
#' }
#'
#' # plotting soil moisture
#' for(i in 1:10){
#'  if(i==1){
#'    plot(soil_moisture_lowshade[,i+3]~soil_moisture_lowshade[,1],xlab = "Date and Time", ylab = "Soil Moisture (% volumetric)"
#'    ,col=i,type = "l",main=paste("soil moisture ",minimum_shade,"% shade",sep=""))
#'  }else{
#'    points(soil_moisture_lowshade[,i+3]~soil_moisture_lowshade[,1],xlab = "Date and Time", ylab = "Soil Moisture
#'     (%)",col=i,type = "l")
#'  }
#' }
micro_uk <- function(
  loc = c(-0.1277583, 51.50735),
  timeinterval = 365,
  ystart = 2015,
  yfinish = 2015,
  nyears = yfinish - ystart + 1,
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
  terrain = 0,
  dailywind = 1,
  wind_multiplier = 1,
  adiab_cor = 1,
  air_temperature_offset = 0,
  spatial = "Q:/CHESS",
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
  horizon_angles = rep(0,24),
  maxima_times=c(1, 1, 0, 0),
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
  snow_model = 1,
  snowtemp = 1.5,
  snowdens = 0.375,
  densfun = c(0.5979, 0.2178, 0.001, 0.0038),
  snowmelt = 1,
  undercatch = 1,
  rainmelt = 0.0125,
  shore = 0,
  tides = matrix(data = 0, nrow = 24 * timeinterval * nyears, ncol = 3),
  hourly = 0,
  rainhourly = 0,
  rainhour = 0,
  rainoff = 0,
  radiation_per_wavelength = 0,
  scattered_uv = 0,
  max_iterations_per_day = 3,
  opendap = 0,
  soilgrids = 0,
  longwave_radiation_model = 0,
  message = 0,
  fail = nyears * 24 * 365,
  runmicro = 1,
  snowcond = 0,
  intercept = max(maximum_shade) / 100 * 0.3,
  save = 0,
  grasshade = 0,
  maxsurf = 85) {

  # error trapping - originally inside the Fortran code, but now checking before executing Fortran
  # --- 1. Input validation ---
  reference_height <- 1.2  # met data reference height (m) for CHESS data
  errors <- 0
  # uk-specific checks
  if(timeinterval < 12 | timeinterval > 365){
    message("ERROR: the variable 'timeinterval' is out of bounds. Please enter a correct value (12 - 365).")
    errors <- 1
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

  if(errors==0){ # continue

    ################## time related variables #################################
    nyears<-yfinish-ystart+1

    doys12<-c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month

    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial max_iterations_per_day day burn in

    daystart<-as.integer(ceiling(365/timeinterval/2))
    start_day <- 1 # start day
    dates<-Sys.time()-60*60*24
    curyear<-as.numeric(format(dates,"%Y"))

    ################## location and terrain #################################
    if (!requireNamespace("terra", quietly = TRUE)) {
      stop("package 'terra' is needed. Please install it.",
        call. = FALSE)
    }
    if (!requireNamespace("ncdf4", quietly = TRUE)) {
      stop("package 'ncdf4' is needed. Please install it.",
        call. = FALSE)
    }
    if (!requireNamespace("sp", quietly = TRUE)) {
      stop("package 'sp' is needed. Please install it.",
        call. = FALSE)
    }
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    requireNamespace("terra")
    requireNamespace("ncdf4")
    # reference meridian for solar noon calculation
    solar_noon_longitude <- get_timezone_alref(lon = x[1], lat = x[2], timezone = timezone)
    hemisphere <- ifelse(x[2]<0,2.,1.) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    latitude_degrees <- abs(trunc(x[2]))
    latitude_minutes <- (abs(x[2])-latitude_degrees)*60
    longitude_degrees <- abs(trunc(x[1]))
    longitude_minutes <- (abs(x[1])-longitude_degrees)*60
    azmuth<-aspect

    ConvertCoordinates <- function(easting,northing) {
      out = cbind(easting,northing)
      mask = !is.na(easting)
      sp <-  sp::spTransform(sp::SpatialPoints(list(easting[mask],northing[mask]),proj4string=sp::CRS(bng)),sp::CRS(wgs84))
      out[mask,]=sp@coords
      out
    }
    # conversion parameters for UTM to lat/lon
    wgs84 = "+init=epsg:4326"
    bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000
+ellps=airy +datum=OSGB36 +units=m +no_defs'
    sp <-  sp::spTransform(sp::SpatialPoints(x,proj4string=sp::CRS(wgs84)),sp::CRS(bng))
    x2 <- sp@coords
    UKDEM <-as.numeric(terra::extract(terra::rast(paste0(spatial,"/terr1000.tif")), x2)[, 2])
    ALTITUDES <- as.numeric(terra::extract(terra::rast(paste0(spatial,"/terr50.tif")), x2)[, 2])
    if(is.na(ALTITUDES)==TRUE){ALTITUDES<-UKDEM}
    if(is.na(elevation) == FALSE){ # check if user-specified elevation
      ALTITUDES <- elevation
    }
    if(terrain==1){
      cat("extracting terrain data \n")

      # now extract terrain data from elevslpasphori.nc
      # get UTM from dec degrees, NZTM
      HORIZONS <- elevslpasphori[4:27]
      SLOPES <- elevslpasphori[2]
      AZMUTHS <- elevslpasphori[3]
      # the horizons have been arranged so that they go from 0 degrees azimuth (north) clockwise - r.horizon starts
      # in the east and goes counter clockwise!
      HORIZONS <- (ifelse(is.na(HORIZONS),0,HORIZONS))/10 # get rid of na and get back to floating point
      HORIZONS <- data.frame(HORIZONS)
      VIEWF_all <- 1-sum(sin(as.data.frame(horizon_angles)*pi/180))/length(horizon_angles) # convert horizon angles to radians and calc view factor(s)
    }else{
      HORIZONS <- horizon_angles
      HORIZONS <- data.frame(HORIZONS)
      VIEWF_all <- 1-sum(sin(as.data.frame(horizon_angles)*pi/180))/length(horizon_angles) # convert horizon angles to radians and calc view factor(s)
      SLOPES<-rep(slope,length(x[,1]))
      AZMUTHS<-rep(aspect,length(x[,1]))
    }
    horizon_angles<-HORIZONS
    row.names(horizon_angles)<-NULL
    horizon_angles<-as.numeric(as.matrix(horizon_angles))
    VIEWF<-VIEWF_all

    if(save != 2){
      if(soilgrids == 1){
        sg <- fetch_soilgrids(x, depths)
        if (!is.null(sg)) { air_entry_water_potential <- sg$air_entry_water_potential; saturated_hydraulic_conductivity <- sg$saturated_hydraulic_conductivity; campbell_b_parameter <- sg$campbell_b_parameter; soil_bulk_density <- sg$soil_bulk_density; bulk_density <- sg$bulk_density }
      }
    }else{
      if(soilgrids == 1){
        cat("loading SoilGrids data from previous run \n")
        load('air_entry_water_potential.Rda')
        load('campbell_b_parameter.Rda')
        load('soil_bulk_density.Rda')
        load('saturated_hydraulic_conductivity.Rda')
        load('bulk_density.Rda')
      }
    }
    if(save == 1 & soilgrids == 1){
      cat("saving SoilGrids data for later \n")
      save(air_entry_water_potential, file = 'air_entry_water_potential.Rda')
      save(campbell_b_parameter, file = 'campbell_b_parameter.Rda')
      save(soil_bulk_density, file = 'soil_bulk_density.Rda')
      save(saturated_hydraulic_conductivity, file = 'saturated_hydraulic_conductivity.Rda')
      save(bulk_density, file = 'bulk_density.Rda')
    }
    delta_elev = UKDEM - ALTITUDES
    adiab_corr_max <- delta_elev * maximum_temperature_lapse_rate
    adiab_corr_min <- delta_elev * minimum_temperature_lapse_rate
    cat("extracting weather data \n")
    yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
    wgs84 = "+init=epsg:4326"
    bng = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000\n+ellps=airy +datum=OSGB36 +units=m +no_defs"

    month <- c("01", "02", "03", "04", "05", "06", "07",
      "08", "09", "10", "11", "12")

    if(save != 2){
      if(opendap != 0){
        if(!require(futile.logger, quietly = TRUE)){
         stop('package "futile.logger" is required for this opendap extraction process. Please install it.')
        }
        require(utils)
        require(futile.logger)
        retry <- function(expr, isError=function(x) "try-error" %in% class(x), maxErrors=20, sleep=0) {
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

        baseurl <- opendap
        nc <- ncdf4::nc_open(paste0(baseurl, "dtrDetail02/chess_dtr_", yearlist[1],
          "12.nc", sep = ""))
        northings <- retry(matrix(ncdf4::ncvar_get(nc, "y")))
        eastings <- retry(matrix(ncdf4::ncvar_get(nc, "x")))
        ncdf4::nc_close(nc)
        easting <- eastings
        northing <- rep(northings[1], length(easting))
        lon <- ConvertCoordinates(easting, northing)[, 1]
        northing <- northings
        easting <- rep(eastings[1], length(northing))
        lat <- ConvertCoordinates(easting, northing)[, 2]
        flat = match(abs(lat - x[2]) < 0.0013,
          1)
        latindex = which(flat %in% 1)
        flon = match(abs(lon - x[1]) < 0.0013,
          1)
        lonindex = which(flon %in% 1)
        start <- c(lonindex, latindex, 1)
        count <- c(1, 1, -1)
        for (j in 1:nyears) {
          for (jj in 1:12) {
            cat(paste("reading weather input for ", yearlist[j],
              " month ", month[jj], " \n", sep = ""))
            nc <- ncdf4::nc_open(paste0(baseurl, "dtrDetail02/chess_dtr_", yearlist[j], month[jj], ".nc"))
            dtr <- retry(as.numeric(ncdf4::ncvar_get(nc, varid = "dtr",start = start, count)))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste0(baseurl, "tasDetail02/chess_tas_", yearlist[j], month[jj], ".nc"))
            tas <- retry(as.numeric(ncdf4::ncvar_get(nc, varid = "tas", start = start, count)))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste0(baseurl, "psurfDetail02/chess_psurf_", yearlist[j], month[jj], ".nc"))
            psurf <- retry(as.numeric(ncdf4::ncvar_get(nc, varid = "psurf", start = start, count)))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste0(baseurl, "hussDetail02/chess_huss_", yearlist[j], month[jj], ".nc"))
            huss <- retry(as.numeric(ncdf4::ncvar_get(nc, varid = "huss",  start = start, count)))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste0(baseurl, "precipDetail02/chess_precip_", yearlist[j], month[jj], ".nc"))
            precip <- retry(as.numeric(ncdf4::ncvar_get(nc, varid = "precip", start = start, count)))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste0(baseurl, "rsdsDetail02/chess_rsds_", yearlist[j], month[jj], ".nc"))
            rsds <- retry(as.numeric(ncdf4::ncvar_get(nc, varid = "rsds", start = start, count)))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste0(baseurl, "sfcWindDetail02/chess_sfcWind_", yearlist[j], month[jj], ".nc"))
            sfcWind <- retry(as.numeric(ncdf4::ncvar_get(nc, varid = "sfcWind", start = start, count)))
            ncdf4::nc_close(nc)
            if (j == 1 & jj == 1) {
              Tmax <- tas - 273.15 + dtr/2
              Tmin <- tas - 273.15 - dtr/2
              Wind <- sfcWind
              Hum <- huss
              Rain <- precip * 3600 * 24
              solar <- rsds
              press <- psurf
            }
            else {
              Tmax <- c(Tmax, tas - 273.15 + dtr/2)
              Tmin <- c(Tmin, tas - 273.15 - dtr/2)
              Wind <- c(Wind, sfcWind)
              Hum <- c(Hum, huss)
              Rain <- c(Rain, precip * 3600 * 24)
              solar <- c(solar, rsds)
              press <- c(press, psurf)
            }
          }
        }
      }else{
        nc <- ncdf4::nc_open(paste(spatial, "/chess_dtr_", yearlist[1],
          "12.nc", sep = ""))
        northings <- matrix(ncdf4::ncvar_get(nc, "y"))
        eastings <- matrix(ncdf4::ncvar_get(nc, "x"))
        ncdf4::nc_close(nc)
        easting <- eastings
        northing <- rep(northings[1], length(easting))
        lon <- ConvertCoordinates(easting, northing)[, 1]
        northing <- northings
        easting <- rep(eastings[1], length(northing))
        lat <- ConvertCoordinates(easting, northing)[, 2]
        lon_1 <- as.numeric(longlat[1])
        lat_1 <- as.numeric(longlat[2])
        dist1 <- abs(lon - lon_1)
        index1 <- which.min(dist1)
        dist2 <- abs(lat - lat_1)
        index2 <- which.min(dist2)
        start <- c(index1, index2, 1)
        count <- c(1, 1, -1)
        for (j in 1:nyears) {
          for (jj in 1:12) {
            cat(paste("reading weather input for ", yearlist[j],
              " month ", month[jj], " \n", sep = ""))
            nc <- ncdf4::nc_open(paste(spatial, "/chess_dtr_", yearlist[j],
              month[jj], ".nc", sep = ""))
            dtr <- as.numeric(ncdf4::ncvar_get(nc, varid = "dtr",
              start = start, count))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste(spatial, "/chess_tas_", yearlist[j],
              month[jj], ".nc", sep = ""))
            tas <- as.numeric(ncdf4::ncvar_get(nc, varid = "tas",
              start = start, count))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste(spatial, "/chess_psurf_",
              yearlist[j], month[jj], ".nc", sep = ""))
            psurf <- as.numeric(ncdf4::ncvar_get(nc, varid = "psurf",
              start = start, count))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste(spatial, "/chess_huss_",
              yearlist[j], month[jj], ".nc", sep = ""))
            huss <- as.numeric(ncdf4::ncvar_get(nc, varid = "huss",
              start = start, count))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste(spatial, "/chess_precip_",
              yearlist[j], month[jj], ".nc", sep = ""))
            precip <- as.numeric(ncdf4::ncvar_get(nc, varid = "precip",
              start = start, count))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste(spatial, "/chess_rsds_",
              yearlist[j], month[jj], ".nc", sep = ""))
            rsds <- as.numeric(ncdf4::ncvar_get(nc, varid = "rsds",
              start = start, count))
            ncdf4::nc_close(nc)
            nc <- ncdf4::nc_open(paste(spatial, "/chess_sfcWind_",
              yearlist[j], month[jj], ".nc", sep = ""))
            sfcWind <- as.numeric(ncdf4::ncvar_get(nc, varid = "sfcWind",
              start = start, count))
            ncdf4::nc_close(nc)
            if (j == 1 & jj == 1) {
              Tmax <- tas - 273.15 + dtr/2
              Tmin <- tas - 273.15 - dtr/2
              Wind <- sfcWind
              Hum <- huss
              Rain <- precip * 3600 * 24
              solar <- rsds
              press <- psurf
            }
            else {
              Tmax <- c(Tmax, tas - 273.15 + dtr/2)
              Tmin <- c(Tmin, tas - 273.15 - dtr/2)
              Wind <- c(Wind, sfcWind)
              Hum <- c(Hum, huss)
              Rain <- c(Rain, precip * 3600 * 24)
              solar <- c(solar, rsds)
              press <- c(press, psurf)
            }
          }
        }
      }


      # compute clear sky solar for the site of interest, for cloud cover computation below
      cat("running micro_global to get clear sky solar \n")
      micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, timeinterval = 365)
      clearskyrad <- micro_clearsky$micromet_lowshade[,c(1, 13)]
      clearsky_mean1 <- aggregate(clearskyrad[,2], by = list(clearskyrad[,1]), FUN = mean)[,2]
      leapyears<-seq(1972,2060,4)
      for(j in 1:nyears){
        if(yearlist[j]%in%leapyears){# add day for leap year if needed
          clearsky_mean<-c(clearsky_mean1[1:59],clearsky_mean1[59],clearsky_mean1[60:365])
        }else{
          clearsky_mean <- clearsky_mean1
        }
        if(j==1){
          allclearsky <- clearsky_mean
        }else{
          allclearsky <- c(allclearsky, clearsky_mean)
        }
      }
      cloud <- (1-solar/allclearsky) * 100
      cloud[cloud<0]<-0
      cloud[cloud>100]<-100
      if(clearsky == 1){
        cloud <- cloud * 0
      }
      cloud_max<-as.numeric(cloud)
      cloud_min<-cloud_max
      cloud_min<-cloud_min*0.5
      cloud_max<-cloud_max*2
      cloud_min[cloud_min>100]<-100
      cloud_max[cloud_max>100]<-100
      Wind<-Wind * wind_multiplier
      Wind[Wind==0]<-0.1
      if(save == 1){
        cat("saving met data for later \n")
        save(cloud_max, file = 'cloud_max.Rda')
        save(cloud_min, file = 'cloud_min.Rda')
        save(Wind, file = 'Wind.Rda')
        save(Tmax, file = 'Tmax.Rda')
        save(Tmin, file = 'Tmin.Rda')
        save(Hum, file = 'Hum.Rda')
        save(Rain, file = 'Rain.Rda')
        save(press, file = 'press.Rda')
      }
    }else{
      cat("loading met data from previous run \n")
      load('cloud_max.Rda')
      load('cloud_min.Rda')
      load('Wind.Rda')
      load('Tmax.Rda')
      load('Tmin.Rda')
      load('hum.Rda')
      load('Rain.Rda')
      load('press.Rda')
    }
    ndays<-length(Tmax)
    doynum<-ndays
    doys<-seq(daystart,doynum,1)
    day_of_year <- subset(doys, doys!=0)
    end_day<-ndays
    start_day <- 1

    tzone<-paste("Etc/GMT-12",sep="") # doing it this way ignores daylight savings!
    ndays<-length(seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days"))
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

    if(is.na(ALTITUDES)!=TRUE){

      aerosol_optical_depth <- compute_tai(longlat, global_aerosol_database, TAI_AUSTRALIA)

      if(adiab_cor==1){
        reference_temperature_max<-as.matrix(Tmax+adiab_corr_max)
        reference_temperature_min<-as.matrix(Tmin+adiab_corr_min)
      }else{
        reference_temperature_max<-as.matrix(Tmax)
        reference_temperature_min<-as.matrix(Tmin)
      }
      if(air_temperature_offset != 0){
        # impose temperature change
        reference_temperature_max<-reference_temperature_max+air_temperature_offset
        reference_temperature_min<-reference_temperature_min+air_temperature_offset
      }
      rainfall<-Rain+rainoff
      rainfall[rainfall < 0] <- 0

      # convert from specific humidity to relative humidity
      es <- WETAIR(db = reference_temperature_max, rh = 100)$esat
      e <- Hum * press / (0.378 * Hum + 0.622)
      reference_humidity_min <- (e / es) * 100
      reference_humidity_min[reference_humidity_min>100]<-100
      reference_humidity_min[reference_humidity_min<0]<-0.01
      es <- WETAIR(db = reference_temperature_min, rh = 100)$esat
      e <- Hum * press / (0.378 * Hum + 0.622)
      reference_humidity_max <- (e / es) * 100
      reference_humidity_max[reference_humidity_max>100]<-100
      reference_humidity_max[reference_humidity_max<0]<-0.01

      ALLMINTEMPS<-reference_temperature_min
      ALLMAXTEMPS<-reference_temperature_max
      ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)

      reference_wind_max <- Wind * 2
      reference_wind_min <- Wind * 0.5
      message('min wind * 0.5 \n')
      message('max wind * 2 \n')

      albedo <- rep(surface_reflectivity, ndays)
      soil_wetness <- rep(soil_wetness, ndays)
      soilwet<-rainfall
      soilwet[soilwet<=rainwet] = 0
      soilwet[soilwet>0] = 90
      soil_wetness<-pmax(soilwet,soil_wetness)


      Intrvls<-rep(0, ndays)
      Intrvls[1] <- 1 # user-supplied last day-of-year in each time interval sequence
      Numtyps <- 1 # number of substrate types
      Numint <- 1  # number of time intervals
      soil_nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
      soil_nodes[1,1] <- 10. # deepest nodes for each substrate type
      solar_noon_longitude <- get_timezone_alref(lon = x[1], lat = x[2], timezone = timezone)

      hemisphere <- ifelse(x[2]<0, 2, 1)
      latitude_degrees <- abs(trunc(x[2]))
      latitude_minutes <- (abs(x[2])-latitude_degrees)*60
      longitude_degrees <- abs(trunc(x[1]))
      longitude_minutes <- (abs(x[1])-longitude_degrees)*60
      if(adiab_cor==1){
        elevation<-ALTITUDES
      }else{
        elevation<-UKDEM
      }
      SLOPE<-SLOPES
      AZMUTH<-AZMUTHS

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

      if(nyears==1){
        avetemp<-(sum(reference_temperature_max)+sum(reference_temperature_min))/(length(reference_temperature_max)*2)
        deep_soil_temperature<-rep(avetemp,ndays)
      }else{
        if(nrow(reference_temperature_max)==1){
          avetemp<-colMeans(cbind(reference_temperature_max, reference_temperature_min), na.rm=TRUE)
        }else{
          avetemp<-rowMeans(cbind(reference_temperature_max, reference_temperature_min), na.rm=TRUE)
        }
        if(length(reference_temperature_max)<365){
          deep_soil_temperature<-rep((sum(reference_temperature_max)+sum(reference_temperature_min))/(length(reference_temperature_max)*2),length(reference_temperature_max))
        }else{
          deep_soil_temperature<-terra::roll(avetemp,n=365,fun=mean,type='to')
          yearone<-rep((sum(reference_temperature_max[1:365])+sum(reference_temperature_min[1:365]))/(365*2),365)
          deep_soil_temperature[1:365]<-yearone
        }
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
      #Smooth, level, grass-covered   0.15
      #Row crops 	0.2
      #Low bushes with a few trees 	0.2
      #Heavy trees 	0.25
      #Several buildings 	0.25
      #Hilly, mountainous terrain 	0.25
      reference_wind_max<-reference_wind_max*(1.2/10)^0.15
      reference_wind_min<-reference_wind_min*(1.2/10)^0.15

      surface_emissivity<-matrix(nrow=ndays,data=0)
      surface_emissivity<-surface_emissivity+surface_emissivity

      moists2<-matrix(nrow=10, ncol = ndays, data=0)
      moists2[1,ndays]<-0.2
      soil_moisture_profile<-moists2

      if(soil_moisture_model==1){
        if(timeinterval==365){
          moists2<-matrix(nrow=10, ncol = ndays, data=0) # set up an empty vector for soil moisture values through time
        }else{
          moists2<-matrix(nrow=10, ncol = timeinterval, data=0) # set up an empty vector for soil moisture values through time
        }
        moists2[1:10,]<-initial_soil_moisture
        soil_moisture_profile<-moists2
      }
      soil_properties<-matrix(data = 0, nrow = 10, ncol = 5)

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

      elevation<-as.numeric(elevation)
      solar_noon_longitude<-as.numeric(solar_noon_longitude)
      longitude_minutes<-as.numeric(longitude_minutes)
      longitude_degrees<-as.numeric(longitude_degrees)
      latitude_minutes<-as.numeric(latitude_minutes)
      latitude_degrees<-as.numeric(latitude_degrees)

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
      if(hourly==0){
        reference_temperature=rep(0,24*ndays)
        reference_humidity=rep(0,24*ndays)
        reference_wind_speed=rep(0,24*ndays)
        cloud_cover=rep(0,24*ndays)
        global_radiation=rep(0,24*ndays)
        zenith_angle_hourly=rep(-1,24*ndays)
        longwave_radiation=rep(-1,24*ndays)
      }
      if(rainhourly==0){
        rainfall_hourly=rep(0,24*ndays)
      }else{
        rainfall_hourly = rainhour
      }

      if(length(leaf_area_index)<ndays){
        leaf_area_index<-rep(leaf_area_index[1],ndays)
      }
      if(shore==0){
        tides<-matrix(data = 0, nrow = 24*ndays, ncol = 3) # make an empty matrix
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
        location<-paste("long",loc[1],"lat",loc[2])
      }else{
        location<-loc
      }
      if(runmicro){
        cat(paste('running microclimate model for',timeinterval,'days by',nyears,'years at site',location,'\n'))
      message('Note: the output column `SOLR` in micromet_lowshade and micromet_highshade is for unshaded horizontal plane solar radiation \n')
      ptm <- proc.time() # Start timing
      microut<-microclimate(micro)
      print(proc.time() - ptm) # Stop the clock

      # --- 5. Process and return results ---
      out <- process_micro_output(microut, soil_moisture_model, snow_model, radiation_per_wavelength)
      if(max(out$micromet_lowshade[,1] == 0)){
        message("ERROR: the model crashed - try a different error tolerance (tolerance) or a different spacing in depths")
      }
      dates <- seq(as.POSIXct(paste0("01/01/",ystart), format = "%d/%m/%Y", tz = 'GMT'),
                   as.POSIXct(paste0("31/12/",yfinish), format = "%d/%m/%Y", tz = 'GMT'),
                   by = 'hours')
      dates2 <- seq(as.POSIXct(paste0("01/01/",ystart), format = "%d/%m/%Y", tz = 'GMT'),
                    as.POSIXct(paste0("31/12/",yfinish), format = "%d/%m/%Y", tz = 'GMT'),
                    by = 'days')
      return(build_micro_return(out, rainfall, ndays, elevation, surface_reflectivity,
        longlat = c(x[1], x[2]), nyears, timeinterval, minimum_shade_daily, maximum_shade_daily,
        depths, dates, dates2, air_entry_water_potential, soil_bulk_density, soil_mineral_density, campbell_b_parameter, saturated_hydraulic_conductivity, dem = NA, diffuse_frac = NA,
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
  } # end error trapping
} # end of micro_UK function
