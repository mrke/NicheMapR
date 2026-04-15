#' NCEP implementation of the microclimate model, with package microclima computing hourly forcing.
#'
#' An implementation of the NicheMapR microclimate model that integrates the RNCEP 6-hourly weather data and the elevatr package for obtaining DEM using downscaling functions from the microclima package as described in more detail in Kearney, M. R., Gillingham, P. K., Bramer, I., Duffy, J. P., & Maclean, I. M. D. (2019). A method for computing hourly, historical, terrain-corrected microclimate anywhere on Earth. Methods in Ecology and Evolution.
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param dstart First day to run, date in format "d/m/Y" e.g. "01/01/2016"
#' @param dfinish Last day to run, date in format "d/m/Y" e.g. "31/12/2016"
#' @param dem A digital elevation model used by microclima for micro-topographic effects, produced by microclima function 'get_dem' via R package 'elevatr' (internally generated via same function based on 'loc' if NA)
#' @param dem2 A digital elevation model used by microclima for meso-climate calculations, produced by microclima function 'get_dem' via R package 'elevatr' (internally generated via same function based on 'loc' if NA)
#' @param dem_resolution Requested resolution of the DEM from elevatr, m
#' @param pixels Number of pixels along one edge of square requested of DEM requested from elevatr, #
#' @param surface_reflectivity Soil solar reflectance, decimal \%
#' @param slope Slope in degrees (if NA, then derived from DEM with package microclima)
#' @param aspect Aspect in degrees (0 = north) (if NA, then derived from DEM with microclima)
#' @param depths Soil depths at which calculations are to be made (cm), must be 10 values starting from 0, and more closely spaced near the surface
#' @param minimum_shade Minimum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param maximum_shade Maximum shade level to use (\%) (can be a single value or a vector of daily values)
#' @param local_height Local height (m) at which air temperature, wind speed and humidity are to be computed for organism of interest
#' @param coastal Compute coastal effects with microclima? T (TRUE) or F (FALSE) (can take a while and may have high memory requirements depending on DEM size)
#' @param hourlydata user input of the hourlydata matrix
#' @param dailyprecip user input of daily rainfall
#' @param weather.elevation optional value indicating the elevation of values in `hourlydata`. Either a numeric value, corresponding to the elevation in (m) of the location from which `hourlydata` were obtained, or `ncep` (default, data derive from NOAA-NCEP reanalysis) project.
#' @param cad.effects optional logical indicating whether to calculate cold air drainage effects (TRUE = Yes, slower. FALSE =  No, quicker)
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
#' @usage micro_ncep(loc = c(-91.415669, -0.287145), dstart = "01/01/2019", dfinish = "31/07/2019",
#' surface_reflectivity = 0.15, slope = 0, aspect = 0, depths = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minimum_shade = 0, maximum_shade = 90,
#' local_height = 0.01, ...)
#' @export
#' @details
#' \strong{ Parameters controlling how the model runs:}\cr\cr
#' \code{runshade}{ = 1, Run the microclimate model twice, once for each shade level (1) or just once for the minimum shade (0)?}\cr\cr
#' \code{clearsky}{ = 0, Run for clear skies (1) or with observed cloud cover (0)}\cr\cr
#' \code{global_aerosol_database}{ = 1, Use the Global Aerosol Database? 1=yes (Fortran version), 2=yes (R version), 0=no}\cr\cr
#' \code{longwave_radiation_model}{ = 0, Clear-sky longwave radiation computed using Campbell and Norman (1998) eq. 10.10 (includes humidity) (0) or Swinbank formula (1) or from NCEP data (2)}\cr\cr
#' \code{solar_model_only}{ = 0, Only run SOLRAD to get solar radiation? 1=yes, 0=no}\cr\cr
#' \code{radiation_per_wavelength}{ = 0, Return wavelength-specific solar radiation output?}\cr\cr
#' \code{scattered_uv}{ = 0, Use gamma function for scattered solar radiation? (computationally intensive)}\cr\cr
#' \code{max_iterations_per_day}{ = 3, iterations of first day to get a steady periodic}\cr\cr
#' \code{initial_soil_temperature}{ = NA, initial soil temperature at each soil node, °C (if NA, will use the mean air temperature to initialise)}\cr\cr
#' \code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr\cr
#' \code{output_to_csv}{ = 0, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{reanalysis}{ = TRUE, Use reanalysis2 NCEP data? TRUE/FALSE}\cr\cr
#' \code{wind_multiplier}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{air_temperature_offset}{ = 0, warming offset vector, °C (negative values mean cooling). Can supply a single value or a vector the length of the number of days to be simulated.}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#' \code{spatial}{ = NA, specify folder with local NCEP data (no trailing forward slash), goes to the web if NA }\cr\cr
#' \code{scenario}{ = 0, climate change scenario from terraclimate (0 is historical, 2 is plus 2 deg C, 4 is plus 4 deg C }\cr\cr
#' \code{terra_source}{ = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data", specify location of terraclimate data, goes to the web by default}\cr\cr
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
#' \code{horizon_angles}{ = rep(NA,24), Horizon angles (degrees), from 0 degrees azimuth (north) clockwise in 15 degree intervals}\cr\cr
#'
#' \strong{ Soil moisture mode parameters:}\cr\cr
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
#' \code{rainhourly}{ = 0, Is hourly rain input being supplied (1 = yes, 0 = no)?}\cr\cr
#' \code{rainhour}{ = 0, Vector of hourly rainfall values - overrides daily NCEP rain if rainhourly = 1}\cr\cr
#' \code{rain_multiplier}{ = 1, Rain multiplier for surface soil moisture (-), used to induce runon}\cr\cr
#' \code{rainoff}{ = 0, Rain offset (mm), used to induce changes in rainfall from NCEP values. Can be a single value or a vector matching the number of days to simulate. If negative values are used, rainfall will be prevented from becomming negative.}\cr\cr
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
#' \code{leaf_area_index}{ = 0.1, leaf area index (can be a single value or a vector of daily values), used to partition traspiration/evaporation from PET in soil moisture model}\cr\cr
#' \code{microclima.leaf_area_index}{ = 0, leaf area index, used by package microclima for radiation calcs}\cr\cr
#' \code{LOR}{ = 1, leaf orientation for package microclima radiation calcs}\cr\cr
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
#' \code{snowcond}{ = 0, effective snow thermal conductivity W/mC (if zero, uses inbuilt function of density)}\cr\cr
#' \code{intercept}{ = max(maximum_shade) / 100 * 0.3, snow interception fraction for when there's shade (0-1)}\cr\cr
#' \code{grasshade}{ = 0, if 1, means shade is removed when snow is present, because shade is cast by grass/low shrubs}\cr\cr
#'
#' \strong{ Intertidal mode parameters:}
#'
#' \code{shore}{ Include tide effects? If 1, the matrix}
#' \code{tides}{ is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0, nrow = length(seq(as.POSIXct(dstart, format = '%d/%m/%Y'), as.POSIXct(dfinish, format = '%d/%m/%Y'), by = 'days')) * 24, ncol = 3), matrix of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#' }
#'
#' \strong{Outputs:}
#'
#' \code{ndays}{ - number of days for which predictions are made}\cr\cr
#' \code{longlat}{ - longitude and latitude for which simulation was run (decimal degrees)}\cr\cr
#' \code{dates}{ - vector of dates (POSIXct, UTC)}\cr\cr
#' \code{nyears}{ - number of years for which predictions are made}\cr\cr
#' \code{rainfall}{ - vector of daily rainfall (mm)}\cr\cr
#' \code{elevation}{ - elevation at point of simulation (m)}\cr\cr
#' \code{minimum_shade}{ - minimum shade for each day of simulation (\%)}\cr\cr
#' \code{maximum_shade}{ - maximum shade for each day of simulation (\%)}\cr\cr
#' \code{dem}{ - digital elevation model obtained via 'get_dev' using package 'elevatr' (m)}\cr\cr
#' \code{depths}{ - vector of depths used (cm)}\cr\cr
#' \code{SLOPE}{ - slope at point of simulation (\%)}\cr\cr
#' \code{ASPECT}{ - aspect at point of simulation (°, 0 is north)}\cr\cr
#' \code{HORIZON}{ - horizon angles at point of simulation (°)}\cr\cr
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
#' \item 13 solar_radiation - solar radiation (W/m2) (unshaded, adjusted for slope, aspect and horizon angle)
#' \item 14 sky_temperature - sky radiant temperature (°C)
#' \item 15 dew - dew fall (mm / h)
#' \item 16 frost - frost presence (0 or 1)
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
#' dstart <- "02/01/2017"
#' dfinish <- "30/12/2017"
#' loc <- c(-91.415669, -0.287145) # Isla Fernandina Galapagos
#' micro<-micro_ncep(loc = loc, dstart = dstart, dfinish = dfinish)
#'
#' micromet_lowshade<-as.data.frame(micro$micromet_lowshade) # above ground microclimatic conditions, min shade
#' soil_temperature_lowshade<-as.data.frame(micro$soil_temperature_lowshade) # soil temperatures, minimum shade
#' soil_moisture_lowshade<-as.data.frame(micro$soil_moisture_lowshade) # soil temperatures, minimum shade
#'
#' # append dates
#' tzone<-paste("Etc/GMT+",0,sep="")
#' dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")
#'
#' micromet_lowshade <- cbind(dates,micromet_lowshade)
#' soil_temperature_lowshade <- cbind(dates,soil_temperature_lowshade)
#' soil_moisture_lowshade <- cbind(dates, soil_moisture_lowshade)
#'
#' # plotting above-ground conditions in minimum shade
#' with(micromet_lowshade,{plot(air_temperature_local ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l",main=paste("air and sky temperature",sep=""), ylim = c(-20, 60))})
#' with(micromet_lowshade,{points(air_temperature_reference ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' , type = "l",lty=2,col='blue')})
#' with(micromet_lowshade,{points(sky_temperature ~ dates,xlab = "Date and Time", ylab = "Temperature (°C)"
#' ,  type = "l",col='light blue',main=paste("sky temperature",sep=""))})
#' with(micromet_lowshade,{plot(relative_humidity_local ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#' , type = "l",ylim=c(0,100),main=paste("humidity",sep=""))})
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
#'    ,col=i,type = "l",main=paste("soil temperature",sep=""))
#'  }else{
#'    points(soil_temperature_lowshade[,i+3]~soil_temperature_lowshade[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)",col=i,type = "l")
#'  }
#' }
#'
#' # plotting soil moisture
#' for(i in 1:10){
#'  if(i==1){
#'    plot(soil_moisture_lowshade[,i+3]*100~soil_moisture_lowshade[,1],xlab = "Date and Time", ylab = "Soil Moisture (% volumetric)"
#'    ,col=i,type = "l",main=paste("soil moisture",sep=""))
#'  }else{
#'    points(soil_moisture_lowshade[,i+3]*100~soil_moisture_lowshade[,1],xlab = "Date and Time", ylab = "Soil Moisture
#'     (%)",col=i,type = "l")
#'  }
#' }
micro_ncep <- function(
  loc = c(-5.3, 50.13),
  dstart = "01/01/2017",
  dfinish = "31/12/2017",
  dem = NA,
  dem2 = dem,
  dem_resolution = 30,
  pixels = 100,
  nyears = as.numeric(substr(dfinish, 7, 10)) - as.numeric(substr(dstart, 7, 10)) + 1,
  surface_reflectivity = 0.15,
  slope = NA,
  aspect = NA,
  depths = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200),
  minimum_shade = 0,
  maximum_shade = 90,
  local_height = 0.01,
  roughness_height_1 = 0,
  roughness_height_2 = 0,
  wind_profile_height_1 = 0,
  wind_profile_height_2 = 0,
  runshade = 1,
  clearsky = 0,
  global_aerosol_database = 1,
  solar_model_only = 0,
  initial_soil_temperature = NA,
  write_input = 0,
  output_to_csv = 0,
  reanalysis = TRUE,
  wind_multiplier = 1,
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
  rainwet = 1.5,
  organic_soil_cap = 1,
  precipitable_water = 1,
  horizon_angles = rep(NA, 24),
  soil_moisture_model = 1,
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
  leaf_resistance = 2e+06,
  stomatal_closure_potential = -1500,
  stomatal_stability_parameter = 10,
  moist_error = 1e-06,
  moist_count = 500,
  leaf_area_index = 0.1,
  microclima.leaf_area_index = 0,
  LOR = 1,
  snow_model = 1,
  snowtemp = 1.5,
  snowdens = 0.375,
  densfun = c(0.5979, 0.2178, 0.001, 0.0038),
  snowmelt = 1,
  undercatch = 1,
  rainmelt = 0.0125,
  shore = 0,
  tides = 0,
  deep_soil_temperature = NA,
  rainhour = 0,
  rainhourly = 0,
  rainoff = 0,
  radiation_per_wavelength = 0,
  scattered_uv = 0,
  max_iterations_per_day = 3,
  soilgrids = 0,
  longwave_radiation_model = 0,
  message = 0,
  fail = nyears * 24 * 365,
  spatial = NA,
  scenario = 0,
  terra_source = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data",
  save = 0,
  runmicro = 1,
  snowcond = 0,
  intercept = max(maximum_shade) / 100 * 0.3,
  grasshade = 0,
  coastal = F,
  hourlydata = NA,
  dailyprecip = NA,
  weather.elevation = 'ncep',
  cad.effects = TRUE,
  maxsurf = 85){ # end function parameters
  # --- 1. Input validation ---
  errors <- 0
  reference_height <- 2 # Reference height (m): NCEP data is at 2 m screen height

  # ncep-specific checks
  if (is.numeric(loc[1]) && (loc[1] > 180 | loc[2] > 90)) {
    message("ERROR: Latitude or longitude (longlat) is out of bounds. Please enter a correct value.")
    errors <- 1
  }
  if (global_aerosol_database == 1) {
    message("If program is crashing, try global_aerosol_database = 2.")
  }
  if (!(write_input %in% c(0, 1))) {
    message("ERROR: the variable 'write_input' must be either 0 or 1. Please correct.")
    errors <- 1
  }
  if (!(scenario %in% c(0, 2, 4))) {
    message("ERROR: Scenario must be either 0 (historical), 2 (plus 2) or 4 (plus 4).")
    errors <- 1
  }

  errors <- errors + validate_micro_inputs(depths, surface_reflectivity,
    slope  = ifelse(is.na(slope),  0, slope),
    aspect = ifelse(is.na(aspect), 0, aspect),
    horizon_angles   = if (is.na(horizon_angles[1])) rep(0, 24) else horizon_angles,
    surface_emissivity, tolerance, roughness_height, zero_plane_displacement, local_height, reference_height, orbital_eccentricity, precipitable_water,
    maxima_times = c(1, 1, 0, 0), minima_times = c(0, 0, 1, 1),
    minimum_shade, maximum_shade, mineral_conductivity, mineral_density, mineral_heat_capacity, bulk_density,
    global_aerosol_database = global_aerosol_database)

  if(errors==0){ # continue
    max.date <- as.Date(paste0("01/", format(Sys.time(), "%m/%Y")), format = "%d/%m/%Y", tz = "UTC")
    if(as.Date(dfinish, format = "%d/%m/%Y") > max.date){
      message(paste0("Cannot simulate past ", max.date, "; reducing timespan accordingly \n"))
      dfinish <- as.character(as.Date(max.date, format = "%Y/%m/%d"))
      dfinish <- paste(substr(dfinish, 9, 10), substr(dfinish, 6, 7), substr(dfinish, 1, 4), sep = "/")
    }
    ystart <- as.numeric(substr(dstart, 7, 10))
    yfinish <- as.numeric(substr(dfinish, 7, 10))
    yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)

    #remove trailing forward slash if necessary
    if(is.na(spatial)==FALSE){
      if(substr(x = spatial, start = nchar(spatial), nchar(spatial))=='/'){
        spatial <- substr(spatial, 1, nchar(spatial)-1)
      }
    }
    ################## time related variables #################################

    # for microclima calculations
    tme <- seq(as.Date(dstart, format = "%d/%m/%Y"), as.Date(dfinish, format = "%d/%m/%Y"), "days")

    day_of_year <- as.numeric(strftime(tme, format = "%j"))
    ndays<-length(day_of_year)
    end_day<-ndays
    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial max_iterations_per_day day burn in
    daystart<-1
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
    start_day <- 1 # start day

    ################## location and terrain #################################
    if (!require("terra", quietly = TRUE)) {
      stop("package 'terra' is needed. Please install it.",
           call. = FALSE)
    }
    if (!require("RNCEP", quietly = TRUE)) {
      stop("package 'RNCEP' is needed. Please install it.",
           call. = FALSE)
    }
    if (!require("RNetCDF", quietly = TRUE)) {
      stop("package 'RNetCDF' is needed. Please install it.",
           call. = FALSE)
    }
    if (!require("microclima", quietly = TRUE)) {
      stop("package 'microclima' is needed. Please install it via command: devtools::install_github('ilyamaclean/microclima').",
           call. = FALSE)
    }
    require("terra")
    require("RNCEP")
    require("RNetCDF")
    require("microclima")
    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    # get the local timezone reference longitude
    solar_noon_longitude <- abs(trunc(x[1]))
    hemisphere <- ifelse(x[2]<0, 2, 1) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    latitude_degrees <- abs(trunc(x[2]))
    latitude_minutes <- (abs(x[2])-latitude_degrees)*60
    longitude_degrees <- abs(trunc(x[1]))
    longitude_minutes <- (abs(x[1])-longitude_degrees)*60
    azmuth<-aspect
    lat <- as.numeric(longlat[2])
    long <- as.numeric(longlat[1])
    loc <- c(long, lat)
    if(class(dem)[1] == "SpatRaster"){
      cat('using DEM provided to function call \n')
    }
    if(save != 2 & class(dem)[1] != "SpatRaster"){
      cat('downloading DEM via package elevatr \n')
      dem <- microclima::get_dem(lat = lat, long = long, resolution = dem_resolution, xdims = pixels, ydims = pixels) # mercator equal area projection
    }
    if(save == 1){
      save(dem, file = 'dem.Rda')
    }
    if(save == 2){
      load('dem.Rda')
    }
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
    if(is.na(horizon_angles[1])){
      horizon_angles<-rep(0, 24)
      VIEWF <- 1 # incorporated already by microclima
    }else{
      VIEWF <- 1-sum(sin(as.data.frame(horizon_angles)*pi/180))/length(horizon_angles) # convert horizon angles to radians and calc view factor(s)
      HORIZON <- spline(x = horizon_angles, n = 36, method =  'periodic')$y
      HORIZON[HORIZON < 0] <- 0
      HORIZON[HORIZON > 90] <- 90
    }
    days <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), by = 'days')
    alldays <- seq(as.POSIXct("01/01/1900", format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), Sys.time()-60*60*24, by = 'days')
    startday <- which(as.character(format(alldays, "%d/%m/%Y")) == format(as.POSIXct(dstart, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), "%d/%m/%Y"))
    endday <- which(as.character(format(alldays, "%d/%m/%Y")) == format(as.POSIXct(dfinish, format = "%d/%m/%Y", origin = "01/01/1900", tz = 'UTC'), "%d/%m/%Y"))
    countday <- endday-startday+1
    dates <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(dfinish, format = "%d/%m/%Y", tz = 'UTC')+23*3600, by = 'hours')
    dates2 <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = 'UTC'), as.POSIXct(dfinish, format = "%d/%m/%Y", tz = 'UTC')+23*3600, by = 'days')
    if(save == 2){
      load('tref.Rda')
      load('SLOPE.Rda')
      load('ASPECT.Rda')
      load('HORIZON.Rda')
      elevation <- tref$elevation[1] # m
      elevation <- elevation
    }

    if(save != 2){
      ncquery <- function(filename, var, start, count, year){
        if (!require("ncdf4", quietly = TRUE)) {
          stop("package 'ncdf4' is needed. Please install it.",
               call. = FALSE)
        }
        if(!file.exists(paste(spatial, "/",filename,year,"_time.nc", sep = ""))){
          stop("looks like you need to generate the reordered version of your locally-stored NCEP data, see instructions on the NicheMapR Google Group under the heading 'speeding up reading of NCEP data'")
        }
        nc <- ncdf4::nc_open(paste(spatial, "/",filename,year,"_time.nc", sep = ""))
        if(nc$ndims == 3){
          start <- start[c(1,3:4)]
          count <- count[c(1,3:4)]
        }
        out <- as.numeric(ncdf4::ncvar_get(nc, varid = var, start = start, count = count))
        ncdf4::nc_close(nc)
        out
      }
      if(is.na(spatial) == FALSE){
        # now getting starting point and count for reading netcdf files
        cat(paste0("extracting weather data locally from ", spatial, " \n"))
        years <- as.numeric(unique(format(tme, "%Y")))
        nc <- RNetCDF::open.nc(paste(spatial, "/air.2m.gauss.", years[1], "_time.nc", sep = ""))
        lon2 <- matrix(RNetCDF::var.get.nc(nc, "lon", unpack = TRUE))
        lat2 <- matrix(RNetCDF::var.get.nc(nc, "lat", unpack = TRUE))
        lon_1 <- long
        if(lon_1 < 0){lon_1 <- 180 - (long*-1) + 180}
        lat_1 <- lat
        dist1 <- abs(lon2 - lon_1)
        index1 <- which.min(dist1)
        lon3 <- lon2[index1]
        dist2 <- abs(lat2 - lat_1)
        index2 <- which.min(dist2)
        lat3 <- lat2[index2]
        RNetCDF::close.nc(nc)
        if(lon3 > 180){lon3 <- lon3 - 180} # ensure longitude is -180 to 180

        if(class(hourlydata) == "logical"){
          prev.year <- 0
          if(format(tme[1], "%d/%m") == "01/01"){ # starting at beginning of year, need to start at end of previous year
            cat(paste("reading weather input for ", years[1]-1," (need a bit of the previous year) \n", sep = ""))
            prev.year <- 1
            start <- c(1460-3, 1, index2, index1) # for year prior to chosen years (getting the last day)
            count <- c(4, 1, 1, 1) # for year prior/year after chosen years (getting the last four or first four values, i.e. hours 0, 6, 12, 18)
            Tkmin1 <- ncquery("tmin.2m.gauss.", "tmin", start, count, years[1]-1)
            Tkmax1 <- ncquery("tmax.2m.gauss.", "tmax", start, count, years[1]-1)
            Tk1 <- ncquery("air.2m.gauss.", "air", start, count, years[1]-1)
            sh1 <- ncquery("shum.2m.gauss.", "shum", start, count, years[1]-1)
            pr1 <- ncquery("pres.sfc.gauss.", "pres", start, count, years[1]-1) # only three dimensions, hence
            tcdc1 <- ncquery("tcdc.eatm.gauss.", "tcdc", start, count, years[1]-1)
            dsw1 <- ncquery("dswrf.sfc.gauss.", "dswrf", start, count, years[1]-1)
            dlw1 <- ncquery("dlwrf.sfc.gauss.", "dlwrf", start, count, years[1]-1)
            ulw1 <- ncquery("ulwrf.sfc.gauss.", "ulwrf", start, count, years[1]-1)
            wu1 <- ncquery("uwnd.10m.gauss.", "uwnd", start, count, years[1]-1)
            wv1 <- ncquery("vwnd.10m.gauss.", "vwnd", start, count, years[1]-1)
            prate1 <- ncquery("prate.sfc.gauss.", "prate", start, count, years[1]-1)
          }

          for(j in 1:length(years)){
            tme3 <- seq(as.POSIXct(paste0("01/01/", years[j]), format = "%d/%m/%Y", tz = "UCT"), as.POSIXct(paste0("31/12/", years[j]), format = "%d/%m/%Y", tz = "UCT") + 18 * 3600, 6 * 3600)
            if(j == 1){
              if(format(tme[1], "%d/%m") == "01/01"){
                start <- c(1, 1, index2, index1) # for chosen years
                count <- c(-1, 1, 1, 1)# for chosen years
              }else{
                strt <- which(format(tme3, '%Y-%m-%d') == tme[1])[1]
                start <- c(strt - 4, 1, index2, index1) # for chosen years
                count <- c(-1, 1, 1, 1)# for chosen years
              }
              if(format(tme[length(tme)], "%d/%m") != "31/12" & length(years) == 1){
                cnt <- which(format(tme3, '%Y-%m-%d') == tme[length(tme)] + 1)[1] - 1
                count <- c(cnt + 5 - start[1], 1, 1, 1)# for chosen years
              }
              cat(paste("reading weather input for ", years[j]," \n", sep = ""))
              Tkmin <- ncquery("tmin.2m.gauss.", "tmin", start, count, years[j])
              Tkmax <- ncquery("tmax.2m.gauss.", "tmax", start, count, years[j])
              Tk <- ncquery("air.2m.gauss.", "air", start, count, years[j])
              sh <- ncquery("shum.2m.gauss.", "shum", start, count, years[j])
              pr <- ncquery("pres.sfc.gauss.", "pres", start, count, years[j]) # only three dimensions, hence
              tcdc <- ncquery("tcdc.eatm.gauss.", "tcdc", start, count, years[j])
              dsw <- ncquery("dswrf.sfc.gauss.", "dswrf", start, count, years[j])
              dlw <- ncquery("dlwrf.sfc.gauss.", "dlwrf", start, count, years[j])
              ulw <- ncquery("ulwrf.sfc.gauss.", "ulwrf", start, count, years[j])
              wu <- ncquery("uwnd.10m.gauss.", "uwnd", start, count, years[j])
              wv <- ncquery("vwnd.10m.gauss.", "vwnd", start, count, years[j])
              prate <- ncquery("prate.sfc.gauss.", "prate", start, count, years[j])
            }else{
              start <- c(1, 1, index2, index1) # for chosen years
              if(format(tme[length(tme)], "%d/%m") == "31/12" | j != length(years)){
                count <- c(-1, 1, 1, 1)# for chosen years
              }else{
                cnt <- which(format(tme3, '%Y-%m-%d') == tme[length(tme)] + 1)[1] - 1
                count <- c(cnt + 5 - start[1], 1, 1, 1)# for chosen years
              }
              cat(paste("reading weather input for ", years[j]," \n", sep = ""))
              Tkmin <- c(Tkmin, ncquery("tmin.2m.gauss.", "tmin", start, count, years[j]))
              Tkmax <- c(Tkmax, ncquery("tmax.2m.gauss.", "tmax", start, count, years[j]))
              Tk <- c(Tk, ncquery("air.2m.gauss.", "air", start, count, years[j]))
              sh <- c(sh, ncquery("shum.2m.gauss.", "shum", start, count, years[j]))
              pr <- c(pr, ncquery("pres.sfc.gauss.", "pres", start, count, years[j]))
              tcdc <- c(tcdc, ncquery("tcdc.eatm.gauss.", "tcdc", start, count, years[j]))
              dsw <- c(dsw, ncquery("dswrf.sfc.gauss.", "dswrf", start, count, years[j]))
              dlw <- c(dlw, ncquery("dlwrf.sfc.gauss.", "dlwrf", start, count, years[j]))
              ulw <- c(ulw, ncquery("ulwrf.sfc.gauss.", "ulwrf", start, count, years[j]))
              wu <- c(wu, ncquery("uwnd.10m.gauss.", "uwnd", start, count, years[j]))
              wv <- c(wv, ncquery("vwnd.10m.gauss.", "vwnd", start, count, years[j]))
              prate <- c(prate, ncquery("prate.sfc.gauss.", "prate", start, count, years[j]))
            }
          }

          next.year <- 0
          if(format(tme[length(tme)], "%d/%m") == "31/12"){ # starting at beginning of year, need to start at end of previous year
            cat(paste("reading weather input for ", years[length(years)]+1," (need a bit of the next year) \n", sep = ""))
            next.year <- 1
            start <- c(1, 1, index2, index1) # for year prior to chosen years (getting the last day)
            count <- c(4, 1, 1, 1) # for year prior/year after chosen years (getting the last four or first four values, i.e. hours 0, 6, 12, 18)
            Tkmin2 <- ncquery("tmin.2m.gauss.", "tmin", start, count, years[length(years)]+1)
            Tkmax2 <- ncquery("tmax.2m.gauss.", "tmax", start, count, years[length(years)]+1)
            Tk2 <- ncquery("air.2m.gauss.", "air", start, count, years[length(years)]+1)
            sh2 <- ncquery("shum.2m.gauss.", "shum", start, count, years[length(years)]+1)
            pr2 <- ncquery("pres.sfc.gauss.", "pres", start, count, years[length(years)]+1) # only three dimensions, hence
            tcdc2 <- ncquery("tcdc.eatm.gauss.", "tcdc", start, count, years[length(years)]+1)
            dsw2 <- ncquery("dswrf.sfc.gauss.", "dswrf", start, count, years[length(years)]+1)
            dlw2 <- ncquery("dlwrf.sfc.gauss.", "dlwrf", start, count, years[length(years)]+1)
            ulw2 <- ncquery("ulwrf.sfc.gauss.", "ulwrf", start, count, years[length(years)]+1)
            wu2 <- ncquery("uwnd.10m.gauss.", "uwnd", start, count, years[length(years)]+1)
            wv2 <- ncquery("vwnd.10m.gauss.", "vwnd", start, count, years[length(years)]+1)
            prate2 <- ncquery("prate.sfc.gauss.", "prate", start, count, years[length(years)]+1)
          }

          if(prev.year == 1){
            Tkmin <- c(Tkmin1, Tkmin)
            Tkmax <- c(Tkmax1, Tkmax)
            Tk <- c(Tk1, Tk)
            sh <- c(sh1, sh)
            pr <- c(pr1, pr)
            tcdc <- c(tcdc1, tcdc)
            dsw <- c(dsw1, dsw)
            dlw <- c(dlw1, dlw)
            ulw <- c(ulw1, ulw)
            wu <- c(wu1, wu)
            wv <- c(wv1, wv)
            prate <- c(prate1, prate)
          }
          if(next.year == 1){
            Tkmin <- c(Tkmin, Tkmin2)
            Tkmax <- c(Tkmax, Tkmax2)
            Tk <- c(Tk, Tk2)
            sh <- c(sh, sh2)
            pr <- c(pr, pr2)
            tcdc <- c(tcdc, tcdc2)
            dsw <- c(dsw, dsw2)
            dlw <- c(dlw, dlw2)
            ulw <- c(ulw, ulw2)
            wu <- c(wu, wu2)
            wv <- c(wv, wv2)
            prate <- c(prate, prate2)
          }
          dsw[dsw < 0] <- 0
          prate[prate < 0] <- 0
          prate <- prate * 3600 * 6 # mm in 6 hrs
          obs_time <- seq(as.POSIXct(dstart, format = "%d/%m/%Y", tz = "UCT") - 3600 * 24, as.POSIXct(dfinish, format = "%d/%m/%Y", tz = "UCT") + 18 * 3600 + 24 * 3600, 6 * 3600)
          ncepdata <- data.frame(obs_time = obs_time, Tk, Tkmin, Tkmax, sh, pr, wu, wv, dlw, ulw, dsw, tcdc) # 6-hourly ncep for chosen period plus a day added either side for interpolation
          hourlydata <- microclima::hourlyNCEP(ncepdata = ncepdata, lat, long, tme, reanalysis) # interpolated to hourly
          cat("computing radiation and elevation effects with package microclima \n")
          microclima.out <- microclima::microclimaforNMR(lat = longlat[2], long = longlat[1], dstart = dstart, dfinish = dfinish, l = mean(microclima.leaf_area_index), x = LOR, coastal = coastal, hourlydata = hourlydata, dailyprecip = prate, dem = dem, demmeso = dem2, albr = 0, resolution = 30, zmin = 0, slope = slope, aspect = aspect, windthresh = 4.5, emthresh = 0.78, reanalysis2 = reanalysis, difani = FALSE, weather.elevation = weather.elevation, cad.effects = cad.effects)
        }else{
          cat("computing radiation and elevation effects with package microclima \n")
          microclima.out <- microclima::microclimaforNMR(lat = longlat[2], long = longlat[1], dstart = dstart, dfinish = dfinish, l = mean(microclima.leaf_area_index), x = LOR, coastal = coastal, hourlydata = hourlydata, dailyprecip = NA, dem = dem, demmeso = dem2, albr = 0, resolution = 30, zmin = 0, slope = slope, aspect = aspect, windthresh = 4.5, emthresh = 0.78, reanalysis2 = reanalysis, difani = FALSE, weather.elevation = weather.elevation, cad.effects = cad.effects)
        }
        if(class(dailyprecip) == "logical"){
          dailyprecip <- microclima.out$dailyprecip[-c(1:4)] # remove extra 4 values from start
          dailyprecip <- dailyprecip[1:(length(dailyprecip) - 4)] # remove extra 5 values from end
          dailyprecip <- aggregate(dailyprecip, by = list(format(hourlydata$obs_time[seq(1, nrow(hourlydata), 6)], "%Y-%m-%d")), sum)$x
        }
      }else{
        if(class(hourlydata) == "logical"){
          cat("downloading weather data with package RNCEP via package microclima \n")
          hourlydata <- microclima::hourlyNCEP(ncepdata = NA, lat, long, tme, reanalysis) # interpolated to hourly
        }
        cat("computing radiation and elevation effects with package microclima \n")
        microclima.out <- microclima::microclimaforNMR(lat = longlat[2], long = longlat[1], dstart = dstart, dfinish = dfinish, l = mean(microclima.leaf_area_index), x = LOR, coastal = coastal, hourlydata = hourlydata, dailyprecip = dailyprecip, dem = dem, demmeso = dem2, albr = 0, resolution = 30, zmin = 0, slope = slope, aspect = aspect, windthresh = 4.5, emthresh = 0.78, reanalysis2 = reanalysis, difani = FALSE, weather.elevation = weather.elevation, cad.effects = cad.effects)
        dailyprecip <- microclima.out$dailyprecip
      }
      hourlyradwind <- microclima.out$hourlyradwind
      SLOPE <- hourlyradwind$slope[1]
      ASPECT <- hourlyradwind$aspect[1]
      tref <- microclima.out$tref
      zenith_angle_hourly <- hourlydata$szenith
      zenith_angle_hourly[zenith_angle_hourly > 90] <- 90
      HORIZON <- horizon_angles
      if(save == 1){
        save(SLOPE, file = 'SLOPE.Rda')
        save(ASPECT, file = 'ASPECT.Rda')
        save(HORIZON, file = 'HORIZON.Rda')
      }
      # NB units for rad = MJ / m^2 / hr (divide by 0.0036 to get to W / m^2)
      # Skyviewfact (time invariant, 1 = complete hemisphere in view)
      # canopyfact (proportion of isotropic radiation blocked out, 1 = no radiation gets in)
      if(save == 1){
        save(tref, file = 'tref.Rda')
      }
      # tref (original hourly NCEP)
      # telev (elevation-corrected temperature)
      # tcad (delta temperature due to cold air drainage)
      elevation <- tref$elevation[1] # m
      elevation <- elevation
      reference_temperature <- tref$telev + tref$tcad # reference Tair corrected for lapse rate and cold air drainage

      if(scenario != 0){
        TAIRhr_orig <- reference_temperature
        yearstodo <- seq(ystart, yfinish)
        nyears <- yfinish - ystart + 1
        if(yfinish > 2025){
          ystart_terra <- 2025 - nyears + 1
          yfinish_terra <- 2025
          message(paste0("terraclimate climate change scenarios are for 1950 to 2025 - using ", ystart, "-", yfinish), '\n')
        }else{
          ystart_terra <- ystart
          yfinish_terra <- yfinish
        }
        leapyears <- seq(1900, 2060, 4)
        message("generate climate change scenario", '\n')
        # diff spline function
        getdiff <- function(diffs){
          diff1 <- (unlist(diffs[1]) + unlist(diffs[length(diffs)])) / 2
          # generate list of days
          leapcount <- 0
          for(ys in 1:nyears){
            if(ys == 1){
              if(nyears == 1){
                if(yearstodo[ys] %in% leapyears){
                  leapcount <- leapcount + 1
                  day<-c(1, 15.5, 45.5, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5, 366)
                }else{
                  day<-c(1, 15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 365)
                }
              }else{
                if(yearstodo[ys] %in% leapyears){
                  leapcount <- leapcount + 1
                  day<-c(1, 15.5, 45.5, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5)
                }else{
                  day<-c(1, 15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5)
                }
              }
            }else{
              if(ys == nyears){
                if(yearstodo[ys] %in% leapyears){
                  leapcount <- leapcount + 1
                  day<-c(15.5, 45.5, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5, 366)
                }else{
                  day<-c(15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 365)
                }
              }else{
                if(yearstodo[ys] %in% leapyears){
                  leapcount <- leapcount + 1
                  day<-c(15.5, 45.5, 75.5, 106, 136.5, 167, 197.5, 228.5, 259, 289.5, 320, 350.5)
                }else{
                  day<-c(15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5)
                }
              }
            }
            if(ys == 1){
              days2 <- day
              days <- day
            }else{
              if(yearstodo[ys] %in% leapyears){
                days2 <- c(days2, (day + 366 * (ys - 1)) + leapcount)
              }else{
                days2 <- c(days2, (day + 365 * (ys - 1)) + leapcount)
              }
              days <- c(days, day)
            }
          }

          diffs3 <- c(diff1, diffs, diff1)
          days_diffs <- data.frame(matrix(NA, nrow = nyears * 12 + 2, ncol = 3))
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

        terra <- as.data.frame(get_terra(x = loc, ystart = ystart_terra, yfinish = yfinish_terra, source = terra_source))
        if(is.infinite(max(terra))){
          message("no TerraClimate data available", '\n')
          stop()
        }
        if(scenario == 4){
          terra_cc <- as.data.frame(get_terra(x = loc, ystart = ystart_terra, yfinish = yfinish_terra, scenario = 4, source = terra_source))
        }
        if(scenario == 2){
          terra_cc <- as.data.frame(get_terra(x = loc, ystart = ystart_terra, yfinish = yfinish_terra, scenario = 2, source = terra_source))
        }

        tmin_diffs <- terra_cc$reference_temperature_min - terra$reference_temperature_min
        tmax_diffs <- terra_cc$reference_temperature_max - terra$reference_temperature_max
        precip_diffs <- terra_cc$rainfall / terra$rainfall
        srad_diffs <- terra_cc$SRAD / terra$SRAD
        vpd_diffs <- terra_cc$VPD / terra$VPD
        srad_diffs[is.na(srad_diffs)] <- 1
        vpd_diffs[is.na(vpd_diffs)] <- 1
        precip_diffs[is.na(precip_diffs)] <- 1
        srad_diffs[is.infinite(srad_diffs)] <- 1
        vpd_diffs[is.infinite(vpd_diffs)] <- 1
        precip_diffs[is.infinite(precip_diffs)] <- 1

        TMINN_diff <- getdiff(tmin_diffs)
        TMAXX_diff <- getdiff(tmax_diffs)
        VPD_diff <- getdiff(vpd_diffs)
        SOLAR_diff <- getdiff(srad_diffs)
        RAIN_diff <- getdiff(precip_diffs)

        TMINN_diff <- rep(TMINN_diff, each=24)[1:length(reference_temperature)]
        TMAXX_diff <- rep(TMAXX_diff, each=24)[1:length(reference_temperature)]
        VPD_diff <- rep(VPD_diff, each=24)[1:length(reference_temperature)]
        SOLAR_diff <- rep(SOLAR_diff, each=24)[1:length(reference_temperature)]


        # code to apply changes in min and max air temperature to hourly air temperature data

        # find times of maxima and minima
        mins <- aggregate(reference_temperature, by = list(format(dates, "%Y-%m-%d")), FUN = min)$x
        maxs <- aggregate(reference_temperature, by = list(format(dates, "%Y-%m-%d")), FUN = max)$x
        mins <- rep(mins, each=24)
        maxs <- rep(maxs, each=24)
        mins2 <- reference_temperature / mins[1:length(reference_temperature)]
        maxs2 <- reference_temperature / maxs
        mins2[mins2 != 1] <- 0
        maxs2[maxs2 != 1] <- 0

        # construct weightings so that changes in min and max air temp can be applied
        minweight <- rep(NA, length(mins2))
        maxweight <- minweight
        minweight[mins2 == 1] <- 1
        minweight[maxs2 == 1] <- 0
        maxweight[mins2 == 1] <- 0
        maxweight[maxs2 == 1] <- 1
        minweight <- zoo::na.approx(minweight, na.rm = FALSE)
        maxweight <- zoo::na.approx(maxweight, na.rm = FALSE)
        minweight[is.na(minweight)] <- 0.5
        maxweight[is.na(maxweight)] <- 0.5
        minweight[minweight > 1] <- 1
        minweight[minweight < 0] <- 0
        maxweight[maxweight > 1] <- 1
        maxweight[maxweight < 0] <- 0

        # construct final weighted correction factor and apply
        diff <- TMINN_diff * minweight + TMAXX_diff * maxweight
        reference_temperature <- reference_temperature + diff
      }
      global_radiation <- hourlyradwind$swrad / 0.0036
      if(scenario != 0){
       global_radiation <- global_radiation * SOLAR_diff
      }
      global_radiation[global_radiation < 0] <- 0
      global_radiation <- zoo::na.approx(global_radiation)
      longwave_radiation <- hourlydata$downlong / 0.0036
      if(longwave_radiation_model != 2){
        longwave_radiation <- longwave_radiation * 0 + -1 # make negative so computed internally in the microclimate model
      }
      cloud_cover <- hourlydata$cloudcover
      if(scenario != 0){
        cloud_cover <- cloud_cover * (1 + 1 - SOLAR_diff)
      }
      cloud_cover[cloud_cover < 0] <- 0
      cloud_cover[cloud_cover > 100] <- 100
      reference_humidity <- suppressWarnings(humidityconvert(h = hourlydata$humidity, intype = 'specific', p = hourlydata$pressure, tc = reference_temperature)$relative)
      reference_humidity[reference_humidity > 100] <- 100
      reference_humidity[reference_humidity < 0] <- 0
      if(scenario != 0){
       VPDhr <- WETAIR(db = TAIRhr_orig, rh = 100)$e - WETAIR(db = TAIRhr_orig, rh = reference_humidity)$e
       VPDhr <- VPDhr - mean(VPD_diff) # note using mean here because otherwise can lead to unrealistically low RH which then stronly affects TSKY
       e <- WETAIR(db = reference_temperature, rh = 100)$e - VPDhr
       es <- WETAIR(db = reference_temperature, rh = 100)$esat
       reference_humidity <- (e / es) * 100
       reference_humidity[reference_humidity > 100] <- 100
       reference_humidity[reference_humidity < 0] <- 0
      }
      reference_wind_speed <- hourlyradwind$windspeed
      reference_wind_speed[is.na(reference_wind_speed)] <- 0.1
      if(rainhourly == 0){
        rainfall_hourly <- rep(0, 24 * ndays)
      }else{
        rainfall_hourly <- rainhour
      }
      if(length(rainhourly) == length(tme)){ # an hourly rainfall vector has been provided
        aggvec <- rep(1:(length(tme)/24), 24)
        aggvec <- aggvec[order(aggvec)]
        rainfall <- aggregate(rainhourly, by = aggvec, FUN = 'sum') # aggregate to daily totals
      }else{
        rainfall <- dailyprecip
      }
      PRESShr <- hourlydata$pressure
      if(scenario != 0){
        rainfall <- rainfall * RAIN_diff[1:length(rainfall)]
      }
      rainfall[rainfall < 0.1] <- 0
      if(rainhourly == 0){ # putting daily rainfall on midnight (account for local solar time)
        midnight <- which(microclima.out$hourlydata$szenith[1:24] == max(microclima.out$hourlydata$szenith[1:24]))
        if(midnight <= 24){
         midnights <- seq(midnight, length(rainfall_hourly / 24) - 1, 24)
        }else{
         midnights <- seq(midnight, length(rainfall_hourly / 24), 24)
        }
        rainfall_hourly[midnights+1] <- rainfall/3
        rainfall_hourly[midnights] <- rainfall/3
        rainfall_hourly[midnights-1] <- rainfall/3
        rainhourly <- 1
      }
      ZENhr2 <- zenith_angle_hourly
      ZENhr2[ZENhr2!=90] <- 0
      dmaxmin <- function(x, fun) {
        dx <- t(matrix(x, nrow = 24))
        apply(dx, 1, fun)
      }
      reference_temperature_max <- dmaxmin(reference_temperature, max)
      reference_temperature_min <- dmaxmin(reference_temperature, min)
      cloud_max <- dmaxmin(cloud_cover, max)
      cloud_min <- dmaxmin(cloud_cover, min)
      reference_humidity_max <- dmaxmin(reference_humidity, max)
      reference_humidity_min <- dmaxmin(reference_humidity, min)
      reference_wind_max <- dmaxmin(reference_wind_speed, max)
      reference_wind_min <- dmaxmin(reference_wind_speed, min)
      PRESS <- dmaxmin(PRESShr, min)

      if(clearsky == 1){
        cat("running micro_global to get clear sky solar \n")
        micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, timeinterval = 365, solar_model_only = 1)
        clearskyrad <- micro_clearsky$micromet_lowshade[, c(1, 13)]
        leapyears <- seq(1900, 2060, 4)
        ystart <- as.numeric(substr(dstart, 7, 10))
        yfinish <- as.numeric(substr(dfinish, 7, 10))
        datestimes <- seq(as.POSIXct(paste0(ystart, '-01-01 00:00:00')), as.POSIXct(paste0(yfinish, '-12-31 23:00:00')), by = 'hours')
        for(year in ystart:yfinish){
          clear <- clearskyrad[, 2]
          if(year %in% leapyears){
            clear <- c(clear[1:1416], clear[1417:(1416 + 24)], clear[(1417 + 1):length(clear)])
          }
          if(year == ystart){
            global_radiation <- clear
          }else{
            global_radiation <- c(global_radiation, clear)
          }
        }
        tzoff <- 12-which(zenith_angle_hourly[1:24] ==min(zenith_angle_hourly[1:24]))+1
        if(tzoff != 12){
          if(tzoff > 12){
            global_radiation <- c(global_radiation[(length(global_radiation)-tzoff+1):(length(global_radiation))], global_radiation[1:(length(global_radiation)-tzoff)])
          }else{
            global_radiation <- c(global_radiation[(tzoff+1):length(global_radiation)], global_radiation[1:tzoff])
          }
        }
      }
      if(save == 1){
        cat("saving met data for later \n")
        save(cloud_max, file = 'cloud_max.Rda')
        save(cloud_min, file = 'cloud_min.Rda')
        save(reference_wind_max, file = 'reference_wind_max.Rda')
        save(reference_wind_min, file = 'reference_wind_min.Rda')
        save(reference_temperature_max, file = 'reference_temperature_max.Rda')
        save(reference_temperature_min, file = 'reference_temperature_min.Rda')
        save(reference_humidity_max, file = 'reference_humidity_max.Rda')
        save(reference_humidity_min, file = 'reference_humidity_min.Rda')
        save(rainfall, file = 'rainfall.Rda')
        save(PRESS, file = 'PRESS.Rda')
        save(cloud_cover, file = 'cloud_cover.Rda')
        save(reference_wind_speed, file = 'reference_wind_speed.Rda')
        save(reference_temperature, file = 'reference_temperature.Rda')
        save(reference_humidity, file = 'reference_humidity.Rda')
        save(rainfall_hourly, file = 'rainfall_hourly.Rda')
        save(global_radiation, file = 'global_radiation.Rda')
        save(zenith_angle_hourly, file = 'zenith_angle_hourly.Rda')
        save(longwave_radiation, file = 'longwave_radiation.Rda')
        save(microclima.out, file = 'microclima.out.Rda')
      }
    }else{
      cat("loading met data from previous run \n")
      load('cloud_max.Rda')
      load('cloud_min.Rda')
      load('reference_wind_max.Rda')
      load('reference_wind_min.Rda')
      load('reference_temperature_max.Rda')
      load('reference_temperature_min.Rda')
      load('reference_humidity_max.Rda')
      load('reference_humidity_min.Rda')
      load('rainfall.Rda')
      load('PRESS.Rda')
      load('cloud_cover.Rda')
      load('reference_wind_speed.Rda')
      load('reference_temperature.Rda')
      load('reference_humidity.Rda')
      load('rainfall_hourly.Rda')
      load('global_radiation.Rda')
      load('zenith_angle_hourly.Rda')
      load('longwave_radiation.Rda')
      load('microclima.out.Rda')
      rainhourly <- 1
    }
    slope <- 0 # already corrected for by microclima
    azmuth <- 0 # already corrected for by microclima

    aerosol_optical_depth <- compute_tai(longlat, global_aerosol_database, TAI_AUSTRALIA)

    if(max(abs(air_temperature_offset)) != 0){
      # impose uniform temperature change
      reference_temperature_max <- reference_temperature_max + air_temperature_offset
      reference_temperature_min <- reference_temperature_min + air_temperature_offset
      if(length(air_temperature_offset) == length(reference_temperature_max)){
        air_temperature_offset.hr <- rep(air_temperature_offset, each = 24)
      }else{
        air_temperature_offset.hr <- air_temperature_offset
      }
      reference_temperature <- reference_temperature + air_temperature_offset.hr
      sigma <- 5.67e-8 #Stefan-Boltzman, W/(m.K)
      if(longwave_radiation_model == 2){
        longwave_radiation <- sigma * ((longwave_radiation / sigma) ^ (1 / 4) + air_temperature_offset.hr) ^ 4 # adjust downward radiation for altered 'sky temperature'
      }
    }
    rainfall <- rainfall + rainoff
    rainfall[rainfall < 0] <- 0
    if(length(rainoff) == length(rainfall)){
      rainoff.hr <- rep(rainoff, each = 24)
    }else{
      rainoff.hr <- rainoff
    }
    rainfall_hourly <- rainfall_hourly + rainoff.hr
    rainfall_hourly[rainfall_hourly < 0] <- 0
    ALLMINTEMPS <- reference_temperature_min
    ALLMAXTEMPS <- reference_temperature_max
    ALLTEMPS <- cbind(ALLMAXTEMPS, ALLMINTEMPS)

    reference_wind_max <- reference_wind_max * wind_multiplier
    reference_wind_min <- reference_wind_min * wind_multiplier
    reference_wind_speed <- reference_wind_speed * wind_multiplier

    albedo <- rep(surface_reflectivity, ndays)
    soil_wetness <- rep(soil_wetness, ndays)
    soilwet <- rainfall
    soilwet[soilwet <= rainwet] <- 0
    soilwet[soilwet > 0] <- 90
    if(ndays < 1){
      soil_wetness <- pmax(soilwet, soil_wetness)
    }

    Intrvls<-rep(0, ndays)
    Intrvls[1] <- 1 # user-supplied last day-of-year in each time interval sequence
    Numtyps <- 10 # number of substrate types
    soil_nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
    soil_nodes[1:10,] <- c(1:10) # deepest nodes for each substrate type
    solar_noon_longitude <- abs(trunc(x[1]))

    hemisphere <- ifelse(x[2] < 0, 2, 1)
    latitude_degrees <- abs(trunc(x[2]))
    latitude_minutes <- (abs(x[2]) - latitude_degrees) * 60
    longitude_degrees <- abs(trunc(x[1]))
    longitude_minutes <- (abs(x[1]) - longitude_degrees) * 60

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

    if(is.na(deep_soil_temperature)){
      if(nyears == 1){
        avetemp <- (sum(reference_temperature_max) + sum(reference_temperature_min)) / (length(reference_temperature_max)*2)
        deep_soil_temperature <-rep(avetemp, ndays)
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

    surface_emissivity <- matrix(nrow = ndays, data = 0)
    surface_emissivity <- surface_emissivity + surface_emissivity

    moists2 <- matrix(nrow = 10, ncol = ndays, data = 0)
    moists2[1, ndays] <- 0.2
    soil_moisture_profile <- moists2

    if(soil_moisture_model == 1){
      moists2 <- matrix(nrow = 10, ncol = ndays, data = 0) # set up an empty vector for soil moisture values through time
      moists2[1:10, ] <- initial_soil_moisture
      soil_moisture_profile <- moists2
    }
    soil_properties <- matrix(data = 0, nrow = 10, ncol = 5)
    soil_properties[,1] <- bulk_density
    soil_properties[,2] <- 1 - bulk_density / mineral_density # not used if soil moisture computed
    soil_properties[soil_properties[,2] < 0.26, 2] <- 0.26
    soil_properties[,3] <- mineral_conductivity
    soil_properties[,4] <- mineral_heat_capacity
    soil_properties[,5] <- mineral_density
    if(organic_soil_cap == 1){
      soil_properties[1:2, 3] <- 0.2
      soil_properties[1:2, 4] <- 1920
    }
    if(organic_soil_cap==2){
      soil_properties[1:2, 3] <- 0.1
      soil_properties[3:4, 3] <- 0.25
      soil_properties[1:4, 4] <- 1920
      soil_properties[1:4, 5] <- 1.3
      soil_properties[1:4, 1] <- 0.7
    }

    elevation <- as.numeric(elevation)
    solar_noon_longitude <- as.numeric(solar_noon_longitude)
    longitude_minutes <- as.numeric(longitude_minutes)
    longitude_degrees <- as.numeric(longitude_degrees)
    latitude_minutes <- as.numeric(latitude_minutes)
    latitude_degrees <-as.numeric(latitude_degrees)
    hourly <- 1

    maxima_times <- c(1, 1, 0, 0)
    minima_times <- c(0, 0, 1, 1)
    # --- 3. Build microclimate model input list ---
    micro_input <- build_microinput(ndays, roughness_height, tolerance, local_height, reference_height, Numtyps,
      roughness_height_1, roughness_height_2, wind_profile_height_1, wind_profile_height_2, start_day, end_day, hemisphere, latitude_degrees, latitude_minutes, longitude_degrees, longitude_minutes,
      solar_noon_longitude, slope, azmuth, elevation, precipitable_water, microdaily, mean_annual_temperature, orbital_eccentricity, VIEWF,
      snowtemp, snowdens, snowmelt, undercatch, rain_multiplier, runshade, soil_moisture_model,
      maximum_pooling_depth, evenrain, snow_model, rainmelt, output_to_csv, densfun, hourly,
      rainhourly, radiation_per_wavelength, scattered_uv, root_resistance, stomatal_closure_potential, leaf_resistance, stomatal_stability_parameter, root_radius, moist_error, moist_count, longwave_radiation_model, message,
      fail, snowcond, intercept, grasshade, solar_model_only, canopy_roughness_height, zero_plane_displacement, maxima_times, minima_times,
      spinup, dewrain = 0, moiststep = 360, maxsurf, max_iterations_per_day)

    if(length(leaf_area_index) < ndays){
      leaf_area_index<-rep(leaf_area_index[1], ndays)
    }
    if(shore == 0){
      tides <- matrix(data = 0, nrow = 24 * ndays, ncol = 3) # make an empty matrix
    }

    micro <- build_micro_list(micro_input, day_of_year, surface_emissivity, depths, soil_nodes, maximum_shade_daily,
      minimum_shade_daily, reference_temperature_max, reference_temperature_min, reference_humidity_max, reference_humidity_min, cloud_max, cloud_min, reference_wind_max, reference_wind_min,
      reference_temperature, reference_humidity, reference_wind_speed, cloud_cover, global_radiation, rainfall_hourly, zenith_angle_hourly, longwave_radiation, albedo, soil_wetness,
      soilinit, horizon_angles, aerosol_optical_depth, soil_properties, soil_moisture_profile, rainfall, deep_soil_temperature = deep_soil_temperature,
      air_entry_water_potential, saturated_hydraulic_conductivity, campbell_b_parameter, soil_bulk_density, soil_mineral_density, root_density, leaf_area_index, tides)

    # --- 4. Write inputs to CSV (optional) ---
    if (write_input == 1) {
      write_micro_csv(micro_input, day_of_year, surface_emissivity, depths, soil_nodes, maximum_shade_daily, minimum_shade_daily,
        maxima_times = c(1, 1, 0, 0), minima_times = c(0, 0, 1, 1),
        reference_temperature_max, reference_temperature_min, reference_humidity_max, reference_humidity_min, cloud_max, cloud_min, reference_wind_max, reference_wind_min,
        albedo, soil_wetness, soilinit, horizon_angles, aerosol_optical_depth, soil_properties, soil_moisture_profile, rainfall,
        deep_soil_temperature, air_entry_water_potential, soil_bulk_density, soil_mineral_density, campbell_b_parameter, saturated_hydraulic_conductivity, root_density, leaf_area_index, tides,
        reference_temperature, reference_humidity, reference_wind_speed, cloud_cover, global_radiation, rainfall_hourly, zenith_angle_hourly, longwave_radiation)
    }
    if(is.numeric(loc[1])){
      location <- paste("long", loc[1], "lat", loc[2])
    }else{
      location <- loc
    }
    if(runmicro){
      cat(paste('running microclimate model for ', ndays, ' days from ', dates[1], ' to ', dates[length(dates)], ' at site ', location, '\n'))
    cat('Note: the output column `SOLR` in micromet_lowshade and micromet_highshade is for unshaded solar radiation adjusted for slope, aspect and horizon angle \n')
    ptm <- proc.time() # Start timing
    microut<-microclimate(micro)
    print(proc.time() - ptm) # Stop the clock

    # --- 5. Process and return results ---
    out <- process_micro_output(microut, soil_moisture_model, snow_model, radiation_per_wavelength)
    if (max(out$micromet_lowshade[, 1] == 0)) {
      message("ERROR: the model crashed - try a different error tolerance (tolerance) or a different spacing in depths")
    }
    return(build_micro_return(out, rainfall, ndays, elevation, surface_reflectivity,
      longlat = longlat, nyears, timeinterval = ndays, minimum_shade_daily, maximum_shade_daily,
      depths, dates, dates2, air_entry_water_potential, soil_bulk_density, soil_mineral_density, campbell_b_parameter, saturated_hydraulic_conductivity, dem = dem, diffuse_frac = NA,
      snow_model, radiation_per_wavelength,
      extra = list(SLOPE = SLOPE, ASPECT = ASPECT, HORIZON = HORIZON,
                   microclima.out = microclima.out)))
  }else{
    return(list(rainfall = rainfall, reference_temperature_max = reference_temperature_max, reference_temperature_min = reference_temperature_min, reference_humidity_max = reference_humidity_max, reference_humidity_min = reference_humidity_min, reference_wind_max = reference_wind_max, reference_wind_min = reference_wind_min, cloud_max = cloud_max, cloud_min = cloud_min, cloud_cover = cloud_cover, reference_wind_speed = reference_wind_speed, reference_temperature = reference_temperature, reference_humidity = reference_humidity, rainfall_hourly = rainfall_hourly, global_radiation = global_radiation, zenith_angle_hourly = zenith_angle_hourly, longwave_radiation = longwave_radiation, dates = dates, dates2 = dates2, air_entry_water_potential=air_entry_water_potential,soil_bulk_density=soil_bulk_density,soil_mineral_density=soil_mineral_density,campbell_b_parameter=campbell_b_parameter,saturated_hydraulic_conductivity=saturated_hydraulic_conductivity))
  }
  } # end error trapping
} # end of micro_ncep function
