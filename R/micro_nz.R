#' New Zealand implementation of the microclimate model.
#'
#' An implementation of the NicheMapR microclimate model that uses the Virtual Climate Station Network (VCSN) for NZ https://www.niwa.co.nz/climate/our-services/virtual-climate-stations
#' @encoding UTF-8
#' @param loc Longitude and latitude (decimal degrees)
#' @param ystart First year to run
#' @param yfinish Last year to run
#' @param surface_reflectivity Soil solar reflectance (decimal \%)
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
#' @usage micro_aust(loc = "Melbourne, Australia", ystart = 1990, yfinish = 1990,
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
#' \code{terrain}{ = 0, Use 250m resolution terrain data? 1=yes, 0=no}\cr\cr
#' \code{dailywind}{ = 1, Make Fortran code write output as csv files? 1=yes, 0=no}\cr\cr
#' \code{wind_multiplier}{ = 1, factor to multiply wind speed by e.g. to simulate forest}\cr\cr
#' \code{adiab_cor}{ = 1, use adiabatic lapse rate correction? 1=yes, 0=no}\cr\cr
#' \code{air_temperature_offset}{ = 0, warming offset vector, °C (negative values mean cooling). Can supply a single value or a vector the length of the number of days to be simulated.}\cr\cr
#' \code{spatial}{ = "c:/Australian Environment/", choose location of terrain data}\cr\cr
#' \code{vlsci}{ = 0, running on the VLSCI system? 1=yes, 0=no}\cr\cr
#' \code{soilgrids}{ = 0, query soilgrids.org database for soil hydraulic properties?}\cr\cr
#' \code{message}{ = 0, allow the Fortran integrator to output warnings? (1) or not (0)}\cr\cr
#' \code{fail}{ = nyears x 24 x 365, how many restarts of the integrator before the Fortran program quits (avoids endless loops when solutions can't be found)}\cr\cr
#' \code{runmicro}{ = 1, call the microclimate model (1) or not (0), if you just want the downscaled input weather data}\cr\cr
#'
#' \strong{ General additional parameters:}\cr\cr
#'
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
#' \code{rainoff}{ = 0, Rain offset (mm), used to induce changes in rainfall from NIWA values. Can be a single value or a vector matching the number of days to simulate. If negative values are used, rainfall will be prevented from becomming negative.}\cr\cr
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
#' \code{snowdens}{ = 0.375, snow density (mg/m3), overridden by densfun}\cr\cr
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
#' \code{tides}
#' { is used to specify tide presence, sea water temperature and presence of wavesplash}\cr\cr
#' \code{tides}{ = matrix(data = 0., nrow = 24*365*nyears, ncol = 3), matrix of 1. tide state (0=out, 1=in), 2. Water temperature (°C) and 3. Wave splash (0=yes, 1=no)}\cr\cr
#' }
#'
#' \strong{Outputs:}
#'
#' \code{ndays}{ - number of days for which predictions are made}\cr\cr
#' \code{longlat}{ - longitude and latitude for which simulation was run (decimal degrees)}\cr\cr
#' \code{dates}{ - vector of dates (hourly, POSIXct, timezone = NZ)}\cr\cr
#' \code{dates2}{ - vector of dates (daily, POSIXct, timezone = NZ)}\cr\cr
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
#'micro<-micro_nz() # run the model with default location (Dunedin) and settings
#'
#'micromet_lowshade<-as.data.frame(micro$micromet_lowshade) # above ground microclimatic conditions, min shade
#'micromet_highshade<-as.data.frame(micro$micromet_highshade) # above ground microclimatic conditions, max shade
#'soil_temperature_lowshade<-as.data.frame(micro$soil_temperature_lowshade) # soil temperatures, minimum shade
#'soil_temperature_highshade<-as.data.frame(micro$soil_temperature_highshade) # soil temperatures, maximum shade
#'
#'# append dates
#'days<-rep(seq(1,12),24)
#'days<-days[order(days)]
#'dates<-days+micromet_lowshade$TIME/60/24-1 # dates for hourly output
#'dates2<-seq(1,12,1) # dates for daily output
#'
#'plotmetout<-cbind(dates,micromet_lowshade)
#'plotsoil<-cbind(dates,soil_temperature_lowshade)
#'plotshadmet<-cbind(dates,micromet_highshade)
#'plotshadsoil<-cbind(dates,soil_temperature_highshade)
#'
#'minimum_shade<-micro$minimum_shade[1]
#'maximum_shade<-micro$maximum_shade[1]
#'
#'# plotting above-ground conditions in minimum shade
#'with(plotmetout,{plot(air_temperature_local ~ dates,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l",main=paste("air temperature, ",minimum_shade,"% shade",sep=""))})
#'with(plotmetout,{points(air_temperature_reference ~ dates,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l",lty=2,col='blue')})
#'with(plotmetout,{plot(relative_humidity_local ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",ylim=c(0,100),main=paste("humidity, ",minimum_shade,"% shade",sep=""))})
#'with(plotmetout,{points(relative_humidity_reference ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",col='blue',lty=2,ylim=c(0,100))})
#'with(plotmetout,{plot(sky_temperature ~ dates,xlab = "Date and Time", ylab = "Sky Temperature (°C)"
#',  type = "l",main=paste("sky temperature, ",minimum_shade,"% shade",sep=""))})
#'with(plotmetout,{plot(wind_speed_reference ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#',  type = "l",main="wind speed",ylim = c(0, 15))})
#'with(plotmetout,{points(wind_speed_local ~ dates,xlab = "Date and Time", ylab = "Wind Speed (m/s)"
#',  type = "l",lty=2,col='blue')})
#'with(plotmetout,{plot(zenith_angle ~ dates,xlab = "Date and Time", ylab = "Zenith Angle of Sun (deg)"
#',  type = "l",main="solar angle, sun")})
#'with(plotmetout,{plot(solar_radiation ~ dates,xlab = "Date and Time", ylab = "Solar Radiation (W/m2)"
#',  type = "l",main="solar radiation")})
#'
#'# plotting soil temperature for minimum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(plotsoil[,i+3]~plotsoil[,1],xlab = "Date and Time", ylab = "Soil Temperature (°C)"
#'    ,col=i,type = "l",main=paste("soil temperature ",minimum_shade,"% shade",sep=""))
#'  }else{
#'    points(plotsoil[,i+3]~plotsoil[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)",col=i,type = "l")
#'  }
#'}
#'
#'# plotting above-ground conditions in maximum shade
#'with(plotshadmet,{plot(air_temperature_local ~ dates,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l",main="air temperature, sun")})
#'with(plotshadmet,{points(air_temperature_reference ~ dates,xlab = "Date and Time", ylab = "Air Temperature (°C)"
#', type = "l",lty=2,col='blue')})
#'with(plotshadmet,{plot(relative_humidity_local ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",ylim=c(0,100),main="humidity, shade")})
#'with(plotshadmet,{points(relative_humidity_reference ~ dates,xlab = "Date and Time", ylab = "Relative Humidity (%)"
#', type = "l",col='blue',lty=2,ylim=c(0,100))})
#'with(plotshadmet,{plot(sky_temperature ~ dates,xlab = "Date and Time", ylab = "Sky Temperature (°C)",
#'  type = "l",main="sky temperature, shade")})
#'
#'# plotting soil temperature for maximum shade
#'for(i in 1:10){
#'  if(i==1){
#'    plot(plotshadsoil[,i+3]~plotshadsoil[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)",col=i,type = "l",main=paste("soil temperature ",maximum_shade,"% shade",sep=""))
#'  }else{
#'    points(plotshadsoil[,i+3]~plotshadsoil[,1],xlab = "Date and Time", ylab = "Soil Temperature
#'     (°C)",col=i,type = "l")
#'  }
#'}
micro_nz <- function(
  loc = c(170.50280, -45.87876),
  ystart = 2000,
  yfinish = 2000,
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
  terrain = 0,
  dailywind = 1,
  wind_multiplier = 1,
  adiab_cor = 1,
  air_temperature_offset = 0,
  spatial = "C:/Spatial_Data/Climate/New Zealand/weather",
  vlsci = 0,
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
  Clay = 20,
  SatWater = rep(0.26, 10),
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
  snowmelt = 0.9,
  undercatch = 1,
  rainmelt = 0.0125,
  shore = 0,
  tides = 0,
  scenario = "",
  year="",
  barcoo="",
  quadrangle = 1,
  hourly = 0,
  rainhourly = 0,
  rainhour = 0,
  rainoff = 0,
  radiation_per_wavelength = 0,
  scattered_uv = 0,
  max_iterations_per_day = 3,
  soilgrids = 0,
  longwave_radiation_model = 0,
  forecast = 0,
  message = 0,
  fail = nyears * 24 * 365,
  runmicro = 1,
  snowcond = 0,
  intercept = max(maximum_shade) / 100 * 0.3,
  grasshade = 0,
  maxsurf = 85) {

  # Reference height (m) at which input air temperature, wind speed and
  # relative humidity are measured
  reference_height <- 1.2

  # --- 1. Input validation --------------------------------------------------
  # Dataset-specific checks first
  errors <- 0
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

  if(errors==0){ # continue

    ################## time related variables #################################
    nyears <- yfinish-ystart+1
    yearlist <- seq(ystart, (ystart + (nyears - 1)), 1)
    tzone <- paste("Etc/GMT+",10,sep="")
    dates <- seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days")
    ndays <- length(dates)
    doys12<-c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month
    microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial max_iterations_per_day day burn in
    daystart<-1
    # Expand scalar or short shade inputs to one value per simulation day
    shades <- setup_shade_vectors(minimum_shade, maximum_shade, ndays)
    minimum_shade_daily <- shades$minimum_shade_daily
    maximum_shade_daily <- shades$maximum_shade_daily
    leapyears<-seq(1972,2060,4)
    for(mm in 1:nyears){
      if(mm == 1){
        currenty <- ystart
      }else{
        currenty <- ystart + mm
      }
      if(currenty %in% leapyears){
        dayoy <- seq(1,366)
      }else{
        dayoy <- seq(1,365)
      }
      if(mm == 1){
        day_of_year <- dayoy
      }else{
        day_of_year <- c(day_of_year, dayoy)
      }
    }
    start_day <- 1 # start day
    end_day<-ndays # end day
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

    longlat <- loc
    x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))

    library(terra)
    library(proj4)
    library(ncdf4)

    cat("running micro_global to get clear sky solar \n")
    if(global_aerosol_database == 0){
      aerosol_optical_depth <- c(0.0670358341290886,0.0662612704779235,0.065497075238002,0.0647431301168489,0.0639993178022531,0.0632655219571553,0.0625416272145492,0.0611230843885423,0.0597427855962549,0.0583998423063099,0.0570933810229656,0.0558225431259535,0.0545864847111214,0.0533843764318805,0.0522154033414562,0.0499736739981675,0.047855059159556,0.0458535417401334,0.0439633201842001,0.0421788036108921,0.0404946070106968,0.0389055464934382,0.0374066345877315,0.0359930755919066,0.0346602609764008,0.0334037648376212,0.0322193394032758,0.0311029105891739,0.0300505736074963,0.0290585886265337,0.0281233764818952,0.0272415144391857,0.0264097320081524,0.0256249068083005,0.0248840604859789,0.0241843546829336,0.0235230870563317,0.0228976873502544,0.0223057135186581,0.0217448478998064,0.0212128934421699,0.0207077699817964,0.0202275105711489,0.0197702578594144,0.0193342605242809,0.0189178697551836,0.0177713140039894,0.0174187914242432,0.0170790495503944,0.0167509836728154,0.0164335684174899,0.0161258546410128,0.0158269663770596,0.0155360978343254,0.0152525104459325,0.0149755299703076,0.0147045436435285,0.0144389973831391,0.0141783930434343,0.0134220329447663,0.0131772403830191,0.0129356456025128,0.0126970313213065,0.0124612184223418,0.0122280636204822,0.01199745718102,0.0115436048739351,0.0110993711778668,0.0108808815754663,0.0106648652077878,0.0104513876347606,0.0102405315676965,0.00982708969547694,0.00962473896278535,0.00903679230300494,0.00884767454432418,0.0083031278398166,0.00796072474935954,0.00755817587626185,0.00718610751850881,0.00704629977586921,0.00684663903049612,0.00654155580333479,0.00642947339729728,0.00627223096874308,0.00603955966866779,0.00580920937536261,0.00568506186880564,0.00563167068287251,0.00556222005081865,0.00550522989971023,0.00547395763028062,0.0054478983436216,0.00541823364504573,0.00539532163908382,0.00539239864119488,0.00541690124712384,0.00551525885358836,0.00564825853509463,0.00577220185074264,0.00584222986640171,0.00581645238345584,0.00566088137411449,0.00535516862329704,0.00489914757707667,0.00432017939770409,0.0036813032251836,0.00309019064543606,0.00270890436501562,0.00276446109239711,0.00356019862584603)
    }else{
      aerosol_optical_depth <- 0
    }
    micro_clearsky <- micro_global(loc = c(x[1], x[2]), clearsky = 1, aerosol_optical_depth = aerosol_optical_depth, timeinterval = 365, solar_model_only = 1)
    clearskyrad <- micro_clearsky$micromet_lowshade[,c(1, 13)]
    clearskysum <- aggregate(clearskyrad[,2], by = list(clearskyrad[,1]), FUN = sum)[,2]

    # Get reference longitude for solar noon correction
    solar_noon_longitude <- get_timezone_alref(lon = x[1], lat = x[2], timezone = timezone)
    hemisphere <- ifelse(x[2]<0,2.,1.) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    latitude_degrees <- abs(trunc(x[2]))
    latitude_minutes <- (abs(x[2])-latitude_degrees)*60
    longitude_degrees <- abs(trunc(x[1]))
    longitude_minutes <- (abs(x[1])-longitude_degrees)*60
    azmuth<-aspect


    soilprop<-cbind(0,0)

    r1<-terra::rast(paste(spatial,'/nz_geo3_km.asc',sep=""))
    NZDEM<-as.numeric(terra::extract(r1,x))*1000

    if(is.na(elevation) == FALSE){ # check if user-specified elevation
      ALTITUDES <- elevation
    }else{
      utm<-project(as.matrix(x), "+proj=tmerc +lat_0=0.0 +lon_0=173.0 +k=0.9996 +x_0=1600000.0 +y_0=10000000.0 +datum=WGS84 +units=m")
      nc<-nc_open(paste(spatial,"/elevslpasphori.nc",sep=""))
      easting<-ncvar_get(nc,"easting")
      northing<-ncvar_get(nc,"northing")
      dist1<-abs(easting-utm[1])
      index1<-which.min(dist1)
      dist2<-abs(northing-utm[2])
      index2<-which.min(dist2)
      start<-c(index1,index2,1)
      count<-c(1,1,-1)
      elevslpasphori<-as.numeric(ncvar_get(nc,varid="variable",start=start,count=count))
      nc_close(nc)

      ALTITUDES <- elevslpasphori[1]
      if(is.na(ALTITUDES)==TRUE){ALTITUDES<-NZDEM}
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
      VIEWF_all <- 1-rowSums(sin(t(HORIZONS)*pi/180))/length(t(HORIZONS)) # convert horizon angles to radians and calc view factor(s)
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

    if(soilgrids == 1){
      sg <- fetch_soilgrids(x, depths)
      if (!is.null(sg)) { air_entry_water_potential <- sg$air_entry_water_potential; saturated_hydraulic_conductivity <- sg$saturated_hydraulic_conductivity; campbell_b_parameter <- sg$campbell_b_parameter; soil_bulk_density <- sg$soil_bulk_density; bulk_density <- sg$bulk_density }
    }

    delta_elev = NZDEM - ALTITUDES
    adiab_corr_max <- delta_elev * maximum_temperature_lapse_rate # Adiabatic temperature correction for elevation (C)
    adiab_corr_min <- delta_elev * minimum_temperature_lapse_rate # Adiabatic temperature correction for elevation (C)

    if(scenario!=""){
      cat("generate climate change scenario", '\n')
      # diff spline function
      getdiff<-function(diffs,grid){
        leapyears<-seq(1972,2060,4)
        yearstodo<-seq(ystart,ystart*nyears)
        diff1<-(unlist(diffs[1])+unlist(diffs[12]))/2

        # generate list of days
        for(ys in 1:nyears){
          if(yearstodo[ys]%in%leapyears){
            day<-c(1,15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 366)
          }else{
            day<-c(1,15.5, 45, 74.5, 105, 135.5, 166, 196.5, 227.5, 258, 288.5, 319, 349.5, 365)
          }
          if(ys==1){
            days2=day
            days=day
          }else{
            if(yearstodo[ys]%in%leapyears){
              days2=c(days2,(day+366*(ys-1)))
            }else{
              days2=c(days2,(day+365*(ys-1)))
            }
            days=c(days,day)
          }
        }

        if(is.na(diffs[1])==TRUE){
          # find the nearest cell with data
          NArem<-grid[[1]]
          NArem<-Which(!is.na(NArem), cells=TRUE)
          dist<-distanceFromPoints(maxTst05[[1]],x)
          distNA<-as.numeric(terra::extract(dist,NArem))
          cellsR<-cbind(distNA,NArem)
          distmin<-which.min(distNA)
          cellrep<-cellsR[distmin,2]
          diffs<-as.numeric(terra::extract(maxTst05,cellrep))
          diff1<-(unlist(diffs[1])+unlist(diffs[12]))/2
        }
        diffs3=rep(c(diff1,diffs,diff1),nyears)
        days_diffs<-data.frame(matrix(NA, nrow = nyears*14, ncol = 3))
        days_diffs[,1]<-days
        days_diffs[,3]<-days2
        days_diffs[,2]<-diffs3
        colnames(days_diffs)<-c("days","diffs","new_day")

        # interpolate monthly differences
        f<-approxfun(x=days_diffs$new_day, y=days_diffs$diffs)
        xx<-seq(1,max(days2),1)
        sp_diff<-f(xx)
        return(sp_diff)
      }

      ########### Max and Min Air Temps ################

      TMEAN<-rast(paste(spatial,"/CC/TMEAN_",year,"_",scenario,".nc",sep="")) # air temperature shift

      diffs<-as.numeric(terra::extract(TMEAN,x))
      TMAXX_diff <- getdiff(diffs,TMEAN)
      TMINN_diff <- TMAXX_diff

      ################ VP ############################

      # We are given a value of % change in absolute humidity per change in degree C of
      # air temperature, for each season. So we need to extract these values, spline them
      # to monthly (making sure we slide the days back by 30 days for the interpolation
      # because the values are centred within each season). Then we get the

      AH<-rast(paste0(spatial,"/CC/AHCC_SUM.nc"),paste0(spatial,"/CC/AHCC_AUT.nc"),paste0(spatial,"/CC/AHCC_WIN.nc"),paste0(spatial,"/CC/AHCC_SPR.nc"))

      diffs<-as.numeric(terra::extract(AH,x)) # extract seasonal values
      diff1<-(diffs[1]+diffs[4])/2 # get mean of first and last
      diffs3=c(diff1,diffs,diff1) # make vector of 6, with the start and end being the mean of the first and last
      day<-c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349) # middle day of each month
      day0<-c(1,45,135,225,315,365) # days for seasonal interpolation
      f<-approxfun(x=day0, y=diffs3) # use approxfun for the linear interpolation
      sp_diff<-f(seq(1,365)) # spline across year
      sp_diff<-cbind(seq(1,365)-30,sp_diff) # add days of year offset by 30 days
      sp_diff[1:30,1]<-336:365 # replace the first 30 days with the last 30 days of year
      sp_diff<-sp_diff[order(sp_diff[,1]),] # reorder per ay
      diffs<-sp_diff[day,2] # extract the value for the middle day of each month
      VP_diff<-getdiff(diffs,AH) # now use the 'getdiffs' function to spline 365 days for however many years

      # change diff to a proportion, multiply by change in air temp to get overall
      # proportional change, and add to one to get the multiplier
      VP_diff2<-1+(VP_diff/100)*TMAXX_diff

      ################ wind ############################

      WIND_diff<-1

      ############# SOLAR/CLOUD COVER ##################

      SOLAR_diff<-1
    }

    cat("extracting weather data \n")
    # read daily weather
    yearlist<-seq(ystart,(ystart+(nyears-1)),1)
    for(j in 1:nyears){ # start loop through years
      cat(paste('reading weather input for ',yearlist[j],' \n',sep=""))
      lon_1<-as.numeric(longlat[1])
      lat_1<-as.numeric(longlat[2])
      lat<-read.csv(paste0(spatial,'/ncdf_lat.csv'))[,2]
      lon<-read.csv(paste0(spatial,'/ncdf_lon.csv'))[,2]
      dist1<-abs(lon-lon_1)
      index1<-which.min(dist1)
      dist2<-abs(lat-lat_1)
      index2<-which.min(dist2)
      start<-c(index1,index2,1)
      count<-c(1,1,-1)
      if(j==1){
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_Tmax.nc',sep=""))
        Tmax<-as.numeric(ncvar_get(nc,varid="variable",start=start,count))
        nc_close(nc)
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_Tmin.nc',sep=""))
        Tmin<-as.numeric(ncvar_get(nc,varid="variable",start=start,count))
        nc_close(nc)
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_VP.nc',sep=""))
        VP<-as.numeric(ncvar_get(nc,varid="variable",start=start,count))
        nc_close(nc)
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_Rain.nc',sep=""))
        Rain<-as.numeric(ncvar_get(nc,varid="variable",start=start,count))
        nc_close(nc)
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_SoilM.nc',sep=""))
        SoilM<-as.numeric(ncvar_get(nc,varid="variable",start=start,count))
        nc_close(nc)
        if(yearlist[j]>=1997){
          nc<-nc_open(paste(spatial,"/",yearlist[j],'_Wind.nc',sep=""))
          Wind<-as.numeric(ncvar_get(nc,varid="variable",start=start,count))
          nc_close(nc)
        }else{
          Wind<-Tmax*0+3
        }
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_Rad.nc',sep=""))
        Rad<-as.numeric(ncvar_get(nc,varid="variable",start=start,count))
        nc_close(nc)
        if(length(Rad)==366){# add day for leap year if needed
          clear<-c(clearskysum[1:59],clearskysum[59],clearskysum[60:365])
        }else{
          clear<-clearskysum
        }
        clear <-  clear * 3600 / 1e6 # to MJ/d
        cloud<-(1-Rad/clear)*100
        cloud[cloud<0]<-0
        cloud[cloud>100]<-100
        cloud_max<-as.numeric(cloud)
      }else{
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_Tmax.nc',sep=""))
        Tmax<-as.numeric(c(Tmax,ncvar_get(nc,varid="variable",start=start,count)))
        nc_close(nc)
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_Tmin.nc',sep=""))
        Tmin<-as.numeric(c(Tmin,ncvar_get(nc,varid="variable",start=start,count)))
        nc_close(nc)
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_VP.nc',sep=""))
        VP<-as.numeric(c(VP,ncvar_get(nc,varid="variable",start=start,count)))
        nc_close(nc)
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_Rain.nc',sep=""))
        Rain<-as.numeric(c(Rain,ncvar_get(nc,varid="variable",start=start,count)))
        nc_close(nc)
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_SoilM.nc',sep=""))
        SoilM<-as.numeric(c(SoilM,ncvar_get(nc,varid="variable",start=start,count)))
        nc_close(nc)
        if(yearlist[j]>=1997){
          nc<-nc_open(paste(spatial,"/",yearlist[j],'_Wind.nc',sep=""))
          Wind<-as.numeric(c(Wind,ncvar_get(nc,varid="variable",start=start,count)))
          nc_close(nc)
        }else{
          Wind<-Tmax*0+3
        }
        nc<-nc_open(paste(spatial,"/",yearlist[j],'_Rad.nc',sep=""))
        Rad<-as.numeric(ncvar_get(nc,varid="variable",start=start,count))
        nc_close(nc)
        if(length(Rad)==366){# add day for leap year if needed
          clear<-c(clearskysum[1:59],clearskysum[59],clearskysum[60:365])
        }else{
          clear<-clearskysum
        }
        cloud<-(1-Rad/clear)*100
        cloud[cloud<0]<-0
        cloud[cloud>100]<-100
        cloud_max<-as.numeric(c(cloud_max,cloud))
      }
    }
    if(clearsky == 1){
      cloud_max <- cloud_max * 0
    }
    cloud_min<-cloud_max
    Wind<-Wind*wind_multiplier
    Wind[Wind==0]<-0.1

    end_day<-ndays
    start_day <- 1 # start month
    # end preliminary test for incomplete year, if simulation includes the present year

    tzone<-paste("Etc/GMT-12",sep="") # doing it this way ignores daylight savings!
    ndays<-length(seq(ISOdate(ystart,1,1,tz=tzone)-3600*12, ISOdate((ystart+nyears),1,1,tz=tzone)-3600*13, by="days"))

    # end preliminary test for incomplete year, if simulation includes the present year

    if(is.na(ALTITUDES)!=TRUE){

      aerosol_optical_depth <- compute_tai(longlat, global_aerosol_database, TAI_AUSTRALIA)

      if(adiab_cor==1){
        reference_temperature_max<-as.matrix(Tmax+adiab_corr_max)
        reference_temperature_min<-as.matrix(Tmin+adiab_corr_min)
      }else{
        reference_temperature_max<-as.matrix(Tmax)
        reference_temperature_min<-as.matrix(Tmin)
      }
      if(scenario!=""){
        reference_temperature_max=reference_temperature_max+TMAXX_diff
        reference_temperature_min=reference_temperature_min+TMINN_diff
        VP<-VP*VP_diff2 # modify the predicted VP by this factor
      }
      rainfall<-Rain+rainoff
      rainfall[rainfall < 0] <- 0

      if(scenario!=""){
        RAIN_current<-as.data.frame(rainfall)
        dates<-seq(ISOdate(ystart,1,1,tz=paste("Etc/GMT-",10,sep=""))-3600*12, ISOdate((ystart+nyears),1,1,tz=paste("Etc/GMT-",10,sep=""))-3600*13, by="days")
        RAINFALL_sum<-aggregate(RAIN_current, by=list(format(dates,"%m-%Y")),FUN=sum)
        dates2<-RAINFALL_sum$Group.1
        RAINFALL_sum<-RAINFALL_sum[order(as.Date(paste("01-",RAINFALL_sum$Group.1,sep=""),"%m-%Y")),2]

        RAIN<-rast(paste(spatial,"/CC/RAIN_",year,"_",scenario,".nc",sep="")) # rainfall shift
        diffs<-rep(as.numeric(terra::extract(RAIN,x)),nyears)

        rainfall_new<-(RAINFALL_sum+RAINFALL_sum*(diffs/100))

        rainfall_new[rainfall_new < 0 ]= 0 # get rid of any negative rainfall values

        #########Now get predicted change in rainfall (could also get this from OzClim or ClimSim layer)#############
        Diff_prop<-rainfall_new/RAINFALL_sum # proportion change
        Diff_prop[Diff_prop=='NaN']= 0
        Diff_prop[Diff_prop=='Inf']= 0 ## If was previously no rainfall and now is rainfall need to alter code so this is simply added

        newRAINFALL=rep(NA,length(rainfall))
        for (k in 1:length(rainfall)){ # loop through each sites applying monthly % changes
          month<-which(dates2==format(dates[k],"%m-%Y"))
          newRAINFALL[k]<-rainfall[k]*Diff_prop[month]
        } # end of loop through each day
        newRAINFALL[newRAINFALL<0.1]<-0
        newRAINFALL[is.na(newRAINFALL)]<-0
        rainfall=newRAINFALL
      }

      VAPRES<-VP*100 # convert from hectopascals to pascals
      TMAXK<-reference_temperature_max+273.15
      loge<-TMAXK
      loge2<-loge
      loge[loge2>273.16]<- -7.90298*(373.16/TMAXK[loge2>273.16]-1.)+5.02808*log10(373.16/TMAXK[loge2>273.16])-1.3816E-07*(10.^(11.344*(1.-TMAXK[loge2>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/TMAXK[loge2>273.16]-1.))-1.)+log10(1013.246)
      loge[loge2<=273.16]<- -9.09718*(273.16/TMAXK[loge2<=273.16]-1.)-3.56654*log10(273.16/TMAXK[loge2<=273.16])+.876793*(1.-TMAXK[loge2<=273.16]/273.16)+log10(6.1071)
      estar<-(10.^loge)*100.
      reference_humidity_min<-(VAPRES/estar)*100
      reference_humidity_min[reference_humidity_min>100]<-100
      reference_humidity_min[reference_humidity_min<0]<-0.01
      #reference_humidity_min
      TMINK<-reference_temperature_min+273.15
      loge<-TMINK
      loge2<-loge
      loge[loge2>273.16]<- -7.90298*(373.16/TMAXK[loge2>273.16]-1.)+5.02808*log10(373.16/TMAXK[loge2>273.16])-1.3816E-07*(10.^(11.344*(1.-TMAXK[loge2>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/TMAXK[loge2>273.16]-1.))-1.)+log10(1013.246)
      loge[loge2<=273.16]<- -9.09718*(273.16/TMAXK[loge2<=273.16]-1.)-3.56654*log10(273.16/TMAXK[loge2<=273.16])+.876793*(1.-TMAXK[loge2<=273.16]/273.16)+log10(6.1071)
      estar<-(10.^loge)*100.
      reference_humidity_max<-(VAPRES/estar)*100
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

      Intrvls<-rep(0,ndays)
      Intrvls[1] <- 1 # user-supplied last day-of-year in each time interval sequence
      Numtyps <- 1 # number of substrate types
      Numint <- 1  # number of time intervals
      soil_nodes <- matrix(data = 0, nrow = 10, ncol = ndays) # deepest nodes for each substrate type
      soil_nodes[1,1] <- 10. # deepest nodes for each substrate type
      solar_noon_longitude <- abs(trunc(x[1]))
      hemisphere <- ifelse(x[2]<0,2.,1.)
      latitude_degrees <- abs(trunc(x[2]))
      latitude_minutes <- (abs(x[2])-latitude_degrees)*60
      longitude_degrees <- abs(trunc(x[1]))
      longitude_minutes <- (abs(x[1])-longitude_degrees)*60
      if(adiab_cor==1){
        elevation<-ALTITUDES
      }else{
        elevation<-NZDEM
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
        deep_soil_temperature<-rep(avetemp,365)
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
      reference_wind_max<-reference_wind_max*(1.2/2)^0.15
      reference_wind_min<-reference_wind_min*(1.2/2)^0.15


      # impose warming
      reference_temperature_max<-reference_temperature_max+air_temperature_offset
      reference_temperature_min<-reference_temperature_min+air_temperature_offset


      surface_emissivity<-matrix(nrow=ndays,data=0)
      surface_emissivity<-surface_emissivity+surface_emissivity

      moists2<-matrix(nrow=10, ncol = ndays, data=0)
      moists2[1,ndays]<-0.2
      soil_moisture_profile<-moists2

      if(soil_moisture_model==1){
        moists2<-matrix(nrow=10, ncol = ndays, data=0) # set up an empty vector for soil moisture values through time
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

      # --- 3. Prepare model inputs ------------------------------------------
      # Ensure numeric types for Fortran
      elevation  <- as.numeric(elevation)
      solar_noon_longitude <- as.numeric(solar_noon_longitude)

      # Hourly input override vectors — unused here (daily VCSN dataset)
      if (hourly == 0) {
        reference_temperature <- rep(0,  24 * ndays)
        reference_humidity   <- rep(0,  24 * ndays)
        reference_wind_speed   <- rep(0,  24 * ndays)
        cloud_cover  <- rep(0,  24 * ndays)
        global_radiation <- rep(0,  24 * ndays)
        zenith_angle_hourly  <- rep(-1, 24 * ndays)
        longwave_radiation  <- rep(-1, 24 * ndays)
      }
      if (rainhourly == 0) {
        rainfall_hourly <- rep(0, 24 * ndays)
      } else {
        rainfall_hourly <- rainhour
      }

      # Build the flat numeric vector consumed by the Fortran microclimate solver
      micro_input <- build_microinput(
        ndays = ndays, roughness_height = roughness_height, tolerance = tolerance, local_height = local_height, reference_height = reference_height,
        Numtyps = Numtyps, roughness_height_1 = roughness_height_1, roughness_height_2 = roughness_height_2, wind_profile_height_1 = wind_profile_height_1, wind_profile_height_2 = wind_profile_height_2,
        start_day = start_day, end_day = end_day, hemisphere = hemisphere, latitude_degrees = latitude_degrees, latitude_minutes = latitude_minutes,
        longitude_degrees = longitude_degrees, longitude_minutes = longitude_minutes, solar_noon_longitude = solar_noon_longitude, slope = slope,
        azmuth = azmuth, elevation = elevation, precipitable_water = precipitable_water, microdaily = microdaily,
        mean_annual_temperature = mean_annual_temperature, orbital_eccentricity = orbital_eccentricity, VIEWF = VIEWF, snowtemp = snowtemp,
        snowdens = snowdens, snowmelt = snowmelt, undercatch = undercatch,
        rain_multiplier = rain_multiplier, runshade = runshade, soil_moisture_model = soil_moisture_model,
        maximum_pooling_depth = maximum_pooling_depth, evenrain = evenrain, snow_model = snow_model,
        rainmelt = rainmelt, output_to_csv = output_to_csv, densfun = densfun,
        hourly = hourly, rainhourly = rainhourly, radiation_per_wavelength = radiation_per_wavelength, scattered_uv = scattered_uv,
        root_resistance = root_resistance, stomatal_closure_potential = stomatal_closure_potential, leaf_resistance = leaf_resistance, stomatal_stability_parameter = stomatal_stability_parameter, root_radius = root_radius, moist_error = moist_error,
        moist_count = moist_count, longwave_radiation_model = longwave_radiation_model, message = message, fail = fail,
        snowcond = snowcond, intercept = intercept, grasshade = grasshade,
        solar_model_only = solar_model_only, canopy_roughness_height = canopy_roughness_height, zero_plane_displacement = zero_plane_displacement, maxima_times = maxima_times, minima_times = minima_times,
        spinup = spinup, maxsurf = maxsurf, max_iterations_per_day = max_iterations_per_day
      )

      if (length(leaf_area_index) < ndays) {
        leaf_area_index <- rep(leaf_area_index[1], ndays)
      }
      if (shore == 0) {
        tides <- matrix(data = 0, nrow = 24 * ndays, ncol = 3)
      }

      # Assemble the complete input list for microclimate()
      micro <- build_micro_list(
        micro_input = micro_input, day_of_year = day_of_year, surface_emissivity = surface_emissivity, depths = depths,
        soil_nodes = soil_nodes, maximum_shade_daily = maximum_shade_daily, minimum_shade_daily = minimum_shade_daily,
        reference_temperature_max = reference_temperature_max, reference_temperature_min = reference_temperature_min, reference_humidity_max = reference_humidity_max, reference_humidity_min = reference_humidity_min,
        cloud_max = cloud_max, cloud_min = cloud_min, reference_wind_max = reference_wind_max, reference_wind_min = reference_wind_min,
        reference_temperature = reference_temperature, reference_humidity = reference_humidity, reference_wind_speed = reference_wind_speed, cloud_cover = cloud_cover,
        global_radiation = global_radiation, rainfall_hourly = rainfall_hourly, zenith_angle_hourly = zenith_angle_hourly, longwave_radiation = longwave_radiation,
        albedo = albedo, soil_wetness = soil_wetness, soilinit = soilinit, horizon_angles = horizon_angles,
        aerosol_optical_depth = aerosol_optical_depth, soil_properties = soil_properties, soil_moisture_profile = soil_moisture_profile, rainfall = rainfall,
        deep_soil_temperature = deep_soil_temperature, air_entry_water_potential = air_entry_water_potential, saturated_hydraulic_conductivity = saturated_hydraulic_conductivity, campbell_b_parameter = campbell_b_parameter, soil_bulk_density = soil_bulk_density, soil_mineral_density = soil_mineral_density,
        root_density = root_density, leaf_area_index = leaf_area_index, tides = tides
      )

      # Optionally dump all inputs to CSV files for debugging
      if (write_input == 1) {
        write_micro_csv(
          micro_input = micro_input, day_of_year = day_of_year, surface_emissivity = surface_emissivity, depths = depths,
          soil_nodes = soil_nodes, maximum_shade_daily = maximum_shade_daily, minimum_shade_daily = minimum_shade_daily,
          maxima_times = maxima_times, minima_times = minima_times,
          reference_temperature_max = reference_temperature_max, reference_temperature_min = reference_temperature_min, reference_humidity_max = reference_humidity_max, reference_humidity_min = reference_humidity_min,
          cloud_max = cloud_max, cloud_min = cloud_min, reference_wind_max = reference_wind_max, reference_wind_min = reference_wind_min,
          albedo = albedo, soil_wetness = soil_wetness, soilinit = soilinit, horizon_angles = horizon_angles,
          aerosol_optical_depth = aerosol_optical_depth, soil_properties = soil_properties, soil_moisture_profile = soil_moisture_profile, rainfall = rainfall,
          deep_soil_temperature = deep_soil_temperature, air_entry_water_potential = air_entry_water_potential, soil_bulk_density = soil_bulk_density, soil_mineral_density = soil_mineral_density, campbell_b_parameter = campbell_b_parameter, saturated_hydraulic_conductivity = saturated_hydraulic_conductivity,
          root_density = root_density, leaf_area_index = leaf_area_index, tides = tides,
          reference_temperature = reference_temperature, reference_humidity = reference_humidity, reference_wind_speed = reference_wind_speed, cloud_cover = cloud_cover,
          global_radiation = global_radiation, rainfall_hourly = rainfall_hourly, zenith_angle_hourly = zenith_angle_hourly, longwave_radiation = longwave_radiation
        )
      }

      # --- 4. Run the microclimate model ------------------------------------
      if (is.numeric(loc[1])) {
        location <- paste("long", loc[1], "lat", loc[2])
      } else {
        location <- loc
      }

      if (runmicro) {
        message(paste('running microclimate model for', ndays, 'days between',
                      ystart, 'and', yfinish, 'at site', location))
        message('Note: the output column `SOLR` in micromet_lowshade and micromet_highshade is for unshaded horizontal plane solar radiation')
        ptm <- proc.time()
        microut <- microclimate(micro)
        message(paste0('runtime ', (proc.time() - ptm)[3], ' seconds'))

        if (max(microut$micromet_lowshade[, 1] == 0)) {
          message("ERROR: the model crashed - try a different error tolerance (tolerance) or a different spacing in depths")
        }

        # NZ-specific date sequences (NZ timezone)
        dates  <- seq(as.POSIXct(paste0("01/01/", ystart),  format = "%d/%m/%Y", tz = 'NZ'),
                      as.POSIXct(paste0("31/12/", yfinish), format = "%d/%m/%Y", tz = 'NZ'),
                      by = 'hours')
        dates2 <- seq(as.POSIXct(paste0("01/01/", ystart),  format = "%d/%m/%Y", tz = 'NZ'),
                      as.POSIXct(paste0("31/12/", yfinish), format = "%d/%m/%Y", tz = 'NZ'),
                      by = 'days')

        # --- 5. Process and return results ----------------------------------
        out <- process_micro_output(microut,
                                    soil_moisture_model  = soil_moisture_model,
                                    snow_model = snow_model,
                                    radiation_per_wavelength      = radiation_per_wavelength)

        return(build_micro_return(
          out          = out,
          rainfall     = rainfall,
          ndays        = ndays,
          elevation         = elevation,
          surface_reflectivity         = surface_reflectivity,
          longlat      = c(x[1], x[2]),
          nyears       = nyears,
          timeinterval = nyears * 365,  # NZ uses daily runs
          minimum_shade_daily    = minimum_shade_daily,
          maximum_shade_daily    = maximum_shade_daily,
          depths          = depths,
          dates        = dates,
          dates2       = dates2,
          air_entry_water_potential           = air_entry_water_potential,
          soil_bulk_density           = soil_bulk_density,
          soil_mineral_density           = soil_mineral_density,
          campbell_b_parameter           = campbell_b_parameter,
          saturated_hydraulic_conductivity           = saturated_hydraulic_conductivity,
          dem          = NA,
          diffuse_frac = NA,
          snow_model    = snow_model,
          radiation_per_wavelength         = radiation_per_wavelength
        ))

      } else {
        # runmicro = FALSE: return only the prepared weather inputs
        return(list(
          rainfall = rainfall, reference_temperature_max = reference_temperature_max, reference_temperature_min = reference_temperature_min,
          reference_humidity_max = reference_humidity_max, reference_humidity_min = reference_humidity_min, reference_wind_max = reference_wind_max, reference_wind_min = reference_wind_min,
          cloud_max = cloud_max, cloud_min = cloud_min, cloud_cover = cloud_cover, reference_wind_speed = reference_wind_speed,
          reference_temperature = reference_temperature, reference_humidity = reference_humidity, rainfall_hourly = rainfall_hourly, global_radiation = global_radiation,
          zenith_angle_hourly = zenith_angle_hourly, longwave_radiation = longwave_radiation, dates = dates2, dates2 = dates2,
          air_entry_water_potential = air_entry_water_potential, soil_bulk_density = soil_bulk_density, soil_mineral_density = soil_mineral_density, campbell_b_parameter = campbell_b_parameter, saturated_hydraulic_conductivity = saturated_hydraulic_conductivity
        ))
      }
    } # end of check for NA sites
  } # end if(errors == 0)
} # end of micro_nz function
