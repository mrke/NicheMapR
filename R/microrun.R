#' microclimate model
#'
#' R wrapper for Fortran binary of Niche Mapper microclimate model
#' @param micro A vector of input variables for the microclimate model
#' @return micromet_lowshade The above ground micrometeorological conditions under the minimum specified shade
#' @return micromet_highshade The above ground micrometeorological conditions under the maximum specified shade
#' @return soil Hourly predictions of the soil temperatures under the minimum specified shade
#' @return soil_temperature_highshade Hourly predictions of the soil temperatures under the maximum specified shade
#' @return soil_moisture_lowshade Hourly predictions of the soil moisture under the minimum specified shade
#' @return soil_moisture_highshade Hourly predictions of the soil moisture under the maximum specified shade
#' @return soil_water_potential_lowshade Hourly predictions of the soil water potential under the minimum specified shade
#' @return soil_water_potential_highshade Hourly predictions of the soil water potential under the maximum specified shade
#' @return humid Hourly predictions of the soil humidity under the minimum specified shade
#' @return soil_humidity_highshade Hourly predictions of the soil humidity under the maximum specified shade
#' @return plant Hourly predictions of plant variables under the minimum specified shade
#' @return plant_output_highshade Hourly predictions of plant variables under the maximum specified shade
#' @return snow_output_lowshade Hourly predictions of the snow temperature under the minimum specified shade
#' @return snow_output_highshade Hourly predictions of the snow temperature under the maximum specified shade
#' @return soil_conductivity_lowshade Hourly predictions of the soil thermal conductivity under the minimum specified shade
#' @return soil_conductivity_highshade Hourly predictions of the soil thermal conductivity under the maximum specified shade
#' @return soil_specific_heat_lowshade Hourly predictions of the soil specific heat capacity under the minimum specified shade
#' @return soil_specific_heat_highshade Hourly predictions of soil specific heat capacity under the maximum specified shade
#' @return soil_density_lowshade Hourly predictions of the soil density under the minimum specified shade
#' @return soil_density_highshade Hourly predictions of the soil density under the maximum specified shade
#' @useDynLib "NicheMapR"
#' @export
microclimate <- function(micro) {
  doynum<-micro$micro_input[1]
  errors <- 0
  if(length(micro$micro_input) != 75){
    message("ERROR: micro_input has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$minimum_shade_daily) != doynum){
    message("ERROR: minimum_shade_daily has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$maximum_shade_daily) != doynum){
    message("ERROR: maximum_shade_daily has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$surface_emissivity) != doynum){
    message("ERROR: surface_emissivity has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$depths) != 10){
    message("ERROR: depths has the wrong number of inputs \n")
    errors <- 1
  }
  if(nrow(micro$soil_moisture_profile) != 10 | ncol(micro$soil_moisture_profile) != doynum){
    message("ERROR: soil_moisture_profile has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$reference_humidity_min) != doynum){
    message("ERROR: reference_humidity_min has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$reference_humidity_max) != doynum){
    message("ERROR: reference_humidity_max has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$cloud_min) != doynum){
    message("ERROR: cloud_min has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$cloud_max) != doynum){
    message("ERROR: cloud_max has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$reference_wind_min) != doynum){
    message("ERROR: reference_wind_min has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$reference_wind_max) != doynum){
    message("ERROR: reference_wind_max has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$reference_temperature_min) != doynum){
    message("ERROR: reference_temperature_min has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$reference_temperature_max) != doynum){
    message("ERROR: reference_temperature_max has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$albedo) != doynum){
    message("ERROR: albedo has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$soil_wetness) != doynum){
    message("ERROR: soil_wetness has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$rainfall) != doynum){
    message("ERROR: rainfall has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$deep_soil_temperature) != doynum){
    message("ERROR: deep_soil_temperature has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$leaf_area_index) != doynum){
    message("ERROR: leaf_area_index has the wrong number of inputs \n")
    errors <- 1
  }
  if(nrow(micro$tides) != doynum * 24){
    message("ERROR: tides has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$reference_temperature) != doynum * 24){
    message("ERROR: reference_temperature has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$reference_humidity) != doynum * 24){
    message("ERROR: reference_humidity has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$reference_wind_speed) != doynum * 24){
    message("ERROR: reference_wind_speed has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$cloud_cover) != doynum * 24){
    message("ERROR: cloud_cover has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$global_radiation) != doynum * 24){
    message("ERROR: global_radiation has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$rainfall_hourly) != doynum * 24){
    message("ERROR: rainfall_hourly has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$zenith_angle_hourly) != doynum * 24){
    message("ERROR: zenith_angle_hourly has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$longwave_radiation) != doynum * 24){
    message("ERROR: longwave_radiation has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$air_entry_water_potential) != 19){
    message("ERROR: air_entry_water_potential has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$campbell_b_parameter) != 19){
    message("ERROR: campbell_b_parameter has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$saturated_hydraulic_conductivity) != 19){
    message("ERROR: saturated_hydraulic_conductivity has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$soil_bulk_density) != 19){
    message("ERROR: soil_bulk_density has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$soil_mineral_density) != 19){
    message("ERROR: soil_mineral_density has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$root_density) != 19){
    message("ERROR: root_density has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$horizon_angles) != 24){
    message("ERROR: horizon_angles has the wrong number of inputs \n")
    errors <- 1
  }
  if(errors == 0){
  a <- .Fortran("microclimate",
                as.integer(doynum),
                as.double(micro$micro_input),
                as.double(micro$day_of_year),
                as.double(micro$surface_emissivity),
                as.double(micro$depths),
                as.double(micro$maximum_shade_daily),
                as.double(micro$minimum_shade_daily),
                as.double(micro$soil_nodes),
                as.double(micro$reference_humidity_max),
                as.double(micro$reference_humidity_min),
                as.double(micro$cloud_max),
                as.double(micro$cloud_min),
                as.double(micro$reference_wind_max),
                as.double(micro$reference_wind_min),
                as.double(micro$reference_temperature_max),
                as.double(micro$reference_temperature_min),
                as.double(micro$albedo),
                as.double(micro$soil_wetness),
                as.double(micro$soilinit),
                as.double(micro$horizon_angles),
                as.double(micro$aerosol_optical_depth),
                as.double(micro$soil_properties),
                as.double(micro$soil_moisture_profile),
                as.double(micro$rainfall),
                as.double(micro$deep_soil_temperature),
                as.double(micro$tides),
                as.double(micro$air_entry_water_potential),
                as.double(micro$saturated_hydraulic_conductivity),
                as.double(micro$campbell_b_parameter),
                as.double(micro$soil_bulk_density),
                as.double(micro$soil_mineral_density),
                as.double(micro$root_density),
                as.double(micro$leaf_area_index),
                as.double(micro$reference_temperature),
                as.double(micro$reference_humidity),
                as.double(micro$reference_wind_speed),
                as.double(micro$cloud_cover),
                as.double(micro$global_radiation),
                as.double(micro$rainfall_hourly),
                as.double(micro$zenith_angle_hourly),
                as.double(micro$longwave_radiation),
                micromet_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 19),
                soil_temperature_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                micromet_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 19),
                soil_temperature_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_moisture_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_moisture_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_humidity_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_humidity_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_water_potential_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_water_potential_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                snow_output_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 11),
                snow_output_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 11),
                plant_output_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 14),
                plant_output_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 14),
                soil_conductivity_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_conductivity_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_specific_heat_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_specific_heat_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_density_lowshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soil_density_highshade=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                direct_solar_spectrum=matrix(data = 0, nrow = 24 * doynum, ncol = 113),
                rayleigh_solar_spectrum=matrix(data = 0, nrow = 24 * doynum, ncol = 113),
                diffuse_solar_spectrum=matrix(data = 0, nrow = 24 * doynum, ncol = 113), PACKAGE = "NicheMapR")
  }
  micromet_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 19)
  micromet_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 19)
  soil_temperature_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_temperature_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_moisture_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_moisture_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_humidity_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_humidity_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_water_potential_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_water_potential_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  snow_output_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 11)
  snow_output_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 11)
  plant_output_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 14)
  plant_output_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 14)
  soil_conductivity_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_conductivity_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_specific_heat_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_specific_heat_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_density_lowshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soil_density_highshade <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  direct_solar_spectrum <- matrix(data = 0, nrow = 24 * doynum, ncol = 113)
  rayleigh_solar_spectrum <- matrix(data = 0, nrow = 24 * doynum, ncol = 113)
  diffuse_solar_spectrum <- matrix(data = 0, nrow = 24 * doynum, ncol = 113)
  storage.mode(micromet_lowshade) <- "double"
  storage.mode(micromet_highshade) <- "double"
  storage.mode(soil_temperature_lowshade) <- "double"
  storage.mode(soil_temperature_highshade) <- "double"
  storage.mode(soil_moisture_lowshade) <- "double"
  storage.mode(soil_moisture_highshade) <- "double"
  storage.mode(soil_humidity_lowshade) <- "double"
  storage.mode(soil_humidity_highshade) <- "double"
  storage.mode(soil_water_potential_lowshade) <- "double"
  storage.mode(soil_water_potential_highshade) <- "double"
  storage.mode(snow_output_lowshade) <- "double"
  storage.mode(snow_output_highshade) <- "double"
  storage.mode(plant_output_lowshade) <- "double"
  storage.mode(plant_output_highshade) <- "double"
  storage.mode(soil_conductivity_lowshade) <- "double"
  storage.mode(soil_conductivity_highshade) <- "double"
  storage.mode(soil_specific_heat_lowshade) <- "double"
  storage.mode(soil_specific_heat_highshade) <- "double"
  storage.mode(soil_density_lowshade) <- "double"
  storage.mode(soil_density_highshade) <- "double"
  storage.mode(direct_solar_spectrum) <- "double"
  storage.mode(rayleigh_solar_spectrum) <- "double"
  storage.mode(diffuse_solar_spectrum) <- "double"
  micromet_lowshade <- a$micromet_lowshade
  micromet_highshade <- a$micromet_highshade
  soil_temperature_lowshade <- a$soil_temperature_lowshade
  soil_temperature_highshade <- a$soil_temperature_highshade
  soil_moisture_lowshade <- a$soil_moisture_lowshade
  soil_moisture_highshade <- a$soil_moisture_highshade
  soil_humidity_lowshade <- a$soil_humidity_lowshade
  soil_humidity_highshade <- a$soil_humidity_highshade
  soil_water_potential_lowshade <- a$soil_water_potential_lowshade
  soil_water_potential_highshade <- a$soil_water_potential_highshade
  snow_output_lowshade <- a$snow_output_lowshade
  snow_output_highshade <- a$snow_output_highshade
  plant_output_lowshade <- a$plant_output_lowshade
  plant_output_highshade <- a$plant_output_highshade
  soil_conductivity_lowshade <- a$soil_conductivity_lowshade
  soil_conductivity_highshade <- a$soil_conductivity_highshade
  soil_specific_heat_lowshade <- a$soil_specific_heat_lowshade
  soil_specific_heat_highshade <- a$soil_specific_heat_highshade
  soil_density_lowshade <- a$soil_density_lowshade
  soil_density_highshade <- a$soil_density_highshade
  direct_solar_spectrum <- a$direct_solar_spectrum
  rayleigh_solar_spectrum <- a$rayleigh_solar_spectrum
  diffuse_solar_spectrum <- a$diffuse_solar_spectrum
  micromet.names <- c("day_of_year", "time", "air_temperature_local", "air_temperature_reference", "relative_humidity_local", "relative_humidity_reference", "wind_speed_local", "wind_speed_reference", "snow_melt", "pool_depth", "soil_wetness", "zenith_angle", "solar_radiation", "sky_temperature", "dew", "frost", "snow_fall", "snow_depth", "snow_density")
  colnames(micromet_lowshade) <- micromet.names
  colnames(micromet_highshade) <- micromet.names
  depth.names <- c("day_of_year", "time", paste0("depth_", micro$depths, "cm"))
  colnames(soil_temperature_lowshade) <- depth.names
  colnames(soil_temperature_highshade) <- depth.names
  colnames(soil_moisture_lowshade) <- depth.names
  colnames(soil_moisture_highshade) <- depth.names
  colnames(soil_humidity_lowshade) <- depth.names
  colnames(soil_humidity_highshade) <- depth.names
  colnames(soil_water_potential_lowshade) <- depth.names
  colnames(soil_water_potential_highshade) <- depth.names
  colnames(soil_conductivity_lowshade) <- depth.names
  colnames(soil_conductivity_highshade) <- depth.names
  colnames(soil_specific_heat_lowshade) <- c("day_of_year", "time", paste0("SP_depth_", micro$depths, "cm"))
  colnames(soil_specific_heat_highshade) <- c("day_of_year", "time", paste0("SP_depth_", micro$depths, "cm"))
  colnames(soil_density_lowshade) <- depth.names
  colnames(soil_density_highshade) <- depth.names
  plant.names <- c("day_of_year", "time", "transpiration", "leaf_water_potential", paste0("depth_", micro$depths, "cm"))
  colnames(plant_output_lowshade) <- plant.names
  colnames(plant_output_highshade) <- plant.names
  snow.names <- c("day_of_year", "time", paste0("SN", c(1, 2, 3, 4, 5, 6, 7, 8, 9)))
  colnames(snow_output_lowshade) <- snow.names
  colnames(snow_output_highshade) <- snow.names
  spectrum.names <- c("day_of_year", "time", "290", "295", "300", "305", "310", "315", "320", "330", "340", "350", "360", "370", "380", "390", "400", "420", "440", "460", "480", "500", "520", "540", "560", "580", "600", "620", "640", "660", "680", "700", "720", "740", "760", "780", "800", "820", "840", "860", "880", "900", "920", "940", "960", "980", "1000", "1020", "1080", "1100", "1120", "1140", "1160", "1180", "1200", "1220", "1240", "1260", "1280", "1300", "1320", "1380", "1400", "1420", "1440", "1460", "1480", "1500", "1540", "1580", "1600", "1620", "1640", "1660", "1700", "1720", "1780", "1800", "1860", "1900", "1950", "2000", "2020", "2050", "2100", "2120", "2150", "2200", "2260", "2300", "2320", "2350", "2380", "2400", "2420", "2450", "2490", "2500", "2600", "2700", "2800", "2900", "3000", "3100", "3200", "3300", "3400", "3500", "3600", "3700", "3800", "3900", "4000")
  colnames(direct_solar_spectrum) <- spectrum.names
  colnames(rayleigh_solar_spectrum) <- spectrum.names
  colnames(diffuse_solar_spectrum) <- spectrum.names
  return (list(micromet_lowshade=micromet_lowshade, soil_temperature_lowshade=soil_temperature_lowshade, micromet_highshade=micromet_highshade, soil_temperature_highshade=soil_temperature_highshade, soil_moisture_lowshade=soil_moisture_lowshade, soil_moisture_highshade=soil_moisture_highshade, soil_humidity_lowshade=soil_humidity_lowshade, soil_humidity_highshade=soil_humidity_highshade, soil_water_potential_lowshade=soil_water_potential_lowshade, soil_water_potential_highshade=soil_water_potential_highshade, plant_output_lowshade=plant_output_lowshade, plant_output_highshade=plant_output_highshade, snow_output_lowshade=snow_output_lowshade, snow_output_highshade=snow_output_highshade, soil_conductivity_lowshade=soil_conductivity_lowshade, soil_conductivity_highshade=soil_conductivity_highshade, soil_specific_heat_lowshade=soil_specific_heat_lowshade, soil_specific_heat_highshade=soil_specific_heat_highshade, soil_density_lowshade=soil_density_lowshade, soil_density_highshade=soil_density_highshade, direct_solar_spectrum=direct_solar_spectrum, rayleigh_solar_spectrum=rayleigh_solar_spectrum, diffuse_solar_spectrum=diffuse_solar_spectrum))
}
