#' Function to extract point data from TerraClimate
#'
#' Extracts variables needed for microclimate modelling from the monthly TerraClimate database
#' described in "Abatzoglou, J. T., Dobrowski, S. Z., Parks, S. A., & Hegewisch, K. C. (2018). TerraClimate,
#' a high-resolution global dataset of monthly climate and climatic water balance from 1958–2015.
#' Scientific Data, 5(1), 170191. https://doi.org/10.1038/sdata.2017.191"
#' This dataset includes climate change scenarios for +2 and +4 deg C warming (via parameter 'scenario')
#' Michael Kearney Feb 2022
#' @param scenario = 0, climate scenario, either 0 (historical), 2 (plus 2 deg C) or 4 (plus 4 deg C)
#' @param x = c(-5.3, 50.13), location expressed as longitude, latitude, decimal degrees
#' @param ystart = start year (min is 1958 for historical, 1985 for climate change)
#' @param yfinish = finish year (max is 'last year' for historical, 2015 for climate change)
#' @param source = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data", location where data is stored - defaults to opendap web location but can be retrieved locally if desired
#' @return TMINN, minimum daily 2m air temperature, deg C
#' @return TMAXX, maximum daily 2m air temperature, deg C
#' @return RAINFALL, daily total precipitation, mm
#' @return VPD, daily vapour pressure deficit, kPa
#' @return SRAD, daily solar radiation, W/m2
#' @return SoilMoist, volumetric soil moisture, m3/m3
#' @return WIND, daily 10m wind speed, m/s
#' @export
get_terra <- function(scenario = 0, x = c(-5.3, 50.13), ystart = 1985, yfinish = 2015, source = "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
  library(ncdf4)
  if(!(scenario %in% c(0, 2, 4))){
    message('ERROR: scenario must be either 0 (historical), 2 (plus 2) or 4 (plus 4) \n')
    break
  }
  maxyear <- as.numeric(substr(Sys.time(), start = 1, stop = 4)) - 2
  if(ystart < 1958 | ystart > maxyear){
    message(paste0('ERROR: data only available for years 1958 to ', maxyear, ' \n'))
    break
  }
  if(scenario %in% c(2, 4) & (ystart < 1985 | yfinish > 2015)){
    message('ERROR: climate change scenarios only available for years 1985 to 2015 \n')
    break
  }
  nyears <- yfinish - ystart + 1
  yearlist <- seq(ystart, yfinish)
  count <- c(1, 1, -1)
  if(scenario == 0){
    for(i in 1:length(yearlist)){
      var <- "tmax"
      ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      if(i == 1){
        lon <- ncvar_get(nc, "lon")
        lat <- ncvar_get(nc, "lat")
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
        start <- c(lonindex, latindex, 1)
      }
      message(paste0('extracting maximum air temperature data from TerraClimate for ', yearlist[i], '\n'))
      if(i == 1){
        TMAXX <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        TMAXX <- c(TMAXX, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- 'tmin'
      message(paste0('extracting minimum air temperature data from TerraClimate for ', yearlist[i], '\n'))
      ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      if(i == 1){
        TMINN <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        TMINN <- c(TMINN, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- 'ppt'
      message(paste0('extracting precipitation data from TerraClimate for ', yearlist[i], '\n'))
      ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      if(i == 1){
        RAINFALL <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        RAINFALL <- c(RAINFALL, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- 'ws'
      message(paste0('extracting wind speed data from TerraClimate for ', yearlist[i], '\n'))
      ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      if(i == 1){
        WIND <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        WIND <- c(WIND, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- 'vpd'
      message(paste0('extracting vapour pressure deficit data from TerraClimate for ', yearlist[i], '\n'))
      ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      if(i == 1){
        VPD <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        VPD <- c(VPD, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- 'srad'
      message(paste0('extracting vapour pressure deficit data from TerraClimate for ', yearlist[i], '\n'))
      ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      if(i == 1){
        SRAD <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        SRAD <- c(SRAD, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- 'soil'
      message(paste0('extracting soil moisture data from TerraClimate for ', yearlist[i], '\n'))
      ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      if(i == 1){
        SoilMoist <- as.numeric(ncvar_get(nc, varid = var, start = start, count)) / 1000 # this is originally in mm/m
      }else{
        SoilMoist <- c(SoilMoist, as.numeric(ncvar_get(nc, varid = var, start = start, count))) / 1000 # this is originally in mm/m
      }
    }
    output <- cbind(TMINN, TMAXX, RAINFALL, VPD, SRAD, SoilMoist, WIND)
  }else{
    for(i in 1:length(yearlist)){

      if(scenario == 2){
        base <- paste0(source, '_plus2C/TerraClimate_2c')
      }else{
        base <- paste0(source, '_plus4C/TerraClimate_4c')
      }

      var <- "tmax"
      ncfile <- paste0(base, "_", var, "_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      if(i == 1){
        lon <- ncvar_get(nc, "lon")
        lat <- ncvar_get(nc, "lat")
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
        start <- c(lonindex, latindex, 1)
      }

      message(paste0('extracting plus ', scenario,' maximum air temperature data from TerraClimate for ', yearlist[i], '\n'))
      if(i == 1){
        TMAXX <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        TMAXX <- c(TMAXX, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- "tmin"
      ncfile <- paste0(base, "_", var, "_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      message(paste0('extracting plus ', scenario,' minimum air temperature data from TerraClimate for ', yearlist[i], '\n'))
      if(i == 1){
        TMINN <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        TMINN <- c(TMINN, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- "ppt"
      ncfile <- paste0(base, "_", var, "_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      message(paste0('extracting plus ', scenario,' precipitation data from TerraClimate for ', yearlist[i], '\n'))
      if(i == 1){
        RAINFALL <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        RAINFALL <- c(RAINFALL, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- "vpd"
      ncfile <- paste0(base, "_", var, "_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      message(paste0('extracting plus ', scenario,' vapour pressure deficit data from TerraClimate for ', yearlist[i], '\n'))
      if(i == 1){
        VPD <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        VPD <- c(VPD, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
      var <- "soil"
      ncfile <- paste0(base, "_", var, "_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      message(paste0('extracting plus ', scenario,' soil moisture data from TerraClimate for ', yearlist[i], '\n'))
      if(i == 1){
        SoilMoist <- as.numeric(ncvar_get(nc, varid = var, start = start, count)) / 1000 # this is originally in mm/m
      }else{
        SoilMoist <- c(SoilMoist, as.numeric(ncvar_get(nc, varid = var, start = start, count)) / 1000)# this is originally in mm/m
      }
      var <- "srad"
      base <- paste0(source, '_plus2C/TerraClimate_2c')
      ncfile <- paste0(base, "_", var, "_", yearlist[i], ".nc")
      nc <- nc_open(ncfile)
      message(paste0('extracting plus ', scenario,' solar radiation data from TerraClimate for ', yearlist[i], '\n'))
      if(i == 1){
        SRAD <- as.numeric(ncvar_get(nc, varid = var, start = start, count))
      }else{
        SRAD <- c(SRAD, as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      }
    }
    output <- cbind(TMINN, TMAXX, RAINFALL, VPD, SRAD, SoilMoist)
  }
  return(output)
}
