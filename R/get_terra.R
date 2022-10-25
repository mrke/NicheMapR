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
  if(!require(futile.logger, quietly = TRUE)){
    stop('package "futile.logger" is required for this opendap extraction process. Please install it.')
  }
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
  find.nearest <- function(varname, indata){
    message(paste0('no data from TerraClimate for ', varname, ' at this site - searching for nearest pixel with data', '\n'))
    counter <- 0
    tryvec1 <- c(1, 0, -1, 0)
    tryvec2 <- c(0, 1, 0, -1)
    data.out <- matrix(data = NA, nrow = 4, ncol = 12)
    lon <- retry(ncvar_get(nc, "lon"))
    lat <- retry(ncvar_get(nc, "lat"))
    while(is.na(max(indata))){
      counter <- counter + 1
      for(ii in 1:length(tryvec1)){
        start <- c(lonindex + tryvec1[ii] * counter, latindex + tryvec2[ii] * counter, 1)
        if(start[1] <= length(lon) & start[2] <= length(lat)){
          #break
        #}else{
        data.out[ii, 1:12] <- as.numeric(ncvar_get(nc, varid = varname, start = start, count))
        }else{
          #data.out[ii, 1:12] <- rep(Inf, 12)
          message(paste0('no terraclimate data available at this site', '\n'))
          stop()
        }
      }
      good.row <- min(which(!is.na(data.out[, 1])))
      if(good.row %in% c(1, 4)){
        indata <- data.out[good.row, ]
      }
    }
    return(indata)
  }
  errors <- 0
  if(!(scenario %in% c(0, 2, 4))){
    message('ERROR: scenario must be either 0 (historical), 2 (plus 2) or 4 (plus 4) \n')
    errors <- 1
  }
  maxyear <- as.numeric(substr(Sys.time(), start = 1, stop = 4)) - 2
  if(ystart < 1958 | ystart > maxyear){
    message(paste0('ERROR: data only available for years 1958 to ', maxyear, ' \n'))
    errors <- 1
  }
  if(scenario %in% c(2, 4) & (ystart < 1985 | yfinish > 2015)){
    message('ERROR: climate change scenarios only available for years 1985 to 2015 \n')
    errors <- 1
  }
  if(errors != 1){
  nyears <- yfinish - ystart + 1
  yearlist <- seq(ystart, yfinish)
  count <- c(1, 1, -1)
  if(scenario == 0){
    for(i in 1:length(yearlist)){
      var <- "tmax"
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
       ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
       ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      if(i == 1){
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
        start <- c(lonindex, latindex, 1)
      }
      message(paste0('extracting maximum air temperature data from TerraClimate for ', yearlist[i], '\n'))
      TMAXX1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(TMAXX1))){
        TMAXX1 <- find.nearest(var, TMAXX1)
      }
      if(i == 1){
        TMAXX <- TMAXX1
      }else{
        TMAXX <- c(TMAXX, TMAXX1)
      }
      nc_close(nc)
      var <- 'tmin'
      message(paste0('extracting minimum air temperature data from TerraClimate for ', yearlist[i], '\n'))
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      TMINN1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(TMINN1))){
        TMINN1 <- find.nearest(var, TMINN1)
      }
      if(i == 1){
        TMINN <- TMINN1
      }else{
        TMINN <- c(TMINN, TMINN1)
      }
      nc_close(nc)
      var <- 'ppt'
      message(paste0('extracting precipitation data from TerraClimate for ', yearlist[i], '\n'))
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      RAINFALL1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(RAINFALL1))){
        RAINFALL1 <- find.nearest(var, RAINFALL1)
      }
      if(i == 1){
        RAINFALL <- RAINFALL1
      }else{
        RAINFALL <- c(RAINFALL, RAINFALL1)
      }
      nc_close(nc)
      var <- 'ws'
      message(paste0('extracting wind speed data from TerraClimate for ', yearlist[i], '\n'))
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      WIND1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(WIND1))){
        WIND1 <- find.nearest(var, WIND1)
      }
      if(i == 1){
        WIND <- WIND1
      }else{
        WIND <- c(WIND, WIND1)
      }
      nc_close(nc)
      var <- 'vpd'
      message(paste0('extracting vapour pressure deficit data from TerraClimate for ', yearlist[i], '\n'))
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      VPD1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(VPD1))){
        VPD1 <- find.nearest(var, VPD1)
      }
      if(i == 1){
        VPD <- VPD1
      }else{
        VPD <- c(VPD, VPD1)
      }
      nc_close(nc)
      var <- 'srad'
      message(paste0('extracting solar radiation data from TerraClimate for ', yearlist[i], '\n'))
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      SRAD1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(SRAD1))){
        SRAD1 <- find.nearest(var, SRAD1)
      }
      if(i == 1){
        SRAD <- SRAD1
      }else{
        SRAD <- c(SRAD, SRAD1)
      }
      nc_close(nc)
      var <- 'soil'
      message(paste0('extracting soil moisture data from TerraClimate for ', yearlist[i], '\n'))
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(source, "/TerraClimate_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      SoilMoist1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))) / 1000 # this is originally in mm/m
      if(is.na(max(SoilMoist1))){
        SoilMoist1 <- find.nearest(var, SoilMoist1) / 1000 # this is originally in mm/m
      }
      if(i == 1){
        SoilMoist <- SoilMoist1
      }else{
        SoilMoist <- c(SoilMoist, SoilMoist1)
      }
      nc_close(nc)
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
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
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
      TMAXX1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(TMAXX1))){
        TMAXX1 <- find.nearest(var, TMAXX1)
      }
      if(i == 1){
        TMAXX <- TMAXX1
      }else{
        TMAXX <- c(TMAXX, TMAXX1)
      }
      nc_close(nc)
      var <- "tmin"
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      message(paste0('extracting plus ', scenario,' minimum air temperature data from TerraClimate for ', yearlist[i], '\n'))
      TMINN1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(TMINN1))){
        TMINN1 <- find.nearest(var, TMINN1)
      }
      if(i == 1){
        TMINN <- TMINN1
      }else{
        TMINN <- c(TMINN, TMINN1)
      }
      nc_close(nc)
      var <- "ppt"
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      message(paste0('extracting plus ', scenario,' precipitation data from TerraClimate for ', yearlist[i], '\n'))
      RAINFALL1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(RAINFALL1))){
        RAINFALL1 <- find.nearest(var, RAINFALL1)
      }
      if(i == 1){
        RAINFALL <- RAINFALL1
      }else{
        RAINFALL <- c(RAINFALL, RAINFALL1)
      }
      var <- "vpd"
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      message(paste0('extracting plus ', scenario,' vapour pressure deficit data from TerraClimate for ', yearlist[i], '\n'))
      VPD1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(VPD1))){
        VPD1 <- find.nearest(var, VPD1)
      }
      if(i == 1){
        VPD <- VPD1
      }else{
        VPD <- c(VPD, VPD1)
      }
      nc_close(nc)
      var <- "soil"
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      message(paste0('extracting plus ', scenario,' soil moisture data from TerraClimate for ', yearlist[i], '\n'))
      SoilMoist1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count))) / 1000 # this is originally in mm/m
      if(is.na(max(SoilMoist1))){
        SoilMoist1 <- find.nearest(var, SoilMoist1) / 1000 # this is originally in mm/m
      }
      if(i == 1){
        SoilMoist <- SoilMoist1
      }else{
        SoilMoist <- c(SoilMoist, SoilMoist1)
      }
      nc_close(nc)
      var <- "srad"
      base <- paste0(source, '_plus2C/TerraClimate_2c')
      if(source == "http://thredds.northwestknowledge.net:8080/thredds/dodsC/TERRACLIMATE_ALL/data"){
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc#fillmismatch")
      }else{
        ncfile <- paste0(base, "_", var,"_", yearlist[i], ".nc")
      }
      nc <- retry(nc_open(ncfile))
      message(paste0('extracting plus ', scenario,' solar radiation data from TerraClimate for ', yearlist[i], '\n'))
      SRAD1 <- retry(as.numeric(ncvar_get(nc, varid = var, start = start, count)))
      if(is.na(max(SRAD1))){
        SRAD1 <- find.nearest(var, SRAD1)
      }
      if(i == 1){
        SRAD <- SRAD1
      }else{
        SRAD <- c(SRAD, SRAD1)
      }
      nc_close(nc)
    }
    output <- cbind(TMINN, TMAXX, RAINFALL, VPD, SRAD, SoilMoist)
  }
  return(output)
  }
}
