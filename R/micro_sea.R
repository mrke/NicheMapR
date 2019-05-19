# get SST time series
# based on https://cran.r-project.org/web/packages/heatwaveR/vignettes/OISST_preparation.html

#install.packages('rerddap') # for workign with ERDDAP servers
library(dplyr)
library(rerddap)
library(ncdf4)

# The information for the NOAA OISST data, using rerddap
info(datasetid = "ncdc_oisst_v2_avhrr_by_time_zlev_lat_lon", url = "https://www.ncei.noaa.gov/erddap/")
info('erdVHNchlamday')
info('erdSW1chla8day')

# This function expects the user to provide it with two values
# that match the time format of the target OISST dataset
oisst_res <- griddap(x = "ncdc_oisst_v2_avhrr_by_time_zlev_lat_lon",
                       url = "https://www.ncei.noaa.gov/erddap/",
                       time = times,
                       depth = c(0, 0),
                       latitude = c(-40, -35),
                       longitude = c(15, 21),
                       fields = "sst")

dstart <- "1997-09-01T00:00:00Z"
dfinish <- "2010-12-15T00:00:00Z"
lon <- 15
lat <- -40
dates <- seq(as.POSIXct(dstart), as.POSIXct(dfinish) + 3600*24, by = "days")
chla <- griddap(x = "erdSW1chla8day",
        time = c(dstart, dfinish),
        latitude = c(lat, lat),
        longitude = c(lon, lon),
        fields = "chlorophyll")

sst <- griddap(x = "ncdc_oisst_v2_avhrr_by_time_zlev_lat_lon",
                  url = "https://www.ncei.noaa.gov/erddap/",
                  time = c(dstart, dfinish),
                  depth = c(0, 0),
                  latitude = c(lat, lat),
                  longitude = c(lon, lon),
                  fields = "sst")

nc <- nc_open(sst$summary$filename)
sst_data <- ncvar_get(nc, varid = 'sst')
plot(dates, sst_data, type = 'l')

nc <- nc_open(chla$summary$filename)
chla_data <- ncvar_get(nc, varid = 'chlorophyll')
plot(dates[seq(1, length(dates), 8)][1:length(chla_data)], chla_data, type = 'l')

OISST_prep <- function(nc_file){

  # Open the NetCDF connection
  nc <- nc_open(nc_file$summary$filename)

  # Extract the SST values and add the lon/lat/time dimension names
  res <- ncvar_get(nc, varid = "sst")
  dimnames(res) <- list(lon = nc$dim$longitude$vals,
                        lat = nc$dim$latitude$vals,
                        t = nc$dim$time$vals)

  # Convert the data into a 'long' dataframe for use in the 'tidyverse' ecosystem
  res <- as.data.frame(reshape2::melt(res, value.name = "temp"), row.names = NULL) %>%
    mutate(t = as.Date(as.POSIXct(t, origin = "1970-01-01 00:00:00")),
           temp = round(temp, 2))

  # Close the NetCDF connection and finish
  nc_close(nc)
  return(res)
}

OISST1 <- OISST_sub(c("1981-09-01T00:00:00Z", "1990-12-31T00:00:00Z"), latitude = c(-40, -40), longitude = c(15, 15))
OISST2 <- OISST_sub(c("1991-01-01T00:00:00Z", "1999-12-31T00:00:00Z"))
OISST3 <- OISST_sub(c("2000-01-01T00:00:00Z", "2008-12-31T00:00:00Z"))
OISST4 <- OISST_sub(c("2009-01-01T00:00:00Z", "2013-12-03T00:00:00Z"))
OISST5 <- OISST_sub(c("2014-01-01T00:00:00Z", "2018-12-03T00:00:00Z"))
