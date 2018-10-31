#' Wrapper for google_geocode from the googleway package.
#'
#' A Google Maps API Key is required for this function, which can be obtained at https://cloud.google.com/maps-platform/ (credit card details need to be provided but it can be used for free for a certain limit of queries)
#' @encoding UTF-8
#' @param place String describing the location, e.g. "Alice Springs, Australia"
#' @param key Your API key
#' @return longlat The longitude and latitude in decimal degrees
#' @usage get.longlat(place = "Alice Springs, Australia", key = 'your key')
get.longlat <- function(place, key){
  if(!require(googleway, quietly = TRUE)){
   stop('package "googleway" is required to do a geocode search and you also need a Google Maps API key. Please install it.')
  }
  library(googleway)
  longlat <- google_geocode(address = place, key = key)$results[[3]][2]
  longlat <- c(as.numeric(longlat$location[2]), as.numeric(longlat$location[1]))
  return(longlat)
}
