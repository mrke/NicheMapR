#' Function to assist with interpolated data in leap years
#'
#' @encoding UTF-8
#' @param indata data needing to be interpolated
#' @param yearlist list of years to do
#' @param mult how many elements per day? (1 = daily, 24 = hourly)
#' @export
leapfix <- function(indata, yearlist, mult = 1){
  leapyears <- seq(1900, 2060, 4)
  for(j in 1:length(yearlist)){
    if(yearlist[j] %in% leapyears){# add day for leap year if needed
      if(mult == 1){
        data <- c(indata[1:59], indata[59], indata[60:365])
      }else{
        data <- c(indata[1:(59 * mult)], indata[(58*mult+1):(59 * mult)], indata[(59 * mult + 1):(365 * mult)])
      }
    }else{
      data <- indata
    }
    if(j == 1){
      alldata <- data
    }else{
      alldata <- c(alldata, data)
    }
  }
  return(alldata)
}
