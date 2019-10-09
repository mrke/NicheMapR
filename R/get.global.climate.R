#' get.global.climate
#'
#' Downloads the global climate, elevation and soil moisture datasets used in the function micro_global. The global
#' climate dataset global.nc is pre-generated (using code in the 'build.global.climate.R' function
#' generated from "New, M., Lister, D., Hulme, M. and Makin, I., 2002: A high-resolution
#' data set of surface climate over global land areas. Climate Research 21:1-25", specically 10' estimates of monthly
#' mean (over 1961-1990) preceipitation, wet-days, air temperature, diurnal temperature range, relative humidity and
#' wind speed. Cloud cover comes from a bilinear interpolation of a lower resolution version of this dataset (New, M.,
#' M. Hulme, and P. D. Jones. 1999. Representing twentieth century space-time climate variability. Part 1: development
#' of a 1961-90 mean monthly terrestrial climatology. Journal of Climate 12:829-856.). The micro_global function
#' optionaly uses a  global monthly soil moisture estimate from NOAA CPC Soil Moisture http://140.172.38.100/psd/
#' thredds/catalog/Datasets/cpcsoil/catalog.html
#' @param folder Path to the folder you want to install the global climate data in
#' @usage get.global.climate(folder)
#' @export
get.global.climate <- function(folder="c:/globalclimate"){
ANSWER<-readline(prompt = "This function downloads and extracts 0.5 GB of data, type 'y' if you want to continue: ")
if(substr(ANSWER, 1, 1) == "y"){
  # check if user put a slash at end, remove if so
  if(substr(folder, nchar(folder), nchar(folder)) == "/"){folder <- substr(folder, 1, nchar(folder)-1)}
  if(substr(folder, nchar(folder)-1, nchar(folder)) == "\\"){folder <- substr(folder, 1, nchar(folder)-2)}
  if(dir.exists(folder)==FALSE){
    dir.create(folder)
  }
  climate.file<-"https://github.com/mrke/NicheMapR/releases/download/v2.0.0/global.climate.zip"
  destin<-paste(folder,"/global climate.zip",sep="")
  cat("downloading 'Global Climate.zip' \n")
  download.file(climate.file, destin, mode="wb")
  cat("unzipping 'Global Climate.zip' \n")
  unzip(destin, exdir = folder)
  file.remove(destin)
  save(folder,file = paste(.libPaths()[1],"/gcfolder.rda",sep=""))
  cat(paste("global climate data is downloaded and extracted into folder ",folder,sep=""))
  cat(paste("folder location has been saved into ",paste(.libPaths()[1],"/gcfolder.rda",sep="")))
}else{
  print("aborting installation of global climate data")
}
}
