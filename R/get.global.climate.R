#' Downloads the global climate, soil moisture and elevation datasets
#' required to run function micro_global
#'
#' the files are:
#' global
#' derived from "New, M., Lister, D., Hulme, M. and Makin, I., 2002: A high-resolution data
#' set of surface climate over global land areas. Climate Research 21:1-25"
#' It also optionally uses a global monthly soil moisture estimate from NOAA CPC Soil Moisture http://140.172.38.100/psd/thredds/catalog/Datasets/cpcsoil/catalog.html
#' Aerosol attenuation can also be computed based on the Global Aerosol Data Set (GADS)
#' @param folder Path to the folder you want to install the global climate data in
#' @usage get.global.climate(folder)
#' @export
get.global.climate <- function(folder="c:/globalclimate"){
ANSWER<-readline(prompt = "This function downloads and extracts 0.93 GB of data, type 'y' if you want to continue: ")
if(substr(ANSWER, 1, 1) == "y"){
  if(dir.exists(folder)==FALSE){
    dir.create(folder)
  }
  climate.file<-"https://www.dropbox.com/s/9eg1aaooid4lbqv/global%20climate.zip?raw=1"
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
