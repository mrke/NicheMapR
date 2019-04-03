#' rungads
#'
#' Function to extract data from the Fortran Global Aerosol Database software,
#' for a specific location, season and relative humidity
#' see http://www.mpimet.mpg.de/fileadmin/publikationen/Reports/MPI-Report_243.pdf
#' @param lat Latitude in decimal degrees
#' @param lon Longitude in decimal degrees
#' @param relhum Integer specifying relative humdity to use, 1:8 corresponds to 0,50,70,80,90,95,98,99 percent relative humidity
#' @param season Season to obtain data for, 0 = summer, 1 = winter
#' @return optdep A vector of wavelength-specific optical depths
#' @useDynLib "NicheMapR"
#' @export
rungads <- function(lat, lon, relhum, season) {
  os = Sys.info()['sysname']
  if (os == "Windows") {
    if (R.Version()$arch=="x86_64") {
      libpath='/NicheMapR/libs/x64/NicheMapR.dll'
    } else {
      libpath='/NicheMapR/libs/i386/NicheMapR.dll'
    }
  } else if (os == "Linux") {
    libpath='/NicheMapR/libs/NicheMapR.so'
  } else if (os == "Darwin") {
    libpath='/NicheMapR/libs/NicheMapR.so'
  }
  if(is.loaded("gads", "NicheMapR", type = "FORTRAN")==FALSE){
    dyn.load(paste(lib.loc = .libPaths()[1],libpath,sep=""))
  }
  curdir<-getwd()
  setwd(path.package("NicheMapR"))
  lat5s<-seq(-90,90,5) #lat range for GADS
  lon5s<-seq(-180,175,5) #long range for GADS
  lat5<-lat5s[which.min(abs(lat5s-lat))] # get nearest latitude square for input location
  lon5<-lon5s[which.min(abs(lon5s-lon))] # get nearest longitude square for input location
    a <- .Fortran("gads",
    as.double(lat5),
    as.double(lon5),
    as.double(relhum),
    as.double(season),
    optdep=matrix(data = 0., nrow = 25, ncol = 2), PACKAGE="NicheMapR")
  setwd(curdir)
  optdep <- matrix(data = 0., nrow = 25, ncol = 2)
  storage.mode(optdep)<-"double"
  optdep<-a$optdep
  optdep.names<-c("LAMBDA","OPTDEPTH")
  colnames(optdep)<-optdep.names
  optdep[,1]<-optdep[,1]*1000
  return (optdep)
}
