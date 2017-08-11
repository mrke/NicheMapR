.onLoad <- function(libname, pkgname) {
  os = Sys.info()['sysname']
  if (os == "Windows") {
      if (R.Version()$arch=="x86_64") {
        microclimate_path='/NicheMapR/libs/win/x64/microclimate.dll'
        ecotherm_path='/NicheMapR/libs/win/x64/ectotherm.dll'
      } else {
        microclimate_path='/NicheMapR/libs/win/i386/microclimate.dll'
        ecotherm_path='/NicheMapR/libs/win/i386/ectotherm.dll'
      }
  } else if (os == "Linux") {
      microclimate_path='/NicheMapR/libs/linux/MICROCLIMATE.so'
      ecotherm_path='/NicheMapR/libs/linux/ECTOTHERM.so'
  } else if (os == "Darwin") {
      microclimate_path='/NicheMapR/libs/mac/MICROCLIMATE.so'
      ectotherm_path='/NicheMapR/libs/mac/ECTOTHERM.so'
  }
  if(is.loaded("microclimate", "MICROCLIMATE", type = "FORTRAN")==FALSE){
    dyn.load(paste(lib.loc = .libPaths()[1],microclimate_path,sep=""))
  }
  if(is.loaded("ectotherm", "ECTOTHERM", type = "FORTRAN")==FALSE){
    dyn.load(paste(lib.loc = .libPaths()[1],ecotherm_path,sep=""))
  }
}
