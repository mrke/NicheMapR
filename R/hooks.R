.onLoad <- function(libname, pkgname) {
  handleLibs("load")
}

.onUnload <- function(libpath) {
  handleLibs("unload")
}

handleLibs <- function(action) {
  # Load/Unload dynamic libraries manually. 
  # This is a workaround as useDynLib() does not easily handle mac and linux 
  # libraries with the same name and extension
  os = Sys.info()["sysname"]
  if (os == "Windows") {
      if (R.Version()$arch == "x86_64") {
          micro_path = "/NicheMapR/libs/win/x64/microclimate.dll"
          ecto_path = "/NicheMapR/libs/win/x64/ectotherm.dll"
      } else {
          micro_path = "/NicheMapR/libs/win/i386/microclimate.dll"
          ecto_path = "/NicheMapR/libs/win/i386/ectotherm.dll"
      }
  } else if (os == "Linux") {
      micro_path = "/NicheMapR/libs/linux/MICROCLIMATE.so"
      ecto_path = "/NicheMapR/libs/linux/ECTOTHERM.so"
  } else if (os == "Darwin") {
      micro_path = "/NicheMapR/libs/mac/MICROCLIMATE.so"
      ecto_path = "/NicheMapR/libs/mac/ECTOTHERM.so"
  }

  micro_lib <- paste(lib.loc = .libPaths()[1], micro_path, sep = "")
  ecto_lib <- paste(lib.loc = .libPaths()[1], ecto_path, sep = "")
  micro_loaded <- is.loaded("microclimate", "MICROCLIMATE", type = "FORTRAN")
  ecto_loaded <- is.loaded("ectotherm", "ECTOTHERM", type = "FORTRAN")

  if (action == "load") {
    if (!micro_loaded) {
        dyn.load(micro_lib)
    }
    if (!ecto_loaded) {
        dyn.load(ecto_lib)
    }
  } else if (action == "unload") {
    if (micro_loaded) {
        dyn.unload(micro_lib)
    }
    if (ecto_loaded) {
        dyn.unload(ecto_lib)
    }
  }
}
