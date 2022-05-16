#' SOLVENDO
#'
#' R wrapper for Fortran binary of SOLVENDO (endotherm model)
#' @encoding UTF-8
#' @param SOLVENDO.input input vector
#' @export
SOLVENDO <- function(SOLVENDO.input){
  a <- .Fortran("SOLVENDO",
                as.double(SOLVENDO.input),
                treg=matrix(data = 0., nrow = 1, ncol = 15),
                morph=matrix(data = 0., nrow = 1, ncol = 20),
                enbal=matrix(data = 0., nrow = 1, ncol = 10),
                masbal=matrix(data = 0., nrow = 1, ncol = 10),
                PACKAGE = "NicheMapR")

  treg <- matrix(data = 0., nrow = 1, ncol = 15)
  morph <- matrix(data = 0., nrow = 1, ncol = 20)
  masbal <- matrix(data = 0., nrow = 1, ncol = 10)
  enbal <- matrix(data = 0., nrow = 1, ncol = 10)

  storage.mode(treg)<-"double"
  storage.mode(morph)<-"double"
  storage.mode(enbal)<-"double"
  storage.mode(masbal)<-"double"
  treg <- a$treg
  morph <- a$morph
  enbal <- a$enbal
  masbal <- a$masbal

  if(max(treg) == 0){
    warning("A solution could not be found and panting/'sweating' options are exhausted; try allowing greater evaporation or allowing higher body maximum body temperature")
  }

  treg.names<-c("TC", "TLUNG", "TSKIN_D", "TSKIN_V", "TFA_D", "TFA_V", "SHAPE_B", "PANT", "PCTWET", "K_FLESH", "K_FUR_EFF", "K_FUR_D", "K_FUR_V", "K_COMPFUR", "Q10")
  morph.names<-c("AREA", "VOLUME", "CHAR_DIM", "MASS_FAT", "FAT_THICK", "FLESH_VOL", "LENGTH", "WIDTH", "HEIGHT", "DIAM_FLESH", "DIAM_FUR", "AREA_SIL", "AREA_SILN", "AREA_ASILP", "AREA_SKIN", "AREA_SKIN_EVAP", "AREA_CONV", "AREA_COND", "F_SKY", "F_GROUND")
  enbal.names<-c("QSOL", "QIRIN", "QGEN", "QEVAP", "QIROUT", "QCONV", "QCOND", "ENB", "NTRY", "SUCCESS")
  masbal.names<-c("AIR_L", "O2_L", "H2OResp_g", "H2OCut_g", "O2_mol_in", "O2_mol_out", "N2_mol_in", "N2_mol_out", "AIR_mol_in", "AIR_mol_out")

  colnames(treg)<-treg.names
  colnames(morph)<-morph.names
  colnames(enbal)<-enbal.names
  colnames(masbal)<-masbal.names
  return (list(treg=treg, morph=morph, enbal=enbal, masbal=masbal))
 }
