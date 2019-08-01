#' ectorun
#'
#' R wrapper for Fortran binary of Niche Mapper microclimate model
#' @param ecto A vector of input variables for the microclimate model
#' @param metout The above ground micrometeorological conditions under the minimum specified shade
#' @param shadmet The above ground micrometeorological conditions under the maximum specified shade
#' @param soil Hourly predictions of the soil temperatures under the minimum specified shade
#' @param shadsoil Hourly predictions of the soil temperatures under the maximum specified shade
#' @param soilmoist Hourly predictions of the soil moisture under the minimum specified shade
#' @param shadmoist Hourly predictions of the soil moisture under the maximum specified shade
#' @param soilpot Hourly predictions of the soil water potential under the minimum specified shade
#' @param shadpot Hourly predictions of the soil water potential under the maximum specified shade
#' @param humid Hourly predictions of the soil humidity under the minimum specified shade
#' @param shadhumid Hourly predictions of the soil humidity under the maximum specified shade
#' @param DEP Depths used for the microclimate model
#' @param rainfall Daily rainfall
#' @param rainhr Hourly rainfall (overwrites rainfall if non-negative
#' @param debmod Dynamic Energy Budget (DEB) model paramters
#' @param deblast Initial DEB state
#' @param foodwaters Food water content (add units)
#' @param foodlevels Food levels (add units)
#' @param wetlandTemps Tempature of water body
#' @param wetlandDepths Depth of water body
#' @param GLMtemps General Lake Model water temperatures
#' @param GLMO2s General Lake Model water O2 partial pressures
#' @param GLMsalts General Lake Model water salinities
#' @param GLMpHs General Lake Model water pHs
#' @param GLMfoods General Lake Model water food densities
#' @param arrhenius Arrhenius thermal response
#' @param thermal_stages Stage-specific thermal physiology
#' @param behav_stages Stage-specific behaviour
#' @param water_stages Stage-specific water related parameters
#' @param nutri_stages Stage-specific nutrition related parameters
#' @param MINSHADES Minimum shade levels
#' @param MAXSHADES Maximum shade levels
#' @param S_instar For the Insect DEB model
#' @useDynLib "NicheMapR"
#' @export
ectorun <- function(ecto) {
  ndays<-ecto$ndays
  a <- .Fortran("ectotherm",
                as.integer(ecto$ndays),
                as.integer(ecto$nstages),
                as.double(ecto$ectoinput),
                as.double(ecto$metout),
                as.double(ecto$shadmet),
                as.double(ecto$soil),
                as.double(ecto$shadsoil),
                as.double(ecto$soilmoist),
                as.double(ecto$shadmoist),
                as.double(ecto$soilpot),
                as.double(ecto$shadpot),
                as.double(ecto$humid),
                as.double(ecto$shadhumid),
                as.double(ecto$DEP),
                as.double(ecto$rainfall),
                as.double(ecto$rainhr),
                as.double(ecto$debmod),
                as.double(ecto$deblast),
                as.double(ecto$foodwaters),
                as.double(ecto$foodlevels),
                as.double(ecto$wetlandTemps),
                as.double(ecto$wetlandDepths),
                as.double(ecto$GLMtemps),
                as.double(ecto$GLMO2s),
                as.double(ecto$GLMsalts),
                as.double(ecto$GLMpHs),
                as.double(ecto$GLMfoods),
                as.double(ecto$arrhenius),
                as.double(ecto$thermal_stages),
                as.double(ecto$behav_stages),
                as.double(ecto$water_stages),
                as.double(ecto$nutri_stages),
                as.double(ecto$minshades),
                as.double(ecto$maxshades),
                as.double(ecto$S_instar),
                environ=matrix(data = 0, nrow = ndays * 24, ncol = 28),
                enbal=matrix(data = 0, nrow = ndays * 24, ncol = 13),
                masbal=matrix(data = 0, nrow = ndays * 24, ncol = 19),
                debout=matrix(data = 0, nrow = ndays * 24, ncol = 21),
                yearout=matrix(data = 0, nrow = 1, ncol = 20),
                yearsout=matrix(data = 0, nrow = ceiling(ndays / 365), ncol = 43),
                PACKAGE = "NicheMapR"
  )

  environ <- matrix(data = 0, nrow = 24 * ndays, ncol = 28)
  enbal <- matrix(data = 0, nrow = 24 * ndays, ncol = 13)
  masbal <- matrix(data = 0, nrow = 24 * ndays, ncol = 19)
  debout <- matrix(data = 0, nrow = 24 * ndays, ncol = 21)
  yearout <- matrix(data = 0, nrow = 1, ncol = 20)
  yearsout <- matrix(data = 0, nrow = ceiling(ndays / 365), ncol = 43)

  storage.mode(environ)<-"double"
  storage.mode(enbal)<-"double"
  storage.mode(masbal)<-"double"
  storage.mode(debout)<-"double"
  storage.mode(yearout)<-"double"
  storage.mode(yearsout)<-"double"
  environ<-a$environ
  enbal<-a$enbal
  masbal<-a$masbal
  debout<-a$debout
  environ[,4] <- environ[,4] - 1 # make first hour midnight
  enbal[,4] <- enbal[,4] - 1 # make first hour midnight
  masbal[,4] <- masbal[,4] - 1 # make first hour midnight
  debout[,4] <- debout[,4] - 1 # make first hour midnight
  yearout<-a$yearout
  yearsout<-a$yearsout
  environ.names<-c("DOY","YEAR","DAY","TIME","TC","SHADE","SOLAR","DEP","ACT","TA","TSUB","TSKY","VEL","RELHUM","ZEN","CONDEP","WATERTEMP","DAYLENGTH","WINGANGLE","WINGTEMP","FLYING","FLYTIME","PO2WATER","SALWATER","ABSAN","PTCOND","POSTURE","PANT")
  enbal.names<-c("DOY","YEAR","DAY","TIME","QSOL","QIRIN","QMET","QEVAP","QIROUT","QCONV","QCOND","ENB","NTRY")
  masbal.names<-c("DOY","YEAR","DAY","TIME","O2_ml","CO2_ml","NWASTE_g","H2OFree_g","H2OMet_g","DryFood_g","WetFood_g","DryFaeces_g","WetFaeces_G","Urine_g","H2OResp_g","H2OCut_g","H2OEye_g","H2OBal_g","H2OCumBal_g")
  debout.names<-c("DOY","YEAR","DAY","TIME","STAGE","V","E","E_H","L_W","WETMASS","WETGONAD","WETGUT","PCT_DESIC","E_R","E_B","BREEDING","PREGNANT","V_BABY","E_BABY","H_S","P_SURV")
  yearout.names<-c("DEVTIME","BIRTHDAY","BIRTHMASS","MONMATURE","MONREPRO","LENREPRO","FECUNDITY","CLUTCHES","ANNUALACT","MINRESERVE","LASTFOOD","TOTFOOD","MINTB","MAXTB","Pct_Des","LifeSpan","GenTime","R0","rmax","Length")
  yearsout.names<-c("YEAR","MaxStg","MaxWgt","MaxLen","Tmax","Tmin","MinRes","MaxDes","MinShade","MaxShade","MinDep","MaxDep","Bsk","Forage","Dist","Food","Drink","NWaste","Faeces","O2","Clutch","Fec","CauseDeath","tLay","tEgg","tStg1","tStg2","tStg3","tStg4","tStg5","tStg6","tStg7","tStg8","mStg1","mStg2","mStg3","mStg4","mStg5","mStg6","mStg7","mStg8","surviv","deathstage")

  colnames(environ)<-environ.names
  colnames(enbal)<-enbal.names
  colnames(masbal)<-masbal.names
  colnames(debout)<-debout.names
  colnames(yearout)<-yearout.names
  colnames(yearsout)<-yearsout.names
  return (list(environ=environ, enbal=enbal, masbal=masbal, debout=debout, yearout=yearout, yearsout=yearsout))
}
