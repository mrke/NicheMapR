#' microclimate model
#'
#' R wrapper for Fortran binary of Niche Mapper microclimate model
#' @param micro A vector of input variables for the microclimate model
#' @return metout The above ground micrometeorological conditions under the minimum specified shade
#' @return shadmet The above ground micrometeorological conditions under the maximum specified shade
#' @return soil Hourly predictions of the soil temperatures under the minimum specified shade
#' @return shadsoil Hourly predictions of the soil temperatures under the maximum specified shade
#' @return soilmoist Hourly predictions of the soil moisture under the minimum specified shade
#' @return shadmoist Hourly predictions of the soil moisture under the maximum specified shade
#' @return soilpot Hourly predictions of the soil water potential under the minimum specified shade
#' @return shadpot Hourly predictions of the soil water potential under the maximum specified shade
#' @return humid Hourly predictions of the soil humidity under the minimum specified shade
#' @return shadhumid Hourly predictions of the soil humidity under the maximum specified shade
#' @return plant Hourly predictions of plant variables under the minimum specified shade
#' @return shadplant Hourly predictions of plant variables under the maximum specified shade
#' @return sunsnow Hourly predictions of the snow temperature under the minimum specified shade
#' @return shdsnow Hourly predictions of the snow temperature under the maximum specified shade
#' @return tcond Hourly predictions of the soil thermal conductivity under the minimum specified shade
#' @return shadtcond Hourly predictions of the soil thermal conductivity under the maximum specified shade
#' @return specheat Hourly predictions of the soil specific heat capacity under the minimum specified shade
#' @return shadspecheat Hourly predictions of soil specific heat capacity under the maximum specified shade
#' @return densit Hourly predictions of the soil density under the minimum specified shade
#' @return shaddensit Hourly predictions of the soil density under the maximum specified shade
#' @useDynLib "NicheMapR"
#' @export
microclimate <- function(micro) {
  doynum<-micro$microinput[1]
  errors <- 0
  if(length(micro$microinput) != 74){
    message("ERROR: microinput has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$MINSHADES) != doynum){
    message("ERROR: MINSHADES has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$MAXSHADES) != doynum){
    message("ERROR: MAXSHADES has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$SLES) != doynum){
    message("ERROR: SLES has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$DEP) != 10){
    message("ERROR: DEP has the wrong number of inputs \n")
    errors <- 1
  }
  if(nrow(micro$moists) != 10 | ncol(micro$moists) != doynum){
    message("ERROR: moists has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$RHMINN) != doynum){
    message("ERROR: RHMINN has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$RHMAXX) != doynum){
    message("ERROR: RHMAXX has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$CCMINN) != doynum){
    message("ERROR: CCMINN has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$CCMAXX) != doynum){
    message("ERROR: CCMAXX has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$WNMINN) != doynum){
    message("ERROR: WNMINN has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$WNMAXX) != doynum){
    message("ERROR: WNMAXX has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$TMINN) != doynum){
    message("ERROR: TMINN has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$TMAXX) != doynum){
    message("ERROR: TMAXX has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$REFLS) != doynum){
    message("ERROR: REFLS has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$PCTWET) != doynum){
    message("ERROR: PCTWET has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$RAINFALL) != doynum){
    message("ERROR: RAINFALL has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$tannulrun) != doynum){
    message("ERROR: tannulrun has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$LAI) != doynum){
    message("ERROR: LAI has the wrong number of inputs \n")
    errors <- 1
  }
  if(nrow(micro$tides) != doynum * 24){
    message("ERROR: tides has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$TAIRhr) != doynum * 24){
    message("ERROR: TAIRhr has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$RHhr) != doynum * 24){
    message("ERROR: RHhr has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$WNhr) != doynum * 24){
    message("ERROR: WNhr has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$CLDhr) != doynum * 24){
    message("ERROR: CLDhr has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$SOLRhr) != doynum * 24){
    message("ERROR: SOLRhr has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$RAINhr) != doynum * 24){
    message("ERROR: RAINhr has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$ZENhr) != doynum * 24){
    message("ERROR: ZENhr has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$IRDhr) != doynum * 24){
    message("ERROR: IRDhr has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$PE) != 19){
    message("ERROR: PE has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$BB) != 19){
    message("ERROR: BB has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$KS) != 19){
    message("ERROR: KS has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$BD) != 19){
    message("ERROR: BD has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$DD) != 19){
    message("ERROR: DD has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$L) != 19){
    message("ERROR: L has the wrong number of inputs \n")
    errors <- 1
  }
  if(length(micro$hori) != 24){
    message("ERROR: hori has the wrong number of inputs \n")
    errors <- 1
  }
  if(errors == 0){
  a <- .Fortran("microclimate",
                as.integer(doynum),
                as.double(micro$microinput),
                as.double(micro$doy),
                as.double(micro$SLES),
                as.double(micro$DEP),
                as.double(micro$MAXSHADES),
                as.double(micro$MINSHADES),
                as.double(micro$Nodes),
                as.double(micro$RHMAXX),
                as.double(micro$RHMINN),
                as.double(micro$CCMAXX),
                as.double(micro$CCMINN),
                as.double(micro$WNMAXX),
                as.double(micro$WNMINN),
                as.double(micro$TMAXX),
                as.double(micro$TMINN),
                as.double(micro$REFLS),
                as.double(micro$PCTWET),
                as.double(micro$soilinit),
                as.double(micro$hori),
                as.double(micro$TAI),
                as.double(micro$soilprops),
                as.double(micro$moists),
                as.double(micro$RAINFALL),
                as.double(micro$tannulrun),
                as.double(micro$tides),
                as.double(micro$PE),
                as.double(micro$KS),
                as.double(micro$BB),
                as.double(micro$BD),
                as.double(micro$DD),
                as.double(micro$L),
                as.double(micro$LAI),
                as.double(micro$TAIRhr),
                as.double(micro$RHhr),
                as.double(micro$WNhr),
                as.double(micro$CLDhr),
                as.double(micro$SOLRhr),
                as.double(micro$RAINhr),
                as.double(micro$ZENhr),
                as.double(micro$IRDhr),
                metout=matrix(data = 0, nrow = 24 * doynum, ncol = 19),
                soil=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                shadmet=matrix(data = 0, nrow = 24 * doynum, ncol = 19),
                shadsoil=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soilmoist=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                shadmoist=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                humid=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                shadhumid=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                soilpot=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                shadpot=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                sunsnow=matrix(data = 0, nrow = 24 * doynum, ncol = 11),
                shdsnow=matrix(data = 0, nrow = 24 * doynum, ncol = 11),
                plant=matrix(data = 0, nrow = 24 * doynum, ncol = 14),
                shadplant=matrix(data = 0, nrow = 24 * doynum, ncol = 14),
                tcond=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                shadtcond=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                specheat=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                shadspecheat=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                densit=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                shaddensit=matrix(data = 0, nrow = 24 * doynum, ncol = 12),
                drlam=matrix(data = 0, nrow = 24 * doynum, ncol = 113),
                drrlam=matrix(data = 0, nrow = 24 * doynum, ncol = 113),
                srlam=matrix(data = 0, nrow = 24 * doynum, ncol = 113), PACKAGE = "NicheMapR")
  }
  metout <- matrix(data = 0, nrow = 24 * doynum, ncol = 19)
  shadmet <- matrix(data = 0, nrow = 24 * doynum, ncol = 19)
  soil <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  shadsoil <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soilmoist <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  shadmoist <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  humid <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  shadhumid <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  soilpot <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  shadpot <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  sunsnow <- matrix(data = 0, nrow = 24 * doynum, ncol = 11)
  shdsnow <- matrix(data = 0, nrow = 24 * doynum, ncol = 11)
  plant <- matrix(data = 0, nrow = 24 * doynum, ncol = 14)
  shadplant <- matrix(data = 0, nrow = 24 * doynum, ncol = 14)
  tcond <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  shadtcond <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  specheat <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  shadspecheat <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  densit <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  shaddensit <- matrix(data = 0, nrow = 24 * doynum, ncol = 12)
  drlam <- matrix(data = 0, nrow = 24 * doynum, ncol = 113)
  drrlam <- matrix(data = 0, nrow = 24 * doynum, ncol = 113)
  srlam <- matrix(data = 0, nrow = 24 * doynum, ncol = 113)
  storage.mode(metout) <- "double"
  storage.mode(shadmet) <- "double"
  storage.mode(soil) <- "double"
  storage.mode(shadsoil) <- "double"
  storage.mode(soilmoist) <- "double"
  storage.mode(shadmoist) <- "double"
  storage.mode(humid) <- "double"
  storage.mode(shadhumid) <- "double"
  storage.mode(soilpot) <- "double"
  storage.mode(shadpot) <- "double"
  storage.mode(sunsnow) <- "double"
  storage.mode(shdsnow) <- "double"
  storage.mode(plant) <- "double"
  storage.mode(shadplant) <- "double"
  storage.mode(tcond) <- "double"
  storage.mode(shadtcond) <- "double"
  storage.mode(specheat) <- "double"
  storage.mode(shadspecheat) <- "double"
  storage.mode(densit) <- "double"
  storage.mode(shaddensit) <- "double"
  storage.mode(drlam) <- "double"
  storage.mode(drrlam) <- "double"
  storage.mode(srlam) <- "double"
  metout <- a$metout
  shadmet <- a$shadmet
  soil <- a$soil
  shadsoil <- a$shadsoil
  soilmoist <- a$soilmoist
  shadmoist <- a$shadmoist
  humid <- a$humid
  shadhumid <- a$shadhumid
  soilpot <- a$soilpot
  shadpot <- a$shadpot
  sunsnow <- a$sunsnow
  shdsnow <- a$shdsnow
  plant <- a$plant
  shadplant <- a$shadplant
  tcond <- a$tcond
  shadtcond <- a$shadtcond
  specheat <- a$specheat
  shadspecheat <- a$shadspecheat
  densit <- a$densit
  shaddensit <- a$shaddensit
  drlam <- a$drlam
  drrlam <- a$drrlam
  srlam <- a$srlam
  metout.names <- c("DOY", "TIME", "TALOC", "TAREF", "RHLOC", "RH", "VLOC", "VREF", "SNOWMELT", "POOLDEP", "PCTWET", "ZEN", "SOLR", "TSKYC", "DEW", "FROST", "SNOWFALL", "SNOWDEP", "SNOWDENS")
  colnames(metout) <- metout.names
  colnames(shadmet) <- metout.names
  soil.names <- c("DOY", "TIME", paste0("D", micro$DEP, "cm"))
  colnames(soil) <- soil.names
  colnames(shadsoil) <- soil.names
  moist.names <- c("DOY", "TIME", paste0("WC", micro$DEP, "cm"))
  humid.names <- c("DOY", "TIME", paste0("RH", micro$DEP, "cm"))
  pot.names <- c("DOY", "TIME", paste0("PT", micro$DEP, "cm"))
  plant.names <- c("DOY", "TIME", "TRANS", "LPT", paste0("RPT", micro$DEP, "cm"))
  tcond.names <- c("DOY", "TIME", paste0("TC", micro$DEP, "cm"))
  spheat.names <- c("DOY", "TIME", paste0("SP", micro$DEP, "cm"))
  denst.names <- c("DOY", "TIME", paste0("DE", micro$DEP, "cm"))
  colnames(soilmoist) <- moist.names
  colnames(shadmoist) <- moist.names
  colnames(humid) <- humid.names
  colnames(shadhumid) <- humid.names
  colnames(soilpot) <- pot.names
  colnames(shadpot) <- pot.names
  colnames(plant) <- plant.names
  colnames(shadplant) <- plant.names
  colnames(tcond) <- tcond.names
  colnames(shadtcond) <- tcond.names
  colnames(specheat) <- spheat.names
  colnames(shadspecheat) <- spheat.names
  colnames(densit) <- denst.names
  colnames(shaddensit) <- denst.names
  snow.names <- c("DOY", "TIME", paste0("SN", c(1, 2, 3, 4, 5, 6, 7, 8, 9)))
  colnames(sunsnow) <- snow.names
  colnames(shdsnow) <- snow.names
  drlam.colnames <- c("DOY", "TIME", "290", "295", "300", "305", "310", "315", "320", "330", "340", "350", "360", "370", "380", "390", "400", "420", "440", "460", "480", "500", "520", "540", "560", "580", "600", "620", "640", "660", "680", "700", "720", "740", "760", "780", "800", "820", "840", "860", "880", "900", "920", "940", "960", "980", "1000", "1020", "1080", "1100", "1120", "1140", "1160", "1180", "1200", "1220", "1240", "1260", "1280", "1300", "1320", "1380", "1400", "1420", "1440", "1460", "1480", "1500", "1540", "1580", "1600", "1620", "1640", "1660", "1700", "1720", "1780", "1800", "1860", "1900", "1950", "2000", "2020", "2050", "2100", "2120", "2150", "2200", "2260", "2300", "2320", "2350", "2380", "2400", "2420", "2450", "2490", "2500", "2600", "2700", "2800", "2900", "3000", "3100", "3200", "3300", "3400", "3500", "3600", "3700", "3800", "3900", "4000")
  colnames(drlam) <- drlam.colnames
  colnames(drrlam) <- drlam.colnames
  colnames(srlam) <- drlam.colnames
  return (list(metout=metout, soil=soil, shadmet=shadmet, shadsoil=shadsoil, soilmoist=soilmoist, shadmoist=shadmoist, humid=humid, shadhumid=shadhumid, soilpot=soilpot, shadpot=shadpot, plant = plant, shadplant = shadplant,  sunsnow=sunsnow, shdsnow=shdsnow, tcond=tcond, shadtcond=shadtcond, specheat=specheat, shadspecheat=shadspecheat, densit=densit, shaddensit=shaddensit, drlam=drlam, drrlam=drrlam, srlam=srlam))
}
