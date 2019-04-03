#' Ellipsoid endotherm model
#'
#' An implementation of the model described in
#' Porter, W. P., and M. Kearney. 2009. Size, shape, and the thermal niche of endotherms.
#' Proceedings of the National Academy of Sciences 106:19666-19672.
#' with additional rough calculations for water loss. Note this model is only relevant for
#' conditions where there is no solar radiation and air/ground/sky temperatures are equal.
#' Originally coded into R from Excel by John Baumgartner.
#' Requires the NicheMapR functions WETAIR.rh, DRYAIR and VAPPRS
#' Modified to allow use with rasters 2 Aug 2018
#' Modified to allow user-specified basal metabolic rates and Q10 effects 2 Aug 2018
#' @param posture = 4.5, Shape, ratio of long to short axis of a prolate ellipsoid
#' @param mass = 0.5, Body Mass (kg)
#' @param coreT = 37, Core temperature (deg C)
#' @param furdepth = 5, Fur depth (mm)
#' @param furcond = 0.04, Conductivity of fur (W/Cm)
#' @param O2eff = 0.2, Oxygen extraction efficiency (decimal \%)
#' @param stress = 0.6, Fraction of basal metabolic rate at which evaporative water loss is required (decimal \%)
#' @param airT = 20, Air temperature (deg C)
#' @param windspd = 1, Wind speed (m/s)
#' @param rh = 50, Relative humidity (\%)
#' @param basal = NA, user specified basal metabolic rate (W)
#' @param basmult = 1, multiplier to adjust mouse-elephant predicted basal metabolic rate
#' @return
#' \itemize{
#' \item 1 AirTemp - Air temperature for calculation (deg C)
#' \item 2 Windspeed - Wind speed for calculation (m/s)
#' \item 3 RelHum - Relative humidity for calculation (\%)
#' \item 4 Tcore - Core temperature for calculation (deg C)
#' \item 5 UCT - Upper Critical Temperature (deg C)
#' \item 6 LCT - Lower Critical Temperature (deg C)
#' \item 7 Qresp_gph - Respiration rate (W)
#' \item 8 Qresp_W - Respiration rate (W)
#' \item 9 Qresp_kjph - Respiration rate (kJ/h)
#' \item 10 Tskin - Skin temperature (deg C)
#' \item 11 Qgen - required heat generation (W)
#' \item 12 QgenFinal - Required heat generation capped at basal (W)
#' \item 13 mlO2ph - Oxygen consumption rate (ml/h)
#' \item 14 PctBasal - Metabolic rate (\% of basal)
#' \item 15 status - 1=cold, 2=happy, 3=hot
#' \item 16 H2Oloss_W - Watts of heat to be dumped
#' \item 17 H2O_gph - Water loss rate required to dump heat based on latent heat of vaporization
#' \item 18 massph_percent - Percent of body mass lost as water per hour
#' \item 19 timetodeath - Time to death (hours) from desiccation (15\% desiccated) if no water to drink
#' }
#' @details
#' Note that the parameter 'stress' is a fudge factor that accounts for physiological thermoregulatory
#' adjustments during heat stress (e.g. core temperature increases, changes in flesh conductivity, changes
#' in posture) that are not captured dynamically in this function (but some of these could be modelled
#' dynamically on the basis of this function).
#' @export
#' @examples
#'
#'# compute metabolic and water loss rate as a function of air temperature
#'endo<-ellipsoid(airT=seq(5,45), windspd = 0.1)
#'endo<-as.data.frame(endo)
#'plot(endo$AirTemp,endo$QgenFinal,type='l', xlab = "Air Temperature (deg C)", ylab = "Metabolic Rate (W)", main="Metabolic Rate vs. Air Temperaure")
#'plot(endo$AirTemp,endo$H2O_gph,type='l', xlab = "Air Temperature (deg C)", ylab = "Water Loss Rate (g/h)", main="Water Loss Rate vs. Air Temperaure")
#'
#'# compute thermoneutral zone as a function of body mass
#'masses<-10^seq(-3,2,0.5) # log sequence of masses
#'endo<-ellipsoid(mass = masses)
#'endo<-as.data.frame(endo)
#'ylims<-c(min(endo$UCT,endo$LCT),max(endo$UCT,endo$LCT))
#'plot(masses,endo$UCT, col='red',type='l', ylim=ylims, xlab="body mass (kg)",ylab="temperature (deg C)", main = "Upper and Lower Critical Temperatures vs Mass")
#'points(masses,endo$LCT, type='l', col='blue')
#'
#'micro<-micro_global(loc = c(139.5, -25.9)) # run the microclimate model at Birdsville with default settings
#'
#'metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#'shadmet<-as.data.frame(micro$shadmet) # above ground microclimatic conditions, max shade
#'soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
#'shadsoil<-as.data.frame(micro$shadsoil) # soil temperatures, maximum shade
#'
#'# append dates
#'days<-rep(seq(1,12),24)
#'days<-days[order(days)]
#'dates<-days+metout$TIME/60/24-1 # dates for hourly output
#'dates2<-seq(1,12,1) # dates for daily output
#'metout<-cbind(dates,metout)
#'soil<-cbind(dates,soil)
#'shadmet<-cbind(dates,shadmet)
#'shadsoil<-cbind(dates,shadsoil)
#'
#'endo<-cbind(metout[,1:3],ellipsoid(airT = shadmet$TALOC, windspd = shadmet$VLOC, rh = shadmet$RHLOC))
#'
#'with(endo,{plot(H2O_gph ~ dates,xlab = "Date and Time", ylab = "Water Loss Rate (g/h)"
#', type = "l",main=paste("Evaporative Water Loss",sep=""))})
#'with(endo,{plot(QgenFinal ~ dates,xlab = "Date and Time", ylab = "Metabolic Rate (W)"
#' , type = "l",main=paste("Metabolic Heat Generation",sep=""))})
#'with(endo,{plot(PctBasal ~ dates,xlab = "Date and Time", ylab = "Metabolic Rate (% of basal)"
#', type = "l",main=paste("Metabolic Heat Generation",sep=""))})
#'with(endo,{plot(Tskin ~ dates,xlab = "Date and Time", ylab = "Skin and Core Temperature (deg C)"
#', type = "l",main=paste("Skin and Core Temperature",sep=""),ylim=c(min(cbind(endo$Tcore,endo$Tskin)),
#'  max(cbind(endo$Tcore,endo$Tskin))))})
#'with(endo,{points(Tcore ~ dates,xlab = "Date and Time",lty=2, type = "l")})
#'with(endo,{plot(timetodeath ~ dates,xlab = "Date and Time", ylab = "Time to Death (h)"
#', type = "l",main=paste("Time to Death by Desiccation",sep=""), ylim=c(0,24))})
ellipsoid <- function(posture = 4.5, mass = 0.5, coreT = 37, furdepth = 5, furcond = 0.04,
  O2eff = 0.2, stress = 0.6, airT = 20, windspd = 1, rh = 50, Q10 = 3, basal = NA, basmult = 1) {
  posture[posture==1]<-1.01 # avoid divide by zero
  if(class(basal)=="logical"){ # this checks if basal is set to 'NA'
   mouseelephant <- 10^(-1.462 + 0.675 * log10(mass * 1000)) * basmult
   basal <- mouseelephant * Q10 ^ ((coreT - 37) / 10) # Q10 correction
  }
  a_coef <- 0.6
  b_coef <- 0.5
  sp_heat_air <- 1005.8
  volume <- mass / 1000
  b <- ((3 * volume) / (4 * 3.14159 * posture))^0.333
  c <- b
  a <- b * posture
  k_body <- 0.5 + (6.14 * b) + 0.439
  numerator <- a^2*b^2*c^2
  denominator <- a^2 * b^2 + a^2 * c^2 + b^2 * c^2
  Fraction <- numerator/denominator
  Rbody <- Fraction / (2 * k_body * volume)
  ao <- b*posture + furdepth/1000
  bo <- b + furdepth/1000
  co <- c + furdepth/1000
  Ecc_outr <- sqrt(ao^2 - co^2) / ao
  Aouter <- 2 * pi * bo^2 + 2 * pi * ((ao*bo)/Ecc_outr) * asin(Ecc_outr)
  Rinsul <- (bo - b)/(furcond * Aouter)
  dryair=DRYAIR(db=airT)
  visc_air <- dryair$visdyn
  k_air <- dryair$thcond
  den_air <- dryair$densty
  volcheck <- (4/3) * 3.14159 * a * b * c
  CharDimens <- volcheck^0.333
  Eccentricity <- sqrt(a^2 - c^2) / a
  area <- 2 * pi * b^2 + 2 * pi * ((a * b) / Eccentricity) * asin(Eccentricity)
  Re_number <- den_air * windspd * CharDimens/visc_air
  Pr_number <- (visc_air * sp_heat_air) / k_air
  q3p_num <- 2 * area * k_body * k_air * (2 + a_coef * (Re_number^b_coef)*Pr_number^0.333)*(coreT - airT)
  q3p_denom <- 2 * k_body * CharDimens * volcheck + area * Fraction * k_air * (2 + (a_coef * Re_number^b_coef) * Pr_number^0.333)
  g_in_air <- q3p_num/q3p_denom
  skinT <- coreT - (g_in_air * Fraction) / (2 * k_body)
  Gr_number <- abs(((den_air^2)*(1/(airT+273.15))*9.80665*(CharDimens^3)*(skinT - airT))/(visc_air^2))
  Nufree <- 2 + 0.6 * ((Gr_number^0.25) * (Pr_number^0.333))
  Nuforced <- 0.37 * Re_number^0.6
  Nutotal <- (Nufree^3 + Nuforced^3)^(1/3)
  hc <- Nutotal * k_air/CharDimens
  Rconv <- 1 / (hc * Aouter)
  Rrad <- 1 / (4 * Aouter * 1 * 0.95 * (5.7 * 10^-8) * (airT + 273.15)^3)
  Rtotal <- Rbody + Rinsul + (Rconv*Rrad)/(Rconv + Rrad)
  upcrit <- coreT - (basal*stress*Rtotal)
  lowcrit <- coreT - basal * Rtotal
  Qgen <- (coreT - airT) / Rtotal
  QgenFinal <- Qgen
  if(length(basal) == 1){
   QgenFinal[QgenFinal<basal]<-basal
  }else{
   QgenFinal[QgenFinal<basal]<-basal[QgenFinal<basal]
  }
  mlO2ph <- QgenFinal / 20.1 * 3600
  esat <- VAPPRS(coreT)
  Qresp_gph <- (mlO2ph / 0.2094 / O2eff) * (WETAIR.rh(db = coreT, rh = 100)$vd - WETAIR.rh(db = airT, rh = rh)$vd) / 1000
  conv_H2O_loss <- 2501200 - 2378.7 * airT
  Qresp_W <- ((Qresp_gph / 3600) * conv_H2O_loss) / 1000
  Qresp_kjph <- Qresp_W / 1000 * 3600
  PctBasal <- QgenFinal / basal * 100
  status<-Qgen
  status[status>basal]<--100000000
  status[status<stress * basal]<--300000000
  status[status<100000000*-1]<--200000000
  status[status==100000000*-1]<-1
  status[status==300000000*-1]<-3
  status[status==200000000*-1]<-2
  H2Oloss_W <- (Qgen * -1) + basal
  H2O_gph <- (((H2Oloss_W) * 1000) / conv_H2O_loss) * 3600
  H2O_gph[H2O_gph<0]<-0
  massph_percent<-H2O_gph
  massph_percent[massph_percent<0]<-0
  timetodeath<-massph_percent
  massph_percent[massph_percent!=0]<-((H2O_gph[massph_percent != 0] / 1000) / mass) * 100
  timetodeath[timetodeath!=0]<-1 / (H2O_gph[timetodeath!=0] / (mass * 0.15 * 1000))

  # check if raster output or table
  if(class(airT)[[1]] == 'RasterLayer'){
   results <- list(upcrit = upcrit, lowcrit = lowcrit, Qresp_gph = Qresp_gph, Qresp_W = Qresp_W,
     Qresp_kjph = Qresp_kjph, skinT = skinT, Qgen = Qgen, QgenFinal = QgenFinal, mlO2ph = mlO2ph,
     PctBasal = PctBasal, status = status, H2Oloss_W = H2Oloss_W, H2O_gph = H2O_gph,
     massph_percent = massph_percent, timetodeath = timetodeath)
  }else{
   results<-cbind(airT, windspd, rh, coreT, upcrit, lowcrit, Qresp_gph, Qresp_W, Qresp_kjph, skinT, Qgen,
           QgenFinal, mlO2ph, PctBasal, status, H2Oloss_W, H2O_gph,
           massph_percent, timetodeath)
  colnames(results)<-c('AirTemp','WindSpeed','RelHum','Tcore','UCT','LCT','Qresp_gph','Qresp_W',
    'Qresp_kjph','Tskin','Qgen','QgenFinal','mlO2ph','PctBasal','status','H2Oloss_W',
    'H2O_gph','massph_percent','timetodeath')
  }
  return(results)
}
