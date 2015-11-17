#' Implementation of the endotherm model .
#'
#' An implementation of described in
#' Porter, W. P., and M. Kearney. 2009. Size, shape, and the thermal niche of endotherms.
#' Proceedings of the National Academy of Sciences 106:19666-19672.
#' with additional rough calculations for water loss. Note this model is only relevant for
#' conditions where there is no solar radiation and air/ground/sky temperatures are equal
#' @param posture = 4.5, Shape, ratio of long to short axis of a prolate ellipsoid
#' @param mass = 0.5, Body Mass (kg)
#' @param coreT = 37, Core temperature (deg C)
#' @param furdepth = 5, Fur depth (mm)
#' @param O2eff = 0.2, Oxygen extraction efficiency (decimal \%)
#' @param furcond = 0.04, Conductivity of fur (W/Cm)
#' @param airT = 20, Air temperature (deg C)
#' @param windspd = 1, Wind speed (m/s)
#' @param rh = 50, Relative humidity (\%)
#' @return AirTemp Air temperature for calculation (deg C)
#' @return Windspeed Wind speed for calculation (m/s)
#' @return RelHum Relative humidity for calculation (\%)
#' @return Tcore Core temperature for calculation (deg C)
#' @return UCT Upper Critical Temperature (deg C)
#' @return LCT Lower Critical Temperature (deg C)
#' @return Qresp_gph Respiration rate (W)
#' @return Qresp_W Respiration rate (W)
#' @return Qresp_kjph Respiration rate (kJ/h)
#' @return Tskin Skin temperature (deg C)
#' @return Qgen required heat generation (W)
#' @return QgenFinal Required heat generation capped at basal (W)
#' @return mlO2ph Oxygen consumption rate (ml/h)
#' @return PctBasal Metabolic rate (\% of basal)
#' @return status 1=cold, 2=happy, 3=hot
#' @return H2Oloss_W Watts of heat to be dumped
#' @return H2O_gph Water loss rate required to dump heat based on latent heat of vaporization
#' @return massph_percent Percent of body mass lost as water per hour
#' @return timetodeath Time to death (hours) from desiccation (15\% desiccated) if no water to drink
#' @export
#' @examples
#'micro<-micro_global(loc = 'Birdsville, Australia') # run the model with default location and settings
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
ellipsoid <- function(posture = 4.5, mass = 0.5, coreT = 37, furdepth = 5, O2eff = 0.2,
  furcond = 0.04, airT = 20, windspd = 1, rh = 50) {
  if(posture==1){posture=1.01} # avoid divide by zero
  mouseelephant <- 10^(-1.462 + 0.675 * log10(mass * 1000))
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
  visc_air <- (0.000018325*((296.16+120)/(airT + 273.15+120)))*(((airT + 273.15)/296.16)^1.5)
  k_air <- 0.02425 + (0.00007038 * airT)
  den_air <- 101325 / (287.04 * (airT + 273.15))
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
  upcrit <- coreT - (mouseelephant*0.6*Rtotal)
  lowcrit <- coreT - mouseelephant * Rtotal
  Qgen <- (coreT - airT) / Rtotal
  QgenFinal <- ifelse(Qgen < mouseelephant, mouseelephant, Qgen)
  mlO2ph <- QgenFinal / 20.1 * 3600
  Qresp_gph <- (mlO2ph / 0.2094 / O2eff * (((10^( - 7.90298 * (373.16 / (coreT + 273.15) - 1) + 5.02808 * log10(373.16 / (coreT + 273.15)) - (1.3816 * 10^ - 7) * (10^(11.344 * (1 - (coreT + 273.15) / 373.16)) - 1) + (8.1328 * 10^ - 3) * (10^( - 3.49149 * (373.16 / (coreT + 273.15) - 1)) - 1) + log10(1013.246))) * 100 * (100 / 100)) * 0.018016 / (0.998 * 8.31434 * (coreT + 273.15)) - ((10^( - 7.90298 * (373.16 / (airT + 273.15) - 1) + 5.02808 * log10(373.16 / (airT + 273.15)) - (1.3816 * 10^ - 7) * (10^(11.344 * (1 - (airT + 273.15) / 373.16)) - 1) + (8.1328 * 10^ - 3) * (10^( - 3.49149 * (373.16 / (airT + 273.15) - 1)) - 1) + log10(1013.246))) * 100 * (rh / 100)) * 0.018016 / (0.998 * 8.31434 * (airT + 273.15)))) / 1000
  conv_H2O_loss <- 2501200 - 2378.7 * airT
  Qresp_W <- ((Qresp_gph / 3600) * conv_H2O_loss) / 1000
  Qresp_kjph <- Qresp_W / 1000 * 3600
  PctBasal <- QgenFinal / mouseelephant * 100
  status <- ifelse(Qgen > mouseelephant, 1, ifelse(Qgen < 0.6 * mouseelephant, 3, 2)) # 1 = cold, 2 = happy, 3 = stress
  H2Oloss_W <- (Qgen * -1) + mouseelephant
  H2O_gph <- (((H2Oloss_W) * 1000) / conv_H2O_loss) * 3600
  H2O_gph[H2O_gph<0]<-0
  massph_percent <- ifelse(H2O_gph < 0, 0, ((H2O_gph / 1000) / mass) * 100)
  timetodeath <- ifelse(H2O_gph < 0, NA, 1 / (H2O_gph / (mass * 0.15 * 1000)))


  results<-cbind(airT, windspd, rh, coreT, upcrit, lowcrit, Qresp_gph, Qresp_W, Qresp_kjph, skinT, Qgen,
           QgenFinal, mlO2ph, PctBasal, status, H2Oloss_W, H2O_gph,
           massph_percent, timetodeath)
  colnames(results)<-c('AirTemp','WindSpeed','RelHum','Tcore','UCT','LCT','Qresp_gph','Qresp_W',
    'Qresp_kjph','Tskin','Qgen','QgenFinal','mlO2ph','PctBasal','status','H2Oloss_W',
    'H2O_gph','massph_percent','timetodeath')
  return(results)
}
