#' Thermoregulatory model using a transient heat budget
#'
#' This model uses the transient heat budget models (i.e. accounting for heat storage and
#' hence lag-effects of body mass) to simulate thermoregulatory behaviour of a diurnally active
#' ectotherm. It uses a set of events to break out of the ordinary differential
#' equation solver of the transient heat budget (onelump_var or twolump functions) to simulate
#' thermoregulation around setpoints, specifically the transition from sitting in the shade
#' overnight to basking in the sun perpendicular to the sun's rays (T_B_min), from basking to
#' foraging in the open in an 'average' posture (T_F_min), from foraging in the open to moving into
#' the shade (T_F_max), and then transitioning to basking in the afternoon and finally retreating
#' to the shade for the evening. It also computes the operative temperature (Te) of a
#' non-thermoregulating animal in the open, i.e. the steady state body temperature the animal would
#' come to if it had zero heat capacity, and the body temperature of a non thermoregulating animal
#' in the open accounting for the effect of thermal mass. The function computes a series of
#' summary statistics about the length of basking and foraging bouts, activity time, time spent
#' above high temperature thresholds (T_F_max and CT_max), and the extremes of body temperature.
#' @encoding UTF-8
#' @param Tc_init = Tairf(1), initial temperature (°C) Organism shape, 0-5, Determines whether standard or custom shapes/surface area/volume relationships are used: 0=plate, 1=cyl, 2=ellips, 3=lizard (desert iguana), 4=frog (leopard frog), 5=custom (see details)
#' @param Ts_init = Tc_init + 0.1, initial shell temperature (°C)
#' @param To_init = Tc_init + 0.2, initial surface temperature (°C)
#' @param Ww_g = 500, weight (g)
#' @param T_F_min = 33, # minimum foraging body temperature threshold (°C)
#' @param T_F_max = 38, # maximum foraging body temperature threshold (°C)
#' @param T_B_min = 25, # basking body temperature threshold (°C)
#' @param CT_max = 43, # critical thermal maximum (°C)
#' @param rho_body = 932, animal density (kg/m3)
#' @param lump = 1, one lump (1) or two lump (2) model?
#' @param x_shell = 0.001, shell thickness, m
#' @param q = 0, metabolic rate (W/m3)
#' @param c_body = 3073, specific heat of flesh (J/kg-°C)
#' @param c_body_inner = c_body, Specific heat of flesh J/(kg-°C)
#' @param c_body_outer = c_body, Specific heat of outer shell J/(kg-°C)
#' @param k_flesh = 0.5, conductivity of flesh (W/m-K)
#' @param k_inner = k_flesh, Thermal conductivity of inner shell (W/m-K, range: 0.412-2.8)
#' @param k_outer = k_flesh, Thermal conductivity of outer shell (W/m-K, range: 0.412-2.8)
#' @param emis = 0.95, emissivity of skin (-)
#' @param alpha = 0.85, animal solar absorptivity (-)
#' @param geom = 2, Organism shape, 0-5, Determines whether standard or custom shapes/surface area/volume relationships are used: 0=plate, 1=cyl, 2=ellips, 3=lizard (desert iguana), 4=frog (leopard frog), 5=custom (see parameter 'shape_coeffs')
#' @param shape_b = 1/5, Proportionality factor (-) for going from volume to area, represents ratio of width:height for a plate, length:diameter for cylinder, b axis:a axis for ellipsoid
#' @param shape_c = 1/5, Proportionality factor (-) for going from volume to area, represents ratio of length:height for a plate, c axis:a axis for ellipsoid
#' @param shape_coefs = c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743), Custom surface area coefficients. Operates if posture = 5, and consists of 4 pairs of values representing the parameters a and b of a relationship AREA=a*Ww_g^b, where AREA is in cm2 and Ww_g is in g. The first pair are a and b for total surface area, then a and b for ventral area, then for sillhouette area normal to the sun, then sillhouette area perpendicular to the sun
#' @param posture = 'n', pointing normal 'n' or parallel 'p' to the sun's rays, or average 'a'?
#' @param orient = 1, does the object orient toward the sun? (0,1)
#' @param fatosk = 0.4, solar configuration factor to sky (-)
#' @param fatosb = 0.4, solar configuration factor to substrate (-)
#' @param alpha_sub = 0.2, substrate solar reflectivity, decimal percent
#' @param pdif = 0.1, proportion of solar energy that is diffuse (rather than direct beam)
#' @param shade = 90, maximum shade level (\%)
#' @param metout = metout, aboveground minimum shade microclimate output table from NicheMapR's microclimate model
#' @param shadmet = shadmet, metout, aboveground maximum shademicroclimate output table from NicheMapR's microclimate model
#' @param soil = soil, minimum shade soil temperature output table from NicheMapR's microclimate model
#' @param shadsoil = shadsoil, maximum shade soil temperature output table from NicheMapR's microclimate model
#' @param press = 101325, air pressure (Pa)
#'
#' \strong{Outputs:}
#'
#' day_results variables:
#' \itemize{
#' \item 1 time - time (hours)
#' \item 2 hour - hour of day (rounded to nearest hour)
#' \item 3 T_air_shd - shaded air temperature at animal height  (°C)
#' \item 4 Tb - body temperature of thermoregulating animal (°C)
#' \item 5 Tb_final - steady state temperature in current environment (°C)
#' \item 6 Tb_open - body temperature of animal staying in the open (°C)
#' \item 7 Te_open  - operative temperature (zero heat capacity) of animal in open (°C)
#' \item 8 time_constant - time constant of animal (minutes)
#' \item 9 dTb_dt - rate of change of body temperature (°C / minute)
#' \item 10 posture - posture of animal ('n' = normal to sun's rays, 'p' = parallel to sun's rays , 'f' = foraging posture)
#' \item 11 active - is the animal inactive (0) or active (1)?
#' \item 12 state - activity state (0 = sleeping, 1 = basking, 2 = foraging in open, 3 = cooling in shade)
#' \item 13 mrate - total metabolic rate for time step (J)
#' }
#'
#' act_window variables:
#' \itemize{
#' \item 1 time - time (hours)
#' \item 2 forage_sun - total foraging time in sun for this hour (minutes)
#' \item 3 max_bout_shd - maximum foraging bout in sun for this hour (minutes)
#' \item 4 forage_shd - total foraging time in shade for this hour (minutes)
#' \item 5 max_bout_shd - maximum foraging bout in shade for this hour (minutes)
#' }
#'
#' sum_stats variables:
#' \itemize{
#' \item 1 Ww_g - wet weight of animal (g)
#' \item 2 T_F_min - minimum foraging temperature (°C)
#' \item 3 T_F_max - maximum foraging temperature (°C)
#' \item 4 max_bout_sun - longest foraging bout in sun (minutes)
#' \item 5 max_bout_shd - longest foraging bout in shade (minutes)
#' \item 6 sum_activity_sun - total activity in sun (minutes)
#' \item 7 sum_activity_shd - total activity in shade (minutes)
#' \item 8 bouts_sun - total number of foraging bouts in sun (#)
#' \item 9 bouts_shd - total number of foraging bouts in shade (#)
#' \item 10 morning_bask - morning basking time (minutes)
#' \item 11 morning_forage_sun - first foraging bout length in sun (minutes)
#' \item 12 morning_forage_shd - first foraging bout length in shade (minutes)
#' \item 13 first_midday_bout_sun - second foraging bout length in sun (minutes)
#' \item 14 first_midday_bout_shd - second foraging bout length in shade (minutes)
#' \item 15 mean_midday_bout_shd - average length of foraging bouts between the first and last foraging bouts in sun (mins)
#' \item 16 mean_midday_bout_sun - average length of foraging bouts between the first and last foraging bouts in shade (mins)
#' \item 17 afternoon_forage_sun - last foraging bout length in sun (minutes)
#' \item 18 afternoon_forage_shd - last foraging bout length in shade (minutes)
#' \item 19 mrate_sum - cumulative metabolic rate (kJ)
#' \item 20 mrate_sum_inactive_sun - cumulative metabolic rate while inactive, sun forager (kJ)
#' \item 21 mrate_sum_inactive_shd - cumulative metabolic rate while inactive, shade forager (kJ)
#' \item 22 mrate_sum_active_sun - cumulative metabolic rate while active in sun (kJ)
#' \item 23 mrate_sum_active_shd - cumulative metabolic rate while active in shade (kJ)
#' \item 24 T_F_max_time_Te - time operative temperature in open spent above maximum foraging temperature (minutes)
#' \item 25 CT_max_time_Te - time operative temperature in open spent above critical thermal maximum (minutes)
#' \item 26 T_F_max_time_Tb_open - time body temperature in open spent above maximum foraging temperature (minutes)
#' \item 27 CT_max_time_Tb_open - time body temperature in open spent above critical thermal maximum (minutes)
#' \item 28 T_F_maxtime - time body temperature of thermoregulating animal spent above maximum foraging temperature (minutes)
#' \item 29 CT_max_time - time body temperature of thermoregulating animal spent above critical thermal maximum (minutes)
#' \item 30 max_Tb_open - maximum body temperature in open (°C)
#' \item 31 min_Tb_open - minimum body temperature in open (°C)
#' \item 32 max_Te - maximum operative temperature in open (°C)
#' \item 33 min_Te - minimum operative temperature in open (°C)
#' \item 34 max_Tb - maximum body temperature of thermoregulating animal (°C)
#' \item 35 min_Tb - minimum body temperature of thermoregulating animal (°C)
#' }
#'
#' @examples
#' library(NicheMapR)
#'
#' # define animal biophysical functional traits
#' Ww_g <- 500 # wet weight (g)
#' Usrhyt <- 0.05 # height of animal (mid-point) above ground (m)
#' alpha <- 0.85 # solar absorptivity (-)
#' T_F_min <- 33 # minimum foraging Tb (deg C)
#' T_F_max <- 43 # maximum foraging Tb (deg C)
#' T_B_min <- 18 # basking Tb, moving from shade to sun (deg C)
#' CT_max <- 48 # critical thermal maximum (deg C)
#' shape_b <- 1/5 # shape coefficient a, -
#' shape_c <- 1/5 # shape coefficient b, -
#' rho_body <- 1000 # animal density, kg/m3
#' c_body <- 3073 # heat capacity (J/kg-C)
#' q <- 0 # metabolic rate, W/m3
#' k_flesh <- 0.5 # thermal conductivity of flesh, W/mK
#' geom <- 2 # shape, -
#'
#' # get microclimate data
#' loc <- c(130, -25)
#' maxshade <- 90
#' micro <- micro_global(loc = loc, maxshade = maxshade, Usrhyt = Usrhyt) # run the model with default location and settings
#' metout <- as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#' soil <- as.data.frame(micro$soil) # soil temperatures, minimum shade
#' shadmet <- as.data.frame(micro$shadmet) # above ground microclimatic conditions, min shade
#' shadsoil <- as.data.frame(micro$shadsoil) # soil temperatures, minimum shade
#' # get air pressure
#' elevation <- micro$elev
#' press <- 101325 * ((1 - (0.0065 * elevation / 288)) ^ (1 / 0.190284))
#'
#' mons <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
#' DOYs <- unique(metout$DOY)
#'
#' # loop through each month and run transient model with behaviour
#' for(i in 1:12){
#'
#'   # subset current month
#'   metout_in <- subset(metout, DOY == DOYs[i])
#'   shadmet_in <- subset(shadmet, DOY == DOYs[i])
#'   soil_in <- subset(soil, DOY == DOYs[i])
#'   shadsoil_in <- subset(shadmet, DOY == DOYs[i])
#'
#'   # run transient behavioural simulation
#'   trans <- trans_behav(Ww_g = Ww_g, alpha = alpha, T_F_min = T_F_min, T_F_max = T_F_max,
#'                        CT_max = CT_max, T_B_min = T_B_min, geom = geom, shape_b = shape_b, shape_c = shape_c,
#'                        rho_body = rho_body, k_flesh = k_flesh, q = q, lump = 1,
#'                        metout = metout_in, shadmet = shadmet_in, soil = soil_in, shadsoil = shadsoil_in,
#'                        press = press, alpha_sub = 1 - micro$REFL, shade = micro$maxshade)
#'
#'   results <- as.data.frame(trans$day_results)
#'   sum_stats <- as.data.frame(trans$sum_stats)
#'   act_window <- as.data.frame(trans$act_window)
#'
#'   # collate
#'   if(i == 1){
#'     all_act_window <- act_window
#'   }else{
#'     all_act_window <- rbind(all_act_window, act_window)
#'   }
#'
#'   results$hours <- results$time / 3600
#'
#'   # plot hourly results for the current day
#'   plot(results$Tb_open ~ results$hours, type = 'l', ylim = c(-10, 80), col = 'grey', xaxs = 'i', ylab = "temperature, deg C", xlab = "time",
#'        main = paste0(if(length(loc) == 2){paste("lon", loc[1], "lat", loc[2])}else{loc}, ", ", mons[i], ", ", Ww_g," g"), xlim = c(0, 23))
#'   grid(nx = 23, ny = 0, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
#'   abline(T_F_max, 0, col = 'red', lty = 2)
#'   abline(T_F_min, 0, col = 'light blue', lty = 2)
#'   abline(CT_max, 0, col = 'red')
#'   points(results$T_air_shd ~ results$hours, type = 'l', col = 'blue')
#'   points(results$Tb ~ results$hours, type = 'l', col = 'orange', lty = 1, lwd = 2)
#'   text(3, 60, paste0("bouts ", round(sum_stats$bouts_sun, 0)), cex = 1)
#'   text(3, 65, paste0("maximum bout ", round(sum_stats$max_foraging_bout_sun / 60, 1), " hrs"), cex = 1)
#'   text(3, 70, paste0("total activity ", round(sum_stats$sum_activity_sun / 60, 1), " hrs"), cex = 1)
#' }
#'
#' # make seasonal activity plot
#' all_act_window$ZEN <- metout$ZEN
#' all_act_window$DOY <- metout$DOY
#' foraging<-subset(all_act_window, forage_sun > 0)
#' night<-subset(all_act_window, ZEN==90)
#' with(night, plot(time ~ DOY, pch=15, cex = 2, xlim = c(1, 365), col = 'dark blue', xlab = 'day of year', ylab = 'hour of day', main = "Seasonal Activity Plot, Sun"))
#' with(foraging, points(time ~ DOY, pch = 15, cex = forage_sun / 30, col = 'orange'))
#' foraging<-subset(all_act_window, forage_shd > 0)
#' with(night, plot(time ~ DOY, pch=15, cex = 2, xlim = c(1, 365), col = 'dark blue', xlab = 'day of year', ylab = 'hour of day', main = "Seasonal Activity Plot, Shade"))
#' with(foraging, points(time ~ DOY, pch = 15, cex = forage_shd / 30, col = 'orange'))
#'
#' mtext(text =  paste0('Seasonal Activity Plot, ', if(length(loc) == 2){paste("lon", loc[1], "lat", loc[2])}else{loc}, " ", Ww_g," g"), outer = TRUE, side = 3, line = 0)
#' @export
trans_behav <- function(Tc_init = rep(20, 60),
                        Ts_init = Tc_init + 0.1,
                        To_init = Tc_init + 0.2,
                        Ww_g = 500,
                        T_F_min = 33,
                        T_F_max = 38,
                        T_B_min = 25,
                        CT_max = 43,
                        rho_body = 932,
                        x_shell = 0.001,
                        lump = 1,
                        q =0,
                        c_body = 3073,
                        c_body_inner = c_body,
                        c_body_outer = c_body,
                        k_flesh = 0.5,
                        k_inner = k_flesh,
                        k_outer = k_flesh,
                        emis = 0.95,
                        alpha = 0.85,
                        geom = 2,
                        shape_b = 1/5,
                        shape_c = 1/5,
                        shape_coefs = c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743),
                        posture = 'n',
                        orient = 1,
                        fatosk = 0.4,
                        fatosb = 0.4,
                        alpha_sub = 0.2,
                        pdif = 0.1,
                        shade = 90,
                        metout = metout,
                        shadmet = shadmet,
                        soil = soil,
                        shadsoil = shadsoil,
                        press = 101325) {

  if (exists("day_results"))
  {
    rm(day_results)
  }  # clear the results, if any already in the memory

  # single site transient analysis using user-input microclimate

  require(deSolve)

  # combine relevant variables
  micro_sun <- cbind(metout[,1:4], metout[,7], soil[,3], metout[,12:14], metout[,5])
  names <- c('DOY', 'TIME', 'TALOC', 'TA2m', 'VLOC', 'TS', 'ZEN', 'SOLR', 'TSKYC', 'RHLOC')
  colnames(micro_sun) <- names
  micro_shd <- cbind(shadmet[,1:4], shadmet[,7], shadsoil[,3], shadmet[,12:14], shadmet[,5])
  colnames(micro_shd) <- names
  time<-seq(0,60*24,60) # 60 minute intervals from microclimate output
  time<-time*60 # minutes to seconds

  Qsolf_sun <- approxfun(time, c(micro_sun[,8],(micro_sun[1,8]+micro_sun[24,8])/2), rule = 2)
  Tradf_sun <- approxfun(time, rowMeans(cbind(c(micro_sun[,6],(micro_sun[1,6]+micro_sun[24,6])/24),c(micro_sun[,9],(micro_sun[1,9]+micro_sun[24,9])/24)),na.rm=TRUE), rule = 2)
  Qsolf_shd <- approxfun(time, c(micro_shd[,8],(micro_shd[1,8]+micro_shd[24,8])/2)*(1-shade/100), rule = 2)
  Tradf_shd <- approxfun(time, rowMeans(cbind(c(micro_shd[,6],(micro_shd[1,6]+micro_shd[24,6])/24),c(micro_shd[,9],(micro_shd[1,9]+micro_shd[24,9])/24)),na.rm=TRUE), rule = 2)
  velf_sun <- approxfun(time, c(micro_sun[,5],(micro_sun[1,5]+micro_sun[24,5])/2), rule = 2)
  velf_shd <- approxfun(time, c(micro_shd[,5],(micro_shd[1,5]+micro_shd[24,5])/2), rule = 2)
  Tairf_sun <- approxfun(time, c(micro_sun[,3],(micro_sun[1,3]+micro_sun[24,3])/2), rule = 2)
  Tairf_shd <- approxfun(time, c(micro_shd[,3],(micro_shd[1,3]+micro_shd[24,3])/2), rule = 2)
  Zenf <- approxfun(time, c(micro_sun[,7],90), rule = 2)

  times_sec <- seq(0, 3600 * 24, 3600)  # hours of day in seconds
  times <- seq(0, 3600 * 24, 10)  # sequence of seconds for a day
  times <- times[1:(length(times) - 1)]
  hours <- times/3600
  times_orig <- times

  # events

  emerge <- function(t, y, pars) {
    # if sun is up and body temperature greater than threshold for basking, then trigger 'emerge' event
    if (Zenf(t) != 90 & y[1] > T_B_min) {
      y[1] <- 0
    }
    return(y)
  }
  retreat <- function(t, y, pars) {
    # if sun is down or body temperature is lower than threshold for basking, then trigger 'retreat' event
    if (Zenf(t) == 90 | y[1] < T_B_min) {
      y[1] <- 0
    }
    return(y)
  }
  shuttle <- function(t, y, pars) {
    # if temperature exceeds voluntary max Tb or is lower than voluntary minimum, trigger 'shuttle' event
    if (y[1] >= T_F_max | y < T_F_min) {
      y[1] <- 0
    }
    return(y)
  }
  forage <- function(t, y, pars) {
    # if Tb exceeds voluntary min foraging temp, or reaches the shaded air temp (plus a bit, one degree, because
    # otherwise they chase a moving target in the afternoon and may never come out) in the case of an animal cooling
    # in the shade when the shaded air temp is higher than the T_F_min, trigger 'forage' event (have to also make sure
    # that shaded air temp isn't approaching the T_F_max)
    if (lump == 1) {
      return(y[1] - max(T_F_min, if (Tairf_shd(t) > T_F_max - 2) {
        0
      } else {
        Tairf_shd(t) + 1
      }))
    }
    if (lump == 2) {
      return(c(y[1] - max(T_F_min, if (Tairf_shd(t) > T_F_max - 2) {
        0
      } else {
        Tairf_shd(t) + 1
      }), y[2:3]))
    }
  }
  eventfun <- function(t, y, pars) {
    if (lump == 1) {
      return(y <- 1)
    }
    if (lump == 2) {
      return(y <- c(1, 1, 1))
    }
  }

  # behaviours

  morning <- function() {
    if (lump == 1) {
      Tbs_ode <- as.data.frame(ode(y = Tc_init, times = subtime, func = onelump_var, parms = indata, events = list(func = eventfun,
                                                                                                                   root = TRUE, terminalroot = 1), rootfun = emerge, method = "lsoda"))
      colnames(Tbs_ode) <- c("time", "Tb", "Tcfinal", "tau", "dTc")
    }
    if (lump == 2) {
      Tbs_ode <- as.data.frame(ode(y = c(Tc_init, Ts_init, To_init), times = subtime, func = twolump, parms = indata,
                                   events = list(func = eventfun, root = TRUE, terminalroot = 1), rootfun = emerge, method = "lsoda"))
      colnames(Tbs_ode) <- c("time", "Tb", "Ts", "To", "Tcf")
    }
    return(Tbs_ode)
  }

  afternoon <- function() {
    if (lump == 1) {
      Tbs_ode <- as.data.frame(ode(y = Tc_init, times = subtime, func = onelump_var, parms = indata, events = list(func = eventfun,
                                                                                                                   root = TRUE, terminalroot = 1), rootfun = retreat, method = "lsoda"))
      colnames(Tbs_ode) <- c("time", "Tb", "Tcfinal", "tau", "dTc")
    }
    if (lump == 2) {
      Tbs_ode <- as.data.frame(ode(y = c(Tc_init, Ts_init, To_init), times = subtime, func = twolump, parms = indata,
                                   events = list(func = eventfun, root = TRUE, terminalroot = 1), rootfun = retreat, method = "lsoda"))
      colnames(Tbs_ode) <- c("time", "Tb", "Ts", "To", "Tcf")
    }
    return(Tbs_ode)
  }

  warming <- function() {
    if (lump == 1) {
      Tbs_ode <- as.data.frame(ode(y = Tc_init, times = subtime, func = onelump_var, parms = indata, events = list(func = eventfun,
                                                                                                                   root = TRUE, terminalroot = 1), rootfun = shuttle, method = "lsoda"))
      colnames(Tbs_ode) <- c("time", "Tb", "Tcfinal", "tau", "dTc")
    }
    if (lump == 2) {
      Tbs_ode <- as.data.frame(ode(y = c(Tc_init, Ts_init, To_init), times = subtime, func = twolump, parms = indata,
                                   events = list(func = eventfun, root = TRUE, terminalroot = 1), rootfun = shuttle, method = "lsoda"))
      colnames(Tbs_ode) <- c("time", "Tb", "Ts", "To", "Tcf")
    }
    return(Tbs_ode)
  }

  cooling <- function() {
    if (lump == 1) {
      Tbs_ode <- as.data.frame(ode(y = Tc_init, times = subtime, func = onelump_var, parms = indata, events = list(func = eventfun,
                                                                                                                   root = TRUE, terminalroot = 1), rootfun = forage, method = "lsoda"))
      colnames(Tbs_ode) <- c("time", "Tb", "Tcfinal", "tau", "dTc")
    }
    if (lump == 2) {
      Tbs_ode <- as.data.frame(ode(y = c(Tc_init, Ts_init, To_init), times = subtime, func = twolump, parms = indata,
                                   events = list(func = eventfun, root = TRUE, terminalroot = 1), rootfun = forage, method = "lsoda"))
      colnames(Tbs_ode) <- c("time", "Tb", "Ts", "To", "Tcf")
    }
    return(Tbs_ode)
  }
  cat("computing transient heat budget \n")

  # initialise
  posture <- "n"  # pointing normal 'n' or parallel 'p' to the sun's rays, or average 'a'?
  Tc_init <- Tairf_shd(0)  # start with Tb at shaded air temp
  Ts_init <- Tc_init
  To_init <- Tc_init
  Tairf <- Tairf_sun
  Qsolf <- Qsolf_sun
  Tradf <- Tradf_sun
  velf <- velf_sun

  # first get Te and Tb in open
  if (lump == 1) {
    indata <- list(alpha = alpha, emis = emis, alpha_sub = alpha_sub, press = press, Ww_g = Ww_g, c_body = c_body, rho_body = rho_body,
                   q = q, k_flesh = k_flesh, geom = geom, posture = posture, orient = orient, shape_b = shape_b, shape_c = shape_c,
                   shape_coefs = shape_coefs, pdif = pdif, fatosk = fatosk, fatosb = fatosb, Tairf = Tairf, velf = velf,
                   Qsolf = Qsolf, Tradf = Tradf, Zenf = Zenf)
    indata$c_body = 0.1  #specific heat of flesh, J/(kg.C)
  }
  if (lump == 2) {
    indata <- list(Ww_g = Ww_g, x_shell = x_shell, geom = geom, k_inner = k_inner, k_outer = k_outer, q = q,
                   c_body_inner = c_body_inner, c_body_outer = c_body_outer, emis = emis, rho_body_body = rho_body_body, alpha = alpha, shape_coefs = shape_coefs,
                   shape_b = shape_b, shape_c = shape_c, posture = posture, orient = orient, fatosk = fatosk, fatosb = fatosb,
                   alpha_sub = alpha_sub, pdif = pdif, press = press)
    indata$c_body_inner = 0.1  #specific heat of flesh, J/(kg.C)
    indata$c_body_outer = 0.1  #specific heat of flesh, J/(kg.C)
  }
  if (lump == 1) {
    Te <- try(ode(y = Tc_init, times = times, func = onelump_var, parms = indata)[, 2])
  }
  if (lump == 2) {
    Te <- try(ode(y = c(Tc_init, Ts_init, To_init), times = times, func = twolump, parms = indata)[, 2])
  }
  if (class(Te) == "try-error") {
    Te <- rep(NA, length(times))
  }
  if (lump == 1) {
    indata$c_body <- 3073  #specific heat of flesh, J/(kg.C)
  }
  if (lump == 2) {
    indata$c_body_inner <- 3073  #specific heat of flesh, J/(kg.C)
    indata$c_body_outer <- 3073  #specific heat of flesh, J/(kg.C)
  }
  if (lump == 1) {
    Tb_open <- ode(y = Tc_init, times = times, func = onelump_var, parms = indata)[, 2]
  }
  if (lump == 2) {
    Tb_open <- ode(y = c(Tc_init, Ts_init, To_init), times = times, func = twolump, parms = indata)[, 2]
  }
  out <- 0  # initial foraging state
  bask <- 1  # initial basking state
  daybreak <- 0  # initialise daybreak flag
  posture <- "n"  # initial postural state
  arvo <- times[(length(times)/2):length(times)]  # second half of day
  zeniths <- as.data.frame(cbind(arvo, Zenf(arvo)))  # afternoon zenith angles
  colnames(zeniths) <- c("time", "zen")
  evening <- subset(zeniths, zen == 90)  # evening times
  if (nrow(evening) == 0) {
    subtime <- times
  } else {
    sunset <- evening[1, 1]  # time of sunset
    times <- times[times < sunset]  # non-sunset times
    subtime <- times  # starting times to work with
  }
  while (length(subtime) > 0) {
    # now go through the non-evening times and check for daybreak start of the simulation, sun is not up yet, keep
    # it in the shade and inactive until sun comes up
    if (daybreak == 0) {
      indata$posture <- "a"
      indata$Tairf <- Tairf_shd  # choose shaded environment
      indata$Tradf <- Tradf_shd
      indata$Qsolf <- Qsolf_shd
      indata$velf <- velf_shd
      Tbs <- morning()  # get Tbs until sun rises and basking threshold is reached
      Tc_last <- Tbs[nrow(Tbs), 2]  # get initial temp for next behavioural phase
      Tc_prev <- Tbs[nrow(Tbs) - 1, 2]
      while (sign(Tc_last) != sign(Tc_prev)) {
        # get past real zeros!
        Tc_init <- 0
        subtime <- subset(times, times > Tbs[nrow(Tbs), 1])  # get times post basking event, for the next behavioural phase
        Tbs2 <- morning()  # get Tbs until sun rises and basking threshold is reached
        Tc_last <- Tbs2[nrow(Tbs2), 2]  # get initial temp for next behavioural phase
        Tc_prev <- Tbs2[nrow(Tbs2) - 1, 2]
        Tbs <- rbind(Tbs, Tbs2)
      }
      Tbs$posture <- 0
      Tbs$active <- 0
      Tbs$state <- 0
      if (exists("day_results")) {
        day_results <- rbind(day_results, Tbs)
      } else {
        day_results <- Tbs
      }
      Tc_init <- Tbs[nrow(Tbs), 2]  # get initial temp for next behavioural phase
      if (lump == 2) {
        Ts_init <- Tbs[nrow(Tbs), 3]
        To_init <- Tbs[nrow(Tbs), 4]
      }
      subtime <- subset(times, times > Tbs[nrow(Tbs), 1])  # get times post basking event, for the next behavioural phase
      daybreak <- 1  # sun has now risen
    }
    while (bask == 1 & length(subtime) > 0) {
      # now in the basking period
      indata$Tairf <- Tairf_sun  # choose full sun environment
      indata$Tradf <- Tradf_sun
      indata$Qsolf <- Qsolf_sun
      indata$velf <- velf_sun
      if (Tc_init < T_F_min) {
        indata$posture <- "n"  # change posture to be normal to the sun - basking
        Tbs <- cooling()  # simulate Tb until it reaches T_F_min - i.e. until it can forage
        Tbs$posture <- 1
        Tbs$active <- 0
        Tbs$state <- 1
        if (exists("day_results")) {
          day_results <- rbind(day_results, Tbs)
        } else {
          day_results <- Tbs
        }
        Tc_init <- Tbs[nrow(Tbs), 2]
        if (lump == 2) {
          Ts_init <- Tbs[nrow(Tbs), 3]
          To_init <- Tbs[nrow(Tbs), 4]
        }
        subtime <- subset(times, times > Tbs[nrow(Tbs), 1])  # exclude basking time from next simulation
        if (length(subtime) == 0)
        {
          break
        }  # stop if got through the rest of the day simply basking
      }
      indata$posture <- "a"  # has now got to foraging temp, change to foraging posture
      if (length(subtime) == 0) {
        break
      }
      Tbs <- warming()
      Tbs$posture <- 0
      Tbs$active <- 1
      Tbs$state <- 2
      if (exists("day_results")) {
        day_results <- rbind(day_results, Tbs)
      } else {
        day_results <- Tbs
      }
      Tc_init <- Tbs[nrow(Tbs), 2]
      if (lump == 2) {
        Ts_init <- Tbs[nrow(Tbs), 3]
        To_init <- Tbs[nrow(Tbs), 4]
      }
      subtime <- subset(times, times > Tbs[nrow(Tbs), 1])
      if (length(subtime) == 0) {
        break
      }
      if (Tc_init > T_F_min) {
        bask <- 0
      }
    }
    if (length(subtime) == 0) {
      break
    }
    # if we got to here, the animal has been out foraging in the sun and has reached the maximum voluntary foraging
    # temp, so needs to go into shade
    indata$Tairf <- Tairf_shd
    indata$Tradf <- Tradf_shd
    indata$Qsolf <- Qsolf_shd
    indata$velf <- velf_shd
    Tbs <- cooling()  # simulate cooling in shade
    Tbs$posture <- 0
    Tbs$active <- 0
    Tbs$state <- 3
    if (exists("day_results")) {
      day_results <- rbind(day_results, Tbs)
    } else {
      day_results <- Tbs
    }
    Tc_init <- Tbs[nrow(Tbs), 2]
    if (lump == 2) {
      Ts_init <- Tbs[nrow(Tbs), 3]
      To_init <- Tbs[nrow(Tbs), 4]
    }
    subtime <- subset(times, times > Tbs[nrow(Tbs), 1])
    indata$Tairf <- Tairf_sun
    indata$Tradf <- Tradf_sun
    indata$Qsolf <- Qsolf_sun
    indata$velf <- velf_sun
    if (length(subtime) == 0) {
      break
    }
    Tbs <- warming()  # now go foraging again in the sun
    Tbs$posture <- 0
    Tbs$active <- 1
    Tbs$state <- 2
    if (exists("day_results")) {
      day_results <- rbind(day_results, Tbs)
    } else {
      day_results <- Tbs
    }
    Tc_init <- Tbs[nrow(Tbs), 2]
    if (lump == 2) {
      Ts_init <- Tbs[nrow(Tbs), 3]
      To_init <- Tbs[nrow(Tbs), 4]
    }
    subtime <- subset(times, times > Tbs[nrow(Tbs), 1])
    if (Tc_init < T_F_min) {
      indata$posture <- "n"
      indata$Tairf <- Tairf_sun
      indata$Tradf <- Tradf_sun
      indata$Qsolf <- Qsolf_sun
      indata$velf <- velf_sun
      Tbs <- afternoon()
      Tc_last <- Tbs[nrow(Tbs), 2]  # get initial temp for next behavioural phase
      Tc_prev <- Tbs[nrow(Tbs) - 1, 2]
      while (sign(Tc_last) != sign(Tc_prev)) {
        # get past real zeros!
        Tc_init <- 0
        if (lump == 2) {
          Ts_init <- 0
          To_init <- 0
        }
        subtime <- subset(times, times > Tbs[nrow(Tbs), 1])  # get times post basking event, for the next behavioural phase
        Tbs2 <- afternoon()
        Tc_last <- Tbs2[nrow(Tbs2), 2]  # get initial temp for next behavioural phase
        Tc_prev <- Tbs2[nrow(Tbs2) - 1, 2]
        Tbs <- rbind(Tbs, Tbs2)
      }
      Tbs$posture <- 1
      Tbs$active <- 0
      Tbs$state <- 1
      if (exists("day_results")) {
        day_results <- rbind(day_results, Tbs)
      } else {
        day_results <- Tbs
      }
      Tc_init <- Tbs[nrow(Tbs), 2]
      if (lump == 2) {
        Ts_init <- Tbs[nrow(Tbs), 3]
        To_init <- Tbs[nrow(Tbs), 4]
      }
      subtime <- subset(times, times > Tbs[nrow(Tbs), 1])
      if (length(subtime) == 0) {
        break
      }
      indata$posture <- "a"  # if got this far, time to retreat to the shade in prep for the evening
      indata$Tairf <- Tairf_shd
      indata$Tradf <- Tradf_shd
      indata$Qsolf <- Qsolf_shd
      indata$velf <- velf_shd
      Tbs <- morning()
      Tc_last <- Tbs[nrow(Tbs), 2]  # get initial temp for next behavioural phase
      Tc_prev <- Tbs[nrow(Tbs) - 1, 2]
      while (sign(Tc_last) != sign(Tc_prev)) {
        # get past real zeros!
        Tc_init <- 0
        if (lump == 2) {
          Ts_init <- 0
          To_init <- 0
        }
        subtime <- subset(times, times > Tbs[nrow(Tbs), 1])  # get times post basking event, for the next behavioural phase
        Tbs2 <- morning()
        Tc_last <- Tbs2[nrow(Tbs2), 2]  # get initial temp for next behavioural phase
        Tc_prev <- Tbs2[nrow(Tbs2) - 1, 2]
        Tbs <- rbind(Tbs, Tbs2)
      }
      Tbs$posture <- 0
      Tbs$active <- 0
      Tbs$state <- 0
      if (exists("day_results")) {
        day_results <- rbind(day_results, Tbs)
      } else {
        day_results <- Tbs
      }
      Tc_init <- Tbs[nrow(Tbs), 2]
      if (lump == 2) {
        Ts_init <- Tbs[nrow(Tbs), 3]
        To_init <- Tbs[nrow(Tbs), 4]
      }
      subtime <- subset(times, times > Tbs[nrow(Tbs), 1])
      if (length(subtime) == 0) {
        break
      }
    }
  }
  if (length(subtime) == 0 & nrow(evening) != 0) {
    # now simulate the evening
    subtime <- evening[, 1]
    indata$posture <- "a"
    indata$Tairf <- Tairf_shd
    indata$Tradf <- Tradf_shd
    indata$Qsolf <- Qsolf_shd
    indata$velf <- velf_shd
    if (lump == 1) {
      Tbs <- as.data.frame(ode(y = Tc_init, times = subtime, func = onelump_var, parms = indata, method = "lsoda"))
      colnames(Tbs) <- c("time", "Tb", "Tcfinal", "tau", "dTc")
    }
    if (lump == 2) {
      Tbs <- as.data.frame(ode(y = c(Tc_init, Ts_init, To_init), times = subtime, func = twolump, parms = indata,
                               method = "lsoda"))
      colnames(Tbs) <- c("time", "Tb", "Ts", "To", "Tcf")
    }
    Tbs$posture <- 0
    Tbs$active <- 0
    Tbs$state <- 0
    if (exists("day_results")) {
      day_results <- rbind(day_results, Tbs)
    } else {
      day_results <- Tbs
    }
  }

  day_results <- subset(day_results, day_results$time %in% times_orig)
  day_results$Te <- Te[1:nrow(day_results)]
  day_results$Tb_open <- Tb_open[1:nrow(day_results)]

  # prevent constant back-and-forth between sun and shade during transition points by imposing a limit on
  # detection of Tb change by the animal
  day_results$state[day_results$Tb < T_F_min + 0.15 & day_results$Tb > T_F_min - 0.15 & day_results$Tb_open > T_F_max] <- 3
  day_results$state[day_results$Tb < T_F_min + 0.15 & day_results$Tb > T_F_min - 0.15 & day_results$state == 2] <- 1
  day_results$active[day_results$Tb < T_F_min + 0.15 & day_results$Tb > T_F_min - 0.15] <- 0

  interval <- length(times_orig) # number of time steps

  # now get metabolic rates MRT (ml O2 per h) <- 0.110 M ^ 0.768 x 10 ^ (Tb - 20) x log10(Q10)/10
  # from Craig White email 11/8/2014
  Q10 <- 2.44
  mrate.reptile <- (0.11 * Ww_g ^ 0.768 * 10 ^ ((day_results$Tb - 20) * log10(Q10) / 10)) * 0.0056 * (24 / interval) * 3600 / 1000  # 0.0056 converts to Watts, then convert to kJ
  day_results <- cbind(day_results, mrate.reptile)
  mrate_sum <- sum(day_results$mrate.reptile) # cumulated metabolic rate
  inactive <- subset(day_results, active == 0) # inactive steps
  active <- subset(day_results, active == 1) # active steps in sun
  active_shd <- subset(day_results, state == 3) # active steps in shade
  mrate_sum_inactive <- sum(inactive$mrate.reptile) # sum of metabolic rate during inactivity, sun active
  mrate_sum_active <- sum(active$mrate.reptile) # sum of metabolic rate during activity, sun active
  mrate_sum_active_shd <- sum(active_shd$mrate.reptile) # sum of metabolic rate during activity, shade active
  mrate_sum_inactive_shd <- mrate_sum - mrate_sum_active_shd # sum of metabolic rate during inactivity, shade active
  T_F_max_time_Te <- length(Te[Te > T_F_max + 0.5])/(interval/24) * 60 # time Te above T_F_max
  CT_max_time_Te <- length(Te[Te > CT_max])/(interval/24) * 60 # time Te above CT_max
  T_F_max_time_Tb_open <- length(Tb_open[Tb_open > T_F_max + 0.5])/(interval/24) * 60 # time Tb in open above T_F_max
  CT_max_time_Tb_open <- length(Tb_open[Tb_open > CT_max])/(interval/24) * 60 # time Tb in open above CT_max
  T_F_maxtime <- length(day_results$Tb[day_results$Tb > T_F_max + 0.5])/(interval/24) * 60 # time Tb thermoregulating animal is above T_F_max
  CT_max_time <- length(day_results$Tb[day_results$Tb > CT_max + 0.5])/(interval/24) * 60 # time Tb of thermoregulating animal is above CT_max
  max_Tb_open <- max(Tb_open) # maximum Tb in open
  min_Tb_open <- min(Tb_open) # minimum Tb in open
  max_Te <- max(Te) # maximum Te in open
  min_Te <- min(Te) # minimum Te in open
  max_Tb <- max(day_results$Tb) # maximum Tb of thermoregulating animal
  min_Tb <- min(day_results$Tb) # minimum Tb of thermoregulating animal

  # now summarize to hourly activity times and max foraging bouts
  Hour <- trunc(day_results$time / 3600) # make a vector of time in units of hours (10 second intervals)
  day_results <- cbind(Hour, day_results) # add hours vector to day_results

  # animal active in sun
  active <- aggregate(day_results$active, by = list(day_results$Hour), sum) # aggregate active vector by whole hours
  active <- active$x / (interval / 24) * 60 # convert from sum of 10 second time intervals to minutes
  y <- rle(day_results$active) # get the run lengths
  max_foraging_bout <- max((y$lengths[y$values == 1]))/(interval/24) * 60 # get the maximum foraging bout
  z <- rle(day_results$state) # get run lengths by state (0-3, inactive, basking, warming, cooling)
  morning.bask <- z$lengths[z$values == 1][1]/(interval/24) * 60 # get first basking event in minutes
  active.bouts <- y$lengths[y$values == 1] # get lengths of all foraging bouts
  total.bouts <- length(active.bouts) # sum total number of foraging bouts
  if (total.bouts > 1) { # check if more than one foraging bout
    morning.bout <- y$lengths[2]/(interval/24) * 60 # length of first foraging bout in the morning
    afternoon.bout <- y$lengths[length(y$lengths) - 1]/(interval/24) * 60 # length of last foraging bout in the afternoon
  } else {
    morning.bout <- max_foraging_bout # only one foraging bout
    afternoon.bout <- max_foraging_bout # only one foraging bout
  }
  if (total.bouts > 2) { # if there is at least one midday bout
    midday.bouts <- active.bouts[2:(length(active.bouts) - 1)]/(interval/24) * 60 # lengths of all foraging bouts between morning and afternoon bouts
    midday.bout1 <- midday.bouts[1] # length of first foraging bout after the morning bout
    mean.midday.bout <- mean(midday.bouts) # mean length of all foraging bouts between morning and afternoon bout
  } else {
    midday.bout1 <- max_foraging_bout # only a morning and afternoon bout - making this equal to the longest bout
    mean.midday.bout <- max_foraging_bout # only a morning and afternoon bout - making this equal to the longest bout
  }

  # check if no foraging at all and set stats to zero
  if (max_foraging_bout == "-Inf") {
    max_foraging_bout <- 0
    morning.bout <- 0
    midday.bout1 <- 0
    mean.midday.bout <- 0
    afternoon.bout <- 0
    morning.bask <- 0
  }
  sum_activity <- sum(active) # sum of all active hours (minutes)
  active_sun <- active

  # animal active in shade
  shade_active <- day_results$state
  shade_active[shade_active != 3] <- 0
  shade_active[shade_active == 3] <- 1
  active <- aggregate(shade_active, by = list(day_results$Hour), sum) # aggregate active vector by whole hours
  active <- active$x / (interval / 24) * 60 # convert from sum of 10 second time intervals to minutes
  y <- rle(shade_active) # get the run lengths
  max_foraging_bout_shd <- max((y$lengths[y$values == 1]))/(interval/24) * 60 # get the maximum foraging bout
  active.bouts <- y$lengths[y$values == 1] # get lengths of all foraging bouts
  total.bouts_shd <- length(active.bouts) # sum total number of foraging bouts
  if (total.bouts_shd > 1) { # check if more than one foraging bout
    morning.bout_shd <- y$lengths[2]/(interval/24) * 60 # length of first foraging bout in the morning
    afternoon.bout_shd <- y$lengths[length(y$lengths) - 1]/(interval/24) * 60 # length of last foraging bout in the afternoon
  } else {
    morning.bout_shd <- max_foraging_bout_shd # only one foraging bout
    afternoon.bout_shd <- max_foraging_bout_shd # only one foraging bout
  }
  if (total.bouts_shd > 2) { # if there is at least one midday bout
    midday.bouts <- active.bouts[2:(length(active.bouts) - 1)]/(interval/24) * 60 # lengths of all foraging bouts between morning and afternoon bouts
    midday.bout1_shd <- midday.bouts[1] # length of first foraging bout after the morning bout
    mean.midday.bout_shd <- mean(midday.bouts) # mean length of all foraging bouts between morning and afternoon bout
  } else {
    midday.bout1_shd <- max_foraging_bout_shd # only a morning and afternoon bout - making this equal to the longest bout
    mean.midday.bout_shd <- max_foraging_bout_shd # only a morning and afternoon bout - making this equal to the longest bout
  }

  # check if no foraging at all and set stats to zero
  if (max_foraging_bout_shd == "-Inf") {
    max_foraging_bout_shd <- 0
    morning.bout_shd <- 0
    midday.bout1_shd <- 0
    mean.midday.bout_shd <- 0
    afternoon.bout_shd <- 0
  }
  sum_activity_shd <- sum(active) # sum of all active hours (minutes)
  active_shd <- active

  # table of summary statistics
  sum_stats <- as.data.frame(cbind(Ww_g, T_F_min, T_F_max, max_foraging_bout, max_foraging_bout_shd, sum_activity, sum_activity_shd, total.bouts, total.bouts_shd, morning.bask, morning.bout, morning.bout_shd, midday.bout1, midday.bout1_shd, mean.midday.bout, mean.midday.bout_shd, afternoon.bout, afternoon.bout_shd, mrate_sum, mrate_sum_inactive, mrate_sum_inactive_shd, mrate_sum_active, mrate_sum_active_shd, T_F_max_time_Te, CT_max_time_Te, T_F_max_time_Tb_open, CT_max_time_Tb_open, T_F_maxtime, CT_max_time, max_Tb_open, min_Tb_open, max_Te, min_Te, max_Tb, min_Tb))


  colnames(sum_stats) <- c("Ww_g", "T_F_min", "T_F_max", "max_bout_sun", "max_bout_shd", "sum_activity_sun", "sum_activity_shd", "bouts_sun", "bouts_shd", "morning_bask", "morning_forage_sun", "morning_forage_shd",
                           "midday_bout1_sun", "midday_bout1_shd", "mean_midday_bout_sun", "mean_midday_bout_shd", "afternoon_forage_sun", "afternoon_forage_shd", "mrate_sum", "mrate_sum_inactive_sun", "mrate_sum_inactive_shd", "mrate_sum_active_sun", "mrate_sum_active_shd", "T_F_max_time_Te", "CT_max_time_Te",
                           "T_F_max_time_Tb_open", "CT_max_time_Tb_open", "T_F_maxtime", "CT_max_time", "max_Tb_open", "min_Tb_open", "max_Te", "min_Te",
                           "max_Tb", "min_Tb")
  # summarise maximum bout length per hour in sun
  for (i in 0:23) {
    run <- subset(day_results, Hour == i)
    y <- rle(run$active)
    run <- max((y$lengths[y$values == 1]))/(interval/24) * 60 # get maximum bout lenght for the current hour
    if (run == "-Inf") {
      run <- 0
    }
    if (i == 0) {
      runs <- run
    } else {
      runs <- c(runs, run)
    }
  }
  runs_sun <- runs
  runs_shd <- runs

  # summarise maximum bout length per hour in sun
  for (i in 0:23) {
    run <- subset(day_results, Hour == i)
    y <- rle(run$state)
    run <- max((y$lengths[y$values == 3]))/(interval/24) * 60 # get maximum bout length for the current hour
    if (run == "-Inf") {
      run <- 0
    }
    if (i == 0) {
      runs <- run
    } else {
      runs <- c(runs, run)
    }
  }
  runs_shd <- runs


  # activity window for plotting hourly activity and maximum bout length per day
  act_window <- cbind(seq(0, 23, 1), active_sun, runs_sun, active_shd, runs_shd)
  act_window <- as.data.frame(act_window)
  Tair_shd <- Tairf_shd(day_results$time) # spline shaded air temperature at local height
  day_results$Tair_shd <- Tair_shd # append shaded air temperature at local height
  day_results <- day_results[, c(2, 1, 13, 3, 4, 11, 10, 5:9, 12)] # rearrange day_results
  colnames(day_results) <- c("time", "hour", "T_air_shd", "Tb", "Tb_final", "Tb_open", "Te_open", "time_constant", "dTb_dt", "posture", "active", "state", "mrate")
  colnames(act_window) <- c("time", "forage_sun", "max_bout_sun", "forage_shd", "max_bout_shd")
  day_results$mrate <- day_results$mrate * 1000 # convert mrate from kJ to J

  return(list(day_results = day_results, sum_stats = sum_stats, act_window = act_window))
}
