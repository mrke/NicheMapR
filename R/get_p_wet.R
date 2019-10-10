#' Function to find effective proportion of surface area that is wet from experimental water loss data
#'
#' A set of equations for computing the proportion of the total surface area of an organism that
#' is acting as a free water surface for evaporation, which is used in the ectotherm model of NicheMapR
#' NicheMapR to compute cutaneous water loss rates. This function requires the WETAIR and DRYAIR
#' functions of NicheMapR.
#' Elia Pirtle and Michael Kearney developed this R function and example in November 2017.
#' @param cut.known = 0, Is cutaneous water loss supplied? 0 = no (invokes estimation by partitioning respiratory), 1 = yes
#' @param resp.known = 0, Is respiratory water loss supplied? 0 = no (invokes estimation by partitioning respiratory), 1 = yes
#' @param eye.known = 0, Is ocular water loss supplied? 0 = no (invokes estimation by ocular surface area paritioning), 1 = yes
#' @param vent.known = 0, Is ventilation rate supplied? 0 = no (rely on estimate from V_O2), 1 = yes
#' @param press.known = 0, Is barometric pressure known (1, uses altitude) or has it been provided as the 'bp' variable (0)
#' @param A_tot.known = 0, Is total surface area known? 0 = no (uses allometric equation for Dipsosaurus dorsalis, the desert iguana), 1 = yes
#' @param E_tot = 0.344 / 60 / 1000, Observed total evaporative water loss rate, g/s
#' @param E_cut = 0.108 / 60 / 1000, Observed cutaneous evaporative water loss rate, g/s
#' @param E_resp = 0.0371 / 60 / 1000, Observed respiratory evaporative water loss rate, g/s
#' @param E_eye = 0.199 / 60 / 1000, Observed ocular evaporative water loss rate, g/s
#' @param mass.g = 26.7, Wet mass of animal, g
#' @param A_tot = 0.00894, Total surface area, m2
#' @param fcond = 0.25, Fraction of total area area touching surfaces, dec\% (therefore no evaporation or convection)
#' @param feye = 0.002, Fraction of total area area comprising the eye
#' @param eye_frac = 0.58, Fraction of time eye is open, dec\%
#' @param T_b = 25, Core body temperature, °C
#' @param T_s = 24.99, Skin (surface) temperature, °C
#' @param T_a = 24.98, Air temperature, °C
#' @param vel = 0.0021, wind speed, m / s
#' @param RHex = 100, relative humidity of exhaled air, \%
#' @param RHin = 30, relative humidity of atmosphere, \%
#' @param E_O2 = 15, Oxygen extraction efficiency, \%
#' @param RQ = 0.766, Respiratory quotient, dec\%
#' @param M_1 = 0.013, Metabolic rate equation parameter 1 V_O2 = M_1 * M ^ M_2 * 10 ^ (M_3 * T_b) based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231
#' @param M_2 = 0.800, Metabolic rate equation parameter 2
#' @param M_3 = 0.038, Metabolic rate equation parameter 3
#' @param V_O2_STP = M_1 * mass.g ^ M_2 * 10 ^ (M_3 * T_b) / 3600 / 1000, Oxygen consumption rate corrected to STP, L O2 / s
#' @param T_ref_O2 = 25, Reference temperature for oxygen consumption rate correction (°C)
#' @param V_air_STP = V_O2_STP / (0.2094 * E_O2 / 100), Ventilation rate at STP, L/s
#' @param T_ref_vent = 25, Reference temperature for ventilation rate correction (°C)
#' @param rho_flesh = 1000, Density of flesh (kg/m3)
#' @param T_A = 8817, Arhhenius temperature (K)
#' @param T_AL = 50000, Arrhenius temperature for decrease below lower boundary of tolerance range T_L (K)
#' @param T_AH = 90000, Arrhenius temperature for decrease above upper boundary of tolerance range T_H (K)
#' @param T_L = 279, Lower boundary (K) of temperature tolerance range for Arrhenius thermal response
#' @param T_H = 306, Upper boundary (K) of temperature tolerance range for Arrhenius thermal response
#' @param alt = 0.03, Altitude (km)
#' @param bp	= 101325, Barometric pressure (pascal)
#' @return E_tot total evaporative water loss, g/s
#' @return E_cut cutaneous evaporative water loss, g/s
#' @return E_resp respiratory evaporative water loss, g/s
#' @return E_eye ocular evaporative water loss, g/s
#' @return R_a boundary layer resitance, s/m
#' @return R_s skin resistance, s/m
#' @return A_tot total surface area, m2
#' @return A_eye ocular surface area, m2
#' @return A_eff effective area for evaporation, m2
#' @return A_wet area acting as a free water surface, m2
#' @return p_wet proportion of surface area acting as a free water surface, dec\%
#' @usage get_p_wet(cut.known = 1, E_tot = 0.344 / 60 / 1000, E_cut = 0.108 / 60 / 1000, T_b = 25, RHin = 30, vel = 0.0021)
#' @export
get_p_wet <- function(cut.known = 0,
                      resp.known = 0,
                      eye.known = 0,
                      vent.known = 0,
                      press.known = 0,
                      A_tot.known = 0,
                      E_tot = 0.344 / 60 / 1000,
                      E_cut = 0.108 / 60 / 1000,
                      E_resp = 0.0371 / 60 / 1000,
                      E_eye = 0.199 / 60 / 1000,
                      mass.g = 26.7,
                      A_tot = 0.00894,
                      fcond = 0.25,
                      feye = 0.002,
                      eye_frac = 0.58,
                      T_b = 25,
                      T_s = T_b - 0.001,
                      T_a = T_b - 0.002,
                      vel = 0.0021,
                      RHex = 100,
                      RHin = 30,
                      E_O2 = 15,
                      RQ = 0.766,
                      M_1 = 0.013,
                      M_2 = 0.8,
                      M_3 = 0.0038,
                      V_O2_STP = M_1 * mass.g ^ M_2 * 10 ^ (M_3 * T_b) / 3600 / 1000,
                      T_ref_O2 = 30,
                      V_air_STP = V_O2_STP / (0.2094 * E_O2 / 100),
                      rho_flesh = 1000,
                      T_ref_vent = 30,
                      T_A = 8817,
                      T_AL = 50000,
                      T_AH = 90000,
                      T_L = 279,
                      T_H = 306,
                      alt = 0,
                      bp = 101325) {
  require(NicheMapR)
  # constants
  fO2_ref <- 0.2095	# ref O2 proportion, dec%
  fO2 <- 0.2095	# percent O2, dec%
  fCO2 <- 0.0003 # percent CO2, dec%
  fN2 <- 0.7902	# percent N2, dec%
  RGC <- 8309.28 # Universal gas constant, (Pa L)/(mol K)
  PO2 <- fO2 * bp # Partial pressure O2, Pa
  PO2_ref <- fO2_ref * 101325 # reference partial pressure O2, Pa
  G <- 9.80665	# acceleration from gravity, m/s
  C_P <- 1.01E+03	# specific heat dry air, J/(kg K)
  BETA <-	1 / (T_a + 273.15) # coeff thermal expansion at constant density, 1/K

  # air properties
  dryair <- DRYAIR(db = T_a, alt = alt)
  if(press.known == 0){bp <- dryair$patmos} # barometric pressure, Pa
  rho_dry <- dryair$densty # density, kg/m3
  mu <- dryair$visdyn # dynamic viscosity, kg/(m s)
  k_air <- dryair$thcond # thermal conductivity, W/(m K)
  D_H2O <- dryair$difvpr # diffusivity of water vapour in air, m2/s

  # vapour densities
  rho_H2O_air <- WETAIR(db = T_a, rh = RHin)$vd # vapour density of air, kg/m3
  rho_H2O_skin <- WETAIR(db = T_s, rh = RHex)$vd # vapour density at skin, kg/m3
  rho_H2O_lung <- WETAIR(db = T_b, rh = RHex)$vd #vapour density in lung, kg/m3

  # Arrhenius temperature correction factor, if ventilation rate measured at different T_b to water loss rate
  T_corr_vent <- exp(T_A * (1 / (273.15 + T_ref_vent) - 1 / (273.15 + T_b))) / (1 + exp(T_AL * (1 / (273.15 + T_b) - 1 / T_L)) + exp(T_AH * (1 / T_H - 1 / (273.15 + T_b))))
  T_corr_O2 <- exp(T_A * (1 / (273.15 + T_ref_O2) - 1 / (273.15 + T_b))) / (1 + exp(T_AL * (1 / (273.15 + T_b) - 1 / T_L)) + exp(T_AH * (1 / T_H - 1 / (273.15 + T_b))))

  # compute molar flux of air through the lungs
  if(vent.known == 1) {
    # compute air mass flux from observed ventilation rate
    V_Tair_STP <- V_air_STP * T_corr_vent # L/s air, corrected to temperature of water loss rate measurement
    J_air_in <-  bp * V_Tair_STP /(RGC * (T_b + 273.15)) # mol/s air
  }else{
    # computing air mass flux from measured oxygen consumption rate
    V_TO2_STP <- V_O2_STP * T_corr_O2 # L/s O2, corrected to temperature of water loss rate measurement
    V_O2 <- (V_TO2_STP * PO2_ref / 273.15) * ((T_b + 273.15) / PO2) # L/s O2 consumption at expermiental temperature and pressure
    J_O2_con <- bp * V_O2 / (RGC * (T_b + 273.15)) # mol/s O2 consumed at expermiental temperature and pressur
    J_O2_in <- J_O2_con / (E_O2 / 100) # mol/s O2 actually flowing through lungs
    Airato <- (fN2 + fO2 + fCO2) / fO2 # air to oxygen ratio
    J_air_in <- J_O2_in * Airato * (fO2_ref / fO2) * (PO2_ref / PO2) # total moles of air at 1 (mol/s)
  }
  #V_air = (J_air_in * RGC * (T_b + 273.15)) / PO2 # L/s, air flow through lungs (note dividing by PO2)

  if(resp.known == 0){
    # compute molar flux of water into the lungs
    e_air <- WETAIR(db = T_a, rh = RHin)$e
    J_H2O_in <- J_air_in * (e_air * (RHin / 100)) / (bp - e_air * (RHin / 100))

    # compute molar flux of water out of the lungs
    if(vent.known == 1) {
      # assume volume/mass that went in came out
      J_air_out <- J_air_in # air flux at exit, mol/s
    }else{
      # adjusting for CO2/O2 ratios between incoming and outgoing breath
      J_O2_out <- J_O2_in - J_O2_con # O2 flux at exit, mol/s
      J_CO2<- J_O2_con * RQ # Co2 flux at exit, mol/s
      J_air_out <- (J_O2_out + J_CO2) * ((fN2 + fO2) / fO2) * (fO2_ref / fO2) * (PO2_ref / PO2) #Moles air at 2 (mol/s)
    }
    e_lung <- WETAIR(db = T_b, rh = 100)$e
    J_H2O_out <- J_air_out*(e_lung/(bp-e_lung)) # water exiting, mol/s

    # total respiratory water loss
    J_H20_resp <- J_H2O_out - J_H2O_in # respiratory evaporative water lost, mol/s
    E_resp <- J_H20_resp * 18 # respiratory evaporative water, g/s (moles lost * gram molecular weight water)
  }

  if(cut.known == 0){
    E_cut <- E_tot - E_resp # cutaneous evaporative water loss, g/s
  }

  # compute convection and mass transfer coefficients using the Colburn analogy
  VOL <- mass.g / 1000 / rho_flesh # animal volume
  D <- VOL ^ (1/3) # m, characteristic dimension for convection
  if(A_tot.known == 0){
    # equation for Dipsosaurus dorsalis, from Porter et al. 1973, Oecologia
    #A_tot <- (10.4713 * mass ^ 0.688) / 10000 # total surface area, m2
    A_tot <- 10 * mass.g ^ (2 / 3) / 100 ^ 2
  }
  A_eye <- A_tot * feye
  A_conv <- A_tot * (1 - fcond)	# area for convection, m2

  # free convection
  PR_free <- C_P * mu / k_air # Prandlt number for free convection, -
  SC_free <- mu / (rho_dry * D_H2O) # Schmidt number for free convection, -
  delta_T <- T_s - T_a # temperature difference between skin and air, °C
  delta_T[abs(delta_T) < 0.01] <- 0.01 * sign(delta_T[abs(delta_T) < 0.01]) # impose lower limit of 0.01 °C
  GR <- abs(((rho_dry ^ 2) * BETA * G * (D ^ 3) * delta_T) / (mu ^ 2)) # Grashof number, -
  RE <- rho_dry * vel * D / mu # Reynold's number, -
  RA <- (GR ^ (1/4)) * (PR_free ^ (1/3)) # Rayleigh number, - (formula for sphere or lizard)
  NU_free <- 2 + 0.60 * RA # Nusslet number for free convection, -
  HC_free <- (NU_free * k_air) / D	# Heat transfer coeff for free convection, W/(m2 K)
  HC_free[HC_free < 5] <- 5 # don't let it go below 5, as this is unrealistic
  NU_free <- HC_free * D / k_air # recalcualting NU_free in case HC_free had gone below lower limit of 5
  SH_free <- NU_free * (SC_free / PR_free) ^ (1/3) # Sherwood number, -
  HD_free <- SH_free * D_H2O / D # mass transfer coeff for free convection, m/s
  #Q_conv_free <- HC_free * convar * (T_s - T_a) # Free convective heat loss at skin, W

  # forced convection
  PR_forc <-	0.72 # Prandlt number for forced convection, -
  SC_forc <-	0.6 # Schmidt number for forced convection, -)
  NU_forc <- 0.35 * RE ^ 0.6	# Nusslet number for forced convection, -
  HC_forc <- NU_forc * k_air / D # Heat transfer coeff for forced convection, W/(m2 K)
  SH_forc <- NU_forc* (SC_forc / PR_forc) ^ (1 / 3) # Sherwood number, -
  HD_forc <-	SH_forc * D_H2O / D	# mass transfer coeff, m/s
  #Q_conv_forc <-	HC_forc * convar * (T_s - T_a) # Forced convective heat loss at skin, W

  # combined heat and mass transfer coefficients
  HC <- HC_free + HC_forc # Final heat transfer coeff, W/(m2 K)
  HD <- HD_free + HD_forc # Final mass transfer coefficient, m/s
  #Q_conv <- Q_conv_free + Q_conv_forc	# Total convective heat loss at skin, W

  # ocular water loss
  if(eye.known  == 0){
    E_eye <- (A_eye * eye_frac) * HD * (rho_H2O_skin - rho_H2O_air) / 1000 # g/s
    E_cut <- E_cut - E_eye # subtract any ocular water loss from the estimate of cutaneous
  }

  # compute actual skin resistance (unlike p_wet, this does not include eyes)
  R_a <- rho_dry * C_P / HC # boundary layer resistance, s/m
  R_s <- (rho_H2O_skin - rho_H2O_air) / ((E_cut / 1000)/(A_tot - A_eye)) - R_a # skin resistance, s/m

  # compute proportion of skin area that is acting like a free-water surface, p_wet, and associated areas
  A_eff <- (E_cut / 1000) / (HD * (rho_H2O_skin - rho_H2O_air)) # effective area for evaporation, m2
  p_wet = A_eff / (A_tot - A_tot * fcond - A_eye * eye_frac + A_eye * (1 - eye_frac))
  A_wet = p_wet * A_tot

  result <- cbind(E_tot, E_cut, E_resp, E_eye, R_a, R_s, A_tot, A_eye, A_eff, A_wet, p_wet)
  colnames(result) <- c('E_tot', 'E_cut', 'E_resp', 'E_eye', 'R_a', 'R_s', 'A_tot', 'A_eye', 'A_eff', 'A_wet', 'p_wet')
  return(result)
}
