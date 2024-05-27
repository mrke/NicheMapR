## ----echo = FALSE, warning=FALSE----------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE, tidy.opts=list(width.cutoff=60), tidy=TRUE  
)
library(knitr)

## ----message=FALSE, warning=FALSE, eval=FALSE---------------------------------
#  library(devtools)
#  install_github('ropensci/rnoaa') # this package is no longer on CRAN
#  install_github('ilyamaclean/microclima') # install the microclima package
#  install_github('mrke/NicheMapR') # install NicheMapR
#  climate_data_path <- "put your path here" # place where you want to put the global climate dataset
#  get.global.climate(folder = climate_data_path) # download and unpack the global climate data
#  install.packages('tealeaves')
#  install.packages('plantecophys')
#  install.packages('futile.logger') # needed for micro_usa

## ----message=FALSE, warning=FALSE---------------------------------------------
library(NicheMapR)
library(tealeaves)
library(plantecophys)

## ----message=FALSE, warning=FALSE---------------------------------------------
# simulate a microclimate
loc <- c(141.861086, -29.050343) # longitude, latitude, decimal degrees
Usrhyt <- 0.01 # leaf height, m
micro <- micro_global(loc = loc,
                      Usrhyt = Usrhyt,
                      microclima = 1
                      ) # set microclima = 1 to get partitioned diffuse and direct solar - needs internet connection for the DEM download

## ----message=FALSE, warning=FALSE---------------------------------------------
metout <- as.data.frame(micro$metout) # microclimate aboveground conditions
soil <- as.data.frame(micro$soil) # microclimate soil conditions
PDIFs <- micro$diffuse_frac # use variable diffuse fraction
#PDIFs <- rep(0.15, nrow(metout)) # assume uniform diffuse solar
P_a <- get_pressure(micro$elev) # air pressure, Pa
P_a <- rep(P_a, nrow(metout)) # create hourly vector
dates <- micro$dates # mock dates

## ----message=FALSE, warning=FALSE---------------------------------------------
# model settings
live <- 0 # don't simulate behaviour or respiration
leaf <- 1 # use the stomatal conductance formulation for evaporation

## ----message=FALSE, warning=FALSE---------------------------------------------
# leaf functional traits
Ww_g <- 1 # wet weight, g
shape <- 2 # 0=plate, 1=cylinder, 2=ellipsoid
shape_b <- 0.0025 # ratio of b axis:a axis for ellipsoid 
shape_c <- 0.1176 # ratio of c axis:a axis for ellipsoid 
g_vs_ab <- 0.2 # leaf vapour conductance, abaxial (bottom of leaf), mol/m2/s
g_vs_ad <- 0.0 # leaf vapour conductance, adaxial (top of leaf), mol/m2/s
g_vs_ab_max <- 0.3 # leaf vapour conductance, abaxial (bottom of leaf), mol/m2/s
g_vs_ad_max <- 0.0 # leaf vapour conductance, adaxial (top of leaf), mol/m2/s
epsilon_sub <- 0.95 # emissivity of the substrate (-)
epsilon_sky <- 1 # emissivity of the sky (-), accounted for already in 'TSKY' output from microclimate model so make it 1
epsilon_L <- 0.97 # emissivity of the leaf (-)
alpha_L <- 0.5 # solar absorptivity of the leaf (-)
fatosk <- 0.5 # radiation configuration factor to sky (-)
fatosb <- 0.5 # radiation configuration factor to substrate (-)
conv_enhance <- 1.4 # convective enhancement factor (-)
pct_cond <- 0 # percent of leaf conducting to the ground (%)
postur <- 3 # leaf oriented horizontally, not tracking sun

## ----message=FALSE, warning=FALSE---------------------------------------------
# set up vapour conductance vectors and simulate stomatal closure at night
f_open <- rep(1, nrow(metout))
f_open[metout$ZEN == 90] <- 0 # close stomata when the sun is set
g_vs_abs <- g_vs_ab * f_open 
g_vs_ads <- g_vs_ad * f_open

## ----message=FALSE, warning=FALSE---------------------------------------------
ptm <- proc.time() # Start timing
ecto <- ectotherm(leaf = leaf, 
                  live = live,
                  Ww_g = Ww_g, 
                  shape = shape, 
                  shape_b = shape_b, 
                  shape_c = shape_c, 
                  g_vs_ab = g_vs_abs, 
                  g_vs_ad = g_vs_ads, 
                  g_vs_ab_max = g_vs_ab_max, 
                  g_vs_ad_max = g_vs_ad_max, 
                  epsilon_sub = epsilon_sub, 
                  epsilon_sky = epsilon_sky,
                  epsilon = epsilon_L, 
                  alpha_max = alpha_L, 
                  alpha_min = alpha_L,
                  M_1 = 0,
                  fatosk = fatosk, 
                  fatosb = fatosb, 
                  conv_enhance = conv_enhance, 
                  pct_cond = pct_cond, 
                  postur = postur,
                  preshr = P_a, 
                  PDIF = PDIFs)
print(proc.time() - ptm) # Stop the clock
environ <- as.data.frame(ecto$environ)
T_leaf_NMR <- environ$TC

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
month <- 1 # choose a month to plot
sub <- which(floor(dates) + 1 == month)
T_air_1.2m <- metout$TAREF
T_air_1cm <- metout$TALOC
plot(seq(0, 23), T_air_1cm[sub], type = 'l', col = 'blue', ylab = 'temperature, deg C', xlab = 'hour of day', ylim = c(15, 50))
points(seq(0, 23), T_air_1.2m[sub], type = 'l', col = 'blue', lty = 2)
points(seq(0, 23), T_leaf_NMR[sub], type = 'l', ylab = 'T_leaf, deg C')

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
TAs <- metout$TALOC
TGRDs <- soil$D0cm
TSKYs <- metout$TSKYC
VELs <- metout$VLOC
RHs <- metout$RHLOC
QSOLRs <- metout$SOLR
Zs <- metout$ZEN
alpha_sub <- (1 - micro$REFL)

n <- nrow(metout)
T_leaf_ectoR <- matrix(nrow = n)
ptm <- proc.time() # Start timing
for(i in 1:n){
  ectoR.out <- ectoR_devel(leaf = leaf,
                           Ww_g = Ww_g,
                           M_1 = 0, # make metabolic rate zero
                           shape = shape,
                           shape_b = shape_b,
                           shape_c = shape_c,
                           g_vs_ab = g_vs_abs[i],
                           g_vs_ad = g_vs_ads[i],
                           epsilon = epsilon_L,
                           alpha = alpha_L,
                           fatosk = fatosk,
                           fatosb = fatosb,
                           conv_enhance = conv_enhance,
                           pct_cond = pct_cond,
                           postur = postur,
                           elev = micro$elev,
                           alpha_sub = alpha_sub,
                           epsilon_sub = epsilon_sub,
                           epsilon_sky = epsilon_sky,
                           TA = TAs[i],
                           TGRD = TGRDs[i],
                           TSKY = TSKYs[i],
                           VEL = VELs[i],
                           RH = RHs[i],
                           QSOLR = QSOLRs[i],
                           Z = Zs[i],
                           pres = P_a[i],
                           PDIF = PDIFs[i])
  T_leaf_ectoR[i] <- ectoR.out$TC
}
print(proc.time() - ptm) # Stop the clock

# plot
plot(seq(0, 23), T_air_1cm[sub], type = 'l', col = 'blue', ylab = 'temperature, deg C', xlab = 'hour of day', ylim = c(15, 50))
points(seq(0, 23), T_air_1.2m[sub], type = 'l', col = 'blue', lty = 2)
points(seq(0, 23), T_leaf_NMR[sub], type = 'l', ylab = 'T_leaf, deg C')
points(seq(0, 23), T_leaf_ectoR[sub], type = 'l', col = 2, lty = 2)

## ----message=FALSE, warning=FALSE---------------------------------------------
# get characteristic dimension and areas using ecto_devel function GEOM_ecto
GEOM.out <- GEOM_ecto(AMASS = Ww_g / 1000, GEOMETRY = shape, SHP = c(1, shape_b, shape_c), PTCOND = pct_cond / 100, PMOUTH = 0, SKINW = 0 / 100)
w <- GEOM.out$AL / 0.7 # leaf width, m
A <- GEOM.out$AREA # total leaf surface area, m2
#A_sil <- (GEOM.out$ASILN + GEOM.out$ASILP) / 2 # leaf silhouette area to direct solar radiation, m2
A_sil <- A / 2 # half the leaf is facing up to direct solar

g_vs_base <- 0.01 # base leaf vapour conductance, mol/m2/s (equivalent to 0.1 (µmol H2O) / (m^2 s Pa) from tealeaves)

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# Campbell and Norman Calculation
tm <- proc.time() # Start timing
T_leaf_Campbell <- unlist(lapply(1:length(TAs), function(x){leaf_temperature(
  w = w, # leaf width, m
  A = A, # leaf surface area, m^2
  A_sil = A_sil, # leaf silhouette area, m^2
  heliotropic = 0, # don't track the sun
  alpha_L = alpha_L, # leaf solar absorptivity, -
  alpha_S = alpha_sub, # ground solar absorptivity, -
  epsilon_L = epsilon_L, # leaf emissivity, -
  epsilon_sub = epsilon_sub, # ground emissivity, -
  epsilon_sky = epsilon_sky, # sky emissivity, -
  g_vs_ab = g_vs_abs[x] + g_vs_base / 2, # vapour conductance, abaxial (bottom of leaf), mol/m2/s
  g_vs_ad = g_vs_ads[x] + g_vs_base / 2, # vapour conductance, adaxial (top of leaf), mol/m2/s
  TA = TAs[x], # air temperature at lizard height from microclimate model, deg C
  TGRD = TGRDs[x], # ground temperature from microclimate model, deg C
  TSKY = TSKYs[x], # sky temperature from microclimate model, deg C
  VEL = VELs[x], # wind speed from microclimate model, m/s
  RH = RHs[x], # relative humidity from microclimate model, %
  QSOLR = QSOLRs[x], # total horizontal plane solar radiation from microclimate model, W/m2
  Z = Zs[x], # solar zenith angle from microclimate model, degrees
  PRESS = P_a[x],
  PDIF = PDIFs[x], # proportion solar radiation that is diffuse, -
  conv_enhance = conv_enhance # convective enhancement factor
)})) # run leaf_temperature across environments
print(proc.time() - ptm) # Stop the clock

# plot
plot(seq(0, 23), T_air_1cm[sub], type = 'l', col = 'blue', ylab = 'temperature, deg C', xlab = 'hour of day', ylim = c(15, 50))
points(seq(0, 23), T_air_1.2m[sub], type = 'l', col = 'blue', lty = 2)
points(seq(0, 23), T_leaf_NMR[sub], type = 'l', ylab = 'T_leaf, deg C')
points(seq(0, 23), T_leaf_Campbell[sub], type = 'l', col = 2, lty = 2)

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# plant ecophys calculation
UMOLPERJ <- 4.57 # mu_mol photons / J
WETAIR_out <- WETAIR(db = TAs, rh = RHs, bp = P_a)
VPD_Pa <- WETAIR_out$esat - WETAIR_out$e
tm <- proc.time() # Start timing
T_leaf_ecophys <- unlist(lapply(1:length(TAs), function(x){FindTleaf(
  Tair = TAs[x],
  gs = g_vs_abs[x] + g_vs_ads[x] + g_vs_base,
  PPFD = QSOLRs[x] * 0.5 * UMOLPERJ, # assuming PAR is 0.5 of global radiation
  VPD = VPD_Pa[x] / 1000,
  Patm = P_a[x] / 1000,
  Wind = VELs[x],
  Wleaf = w,
  StomatalRatio = 1,
  LeafAbs = alpha_L)})) # run FindTleaf across environments
print(proc.time() - ptm) # Stop the clock

# plot
plot(seq(0, 23), T_air_1cm[sub], type = 'l', col = 'blue', ylab = 'temperature, deg C', xlab = 'hour of day', ylim = c(15, 50))
points(seq(0, 23), T_air_1.2m[sub], type = 'l', col = 'blue', lty = 2)
points(seq(0, 23), T_leaf_NMR[sub], type = 'l', ylab = 'T_leaf, deg C')
points(seq(0, 23), T_leaf_ecophys[sub], type = 'l', col = 2, lty = 2)

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# Leaving the make_* functions empty will automatically default to defaults
# parameters.
leaf_par   <- make_leafpar()   # leaf parameters
enviro_par <- make_enviropar() # environmental parameters
constants  <- make_constants() # physical constants
T_leaf <- tleaf(leaf_par, enviro_par, constants, quiet = TRUE)
T_leaf_tealeaves <- matrix(data = NA, nrow = length(TAs), ncol = 12)

pctdones <- seq(0, length(TAs), round(0.01*length(TAs),0))
ptm <- proc.time() # Start timing
for(i in 1:length(TAs)){
  enviro_par <- make_enviropar(
    replace = list(
      RH = set_units(RHs[i] / 100, "1"),
      S_sw = set_units(QSOLRs[i], "W/m^2"),
      T_air = set_units(TAs[i] + 273.15, "K"),
      T_sky = set_units(TSKYs[i] + 273.15, "K"),
      wind = set_units(VELs[i], "m/s"),
      P = set_units(P_a[i] / 1000, "kPa"),
      r = set_units(micro$REF, "1")
    )
  )
  leaf_par <- make_leafpar(
    replace = list(
      g_sw = set_units((g_vs_abs[i] + g_vs_ads[i]) * 1e6 / P_a[i], "umol/m^2/s/Pa"),
      g_uw = set_units(g_vs_base * 1e6 / P_a[i], "umol/m^2/s/Pa"),
      abs_l = set_units(epsilon_L, "1"),
      abs_s = set_units(alpha_L, "1"),
      leafsize = set_units(GEOM.out$AL, "m")
    )
  )
  T_leaf <- tleaf(leaf_par, enviro_par, constants, quiet = TRUE)
  T_leaf_tealeaves[i, ] <- unlist(T_leaf)
  pctdone <- round(i/length(TAs)*100,1)
  if(i %in% pctdones){
    #cat(round(i/length(TAs)*100,1), 'percent done \n')
  }
}
print(proc.time() - ptm) # Stop the clock
colnames(T_leaf_tealeaves) <- colnames(T_leaf)
T_leaf_tealeaves <- as.data.frame(T_leaf_tealeaves)$T_leaf - 273.15

# plot
plot(seq(0, 23), T_air_1cm[sub], type = 'l', col = 'blue', ylab = 'temperature, deg C', xlab = 'hour of day', ylim = c(15, 50))
points(seq(0, 23), T_air_1.2m[sub], type = 'l', col = 'blue', lty = 2)
points(seq(0, 23), T_leaf_NMR[sub], type = 'l', ylab = 'T_leaf, deg C')
points(seq(0, 23), T_leaf_tealeaves[sub], type = 'l', col = 2, lty = 2)

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
par(mfrow = c(1, 1))
month <- 1 # choose a month to plot
sub <- which(floor(dates) + 1 == month)
plot(seq(0, 23), T_air_1cm[sub], type = 'l', col = 'blue', ylab = 'temperature, deg C', xlab = 'hour of day', ylim = c(15, 50))
points(seq(0, 23), T_air_1.2m[sub], type = 'l', col = 'blue', lty = 2)
points(seq(0, 23), T_leaf_NMR[sub], type = 'l')
points(seq(0, 23), T_leaf_Campbell[sub], col = 'darkgreen', lty = 2, type = 'l')
points(seq(0, 23), T_leaf_ecophys[sub], col = 'brown', lty = 1, type = 'l')
points(seq(0, 23), T_leaf_tealeaves[sub], type = 'l', col = 'darkgreen')
legend(0, 50, legend = c('T_air 1.2m', 'T_air 1cm', 'T_leaf NMR', 'T_leaf C&N', 'T_leaf ecophys', 'T_leaf TL'),
       lty = c(2, 1, 1, 2, 1, 1), col = c("blue", "blue", "black", "darkgreen", "brown", "darkgreen"), bty = 'n')

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
months <- c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
par(mfrow = c(1, 1))
for(i in seq(1, 12, 2)){
month <- i # choose a month to plot
sub <- which(floor(dates) + 1 == month)
plot(seq(0, 23), T_air_1cm[sub], type = 'l', col = 'blue', ylab = 'temperature, deg C', xlab = 'hour of day', ylim = c(0, 50), main = months[i])
points(seq(0, 23), T_air_1.2m[sub], type = 'l', col = 'blue', lty = 2)
points(seq(0, 23), T_leaf_NMR[sub], type = 'l')
points(seq(0, 23), T_leaf_Campbell[sub], col = 'darkgreen', lty = 2, type = 'l')
points(seq(0, 23), T_leaf_ecophys[sub], col = 'brown', lty = 1, type = 'l')
points(seq(0, 23), T_leaf_tealeaves[sub], type = 'l', col = 'darkgreen')
legend(0, 50, legend = c('T_air 1.2m', 'T_air 1cm', 'T_leaf NMR', 'T_leaf C&N', 'T_leaf ecophys', 'T_leaf TL'),
       lty = c(2, 1, 1, 2, 1, 1), col = c("blue", "blue", "black", "darkgreen", "brown", "darkgreen"), bty = 'n')
}

## ----message=FALSE, warning=FALSE---------------------------------------------
# now include simulated stomatal conductance from plantecophys
photo_out <- Photosyn(Tleaf = T_leaf_NMR,
                      g1 = 5.2, # Parameter of Ball-Berry type stomatal conductance models
                      Ca = 400, # Atmospheric CO2 concentration (ppm)
                      VPD = VPD_Pa/1000, # Vapour pressure deficit (kPa)
                      vpdmin = 0.5, # lower threshold on VPD
                      PPFD = QSOLRs * UMOLPERJ * 1, # Photosynthetic photon flux density ('PAR') (mu mol m-2 s-1)
                      Patm = P_a/1000 # Atmospheric pressure (kPa)
                      )
g_vs_abs_plantecophys <- photo_out$GS

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
month <- 1 # choose a month to plot
sub <- which(floor(dates) + 1 == month)
par(mfrow = c(1, 1))
plot(seq(0, 23), g_vs_abs_plantecophys[sub], type = 'l', ylim = c(0, g_vs_ab_max), ylab = 'stomatal conductance, mol/m2/s', xlab = 'hour of day')
abline(h = g_vs_ab_max, lty = 2, col = 3)

## ----message=FALSE, warning=FALSE---------------------------------------------
ptm <- proc.time() # Start timing
ecto <- ectotherm(leaf = leaf, 
                  live = live,
                  Ww_g = Ww_g, 
                  shape = shape, 
                  shape_b = shape_b, 
                  shape_c = shape_c, 
                  g_vs_ab = g_vs_abs_plantecophys, 
                  g_vs_ad = g_vs_ads, 
                  g_vs_ab_max = g_vs_ab_max, 
                  g_vs_ad_max = g_vs_ad_max, 
                  epsilon_sub = epsilon_sub, 
                  epsilon_sky = epsilon_sky,
                  epsilon = epsilon_L, 
                  alpha_max = alpha_L, 
                  alpha_min = alpha_L,
                  M_1 = 0,
                  fatosk = fatosk, 
                  fatosb = fatosb, 
                  conv_enhance = conv_enhance, 
                  pct_cond = pct_cond, 
                  postur = postur,
                  preshr = P_a, 
                  PDIF = PDIFs)
print(proc.time() - ptm) # Stop the clock
environ <- as.data.frame(ecto$environ)
T_leaf_NMR_plantecophys <- environ$TC

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
plot(seq(0, 23), T_leaf_NMR_plantecophys[sub], type = 'l', lty = 2, ylab = 'temperature, deg C', xlab = 'hour of day', ylim = c(15, 50))
points(seq(0, 23), T_leaf_NMR[sub], type = 'l')
abline(h = 45, lty = 2, col = 2)
abline(h = 50, lty = 1, col = 2)

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# now do a loop of n iterations to see when convergence occurs

n <- 6 # iterations to do to find steady state stomatal conductance
par(mfrow = c(2, 3))
g_vs_abs_iterate <- g_vs_abs
T_leaf_NMR_iterate <- NULL
for(i in 1:n){
  ecto <- ectotherm(leaf = leaf, 
                  live = live,
                  Ww_g = Ww_g, 
                  shape = shape, 
                  shape_b = shape_b, 
                  shape_c = shape_c, 
                  g_vs_ab = g_vs_abs_iterate, 
                  g_vs_ad = g_vs_ads, 
                  g_vs_ab_max = g_vs_ab_max, 
                  g_vs_ad_max = g_vs_ad_max, 
                  epsilon_sub = epsilon_sub, 
                  epsilon_sky = epsilon_sky,
                  epsilon = epsilon_L, 
                  alpha_max = alpha_L, 
                  alpha_min = alpha_L,
                  M_1 = 0,
                  fatosk = fatosk, 
                  fatosb = fatosb, 
                  conv_enhance = conv_enhance, 
                  pct_cond = pct_cond, 
                  postur = postur,
                  preshr = P_a, 
                  PDIF = PDIFs)
  environ <- as.data.frame(ecto$environ)
  if(i > 1){
   T_leaf_prev <- T_leaf_NMR_iterate
  }
  T_leaf_NMR_iterate <- environ$TC
  photo_out <- Photosyn(Tleaf = T_leaf_NMR_iterate,
                        g1 = 5.2,
                        Ca = 400,
                        VPD = VPD_Pa/1000,
                        vpdmin = 0.5,
                        PPFD = QSOLRs * UMOLPERJ * 1,
                        Patm = P_a/1000)
  g_vs_abs_iterate <- photo_out$GS
  plot(seq(0, 23), T_leaf_NMR_iterate[sub], type = 'l', col = 1, ylim = c(15, 55), main = i, ylab = 'leaf temperature, deg C', xlab = 'hour of day')
  points(seq(0, 23), T_air_1cm[sub], type = 'l', col = 'blue')
  if(i > 1){
     points(seq(0, 23), T_leaf_prev[sub], type = 'l', lty = 2, col = 'grey')
  }
  abline(h = 45, lty = 2, col = 2)
  abline(h = 50, lty = 1, col = 2)
}


## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# plot photosynthesis outputs
par(mfrow = c(1, 1))
A_j <- photo_out$Aj
A_c <- photo_out$Ac
R_d <- photo_out$Rd
A_n <- photo_out$ALEAF # A_n = pmin(A_c, A_j) - R_d
plot(seq(0, 23), A_j[sub], type = 'l', ylab = 'mu mol m-2 s-1', main = 'computed leaf temperature', xlab = 'hour of day', ylim = c(0, 25), col = 'grey')
points(seq(0, 23), A_c[sub], type = 'l', lty = 2)
points(seq(0, 23), R_d[sub], type = 'l', ylab = 'mu mol m-2 s-1', col = 'brown', lty = 2)
points(seq(0, 23), A_n[sub], type = 'l', col = 'red')
legend(2, 25, c('Rubisco-limited', 'RUBP-limited', 'Net', 'respiration'), col = c("grey", "black", "red", "brown"), lty = c(1, 2, 1, 2), bty = 'n')

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# redo with T_leaf = T_air
photo_out_T_air <- Photosyn(Tleaf = T_air_1cm,
                      g1 = 5.2,
                      Ca = 400,
                      VPD = VPD_Pa/1000,
                      vpdmin = 0.5,
                      PPFD = QSOLRs * UMOLPERJ * 1,
                      Patm = P_a/1000)

# plot photosynthesis outputs
A_j2 <- photo_out_T_air$Aj
A_c2 <- photo_out_T_air$Ac
R_d2 <- photo_out_T_air$Rd
A_n2 <- photo_out_T_air$ALEAF
plot(seq(0, 23), A_j2[sub], type = 'l', ylab = 'mu mol m-2 s-1', main = 'leaf temperature = air temperature', xlab = 'hour of day', ylim = c(0, 25), col = 'grey')
points(seq(0, 23), A_c2[sub], type = 'l', lty = 2)
points(seq(0, 23), R_d2[sub], type = 'l', ylab = 'mu mol m-2 s-1', col = 'brown', lty = 2)
points(seq(0, 23), A_n2[sub], type = 'l', col = 'red')
legend(2, 25, c('Rubisco-limited', 'RUBP-limited', 'Net', 'respiration'), col = c("grey", "black", "red", "brown"), lty = c(1, 2, 1, 2), bty = 'n')

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# compare net photosynthesis using calculated leaf temperature or air temperature
plot(seq(0, 23), A_n[sub], type = 'l', ylab = 'mu mol m-2 s-1', main = 'leaf vs air temperature', xlab = 'hour of day', ylim = c(-2, 15))
points(seq(0, 23), A_n2[sub], type = 'l', col = 'red', lty = 2)
legend(2, 12, c('calculated T_leaf', 'T_leaf = T_air'), col = c("black", "red"), lty = c(1, 2), bty = 'n')

## ----message=FALSE, warning=FALSE---------------------------------------------
loc <- c(-116.53, 33.83)
Usrhyt <- 0.30 # leaf height
PC <- -1500 # critical leaf water potential for stomatal closure, J/kg
SP <- 10 # stability parameter for stomatal closure, -

micro <- micro_usa(loc = loc,
                   Usrhyt = Usrhyt,
                   PC = PC,
                   SP = SP,
                   dstart = "01/01/2016",
                   dfinish = "31/12/2017"
                   )
# define time zone
tz <- paste0("Etc/GMT-", floor(micro$longlat[1] / 15 * -1))
dates <- micro$dates
attr(dates, "tzone") <- tz
dates2 <- micro$dates2

## ----message=FALSE, warning=FALSE---------------------------------------------
metout <- as.data.frame(micro$metout)
soilpot <- as.data.frame(micro$soilpot)
plant <- as.data.frame(micro$plant)

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# calculate stomata opennes as a function of leaf water potential
f_open <- 1 / (1 + (plant$LPT / PC) ^ SP)
plot(dates, f_open, type = 'l', ylab = 'fraction open', xlab = '', main = 'stomatal closure due to water stress')

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
g_vs_ab <- 0.4 # leaf vapour conductance, abaxial (bottom of leaf), mol/m2/s
g_vs_ab_max <- 0.4 # leaf vapour conductance, abaxial (bottom of leaf), mol/m2/s

# use Photosyn to get stomatal dynamics
P_a <- get_pressure(micro$elev) # air pressure, Pa
P_a <- rep(P_a, nrow(metout)) # create hourly vector
TAs <- metout$TALOC
RHs <- metout$RHLOC
QSOLRs <- metout$SOLR
WETAIR_out <- WETAIR(db = TAs, rh = RHs, bp = P_a)
VPD_Pa <- WETAIR_out$esat - WETAIR_out$e
UMOLPERJ <- 4.57 # mu_mol photons / J

photo_out <- Photosyn(Tleaf = TAs,
                      g1 = 5.2,
                      Ca = 400,
                      VPD = VPD_Pa/1000,
                      vpdmin = 0.5,
                      PPFD = QSOLRs * UMOLPERJ * 1,
                      Patm = P_a/1000)
if(g_vs_ab > 0){
  g_vs_abs <- photo_out$GS
}else{
  g_vs_abs <- photo_out$GS * 0
}
if(g_vs_ad > 0){
  g_vs_ads <- photo_out$GS
}else{
  g_vs_ads <- photo_out$GS * 0
}
g_vs_abs_photosyn <- g_vs_abs # save for comparison plot
sub <- which(as.numeric(format(dates, '%Y')) >= 2017)
plot(dates[sub], g_vs_abs[sub], type = 'l', ylab = 'stomatal conductance, mol/m2/s', xlab = 'time')

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
ptm <- proc.time() # Start timing
ecto <- ectotherm(leaf = leaf, 
                  live = live,
                  Ww_g = Ww_g, 
                  shape = shape, 
                  shape_b = shape_b, 
                  shape_c = shape_c, 
                  g_vs_ab = g_vs_abs, 
                  g_vs_ad = g_vs_ads, 
                  g_vs_ab_max = g_vs_ab_max, 
                  g_vs_ad_max = g_vs_ad_max, 
                  epsilon_sub = epsilon_sub, 
                  epsilon_sky = epsilon_sky,
                  epsilon = epsilon_L, 
                  alpha_max = alpha_L, 
                  alpha_min = alpha_L,
                  M_1 = 0,
                  fatosk = fatosk, 
                  fatosb = fatosb, 
                  conv_enhance = conv_enhance, 
                  pct_cond = pct_cond, 
                  postur = postur,
                  preshr = P_a, 
                  PDIF = PDIFs)
print(proc.time() - ptm) # Stop the clock
environ <- as.data.frame(ecto$environ)
T_leaf_photosyn <- environ$TC
plot(dates[sub], T_leaf_photosyn[sub], type = 'l', ylim = c(-5, 70), ylab = 'leaf temperature, deg C', xlab = 'time', col = 'blue')

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
# impose water stress effects
g_vs_abs[g_vs_ab * f_open < g_vs_abs] <- g_vs_ab * f_open[g_vs_ab * f_open < g_vs_abs]
g_vs_ads[g_vs_ad * f_open < g_vs_ads] <- g_vs_ad * f_open[g_vs_ad * f_open < g_vs_ads]

plot(dates[sub], g_vs_abs_photosyn[sub], type = 'l', ylab = 'stomatal conductance, mol/m2/s', xlab = 'time')
points(dates[sub], g_vs_abs[sub], type = 'l', ylab = 'stomatal conductance, mol/m2/s', xlab = 'time', col = 'grey')

## ----message=FALSE, warning=FALSE---------------------------------------------
ptm <- proc.time() # Start timing
ecto <- ectotherm(leaf = leaf, 
                  live = live,
                  Ww_g = Ww_g, 
                  shape = shape, 
                  shape_b = shape_b, 
                  shape_c = shape_c, 
                  g_vs_ab = g_vs_abs, 
                  g_vs_ad = g_vs_ads, 
                  g_vs_ab_max = g_vs_ab_max, 
                  g_vs_ad_max = g_vs_ad_max, 
                  epsilon_sub = epsilon_sub, 
                  epsilon_sky = epsilon_sky,
                  epsilon = epsilon_L, 
                  alpha_max = alpha_L, 
                  alpha_min = alpha_L,
                  M_1 = 0,
                  fatosk = fatosk, 
                  fatosb = fatosb, 
                  conv_enhance = conv_enhance, 
                  pct_cond = pct_cond, 
                  postur = postur,
                  preshr = P_a, 
                  PDIF = PDIFs)
print(proc.time() - ptm) # Stop the clock
environ <- as.data.frame(ecto$environ)
T_leaf <- environ$TC
#plot(dates[sub], T_leaf[sub], type = 'l', ylim = c(-5, 70), ylab = 'leaf #temperature, deg C', xlab = 'time', col = 'red')
#points(dates[sub], T_leaf_photosyn[sub], type = 'l', col = 'blue')

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
photo_out <- Photosyn(Tleaf = T_leaf,
                      g1 = 5.2,
                      Ca = 400,
                      VPD = VPD_Pa/1000,
                      vpdmin = 0.5,
                      PPFD = QSOLRs * UMOLPERJ * 1,
                      Patm = P_a/1000)
if(g_vs_ab > 0){
  g_vs_abs <- photo_out$GS
}else{
  g_vs_abs <- photo_out$GS * 0
}
if(g_vs_ad > 0){
  g_vs_ads <- photo_out$GS
}else{
  g_vs_ads <- photo_out$GS * 0
}
g_vs_abs[g_vs_abs * f_open < g_vs_abs] <- g_vs_abs[g_vs_abs * f_open < g_vs_abs] * f_open[g_vs_abs * f_open < g_vs_abs]
g_vs_ads[g_vs_ads * f_open < g_vs_ads] <- g_vs_ads[g_vs_ads * f_open < g_vs_ads] * f_open[g_vs_ads * f_open < g_vs_ads]

ptm <- proc.time() # Start timing
ecto <- ectotherm(leaf = leaf, 
                  live = live,
                  Ww_g = Ww_g, 
                  shape = shape, 
                  shape_b = shape_b, 
                  shape_c = shape_c, 
                  g_vs_ab = g_vs_abs, 
                  g_vs_ad = g_vs_ads, 
                  g_vs_ab_max = g_vs_ab_max, 
                  g_vs_ad_max = g_vs_ad_max, 
                  epsilon_sub = epsilon_sub, 
                  epsilon_sky = epsilon_sky,
                  epsilon = epsilon_L, 
                  alpha_max = alpha_L, 
                  alpha_min = alpha_L,
                  M_1 = 0,
                  fatosk = fatosk, 
                  fatosb = fatosb, 
                  conv_enhance = conv_enhance, 
                  pct_cond = pct_cond, 
                  postur = postur,
                  preshr = P_a, 
                  PDIF = PDIFs)
print(proc.time() - ptm) # Stop the clock
environ <- as.data.frame(ecto$environ)
T_leaf <- environ$TC
plot(dates[sub], T_leaf[sub], type = 'l', ylim = c(-5, 70), ylab = 'leaf temperature, deg C', xlab = 'time', col = 'red')
points(dates[sub], T_leaf_photosyn[sub], type = 'l', col = 'blue')

## ----fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE----
plot(dates[sub], T_leaf[sub] - T_leaf_photosyn[sub], type = 'l', ylab = 'difference, deg C', xlab = 'time')

