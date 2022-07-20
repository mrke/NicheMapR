#' Dynamic Energy Budget model
#'
#' Implementation of the Standard Dynamic Energy Budget model of Kooijman
#' Note that this uses the deSolve package 'ode' function. The older version
#' that uses Euler integration is now called DEB_euler (and is faster and
#' may be preferable in some cases, though accuracy of the latter will depend
#' on the step size chosen)
#' Michael Kearney Dec 2015, updated to include ODE solver Feb 2019
#' @param step = 1/24, step size (days)
#' @param z = 7.997, Zoom factor (cm)
#' @param del_M =  0.242, Shape coefficient (-)
#' @param p_Xm = 13290*step, Surface area-specific maximum feeding rate J/cm2/h
#' @param kap_X = 0.85, Digestive efficiency (decimal \%)
#' @param v = 0.065*step, Energy conductance (cm/h)
#' @param kap = 0.886, fraction of mobilised reserve allocated to soma
#' @param p_M = 32*step, Volume-specific somatic maintenance (J/cm3/h)
#' @param p_T = 0, (Structural-)Surface-area-specific heating cost (J/cm2/h)
#' @param E_G = 7767, Cost of structure (J/cm3)
#' @param kap_R = 0.95, Fraction of reproduction energy fixed in eggs
#' @param k_J = 0.002*step, Maturity maintenance rate coefficient (1/h)
#' @param E_Hb = 7.359e+04, Maturity at birth (J)
#' @param E_Hj = E_Hb, Maturity at metamorphosis (J)
#' @param E_Hp = 1.865e+05, Maturity at puberty
#' @param h_a = 2.16e-11*(step^2), Weibull ageing acceleration (1/h2)
#' @param s_G = 0.01, Gompertz stress coefficient (-)
#' @param E_0 = 1.04e+06, Energy content of the egg (derived from core parameters) (J)
#' @param T_REF = 20+273.15, Reference temperature for rate correction (deg C)
#' @param T_A = 8085 Arrhenius temperature
#' @param T_AL = 18721, Arrhenius temperature for decrease below lower boundary of tolerance range \code{T_L}
#' @param T_AH = 90000, Arrhenius temperature for decrease above upper boundary of tolerance range \code{T_H}
#' @param T_L = 288, Lower boundary (K) of temperature tolerance range for Arrhenius thermal response
#' @param T_H = 315, Upper boundary (K) of temperature tolerance range for Arrhenius thermal response
#' @param T_A2 = 8085 Arrhenius temperature for maturity maintenance (causes 'Temperature Size Rule' effect)
#' @param T_AL2 = 18721, Arrhenius temperature for decrease below lower boundary of tolerance range \code{T_L}  for maturity maintenance (causes 'Temperature Size Rule' effect)
#' @param T_AH2 = 90000, Arrhenius temperature for decrease above upper boundary of tolerance range \code{T_H} for maturity maintenance (causes 'Temperature Size Rule' effect)
#' @param T_L2 = 288, Lower boundary (K) of temperature tolerance range for Arrhenius thermal response for maturity maintenance (causes 'Temperature Size Rule' effect)
#' @param T_H2 = 315, Upper boundary (K) of temperature tolerance range for Arrhenius thermal response for maturity maintenance (causes 'Temperature Size Rule' effect)
#' @param f = 1, functional response (-), usually kept at 1 because gut model controls food availability such that f=0 when gut empty
#' @param E_sm = 1116, Maximum volume-specific energy density of stomach (J/cm3)
#' @param K = 500, Half saturation constant (#/cm2)
#' @param X = 11, Food density (J/cm2)
#' @param andens_deb = 1, Animal density (g/cm3)
#' @param d_V = 0.3, Dry mass fraction of structure
#' @param d_E = 0.3, Dry mass fraction of reserve
#' @param d_Egg = 0.3, Dry mass fraction of egg
#' @param stoich_mode = 0, adjust chemical indices to chemical potentials (0) or vice versa (1), or leave as is (2)
#' @param mu_X = 525000, Molar Gibbs energy (chemical potential) of food (J/mol)
#' @param mu_E = 585000, Molar Gibbs energy (chemical potential) of reserve (J/mol)
#' @param mu_V = 500000, Molar Gibbs energy (chemical potential) of structure (J/mol)
#' @param mu_P = 480000, Molar Gibbs energy (chemical potential) of faeces (J/mol)
#' @param mu_N = 244e3/5, Molar Gibbs energy (chemical potential) of nitrogenous waste (J/mol), synthesis from NH3, Withers page 119
#' @param kap_X_P = 0.1, Faecation efficiency of food to faeces (-)
#' @param n_X = c(1,1.8,0.5,.15), Chem. indices of C, O, H and N in food
#' @param n_E = c(1,1.8,0.5,.15), Chem. indices of C, O, H and N in reserve
#' @param n_V = c(1,1.8,0.5,.15), Chem. indices of C, O, H and N in structure
#' @param n_P = c(1,1.8,0.5,.15), Chem. indices of C, O, H and N in faeces
#' @param fdry = 0.3, Dry mass fraction of food
#' @param n_M_nitro = c(1,4/5,3/5,4/5), Chem. indices of C, O, H and N in nitrogenous waste
#' @param h_N = 384238, molar enthalpy of nitrogenous waste (combustion frame of reference) (J/mol), overridden if n_M_nitro specified as urea, uric acid or ammonia
#' @param stages = 3, how many life stages?
#' @param stage = 0, Initial stage (0=embryo, for STD 1=juvenile, 2=mature but not yet reproducing, 3=beyond first reproduction, for ABP 1-(stages-1) = instars, stages = adult)
#' @param S_instar = rep(1.6, stages), stress at instar n: L_n^2/ L_n-1^2 (-)
#' @param clutchsize = 2, Clutch size (#), overridden by \code{clutch_ab}
#' @param clutch_ab = c(0,0), paramters for relationship between length (cm) and clutch size: clutch size = a*L_w-b, make a and b zero if fixed clutch size
#' @param minclutch = 0, Minimum clutch size if not enough in reproduction buffer for clutch size predicted by \code{clutch_ab} - if zero, will not operate
#' @param batch = 1, Invoke Pequerie et al.'s batch laying model?
#' @param lambda = 1/2
#' @param VTMIN = 26, Voluntary thermal maximum, degrees C, controls whether repro event can occur at a given time
#' @param VTMAX = 39, Voluntary thermal maximum, degrees C, controls whether repro event can occur at a given time
#' @param arrhenius = matrix(data = matrix(data = c(rep(T_A,8),rep(T_AL,8),rep(T_AH,8),rep(T_L,8),rep(T_H,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5), Stage-specific 5-parameter Arrhenius thermal response for DEB model (T_A,T_AL,T_AH,T_L,T_H)
#' @param arrhenius2 = matrix(data = matrix(data = c(rep(T_A2,8),rep(T_AL2,8),rep(T_AH2,8),rep(T_L2,8),rep(T_H2,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5), Stage-specific 5-parameter Arrhenius thermal response for maturity maintenance (causes 'Temperature Size Rule' effect) DEB model (T_A2,T_AL2,T_AH2,T_L2,T_H2)
#' @param acthr = 1
#' @param X = 11
#' @param E_pres = 6011.93
#' @param V_pres = 3.9752^3
#' @param E_H_pres = 73592
#' @param q_pres =0
#' @param hs_pres =0
#' @param p_surv = 1
#' @param E_s_pres = 0
#' @param E_R = 0
#' @param E_B = 0
#' @param stage = 0
#' @param breeding = 0
#' @param Tb = 33
#' @return stage Life cycle stage, -
#' @return V Structure, cm^3
#' @return E Reserve density, J/cm^3
#' @return E_H Maturity, J
#' @return E_s Stomach energy content, J
#' @return E_R Reproduction buffer energy, J
#' @return E_B Reproduction batch energy, J
#' @return q Aging acceleration
#' @return hs Hazard rate
#' @return length Physical length, cm
#' @return wetmass Total wet mass, g
#' @return wetgonad Wet mass of gonad, g
#' @return wetgut Wet mass of food in gut, g
#' @return wetstorage Wet mass of reserve, g
#' @return p_surv Survival probability, -
#' @return fecundity Eggs produced at a given time point, #
#' @return clutches Clutches produced at a given time point
#' @return JMO2 Oxygen flux, mol/time
#' @return JMCO2 Carbon dioxide flux, mol/time
#' @return JMH2O metabolic water flux, mol/time
#' @return JMNWASTE nitrogenous waste flux, mol/time
#' @return O2ML Oxgen consumption rate, ml/hour
#' @return CO2ML Carbon dioxide production rate, ml/time
#' @return GH2OML Metabolic water production rate, ml/time
#' @return DEBQMET Metabolic heat generation, J/time
#' @return GDRYFOOD Dry food intake, g/time
#' @return GFAECES Faeces production, g/time
#' @return GNWASTE Nitrogenous waste production, g/time
#' @return p_A Assimilation power, J/time
#' @return p_C Catabolic power, J/time
#' @return p_M Somatic maintenance power, J/time
#' @return p_G Growth power, J/time
#' @return p_D Dissipation power, J/time
#' @return p_J Maturity power, J/time
#' @return p_R Reproduction power, J/time
#' @return p_B Reproduction batch power, J/time
#' @return L_b Structural length at birth, cm
#' @return L_j Structural length at end of metabolic acceleration (if occurring), cm
#' @examples
#' # simulate growth and reproduction at different constant body temperatures at
#' # constant food for a lizard (Tiliqua rugosa - default parameter values, starting
#' # as an egg)
#'
#' n <- 3000 # time steps
#' step <- 1 # step size (days)
#'
#' Tbs=seq(25, 35, 5) # sequence of body temperatures to use
#'
#' for(j in 1:length(Tbs)){
#'   debout<-matrix(data = 0, nrow = n, ncol = 38)
#'   deb.names <- c("stage", "V", "E", "E_H", "E_s", "E_R", "E_B", "q", "hs", "length", "wetmass", "wetgonad", "wetgut", "wetstorage", "p_surv", "fecundity", "clutches", "JMO2", "JMCO2", "JMH2O", "JMNWASTE", "O2ML", "CO2ML", "GH2OMET", "DEBQMETW", "GDRYFOOD", "GFAECES", "GNWASTE", "p_A", "p_C", "p_M", "p_G", "p_D", "p_J", "p_R", "p_B", "L_b", "L_j")
#'   colnames(debout)<-deb.names
#'
#'   # initialise
#'   debout[1,]<-DEB(Tb = Tbs[j], step = step)
#'
#'   for(i in 2:n){
#'     debout[i,] <- DEB(Tb = Tbs[j], breeding = 1, step = step, E_pres = debout[(i - 1), 3], V_pres = debout[(i - 1), 2], E_H_pres = debout[(i - 1), 4], q = debout[(i - 1), 8], hs = debout[(i - 1), 9], p_surv = debout[(i - 1), 15], E_s_pres = debout[(i - 1), 5], E_R = debout[(i - 1), 6], E_B = debout[(i - 1), 7])
#'   }
#'
#'   if(j == 1){
#'     plot((seq(1, n) / 365), debout[, 11], ylim = c(100, 1500), type = 'l', xlab = 'years', ylab = 'wet mass, g', col = j)
#'   }else{
#'     points((seq(1,n) / 365), debout[, 11], ylim = c(100, 1500), type = 'l', xlab = 'years', ylab = 'wet mass, g',col = j)
#'   }
#'
#' } #end loop through body temperatures
#' @export
DEB<-function(
  step=1/24,
  z=7.997,
  del_M=0.242,
  p_Xm=13290*step,
  kap_X=0.85,
  v=0.065*step,
  kap=0.886,
  p_M=32*step,
  p_T=0,
  E_G=7767,
  kap_R=0.95,
  k_J=0.002*step,
  E_Hb=7.359e+04,
  E_Hj=E_Hb,
  E_Hp=1.865e+05,
  h_a=2.16e-11/(step^2),
  s_G=0.01,
  T_REF=20+273.15,
  T_A=8085,
  T_AL=18721,
  T_AH=9.0E+04,
  T_L=288,
  T_H=315,
  T_A2=T_A,
  T_AL2=T_AL,
  T_AH2=T_AH,
  T_L2=T_L,
  T_H2=T_H,
  E_0=1.04e+06,
  f=1,
  E_sm=1116,
  K=1,
  andens_deb=1,
  d_V=0.3,
  d_E=0.3,
  d_Egg=0.3,
  stoich_mode=0,
  mu_X=525000,
  mu_E=585000,
  mu_V=500000,
  mu_P=480000,
  mu_N=244e3/5,
  kap_X_P=0.1,
  n_X=c(1,1.8,0.5,0.15),
  n_E=c(1,1.8,0.5,0.15),
  n_V=c(1,1.8,0.5,0.15),
  n_P=c(1,1.8,0.5,0.15),
  n_M_nitro=c(1,4/5,3/5,4/5),
  h_N = 384238,
  clutchsize=2,
  clutch_ab=c(0.085,0.7),
  minclutch=0,
  batch=1,
  lambda=1/2,
  VTMIN=26,
  VTMAX=39,
  ma=1e-4,
  mi=0,
  mh=0.5,
  arrhenius=matrix(data = matrix(data = c(rep(T_A,8),rep(T_AL,8),rep(T_AH,8),rep(T_L,8),rep(T_H,8)),nrow = 8, ncol = 5), nrow = 8, ncol = 5),
  arrhenius2=matrix(data = matrix(data = c(rep(T_A2,8),rep(T_AL2,8),rep(T_AH2,8),rep(T_L2,8),rep(T_H2,8)),nrow = 8, ncol = 5), nrow = 8, ncol = 5),
  acthr=1,
  X=10,
  E_pres=E_0/3e-9,
  V_pres=3e-9,
  E_H_pres=0,
  q_pres=0,
  hs_pres=0,
  p_surv_pres=1,
  E_s_pres=0,
  p_B_pres=0,
  E_R=0,
  E_B=0,
  stages=7,
  stage=0,
  breeding=0,
  Tb=33,
  fdry=0.3,
  L_b=0.42,
  L_j=1.376,
  S_instar=rep(1.618, stages),
  spawnday=1,
  day=1,
  metab_mode=0,
  age=0){

  if (!require("deSolve", quietly = TRUE)) {
    stop("package 'deSolve' is needed. Please install it.",
         call. = FALSE)
  }

  # initialise for reproduction and starvation
  if(clutch_ab[1] > 0){
    clutchsize <- floor(clutch_ab[1] * (V_pres ^ (1 / 3) / del_M) - clutch_ab[2])
    clutchsize[clutchsize < 0] <- 0
  }
  orig_clutchsize <- clutchsize
  fecundity <- 0
  clutches <- 0
  clutchenergy <- E_0 * clutchsize
  starve <- 0
  p_B <- 0

  #DEB mass balance-related calculations
  if(stoich_mode == 0){
    # match H fraction in organics to stated chemical potentials (needed later for heat production)
    n_X[2] <- ((mu_X / 10 ^ 5) - 4.3842 * n_X[1] - (-1.8176) * n_X[3] - (0.0593) * n_X[4]) / 0.9823
    n_V[2] <- ((mu_V / 10 ^ 5) - 4.3842 * n_V[1] - (-1.8176) * n_V[3] - (0.0593) * n_V[4]) / 0.9823
    n_E[2] <- ((mu_E / 10 ^ 5) - 4.3842 * n_E[1] - (-1.8176) * n_E[3] - (0.0593) * n_E[4]) / 0.9823
    n_P[2] <- ((mu_P / 10 ^ 5) - 4.3842 * n_P[1] - (-1.8176) * n_P[3] - (0.0593) * n_P[4]) / 0.9823
  }else{
    if(stoich_mode == 1){
    # match stated chemical potentials to H fraction in organics
    mu_X <- (n_X[2] * 0.9823 + 4.3842 * n_X[1] + (-1.8176) * n_X[3] + (0.0593) * n_X[4]) * 10 ^ 5
    mu_V <- (n_V[2] * 0.9823 + 4.3842 * n_V[1] + (-1.8176) * n_V[3] + (0.0593) * n_V[4]) * 10 ^ 5
    mu_E <- (n_E[2] * 0.9823 + 4.3842 * n_E[1] + (-1.8176) * n_E[3] + (0.0593) * n_E[4]) * 10 ^ 5
    mu_P <- (n_P[2] * 0.9823 + 4.3842 * n_P[1] + (-1.8176) * n_P[3] + (0.0593) * n_P[4]) * 10 ^ 5
    }
  }
  # enthalpies (combustion frame)
  h_X <- 10^5 * (4.3284 * n_X[1] + 1.0994 * n_X[2] + (-2.0915) * n_X[3] + (-0.1510) * n_X[4]) #J mol^(-1)
  h_V <- 10^5 * (4.3284 * n_V[1] + 1.0994 * n_V[2] + (-2.0915) * n_V[3] + (-0.1510) * n_V[4]) #J mol^(-1)
  h_E <- 10^5 * (4.3284 * n_E[1] + 1.0994 * n_E[2] + (-2.0915) * n_E[3] + (-0.1510) * n_E[4]) #J mol^(-1)
  h_P <- 10^5 * (4.3284 * n_P[1] + 1.0994 * n_P[2] + (-2.0915) * n_P[3] + (-0.1510) * n_P[4]) #J mol^(-1)
  h_CO2 <- 0 #J mol^(-1)
  h_O2 <- 0 #J mol^(-1)
  h_H2O <- 0 #J mol^(-1)
  if(all(n_M_nitro == c(0, 3, 0, 1))){ # ammonia
    h_N <- 382805
    mu_N <- 0
  }
  if(all(n_M_nitro == c(1.0, 0.8, 0.6, 0.8))){ # uric acid
    h_N <- 384238
    mu_N <- 244e3/5
  }
  if(all(n_M_nitro == c(1, 2, 1, 2))){ # urea
    h_N <- 631890
    mu_N <- 122e3
  }
  h_O <- c(h_X, h_V, h_E, h_P)
  h_M <- c(h_CO2, h_H2O, h_O2, h_N)
  n_O <- cbind(n_X, n_V, n_E, n_P) # matrix of composition of organics, i.e. food, structure, reserve and faeces
  CHON <- c(12, 1, 16, 14)
  wO <- CHON %*% n_O
  w_V <- wO[2] # molar mass of structure
  M_V <- d_V / w_V # cmoles structure per volume structure
  y_EX <- kap_X * mu_X / mu_E # yield of reserve on food
  y_XE <- 1 / y_EX # yield of food on reserve
  y_VE <- mu_E * M_V / E_G  # yield of structure on reserve
  y_PX <- kap_X_P * mu_X / mu_P # yield of faeces on food
  y_PE <- y_PX / y_EX # yield of faeces on reserve
  nM <- matrix(c(1, 0, 2, 0, 0, 2, 1, 0, 0, 0, 2, 0, n_M_nitro), nrow = 4)
  n_M_nitro_inv <- c(-1 * n_M_nitro[1] / n_M_nitro[4], (-1 * n_M_nitro[2]) / (2 * n_M_nitro[4]), (4 * n_M_nitro[1] + n_M_nitro[2] - 2 * n_M_nitro[3]) / (4 * n_M_nitro[4]), 1 /n_M_nitro[4])
  n_M_inv <- matrix(c(1, 0, -1, 0, 0, 1 / 2, -1 / 4, 0, 0, 0, 1 / 2, 0, n_M_nitro_inv), nrow = 4)
  JM_JO <- -1 * n_M_inv %*% n_O
  etaO <- matrix(c(y_XE / mu_E * -1, 0, 1 / mu_E, y_PE / mu_E, 0, 0, -1 / mu_E, 0, 0, y_VE / mu_E, -1 / mu_E, 0), nrow = 4)
  w_N <- CHON %*% n_M_nitro

  # Arrhenius temperature correction factor
  #Tcorr <- exp(T_A * (1 / (273.15 + T_REF) - 1 / (273.15 + Tb))) / (1 + exp(T_AL * (1 / (273.15 + Tb) - 1 / T_L)) + exp(T_AH * (1 / T_H - 1 / (273.15 + Tb))))
  Tcorr <- exp(T_A / T_REF - T_A / (273.15 + Tb)) * (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)))
  Tcorr2 <- exp(T_A2 / T_REF - T_A2 / (273.15 + Tb)) * (1 + exp(T_AL2 / T_REF - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / T_REF)) / (1 + exp(T_AL2 / (273.15 + Tb) - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / (273.15 + Tb)))

  # metabolic acceleration if present
  s_M <- 1 # -, multiplication factor for v and {p_Am} under metabolic acceleration
  if(E_Hj != E_Hb){
    if(E_H_pres < E_Hb){
      s_M <- 1
    }else{
      if(E_H_pres < E_Hj){
        s_M <- V_pres ^ (1 / 3) / L_b
      }else{
        s_M <- L_j / L_b
      }
    }
  }

  # temperature corrections and compound parameters
  M_V <- d_V / w_V
  p_MT <- p_M * Tcorr
  k_M <- p_MT / E_G
  k_JT <- k_J * Tcorr2
  vT <- v * Tcorr * s_M
  p_AmT <- p_MT * z / kap * s_M
  p_XmT <- p_Xm * Tcorr * s_M
  h_aT <- h_a * Tcorr
  E_m <- p_AmT / vT
  g <- E_G / (kap * E_m) # energy investment ratio
  e <- E_pres / E_m # scaled reserve density
  V_m <- (kap * p_AmT / p_MT) ^ 3 # maximum structural volume
  L_T <- p_T / p_MT # heating length
  L_pres <- V_pres ^ (1 / 3)
  L_m <- V_m ^ (1 / 3)
  scaled_l <- L_pres / L_m
  kap_G <- (d_V * mu_V) / (w_V * E_G)
  yEX <- kap_X * mu_X / mu_E
  yXE <- 1 / yEX
  yPX <- kap_X_P * mu_X / mu_P
  mu_AX <- mu_E / yXE
  eta_PA <- yPX / mu_AX
  w_X <- wO[1]
  w_E <- wO[3]
  w_V <- wO[2]
  w_P <- wO[4]

  # initial conditions for solver
  init <- c(V_pres, E_pres, E_H_pres, E_s_pres, starve, q_pres, hs_pres, E_R, E_B)

  # parameters
  indata <- list(k_J = k_JT, p_Am = p_AmT, k_M = k_M, p_M = p_MT,
                 p_Xm = p_XmT, v = vT, E_m = E_m, L_m = L_m, L_T = L_T,
                 kap = kap, g = g, M_V = M_V, mu_E = mu_E,
                 mu_V = mu_V, d_V = d_V, w_V = w_V, acthr = acthr,
                 X = X, K = K, E_Hp = E_Hp, E_Hb = E_Hb, E_Hj = E_Hj, s_G = s_G, h_a = h_aT,
                 batch = batch, kap_R = kap_R, lambda = lambda,
                 breeding = breeding, kap_X = kap_X, f = f, E_sm = E_sm, s_M = s_M,
                 L_j = L_j, metab_mode = metab_mode)

  # function for solver (running for one time step)
  dget_DEB <- function(t, y, indata){
    with(as.list(c(indata, y)), {

      # unpack variables
      V <- y[1]# cm^3, structural volume
      E <- y[2]# J/cm3, reserve density
      H <- y[3]# J, maturity
      E_s <- max(0, y[4])# J, stomach energy
      S <- y[5]# J, starvation energy
      q <- y[6]# -, aging acceleration
      hs <- y[7]# -, hazard rate
      R <- y[8]# J, reproduction buffer energy
      B <- y[9]# J, egg batch energy

      L <- V ^ (1/3) # cm, structural length
      V_m <- L_m ^ 3 # cm ^ 3, maximum structural volume
      e <- E / E_m  # -, scaled reserve density
      r <- v * (e / L - (1 + L_T / L) / L_m) / (e + g) # specific growth rate
      p_C <- E * V * (v / L - r) # J / t, mobilisation rate, equation 2.12 DEB3
      if(metab_mode == 1 & H >= E_Hj){
        r <- min(0, r) # no growth in abp after puberty, but could still be negative because starving
        p_C <- E * V * v / L
      }
      dV <- V * r # cm^3 / t, change in structure

      if(H < E_Hb){ # embryo
        # structure
        dE <- (- 1 *  E * v) / L
        dH <- (1 - kap) * p_C - k_J * H # J/d, change in maturity
        p_J <- k_J * H
        p_R <- (1 - kap) * p_C - p_J
        p_B <- 0

        # no aging or stomach in embryo
        dS <- 0
        dEs <- 0
        dq <- 0
        dhs <- 0
        dR <- p_R
        dB <- 0
      }else{ # post-embryo

        # structure and starvation
        if(V * r < 0){
          dS <- V * r * -1 * mu_V * d_V / w_V # J / t, starvation energy to be subtracted from reproduction buffer if necessary
          dV <- 0
          if(B + R < dS){ # reproduction and batch buffer has run out so draw from structure
            dV <- V * r
            dS <- 0
          }
        }else{
          dS <- 0
        }

        # assimilation
        p_A <- p_Am * f * L ^ 2

        # reserve
        if(E_s > p_A){
          dE <- p_A / L ^ 3 - (E * v) / L
        }else{
          dE <- max(0, E_s / L ^ 3) - (E * v) / L
        }

        if(metab_mode == 1 & H >= E_Hj){
         p_C <- p_A - dE * V
        }

        # maturation
        p_J <- k_J * H
        if(H < E_Hp){
          dH <- (1 - kap) * p_C - p_J
        }else{
          dH <- 0
        }

        # feeding
        if(acthr > 0){
          # Regulates X dynamics
          p_X <- f * p_Xm * ((X / K) / (1 + X / K)) * V ^ (2 / 3)
        }else{
          p_X <- 0
        }
        dEs <- p_X - (p_Am / kap_X) * V ^ (2 / 3)

        if(metab_mode == 1 & H >= E_Hj){
         r <- 0 # no growth in abp after puberty - not setting this to zero messes up ageing calculation
        }

        # ageing (equation 6.2 in Kooijman 2010 (DEB3)
        dq <- (q * (V / V_m) * s_G + h_a) * e * ((v / L) - r) - r * q # ageing acceleration
        dhs <- q - r * hs # hazard

        # reproduction
        if(metab_mode == 1 & H >= E_Hj){
          p_R <- (1 - kap) * p_A - p_J
        }else{
          p_R <- (1 - kap) * p_C - p_J
        }

        if(R <= 0 & B <= 0 & S > 0 &  p_R < S){
         dV <- -1 * abs(p_R) * w_V / (mu_V * d_V)  # subtract from structure since not enough flow to reproduction to pay for pay for somatic maintenance
         p_R <- 0
        }
        if(H < E_Hp){
          p_B <- 0
        }else{
          if(batch == 1){
             batchprep <- (kap_R / lambda) * ((1 - kap) * (E_m * (v * V ^ (2 / 3) + k_M * V) / (1 + (1 / g))) - p_J)
            if(breeding == 0){
              p_B <- 0
            }else{
              #if the repro buffer is lower than what p_B would be (see below), p_B is p_R
              if(R < batchprep){
                p_B <- p_R
              }else{
                #otherwise it is a faster rate, as specified in Pecquerie et. al JSR 2009 Anchovy paper,
                #with lambda (the fraction of the year the animals breed if food/temperature not limiting) = 0.583 or 7 months of the year
                p_B <- max(batchprep, p_R)
              }
            }
          }else{
            p_B <- p_R
          }#end check for whether batch mode is operating
        }#end check for immature or mature
        p_R <- max(0, p_R - p_B) # take finalised value of p_B from p_R

        # draw from reproduction and then batch buffers under starvation
        if(dS > 0 & R > dS){
         p_R <- p_R - dS
         dS <- 0
        }
        if(dS > 0 & B > dS){
         p_B <- p_B - dS
         dS <- 0
        }
        #accumulate energy/matter in reproduction and batch buffers
        dR <- p_R
        dB <- p_B * kap_R
      }

      y = list(c(dV, dE, dH, dEs, dS, dq, dhs, dR, dB))
    })
  }

  DEB.state <- as.data.frame(deSolve::ode(y = init, times = c(0, 1), func = dget_DEB, parms = indata, method = "ode45"))[2,2:10]
  colnames(DEB.state) <- c("V", "E", "H", "E_s", "S", "q", "hs", "R", "B")
  V <- max(DEB.state$V, 0)
  E <- max(DEB.state$E, 0)
  E_H <- max(DEB.state$H, 0)
  E_s <- max(DEB.state$E_s, 0)
  if(E_s > (E_sm * V)){
   resid <- E_s - E_sm * V # excess food intake to stomach capacity
   E_s <- E_sm * V
  }else{
   resid <- 0
  }
  starve <- max(DEB.state$S, 0)
  q <- max(DEB.state$q, 0)
  hs <- max(DEB.state$hs, 0)
  p_R <- max(DEB.state$R - E_R, 0)
  p_B <- max(DEB.state$B - E_B, 0) / kap_R
  if(E_H >= E_Hp){
   E_R <- max(DEB.state$R, 0)
   E_B <- max(DEB.state$B, 0)
  }
  e <- E / E_m # use new value of e
  L_w = V ^ (1 / 3) / del_M * 10 # length in mm
  if(E_H >= E_Hb & E_H_pres < E_Hb){
    L_b <- V ^ (1 / 3)
  }
  if(E_H >= E_Hj & E_H_pres < E_Hj){
    L_j <- V ^ (1 / 3)
  }
  # some powers
  p_M2 <- max(0, p_MT * V + p_T * V ^ (2 / 3))

  p_J <- k_JT * E_H
  if(E_s > V ^ (2 / 3) * p_AmT * f){
    p_A <- V ^ (2 / 3) * p_AmT * f
  }else{
    p_A <- E_s
  }
  r <- vT * (e / V ^ (1 / 3) - (1 + L_T / V ^ (1 / 3)) / L_m) / (e + g)
  p_C <- E * (vT / V ^ (1 / 3) - r) * V # J / t, mobilisation rate, equation 2.12 DEB3
  if(metab_mode == 1){
    if(E_H >= E_Hj){
      p_A <- p_R + p_B + p_M2 + p_J + (E_pres - E) * V
      p_C <- p_A - (E_pres - E) * V
    }
  }
  p_A <- max(0, p_A)
  p_C <- max(0, p_C)

  if(E_H >= E_Hp){
    p_D <- p_M2 + p_J + (1 - kap_R) * p_B
  }else{
    p_D <- p_M2 + p_J + p_R
  }
  if(metab_mode == 1 & E_H >= E_Hj){
    p_G <- 0
  }else{
    p_G <- p_C - p_M2 - p_J - p_R - p_B
  }

  # J food eaten per hour
  if(acthr > 1){
   p_X <- f * p_XmT * V ^ (2 / 3) * (X / K) / (1 + X / K) - resid
  }else{
   p_X <- 0
  }

  testclutch <- floor((E_R + E_B) / E_0)
  # FOR VARIABLE CLUTCH SIZE FROM REPRO AND BATCH BUFFERS
  if(minclutch > 0 & floor(E_R + E_B) / E_0 > minclutch){
    if(testclutch <= orig_clutchsize){# ! MAKE SMALLEST CLUTCH ALLOWABLE FOR THIS REPRO EVENT
      clutchsize <- minclutch
      clutchenergy <- clutchsize * E_0
    }
  }

  # determine stages

  # STD MODEL
  if(metab_mode == 0 & E_Hb == E_Hj){
    if(E_H < E_Hb){
      stage <- 0
    }else{
      if(E_H < E_Hp){
        stage <- 1
      }else{
        stage <- 2
      }
    }
    if(E_B > 0){
      if(E_H >= E_Hp){
        stage <- 3
      }else{
        stage <- stage
      }
    }
  }

  # ABJ MODEL
  if(metab_mode == 0 & E_Hb != E_Hj){
    if(E_H < E_Hb){
      stage <- 0
    }else{
      if(E_H < E_Hj){
        stage <- 1
      }
      if(E_H >= E_Hj){
        stage <- 2
      }
      if(E_H >= E_Hp){
        stage <- 3
      }
    }
    if(E_B > 0){
      if(E_H >= E_Hp){
        stage <- 4
      }else{
        stage <- stage
      }
    }
  }

  # ABP acceleration model
  if(metab_mode == 1){
    L_instar <- rep(0, stages)
    L_instar[1] <- S_instar[1] ^ 0.5 * L_b
    for(j in 2:stages){
      L_instar[j] <- S_instar[j] ^ 0.5 * L_instar[j - 1]
    }
    L_thresh <- L_instar[stage]
    if(stage == 0){
      if(E_H >= E_Hb){
        stage <- stage + 1
      }
    }else{
      if(stage < stages - 1){
        if(V^(1/3) > L_thresh){
          stage <- stage + 1
        }
      }
      if(E_H >= E_Hp){
        stage <- stages
      }
    }
  }


  if(E_B>clutchenergy){
    if((Tb >= VTMIN)  |  (Tb <= VTMAX)){
      if(day == spawnday & spawnday != 0){
        testclutch <- floor(E_B / E_0)
        if(testclutch > clutchsize){
          clutchsize <- testclutch
          clutchenergy <- clutchsize * E_0
        }
        if(spawnday > 0){
          clutchsize <- testclutch
          clutchenergy <- clutchsize * E_0
        }
        E_B <- E_B - clutchenergy
        repro <- 1
        fecundity <- clutchsize
        clutches <- 1
      }else{
        if(spawnday == 0){
          testclutch <- floor(E_B / E_0)
          if(testclutch > clutchsize){
            clutchsize <- testclutch
            clutchenergy <- clutchsize * E_0
          }
          E_B <- E_B - clutchenergy
          repro <- 1
          fecundity <- clutchsize
          clutches <- 1
        }
      }
    }
  }

  #mass balance
  JOJx <- p_A * etaO[1,1] + p_D * etaO[1,2] + p_G * etaO[1,3] # molar flux of food (mol/time step)
  JOJv <- p_A * etaO[2,1] + p_D * etaO[2,2] + p_G * etaO[2,3] # molar flux of reserve (mol/time step)
  JOJe <- p_A * etaO[3,1] + p_D * etaO[3,2] + p_G * etaO[3,3] # molar flux of structure (mol/time step)
  JOJp <- p_A * etaO[4,1] + p_D * etaO[4,2] + p_G * etaO[4,3] # molar flux of faeces (mol/time step)

  JOJx_GM <- p_D * etaO[1,2] + p_G * etaO[1,3] # non-assimilation (i.e. growth and maintenance) molar flux of food (mol/time step)
  JOJv_GM <- p_D * etaO[2,2] + p_G * etaO[2,3] # non-assimilation (i.e. growth and maintenance) molar flux of reserve (mol/time step)
  JOJe_GM <- p_D * etaO[3,2] + p_G * etaO[3,3] # non-assimilation (i.e. growth and maintenance) molar flux of structure (mol/time step)
  JOJp_GM <- p_D * etaO[4,2] + p_G * etaO[4,3] # non-assimilation (i.e. growth and maintenance) molar flux of faeces (mol/time step)

  JMCO2 <- JOJx * JM_JO[1,1] + JOJv * JM_JO[1,2] + JOJe * JM_JO[1,3] + JOJp * JM_JO[1,4] # molar flux of CO2 (mol/time step)
  JMH2O <- JOJx * JM_JO[2,1] + JOJv * JM_JO[2,2] + JOJe * JM_JO[2,3] + JOJp * JM_JO[2,4] # molar flux of H2O (mol/time step)
  JMO2 <- JOJx * JM_JO[3,1] + JOJv * JM_JO[3,2] + JOJe * JM_JO[3,3] + JOJp * JM_JO[3,4] # molar flux of O2 (mol/time step)
  JMNWASTE <- JOJx * JM_JO[4,1] + JOJv * JM_JO[4,2] + JOJe * JM_JO[4,3] + JOJp * JM_JO[4,4] # molar flux of nitrogenous waste (mol/time step)

  JMCO2_GM <- JOJx_GM * JM_JO[1,1] + JOJv_GM * JM_JO[1,2] + JOJe_GM * JM_JO[1,3] + JOJp_GM * JM_JO[1,4]
  JMH2O_GM <- JOJx_GM * JM_JO[2,1] + JOJv_GM * JM_JO[2,2] + JOJe_GM * JM_JO[2,3] + JOJp_GM * JM_JO[2,4]
  JMO2_GM <- JOJx_GM * JM_JO[3,1] + JOJv_GM * JM_JO[3,2] + JOJe_GM * JM_JO[3,3] + JOJp_GM * JM_JO[3,4]
  JMNWASTE_GM <- JOJx_GM * JM_JO[4,1] + JOJv_GM * JM_JO[4,2] + JOJe_GM * JM_JO[4,3] + JOJp_GM * JM_JO[4,4]

  #RQ <- JMCO2 / JMO2 # respiratory quotient

  #PV=nRT
  #T=273.15 #K
  #R=0.082058 #L*atm/mol*K
  #n=1 #mole
  #P=1 #atm
  #V=nRT/P=1*0.082058*273.15=22.41414
  #T=293.15
  #V=nRT/P=1*0.082058*293.15/1=24.0553
  P_atm <- 1
  R_const <- 0.082058
  gas_cor <- R_const * T_REF / P_atm * (Tb + 273.15) / T_REF * 1000 # 1 mole to ml/time at Tb and atmospheric pressure
  O2ML <- -1 * JMO2 * gas_cor # mlO2/time, temperature corrected (including SDA)
  CO2ML <- JMCO2 * gas_cor # mlCO2/time, temperature corrected (including SDA)
  GH2OMET <- JMH2O * 18.01528 # g metabolic water/time

  #metabolic heat production (Watts) - growth overhead plus dissipation power (maintenance, maturity maintenance,
  #maturation/repro overheads) plus assimilation overheads
  #DEBQMETW <- ((1 - kap_G) * p_G + p_D + (p_A / kap_X - p_A - p_A * mu_P * eta_PA)) / 3600 / Tcorr
  mu_O <- c(mu_X, mu_V, mu_E, mu_P) # J/mol, chemical potentials of organics
  mu_M <- c(0, 0, 0, mu_N)          # J/mol, chemical potentials of minerals C: CO2, H: H2O, O: O2, N: nitrogenous waste
  J_O <- c(JOJx, JOJv, JOJe, JOJp) # eta_O * diag(p_ADG(2,:)); # mol/d, J_X, J_V, J_E, J_P in rows, A, D, G in cols
  J_M <- c(JMCO2, JMH2O, JMO2, JMNWASTE) # - (n_M\n_O) * J_O;        # mol/d, J_C, J_H, J_O, J_N in rows, A, D, G in cols

  # compute heat production
  p_T <- sum(-1 * J_O * h_O -J_M * h_M) / 3600
  DEBQMETW <- p_T

  GDRYFOOD <- -1 * JOJx * w_X
  GFAECES <- JOJp * w_P
  GNWASTE <- JMNWASTE * w_N
  wetgonad <- ((E_R / mu_E) * w_E) / d_Egg + ((E_B / mu_E) * w_E) / d_Egg
  wetstorage <- ((V * E / mu_E) * w_E) / d_E
  wetgut <- ((E_s / mu_E) * w_E) / fdry
  wetmass <- V * andens_deb + wetgonad + wetstorage + wetgut

  dsurvdt <- -1 * p_surv_pres * hs
  p_surv <- p_surv_pres + dsurvdt

  # new states
  E_pres <- E
  V_pres <- V
  E_H_pres <- E_H
  q_pres <- q
  hs_pres <- hs
  p_surv_pres <- p_surv
  E_s_pres <- E_s

  deb.names <- c("stage", "V", "E", "E_H", "E_s", "E_R", "E_B", "q", "hs", "length", "wetmass", "wetgonad", "wetgut", "wetstorage", "p_surv", "fecundity", "clutches", "JMO2", "JMCO2", "JMH2O", "JMNWASTE", "O2ML", "CO2ML", "GH2OMET", "DEBQMETW", "GDRYFOOD", "GFAECES", "GNWASTE", "p_A", "p_C", "p_M", "p_G", "p_D", "p_J", "p_R", "p_B", "L_b", "L_j")
  results.deb <- c(stage, V_pres, E_pres, E_H_pres, E_s_pres, E_R, E_B, q_pres, hs_pres, L_w, wetmass, wetgonad, wetgut, wetstorage, p_surv, fecundity, clutches, JMO2, JMCO2, JMH2O, JMNWASTE, O2ML, CO2ML, GH2OMET, DEBQMETW, GDRYFOOD, GFAECES, GNWASTE, p_A, p_C, p_M2, p_G, p_D, p_J, p_R, p_B, L_b, L_j)
  names(results.deb) <- deb.names
  return(results.deb)
}
