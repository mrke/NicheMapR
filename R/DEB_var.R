#' Dynamic Energy Budget model
#'
#' Implementation of the Standard Dynamic Energy Budget model of Kooijman
#' Note that this uses the deSolve package 'ode' function with events and
#' can handle variable food and variable temperature. It runs faster than the 'DEB'
#' function.
#' Michael Kearney May 2021
#' @param ndays = 365, length of simulation (days)
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
#' @param acthr = 1
#' @param E_init = 6011.93
#' @param V_init = 3.9752^3
#' @param E_H_init = 73592
#' @param q_init = 0
#' @param hs_init = 0
#' @param p_surv_init = 1
#' @param E_s_init = 0
#' @param E_R_init = 0
#' @param E_B_init = 0
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
#' # simulate growth and reproduction at fluctuating body temperatures (Tb = 5 cm air temperature) at
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
DEB_var<-function(
  ndays=365,
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
  E_He=E_Hp,
  h_a=2.16e-11*(step^2),
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
  batch=1,
  lambda=1/2,
  E_init=E_0/3e-9,
  V_init=3e-9,
  E_H_init=0,
  q_init=0,
  hs_init=0,
  p_surv_init=1,
  E_s_init=0,
  p_B_init=0,
  E_R_init=0,
  E_B_init=0,
  stages=7,
  fdry=0.3,
  L_b=0.07101934,
  L_j=1.376,
  s_j=0.9977785,
  k_EV=0.002948759*step,
  kap_V=0.8,
  S_instar=rep(1.618, stages),
  day=1,
  metab_mode=0,
  age=0,
  foetal=0){

  if (!require("deSolve", quietly = TRUE)) {
    stop("package 'deSolve' is needed. Please install it.",
         call. = FALSE)
  }
  L.b <- L_b
  L.j <- L_j
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
  #Tcorr <- exp(T_A / T_REF - T_A / (273.15 + Tb)) * (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)))
  #Tcorr2 <- exp(T_A2 / T_REF - T_A2 / (273.15 + Tb)) * (1 + exp(T_AL2 / T_REF - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / T_REF)) / (1 + exp(T_AL2 / (273.15 + Tb) - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / (273.15 + Tb)))

  # temperature corrections and compound parameters
  M_V <- d_V / w_V
  p_Am <- p_M * z / kap
  k_M <- p_M / E_G
  E_m <- p_Am / v
  g <- E_G / (kap * E_m) # energy investment ratio
  e <- E_init / E_m # scaled reserve density
  V_m <- (kap * p_Am / p_M) ^ 3 # maximum structural volume
  L_T <- p_T / p_M # heating length
  L_init <- V_init ^ (1 / 3)
  L_m <- V_m ^ (1 / 3)
  scaled_l <- L_init / L_m
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


  breeding <- 1

  # general DEB function for solver (std, abj, abp)
  dget_DEB <- function(t, y, indata){
    with(as.list(c(indata, y)), {

      # temperature correction
      Tb <- Tbf(t)
      X <- Xf(t)
      Tcorr <- exp(T_A / T_REF - T_A / (273.15 + Tb)) * (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)))
      Tcorr2 <- exp(T_A2 / T_REF - T_A2 / (273.15 + Tb)) * (1 + exp(T_AL2 / T_REF - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / T_REF)) / (1 + exp(T_AL2 / (273.15 + Tb) - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / (273.15 + Tb)))

      # unpack variables
      V <- max(0, y[1])# cm^3, structural volume
      E <- max(0, y[2])# J/cm3, reserve density
      H <- max(0, y[3])# J, maturity
      E_s <- max(0, y[4])# J, stomach energy
      S <- max(0, y[5])# J, starvation energy
      q <- y[6]# -, aging acceleration
      hs <- y[7]# -, hazard rate
      R <-  max(0, y[8])# J, reproduction buffer energy
      B <-  max(0, y[9])# J, egg batch energy

      s_M <- 1

      if(E_Hj != E_Hb){
        if(H < E_Hb){
          s_M <- 1
        }else{
          if(H < E_Hj){
            s_M <- V ^ (1 / 3) / L_b
          }else{
            s_M <- L_j / L_b
          }
        }
      }else{
        s_M <- 1
      }

      L <- V ^ (1 / 3) # cm, structural length
      V_m <- (kap * (p_Am * s_M) / p_M) ^ 3 # cm ^ 3, maximum structural volume
      L_m <- V_m ^ (1 / 3)
      e <- E / E_m  # -, scaled reserve density
      r <- ((v * s_M) * (e / L - (1 + L_T / L) / L_m) / (e + g)) * Tcorr # specific growth rate
      p_C <- E * V * ((v * Tcorr * s_M) / L - r) # J / t, mobilisation rate, equation 2.12 DEB3
      if(metab_mode == 1 & H >= E_Hj){
        r <- min(0, r) # no growth in abp after puberty, but could still be negative because starving
        p_C <- E * V * (v * Tcorr * s_M) / L
      }
      dV <- V * r # cm^3 / t, change in structure

      if(H < E_Hb){ # embryo
        # structure
        dE <- (- 1 *  E * v * Tcorr) / L
        dH <- (1 - kap) * p_C - k_J * H * Tcorr2 # J/d, change in maturity
        p_J <- k_J * H * Tcorr2
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
        p_A <- (p_Am * Tcorr * s_M) * f * L ^ 2

        # reserve
        if(E_s > p_A){
          dE <- max(0, p_A / L ^ 3) - (E * (v * Tcorr * s_M)) / L
        }else{
          dE <- max(0, E_s / L ^ 3) - (E * (v * Tcorr * s_M)) / L
        }

        if(metab_mode == 1 & H >= E_Hj){
          p_C <- p_A - dE * V
        }
        p_C <- max(0, p_C)

        # maturation
        p_J <- k_J * H * Tcorr2
        if(H < E_Hp){
          dH <- (1 - kap) * p_C - p_J
        }else{
          dH <- 0
        }
        acthr <- 1
        # feeding
        if(acthr > 0){
          # Regulates X dynamics
          p_X <- f * (p_Xm * Tcorr * s_M) * ((X / K) / (1 + X / K)) * V ^ (2 / 3)
        }else{
          p_X <- 0
        }
        dEs <- p_X - ((p_Am * Tcorr * s_M) / kap_X) * V ^ (2 / 3)

        if(metab_mode == 1 & H >= E_Hj){
          r <- 0 # no growth in abp after puberty - not setting this to zero messes up ageing calculation
        }

        # ageing (equation 6.2 in Kooijman 2010 (DEB3)
        dq <- (q * (V / V_m) * s_G + h_a * Tcorr) * e * (((v * Tcorr * s_M) / L) - r) - r * q # ageing acceleration
        dhs <- q - r * hs * Tcorr # hazard

        # reproduction
        if(metab_mode == 1 & H >= E_Hj){
          p_R <- (1 - kap) * p_A - p_J
        }else{
          p_R <- (1 - kap) * p_C - p_J
        }
        p_R <- max(0, p_R)
        if(R <= 0 & B <= 0 & S > 0 & p_R < S){
          dV <- -1 * abs(p_R) * w_V / (mu_V * d_V)  # subtract from structure since not enough flow to reproduction to pay for pay for somatic maintenance
          p_R <- 0
        }

        if(H < E_Hp){
          p_B <- 0
        }else{
          if(batch == 1){
            batchprep <- max(0, (kap_R / lambda) * ((1 - kap) * (E_m * ((v * Tcorr * s_M) * V ^ (2 / 3) + k_M * Tcorr * V) / (1 + (1 / g))) - p_J))
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
        p_B <- max(0, p_B)
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
        dB <- p_B
      }

      y = list(c(dV, dE, dH, dEs, dS, dq, dhs, dR, dB))
    })
  }

  # embryo for hex model
  dget_ELH <- function(t, y, indata){
    with(as.list(c(indata, y)), {

      # temperature correction
      Tb <- Tbf(t)
      Tcorr <- exp(T_A / T_REF - T_A / (273.15 + Tb)) * (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)))
      Tcorr2 <- exp(T_A2 / T_REF - T_A2 / (273.15 + Tb)) * (1 + exp(T_AL2 / T_REF - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / T_REF)) / (1 + exp(T_AL2 / (273.15 + Tb) - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / (273.15 + Tb)))

      # unpack variables
      E <- y[1] # J, RESERVE
      L <- y[2] # CM, STRUCTURAL LENGTH
      H <- y[3] # J, MATURITY

      # use embryo equation for length, from Kooijman 2009 eq. 2
      V <- L ^ 3                                                # CM^3, STRUCTURAL VOLUME
      E_s <- E / V / E_m                                        # -, SCALED RESERVE DENSITY
      dL <- (v * E_s - k_M * Tcorr * g * L) / (3 * (E_s + g))   # CM/TIME, CHANGE IN LENGTH
      #r <- v * (E_s / L - 1 / L_m) / (E_s + g)                 # 1/TIME, GROWTH RATE
      SC <- L ^ 2 * (g * E_s)/(g + E_s) * (1 + (k_M * Tcorr * L) / v * Tcorr)
      dE <- -1 * SC * p_Am * Tcorr                              # J/TIME, CHANGE IN RESERVE
      U_H <- H / p_Am * Tcorr                                   # SCALED MATURITY
      dH <- ((1 - kap) * SC - k_J * Tcorr2 * U_H) * p_Am * Tcorr# J/TIME, CHANGE IN MATURITY

      y <- list(c(dE, dL, dH))
    })
  }

  # larva for hex model
  dget_AELES <- function(t, y, indata){
    with(as.list(c(indata, y)), {

      # temperature correction
      Tb <- Tbf(t)
      Tcorr <- exp(T_A / T_REF - T_A / (273.15 + Tb)) * (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)))
      Tcorr2 <- exp(T_A2 / T_REF - T_A2 / (273.15 + Tb)) * (1 + exp(T_AL2 / T_REF - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / T_REF)) / (1 + exp(T_AL2 / (273.15 + Tb) - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / (273.15 + Tb)))

      # food
      X <- Xf(t)

      # unpack variables
      E <- y[1] # J, RESERVE
      L <- y[2] # CM, STRUCTURAL LENGTH
      E_R <- y[3] #J, REPRODUCTION BUFFER
      E_S <- min(y[4], E_sm * (L ^ 3)) #J, STOMACH ENERGY

      V <- L ^ 3                                         # CM^3, STRUCTURAL VOLUME
      e <- E / V / E_m                                   # -, SCALED RESERVE DENSITY
      r <- (e * k_E * Tcorr - g * k_M * Tcorr)/ (e + g)  # 1/TIME, SPECIFIC GROWTH RATE
      p_C <- E * (k_E * Tcorr - r)                       # J/TIME, MOBILISATION RATE
      p_A <- f2 * p_Am * Tcorr * V                       # J/TIME, ASSIMILATION RATE, NOTE MULTPLYING BY V SINCE ALREADY DIVIDED BY L_B WHICH IS THE SAME AS MULT BY V^2/3 AND BY L/L_b
      p_X <- p_Xm * Tcorr * ((X / K) / (f2 + X / K)) * V # J/TIME, FOOD ENERGY INTAKE RATE, NOTE MULTPLYING BY V SINCE ALREADY DIVIDED BY L_B WHICH IS THE SAME AS MULT BY V^2/3 AND BY L/L_b
      #p_X <- p_A / kap_X                                # J/TIME, FOOD ENERGY INTAKE RATE, NOTE MULTPLYING BY V SINCE ALREADY DIVIDED BY L_B WHICH IS THE SAME AS MULT BY V^2/3 AND BY L/L_b
      if(E_S < p_A){                                     # NO ASSIMILATION IF STOMACH TOO EMPTY
        dE <- max(0, E_S) - p_C                                  # J/TIME, CHANGE IN RESERVE
      }else{
        dE <- f * p_Am * Tcorr * V - p_C                 # J/TIME, CHANGE IN RESERVE, NOTE MULTPLYING BY V SINCE ALREADY DIVIDED BY L_B WHICH IS THE SAME AS MULT BY V^2/3 AND BY L/L_b
      }
      dL <- r * L / 3                                    # CM/TIME, CHANGE IN LENGTH
      dER <- (1 - kap) * p_C - p_J * Tcorr2              # J/TIME, CHANGE IN REPROD BUFFER
      dEs <- p_X * Tcorr - f * (p_Am * Tcorr / kap_X) * V# J/TIME, CHANGE IN STOMACH ENERGY, NOTE MULTPLYING BY V SINCE ALREADY DIVIDED BY L_B WHICH IS THE SAME AS MULT BY V^2/3 AND BY L/L_b

      y <- list(c(dE, dL, dER, dEs))
    })
  }

  # pupa for hex model
  dget_AVELHS <- function(t, y, indata){
    with(as.list(c(indata, y)), {

      # temperature correction
      Tb <- Tbf(t)
      Tcorr <- exp(T_A / T_REF - T_A / (273.15 + Tb)) * (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)))
      Tcorr2 <- exp(T_A2 / T_REF - T_A2 / (273.15 + Tb)) * (1 + exp(T_AL2 / T_REF - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / T_REF)) / (1 + exp(T_AL2 / (273.15 + Tb) - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / (273.15 + Tb)))

      # unpack variables
      V <- y[1] # CM3, STRUCTURAL VOLUME OF LARVA
      E <- y[2] # J, RESERVE OF LARVA
      L <- y[3] # CM, STRUCTURAL LENGTH OF IMAGO
      H <- y[4] # J, MATURITY

      dV <- -1 * V * k_E * Tcorr                    # CM^3/TIME, CHANGE IN LARVAL STRUCTURAL VOLUME
      e <- E / L ^ 3 / E_m                          # -, SCALED RESERVE DENSITY
      r <- v_j * Tcorr * (e / L - 1 / L_m)/ (e + g) # 1/TIME, SPECIFIC GROWTH RATE
      p_C <- E * (v_j * Tcorr / L - r)              # J/TIME, MOBILISATION RATE
      dE <-  dV * g * E_m * kap * kap_V - p_C       # J/TIME, CHANGE IN RESERVE
      dL <- r * L / 3                               # CM/TIME, CHANGE IN LENGTH OF IMAGO
      dH <- (1 - kap) * p_C - k_J * Tcorr2 * H       # J/TIME, CHANGE IN MATURITY

      y <- list(c(dV, dE, dL, dH))
    })
  }

  # imago for hex model
  dget_EEES <- function(t, y, indata){
    with(as.list(c(indata, y)), {

      # temperature correction
      Tb <- Tbf(t)
      Tcorr <- exp(T_A / T_REF - T_A / (273.15 + Tb)) * (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)))
      Tcorr2 <- exp(T_A2 / T_REF - T_A2 / (273.15 + Tb)) * (1 + exp(T_AL2 / T_REF - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / T_REF)) / (1 + exp(T_AL2 / (273.15 + Tb) - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / (273.15 + Tb)))

      # food
      X <- Xf(t)

      # unpack variables
      E <- y[1]                           # J, RESERVE
      E_R <- y[2]                         # J, REPRODUCTION BUFFER
      E_B <- y[3]                         # J, EGG BATCH BUFFER
      E_s <- min(y[4], E_sm * (L ^ 3))    # J, ENERGY OF THE STOMACH
      q <- y[5]                           # -, aging acceleration
      hs <- y[6]                          # -, hazard rate
      V <- L ^ 3
      E_m <- p_Am / v * Tcorr             # -, SCALED RESERVE DENSITY
      e <- E / V / E_m
      p_C <- E * k_E * Tcorr              # J/TIME, RESERVE MOBILISATION
      p_R <- p_C - p_M - p_J * Tcorr2     # J/TIME, ENERGY ALLOCATION FROM RESERVE TO E_R
      p_X <- p_Xm * Tcorr * ((X / K) / (f2 + X / K)) * V ^ (2 / 3) # J/TIME, FOOD ENERGY INTAKE RATE
      if(breeding == 1){
        p_CR <- kap_R * p_R               # J/H, DRAIN FROM E_R TO EGGS
      }else{
        p_CR <- 0
      }
      if(E_s < p_A){      # NO ASSIMILATION IF STOMACH TOO EMPTY
        dE <- max(0, E_s) - p_C           # J/TIME, CHANGE IN RESERVE
      }else{
        dE <- f * p_Am * Tcorr * V - p_C  # J/TIME, CHANGE IN RESERVE
      }
      dER <- max(0, p_R - p_CR) # J/TIME, CHANGE IN REPROD BUFFER
      dEB <- p_CR                # J/TIME, CHANGE IN EGG BUFFER
      dES <- p_X - f * (p_Am / kap_X) * V # J/TIME, CHANGE IN STOMACH ENERGY
      r <- 0
      # ageing (equation 6.2 in Kooijman 2010 (DEB3)
      dq <- (q * (V / L_m ^ 3) * s_G + h_a * Tcorr) * e * (((v * Tcorr * s_M) / L) - r) - r * q # ageing acceleration
      dhs <- q - r * hs # hazard

      y <- list(c(dE, dER, dEB, dES, dq, dhs))
    })
  }

  # events

  eventfun <- function(t, y, pars) y

  birth <- function(t, y, pars) {
    if(y[3] > E_Hb) {
      y[3] <- 0
    }
    return(y)
  }
  puberty <- function(t, y, pars) {
    if(y[3] > E_Hp) {
      y[3] <- 0
    }
    return(y)
  }
  birth_stx <- function(t, y, pars) {
    if(y[1] > vHb) {
      y[1] <- 0
    }
    return(y)
  }
  metamorphosis <- function(t, y, pars) {
    if(y[3] > E_Hj) {
      y[3] <- 0
    }
    return(y)
  }
  pupate <- function(t, y, pars) {
    if(y[3] / (y[2] ^ 3) > E_RJ) {
      y[3] <- 0
    }
    return(y)
  }
  eclose <- function(t, y, pars) {
    if(y[4] > E_He) {
      y[4] <- 0
    }
    return(y)
  }

  if(metab_mode < 2){

    # parameters
    indata <- list(k_J = k_J,
                   p_Am = p_Am,
                   k_M = k_M,
                   p_M = p_M,
                   p_Xm = p_Xm,
                   v = v,
                   E_m = E_m,
                   L_m = L_m,
                   L_T = L_T,
                   kap = kap,
                   g = g,
                   M_V = M_V,
                   mu_E = mu_E,
                   mu_V = mu_V,
                   d_V = d_V,
                   w_V = w_V,
                   K = K,
                   E_Hp = E_Hp,
                   E_Hb = E_Hb,
                   E_Hj = E_Hj,
                   s_G = s_G,
                   h_a = h_a,
                   batch = batch,
                   kap_R = kap_R,
                   lambda = lambda,
                   breeding = breeding,
                   kap_X = kap_X,
                   f = f,
                   E_sm = E_sm,
                   L_b = L_b,
                   L_j = L_j,
                   metab_mode = metab_mode,
                   T_A = T_A,
                   T_L = T_L,
                   T_H = T_H,
                   T_AL = T_AL,
                   T_AH = T_AH,
                   T_REF = T_REF,
                   T_A2 = T_A2,
                   T_L2 = T_L2,
                   T_H2 = T_H2,
                   T_AL2 = T_AL2,
                   T_AH2 = T_AH2)

    times <- seq(0, (1 / step) * ndays)

    # initial conditions for solver
    init <- c(V_init, E_init, E_H_init + 1e-10, E_s_init + 1e-10, hs_init + 1e-10, q_init + 1e-10, hs_init + 1e-10, E_R_init + 1e-10, E_B_init + 1e-10)

    #"egg", "hatchling", "puberty", "adult"

    if(E_H_init < E_Hb){
      if(foetal == 1){
        indata$f
      }
      # to birth
      DEB.state.birth <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_DEB, parms = indata, method = "lsodes", events = list(func = eventfun, root = TRUE, terminalroot = 3), rootfunc = birth))[, 2:10]
      colnames(DEB.state.birth) <- c("V", "E", "H", "E_s", "S", "q", "hs", "R", "B")
      L.b <- max(DEB.state.birth$V ^ (1 / 3))
      t_birth <- which(DEB.state.birth$V ^ (1 / 3) == L.b)[1]
      DEB.state.birth <- DEB.state.birth[1:t_birth, ]
      init <- as.numeric(DEB.state.birth[nrow(DEB.state.birth), ])
      # if(nrow(DEB.state.birth) > 1){
      #   DEB.state.birth <- head(DEB.state.birth, -1)
      # }
      #times <- seq(0, (1 / step) * ndays - t_birth)
      times <- seq(t_birth, (1 / step) * ndays)

      # parameters
      indata <- list(k_J = k_J,
                     p_Am = p_Am,
                     k_M = k_M,
                     p_M = p_M,
                     p_Xm = p_Xm,
                     v = v,
                     E_m = E_m,
                     L_m = L_m,
                     L_T = L_T,
                     kap = kap,
                     g = g,
                     M_V = M_V,
                     mu_E = mu_E,
                     mu_V = mu_V,
                     d_V = d_V,
                     w_V = w_V,
                     K = K,
                     E_Hp = E_Hp,
                     E_Hb = E_Hb,
                     E_Hj = E_Hj,
                     s_G = s_G,
                     h_a = h_a,
                     batch = batch,
                     kap_R = kap_R,
                     lambda = lambda,
                     breeding = breeding,
                     kap_X = kap_X,
                     f = f,
                     E_sm = E_sm,
                     L_b = L.b,
                     L_j = L_j,
                     metab_mode = metab_mode,
                     T_A = T_A,
                     T_L = T_L,
                     T_H = T_H,
                     T_AL = T_AL,
                     T_AH = T_AH,
                     T_REF = T_REF,
                     T_A2 = T_A2,
                     T_L2 = T_L2,
                     T_H2 = T_H2,
                     T_AL2 = T_AL2,
                     T_AH2 = T_AH2)
    }
    if(E_H_init < E_Hj){
      if(E_Hb != E_Hj){
        # to metamorphosis
        DEB.state.meta <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_DEB, parms = indata, method = "lsodes", events = list(func = eventfun, root = TRUE, terminalroot = 3), rootfun = metamorphosis))[, 2:10]
        colnames(DEB.state.meta) <- c("V", "E", "H", "E_s", "S", "q", "hs", "R", "B")
        L.j <- max(DEB.state.meta$V ^ (1 / 3))
        t_meta <- which(DEB.state.meta$V ^ (1 / 3) == L.j)
        DEB.state.meta <- DEB.state.meta[1:t_meta, ]
        init <- as.numeric(DEB.state.meta[nrow(DEB.state.meta), ])
        if(nrow(DEB.state.meta) > 1){
          DEB.state.meta <- head(DEB.state.meta, -1)
        }
        indata <- list(k_J = k_J,
                       p_Am = p_Am,
                       k_M = k_M,
                       p_M = p_M,
                       p_Xm = p_Xm,
                       v = v,
                       E_m = E_m,
                       L_m = L_m,
                       L_T = L_T,
                       kap = kap,
                       g = g,
                       M_V = M_V,
                       mu_E = mu_E,
                       mu_V = mu_V,
                       d_V = d_V,
                       w_V = w_V,
                       K = K,
                       E_Hp = E_Hp,
                       E_Hb = E_Hb,
                       E_Hj = E_Hj,
                       s_G = s_G,
                       h_a = h_a,
                       batch = batch,
                       kap_R = kap_R,
                       lambda = lambda,
                       breeding = breeding,
                       kap_X = kap_X,
                       f = f,
                       E_sm = E_sm,
                       L_b = L.b,
                       L_j = L.j,
                       metab_mode = metab_mode,
                       T_A = T_A,
                       T_L = T_L,
                       T_H = T_H,
                       T_AL = T_AL,
                       T_AH = T_AH,
                       T_REF = T_REF,
                       T_A2 = T_A2,
                       T_L2 = T_L2,
                       T_H2 = T_H2,
                       T_AL2 = T_AL2,
                       T_AH2 = T_AH2)
        if(E_H_init < E_Hb){
          #times <- seq(0, (1 / step) * ndays - t_birth - t_meta)
          times <- seq(t_birth + t_birth, (1 / step) * ndays)
        }else{
          #times <- seq(0, (1 / step) * ndays - t_meta)
          times <- seq(t_meta, (1 / step) * ndays)
        }
      }
    }

    if(E_H_init < E_Hp & metab_mode != 1){

      # to puberty

      #DEB.state.pub <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_DEB, parms = indata, method = "lsodes", events = list(func = eventfun, root = TRUE, terminalroot = 3), rootfun = puberty))[, 2:10]
      DEB.state.pub <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_DEB, parms = indata, method = "lsoda", events = list(func = eventfun, root = TRUE, terminalroot = 3), rootfun = puberty))[, 2:10]
      colnames(DEB.state.pub) <- c("V", "E", "H", "E_s", "S", "q", "hs", "R", "B")
      L.p <- max(DEB.state.pub$V ^ (1 / 3))
      t_puberty <- which(DEB.state.pub$V ^ (1 / 3) == L.p)[1]
      DEB.state.pub <- DEB.state.pub[1:t_puberty, ]
      init <- as.numeric(DEB.state.pub[nrow(DEB.state.pub), ])
      if(nrow(DEB.state.pub) > 1){
        DEB.state.pub <- head(DEB.state.pub, -1)
      }
      if(E_Hb != E_Hj){

        if(E_H_init < E_Hb){
          #times <- seq(0, (1 / step) * ndays - t_birth - t_meta - t_puberty)
          times <- seq(t_birth + t_meta + t_puberty, (1 / step) * ndays)
        }else{
          #times <- seq(0, (1 / step) * ndays - t_meta - t_puberty)
          times <- seq(t_meta + t_puberty, (1 / step) * ndays)
        }
      }else{
        if(E_H_init < E_Hb){
          #times <- seq(0, (1 / step) * ndays - t_birth - t_puberty)
          times <- seq(t_birth + t_puberty, (1 / step) * ndays)
        }else{
          #times <- seq(0, (1 / step) * ndays - t_puberty)
          times <- seq(t_puberty, (1 / step) * ndays)
        }
      }
    }
    # to end of time sequence

    DEB.state.end <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_DEB, parms = indata, method = "lsoda"))[, 2:10]
    colnames(DEB.state.end) <- c("V", "E", "H", "E_s", "S", "q", "hs", "R", "B")

    if(E_Hb != E_Hj){
      if(metab_mode == 1){
        if(E_H_init < E_Hb){
          DEB.state <- rbind(DEB.state.birth, DEB.state.meta, DEB.state.end)
        }else{
          if(E_H_init < E_Hj){
            DEB.state <- rbind(DEB.state.meta, DEB.state.end)
          }else{
            if(E_H_init < E_Hp){
              DEB.state <- rbind(DEB.state.pub)
            }
          }
        }
      }else{
        if(E_H_init < E_Hb){
          DEB.state <- rbind(DEB.state.birth, DEB.state.meta, DEB.state.pub, DEB.state.end)
        }else{
          if(E_H_init < E_Hj){
            DEB.state <- rbind(DEB.state.meta, DEB.state.pub, DEB.state.end)
          }else{
            if(E_H_init < E_Hp){
              DEB.state <- rbind(DEB.state.pub, DEB.state.end)
            }else{
              DEB.state <- DEB.state.end
            }
          }
        }
      }
    }else{
      if(E_H_init < E_Hb){
        DEB.state <- rbind(DEB.state.birth, DEB.state.pub, DEB.state.end)
      }else{
        if(E_H_init < E_Hp){
          DEB.state <- rbind(DEB.state.pub, DEB.state.end)
        }else{
          DEB.state <- DEB.state.end
        }
      }
    }
  }else{
    # embryo

    # parameters
    indata <- list(f = f,
                   k_M = k_M,
                   v = v,
                   k_J = k_J,
                   E_m = E_m,
                   g = g,
                   kap = kap,
                   E_G = E_G,
                   p_M = p_M,
                   p_Am = p_Am,
                   L_m = L_m,
                   T_A = T_A,
                   T_L = T_L,
                   T_H = T_H,
                   T_AL = T_AL,
                   T_AH = T_AH,
                   T_REF = T_REF,
                   T_A2 = T_A2,
                   T_L2 = T_L2,
                   T_H2 = T_H2,
                   T_AL2 = T_AL2,
                   T_AH2 = T_AH2)

    times <- seq(0, (1 / step) * ndays)

    #if(E_H_init < E_Hb){
    # to birth

    # initial conditions for solver
    init <- c(E_init * V_init, V_init ^ (1 / 3), E_H_init + 1e-10)

    DEB.state.embryo <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_ELH, parms = indata, method = "lsodes", events = list(func = eventfun, root = TRUE, terminalroot = 3), rootfunc = birth))[, 2:4]
    colnames(DEB.state.embryo) <- c("E", "L", "H")
    L.b <- max(DEB.state.embryo$L)
    t_birth <- which(DEB.state.embryo$L == L.b)[1]
    DEB.state.embryo <- DEB.state.embryo[1:t_birth, ]
    init <- as.numeric(DEB.state.embryo[nrow(DEB.state.embryo), ])
    E_pres <- init[1]
    L_pres <- init[2]
    E_H_pres <- init[3]
    if(nrow(DEB.state.embryo) > 1){
      DEB.state.embryo <- head(DEB.state.embryo, -1)
    }
    #times <- seq(0, (1 / step) * ndays - t_birth)
    times <- seq(t_birth, (1 / step) * ndays)
    #}

    # larva

    p_J <- E_H_pres * k_J
    k_E <- v / L.b
    E_RJ <- s_j * ((1 - kap) * E_m * g *(((v / L.b) + (p_M / E_G))/ ((v / L.b) - g * (p_M / E_G))))

    # parameters
    indata <- list(f = f,
                   k_m = k_M,
                   k_E = k_E,
                   p_J = p_J,
                   p_Am = (p_M * z / kap) / L.b,
                   E_m = E_m,
                   g = g,
                   kap = kap,
                   p_Xm = p_Xm / L.b,
                   K = K,
                   f2 = 1,
                   kap_X = kap_X,
                   E_sm = E_sm,
                   E_RJ = E_RJ,
                   T_A = T_A,
                   T_L = T_L,
                   T_H = T_H,
                   T_AL = T_AL,
                   T_AH = T_AH,
                   T_REF = T_REF,
                   T_A2 = T_A2,
                   T_L2 = T_L2,
                   T_H2 = T_H2,
                   T_AL2 = T_AL2,
                   T_AH2 = T_AH2)

    # initial conditions for solver
    init <- c(E_pres, L_pres, E_R_init + 1e-10, E_s_init + 1e-10)

    # to pupation
    DEB.state.larva <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_AELES, parms = indata, method = "lsodes", events = list(func = eventfun, root = TRUE, terminalroot = 3), rootfunc = pupate))[, 2:5]
    colnames(DEB.state.larva) <- c("E", "L", "E_R", "E_s")
    L.j <- max(DEB.state.larva$L)
    t_pupate <- which(DEB.state.larva$L == L.j)
    DEB.state.larva <- DEB.state.larva[1:t_pupate, ]
    DEB.state.larva$E_s[DEB.state.larva$E_s > E_sm * DEB.state.larva$L ^ 3] <- E_sm * DEB.state.larva$L[DEB.state.larva$E_s > E_sm * DEB.state.larva$L ^ 3] ^ 3
    init <- as.numeric(DEB.state.larva[nrow(DEB.state.larva), ])
    if(nrow(DEB.state.larva) > 1){
      DEB.state.larva <- head(DEB.state.larva, -1)
    }
    E_pres <- init[1]
    L_pres <- init[2]
    E_R_pres <- init[3]
    E_s_pres <- init[4]
    #if(E_H_init < E_Hb){
    #times <- seq(0, (1 / step) * ndays - t_birth - t_pupate)
    times <- seq(t_birth + t_pupate, (1 / step) * ndays)
    #}else{
    #  times <- seq(0, (1 / step) * ndays - t_pupate)
    #}

    # pupa

    s_M <- L.j / L.b
    v_j <- v * s_M
    #k_EV <- k_EV * Tcorr#v_jT / L.j#k_EV * Tcorr
    L_m <- kap * p_Am * s_M / p_M
    # parameters
    indata <- list(f = f,
                   k_E = k_EV,
                   v_j = v_j,
                   E_m = E_m,
                   g = g,
                   kap = kap,
                   kap_V = kap_V,
                   k_J = k_J,
                   L_m = L_m,
                   E_He = E_He,
                   T_A = T_A,
                   T_L = T_L,
                   T_H = T_H,
                   T_AL = T_AL,
                   T_AH = T_AH,
                   T_REF = T_REF,
                   T_A2 = T_A2,
                   T_L2 = T_L2,
                   T_H2 = T_H2,
                   T_AL2 = T_AL2,
                   T_AH2 = T_AH2)
    # initial conditions for solver
    init <- c(L_pres ^ 3, E_pres, 1e-4, E_H_init + 1e-10)
    # to eclosion
    DEB.state.pupa <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_AVELHS, parms = indata, method = "lsodes", events = list(func = eventfun, root = TRUE, terminalroot = 4), rootfunc = eclose))[, 2:5]
    colnames(DEB.state.pupa) <- c("V", "E", "L", "H")
    # plot(DEB.state.pupa$V + DEB.state.pupa$L^3,type='l',ylim=c(0, max(DEB.state.pupa$V + DEB.state.pupa$L^3)))
    # points(DEB.state.pupa$V,type='l')
    # points(DEB.state.pupa$L^3,type='l')
    L.e <- max(DEB.state.pupa$L)
    t_eclose <- which(DEB.state.pupa$L == L.e)
    DEB.state.pupa <- DEB.state.pupa[1:t_eclose, ]
    init <- as.numeric(DEB.state.pupa[nrow(DEB.state.pupa), ])
    if(nrow(DEB.state.pupa) > 1){
      DEB.state.pupa <- head(DEB.state.pupa, -1)
    }
    V_old_pres <- init[1]
    E_pres <- init[2]
    L_pres <- init[3]
    H_pres <- init[4]
    #if(E_H_init < E_Hb){
    #times <- seq(0, (1 / step) * ndays - t_birth - t_pupate - t_eclose)
    times <- seq(t_birth + t_pupate + t_eclose, (1 / step) * ndays)
    #}else{
    #  times <- seq(0, (1 / step) * ndays - t_pupate - t_eclose)
    #}

    # imago
    V_pres <- L_pres ^ 3
    p_J <- E_H_pres * k_J
    p_A <- V_pres ^ (2 / 3) * p_Am * f
    # parameters
    indata <- list(f = f,
                   k_E = v / L.e,
                   p_J = p_J,
                   p_Am = p_Am * s_M,
                   p_M = V_pres * p_M,
                   p_A = p_A,
                   kap = kap,
                   p_Xm = p_Xm * s_M,
                   K = K,
                   f2 = 1,
                   kap_X = kap_X,
                   kap_R = kap_R,
                   v = v,
                   E_sm = E_sm,
                   s_G = s_G,
                   h_a = h_a,
                   L = L_pres,
                   breeding = breeding,
                   s_M = s_M,
                   E_m = E_m,
                   L = L.e,
                   L_m = L_m,
                   T_A = T_A,
                   T_L = T_L,
                   T_H = T_H,
                   T_AL = T_AL,
                   T_AH = T_AH,
                   T_REF = T_REF,
                   T_A2 = T_A2,
                   T_L2 = T_L2,
                   T_H2 = T_H2,
                   T_AL2 = T_AL2,
                   T_AH2 = T_AH2)

    init <- c(E_pres, E_R_pres, 1e-10, 1e-10, 1e-10, 1e-10)

    # to end
    DEB.state.imago <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_EEES, parms = indata, method = "lsodes"))[, 2:7]
    colnames(DEB.state.imago) <- c("E", "E_R", "E_B", "E_s", "q", "hs")
    DEB.state.imago$E_s[DEB.state.imago$E_s > E_sm * L.j ^ 3] <- E_sm * L.j ^ 3
    DEB.state.imago$E[DEB.state.imago$E / L.e ^ 3 > E_m] <- E_m * L.e ^ 3

    V <- c(DEB.state.embryo$L ^ 3, DEB.state.larva$L ^ 3, DEB.state.pupa$V + DEB.state.pupa$L ^ 3, DEB.state.imago$E_s * 0 + L.e ^ 3)
    V2 <- c(DEB.state.embryo$L ^ 3, DEB.state.larva$L ^ 3, DEB.state.pupa$L ^ 3, DEB.state.imago$E_s * 0 + L.e ^ 3)
    E2 <- c(DEB.state.embryo$E / DEB.state.embryo$L ^ 3, DEB.state.larva$E / DEB.state.larva$L ^ 3, DEB.state.pupa$E / DEB.state.pupa$L ^ 3, DEB.state.imago$E / L.e ^ 3)
    E <- c(DEB.state.embryo$E / DEB.state.embryo$L ^ 3, DEB.state.larva$E / DEB.state.larva$L ^ 3, DEB.state.pupa$E / (DEB.state.pupa$V + DEB.state.pupa$L ^ 3), DEB.state.imago$E / L.e ^ 3)
    #E2 <- c(DEB.state.embryo$E, DEB.state.larva$E, DEB.state.pupa$E, DEB.state.imago$E)
    H <- c(DEB.state.embryo$H, DEB.state.larva$E * 0 + tail(DEB.state.embryo$H, 1), DEB.state.pupa$H, DEB.state.imago$E * 0 + tail(DEB.state.pupa$H, 1))
    E_s <- c(DEB.state.embryo$L * 0, DEB.state.larva$E_s, DEB.state.pupa$V * 0, DEB.state.imago$E_s)
    E_s[E_s > V * E_sm] <- V[E_s > V * E_sm] * E_sm
    q <- c(DEB.state.embryo$L * 0, DEB.state.larva$L * 0, DEB.state.pupa$V * 0, DEB.state.imago$q)
    hs <- c(DEB.state.embryo$L * 0, DEB.state.larva$L * 0, DEB.state.pupa$V * 0, DEB.state.imago$hs)
    S <- hs * 0
    R <- c(DEB.state.embryo$L * 0, DEB.state.larva$E_R, DEB.state.pupa$V * 0 + tail(DEB.state.larva$E_R, 1), DEB.state.imago$E_R)
    B <- 0
    DEB.state <- as.data.frame(cbind(V, E, H, E_s, S, q, hs, R, B))
    E_Hp <- max(H)
    E_Hj <- E_Hp
    E_He <- E_Hp
  }

  # temperature correction
  Tb <- Tbf(seq(0, (nrow(DEB.state)-1)))
  Tcorr <- exp(T_A / T_REF - T_A / (273.15 + Tb)) * (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)))
  Tcorr2 <- exp(T_A2 / T_REF - T_A2 / (273.15 + Tb)) * (1 + exp(T_AL2 / T_REF - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / T_REF)) / (1 + exp(T_AL2 / (273.15 + Tb) - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / (273.15 + Tb)))
  X <- Xf(seq(0, nrow(DEB.state)-1))

  DEB.state$R[DEB.state$R < 0] <- 0
  DEB.state$B[DEB.state$B < 0] <- 0
  DEB.state$E_s[DEB.state$E_s < 0] <- 0
  V <- DEB.state$V
  if(metab_mode !=2){
    V2 <- V
  }
  E <- DEB.state$E
  E_H <- DEB.state$H
  E_s <- DEB.state$E_s
  resid <- E_s - E_sm * V2 # excess food intake to stomach capacity
  E_s[E_s > E_sm * V2] <- E_sm * V2[E_s > E_sm * V2]
  E_s[E_s < 0] <- 0
  resid[resid < 0] <- 0

  starve <- DEB.state$S
  q <- DEB.state$q
  hs <- DEB.state$hs
  p_R <- c(0, DEB.state$R[2:(length(DEB.state$R))] - DEB.state$R[1:((length(DEB.state$R)-1))])
  #p_R <- c(DEB.state$R[2:(length(DEB.state$R))] - DEB.state$R[1:((length(DEB.state$R)-1))], 0)
  if(metab_mode == 2){ # add in p_R due to maturation
    #dEH <- c(DEB.state$H[2:(length(DEB.state$H))] - DEB.state$H[1:((length(DEB.state$H)-1))], 0)
    dEH <- c(0, DEB.state$H[2:(length(DEB.state$H))] - DEB.state$H[1:((length(DEB.state$H)-1))])
    dEH[dEH < 0] <- 0
    p_R[1:t_birth] <- dEH[1:t_birth]
    p_R[(t_birth + t_pupate):(t_birth + t_pupate + t_eclose)] <- dEH[(t_birth + t_pupate):(t_birth + t_pupate + t_eclose)]
  }
  p_B <- c(0, DEB.state$B[2:(length(DEB.state$B))] - DEB.state$B[1:((length(DEB.state$B)-1))])
  #p_B <- c(DEB.state$B[2:(length(DEB.state$B))] - DEB.state$B[1:((length(DEB.state$B)-1))], 0)
  if(metab_mode == 2){
    p_R[which(hs == init[5])-1] <- 0 # fixing spike caused by transition between pupa and imago
  }
  #p_B <- p_R + p_B
  p_B[E_H < E_Hp] <- 0
  p_R[p_R < 0] <- 0
  p_B[p_B < 0] <- 0
  # if(max(E_H) > E_Hp){ # temporary workaround for weird p_G driven by low p_R at transition to maturity
  #   suppressWarnings(p_R_fix <- which(p_R == p_R[E_H > E_Hp])[1] - 1)
  #   p_R[p_R_fix] <- (p_R[p_R_fix - 1] + p_B[p_R_fix + 1]) / 2
  # }
  E_R <- p_R * 0
  E_B <- E_R
  if(metab_mode != 2){
    E_R[E_H >= E_Hp] <- cumsum(p_R[E_H >= E_Hp])
  }else{
    E_R <- cumsum(p_R)
  }
  E_B[E_H >= E_Hp] <- cumsum(p_B[E_H >= E_Hp] * kap_R)
  s_M <- rep(1, length(E_H))
  if(metab_mode == 2){
    s_M <- V ^ (1 / 3) / L.b
    s_M[E_H < E_Hb & p_A == 0] <- 1
    s_M[DEB.state$R / DEB.state$V > E_RJ] <- L.j / L.b
  }else{
    if(E_Hb != E_Hj){
      s_M[E_H > E_Hb] <- V[E_H > E_Hb] ^ (1 / 3) / L.b
      s_M[E_H > E_Hj] <- L.j / L.b
    }
  }
  e <- E / E_m # use new value of e
  L_w <- V ^ (1 / 3) / del_M * 10 # length in mm
  L_b <- V[E_H >= E_Hb][1] ^ (1 / 3)
  if(metab_mode < 2){
    L_j <- V[E_H >= E_Hj][1] ^ (1 / 3)
  }else{
    L_j <- V[DEB.state$R/DEB.state$V >= E_RJ][1] ^ (1 / 3)
  }
  V_m <- (kap * (p_Am * s_M) / p_M) ^ 3 # cm ^ 3, maximum structural volume
  L_m <- V_m ^ (1 / 3)

  # some powers
  p_M2 <- p_M * Tcorr * V + p_T * V ^ (2 / 3)
  p_M2[p_M2 < 0] <- 0
  p_J <- k_J * E_H * Tcorr2
  p_J[p_J < 0] <- 0
  p_A <- V ^ (2 / 3) * p_Am * Tcorr * s_M * f
  p_A[E_s == 0] <- 0
  p_A[p_A < 0] <- 0
  #p_A[E_s > V ^ (2 / 3) * p_Am * Tcorr * s_M * f] <- V[E_s > V ^ (2 / 3) * p_Am * Tcorr * s_M * f] ^ (2 / 3) * p_Am * s_M[E_s > V ^ (2 / 3) * p_Am * Tcorr * s_M * f] * f
  r <- v * Tcorr * s_M * (e / V ^ (1 / 3) - (1 + L_T / V ^ (1 / 3)) / L_m) / (e + g)
  p_C <- E * ((v * Tcorr * s_M) / V ^ (1 / 3) - r) * V # J / t, mobilisation rate, equation 2.12 DEB3

  if(metab_mode == 1){
    dE <- c(DEB.state$E[2:(length(DEB.state$E))] - DEB.state$E[1:((length(DEB.state$E)-1))], 0)
    p_A[E_H >= E_Hj] <- p_R[E_H >= E_Hj] + p_B[E_H >= E_Hj] + p_M2[E_H >= E_Hj] + p_J[E_H >= E_Hj] + dE[E_H >= E_Hj] * V[E_H >= E_Hj]
    p_C[E_H >= E_Hj] <- p_A[E_H >= E_Hj] - dE[E_H >= E_Hj] * V[E_H >= E_Hj]
    p_A[p_A < 0] <- 0
  }
  p_C[p_C < 0] <- 0
  p_D <- p_M2 + p_J + p_R
  p_D[E_H >= E_Hp] <- p_M2[E_H >= E_Hp] + p_J[E_H >= E_Hp] + (1 - kap_R) * p_B[E_H >= E_Hp]

  if(metab_mode == 2){
    p_G <- p_C - p_M2 - p_J
    p_G[E_H < E_He & p_A > 0] <- p_G[E_H < E_He & p_A > 0] - p_R[E_H < E_He & p_A > 0]
  }else{
    p_G <- p_C - p_M2 - p_J - p_R - p_B
  }
  if(metab_mode >= 1){
    p_G[E_H >= E_Hj] <- 0
  }
  p_G[p_G < 0] <- 0
  p_X <- max(0, f * p_Xm * s_M * V ^ (2 / 3) * (X / K) / (1 + X / K) - resid)
  if(metab_mode < 2){
    p_X[E_H < E_Hb] <- 0
  }else{
    p_X[E_s == 0] <- 0
  }

  # initialise for reproduction and starvation
  if(clutch_ab[1] > 0){
    clutchsize <- floor(clutch_ab[1] * (V ^ (1 / 3) / del_M) - clutch_ab[2])
    clutchsize[clutchsize < 0] <- 0
  }
  clutchenergy <- E_0 * clutchsize
  clutches <- rep(0, length(E_B))
  fecundity <- clutches
  firstlay <- 0
  for(i in 2:length(E_B)){
    if(clutch_ab[1] > 0){
      clutchnrg <- clutchenergy[i]
      sizeclutch <- clutchsize[i]
    }else{
      clutchnrg <- clutchenergy
      sizeclutch <- clutchsize
    }
    if(E_B[i - 1] > clutchnrg){
      E_B[i:length(E_B)] <- E_B[i:length(E_B)] - clutchnrg
      fecundity[i] <- sizeclutch
      clutches[i] <- 1
      if(firstlay == 0){
        firstlay <- i
      }
    }
  }

  # determine stages

  stage <- rep(0, length(E_H))

  # STD MODEL
  if(metab_mode == 0 & E_Hb == E_Hj){
    stage[E_H > E_Hb] <- 1
    stage[E_H > E_Hp] <- 2
    stage[firstlay:length(stage)] <- 3
  }

  # ABJ MODEL
  if(metab_mode == 0 & E_Hb != E_Hj){
    stage[E_H > E_Hb] <- 1
    stage[E_H > E_Hj] <- 2
    stage[E_H > E_Hp] <- 3
    stage[firstlay:length(stage)] <- 4
  }

  # ABP acceleration model
  if(metab_mode == 1){
    L_instar <- rep(0, stages)
    L_instar[1] <- S_instar[1] ^ 0.5 * L_b
    stage[E_H > E_Hb] <- 1
    for(j in 2:stages){
      L_instar[j] <- S_instar[j] ^ 0.5 * L_instar[j - 1]
      L_thresh <- L_instar[j]
      stage[V^(1/3) > L_thresh] <- j
    }
    stage[E_H > E_Hp] <- stages
    stage[firstlay:length(stage)] <- stages + 1
  }

  # hex model
  if(metab_mode == 2){
    L_instar <- rep(0, stages)
    L_instar[1] <- S_instar[1] ^ 0.5 * L_b
    stage[E_H > E_Hb] <- 1
    for(j in 2:stages){
      L_instar[j] <- S_instar[j] ^ 0.5 * L_instar[j - 1]
      L_thresh <- L_instar[j]
      stage[V ^ (1 / 3) > L_thresh] <- j
    }
    stage[E_R >= E_RJ] <- stages - 2
    stage[E_H >= E_He] <- stages - 1
    stage[firstlay:length(stage)] <- stages
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

  # metabolic heat production (Watts) - growth overhead plus dissipation power (maintenance, maturity maintenance,
  # maturation/repro overheads) plus assimilation overheads
  # DEBQMETW <- ((1 - kap_G) * p_G + p_D + (p_A / kap_X - p_A - p_A * mu_P * eta_PA)) / 3600 / Tcorr
  mu_O <- c(mu_X, mu_V, mu_E, mu_P) # J/mol, chemical potentials of organics
  mu_M <- c(0, 0, 0, mu_N)          # J/mol, chemical potentials of minerals C: CO2, H: H2O, O: O2, N: nitrogenous waste
  J_O <- c(JOJx, JOJv, JOJe, JOJp) # eta_O * diag(p_ADG(2,:)); # mol/d, J_X, J_V, J_E, J_P in rows, A, D, G in cols
  J_M <- c(JMCO2, JMH2O, JMO2, JMNWASTE) # - (n_M\n_O) * J_O;        # mol/d, J_C, J_H, J_O, J_N in rows, A, D, G in cols

  # compute heat production
  p_T <- sum(-1 * J_O * h_O -J_M * h_M) / 3600
  DEBQMETW <- p_T

  E_R[E_R < 0] <- 0
  E_B[E_B < 0] <- 0
  GDRYFOOD <- -1 * JOJx * w_X
  GFAECES <- JOJp * w_P
  GNWASTE <- JMNWASTE * w_N[1]
  wetgonad <- ((E_R / mu_E) * w_E) / d_Egg + ((E_B / mu_E) * w_E) / d_Egg
  wetstorage <- ((V * E / mu_E) * w_E) / d_E
  wetgut <- ((E_s / mu_E) * w_E) / fdry
  wetmass <- V * andens_deb + wetgonad + wetstorage + wetgut

  p_surv <- hs * 0
  p_surv[1] <- p_surv_init
  for(i in 2:length(p_surv)){
    p_surv[i] <- p_surv[i - 1] + p_surv[i - 1] * -1 * hs[i]
  }
  p_M <- p_M2
  deb.names <- c("stage", "V", "E", "E_H", "E_s", "E_R", "E_B", "q", "hs", "length", "wetmass", "wetgonad", "wetgut", "wetstorage", "p_surv", "fecundity", "clutches", "JMO2", "JMCO2", "JMH2O", "JMNWASTE", "O2ML", "CO2ML", "GH2OMET", "DEBQMETW", "GDRYFOOD", "GFAECES", "GNWASTE", "p_A", "p_C", "p_M", "p_G", "p_D", "p_J", "p_R", "p_B", "L_b", "L_j", "Tb")
  results.deb <- cbind(stage, V, E, E_H, E_s, E_R, E_B, q, hs, L_w, wetmass, wetgonad, wetgut, wetstorage, p_surv, fecundity, clutches, JMO2, JMCO2, JMH2O, JMNWASTE, O2ML, CO2ML, GH2OMET, DEBQMETW, GDRYFOOD, GFAECES, GNWASTE, p_A, p_C, p_M, p_G, p_D, p_J, p_R, p_B, L_b, L_j, Tb)
  names(results.deb) <- deb.names
  return(results.deb)
}
