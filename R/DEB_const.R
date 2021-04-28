#' Dynamic Energy Budget model
#'
#' Implementation of the Standard Dynamic Energy Budget model of Kooijman
#' Note that this uses the deSolve package 'ode' function with events and
#' can only handle constant food and temperature. It runs faster than the 'DEB'
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
#' @param X = 11, Food density (J/cm2)
#' @param andens_deb = 1, Animal density (g/cm3)
#' @param d_V = 0.3, Dry mass fraction of structure
#' @param d_E = 0.3, Dry mass fraction of reserve
#' @param d_Egg = 0.3, Dry mass fraction of egg
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
#' @param stages = 3, how many life stages?
#' @param stage = 0, Initial stage (0=embryo, for STD 1=juvenile, 2=mature but not yet reproducing, 3=beyond first reproduction, for ABP 1-(stages-1) = instars, stages = adult)
#' @param S_instar = rep(1.6, stages), stress at instar n: L_n^2/ L_n-1^2 (-)
#' @param clutchsize = 2, Clutch size (#), overridden by \code{clutch_ab}
#' @param clutch_ab = c(0,0), paramters for relationship between length (cm) and clutch size: clutch size = a*L_w-b, make a and b zero if fixed clutch size
#' @param minclutch = 0, Minimum clutch size if not enough in reproduction buffer for clutch size predicted by \code{clutch_ab} - if zero, will not operate
#' @param batch = 1, Invoke Pequerie et al.'s batch laying model?
#' @param lambda = 1/2
#' @param acthr = 1
#' @param X = 11
#' @param E_init = 6011.93
#' @param V_init = 3.9752^3
#' @param E_H_init = 73592
#' @param q_init = 0
#' @param hs_init = 0
#' @param p_surv_init = 1
#' @param E_s_init = 0
#' @param E_R_init = 0
#' @param E_B_init = 0
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
DEB_const<-function(
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
  clutchsize=2,
  clutch_ab=c(0.085,0.7),
  batch=1,
  lambda=1/2,
  X=10,
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
  Tb=33,
  fdry=0.3,
  L_b=0.42,
  L_j=1.376,
  S_instar=rep(1.618, stages),
  day=1,
  metab_mode=0,
  age=0){

  if (!require("deSolve", quietly = TRUE)) {
    stop("package 'deSolve' is needed. Please install it.",
         call. = FALSE)
  }

  #DEB mass balance-related calculations
  n_O <- cbind(n_X, n_V, n_E, n_P) # matrix of composition of organics, i.e. food, structure, reserve and faeces
  CHON <- c(12, 1, 16, 14)
  wO <- CHON %*% n_O
  w_V <- wO[3]
  M_V <- d_V / w_V
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

  # temperature corrections and compound parameters
  M_V <- d_V / w_V
  p_MT <- p_M * Tcorr
  k_M <- p_MT / E_G
  k_JT <- k_J * Tcorr2
  vT <- v * Tcorr
  p_AmT <- p_MT * z / kap
  p_XmT <- p_Xm * Tcorr
  h_aT <- h_a * Tcorr
  E_m <- p_AmT / vT
  g <- E_G / (kap * E_m) # energy investment ratio
  e <- E_init / E_m # scaled reserve density
  V_m <- (kap * p_AmT / p_MT) ^ 3 # maximum structural volume
  L_T <- p_T / p_MT # heating length
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

  # function for solver (running for one time step)
  dget_DEB <- function(t, y, indata){
    with(as.list(c(indata, y)), {

      # unpack variables
      V <- y[1]# cm^3, structural volume
      E <- y[2]# J/cm3, reserve density
      H <- y[3]# J, maturity
      E_s <- y[4]# J, stomach energy
      S <- y[5]# J, starvation energy
      q <- y[6]# -, aging acceleration
      hs <- y[7]# -, hazard rate
      R <- y[8]# J, reproduction buffer energy
      B <- y[9]# J, egg batch energy

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
      r <- (v * s_M) * (e / L - (1 + L_T / L) / L_m) / (e + g) # specific growth rate
      p_C <- E * V * ((v * s_M) / L - r) # J / t, mobilisation rate, equation 2.12 DEB3
      if(metab_mode == 1 & H >= E_Hj){
        r <- min(0, r) # no growth in abp after puberty, but could still be negative because starving
        p_C <- E * V * (v * s_M) / L
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
          if(B < dS){ # batch buffer has run out so draw from structure
            dV <- V * r
            dS <- 0
          }
        }else{
          dS <- 0
        }

        # assimilation
        p_A <- (p_Am * s_M) * f * L ^ 2

        # reserve
        if(E_s > p_A){
          dE <- p_A / L ^ 3 - (E * (v * s_M)) / L
        }else{
          dE <- E_s / L ^ 3 - (E * (v * s_M)) / L
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
        acthr <- 1
        # feeding
        if(acthr > 0){
          # Regulates X dynamics
          p_X <- f * (p_Xm * s_M) * ((X / K) / (1 + X / K)) * V ^ (2 / 3)
        }else{
          p_X <- 0
        }
        dEs <- p_X - ((p_Am * s_M) / kap_X) * V ^ (2 / 3)

        if(metab_mode == 1 & H >= E_Hj){
          r <- 0 # no growth in abp after puberty - not setting this to zero messes up ageing calculation
        }

        # ageing (equation 6.2 in Kooijman 2010 (DEB3)
        dq <- (q * (V / V_m) * s_G + h_a) * e * (((v * s_M) / L) - r) - r * q # ageing acceleration
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
            batchprep <- (kap_R / lambda) * ((1 - kap) * (E_m * ((v * s_M) * V ^ (2 / 3) + k_M * V) / (1 + (1 / g))) - p_J)
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
        p_R <- p_R - p_B # take finalised value of p_B from p_R

        # draw from reproduction and then batch buffers under starvation
        if(dS > 0 & p_R > dS){
          p_R <- p_R - dS
          dS <- 0
        }
        if(dS > 0 & p_B > dS){
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
  metamorphosis <- function(t, y, pars) {
    if(y[3] > E_Hj) {
      y[3] <- 0
    }
    return(y)
  }

  # parameters
  indata <- list(k_J = k_JT,
                 p_Am = p_AmT,
                 k_M = k_M,
                 p_M = p_MT,
                 p_Xm = p_XmT,
                 v = vT,
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
                 X = X,
                 K = K,
                 E_Hp = E_Hp,
                 E_Hb = E_Hb,
                 E_Hj = E_Hj,
                 s_G = s_G,
                 h_a = h_aT,
                 batch = batch,
                 kap_R = kap_R,
                 lambda = lambda,
                 breeding = breeding,
                 kap_X = kap_X,
                 f = f,
                 E_sm = E_sm,
                 L_b = L_b,
                 L_j = L_j,
                 metab_mode = metab_mode)

  times <- seq(0, (1 / step) * ndays)

  # initial conditions for solver
  init <- c(V_init, E_init, E_H_init + 1e-10, E_s_init + 1e-10, hs_init + 1e-10, q_init + 1e-10, hs_init + 1e-10, E_R_init + 1e-10, E_B_init + 1e-10)

  #"egg", "hatchling", "puberty", "adult"

  if(E_H_init < E_Hb){
    # to birth
    DEB.state.birth <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_DEB, parms = indata, method = "lsodes", events = list(func = eventfun, root = TRUE, terminalroot = 3), rootfunc = birth))[, 2:10]
    colnames(DEB.state.birth) <- c("V", "E", "H", "E_s", "S", "q", "hs", "R", "B")
    L.b <- max(DEB.state.birth$V ^ (1 / 3))
    t_birth <- which(DEB.state.birth$V ^ (1 / 3) == L.b)[1]
    DEB.state.birth <- DEB.state.birth[1:t_birth, ]
    init <- as.numeric(DEB.state.birth[nrow(DEB.state.birth), ])
    if(nrow(DEB.state.birth) > 1){
      DEB.state.birth <- head(DEB.state.birth, -1)
    }
    times <- seq(0, (1 / step) * ndays - t_birth)

    # parameters
    indata <- list(k_J = k_JT, p_Am = p_AmT, k_M = k_M, p_M = p_MT,
                   p_Xm = p_XmT, v = vT, E_m = E_m, L_m = L_m, L_T = L_T,
                   kap = kap, g = g, M_V = M_V, mu_E = mu_E,
                   mu_V = mu_V, d_V = d_V, w_V = w_V,
                   X = X, K = K, E_Hp = E_Hp, E_Hb = E_Hb, E_Hj = E_Hj, s_G = s_G, h_a = h_aT,
                   batch = batch, kap_R = kap_R, lambda = lambda,
                   breeding = breeding, kap_X = kap_X, f = f, E_sm = E_sm, L_b = L.b,
                   L_j = L_j, metab_mode = metab_mode)
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
      indata <- list(k_J = k_JT, p_Am = p_AmT, k_M = k_M, p_M = p_MT,
                     p_Xm = p_XmT, v = vT, E_m = E_m, L_m = L_m, L_T = L_T,
                     kap = kap, g = g, M_V = M_V, mu_E = mu_E,
                     mu_V = mu_V, d_V = d_V, w_V = w_V,
                     X = X, K = K, E_Hp = E_Hp, E_Hb = E_Hb, E_Hj = E_Hj, s_G = s_G, h_a = h_aT,
                     batch = batch, kap_R = kap_R, lambda = lambda,
                     breeding = breeding, kap_X = kap_X, f = f, E_sm = E_sm, L_b = L.b,
                     L_j = L.j, metab_mode = metab_mode)
      if(E_H_init < E_Hb){
        times <- seq(0, (1 / step) * ndays - t_birth - t_meta)
      }else{
        times <- seq(0, (1 / step) * ndays - t_meta)
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
        times <- seq(0, (1 / step) * ndays - t_birth - t_meta - t_puberty)
      }else{
        times <- seq(0, (1 / step) * ndays - t_meta - t_puberty)
      }
    }else{
      if(E_H_init < E_Hb){
        times <- seq(0, (1 / step) * ndays - t_birth - t_puberty)
      }else{
        times <- seq(0, (1 / step) * ndays - t_puberty)
      }
    }
  }
  # to end of time sequence

  DEB.state.end <- as.data.frame(deSolve::ode(y = init, times = times, func = dget_DEB, parms = indata, method = "lsodes"))[, 2:10]
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

  V <- DEB.state$V
  E <- DEB.state$E
  E_H <- DEB.state$H
  E_s <- DEB.state$E_s
  resid <- E_s - E_sm * V # excess food intake to stomach capacity
  E_s[E_s > E_sm * V] <- E_sm * V[E_s > E_sm * V]

  starve <- DEB.state$S
  q <- DEB.state$q
  hs <- DEB.state$hs
  p_R <- c(DEB.state$R[2:(length(DEB.state$R))] - DEB.state$R[1:((length(DEB.state$R)-1))], 0)
  p_B <- c(DEB.state$B[2:(length(DEB.state$B))] - DEB.state$B[1:((length(DEB.state$B)-1))], 0)
  p_B <- p_R + p_B
  p_B[E_H < E_Hp] <- 0
  p_R[p_R < 0] <- 0
  if(max(E_H) > E_Hp){ # temporary workaround for weird p_G driven by low p_R at transition to maturity
    suppressWarnings(p_R_fix <- which(p_R == p_R[E_H > E_Hp])[1] - 1)
    p_R[p_R_fix] <- (p_R[p_R_fix - 1] + p_B[p_R_fix + 1]) / 2
  }
  E_R <- p_R * 0
  E_B <- E_R
  E_R[E_H >= E_Hp] <- cumsum(p_R[E_H >= E_Hp])
  E_B[E_H >= E_Hp] <- cumsum(p_B[E_H >= E_Hp] * kap_R)
  s_M <- rep(1, length(E_H))
  if(E_Hb != E_Hj){
    s_M[E_H > E_Hb] <- V[E_H > E_Hb] ^ (1 / 3) / L.b
    s_M[E_H > E_Hj] <- L.j / L.b
  }
  e <- E / E_m # use new value of e
  L_w <- V ^ (1 / 3) / del_M * 10 # length in mm
  L_b <- V[E_H >= E_Hb][1] ^ (1 / 3)
  L_j <- V[E_H >= E_Hj][1] ^ (1 / 3)
  V_m <- (kap * (p_AmT * s_M) / p_MT) ^ 3 # cm ^ 3, maximum structural volume
  L_m <- V_m ^ (1 / 3)

  # some powers
  p_M2 <- p_MT * V + p_T * V ^ (2 / 3)
  p_J <- k_JT * E_H - starve
  p_A <- E_s
  p_A[E_s > V ^ (2 / 3) * p_AmT * s_M * f] <- V[E_s > V ^ (2 / 3) * p_AmT * s_M * f] ^ (2 / 3) * p_AmT * s_M[E_s > V ^ (2 / 3) * p_AmT * s_M * f] * f
  r <- vT * (e / V ^ (1 / 3) - (1 + L_T / V ^ (1 / 3)) / L_m) / (e + g)
  p_C <- E * ((vT * s_M) / V ^ (1 / 3) - r) * V # J / t, mobilisation rate, equation 2.12 DEB3
  if(metab_mode == 1){
    dE <- c(DEB.state$E[2:(length(DEB.state$E))] - DEB.state$E[1:((length(DEB.state$E)-1))], 0)
    p_A[E_H >= E_Hj] <- p_R[E_H >= E_Hj] + p_B[E_H >= E_Hj] + p_M2[E_H >= E_Hj] + p_J[E_H >= E_Hj] + dE[E_H >= E_Hj] * V[E_H >= E_Hj]
    p_C[E_H >= E_Hj] <- p_A[E_H >= E_Hj] - dE[E_H >= E_Hj] * V[E_H >= E_Hj]
  }

  p_D <- p_M2 + p_J + p_R
  p_D[E_H >= E_Hp] <- p_M2[E_H >= E_Hp] + p_J[E_H >= E_Hp] + (1 - kap_R) * p_B[E_H >= E_Hp]

  p_G <- p_C - p_M2 - p_J - p_R - p_B

  if(metab_mode == 1){
    p_G[E_H >= E_Hj] <- 0
  }

  p_X <- f * p_XmT * s_M * V ^ (2 / 3) * (X / K) / (1 + X / K) - resid
  p_X[E_H < E_Hb] <- 0

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
  # maturation/repro overheads) plus assimilation overheads - correct to 20 degrees so it can be temperature corrected
  # in MET.f for the new guessed Tb
  # DEBQMETW <- ((1 - kap_G) * p_G + p_D + (p_A / kap_X - p_A - p_A * mu_P * eta_PA)) / 3600 / Tcorr
  mu_O <- c(mu_X, mu_V, mu_E, mu_P) # J/mol, chemical potentials of organics
  mu_M <- c(0, 0, 0, mu_N)          # J/mol, chemical potentials of minerals C: CO2, H: H2O, O: O2, N: nitrogenous waste
  J_O <- c(JOJx, JOJv, JOJe, JOJp) # eta_O * diag(p_ADG(2,:)); # mol/d, J_X, J_V, J_E, J_P in rows, A, D, G in cols
  J_M <- c(JMCO2, JMH2O, JMO2, JMNWASTE) # - (n_M\n_O) * J_O;        # mol/d, J_C, J_H, J_O, J_N in rows, A, D, G in cols
  p_T <- sum(-J_O * mu_O -J_M * mu_M) / 3600 / Tcorr # W
  DEBQMETW <- p_T

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
  deb.names <- c("stage", "V", "E", "E_H", "E_s", "E_R", "E_B", "q", "hs", "length", "wetmass", "wetgonad", "wetgut", "wetstorage", "p_surv", "fecundity", "clutches", "JMO2", "JMCO2", "JMH2O", "JMNWASTE", "O2ML", "CO2ML", "GH2OMET", "DEBQMETW", "GDRYFOOD", "GFAECES", "GNWASTE", "p_A", "p_C", "p_M", "p_G", "p_D", "p_J", "p_R", "p_B", "L_b", "L_j")
  results.deb <- cbind(stage, V, E, E_H, E_s, E_R, E_B, q, hs, L_w, wetmass, wetgonad, wetgut, wetstorage, p_surv, fecundity, clutches, JMO2, JMCO2, JMH2O, JMNWASTE, O2ML, CO2ML, GH2OMET, DEBQMETW, GDRYFOOD, GFAECES, GNWASTE, p_A, p_C, p_M, p_G, p_D, p_J, p_R, p_B, L_b, L_j)
  names(results.deb) <- deb.names
  return(results.deb)
}
