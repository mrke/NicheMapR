#' Dynamic Energy Budget model
#'
#' Implementation of the Standarad Dynamic Energy Budget model of Kooijman
#' Note that this uses the deSolve package 'ode' function. The older version
#' that uses Euler integration is now called DEB_euler (and is faster and
#' may be preferable in some cases, though accuracy of the latter will depend
#' on the step size chosen)
#' Michael Kearney Dec 2015, updated to include ODE solver Feb 2019
#' @param step = 1/24, step size (days)
#' @param z = 7.997, Zoom factor (cm)
#' @param del_M =  0.242, Shape coefficient (-)
#' @param F_m = 13290*step, Surface area-specific maximum feeding rate J/cm2/h
#' @param kap_X = 0.85, Digestive efficiency (decimal \%)
#' @param v = 0.065*step, Energy conductance (cm/h)
#' @param kap = 0.886, fraction of mobilised reserve allocated to soma
#' @param p_M = 32*step, Volume-specific somatic maintenance (J/cm3/h)
#' @param E_G = 7767, Cost of structure (J/cm3)
#' @param kap_R = 0.95, Fraction of reproduction energy fixed in eggs
#' @param k_J = 0.002*step, Maturity maintenance rate coefficient (1/h)
#' @param E_Hb = 7.359e+04, Maturity at birth (J)
#' @param E_Hj = E_Hb, Maturity at metamorphosis (J)
#' @param E_Hp = 1.865e+05, Maturity at puberty
#' @param h_a = 2.16e-11*(step^2), Weibull ageing acceleration (1/h2)
#' @param s_G = 0.01, Gompertz stress coefficient (-)
#' @param E_0 = 1.04e+06, Energy content of the egg (derived from core parameters) (J)
#' @param T_REF = 20, Reference temperature for rate correction (deg C)
#' @param T_A = 8085 Arhhenius temperature
#' @param T_AL = 18721, Arrhenius temperature for decrease below lower boundary of tolerance range \code{T_L}
#' @param T_AH = 90000, Arrhenius temperature for decrease above upper boundary of tolerance range \code{T_H}
#' @param T_L = 288, Lower boundary (K) of temperature tolerance range for Arrhenius thermal response
#' @param T_H = 315, Upper boundary (K) of temperature tolerance range for Arrhenius thermal response
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
#' @param viviparous = 0, Viviparous reproduction? 1=yes, 0=no (if yes, animal will be held in adult-sided female's body for duration of development and will experience her body temperature
#' @param minclutch = 0, Minimum clutch size if not enough in reproduction buffer for clutch size predicted by \code{clutch_ab} - if zero, will not operate
#' @param batch = 1, Invoke Pequerie et al.'s batch laying model?
#' @param lambda = 1/2
#' @param VTMIN = 26, Voluntary thermal maximum, degrees C, controls whether repro event can occur at a given time
#' @param VTMAX = 39, Voluntary thermal maximum, degrees C, controls whether repro event can occur at a given time
#' @param arrhenius = matrix(data = matrix(data = c(rep(T_A,8),rep(T_AL,8),rep(T_AH,8),rep(T_L,8),rep(T_H,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5), Stage-specific 5-parameter Arrhenius thermal response for DEB model (T_A,T_AL,T_AH,T_L,T_H)
#' @param acthr = 1
#' @param X = 11
#' @param E_pres = 6011.93
#' @param V_pres = 3.9752^3
#' @param E_H_pres = 73592
#' @param q_pres =0
#' @param hs_pres =0
#' @param surviv_pres = 1
#' @param Es_pres = 0
#' @param cumrepro = 0
#' @param cumbatch = 0
#' @param stage = 1
#' @param breeding = 0
#' @param pregnant = 0
#' @param Tb = 33
#' @return E_pres
#' @return V_pres
#' @return E_H_pres
#' @return q_pres
#' @return hs_pres
#' @return surviv_pres
#' @return Es_pres
#' @return cumrepro
#' @return cumbatch
#' @return O2FLUX
#' @return CO2FLUX
#' @return MLO2
#' @return GH2OMET
#' @return DEBQMET
#' @return DRYFOOD,
#' @return FAECES
#' @return NWASTE
#' @return wetgonad
#' @return wetstorage
#' @return wetfood
#' @return wetmass
#' @return gutfreemass
#' @return gutfull
#' @return fecundity
#' @return clutches
#' @examples
#' # simulate growth and reproduction at different constant body temperatures at constant food for a lizard (Eulamprus quoyii - default parameter values, starting as a hatchling)
#'
#' n<-3000 # time steps
#' step<-1 # step size (days)
#'
#' Tbs=seq(25,35,2.5) # sequence of body temperatures to use
#'
#' for(j in 1:length(Tbs)){
#' debout<-matrix(data = 0, nrow = n, ncol=27)
#' deb.names<-c("E_pres","V_pres","E_H_pres","q_pres","hs_pres","surviv_pres","Es_pres","cumrepro","cumbatch","O2FLUX","CO2FLUX","MLO2","GH2OMET","DEBQMET","DRYFOOD","FAECES","NWASTE","wetgonad","wetstorage","wetfood","wetmass","gutfreemass","gutfull","fecundity","clutches","potfreemass")
#' colnames(debout)<-deb.names
#'
#' # initialise
#' debout[1,]<-DEB(Tb = Tbs[j], step = step)
#'
#' for(i in 2:n){
#' debout[i,]<-DEB(Tb = Tbs[j], breeding = 1, step = step,E_pres=debout[(i-1),1],V_pres=debout[(i-1),2],E_H_pres=debout[(i-1),3],q_pres=debout[(i-1),4],hs_pres=debout[(i-1),5],surviv_pres=debout[(i-1),6],Es_pres=debout[(i-1),7],cumrepro=debout[(i-1),8],cumbatch=debout[(i-1),9])
#' }
#'
#' if(j==1){
#'   plot((seq(1,n)/365),debout[,21],ylim=c(100,1500),type='l',xlab='years',ylab='wet mass, g', col=j)
#' }else{
#'   points((seq(1,n)/365),debout[,21],ylim=c(100,1500),type='l',xlab='years',ylab='wet mass, g',col=j)
#' }
#'
#' } #end loop through body temperatures
#'
#'
#'
#' @export
DEB<-function(
  step=1/24,
  z=7.997,
  del_M=0.242,
  F_m=13290*step,
  kap_X=0.85,
  v=0.065*step,
  kap=0.886,
  p_M=32*step,
  E_G=7767,
  kap_R=0.95,
  k_J=0.002*step,
  E_Hb=7.359e+04,
  E_Hj=E_Hb,
  E_Hp=1.865e+05,
  h_a=2.16e-11/(step^2),
  s_G=0.01,
  T_REF=20,
  T_A=8085,
  T_AL=18721,
  T_AH=9.0E+04,
  T_L=288,
  T_H=315,
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
  kap_X_P=0.1,
  n_X=c(1,1.8,0.5,.15),
  n_E=c(1,1.8,0.5,0.15),
  n_V=c(1,1.8,0.5,.15),
  n_P=c(1,1.8,0.5,.15),
  n_M_nitro=c(1,4/5,3/5,4/5),
  clutchsize=2,
  clutch_ab=c(0.085,0.7),
  viviparous=0,
  minclutch=0,
  batch=1,
  lambda=1/2,
  VTMIN=26,
  VTMAX=39,
  ma=1e-4,
  mi=0,
  mh=0.5,
  arrhenius=matrix(data = matrix(data = c(rep(T_A,8),rep(T_AL,8),rep(T_AH,8),rep(T_L,8),rep(T_H,8)),nrow = 8, ncol = 5), nrow = 8, ncol = 5),
  acthr=1,
  X=10,
  E_pres=6011.93,
  V_pres=3.9752^3,
  E_H_pres=73592,
  q_pres=0,
  hs_pres=0,
  surviv_pres=1,
  Es_pres=0,
  cumrepro=0,
  cumbatch=0,
  stages=3,
  stage=1,
  breeding=0,
  pregnant=0,
  Tb=33,
  fdry=0.3,
  L_b=0.42,
  L_j=1.376,
  S_instar=rep(1.6, stages),
  spawnday=1,
  day=1,
  metab_mode=0){

  if (!require("deSolve", quietly = TRUE)) {
    stop("package 'deSolve' is needed. Please install it.",
         call. = FALSE)
  }

  # initialise for reproduction and starvation
  if(clutch_ab[1] > 0){
    clutchsize <- floor(clutch_ab[1] * (V_pres ^ (1 / 3) / del_M) - clutch_ab[2])
  }
  orig_clutchsize <- clutchsize
  fecundity <- 0
  clutches <- 0
  clutchenergy <- E_0 * clutchsize
  starve <- 0

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
  n_M_inv <- matrix(c(1, 0, -1, 0, 0, 1/2, -1/4, 0, 0, 0, 1/2, 0, n_M_nitro_inv), nrow=4)
  JM_JO <- -1 * n_M_inv %*% n_O
  etaO <- matrix(c(y_XE / mu_E * -1, 0, 1 / mu_E, y_PE / mu_E, 0, 0, -1 / mu_E, 0, 0, y_VE / mu_E, -1 / mu_E, 0), nrow = 4)
  w_N <- CHON %*% n_M_nitro

  # Arrhenius temperature correction factor
  Tcorr <- exp(T_A * (1 / (273.15 + T_REF) - 1 / (273.15 + Tb))) / (1 + exp(T_AL * (1 / (273.15 + Tb) - 1 / T_L)) + exp(T_AH * (1 / T_H - 1 / (273.15 + Tb))))
  Tcorr <- exp(T_A / (273.15 + T_REF) - T_A / (273.15 + Tb)) * (1 + exp(T_AL / (273.15 + T_REF) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + T_REF))) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)))

  # metabolic acceleration if present
  s_M <- 1 # -, multiplication factor for v and {p_Am} under metabolic acceleration
  if(E_Hj != E_Hb){
    if(E_H_pres < E_Hb){
      s_M <- 1 # -, multiplication factor for v and {p_Am}
    }else{
      if(E_H_pres < E_Hj){
        s_M <- V_pres ^ (1 / 3) / L_b
      }else{
        s_M <- L_j / L_b
      }
    }
  }

  # temperature corrections and compound parameters
  p_MT <- p_M * Tcorr
  k_M <- p_MT / E_G
  k_JT <- k_J * Tcorr
  p_AmT <- p_MT * z / kap * s_M
  vT <- v * Tcorr * s_M
  E_m <- p_AmT / vT
  F_mT <- F_m * Tcorr * s_M
  g <- E_G / (kap * E_m) # energy investment ratio
  e <- E_pres / E_m # scaled reserve density
  V_m <- (kap * p_AmT / p_MT) ^ (3) # maximum structural volume
  h_aT <- h_a * Tcorr
  L_T <- 0 # heating length - not used for now
  L_pres <- V_pres ^ (1 / 3)
  L_m <- V_m ^ (1 / 3)
  scaled_l <- L_pres / L_m
  kappa_G <- (d_V * mu_V) / (w_V * E_G)
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
  init <- c(V_pres, E_pres, E_H_pres, Es_pres, starve, q_pres, hs_pres, cumrepro, cumbatch)

  # parameters
  indata <- list(k_J = k_JT, p_Am = p_AmT, k_M = k_M, p_M = p_MT,
                 F_m = F_mT, v = vT, E_m = E_m, L_m = L_m, L_T = L_T,
                 kap = kap, g = g, M_V = M_V, mu_E = mu_E,
                 mu_V = mu_V, d_V = d_V, w_V = w_V, acthr = acthr,
                 X = X, K = K, E_Hp = E_Hp, E_Hb = E_Hb, E_Hj = E_Hj, s_G = s_G, h_a = h_aT,
                 pregnant = pregnant, batch = batch, kap_R = kap_R, lambda = lambda,
                 breeding = breeding, kap_X = kap_X, f = f, E_sm = E_sm, s_M = s_M,
                 L_j = L_j, metab_mode = metab_mode)

  # function for solver (running for one time step)
  dget_DEB <- function(t, y, indata){
    with(as.list(c(indata, y)), {

      # unpack variables
      V <- y[1]# cm^3, structural volume
      E <- y[2]# J/cm3, reserve density
      H <- y[3]# J, maturity
      Es <- y[4]# J, gut energy
      S <- y[5]# J, starvation energy
      q <- y[6]# -, aging acceleration
      hs <- y[7]# -, hazard rate
      R <- y[8]# J, reproduction buffer energy
      B <- y[9]# J, egg batch energy

      L <- V ^ (1/3) # cm, structural length
      V_m <- L_m ^ 3 # cm ^ 3, maximum structural volume
      e <- E / E_m  # -, scaled reserve density
      r <- v * (e / L - (1 + L_T / L) / L_m) / (e + g) # specific growth rate
      dV <- V * r                 # cm^3 / t, change in structure
      p_C <- (E_m * (v / L + k_M * (1 + L_T / L)) * (e * g) / (e + g)) * V # J / t, mobilisation rate, equation 2.20 DEB3

      if(H < E_Hb){ # embryo
        # structure
        dE <- (- 1 *  E * v) / L
        dH <- (1 - kap) * p_C - k_J * H # J/d, change in maturity
        # no aging or stomach in embryo
        dS <- 0
        dEs <- 0
        dq <- 0
        dhs <- 0
        dR <- 0
        dB <- 0
      }else{ # post-embryo
        # structure
        p_J <- k_J * H + S # adding starvation costs to P_J so maturation time (immature) or reproduction (mature) can be sacrificed to pay somatic maintenance
        if(dV < 0){
          dS <- dV * -1 * mu_V * d_V / w_V # J / t, starvation energy to be subtracted from reproduction buffer if necessary
          dV <- 0
        }else{
          dS <- 0
        }
        if(metab_mode == 0){
          p_R <- (1 - kap) * p_C - p_J
        }
        if(metab_mode == 1){
          if(H > E_Hj){
            p_R <- p_C - p_M * V - p_J # no kappa-rule - absolute reserve amount never reaches steady state so reproduction gets all of what would otherwise have gone to growth
            dV <- 0
          }else{
            p_R <- (1 - kap) * p_C - p_J
          }
        }
        if(p_R < 0 & S > 0){
          dV <- abs(p_R) * w_V / (mu_V * d_V) * -1 # subtract from structure since not enough flow to reproduction to pay for somatic maintenance
          p_R <- 0
        }
        p_A <- p_Am * f * L ^ 2
        if(metab_mode == 1){
          if(p_A > p_C & E == E_m){
            p_A <- p_C
          }
        }
        if(Es > 0){
          dE <- p_A / L ^ 3 - (E * v) / L
        }else{
          dE <- (-E * v) / L
        }
        if(acthr > 0){
          # Regulates X dynamics
          J_X <- F_m * ((X / K) / (1 + X / K)) * V ^ (2 / 3)
          dEs <- J_X * f - (p_Am / kap_X) * V ^ (2 / 3)
        }else{
          dEs <- -1 * J_X * (p_Am / kap_X) * V ^(2 / 3)
        }
        if(H < E_Hp){
          dH <- (1 - kap) * p_C - p_J
        }else{
          dH <- 0
        }

        dq <- (q * (V / V_m) * s_G + h_a) * e * ((v / L) - r) - r * q # aging acceleration
        dhs <- q - r * hs # hazard

        # reproduction
        if((H <= E_Hp) | (pregnant==1)){
          p_B <- 0
        }else{
          if(batch == 1){
            if(metab_mode == 0){
              batchprep <- (kap_R / lambda) * ((1 - kap) * (E_m * (v * V ^ (2 / 3) + k_M * V) / (1 + (1 / g))) - p_J)
            }else{
              batchprep <- (kap_R / lambda) * ((E_m * (v * V ^ (2 / 3) + k_M * V) / (1 + (1 / g))) - p_J - p_M)
            }
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

        #accumulate energy/matter in reproduction and batch buffers
        if(H > E_Hp){
          dR <- p_R * kap_R - p_B
          dB <- p_B
        }else{
          dR <- 0
          dB <- 0
        }
      }

      y = list(c(dV, dE, dH, dEs, dS, dq, dhs, dR, dB))
    })
  }

  DEB.state <- as.data.frame(ode(y = init, times = c(0, 1), func = dget_DEB, parms = indata, method = "ode45"))[2,2:10]
  colnames(DEB.state) <- c("V", "E", "H", "Es", "S", "q", "hs", "R", "B")
  V <- max(DEB.state$V, 0)
  E <- max(DEB.state$E, 0)
  E_H <- max(DEB.state$H, 0)
  Es <- max(DEB.state$Es, 0)
  starve <- max(DEB.state$S, 0)
  q <- max(DEB.state$q, 0)
  hs <- max(DEB.state$hs, 0)
  cumrepro <- max(DEB.state$R, 0)
  cumbatch <- max(DEB.state$B, 0)

  L_w = V ^ (1 / 3) / del_M * 10 # length in mm
  if(Es > E_sm * V){
    Es <- E_sm * V
  }
  gutfull <- Es / (E_sm * V)
  if(gutfull > 1){
    gutfull <- 1
  }

  # some powers
  p_M2 <- p_MT * V
  p_J <- k_JT * E_H - starve
  if(Es > 0){
    p_A = V ^ (2 / 3) * p_AmT * f
  }else{
    p_A = 0
  }
  p_C <- (E_m * (vT / V^(1/3) + k_M * (1 + L_T / V^(1/3))) * (e * g) / (e + g)) * V #equation 2.20 DEB3

  if(metab_mode == 0){
    p_R <- (1 - kap) * p_C - p_J
  }
  if(metab_mode == 1){
    if(E_H_pres > E_Hj){
      p_R <- p_C - p_J - p_M2
    }else{
      p_R <- (1 - kap) * p_C - p_J
    }
  }

  if(metab_mode == 1){
    #if(E_H_pres > E_Hj){
      if(p_A > p_C & E == E_m){
       p_A <- p_C
      }
    #}
  }

  p_X <- p_A / kap_X #J food eaten per hour
  if(E_H_pres >= E_Hp){
    p_D = p_M2 + p_J + (1 - kap_R) * p_R
  }else{
    p_D = p_M2 + p_J + p_R
  }
  p_G = p_C - p_M2 - p_J - p_R

  testclutch <- floor((cumrepro + cumbatch) / E_0)
  # FOR VARIABLE CLUTCH SIZE FROM REPRO AND BATCH BUFFERS
  if(minclutch > 0 & floor(cumrepro + cumbatch) / E_0 > minclutch){
    if(testclutch <= orig_clutchsize){# ! MAKE SMALLEST CLUTCH ALLOWABLE FOR THIS REPRO EVENT
      clutchsize <- minclutch
      clutchenergy <- clutchsize * E_0
    }
  }

  # determine stages

  # STD MODEL
  if(metab_mode == 0){
    if(stage == 2){
      if(cumbatch < 0.1 * clutchenergy){
        stage <- 3
      }
    }
    if(E_H <= E_Hb){
      stage <- 0
    }else{
      if(E_H < E_Hj){
        stage <- 1
      }else{
        if(E_H < E_Hp){
          stage <- 2
        }else{
          stage <- 3
        }
      }
    }
    if(cumbatch > 0){
      if(E_H > E_Hp){
        stage <- 4
      }else{
        stage <- stage
      }
    }
  }

  if(metab_mode == 1){
    L_instar <- rep(0, stages)
    L_instar[1] <- S_instar[1] ^ 0.5 * L_b
    for(j in 2:stages){
      L_instar[j] <- S_instar[j] ^ 0.5 * L_instar[j - 1]
    }
    if(stage > 0 & stage < stages - 1){
      L_thresh <- L_instar[stage]
    }
    if(stage == 0){
      if(E_H > E_Hb){
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


  if((cumbatch>clutchenergy) | (pregnant==1)){
    if(viviparous == 1){
      if((pregnant == 0) & (breeding == 1)){
        v_baby <- v_init_baby
        e_baby <- e_init_baby
        EH_baby <- 0
        pregnant <- 1
        testclutch <- floor(cumbatch / E_0)
        if(testclutch > clutchsize){
          clutchsize <- testclutch
          clutchenergy <- E_0*clutchsize
        }
        # for variable clutch size from repro and batch buffers
        if(cumbatch<clutchenergy){
          # needs to draw from repro buffer - temporarily store current repro as cumrepro_temp,
          # { remove what is needed from the repro buffer and add it to the batch buffer
          cumrepro_temp <- cumrepro
          cumrepro <- cumrepro+cumbatch-clutchenergy
          cumbatch <- cumbatch+cumrepro_temp-cumrepro
        }
      }
      if(hour==1){
        v_baby <- v_baby_init
        e_baby <- e_baby_init
        EH_baby <- EH_baby_init
      }
      #if(pregnant==1){
      #call deb_baby
      #}
      if(EH_baby>E_Hb){
        if((Tb < VTMIN)  |  (Tb > VTMAX)){
          #goto 898
        }
        cumbatch(hour) <- cumbatch(hour) - clutchenergy
        repro(hour) <- 1
        pregnant <- 0
        v_baby <- v_init_baby
        e_baby <- e_init_baby
        EH_baby <- 0
        newclutch <- clutchsize
        fecundity <- clutchsize
        clutches <- 1
        pregnant <- 0
      }
    }else{
      #not viviparous, so lay the eggs at next period of activity
      if((Tb >= VTMIN)  |  (Tb <= VTMAX)){
        #    change below to active or not active rather than depth-based, in case of fossorial
        #if((Tb < VTMIN)  |  (Tb > VTMAX)){
        #}
        if(day == spawnday & spawnday != 0){
          testclutch <- floor(cumbatch / E_0)
          if(testclutch > clutchsize){
            clutchsize <- testclutch
            clutchenergy <- clutchsize * E_0
          }
          if(spawnday > 0){
            clutchsize <- testclutch
            clutchenergy <- clutchsize * E_0
          }
          cumbatch <- cumbatch - clutchenergy
          repro <- 1
          fecundity <- clutchsize
          clutches <- 1
        }else{
          if(spawnday == 0){
            testclutch <- floor(cumbatch / E_0)
            if(testclutch > clutchsize){
              clutchsize <- testclutch
              clutchenergy <- clutchsize * E_0
            }
            cumbatch <- cumbatch - clutchenergy
            repro <- 1
            fecundity <- clutchsize
            clutches <- 1
          }
        }
      }
    }
  }

  # feeding (gut) model (for output of food in - Es computed internally in dget_DEB above)
  if(E_H_pres > E_Hb){
    if(acthr > 0){
      # Regulates X dynamics
      J_X <- F_mT * ((X / K) / (1 + X / K)) * V_pres ^ (2 / 3)
      dEsdt <- J_X * f - (p_AmT / kap_X) * V_pres ^ (2 / 3)
    }else{
      dEsdt <- -1 * (p_AmT / kap_X) * V_pres ^ (2 / 3)
    }
  }else{
    dEsdt <- -1 * (p_AmT / kap_X) * V_pres ^ (2 / 3)
  }

  #mass balance
  JOJx <- p_A * etaO[1,1] + p_D * etaO[1,2] + p_G * etaO[1,3]
  JOJv <- p_A * etaO[2,1] + p_D * etaO[2,2] + p_G * etaO[2,3]
  JOJe <- p_A * etaO[3,1] + p_D * etaO[3,2] + p_G * etaO[3,3]
  JOJp <- p_A * etaO[4,1] + p_D * etaO[4,2] + p_G * etaO[4,3]

  JOJx_GM <- p_D * etaO[1,2] + p_G * etaO[1,3]
  JOJv_GM <- p_D * etaO[2,2] + p_G * etaO[2,3]
  JOJe_GM <- p_D * etaO[3,2] + p_G * etaO[3,3]
  JOJp_GM <- p_D * etaO[4,2] + p_G * etaO[4,3]

  JMCO2 <- JOJx * JM_JO[1,1] + JOJv * JM_JO[1,2] + JOJe * JM_JO[1,3] + JOJp * JM_JO[1,4]
  JMH2O <- JOJx * JM_JO[2,1] + JOJv * JM_JO[2,2] + JOJe * JM_JO[2,3] + JOJp * JM_JO[2,4]
  JMO2 <- JOJx * JM_JO[3,1] + JOJv * JM_JO[3,2] + JOJe * JM_JO[3,3] + JOJp * JM_JO[3,4]
  JMNWASTE <- JOJx * JM_JO[4,1] + JOJv * JM_JO[4,2] + JOJe * JM_JO[4,3] + JOJp * JM_JO[4,4]

  JMCO2_GM <- JOJx_GM * JM_JO[1,1] + JOJv_GM * JM_JO[1,2] + JOJe_GM * JM_JO[1,3] + JOJp_GM * JM_JO[1,4]
  JMH2O_GM <- JOJx_GM * JM_JO[2,1] + JOJv_GM * JM_JO[2,2] + JOJe_GM * JM_JO[2,3] + JOJp_GM * JM_JO[2,4]
  JMO2_GM <- JOJx_GM * JM_JO[3,1] + JOJv_GM * JM_JO[3,2] + JOJe_GM * JM_JO[3,3] + JOJp_GM * JM_JO[3,4]
  JMNWASTE_GM <- JOJx_GM * JM_JO[4,1] + JOJv_GM * JM_JO[4,2] + JOJe_GM * JM_JO[4,3] + JOJp_GM * JM_JO[4,4]

  #RQ <- JMCO2/JMO2

  O2FLUX <- -1 * JMO2/(T_REF / Tb / 24.4) * 1000 #mlO2/h, temperature corrected (including SDA)
  CO2FLUX <- JMCO2 / (T_REF / Tb / 24.4) * 1000
  MLO2 <- (-1 * JMO2 * (0.082058 * (Tb + 273.15)) / (0.082058 * 293.15)) * 24.06 * 1000 #mlO2/h, stp
  GH2OMET <- JMH2O * 18.01528 #g metabolic water/h
  #metabolic heat production (Watts) - growth overhead plus dissipation power (maintenance, maturity maintenance,
  #maturation/repro overheads) plus assimilation overheads - correct to 20 degrees so it can be temperature corrected
  #in MET.f for the new guessed Tb
  DEBQMET <- ((1 - kappa_G) * p_G + p_D + (p_X - p_A - p_A * mu_P * eta_PA)) / 3600 / Tcorr

  DRYFOOD <- -1 * JOJx * w_X
  FAECES <- JOJp * w_P
  NWASTE <- JMNWASTE * w_N
  if(pregnant==1){
    wetgonad <- ((cumrepro / mu_E) * w_E) / d_Egg + ((((v_baby * e_baby) / mu_E) * w_E) / d_V + v_baby) * clutchsize
  }else{
    wetgonad <- ((cumrepro/mu_E) * w_E) / d_Egg + ((cumbatch / mu_E) * w_E) / d_Egg
  }
  wetstorage <- ((V * E / mu_E) *w_E) / d_V
  wetfood <- ((Es / mu_E) * w_E) / fdry
  foodin <- ((dEsdt / mu_E) * w_E) / fdry
  #wetfood <- Es / 21525.37 / fdry
  wetmass <- V * andens_deb + wetgonad + wetstorage + wetfood
  gutfreemass <- V * andens_deb + wetgonad + wetstorage
  potfreemass <- V * andens_deb + (((V * E_m) / mu_E) * w_E) / d_V # this is the max potential mass if reserve density is at max value

  dsurvdt <- -1*surviv_pres * hs
  surviv <- surviv_pres + dsurvdt

  # new states
  E_pres <- E
  V_pres <- V
  E_H_pres <- E_H
  q_pres <- q
  hs_pres <- hs
  surviv_pres <- surviv
  Es_pres <- Es

  deb.names <- c("E_pres", "V_pres", "E_H_pres", "q_pres", "hs_pres" ,"surviv_pres", "Es_pres", "cumrepro", "cumbatch", "O2FLUX", "CO2FLUX", "MLO2", "GH2OMET", "DEBQMET", "DRYFOOD", "FAECES", "NWASTE", "wetgonad", "wetstorage", "wetfood", "wetmass", "gutfreemass", "gutfull", "fecundity", "clutches", "potfreemass", "length", "p.R", "foodin", "stage", "p_G", "p_M2", "p_D", "p_J","p_C", "p_A")
  results_deb <- c(E_pres ,V_pres ,E_H_pres, q_pres, hs_pres, surviv_pres, Es_pres, cumrepro, cumbatch, O2FLUX, CO2FLUX, MLO2, GH2OMET, DEBQMET, DRYFOOD, FAECES, NWASTE, wetgonad, wetstorage, wetfood, wetmass, gutfreemass, gutfull, fecundity, clutches, potfreemass, L_w, p_R, foodin, stage, p_G, p_M2, p_D, p_J, p_C, p_A)
  names(results_deb)<-deb.names
  return(results_deb)
}
