#' Dynamic Energy Budget model
#'
#' Implementation of the Standarad Dynamic Energy Budget model of Kooijman
#' Note does not do a proper integration, i.e. rate of change is the difference between two fixed time periods
#' which is ok as long as the timestep is very small relative to the dynamics - hourly time steps or shorter
#' are ok for long-lived (months to years) animals
#' Michael Kearney Dec 2015
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
#' @param TA = 8085 Arhhenius temperature
#' @param TAL = 18721, Arrhenius temperature for decrease below lower boundary of tolerance range \code{TL}
#' @param TAH = 90000, Arrhenius temperature for decrease above upper boundary of tolerance range \code{TH}
#' @param TL = 288, Lower boundary (K) of temperature tolerance range for Arrhenius thermal response
#' @param TH = 315, Upper boundary (K) of temperature tolerance range for Arrhenius thermal response
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
#' @param stage = 0, Initial stage (0=embryo, 1=juvenile, 2=mature but not yet reproducing, 3=beyond first reproduction)
#' @param clutchsize = 2, Clutch size (#), overridden by \code{clutch_ab}
#' @param clutch_ab = c(0,0), paramters for relationship between length (cm) and clutch size: clutch size = a*SVL-b, make a and b zero if fixed clutch size
#' @param viviparous = 0, Viviparous reproduction? 1=yes, 0=no (if yes, animal will be held in adult-sided female's body for duration of development and will experience her body temperature
#' @param minclutch = 0, Minimum clutch size if not enough in reproduction buffer for clutch size predicted by \code{clutch_ab} - if zero, will not operate
#' @param batch = 1, Invoke Pequerie et al.'s batch laying model?
#' @param lambda = 1/2
#' @param VTMIN = 26, Voluntary thermal maximum, degrees C, controls whether repro event can occur at a given time
#' @param VTMAX = 39, Voluntary thermal maximum, degrees C, controls whether repro event can occur at a given time
#' @param arrhenius = matrix(data = matrix(data = c(rep(TA,8),rep(TAL,8),rep(TAH,8),rep(TL,8),rep(TH,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5), Stage-specific 5-parameter Arrhenius thermal response for DEB model (TA,TAL,TAH,TL,TH)
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
#' @param p_B_past = 0
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
#' @return p_B_past
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
#' # simulate growth and reproduction at different constant body temperatures at constant food
#'
#' n<-3000 # time steps
#' step<-1 # step size (days)
#'
#' Tbs=seq(25,35,2.5) # sequence of body temperatures to use
#'
#' for(j in 1:length(Tbs)){
#' debout<-matrix(data = 0, nrow = n, ncol=27)
#' deb.names<-c("E_pres","V_pres","E_H_pres","q_pres","hs_pres","surviv_pres","Es_pres","cumrepro","cumbatch","p_B_past","O2FLUX","CO2FLUX","MLO2","GH2OMET","DEBQMET","DRYFOOD","FAECES","NWASTE","wetgonad","wetstorage","wetfood","wetmass","gutfreemass","gutfull","fecundity","clutches","potfreemass")
#' colnames(debout)<-deb.names
#'
#' # initialise
#' debout[1,]<-DEB(Tb = Tbs[j], step = step)
#'
#' for(i in 2:n){
#' debout[i,]<-DEB(Tb = Tbs[j], breeding = 1, step = step,E_pres=debout[(i-1),1],V_pres=debout[(i-1),2],E_H_pres=debout[(i-1),3],q_pres=debout[(i-1),4],hs_pres=debout[(i-1),5],surviv_pres=debout[(i-1),6],Es_pres=debout[(i-1),7],cumrepro=debout[(i-1),8],cumbatch=debout[(i-1),9],p_B_past=debout[(i-1),10])
#' }
#'
#' if(j==1){
#'   plot((seq(1,n)/365),debout[,23],ylim=c(100,1500),type='l',xlab='years',ylab='wet mass, g', col=j)
#' }else{
#'   points((seq(1,n)/365),debout[,23],ylim=c(100,1500),type='l',xlab='years',ylab='wet mass, g',col=j)
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
  TA=8085,
  TAL=18721,
  TAH=9.0E+04,
  TL=288,
  TH=315,
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
  arrhenius=matrix(data = matrix(data = c(rep(TA,8),rep(TAL,8),rep(TAH,8),rep(TL,8),rep(TH,8)),nrow = 8, ncol = 5), nrow = 8, ncol = 5),
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
  p_B_past=0,
  stage=1,
  breeding=0,
  pregnant=0,
  Tb=33,
  fdry=0.3,
  L_b=0.42,
  L_j=1.376){

  if(clutch_ab[1]>0){
  clutchsize=floor(clutch_ab[1]*(V_pres^(1/3)/del_M)-clutch_ab[2])
  }
  orig_clutchsize <- clutchsize

  q_init<-q_pres
  E_H_init<-E_H_pres
  hs_init<-hs_pres
  fecundity<-0
  clutches<-0
  clutchenergy <- E_0*clutchsize
  starve <- 0
  #DEB mass balance-related calculations
  n_O<-cbind(n_X,n_V,n_E,n_P) # matrix of composition of organics, i.e. food, structure, reserve and faeces
  CHON<-c(12,1,16,14)
  wO<-CHON%*%n_O
  w_V=wO[3]
  M_V<-d_V/w_V
  y_EX<-kap_X*mu_X/mu_E # yield of reserve on food
  y_XE<-1/y_EX # yield of food on reserve
  y_VE<-mu_E*M_V/E_G  # yield of structure on reserve
  y_PX<-kap_X_P*mu_X/mu_P # yield of faeces on food
  y_PE<-y_PX/y_EX # yield of faeces on reserve  0.143382353
  nM<-matrix(c(1,0,2,0,0,2,1,0,0,0,2,0,n_M_nitro),nrow=4)
  n_M_nitro_inv<-c(-1*n_M_nitro[1]/n_M_nitro[4],(-1*n_M_nitro[2])/(2*n_M_nitro[4]),(4*n_M_nitro[1]+n_M_nitro[2]-2*n_M_nitro[3])/(4*n_M_nitro[4]),1/n_M_nitro[4])
  n_M_inv<-matrix(c(1,0,-1,0,0,1/2,-1/4,0,0,0,1/2,0,n_M_nitro_inv),nrow=4)
  JM_JO<--1*n_M_inv%*%n_O
  etaO<-matrix(c(y_XE/mu_E*-1,0,1/mu_E,y_PE/mu_E,0,0,-1/mu_E,0,0,y_VE/mu_E,-1/mu_E,0),nrow=4)
  w_N<-CHON%*%n_M_nitro

  # Arrhenius temperature correction factor
  Tcorr = exp(TA*(1/(273+T_REF)-1/(273+Tb)))/(1+exp(TAL*(1/(273+Tb)-1/TL))+exp(TAH*(1/TH-1/(273+Tb))))

  s <- 1
  if(E_Hj != E_Hb){
    if(E_H_pres < E_Hb){
      s <- 1# % -, multiplication factor for v and {p_Am}
    }else{
      if(E_H_pres < E_Hj){
        s <- V_pres^(1/3) / L_b#
      }else{
        s <- L_j / L_b#
      }
    }
  }

  M_V = d_V/w_V
  p_MT = p_M*Tcorr
  k_Mdot = p_MT/E_G
  k_JT = k_J*Tcorr
  p_AmT = p_MT*z/kap*s
  vT = v*Tcorr*s
  E_m = p_AmT/vT
  F_mT = F_m*Tcorr*s
  g = E_G/(kap*E_m)
  E_scaled=E_pres/E_m
  V_max=(kap*p_AmT/p_MT)^(3.)
  h_aT = h_a*Tcorr
  L_T = 0.
  L_pres = V_pres^(1./3.)
  L_max = V_max^(1./3.)
  scaled_l = L_pres/L_max
  kappa_G = (d_V*mu_V)/(w_V*E_G)
  yEX=kap_X*mu_X/mu_E
  yXE=1/yEX
  yPX=kap_X_P*mu_X/mu_P
  mu_AX=mu_E/yXE
  eta_PA=yPX/mu_AX
  w_X=wO[1]
  w_E=wO[3]
  w_V=wO[2]
  w_P=wO[4]

  # growth of structure
  if(E_H_pres<=E_Hb){
    #use embryo equation for length, from Kooijman 2009 eq. 2
    dLdt=(vT*E_scaled-k_Mdot*g*V_pres^(1./3.))/(3*(E_scaled+g))
    V_temp=(V_pres^(1./3.)+dLdt)^3
    dVdt = V_temp-V_pres
    rdot=vT*(E_scaled/L_pres-(1+L_T/L_pres)/L_max)/(E_scaled+g)
  }else{
    #equation 2.21 from DEB3
    rdot=vT*(E_scaled/L_pres-(1+L_T/L_pres)/L_max)/(E_scaled+g)
    dVdt = V_pres*rdot
    if(dVdt<0){
      dVdt=0
      starve=dVdt*-1*mu_V*d_V/w_V
    }else{
      starve=0
    }
  }
  V=V_pres+dVdt
  if(V<0){
    V=0
  }
  #svl in mm
  svl = V^(0.3333333333333)/del_M*10

  # reserve dynamics
  if(E_H_pres<=E_Hb){
    #use embryo equation for scaled reserve, U_E, from Kooijman 2009 eq. 1
    Sc = L_pres^2*(g*E_scaled)/(g+E_scaled)*(1+((k_Mdot*L_pres)/vT))
    dUEdt = -1*Sc
    E_temp=((E_pres*V_pres/p_AmT)+dUEdt)*p_AmT/(V_pres+dVdt)
    dEdt=E_temp-E_pres
  }else{
    if(Es_pres>0.0000001*E_sm*V_pres){
      dEdt = (p_AmT*f-E_pres*vT)/L_pres
    }else{
      dEdt = (p_AmT*0-E_pres*vT)/L_pres
    }
  }
  E = E_pres+dEdt
  if(E<0){ #make sure E doesn't go below zero
    E=0
  }

  # some powers
  p_M2 = p_MT*V_pres
  p_J = k_JT*E_H_pres
  if(Es_pres>0.0000001*E_sm*V_pres){
    p_A = V_pres^(2./3.)*p_AmT*f
  }else{
    p_A = 0
  }
  p_X = p_A/kap_X #J food eaten per hour
  p_C = (E_m*(vT/L_pres+k_Mdot*(1+L_T/L_pres))*(E_scaled*g)/(E_scaled+g))*V_pres #equation 2.20 DEB3
  p_R = (1-kap)*p_C-p_J-starve

  #maturation
  if(E_H_pres<E_Hp){
    if(E_H_pres<=E_Hb){
      #use embryo equation for scaled maturity, U_H, from Kooijman 2009 eq. 3
      U_H_pres=E_H_pres/p_AmT
      dUHdt=(1-kap)*Sc-k_JT*U_H_pres
      dE_Hdt=dUHdt*p_AmT
    }else{
      dE_Hdt = (1-kap)*p_C-p_J
    }
  }else{
    dE_Hdt = 0
  }
  E_H = E_H_init + dE_Hdt

  # more powers
  if(E_H_pres>=E_Hp){
    p_D = p_M2+p_J+(1-kap_R)*p_R
  }else{
    p_D = p_M2+p_J+p_R
  }
  p_G = p_C-p_M2-p_J-p_R

  # reproduction
  if((E_H_pres<=E_Hp) | (pregnant==1)){
    p_B = 0.
  }else{
    if(batch==1){
      batchprep=(kap_R/lambda)*((1-kap)*(E_m*(vT*V_pres^(2./3.)+k_Mdot*V_pres)/(1+(1/g)))-p_J)
      if(breeding==0){
        p_B = 0
      }else{
        #if the repro buffer is lower than what p_B would be(see below), p_B is p_R
        if(cumrepro<batchprep){
          p_B = p_R
        }else{
          #otherwise it is a faster rate, as specified in Pecquerie et. al JSR 2009 Anchovy paper,
          #with lambda (the fraction of the year the animals breed if food/temperature not limiting) = 0.583 or 7 months of the year
          p_B = batchprep
        }
      }
    }else{
      p_B=p_R
    }#end check for whether batch mode is operating
  }#end check for immature or mature

  #accumulate energy/matter in reproduction buffer
  if(E_H_pres>E_Hp){
    # if buffer ran on previous day
    if(cumrepro<0){
      #       keep it empty, bring it up to 0
      cumrepro=0
    }else{
      cumrepro = cumrepro+p_R*kap_R-p_B_past
    }
  }
  cumbatch = cumbatch+p_B #accumulate energy/matter in egg batch buffer

  # determine stages
  if(stage==2){
    if(cumbatch<0.1*clutchenergy){
      stage=3
    }
  }
  if(E_H<=E_Hb){
    stage=0
  }else{
    if(E_H<E_Hj){
      stage=1
    }else{
      if(E_H<E_Hp){
        stage=2
      }else{
        stage=3
      }
    }
  }
  if(cumbatch>0){
    if(E_H>E_Hp){
      stage=4
    }else{
      stage=stage
    }
  }
      testclutch=floor((cumrepro+cumbatch)/E_0)
#     FOR VARIABLE CLUTCH SIZE FROM REPRO AND BATCH BUFFERS
      if(minclutch > 0 & floor(cumrepro+cumbatch)/E_0 > minclutch){
       if(testclutch <= orig_clutchsize){# ! MAKE SMALLEST CLUTCH ALLOWABLE FOR THIS REPRO EVENT
        clutchsize <- minclutch
        clutchenergy <- clutchsize*E_0
       }
      }

  if((cumbatch>clutchenergy) | (pregnant==1)){
    #for variable clutch size from repro and batch buffers
    #if((cumbatch(hour)>clutchenergy) | (pregnant==1).or
    #&.((viviparous==1) & (cumbatch(hour)+cumrepro(hour)>
    #&clutchenergy))){
    #batch is ready so if viviparous, start gestation, }else{ dump it
    if(viviparous==1){
      if((pregnant==0) & (breeding==1)){
        v_baby=v_init_baby
        e_baby=e_init_baby
        EH_baby=0.
        pregnant=1
        testclutch=floor(cumbatch/E_0)
        #for variable clutch size from repro and batch buffers
        #testclutch=floor((cumbatch(hour)+cumrepro(hour))/E_0)
        #testclutch=real(testclutch)
        if(testclutch>clutchsize){
          clutchsize=testclutch
          clutchenergy = E_0*clutchsize
        }
        #       for variable clutch size from repro and batch buffers
        if(cumbatch<clutchenergy){
          #        needs to draw from repro buffer - temporarily store current repro as cumrepro_temp,
          #        { remove what is needed from the repro buffer and add it to the batch buffer
          cumrepro_temp=cumrepro
          cumrepro=cumrepro+cumbatch-clutchenergy
          cumbatch=cumbatch+cumrepro_temp-cumrepro
        }
      }
      if(hour==1){
        v_baby=v_baby_init
        e_baby=e_baby_init
        EH_baby=EH_baby_init
      }
      #if(pregnant==1){
      #call deb_baby
      #}
      if(EH_baby>E_Hb){
        if((Tb < VTMIN)  |  (Tb > VTMAX)){
          #goto 898
        }
        cumbatch(hour) = cumbatch(hour)-clutchenergy
        repro(hour)=1
        pregnant=0
        v_baby=v_init_baby
        e_baby=e_init_baby
        EH_baby=0
        newclutch=clutchsize
        fecundity=clutchsize
        clutches=1
        pregnant=0
      }
    }else{
      #not viviparous, so lay the eggs at next period of activity
      if((Tb < VTMIN)  |  (Tb > VTMAX)){
      }
      #    change below to active or not active rather than depth-based, in case of fossorial
      if((Tb < VTMIN)  |  (Tb > VTMAX)){
      }
      testclutch=floor(cumbatch/E_0)
      if(testclutch>clutchsize){
        clutchsize=testclutch
        clutchenergy <- clutchsize*E_0
      }
      cumbatch = cumbatch-clutchenergy
      repro=1
      fecundity=clutchsize
      clutches=1
    }
  }

  # feeding (gut) model
  if(E_H_pres>E_Hb){
    if(acthr > 0){
      # Regulates X dynamics
      dEsdt = F_mT*(X/(K+X))*V_pres^(2./3.)*f-1.*(p_AmT/kap_X)*V_pres^(2./3.)
    }else{
      dEsdt = -1.*(p_AmT/kap_X)*V_pres^(2./3.)
    }
  }else{
    dEsdt = -1.*(p_AmT/kap_X)*V_pres^(2./3.)
  }

  if(V_pres==0){
    dEsdt=0
  }
  Es = Es_pres+dEsdt
  if(Es<0){
    Es=0
  }
  if(Es>E_sm*V_pres){
    Es=E_sm*V_pres
  }
  gutfull=Es/(E_sm*V_pres)
  if(gutfull>1){
    gutfull=1
  }

  #mass balance
  JOJx=p_A*etaO[1,1]+p_D*etaO[1,2]+p_G*etaO[1,3]
  JOJv=p_A*etaO[2,1]+p_D*etaO[2,2]+p_G*etaO[2,3]
  JOJe=p_A*etaO[3,1]+p_D*etaO[3,2]+p_G*etaO[3,3]
  JOJp=p_A*etaO[4,1]+p_D*etaO[4,2]+p_G*etaO[4,3]

  JOJx_GM=p_D*etaO[1,2]+p_G*etaO[1,3]
  JOJv_GM=p_D*etaO[2,2]+p_G*etaO[2,3]
  JOJe_GM=p_D*etaO[3,2]+p_G*etaO[3,3]
  JOJp_GM=p_D*etaO[4,2]+p_G*etaO[4,3]

  JMCO2=JOJx*JM_JO[1,1]+JOJv*JM_JO[1,2]+JOJe*JM_JO[1,3]+JOJp*JM_JO[1,4]
  JMH2O=JOJx*JM_JO[2,1]+JOJv*JM_JO[2,2]+JOJe*JM_JO[2,3]+JOJp*JM_JO[2,4]
  JMO2=JOJx*JM_JO[3,1]+JOJv*JM_JO[3,2]+JOJe*JM_JO[3,3]+JOJp*JM_JO[3,4]
  JMNWASTE=JOJx*JM_JO[4,1]+JOJv*JM_JO[4,2]+JOJe*JM_JO[4,3]+JOJp*JM_JO[4,4]

  JMCO2_GM=JOJx_GM*JM_JO[1,1]+JOJv_GM*JM_JO[1,2]+JOJe_GM*JM_JO[1,3]+JOJp_GM*JM_JO[1,4]
  JMH2O_GM=JOJx_GM*JM_JO[2,1]+JOJv_GM*JM_JO[2,2]+JOJe_GM*JM_JO[2,3]+JOJp_GM*JM_JO[2,4]
  JMO2_GM=JOJx_GM*JM_JO[3,1]+JOJv_GM*JM_JO[3,2]+JOJe_GM*JM_JO[3,3]+JOJp_GM*JM_JO[3,4]
  JMNWASTE_GM=JOJx_GM*JM_JO[4,1]+JOJv_GM*JM_JO[4,2]+JOJe_GM*JM_JO[4,3]+JOJp_GM*JM_JO[4,4]

  #RQ = JMCO2/JMO2

  O2FLUX = -1*JMO2/(T_REF/Tb/24.4)*1000 #mlO2/h, temperature corrected (including SDA)
  CO2FLUX = JMCO2/(T_REF/Tb/24.4)*1000
  MLO2 = (-1*JMO2*(0.082058*(Tb+273.15))/(0.082058*293.15))*24.06*1000 #mlO2/h, stp
  GH2OMET = JMH2O*18.01528 #g metabolic water/h
  #metabolic heat production (Watts) - growth overhead plus dissipation power (maintenance, maturity maintenance,
  #maturation/repro overheads) plus assimilation overheads - correct to 20 degrees so it can be temperature corrected
  #in MET.f for the new guessed Tb
  DEBQMET = ((1-kappa_G)*p_G+p_D+(p_X-p_A-p_A*mu_P*eta_PA))/3600/Tcorr

  DRYFOOD=-1*JOJx*w_X
  FAECES=JOJp*w_P
  NWASTE=JMNWASTE*w_N
  if(pregnant==1){
    wetgonad = ((cumrepro/mu_E)*w_E)/d_Egg+((((v_baby*e_baby)/mu_E)*w_E)/d_V + v_baby)*clutchsize
  }else{
    wetgonad = ((cumrepro/mu_E)*w_E)/d_Egg+((cumbatch/mu_E)*w_E)/d_Egg
  }
  wetstorage = ((V*E/mu_E)*w_E)/d_V
  #    wetfood(hour) = ((Es(hour)/mu_E)*w_E)/d_V
  wetfood = Es / 21525.37 / fdry
  wetmass = V*andens_deb+wetgonad+wetstorage+wetfood
  gutfreemass=V*andens_deb+wetgonad+wetstorage
  potfreemass=V*andens_deb+(((V*E_m)/mu_E)*w_E)/d_V # this is the max potential mass if reserve density is at max value

  #aging
  dqdt = (q_pres*(V_pres/V_max)*s_G+h_aT)*(E_pres/E_m)*((vT/L_pres)-rdot)-rdot*q_pres
  if(E_H_pres>E_Hb){
    q = q_init + dqdt
  }else{
    q = 0
  }
  dhsds = q_pres-rdot*hs_pres
  if(E_H_pres>E_Hb){
    hs = hs_init + dhsds
  }else{
    hs = 0
  }
  h_w = ((h_aT*(E_pres/E_m)*vT)/(6*V_pres^(1./3.)))^(1./3.)
  dsurvdt = -1*surviv_pres*hs
  surviv = surviv_pres+dsurvdt

  # new states
  p_B_past=p_B
  E_pres=E
  V_pres=V
  E_H_pres=E_H
  q_pres=q
  hs_pres=hs
  suriv_pres=surviv_pres
  Es_pres=Es

  deb.names<-c("E_pres","V_pres","E_H_pres","q_pres","hs_pres","surviv_pres","Es_pres","cumrepro","cumbatch","p_B_past","O2FLUX","CO2FLUX","MLO2","GH2OMET","DEBQMET","DRYFOOD","FAECES","NWASTE","wetgonad","wetstorage","wetfood","wetmass","gutfreemass","gutfull","fecundity","clutches","potfreemass","length","p.R")
  results_deb<-c(E_pres,V_pres,E_H_pres,q_pres,hs_pres,surviv_pres,Es_pres,cumrepro,cumbatch,p_B_past,O2FLUX,CO2FLUX,MLO2,GH2OMET,DEBQMET,DRYFOOD,FAECES,NWASTE,wetgonad,wetstorage,wetfood,wetmass,gutfreemass,gutfull,fecundity,clutches,potfreemass,svl,p_R)
  names(results_deb)<-deb.names
  return(results_deb)
}
