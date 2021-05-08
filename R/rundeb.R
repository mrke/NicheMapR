#' rundeb
#'
#' Function to simulate the development, growth and reproduction trajectory of
#' an organism using DEB theory, drawing parameters from the 'AmP' parameter
#' database at https://github.com/add-my-pet/AmPtool. It requires the
#' 'allStat.mat' file to have been converted to 'allStat.Rda' via the R.matlab
#' package (i.e. allStat <- readMat('allStat.mat' then
#' save(allStat, file = 'allstat.Rda'))).
#' @param allstat = allstat, the allstat data set
#' @param species = 'Daphnia.magna', a species in the allstat file
#' @param Euler = 0, use Euler integration? (faster, but less accurate), 0 or 1
#' @param start.stage = 0, stage at which simulation starts, 0=embryo, 1=birth, 2=puberty
#' @param stages = 6, how many life cycle stages? (e.g. 6 for a abp insect model), #
#' @param S_instar = rep(1.618, stages), 'stress' factor determining body lengths at which molts occur (for abp insect model), -
#' @param ndays = 50, days to run the simulation for
#' @param div = 24, value to divide the default step size (days) by, which determines output frequency (and integration frequency if Euler = 1)
#' @param Tbs = rep(20, ndays*div), vector of body temperatures for each time step, °C
#' @param Xs = rep(100, ndays*div), vector of food densities,J/cm2 or J/cm3
#' @param E_sm = 350, maximum stomach energy density per structure, J/cm3
#' @param fdry = 0.3, dry fraction of food, -
#' @param clutchsize = 5, eggs in a clutch, #
#' @param kap.mult = 1, multiplier on original value for kappa, -
#' @param z.mult = 1, multiplier on original value of z, to explore implications of DEB covariation rules for body size scaling, -
#' @param v.mult = 1, multiplier on original value of energy conductance, v, -
#' @param p.M.mult = 1, multiplier on original value of somatic maintenance p.M, -
#' @param E.0.mult = 1, multiplier on original value of egg energy E.0, -
#' @param plot = 1, produce example plots? 0=no, 1=yes
#' @param mass.unit = 'g', mass unit for the plots, 'mg', 'g' or 'kg'
#' @param length.unit = 'cm', length unit for the plots, 'mm', 'cm' or 'm'
#' @param ageing = 1, impose ageing? 0=immortal, 1=aging according to parameter h.a
#' @examples
#' #library(R.matlab)
#' #allStat<-readMat('allStat.mat') # this will take a few minutes
#' #save(allStat, file = 'allstat.Rda') # save it as an R data file for faster future loading
#' load('allStat.Rda')
#' species <- "Daphnia.magna" # must be in the AmP collection - see allDEB.species list
#' Euler <- 0 # use Euler integration (faster but less accurate) (1) or deSolve's ODE solve (0)
#' start.stage <- 0 # stage in life cycle to start (0 = egg, 1 = juvenile, 2 = puberty)
#' ndays <- 50 # number days to run the simulation for
#' div <- 24 # time step divider (1 = days, 24 = hours, etc.) - keep small if using Euler method for integration
#' Tbs <- rep(20, ndays * div) # deg C, body temperature
#' starvetime <- 0 # length of low food period when simulating starvation
#' X <- 100 # J/cm2 base food density
#' Xs <- c(rep(X, ndays * div/2), rep(0.000005, starvetime), rep(X, ndays * div/2-starvetime + 1)) # body temperature
#' E_sm <- 350 # J/cm3, volume-specific stomach energy
#' fdry <- 0.3 # J/cm3, volume-specific stomach energy
#' clutchsize <- 50 # -, clutch size
#' kap.mult <- 1 # -, factor to multiply kappa by (caps at 1)
#' p.M.mult <- 1 # -, factor to multiply p.M by (caps at 1)
#' v.mult <- 1 # -, factor to multiply energy conductance by
#' z.mult <- 1 # -, factor to multiply zoom factor by
#' E.0.mult <- 1 # -, factor by which to multiply initial egg energy
#' plot <- 1 # plot results?
#' mass.unit <- 'g'
#' length.unit <- 'mm'
#'
#' deb <- rundeb(species = species, ndays = ndays, div = div, Tbs = Tbs,
#'               clutchsize = clutchsize, kap.mult = kap.mult, v.mult = v.mult,
#'               p.M.mult = p.M.mult, Xs = Xs, z.mult = z.mult, E.0.mult = E.0.mult,
#'               mass.unit = mass.unit, length.unit = length.unit, start.stage = start.stage,
#'               E_sm = E_sm, Euler = Euler, plot = plot)
#'
#' debout <- as.data.frame(deb$debout) # retrieve the output
#' parameters <- deb$pars # retrieve the extracted parameters
#' @export
rundeb <- function(
  allstat = allstat,
  species = 'Daphnia.magna',
  Euler = 0,
  start.stage = 0,
  stages = 6,
  S_instar = rep(1.618, stages),
  ndays = 50,
  div = 24,
  Tbs = rep(20, ndays*div),
  Xs = rep(100, ndays*div),
  E_sm = 350,
  fdry = 0.3,
  clutchsize = 5,
  kap.mult = 1,
  z.mult = 1,
  v.mult = 1,
  p.M.mult = 1,
  E.0.mult = 1,
  plot = 1,
  mass.unit = 'g',
  length.unit = 'cm',
  ageing = 1){ # end function parameters

  n <- div * ndays # time steps
  step<-1/div # step size (hours)

  allDEB.species <- unlist(labels(allStat$allStat))
  if(exists("E.Hj") == TRUE){rm(E.Hj)} # remove previous value from memory if present from prior run
  if(exists("L.j") == TRUE){rm(L.j)} # remove previous value from memory if present from prior run
  species.slot <- which(allDEB.species == species) # find the slot with the species of interest
  par.names <- unlist(labels(allStat$allStat[[species.slot]])) # get parameter names
  for(i in 1:length(par.names)){ # pull parameters into memory
    assign(par.names[i], unlist(allStat$allStat[[species.slot]][i]))
  }
  if(exists("E.Hj")==FALSE){E.Hj <- E.Hb} # if no acceleration, give E.Hj the value at birth
  if(exists("L.j")==FALSE){L.j <- L.b} # if no acceleration, give L.j the value at birth
  if(exists("del.M")==FALSE){ # if no del.M, assume 0.2 and warn
    del.M <- 0.2
    cat("NOTE: No length data in parameter estimation, assuming del.M of 0.2, predictions of physical lengths are suspect.", '\n')
  }
  metab_mode <- 0 # default (standard model or abj)
  if(abs(E.Hj - E.Hp) < .001){ # check if abp mode
    metab_mode <- 1
  }

  # checking if 5 parameter Arrhenius function has been defined
  if(exists("T.L")==FALSE){T_L <- 173.15}else{T_L <- T.L}
  if(exists("T.H")==FALSE){T_H <- 373.15}else{T_H <- T.H}
  if(exists("T.AL")==FALSE){T_AL <- 5E04}else{T_AL <- T.AL}
  if(exists("T.AH")==FALSE){T_AH <- 9E04}else{T_AH <- T.AH}
  if(exists("T.ref")==FALSE){T_REF <- 20 + 273.15}else{T_REF <- T.ref}

  # user-specified parameter changes
  kap.orig <- kap # keep original value
  p.M.orig <- p.M # keep original value
  z.orig <- z # keep original value
  kap <- kap * kap.mult # multiply by chosen value
  repro.count <- 0
  if(kap > 1){
    kap <- 1
  }
  v <- v * v.mult # multiply by chosen value
  p.M <- p.M * p.M.mult # multiply by chosen value
  p.Am <- p.M.orig * z / kap.orig # get new p.Am if kap.mult != 1
  z <- p.Am * kap / p.M # get new z if p.M.mult != 1
  # covariation rules for scaling
  z <- z * z.mult # multiply by chosen value
  E.0 <- E.0 * z.mult ^ 4 * E.0.mult # multiply by chosen value
  L.b <- L.b * z.mult # multiply by chosen value
  L.j <- L.j * z.mult # multiply by chosen value
  L.p <- L.p * z.mult # multiply by chosen value
  E.Hb <- E.Hb * z.mult ^ 3 # multiply by chosen value
  E.Hj <- E.Hj * z.mult ^ 3 # multiply by chosen value
  E.Hp <- E.Hp * z.mult ^ 3 # multiply by chosen value
  p.Xm <- p.Xm * z.mult ^ 3 # multiply by chosen value
  h.a <- h.a * z.mult # multiply by chosen value
  p.Am <- p.M * z / kap # recompute new p.Am
  E.m <- p.Am / v # recompute new E.m
  p.Xm <- p.Am / kap.X * 5 # food intake a large value - rapid stomach fill
  if(ageing == 0){
    h.a <- 0
  }
  # save parameters for output
  pars <- c(p.Xm = p.Xm, kap.X = kap.X, kap.P = kap.P, p.Am = p.Am, p.M = p.M, k.J = k.J, v = v, E.m = E.m, E.G = E.G, kap = kap, kap.R = kap.R, E.Hb = E.Hb, E.Hj = E.Hj, E.Hp = E.Hp, h.a = h.a, s.G = s.G, T.A = T.A)

  # initialise DEB output matrix
  deb.names <- c("stage", "V", "E", "E_H", "E_s", "E_R", "E_B", "q", "hs", "length", "wetmass", "wetgonad", "wetgut", "wetstorage", "p_surv", "fecundity", "clutches", "JMO2", "JMCO2", "JMH2O", "JMNWASTE", "O2ML", "CO2ML", "GH2OMET", "DEBQMETW", "GDRYFOOD", "GFAECES", "GNWASTE", "p_A", "p_C", "p_M", "p_G", "p_D", "p_J", "p_R", "p_B", "L.b", "L.j")
  debout <- matrix(data = 0, nrow = n, ncol=38)
  colnames(debout) <- deb.names

  # initial conditions
  i <- 1
  age <- 0
  if(start.stage == 0){ # egg
    V_pres <- 3e-9
    E_pres <- E.0 / V_pres
    E_H_pres <- 0
    E_s <- 0
  }
  if(start.stage == 1){ # birth
    V_pres <- L.b^3
    E_pres <- E.m
    E_H_pres <- E.Hb
    E_s <- 0
  }
  if(start.stage == 2){ # puberty
    V_pres <- L.p^3
    E_pres <- E.m
    E_H_pres <- E.Hp
    E_s <- 0
  }
  if(Euler == 1){
    debout[1,]<-DEB_euler(step=step,
                          z=z,
                          del_M=del.M,
                          p_Xm=p.Xm*step,
                          kap_X=kap.X,
                          v=v*step,
                          kap=kap,
                          p_M=p.M*step,
                          E_G=E.G,
                          kap_R=kap.R,
                          k_J=k.J*step,
                          E_Hb=E.Hb,
                          E_Hj=E.Hj,
                          E_Hp=E.Hp,
                          h_a=h.a*(step^2),
                          s_G=s.G,
                          T_REF=T.ref,
                          T_A=T.A,
                          T_AL=T_AL,
                          T_AH=T_AH,
                          T_L=T_L,
                          T_H=T_H,
                          E_0=E.0,
                          f=f,
                          E_sm=E_sm,
                          K=1,
                          andens_deb=1,
                          d_V=d.V,
                          d_E=d.E,
                          d_Egg=d.E,
                          mu_X=mu.X,
                          mu_E=mu.E,
                          mu_V=mu.V,
                          mu_P=mu.P,
                          kap_X_P=kap.P,
                          n_X=c(1,1.8,0.5,0.15),
                          n_E=c(1,1.8,0.5,0.15),
                          n_V=c(1,1.8,0.5,0.15),
                          n_P=c(1,1.8,0.5,0.15),
                          n_M_nitro=c(1,4/5,3/5,4/5),
                          clutchsize=clutchsize,
                          clutch_ab=c(0,0),
                          minclutch=clutchsize,
                          batch=0,
                          lambda=1,
                          VTMIN=-100,
                          VTMAX=100,
                          ma=1e-4,
                          mi=0,
                          mh=0.5,
                          arrhenius=matrix(data = matrix(data = c(rep(T_A,8),rep(T_AL,8),rep(T_AH,8),rep(T_L,8),rep(T_H,8)),nrow = 8, ncol = 5), nrow = 8, ncol = 5),
                          acthr=1,
                          X=Xs[i],
                          E_pres=E_pres,
                          V_pres=V_pres,
                          E_H_pres=E_H_pres,
                          q=0,
                          hs=0,
                          p_surv=1,
                          E_s=E_s,
                          p_B_pres=0,
                          E_R=0,
                          E_B=0,
                          stage=start.stage,
                          breeding=0,
                          Tb=Tbs[i],
                          fdry=d.V,
                          L_b=L.b,
                          L_j=L.j,
                          spawnday=spawnday,
                          day=day,
                          metab_mode=metab_mode,
                          S_instar=S_instar,
                          stages=stages,
                          age=age)
  }else{
    debout[1,]<-DEB(step=step,
                    z=z,
                    del_M=del.M,
                    p_Xm=p.Xm*step,
                    kap_X=kap.X,
                    v=v*step,
                    kap=kap,
                    p_M=p.M*step,
                    E_G=E.G,
                    kap_R=kap.R,
                    k_J=k.J*step,
                    E_Hb=E.Hb,
                    E_Hj=E.Hj,
                    E_Hp=E.Hp,
                    h_a=h.a*(step^2),
                    s_G=s.G,
                    T_REF=T.ref,
                    T_A=T.A,
                    T_AL=T_AL,
                    T_AH=T_AH,
                    T_L=T_L,
                    T_H=T_H,
                    E_0=E.0,
                    f=f,
                    E_sm=E_sm,
                    K=1,
                    andens_deb=1,
                    d_V=d.V,
                    d_E=d.E,
                    d_Egg=d.E,
                    mu_X=mu.X,
                    mu_E=mu.E,
                    mu_V=mu.V,
                    mu_P=mu.P,
                    kap_X_P=kap.P,
                    n_X=c(1,1.8,0.5,0.15),
                    n_E=c(1,1.8,0.5,0.15),
                    n_V=c(1,1.8,0.5,0.15),
                    n_P=c(1,1.8,0.5,0.15),
                    n_M_nitro=c(1,4/5,3/5,4/5),
                    clutchsize=clutchsize,
                    clutch_ab=c(0,0),
                    minclutch=clutchsize,
                    batch=0,
                    lambda=1,
                    VTMIN=-100,
                    VTMAX=100,
                    ma=1e-4,
                    mi=0,
                    mh=0.5,
                    arrhenius=matrix(data = matrix(data = c(rep(T_A,8),rep(T_AL,8),rep(T_AH,8),rep(T_L,8),rep(T_H,8)),nrow = 8, ncol = 5), nrow = 8, ncol = 5),
                    acthr=1,
                    X=Xs[i],
                    E_pres=E_pres,
                    V_pres=V_pres,
                    E_H_pres=E_H_pres,
                    q=0,
                    hs=0,
                    p_surv=1,
                    E_s=E_s,
                    E_R=0,
                    E_B=0,
                    stage=start.stage,
                    breeding=0,
                    Tb=Tbs[i],
                    fdry=d.V,
                    L_b=L.b,
                    L_j=L.j,
                    spawnday=spawnday,
                    day=day,
                    metab_mode=metab_mode,
                    S_instar=S_instar,
                    stages=stages,
                    age=age)
  }

  # run through all remaining time steps
  day <- 1
  breeding <- 1
  for(i in 2:n){
    day <- day + step
    if(day > 365){
      day <- 1
    }
    if(debout[(i-1), 14] > E.Hb){
      age <- age + 1
    }
    spawnday <- day
    if(debout[(i-1),15] > 0.5){ #check if hasn't died of old age
      if(Euler == 1){
        debout[i,]<-DEB_euler(step=step,
                              z=z,
                              del_M=del.M,
                              p_Xm=p.Xm*step,
                              kap_X=kap.X,
                              v=v*step,
                              kap=kap,
                              p_M=p.M*step,
                              E_G=E.G,
                              kap_R=kap.R,
                              k_J=k.J*step,
                              E_Hb=E.Hb,
                              E_Hj=E.Hj,
                              E_Hp=E.Hp,
                              h_a=h.a*(step^2),
                              s_G=s.G,
                              T_REF=T.ref,
                              T_A=T.A,
                              T_AL=T_AL,
                              T_AH=T_AH,
                              T_L=T_L,
                              T_H=T_H,
                              E_0=E.0,
                              f=f,
                              E_sm=E_sm,
                              K=1,
                              andens_deb=1,
                              d_V=d.V,
                              d_E=d.E,
                              d_Egg=d.E,
                              mu_X=mu.X,
                              mu_E=mu.E,
                              mu_V=mu.V,
                              mu_P=mu.P,
                              kap_X_P=kap.P,
                              n_X=c(1,1.8,0.5,0.15),
                              n_E=c(1,1.8,0.5,0.15),
                              n_V=c(1,1.8,0.5,0.15),
                              n_P=c(1,1.8,0.5,0.15),
                              n_M_nitro=c(1,4/5,3/5,4/5),
                              clutchsize=clutchsize,
                              clutch_ab=c(0,0),
                              minclutch=clutchsize,
                              batch=1,
                              lambda=1,
                              VTMIN=-100,
                              VTMAX=100,
                              ma=1e-4,
                              mi=0,
                              mh=0.5,
                              arrhenius=matrix(data = matrix(data = c(rep(T_A,8), rep(T_AL,8), rep(T_AH,8), rep(T_L,8), rep(T_H,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5),
                              acthr=1,
                              X=Xs[i],
                              E_pres=debout[(i-1),3],
                              V_pres=debout[(i-1),2],
                              E_H_pres=debout[(i-1),4],
                              q=debout[(i-1),8],
                              hs=debout[(i-1),9],
                              p_surv=debout[(i-1),15],
                              E_s=debout[(i-1),5],
                              p_B_pres=debout[(i-1),36],
                              E_R=debout[(i-1),6],
                              E_B=debout[(i-1),7],
                              stage=debout[(i-1),1],
                              breeding=breeding,
                              Tb=Tbs[i],
                              fdry=d.V,
                              L_b=L.b,
                              L_j=L.j,
                              spawnday=spawnday,
                              day=day,
                              metab_mode=metab_mode,
                              S_instar=S_instar,
                              stages=stages,
                              age=age)
      }else{
        debout[i,]<-DEB(step=step,
                        z=z,
                        del_M=del.M,
                        p_Xm=p.Xm*step,
                        kap_X=kap.X,
                        v=v*step,
                        kap=kap,
                        p_M=p.M*step,
                        E_G=E.G,
                        kap_R=kap.R,
                        k_J=k.J*step,
                        E_Hb=E.Hb,
                        E_Hj=E.Hj,
                        E_Hp=E.Hp,
                        h_a=h.a*(step^2),
                        s_G=s.G,
                        T_REF=T.ref,
                        T_A=T.A,
                        T_AL=T_AL,
                        T_AH=T_AH,
                        T_L=T_L,
                        T_H=T_H,
                        E_0=E.0,
                        f=f,
                        E_sm=E_sm,
                        K=1,
                        andens_deb=1,
                        d_V=d.V,
                        d_E=d.E,
                        d_Egg=d.E,
                        mu_X=mu.X,
                        mu_E=mu.E,
                        mu_V=mu.V,
                        mu_P=mu.P,
                        kap_X_P=kap.P,
                        n_X=c(1,1.8,0.5,0.15),
                        n_E=c(1,1.8,0.5,0.15),
                        n_V=c(1,1.8,0.5,0.15),
                        n_P=c(1,1.8,0.5,0.15),
                        n_M_nitro=c(1,4/5,3/5,4/5),
                        clutchsize=clutchsize,
                        clutch_ab=c(0,0),
                        minclutch=clutchsize,
                        batch=1,
                        lambda=1,
                        VTMIN=-100,
                        VTMAX=100,
                        ma=1e-4,
                        mi=0,
                        mh=0.5,
                        arrhenius=matrix(data = matrix(data = c(rep(T_A,8), rep(T_AL,8), rep(T_AH,8), rep(T_L,8), rep(T_H,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5),
                        acthr=1,
                        X=Xs[i],
                        E_pres=debout[(i-1),3],
                        V_pres=debout[(i-1),2],
                        E_H_pres=debout[(i-1),4],
                        q=debout[(i-1),8],
                        hs=debout[(i-1),9],
                        p_surv=debout[(i-1),15],
                        E_s=debout[(i-1),5],
                        p_B_pres=debout[(i-1),36],
                        E_R=debout[(i-1),6],
                        E_B=debout[(i-1),7],
                        stage=debout[(i-1),1],
                        breeding=breeding,
                        Tb=Tbs[i],
                        fdry=d.V,
                        L_b=L.b,
                        L_j=L.j,
                        spawnday=spawnday,
                        day=day,
                        metab_mode=metab_mode,
                        S_instar=S_instar,
                        stages=stages,
                        age=age)
      }
    }
  }
  if(plot == 1){ # make some plots
    par(mfrow = c(2, 2))
    par(oma = c(2, 1, 2, 2) + 0.1) # margin spacing
    par(mar = c(3, 3, 1.5, 1) + 0.1) # margin spacing
    par(mgp = c(2, 1, 0)) # margin spacing
    xmax <- n / div
    if(mass.unit == 'mg'){
      mass.mult <- 1000
    }
    if(mass.unit == 'g'){
      mass.mult <- 1
    }
    if(mass.unit == 'kg'){
      mass.mult <- 0.001
    }
    if(length.unit == 'cm'){
      length.mult <- .1
    }
    if(length.unit == 'mm'){
      length.mult <- 1
    }
    debout.df <- as.data.frame(debout)
    plot(seq(1, n) / div, debout.df$V * mass.mult, type = 'l', xlab = 'Age (days)', ylab = paste0('wet mass (', mass.unit, ')'), col = 'dark green', lwd = 2, xlim = c(0, xmax), ylim = c(0, max(debout.df$wetmass) * mass.mult), main = "growth, weight")
    points(seq(1, n) / div, (debout.df$V + debout.df$wetstorage + debout.df$wetgonad + debout.df$wetgut) * mass.mult, type = 'l', lwd = 2, col = 'pink')
    points(seq(1, n) / div, (debout.df$V + debout.df$wetstorage + debout.df$wetgut) * mass.mult, type = 'l', col = 'brown', lwd =2)
    points(seq(1, n) / div, (debout.df$V + debout.df$wetstorage) * mass.mult, type = 'l', col = 'grey', lwd =2)
    abline(v = which(debout.df$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
    abline(v = which(debout.df$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
    legend(-(n/div/20), max(debout.df$wetmass) * 1.05, c('repro. buffer', 'food in gut', 'reserve', 'structure'), lty = c(1, 1, 1, 1), col = c("pink", "brown", "grey", "dark green"), bty = 'n')
    plot(seq(1, n) / div, debout.df$length * length.mult, type = 'l', xlab = 'Age (days)', ylab = paste0('Length (', length.unit, ')'), col = 'black', lwd = 2, xlim = c(0, xmax), main = "growth, length")
    abline(v = which(debout.df$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
    abline(v = which(debout.df$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
    debout.df$cumfood <- cumsum(debout.df$GDRYFOOD) / d.V
    plot(seq(1, n) / div, debout.df$cumfood * mass.mult, type = 'l', xlab = 'Age (days)', ylab = paste0('cum. wet food (', mass.unit, ')'), col = 'black', lwd = 2, xlim = c(0, xmax), main = "food intake")
    abline(v = which(debout.df$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
    abline(v = which(debout.df$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
    debout.postembryo <- subset(debout.df, E_H > E.Hb)
    nonrepro.mass <- (debout.postembryo$wetmass - debout.postembryo$wetgonad) * mass.mult
    slope <- lm(log10(debout.postembryo$O2ML) ~ log10((nonrepro.mass)))$coefficients[2]
    plot(nonrepro.mass, debout.postembryo$O2ML, type = 'l', xlab = paste0('wet mass (', mass.unit, ')'), ylab = "O2 consumption, ml/hr", col = 'black', lwd = 2, main = paste0('respiration, allometric exponent = ',round(slope, 3)))
    abline(v = debout.df$wetmass[which(debout.df$E_H > E.Hb)[1]] * mass.mult, lty = 2, col = 'grey')
    abline(v = debout.df$wetmass[which(debout.df$E_H > E.Hp)[1]] * mass.mult, lty = 2, col = 'grey')
  }
  return(list(debout = debout, pars = pars))
}
