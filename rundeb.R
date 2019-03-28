rundeb <- function(species = 'Daphnia.magna', ndays = 50, div = 24, Tbs = rep(20, ndays*div),
                   clutchsize = 5, kap.mult = 1, z.mult = 1, v.mult = 1, p.M.mult = 1, mass.mult = 1,
                   F.m.mult = 1000, Xs = rep(100, ndays*div), E.0.mult = E.0.mult,
                   mass.unit = 'g', length.unit = 'cm', plot = 1, start.stage = 0, E_sm = 1000,
                   Euler = 0, stages = 4, S_instar = 1.618){
  library(NicheMapR)
  allDEB.species<-unlist(labels(allStat$allStat))
  if(exists("E.Hj")==TRUE){rm(E.Hj)}
  if(exists("L.j")==TRUE){rm(L.j)}
  species.slot <- which(allDEB.species == species)
  par.names <- unlist(labels(allStat$allStat[[species.slot]]))
  for(i in 1:length(par.names)){
    assign(par.names[i], unlist(allStat$allStat[[species.slot]][i]))
  }
  if(exists("E.Hj")==FALSE){E.Hj <- E.Hb}
  if(exists("L.j")==FALSE){L.j <- L.b}
  metab_mode <- 0 # default (standard model or abj)
  if(abs(E.Hj - E.Hp) < .001){ # abp mode
    metab_mode <- 1
  }

  nyears <- ndays/365
  n<-365*div*nyears # time steps
  step<-1/div # step size (hours)
  T_F_min <- 0
  T_F_max <- 35
  T_AL <- 9.0E04
  T_AH <- 9.0E04
  T_L <- 213
  T_H <- 355
  T_REF <- 20 + 273.15

  # user-specified parameter changes
  kap.orig <- kap
  p.M.orig <- p.M
  z.orig <- z
  kap <- kap * kap.mult
  repro.count <- 0
  if(kap > 1){
    kap <- 1
  }
  F.m <- F.m * F.m.mult
  p.M <- p.M * p.M.mult
  p.Am <- p.M.orig * z / kap.orig
  z <- p.Am * kap / p.M
  v <- v * v.mult

  # covariation rules for scaling
  z <- z * z.mult
  E.0 <- E.0*z.mult^4*E.0.mult
  L.b <- L.b*z.mult
  L.j <- L.j*z.mult
  L.p <- L.p*z.mult
  E.Hb <- E.Hb*z.mult^3
  E.Hj <- E.Hj*z.mult^3
  E.Hp <- E.Hp*z.mult^3
  F.m <- F.m*z.mult^3
  h.a <- h.a*z.mult
  p.Am <- p.M * z / kap # recompute new p.Am
  E.m <- p.Am / v # recompute new E.m

  #sk.J <- p.M/E.G
  # save parameters for output
  pars <- c(F.m = F.m, kap.X = kap.X, kap.P = kap.P, p.Am = p.Am, p.M = p.M, k.J = k.J, v = v, E.m = E.m, E.G = E.G, kap = kap, kap.R = kap.R, E.Hb = E.Hb, E.Hj = E.Hj, E.Hp = E.Hp, h.a = h.a, s.G = s.G, T.A = T.A)

  # initialise DEB output matrix
  if(Euler == 1){
    debout<-matrix(data = 0, nrow = n, ncol=37)
    deb.names<-c("E_pres","V_pres","E_H_pres","q_pres","hs_pres","surviv_pres","Es_pres","cumrepro","cumbatch","p_B_past","O2FLUX","CO2FLUX","MLO2","GH2OMET","DEBQMET","DRYFOOD","FAECES","NWASTE","wetgonad","wetstorage","wetfood","wetmass","gutfreemass","gutfull","fecundity","clutches","potfreemass","length","p.R","foodin","stage", "p_G", "p_M2", "p_D", "p_J", "p_C", "p_A")
  }else{
    debout<-matrix(data = 0, nrow = n, ncol=36)
    deb.names<-c("E_pres","V_pres","E_H_pres","q_pres","hs_pres","surviv_pres","Es_pres","cumrepro","cumbatch","O2FLUX","CO2FLUX","MLO2","GH2OMET","DEBQMET","DRYFOOD","FAECES","NWASTE","wetgonad","wetstorage","wetfood","wetmass","gutfreemass","gutfull","fecundity","clutches","potfreemass","length","p.R","foodin","stage", "p_G", "p_M2", "p_D", "p_J", "p_C", "p_A")
  }
  colnames(debout)<-deb.names

  # initial conditions
  i <- 1
  age <- 0
  if(start.stage == 0){
    V_pres <- 1e-9
    E_pres <- E.0 / V_pres
    E_H_pres <- 0
    Es_pres <- 0
  }
  if(start.stage == 1){
    V_pres <- L.b^3
    E_pres <- E.m
    E_H_pres <- E.Hb
    Es_pres <- 0
  }
  if(start.stage == 2){
    V_pres <- L.p^3
    E_pres <- E.m
    E_H_pres <- E.Hp
    Es_pres <- 0
  }
  if(Euler == 1){
    debout[1,]<-DEB_euler(step=step, z=z, del_M=del.M, F_m=F.m*step, kap_X=kap.X, v=v*step, kap=kap, p_M=p.M*step, E_G=E.G,
                    kap_R=kap.R, k_J=k.J*step, E_Hb=E.Hb, E_Hj=E.Hj, E_Hp=E.Hp, h_a=h.a*(step^2), s_G=s.G, T_REF=T.ref-273.15, T_A=T.A,
                    T_AL=T_AL, T_AH=T_AH, T_L=T_L, T_H=T_H, E_0=E.0, f=f, E_sm=E_sm, K=1, andens_deb=1, d_V=d.V, d_E=d.E, d_Egg=d.E,
                    mu_X=mu.X, mu_E=mu.E, mu_V=mu.V, mu_P=mu.P, kap_X_P=kap.P, n_X=c(1,1.8,0.5,.15), n_E=c(1,1.8,0.5,0.15),
                    n_V=c(1,1.8,0.5,.15), n_P=c(1,1.8,0.5,.15), n_M_nitro=c(1,4/5,3/5,4/5), clutchsize=clutchsize, clutch_ab=c(0,0),
                    viviparous=0, minclutch=clutchsize, batch=0, lambda=1, VTMIN=T_F_min, VTMAX=T_F_max, ma=1e-4, mi=0, mh=0.5,
                    arrhenius=matrix(data = matrix(data = c(rep(T_A,8),rep(T_AL,8),rep(T_AH,8),rep(T_L,8),rep(T_H,8)),nrow = 8, ncol = 5), nrow = 8, ncol = 5), acthr=1, X=Xs[i], E_pres=E_pres, V_pres=V_pres, E_H_pres=E_H_pres, q_pres=0, hs_pres=0,
                    surviv_pres=1, Es_pres=Es_pres, cumrepro=0, cumbatch=0, p_B_past=0, stage=1,breeding=0, pregnant=0, Tb=Tbs[i],fdry=d.V,
                    L_b=L.b, L_j=L.j, spawnday=spawnday, day=day, metab_mode=metab_mode, S_instar=S_instar, stages=stages, age=age)
  }else{
  debout[1,]<-DEB(step=step, z=z, del_M=del.M, F_m=F.m*step, kap_X=kap.X, v=v*step, kap=kap, p_M=p.M*step, E_G=E.G,
                  kap_R=kap.R, k_J=k.J*step, E_Hb=E.Hb, E_Hj=E.Hj, E_Hp=E.Hp, h_a=h.a*(step^2), s_G=s.G, T_REF=T.ref-273.15, T_A=T.A,
                  T_AL=T_AL, T_AH=T_AH, T_L=T_L, T_H=T_H, E_0=E.0, f=f, E_sm=E_sm, K=1, andens_deb=1, d_V=d.V, d_E=d.E, d_Egg=d.E,
                  mu_X=mu.X, mu_E=mu.E, mu_V=mu.V, mu_P=mu.P, kap_X_P=kap.P, n_X=c(1,1.8,0.5,.15), n_E=c(1,1.8,0.5,0.15),
                  n_V=c(1,1.8,0.5,.15), n_P=c(1,1.8,0.5,.15), n_M_nitro=c(1,4/5,3/5,4/5), clutchsize=clutchsize, clutch_ab=c(0,0),
                  viviparous=0, minclutch=clutchsize, batch=0, lambda=1, VTMIN=T_F_min, VTMAX=T_F_max, ma=1e-4, mi=0, mh=0.5,
                  arrhenius=matrix(data = matrix(data = c(rep(T_A,8),rep(T_AL,8),rep(T_AH,8),rep(T_L,8),rep(T_H,8)),nrow = 8, ncol = 5), nrow = 8, ncol = 5), acthr=1, X=Xs[i], E_pres=E_pres, V_pres=V_pres, E_H_pres=E_H_pres, q_pres=0, hs_pres=0,
                  surviv_pres=1, Es_pres=Es_pres, cumrepro=0, cumbatch=0, stage=1,breeding=0, pregnant=0, Tb=Tbs[i],fdry=d.V,
                  L_b=L.b, L_j=L.j, spawnday=spawnday, day=day, metab_mode=metab_mode, S_instar=S_instar, stages=stages, age=age)
  }

  # run through all remaining time steps
  day<-1
  breeding <- 1
  for(i in 2:n){
    day <- day + step
    if(day > 365){
      day <- 1
    }
    if(debout[(i-1),3] > E.Hb){
      age <- age + 1
    }
    spawnday <- day
    if(debout[(i-1),6] > 0.5){ #check if hasn't died of old age
      if(Euler == 1){
      debout[i,]<-DEB_euler(step=step, z=z, del_M=del.M, F_m=F.m*step, kap_X=kap.X, v=v*step, kap=kap, p_M=p.M*step, E_G=E.G,
                      kap_R=kap.R, k_J=k.J*step, E_Hb=E.Hb, E_Hj=E.Hj, E_Hp=E.Hp, h_a=h.a*(step^2), s_G=s.G, T_REF=T.ref-273.15,
                      T_A=T.A, T_AL=T_AL, T_AH=T_AH, T_L=T_L, T_H=T_H, E_0=E.0, f=f, E_sm=E_sm, K=1, andens_deb=1, d_V=d.V, d_E=d.E,
                      d_Egg=d.E, mu_X=mu.X, mu_E=mu.E, mu_V=mu.V, mu_P=mu.P, kap_X_P=kap.P, n_X=c(1,1.8,0.5,.15),
                      n_E=c(1,1.8,0.5,0.15), n_V=c(1,1.8,0.5,.15), n_P=c(1,1.8,0.5,.15), n_M_nitro=c(1,4/5,3/5,4/5),
                      clutchsize=clutchsize, clutch_ab=c(0,0), viviparous=0, minclutch=clutchsize, batch=1,lambda=1,
                      VTMIN=T_F_min, VTMAX=T_F_max, ma=1e-4, mi=0, mh=0.5, arrhenius=matrix(data = matrix(data = c(rep(T_A,8),
                      rep(T_AL,8), rep(T_AH,8), rep(T_L,8), rep(T_H,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5), acthr=1, X=Xs[i],
                      E_pres=debout[(i-1),1], V_pres=debout[(i-1),2], E_H_pres=debout[(i-1),3], q_pres=debout[(i-1),4],
                      hs_pres=debout[(i-1),5], surviv_pres=debout[(i-1),6], Es_pres=debout[(i-1),7],cumrepro=debout[(i-1),8],
                      cumbatch=debout[(i-1),9], p_B_past=debout[(i-1),10], stage=1, breeding=breeding, pregnant=0, Tb=Tbs[i],
                      fdry=d.V, L_b=L.b, L_j=L.j, spawnday=spawnday, day=day, metab_mode=metab_mode, S_instar=S_instar, stages=stages, age=age)
      }else{
      debout[i,]<-DEB(step=step, z=z, del_M=del.M, F_m=F.m*step, kap_X=kap.X, v=v*step, kap=kap, p_M=p.M*step, E_G=E.G,
                      kap_R=kap.R, k_J=k.J*step, E_Hb=E.Hb, E_Hj=E.Hj, E_Hp=E.Hp, h_a=h.a*(step^2), s_G=s.G, T_REF=T.ref-273.15,
                      T_A=T.A, T_AL=T_AL, T_AH=T_AH, T_L=T_L, T_H=T_H, E_0=E.0, f=f, E_sm=E_sm, K=1, andens_deb=1, d_V=d.V, d_E=d.E,
                      d_Egg=d.E, mu_X=mu.X, mu_E=mu.E, mu_V=mu.V, mu_P=mu.P, kap_X_P=kap.P, n_X=c(1,1.8,0.5,.15),
                      n_E=c(1,1.8,0.5,0.15), n_V=c(1,1.8,0.5,.15), n_P=c(1,1.8,0.5,.15), n_M_nitro=c(1,4/5,3/5,4/5),
                      clutchsize=clutchsize, clutch_ab=c(0,0), viviparous=0, minclutch=clutchsize, batch=1,lambda=1,
                      VTMIN=T_F_min, VTMAX=T_F_max, ma=1e-4, mi=0, mh=0.5, arrhenius=matrix(data = matrix(data = c(rep(T_A,8),
                      rep(T_AL,8), rep(T_AH,8), rep(T_L,8), rep(T_H,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5), acthr=1, X=Xs[i],
                      E_pres=debout[(i-1),1], V_pres=debout[(i-1),2], E_H_pres=debout[(i-1),3], q_pres=debout[(i-1),4],
                      hs_pres=debout[(i-1),5], surviv_pres=debout[(i-1),6], Es_pres=debout[(i-1),7],cumrepro=debout[(i-1),8],
                      cumbatch=debout[(i-1),9], stage=1, breeding=breeding, pregnant=0, Tb=Tbs[i],
                      fdry=d.V, L_b=L.b, L_j=L.j, spawnday=spawnday, day=day, metab_mode=metab_mode, S_instar=S_instar, stages=stages, age=age)
      }
    }
  }
  if(plot == 1){
  par(mfrow = c(2,2))
  par(oma = c(2,1,2,2) + 0.1) # margin spacing stuff
  par(mar = c(3,3,1.5,1) + 0.1) # margin spacing stuff
  par(mgp = c(2,1,0)) # margin spacing stuff
  xmax <- n/div
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
  plot(seq(1,n)/div, debout.df$V * mass.mult, type = 'l', xlab = 'Age (days)', ylab = paste0('wet mass (',mass.unit,')'), col = 'dark green', lwd = 2, xlim = c(0, xmax), ylim = c(0, max(debout.df$wetmass)), main = "growth, weight")
  points(seq(1,n)/div, (debout.df$V + debout.df$wetstorage + debout.df$wetgonad)*mass.mult, type = 'l', lwd = 2, col = 'pink')
  points(seq(1,n)/div, (debout.df$V + debout.df$wetstorage) * mass.mult, type = 'l', col = 'grey', lwd =2)
  abline(v = which(debout.df$E_H_pres > E.Hb)[1]/div, lty = 2, col = 'grey')
  abline(v = which(debout.df$E_H_pres > E.Hp)[1]/div, lty = 2, col = 'grey')
  legend(-(n/div/20), max(debout.df$wetmass)*1.1, c('repro. buffer', 'reserve', 'structure'), lty = c(1, 1, 1), col = c("pink", "grey", "dark green"), bty = 'n')
  plot(seq(1,n)/div, debout.df$length * length.mult, type = 'l', xlab = 'Age (days)', ylab = paste0('Length (',length.unit,')'), col = 'black', lwd = 2, xlim = c(0, xmax), main = "growth, length")
  abline(v = which(debout.df$E_H_pres>E.Hb)[1]/div, lty = 2, col = 'grey')
  abline(v = which(debout.df$E_H_pres>E.Hp)[1]/div, lty = 2, col = 'grey')
  #debout.df$cumeggs <- cumsum(debout.df$fecundity)
  #plot(seq(1,n)/div, debout.df$cumeggs, type = 'l', xlab = 'Age (days)', ylab = 'cum. no. of young', col = 'black', lwd = 2, xlim = c(0, xmax), main = "reproduction")
  debout.df$cumfood <- cumsum(debout.df$DRYFOOD)/div
  plot(seq(1,n)/div, debout.df$cumfood*mass.mult, type = 'l', xlab = 'Age (days)', ylab = paste0('cum. dry food (',mass.unit,')'), col = 'black', lwd = 2, xlim = c(0, xmax), main = "food intake")
  abline(v = which(debout.df$E_H_pres>E.Hb)[1]/div, lty = 2, col = 'grey')
  abline(v = which(debout.df$E_H_pres>E.Hp)[1]/div, lty = 2, col = 'grey')
  debout.postembryo <- subset(debout.df, E_H_pres > E.Hb)
  nonrepro.mass <- (debout.postembryo$wetmass-debout.postembryo$wetgonad)*mass.mult
  slope <- lm(log10(debout.postembryo$MLO2) ~ log10((nonrepro.mass)))$coefficients[2]
  plot(log10(nonrepro.mass), log10(debout.postembryo$MLO2), type = 'l', xlab = paste0('log10 wet mass (',mass.unit,')'), ylab = "log10 O2 consumption, ml/hr", col = 'black', lwd = 2, main = paste0('respiration, allometric exponent = ',round(slope, 3)))
  abline(v = log10(nonrepro.mass[which(debout.df$E_H_pres>E.Hb)[1]]), lty = 2, col = 'grey')
  abline(v = log10(nonrepro.mass[which(debout.df$E_H_pres>E.Hp)[1]]), lty = 2, col = 'grey')
  }
  return(list(debout = debout, pars = pars))
}
