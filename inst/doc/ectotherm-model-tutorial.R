## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ---- message=FALSE, warnings=FALSE--------------------------------------
library(NicheMapR)
longlat <- c(146.77, -19.29) # Townsville, northern Australia
micro <- micro_global(loc = longlat)
ecto <- ectotherm() # uses default settings (the Eastern Water Skink)

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE----------
knitr::kable(head(ecto$environ[,c(4:15,25)], 13), digits = 2)
knitr::kable(head(ecto$environ[,c(26:27)], 13), digits = 2)

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE----------
knitr::kable(head(ecto$enbal, 13), digits = 2)

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE----------
knitr::kable(head(ecto$masbal[,c(1:5, 15:17)], 12))

## ---- warning=FALSE, message=FALSE---------------------------------------
Ww_g <- 40        # wet weight of animal (g)
pct_wet <- 0.2    # % of surface area acting as a free-water exchanger
alpha_min <-0.85  # minimum solar absorbtivity (dec %)
alpha_max <- 0.85 # maximum solar absorbtivity (dec %)
shape <- 3        # lizard shape
T_RB_min <- 17.5  # min Tb at which they will attempt to leave retreat
T_B_min <- 17.5   # min Tb at which leaves retreat to bask
T_F_min <- 24     # minimum Tb at which activity occurs
T_F_max <- 34     # maximum Tb at which activity occurs
T_pref <- 30      # preferred Tb (will try and regulate to this)
CT_max <- 40      # critical thermal minimum (affects choice of retreat)
CT_min <- 6       # critical thermal maximum (affects choice of retreat)
mindepth <- 2     # min depth (node, 1-10) allowed
maxdepth <- 10    # max depth (node, 1-10) allowed
shade_seek <- 1   # shade seeking?
burrow <- 1       # can it burrow?
climb <- 0        # can it climb to thermoregulate?
nocturn <- 0      # nocturnal activity
crepus <- 0       # crepuscular activity
diurn <- 1        # diurnal activity
minshades <- rep(0, 12)   # min available shade?
maxshades <- micro$MAXSHADES # max available shade?

## ---- warning=FALSE, message=FALSE---------------------------------------
ecto <- ectotherm(Ww_g = Ww_g, alpha_max = alpha_max, alpha_min = alpha_min, shape = shape, pct_wet = pct_wet,
                  T_F_max = T_F_max, T_F_min = T_F_min, T_B_min = T_B_min, T_RB_min = T_RB_min,
                  CT_max = CT_max, CT_min = CT_min, T_pref = T_pref, mindepth = mindepth, 
                  maxdepth = maxdepth, shade_seek = shade_seek, burrow = burrow, climb = climb,
                  minshades = minshades, nocturn = nocturn, diurn = diurn, crepus = crepus,
                  maxshades = maxshades)

# retrieve output
environ <- as.data.frame(ecto$environ) # behaviour, Tb and environment
enbal <- as.data.frame(ecto$enbal) # heat balance outputs
masbal <- as.data.frame(ecto$masbal) # mass balance outputs
metout <- as.data.frame(micro$metout) # above ground microclimate
environ <- cbind(environ,metout$SOLR) # add solar radiation for activity window plots
colnames(environ)[ncol(environ)] <- "Solar"

# append dates
days <- rep(seq(1,12),24)
days <- days[order(days)]
dates <- days+metout$TIME/60/24-1 # dates for hourly output
dates2 <- seq(1,12,1) # dates for daily output
metout <- cbind(dates,metout)
environ <- cbind(dates,environ)
masbal <- cbind(dates,masbal)
enbal <- cbind(dates,enbal)

## ---- fig.width=7, fig.height=7, fig.show = "hold", message=FALSE, warnings=FALSE, fig.cap="**Body temperature, depth, shade and activity of for the lizard *Eulamprus quoyii* with shade options ranging from 0% to 90%**"----
with(environ, plot(TC ~ dates, ylab = "", xlab="month of year", col = 'black', xlim = c(-0.25, 12), ylim = c(-20, 40), type = "l", yaxt = 'n'))
with(environ, points(ACT * 2 + 7 ~ dates, type = "p", pch = 16, col = "orange"))
with(environ, points(SHADE / 10 - 6 ~ dates, type = "l", col = "dark green"))
with(environ, points(DEP - 10 ~ dates, type = "l", col = "brown"))
abline(ecto$T_F_min, 0, lty = 2, col = 'blue')
abline(ecto$T_F_max, 0, lty = 2, col = 'red')
ytick<-seq(15, 40, by=5)
axis(side=2, at=ytick, labels = TRUE)
mtext(text = c('A', 'B', 'I'), side = 2, line = 1, at = c(11, 9, 7))
ytick<-seq(-6, 4, by=2)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = seq(0, 100, 20), side = 2, line = 1, at = seq(-6, 4, 2), las = 2)
ytick<-seq(-20, -10, by=2)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = rev(seq(0, 100, 20)), side = 2, line = 1, at = seq(-20, -10, 2), las = 2)
abline(h = -10, lty = 2, col = 'grey')
mtext(text = c('body temperature (°C)', 'activity', 'shade (%)', 'depth (cm)'), side = 2, line = 2.5, at = c(30, 9, 0, -15))
text(0.1, c(ecto$T_F_max + 1, ecto$T_F_min + 1), c('T_F_max', 'T_F_min'), col = c('red', 'blue'), cex = 0.75)

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE, fig.cap="**Annual activity window for the lizard *Eulamprus quoyii* with shade options ranging from 0% to 90%**"----
# seasonal activity plot (dark blue = night, light blue = basking, orange = foraging)
forage <- subset(environ, ACT == 2) # get foraging hours
bask <- subset(environ, ACT == 1) # get basking hours
night <- subset(environ, Solar == 0) # get night hours
with(night, plot(TIME ~ DOY, ylab = "Hour of Day", xlab = "Day of Year", pch = 15, cex = 2, 
                 col = 'dark blue')) # nighttime hours
with(forage, points(TIME ~ DOY, pch = 15, cex = 2, col = 'orange')) # foraging Tbs
with(bask, points(TIME ~ DOY, pch = 15, cex = 2, col = 'light blue')) # basking Tbs

## ---- echo=FALSE, fig.width=7, fig.height=7, fig.show = "hold", message=FALSE, warnings=FALSE, fig.cap="**Body temperature, depth, shade and activity of for the lizard *Eulamprus quoyii* with shade options ranging from 0% to 10%**"----
micro <- micro_global(loc = longlat, maxshade = 10)
maxshades <- micro$MAXSHADES

ecto <- ectotherm(Ww_g = Ww_g, alpha_max = alpha_max, alpha_min = alpha_min, pct_wet = pct_wet, T_F_max = T_F_max, T_F_min = T_F_min, T_B_min = T_B_min, T_RB_min = T_RB_min, CT_max = CT_max, CT_min = CT_min, T_pref = T_pref, mindepth = mindepth, maxdepth = maxdepth, shade_seek = shade_seek, burrow = burrow, climb = climb, minshades = minshades, nocturn = nocturn, diurn = diurn, crepus = crepus, maxshades = maxshades)

# retrieve output
environ <- as.data.frame(ecto$environ) # behaviour, Tb and environment
enbal <- as.data.frame(ecto$enbal) # heat balance outputs
masbal <- as.data.frame(ecto$masbal) # mass balance outputs
metout <- as.data.frame(micro$metout) # above ground microclimate
environ <- cbind(environ,metout$SOLR) # add solar radiation for activity window plots
colnames(environ)[ncol(environ)] <- "Solar"

# append dates
days <- rep(seq(1,12),24)
days <- days[order(days)]
dates <- days+metout$TIME/60/24-1 # dates for hourly output
dates2 <- seq(1,12,1) # dates for daily output
metout <- cbind(dates,metout)
environ <- cbind(dates,environ)
masbal <- cbind(dates,masbal)
enbal <- cbind(dates,enbal)

with(environ, plot(TC ~ dates, ylab = "", xlab="month of year", col = 'black', xlim = c(-0.25, 12), ylim = c(-20, 40), type = "l", yaxt = 'n'))
with(environ, points(ACT * 2 + 7 ~ dates, type = "p", pch = 16, col = "orange"))
with(environ, points(SHADE / 10 - 6 ~ dates, type = "l", col = "dark green"))
with(environ, points(DEP - 10 ~ dates, type = "l", col = "brown"))
abline(ecto$T_F_min, 0, lty = 2, col = 'blue')
abline(ecto$T_F_max, 0, lty = 2, col = 'red')
ytick<-seq(15, 40, by=5)
axis(side=2, at=ytick, labels = TRUE)
mtext(text = c('A', 'B', 'I'), side = 2, line = 1, at = c(11, 9, 7))
ytick<-seq(-6, 4, by=2)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = seq(0, 100, 20), side = 2, line = 1, at = seq(-6, 4, 2), las = 2)
ytick<-seq(-20, -10, by=2)
axis(side=2, at=ytick, labels = FALSE)
mtext(text = rev(seq(0, 100, 20)), side = 2, line = 1, at = seq(-20, -10, 2), las = 2)
abline(h = -10, lty = 2, col = 'grey')
mtext(text = c('body temperature (°C)', 'activity', 'shade (%)', 'depth (cm)'), side = 2, line = 2.5, at = c(30, 9, 0, -15))
text(0.1, c(ecto$T_F_max + 1, ecto$T_F_min + 1), c('T_F_max', 'T_F_min'), col = c('red', 'blue'), cex = 0.75)

## ---- fig.width=7, fig.height=5, fig.show = "hold", message=FALSE, warnings=FALSE, fig.cap="**Annual activity window for the lizard *Eulamprus quoyii* with shade options ranging from 0% to 10%**"----
forage<-subset(environ,ACT==2)
bask<-subset(environ,ACT==1)
night<-subset(environ,Solar==0)
day<-subset(environ,Solar==0)
with(night,plot(TIME ~ DOY,ylab="Hour of Day",xlab="Day of Year",pch=15,cex=2,col=
    'dark blue'))
# nighttime hours
with(forage,points(TIME~DOY,pch=15,cex=2,col='orange')) # foraging Tbs
with(bask,points(TIME~DOY,pch=15,cex=2,col='light blue')) # basking Tbs

## ---- warning=FALSE, message=FALSE, echo=FALSE, fig.width=8, fig.height=5----
rm(list=ls())
longlat <- c(146.77, -19.29) # Townsville, northern Australia
nyear <- 5
micro <- micro_global(loc = longlat, timeinterval = 365, nyear = nyear)

load('allstat.Rda') # load the allstat file

species <- "Eulamprus.quoyii"

allDEB.species<-unlist(labels(allStat$allStat)) # get all the species names
allDEB.species<-allDEB.species[1:(length(allDEB.species)-2)] # last two elements are not species
species.slot <- which(allDEB.species == species)
par.names <- unlist(labels(allStat$allStat[[species.slot]]))
# clear possible missing parameters
if(exists("E.Hj")==TRUE){rm(E.Hj)}
if(exists("E.He")==TRUE){rm(E.He)}
if(exists("L.j")==TRUE){rm(L.j)}
if(exists("T.L")==TRUE){rm(T.L)}
if(exists("T.H")==TRUE){rm(T.H)}
if(exists("T.AL")==TRUE){rm(T.AL)}
if(exists("T.AH")==TRUE){rm(T.AH)}

for(i in 1:length(par.names)){
  assign(par.names[i], unlist(allStat$allStat[[species.slot]][i]))
}
# assign possible missing parameters
if(exists("E.Hj")==FALSE){E.Hj <- E.Hb}
if(exists("E.He")==FALSE){E.He <- E.Hb}
if(exists("L.j")==FALSE){L.j <- L.b}
F.m <- p.Am / kap.X # redefining F.m to max possible value
z.mult <- 1         # DEB body size scaling parameter

# assign missing 5-par thermal response curve parameters if necessary
if(exists("T.L")==FALSE){T.L <- CT_min + 273.15}
if(exists("T.H")==FALSE){T.H <- CT_max + 273.15}
if(exists("T.AL")==FALSE){T.AL <- 5E04}
if(exists("T.AH")==FALSE){T.AH <- 9E04}

# overwrite nitrogenous waste indices with those of uric acid (currently ammonia by default)
n.NC <- 1
n.NH <- 4/5
n.NO <- 3/5
n.NN <- 4/5

# morph, behav and water loss
pct_wet <- 0.2    # % of surface area acting as a free-water exchanger
alpha_max <- 0.85 # maximum solar absorptivity
alpha_min <- 0.85 # minimum solar absorptivity
shape <- 3        # animal shape - 3 = lizard
T_RB_min <- 17.5  # min Tb at which they will attempt to leave retreat
T_B_min <- 17.5   # min Tb at which leaves retreat to bask
T_F_min <- 24     # minimum Tb at which activity occurs
T_F_max <- 34     # maximum Tb at which activity occurs
T_pref <- 30      # preferred Tb (will try and regulate to this)
CT_max <- 40      # critical thermal minimum (affects choice of retreat)
CT_min <- 6       # critical thermal maximum (affects choice of retreat)
mindepth <- 2     # min depth (node, 1-10) allowed
maxdepth <- 10    # max depth (node, 1-10) allowed
shade_seek <- 1   # shade seeking?
burrow <- 1       # can it burrow?
climb <- 0        # can it climb to thermoregulate?
minshade <- 0     # min available shade?
maxshade <- 90    # min available shade?
nocturn <- 0      # nocturnal activity
crepus <- 0       # crepuscular activity
diurn <- 1        # diurnal activity
minshades <- rep(0, nyear * 365)   # min available shade?
maxshades <- micro$MAXSHADES # max available shade?

# DEB initial state

# egg
V_init <- 3e-9
E_init <- E.0 / V_init
E_H_init <- 0
stage <- 0

# hatchling
# V_init <- L.b^3
# E_init <- E.m
# E_H_init <- E.Hb
# stage <- 1

# mature
# V_init <- L.p^3
# E_init <- E.m
# E_H_init <- E.Hp+2
# stage <- 2

# reproduction parameters
viviparous <- 1 # live bearing (1) or egg laying (0)
clutchsize <- 5 # how many eggs per clutch?
photostart <- 3 # winter solstice is the start of the reproduction cycle
photofinish <- 2 # autumnal equinox is the end of the reproduction cycle

# run the ectotherm model
ecto<-ectotherm(DEB=1,
                viviparous=viviparous,
                clutchsize=clutchsize,
                z.mult=z.mult,
                shape=shape,
                alpha_max=alpha_max,
                alpha_min=alpha_min,
                T_F_min=T_F_min,
                T_F_max=T_F_max,
                T_B_min=T_B_min,
                T_RB_min=T_RB_min,
                T_pref=T_pref,
                CT_max=CT_max,
                CT_min=CT_min,
                diurn=diurn,
                nocturn=nocturn,
                crepus=crepus,
                shade_seek=shade_seek,
                burrow=burrow,
                climb=climb,
                mindepth=mindepth,
                maxdepth=maxdepth,
                minshade=minshade,
                pct_wet=pct_wet,
                z=z*z.mult,
                del_M=del.M,
                F_m=F.m,
                kap_X=kap.X,
                v=v/24,
                kap=kap,
                p_M=p.M/24,
                E_G=E.G,
                kap_R=kap.R,
                k_J=k.J/24,
                E_Hb=E.Hb*z.mult^3,
                E_Hj=E.Hj*z.mult^3,
                E_Hp=E.Hp*z.mult^3,
                E_He=E.He*z.mult^3,
                h_a=h.a/(24^2),
                s_G=s.G,
                T_REF=T.ref,
                T_A=T.A,
                T_AL=T.AL,
                T_AH=T.AH,
                T_L=T.L,
                T_H=T.H,
                E_0=E.0*z.mult^4,
                f=f,
                d_V=d.V,
                d_E=d.E,
                d_Egg=d.E,
                mu_X=mu.X,
                mu_E=mu.E,
                mu_V=mu.V,
                mu_P=mu.P,
                kap_X_P=kap.P,
                n_X=c(n.CX,n.HX,n.OX,n.ON),
                n_E=c(n.CE,n.HE,n.OE,n.OE),
                n_V=c(n.CV,n.HV,n.OV,n.OV),
                n_P=c(n.CP,n.HP,n.OP,n.OP),
                n_M_nitro=c(n.NC,n.NH,n.NO,n.NN),
                L_b=L.b,
                V_init=V_init,
                E_init=E_init,
                E_H_init=E_H_init,
                stage=stage,
                photostart = photostart,
                photofinish = photofinish)

# retrieve output
environ<-as.data.frame(ecto$environ) # behaviour, Tb and environment
enbal<-as.data.frame(ecto$enbal) # heat balance outputs
masbal<-as.data.frame(ecto$masbal) # mass balance outputs
debout<-as.data.frame(ecto$debout) # DEB model outputs
yearout <- as.data.frame(ecto$yearout) # whole life cycle summary
yearsout <- as.data.frame(ecto$yearsout) # annual summaries
ndays <- nyear * 365
par(mfrow = c(1,1))
plot(seq(1, ndays * 24) / 24, debout$WETMASS, type = 'l', xlab = 'date', 
     ylab = paste0('wet mass (g)'), col = 'pink', lwd = 2, 
     ylim = c(0, max(debout$WETMASS)))
points(seq(1, ndays * 24) / 24, debout$V, type = 'l', xlab = 'date', 
       ylab = paste0('wet mass (g)'), col = 'dark green', lwd = 2)
points(seq(1, ndays * 24) / 24, debout$WETMASS-debout$WETGONAD, type = 'l', 
       lwd = 2, col = 'brown')
points(seq(1, ndays * 24) / 24, debout$WETMASS-debout$WETGONAD-debout$WETGUT,
       type = 'l', lwd = 2, col = 'grey')
abline(v = (seq(1, ndays * 24) / 24)[which(debout$E_H>E.Hb)[1]], lty = 2, col = 'grey')
abline(v = (seq(1, ndays * 24) / 24)[which(debout$E_H>E.Hp)[1]], lty = 2, col = 'grey')
legend(1200, max(debout$WETMASS) * 0.3, 
       c('repro. buffer', 'food in gut', 'reserve', 'structure'), lty = rep(1, 4), 
       col = c("pink", "brown", "grey", "dark green"), bty = 'n')
text(0, max(debout$WETMASS) * 1, labels = "embryo", cex = 0.85)
text((which(debout$E_H > E.Hp)[1] - which(debout$E_H > E.Hp)[1] * .5) / 24 , max(debout$WETMASS) * 1, labels = "immature", cex = 0.85)
text(which(debout$E_H > E.Hp)[1] * 1.2 / 24, max(debout$WETMASS) * 1, labels = "adult", cex = 0.85)

