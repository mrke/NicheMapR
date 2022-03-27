## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  eval = TRUE
)

## ----fig1, warning=FALSE, message=FALSE, echo=FALSE, fig.width=8, fig.height=5, eval=TRUE, fig.cap="DEB model prediction of length through time of the lizard *Eulamprus quoyii* growing under constant, *ad libitum* food at 25 &deg;C. Photo by the M. Kearney"----
library(NicheMapR)
library(knitr) # this packages has a function for producing formatted tables.

load('allStat.Rda')
species <- "Eulamprus.quoyii" # must be in the AmP collection - see allDEB.species list
species.name <- gsub(pattern = "[.]", " ", species)
ndays <- 365*5 # number days to run the simulation for
div <- 1 # time step divider (1 = days, 24 = hours, etc.) - keep small if using Euler method
Tbs <- rep(25, ndays * div) # °C, body temperature
starvetime <- 0 # length of low food period when simulating starvation
X <- 100 # J/cm2 base food density
Xs <- c(rep(X, ndays * div / 2), rep(0.000005, starvetime), rep(X, ndays * div / 2 - starvetime + 1)) # food density (J/cm2 or J/cm3)
E_sm <- 350 # J/cm3, volume-specific stomach energy
clutchsize <- 5 # -, clutch size

mass.unit <- 'g'
length.unit <- 'mm'
Euler <- 0 # use Euler integration (faster but less accurate) (1) or deSolve's ODE solve (0)
plot <- 0 # plot results?
start.stage <- 0 # stage in life cycle to start (0 = egg, 1 = juvenile, 2 = puberty)

deb <- rundeb(species = species, ndays = ndays, div = div, Tbs = Tbs, clutchsize = clutchsize, Xs = Xs, mass.unit = mass.unit, length.unit = length.unit, start.stage = start.stage, E_sm = E_sm, Euler = Euler, plot = plot)

debout <- as.data.frame(deb$debout)
debpars <- as.data.frame(deb$pars)
for(i in 1:length(deb$pars)){
assign(names(deb$pars[i]), deb$pars[i])
}

plot(seq(1, ndays) / div, debout$length, type = 'l', xlab = 'Age (days)', ylab = paste0('Length (', length.unit, ')'), col = 'black', lwd = 2)
abline(v = which(debout$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
abline(v = which(debout$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
text(0, max(debout$length) * 1, labels = "embryo", cex = 0.85)
text(which(debout$E_H > E.Hp)[1] / div - which(debout$E_H > E.Hp)[1] / div * .5, max(debout$length) * 1, labels = "immature", cex = 0.85)
text(which(debout$E_H > E.Hp)[1] / div * 1.2, max(debout$length) * 1, labels = "adult", cex = 0.85)
library(grid)
library(png)
img<-readPNG("Eulamprus_tympanum.png")
grid.raster(img, width = 0.38, height = 0.38, x = .7)

## ----fig2, warning=FALSE, message=FALSE, echo=FALSE, fig.width=8, fig.height=5, eval=TRUE, fig.cap="DEB model prediction of wet weight through time of the lizard *Eulamprus quoyii* growing under constant *ad libitum* food at 25 &deg;C, partitioning total biomass mass into structure, reserve, reproduction buffer and gut contents."----
plot(seq(1, ndays) / div, debout$V, type = 'l', xlab = 'Age (days)', ylab = paste0('wet mass (', mass.unit, ')'), col = 'dark green', lwd = 2, ylim = c(0, max(debout$wetmass)))
points(seq(1, ndays) / div, debout$V + debout$wetstorage + debout$wetgonad + debout$wetgut, type = 'l', lwd = 2, col = 'pink')
points(seq(1, ndays) / div, debout$V + debout$wetstorage + debout$wetgut, type = 'l', col = 'brown', lwd = 2)
points(seq(1, ndays) / div, debout$V + debout$wetstorage, type = 'l', col = 'grey', lwd =2)
abline(v = which(debout$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
abline(v = which(debout$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
legend(ndays/div/1.5, max(debout$wetmass) * 0.3, c('repro. buffer', 'food in gut', 'reserve', 'structure'), lty = c(1, 1, 1, 1), col = c("pink", "brown", "grey", "dark green"), bty = 'n', lwd = 1.5)
text(0, max(debout$wetmass) * 1, labels = "embryo", cex = 0.85)
text(which(debout$E_H > E.Hp)[1] / div - which(debout$E_H > E.Hp)[1] / div * .5, max(debout$wetmass) * 1, labels = "immature", cex = 0.85)
text(which(debout$E_H > E.Hp)[1] / div * 1.2, max(debout$wetmass) * 1, labels = "adult", cex = 0.85)

## ----fig3, fig.width=8, fig.height=4,echo=FALSE, results='asis', fig.align='center', fig.cap='Stucture of the standard model of Dynamic Energy Budget theory, from Marques et al. (2018). Boxes represent state variables, arrows represent processes.', fig.pos='!h', message=FALSE, warnings=FALSE, eval=TRUE----
library(png)
library(grid)
img <- readPNG("DEB_scheme.png")
 grid.raster(img)

## ----fig4, fig.width=8, fig.height=4,echo=FALSE, results='asis', fig.align='center', fig.cap='Stoichiometric equations of the standard model of Dynamic Energy Budget theory, taken from Kooijman (2010).', fig.pos='!h', message=FALSE, warnings=FALSE, eval=TRUE----
library(png)
library(grid)
img <- readPNG("DEB_mass_budget.png")
 grid.raster(img)

## ----fig5, warning=FALSE, message=FALSE, echo=FALSE, fig.width=8, fig.height=5, eval=TRUE, fig.cap="DEB predictions of some mass fluxes for the lizard *Eulamprus quoyii* growing under constant *ad libitum* food at 25 &deg;C. In panel a, the red line is the prediction from a general allometric relation for lizards (Andrews and Pough 1985). 'RQ' is the respiratory quotient, the ratio of $CO_2$ to $O_2$ production."----
par(mfrow = c(2,3))
plot(seq(1, ndays) / div, debout$O2ML, type = 'l', xlab = 'Age (days)', ylab = 'ml O2 / hour', col = 'dark green', lwd = 2, ylim = c(0, max(debout$O2ML)), main = "a. oxygen consumption")
points(seq(1, ndays) / div, 10^(0.038*Tbs)*0.013*(debout$wetmass-debout$wetgonad-debout$wetgut)**0.800*(debout$wetmass-debout$wetgonad-debout$wetgut), type = 'l', col = 'red') 
abline(v = which(debout$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
abline(v = which(debout$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
RQ <- max(debout$CO2ML/debout$O2ML)
text(ndays*.75, max(debout$O2ML)/2, paste0('RQ = ', round(RQ, 2)))
plot(seq(1, ndays) / div, debout$CO2ML, type = 'l', xlab = 'Age (days)', ylab = 'ml CO2 / hour', col = 'dark green', lwd = 2, ylim = c(0, max(debout$CO2ML)), main = "b. carbon dioxide production")
abline(v = which(debout$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
abline(v = which(debout$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
plot(seq(1, ndays) / div, debout$GH2OMET * 1000, type = 'l', xlab = 'Age (days)', ylab = 'mg metabolic water / hour', col = 'dark green', lwd = 2, ylim = c(0, max(debout$GH2OMET * 1000)), main = "c. metabolic water production")
abline(v = which(debout$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
abline(v = which(debout$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
plot(seq(1, ndays) / div, debout$GNWASTE * 1000, type = 'l', xlab = 'Age (days)', ylab = 'mg uric acid / hour', col = 'dark green', lwd = 2, ylim = c(0, max(debout$GNWASTE * 1000)), main = "d. nitrogenous waste production")
abline(v = which(debout$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
abline(v = which(debout$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
plot(seq(1, ndays) / div, debout$GDRYFOOD * 1000, type = 'l', xlab = 'Age (days)', ylab = 'mg dry food / hour', col = 'dark green', lwd = 2, ylim = c(0, max(debout$GDRYFOOD * 1000)), main = "e. food intake")
abline(v = which(debout$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
abline(v = which(debout$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')
plot(seq(1, ndays) / div, debout$GFAECES * 1000, type = 'l', xlab = 'Age (days)', ylab = 'mg faeces / hour', col = 'dark green', lwd = 2, ylim = c(0, max(debout$GFAECES * 1000)), main = "f. faeces output")
abline(v = which(debout$E_H > E.Hb)[1] / div, lty = 2, col = 'grey')
abline(v = which(debout$E_H > E.Hp)[1] / div, lty = 2, col = 'grey')

## ---- echo=FALSE, message=FALSE, warnings=FALSE, results='asis'---------------
tabl <- "
*Parameter*                         |	*Parameter Symbol*    |	*Unit*
----------------------------------- | --------------------- | ------------------
max surface-specific assimilation rate | $\\{\\dot{p}_{Am}\\}$ | $J\\:cm^{-2}\\:d^{-1}$
energy conductance                  | $\\dot{v}$            | $d^{-1}$
allocation fraction to soma         | $\\kappa$             | -
volume-specific somatic maintenance cost | $[\\dot{p}_{M}]$ | $J\\:cm^{-3}\\:d^{-1}$
volume-specific cost for structure         | $[E_g]$               | $J\\:cm^{-3}$
maturity at birth                   | $E_{H}^{b}$           | $J$
maturity at puberty                 | $E_{H}^{p}$           | $J$

Table: **Table 1. Minimal required core parameters of the standard DEB model.**
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ---- echo=FALSE, message=FALSE, warnings=FALSE, results='asis'---------------
tabl <- "
*Parameter*                         |	*Parameter Symbol*    |	*Unit*
----------------------------------- | --------------------- | ------------------
Weibull aging acceleration          | $\\ddot{h}_a$         | $d^{-2}$
Arhhenius temperature               | $T_A$                 | Kelvin
Arrhenius lower boundary            | $T_L$                 | Kelvin
Arrhenius upper boundary            | $T_H$                 | Kelvin
Arrhenius temperature at lower boundary | $T_{AL}$          | Kelvin
Arrhenius temperature at lower boundary | $T_{AH}$          | Kelvin

Table: **Table 2. Parameters required to model aging and the thermal response (one may also need estimates of the parameter $s_G$, the Gompertz stress coefficient, to model aging).**
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ---- echo=FALSE, message=FALSE, warnings=FALSE, results='asis'---------------
tabl <- "
*Parameter*                         |	*Parameter Symbol*    |	*Unit*
----------------------------------- | --------------------- | ------------------
shape correction factor             | $\\delta_M$           | -
density of structure                | $d_V$                 | $g \\: cm^{-3}$
density of reserve                  | $d_E$                 | $g \\: cm^{-3}$
molecular weight of reserve         | $w_E$                 | $g \\: mol^{-1}$
chemical potential of reserve       | $\\mu_E$              | $J \\: mol^{-1}$

Table: **Table 3. Parameters required to convert structure, reserve and reproduction buffer (same composition as reserve) into wet mass.**
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ---- echo=FALSE, message=FALSE, warnings=FALSE, results='asis'---------------
tabl <- "
*Parameter*                         |	*Parameter Symbol*    |	*Unit*
----------------------------------- | --------------------- | ------------------
max spec searching rate             | $\\{\\dot{F}_{m}\\}$  | $d^{-2}$
digestion efficiency of food to reserve  | $\\kappa_X$      | $t^{-1}$
faecation efficiency of food to faeces | $\\kappa_P$        | -
maximum specific stomach energy     | $E^S_m$              | $J \\: cm^{-3}$

Table: Table 4. **Parameters required to compute food-density-specific feeding rate.**
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages('R.matlab')
#  library(R.matlab)
#  allStat <- readMat('allStat.mat') # this will take a few minutes
#  save(allStat, file = 'allstat.Rda') # save it as an R data file for faster future loading

## ---- warning=FALSE, message=FALSE, eval=TRUE---------------------------------
library(knitr) # this packages has a function for producing formatted tables.
load('allStat.Rda')

allDEB.species<-unlist(labels(allStat$allStat)) # get all the species names
allDEB.species<-allDEB.species[1:(length(allDEB.species)-2)] # last two elements are not species names
kable(head(allDEB.species))

Nspecies <- length(allStat$allStat)
Nspecies

## ---- warning=FALSE, message=FALSE, eval=TRUE---------------------------------
species <- "Eulamprus.quoyii"
species.slot <- which(allDEB.species == species)
par.names <- unlist(labels(allStat$allStat[[species.slot]]))

for(i in 1:length(par.names)){
 assign(par.names[i], unlist(allStat$allStat[[species.slot]][i]))
}

## ----fig6, warning=FALSE, message=FALSE, fig.height=5, fig.width=8, eval=TRUE----
library(NicheMapR)
species <- "Eulamprus.quoyii" # must be in the AmP collection - see allDEB.species list
ndays <- 365*5 # number days to run the simulation for
div <- 1 # time step divider (1 = days, 24 = hours, etc.) - keep small if using Euler method
Tbs <- rep(25, ndays * div) # °C, body temperature
starvetime <- 0 # length of low food period when simulating starvation
X <- 100 # J/cm2 base food density
Xs <- c(rep(X, ndays * div / 2), rep(0.000005, starvetime), 
        rep(X, ndays * div / 2 - starvetime + 1)) # food density (J/cm2 or J/cm3)
E_sm <- 350 # J/cm3, volume-specific stomach energy
clutchsize <- 5 # -, clutch size

mass.unit <- 'g'
length.unit <- 'mm'
plot <- 1 # plot results?
start.stage <- 1 # stage in life cycle to start (0 = egg, 1 = juvenile, 2 = puberty)

deb <- rundeb(species = species, ndays = ndays, div = div, Tbs = Tbs, clutchsize = clutchsize, Xs = Xs, mass.unit = mass.unit, length.unit = length.unit, start.stage = start.stage, E_sm = E_sm, plot = plot)

## ---- warning=FALSE, message=FALSE, eval=TRUE---------------------------------
longlat <- c(146.77, -19.29) # Townsville, northern Australia
nyear <- 5
micro <- micro_global(loc = longlat, timeinterval = 365, nyear = nyear)

## ---- warning=FALSE, message=FALSE, eval=TRUE---------------------------------
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
p.Xm <- p.Am / kap.X * 10 # redefining p.Xm to a large value relative to p.Am because 
# running a stomach model
z.mult <- 1 # DEB body size scaling parameter

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

## ---- warning=FALSE, message=FALSE, eval=TRUE---------------------------------
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
nocturn <- 0      # nocturnal activity
crepus <- 0       # crepuscular activity
diurn <- 1        # diurnal activity

## ---- warning=FALSE, message=FALSE, eval=TRUE---------------------------------
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

## ---- warning=FALSE, message=FALSE, eval=TRUE---------------------------------
# reproduction parameters
viviparous <- 1 # live bearing (1) or egg laying (0)
clutchsize <- 5 # how many eggs per clutch?
photostart <- 3 # winter solstice is the start of the reproduction cycle
photofinish <- 2 # autumnal equinox is the end of the reproduction cycle

## ---- warning=FALSE, message=FALSE, eval=TRUE---------------------------------
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
                pct_wet=pct_wet,
                z=z*z.mult,
                del_M=del.M,
                p_Xm=p.Xm,
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
                n_X=c(n.CX,n.HX,n.OX,n.NX),
                n_E=c(n.CE,n.HE,n.OE,n.NE),
                n_V=c(n.CV,n.HV,n.OV,n.NV),
                n_P=c(n.CP,n.HP,n.OP,n.NP),
                n_M_nitro=c(n.CN,n.HN,n.ON,n.NN),
                L_b=L.b,
                V_init=V_init,
                E_init=E_init,
                E_H_init=E_H_init,
                stage=stage,
                photostart = photostart,
                photofinish = photofinish)

# retrieve output
environ <- as.data.frame(ecto$environ) # behaviour, Tb and environment
enbal <- as.data.frame(ecto$enbal) # heat balance outputs
masbal <- as.data.frame(ecto$masbal) # mass balance outputs
debout <- as.data.frame(ecto$debout) # DEB model outputs
yearout <- as.data.frame(ecto$yearout) # whole life cycle summary
yearsout <- as.data.frame(ecto$yearsout) # annual summaries

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE---------------
knitr::kable(tail(masbal[, c(1:9)], 12))
knitr::kable(tail(masbal[, c(10:15)], 12))
knitr::kable(tail(masbal[, c(16:19)], 12))

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE---------------
knitr::kable(tail(debout[, c(1:10)], 12))
knitr::kable(tail(debout[, c(11:16)], 12))
knitr::kable(tail(debout[, c(17:21)], 12))

## ---- warning=FALSE, message=FALSE, eval=TRUE, fig.width=8, fig.height=5, fig.cap="DEB model prediction of wet weight through time of the lizard *Eulamprus quoyii* growing under constant *ad libitum* under field conditions in Townsville, Queensland Australia, partitioning total biomass mass into structure, reserve, reproduction buffer and gut contents."----
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
text((which(debout$E_H > E.Hp)[1] - which(debout$E_H > E.Hp)[1] * .5) / 24 ,
     max(debout$WETMASS) * 1, labels = "immature", cex = 0.85)
text(which(debout$E_H > E.Hp)[1] * 1.2 / 24, max(debout$WETMASS) * 1,
     labels = "adult", cex = 0.85)

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE---------------
knitr::kable(yearsout[, c(1:7)])
knitr::kable(yearsout[, c(8:14)])
knitr::kable(yearsout[, c(15:22)])
knitr::kable(yearsout[, c(23:28, 34:36, 42:43)])

## ---- echo=FALSE, results='asis', message=FALSE, warnings=FALSE---------------
knitr::kable(yearout[, c(1:7)])
knitr::kable(yearout[, c(8:14)])
knitr::kable(yearout[, c(15:20)])

