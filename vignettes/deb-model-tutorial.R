## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  eval = TRUE
)

## ----fig1, warning=FALSE, message=FALSE, echo=FALSE, fig.width=8, fig.height=5, eval=TRUE, fig.cap="DEB model prediction of length through time of the lizard *Eulamprus quoyii* growing under constant, *ad libitum* food at 25 &deg;C. Photo by the author."----
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

## ---- echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----------
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

## ---- echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----------
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

## ---- echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----------
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

## ---- echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----------
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

## ---- eval=FALSE---------------------------------------------------------
#  install.packages('R.matlab')
#  library(R.matlab)
#  allStat <- readMat('allStat.mat') # this will take a few minutes
#  save(allStat, file = 'allstat.Rda') # save it as an R data file for faster future loading

## ---- warning=FALSE, message=FALSE, eval=TRUE----------------------------
library(knitr) # this packages has a function for producing formatted tables.
load('allStat.Rda')

allDEB.species<-unlist(labels(allStat$allStat)) # get all the species names
allDEB.species<-allDEB.species[1:(length(allDEB.species)-2)] # last two elements are not species names
kable(head(allDEB.species))

Nspecies <- length(allStat$allStat)
Nspecies

## ---- warning=FALSE, message=FALSE, eval=TRUE----------------------------
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

deb <- rundeb(species = species, ndays = ndays, div = div, Tbs = Tbs, clutchsize = clutchsize, 
              Xs = Xs, mass.unit = mass.unit, length.unit = length.unit, 
              start.stage = start.stage, E_sm = E_sm, plot = plot)

