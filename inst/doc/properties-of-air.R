## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ----table1, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Physical Quantity*       |	*Quantity Symbol* |	*Unit*   |	*Unit Symbol*
------------------------- | ----------------- | -------- | --------------
length                    | *l*               | metre    | m
mass                      | *m*               | kilogram | kg
time                      | *t*               | second   | s
thermodynamic temperature | *T*               | kelvin   | K
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table2, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Physical Quantity*  | *Quantity Symbol*   | *Unit*   | *Unit Symbol*   | *Unit Definition*
-------------------- | ------------------- | -------- | --------------- | -----------------
force                | *F*                 | newton   | N               | $J$ $m^{-1}$
energy               | *E*                 | joule    | J               | $N$ $m$
pressure             | *p*                 | pascal   | Pa              | $N$ $m^{-2}$, $J$ $m^{-3}$
power                | *P*                 | watt     | W               | $J$ $s^{-1}$
volume               | *t*                 | cubic metre | $m^3$        | $m^3$
Celsius temperature$^*$ | $t$              | degree Celsius   | &deg;C   | $K$
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table3, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Multiple* | *Prefix*	| *Symbol* |	*Multiple* | *Prefix*	 | *Symbol*
---------- | -------- | -------- | ----------- | --------- | -------
$10^{-1}$	 | deci	    | *d*      | $10^{1}$	   | deca	     | da
$10^{-2}$	 | centi    | *c*      | $10^{2}$	   | hecto     | h
$10^{-3}$	 | milli    | *m*      | $10^{3}$	   | kilo	     | k
$10^{-6}$	 | micro    | $\\mu$      | $10^{6}$	   | mega	     | M
$10^{-9}$	 | nano	    | *n*      | $10^{9}$	   | giga	     | G
$10^{-12}$ | pico	    | *p*      | $10^{12}$   | tera	     | T
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table4, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
**Force**   
&nbsp;&nbsp;&nbsp;&nbsp;dyne = $10^5$ N$^*$  
&nbsp;&nbsp;&nbsp;&nbsp;kilogram force = $9.80665$ N$^*$  
&nbsp;&nbsp;&nbsp;&nbsp;ounce force = $2.7801385$ x $10^{-1}$ N  
&nbsp;&nbsp;&nbsp;&nbsp;pound force = $4.4482216152605$ N$^*$  
   
**Energy**   
&nbsp;&nbsp;&nbsp;&nbsp;British thermal unit (thermochemical) = $1.054350$ x $10^{3}$ J  
&nbsp;&nbsp;&nbsp;&nbsp;calorie (thermochemical) = $4.184$ J  
&nbsp;&nbsp;&nbsp;&nbsp;what hour = $3.60$ x $10^{3}$ J$^*$  
   
**Pressure**   
&nbsp;&nbsp;&nbsp;&nbsp;bar = $10^{5}$ $Pa^*$  
&nbsp;&nbsp;&nbsp;&nbsp;centimetre of water (4 &deg;C) = $9.80638$ x $10^{1}$ Pa   
&nbsp;&nbsp;&nbsp;&nbsp;millimetre of mercury (0 &deg;C) = $1.333224$ x $10^{2}$ Pa   
&nbsp;&nbsp;&nbsp;&nbsp;pound per square inch = $6.894757$ x $10^{3}$ Pa   
&nbsp;&nbsp;&nbsp;&nbsp;standard atmosphere = $1.01325$ x $10^{5}$ Pa   
&nbsp;&nbsp;&nbsp;&nbsp;torr = $1.333224$ x $10^{2}$ Pa   
   
**Power**   
&nbsp;&nbsp;&nbsp;&nbsp;British thermal unit (thermochemical) = $1.054350$ x $10^{3}$ W   
&nbsp;&nbsp;&nbsp;&nbsp;calorie (thermochemical) per second = $4.184$ W$^*$   
&nbsp;&nbsp;&nbsp;&nbsp;calorie (thermochemical) per minute = $6.973333$ x $10^{-2}$ W   
&nbsp;&nbsp;&nbsp;&nbsp;kilocalorie (thermochemical) per hour = $8.604207$ x $10^{-1}$ W   
   
**Volume**   
&nbsp;&nbsp;&nbsp;&nbsp;cubic foot = $2.8316592$ m$^{3*}$   
&nbsp;&nbsp;&nbsp;&nbsp;gallon = $3.785411784$ x $10^{-3}$ m$^{-3*}$    
&nbsp;&nbsp;&nbsp;&nbsp;litre = $10^{-3}$ m$^{-3*}$   
&nbsp;&nbsp;&nbsp;&nbsp;pint (U.S. liquid) = $4.73176473$ x $10^{-4}$  m$^{3*}$   
&nbsp;&nbsp;&nbsp;&nbsp;quart (U.S. liquid) = $9.4635295$ x $10^{-4}$  m$^{3*}$  
   
**Length**   
&nbsp;&nbsp;&nbsp;&nbsp;foot = $3.048$ x $10^{-1}$ m$^{*}$   
&nbsp;&nbsp;&nbsp;&nbsp;inch = $2.54$ x $10^{-2}$ m$^{*}$   
&nbsp;&nbsp;&nbsp;&nbsp;mile (U.S. statute) = $1.609344$ x $10^{3}$ m$^{*}$   
&nbsp;&nbsp;&nbsp;&nbsp;mile (U.S. nautical) = $1.852$ x $10^{3}$ m$^{*}$   
&nbsp;&nbsp;&nbsp;&nbsp;yard = $9.144$ x $10^{-1}$ m$^{*}$   
   
**Mass**
&nbsp;&nbsp;&nbsp;&nbsp;ounce = $2.8349523125$ x $10^{-2}$ kg$^{*}$   
&nbsp;&nbsp;&nbsp;&nbsp;pound = $4.5359237$ x $10^{-1}$ kg$^{*}$   
&nbsp;&nbsp;&nbsp;&nbsp;ton (short, 2000 pound) = $9.0718474$ x $10^{2}$ kg$^{*}$   
   
$^{*}$exact, by definition
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table5, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Quantity*                                           | *Symbol*	| *Value*  |	*Units*
---------------------------------------------------- | -------- | -------- | -----------
Density of dry air at 25 &deg;C and 101325 Pa$^{*}$	 | $\\rho$	    | $1.184$    | kg m$^3$	   
Diffusivity of water vapor in air at 25 &deg;C and 101325 Pa$^{*}$	 | *D*	    | $2.614$ x $10^{-5}$    | m$^2$ s$^{-1}$	   
Dynamic viscosity of dry air at 25 &deg;C and 101325 Pa$^{*}$	 | $\\mu$	    | $1.806$ x $10^{-5}$    | kg m$^{-1}$ s$^{-1}$	   
Gas constant for ideal gas	 | *R*	    | $8.31434$    | J mol$^{-1}$ K$^{-1}$   
Gas constant for water vapor	 | *R*$_v$	    | $4.6150$ x $10^{2}$    | Pa m$^3$ kg$^{-1}$ K$^{-1}$	   
Gravitational constant	 | *R*	    | $6.6732$ x $10^{-11}$    | N m$^2$ kg$^{-2}$	   
Kinematic viscosity of dry air at 25 &deg;C and 101325 Pa$^{*}$	 | $\\nu$	    | $1.525$ x $10^{-5}$    | m$^3$ s$^{-1}$	
Latent heat of vaporization of water at 25 &deg;C	 | *L*	    | $2.442$ x $10^{6}$    | J kg$^{-1}$	   
Molecular weight of dry air	 | *M*$_a$	    | $2.8966$ x $10^{-2}$    | kg mol$^{-1}$	   
Molecular weight of water | *M*$_v$	    | $1.8016$ x $10^{-2}$    | kg mol$^{-1}$	 
Specific heat of dry air	 | *c*$_p$	    | $1.00484$ x $10^{3}$    | K kg$^{-1}$ K$^{-1}$	   
Standard pressure	 | *P*$_0$ 	    | $1.01325$ x $10^{5}$    | Pa	   
Standard temperature	 | *T*$_0$ 	    | $2.7315$ x $10^{2}$    | K	   
Stefan-Boltzmann constant	 | $\\sigma$	    | $5.67032$ x $10^{-8}$    | W m$^{-2}$ K$^{-4}$	   
Thermal conductivity of dry air at 25 &deg;C and 101325 Pa$^{*}$	 | *k*	    | $2.601$ x $10^{-2}$    | W m$^{-1}$ K$^{-1}$	   
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table6, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "

*Compound* | *10&deg;C*	| *15&deg;C*  |	 *20&deg;C* | *25&deg;C*	| *30&deg;C*	| *35&deg;C*	
--- | --- | --- | --- | --- | --- | --- 
NH$_4$Cl|79.0 | 79.5 | 79.5 | 78.0 | 77.5 | -
NH$_4$HPO$_4$ (monophosphate)|- | - | 93.0 | 93.0 | 92.0 | -
NH$_4$NO$_3$|75.0 | 70.0 | 65.5 | 62.5 | 59.5 | 55.0
NH$_4$NO$_3$ + NaNO$_3$|58.0 | 55.0 | 52.0 | 50.0 | 47.0 | 44.5
NH$_4$NO$_3$ + AgNO$_3$|70.5 | 68.0 | 65.1 | 61.5 | 58.0 | 55.0
(NH$_4$)$_2$ SO$_4$|80.5 | 81.0 | 80.5 | 80.0 | 80.0 | 79.5
CaCl$_2$&middot;6H$_2$O$^2$|38.0 | 35.0 | 32.5 | 29.5 | - | -
CaHPO$_4$&middot;2H$_2$O$^3$|97.5 | 99.5 | 95.0 | 97.0 | 95.0 | -
CaH$_4$(PO$_4$)&middot;H$_2$O$^4$|98.0 | 99.0 | 94.0 | 96.0 | 93.5 | -
Ca(NO$_3$)$_2$&middot;4H$_2$O$^5$|- | 56.0 | 55.5 | 50.5 | 47.0 | -
C0Cl$_2$$^6$|- | 72.5 (18 &deg;C) | 67.0 | - | 62.0 | -
Cr$_2$O$_3$$^7$|- | 45.5 | (39.0) | - | 44.5 | -
Glucose|57 (12 &deg;C) | - | 55.0 | 55.0 | - | 53.0
LiCl&middot;H$_2$O$^8$|13.5 | 13.0 | 12.5 | 12.0 | 11.5 | 11.5
Pb(NO$_3$)$_2$$^9$|98.0 | 97.0 | 97.0 | 95.5 | 95.0 | 94.5
Pb(NO$_3$)$_2$+NH$_4$NO$_3$|66.5 | 62.0 | 58.0 | 55.0 | 52.5 | 49.5
MgCl$_2$&middot;6H$_2$O|34.0 | 34.0 | 33.0 | 32.5 | 32.5 | 32.5
Mg(NO$_3$)$_2$&middot;6H$_2$O|58.0 | 56.0 | 55.0 | 53.0 | 52.0 | 50.5
P$_2$O$_5$|0.0 | 0.0 | 0.0 | 0.0 | 0.0 | -
KAc|(21.0) | - | 20.0 | 22.5 | 22.0 | -
KBr|86.0 | - | 84.0 | 80.0 | 82.0 | -
K$_2$CO$_3$&middot;2H$_2$O$^{10}$|47.0 | 44.0 | 44.0 | 43.0 | 43.0 | -
KCl|88.0 | 86.5 | 85.0 | 85.0 | 84.5 | 83.0
K$_2$CrO$_4$|- | - | 88.0 | - | 86.5 | -
K$_2$Cr$_2$O$_7$|- | - | 98.0 | 98.0 | 97.5 | 96.5
KH$_2$PO$_4$|98.0 | 99.0 | 96.5 | 96.0 | 93.5 | -
KF&middot;2H$_2$O$^{11}$|- | - | - | 30.5 | 27.4 | -
KNO$_3$$^{12}$|96.0 | 95.5 | 93.5 | 92.5 | 91.0 | 89.5
KNO$_2$|- | - | 48.5 | - | 47.0 | -
K$_2$SO$_4$|98.5 | 99.0 | 98.0 | 97.5 | 96.5 | 96.0
KCNS$^{13}$|52.0 | 50.0 | 47.0 | 46.5 | 43.5 | 41.5
K tartrate$^{14}$|- | 75.0 | 75.0 | 75.0 | 74.0 | -
KNa tartrate$^{15}$|87.5 | - | 87.0 | 87.0 | 87.0 | -
Pyrocatechol|- | 99.0 | 95.5 | 93.5 | 92.5 | 90.5
Resorinol|- | 95.0 | 87.0 | 85.0 | 82.5 | 79.5
AgNO$_3$|88.0 | 86.0 | 84.0 | 82.0 | 80.0 | 78.0
AgNO$_3$+Pb(NO$_3$)$_2$|84.0 | 83.0 | 80.5 | 78.5 | 76.5 | 75.0
NaAc&middot;3H$_2$O|- | - | 75.0 | 73.0 | 71.0 | -
NaBr&middot;2H$_2$O$^{16}$|63.0 | 61.0 | 59.0 | 57.5 | 56.0 | 54.5
Na$_2$CO$_3$&middot;10H$_2$O$^{17}$|(99.0) | - | 92.0 | 87.0 | 87.0 | -
NaCl|76.5 | 76.0 | 76.0 | 75.5 | 75.5 | 75.5
NaCl+KCl|- | - | 70.0 | 71.5 | 71.0 | -
Na$_2$CrO$_4$&middot;4H$_2$O|- | - | - | - | 64.5 | -
Na$_2$Cr$_2$O$_7$&middot;H$_2$O$^{18}$|60.0 | 56.5 | 54.5 | 53.0 | 52.5 | 51.0
NaOH|(5.5) | - | 5.5 | 7.0 | (4.0) | -
NaI|- | - | 38.0 | 38.0 | 36.0 | -
NaNO$_3$|77.5 | 76.5 | 76.0 | 74.0 | 72.5 | 71.0
NaNO$_2$|- | - | 65.5 | 64.0 | 63.0 | -
Na$_2$SO$_4$&middot;10H$_2$O|- | - | (93.0) | (93.0) | - | (87.0)
Na tartrate$^{19}$|- | 94.0 | 92.0 | 92.0 | 92.0 | -
Urea|81.5 | 80.0 | 80.0 | 76.0 | 73.3 | -
ZnCl$_2$&middot;1.5H$_2$O$^{20}$|(10.0) | (10.0) | 10.0 | - | (10.0) | -
ZnSO$_4$&middot;7H$_2$O$^{21}$|- | - | - | 88.5 | - | -
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table7, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "

*Compound* | *Relative Humidity*	| *Compound*  |	 *Relative Humidity*	
--- | --- | --- | --- 
BaI$_2$&middot;2H$_2$O|44.0 | NaHSO$_4$ | 52.0
CaSO$_4$&middot;2H$_2$O|98.0 | NaBrO$_3$ | 92.0
NH$_4$&middot;4H$_2$O|65.0 | NaNO$_3$+KNO$_3$| 71.0
NH$_4$&middot;7H$_2$O|90.0 | Na$_2$HPO$_4$&middot;7H$_2$O | 95.0
NH$_4$&middot;4H$_2$O|57.5 | Na$_2$SO$_3$$^*$ | 95.0
NH$_4$&middot;6H$_2$O|55.5 | NaCNS | 35.5
H$_3$PO$_4$|9.0 | Na$_2$S$_2$O$_3$ | 78.0
KHSO$_4$|86.0 | 2n(NO$_3$)$_2$ | 42.0
K$_2$HPO$_4$|44.5 | - | -
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table8, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "

*Compound* | *Relative Humidity*	| *Compound*  |	 *Relative Humidity*	
--- | --- | --- | --- 
NH$_4$Br$_2$|75.0 | MgSO$_4$&middot;7H$_2$O | 89.0
BaBr$_2$|74.5 | Mg(CNS)$_2$ | 47.5
BaCl$_2$|90.0 | MnBr$_2$&middot;6H$_2$O | 34.5
BaI$_2$&middot;2H$_2$O|43.0 | MnCl$_2$&middot;4H$_2$O | 56.0
Ba(CNS)$_2$&middot;2H$_2$O|54.5 | NiBr$_2$&middot;3H$_2$O | 27.0
CaAc$_2$&middot;H$_2$O|17.0 | NiCl$_2$&middot;6H$_2$O | 53.0
CaAc$_2$&middot;H$_2$O + sucrose|13.0 | H$_3$PO$_4$ | 9.0
CaBr$_2$|16.5 | KBr + sucrose | 63.0
CaI$_2$&middot;6H$_2$O|11.5 | KBr + Urea | 51.0
Ca methane sulphonate|72.5 | KClO$_3$ | 98.0
Ca(MnO$_4$)$_2$&middot;4H$_2$O|37.5 | KCl + KClO$_3$ | 85.0
Ca(CNS)$_2$&middot;3H$_2$O|17.5 | KHCO$_4$ | 21.5
CaZnCl$_2$|20.0 | KOH | 8.0
CeCl$_3$ (tech.)|45.5 | K$_4$P$_2$O$_7$&middot;3H$_2$O | 49.5
CrCl$_3$|42.5 | K$_2$S$_2$O$_3$ | 66.0
CoBr$_2$|41.5 | NaAc + sucrose | 37.6
Co(NO$_3$)$_2$&middot;6H$_2$O|(49.0) | NaCl + sucrose | 63.0
CuBr$_2$|62.5 | NaCl + NaNO$_3$ | 69.0
CuCl$_2$&middot;2H$_2$O|67.0 | NaCl + Na$_2$SO$_4$&middot;7H$_2$O | 74.0
Cu(NO$_3$)$_2$&middot;6H$_2$O|(35.0) | NaHSO$_4$ | 81.0
Ethanolamine sulphate|12.5 | NaH$_2$PO$_4$ | 70.0
FeBr$_2$&middot;6H$_2$O|39.0 | Na methane sulphonate | 61.5
FeCl$_2$&middot;4H$_2$O|60.0 | SrBr$_2$&middot;6H$_2$O | 58.5
LiBr&middot;2H$_2$O|7.0 | SrI$_2$&middot;6H$_2$O | 33.0
LiI&middot;3H$_2$O|18.0 | Sr(CNS)$_2$&middot;3H$_2$O | 31.5
LiNo$_3$&middot;3H$_2$O|47.0 | Sucrose | 85.0
LiCNS$^2$|12.5 | Surcose + urea | 52.0
MgBr$_2$&middot;6H$_2$O|31.5 | ZnBr$_2$ | 8.5
MgI$_2$|27.0 | ZbI$_2$$^2$ | 20.0
Mg(ClO$_4$)$_2$&middot;6H$_2$O|(41.0) | Zn(MnO$_4$)$_2$&middot;6H$_2$O | 51.0
Mg silicoflouride|91.5 | Zn(CNS)$_2$ | 80.5
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table9, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "

*Compound* | *Mass*, g m$^{-3}$ (STPD)	| *Relative Humidity*  	
--- | --- | --- | --- 
P$_4$O$_5$|<2 x 10$^{-3}$ | <0.1
Mg(ClO$_4$)$_2$|<5 x 10$^{-4}$ | <0.1
Mg(ClO$_4$)$_2$&middot;3H$_2$O|<2 x 10$^{-5}$ | <0.1
KOH (fused)|0.002 | <0.1
Al$_2$O$_3$|0.003 | <0.1
H$_2$SO$_4$|0.003 | <0.1
MgO|0.008 | <0.1
NaOH (fused)|0.16 | 0.6
CaBr$_2$|0.2 | 0.8
CaO|0.2 | 0.8
CaCl$_2$|0.14-0.25 | 0.6-1.0
H$_2$SO$_4$ (95%)|0.3 | 1.2
CaCl$_2$ (fused)|0.36 | 1.4
ZnCl$_2$|0.8 | 3.2
ZnBr$_2$|1.1 | 4.4
CuSO$_4$|1.4 | 5.6
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table10, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
library(NicheMapR)
knitr::kable(PropAirTable10) # output the table in a format good for HTML/PDF/docx conversion

## ----table11, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
library(NicheMapR)
knitr::kable(PropAirTable11) # output the table in a format good for HTML/PDF/docx conversion

## ----table12, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
library(NicheMapR)
knitr::kable(PropAirTable12) # output the table in a format good for HTML/PDF/docx conversion

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylim=c(0,1000)
ylab = expression(paste("BLACK BODY EMITTANCE ( ", italic(phi), ") W m"^{-2}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs)$bbemit,type='n',ylim=ylim,ylab=ylab, xlab=xlab,cex.lab=1.5,cex.axis=1.5,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
minor.ticks<-function (nx = 2, ny = 2, tick.ratio = 0.5)
{
    ax <- function(w, n, tick.ratio) {
        range <- par("usr")[if (w == "x")
            1:2
        else 3:4]
        tick.pos <- if (w == "x")
            par("xaxp")
        else par("yaxp")
        distance.between.minor <- (tick.pos[2] - tick.pos[1])/tick.pos[3]/n
        possible.minors <- tick.pos[1] - (0:100) * distance.between.minor
        low.candidates <- possible.minors >= range[1]
        low.minor <- if (any(low.candidates)) 
            min(possible.minors[low.candidates])
        else tick.pos[1]
        possible.minors <- tick.pos[2] + (0:100) * distance.between.minor
        hi.candidates <- possible.minors <= range[2]
        hi.minor <- if (any(hi.candidates)) 
            max(possible.minors[hi.candidates])
        else tick.pos[2]
        axis(if (w == "x") 
            1
        else 2, seq(low.minor, hi.minor, by = distance.between.minor), 
            labels = FALSE, tcl = par("tcl") * tick.ratio)
    }
    if (nx > 1)
        xticks=ax("x", ny, tick.ratio = tick.ratio)
    if (ny > 1)
        yticks=ax("y", ny, tick.ratio = tick.ratio)
    return(list(xticks=xticks, yticks=yticks))
}
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs)$bbemit,type='l',ylim=ylim,ylab=ylab, xlab=xlab,cex.lab=1.5,cex.axis=1.5,lwd=2)
plotrix::boxed.labels(10, 550, expression(paste(phi," = ",5.67032," x 10"^{-8},"(",italic(t)," + 273.15)"^{4})),cex=1.5,bg="white",border=NA)

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("DENSITY OF DRY AIR ( ", italic(rho), ") kg m"^{-3}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs,bp = 85000)$densty,type='n',ylim=c(1,1.5),ylab=ylab, xlab=xlab,cex.lab=1.5,cex.axis=1.5,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs,bp = 85000)$densty,type='l',ylim=c(1,1.5),ylab=ylab, xlab=xlab,cex.lab=1.5,cex.axis=1.5,lwd=2)
points(tairs,DRYAIR(db=tairs,bp = 100000)$densty,type='l',lwd=2)
plotrix::boxed.labels(30, 1.4, expression(paste(italic(rho)," = ",frac(italic(P), paste(287.04,"(",italic(t)," + 273.15)",sep="")))),cex=1.5,bg="white",border=NA)
plotrix::boxed.labels(10,1.3,"100 000 Pa",cex=1.5,bg="white",border=NA)
plotrix::boxed.labels(0,1.14,"85 000 Pa",cex=1.5,bg="white",border=NA)

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(2,1)) # set up for 2 plots in 1 columns
dps=seq(0,50,5)
par(oma = c(5,5,2,2) + 0.1)
par(mar = c(0,0,1,0) + 0.1)
par(xaxs="i")
par(yaxs="i")
xlab = expression(paste("SATURATION VAPOR PRESSURE AT ",italic(t)[{d}]," OVER WATER (", italic(e)[{d}],") Pa"))
ylab = expression(paste("DEW POINT TEMPERATURE (",italic(t)[{d}],") ",degree*C))
plot(WETAIR(db=dps,rh=100)$esat,dps,type='n',xlim=c(0,14000),xaxt='n',yaxt='n',ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(WETAIR(db=dps,rh=100)$esat,dps,type='l',xlim=c(0,14000),xaxt='n',yaxt='n',ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,lwd=2)
axis(side=3, at=seq(0,14,2)*1000)
mtext(xlab, side=3, line=2,cex=1.25)
axis(side=2, at=seq(0,50,10))
x=0
y=10
plotrix::boxed.labels(750+x,25+y, expression(paste(italic(alpha)," = ",7.5)),cex=.8,bg="white",border=NA) 
plotrix::boxed.labels(900+x,22+y, expression(paste(italic(beta)," = ",237.3)),cex=.8,bg="white",border=NA) 
plotrix::boxed.labels(1800+x,23.5+y, "}",cex=1.5,bg="white",border=NA)
plotrix::boxed.labels(3400+x,23.5+y, "OVER WATER",cex=.8,bg="white",border=NA)
dps=seq(-50,0,5)
xlab = expression(paste("SATURATION VAPOR PRESSURE AT ",italic(t)[{d}]," OVER ICE (", italic(e)[{d}],") Pa"))
ylab = expression(paste("DEW POINT TEMPERATURE (  ",italic(t)[{d}]," ) ",degree*C))
plot(WETAIR(db=dps,rh=100)$esat,dps,type='n',xlim=c(0,700),xaxt='n',yaxt='n',ylab="", xlab="",cex.lab=1.2,cex.axis=1.5,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(WETAIR(db=dps,rh=100)$esat,dps,type='l',xlim=c(0,700),xaxt='n',yaxt='n',ylab="", xlab="",cex.lab=1.5,cex.axis=1.5,lwd=2)
axis(side=1, at=seq(0,7,1)*100)
mtext(xlab, side=1, line=2.3, cex = 1.25)
axis(side=2, at=seq(-50,0,10))
par(cex=1.5)
title(ylab = ylab, outer = TRUE,line=2)
plotrix::boxed.labels(400,-35, expression(paste(italic(e)[d]^{"*"}," = 100 x 10"^{bgroup("[",0.7858+frac(italic(t)[d]*alpha,italic(t)[d]+beta),"]")})),cex=1.5,bg="white",border=NA)
x=-170
y=15
plotrix::boxed.labels(300+x,-20+y, expression(paste(italic(alpha)," = ",9.5)),cex=.8,bg="white",border=NA) 
plotrix::boxed.labels(310+x,-23.5+y, expression(paste(italic(beta)," = ",265.5)),cex=.8,bg="white",border=NA) 
plotrix::boxed.labels(360+x,-21.5+y, "}",cex=1.5,bg="white",border=NA)
plotrix::boxed.labels(420+x,-21.5+y, "OVER ICE",cex=.8,bg="white",border=NA)

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("DIFFUSIVITY OF WATER VAPOR IN AIR (  ", italic(D), ") m"^{2},"s"^{-1}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs,bp = 70000)$difvpr,type='n',ylim=c(2E-5,3E-5),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs,bp = 70000)$difvpr,type='l',ylim=c(2E-5,3E-5),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
points(tairs,DRYAIR(db=tairs,bp = 85000)$difvpr,type='l',lwd=2)
points(tairs,DRYAIR(db=tairs,bp = 100000)$difvpr,type='l',lwd=2)
rect(xleft = 26.1,xright = 48, ytop = 2.56E-5, ybottom = 2.04E-5,col="white",border = NA)
plotrix::boxed.labels(17,2.35E-5,"100 000 Pa",cex=1.25,bg="white",border=NA)
plotrix::boxed.labels(4,2.55E-5,"85 000 Pa",cex=1.25,bg="white",border=NA)
plotrix::boxed.labels(-9,2.85E-5,"70 000 Pa",cex=1.25,bg="white",border=NA)
plotrix::boxed.labels(37, 2.5E-5, expression(paste(italic(D)," = ",italic(D)[0],(frac(italic(T),italic(T)[0]))^italic(n),frac(italic(rho),italic(rho)[0]))),cex=1.15,bg="white",border=NA) 
plotrix::boxed.labels(37, 2.40E-5, expression(paste(italic(D)[0]," = ",2.26," x ",10^{-5},sep="")),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(37, 2.30E-5, expression(paste(italic(T)," = ",italic(t)+273.15)),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(37, 2.2E-5, expression(paste(italic(n)," = ",1.81)),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(37, 2.1E-5, expression(paste(italic(rho)[0]," = ",2.26," x ",10^{5})),cex=1.05,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("DYNAMIC VISCOSITY OF AIR (  ", italic(mu), ") kg m"^{-1},"s"^{-1}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs)$visdyn,type='n',ylim=c(.4E-5,2.4E-5),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs)$visdyn,type='l',ylim=c(.4E-5,2.4E-5),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
rect(xleft = -4,xright = 24, ytop = 1.65E-5, ybottom = 0.7E-5,col="white",border = NA)
plotrix::boxed.labels(10, 1.5E-5, expression(paste(italic(mu)," = ",italic(mu)[0],bgroup("[",frac(italic(T)[0]+C,italic(T)+C)*bgroup("(",frac(italic(T),italic(T)[0]),")")^italic(1.5),"]"))),cex=1.15,bg="white",border=NA) 
plotrix::boxed.labels(10, 1.2E-5, expression(paste(italic(mu)[0]," = ",1.8325," x ",10^{-5},sep="")),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(10, 1.04E-5, expression(paste(italic(T)[0]," = ",296.16)),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(10, 0.94E-5, expression(paste(italic(C)," = ",120)),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(10, 0.81E-5, expression(paste(italic(rho)[0]," = ",italic(t)+273.15)),cex=1.05,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("GROUP OF VARIABLES IN GRASHOF NUMBER (  ", italic(gamma), ") m"^{-3},"K"^{-1}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs)$ggroup,type='n',ylim=c(9E5,10E6),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs)$ggroup,type='l',ylim=c(9E5,10E6),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
plotrix::boxed.labels(10, 8E6, expression(paste(italic(gamma)," = ",frac(0.0980618*italic(beta),italic(nu)^{2}))),cex=1.15,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
par(xaxs="i")
par(yaxs="i")
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
ylab = expression(paste("KINEMATIC VISCOSITY OF AIR AT 10 000 Pa (  ", italic(nu), ") m"^{2},"K"^{-1}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs)$viskin,type='n',ylim=c(0.2E-5,2.2E-5),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs)$viskin,type='l',ylim=c(0.2E-5,2.2E-5),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
plotrix::boxed.labels(-10, 2E-5, expression(paste(italic(nu)," = ",frac(italic(mu),italic(rho)))),cex=1.5,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
par(xaxs="i")
par(yaxs="i")
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
ylab = expression(paste("LATENT HEAT OF VAPORIZATION OF WATER (  ", italic(L), ") J kg"^{-1}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs)$htovpr,type='n',ylim=c(2.35E6,2.6E6),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs)$htovpr,type='l',ylim=c(2.35E6,2.6E6),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
plotrix::boxed.labels(20, 2.55E6, expression(paste(italic(L)," = ",2.5012," x 10"^6-2378.7*italic(t))),cex=1.25,bg="white",border=NA) 
plotrix::boxed.labels(20, 2.53E6, expression(-20<paste(italic(L)<60)),cex=1.25,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("MIXING RATIO OVER WATER AT 100 000 Pa (  ", italic(r)[w], ") kg kg"^{-1}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,WETAIR(db=tairs, rh=100)$rw,type='n',ylim=c(.9E-3,10E-2),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,WETAIR(db=tairs, rh=100)$rw,type='l',ylim=c(.9E-3,10E-2),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
points(tairs,WETAIR(db=tairs, rh=50)$rw,type='l',cex.lab=1.25,cex.axis=1.25,lwd=2)
points(tairs,WETAIR(db=tairs, rh=25)$rw,type='l',cex.lab=1.25,cex.axis=1.25,lwd=2)
plotrix::boxed.labels(0, 0.08, expression(paste(italic(r)[w]," = ",frac(0.62570*italic(e),italic(p)-1.0060*italic(e)))),cex=1.15,bg="white",border=NA) 
plotrix::boxed.labels(28, 0.035, expression(paste(100,"% ",italic(rh))),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(35, 0.025, expression(paste(50,"% ",italic(rh))),cex=1.05,bg="white",border=NA)
plotrix::boxed.labels(39, 0.016, expression(paste(25,"% ",italic(rh))),cex=1.05,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("SPECIFIC HEAT OF AIR AT 100 000 Pa (  ", italic(c)[p], ") J kg"^{-1}, "K"^{-1}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,WETAIR(db=tairs, rh=100)$cp,type='n',ylim=c(1E3,1.1E3),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,WETAIR(db=tairs, rh=100)$cp,type='l',ylim=c(1E3,1.1E3),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
points(tairs,WETAIR(db=tairs, rh=50)$cp,type='l',cex.lab=1.25,cex.axis=1.25,lwd=2)
points(tairs,WETAIR(db=tairs, rh=25)$cp,type='l',cex.lab=1.25,cex.axis=1.25,lwd=2)
plotrix::boxed.labels(0, 1080, expression(paste(italic(c)[p]," = ",frac(1004.84+(1864.40*italic(r)[w]),1+italic(r)[w]))),cex=1.15,bg="white",border=NA) 
plotrix::boxed.labels(28, 1040, expression(paste(100,"% ",italic(rh))),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(35, 1025, expression(paste(50,"% ",italic(rh))),cex=1.05,bg="white",border=NA)
plotrix::boxed.labels(39, 1018, expression(paste(25,"% ",italic(rh))),cex=1.05,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
par(xaxs="i")
par(yaxs="i")
alts=seq(-500,3000,250)
par(mar = c(5,5,4,2) + 0.1)
ylab = expression(paste("STANDARD ATMOSPHERIC PRESSURE (  ", italic(p), ") Pa"))
xlab = expression(paste("ALTITUDE (",italic(Z),") m"))
plot(alts,DRYAIR(db=20, alt=alts)$patmos,type='n',ylim=c(0.7E5,1.2E5),ylab=ylab,xlab=xlab,cex.lab=1.1,cex.axis=1.15,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(alts,DRYAIR(db=20, alt=alts)$patmos,type='l',ylim=c(0.7E5,1.2E5),ylab=ylab,xlab=xlab,cex.lab=1.1,cex.axis=1.15,lwd=2)
plotrix::boxed.labels(1500, 1.1E5, expression(paste(italic(p)," = ",101325,bgroup("[",1-(2.2569*10^{-5}*Z),"]")^5.2553)),cex=1.25,bg="white",border=NA) 
plotrix::boxed.labels(1500, 1.05E5, expression(-1000<paste(italic(Z)<20000)),cex=1.25,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("TEMPERATURE (", italic(t)[F],")",degree*F))
xlab = expression(paste("TEMPERATURE (", italic(t)[C],")",degree*C))
plot(tairs,((9/5)*tairs)+32,type='n',ylim=c(10,110),ylab=ylab,xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)
tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,((9/5)*tairs)+32,type='l',ylim=c(10,110),ylab=ylab,xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
tairs=seq(-40,-20,5)
Hmisc::subplot(fun=plot(tairs,((9/5)*tairs)+32,type='l',xlim=c(-40,-20),ylim=c(-40,0),ylab="",xlab="",cex.lab=.9,cex.axis=0.98,lwd=2),x=c(30,50),y=c(20,60)) 
plotrix::boxed.labels(-5, 95, expression(paste(italic(t)[f]," = ",bgroup("[",bgroup("(",frac(9,5),")")*italic(t)[c],"]")+32)),cex=1.25,bg="white",border=NA) 
plotrix::boxed.labels(-5, 75, expression(paste(italic(t)[c]," = ",bgroup("(",frac(5,9),")")*(italic(t)[f]-32))),cex=1.25,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
par(xaxs="i")
par(yaxs="i")
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
ylab = expression(paste("TEMPERATURE COEFFICIENT OF VOLUME EXPANSION FOR DRY AIR (  ", italic(beta), ") K"^{-1}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs)$tcoeff,type='n',ylim=c(3E-3,4E-3),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)

tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs)$tcoeff,type='l',ylim=c(3E-3,4E-3),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
plotrix::boxed.labels(15,3.8E-3, expression(paste(italic(beta)," = ",frac(1,italic(t)+273.15))),cex=1.25,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
par(xaxs="i")
par(yaxs="i")
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
ylab = expression(paste("THERMAL CONDUCTIVITY OF AIR (  ", italic(k), ") W m"^{-1}, " K"^{-1}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs)$thcond,type='n',ylim=c(1E-3,3E-2),ylab=ylab, xlab=xlab,cex.lab=1.2,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)

tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs)$thcond,type='l',ylim=c(1E-3,3E-2),ylab=ylab, xlab=xlab,cex.lab=1.2,cex.axis=1.25,lwd=2)
plotrix::boxed.labels(20, 1.3E-2, expression(paste(italic(k)," = ",(2.425*10^{2})+(7.038*10^{-5}*t))),cex=1.25,bg="white",border=NA) 
plotrix::boxed.labels(20, 1.1E-2, expression(-20<paste(italic(t)<40)),cex=1.25,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=7.5----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(0,60,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("VAPOR DENSITY OVER WATER (", italic(rho)[v], ") kg m"^{-3}))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,WETAIR(db=tairs, rh=100)$vd,type='l',ylim=c(0E-1,9E-2),ylab=ylab,xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)

tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,WETAIR(db=tairs, rh=100)$vd,type='n',ylim=c(0E-1,9E-2),ylab=ylab,xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
for(rh in c(0,10,20,30,40,50,60,80,100)){
points(tairs,WETAIR(db=tairs, rh=rh)$vd,type='l',lwd=2)
}
for(wb in seq(0,50,2)){
points(c(wb,wb+WETAIR(db=wb, wb=wb, rh=100)$esat/(0.000660 * (1.0 + 0.00115 * wb) * 101325)),c(WETAIR(db=wb, wb=wb, rh=100)$vd,0),type='l',lwd=2)
}
plotrix::boxed.labels(30, 0.08, expression(paste(italic(rho)[v]," = ",frac(italic(e),461.5*(italic(t)+273.15)))),cex=1.25,bg="white",border=NA) 
xlab = expression(paste("RELATIVE HUMIDITY (",italic(rh),"), %"))
mtext(xlab, side=4, line=1.25,cex=1.25)
text(27,.044,expression(paste("WET BULB TEMPERATURE ( ",italic(t)[w],") ",degree*C)),cex=1.25, srt=45)
for(rh in seq(0,60,10)){
  mtext(at=c(60,WETAIR(db=60, rh=rh)$vd),text=rh,cex=1.25,side=4,srt=90)
}
for(wb in seq(10,50,10)){
  plotrix::boxed.labels(wb-1,WETAIR(db=wb, rh=100)$vd+.002,wb,cex=1.25,bg="white",border=NA)
}

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("VAPOR PRESSURE OVER WATER (", italic(e), "), Pa"))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,WETAIR(db=tairs, rh=100)$e,type='n',ylim=c(0,10E3),ylab=ylab,xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)

tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,WETAIR(db=tairs, rh=100)$e,type='l',ylim=c(0,10E3),ylab=ylab,xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
points(tairs,WETAIR(db=tairs, rh=50)$e,type='l',lwd=2)
points(tairs,WETAIR(db=tairs, rh=25)$e,type='l',lwd=2)
plotrix::boxed.labels(0,8000, expression(paste(italic(e)," = ",461.50*italic(rho)[v]*(t+273.15))),cex=1.25,bg="white",border=NA) 
plotrix::boxed.labels(26, 4200, expression(paste(100,"% ",italic(rh))),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(32, 3000, expression(paste(50,"% ",italic(rh))),cex=1.05,bg="white",border=NA)
plotrix::boxed.labels(37, 2000, expression(paste(25,"% ",italic(rh))),cex=1.05,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("VIRTUAL TEMPERATURE INCREMENT AT 100 000 Pa ( ", Delta*italic(T)[v], "), K"))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,WETAIR(db=tairs, rh=100)$tvinc,type='n',ylim=c(0,20),ylab=ylab,xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)

tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,WETAIR(db=tairs, rh=100)$tvinc,type='l',ylim=c(0,20),ylab=ylab,xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
points(tairs,WETAIR(db=tairs, rh=50)$tvinc,type='l',lwd=2)
points(tairs,WETAIR(db=tairs, rh=25)$tvinc,type='l',lwd=2)
plotrix::boxed.labels(0,14, expression(paste(Delta*italic(T)[v]," = ",italic(T)*bgroup("[",frac(1+(frac(italic(r)[w],0.622)),1+italic(r)[w]),"]")-italic(T))),cex=1.25,bg="white",border=NA) 
plotrix::boxed.labels(0,11, expression(paste(italic(T)," = ",italic(t)+273.15)),cex=1.25,bg="white",border=NA) 
plotrix::boxed.labels(26, 5.5, expression(paste(100,"% ",italic(rh))),cex=1.05,bg="white",border=NA) 
plotrix::boxed.labels(32, 4, expression(paste(50,"% ",italic(rh))),cex=1.05,bg="white",border=NA)
plotrix::boxed.labels(37, 2.5, expression(paste(25,"% ",italic(rh))),cex=1.05,bg="white",border=NA) 

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(2,1)) # set up for 2 plots in 1 columns
tairs=seq(0,50,5)
par(oma = c(5,5,2,2) + 0.1)
par(mar = c(0,0,1,0) + 0.1)
par(xaxs="i")
par(yaxs="i")
ylab = expression(paste("-1* WATER POTENTIAL ( ", italic(phi), "), PA"))
xlab = expression(paste("RELATIVE HUMIDITY (%)"))

rhs=seq(65,100,5)
plot(rhs,WETAIR(rh=rhs,db=0)$wtrpot*-1,type='n',ylab="", xlab="",xaxt='n',yaxt='n',cex.lab=1.2,cex.axis=1.5,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)

tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(rhs,WETAIR(rh=rhs,db=0)$wtrpot*-1,type='l',ylab="", xlab="",xaxt='n',yaxt='n',cex.lab=1.2,cex.axis=1.5,lwd=2)
points(rhs,WETAIR(rh=rhs,db=50)$wtrpot*-1,type='l',lwd=2)
axis(side=3, at=seq(65,100,10))
axis(side=2, at=seq(0,5,1)*1E7)
plotrix::boxed.labels(80,0.5E7, expression(paste(psi," = ",(4.615*10^5)*(italic(t)+273.15)*ln(italic(rh)/100))),cex=1.25,bg="white",border=NA) 
plotrix::boxed.labels(76,2.8E7,expression(paste(0,degree*C)),cex=1.25,bg="white",border=NA)
plotrix::boxed.labels(87,2.8E7,expression(paste(50,degree*C)),cex=1.25,bg="white",border=NA)

rhs=seq(99.3,100,0.1)
plot(rhs,WETAIR(rh=rhs,db=0)$wtrpot*-1,type='n',ylab="", xlab="",xaxt='n',yaxt='n',cex.lab=1.2,cex.axis=1.5,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)

tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(rhs,WETAIR(rh=rhs,db=0)$wtrpot*-1,type='l',ylab="", xlab="",xaxt='n',yaxt='n',cex.lab=1.2,cex.axis=1.5,lwd=2)
points(rhs,WETAIR(rh=rhs,db=50)$wtrpot*-1,type='l',lwd=2)
axis(side=1, at=seq(99.4,100,0.2))
axis(side=2, at=seq(0,9,1)*1E5)
mtext(text=ylab, side=2, outer=TRUE,cex=1.25,line=2)
mtext(text=xlab, side=1, outer=TRUE,cex=1.25,line=2)
plotrix::boxed.labels(99.5,5E5,expression(paste(0,degree*C)),cex=1.25,bg="white",border=NA)
plotrix::boxed.labels(99.75,5E5,expression(paste(50,degree*C)),cex=1.25,bg="white",border=NA)

## ----echo=FALSE, message=FALSE, results='asis', fig.width=7, fig.height=8----
par(mfrow = c(1,1)) # set up for 2 plots in 1 columns
par(xaxs="i")
par(yaxs="i")
tairs=seq(-20,50,5)
par(mar = c(5,5,4,2) + 0.1)
ylab = expression(paste("WAVELENGTH OF MAXIMUM EMITTANCE FROM BLACK-BODY (  ", italic(lambda)[m], ") m"))
xlab = expression(paste("AIR TEMPERATURE (",italic(t),") ",degree*C))
plot(tairs,DRYAIR(db=tairs)$emtmax,type='n',ylim=c(0.7E-5,1.2E-5),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
Hmisc::minor.tick(nx=10,ny=10)

tix<-minor.ticks(nx=10,ny=10)
abline(h = tix$yticks, col = 'grey', lty = 1)
abline(v = tix$xticks, col = 'grey', lty = 1)
grid(lty=1,col=1)
points(tairs,DRYAIR(db=tairs)$emtmax,type='l',ylim=c(0.7E-5,1.2E-5),ylab=ylab, xlab=xlab,cex.lab=1.25,cex.axis=1.25,lwd=2)
plotrix::boxed.labels(0, 0.9E-5,expression(paste(italic(lambda)[m]," = ",2.897*10^{3}/(italic(t)+273.15))),cex=1.25,bg="white",border=NA) 

## ---- results='asis', fig.width=7, fig.height=8, warning=FALSE, message=FALSE----
library(NicheMapR) # load the NicheMapR package
library(raster) # package for working with rasters
library(ncdf4) # package for dealing with netcdf files (a kind of layered raster)
# read the global_climate dataset
global_climate<-brick(paste("c:/globalclimate/global_climate.nc",sep="")) 
 # 38th layer of global_climate is min January air temperature*10
Tair_min_january=global_climate[[38]]/10
# 50th layer of global_climate is max January air temperature*10
Tair_max_january=global_climate[[50]]/10 
# 62nd layer of global_climate is max January relative humidity*10
RH_min_january=global_climate[[62]]/10 
# 74th layer of global_climate is max January relative humidity*10
RH_max_january=global_climate[[74]]/10 
# use WETAIR.rh to get the vapor pressure for January based on min
# relative humidty and max air temperature
e=WETAIR.rh(rh=RH_min_january,db=Tair_max_january)$e 
# use the VAPPRS function to get the saturation vapor pressure for
# the new temperature, Tmin January
esat=VAPPRS(Tair_min_january) 
# compute new relative humidity for minimum air temperature
RH_max_january=(e/esat)*100 
# conditional replace of any values >100 with 100
values(RH_max_january) <- ifelse(values(RH_max_january > 100), 100, values(RH_max_january)) 
# now plot the results
par(mfrow = c(3,2)) # set up for 6 plots in 2 columns
# plot the January max air temperature
plot(Tair_max_january,zlim=c(-40,50),main="Max Air Temperature, January")
# plot the January min relative humidity
plot(RH_min_january,zlim=c(0,100),main="Min Relative Humidity, January") 
plot(e, main="vapor pressure, January")
plot(esat, main="sat. vapor pressure at Tmin")
# plot the January max air temperature
plot(Tair_min_january,zlim=c(-40,50),main="Min Air Temperature, January") 
# plot the January min relative humidity
plot(RH_max_january,zlim=c(0,100),main="Max Relative Humidity, January") 

## ------------------------------------------------------------------------
#' DRYAIR
#'
#' Calculates several properties of dry air and related characteristics shown
#' as output variables below. The program is based on equations from List, R. J. 1971.
#' Smithsonian Meteorological Tables. Smithsonian Institution Press. Washington, DC.
#' WETAIR must be used in conjunction with function VAPPRS.
#'
#' The user must supply values for the input variables (db, bp and alt).
#' If alt is known (-1000 < alt < 20000) but not BP, then set BP=0
#' @param db Dry bulb temperature (degrees C)
#' @param bp Barometric pressure (pascal)
#' @param alt Altitude (m)
#' @return patmos Standard atmospheric pressure (P)
#' @return densty Density (kg m-3)
#' @return visdyn Dynamic viscosity (kg m-1 s-1)
#' @return viskin Kinematic viscosity (m2 s-1)
#' @return difvpr Diffusivity of water vapour in air (m2 s-1)
#' @return thcond Thermal conductivity (W m-1 K-1)
#' @return htovpr Latent heat of vapourisation of water (J kg-1)
#' @return tcoeff Temperature coefficient of volume expansion (K-1)
#' @return ggroup Group of variables in Grashof number (1-m3 -K)
#' @return bbemit black body emittance (W m-2)
#' @return emtmax Wave length of maximum emittance (m)
#' @export
DRYAIR <- function(db=db, bp=0, alt=0){
  tstd=273.15
  pstd=101325.
  patmos=pstd*((1-(0.0065*alt/288))^(1/0.190284))
  bp=rep(bp,length(patmos))
  bp[bp<=0]=patmos[bp<=0]
  densty=bp/(287.04*(db+tstd))
  visnot=1.8325E-5
  tnot=296.16
  c=120.
  visdyn=(visnot*((tnot+c)/(db+tstd+c)))*(((db+tstd)/tnot)^1.5)
  viskin=visdyn/densty
  difvpr=2.26E-5*(((db+tstd)/tstd)^1.81)*(1.E5/bp)
  thcond=0.02425+(7.038E-5*db)
  htovpr=2.5012E6-2.3787E3*db
  tcoeff=1./(db+tstd)
  ggroup=0.0980616*tcoeff/(viskin*viskin)
  bbemit=5.670367E-8*((db+tstd)^4)
  emtmax=2.897E-3/(db+tstd)
  return(list(patmos=patmos, densty=densty, visdyn=visdyn, viskin=viskin, difvpr=difvpr,
    thcond=thcond, htovpr=htovpr, tcoeff=tcoeff, ggroup=ggroup, bbemit=bbemit, emtmax=emtmax))
}

## ------------------------------------------------------------------------
#' WETAIR
#'
#' Calculates several properties of humid air as output variables below. The program
#' is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
#' Institution Press. Washington, DC. WETAIR must be used in conjunction with function VAPPRS.
#'
#' Input variables are shown below. The user must supply known values for DB and BP 
#' (BP at one standard atmosphere is 101 325 pascals). Values for the remaining variables
#' are determined by whether the user has either (1) psychrometric data (WB or RH),
#' or (2) hygrometric data (DP)
#'
#' (1) Psychrometric data:
#' If WB is known but not RH, then set RH=-1 and DP=999
#' If RH is known but not WB then set WB=0 and DP=999
#'
#' (2) Hygrometric data:
#' If DP is known then set WB = 0 and RH = 0.
#' @param db Dry bulb temperature (degrees C)
#' @param wb Wet bulb temperature (degrees C)
#' @param rh Relative humidity (\%)
#' @param dp Dew point temperature (degrees C)
#' @param bp Barometric pressure (pascal)
#' @return e Vapour pressure (P)
#' @return esat Saturation vapour pressure (P)
#' @return vd Vapour density (kg m-3)
#' @return rw Mixing ration (kg kg-1)
#' @return tvir Virtual temperature (K)
#' @return tvinc Virtual temperature increment (K)
#' @return denair Hourly predictions of the soil moisture under the maximum specified shade
#' @return cp Specific heat of air at constant pressure (J kg-1 K-1)
#' @return wtrpot Water potential (P)
#' @return Relative humidity (\%)
#' @export
WETAIR <- function(db=db, wb=db, rh=0, dp=999, bp=101325){

  tk = db + 273.15
  esat = VAPPRS(db)
  if(dp < 999.0){
    e = VAPPRS(dp)
    rh = (e / esat) * 100
  }else{
    if(min(rh) > -1){
      e = esat * rh / 100
    }else{
      wbd = db - wb
      wbsat = VAPPRS(wb)
      dltae = 0.000660 * (1.0 + 0.00115 * wb) * bp * wbd
      e = wbsat - dltae
      rh = (e / esat) * 100
    }
  }
  rw = ((0.62197 * 1.0053 * e) / (bp - 1.0053 * e))
  vd = e * 0.018016 / (0.998 * 8.31434 * tk)
  tvir = tk * ((1.0 + rw / (18.016 / 28.966)) / (1.0 + rw))
  tvinc = tvir - tk
  denair = 0.0034838 * bp / (0.999 * tvir)
  cp = (1004.84 + (rw * 1846.40)) / (1.0 + rw)
  if (min(rh)<=0.0){
    wtrpot = -999
  }else{
    wtrpot = 4.615e+5 * tk * log(rh / 100.0)
  }
  return(list(e=e, esat=esat, vd=vd, rw=rw, tvinc=tvinc, denair=denair, cp=cp,
  wtrpot=wtrpot, rh=rh))
}

## ------------------------------------------------------------------------
#' VAPPRS
#'
#' Calculates saturation vapour pressure for a given air temperature.
#' @param db Dry bulb temperature (degrees C)
#' @return esat Saturation vapour pressure (P)
#' @export
VAPPRS <- function(db=db){
  t=db+273.16
  loge=t
  loge[t<=273.16]=-9.09718*(273.16/t[t<=273.16]-1.)-3.56654*log10(273.16/t[t<=273.16])+
  .876793*(1.-t[t<=273.16]/273.16)+log10(6.1071)
  loge[t>273.16]=-7.90298*(373.16/t[t>273.16]-1.)+5.02808*log10(373.16/t[t>273.16])-
  1.3816E-07*(10.^(11.344*(1.-t[t>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16
  /t[t>273.16]-1.))-1.)+log10(1013.246)
  esat=(10.^loge)*100
  return(esat)
}

