## ---- echo = FALSE, warning = FALSE, message = FALSE---------------------
knitr::opts_chunk$set(
 eval = TRUE, tidy.opts=list(width.cutoff=60), tidy=TRUE  
)
library(jpeg)
library(knitr)

## ----table, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "

*Symbols* | *Definition*	| *Units*  	
--- | --- | --- | --- 
$a$|constant | -
$A$|area | $m^2$
$b$| constant | -
$C$|thermal capacitance, $mc_p$ | kg-J/(kg-°C)
$C_1$,$C_2$|constants | none
$c_p$|specific heat | J/kg-°C
$E$|voltage | volts
$F$|configuration factor | -
$h$|heat transfer coefficient | W/(m$^2$-°C)
$I$|current | amps
$j$|constant | C
$k$|thermal conductivity| W/(m-°C)
$L$|characteristic dimension | m
$m$|mass | kg
$\\dot{m}$ | mass flow rate | kg/s
$P_1$,$P_2$|constants | none
$P$|d/dt | place holder to ease algebraic operations
$q^{'''}$|volume-specific heat generation | W/m$^3$
$Q$|heat flux | W
$R$|thermal resistance | 1/(W-C)
$S^2$|ellipsoid geometric fraction | m$^2$
$t$|time | s
$T$|temperature | °C
$v$|wind speed | m/s
$V$|volume | m$^3$
$x$|length | m
$\\epsilon$ | emissivity | -
$\\mu$ | dynamic viscosity | kg/(s-m)
$\\rho$ |density | kg/m$^3$
$\\sigma$ | Stefan-Boltzmann constant | W/(m$^2$-K$^4$)
$\\tau$|time constant | s
$\\theta$| inverse of $\\tau$| 1/s
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

tabl <- "

*Subscripts* | *Definition* 	
--- | --- | --- 
$a$|air 
$b$|body 
$bl$|blood 
$c$|core 
$conv$|convection 
$e$|environment 
$evap$|evaporation 
$f$|final 
$gen$|generated 
$i$|initial 
$met$|metabolism
$o$|object
$rad$|radiation 
$s$|surface 
$sk$|skin 
$sol$|solar 
$st$|stored 
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ---- echo=FALSE, out.width = "40%", fig.align="center", fig.cap = "A cross section system diagram for heat flow in an animal", fig.pos='!h'----
include_graphics("Lizard_System_Diagram1.jpg")

## ---- echo=FALSE, out.width = "70%", fig.align="center", fig.cap = "An equivalent thermal circuit diagram for the system diagram in Figure 1. This assumes a point source of heat at the centre and a slab approximation of the geometry", fig.pos='!h'----
include_graphics("Lizard_System_Diagram2.jpg")

## ---- echo=FALSE, fig.cap = "The ratio of internal and external resistance (Biot number) for desert iguana-shaped lizards of different body mass at 0.45 m/s wind speed, 24 °C air and radiant temperature, 26 °C skin temperature, and flesh thermal conductivity of 0.5 W/mC."----
masses <- seq(40, 5000, 100) # sequence of wet masses (g)
Tair <- 24 # deg C
vel <- 0.45 # m/s
T_s <- 26 + 273.15 # surface temp converted to Kelvin
T_rad <- Tair + 273.15 # radiant temperature, converted to Kelvin
press <- 101325 # air pressure (Pa)
kflesh <- 0.5 # W/mC
Biot <- vector(length = length(masses))
for(i in 1:length(masses)){
  mass<-masses[i]
  DENSTY <- press / (287.04 * (Tair + 273.15)) # air density, kg/m3
  THCOND <- 0.02425 + (7.038 * 10 ^ -5 * Tair) # air thermal conductivity, W/(m.K)
  VISDYN <- (1.8325 * 10 ^ -5 * ((296.16 + 120) / ((Tair + 273.15) + 120))) * (((Tair + 273.15) / 296.16) ^ 1.5) # dynamic viscosity of air, kg/(m.s)
  spheat <- 3474 # specific heat of flesh (J/kgC)
  rho <- 932 # density of flesh (kg/m3)
  ATOT <- (10.4713 * mass ^ .688) / 10000.
  m <- mass / 1000 # convert mass to kg
  C <- m * spheat # thermal capacitance, J/C
  V <- m / rho # volume, m3
  L <- V ^ (1 / 3) # characteristic dimension, m
  x_sk <- (L / ((0.04 / rho) ^ (1 / 3)) ) * 0.001 # scale skin depth with body length relative to 40g desert iguana
  #x_sk <- 0.001
  Re <- DENSTY * vel * L / VISDYN # Reynolds number
  PR <- 1005.7 * VISDYN / THCOND # Prandlt number    
  NUfor <- 0.35 * Re ^ 0.6 # Nusselt number
  hc_forced <- NUfor * THCOND / L # convection coefficent, forced
  GR <- abs(DENSTY ^ 2 * (1 / (Tair + 273.15)) * 9.80665 * L ^ 3 * (T_s - Tair) / VISDYN ^ 2) # Grashof number
  Raylei <- GR * PR # Rayleigh number
  if (Raylei < 1.0e-05) {
    NUfre = 0.4
  } else{
    if (Raylei < 0.1) {
      NUfre = 0.976 * Raylei ^ 0.0784
    } else{
      if (Raylei < 100) {
        NUfre = 1.1173 * Raylei ^ 0.1344
      } else{
        if (Raylei < 10000.) {
          NUfre = 0.7455 * Raylei ^ 0.2167
        } else{
          if (Raylei < 1.0e+09) {
            NUfre = 0.5168 * Raylei ^ 0.2501
          } else{
            if (Raylei < 1.0e+12) {
              NUfre = 0.5168 * Raylei ^ 0.2501
            } else{
              NUfre = 0.5168 * Raylei ^ 0.2501
            }
          }
        }
      }
    }
  }
  hc_free <- NUfre * THCOND / L # convection coefficent, forced
  hc_comb <- hc_free + hc_forced # convection coefficent, free plus forced
  R_conv <- 1 / (hc_comb * ATOT)
  
  R_sk <- (x_sk / 2) / (kflesh * ATOT)
  sigma <- 5.67E-8
  emis <- 1
  F_o_e <- 0.8
  spheat <- 3474
  m_dot_bl <- 3.5 * 0.4 / 1000 / 60 # blood flow g/(min-100g) to kg/s for 40 g lizard
  R_bl <- 1 / (m_dot_bl * spheat)
  R_rad <- (T_s - T_rad)/(emis * sigma * F_o_e * ATOT * (T_s^4 - T_rad^4))
  R_b <- (R_sk * R_bl)/(R_sk + R_bl)
  R_e <- (R_conv * R_rad)/(R_conv + R_rad)
  Biot[i] <- R_b/R_e
}

plot(Biot ~ masses, ylab = 'Biot number', xlab = 'mass (g)', ylim = c(0, .25))
abline(h=0.1)

## ---- echo=FALSE, out.width = "40%", fig.align="center", fig.cap = "A one-lump model of heat loss assuming a point (line) source of heat in the centre of the cylindrical geometry", fig.pos='!h'----
include_graphics("C:/Users/mrke/Dropbox/Current Research Projects/global transient/manuscript/Lizard_System_Diagram3.jpg")

## ---- echo=FALSE, out.width = "40%", fig.align="center", fig.cap = "A one lump model of heat exchange assuming distributed metabolic heat generation", fig.pos='!h'----
include_graphics("Bat_System_Diagram.jpg")

