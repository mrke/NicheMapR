% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_p_wet.R
\name{get_p_wet}
\alias{get_p_wet}
\title{Function to find effective proportion of surface area that is wet from experimental water loss data}
\usage{
get_p_wet(cut.known = 1, E_tot = 0.344 / 60 / 1000, E_cut = 0.108 / 60 / 1000, T_b = 25, RHin = 30, vel = 0.0021)
}
\arguments{
\item{cut.known}{= 0, Is cutaneous water loss supplied? 0 = no (invokes estimation by partitioning respiratory), 1 = yes}

\item{resp.known}{= 0, Is respiratory water loss supplied? 0 = no (invokes estimation by partitioning respiratory), 1 = yes}

\item{eye.known}{= 0, Is ocular water loss supplied? 0 = no (invokes estimation by ocular surface area paritioning), 1 = yes}

\item{vent.known}{= 0, Is ventilation rate supplied? 0 = no (rely on estimate from V_O2), 1 = yes}

\item{press.known}{= 0, Is barometric pressure known (1, uses altitude) or has it been provided as the 'bp' variable (0)}

\item{A_tot.known}{= 0, Is total surface area known? 0 = no (uses Meeh's formula A_tot = 10 x mass.g ^ (2/3) / 10000), 1 = yes}

\item{E_tot}{= 0.344 / 60 / 1000, Observed total evaporative water loss rate, g/s}

\item{E_cut}{= 0.108 / 60 / 1000, Observed cutaneous evaporative water loss rate, g/s}

\item{E_resp}{= 0.0371 / 60 / 1000, Observed respiratory evaporative water loss rate, g/s}

\item{E_eye}{= 0.199 / 60 / 1000, Observed ocular evaporative water loss rate, g/s}

\item{mass.g}{= 26.7, Wet mass of animal, g}

\item{A_tot}{= 0.00894, Total surface area, m2}

\item{fcond}{= 0.25, Fraction of total area area touching surfaces, dec\% (therefore no evaporation or convection)}

\item{feye}{= 0.002, Fraction of total area area comprising both eyes}

\item{eye_frac}{= 0.58, Fraction of time eye is open, 1 meaning both eyes fully open, dec\%}

\item{T_b}{= 25, Core body temperature, °C}

\item{T_s}{= 24.99, Skin (surface) temperature, °C}

\item{T_a}{= 24.98, Air temperature, °C}

\item{vel}{= 0.0021, wind speed, m / s}

\item{RHex}{= 100, relative humidity of exhaled air, \%}

\item{RHin}{= 30, relative humidity of atmosphere, \%}

\item{E_O2}{= 15, Oxygen extraction efficiency, \%}

\item{RQ}{= 0.766, Respiratory quotient, dec\%}

\item{M_1}{= 0.013, Metabolic rate equation parameter 1 V_O2 = M_1 * M ^ M_2 * 10 ^ (M_3 * T_b) based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231}

\item{M_2}{= 0.800, Metabolic rate equation parameter 2}

\item{M_3}{= 0.038, Metabolic rate equation parameter 3}

\item{V_O2_STP}{= M_1 * mass.g ^ M_2 * 10 ^ (M_3 * T_b) / 3600 / 1000, Oxygen consumption rate corrected to STP, L O2 / s}

\item{T_ref_O2}{= 25, Reference temperature for oxygen consumption rate correction (°C)}

\item{V_air_STP}{= V_O2_STP / (0.2094 * E_O2 / 100), Ventilation rate at STP, L/s}

\item{rho_flesh}{= 1000, Density of flesh (kg/m3)}

\item{T_ref_vent}{= 25, Reference temperature for ventilation rate correction (°C)}

\item{T_A}{= 8817, Arhhenius temperature (K)}

\item{T_AL}{= 50000, Arrhenius temperature for decrease below lower boundary of tolerance range T_L (K)}

\item{T_AH}{= 90000, Arrhenius temperature for decrease above upper boundary of tolerance range T_H (K)}

\item{T_L}{= 279, Lower boundary (K) of temperature tolerance range for Arrhenius thermal response}

\item{T_H}{= 306, Upper boundary (K) of temperature tolerance range for Arrhenius thermal response}

\item{alt}{= 0.03, Altitude (km)}

\item{bp}{= 101325, Barometric pressure (pascal)}
}
\value{
E_tot total evaporative water loss, g/s

E_cut cutaneous evaporative water loss, g/s

E_resp respiratory evaporative water loss, g/s

E_eye ocular evaporative water loss, g/s

R_a boundary layer resitance, s/m

R_s skin resistance, s/m

A_tot total surface area, m2

A_eye ocular surface area, m2

A_eff effective area for evaporation, m2

A_wet area acting as a free water surface, m2

p_wet proportion of surface area acting as a free water surface, dec\%
}
\description{
A set of equations for computing the proportion of the total surface area of an organism that
is acting as a free water surface for evaporation, which is used in the ectotherm model of NicheMapR
NicheMapR to compute cutaneous water loss rates. This function requires the WETAIR and DRYAIR
functions of NicheMapR.
Elia Pirtle and Michael Kearney developed this R function and example in November 2017.
}
