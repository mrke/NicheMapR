#' Ectotherm model
#'
#' An implementation of the Niche Mapper ectotherm model that computes body temperature, water loss,
#' activity and microhabitat selection. It optionally runs the Dynamic Energy Budget (DEB) model for
#' computing mass budgets (inc. water budgets) and growth, development, reproduction trajectories
#' as constrained by food, activity and temperature (see Details). When not running the DEB model
#' a user-specified mass is used as well as a allometric (mass and body temperature) function to
#' compute metabolic rate. \cr\cr NOTE: The microclimate model, e.g. \code{\link{micro_global}}, must be run prior to running the ectotherm model
#'
#' @encoding UTF-8
#' @param Ww_g = 40, Wet weight of animal (g), note this model is 'steady state' so no lags in heating/cooling due to mass
#' @param shape = 3, Organism shape, 0-5, Determines whether standard or custom shapes/surface area/volume relationships are used: 0=plate, 1=cyl, 2=ellips, 3=lizard (desert iguana), 4=frog (leopard frog), 5=custom (see details)
#' @param alpha_max = 0.85, Maximum solar absorptivity, 0-1
#' @param alpha_min = 0.85, Maximum solar absorptivity, 0-1
#' @param T_F_min = 24, Minimum foraging temperature, °C (also affects burrow depth selection)
#' @param T_F_max = 34, Maximum foraging temperature, °C
#' @param T_B_min = 17.5, Minimum basking temperature, °C
#' @param T_RB_min = 17.5, Minimum temperature at which animal will move from retreat to basking site, °C
#' @param T_pref = 30, Preferred body temperature, °C
#' @param CT_max = 40, Critical thermal maximum, °C (affects burrow depth selection and may be used to impose death from heat stress)
#' @param CT_min = 6, Critical thermal minimum, °C (affects burrow depth selection and may be used to impose death from cold stress)
#' @param diurn = 1, Diurnal activity allowed?  1=yes, 0=no
#' @param nocturn = 0, Nocturnal activity allowed?  1=yes, 0=no
#' @param crepus = 0, Crepuscular activity allowed?  1=yes, 0=no
#' @param shade_seek = 1, Shade seeking allowed?  1=yes, 0=no
#' @param burrow = 1 Shelter in burrow allowed?  1=yes, 0=no
#' @param climb = 0, Climbing to seek cooler habitats allowed?  1=yes, 0=no
#' @param shdburrow = 0, Choose if the animal's retreat is in the open (0), in the shade when above or below CTmin in sun (1) or in shade always (2)
#' @param mindepth = 2, Minimum depth (soil node #) to which animal can retreat if burrowing
#' @param maxdepth = 10, Maximum depth (soil node #) to which animal can retreat if burrowing
#' @param aestdepth = 10, Depth (soil node #) to which animal retreats if burrowing and aestivating due to desiccation
#' @param M_1 = 0.013, Metabolic rate parameter 1 V_O2=M_1*M^M_2*10^(M_3*Tb), in ml O2 / h, default parameters for lizards based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231
#' @param M_2 = 0.800, Metabolic rate parameter 2
#' @param M_3 = 0.038, Metabolic rate parameter 3
#' @param pct_wet = 0.2, \% of surface area acting as a free-water exchanger, for computing cutaneous water loss
#' @param pct_eyes = 0.03, \% of surface area taken up by open eyes, for computing ocular water loss (only when active)
#' @param pct_mouth = 5, \% of surface area taken up by open mouth, for computing panting water loss
#' @param pantmax = 5, maximum multiplier on breathing rate, for respiratory water loss via panting (value of 1 prevents panting)
#' @param F_O2 = 20, \% oxygen extraction efficiency, for respiratory water loss
#' @param delta_air = 0.1, Temperature difference (°C) between expired and inspired air, for computing respiratory water loss
#' @usage ectotherm(Ww_g, shape, alpha_max, alpha_min, T_F_min, T_F_max, T_B_min, T_RB_min, CT_max, CT_min,
#'  T_pref, diurn, nocturn, crepus, shade_seek, burrow, climb, shdburrow, mindepth, maxdepth,
#'  M_1, M_2, M_3, pct_wet, F_O2, delta_air, ...)
#' @details
#' \strong{ Parameters controling how the model runs:}
#' \itemize{
#' \item{\code{nyears}{ = micro$nyears, Number of years the simulation runs for - must be consistent with dimensions of environmental input data}\cr}
#' \item{\code{enberr}{ = 0.01, Factor by which the mass is multiplied to obtain a tolerance level for the heat budget solution}\cr}
#' \item{\code{live}{ = 1, Live (metabolism/behaviour) or dead animal?}\cr}
#' \item{\code{transient}{ = 0, Run a transient (i.e. include heat storage) simulation (1=yes, 0=no)? No behaviour yet - assums full sun}\cr}
#' \item{\code{delta_shade}{ = 3, Percent shade increment step, 0-100\%, allowing different thermoregulatory precision (smaller values increase run time)}\cr}
#' \item{\code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1 = yes, 0 = no}\cr}
#' \item{\code{write_csv}{ = 0, Write csv files of final output? 0 = no, 1 = just yearout and yearsout, 2 = all output}\cr}
#' \item{\code{startday}{ = 1, Day of year at which simulation starts}\cr}
#'}
#' \strong{ Environmental inputs:}
#'
#' \itemize{
#' \item{\code{minshades}{ = micro$minshade, Vector of daily minimum shade values - can be different to value used in microclimate model (e.g. to simulate sunspot tracking on a forest floor) (\%)}\cr}
#' \item{\code{maxshades}{ = micro$maxshade, Vector of daily maximum shade values - can be different to value used in microclimate model (e.g. to simulate use of fine-scale shade in a generally unshaded habitat) (\%)}\cr}
#' \item{\code{fluid}{ = 0, Fluid type 0=air, 1=water }\cr}
#' \item{\code{O2gas}{ = 20.95, oxygen concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)}\cr}
#' \item{\code{N2gas}{ = 79.02, nitrogen concetration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)}\cr}
#' \item{\code{CO2gas}{ = 0.0412, carbon dioxide concentration of air, to account for non-atmospheric concentrations e.g. in burrows (\%)}\cr}
#' \item{\code{alpha_sub}{ = 1 - micro$REFL, Vector of daily substrate reflectances (0-1)}\cr}
#' \item{\code{epsilon_sub}{ = 1, Emissivity of substrate (0-1)}\cr}
#' \item{\code{epsilon_sky}{ = 1, Emissivity of sky (0-1)}\cr}
#' \item{\code{PDIF}{ = 0.1, Fraction of total solar radiation that is diffuse (0-1), will ultimately be made an optionally hourly vector}\cr}
#' \item{\code{DEP}{ = micro$DEP, Depths available from the microclimate model simulation}\cr}
#' \item{\code{KS}{ = micro$KS[seq(1, 19, 2)], Depth-specific saturated hydraulic conductivity (kg s/m3) from the microclimate model simulation, for modelling liquid exchange with substrate}\cr}
#' \item{\code{b}{ = micro$BB[seq(1, 19, 2)], Depth-specific Campbell's b parameter (-) from the microclimate model simulation, for modelling liquid exchange with substrate}\cr}
#' \item{\code{PE}{ = micro$PE[seq(1, 19, 2)], Depth-specific air-entry water potential (J / kg) from the microclimate model simulation, for modelling liquid exchange with substrate}\cr}
#' \item{\code{metout}{ = micro$metout, Microclimate model output for above ground, minimum shade conditions}\cr}
#' \item{\code{shadmet}{ = micro$shadmet, Microclimate model output for above ground, maximum shade conditions}\cr}
#' \item{\code{soil}{ = micro$soil, Microclimate model output for soil temperature, minimum shade conditions}\cr}
#' \item{\code{shadsoil}{ = micro$shadsoil, Microclimate model output for soil temperature, maximum shade conditions}\cr}
#' \item{\code{soilmoist}{ = micro$soilmoist, Microclimate model output for soil moisture, minimum shade conditions}\cr}
#' \item{\code{shadmoist}{ = micro$shadmoist, Microclimate model output for soil moisture, maximum shade conditions}\cr}
#' \item{\code{humid}{ = micro$humid, Microclimate model output for soil humidity, minimum shade conditions}\cr}
#' \item{\code{shadhumid}{ = micro$shadhumid, Microclimate model output for soil humidity, maximum shade conditions}\cr}
#' \item{\code{soilpot}{ = micro$soilpot, Microclimate model output for soil water potential, minimum shade conditions}\cr}
#' \item{\code{shadpot}{ = micro$shadpot, Microclimate model output for soil water potential, maximum shade conditions}\cr}
#' \item{\code{tcond}{ = micro$tcond, Microclimate model output for soil thermal conductivity, minimum shade conditions}\cr}
#' \item{\code{shadtcond}{ = micro$shadtcond, Microclimate model output for soil thermal conductivity, maximum shade conditions}\cr}
#' \item{\code{rainfall}{ = micro$RAINFALL, Vector of daily rainfall (mm)}\cr}
#' \item{\code{rainhr}{ = rep(-1,nrow(metout)), Vector of hourly rainfall (mm), overwrites rainfall if not negative}\cr}
#' \item{\code{preshr}{ = rep(101325 * ((1 - (0.0065 * as.numeric(micro$elev) / 288)) ^ (1/0.190284)), nrow(metout)), Vector of hourly atmospheric pressure (Pa), defaulting to elevation-adjusted values}\cr}
#' \item{\code{elev}{ = as.numeric(micro$elev), Elevation of simulation (m), obtained from microclimate model output by default}\cr}
#' \item{\code{longitude}{ = micro$longlat[1], Longitude (decimal degrees), obtained from microclimate model output by default}\cr}
#' \item{\code{latitude}{ = micro$longlat[2], Latitude (decimal degrees), obtained from microclimate model output by default}\cr}
#'}
#' \strong{ Morphological parameters:}
#'
#' \itemize{
#' \item{\code{custom_shape}{ = c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743), Custom shape coefficients. Operates if shape=5, and consists of 4 pairs of values representing the parameters a and b of a relationship AREA=a*mass^b, where AREA is in cm2 and mass is in g. The first pair are a and b for total surface area, then a and b for ventral area, then for sillhouette area normal to the sun, then sillhouette area perpendicular to the sun}\cr}
#' \item{\code{shape_a}{ = 1, Proportionality factor (-) for going from volume to area, keep this 1 (redundant parameter that should be removed)}\cr}
#' \item{\code{shape_b}{ = 3, Proportionality factor (-) for going from volume to area, represents ratio of width:height for a plate, length:diameter for cylinder, b axis:a axis for ellipsoid }\cr}
#' \item{\code{shape_c}{ = 0.6666666667, Proportionality factor (-) for going from volume to area, represents ratio of length:height for a plate, c axis:a axis for ellipsoid}\cr}
#' \item{\code{fatosk}{ = 0.4, Configuration factor to sky (-) for infrared calculations}\cr}
#' \item{\code{fatosb}{ = 0.4, Configuration factor to subsrate for infrared calculations}\cr}
#' \item{\code{rinsul}{ = 0, Insulative fat layer thickness (not yet functional) (m)}\cr}
#' \item{\code{pct_cond}{ = 10, Percentage of animal surface contacting the substrate (\%)}\cr}
#' \item{\code{pct_touch}{ = 0, Percentage of animal surface area contacting another animal of same temperature (\%)}\cr}
#' \item{\code{c_body}{ = 3073, Specific heat of flesh J/(kg-K)}\cr}
#' \item{\code{k_flesh}{ = 0.5, Thermal conductivity of flesh (W/mC, range: 0.412-2.8)}\cr}
#' \item{\code{rho_body}{ = 1000, Density of flesh (kg/m3)}\cr}
#' \item{\code{epsilon}{ = 0.95, Emissivity of animal (0-1)}\cr}
#' \item{\code{eggshape_a}{ = 1, Proportionality factor (-) for going from volume to area, keep this 1 (redundant parameter that should be removed)}\cr}
#' \item{\code{eggshape_b}{ = 0.6666666667, Proportionality factor (-) for going from volume to area, represents ratio of width:height for a plate, length:diameter for cylinder, b axis:a axis for ellipsoid }\cr}
#' \item{\code{eggshape_c}{ = 0.6666666667, Proportionality factor (-) for going from volume to area, represents ratio of length:height for a plate, c axis:a axis for ellipsoid}\cr}
#' \item{\code{eggmult }{ = 1 # multiply egg mass by clutch size for heat and water exchange calculations?}\cr}
#' \item{\code{pct_cond_egg}{ = 50, Percentage of egg surface contacting the substrate (\%)}\cr}
#'}
#' \strong{ Behavioural parameters:}
#'
#' \itemize{
#' \item{\code{postur}{ = 1, postural orientation to sun, 1 = perpendicular, 2 = parallel, 0 = half way between, relevant if live = 0}\cr}
#' \item{\code{warmsig}{ = 0, Warming signal for emergence? °C/h (if in burrow deeper than node 2, change in burrow temp must be exceed warmsig)}\cr}
#' \item{\code{fossorial}{ = 0, Fossorial activity? 1=yes, 0=no (this option hasn't been properly implemented)}\cr}
#' \item{\code{rainact}{ = 0, Activity is limited by rainfall? 1=yes, 0=no, threshold rainfall for activity set by \code{actrainthresh}}\cr}
#' \item{\code{actrainthresh}{ = 0.1, Threshold (mm) of rain causing activity if \code{rainact}=1}\cr}
#' \item{\code{soilnode}{ = 4, Soil node (1-10, corresponding to values in \code{DEP}) at which eggs are laid (overridden if \code{amphibreed}=1)}\cr}
#' \item{\code{eggshade}{ = 0, are eggs laid in shade? 0=no, 1=yes}\cr}
#' \item{\code{aquabask}{ = 0, If aquatic, does it bask? 0=no, stay at water temp, 1=yes, when not hungry, 2=all the time}\cr}
#'}
#' \strong{ Thermal physiological parameters:}
#'
#' \itemize{
#' \item{\code{CT_minthresh}{ = 12, Number of consecutive hours below CT_min that leads to death - simulation will terminate beyond this threshold if \code{CT_kill}=1}\cr}
#' \item{\code{CT_kill}{ = 0, Animal dies when it hits critical thermal limits? 1=yes, 0=no}\cr}
#'}
#' \strong{ Water and food budget parameters (first eight only relevant if \code{DEB}=1):}
#'
#' \itemize{
#' \item{\code{pct_H_P}{ = 73, Water in faeces (product) (\%)}\cr}
#' \item{\code{pct_H_N}{ = 0, Water in excreted nitrogenous waste (\%)}\cr}
#' \item{\code{pct_H_X}{ = 82, Water content of food (\%)}\cr}
#' \item{\code{pct_H_R}{ = 15, Minimum tolerated dehydration (\% of wet mass) - prohibits foraging if greater than this}\cr}
#' \item{\code{pct_H_death}{ = 35, Maximum tolerated dehydration (\% of wet mass) - causes death if greater than this}\cr}
#' \item{\code{gutfill}{ = 75, Gut fill (\%) at which satiation occurs - if greater than 100\%, animal always tries to forage}\cr}
#' \item{\code{raindrink}{ = 0, Rainfall level at which rehydration from drinking occurs - if 0 animal can always drink}\cr}
#' \item{\code{foodlim}{ = 1, Is the animal food limited - if 0 animal can always find food (useful for making different life stages dependent on soil moisture-based food estimates}\cr}
#' \item{\code{RQ}{ = 0.8, respiratory quotient (0-1), computed from first principles if DEB model running}\cr}
#' \item{\code{K_skin}{ = 2.8e-09, - Hydraulic conductivity of skin (kg/(m s (J/kg)) - drives liquid water exchange with substrate}\cr}
#' \item{\code{spec_hyd_body}{ = 0.000304, Specific hydration of body (m3 / (m3 (J/kg))) - drives liquid water exchange with substrate if K_skin > 0 }\cr}
#' \item{\code{psi_body}{ = -707, Water potential of body (J/kg) - drives liquid water exchange with substrate if K_skin > 0 and will also affect skin humidity for water vapour exchange}\cr}
#' \item{\code{K_egg}{ = 2.8e-09, Hydraulic conductivity of egg shell (kg/(m s (J/kg)) - drives liquid water exchange with substrate}\cr}
#' \item{\code{spec_hyd_egg}{ = 0.000304, Specific hydration of egg (m3 / (m3 (J/kg))) - drives liquid water exchange with substrate if K_skin > 0 }\cr}
#' \item{\code{psi_egg}{ = -707, Water potential of egg (J/kg) - drives liquid water exchange with substrate if K_skin > 0}\cr}
#' }
#' \strong{ Dynamic Energy Budget (DEB) model parameters:}
#' \itemize{
#' \item{\code{DEB}{ = 0, Run the DEB model (1) or just heat balance (0). Latter uses allometrically predicted respiration base on \code{M_1}, \code{M_2} and \code{M_3}}\cr}
#' \item{\code{intmethod}{ = 1, Use Euler (0) or DOPRI (1) method for integrating non-insect DEB model. Latter will be more accurate but slower}\cr}
#' \item{\code{metab_mode}{ = 0, Run insect model? 0 = no, 1 = hemimetabolus model (abp DEB model), 2 = holometabolous model (hex DEB model)}\cr}
#' \item{\code{z_mult}{ = 1, Scaling factor for DEB body-size covariation relationships - use it to make a metabolically scaled larger or smaller version of your animal}\cr}
#'}
#' \strong{ Core DEB parameters:}
#' \itemize{
#' \item{\code{z}{ = 2.825*z_mult, Zoom factor (cm)}\cr}
#' \item{\code{del_M}{ =  0.2144, Shape coefficient (-)}\cr}
#' \item{\code{p_Xm}{ = 12420, Surface area-specific maximum feeding rate J/cm2/h}\cr}
#' \item{\code{kap_X}{ = 0.85, Digestive efficiency (0-1)}\cr}
#' \item{\code{v}{ = 0.02795/24, Energy conductance (cm/h)}\cr}
#' \item{\code{kap}{ = 0.8206, Fraction of mobilised reserve allocated to soma}\cr}
#' \item{\code{p_M}{ = 48.81/24, Volume-specific somatic maintenance (J/cm3/h)}\cr}
#' \item{\code{E_G}{ = 7512, Cost of structure, including overheads (J/cm3)}\cr}
#' \item{\code{kap_R}{ = 0.95, Fraction of reproduction energy fixed in eggs}\cr}
#' \item{\code{k_J}{ = 0.006498/24, Maturity maintenance rate coefficient (1/h)}\cr}
#' \item{\code{E_Hb}{ = 866.6*z_mult^3, Maturity at birth (J)}\cr}
#' \item{\code{E_Hj}{ = E_Hb*z_mult^3, Maturity at metamorphosis (if different to E_Hb, triggers metabolic acceleration) (J)}\cr}
#' \item{\code{E_Hp}{ = 1.019e+04*z.mult^3, Maturity at puberty}\cr}
#' \item{\code{E_He}{ = E_He*z_mult^3, Maturity at eclosion (J) (relevant only for holometabolous model)}\cr}
#' \item{\code{h_a}{ = 1.051e-08*z_mult/(24^2), Weibull ageing acceleration (1/h2)}\cr}
#' \item{\code{s_G}{ = 0.01, Gompertz stress coefficient (-)}\cr}
#' \item{\code{E_0}{ = 9220*z_mult^4, Energy content of the egg (derived from core parameters via backwards integation from birth to time zero) (J)}\cr}
#'}
#' \strong{ Thermal DEB parameters:}
#' \itemize{
#' \item{\code{arrhen_mode}{ = 1, equation used for Arrhenius rate correction - 0 is original Sharpe-Schoolfield formulation, 1 is version used in DEB tempcorr function from AmPtool}\cr}
#' \item{\code{T_REF}{ = 20 + 273.15, Reference temperature for rate correction (°C)}\cr}
#' \item{\code{T_A}{ = 8817 Arhhenius temperature}\cr}
#' \item{\code{T_AL}{ = 50000, Arrhenius temperature for decrease below lower boundary of tolerance range \code{T_L}}\cr}
#' \item{\code{T_AH}{ = 90000, Arrhenius temperature for decrease above upper boundary of tolerance range \code{T_H}}\cr}
#' \item{\code{T_L}{ = 6 + 273.15, Lower boundary (K) of temperature tolerance range for Arrhenius thermal response}\cr}
#' \item{\code{T_H}{ = 33 + 273.15, Upper boundary (K) of temperature tolerance range for Arrhenius thermal response}\cr}
#' \item{\code{T_A2}{ = 8817 Arhhenius temperature}\cr} for maturity maintenance (causes 'Temperature Size Rule' effect)
#' \item{\code{T_AL2}{ = 50000, Arrhenius temperature for decrease below lower boundary of tolerance range \code{T_L}}\cr} for maturity maintenance (causes 'Temperature Size Rule' effect)
#' \item{\code{T_AH2}{ = 90000, Arrhenius temperature for decrease above upper boundary of tolerance range \code{T_H}}\cr} for maturity maintenance (causes 'Temperature Size Rule' effect)
#' \item{\code{T_2L}{ = 6 + 273.15, Lower boundary (K) of temperature tolerance range for Arrhenius thermal response}\cr} for maturity maintenance (causes 'Temperature Size Rule' effect)
#' \item{\code{T_H2}{ = 33 + 273.15, Upper boundary (K) of temperature tolerance range for Arrhenius thermal response}\cr} for maturity maintenance (causes 'Temperature Size Rule' effect)
#'}
#' \strong{ Compound/derived DEB parameters:}
#' \itemize{
#' \item{\code{E_m}{ = (p_M * z / kap) / v, Maximum reserve density (J/cm3)}\cr}
#'}
#' \strong{ Food-related axilliary DEB parameters:}
#' \itemize{
#' \item{\code{f}{ = 1, functional response (-), usually kept at 1 because gut model controls food availability such that f=0 when gut empty}\cr}
#' \item{\code{E_sm}{ = 350, Maximum structure-specific energy density of stomach (J/cm3)}\cr}
#' \item{\code{K}{ = 1, Half saturation constant (J/cm2)}\cr}
#' \item{\code{X}{ = 10, Food density (J/cm2)}\cr}
#'}
#' \strong{ Composition-related axilliary DEB parameters:}
#' \itemize{
#' \item{\code{rho_body_deb}{ = rho_body/1000, Animal density (g/cm3)}\cr}
#' \item{\code{d_V}{ = 0.3, Dry mass fraction of structure (0-1)}\cr}
#' \item{\code{d_E}{ = 0.3, Dry mass fraction of reserve (0-1)}\cr}
#' \item{\code{d_Egg}{ = 0.3, Dry mass fraction of egg (0-1)}\cr}
#' \item{\code{stoich_mode}{ = 0, adjust chemical indices to chemical potentials (0) or vice versa (1) or leave as is (2)}\cr}
#' \item{\code{mu_X}{ = 525000, Molar Gibbs energy (chemical potential) of food (J/mol)}\cr}
#' \item{\code{mu_E}{ = 585000, Molar Gibbs energy (chemical potential) of reserve (J/mol)}\cr}
#' \item{\code{mu_V}{ = 500000, Molar Gibbs energy (chemical potential) of structure (J/mol)}\cr}
#' \item{\code{mu_P}{ = 480000, Molar Gibbs energy (chemical potential) of faeces (J/mol)}\cr}
#' \item{\code{mu_N}{ = 244e3/5, Molar Gibbs energy (chemical potential) of nitrogenous waste (J/mol), synthesis from NH3, Withers page 119}\cr}
#' \item{\code{kap_X_P}{ = 0.1, Faecation efficiency of food to faeces (-)}\cr}
#' \item{\code{n_X}{ = c(1, 1.8, 0.5, 0.15), chem. indices of C, O, H and N in food}\cr}
#' \item{\code{n_E}{ = c(1, 1.8, 0.5, 0.15), chem. indices of C, O, H and N in reserve}\cr}
#' \item{\code{n_V}{ = c(1, 1.8, 0.5, 0.15), chem. indices of C, O, H and N in structure}\cr}
#' \item{\code{n_P}{ = c(1, 1.8, 0.5, 0.15), chem. indices of C, O, H and N in faeces}\cr}
#' \item{\code{n_M_nitro}{ = c(1, 4/5, 3/5, 4/5), chem. indices of C, O, H and N in nitrogenous waste}\cr}
#' \item{\code{h_N}{ = 384238, molar enthalpy of nitrogenous waste (combustion frame of reference) (J/mol), overridden if n_M_nitro specified as urea, uric acid or ammonia}\cr}
#'}
#' \strong{ Holometabolous insect DEB model parameters:}
#' \itemize{
#' \item{\code{stages}{ = 8, number of stages = number of instars plus 1 for egg + 1 for pupa + 1 for imago}\cr}
#' \item{\code{kap_V}{ = 0.8, conversion efficient E -> V -> E (-)}\cr}
#' \item{\code{k_Ee}{ = 0.06293 / 24, reproduction buffer turnover of imago (1/h)}\cr}
#' \item{\code{k_EV}{ = 0.03111 / 24, spec decay rate of larval structure in pupa (1/h)}\cr}
#' \item{\code{S_instar}{ = rep(2.049137, stages), stress at instar n: L_n^2/ L_n-1^2 (-)}\cr}
#' \item{\code{s_j}{ = 0.9985855, Reprod buffer/structure ratio at pupation as fraction of max}\cr}
#' \item{\code{L_b}{ = 0.0734, Structural length at birth (cm)}\cr}
#'}
#' \strong{ Inital conditions for DEB model:}
#' \itemize{
#' \item{\code{V_init}{ = 3e-9, Initial structural volume (cm3)}\cr}
#' \item{\code{E_init}{ = E_0/V_init, Initial reserve density (J/cm3)}\cr}
#' \item{\code{E_H_init}{ = 0, Initial maturity (J)}\cr}
#' \item{\code{stage}{ = 0, Initial stage (STD model: 0=embryo, 1=juvenile, 2=mature but not yet reproducing, 3=beyond first reproduction, Insect models: 0=embryo, 1-(stages-1)=instar, with pupa as penultimate, stages=adult}\cr}
#'}
#' \strong{ Metabolic depression parameters:}
#' \itemize{
#' \item{\code{aestivate}{ = 0, Does the animal aestivate/go into torpor? 1=yes, 0=no}\cr}
#' \item{\code{depress}{ = 1, Fraction by which \code{p_M}, \code{k_J} and \code{v} are reduced during torpor}\cr}
#'}
#' \strong{ Reproductive phenology model parameters:}
#' \itemize{
#' \item{\code{clutchsize}{ = 5, Clutch size (#), overridden by \code{clutch_ab}}\cr}
#' \item{\code{clutch_ab}{ = c(0,0), # paramters for relationship between length (cm) and clutch size: clutch size = a*SVL-b, make a and b zero if fixed clutch size}\cr}
#' \item{\code{viviparous}{ = 1, Viviparous reproduction? 1=yes, 0=no (if yes, animal will be held in adult-sided female's body for duration of development and will experience her body temperature}\cr}
#' \item{\code{minclutch}{ = 0, Minimum clutch size if not enough in reproduction buffer for clutch size predicted by \code{clutch_ab} - if zero, will not operate}\cr}
#' \item{\code{batch}{ = 1, Invoke Pequerie et al.'s batch laying model?}\cr}
#' \item{\code{photostart}{ = 3, Photoperiod response triggering ovulation, none (0), summer solstice (1), autumnal equinox (2), winter solstice (3), vernal equinox (4), specified daylength thresholds (5 - uses \code{daylengthstart} and \code{daylengthfinish})}\cr}
#' \item{\code{photofinish}{ = 1, Photoperiod terminating ovulation, none (0), summer solstice (1), autumnal equinox (2), winter solstice (3), vernal equinox (4), specified daylength thresholds (5 - uses \code{daylengthstart} and \code{daylengthfinish})}\cr}
#' \item{\code{daylengthstart}{ = 12.5, Threshold daylength (h) for initiating breeding}\cr}
#' \item{\code{daylengthfinish}{ = 13, Threshold daylength (h) for terminating breeding}\cr}
#' \item{\code{photodirs }{ =  1, Is the start daylength trigger during a decrease (0) or increase (1) in day length?}\cr}
#' \item{\code{photodirf }{ =  0, Is the finish daylength trigger during a decrease (0) or increase (1) in day length?}\cr}
#' \item{\code{amphibreed}{ = 0, Amphibious animal breeding mode: 0 is off, 1 is exotrophic aquatic (eggs start when water present in container and within breeding season), 2 is exotrophic terrestrial/aquatic (eggs start at specified soil node within breeding season, diapause at birth threshold, start larval phase if water present in container), 3 endotrophic terrestrial (eggs start at specified soil node within breeding season and continue to metamorphosis on land), 4 turtle mode (eggs start at specified soil node within breeding season, hatch and animals enter water and stay there for the rest of their life, but leave the water if no water is present)}\cr}
#' \item{\code{amphistage}{ = 0, Life cycle control for amphibious animal: 0 is whole life cycle, 1 is just to metamorphosis (then reset and start again)}\cr}
#' \item{\code{reset}{ = 0, Life cycle reset options, 0=quit simulation upon death, 1=restart at emergence, 2=restart at first egg laid, 3=restart at end of breeding season, 4=reset at death}\cr}
#' \item{\code{act_breed}{ = 0, Threshold numbers of hours active after start of breeding season before eggs can be laid (simulating movement to the breeding site)}\cr}
#' \item{\code{rain_breed}{ = 0, Rain dependent breeder? 0 means no, otherwise enter rainfall threshold in mm}\cr}
#' \item{\code{Tb_breed}{ = 200, Body temperature threshold below which breeding will occur}\cr}
#' \item{\code{Tb_breed_hrs}{ = 24*7, Cumulative time below temperature threshold for breeding, hrs \code{Tb_breed} that will trigger breeding}\cr}
#'}
#' \strong{ Mortality rate parameters:}
#' \itemize{
#' \item{\code{m_a}{ = 1e-4, Hourly active mortality rate (probability of mortality per hour)}\cr}
#' \item{\code{m_i}{ = 0, Hourly inactive mortality rate (probability of mortality per hour)}\cr}
#' \item{\code{m_h}{ = 0.5, Survivorship of hatchling in first year}\cr}
#'}
#' \strong{ Water body model parameters (not ready yet):}
#' \itemize{
#' \item{\code{container}{ = 0, Run the water body/container model? (aquatic start of life cycle, e.g. frog or mosquito)}\cr}
#' \item{\code{wetmod}{ = 0, Use precomputed wetland temperature \code{wetlandTemps} and depths \code{wetlandDepthds}?}\cr}
#' \item{\code{conth}{ = 100, Cylindrical container/pond height (mm)}\cr}
#' \item{\code{contw}{ = 1000, Cylindrical container/pond diameter (mm)}\cr}
#' \item{\code{contype}{ = 1, Is 'containter' sitting on the surface, like a bucket (0) or sunk into the ground like a pond (1)}\cr}
#' \item{\code{rainmult}{ = 1, Rainfall multiplier to reflect catchment (don't make this zero unless you want a drought!)}\cr}
#' \item{\code{continit}{ = 0, Initial container water level (cm)}\cr}
#' \item{\code{conthole}{ = 0, Daily (or hourly if rainhr vector filled) loss of height (mm) due to 'hole' in container (e.g. infiltration to soil, drawdown from water tank)}\cr}
#' \item{\code{contonly}{ = 1, Just run the container model and quit?}\cr}
#' \item{\code{contwet}{ = 80, \% of container surface acting as a free water exchanger}\cr}
#' \item{\code{wetlandTemps}{ = matrix(data = 0, nrow = 24 * ndays, ncol = 1), Matrix of hourly wetland temperaures (°C)}\cr}
#' \item{\code{wetlandDepths}{ = matrix(data = 0, nrow = 24 * ndays, ncol = 1), Matrix of hourly wetland depths (cm)}\cr}
#' \item{\code{GLMtemps}{ = matrix(data = 0, nrow = 24 * ndays, ncol = 20), Matrix of hourly wetland temperatures (C) with depth}\cr}
#' \item{\code{GLMO2s}{ = matrix(data = 0, nrow = 24 * ndays, ncol = 20), Matrix of hourly wetland PO2 (kPa) with depth}\cr}
#' \item{\code{GLMsalts}{ = matrix(data = 0, nrow = 24 * ndays, ncol = 20), Matrix of hourly wetland salinity (ppm) with depth}\cr}
#' \item{\code{GLMpHs}{ = matrix(data = 0, nrow = 24 * ndays, ncol = 20), Matrix of hourly wetland pH with depth}\cr}
#' \item{\code{GLMfoods}{ = matrix(data = 0, nrow = 24 * ndays, ncol = 20), Matrix of hourly wetland food density (J/cm3) with depth}\cr}
#' \item{\code{pO2thresh}{ = 10, Oxygen partial pressure tolerance threshold}\cr}
#'}
#' \strong{ Life stage-specific parameter allocation:}
#' \itemize{
#' \item{\code{thermal_stages}{ = matrix(data = c(rep(CT_min, stages), rep(CT_max, stages), rep(T_F_min, stages), rep(T_F_max, stages), rep(T_B_min, stages), rep(T_pref, stages)),  nrow = stages,  ncol = 6), Stage specific thermal thresholds (CT_min, CT_max, T_F_min, T_F_max, T_B_min, T_pref)}\cr}
#' \item{\code{behav_stages}{ = matrix(data = c(rep(diurn, stages), rep(nocturn, stages), rep(crepus, stages), rep(burrow, stages), rep(shdburrow, stages), rep(mindepth, stages), rep(maxdepth, stages), rep(shade_seek, stages), rep(climb, stages), rep(fossorial, stages), rep(rainact, stages), rep(actrainthresh, stages), rep(act_breed, stages), rep(flyer, stages), rep(aquabask, stages)), nrow = stages, ncol = 15), Stage specific behaviour diurn, nocturn, crepus, burrow, shdburrow, mindepth, maxdepth, shade_seek, climb, fossorial, rainact, actrainthresh, act_breed, flyer, aquabask)}\cr}
#' \item{\code{water_stages}{ = matrix(data = c(rep(pct_wet, stages), rep(F_O2, stages), rep(pct_H_P, stages), rep(pct_H_N, stages), rep(pct_H_X, stages), rep(pct_H_R, stages), rep(raindrink, stages), rep(gutfill, stages)),  nrow = stages,  ncol = 8), Stage-specific water budget parameters (pct_wet, F_O2, pct_H_P, pct_H_N, pct_H_X, pct_H_R, raindrink, gutfill)}\cr}
#' \item{\code{nutri_stages}{ = matrix(data = c(rep(foodlim, stages), rep(0, stages)),  nrow = stages,  ncol = 1),  Stage-specific nutritional parameters (foodlim)}\cr}
#' \item{\code{arrhenius}{ = matrix(data = matrix(data = c(rep(T_A, stages), rep(T_AL, stages), rep(T_AH, stages), rep(T_L, stages), rep(T_H, stages)),  nrow = stages,  ncol = 5),  nrow = stages,  ncol = 5), Stage-specific 5-parameter Arrhenius thermal response for DEB model (T_A, T_AL, T_AH, T_L, T_H)}\cr}
#' \item{\code{arrhenius2}{ = matrix(data = matrix(data = c(rep(T_A2, stages), rep(T_AL2, stages), rep(T_AH2, stages), rep(T_L2, stages), rep(T_H2, stages)),  nrow = stages,  ncol = 5),  nrow = stages,  ncol = 5), Stage-specific 5-parameter Arrhenius thermal response for DEB model (T_A, T_AL, T_AH, T_L, T_H) for maturity maintenance (causes 'Temperature Size Rule' effect)}\cr}
#'}
#' \strong{ Butterfly model parameters (not yet tested):}
#' \itemize{
#' \item{\code{wings}{ = 0, Turn wing model on? 1=yes, 0=no}\cr}
#' \item{\code{rho1_3}{ = 0.2, Wing reflectance (0-1)}\cr}
#' \item{\code{trans1}{ = 0.00, Wing transmissivity (0-1)}\cr}
#' \item{\code{aref}{ = 0.26, Width of surface #2 (cm) (back or horizontal or reference surface)}\cr}
#' \item{\code{bref}{ = 2.04, Common length (cm) where the two rectangles join}\cr}
#' \item{\code{cref}{ = 1.47, Width of surface #1 (cm) (wing)}\cr}
#' \item{\code{phi}{ = 179., Initial wing angle (degrees) (90 = vertical relative to body)}\cr}
#' \item{\code{phimax}{ =  phi, Maximum wing angle (degrees) (90 = vertical relative to body)}\cr}
#' \item{\code{phimin}{ =  phi, Minimum wing angle (degrees) (90 = vertical relative to body)}\cr}
#' \item{\code{flyer}{ = 0, Does the animal fly? 1=yes, 0=no}\cr}
#' \item{\code{flyspeed}{ = 5, Flying speed (m/s)}\cr}
#' \item{\code{flymetab}{ = 0.1035, Flight metabolic excess (W/g)}\cr}
#'}
#' \strong{Outputs:}
#'
#' environ variables:
#' \itemize{
#' \item 1 DOY - Day of year
#' \item 2 YEAR - Year of simulation
#' \item 3 DAY - Day of simulation
#' \item 4 TIME - Time of day (hours)
#' \item 5 TC - Body temperature (°C)
#' \item 6 SHADE - Shade selected (\%)
#' \item 7 SOLAR  - Solar radiation (W/m2) at animal location
#' \item 8 DEP - Depth below ground (cm)
#' \item 9 ACT - Activity state (0=inactive, 1=basking, 2=foraging)
#' \item 10 TA - Air temperature (°C) at animal location
#' \item 11 TSUB - air temperature (&deg;C) at animal location
#' \item 12 TSKY - air temperature (&deg;C) at animal location
#' \item 13 VEL - Wind speed (m/s) at animal location
#' \item 14 RELHUM - Relative humidity (\%) at animal location
#' \item 15 ZEN - Zenith angle of sun (degrees - 90 = below the horizon)
#' \item 16 CONDEP - Depth of water body (mm) (may not be simulated or supplied)
#' \item 17 WATERTEMP - Temperature of water body (°C) (may not be simulated or supplied)
#' \item 18 DAYLENGTH - Day length (hours)
#' \item 19 WINGANGLE - Wing angle (degrees) for butterfly model
#' \item 20 WINGTEMP - Wing temperature (°C) for butterfly model
#' \item 21 FLYING - Flying state (1=flying, 0=not flying) for butterfly model
#' \item 22 FLYTIME - Flying time (hours) for butterfly model
#' \item 23 PO2WATER - dissolved oxygen in water, if running GLM water body model
#' \item 24 SALWATER - salinity of water, if running GLM water body model
#' \item 25 ABSAN - solar absorptivity (fractional)
#' \item 26 PCOND - proportion of animal's surface in contact with ground (fractional)
#' \item 27 POSTURE - postural orientation (1=perpendicular to sun, 2=parallel, 0=in-between)
#' }
#' enbal variables:
#' \itemize{
#' \item 1 DOY - Day of year
#' \item 2 YEAR - Year of simulation
#' \item 3 DAY - Day of simulation
#' \item 4 TIME - Time of day (hours)
#' \item 5 QSOL - Solar radiation absorbed (W)
#' \item 6 QIRIN - Infrared radiation absorbed (W)
#' \item 7 QMET - Metabolic heat production (W)
#' \item 8 QEVAP - Evaporative heat loss (W)
#' \item 9 QIROUT - Infrared radiation lost (W)
#' \item 10 QCONV - Heat lost by convection (W)
#' \item 11 QCOND - Heat lost by conduction (W)
#' \item 12 ENB - Energy balance (W)
#' \item 13 NTRY - Iterations that were required for solution to heat balance equation
#'}
#' masbal variables:
#' \itemize{
#' \item 1 DOY - Day of year
#' \item 2 YEAR - Year of simulation
#' \item 3 DAY - Day of simulation
#' \item 4 TIME - Time of day (hours)
#' \item 5 O2_ml - Oxygen consumption rate (ml/h)
#' \item 6 CO2_ml - Carbon dioxide production rate (ml/h)
#' \item 7 NWASTE_g - Nitrogenous waste production (g/h)
#' \item 8 H2OFree_g - Water from food (g/h)
#' \item 9 H2OMet_g - Metabolic water production (g/h)
#' \item 10 DryFood_g - Dry food intake (g/h)
#' \item 11 WetFood_g - Wet foood intake (g/h)
#' \item 12 DryFaeces_g - Dry faeces production (g/h)
#' \item 13 WetFaeces_g - Wet faeces production (g/h)
#' \item 14 Urine_g - Urine production (g/h)
#' \item 15 H2OResp_g - Respiratory water loss (g/h)
#' \item 16 H2OCut_g - Cutaneous water loss (g/h)
#' \item 17 H2OEye_g - Ocular water loss (g/h)
#' \item 18 H2OBal_g - Instantaneous water balance (g/h)
#' \item 19 H2OCumBal_g - Cumulative water balance (g)
#' \item 20 H2OLiq_g - Change in liquid water exchange with substrate (g)
#' \item 21 PSI_kPa - Organism water potential (J/kg = kPa)
#'}
#' debout variables:
#' \itemize{
#' \item 1 DOY - Day of year
#' \item 2 YEAR - Year of simulation
#' \item 3 DAY - Day of simulation
#' \item 4 TIME - Time of day (hours)
#' \item 5 Stage - Life cycle stage (0=embryo, 1=juvenile, 2=mature but not yet reproducing, 3=beyond first reproduction)
#' \item 6 V - Structural volume (cm3)
#' \item 7 E - Reserve density (J/cm3)
#' \item 8 E_H - Maturity state (J)
#' \item 9 L_W - Physical length (mm) (what this represents depends on choice of length measure for DEB paramter fitting, e.g. snout-vent length, head length, etc.)
#' \item 10 WETMASS - Wet mass total (reserve, structure, reproduction buffer, stomach contents) (g)
#' \item 11 WETGONAD - Wet mass of gonad (batch and reproduction buffers) (g)
#' \item 12 WETGUT - Wet mass of food in gut (g)
#' \item 13 PCT_DESIC - \% desiccated
#' \item 14 CUMREPRO - Energy in reproduction buffer (J)
#' \item 15 CUMBATCH - Energy in batch for egg production (J)
#' \item 16 BREEDING - Breeding state (1=breeding, 0=not breeding)
#' \item 17 PREGNANT - Pregnant? (only if viviparous) (0 or 1)
#' \item 18 V_BABY - Structure of baby (cm3) (only if viviparous and pregnant)
#' \item 19 E_BABY - Reserve density of baby (J/cm3) (only if viviparous and pregnant)
#' \item 20 H_S - Hazard rate (1/h)
#' \item 21 P_SURV - Survival probability due to joint influence of ageing and mortality rates
#' \item 22 P_A - assimilation flux, J/h
#' \item 23 P_C - mobilisation flux, J/h
#' \item 24 P_M - maintenance flux, J/h
#' \item 25 P_G - growth flux, J/h
#' \item 26 P_D - dissipation flux, J/h
#' \item 27 P_J - maturity maintenance flux, J/h
#' \item 28 P_R - reproduction/maturation flux, J/h
#' \item 29 P_D - egg flux, J/h
#'}
#' yearout variables:
#' \itemize{
#' \item 1 DEVTIME - Development time (days)
#' \item 2 BIRTHDAY - Birth day (day of year)
#' \item 3 BIRTHMASS - Mass at birth (g)
#' \item 4 MONMATURE - Months to maturity
#' \item 5 LENREPRO - Length (mm) at first reproduction
#' \item 6 FECUNDITY  - Total fecundity
#' \item 7 CLUTCHES - Total clutches
#' \item 8 MINRESERVE - Minimum reserve density (J/cm3)
#' \item 9 LASTFOOD - Food in last year (kg)
#' \item 10 TOTFOOD - Total food eaten (kg)
#' \item 11 MINTB - Minimum body temperature (°C)
#' \item 12 MAXTB - Maxmium body temperature (°C)
#' \item 13 Pct_Des - Maximum level of desiccation (/%)
#' \item 14 LifeSpan - Maximum life span (days)
#' \item 15 GenTime - Generation time (years)
#' \item 16 R0 - Net reproductive rate
#' \item 17 rmax - Intrinsic rate of increase
#' \item 18 LENGTH - Maximum length (mm)
#'}
#'
#' yearsout variables:
#' \itemize{
#' \item 1 YEAR - Year of simulation
#' \item 2 MaxStg - Maximum stage reached in the year
#' \item 3 MaxWgt - Maximum weight reached in the year (g)
#' \item 4 MaxLen - Maximum length in the year (mm)
#' \item 5 Tmax - Maximum annual body temperature (°C)
#' \item 6 Tmin  - Minimum annual body temperature (°C)
#' \item 7 MinRes - Minimum annual reserve density (J/cm3)
#' \item 8 MaxDes - Maximum annual desiccation level (% of normal wet body mass)
#' \item 9 MinShade - Minimum annual shade selected
#' \item 10 MaxShade - Maximum annual shade selected
#' \item 11 MinDep - Minimum annual depth selected (cm)
#' \item 12 MaxDep - Maximum annual depth selected (cm)
#' \item 13 Bsk - Annual basking hours
#' \item 14 Forage - Annual foraging hours
#' \item 15 Dist - Annual distance travelled (flying insect) (km)
#' \item 16 Food - Annual food eaten (g, dry)
#' \item 17 Drink - Annual water drunk (g)
#' \item 18 NWaste - Annual nitrogenous waste (g)
#' \item 19 Faeces - Annual faeces production (g, dry)
#' \item 20 O2 - Annual O2 production (ml)
#' \item 21 Clutch - Annual clutches (#)
#' \item 22 Fec - Annual fecundity (#)
#' \item 23 CauseDeath - Cause of death, 0=no death, 1=cold, 2=heat, 3=desiccation, 4=starvation, 5=ageing
#' \item 24 tLay - Day of year at which eggs laid
#' \item 25 tEgg - Day of year entering egg stage
#' \item 26-33 tStg1-tStg8 - Day of year entering life cycle stages 1-8
#' \item 34-41 mStg1-mStg8 - Body weight upon entering life cycle stages 1-8 (g, wet)
#' \item 42 surviv - Survival probability at end of given year
#' \item 43 deathstage - Life stage at which death occurred
#'}
#' @examples
#'# run the microclimate model
#'micro <- micro_global(loc = c(145.620, -16.821)) #Kuranda, Queensland
#'
#'# retrieve output
#'metout <- as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#'shadmet <- as.data.frame(micro$shadmet) # above ground microclimatic conditions, max shade
#'soil <- as.data.frame(micro$soil) # soil temperatures, minimum shade
#'shadsoil <- as.data.frame(micro$shadsoil) # soil temperatures, maximum shade
#'
#'# append dates
#'dates <- micro$dates
#'metout <- cbind(dates, metout)
#'soil <- cbind(dates, soil)
#'shadmet <- cbind(dates, shadmet)
#'shadsoil <- cbind(dates, shadsoil)
#'
#'# run the ectotherm model
#'ecto <- ectotherm(T_F_min = 30, T_F_max = 35, T_pref = 33, T_B_min = 20, T_RB_min = 10)
#'
#'# retrieve output
#'environ <- as.data.frame(ecto$environ) # activity, Tb and environment
#'enbal <- as.data.frame(ecto$enbal) # energy balance values
#'masbal <- as.data.frame(ecto$masbal) # mass balance value (note most missing if DEB model not running)
#'
#'# append dates
#'environ <- cbind(dates, environ)
#'masbal <- cbind(dates, masbal)
#'enbal <- cbind(dates, enbal)
#'
#'############### plot results ######################
#'
#'# Hourly Tb (black), activity (orange, 5 = bask, 10 = forage), depth (brown, m) and shade (green, %/10)
#'with(environ, plot(TC ~ dates, ylab = "", xlab="month of year", col = 'black', xlim = c(-0.25, 12), ylim = c(-20, 40), type = "l", yaxt = 'n'))
#'with(environ, points(ACT * 2 + 7 ~ dates, type = "p", pch = 16, col = "orange"))
#'with(environ, points(SHADE / 10 - 6 ~ dates, type = "l", col = "dark green"))
#'with(environ, points(DEP - 10 ~ dates, type = "l", col = "brown"))
#'abline(ecto$T_F_min, 0, lty = 2, col = 'blue')
#'abline(ecto$T_F_max, 0, lty = 2, col = 'red')
#'ytick<-seq(15, 40, by=5)
#'axis(side=2, at=ytick, labels = TRUE)
#'mtext(text = c('A', 'B', 'I'), side = 2, line = 1, at = c(11, 9, 7))
#'ytick<-seq(-6, 4, by=2)
#'axis(side=2, at=ytick, labels = FALSE)
#'mtext(text = seq(0, 100, 20), side = 2, line = 1, at = seq(-6, 4, 2), las = 2)
#'ytick<-seq(-20, -10, by=2)
#'axis(side=2, at=ytick, labels = FALSE)
#'mtext(text = rev(seq(0, 100, 20)), side = 2, line = 1, at = seq(-20, -10, 2), las = 2)
#'abline(h = -10, lty = 2, col = 'grey')
#'mtext(text = c('body temperature (°C)', 'activity', 'shade (%)', 'depth (cm)'), side = 2, line = 2.5, at = c(30, 9, 0, -15))
#'text(-0.2, c(ecto$T_F_max + 1, ecto$T_F_min + 1), c('T_F_max', 'T_F_min'), col = c('red', 'blue'), cex = 0.75)
#'
#'# seasonal activity plot (dark blue = night, light blue = basking, orange = foraging)
#'forage <- subset(environ, ACT == 2)
#'bask <- subset(environ, ACT == 1)
#'night <- subset(metout, ZEN == 90)
#'day <- subset(metout, ZEN != 90)
#'with(night, plot(TIME / 60 ~ DOY, ylab = "Hour of Day", xlab = "Day of Year", pch = 15, cex = 2, col = 'dark blue'))
#' # nighttime hours
#'with(forage, points(TIME ~ DOY, pch = 15, cex = 2, col = 'orange')) # foraging Tbs
#'with(bask, points(TIME ~ DOY, pch = 15, cex = 2, col = 'light blue')) # basking Tbs
#' @export
ectotherm <- function(
  Ww_g = 40,
  shape = 3,
  alpha_max = 0.85,
  alpha_min = 0.85,
  T_F_min = 24,
  T_F_max = 34,
  T_B_min = 17.5,
  T_RB_min = 17.5,
  T_pref = 30,
  CT_max = 40,
  CT_min = 6,
  diurn = 1,
  nocturn = 0,
  crepus = 0,
  shade_seek = 1,
  burrow = 1,
  postur = 1,
  climb = 0,
  shdburrow = 0,
  mindepth = 2,
  maxdepth = 10,
  aquabask = 0,
  M_1 = 0.013,
  M_2 = 0.8,
  M_3 = 0.038,
  pct_wet = 0.1,
  pct_eyes = 0.03,
  pct_mouth = 5,
  pantmax = 1,
  F_O2 = 20,
  delta_air = 0.1,
  RQ = 0.8,
  K_skin = 2.8e-09,
  psi_body = -707,
  spec_hyd_body =0.000304,
  K_egg = 2.8e-09,
  psi_egg = -707,
  spec_hyd_egg = 0.000304,
  nyears = micro$nyears,
  enberr = 0.01,
  live = 1,
  write_input = 0,
  transient = 0,
  delta_shade = 3,
  startday = 1,
  minshades = micro$minshade,
  maxshades = micro$maxshade,
  fluid = 0,
  pct_touch = 0,
  O2gas = 20.95,
  CO2gas = 0.03,
  N2gas = 79.02,
  alpha_sub = (1 - micro$REFL),
  PDIF = 0.1,
  DEP = micro$DEP,
  KS = micro$KS,
  b = micro$BB,
  PE = micro$PE,
  metout = micro$metout,
  shadmet = micro$shadmet,
  soil = micro$soil,
  shadsoil = micro$shadsoil,
  soilmoist = micro$soilmoist,
  shadmoist = micro$shadmoist,
  humid = micro$humid,
  shadhumid = micro$shadhumid,
  soilpot = micro$soilpot,
  shadpot = micro$shadpot,
  tcond = micro$tcond,
  shadtcond = micro$shadtcond,
  rainfall = micro$RAINFALL,
  rainhr = rep(-1,nrow(metout)),
  preshr = rep(101325 * ((1 - (0.0065 * as.numeric(micro$elev) / 288)) ^ (1/0.190284)), nrow(metout)),
  elev = as.numeric(micro$elev),
  longitude = as.numeric(micro$longlat[1]),
  latitude = as.numeric(micro$longlat[2]),
  custom_shape = c(10.4713, 0.688, 0.425, 0.85, 3.798, 0.683, 0.694, 0.743),
  shape_a = 1,
  shape_b = 3,
  shape_c = 2 / 3,
  eggshape_a = 1,
  eggshape_b = 2 / 3,
  eggshape_c = 2 / 3,
  fatosk = 0.4,
  fatosb = 0.4,
  rinsul = 0,
  pct_cond = 10,
  pct_cond_egg = 10,
  c_body = 3073,
  k_flesh = 0.5,
  rho_body = 1000,
  epsilon = 0.95,
  epsilon_sub = 1,
  epsilon_sky = 1,
  warmsig = 0,
  fossorial = 0,
  rainact = 0,
  actrainthresh = 0.1,
  soilnode = 4,
  eggshade = 0,
  pO2thresh = 10,
  CT_minthresh = 12,
  CT_kill = 0,
  pct_H_P = 73,
  pct_H_N = 0,
  pct_H_X = 82,
  pct_H_R = 15,
  gutfill = 75,
  raindrink = 0,
  foodlim = 1,
  DEB = 0,
  intmethod = 1,
  z.mult = 1,
  z = 2.825 * z.mult,
  del_M = 0.2144,
  p_Xm = 12420,
  kap_X = 0.85,
  v = 0.02795 / 24,
  kap = 0.8206,
  p_M = 48.81 / 24,
  E_G = 7512,
  kap_R = 0.95,
  k_J = 0.00628 / 24,
  E_Hb = 866.6 * z.mult ^ 3,
  E_Hj = E_Hb * z.mult ^ 3,
  E_Hp = 1.019e+04 * z.mult ^ 3,
  E_He = 1.019e+04 * z.mult ^ 3,
  h_a = 1.051e-08 / (24 ^ 2),
  s_G = 0.01,
  arrhen_mode = 1,
  T_REF = 20 + 273.15,
  T_A = 8817,
  T_AL = 5.0e+04,
  T_AH = 9.0e+04,
  T_L = 6 + 273.15,
  T_H = 33 + 273.15,
  T_A2 = T_A,
  T_AL2 = T_AL,
  T_AH2 = T_AH,
  T_L2 = T_L,
  T_H2 = T_H,
  E_0 = 9220 * z.mult ^ 4,
  f = 1,
  E_sm = 350,
  K = 1,
  X = 10,
  rho_body_deb = rho_body / 1000,
  d_V = 0.3,
  d_E = 0.3,
  d_Egg = 0.3,
  stoich_mode=0,
  mu_X = 525000,
  mu_E = 585000,
  mu_V = 500000,
  mu_P = 480000,
  mu_N = 244e3/5,
  kap_X_P = 0.1,
  n_X = c(1, 1.8,0.5,0.15),
  n_E = c(1, 1.8,0.5,0.15),
  n_V = c(1, 1.8,0.5,0.15),
  n_P = c(1, 1.8,0.5,0.15),
  n_M_nitro = c(1, 4 / 5, 3 / 5, 4 / 5),
  h_N = 384238,
  metab_mode = 0,
  stages = 8,
  S_instar = rep(2.660, stages),
  s_j = 0.999,
  L_b = 0.06148,
  kap_V = 0.8,
  k_Ee = 0.005832307 / 24,
  k_EV = 0.07077021 / 24,
  V_init = 3e-9,
  E_init = E_0 / V_init,
  E_H_init = 0,
  stage = 0,
  aestivate = 0,
  depress = 1,
  clutchsize = 5,
  clutch_ab = c(0, 0),
  eggmult = 0,
  viviparous = 0,
  minclutch = 0,
  batch = 1,
  photostart = 3,
  photofinish = 1,
  daylengthstart = 12.5,
  daylengthfinish = 13,
  photodirs = 1,
  photodirf = 0,
  amphibreed = 0,
  amphistage = 0,
  reset = 0,
  act_breed = 0,
  rain_breed = 0,
  Tb_breed = 200,
  Tb_breed_hrs = 24 * 7,
  m_a = 1e-4,
  m_i = 0,
  m_h = 0.5,
  container = 0,
  wetmod = 0,
  conth = 100,
  contw = 1000,
  contype = 1,
  rainmult = 1,
  continit = 0,
  conthole = 0,
  contonly = 1,
  contwet = 80,
  wetlandTemps = matrix(data = 0, nrow = 24 * ndays, ncol = 1),
  wetlandDepths = matrix(data = 0, nrow = 24 * ndays, ncol = 1),
  GLMtemps = matrix(data = 0, nrow = 24 * ndays, ncol = 20),
  GLMO2s = matrix(data = 10, nrow = 24 * ndays, ncol = 20),
  GLMsalts = matrix(data = 0, nrow = 24 * ndays, ncol = 20),
  GLMpHs = matrix(data = 7, nrow = 24 * ndays, ncol = 20),
  GLMfoods = matrix(data = 10, nrow = 24 * ndays, ncol = 20),
  thermal_stages = matrix(data = c(rep(CT_min, stages), rep(CT_max, stages), rep(T_F_min, stages), rep(T_F_max, stages), rep(T_B_min, stages),
                                 rep(T_pref, stages)), nrow = stages, ncol = 6),
  behav_stages = matrix(data = c(rep(diurn, stages), rep(nocturn, stages), rep(crepus, stages), rep(burrow, stages),
                               rep(shdburrow, stages), rep(mindepth, stages), rep(maxdepth, stages), rep(shade_seek, stages), rep(climb, stages), rep(fossorial, stages),
                               rep(rainact, stages), rep(actrainthresh, stages), rep(act_breed, stages), rep(flyer, stages), rep(aquabask, stages)), nrow = stages, ncol = 15),
  water_stages  =  matrix(data = c(rep(pct_wet, stages), rep(F_O2, stages), rep(pct_H_P, stages), rep(pct_H_N, stages),
                               rep(pct_H_X[1], stages), rep(pct_H_R, stages), rep(raindrink, stages), rep(gutfill, stages)), nrow = stages, ncol = 8),
  nutri_stages = matrix(data = c(rep(foodlim, stages)),  nrow = stages,  ncol = 1),
  arrhenius = matrix(data = matrix(data = c(rep(T_A, stages), rep(T_AL, stages), rep(T_AH, stages), rep(T_L, stages), rep(T_H, stages)),
                                 nrow = stages, ncol = 5), nrow = stages, ncol = 5),
  arrhenius2 = matrix(data = matrix(data = c(rep(T_A2, stages), rep(T_AL2, stages), rep(T_AH2, stages), rep(T_L2, stages), rep(T_H2, stages)),
                                   nrow = stages, ncol = 5), nrow = stages, ncol = 5),
  wings = 0,
  rho1_3 = 0.2,
  trans1 = 0,
  aref = 0.26,
  bref = 2.04,
  cref = 1.47,
  phi = 179,
  phimax = phi,
  phimin = phi,
  flyer = 0,
  flyspeed = 5,
  flymetab = 0.1035,
  pct_H_death = 35,
  write_csv = 0,
  aestdepth = 7){ # end function parameters

  errors <- 0
  ndays <- length(rainfall) # get number of days of simulation

  # error trapping
  if(shape < 0 | shape > 5  | shape%%1 != 0){
    message("error: shape can only be an integer from 0 to 5 \n")
    errors<-1
  }
  if(alpha_max < 0 | alpha_max > 1){
    message("error: alpha_max can only be from 0 to 1 \n")
    errors<-1
  }
  if(alpha_min < 0 | alpha_min > 1){
    message("error: alpha_min can only be from 0 to 1 \n")
    errors<-1
  }
  if(T_F_min > T_F_max){
    message("error: T_F_min must be less than T_F_max \n")
    errors<-1
  }
  if(T_B_min > T_F_min){
    message("error: T_B_min must be less than or equal to T_F_min \n")
    errors<-1
  }
  if(T_RB_min > T_F_min){
    message("error: T_RB_min must be less than or equal to T_F_min \n")
    errors<-1
  }
  if(T_RB_min > T_B_min){
    message("error: T_RB_min must be less than or equal to T_B_min \n")
    errors<-1
  }
  if(!diurn %in% c(0,1)){
    message("error: diurn must be 0 or 1 \n")
    errors<-1
  }
  if(!nocturn %in% c(0,1)){
    message("error: nocturn must be 0 or 1 \n")
    errors<-1
  }
  if(!crepus %in% c(0,1)){
    message("error: crepus must be 0 or 1 \n")
    errors<-1
  }
  if(!shade_seek %in% c(0,1)){
    message("error: shade_seek must be 0 or 1 \n")
    errors<-1
  }
  if(!burrow %in% c(0,1)){
    message("error: burrow must be 0 or 1 \n")
    errors<-1
  }
  if(!postur %in% c(0,1,2)){
    message("error: postur must be 0, 1 or 2 \n")
    errors<-1
  }
  if(!climb %in% c(0,1)){
    message("error: climb must be 0 or 1 \n")
    errors<-1
  }
  if(!shdburrow %in% c(0,1,2)){
    message("error: shdburrow must be 0, 1 or 2 \n")
    errors<-1
  }
  if(!mindepth %in% seq(1,10)){
    message("error: mindepth must be an integer between 1 and 10 \n")
    errors<-1
  }
  if(!maxdepth %in% seq(1,10)){
    message("error: maxdepth must be an integer between 1 and 10 \n")
    errors<-1
  }
  if(mindepth > maxdepth){
    message("error: mindepth must be less than or equal to maxdepth \n")
    errors<-1
  }
  if(!aquabask %in% c(0,1,2)){
    message("error: aquabask must be 0, 1 or 2 \n")
    errors<-1
  }
  if(pct_wet < 0 | pct_wet > 100){
    message("error: pct_wet can only be from 0 to 100 \n")
    errors<-1
  }
  if(pct_eyes < 0 | pct_eyes > 100){
    message("error: pct_eyes can only be from 0 to 100 \n")
    errors<-1
  }
  if(pct_mouth < 0 | pct_mouth > 100){
    message("error: pct_mouth can only be from 0 to 100 \n")
    errors<-1
  }
  if((pantmax < 0) | (pantmax > 0 & pantmax < 1)){
    message("error: pantmax should be greater than or equal to 1, or zero if you want to simluate the effect of no respiratory water loss\n")
    errors<-1
  }
  if(F_O2 < 0 | F_O2 > 100){
    message("error: F_O2 can only be from 0 to 100 \n")
    errors<-1
  }
  if(nyears < 1){
    message("error: nyears must be greater than or equal to 1 \n")
    errors<-1
  }
  if(enberr <= 0){
    message("error: enberr must be greater than 0 \n")
    errors<-1
  }
  if(!live %in% c(0,1)){
    message("error: live must be 0 or 1 \n")
    errors<-1
  }
  if(!write_input %in% c(0,1,2)){
    message("error: write_input must be 0, 1 or 2 \n")
    errors<-1
  }
  if(!transient %in% c(0,1)){
    message("error: transient must be 0 or 1 \n")
    errors<-1
  }
  if(delta_shade <= 0){
    message("error: delta_shade must be greater than 0 \n")
    errors<-1
  }
  if(startday < 1){
    message("error: startday must be greater than or equal to 1 \n")
    errors<-1
  }
  if(min(minshades) < 0 | max(minshades) > 100){
    message("error: minshades can only be from 0 to 100 \n")
    errors<-1
  }
  if(min(maxshades) < 0 | max(maxshades) > 100){
    message("error: maxshades can only be from 0 to 100 \n")
    errors<-1
  }
  if(length(maxshades) != ndays){
    message("error: maxshades must be a vector with a length equal to the number of days simulated \n")
    errors<-1
  }
  if(length(minshades) != ndays){
    message("error: minshades must be a vector with a length equal to the number of days simulated \n")
    errors<-1
  }
  if(!fluid %in% c(0,1)){
    message("error: fluid must be 0 or 1 \n")
    errors<-1
  }
  if(alpha_sub < 0 | alpha_sub > 1){
    message("error: alpha_sub can only be from 0 to 1 \n")
    errors<-1
  }
  if(length(DEP) != 10){
    message("error: DEP must be a vector of length 10 \n")
    errors<-1
  }
  if(ncol(metout) < 19){
    message("error: metout must have 19 columns \n")
    errors<-1
  }
  if(ncol(shadmet) < 19){
    message("error: shadmet must have 19 columns \n")
    errors<-1
  }
  if(ncol(soil) < 12){
    message("error: soil must have 12 columns \n")
    errors<-1
  }
  if(ncol(shadsoil) < 12){
    message("error: shadsoil must have 12 columns \n")
    errors<-1
  }
  if(ncol(soilmoist) < 12){
    message("error: soilmoist must have 12 columns \n")
    errors<-1
  }
  if(ncol(shadmoist) < 12){
    message("error: shadmoist must have 12 columns \n")
    errors<-1
  }
  if(ncol(humid) < 12){
    message("error: humid must have 12 columns \n")
    errors<-1
  }
  if(ncol(shadhumid) < 12){
    message("error: shadhumid must have 12 columns \n")
    errors<-1
  }
  if(ncol(soilpot) < 12){
    message("error: soilpot must have 12 columns \n")
    errors<-1
  }
  if(ncol(shadpot) < 12){
    message("error: shadpot must have 12 columns \n")
    errors<-1
  }
  if(ncol(tcond) < 12){
    message("error: tcond must have 12 columns \n")
    errors<-1
  }
  if(ncol(shadtcond) < 12){
    message("error: shadtcond must have 12 columns \n")
    errors<-1
  }
  if(min(rainfall) < 0){
    message("error: rainfall contains some negative values \n")
    errors<-1
  }
  if(max(rainhr) > 0 & min(rainhr) < 0){
    message("error: rainhr contains some negative values \n")
    errors<-1
  }
  if(longitude < -180 | longitude > 180){
    message("error: longitude must be between -180 and 180 \n")
    errors<-1
  }
  if(latitude < -90 | latitude > 90){
    message("error: latitude must be between -90 and 90 \n")
    errors<-1
  }
  if(shape_a < 0){
    message("error: shape_a can't be negative \n")
    errors<-1
  }
  if(shape_b < 0){
    message("error: shape_b can't be negative \n")
    errors<-1
  }
  if(shape_c < 0){
    message("error: shape_c can't be negative \n")
    errors<-1
  }
  if(eggshape_a < 0){
    message("error: shape_a can't be negative \n")
    errors<-1
  }
  if(eggshape_b < 0){
    message("error: shape_b can't be negative \n")
    errors<-1
  }
  if(eggshape_c < 0){
    message("error: shape_c can't be negative \n")
    errors<-1
  }
  if(fatosk < 0 | fatosk > 1){
    message("error: fatosk can only be from 0 to 1 \n")
    errors<-1
  }
  if(fatosb < 0 | fatosb > 1){
    message("error: fatosb can only be from 0 to 1 \n")
    errors<-1
  }
  if(fatosk + fatosb > 1){
    message("error: the sum of fatosb and fatosk can't exceed 1 \n")
    errors<-1
  }
  if(rinsul < 0){
    message("error: rinsul can't be negative \n")
    errors<-1
  }
  if(pct_cond < 0 | pct_cond > 100){
    message("error: pct_cond can only be from 0 to 100 \n")
    errors<-1
  }
  if(pct_cond_egg < 0 | pct_cond_egg > 100){
    message("error: pct_cond_egg can only be from 0 to 100 \n")
    errors<-1
  }
  if(c_body < 0){
    message("error: c_body can't be negative \n")
    errors<-1
  }
  if(k_flesh < 0){
    message("error: k_flesh can't be negative \n")
    errors<-1
  }
  if(rho_body < 0){
    message("error: rho_body can't be negative \n")
    errors<-1
  }
  if(epsilon < 0 | epsilon > 1){
    message("error: epsilon can only be from 0 to 1 \n")
    errors<-1
  }
  if(epsilon < 0.9){
    message("warning: epsilon is rarely below 0.9 for living things \n")
    errors<-0
  }
  if(!fossorial %in% c(0,1)){
    message("error: fossorial must be 0 or 1 \n")
    errors<-1
  }
  if(!rainact %in% c(0,1)){
    message("error: rainact must be 0 or 1 \n")
    errors<-1
  }
  if(!soilnode %in% seq(1,10)){
    message("error: soilnode must be an integer between 1 and 10 \n")
    errors<-1
  }
  if(!eggshade %in% c(0,1)){
    message("error: eggshade must be 0 or 1 \n")
    errors<-1
  }
  if(pO2thresh < 0){
    message("error: pO2thresh can't be negative \n")
    errors<-1
  }
  if(!eggshade %in% c(0,1)){
    message("error: eggshade must be 0 or 1 \n")
    errors<-1
  }
  if(!CT_kill %in% c(0,1)){
    message("error: CT_kill must be 0 or 1 \n")
    errors<-1
  }
  if(pct_H_P < 0 | pct_H_P > 100){
    message("error: pct_H_P can only be from 0 to 100 \n")
    errors<-1
  }
  if(pct_H_N < 0 | pct_H_N > 100){
    message("error: pct_H_N can only be from 0 to 100 \n")
    errors<-1
  }
  if(min(pct_H_X) < 0 | max(pct_H_X) > 100){
    message("error: pct_H_X can only be from 0 to 100 \n")
    errors<-1
  }
  if(pct_H_R < 0 | pct_H_R > 100){
    message("error: pct_H_R can only be from 0 to 100 \n")
    errors<-1
  }
  if(gutfill < 0 | gutfill > 100){
    message("error: gutfill can only be from 0 to 100 \n")
    errors<-1
  }
  if(!foodlim %in% c(0,1)){
    message("error: foodlim must be 0 or 1 \n")
    errors<-1
  }
  if(!DEB %in% c(0,1)){
    message("error: DEB must be 0 or 1 \n")
    errors<-1
  }
  if(!metab_mode %in% c(0,1,2)){
    message("error: metab_mode must be 0, 1 or 2 \n")
    errors<-1
  }
  if(!aestivate %in% c(0,1)){
    message("error: aestivate must be 0 or 1 \n")
    errors<-1
  }
  if(!viviparous %in% c(0,1)){
    message("error: viviparous must be 0 or 1 \n")
    errors<-1
  }
  if(!batch %in% c(0,1)){
    message("error: batch must be 0 or 1 \n")
    errors<-1
  }
  if(!photostart %in% seq(0,5)){
    message("error: photostart must be an integer between 0 and 5 \n")
    errors<-1
  }
  if(!photofinish %in% seq(0,5)){
    message("error: photofinish must be an integer between 0 and 5 \n")
    errors<-1
  }
  if(!photodirs %in% c(0,1)){
    message("error: photodirs must be 0 or 1 \n")
    errors<-1
  }
  if(!amphibreed %in% c(0,1)){
    message("error: amphibreed must be 0 or 1 \n")
    errors<-1
  }
  if(!amphistage %in% c(0,1)){
    message("error: amphistage must be 0 or 1 \n")
    errors<-1
  }
  if(!wetmod %in% c(0,1)){
    message("error: wetmod must be 0 or 1 \n")
    errors<-1
  }
  if(!contype %in% c(0,1)){
    message("error: contype must be 0 or 1 \n")
    errors<-1
  }
  if(conthole < 0){
    message("error: conthole must be >= 0 \n")
    errors<-1
  }
  if(!container %in% c(0,1)){
    message("error: container must be 0 or 1 \n")
    errors<-1
  }
  if(!contonly %in% c(0,1)){
    message("error: contonly must be 0 or 1 \n")
    errors<-1
  }
  if(rainmult < 1){
    message("warning: rainfall is being reduced from original values because rainmult < 1")
  }
  if(contwet < 0 | contwet > 100){
    message("error: contwet can only be from 0 to 100 \n")
    errors<-1
  }
  if(!wings %in% c(0,1)){
    message("error: wings must be 0 or 1 \n")
    errors<-1
  }
  if(phi < 0 | phi > 180){
    message("error: phi can only be from 0 to 180 \n")
    errors<-1
  }
  if(phimax < 0 | phimax > 180){
    message("error: phimax can only be from 0 to 180 \n")
    errors<-1
  }
  if(phimin < 0 | phimin > 180){
    message("error: phimin can only be from 0 to 180 \n")
    errors<-1
  }
  if(phimin < phimax){
    message("error: phimin must be less than or equal to phimax \n")
    errors<-1
  }
  if(!flyer %in% c(0,1)){
    message("error: flyer must be 0 or 1 \n")
    errors<-1
  }
  if(!write_csv %in% c(0,1,2)){
    message("error: write_csv must be 0, 1 or 2 \n")
    errors<-1
  }
  if(pct_H_death < 0 | pct_H_death > 100){
    message("error: pct_H_death can only be from 0 to 100 \n")
    errors<-1
  }
  if(!aestdepth %in% seq(1,10)){
    message("error: aestdepth must be an integer between 1 and 10 \n")
    errors<-1
  }

  if(shape == 3){ # lizard proportions
    shape_a <- 1
    shape_b <- 1
    shape_c <- 4
  }

  if(shape == 4){ # frog proportions
    shape_a <- 1
    shape_b <- 1
    shape_c <- 0.5
  }

  # turn on container model if aquatic egg/larval phase
  if(amphibreed == 1 | amphibreed == 2){
    container <- 1
  }
  if(amphibreed == 3){
    container <- 0
  }
  if(errors == 0){

    #initializing

    DOYstart <- metout[1, 2] # starting day of year
    DOY <- 1 # day of year at start
    iyear <- 0 # initializing year counter
    countday <- 1 # initializing day counter
    # container/pond initial conditions
    contlast <- 0 # last container depth, cm
    templast <- 7 # last container temperature, deg C

    tannul <- as.numeric(mean(soil[, 12])) # annual mean temperature, deg C
    tcinit <- metout[1, "TALOC"] # initial temperature for transient heat budget

    # parameter name translations

    lat <- latitude # latitude

    # habitat
    ALT <- elev # altitude (m)
    EMISSK <- epsilon_sky # emissivity of the sky (0-1)
    EMISSB <- epsilon_sub # emissivity of the substrate (0-1)
    ABSSB <- alpha_sub # solar absorbtivity of the substrate (0-1)

    # animal properties
    Ww_kg <- Ww_g / 1000 # animal wet weight (kg)
    absan <- alpha_max # animal solar absorbtivity
    SKINW <- pct_wet # skin wetness %
    o2max <- F_O2 # O2 extraction efficiency

    # conversions from percent to proportion
    skint <- pct_touch / 100
    PTUREA1 <- pct_H_N / 100
    PFEWAT1 <- pct_H_P / 100
    pct_H_X <- pct_H_X / 100
    water_stages[, 5] <- water_stages[, 5] / 100 # pct_H_X
    FoodWater1 <- pct_H_X[1]
    water_stages[,3] <- water_stages[, 3] / 100
    water_stages[,4] <- water_stages[, 4] / 100
    water_stages[,5] <- water_stages[, 5] / 100

    # DEB mass/stoichiometry and entropy/heat calculations
    if(stoich_mode == 0){
      # match H fraction in organics to stated chemical potentials (needed later for heat production)
      n_X[2] <- ((mu_X / 10 ^ 5) - 4.3842 * n_X[1] - (-1.8176) * n_X[3] - (0.0593) * n_X[4]) / 0.9823
      n_V[2] <- ((mu_V / 10 ^ 5) - 4.3842 * n_V[1] - (-1.8176) * n_V[3] - (0.0593) * n_V[4]) / 0.9823
      n_E[2] <- ((mu_E / 10 ^ 5) - 4.3842 * n_E[1] - (-1.8176) * n_E[3] - (0.0593) * n_E[4]) / 0.9823
      n_P[2] <- ((mu_P / 10 ^ 5) - 4.3842 * n_P[1] - (-1.8176) * n_P[3] - (0.0593) * n_P[4]) / 0.9823
    }else{
      if(stoich_mode == 1){
        # match stated chemical potentials to H fraction in organics
        mu_X <- (n_X[2] * 0.9823 + 4.3842 * n_X[1] + (-1.8176) * n_X[3] + (0.0593) * n_X[4]) * 10 ^ 5
        mu_V <- (n_V[2] * 0.9823 + 4.3842 * n_V[1] + (-1.8176) * n_V[3] + (0.0593) * n_V[4]) * 10 ^ 5
        mu_E <- (n_E[2] * 0.9823 + 4.3842 * n_E[1] + (-1.8176) * n_E[3] + (0.0593) * n_E[4]) * 10 ^ 5
        mu_P <- (n_P[2] * 0.9823 + 4.3842 * n_P[1] + (-1.8176) * n_P[3] + (0.0593) * n_P[4]) * 10 ^ 5
      }
    }
    # enthalpies (combustion frame)
    h_X <- 10^5 * (4.3284 * n_X[1] + 1.0994 * n_X[2] + (-2.0915) * n_X[3] + (-0.1510) * n_X[4]) #J mol^(-1)
    h_V <- 10^5 * (4.3284 * n_V[1] + 1.0994 * n_V[2] + (-2.0915) * n_V[3] + (-0.1510) * n_V[4]) #J mol^(-1)
    h_E <- 10^5 * (4.3284 * n_E[1] + 1.0994 * n_E[2] + (-2.0915) * n_E[3] + (-0.1510) * n_E[4]) #J mol^(-1)
    h_P <- 10^5 * (4.3284 * n_P[1] + 1.0994 * n_P[2] + (-2.0915) * n_P[3] + (-0.1510) * n_P[4]) #J mol^(-1)
    h_CO2 <- 0 #J mol^(-1)
    h_O2 <- 0 #J mol^(-1)
    h_H2O <- 0 #J mol^(-1)
    if(all(n_M_nitro == c(0, 3, 0, 1))){ # ammonia
      h_N <- 382805
      mu_N <- 0
    }
    if(all(n_M_nitro == c(1.0, 0.8, 0.6, 0.8))){ # uric acid
      h_N <- 384238
      mu_N <- 244e3/5
    }
    if(all(n_M_nitro == c(1, 2, 1, 2))){ # urea
      h_N <- 631890
      mu_N <- 122e3
    }
    h_O <- c(h_X, h_V, h_E, h_P)
    h_M <- c(h_CO2, h_H2O, h_O2, h_N)
    n_O <- cbind(n_X, n_V, n_E, n_P) # matrix of C-mole composition of organics, i.e. food, structure, reserve and faeces
    CHON <- c(12, 1, 16, 14) # molar masses of carbon, hydrogen, oxygen and nitrogen, g/mol
    wO <- CHON %*% n_O # molar weight of organics, g/mol
    w_V <- wO[3] # molar weight of structure, g/mol
    M_V <- d_V / w_V # molar mass of structure, mol/cm3
    y_EX <- kap_X * mu_X / mu_E # yield of reserve on food, mol/mol
    y_XE <- 1 / y_EX # yield of food on reserve, mol/mol
    y_VE <- mu_E * M_V / E_G  # yield of structure on reserve, mol/mol
    y_PX <- kap_X_P * mu_X / mu_P # yield of faeces on food, mol/mol
    y_PE <- y_PX / y_EX # yield of faeces on reserve, mol/mol
    nM <- matrix(c(1, 0, 2, 0, 0, 2, 1, 0, 0, 0, 2, 0, n_M_nitro), nrow = 4)
    n_M_nitro_inv <- c(-1 * n_M_nitro[1] / n_M_nitro[4], (-1 * n_M_nitro[2]) / (2 * n_M_nitro[4]), (4 * n_M_nitro[1] + n_M_nitro[2] - 2 * n_M_nitro[3]) / (4 * n_M_nitro[4]), 1 / n_M_nitro[4])
    n_M_inv <- matrix(c(1, 0, -1, 0, 0, 1 / 2, -1 / 4, 0, 0, 0, 1 / 2, 0, n_M_nitro_inv), nrow = 4)
    JM_JO <- -1 * n_M_inv %*% n_O
    eta_O <- matrix(c(y_XE / mu_E * -1, 0, 1 / mu_E, y_PE / mu_E, 0, 0, -1 / mu_E, 0, 0, y_VE / mu_E, -1 / mu_E, 0), nrow = 4)
    w_N <- CHON %*% n_M_nitro
    E_m <- (p_M * z / kap) / v # maximum reserve density, J/cm3

    # DEB model initial conditions
    V_init_baby <- 3e-9 # initial structure, cm3
    E_init_baby <- E_0 / V_init_baby # initial reserve density, J/cm3
    E_baby_init <- E_init_baby #
    V_baby_init <- V_init_baby
    ES_init <- 0 # intial stomach energy, J
    cumrepro_init <- 0 # initial reproductive energy, J
    cumbatch_init <- 0 #initial reproduction batch energy
    q_init <- 0 # initial surivival probability
    hs_init <- 0 # specific death probability, 1/t
    pregnant <- 0 # initial pregnancy state

    # food and food water levels
    if(length(X) == 1){ # no day-specific food levels given
      foodlevels <- rep(X, nrow(metout) / 24)
    }else{
      foodlevels <- X
    }
    if(length(pct_H_X) == 1){ # no day-specific food water levels given
      foodwaters <- rep(pct_H_X, nrow(metout) / 24)
    }else{
      foodwaters <- pct_H_X
    }

    # unused /spare parameters
    tester <- 0 # unused
    microyear <- 1 # extraneous, not used
    nodnum <- 10 # depth at which foraging occurs in fossorial species, probably not working properly, may not need it
    OBJDIS <- 1.0 # currently unused - distance (m) from nearby object of different temp to sky and ground (e.g. warm rock, fire)
    OBJL <- 0.0001 # currently unused - diameter (m) of nearby object of different temp to sky and ground (e.g. warm rock, fire)
    FATOBJ <- 0 # configuration factor to nearby object of different temp to sky and ground (e.g. warm rock, fire)
    SPARE4 <- 1 # spare input
    SPARE2 <- 1 # spare input
    SPARE3 <- 0 # spare input

    # collate parameters
    gas <- c(O2gas, CO2gas, N2gas) # gas vector
    behav <- c(diurn, nocturn, crepus, rainact, burrow, shade_seek, climb, fossorial, SPARE3) # behaviour vector
    ectoinput <- as.matrix(c(ALT, fluid, OBJDIS, OBJL, PDIF, EMISSK, EMISSB, ABSSB, K_skin, enberr, Ww_kg, epsilon, absan, RQ, rinsul, shape, live, pantmax, k_flesh, c_body, rho_body, alpha_max, alpha_min, fatosk, fatosb, FATOBJ, T_F_max, T_F_min, delta_air, SKINW, pct_eyes, pct_mouth, F_O2, T_pref, pct_cond, skint, gas, transient, soilnode, o2max, SPARE4, tannul, nodnum, postur, psi_body, spec_hyd_body, CT_max, CT_min, behav, DOY, actrainthresh, viviparous, pregnant, conth, contw, contlast, arrhen_mode, tcinit, nyears, lat, rainmult, DOYstart, delta_shade, custom_shape, M_1, M_2, M_3, DEB, tester, rho1_3, trans1, aref, bref, cref, phi, wings, phimax, phimin, shape_a, shape_b, shape_c, pct_H_R, microyear, container, flyer, flyspeed, ndays, maxdepth, CT_minthresh, CT_kill, gutfill, mindepth, T_B_min, T_RB_min, p_Xm, eggmult, flymetab, continit, wetmod, contonly, conthole, contype, shdburrow, Tb_breed, Tb_breed_hrs, contwet, warmsig, aquabask, pct_H_death, write_csv, aestdepth, eggshade, pO2thresh, intmethod, eggshape_a, eggshape_b, eggshape_c, pct_cond_egg, K_egg, psi_egg, spec_hyd_egg, b, KS, PE))
    debmod <- c(clutchsize, rho_body_deb, d_V, d_Egg, mu_X, mu_E, mu_V, mu_P, T_REF - 273.15, z, kap, kap_X, p_M, v, E_G, kap_R, E_sm, del_M, h_a, V_init_baby, E_init_baby, k_J, E_Hb, E_Hj, E_Hp, clutch_ab[2], batch, rain_breed, photostart, photofinish, daylengthstart, daylengthfinish, photodirs, photodirf, clutch_ab[1], amphibreed, amphistage, eta_O, JM_JO, E_0, kap_X_P, PTUREA1, PFEWAT1, wO, w_N, FoodWater1, f, s_G, K, X[1], metab_mode, stages, kap_V, s_j, startday, raindrink, reset, m_a, m_i, m_h, aestivate, depress, minclutch, L_b, E_He, k_Ee, k_EV, mu_N, h_O, h_M[4])
    deblast <- c(iyear, countday, V_init, E_init, ES_init, cumrepro_init, q_init, hs_init, cumbatch_init, V_baby_init, E_baby_init, E_H_init, stage)

    # code to determine wet periods for activity in a pond
    if(wetmod==1){
      wet_thresh <- 10 * 24 # threshold pond duration
      wet_depth <- 100 # threshold pond depth (mm)
      wet_temp <- 28 # threshold exit temp (°C)
      b <- cbind(as.data.frame(wetlandDepths), as.data.frame(wetlandTemps))
      colnames(b) <- c('depth', 'temp')
      b$depth[b$temp > wet_temp] <- 0
      b <- b$depth
      b[b >= wet_depth] <- 1
      b[b != 1] <- 0
      bb <- rle(b)
      bb$values[bb$lengths < wet_thresh] <- 0
      c <- b * 0
      values <- bb$values
      lengths <- bb$lengths
      for(k in 1:length(bb$values)){
        d <- c(rep(values[k], lengths[k]))
        if(k == 1){
          e <- d
        }else{
          e <- c(e, d)
        }
      }
      wetlandDepths <- wetlandDepths * e
    }

    if(write_input == 1){ # write out input as csv files for debugging
      if(dir.exists("ecto csv input") == FALSE){
        dir.create("ecto csv input")
      }
      message('writing input csv files \n')
      write.csv(ectoinput, file = "ecto csv input/ectoinput.csv")
      write.csv(debmod, file = "ecto csv input/debmod.csv")
      write.csv(deblast, file = "ecto csv input/deblast.csv")
      write.csv(rainfall, file = "ecto csv input/rainfall.csv")
      write.csv(rainhr, file = "ecto csv input/rainhr.csv")
      write.csv(preshr, file = "ecto csv input/preshr.csv")
      write.csv(DEP, file = "ecto csv input/dep.csv")
      write.csv(foodwaters, file = "ecto csv input/foodwaters.csv")
      write.csv(foodlevels, file = "ecto csv input/foodlevels.csv")
      write.csv(wetlandTemps, file = "ecto csv input/wetlandTemps.csv")
      write.csv(wetlandDepths, file = "ecto csv input/wetlandDepths.csv")
      write.csv(GLMtemps, file = "ecto csv input/GLMtemps.csv", row.names = F)
      write.csv(GLMO2s, file = "ecto csv input/GLMO2s.csv", row.names = F)
      write.csv(GLMsalts, file = "ecto csv input/GLMsalts.csv", row.names = F)
      write.csv(GLMpHs, file = "ecto csv input/GLMpHs.csv", row.names = F)
      write.csv(GLMfoods, file = "ecto csv input/GLMfoods.csv", row.names = F)
      write.csv(arrhenius, file = "ecto csv input/arrhenius.csv")
      write.csv(arrhenius, file = "ecto csv input/arrhenius2.csv")
      write.csv(thermal_stages, file = "ecto csv input/thermal_stages.csv")
      write.csv(behav_stages, file = "ecto csv input/behav_stages.csv")
      write.csv(water_stages, file = "ecto csv input/water_stages.csv")
      write.csv(nutri_stages, file = "ecto csv input/nutri_stages.csv")
      write.csv(minshades, file = "ecto csv input/Minshades.csv")
      write.csv(maxshades, file = "ecto csv input/Maxshades.csv")
      write.csv(S_instar, file = "ecto csv input/S_instar.csv")
      write.table(metout[(seq(1, ndays * 24)), ], file = "ecto csv input/metout.csv", sep = ",", row.names = FALSE)
      write.table(shadmet[(seq(1, ndays * 24)), ], file = "ecto csv input/shadmet.csv", sep = ",", row.names = FALSE)
      write.table(soil[(seq(1, ndays * 24)), ], file = "ecto csv input/soil.csv", sep = ",", row.names = FALSE)
      write.table(shadsoil[(seq(1, ndays * 24)), ], file = "ecto csv input/shadsoil.csv", sep = ",", row.names = FALSE)
      write.table(soilmoist[(seq(1, ndays * 24)), ], file = "ecto csv input/soilmoist.csv", sep = ",", row.names = FALSE)
      write.table(shadmoist[(seq(1, ndays * 24)), ], file = "ecto csv input/shadmoist.csv", sep = ",", row.names = FALSE)
      write.table(soilpot[(seq(1, ndays * 24)), ], file = "ecto csv input/soilpot.csv", sep = ",", row.names = FALSE)
      write.table(shadpot[(seq(1, ndays * 24)), ], file = "ecto csv input/shadpot.csv", sep = ",", row.names = FALSE)
      write.table(humid[(seq(1, ndays * 24)), ], file = "ecto csv input/humid.csv", sep = ",", row.names = FALSE)
      write.table(shadhumid[(seq(1, ndays * 24)), ], file = "ecto csv input/shadhumid.csv", sep = ",", row.names = FALSE)
      write.table(tcond[(seq(1, ndays * 24)), ], file = "ecto csv input/tcond.csv", sep = ",", row.names = FALSE)
      write.table(shadtcond[(seq(1, ndays * 24)), ], file = "ecto csv input/shadtcond.csv", sep = ",", row.names = FALSE)
    }
    # final input list
    ecto <- list(ndays = ndays, nstages = stages, ectoinput = ectoinput, metout = metout[, 1:18], shadmet = shadmet[, 1:18], soil = soil, shadsoil = shadsoil, soilmoist = soilmoist, shadmoist = shadmoist, soilpot = soilpot, shadpot = shadpot, humid = humid, shadhumid = shadhumid, tcond = tcond, shadtcond = shadtcond, DEP = DEP, rainfall = rainfall, rainhr = rainhr, preshr = preshr, iyear = iyear, countday = countday, debmod = debmod, deblast = deblast, foodwaters = foodwaters, foodlevels = foodlevels, wetlandTemps = wetlandTemps, wetlandDepths = wetlandDepths, GLMtemps = GLMtemps, GLMO2s = GLMO2s, GLMsalts = GLMsalts, GLMpHs = GLMpHs, GLMfoods = GLMfoods, arrhenius = arrhenius, arrhenius2 = arrhenius2, thermal_stages = thermal_stages, behav_stages = behav_stages, water_stages = water_stages, nutri_stages = nutri_stages, minshades = minshades, maxshades = maxshades, S_instar = S_instar)

    message('running ectotherm model ... \n')

    ptm <- proc.time() # Start timing
    ectout <- ectorun(ecto) # call Fortran
    message(paste0('runtime ', (proc.time() - ptm)[3], ' seconds')) # Stop the clock

    environ <- ectout$environ[1:(ndays * 24), ]
    enbal <- ectout$enbal[1:(ndays * 24), ]
    masbal <- ectout$masbal[1:(ndays * 24), ]
    debout <- ectout$debout[1:(ndays * 24), ]
    yearout <- ectout$yearout
    yearsout <- ectout$yearsout

    if(DEB==0){
      return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,soilpot=soilpot,shadpot=shadpot,humid=humid,shadhumid=shadhumid,tcond=tcond,shadtcond=shadtcond,rainfall=rainfall,rainhr=rainhr,enbal=enbal,environ=environ,masbal=masbal,yearout=yearout,yearsout=yearsout,foodwaters=foodwaters,foodlevels=foodlevels,T_F_min=T_F_min,T_F_max=T_F_max,CT_max=CT_max,CT_min=CT_min,T_B_min=T_B_min,T_RB_min=T_RB_min))
    }else{
      return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,soilpot=soilpot,shadpot=shadpot,humid=humid,shadhumid=shadhumid,tcond=tcond,shadtcond=shadtcond,rainfall=rainfall,rainhr=rainhr,enbal=enbal,masbal=masbal,environ=environ,debout=debout,yearout=yearout,yearsout=yearsout,foodwaters=foodwaters,foodlevels=foodlevels,T_F_min=T_F_min,T_F_max=T_F_max,CT_max=CT_max,CT_min=CT_min,T_B_min=T_B_min,T_RB_min=T_RB_min))
    }
  } # end error check
}
