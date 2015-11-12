#' Ectotherm model.
#'
#' An implementation of the Niche Mapper ectotherm model that computes body temperature, water loss,
#' activity and microhabitat selection. It optionally runs the Dynamic Energy Budget (DEB) model for
#' computing mass budgets (inc. water budgets) and growth, development, reproduction trajectories
#' as constrained by food, activity and temperature (see Details). When not running the DEB model
#' a user-specified mass is used as well as a allometric (mass and body temperature) function to
#' compute metabolic rate. \cr\cr NOTE: The microclimate model, e.g. micro_global(), must be run prior to running the ectotherm model
#'
#' @param amass Mass of animal (g), note this model is 'steady state' so no lags in heating/cooling due to mass
#' @param lometry Organism shape, 0-5, Determines whether standard or custom shapes/surface area/volume relationships are used: 0=plate, 1=cyl, 2=ellips, 3=lizard (desert iguana), 4=frog (leopard frog), 5=custom (see details)
#' @param ABSMAX Maximum solar absorptivity, decimal percent
#' @param ABSMIN Maximum solar absorptivity, decimal percent
#' @param TMAXPR Voluntary thermal maximum, degrees C (upper body temperature for foraging and also affects burrow depth selection)
#' @param TMINPR Voluntary thermal minimum, degrees C (lower body temperature for foraging)
#' @param TBASK Minimum basking temperature, degrees C
#' @param TEMERGE Temperature at which animal will move to a basking site, degrees C
#' @param TPREF Preferred body temperature, degrees C
#' @param ctmax Critical thermal maximum, degrees C (affects burrow depth selection)
#' @param ctmin Critical thermal minimum, degrees C (affects burrow depth selection)
#' @param dayact Diurnal activity allowed?  1=yes, 0=no
#' @param nocturn Nocturnal activity allowed?  1=yes, 0=no
#' @param crepus Crepuscular activity allowed?  1=yes, 0=no
#' @param CkGrShad shade seeking allowed?  1=yes, 0=no
#' @param burrow Shelter in burrow allowed?  1=yes, 0=no
#' @param climb climbing to seek cooler habitats allowed?  1=yes, 0=no
#' @param shdburrow choose if the animal's retreat is in the shade (1) or in the open (0)
#' @param mindepth Minimum depth (soil node #) to which animal can retreat if burrowing
#' @param maxdepth Maximum depth (soil node #) to which animal can retreat if burrowing
#' @param MR_1 Metabolic rate parameter MR=MR_1*M^MR_2*10^(MR_3*Tb) based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231
#' @param MR_2 Metabolic rate parameter
#' @param MR_3 Metabolic rate parameter
#' @param skinwet \% of surface area acting as a free-water exchanger, for computing cutaneous water loss
#' @param extref \% oxygen extraction efficiency, for respiratory water loss
#' @param DELTAR Temperature difference (deg C) between expired and inspired air, for computing respiratory water loss
#' @usage ectotherm(amass, lometry, ABSMAX, ABSMIN, TMAXPR, TMINPR, TBASK, TEMERGE, ctmax, ctmin,
#'  TPREF, dayact, nocturn, crepus, CkGrShad, burrow, climb, shdburrow, mindepth, maxdepth,
#'  MR_1, MR_2, MR_3, skinwet, extref, DELTAR, ...)
#' @details
#' \strong{ Parameters controling how the model runs:}
#' \itemize{
#' \item{\code{microin}{ = "none", Directory where precomputed microclimate outputs (.csv files) are saved - if "none", assumes results are in list 'micro'}\cr}
#' \item{\code{nyears}{ = micro$nyears, Number of years the simulation runs for - must be consistent with dimensions of environmental input data}\cr}
#' \item{\code{ystrt}{ = 0, Year of the simulation to start at - will loop back over earlier periods at the end, useful for exploring cohort effects when running a life cycle with the DEB model}\cr}
#' \item{\code{enberr}{ = 0.0002, Factor by which the mass is multiplied to obtain a tolerance level for the heat budget solution}\cr}
#' \item{\code{live=}{ = 1, Live (metabolism/behaviour) or dead animal?}\cr}
#' \item{\code{write_input}{ = 0, Write csv files of final input to folder 'csv input' in working directory? 1=yes, 0=no}\cr}
#' \item{\code{startday}{ = 1, Julian day of year at which simulation starts}\cr}
#'}
#' \strong{ Environmental inputs:}
#'
#' \itemize{
#' \item{\code{minshade}{ = 0., Minimum shade (\%) available to the animal}\cr}
#' \item{\code{maxshade}{ = micro$MAXSHADES[1], Maximum shade (\%) available to the animal}\cr}
#' \item{\code{FLTYPE}{ = 0.0, Fluid type 0=air, 1=water }\cr}
#' \item{\code{SUBTK}{ = 2.79, Substrate thermal conductivity (W/mC)}\cr}
#' \item{\code{REFL}{ = micro$REFL, Vector of daily substrate reflectances (\%)}\cr}
#' \item{\code{DEP}{ = micro$DEP,}\cr}
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
#' \item{\code{RAINFALL}{ = micro$RAINFALL, Vector of daily rainfall}\cr}
#' \item{\code{MAXSHADES}{ = micro$MAXSHADES, Vector of daily maximum shade values}\cr}
#' \item{\code{ectoin}{ = rbind(as.numeric(micro$ALTT),as.numeric(micro$REFL)[1],micro$longlat[1],micro$longlat[2],0,0,1990,1990), Other items needed by the model - this needs to be tidied up}\cr}
#'}
#' \strong{ Morphological parameters:}
#'
#' \itemize{
#' \item{\code{customallom}{ = c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743), Custom allometry coefficients. Operates if lometry=5, and consists of 4 pairs of values representing the parameters a and b of a relationship AREA=a*mass^b, where AREA is in cm2 and mass is in g. The first pair are a and b for total surface area, then a and b for ventral area, then for sillhouette area normal to the sun, then sillhouette area perpendicular to the sun}\cr}
#' \item{\code{shape_a}{ = 1., Proportionality factor (-) for going from volume to area, keep this 1 (redundant parameter that should be removed)}\cr}
#' \item{\code{shape_b}{ = 3, Proportionality factor (-) for going from volume to area, represents ratio of width:height for a plate, length:diameter for cylinder, b axis:a axis for ellipsoid }\cr}
#' \item{\code{shape_c}{ = 0.6666666667, Proportionality factor (-) for going from volume to area, represents ratio of length:height for a plate, c axis:a axis for ellipsoid}\cr}
#' \item{\code{FATOSK}{ = 0.4, Configuration factor to sky (-) for infrared calculations}\cr}
#' \item{\code{FATOSB}{ = 0.4, Configuration factor to subsrate for infrared calculations}\cr}
#' \item{\code{rinsul}{ = 0., Insulative fat layer thickness (m)}\cr}
#' \item{\code{ptcond}{ = 0.25, Fraction of surface contacting the substrate}\cr}
#' \item{\code{Spheat}{ = 4185, Specific heat of flesh J/(kg-K)}\cr}
#' \item{\code{Flshcond}{ = 0.5, Thermal conductivity of flesh (W/mC, range: 0.412-2.8)}\cr}
#' \item{\code{Andens}{ = 1000, Density of flesh (kg/m3)}\cr}
#' \item{\code{EMISAN}{ = 0.95, Emissivity of animal (0-1)}\cr}
#'}
#' \strong{ Behavioural parameters:}
#'
#' \itemize{
#' \item{\code{fosorial}{ = 0, Fossorial activity? 1=yes, 0=no (this option hasn't been tested)}\cr}
#' \item{\code{rainact}{ = 0, Activity is limited by rainfall? 1=yes, 0=no, threshold rainfall for activity set by \code{actrainthresh}}\cr}
#' \item{\code{actrainthresh}{ = 0.1, Threshold (mm) of rain causing activity if \code{rainact}=1}\cr}
#' \item{\code{soilnode}{ = 4., Soil node (1-10, corresponding to values in \code{DEP}) at which eggs are laid (overridden if \code{frogbreed}=1)}\cr}
#'}
#' \strong{ Thermal physiological parameters:}
#'
#' \itemize{
#' \item{\code{ctminthresh}{ = 12, Number of consecutive hours below CTmin that leads to death - simulation will terminate beyond this threshold if \code{ctkill}=1}\cr}
#' \item{\code{ctkill}{ = 0, Animal dies when it hits critical thermal limits? 1=yes, 0=no}\cr}
#'}
#' \strong{ Water and food budget parameters (only relevant if \code{DEB}=1):}
#'
#' \itemize{
#' \item{\code{PFEWAT}{ = 73, Fecal water content (\%)}\cr}
#' \item{\code{PTUREA}{ = 0, Water in excreted nitrogenous waste (\%)}\cr}
#' \item{\code{FoodWater}{ = 82, Wwater content of food (\%)}\cr}
#' \item{\code{minwater}{ = 15, Minimum tolerated dehydration (\% of wet mass) - prohibits foraging if greater than this}\cr}
#' \item{\code{gutfill}{ = 75.,}\cr}
#' \item{\code{soilmoisture1}{ = 0,}\cr}
#' \item{\code{raindrink}{ = 0., Gut fill (\%) at which satiation occurs - if greater than 100\%, animal always tries to forage}\cr}
#'}
#' \strong{ Dynamic Energy Budget (DEB) model parameters:}
#' \itemize{
#' \item{\code{DEB}{ = 0, Run the DEB model (1) or just heat balance (0). Latter uses allometrically predicted respiration base on \code{MR_1}, \code{MR_2} and \code{MR_3}}\cr}
#' \item{\code{fract}{ = 1, Scaling factor for DEB body-size covariation relationships - use it to make a metabolically scaled larger or smaller version of your animal}\cr}
#'}
#' \strong{ Core DEB parameters:}
#' \itemize{
#' \item{\code{z}{ = 7.174*fract, Zoom factor (cm)}\cr}
#' \item{\code{delta}{ =  0.217, Shape coefficient (-)}\cr}
#' \item{\code{p_Xm}{ = 13290, Surface area-specific maximum feeding rate J/cm2/h}\cr}
#' \item{\code{kappa_X}{ = 0.85, Digestive efficiency (decimal \%)}\cr}
#' \item{\code{v_dotref}{ = 0.05591/24., Energy conductance (cm/h)}\cr}
#' \item{\code{kappa}{ = 0.8501, fraction of mobilised reserve allocated to soma}\cr}
#' \item{\code{p_Mref}{ = 45.14/24., Volume-specific somatic maintenance (J/cm3/h)}\cr}
#' \item{\code{E_G}{ = 7189, Cost of structure (J/cm3)}\cr}
#' \item{\code{k_R}{ = 0.95, Fraction of reproduction energy fixed in eggs}\cr}
#' \item{\code{k_J}{ = 0.00628/24., Maturity maintenance rate coefficient (1/h)}\cr}
#' \item{\code{E_Hb}{ = 6.533e+04*fract^3, Maturity at birth (J)}\cr}
#' \item{\code{E_Hj}{ = E_Hb*fract^3, Maturity at metamorphosis (J)}\cr}
#' \item{\code{E_Hp}{ = 1.375e+05*fract^3, Maturity at puberty}\cr}
#' \item{\code{h_aref}{ = 3.61e-13/(24.^2), Weibull ageing acceleration (1/h2)}\cr}
#' \item{\code{s_G}{ = 0.01, Gompertz stress coefficient (-)}\cr}
#'}
#' \strong{ Thermal DEB parameters:}
#' \itemize{
#' \item{\code{T_REF}{ = 20, Reference temperature for rate correction (deg C)}\cr}
#' \item{\code{TA}{ = 7130, Arhhenius temperature}\cr}
#' \item{\code{TAL}{ = 5.305e+04, Arrhenius temperature for decrease below lower boundary of tolerance range \code{TL}}\cr}
#' \item{\code{TAH}{ = 9.076e+04, Arrhenius temperature for decrease above upper boundary of tolerance range \code{TH}}\cr}
#' \item{\code{TL}{ = 288., Lower boundary (K) of temperature tolerance range for Arrhenius thermal response}\cr}
#' \item{\code{TH}{ = 315., Upper boundary (K) of temperature tolerance range for Arrhenius thermal response}\cr}
#'}
#' \strong{ Compound/derived DEB parameters:}
#' \itemize{
#' \item{\code{E_Egg}{ = 1.04e+06*fract^4, Energy content of the egg (J)}\cr}
#' \item{\code{E_m}{ = (p_Mref*z/kappa)/v_dotref, Maximum reserve density (J/cm3)}\cr}
#'}
#' \strong{ Food-related axilliary DEB parameters:}
#' \itemize{
#' \item{\code{f}{ = 1, functional response (-), usually kept at 1 because gut model controls food availability such that f=0 when gut empty}\cr}
#' \item{\code{MsM}{ = 186.03*6, Maximum volume-specific energy density of stomach (J/cm3)}\cr}
#' \item{\code{K}{ = 1, Half saturation constant (#/cm2)}\cr}
#' \item{\code{X}{ = 10, Food density (J/cm2)}\cr}
#'}
#' \strong{ Composition-related axilliary DEB parameters:}
#' \itemize{
#' \item{\code{andens_deb}{ = Andens/1000, Animal density (g/cm3)}\cr}
#' \item{\code{d_V}{ = 0.3, Dry mass fraction of structure}\cr}
#' \item{\code{d_E}{ = 0.3, Dry mass fraction of reserve}\cr}
#' \item{\code{eggdryfrac}{ = 0.3, Dry mass fraction of egg}\cr}
#' \item{\code{mu_X}{ = 525000, Molar Gibbs energy (chemical potential) of food (J/mol)}\cr}
#' \item{\code{mu_E}{ = 585000, Molar Gibbs energy (chemical potential) of reserve (J/mol)}\cr}
#' \item{\code{mu_V}{ = 500000, Molar Gibbs energy (chemical potential) of structure (J/mol)}\cr}
#' \item{\code{mu_P}{ = 480000, Molar Gibbs energy (chemical potential) of faeces (J/mol)}\cr}
#' \item{\code{kappa_X_P}{ = 0.1, Faecation efficiency of food to faeces (-)}\cr}
#' \item{\code{nX}{ = c(1,1.8,0.5,.15), chem. indices of C, O, H and N in food}\cr}
#' \item{\code{nE}{ = c(1,1.8,0.5,.15), chem. indices of C, O, H and N in reserve}\cr}
#' \item{\code{nV}{ = c(1,1.8,0.5,.15), chem. indices of C, O, H and N in structure}\cr}
#' \item{\code{nP}{ = c(1,1.8,0.5,.15), chem. indices of C, O, H and N in faeces}\cr}
#' \item{\code{N_waste}{ = c(1,4/5,3/5,4/5), chem. indices of C, O, H and N in nitrogenous waste}\cr}
#'}
#' \strong{ Insect DEB model parameters (not yet in operation):}
#' \itemize{
#' \item{\code{metab_mode}{ = 0, Run insect model? 0 = no, 1 = hemimetabolus model (to do), 2 = holometabolous model}\cr}
#' \item{\code{stages}{ = 7, number of stages (max = 8) = number of instars plus 1 for egg + 1 for pupa + 1 for imago}\cr}
#' \item{\code{y_EV_l}{ = 0.95, yield of imago reserve on larval structure (mol/mol)}\cr}
#' \item{\code{S_instar}{ = c(2.660,2.310,1.916,0), stress at instar n: L_n^2/ L_n-1^2 (-)}\cr}
#' \item{\code{s_j}{ = 0.999, Reprod buffer/structure ratio at pupation as fraction of max}\cr}
#'}
#' \strong{ Inital conditions for DEB model:}
#' \itemize{
#' \item{\code{v_init}{ = 3e-9, Initial structural volume (cm3)}\cr}
#' \item{\code{E_init}{ = E_Egg/v_init, Initial reserve density (J/cm3)}\cr}
#' \item{\code{E_H_init}{ = 0, Initial maturity (J)}\cr}
#' \item{\code{stage}{ = 0, Initial stage (0=embryo, 1=juvenile, 2=mature but not yet reproducing, 3=beyond first reproduction)}\cr}
#'}
#' \strong{ Metabolic depression parameters (not yet functional):}
#' \itemize{
#' \item{\code{aestivate}{ = 0, Does the animal aestivate/go into torpor? 1=yes, 0=no}\cr}
#' \item{\code{depress}{ = 0.3, Fraction by which \code{p_Mref}, \code{k_J} and \code{v_dotref} are reduced during torpor}\cr}
#'}
#' \strong{ Reproductive phenology model parameters:}
#' \itemize{
#' \item{\code{clutchsize}{ = 2., Clutch size (#), overridden by \code{clutch_ab}}\cr}
#' \item{\code{clutch_ab}{ = c(0,0), # paramters for relationship between length (cm) and clutch size: clutch size = a*SVL-b, make a and b zero if fixed clutch size}\cr}
#' \item{\code{viviparous}{ = 0, Viviparous reproduction? 1=yes, 0=no (if yes, animal will be held in adult-sided female's body for duration of development and will experience her body temperature}\cr}
#' \item{\code{minclutch}{ = 0, Minimum clutch size if not enough in reproduction buffer for clutch size predicted by \code{clutch_ab} - if zero, will not operate}\cr}
#' \item{\code{batch}{ = 1, Invoke Pequerie et al.'s batch laying model?}\cr}
#' \item{\code{photostart}{ = 3, Photoperiod response triggering ovulation, none (0), summer solstice (1), autumnal equinox (2), winter solstice (3), vernal equinox (4), specified daylength thresholds (5 - uses \code{daylengthstart} and \code{daylengthfinish})}\cr}
#' \item{\code{photofinish}{ = 1, Photoperiod terminating ovulation, none (0), summer solstice (1), autumnal equinox (2), winter solstice (3), vernal equinox (4), specified daylength thresholds (5 - uses \code{daylengthstart} and \code{daylengthfinish})}\cr}
#' \item{\code{daylengthstart}{ = 12.5, Threshold daylength (h) for initiating breeding}\cr}
#' \item{\code{daylengthfinish}{ = 13., Threshold daylength (h) for terminating breeding}\cr}
#' \item{\code{photodirs }{ =  1, Is the start daylength trigger during a decrease (0) or increase (1) in day length?}\cr}
#' \item{\code{photodirf }{ =  0, Is the finish daylength trigger during a decrease (0) or increase (1) in day length?}\cr}
#' \item{\code{frogbreed}{ = 0, Amphibious animal breeding mode: 0 is off, 1 is exotrophic aquatic (eggs start when water present in container and within breeding season), 2 is exotrophic terrestrial/aquatic (eggs start at specified soil node within breeding season, diapause at birth threshold, start larval phase if water present in container), 3 endotrophic terrestrial (eggs start at specified soil node within breeding season and continue to metamorphosis on land), 4 turtle mode (eggs start at specified soil node within breeding season, hatch and animals enter water and stay there for the rest of their life, but leave the water if no water is present)}\cr}
#' \item{\code{frogstage}{ = 0, Life cycle control for amphibious animal: 0 is whole life cycle, 1 is just to metamorphosis (then reset and start again)}\cr}
#' \item{\code{reset}{ = 0, Life cycle reset options, 0=quit simulation upon death, 1=restart at emergence, 2=restart at first egg laid, 3=restart at end of breeding season, 4=reset at death}\cr}
#' \item{\code{breedactthresh}{ = 1, Threshold numbers of hours active after start of breeding season before eggs can be laid (simulating movement to the breeding site)}\cr}
#' \item{\code{breedrainthresh}{ = 0, Rain dependent breeder? 0 means no, otherwise enter rainfall threshold in mm}\cr}
#' \item{\code{breedtempthresh}{ = 200, Body temperature threshold below which breeding will occur}\cr}
#' \item{\code{breedtempcum}{ = 24*7, Cumulative time below temperature threshold for breeding \code{breedtempthresh} that will trigger breeding}\cr}
#'}
#' \strong{ Mortality rate parameters:}
#' \itemize{
#' \item{\code{ma}{ = 1e-4, Hourly active mortality rate (probability of mortality per hour)}\cr}
#' \item{\code{mi}{ = 0, Hourly inactive mortality rate (probability of mortality per hour)}\cr}
#' \item{\code{mh}{ = 0.5, Survivorship of hatchling in first year}\cr}
#'}
#' \strong{ Water body model parameters:}
#' \itemize{
#' \item{\code{container}{ = 0, Run the water body/container model? (aquatic start of life cycle, e.g. frog or mosquito)}\cr}
#' \item{\code{wetmod}{ = 0, Use precomputed wetland temperature \code{wetlandTemps} and depths \code{wetlandDepthds}?}\cr}
#' \item{\code{conth}{ = 10, Cylindrical container/pond height (cm)}\cr}
#' \item{\code{contw}{ = 100., Cylindrical container/pond diameter (cm)}\cr}
#' \item{\code{contype}{ = 1, Is 'containter' sitting on the surface, like a bucket (0) or sunk into the ground like a pond (1)}\cr}
#' \item{\code{rainmult}{ = 1, Rainfall multiplier to reflect catchment (don't make this zero unless you want a drought!)}\cr}
#' \item{\code{continit}{ = 0, Initial container water level (cm)}\cr}
#' \item{\code{conthole}{ = 0, Daily loss of height (mm) due to 'hole' in container (e.g. infiltration to soil, drawdown from water tank)}\cr}
#' \item{\code{contonly}{ = 1, Just run the container model and quit?}\cr}
#' \item{\code{contwet}{ = 80, \% of container surface acting as a free water exchanger}\cr}
#' \item{\code{wetlandTemps}{ = matrix(data = 0., nrow = 24*dim, ncol = 1), Matrix of hourly wetland temperaures (deg C)}\cr}
#' \item{\code{wetlandDepths}{ = matrix(data = 0., nrow = 24*dim, ncol = 1), Matrix of hourly wetland depths (cm)}\cr}
#'}
#' \strong{ Life stage-specific parameter allocation:}
#' \itemize{
#' \item{\code{thermal_stages}{ = matrix(data = c(rep(ctmin,8),rep(ctmax,8),rep(TMINPR,8),rep(TMAXPR,8),rep(TBASK,8),rep(TPREF,8)), nrow = 8, ncol = 6), Stage specific thermal thresholds (ctmin,ctmax,TMINPR,TMAXPR,TBASK,TPREF)}\cr}
#' \item{\code{behav_stages}{ = matrix(data = c(rep(dayact,8),rep(nocturn,8),rep(crepus,8),rep(burrow,8),rep(shdburrow,8),rep(mindepth,8),rep(maxdepth,8),rep(CkGrShad,8),rep(climb,8),rep(fosorial,8),rep(rainact,8),rep(actrainthresh,8),rep(breedactthresh,8),rep(flyer,8)), nrow = 8, ncol = 14), Stage specific behaviour dayact,nocturn,crepus,burrow,shdburrow,mindepth,maxdepth,CkGrShad,climb,fosorial,rainact,actrainthresh,breedactthresh,flyer)}\cr}
#' \item{\code{water_stages}{ = matrix(data = c(rep(skinwet,8),rep(extref,8),rep(PFEWAT,8),rep(PTUREA,8),rep(FoodWater,8),rep(minwater,8),rep(raindrink,8),rep(gutfill,8)), nrow = 8, ncol = 8), Stage-specific water budget parameters (skinwet,extref,PFEWAT,PTUREA,FoodWater,minwater,raindrink,gutfill)}\cr}
#' \item{\code{arrhenius}{ = matrix(data = matrix(data = c(rep(TA,8),rep(TAL,8),rep(TAH,8),rep(TL,8),rep(TH,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5), Stage-specific 5-parameter Arrhenius thermal response for DEB model (TA,TAL,TAH,TL,TH)}\cr}
#'}
#' \strong{ Butterfly model parameters (not yet tested):}
#' \itemize{
#' \item{\code{wings}{ = 0, Turn wing model on? 1=yes, 0=no}\cr}
#' \item{\code{rho1_3}{ = 0.2, Wing reflectance (dec \%)}\cr}
#' \item{\code{trans1}{ = 0.00, Wing transmissivity (dec \%)}\cr}
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
#' environ variables: !add solar
#' \itemize{
#' \item 1 JULDAY - day of year
#' \item 2 YEAR - year of simulation
#' \item 3 DAY - day of simulation
#' \item 4 TIME - time of day (hours)
#' \item 5 TC - body temperature (deg C)
#' \item 6 SHADE - shade selected (\%)
#' \item 7 ORIENT  - orientation to sun (1=normal, 0=parallel) (currently fixed at 1 in code)
#' \item 8 DEP - depth below ground (cm)
#' \item 9 ACT - activity state (0=inactive, 1=basking, 2=foraging)
#' \item 10 TA - air temperature (deg C) at animal location
#' \item 11 VEL - wind speed (m/s) at animal location
#' \item 12 RELHUM - relative humidity (\%) at animal location
#' \item 13 ZEN - zenith angle of sun (degrees - 90 = below the horizon)
#' \item 14 CONDEP - depth of water body (cm) (may not be simulated or supplied)
#' \item 15 WATERTEMP - temperature of water body (deg C) (may not be simulated or supplied)
#' \item 16 DAYLENGTH - day length (hours)
#' \item 17 WINGANGLE - wing angle (degrees) for butterfly model
#' \item 18 WINGTEMP - wing temperature (degrees C) for butterfly model
#' \item 19 FLYING - flying state (1=flying, 0=not flying) for butterfly model
#' \item 20 FLYTIME - flying time (hours) for butterfly model
#'}
#' enbal variables:
#' \itemize{
#' \item 1 JULDAY - day of year
#' \item 2 YEAR - year of simulation
#' \item 3 DAY - day of simulation
#' \item 4 TIME - time of day (hours)
#' \item 5 TC - body temperature (deg C) ! redundant
#' \item 6 QSOL - solar radiation absorbed (W)
#' \item 7 QIRIN  - infrared radiation absorbed (W)
#' \item 8 QMET - metabolic heat production (W)
#' \item 9 QEVAP - evaporative heat loss (W)
#' \item 10 QIROUT - infrared radiation lost (W)
#' \item 11 QCONV - heat lost by convection (W)
#' \item 12 QCOND - heat lost by conduction (W)
#' \item 13 ENB - energy balance (degrees C)
#' \item 14 NTRY - iterations required for solution to heat balance equation
#'}
#' masbal variables:
#' \itemize{
#' \item 1 JULDAY - day of year
#' \item 2 YEAR - year of simulation
#' \item 3 DAY - day of simulation
#' \item 4 TIME - time of day (hours)
#' \item 5 TC - body temperature (deg C) ! redundant
#' \item 6 O2_ml - oxygen consumption rate (ml/h)
#' \item 7 CO2_ml - carbon dioxide production rate (ml/h)
#' \item 8 NWASTE_g - nitrogenous waste production (g/h)
#' \item 9 H2OFree_g - water from food (g/h)
#' \item 10 H2OMet_g - metabolic water production (g/h)
#' \item 11 DryFood_g - dry food intake (g/h)
#' \item 12 WetFood_g - wet foood intake (g/h)
#' \item 13 DryFaeces_g - dry faeces production (g/h)
#' \item 14 WetFaeces_G - wet faeces production (g/h)
#' \item 15 Urine_g - urine production (g/h)
#' \item 16 H2OResp_g - respiratory water loss (g/h)
#' \item 17 H2OCut_g - cutaneous water loss (g/h)
#' \item 18 H2OEvap_g - evaporative water loss (cutaneous plus respiratory) (g/h) !redundant
#' \item 19 H2OBal_g - instantaneous water balance (g/h)
#' \item 20 H2OCumBal_g - cumulative water balance (g)
#' \item 21 GutFreeMass_g - gut-free wet mass (g)
#'}
#' debout variables:
#' \itemize{
#' \item 1 JULDAY - day of year
#' \item 2 YEAR - year of simulation
#' \item 3 DAY - day of simulation
#' \item 4 TIME - time of day (hours)
#' \item 5 WETMASS - wet mass (g)
#' \item 6 RESERVE_DENS  - reserve density (J/cm3)
#' \item 7 CUMREPRO - energy in reproduction buffer (J)
#' \item 8 HS - hazard rate (1/h)
#' \item 9 MASS_GUT - wet gut contents (g)
#' \item 10 SVL - length (mm) (might be Snout Vent Length but depends on choice of length measure for DEB paramter fitting)
#' \item 11 V - structural volume (cm3)
#' \item 12 E_H - maturity state (J)
#' \item 13 CUMBATCH - energy in batch for egg production (J)
#' \item 14 V_baby - structure of baby (cm3) (only if viviparous and pregnant
#' \item 15 E_baby - reserve density of baby (J/cm3) (only if viviparous and pregnant)
#' \item 16 Pregnant - pregnant? (only if viviparous) (0 or 1)
#' \item 17 Stage - life cycle stage (0=embryo, 1=juvenile, 2=mature but not yet reproducing, 3=beyond first reproduction)
#' \item 18 WETMASS_STD - wet mass (g) excluding stomach contents
#' \item 19 Body_cond - \% desiccated
#' \item 20 Surviv_Prob - survival probability due to joint influence of ageing and mortality rates
#' \item 21 Breeding - breeding state (1=breeding, 0=not breeding)
#'}
#' yearout variables:
#' \itemize{
#' \item 1 DEVTIME - development time (days)
#' \item 2 BIRTHDAY - birth day (day of year)
#' \item 3 BIRTHMASS - mass at birth (g)
#' \item 4 MONMATURE - months to maturity
#' \item 5 SVLREPRO - length (mm) at first reproduction
#' \item 6 FECUNDITY  - total fecundity
#' \item 7 CLUTCHES - total clutches
#' \item 8 MINRESERVE - minimum reserve density (J/cm3)
#' \item 9 LASTFOOD - food in last year (kg)
#' \item 10 TOTFOOD - total food eaten (kg)
#' \item 11 MINTB - minimum body temperature (deg C)
#' \item 12 MAXTB - maxmium body temperature (deg C)
#' \item 13 Pct_Dess - maximum level of desiccation (/%)
#' \item 14 LifeSpan - maximum life span (days)
#' \item 15 GenTime - generation time (years)
#' \item 16 R0 - net reproductive rate
#' \item 17 rmax - intrinsic rate of increase
#' \item 18 SVL - maximum length (mm)
#'}
#' yearsout variables:
#' \itemize{
#' \item 1 YEAR - year of simulation
#' \item 2 MaxStg - maximum stage reached in the year
#' \item 3 MaxWgt - maximum weight reached in the year (g)
#' \item 4 MaxLen - maximum length in the year (mm)
#' \item 5 Tmax - maximum annual body temperature (deg C)
#' \item 6 Tmin  - minimum annual body temperature (deg C)
#' \item 7 MinRes - minimum annual reserve density (J/cm3)
#' \item 8 MaxDess - maximum annual desiccation level (% of normal wet body mass)
#' \item 9 MinShade - minimum annual shade selected
#' \item 10 MaxShade - maximum annual shade selected
#' \item 11 MinDep - minimum annual depth selected (cm)
#' \item 12 MaxDep - maximum annual depth selected (cm)
#' \item 13 Bsk - annual basking hours
#' \item 14 Forage - annual foraging hours
#' \item 15 Dist - annual distance travelled (km)
#' \item 16 Food - annual food eaten (g)
#' \item 17 Drink - annual water drunk (g)
#' \item 18 NWaste - annual nitrogenous waste (g)
#' \item 19 O2 - annual O2 production (ml)
#' \item 20 Clutch - annual clutches (#)
#' \item 21 Fec - annual fecundity (#)
#' \item 22 CauseDeath - cause of death, 0=no death, 1=cold, 2=heat, 3=desiccation, 4=starvation, 5=ageing
#' \item 23 tLay - day of year at which eggs laid
#' \item 24 tEgg - day of year entering egg stage
#' \item 25 tStg1 - day of year entering life cycle stage 1
#' \item 26 tStg2 - day of year entering life cycle stage 2
#' \item 27 tStg3 - day of year entering life cycle stage 3
#' \item 28 tStg4 - day of year entering life cycle stage 4
#' \item 29 tStg5 - day of year entering life cycle stage 5
#' \item 30 tStg6 - day of year entering life cycle stage 6
#' \item 31 tStg7 - day of year entering life cycle stage 7
#' \item 32 tStg8 - day of year entering life cycle stage 8
#' \item 33 mStg1 - mass entering life cycle stage 1 (g)
#' \item 34 mStg2 - mass entering life cycle stage 2 (g)
#' \item 35 mStg3 - mass entering life cycle stage 3 (g)
#' \item 36 mStg4 - mass entering life cycle stage 4 (g)
#' \item 37 mStg5 - mass entering life cycle stage 5 (g)
#' \item 38 mStg6 - mass entering life cycle stage 6 (g)
#' \item 39 mStg7 - mass entering life cycle stage 7 (g)
#' \item 40 mStg8 - mass entering life cycle stage 8 (g)
#' \item 41 surviv - survival probability at end of given year
#' \item 42 ovip_surviv - breeding state (1=breeding, 0=not breeding)
#' \item 43 fitness - fecundity by survival probability at end of given year
#' \item 44 deathstage - life stage at which death occurred
#'}
#' @examples
#'# run the microclimate model
#'micro<-micro_global(loc="Kuranda, Queensland")
#'
#'# retrieve output
#'metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#'shadmet<-as.data.frame(micro$shadmet) # above ground microclimatic conditions, max shade
#'soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
#'shadsoil<-as.data.frame(micro$shadsoil) # soil temperatures, maximum shade
#'
#'# append dates
#'days<-rep(seq(1,12),24)
#'days<-days[order(days)]
#'dates<-days+metout$TIME/60/24-1 # dates for hourly output
#'dates2<-seq(1,12,1) # dates for daily output
#'metout<-cbind(dates,metout)
#'soil<-cbind(dates,soil)
#'shadmet<-cbind(dates,shadmet)
#'shadsoil<-cbind(dates,shadsoil)
#'
#'# run the ectotherm model
#'ecto<-ectotherm(TMAXPR=35,TMINPR=30,TPREF=33,TBASK=20,TEMERGE=10)
#'
#'# retrieve output
#'environ<-as.data.frame(ecto$environ) # activity, Tb and environment
#'enbal<-as.data.frame(ecto$enbal) # energy balance values
#'masbal<-as.data.frame(ecto$masbal) # mass balance value (note most missing if DEB model not running)
#'
#'# append dates
#'environ<-cbind(dates,environ)
#'masbal<-cbind(dates,masbal)
#'enbal<-cbind(dates,enbal)
#'
#'############### plot results ######################
#'
#'# Hourly Tb (black), activity (orange, 5=bask, 10=forage), depth (brown, m) and shade (green, %/10)
#'with(environ, plot(TC~dates,ylab="Tb, depth, activity and shade", xlab="month of year",
#' ylim=c(-20,70),type = "l"))
#'with(environ, points(ACT*5~dates,type = "l",col="orange"))
#'with(environ, points(SHADE/10~dates,type = "l",col="green"))
#'with(environ, points(DEP/10~dates,type = "l",col="brown"))
#'#with(metout, points(TAREF~dates,type = "l",col="light blue"))
#'abline(ecto$TMAXPR,0,lty=2,col='red')
#'abline(ecto$TMINPR,0,lty=2,col='blue')
#'
#'# seasonal activity plot (dark blue = night, light blue = basking, orange = foraging)
#'forage<-subset(environ,ACT==2)
#'bask<-subset(environ,ACT==1)
#'night<-subset(metout,ZEN==90)
#'day<-subset(metout,ZEN!=90)
#'with(night,plot(TIME/60~JULDAY,ylab="Hour of Day",xlab="Day of Year",pch=15,cex=2,col='dark blue'))
#' # nighttime hours
#'with(forage,points((TIME-1)~JULDAY,pch=15,cex=2,col='orange')) # foraging Tbs
#'with(bask,points((TIME-1)~JULDAY,pch=15,cex=2,col='light blue')) # basking Tbs
#' @export
ectotherm<-function(amass=5,lometry=3,ABSMAX=0.85,ABSMIN=0.85,TMAXPR=35,TMINPR=25,TBASK=20,TEMERGE=10,
TPREF=30,ctmax=40,ctmin=5,dayact=1,nocturn=0,crepus=0,CkGrShad=1,burrow=1,climb=0,shdburrow=0,
mindepth=2,maxdepth=10,MR_1=0.013,MR_2=0.8,MR_3=0.038,skinwet=0.2,extref=20.,DELTAR=0.1,
microin="none",nyears=micro$nyears,ystrt=0,enberr=0.0002,live=1,write_input=0,
startday=1,minshade=0.,maxshade=micro$MAXSHADES[1],FLTYPE=0.0,SUBTK=2.79,REFL=micro$REFL,
DEP=micro$DEP,metout=micro$metout,shadmet=micro$shadmet,soil=micro$soil,shadsoil=micro$shadsoil,
soilmoist=micro$soilmoist,shadmoist=micro$shadmoist,humid=micro$humid,shadhumid=micro$shadhumid,
soilpot=micro$soilpot,shadpot=micro$shadpot,RAINFALL=micro$RAINFALL,MAXSHADES=micro$MAXSHADES,
ectoin=rbind(as.numeric(micro$ALTT),as.numeric(micro$REFL)[1],micro$longlat[1],micro$longlat[2]
  ,0,0,1990,1990),customallom=c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743),
shape_a=1.,shape_b=3,shape_c=0.6666666667,FATOSK=0.4,FATOSB=0.4,rinsul=0.,ptcond=0.25,
Spheat=4185,Flshcond=0.5,Andens=1000,EMISAN=0.95,
fosorial=0,rainact=0,actrainthresh=0.1,soilnode=4.,ctminthresh=12,ctkill=0,
PFEWAT=73,PTUREA=0,FoodWater=82,minwater=15,gutfill=75.,soilmoisture1=0,raindrink=0.,
DEB=0,fract=1,z=7.174*fract,delta= 0.217,p_Xm=13290,kappa_X=0.85,v_dotref=0.05591/24.,
kappa=0.8501,p_Mref=45.14/24.,E_G=7189,k_R=0.95,k_J=0.00628/24.,E_Hb=6.533e+04*fract^3,
E_Hj=E_Hb*fract^3,E_Hp=1.375e+05*fract^3,h_aref=3.61e-13/(24.^2),s_G=0.01,
T_REF=20,TA=7130,TAL=5.305e+04,TAH=9.076e+04,TL=288.,TH=315.,
E_Egg=1.04e+06*fract^4,E_m=(p_Mref*z/kappa)/v_dotref,f=1.,MsM=186.03*6.,K=1,X=10,andens_deb=
    Andens/1000,d_V=0.3,d_E=0.3,eggdryfrac=0.3,mu_X=525000,mu_E=585000,mu_V=500000,mu_P=480000,
kappa_X_P=0.1,nX=c(1,1.8,0.5,.15),nE=c(1,1.8,0.5,.15),nV=c(1,1.8,0.5,.15),nP=c(1,1.8,0.5,.15),
N_waste=c(1,4/5,3/5,4/5),metab_mode=0,stages=7,y_EV_l=0.95,S_instar=c(2.660,2.310,1.916,0),
  s_j=0.999,v_init=3e-9,E_init=E_Egg/v_init,E_H_init=0,stage=0,aestivate=0,depress=0.3,
clutchsize=2.,clutch_ab=c(0,0),viviparous=0,minclutch=0,batch=1,photostart=3,photofinish=1,
daylengthstart=12.5,daylengthfinish=13.,photodirs=1,photodirf=0,frogbreed=0,frogstage=0,
reset=0,breedactthresh=1,breedrainthresh=0,breedtempthresh=200,breedtempcum=24*7,ma=1e-4,mi=0,mh=0.5,
container=0,wetmod=0,conth=10,contw=100,contype=1,rainmult=1,continit=0,conthole=0,contonly=1,
contwet=80,wetlandTemps=matrix(data = 0., nrow = 24*dim, ncol = 1),
wetlandDepths=matrix(data = 0., nrow = 24*dim, ncol = 1),
thermal_stages=matrix(data = c(rep(ctmin,8),rep(ctmax,8),rep(TMINPR,8),rep(TMAXPR,8),rep(TBASK,8),
  rep(TPREF,8)), nrow = 8, ncol = 6),
behav_stages=matrix(data = c(rep(dayact,8),rep(nocturn,8),rep(crepus,8),rep(burrow,8),
  rep(shdburrow,8),rep(mindepth,8),rep(maxdepth,8),rep(CkGrShad,8),rep(climb,8),rep(fosorial,8),
  rep(rainact,8),rep(actrainthresh,8),rep(breedactthresh,8),rep(flyer,8)), nrow = 8, ncol = 14),
water_stages=matrix(data = c(rep(skinwet,8),rep(extref,8),rep(PFEWAT,8),rep(PTUREA,8),
  rep(FoodWater,8),rep(minwater,8),rep(raindrink,8),rep(gutfill,8)), nrow = 8, ncol = 8),
arrhenius=matrix(data = matrix(data = c(rep(TA,8),rep(TAL,8),rep(TAH,8),rep(TL,8),rep(TH,8)),
  nrow = 8, ncol = 5), nrow = 8, ncol = 5),
wings=0,rho1_3=0.2,trans1=0.00,aref=0.26,bref=2.04,cref=1.47,phi=179.,phimax=phi,phimin=phi,
flyer=0,flyspeed=5,flymetab=0.1035){
#

  if(lometry==3){
    shape_a<-1.
    shape_b<-1.
    shape_c<-4.
  }
  if(lometry==4){
    shape_a<-1.
    shape_b<-1.
    shape_c<-0.5
  }

  #turn on container model if aquatic egg/larval phase
  if(frogbreed==1 | frogbreed==2){
    container<-1
  }
  if(frogbreed==3){
    container<-0
  }

  # container/pond initial conditons
  contlast<-0.
  templast<-7.

  iyear<-0 #initializing year counter
  countday<-1 #initializing day counter

  if(microin!="none"){
  cat('reading microclimate input \n')
  RAINFALL<-as.matrix(read.csv(file=paste(microin,'rainfall.csv',sep=""),sep=","))[,2]
  dim=length(RAINFALL)
  metout<-read.csv(file=paste(microin,'metout.csv',sep=""),sep=",")[,-1]
  shadmet<-read.csv(file=paste(microin,'shadmet.csv',sep=""),sep=",")[,-1]
  soil<-read.csv(file=paste(microin,'soil.csv',sep=""),sep=",")[,-1]
  shadsoil<-read.csv(file=paste(microin,'shadsoil.csv',sep=""),sep=",")[,-1]
  if(file.exists(paste(microin,'wetlandTemps.csv',sep=""))){
    wetlandTemps<-read.csv(file=paste(microin,'wetlandTemps.csv',sep=""),sep=",")[,-1]
    wetlandDepths<-read.csv(file=paste(microin,'wetlandDepths.csv',sep=""),sep=",")[,-1]
  }else{
    wetlandTemps=matrix(data = 0., nrow = 24*dim, ncol = 1)
    wetlandDepths=matrix(data = 0., nrow = 24*dim, ncol = 1)
  }
    if(file.exists(paste(microin,'soilpot.csv',sep=""))){
      soilpot<-read.csv(file=paste(microin,'soilpot.csv',sep=""),sep=",")[,-1]
      soilmoist<-read.csv(file=paste(microin,'soilmoist.csv',sep=""),sep=",")[,-1]
      shadpot<-read.csv(file=paste(microin,'shadpot.csv',sep=""),sep=",")[,-1]
      shadmoist<-read.csv(file=paste(microin,'shadmoist.csv',sep=""),sep=",")[,-1]
      humid<-read.csv(file=paste(microin,'humid.csv',sep=""),sep=",")[,-1]
      shadhumid<-read.csv(file=paste(microin,'shadhumid.csv',sep=""),sep=",")[,-1]
    }else{
      soilpot<-soil
      soilmoist<-soil
      shadpot<-soil
      shadmoist<-soil
      humid<-soil
      shadhumid<-soil
      soilpot[,3:12]<-0
      soilmoist[,3:12]<-0.5
      shadpot[,3:12]<-0
      shadmoist[,3:12]<-0.5
      humid[,3:12]<-0.99
      shadhumid[,3:12]<-0.99
    }
  metout<-as.matrix(metout)
  shadmet<-as.matrix(shadmet)
  shadsoil<-as.matrix(shadsoil)
  soil<-as.matrix(soil)
  soilmoist<-as.matrix(soilmoist)
  shadmoist<-as.matrix(shadmoist)
  soilpot<-as.matrix(soilpot)
  shadpot<-as.matrix(shadpot)
  humid<-as.matrix(humid)
  shadhumid<-as.matrix(shadhumid)
  ectoin<-read.csv(file=paste(microin,'ectoin.csv',sep=""),sep=",")[,-1]
  DEP<-as.matrix(read.csv(file=paste(microin,'DEP.csv',sep=""),sep=","))[,2]
  MAXSHADES<-as.matrix(read.csv(file=paste(microin,'MAXSHADES.csv',sep=""),sep=","))[,2]
  metout2=matrix(data = 0., nrow = 24*dim, ncol = 18)
  soil2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  shadmet2=matrix(data = 0., nrow = 24*dim, ncol = 18)
  shadsoil2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  soilmoist2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  shadmoist2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  soilpot2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  shadpot2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  humid2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  shadhumid2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  wetlandTemps2=matrix(data = 0., nrow = 24*dim, ncol = 1)
  wetlandDepths2=matrix(data = 0., nrow = 24*dim, ncol = 1)
  metout2[1:nrow(metout),]<-metout
  shadmet2[1:nrow(metout),]<-shadmet
  soil2[1:nrow(metout),]<-soil
  shadsoil2[1:nrow(metout),]<-shadsoil
  soilmoist2[1:nrow(metout),]<-soilmoist
  shadmoist2[1:nrow(metout),]<-shadmoist
  soilpot2[1:nrow(metout),]<-soilpot
  shadpot2[1:nrow(metout),]<-shadpot
  humid2[1:nrow(metout),]<-humid
  shadhumid2[1:nrow(metout),]<-shadhumid
  wetlandTemps2[1:nrow(metout)]<-wetlandTemps
  wetlandDepths2[1:nrow(metout)]<-wetlandDepths
  metout<-metout2
  shadmet<-shadmet2
  soil<-soil2
  shadsoil<-shadsoil2
  soilmoist<-soilmoist2
  shadmoist<-shadmoist2
  soilpot<-soilpot2
  shadpot<-shadpot2
  humid<-humid2
  shadhumid<-shadhumid2
  wetlandTemps<-wetlandTemps2
  wetlandDepths<-wetlandDepths2
  metout.names<-c("JULDAY","TIME","TALOC","TAREF","RHLOC","RH","VLOC","VREF","SOILMOIST3","POOLDEP","TDEEP","ZEN","SOLR","TSKYC","DEW","FROST","SNOWFALL","SNOWDEP")
  colnames(metout)<-metout.names
  colnames(shadmet)<-metout.names
  soil.names<-c("JULDAY","TIME",paste("D",DEP,"cm", sep = ""))
  colnames(soil)<-soil.names
  colnames(shadsoil)<-soil.names
  moist.names<-c("JULDAY","TIME",paste("WC",DEP,"cm", sep = ""))
  colnames(soilmoist)<-moist.names
  colnames(shadmoist)<-moist.names
  pot.names<-c("JULDAY","TIME",paste("PT",DEP,"cm", sep = ""))
  colnames(soilpot)<-pot.names
  colnames(shadpot)<-pot.names
  hum.names<-c("JULDAY","TIME",paste("RH",DEP,"cm", sep = ""))
  colnames(humid)<-hum.names
  colnames(shadhumid)<-hum.names
  }else{
    dim=length(RAINFALL)
  }


  if(soilmoisture1==1){
  soilpotb<-soilpot
  soilmoistb<-soilmoist
  }
  # habitat
  ALT<-ectoin[1] # altitude (m)
  OBJDIS<-1.0 # distance from object (e.g. bush)
  OBJL<-0.0001
  PCTDIF<-0.1 # percent of sunlight that is diffuse (decimal %)
  EMISSK<-1.0 # emissivity of the sky (decimal %)
  EMISSB<-1.0 # emissivity of the substrate (decimal %)
  ABSSB<-1-ectoin[2] # solar absorbtivity of the substrate (decimal %)
  shade<-minshade # shade (%)

  # animal properties
  AMASS<-amass/1000 # animal mass (kg)
  absan<-ABSMAX # animal solar absorbtivity
  RQ<-0.8 # respiratory quotient

  FATOBJ<-0.
  #  if(container==1){
  #    live<-0}else{live<-1
  #  }
  #live<-1
  TIMBAS<-1.
  #  if(container==1){
  #    SKINW<-100.}else{
  SKINW<-skinwet
  #    }
  skint<-0.
  O2gas<-20.95
  CO2gas<-0.03
  N2gas<-79.02
  gas<-c(O2gas,CO2gas,N2gas)
  #  if(container==1){
  #    transt<-1
  #  }else{
  transt<-0
  #  }
  tranin<-1
  tcinit<-metout[1,"TALOC"]

  ACTLVL<-1
  nodnum<-10
  spec<-0. # spectacle covering eye surface? (adds to water loss for lizard/frog/turtle geometry)
  xbas<-1.
  nofood<-0
  tdigpr<-TPREF
  o2max<-extref
  #  if(container==1){
  #  maxshd<-1.
  #  minshd<-0.
  #  }else{
  maxshd<-maxshade
  minshd<-minshade
  #  }
  behav=c(dayact,nocturn,crepus,rainact,burrow,CkGrShad,climb,fosorial,nofood)
  julday<-1

  # DEB model initial conditions
  V_init_baby<-3e-9
  E_init_baby<-E_Egg/V_init_baby
  E_baby_init<-E_init_baby
  V_baby_init<-V_init_baby
  ms_init<-0.
  cumrepro_init<-0.
  q_init<-0.
  hs_init<-0.
  cumbatch_init<-0.
  pregnant<-0
  E_m<-(p_Mref*z/kappa)/v_dotref

  # conversions from percent to proportion
  PTUREA1<-PTUREA/100
  PFEWAT1<-PFEWAT/100
  FoodWater1<-FoodWater/100
  water_stages[,3]<-water_stages[,3]/100
  water_stages[,4]<-water_stages[,4]/100
  water_stages[,5]<-water_stages[,5]/100
  eggmass<-0 # initial dry mass of an egg (g) - no longer used so delete

  #DEB mass balance calculations
  nO<-cbind(nX,nV,nE,nP) # matrix of composition of organics, i.e. food, structure, reserve and faeces
  CHON<-c(12,1,16,14)
  wO<-CHON%*%nO
  w_V=wO[3]
  M_V<-d_V/w_V
  yEX<-kappa_X*mu_X/mu_E # yield of reserve on food
  yXE<-1/yEX # yield of food on reserve
  yVE<-mu_E*M_V/E_G  # yield of structure on reserve
  yPX<-kappa_X_P*mu_X/mu_P # yield of faeces on food
  yXP<-1/yPX # yield of food on faeces
  yPE<-yPX/yEX # yield of faeces on reserve  0.143382353
  nM<-matrix(c(1,0,2,0,0,2,1,0,0,0,2,0,N_waste),nrow=4)
  N_waste_inv<-c(-1*N_waste[1]/N_waste[4],(-1*N_waste[2])/(2*N_waste[4]),(4*N_waste[1]+N_waste[2]-2*N_waste[3])/(4*N_waste[4]),1/N_waste[4])
  nM_inv<-matrix(c(1,0,-1,0,0,1/2,-1/4,0,0,0,1/2,0,N_waste_inv),nrow=4)
  JM_JO<--1*nM_inv%*%nO
  etaO<-matrix(c(yXE/mu_E*-1,0,1/mu_E,yPE/mu_E,0,0,-1/mu_E,0,0,yVE/mu_E,-1/mu_E,0),nrow=4)
  w_N<-CHON%*%N_waste

  lat<-ectoin[4]
  if(soilmoisture1==1){
  humid[,3:9]<-metout[,5]/100 # assume ambient humidity down to 30cm
  shadhumid[,3:9]<-shadmet[,5]/100 # assume ambient humidity down to 30cm
  humid[,7:12]<-0.8 # assume higher humidity in burrow, 60cm and lower
  shadhumid[,7:12]<-0.8 # assume higher humidity in burrow, 60cm and lower


  grassgrowths<-as.data.frame(soilpotb)
  soilmoist2b<-as.data.frame(soilmoistb)
  soilmoist2b<-subset(soilmoist2b,soilmoist2b$TIME==720)
  grassgrowths<-subset(grassgrowths,soilmoist2b$TIME==720)
  grassgrowths<-grassgrowths$PT5cm # assume plant growth driven by 5cm depth

    grow<-grassgrowths
    grow[grow>-1500]<-1 # find times when below permanent wilting point
    grow[grow<=-1500]<-0
    counter<-0
    grow2<-grow*0
      for(j in 1:length(grow)){
        if(j==1){
            if(grow[j]==1){
            counter<-counter+1
            }
          grow2[j]<-counter
        }else{
          if(grow[j-1]>0 & grow[j]==1){
            counter<-counter+1
          }else{
            counter<-0
          }
          grow2[j]<-counter
        }
      }
     grow3<-grow2
     grow3[grow3<7]<-0 # use one week in a row as time required for plats to come back after PWP has been hit
     grow3[grow3>0]<-1 # make vector of 0 and 1 where 1 means plants could have come back from drought

  soilmoist2b<-soilmoist2b$WC5cm
  grassgrowths<-as.data.frame(cbind(grassgrowths,soilmoist2b))
  colnames(grassgrowths)<-c('pot','moist')
  grassgrowths$pot[grassgrowths$pot>-200]<-FoodWater # assume plants start wilting at about 2 bar, but above this they are at max water content
  grassgrowths$moist<-grassgrowths$moist*100 # convert to percent
  potmult<-grassgrowths$pot
  potmult[potmult!=82]<-0
  potmult[potmult!=0]<-1
  wilting<-subset(grassgrowths,grassgrowths$pot==FoodWater) # find soil moisture range corresponding to values above the wilting point
  wilting<-min(wilting$moist) # get the min soil moisture at which plants aren't wilting
  grassgrowths<-grassgrowths$moist
  grassgrowths[grassgrowths>wilting]<-FoodWater # now have vector of either max plant water content or soil moisture content - need to convert the latter into a smooth decline to zero from max value
  minmoist<-0
  grassgrowths[grassgrowths<FoodWater]<-(grassgrowths[grassgrowths<FoodWater]-minmoist)/(wilting-minmoist)*FoodWater # for just the values less than max water content, make them equal to the
  grassgrowths<-grassgrowths/100*grow3
  grasstsdms<-grassgrowths
  }else{
    grassgrowths<-rep(FoodWater,nrow(metout))
    grasstsdms<-grassgrowths
  }
  julstart<-metout[1,2]
  tannul<-as.numeric(mean(soil[,12]))
  monthly<-0
  tester<-0
  microyear<-1

  # bucket model for soil moisture
  fieldcap<-ectoin[5]# %vol, water content at 0.1ba = 10kPa
  fieldcap<-30 # field capacity, m3/m3*100
  if(soilmoisture1==1){
    conth<-fieldcap/10 # containter height, cm
    contw<-100
    contype<-1 # is 'containter' sitting on the surface, like a bucket (0) or sunk into the ground like a pond (1)
    rainmult<-0.3 # !!!!!!!!!!!!!!rainfall multiplier to reflect catchment (don't make this zero unless you want a drought!)
    continit<-0 # initial container water level (cm)
    conthole<-0#2.8 # daily loss of height (mm) due to 'hole' in container (e.g. infiltration to soil, drawdown from water tank)
    contwet<- 2 # percent wet value for container
  }
  ectoinput<-as.matrix(c(ALT,FLTYPE,OBJDIS,OBJL,PCTDIF,EMISSK,EMISSB,ABSSB,shade,enberr,AMASS,EMISAN,absan,RQ,rinsul,lometry,live,TIMBAS,Flshcond,Spheat,Andens,ABSMAX,ABSMIN,FATOSK,FATOSB,FATOBJ,TMAXPR,TMINPR,DELTAR,SKINW,spec,xbas,extref,TPREF,ptcond,skint,gas,transt,soilnode,o2max,ACTLVL,tannul,nodnum,tdigpr,maxshd,minshd,ctmax,ctmin,behav,julday,actrainthresh,viviparous,pregnant,conth,contw,contlast,tranin,tcinit,nyears,lat,rainmult,julstart,monthly,customallom,MR_1,MR_2,MR_3,DEB,tester,rho1_3,trans1,aref,bref,cref,phi,wings,phimax,phimin,shape_a,shape_b,shape_c,minwater,microyear,container,flyer,flyspeed,dim,maxdepth,ctminthresh,ctkill,gutfill,mindepth,TBASK,TEMERGE,p_Xm,SUBTK,flymetab,continit,wetmod,contonly,conthole,contype,shdburrow,breedtempthresh,breedtempcum,contwet))
  debmod<-c(clutchsize,andens_deb,d_V,eggdryfrac,mu_X,mu_E,mu_V,mu_P,T_REF,z,kappa,kappa_X,p_Mref,v_dotref,E_G,k_R,MsM,delta,h_aref,V_init_baby,E_init_baby,k_J,E_Hb,E_Hj,E_Hp,clutch_ab[2],batch,breedrainthresh,photostart,photofinish,daylengthstart,daylengthfinish,photodirs,photodirf,clutch_ab[1],frogbreed,frogstage,etaO,JM_JO,E_Egg,kappa_X_P,PTUREA1,PFEWAT1,wO,w_N,FoodWater1,f,s_G,K,X,metab_mode,stages,y_EV_l,s_j,startday,raindrink,reset,ma,mi,mh,aestivate,depress,minclutch)
  deblast<-c(iyear,countday,v_init,E_init,ms_init,cumrepro_init,q_init,hs_init,cumbatch_init,V_baby_init,E_baby_init,E_H_init,stage)

  origjulday<-metout[,1]
  if(ystrt>0){
    metout<-rbind(metout[((ystrt)*365*24+1):(dim*24),],metout[1:((ystrt)*365*24),])
    shadmet<-rbind(shadmet[((ystrt)*365*24+1):(dim*24),],shadmet[1:((ystrt)*365*24),])
    soil<-rbind(soil[((ystrt)*365*24+1):(dim*24),],soil[1:((ystrt)*365*24),])
    shadsoil<-rbind(shadsoil[((ystrt)*365*24+1):(dim*24),],shadsoil[1:((ystrt)*365*24),])
    soilmoist<-rbind(soilmoist[((ystrt)*365*24+1):(dim*24),],soilmoist[1:((ystrt)*365*24),])
    shadmoist<-rbind(shadmoist[((ystrt)*365*24+1):(dim*24),],shadmoist[1:((ystrt)*365*24),])
    soilpot<-rbind(soilpot[((ystrt)*365*24+1):(dim*24),],soilpot[1:((ystrt)*365*24),])
    shadpot<-rbind(shadpot[((ystrt)*365*24+1):(dim*24),],shadpot[1:((ystrt)*365*24),])
    humid<-rbind(humid[((ystrt)*365*24+1):(dim*24),],humid[1:((ystrt)*365*24),])
    shadhumid<-rbind(shadhumid[((ystrt)*365*24+1):(dim*24),],shadhumid[1:((ystrt)*365*24),])
    wetlandDepths<-c(wetlandDepths[((ystrt)*365*24+1):(dim*24)],wetlandDepths[1:((ystrt)*365*24)])
    wetlandTemps<-c(wetlandTemps[((ystrt)*365*24+1):(dim*24)],wetlandTemps[1:((ystrt)*365*24)])
    MAXSHADES<-c(MAXSHADES[((ystrt)*365+1):(dim)],MAXSHADES[1:((ystrt)*365)])
    RAINFALL<-c(RAINFALL[((ystrt)*365+1):(dim)],RAINFALL[1:((ystrt)*365)])
    grassgrowths<-c(grassgrowths[((ystrt)*365+1):(dim)],grassgrowths[1:((ystrt)*365)])
    metout[,1]<-origjulday
    shadmet[,1]<-origjulday
    soil[,1]<-origjulday
    shadsoil[,1]<-origjulday
    soilmoist[,1]<-origjulday
    shadmoist[,1]<-origjulday
    soilpot[,1]<-origjulday
    shadpot[,1]<-origjulday
    humid[,1]<-origjulday
    shadhumid[,1]<-origjulday
  }


  # code to determine wet periods for activity in a pond

  if(wetmod==1){
  wet_thresh<-10*24 # threshold pond duration
  wet_depth<-100 # threshold pond depth (mm)
  wet_temp<-28 # threshold exit temp (deg C)
  b<-cbind(as.data.frame(wetlandDepths),as.data.frame(wetlandTemps))
  colnames(b)<-c('depth','temp')
  b$depth[b$temp>wet_temp]<-0
  b<-b$depth
  b[b>=wet_depth]<-1
  b[b!=1]<-0
  bb<-rle(b)
  bb$values[bb$lengths<wet_thresh]<-0
  c<-b*0
  values<-bb$values
  lengths<-bb$lengths
  for(k in 1:length(bb$values)){
    d<-c(rep(values[k],lengths[k]))
    if(k==1){
    e<-d
    }else{
    e<-c(e,d)
    }
  }
  wetlandDepths<-wetlandDepths*e
  }

  if(write_input==1){
      if(dir.exists("ecto csv input")==FALSE){
        dir.create("ecto csv input")
      }
    cat('writing input csv files \n')
    write.csv(ectoinput, file = "ecto csv input/ectoinput.csv")
    write.csv(debmod, file = "ecto csv input/debmod.csv")
    write.csv(deblast, file = "ecto csv input/deblast.csv")
    write.csv(RAINFALL, file = "ecto csv input/rainfall.csv")
    write.csv(DEP, file = "ecto csv input/dep.csv")
    write.csv(grassgrowths, file = "ecto csv input/grassgrowths.csv")
    write.csv(grasstsdms, file = "ecto csv input/grasstsdms.csv")
    write.csv(wetlandTemps, file = "ecto csv input/wetlandTemps.csv")
    write.csv(wetlandDepths, file = "ecto csv input/wetlandDepths.csv")
    write.csv(arrhenius, file = "ecto csv input/arrhenius.csv")
    write.csv(thermal_stages, file = "ecto csv input/thermal_stages.csv")
    write.csv(behav_stages, file = "ecto csv input/behav_stages.csv")
    write.csv(water_stages, file = "ecto csv input/water_stages.csv")
    write.csv(MAXSHADES, file = "ecto csv input/Maxshades.csv")
    write.csv(S_instar, file = "ecto csv input/S_instar.csv")
    write.table(metout[(seq(1,dim*24)),], file = "ecto csv input/metout.csv",sep=",",row.names=FALSE)
    write.table(shadmet[(seq(1,dim*24)),], file = "ecto csv input/shadmet.csv",sep=",",row.names=FALSE)
    write.table(soil[(seq(1,dim*24)),], file = "ecto csv input/soil.csv",sep=",",row.names=FALSE)
    write.table(shadsoil[(seq(1,dim*24)),], file = "ecto csv input/shadsoil.csv",sep=",",row.names=FALSE)
    write.table(soilmoist[(seq(1,dim*24)),], file = "ecto csv input/soilmoist.csv",sep=",",row.names=FALSE)
    write.table(shadmoist[(seq(1,dim*24)),], file = "ecto csv input/shadmoist.csv",sep=",",row.names=FALSE)
    write.table(soilpot[(seq(1,dim*24)),], file = "ecto csv input/soilpot.csv",sep=",",row.names=FALSE)
    write.table(shadpot[(seq(1,dim*24)),], file = "ecto csv input/shadpot.csv",sep=",",row.names=FALSE)
    write.table(humid[(seq(1,dim*24)),], file = "ecto csv input/humid.csv",sep=",",row.names=FALSE)
    write.table(shadhumid[(seq(1,dim*24)),], file = "ecto csv input/shadhumid.csv",sep=",",row.names=FALSE)
  }
  ecto<-list(dim=dim,ectoinput=ectoinput,metout=metout,shadmet=shadmet,soil=soil,shadsoil=shadsoil,soilmoist=soilmoist,shadmoist=shadmoist,soilpot=soilpot,shadpot=shadpot,humid=humid,shadhumid=shadhumid,DEP=DEP,RAINFALL=RAINFALL,iyear=iyear,countday=countday,debmod=debmod,deblast=deblast,grassgrowths=grassgrowths,grasstsdms=grasstsdms,wetlandTemps=wetlandTemps,wetlandDepths=wetlandDepths,arrhenius=arrhenius,thermal_stages=thermal_stages,behav_stages=behav_stages,water_stages=water_stages,MAXSHADES=MAXSHADES,S_instar=S_instar)

  cat('running ectotherm model ... \n')

  ptm <- proc.time() # Start timing
  ectout<-ectorun(ecto)
  print(proc.time() - ptm) # Stop the clock


  environ<-ectout$environ[1:(dim*24),]
  enbal<-ectout$enbal[1:(dim*24),]
  masbal<-ectout$masbal[1:(dim*24),]
  debout<-ectout$debout[1:(dim*24),]
  yearout<-ectout$yearout
  yearsout<-ectout$yearsout[1:nyears,]

  if(DEB==0){
    return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,soilpot=soilpot,shadpot=shadpot,humid=humid,shadhumid=shadhumid,RAINFALL=RAINFALL,enbal=enbal,environ=environ,masbal=masbal,yearout=yearout,yearsout=yearsout,grassgrowths=grassgrowths,grasstsdms=grasstsdms,TMAXPR=TMAXPR,TMINPR=TMINPR,ctmax=ctmax,ctmin=ctmin,TBASK=TBASK,TEMERGE=TEMERGE))
  }else{
    return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,soilpot=soilpot,shadpot=shadpot,humid=humid,shadhumid=shadhumid,RAINFALL=RAINFALL,enbal=enbal,masbal=masbal,environ=environ,debout=debout,yearout=yearout,yearsout=yearsout,grassgrowths=grassgrowths,grasstsdms=grasstsdms,TMAXPR=TMAXPR,TMINPR=TMINPR,ctmax=ctmax,ctmin=ctmin,TBASK=TBASK,TEMERGE=TEMERGE))
  }

}
