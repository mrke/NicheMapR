## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ----table1, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
writecsv                  | -                 | 0 (off) or 1 (on) | make Fortran program write output as csv files  
microdaily                | -                 | 0 (off) or 1 (on) | run in daily mode (initial conditions from previous day)
runshade                  | -                 | 0 (off) or 1 (on) | run the model twice, once for each shade level
runmoist                  | -                 | 0 (off) or 1 (on) | run soil moisture model
snowmodel                 | -                 | 0 (off) or 1 (on) | run the snow model
hourly                    | -                 | 0 (off) or 1 (on) | run the model from hourly weather inputs
IR                        | -                 | 0 or 1            | longwave radiation algorithm
message                   | -                 | 0 or 1            | integrator messages
fail                      | -                 | integer            | integrator failure count before quitting
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table2, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
doynum                    | days              | positive integer  | number of days to run the model
doy                       | day-of-year       | 1-365             | vector of days of year (length must equal doynum)
idayst                    | -                 | 1-doynum          | start day (usually 1)
ida                       | -                 | 1-doynum          | end day (usually value of doynum)
HEMIS                     | -                 | 1 (N) or 2 (S)    | hemisphere to run
ALAT                      | degrees           | 0-90              | latitude (degrees)
AMINUT                    | dec. minutes      | 0-60              | latitude (minutes)
ALONG                     | degrees           | 0-180             | longitude (degrees)
ALMINT                    | dec. minutes      | 0-60              | longitude (minutes)
ALREF                     | degrees           | 0-180             | reference longitude (degrees) for time zone
EC                        | -                 | 0.0034 to 0.058   | eccenricity of the earth's orbit (presently 0.0167238)
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table3, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
RUF                       | m                 | 0.01-200          | roughness height
Refhyt                    | m                 | 0.50-10           | reference height for air temp, wind speed and humidity input data
Usrhyt                    | m                 | > 0.005, < Refhyt | local height at which to compute air temp, wind speed and humidity 
ZH                       | m                 | > 0 (or else not used)   | heat transfer roughness height (m)
D0                       | m                 | > 0                      | zero plane displacement correction factor (m) used with ZH
Z01                       | m                 | 0 or (> Z02, < ZH2)   | 1st segment, roughness height$^1$
Z02                       | m                 | 0 or (> RUF, < Z01)   | 2nd segment, roughness height$^1$
ZH1                       | m                 | 0 or (> ZH2, < Refhyt)| 1st segment, height above surface$^1$
ZH2                       | m                 | 0 or (> RUF, < ZH1)   | 2nd segment, height above surface$^1$
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table4, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
SLES                      | -                 | 0-1               | substrate longwave IR emissivity$^1$ 
REFLS                      | -                 | 0-1               | substrate solar reflectivity$^1$ 
CMH2O                     | cm                | 0.1-2             | precipitable cm H$_2$O in air column
TAI                       | -                 | 1-doynum          | aerosol optical extinction coefficients$^{1,2}$ 
solonly                   | -                | 0-1             | flag to do only solar calcs
lamb                     | -                | 0-1             | flag to return wavelength-specific solar radiation output
IUV                      | -                | 0-1             | flag to use gamma function for scattered solar radiation$^{3}$
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table5, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
ALTT                      | m                 | 0-10,000          | elevation
slope                     | decimal degrees   | 0-90              | slope
azmuth                    | decimal degrees   | 0-360             | aspect (0 is north)
hori                      | decimal degrees   | 0-90              | horizon angles$^1$  
VIEWF                     | -                 | 0-1               | view factor to sky$^2$ 
MINSHADES                 | %                 | 0-100             | minimum shade$^3$  
MAXSHADES                 | %                 | 0-100             | maximum shade$^{3,4}$ 
PCTWETS                   | %                 | 0-100             | percentage of substrate unit area that is wet$^{3,5}$  
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table6, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
DEP                       | cm                | 0-1000            | vector of 10 depths$^1$
ERR                       | -                 | >0                | integrator error (typically 1.5-2)
tannul                    | &deg;C            | -80 - +60         | annual mean temperature
soilinit                  | &deg;C            | -80 - +60         | initial substrate temperature profile$^1$
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table7, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
TIMINS                    | h                 | 0-23              | time of minima for air temp, wind, humidity and cloud cover$^1$
TIMAXS                    | h                 | 0-23              | time of maxima for air temp, wind, humidity and cloud cover$^2$
TMINN                     | &deg;C            | -80 - +60         | minimum air temperature (at reference height, Refhyt)$^3$
TMAXX                     | &deg;C            | -80 - +60         | maximum air temperature (at reference height, Refhyt)$^3$
RHMINN                    | %                 | 0-100             | minimum relative humidity (at reference height, Refhyt)$^3$
RHMAXX                    | %                 | 0-100             | maximum relative humidity (at reference height, Refhyt)$^3$
WNMINN                    | m s$^{-1}$        | 0-100             | minimum wind speed (at reference height, Refhyt)$^3$
WNMAXX                    | m s$^{-1}$        | 0-100             | maximum wind speed (at reference height, Refhyt)$^3$
CCMINN                    | %                 | 0-100             | minimum cloud cover; vector of length = doynum
CCMAXX                    | %                 | 0-100             | maximum cloud cover; vector of length = doynum
RAINFALL                  | mm                | 0-2000            | daily total rainfall; vector of length = doynum
tannulrun                 | &deg;C            | -80 - +60         | daily deep soil temperature; vector of length = doynum
moists                    | decimal %         | 0-1               | predefined soil daily moisture profile through time$^4$
TAIRhr                    | &deg;C            | -80 - +60         | hourly air temperature (at reference height, Refhyt)$^5$
RHhr                      | %                 | 0-100             | hourly relative humidity (at reference height, Refhyt)$^5$
WNhr                      | m s$^{-1}$        | 0-100             | hourly wind speed (at reference height, Refhyt)$^5$
CLDhr                     | %                 | 0-100             | hourly cloud cover$^5$
SOLRhr                    | W m$^{2}$         | 0-1367            | hourly solar radiation$^5$
RAINhr                    | mm                | 0-2000            | hourly rainfall$^5$
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table8, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *             |	*Allowed Range*   |	*Description*
------------------------- | -------------------- | ------------------ | --------------
Numtyps                   | -                    | 1-10               | number of substrate types
Nodes$^{1,2}$             | -                    | 1-10               | nodes where substrate type transitions occur
soilprops[,1]$^{3}$       | Mg m$^{-3}$          | >0                 | bulk density (must not exceed mineral density)
soilprops[,2]$^{3}$       | m$^{3}$ m$^{-3}$     | >0                 | volumetric water content at saturation$^{4}$ 
soilprops[,3]$^{3}$       | W m$^{-1}$ K$^{-1}$  | >0                 | thermal conductivity
soilprops[,4]$^{3}$       | J kg$^{-1}$ K$^{-1}$ | >0                 | specific heat capacity
soilprops[,5]$^{3}$       | Mg m$^{-3}$          | >0                 | mineral density
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table9, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *             |	*Allowed Range*  |	*Description*
------------------------- | -------------------- | ----------------- | --------------
PE                        | J kg$^{-1}$          | 0.7-3.7           | air entry water potential; vector of 19 values$^{1}$
KS                        | kg s m$^{-3}$        | 1 x $10^{-5}$     | saturated conductivity; vector of 19 values$^{1}$
BB                        | -                    | 1.7-7.6           | Campbell's 'b' parameter; vector of 19 values$^{1}$
BD                        | Mg m$^{-3}$          | >0                | bulk density (should be matched to same variable in the soilprops matrix)
DD                        | Mg m$^{-3}$          | >0                | soil mineral density (should be matched to same variable in the soilprops matrix)
L                         | m$^{3}$ m$^{-3}$     | 0-90              | root density; vector of 19 values$^{1}$  
R1                        | m                    | 0-90              | root radius  
RW                        | m$^{3}$ kg$^{-1}$ s$^{-1}$     | 0-90              | resistance per unit length of root  
RL                        | m$^{3}$ kg$^{-1}$ s$^{-1}$     | 0-90              | resistance per unit length of leaf  
PC                        | J kg$^{-1}$          | 0-90              | critical leaf water potential for stomatal closure  
SP                        | -                    | >0                | stability parameter for stomatal closure equation
IM                        | kg                   | >0                | maximum allowable mass balance error
MAXCOUNT                  | -                    | >0                | maximum iterations for mass balance
LAI                       | -                    | 0-90              | leaf area index, used to partition traspiration/evaporation from PET 
rainmult                  | -                    | >0                | rainfall multiplier to impose a catchment 
maxpool                   | mm                   | >0                | max depth for surface water pooling$^{2}$  
evenrain                  | -                    | 1 or 2            | even rainfall over 24hrs (1) or one event at midnight (2)
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table10, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *             |	*Allowed Range*   |	*Description*
------------------------- | -------------------- | ------------------ | --------------
snowtemp                  | &deg;C               | -80 - 80           | temperature at which precipitation falls as snow
snowdens                  | Mg m$^{-3}$          | 0.05-1             | snow density
densfun                   | -                    | -                  | parameters of snow density function$^{1}$ 
snowmelt                  | decimal %            | 0-1                | proportion of calculated snowmelt that doesn't refreeze  
undercatch                | -                    | >=1                | undercatch multipier for converting rainfall to snow 
rainmelt                  | -                    | 0-2                | parameter for rain melting snow$^{2}$ 
snowcond                  | -                    | >0                 | effective thermal conductivity of snow $W/mK$ 
intercept                 | -                    | 0-1                | proportion of snow intercepted by vegetation under shaded conditions 
grasshade                 | -                    | 0 or 1             | if 1, shade is not present when snow is present (because shade is cast by grass/low veg) 
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table11, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *             |	*Allowed Range*   |	*Description*
------------------------- | -------------------- | ------------------ | --------------
tides[1:doynum,1]         | -                    | 0 or 1             | tide in (1) or out(0); vector of length doynum
tides[1:doynum,2]         | &deg;C               | -80 - 80           | sea water temperature; vector of length doynum
tides[1:doynum,3]         | %                    | -                  | % surface wetness from wave splash; vector of length doynum
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

