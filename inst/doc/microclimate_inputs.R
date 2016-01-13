## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
 eval = TRUE
)

## ----table1, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
writecsv                  | *-*               | 0 (off) or 1 (on) | make Fortran program write output as csv files  
microdaily                | *-*               | 0 (off) or 1 (on) | run in daily mode (initial conditions from previous day)
runshade                  | *-*               | 0 (off) or 1 (on) | run the model twice, once for each shade level
runmoist                  | *-*               | 0 (off) or 1 (on) | run soil moisture model
snowmodel                 | *-*               | 0 (off) or 1 (on) | run the snow model
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table2, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
julnum                    | *days*            | positive integer  | number of days to run the model
julday                    | *day of year*     | 1-365             | vector of julian days (length must equal julnum)
idayst                    | *-*               | 1-julnum          | start day (usually 1)
ida                       | *-*               | 1-julnum          | end day (usually value of julnum)
HEMIS                     | *-*               | 1 (N) or 2 (S)    | hemisphere to run
ALAT                      | *degrees*         | 0-90              | latitude (degrees)
AMINUT                    | *dec. minutes*    | 0-60              | latitude (minutes)
ALONG                     | *degrees*         | 0-180             | longitude (degrees)
ALMINT                    | *dec. minutes*    | 0-60              | longitude (minutes)
ALREF                     | *degrees*         | 0-180             | reference longitude (degrees) for time zone
EC                        | *-*               | 0.0034 to 0.058   | eccenricity of the earth's orbit (presently 0.0167238)
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----table3, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
*Name*                    | *Units *          |	*Allowed Range*   |	*Description*
------------------------- | ----------------- | ----------------- | --------------
RUF                       | *m*               | 0.01-200          | roughness height
Refhyt                    | *m*               | 0.50-10           | reference measurement height for air temp, wind speed and humidity input data
Usrhyt                    | *m*               | > 0.005, < Refhyt | local height at which to compute air temp, wind speed and humidity
Z01                       | *m*               | > Z02, < ZH2      | end day (usually value of julnum)
Z02                       | *m*               | > RUF, < Z01      | hemisphere to run
ZH1                       | *m*               | > ZH2, < Refhyt   | experimental reference height 2 (top)
ZH2                       | *m*               | > RUF, < ZH1      | latitude (minutes)
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

