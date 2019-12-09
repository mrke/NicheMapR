#' WETAIR.dp
#'
#' Calculates several properties of humid air as output variables below. The program
#' is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
#' Institution Press. Washington, DC. WETAIR must be used in conjunction with function VAPPRS.
#'
#' Input variables are shown below. The user must supply known values for dp, db and pb (bp at one standard
#' atmosphere is 101 325 pascals).
#' @param db Dry bulb temperature (degrees C)
#' @param dp Dew point temperature (degrees C)
#' @param bp Barometric pressure (pascal)
#' @return e Vapour pressure (Pa)
#' @return esat Saturation vapour pressure (Pa)
#' @return vd Vapour density (kg m-3)
#' @return rw Mixing ratio (kg kg-1)
#' @return tvir Virtual temperature (K)
#' @return tvinc Virtual temperature increment (K)
#' @return denair Density of the air (kg m-3)
#' @return cp Specific heat of air at constant pressure (J kg-1 K-1)
#' @return wtrpot Water potential (Pa)
#' @return rh Relative humidity (\%)
#' @export
WETAIR.dp <- function(db=db, dp=dp, bp=101325){
  tk = db + 273.15
  esat = VAPPRS(db)
  e = VAPPRS(dp)
  rh = (e / esat) * 100
  rw = ((0.62197 * 1.0053 * e) / (bp - 1.0053 * e))
  vd = e * 0.018016 / (0.998 * 8.31434 * tk)
  tvir = tk * ((1.0 + rw / (18.016 / 28.966)) / (1.0 + rw))
  tvinc = tvir - tk
  denair = 0.0034838 * bp / (0.999 * tvir)
  cp = (1004.84 + (rw * 1846.40)) / (1.0 + rw)
  wtrpot = 4.615e+5 * tk * log(rh / 100.0)
  return(list(e=e, esat=esat, vd=vd, rw=rw, tvinc=tvinc, denair=denair, cp=cp, wtrpot=wtrpot, rh=rh))
}
