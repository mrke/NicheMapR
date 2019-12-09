#' WETAIR
#'
#' Calculates several properties of humid air as output variables below. The program
#' is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
#' Institution Press. Washington, DC. WETAIR must be used in conjunction with function VAPPRS.
#'
#' Input variables are shown below. The user must supply known values for db and bp (bp at one standard
#' atmosphere is 101 325 pascals). Values for the remaining variables are determined by whether the user has
#' either (1) psychrometric data (wb or rh), or (2) hygrometric data (dp)
#'
#' (1) Psychrometric data:
#' If wb is known but not rh, then set rh=-1 and dp=999
#' If rh is known but not wb then set wb=0 and dp=999
#'
#' (2) Hygrometric data:
#' If dp is known then set wb = 0 and rh = 0.
#' @param db Dry bulb temperature (degrees C)
#' @param wb Wet bulb temperature (degrees C)
#' @param rh Relative humidity (\%)
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
  return(list(e=e, esat=esat, vd=vd, rw=rw, tvinc=tvinc, denair=denair, cp=cp, wtrpot=wtrpot, rh=rh))
}
