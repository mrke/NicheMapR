#' WETAIR
#'
#' Calculates several properties of humid air as output variables below. The program
#' is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
#' Institution Press. Washington, DC. WETAIR must be used in conjunction with function VAPPRS.
#'
#' Input variables are shown below. The user must supply known values for DB and BP (BP at one standard
#' atmosphere is 101 325 pascals). Values for the remaining variables are determined by whether the user has
#' either (1) psychrometric data (WB or RH), or (2) hygrometric data (DP)
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
  return(list(e=e, esat=esat, vd=vd, rw=rw, tvinc=tvinc, denair=denair, cp=cp, wtrpot=wtrpot, rh=rh))
}
