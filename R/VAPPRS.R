#' VAPPRS
#'
#' Calculates saturation vapour pressure for a given air temperature.
#' @param db Dry bulb temperature (degrees C)
#' @return esat Saturation vapour pressure (Pa)
#' @export
VAPPRS <- function(db=db){
  t=db+273.16
  loge=t
  loge[t<=273.16]=-9.09718*(273.16/t[t<=273.16]-1.)-3.56654*log10(273.16/t[t<=273.16])+.876793*(1.-t[t<=273.16]/273.16)+log10(6.1071)
  loge[t>273.16]=-7.90298*(373.16/t[t>273.16]-1.)+5.02808*log10(373.16/t[t>273.16])-1.3816E-07*(10.^(11.344*(1.-t[t>273.16]/373.16))-1.)+8.1328E-03*(10.^(-3.49149*(373.16/t[t>273.16]-1.))-1.)+log10(1013.246)
  esat=(10.^loge)*100
  return(esat)
}
