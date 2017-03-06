#' DRYAIR
#'
#' Calculates several properties of dry air and related characteristics shownas output variables below. The program
#' is based on equations from List, R. J. 1971. Smithsonian Meteorological Tables. Smithsonian
#' Institution Press. Washington, DC. WETAIR must be used in conjunction with function VAPPRS.
#'
#' The user must supply values for the input variables (db, bp and alt).
#' If alt is known (-1000 < alt < 20000) but not BP, then set BP=0
#' @param db Dry bulb temperature (degrees C)
#' @param bp Barometric pressure (pascal)
#' @param alt Altitude (m)
#' @return patmos Standard atmospheric pressure (Pa)
#' @return densty Density (kg m-3)
#' @return visdyn Dynamic viscosity (kg m-1 s-1)
#' @return viskin Kinematic viscosity (m2 s-1)
#' @return difvpr Diffusivity of water vapour in air (m2 s-1)
#' @return thcond Thermal conductivity (W m-1 K-1)
#' @return htovpr Latent heat of vapourisation of water (J kg-1)
#' @return tcoeff Temperature coefficient of volume expansion (K-1)
#' @return ggroup Group of variables in Grashof number (1-m3 -K)
#' @return bbemit black body emittance (W m-2)
#' @return emtmax Wave length of maximum emittance (m)
#' @export
DRYAIR <- function(db=db, bp=0, alt=0){
  tstd=273.15
  pstd=101325.
  patmos=pstd*((1-(0.0065*alt/288))^(1/0.190284))
  bp=rep(bp,length(patmos))
  bp[bp<=0]=patmos[bp<=0]
  densty=bp/(287.04*(db+tstd))
  visnot=1.8325E-5
  tnot=296.16
  c=120.
  visdyn=(visnot*((tnot+c)/(db+tstd+c)))*(((db+tstd)/tnot)^1.5)
  viskin=visdyn/densty
  difvpr=2.26E-5*(((db+tstd)/tstd)^1.81)*(1.E5/bp)
  thcond=0.02425+(7.038E-5*db)
  htovpr=2.5012E6-2.3787E3*db
  tcoeff=1./(db+tstd)
  ggroup=0.0980616*tcoeff/(viskin*viskin)
  bbemit=5.670367E-8*((db+tstd)^4)
  emtmax=2.897E-3/(db+tstd)
  return(list(patmos=patmos, densty=densty, visdyn=visdyn, viskin=viskin, difvpr=difvpr,
    thcond=thcond, htovpr=htovpr, tcoeff=tcoeff, ggroup=ggroup, bbemit=bbemit, emtmax=emtmax))
}
