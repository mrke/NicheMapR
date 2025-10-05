#' Atmospheric pressure from reference elevation
#'
#' Function to get atmospheric pressure at an elevation of interest (h), given
#' reference level pressure, temperature, lapse rate and molar mass of air
#' @param h_ref = 0, reference elevation, m
#' @param P_ref = 101325, static pressure at sea level, Pa
#' @param L_ref = -0.0065, reference temperature lapse rate, K/m
#' @param T_ref = 288, reference temperature (at sea level), K
#' @param g_0 = 9.80665, gravitational acceleration constant, m/s^2
#' @param M = 0.0289644, molar mass of air, kg/mol
#' @return pressure at focal height, Pa
#' @usage get_pressure(h = 1000)
#' @export
get_pressure <- function(h,
                         h_ref = 0,
                         P_ref = 101325,
                         L_ref = -0.0065,
                         T_ref = 288,
                         g_0 = 9.80665,
                         M = 0.0289644){

  # constants
  R <- 8.31446261815324 # universal gas constant, (N m)/(mol K)
  P_a <- P_ref * (1 + (L_ref / T_ref) * h) ^ ((-g_0 * M) / (R * L_ref))
  return(P_a)
}
