#' plantgro
#'
#' Function to compute plant water content given threshold values of soil water potential at which wilting point and permanent wilting point occurs
#' @param soilpot = micro$soilpot, matrix of soil water potential from NicheMapR microclimate model to use for calculations
#' @param soilmoist = micro$soilmoist, matrix of soil moisture to use for calculations
#' @param plant_temp = micro$soil[,6], vector of plant temperature to use for calculations
#' @param root_shallow = 4, shallowest soil node to use when getting average soil moisture/water potential for calculating plant moisture/growth
#' @param root_deep = 4, deepest soil node to use when getting average soil moisture/water potential for calculating plant moisture/growth
#' @param growth_delay = 0, time required for plants to recover after hitting permanent wilting point, days
#' @param temp_thresh = 15, growth temperature threshold for degree day calculation, °C
#' @param wilting_point = -200, soil water potential at wilting point, J/kg
#' @param permanent_wilting_point = -1500, soil water potential at permanent wilting point, J/kg
#' @param foodwater = 82, Maximum water conten of plant, \%
#' @export
plantgro<-function(soilpot=micro$soilpot, soilmoist=micro$soilmoist, plant_temp=micro$soil[,6], root_shallow=4, root_deep=8, temp_thresh=15, growth_delay=1, wilting_point=-200, permanent_wilting_point=-1500, foodwater=82){

  if(root_shallow == root_deep){
    meanpot <- soilpot[, ((root_shallow + 2):(root_deep + 2))]
    meanmoist <- soilmoist[, ((root_shallow + 2):(root_deep + 2))]  # get range of soil water moisture depths to take mean of
  }else{
    meanpot <- soilpot[, ((root_shallow + 2):(root_deep + 2))] # get range of soil water potential depths to take mean of
    meanpot <- apply(meanpot, 1, mean) # get average soil water potential across chosen depth range
    meanmoist <- soilmoist[, ((root_shallow + 2):(root_deep + 2))]  # get range of soil water moisture depths to take mean of
    meanmoist <- apply(meanmoist, 1, mean)
  }

  # compute plant presence
  pres <- meanpot
  pres[meanpot > permanent_wilting_point] <- 1 # find times above the PWP
  pres[meanpot <= permanent_wilting_point] <- 0 # find times when below the PWP (plant dead)
  plant.pres <- ave(pres, cumsum(pres == 0), FUN = cumsum) - growth_delay * 24 # cumulate hours but take away growth delay
  plant.pres[plant.pres > 0] <- 1 # converting to presence/absence
  plant.pres[plant.pres < 0] <- 0 # converting to presence/absence

  # compute degree day growth
  grow <- plant.pres
  grow[meanpot > wilting_point] <- 1 # growth above wilting point
  grow[meanpot <= wilting_point] <- 0 # no growth below wilting point
  plant_temp[plant_temp > 45] <- temp_thresh # cap plant temperature at 45 deg C
  degdays <- (plant_temp - temp_thresh) / 24 # degree days (divide by 24 because hourly)
  degdays[degdays < 0] <- 1e-10 # make minimum degree days not zero so that cumsum operation in two lines doesn't reset plant growth when thermally unsuitable
  degdays <- degdays * plant.pres # knock out times when below permanent wilting point
  plant.growth <- ave(degdays, cumsum(degdays == 0), FUN = cumsum) # accumulate degree days

  # compute % plant water content
  pct.water <- meanpot
  pct.water[meanpot > wilting_point] <- foodwater # assume plants start wilting at about 2 bar, but above this they are at max water content
  pct.water[meanpot < permanent_wilting_point] <- 0 # assume plants start wilting at about 2 bar, but above this they are at max water content
  pct.water[pct.water < 0]<-abs((abs(pct.water[pct.water < 0]) + permanent_wilting_point)) / abs(permanent_wilting_point - wilting_point) * foodwater # scale % water content from max specified value to zero
  pct.water <- pct.water * plant.pres

  plantgro<-cbind(meanmoist, meanpot, pct.water, plant.pres, plant.growth)

  return(plantgro)
}
