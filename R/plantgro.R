#' plantgro
#'
#' Function to compute plant water content given threshold values of soil water potential at which wilting point and permanent wilting point occurs
#' @param soilpot = soilpot, vector of soil water potential to use for calculations
#' @param soilmoist = soilmoist, vector of soil moisture to use for calculations
#' @param root_shallow = 4shallowest soil node to use when getting average soil moisture/water potential for calculating plant moisture/growth
#' @param root_deep = 4, deepest soil node to use when getting average soil moisture/water potential for calculating plant moisture/growth
#' @param growth_delay = 0, time required for plants to recover after hitting permanent wilting point, days
#' @param wilting_thresh = -200, soil water potential at wilting point, J/kg
#' @param permanent_wilting_point = -1500, soil water potential at permanent wilting point, J/kg
#' @param FoodWater = 82, Maximum water conten of plant, \%
#' @export
plantgro<-function(soilpot=soilpot,soilmoist=soilmoist, root_shallow=4, root_deep=8, growth_delay=1, wilting_thresh=-200, permanent_wilting_point=-1500, FoodWater=82){

  soilmoist<-subset(soilmoist,soilmoist[,"TIME"]==720) # just use midday value
  soilpot<-subset(soilpot,soilpot[,"TIME"]==720) # just use midday
  if(root_shallow==root_deep){
   meanpot<-as.data.frame(soilpot)[,((root_shallow+2):(root_deep+2))]
  }else{
   meanpot<-as.data.frame(soilpot)[,((root_shallow+2):(root_deep+2))] # get range of soil water potential depths to take mean of
   meanpot<-apply(meanpot, 1, mean) # get average soil water potential across chosen depth range
  }
  grow<-meanpot
  grow[grow>permanent_wilting_point]<-1 # find times above the PWP (growth possible)
  grow[grow<=permanent_wilting_point]<-0 # find times when below the PWP (plant dead)
  counter<-0
  grow2<-grow*0 # create empty vector
  for(j in 1:length(grow)){ # accumulate runs of days above the permanent wilting point (growth possible)
    if(j==1){ # check, if first hour, whether we're starting with a growth day or not
      if(grow[j]==1){
        counter<-counter+1 # if so, increment the counter
      }
      grow2[j]<-counter
    }else{ # otherwised, only increment the counter if growth on both present hour and the hour just past
      if(grow[j-1]>0 & grow[j]==1){
        counter<-counter+1
      }else{
        counter<-0
      }
      grow2[j]<-counter
    }
  }
  grow3<-grow2
  grow3[grow3<growth_delay*24]<-0 # apply growth delay specified by the user for time required for plants to come back after PWP has been hit
  grow3[grow3>0]<-1 # make vector of 0 and 1 where 1 means plants could have come back from drought

  if(root_shallow==root_deep){
   meanmoist<-as.data.frame(soilmoist)[,((root_shallow+2):(root_deep+2))]  # get range of soil water moisture depths to take mean of
  }else{
   meanmoist<-as.data.frame(soilmoist)[,((root_shallow+2):(root_deep+2))]  # get range of soil water moisture depths to take mean of
   meanmoist<-apply(meanmoist, 1, mean)
  }
  mean.moist.pot<-as.data.frame(cbind(meanpot,meanmoist))
  colnames(mean.moist.pot)<-c('pot','moist')
  mean.moist.pot$pot[mean.moist.pot$pot>wilting_thresh]<-FoodWater # assume plants start wilting at about 2 bar, but above this they are at max water content
  mean.moist.pot$moist<-mean.moist.pot$moist*100 # convert to percent
  potmult<-mean.moist.pot$pot
  potmult[potmult!=FoodWater]<-0
  potmult[potmult!=0]<-1
  wilting<-subset(mean.moist.pot,pot==FoodWater) # find soil moisture range corresponding to values above the wilting point
  wilting<-min(wilting$moist) # get the min soil moisture at which plants aren't wilting
  pct.water<-mean.moist.pot$moist
  pct.water[pct.water>wilting]<-FoodWater # now have vector of either max plant water content or soil moisture content - need to convert the latter into a smooth decline to zero from max value
  minmoist<-0
  pct.water[pct.water<FoodWater]<-(pct.water[pct.water<FoodWater]-minmoist)/(wilting-minmoist)*FoodWater # for just the values less than max water content, make them equal to the ratio of moisture level and wilting point moisture and then multiplied by food water content so have moisture content ranging from max level down to min moisture possible
  pct.water<-pct.water/100*grow3 # now convert to proportion and cut out times below PWP (including regrowth penalty)

  plantmoist<-mean.moist.pot$moist
  moist.index<-plantmoist/max(plantmoist)*11 # put in units scaling from 0-11
  # next four lines spread the values out more evenly over 11 categories
  minval<-min(moist.index[plantmoist!=0])
  moist.index[plantmoist==0]<-minval
  moist.index<-moist.index-minval
  moist.index<-moist.index/max(moist.index)*11 # put in units scaling from 0-11

  plant.pres<-pct.water
  plant.pres[plant.pres>.6]<-1
  plant.pres[plant.pres<1]<-0

  plantgro<-cbind(plantmoist,moist.index,pct.water,plant.pres)

  return(plantgro)
}
