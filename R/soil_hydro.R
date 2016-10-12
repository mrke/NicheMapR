#' Soil hydrological properties calculator
#'
#' A function to compute soil hydrological properties from information on bulk density and soil texture (clay/silt/sand composition) at particular depths, with capacity to spline results to other depths (used in NicheMapR). Calculations are based on equations in Campbell, G. S. 1985. Soil Physics with Basic: Transport Models for Soil-Plant Systems. Elsevier, Amsterdam, and Rab, M. A., S. Chandra, P. D. Fisher, N. J. Robinson, M. Kitching, C. D. Aumann, and M. Imhof. 2011. Modelling and prediction of soil water contents at field capacity and permanent wilting point of dryland cropping soils. Soil Research 49:389-407.
#' @param soilpro Matrix of n x 5 matrix of soil composition with the following columns 1. depth (cm), 2. bulk density (Mg/m3), 3. clay (%), 4. silt (%), 5. clay (%)
#' @return PE air entry water potential (J/kg), Campbell (1985) eq. 5.12, p. 46
#' @return BB Campbell's b parameter, Campbell (1985) eq. 5.11, p. 45
#' @return BD bulk density, Mg/m3
#' @return KS saturated hydraulic conductivity (kg s / m3), Campbell (1985) eq. 6.12, p. 54
#' @return FC Field capacity (m3/m3, %) Based on model 6 in Table 6 of Rab, M. A., S. Chandra, P. D. Fisher, N. J. Robinson, M. Kitching, C. D. Aumann, and M. Imhof. 2011. Modelling and prediction of soil water contents at field capacity and permanent wilting point of dryland cropping soils. Soil Research 49:389-407.
#' @return PWP Permanent Wilting Point (m3/m3, %) Based on model 2 in Table 7 of Rab et al. 2011 (cited above)
#' @usage soil.hydro(soilpro)
#' @export
soil.hydro<-function(soilpro = as.data.frame(soilpro), DEP = soilpro[,1]){

  # field capacity, model 6 from Table 6 of Rab et al. (2011).
  FC<-(7.561+1.176*soilpro$clay-0.009843*soilpro$clay^2+0.2132*soilpro$silt)/100

  # permanent wilting point, model 2, Table 7 of Rab et al. (2011)
  PWP<-(-1.304+1.117*soilpro$clay-0.009309*soilpro$clay^2)/100

  # depths at which soil texture is known
  soil_depths<-soilpro$depth
  nodeout = length(DEP)*2-2

  # find half-way points between given depths
  DEP2<-rep(0,nodeout)
  j<-1
  for(i in 1:length(DEP2)){ # loop through all depths
    if(i%%2==0){ # alternates between nodes to get half way points
      DEP2[i]<-DEP2[i-1]+(DEP[j]-DEP2[i-1])/2
    }else{
      DEP2[i]<-DEP[j]
      j<-j+1
    }
  }
  DEP2<-as.data.frame(floor(DEP2))
  colnames(DEP2)<-"DEPTH"

  # particle diameters from Campbell (1985) p. 10
  dclay<-0.001 #mm
  dsilt<-0.026 #mm
  dsand<-1.05 #mm

  # Campbell (1985) eq. 2.17, p. 9
  a<-(soilpro$clay/100)*log(dclay) + (soilpro$sand/100)*log(dsand) + (soilpro$silt/100)*log(dsilt)

  # Campbell (1985) eq. 2.18, p. 9
  b.1<-(((soilpro$clay/100)*log(dclay)^2+(soilpro$sand/100)*log(dsand)^2+(soilpro$silt/100)*log(dsilt)^2)-a^2)^(1/2)

  # Campbell (1985) eq. 2.15, p. 9
  dg<-exp(a) # geometric mean of particle diameter

  # Campbell (1985) eq. 2.16, p. 9
  sigma_g<-exp(b.1) # geometric standard deviation of particle diameter

  # Reference air entry water potential Campbell (1985) eq. 5.10, p. 45
  PES<-(0.5*dg^(-1/2))*-1 # air entry water potential reference value at bulk density of 1.3 Mg/m3

  # Campbell's b parameter, Campbell (1985) eq. 5.11, p. 45
  b<--2*PES+0.2*sigma_g # slope of ln of water potential against ln of volumetric water content

  # air entry water potential (J/kg), Campbell (1985) eq. 5.12, p. 46
  PE<-PES*(soilpro$blkdens/1.3)^(0.67*b) #

  # saturated hydraulic conductivity (kg s / m3), Campbell (1985) eq. 6.12, p. 54
  KS<-0.004*(1.3/soilpro$blkdens)^(1.3*b)*exp(-6.9*soilpro$clay/100-3.7*soilpro$silt/100)

  # bulk density, Mg/m3
  BD<-soilpro$blkdens

  # spline to new depths if needed
  if(nodeout != length(DEP)){

    KS_spline <-spline(soil_depths,KS,n=201,xmin=0,xmax=200,method='natural')
    KS_spline<-as.data.frame(cbind(KS_spline$x,KS_spline$y))
    colnames(KS_spline)<-c('DEPTH','VALUE')
    KS<-merge(DEP2,KS_spline)
    KS<-c(KS[1,2],KS[,2])
    KS[KS<0.000017]<-0.000017

    PE_spline <-spline(soil_depths,PE,n=201,xmin=0,xmax=200,method='natural')
    PE_spline<-as.data.frame(cbind(PE_spline$x,PE_spline$y))
    colnames(PE_spline)<-c('DEPTH','VALUE')
    PE<-merge(DEP2,PE_spline)
    PE<-c(-1*PE[1,2],-1*PE[,2])

    b_spline <-spline(soil_depths,b,n=201,xmin=0,xmax=200,method='natural')
    b_spline<-as.data.frame(cbind(b_spline$x,b_spline$y))
    colnames(b_spline)<-c('DEPTH','VALUE')
    b<-merge(DEP2,b_spline)
    BB<-c(b[1,2],b[,2])

    BD_spline <-spline(soil_depths,BD,n=201,xmin=0,xmax=200,method='natural')
    BD_spline<-as.data.frame(cbind(BD_spline$x,BD_spline$y))
    colnames(BD_spline)<-c('DEPTH','VALUE')
    BD<-merge(DEP2,BD_spline)
    BD<-c(BD[1,2],BD[,2])

    FC_spline <-spline(soil_depths,FC,n=201,xmin=0,xmax=200,method='natural')
    FC_spline<-as.data.frame(cbind(FC_spline$x,FC_spline$y))
    colnames(FC_spline)<-c('DEPTH','VALUE')
    FC<-merge(DEP2,FC_spline)
    FC<-c(FC[1,2],FC[,2])

    PWP_spline <-spline(soil_depths,PWP,n=201,xmin=0,xmax=200,method='natural')
    PWP_spline<-as.data.frame(cbind(PWP_spline$x,PWP_spline$y))
    colnames(PWP_spline)<-c('DEPTH','VALUE')
    PWP<-merge(DEP2,PWP_spline)
    PWP<-c(PWP[1,2],PWP[,2])
  } # end check if splining

  return<-list(BD = BD, BB = BB, KS = KS, PE = PE, FC = FC, PWP = PWP)
}
