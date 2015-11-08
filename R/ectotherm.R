#' Ectotherm model.
#'
#' An implementation of the Niche Mapper ectotherm model that computes body temperature, water loss,
#' activity and microhabitat selection. It optionally runs the Dynamic Energy Budget (DEB) model for
#' computing mass budgets (inc. water budgets) and growth, development, reproduction trajectories
#' as constrained by food, activity and temperature (see Details). When not running the DEB model
#' a user-specified mass is used as well as a allometric (mass and body temperature) function to
#' compute metabolic rate.
#'
#' The microclimate model, e.g. micro_global(), must be run prior to running the ectotherm model
#' @param amass Mass of animal (g), note this model is 'steady state' so no lags in heating/cooling due to mass
#' @param lometry Organism shape, 0-4 (see details)
#' @param ABSMAX Maximum solar absorptivity, decimal percent
#' @param ABSMIN Maximum solar absorptivity, decimal percent
#' @param TMAXPR Voluntary thermal maximum, degrees C (upper body temperature for foraging and also affects burrow depth selection)
#' @param TMINPR Voluntary thermal minimum, degrees C (lower body temperature for foraging)
#' @param TBASK Minimum basking temperature, degrees C
#' @param TEMERGE Temperature at which animal will move to a basking site, degrees C
#' @param ctmax Critical thermal maximum, degrees C (affects burrow selection)
#' @param ctmin Critical thermal minimum, degrees C (affects burrow selection)
#' @param tpref Preferred body temperature, degrees C
#' @param dayact Diurnal activity allowed (1) or not (0)?
#' @param nocturn Nocturnal activity allowed (1) or not (0)?
#' @param crepus Crepuscular activity allowed (1) or not (0)?
#' @param CkGrShad shade seeking allowed (1) or not (0)?
#' @param burrow Shelter in burrow allowed (1) or not (0)?
#' @param climb climbing to seek cooler habitats allowed (1) or not (0)?
#' @param shdburrow choose if the animal's retreat is in the shade (1) or in the open (0)
#' @param mindepth Minimum depth (soil node #) to which animal can retreat if burrowing
#' @param maxdepth Maximum depth (soil node #) to which animal can retreat if burrowing
#' @param MR_1 Metabolic rate parameter MR=MR_1*M^MR_2*10^(MR_3*Tb) based on Eq. 2 from Andrews & Pough 1985. Physiol. Zool. 58:214-231
#' @param MR_2 Metabolic rate parameter
#' @param MR_3 Metabolic rate parameter
#' @param skinwet Percentage of surface area acting as a free-water exchanger, for computing cutaneous water loss
#' @param extref Percent oxygen extraction efficiency, for respiratory water loss
#' @param DELTAR Temperature difference (deg C) between expired and inspired air, , for respiratory water loss
#' @usage ectotherm(amass, lometry, ABSMAX, ABSMIN, TMAXPR, TMINPR, TBASK, TEMERGE, ctmax, ctmin,
#'  tpref, dayact, nocturn, crepus, CkGrShad, burrow, climb, shdburrow, mindepth, maxdepth,
#'  MR_1, MR_2, MR_3, ...)
#' @examples
#' # run the microclimate model
#' micro<-micro_global(loc="Paluma, Queensland")
#'
#' # run the ectotherm model
#' ecto<-ectotherm(TMAXPR=35,TMINPR=30,TPREF=33,TBASK=20,TEMERGE=10)
#'
#' # retrieve output
#' metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
#' shadmet<-as.data.frame(micro$shadmet) # above ground microclimatic conditions, max shade
#' soil<-as.data.frame(micro$soil) # soil temperatures, minimum shade
#' shadsoil<-as.data.frame(micro$shadsoil) # soil temperatures, maximum shade
#' environ<-as.data.frame(ecto$environ) # activity, Tb and environment
#' enbal<-as.data.frame(ecto$enbal) # energy balance values
#' masbal<-as.data.frame(ecto$masbal) # mass balance value (note most missing if DEB model not running)
#'
#' # append dates
#' days<-rep(seq(1,12),24)
#'days<-days[order(days)]
#'dates<-days+metout$TIME/60/24-1 # dates for hourly output
#'dates2<-seq(1,12,1) # dates for daily output
#'metout<-cbind(dates,metout)
#'soil<-cbind(dates,soil)
#'shadmet<-cbind(dates,shadmet)
#'shadsoil<-cbind(dates,shadsoil)
#'environ<-cbind(dates,environ)
#'masbal<-cbind(dates,masbal)
#'enbal<-cbind(dates,enbal)
#'
#'############### plot results ######################
#'
#'# Hourly Tb (black), activity (orange, 5=bask, 10=forage), depth (brown, m) and shade (green, %/10)
#'with(environ, plot(TC~dates,ylab="Tb, depth, activity and shade", xlab="month of year",ylim=c(-20,70),type = "l"))
#'with(environ, points(ACT*5~dates,type = "l",col="orange"))
#'with(environ, points(SHADE/10~dates,type = "l",col="green"))
#'with(environ, points(DEP/10~dates,type = "l",col="brown"))
#'#with(metout, points(TAREF~dates,type = "l",col="light blue"))
#'abline(ecto$TMAXPR,0,lty=2,col='red')
#'abline(ecto$TMINPR,0,lty=2,col='blue')
#'
#'# seasonal activity plot (dark blue = night, light blue = basking, orange = foraging)
#'forage<-subset(environ,ACT==2)
#'bask<-subset(environ,ACT==1)
#'night<-subset(metout,ZEN==90)
#'day<-subset(metout,ZEN!=90)
#'with(night,plot(TIME/60~JULDAY,ylab="Hour of Day",xlab="Day of Year",pch=15,cex=2,col='dark blue')) # nighttime hours
#'with(forage,points((TIME-1)~JULDAY,pch=15,cex=2,col='orange')) # foraging Tbs
#'with(bask,points((TIME-1)~JULDAY,pch=15,cex=2,col='light blue')) # basking Tbs
#' @export
ectotherm<-function(amass=5,lometry=3,ABSMAX=0.85,ABSMIN=0.85,
TMAXPR=35,TMINPR=25,TBASK=20,TEMERGE=10,ctmax=40,ctmin=5,TPREF=30,
dayact=1,nocturn=0,crepus=0,CkGrShad=1,burrow=1,climb=0,shdburrow=0,mindepth=2,maxdepth=10,
MR_1=0.013,MR_2=0.8,MR_3=0.038,skinwet=0.2,extref=20.,DELTAR=0.1,
microin="none",write_input=0,enberr=0.0002,minshade=0.,maxshade=micro$MAXSHADES[1],
timeinterval=micro$timeinterval,nyears=length(RAINFALL)/timeinterval,
live=1,ctminthresh=12,ctkill=0,FLTYPE=0.0,SUBTK=2.79,soilnode=4.,REFL=micro$REFL,rinsul=0.,
ptcond=0.25,customallom=c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743),
shape_a=1.,shape_b=3,shape_c=0.6666666667,Flshcond=0.5,Spheat=4185,Andens=1000,
EMISAN=0.95,FATOSK=0.4,FATOSB=0.4,fosorial=0,rainact=0,actrainthresh=0.1,
wings=0,rho1_3=0.2,trans1=0.00,aref=0.26,bref=2.04,cref=1.47,phi=179.,phimax= phi,phimin= phi,
flyer=0,flyspeed=5,flymetab=0.1035,
container=0,conth=10,contw=100.,contype=1,rainmult=1,continit=0,conthole=0,contonly=1,contwet=80,
wetmod=0,soilmoisture1=0,breedactthresh=1,
PFEWAT=73,PTUREA=0,FoodWater=82,minwater=15,raindrink=0.,gutfill=75.,
thermal_stages=matrix(data = c(rep(ctmin,8),rep(ctmax,8),rep(TMINPR,8),rep(TMAXPR,8),rep(TBASK,8),rep(TPREF,8)), nrow = 8, ncol = 6),
behav_stages=matrix(data = c(rep(dayact,8),rep(nocturn,8),rep(crepus,8),rep(burrow,8),rep(shdburrow,8),rep(mindepth,8),rep(maxdepth,8),rep(CkGrShad,8),rep(climb,8),rep(fosorial,8),rep(rainact,8),rep(actrainthresh,8),rep(breedactthresh,8),rep(flyer,8)), nrow = 8, ncol = 14),
water_stages=matrix(data = c(rep(skinwet,8),rep(extref,8),rep(PFEWAT,8),rep(PTUREA,8),rep(FoodWater,8),rep(minwater,8),rep(raindrink,8),rep(gutfill,8)), nrow = 8, ncol = 8),
DEB=0,fract=1,f=1.,MsM=186.03*6.,z=7.174*fract,delta= 0.217,kappa_X=0.85,v_dotref=0.05591/24.,
kappa=0.8501,p_Mref=45.14/24.,E_G=7189,k_R=0.95,k_J=0.00628/24.,E_Hb=6.533e+04*fract^3,E_Hj=E_Hb*fract^3,
E_Hp=1.375e+05*fract^3,h_aref=3.61e-13/(24.^2),s_G=0.01,E_Egg=1.04e+06*fract^3,E_m=(p_Mref*z/kappa)/v_dotref,
p_Xm=13290,K=1,X=10,metab_mode=2,stages=7,y_EV_l=0.95,S_instar=c(2.660,2.310,1.916,0),
s_j=0.999,T_REF=20,TA=7130,TAL=5.305e+04,TAH=9.076e+04,TL=288.,TH=315.,
arrhenius=matrix(data = matrix(data = c(rep(TA,8),rep(TAL,8),rep(TAH,8),rep(TL,8),rep(TH,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5),
andens_deb=1.,d_V=0.3,d_E=0.3,eggdryfrac=0.3,mu_X=525000,mu_E=585000,mu_V=500000,mu_P=480000,
kappa_X_P=0.1,nX=c(1,1.8,0.5,.15),nE=c(1,1.8,0.5,.15),nV=c(1,1.8,0.5,.15),nP=c(1,1.8,0.5,.15),
N_waste=c(1,4/5,3/5,4/5),clutchsize=2.,clutch_ab=c(0,0),viviparous=0,minclutch=0,batch=1,breedrainthresh=0,
photostart=3,photofinish=1,daylengthstart=12.5,daylengthfinish=13.,photodirs = 1,photodirf = 0,
startday=1,breedtempthresh=200,breedtempcum=24*7,reset=0,frogbreed=0,frogstage=0,aestivate=0,depress=0.3,
v_init=(7.063^3)*fract^3*0.85,E_init=E_m,E_H_init=E_Hp+1,stage=3,ma=1e-4,mi=0,mh=0.5,wilting=1,ystrt=0,grasshade=0,
wetlandTemps=matrix(data = 0., nrow = 24*dim, ncol = 1),wetlandDepths=matrix(data = 0., nrow = 24*dim, ncol = 1),
DEP=micro$DEP,ectoin=rbind(as.numeric(micro$ALTT),as.numeric(micro$REFL)[1],micro$longlat[1],micro$longlat[2],0,0,1990,1990),RAINFALL=micro$RAINFALL,
metout=micro$metout,shadmet=micro$shadmet,soil=micro$soil,shadsoil=micro$shadsoil,soilmoist=micro$soilmoist,
shadmoist=micro$shadmoist,humid=micro$humid,shadhumid=micro$shadhumid,soilpot=micro$soilpot,
shadpot=micro$shadpot,MAXSHADES=micro$MAXSHADES){

# amass=5
# lometry=3
# ABSMAX=0.85
# ABSMIN=0.85
# TMAXPR=35
# TMINPR=25
# TBASK=20
# TEMERGE=10
# ctmax=40
# ctmin=5
# TPREF=30
#
# dayact=1
# nocturn=0
# crepus=0
# CkGrShad=1
# burrow=1
# climb=0
# shdburrow=0
# mindepth=2
# maxdepth=10
#
# MR_1=0.013
# MR_2=0.8
# MR_3=0.038
#
# microin="none"
# write_input=0
# timeinterval=micro$timeinterval
# nyears=length(RAINFALL)/timeinterval
# enberr=0.0002
# minshade=0.
# maxshade=micro$MAXSHADES[1]
#
# live=1
# ctminthresh=12
# ctkill=0
#
#
# grasshade=0
# FLTYPE=0.0
# SUBTK=2.79
# soilnode=4.
#
# REFL=micro$REFL
# rinsul=0.
#
# ptcond=0.25
# customallom=c(10.4713,.688,0.425,0.85,3.798,.683,0.694,.743)
# shape_a=1.
# shape_b=3
# shape_c=0.6666666667
# Flshcond=0.5
# Spheat=4185
# Andens=1000
#
# EMISAN=0.95
# ptcond=0.25
# FATOSK=0.4
# FATOSB=0.4
#
# wings=0
# rho1_3=0.2
# trans1=0.00
# aref=0.26
# bref=2.04
# cref=1.47
# phi=179.
# phimax= phi
# phimin= phi
# fosorial=0
# rainact=0
# actrainthresh=0.1
#
#
# flyer=0
# flyspeed=5
# flymetab=0.1035
#
# container=0
# conth=10
# contw=100.
# contype=1
# rainmult=1
# continit=0
# conthole= 0
# contonly=1
# contwet=80
# wetmod=0
# soilmoisture1=0
# breedactthresh=1
#
# skinwet=0.229
# extref=20.
# PFEWAT=73
# PTUREA=0
# FoodWater=82
# minwater=15
# raindrink=0.
# gutfill=75.
# DELTAR=0.1
#
# thermal_stages=matrix(data = c(rep(ctmin,8),rep(ctmax,8),rep(TMINPR,8),rep(TMAXPR,8),rep(TBASK,8),rep(TPREF,8)), nrow = 8, ncol = 6)
# behav_stages=matrix(data = c(rep(dayact,8),rep(nocturn,8),rep(crepus,8),rep(burrow,8),rep(shdburrow,8),rep(mindepth,8),rep(maxdepth,8),rep(CkGrShad,8),rep(climb,8),rep(fosorial,8),rep(rainact,8),rep(actrainthresh,8),rep(breedactthresh,8),rep(flyer,8)), nrow = 8, ncol = 14)
# water_stages=matrix(data = c(rep(skinwet,8),rep(extref,8),rep(PFEWAT,8),rep(PTUREA,8),rep(FoodWater,8),rep(minwater,8),rep(raindrink,8),rep(gutfill,8)), nrow = 8, ncol = 8)
#
# DEB=0
# fract=1
# f=1.
# MsM=186.03*6.
# z=7.174*fract
# delta= 0.217
# kappa_X=0.85
# v_dotref=0.05591/24.
# kappa=0.8501
# p_Mref=45.14/24.
# E_G=7189
# k_R=0.95
# k_J=0.00628/24.
# E_Hb=6.533e+04*fract^3
# E_Hj=E_Hb*fract^3
# E_Hp=1.375e+05*fract^3
# h_aref=3.61e-13/(24.^2)
# s_G=0.01
# E_Egg=1.04e+06*fract^3
# E_m=(p_Mref*z/kappa)/v_dotref
# p_Xm=13290
# K=1
# X=10
# metab_mode=2
# stages=7
# y_EV_l=0.95
# S_instar=c(2.660,2.310,1.916,0)
# s_j=0.999
# T_REF=20
# TA=7130
# TAL=5.305e+04
# TAH=9.076e+04
# TL=288.
# TH=315.
# arrhenius=matrix(data = matrix(data = c(rep(TA,8),rep(TAL,8),rep(TAH,8),rep(TL,8),rep(TH,8)), nrow = 8, ncol = 5), nrow = 8, ncol = 5)
# andens_deb=1.
# d_V=0.3
# d_E=0.3
# eggdryfrac=0.3
# mu_X=525000
# mu_E=585000
# mu_V=500000
# mu_P=480000
# kappa_X_P=0.1
# nX=c(1,1.8,0.5,.15)
# nE=c(1,1.8,0.5,.15)
# nV=c(1,1.8,0.5,.15)
# nP=c(1,1.8,0.5,.15)
# N_waste=c(1,4/5,3/5,4/5)
# clutchsize=2.
# clutch_ab=c(0,0)
# viviparous=0
# minclutch=0
# batch=1
# breedrainthresh=0
# photostart= 3
# photofinish= 1
# daylengthstart= 12.5
# daylengthfinish= 13.
# photodirs = 1
# photodirf = 0
# startday=1
# breedtempthresh=200
# breedtempcum=24*7
# reset=0
# frogbreed=0
# frogstage=0
# aestivate=0
# depress=0.3
# v_init=(7.063^3)*fract^3*0.85
# E_init=E_m
# E_H_init=E_Hp+1
# stage=3
# ma=1e-4
# mi=0
# mh=0.5
# wilting=1
# ystrt=0

  if(lometry==3){
    shape_a<-1.
    shape_b<-1.
    shape_c<-4.
  }
  if(lometry==4){
    shape_a<-1.
    shape_b<-1.
    shape_c<-0.5
  }

  #turn on container model if aquatic egg/larval phase
  if(frogbreed==1 | frogbreed==2){
    container<-1
  }
  if(frogbreed==3){
    container<-0
  }

  # container/pond initial conditons
  contlast<-0.
  templast<-7.

  iyear<-0 #initializing year counter
  countday<-1 #initializing day counter

  if(microin!="none"){
  cat('reading microclimate input \n')
  RAINFALL<-as.matrix(read.csv(file=paste(microin,'rainfall.csv',sep=""),sep=","))[,2]
  dim=length(RAINFALL)
  metout<-read.csv(file=paste(microin,'metout.csv',sep=""),sep=",")[,-1]
  shadmet<-read.csv(file=paste(microin,'shadmet.csv',sep=""),sep=",")[,-1]
  soil<-read.csv(file=paste(microin,'soil.csv',sep=""),sep=",")[,-1]
  shadsoil<-read.csv(file=paste(microin,'shadsoil.csv',sep=""),sep=",")[,-1]
  if(file.exists(paste(microin,'wetlandTemps.csv',sep=""))){
    wetlandTemps<-read.csv(file=paste(microin,'wetlandTemps.csv',sep=""),sep=",")[,-1]
    wetlandDepths<-read.csv(file=paste(microin,'wetlandDepths.csv',sep=""),sep=",")[,-1]
  }else{
    wetlandTemps=matrix(data = 0., nrow = 24*dim, ncol = 1)
    wetlandDepths=matrix(data = 0., nrow = 24*dim, ncol = 1)
  }
    if(file.exists(paste(microin,'soilpot.csv',sep=""))){
      soilpot<-read.csv(file=paste(microin,'soilpot.csv',sep=""),sep=",")[,-1]
      soilmoist<-read.csv(file=paste(microin,'soilmoist.csv',sep=""),sep=",")[,-1]
      shadpot<-read.csv(file=paste(microin,'shadpot.csv',sep=""),sep=",")[,-1]
      shadmoist<-read.csv(file=paste(microin,'shadmoist.csv',sep=""),sep=",")[,-1]
      humid<-read.csv(file=paste(microin,'humid.csv',sep=""),sep=",")[,-1]
      shadhumid<-read.csv(file=paste(microin,'shadhumid.csv',sep=""),sep=",")[,-1]
    }else{
      soilpot<-soil
      soilmoist<-soil
      shadpot<-soil
      shadmoist<-soil
      humid<-soil
      shadhumid<-soil
      soilpot[,3:12]<-0
      soilmoist[,3:12]<-0.5
      shadpot[,3:12]<-0
      shadmoist[,3:12]<-0.5
      humid[,3:12]<-0.99
      shadhumid[,3:12]<-0.99
    }
  metout<-as.matrix(metout)
  shadmet<-as.matrix(shadmet)
  shadsoil<-as.matrix(shadsoil)
  soil<-as.matrix(soil)
  soilmoist<-as.matrix(soilmoist)
  shadmoist<-as.matrix(shadmoist)
  soilpot<-as.matrix(soilpot)
  shadpot<-as.matrix(shadpot)
  humid<-as.matrix(humid)
  shadhumid<-as.matrix(shadhumid)
  ectoin<-read.csv(file=paste(microin,'ectoin.csv',sep=""),sep=",")[,-1]
  DEP<-as.matrix(read.csv(file=paste(microin,'DEP.csv',sep=""),sep=","))[,2]
  MAXSHADES<-as.matrix(read.csv(file=paste(microin,'MAXSHADES.csv',sep=""),sep=","))[,2]
  metout2=matrix(data = 0., nrow = 24*dim, ncol = 18)
  soil2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  shadmet2=matrix(data = 0., nrow = 24*dim, ncol = 18)
  shadsoil2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  soilmoist2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  shadmoist2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  soilpot2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  shadpot2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  humid2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  shadhumid2=matrix(data = 0., nrow = 24*dim, ncol = 12)
  wetlandTemps2=matrix(data = 0., nrow = 24*dim, ncol = 1)
  wetlandDepths2=matrix(data = 0., nrow = 24*dim, ncol = 1)
  metout2[1:nrow(metout),]<-metout
  shadmet2[1:nrow(metout),]<-shadmet
  soil2[1:nrow(metout),]<-soil
  shadsoil2[1:nrow(metout),]<-shadsoil
  soilmoist2[1:nrow(metout),]<-soilmoist
  shadmoist2[1:nrow(metout),]<-shadmoist
  soilpot2[1:nrow(metout),]<-soilpot
  shadpot2[1:nrow(metout),]<-shadpot
  humid2[1:nrow(metout),]<-humid
  shadhumid2[1:nrow(metout),]<-shadhumid
  wetlandTemps2[1:nrow(metout)]<-wetlandTemps
  wetlandDepths2[1:nrow(metout)]<-wetlandDepths
  metout<-metout2
  shadmet<-shadmet2
  soil<-soil2
  shadsoil<-shadsoil2
  soilmoist<-soilmoist2
  shadmoist<-shadmoist2
  soilpot<-soilpot2
  shadpot<-shadpot2
  humid<-humid2
  shadhumid<-shadhumid2
  wetlandTemps<-wetlandTemps2
  wetlandDepths<-wetlandDepths2
  metout.names<-c("JULDAY","TIME","TALOC","TAREF","RHLOC","RH","VLOC","VREF","SOILMOIST3","POOLDEP","TDEEP","ZEN","SOLR","TSKYC","DEW","FROST","SNOWFALL","SNOWDEP")
  colnames(metout)<-metout.names
  colnames(shadmet)<-metout.names
  soil.names<-c("JULDAY","TIME",paste("D",DEP,"cm", sep = ""))
  colnames(soil)<-soil.names
  colnames(shadsoil)<-soil.names
  moist.names<-c("JULDAY","TIME",paste("WC",DEP,"cm", sep = ""))
  colnames(soilmoist)<-moist.names
  colnames(shadmoist)<-moist.names
  pot.names<-c("JULDAY","TIME",paste("PT",DEP,"cm", sep = ""))
  colnames(soilpot)<-pot.names
  colnames(shadpot)<-pot.names
  hum.names<-c("JULDAY","TIME",paste("RH",DEP,"cm", sep = ""))
  colnames(humid)<-hum.names
  colnames(shadhumid)<-hum.names
  }else{
    dim=length(RAINFALL)
  }


  if(soilmoisture1==1){
  soilpotb<-soilpot
  soilmoistb<-soilmoist
  }
  # habitat
  ALT<-ectoin[1] # altitude (m)
  OBJDIS<-1.0 # distance from object (e.g. bush)
  OBJL<-0.0001
  PCTDIF<-0.1 # percent of sunlight that is diffuse (decimal %)
  EMISSK<-1.0 # emissivity of the sky (decimal %)
  EMISSB<-1.0 # emissivity of the substrate (decimal %)
  ABSSB<-1-ectoin[2] # solar absorbtivity of the substrate (decimal %)
  shade<-minshade # shade (%)

  # animal properties
  AMASS<-amass/1000 # animal mass (kg)
  absan<-ABSMAX # animal solar absorbtivity
  RQ<-0.8 # respiratory quotient

  FATOBJ<-0.
  #  if(container==1){
  #    live<-0}else{live<-1
  #  }
  #live<-1
  TIMBAS<-1.
  #  if(container==1){
  #    SKINW<-100.}else{
  SKINW<-skinwet
  #    }
  skint<-0.
  O2gas<-20.95
  CO2gas<-0.03
  N2gas<-79.02
  gas<-c(O2gas,CO2gas,N2gas)
  #  if(container==1){
  #    transt<-1
  #  }else{
  transt<-0
  #  }
  tranin<-1
  tcinit<-metout[1,"TALOC"]

  ACTLVL<-1
  nodnum<-10
  spec<-0. # spectacle covering eye surface? (adds to water loss for lizard/frog/turtle geometry)
  xbas<-1.
  nofood<-0
  tdigpr<-TPREF
  o2max<-extref
  #  if(container==1){
  #  maxshd<-1.
  #  minshd<-0.
  #  }else{
  maxshd<-maxshade
  minshd<-minshade
  #  }
  behav=c(dayact,nocturn,crepus,rainact,burrow,CkGrShad,climb,fosorial,nofood)
  julday<-1

  # DEB model initial conditions
  V_init_baby<-3e-9
  E_init_baby<-E_Egg/V_init_baby
  E_baby_init<-E_init_baby
  V_baby_init<-V_init_baby
  ms_init<-0.
  cumrepro_init<-0.
  q_init<-0.
  hs_init<-0.
  cumbatch_init<-0.
  pregnant<-0
  E_m<-(p_Mref*z/kappa)/v_dotref

  # conversions from percent to proportion
  PTUREA1<-PTUREA/100
  PFEWAT1<-PFEWAT/100
  FoodWater1<-FoodWater/100
  water_stages[,3]<-water_stages[,3]/100
  water_stages[,4]<-water_stages[,4]/100
  water_stages[,5]<-water_stages[,5]/100
  eggmass<-0 # initial dry mass of an egg (g) - no longer used so delete

  #DEB mass balance calculations
  nO<-cbind(nX,nV,nE,nP) # matrix of composition of organics, i.e. food, structure, reserve and faeces
  CHON<-c(12,1,16,14)
  wO<-CHON%*%nO
  w_V=wO[3]
  M_V<-d_V/w_V
  yEX<-kappa_X*mu_X/mu_E # yield of reserve on food
  yXE<-1/yEX # yield of food on reserve
  yVE<-mu_E*M_V/E_G  # yield of structure on reserve
  yPX<-kappa_X_P*mu_X/mu_P # yield of faeces on food
  yXP<-1/yPX # yield of food on faeces
  yPE<-yPX/yEX # yield of faeces on reserve  0.143382353
  nM<-matrix(c(1,0,2,0,0,2,1,0,0,0,2,0,N_waste),nrow=4)
  N_waste_inv<-c(-1*N_waste[1]/N_waste[4],(-1*N_waste[2])/(2*N_waste[4]),(4*N_waste[1]+N_waste[2]-2*N_waste[3])/(4*N_waste[4]),1/N_waste[4])
  nM_inv<-matrix(c(1,0,-1,0,0,1/2,-1/4,0,0,0,1/2,0,N_waste_inv),nrow=4)
  JM_JO<--1*nM_inv%*%nO
  etaO<-matrix(c(yXE/mu_E*-1,0,1/mu_E,yPE/mu_E,0,0,-1/mu_E,0,0,yVE/mu_E,-1/mu_E,0),nrow=4)
  w_N<-CHON%*%N_waste

    #wilting<-ectoin[6] # %vol, water content at 15ba = 1500kPa (wiki for thresholds)

  lat<-ectoin[4]
  if(soilmoisture1==1){
  humid[,3:9]<-metout[,5]/100 # assume ambient humidity down to 30cm
  shadhumid[,3:9]<-shadmet[,5]/100 # assume ambient humidity down to 30cm
  humid[,7:12]<-0.8 # assume higher humidity in burrow, 60cm and lower
  shadhumid[,7:12]<-0.8 # assume higher humidity in burrow, 60cm and lower


  grassgrowths<-as.data.frame(soilpotb)
  soilmoist2b<-as.data.frame(soilmoistb)
  soilmoist2b<-subset(soilmoist2b,soilmoist2b$TIME==720)
  grassgrowths<-subset(grassgrowths,soilmoist2b$TIME==720)
  grassgrowths<-grassgrowths$PT5cm # assume plant growth driven by 5cm depth

    grow<-grassgrowths
    grow[grow>-1500]<-1 # find times when below permanent wilting point
    grow[grow<=-1500]<-0
    counter<-0
    grow2<-grow*0
      for(j in 1:length(grow)){
        if(j==1){
            if(grow[j]==1){
            counter<-counter+1
            }
          grow2[j]<-counter
        }else{
          if(grow[j-1]>0 & grow[j]==1){
            counter<-counter+1
          }else{
            counter<-0
          }
          grow2[j]<-counter
        }
      }
     grow3<-grow2
     grow3[grow3<7]<-0 # use one week in a row as time required for plats to come back after PWP has been hit
     grow3[grow3>0]<-1 # make vector of 0 and 1 where 1 means plants could have come back from drought

  soilmoist2b<-soilmoist2b$WC5cm
  grassgrowths<-as.data.frame(cbind(grassgrowths,soilmoist2b))
  colnames(grassgrowths)<-c('pot','moist')
  grassgrowths$pot[grassgrowths$pot>-200]<-FoodWater # assume plants start wilting at about 2 bar, but above this they are at max water content
  grassgrowths$moist<-grassgrowths$moist*100 # convert to percent
  potmult<-grassgrowths$pot
  potmult[potmult!=82]<-0
  potmult[potmult!=0]<-1
  wilting<-subset(grassgrowths,grassgrowths$pot==FoodWater) # find soil moisture range corresponding to values above the wilting point
  wilting<-min(wilting$moist) # get the min soil moisture at which plants aren't wilting
  grassgrowths<-grassgrowths$moist
  grassgrowths[grassgrowths>wilting]<-FoodWater # now have vector of either max plant water content or soil moisture content - need to convert the latter into a smooth decline to zero from max value
  #minmoist<-min(grassgrowths[grassgrowths<FoodWater])
  minmoist<-0
  grassgrowths[grassgrowths<FoodWater]<-(grassgrowths[grassgrowths<FoodWater]-minmoist)/(wilting-minmoist)*FoodWater # for just the values less than max water content, make them equal to the
  grassgrowths<-grassgrowths/100*grow3
  #grassgrowths<-c(grassgrowths,rep(0,nrow(metout)/24-length(grassgrowths))) # this was to ensure a vector of 20 years x 24 hours, but not needed anymore
  #grassgrowths[grassgrowths<FoodWater/100]<-0
  grasstsdms<-grassgrowths
  #minmoist<-min(grassgrowths)
  #grassgrowths<-(as.numeric(grassgrowths)-minmoist)
  #maxmoist<-max(grassgrowths)
  #grassgrowths[grassgrowths<wilting]<-0
  #grassgrowths<-grassgrowths/maxmoist*X#rep(X,timeinterval*nyears)#
  #grasstsdms<-grassgrowths/maxmoist*X#rep(X,timeinterval*nyears)#
  }else{
    grassgrowths<-rep(FoodWater,nrow(metout))
    grasstsdms<-grassgrowths
  }
  julstart<-metout[1,2]
  tannul<-as.numeric(mean(soil[,12]))
  monthly<-0
  tester<-0
  microyear<-1

  # bucket model for soil moisture
  fieldcap<-ectoin[5]# %vol, water content at 0.1ba = 10kPa
  fieldcap<-30 # field capacity, m3/m3*100
  if(soilmoisture1==1){
    conth<-fieldcap/10 # containter height, cm
    contw<-100
    contype<-1 # is 'containter' sitting on the surface, like a bucket (0) or sunk into the ground like a pond (1)
    rainmult<-0.3 # !!!!!!!!!!!!!!rainfall multiplier to reflect catchment (don't make this zero unless you want a drought!)
    continit<-0 # initial container water level (cm)
    conthole<-0#2.8 # daily loss of height (mm) due to 'hole' in container (e.g. infiltration to soil, drawdown from water tank)
    contwet<- 2 # percent wet value for container
  }

  ectoinput<-as.matrix(c(ALT,FLTYPE,OBJDIS,OBJL,PCTDIF,EMISSK,EMISSB,ABSSB,shade,enberr,AMASS,EMISAN,absan,RQ,rinsul,lometry,live,TIMBAS,Flshcond,Spheat,Andens,ABSMAX,ABSMIN,FATOSK,FATOSB,FATOBJ,TMAXPR,TMINPR,DELTAR,SKINW,spec,xbas,extref,TPREF,ptcond,skint,gas,transt,soilnode,o2max,ACTLVL,tannul,nodnum,tdigpr,maxshd,minshd,ctmax,ctmin,behav,julday,actrainthresh,viviparous,pregnant,conth,contw,contlast,tranin,tcinit,nyears,lat,rainmult,julstart,monthly,customallom,MR_1,MR_2,MR_3,DEB,tester,rho1_3,trans1,aref,bref,cref,phi,wings,phimax,phimin,shape_a,shape_b,shape_c,minwater,microyear,container,flyer,flyspeed,dim,maxdepth,ctminthresh,ctkill,gutfill,mindepth,TBASK,TEMERGE,p_Xm,SUBTK,flymetab,continit,wetmod,contonly,conthole,contype,shdburrow,breedtempthresh,breedtempcum,contwet,fieldcap,wilting,soilmoisture1,grasshade))
  debmod<-c(clutchsize,andens_deb,d_V,eggdryfrac,mu_X,mu_E,mu_V,mu_P,T_REF,z,kappa,kappa_X,p_Mref,v_dotref,E_G,k_R,MsM,delta,h_aref,V_init_baby,E_init_baby,k_J,E_Hb,E_Hj,E_Hp,clutch_ab[2],batch,breedrainthresh,photostart,photofinish,daylengthstart,daylengthfinish,photodirs,photodirf,clutch_ab[1],frogbreed,frogstage,etaO,JM_JO,E_Egg,kappa_X_P,PTUREA1,PFEWAT1,wO,w_N,FoodWater1,f,s_G,K,X,metab_mode,stages,y_EV_l,s_j,startday,raindrink,reset,ma,mi,mh,aestivate,depress,minclutch)
  deblast<-c(iyear,countday,v_init,E_init,ms_init,cumrepro_init,q_init,hs_init,cumbatch_init,V_baby_init,E_baby_init,E_H_init,stage)

  origjulday<-metout[,1]
  if(ystrt>0){
    metout<-rbind(metout[((ystrt)*365*24+1):(dim*24),],metout[1:((ystrt)*365*24),])
    shadmet<-rbind(shadmet[((ystrt)*365*24+1):(dim*24),],shadmet[1:((ystrt)*365*24),])
    soil<-rbind(soil[((ystrt)*365*24+1):(dim*24),],soil[1:((ystrt)*365*24),])
    shadsoil<-rbind(shadsoil[((ystrt)*365*24+1):(dim*24),],shadsoil[1:((ystrt)*365*24),])
    soilmoist<-rbind(soilmoist[((ystrt)*365*24+1):(dim*24),],soilmoist[1:((ystrt)*365*24),])
    shadmoist<-rbind(shadmoist[((ystrt)*365*24+1):(dim*24),],shadmoist[1:((ystrt)*365*24),])
    soilpot<-rbind(soilpot[((ystrt)*365*24+1):(dim*24),],soilpot[1:((ystrt)*365*24),])
    shadpot<-rbind(shadpot[((ystrt)*365*24+1):(dim*24),],shadpot[1:((ystrt)*365*24),])
    humid<-rbind(humid[((ystrt)*365*24+1):(dim*24),],humid[1:((ystrt)*365*24),])
    shadhumid<-rbind(shadhumid[((ystrt)*365*24+1):(dim*24),],shadhumid[1:((ystrt)*365*24),])
    wetlandDepths<-c(wetlandDepths[((ystrt)*365*24+1):(dim*24)],wetlandDepths[1:((ystrt)*365*24)])
    wetlandTemps<-c(wetlandTemps[((ystrt)*365*24+1):(dim*24)],wetlandTemps[1:((ystrt)*365*24)])
    MAXSHADES<-c(MAXSHADES[((ystrt)*365+1):(dim)],MAXSHADES[1:((ystrt)*365)])
    RAINFALL<-c(RAINFALL[((ystrt)*365+1):(dim)],RAINFALL[1:((ystrt)*365)])
    grassgrowths<-c(grassgrowths[((ystrt)*365+1):(dim)],grassgrowths[1:((ystrt)*365)])
    metout[,1]<-origjulday
    shadmet[,1]<-origjulday
    soil[,1]<-origjulday
    shadsoil[,1]<-origjulday
    soilmoist[,1]<-origjulday
    shadmoist[,1]<-origjulday
    soilpot[,1]<-origjulday
    shadpot[,1]<-origjulday
    humid[,1]<-origjulday
    shadhumid[,1]<-origjulday
  }


  # code to determine wet periods for activity in a pond

  if(wetmod==1){
  wet_thresh<-10*24 # threshold pond duration
  wet_depth<-100 # threshold pond depth (mm)
  wet_temp<-28 # threshold exit temp (deg C)
  b<-cbind(as.data.frame(wetlandDepths),as.data.frame(wetlandTemps))
  colnames(b)<-c('depth','temp')
  b$depth[b$temp>wet_temp]<-0
  b<-b$depth
  b[b>=wet_depth]<-1
  b[b!=1]<-0
  bb<-rle(b)
  bb$values[bb$lengths<wet_thresh]<-0
  c<-b*0
  values<-bb$values
  lengths<-bb$lengths
  for(k in 1:length(bb$values)){
    d<-c(rep(values[k],lengths[k]))
    if(k==1){
    e<-d
    }else{
    e<-c(e,d)
    }
  }
  wetlandDepths<-wetlandDepths*e
  }

  if(write_input==1){
      if(dir.exists("ecto csv input")==FALSE){
        dir.create("ecto csv input")
      }
    cat('writing input csv files \n')
    write.csv(ectoinput, file = "ecto csv input/ectoinput.csv")
    write.csv(debmod, file = "ecto csv input/debmod.csv")
    write.csv(deblast, file = "ecto csv input/deblast.csv")
    write.csv(RAINFALL, file = "ecto csv input/rainfall.csv")
    write.csv(DEP, file = "ecto csv input/dep.csv")
    write.csv(grassgrowths, file = "ecto csv input/grassgrowths.csv")
    write.csv(grasstsdms, file = "ecto csv input/grasstsdms.csv")
    write.csv(wetlandTemps, file = "ecto csv input/wetlandTemps.csv")
    write.csv(wetlandDepths, file = "ecto csv input/wetlandDepths.csv")
    write.csv(arrhenius, file = "ecto csv input/arrhenius.csv")
    write.csv(thermal_stages, file = "ecto csv input/thermal_stages.csv")
    write.csv(behav_stages, file = "ecto csv input/behav_stages.csv")
    write.csv(water_stages, file = "ecto csv input/water_stages.csv")
    write.csv(MAXSHADES, file = "ecto csv input/Maxshades.csv")
    write.csv(S_instar, file = "ecto csv input/S_instar.csv")
    write.table(metout[(seq(1,dim*24)),], file = "ecto csv input/metout.csv",sep=",",row.names=FALSE)
    write.table(shadmet[(seq(1,dim*24)),], file = "ecto csv input/shadmet.csv",sep=",",row.names=FALSE)
    write.table(soil[(seq(1,dim*24)),], file = "ecto csv input/soil.csv",sep=",",row.names=FALSE)
    write.table(shadsoil[(seq(1,dim*24)),], file = "ecto csv input/shadsoil.csv",sep=",",row.names=FALSE)
    write.table(soilmoist[(seq(1,dim*24)),], file = "ecto csv input/soilmoist.csv",sep=",",row.names=FALSE)
    write.table(shadmoist[(seq(1,dim*24)),], file = "ecto csv input/shadmoist.csv",sep=",",row.names=FALSE)
    write.table(soilpot[(seq(1,dim*24)),], file = "ecto csv input/soilpot.csv",sep=",",row.names=FALSE)
    write.table(shadpot[(seq(1,dim*24)),], file = "ecto csv input/shadpot.csv",sep=",",row.names=FALSE)
    write.table(humid[(seq(1,dim*24)),], file = "ecto csv input/humid.csv",sep=",",row.names=FALSE)
    write.table(shadhumid[(seq(1,dim*24)),], file = "ecto csv input/shadhumid.csv",sep=",",row.names=FALSE)
  }
  ecto<-list(dim=dim,ectoinput=ectoinput,metout=metout,shadmet=shadmet,soil=soil,shadsoil=shadsoil,soilmoist=soilmoist,shadmoist=shadmoist,soilpot=soilpot,shadpot=shadpot,humid=humid,shadhumid=shadhumid,DEP=DEP,RAINFALL=RAINFALL,iyear=iyear,countday=countday,debmod=debmod,deblast=deblast,grassgrowths=grassgrowths,grasstsdms=grasstsdms,wetlandTemps=wetlandTemps,wetlandDepths=wetlandDepths,arrhenius=arrhenius,thermal_stages=thermal_stages,behav_stages=behav_stages,water_stages=water_stages,MAXSHADES=MAXSHADES,S_instar=S_instar)

  cat('running ectotherm model ... \n')

  ptm <- proc.time() # Start timing
  ectout<-ectorun(ecto)
  print(proc.time() - ptm) # Stop the clock


  environ<-ectout$environ[1:(dim*24),]
  enbal<-ectout$enbal[1:(dim*24),]
  masbal<-ectout$masbal[1:(dim*24),]
  debout<-ectout$debout[1:(dim*24),]
  yearout<-ectout$yearout
  yearsout<-ectout$yearsout[1:nyears,]

  if(DEB==0){
    return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,soilpot=soilpot,shadpot=shadpot,humid=humid,shadhumid=shadhumid,RAINFALL=RAINFALL,enbal=enbal,environ=environ,masbal=masbal,yearout=yearout,yearsout=yearsout,grassgrowths=grassgrowths,grasstsdms=grasstsdms,TMAXPR=TMAXPR,TMINPR=TMINPR,ctmax=ctmax,ctmin=ctmin,TBASK=TBASK,TEMERGE=TEMERGE))
  }else{
    return(list(soil=soil,shadsoil=shadsoil,metout=metout,shadmet=shadmet,soilmoist=soilmoist,shadmoist=shadmoist,soilpot=soilpot,shadpot=shadpot,humid=humid,shadhumid=shadhumid,RAINFALL=RAINFALL,enbal=enbal,masbal=masbal,environ=environ,debout=debout,yearout=yearout,yearsout=yearsout,grassgrowths=grassgrowths,grasstsdms=grasstsdms,TMAXPR=TMAXPR,TMINPR=TMINPR,ctmax=ctmax,ctmin=ctmin,TBASK=TBASK,TEMERGE=TEMERGE))
  }

}
