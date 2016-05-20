manmo<-function(Body.Mass = 58, Height = 162, Activity.MR = 121.6, Maximum.SR = 9,
  Temp.Increase = 5, Area.Clothed = 0.40, Albedo.body = 0.30, Tskin = 36, Speed = 0,
  Evap.Eff = 0.85, sigma = 0.0000000567, emiss = 0.95, Tair = 25, Tsky = 25, Tsub = 25,
  RH = 20, WS = 1,  Solar.Rad = 0, Solar.Load = -1,  Zen = 0, diff = 0.15,  refl = 0.15){

#   Body.Mass = 58.00 #kg
#   Height = 162.00 #cm
#   Activity.MR = 121.60 #watts/m2
#   Maximum.SR = 9.00 #g/min.m2
#   Temp.Increase = 5.00
#   Area.Clothed = 0.40 #proportion
#   Albedo.body = 0.30
#   Tskin = 36.00
#   Speed = 0.00
#   Evap.Eff = 0.85
#   sigma = 0.0000000567
#   emiss = 0.95
#   Tair = 25
#   Tsky = 25
#   Tsub = 25
#   RH = 20
#   WS = 1
#   Solar.Rad = 0 #W/m2
#   Solar.Load = -1 #W
#   diff = 0.15 #% diffuse sunlight
#   refl = 0.15 #% reflected sunlight
#   Zen=90

  # get max dimension of problem
  dim = max(length(Body.Mass), length(Height), length(Activity.MR), length(Maximum.SR),
    length(Temp.Increase), length(Area.Clothed), length(Albedo.body), length(Tskin), length(Speed),
    length(Evap.Eff), length(sigma), length(emiss), length(Tair), length(Tsky), length(Tsub), length(RH),
    length(WS), length(Solar.Rad), length(Solar.Load), length(Zen))

  environs<-list(Tair,Tsky,Tsub,RH,WS,Solar.Rad,Solar.Load,Zen)
  datatypes=sapply(environs, FUN=class)
  if("RasterLayer" %in% datatypes==TRUE){
  rasterdata=1
  # initialise rasters involving if statements
   Sweat = environs[[which(datatypes=="RasterLayer")[1]]]*0
   Sweat.limited = Sweat
   for(i in 1:length(environs)){
     if(datatypes[i]!="RasterLayer"){
       environs[i]=Sweat+environs[[i]]
     }
   }
   Tair=environs[[1]]
   Tsky=environs[[2]]
   Tsub=environs[[3]]
   RH=environs[[4]]
   WS=environs[[5]]
   Solar.Rad=environs[[6]]
   if(Solar.Load>=0){
   Solar.Load=environs[[7]]
   }
   Zenith=environs[[8]]*pi/180
   hc = WS
  }else{
  rasterdata=0
  # initialise arrays involving if statements
   Sweat = rep(0,dim)
   Sweat.limited = rep(0,dim)
   hc = rep(0,dim)
   Zenith=Zen*pi/180
  }
  sigma.emiss = sigma * emiss

  SA = 0.007184*(Body.Mass^0.425)*(Height^0.725) #m2
  Metab = Activity.MR*SA #watts
  Max.Sweat.Rate = Maximum.SR*SA*60 #g/h
  Max.Evap.Rate = ((Max.Sweat.Rate*2400)/3600)*Evap.Eff #watts

  WS.effective = WS + Speed

  Qnorm <- (Solar.Rad / cos(Zenith))
  Qnorm[Zen>=90]<-0
  Qnorm[Qnorm>1367]<-1367 #making sure that low sun angles don't lead to solar values greater than the solar constant

  if(Solar.Load<0){
    # equation below based on section 9 of Klein Klein, W. H. 1948. Calculation of solar radiation and the
    # solar heat load on man. — Journal of Meteorology 5: 121-130.
    Solar.Load<-(Qnorm*(1-diff)*tan(Zenith)*.51+Solar.Rad*diff*0.85+Solar.Rad*refl*0.85)
  }
  maxload=Qnorm*(1-diff)*.51
  Solar.Load[Solar.Load>maxload]<-maxload[Solar.Load>maxload]
  Total.Solar.Load.Incident.Absorbed = Solar.Load*(1-Albedo.body)
  Total.Heat.Load = Metab + Total.Solar.Load.Incident.Absorbed
  hc[WS>4] = 5.38*WS[WS>4]^0.9
  hc[WS>2] = 5.83*WS[WS>2]^0.805
  hc[WS>0.5] = 6.56*WS[WS>0.5]^0.618
  hc[WS<=0.5] = 2.3+5.6*WS[WS<=0.5]^0.67

  hc = hc * 1.16
  Clothing.EVAP = 0.858-0.246*WS.effective+0.0372*WS.effective^2-0.0019*WS.effective^3
  Clothing.CONV = 0.577-0.0886*WS.effective+0.00516*WS.effective^2
  Clothing.RESP = 0.001173*(37-Tair)*Metab
  Conv.clothed = (hc*(Tskin-Tair)*Clothing.CONV*(Area.Clothed*SA))*1.08
  Conv.unclothed = hc*(Tskin-Tair)*(SA*(1-Area.Clothed))
  Convection = Clothing.RESP + Conv.clothed + Conv.unclothed
  Tclothing = Tair-21.9+Tskin
  Trad=(Tsky+Tsub)/2
  IR.Heat.loss = -1*(((sigma.emiss*(Trad+273.16)^4)-(sigma.emiss*(Tclothing+273.16)^4*Area.Clothed)-(sigma.emiss*(Tskin+273.16)^4*(1-Area.Clothed)))*(SA*0.725))
  #Evap.unclothed = (2.2*hc*((1.91*Tskin-25.33)-(RH/100)*(1.372+0.851*Tair-0.0161*Tair^2+0.00071*Tair^3)))*(SA*(1-Area.Clothed))
  Evap.unclothed = 2.2*hc*(WETAIR.rh(db=36, rh = 100)$vd*1000-WETAIR.rh(db=Tair, rh = RH)$vd*1000)*(SA*(1-Area.Clothed))
  #Evap.clothed = (2.2*hc*((1.91*Tskin-25.33)-(RH/100)*(1.372+0.851*Tair-0.0161*Tair^2+0.00071*Tair^3)))*(SA*(Area.Clothed))*Clothing.EVAP
  Evap.clothed = (2.2*hc*(WETAIR.rh(db=36, rh = 100)$vd*1000-WETAIR.rh(db=Tair, rh = RH)$vd*1000))*(SA*(Area.Clothed))*Clothing.EVAP
  Evap = Evap.unclothed+Evap.clothed
  Evap.unlimited = Evap
  Evap[Evap > Max.Evap.Rate] = Max.Evap.Rate
  #Evap.resp = Metab*0.0023*(44-(RH/100)*(1.372+0.851*Tair-0.0161*Tair^2+0.00071*Tair^3))
  Evap.resp = Metab*0.0023*(44-WETAIR.rh(db=Tair, rh = RH)$vd*1000)
  Evap.total = Evap + Evap.resp
  Storage = Total.Heat.Load - Convection - IR.Heat.loss - Evap.total
  Storage.uncapped = Storage
  Storage[Storage < 0] = 0
  Sweat[Total.Heat.Load > Convection + IR.Heat.loss] = 1
  Sweat[Total.Heat.Load <= Convection + IR.Heat.loss] = 0

  Sweat.limited[Evap.unclothed + Evap.clothed > Evap] = 1
  Sweat.limited[Storage == 0] = 0
  if(rasterdata==1){
   output = list(Body.Mass=Body.Mass, Height=Height, SA=SA, Activity.MR=Activity.MR, Maximum.SR=Maximum.SR, Temp.Increase=Temp.Increase, Area.Clothed=Area.Clothed, Albedo.body=Albedo.body, Tskin=Tskin, Speed=Speed, Evap.Eff=Evap.Eff, Tair=Tair, RH=RH, WS=WS, WS.effective=WS.effective, Solar.Rad=Solar.Rad, Solar.Load=Solar.Load, Total.Solar.Load.Incident.Absorbed=Total.Solar.Load.Incident.Absorbed, Metab=Metab, Total.Heat.Load=Total.Heat.Load, hc=hc, Clothing.EVAP=Clothing.EVAP, Clothing.CONV=Clothing.CONV, Clothing.RESP=Clothing.RESP, Conv.clothed=Conv.clothed, Conv.unclothed=Conv.unclothed, Convection=Convection, Tclothing=Tclothing, IR.Heat.loss=IR.Heat.loss, Evap.unclothed=Evap.unclothed, Evap.clothed=Evap.clothed, Evap=Evap, Evap.unlimited=Evap.unlimited, Evap.resp=Evap.resp, Evap.total=Evap.total, Storage=Storage, Storage.uncapped=Storage.uncapped, Max.Sweat.Rate=Max.Sweat.Rate, Max.Evap.Rate=Max.Evap.Rate, Sweat=Sweat, Sweat.limited=Sweat.limited)
  }else{
   output = cbind(Body.Mass, Height, SA, Activity.MR, Maximum.SR, Temp.Increase, Area.Clothed, Albedo.body, Tskin, Speed, Evap.Eff, Tair, RH, WS, WS.effective, Solar.Rad, Solar.Load, Total.Solar.Load.Incident.Absorbed, Metab, Total.Heat.Load, hc, Clothing.EVAP, Clothing.CONV, Clothing.RESP, Conv.clothed, Conv.unclothed, Convection, Tclothing, IR.Heat.loss, Evap.unclothed, Evap.clothed, Evap, Evap.unlimited, Evap.resp, Evap.total, Storage, Storage.uncapped, Max.Sweat.Rate, Max.Evap.Rate, Sweat, Sweat.limited)
  }
  return(output)
}

#output<-manmo(Tair = seq(0,50))

#plot(Sweat ~ Tair, data = output, type = 'l')

#write.csv(output,'manmo_out.csv')

