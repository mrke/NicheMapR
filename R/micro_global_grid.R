micro_global_grid <- function(loc="Madison, Wisconsin USA",timeinterval=12,nyears=1,soiltype=4,
  REFL=0.15,slope=0,aspect=0,
  DEP=c(0., 2.5,  5.,  10.,  15.,  20.,  30.,  50.,  100.,  200.),
  minshade=0,maxshade=90,Usrhyt=1,
  runshade=1,rungads=1,write_input=0,writecsv=0,
  ERR=1.5,RUF=0.004,EC=0.0167238,SLE=0.95,Thcond=2.5,Density=2560,SpecHeat=870,BulkDensity=1300,
  PCTWET=0,cap=1,CMH2O=1.,hori=rep(0,24),
  TIMAXS=c(1.0, 1.0, 0.0, 0.0),TIMINS=c(0.0, 0.0, 1.0, 1.0),timezone=0,
  runmoist=0,PE=rep(1.1,19),KS=rep(0.0037,19),BB=rep(4.5,19),BD=rep(1.3,19),Clay=20,
  maxpool=10000,rainmult=1,evenrain=0,
  SoilMoist_Init=c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3),
  L=c(0,0,8.18990859,7.991299442,7.796891252,7.420411664,7.059944542,6.385001059,5.768074989,
    4.816673431,4.0121088,1.833554792,0.946862989,0.635260544,0.804575,0.43525621,0.366052856,
    0,0)*10000,
  LAI=0.1,
  snowmodel=0,snowtemp=1.5,snowdens=0.375,densfun=c(0,0),snowmelt=0.9,undercatch=1,rainmelt=0.0125,
  rainfrac=0.5,
  shore=0,tides=matrix(data = 0., nrow = 24*timeinterval*nyears, ncol = 3)) {

  errors<-0

  # error trapping - originally inside the Fortran code, but now checking before executing Fortran
  if(DEP[2]-DEP[1]>3 | DEP[3]-DEP[2]>3){
    cat("warning, nodes might be too far apart near the surface, try a different spacing if the program is crashing \n")
  }
  if(timeinterval<12 | timeinterval > 365){
    cat("ERROR: the variable 'timeinterval' is out of bounds.
        Please enter a correct value (12 - 365).", '\n')
    errors<-1
  }
  if(is.numeric(loc[1])){
    if(loc[1]>180 | loc[2] > 90){
      cat("ERROR: Latitude or longitude (longlat) is out of bounds.
        Please enter a correct value.", '\n')
      errors<-1
    }
  }
  if(timezone%in%c(0,1)==FALSE){
    cat("ERROR: the variable 'timezone' be either 0 or 1.
      Please correct.", '\n')
    errors<-1
  }
  if(rungads%in%c(0,1)==FALSE){
    cat("ERROR: the variable 'rungads' be either 0 or 1.
      Please correct.", '\n')
    errors<-1
  }
  if(Clay<0){
    cat("ERROR: Clay density value (Clay) is negative.
      Please input a positive value.", '\n')
    errors<-1
  }
  if(write_input%in%c(0,1)==FALSE){
    cat("ERROR: the variable 'write_input' be either 0 or 1.
      Please correct.", '\n')
    errors<-1
  }
  if(EC<0.0034 | EC > 0.058){
    cat("ERROR: the eccentricity variable (EC) is out of bounds.
        Please enter a correct value (0.0034 - 0.058).", '\n')
    errors<-1
  }
  if(RUF<0.0001){
    cat("ERROR: The roughness height (RUF) is too small ( < 0.0001).
        Please enter a larger value.", '\n')
    errors<-1
  }
  if(RUF>2){
    cat("ERROR: The roughness height (RUF) is too large ( > 2).
        Please enter a smaller value.", '\n')
    errors<-1
  }
  if(DEP[1]!=0){
    cat("ERROR: First soil node (DEP[1]) must = 0 cm.
        Please correct", '\n')
    errors<-1
  }
  if(length(DEP)!=10){
    cat("ERROR: You must enter 10 different soil depths.", '\n')
    errors<-1
  }
  for(i in 1:9){
    if(DEP[i+1]<=DEP[i]){
      cat("ERROR: Soil depth (DEP array) is not in ascending size", '\n')
      errors<-1
    }
  }
  if(DEP[10]>500){
    cat("ERROR: Deepest soil depth (DEP array) is too large (<=500 cm)", '\n')
    errors<-1
  }
  if(Thcond<0){
    cat("ERROR: Thermal variable conductivity (THCOND) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(Density<0){
    cat("ERROR: Density variable (Density) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(SpecHeat<0){
    cat("ERROR: Specific heat variable (SpecHeat) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(BulkDensity<0){
    cat("ERROR: Bulk density value (BulkDensity) is negative.
        Please input a positive value.", '\n')
    errors<-1
  }
  if(REFL<0 | REFL>1){
    cat("ERROR: Soil reflectivity value (REFL) is out of bounds.
        Please input a value between 0 and 1.", '\n')
    errors<-1
  }
  if(slope<0 | slope>90){
    cat("ERROR: Slope value (slope) is out of bounds.
        Please input a value between 0 and 90.", '\n')
    errors<-1
  }
  if(aspect<0 | aspect>365){
    cat("ERROR: Aspect value (aspect) is out of bounds.
        Please input a value between 0 and 365.", '\n')
    errors<-1
  }
  if(max(hori)>90 | min(hori)<0){
    cat("ERROR: At least one of your horizon angles (hori) is out of bounds.
        Please input a value between 0 and 90", '\n')
    errors<-1
  }
  if(length(hori)!=24){
    cat("ERROR: You must enter 24 horizon angle values.", '\n')
    errors<-1
  }
  if(SLE<0.05 | SLE > 1){
    cat("ERROR: Emissivity (SLE) is out of bounds.
        Please enter a correct value (0.05 - 1.00).", '\n')
    errors<-1
  }
  if(ERR<0){
    cat("ERROR: Error bound (ERR) is too small.
        Please enter a correct value (> 0.00).", '\n')
    errors<-1
  }
  if(Usrhyt<RUF){
    cat("ERROR: Reference height (Usrhyt) smaller than roughness height (RUF).
        Please use a larger height above the surface.", '\n')
    errors<-1
  }
  if(Usrhyt<0.5 | Usrhyt>120){
    cat("ERROR: Reference height (Usrhyt) is out of bounds.
        Please enter a correct value (0.05 - 120).", '\n')
    errors<-1
  }
  if(CMH2O<0.5 | CMH2O>120){
    cat("ERROR: Preciptable water in air column (CMH2O) is out of bounds.
        Please enter a correct value (0.1 - 2).", '\n')
    errors<-1
  }
  if(max(TIMAXS)>24 | min(TIMAXS)<0){
    cat("ERROR: At least one of your times of weather maxima (TIMAXS) is out of bounds.
        Please input a value between 0 and 24", '\n')
    errors<-1
  }
  if(max(TIMINS)>24 | min(TIMINS)<0){
    cat("ERROR: At least one of your times of weather minima (TIMINS) is out of bounds.
        Please input a value between 0 and 24", '\n')
    errors<-1
  }
  if(minshade>maxshade | minshade==maxshade){
    cat("ERROR: Your value for minimum shade (minshade) is greater than or equal to the maximum shade (maxshade).
        Please correct this.", '\n')
    errors<-1
  }
  if(minshade>100 | minshade<0){
    cat("ERROR: Your value for minimum shade (minshade) is out of bounds.
        Please input a value between 0 and 100.", '\n')
    errors<-1
  }
  if(maxshade>100 | maxshade<0){
    cat("ERROR: Your value for maximum shade (maxshade) is out of bounds.
        Please input a value between 0 and 100.", '\n')
    errors<-1
  }
  if(soiltype<0 | soiltype>11){
    cat("ERROR: the soil type must range between 1 and 11.
      Please correct.", '\n')
    errors<-1
  }
  # end error trapping

  if(errors==0){ # continue
    if(rungads==1){
      if(library(GADS,quietly = TRUE,logical.return = TRUE)){
      }else{
        if(library(devtools, quietly = TRUE, logical.return = TRUE)){
          devtools::install_github('mrke/GADS', args="--no-multiarch")
          if(library(GADS,quietly = TRUE,logical.return = TRUE)){
          }else{
            stop("could not install GADS")
          }
        }else{
          install.packages("devtools")
          if(library(devtools, quietly = TRUE, logical.return = TRUE)){
            devtools::install_github('mrke/GADS', args="--no-multiarch")
            if(library(GADS,quietly = TRUE,logical.return = TRUE)){
            }else{
              stop("could not install GADS")
            }
          }else{
            stop("could not install devtools or GADS")
          }
        }
      }
    }
    ################## time related variables #################################

    DOYs12<-c(15.,46.,74.,105.,135.,166.,196.,227.,258.,288.,319.,349.) # middle day of each month
    DOYsn<-DOYs12 # variable of DOYs for when doing multiple years
    if(nyears>1 & timeinterval==365){ # create sequence of days for splining across multiple years
      for(i in 1:(nyears-1)){
        DOYsn<-c(DOYsn,(DOYs12+365*i))
      }
    }

    if(timeinterval<365){
      microdaily<-0 # run microclimate model as normal, where each day is iterated 3 times starting with the initial condition of uniform soil temp at mean monthly temperature
    }else{
      microdaily<-1 # run microclimate model where one iteration of each day occurs and last day gives initial conditions for present day with an initial 3 day burn in
    }

    # now check if doing something other than middle day of each month, and create appropriate vector of days of year
    daystart<-as.integer(ceiling(365/timeinterval/2))
    if(timeinterval!=12){
      DOYs<-seq(daystart,365,as.integer(floor(365/timeinterval)))
    }else{
      DOYs<-DOYsn
    }
    DOYnum <- timeinterval*nyears # total days to do
    DOY <- subset(DOYs, DOYs!=0) # final vector of days of year
    DOY<-rep(DOY,nyears)
    idayst <- 1 # start day
    ida<-timeinterval*nyears # end day

    ################## location and terrain #################################

    if(is.numeric(loc)==FALSE){ # use geocode to get location from site name via googlemaps
      if (!requireNamespace("dismo", quietly = TRUE)) {
        stop("dismo needed for the place name geocode function to work. Please install it.",
          call. = FALSE)
      }
      if (!requireNamespace("XML", quietly = TRUE)) {
        stop("XML needed for the place name geocode function to work. Please install it.",
          call. = FALSE)
      }
      longlat <- dismo::geocode(loc)[3:4] # assumes first geocode match is correct
      if(nrow(longlat>1)){longlat<-longlat[1,]}
      x <- t(as.matrix(as.numeric(c(longlat[1,1],longlat[1,2]))))
    }else{
      longlat <- loc
      x <- t(as.matrix(as.numeric(c(loc[1],loc[2]))))
    }

    # get the local timezone reference longitude
    if(timezone==1){ # this now requires registration
      if(!require(geonames, quietly = TRUE)){
        stop('package "geonames" is required to do a specific time zone (timezone=1). Please install it.')
      }
      ALREF<-(geonames::GNtimezone(longlat[2],longlat[1])[4])*-15
    }else{  # just use local solar noon
      ALREF <- abs(trunc(x[1]))
    }
    HEMIS <- ifelse(x[2]<0,2.,1.) # 1 is northern hemisphere
    # break decimal degree lat/lon into deg and min
    ALAT <- abs(trunc(x[2]))
    AMINUT <- (abs(x[2])-ALAT)*60
    ALONG <- abs(trunc(x[1]))
    ALMINT <- (abs(x[1])-ALONG)*60
    azmuth<-aspect

    hori<-as.matrix(hori) #horizon angles
    VIEWF <- 1-sum(sin(as.data.frame(hori)*pi/180))/length(hori) # convert horizon angles to radians and calc view factor(s)
    SLES <- rep(SLE,timeinterval*nyears)
    # creating the shade array
    MAXSHADES <- rep(0,(timeinterval*nyears))+maxshade # daily max shade (%)
    MINSHADES <- rep(0,(timeinterval*nyears))+minshade # daily min shade (%)

    # load global climate files
    gcfolder<-paste(.libPaths()[1],"/gcfolder.rda",sep="")
    if(file.exists(gcfolder)==FALSE){
      cat("You don't appear to have the global climate data set - \n run function get.global.climate(folder = 'folder you want to put it in') .....\n exiting function micro_global")
      opt <- options(show.error.messages=FALSE)
      on.exit(options(opt))
      stop()
    }
    load(gcfolder)
    global_climate<-raster::brick(paste(folder,"/global_climate.nc",sep=""))
    elev<-raster::raster(paste(folder,"/elev.nc",sep=""))
    soilmoisture<-suppressWarnings(raster::brick(paste(folder,"/soilw.mon.ltm.v2.nc",sep="")))

    ALTT<-as.numeric(raster::extract(elev,x)*1000) # convert from km to m
    cat("extracting climate data", '\n')
    CLIMATE <- raster::extract(global_climate,x)
    RAINFALL <- CLIMATE[,1:12]
    WNMAXX <- CLIMATE[,13:24]
    WNMINN<-WNMAXX*0.1 # impose diurnal cycle
    TMINN <- CLIMATE[,25:36]
    TMAXX <- CLIMATE[,37:48]
    ALLMINTEMPS<-TMINN
    ALLMAXTEMPS<-TMAXX
    ALLTEMPS <- cbind(ALLMAXTEMPS,ALLMINTEMPS)
    CCMAXX <- 100-CLIMATE[,49:60]
    CCMINN<-CCMAXX
    RAINYDAYS <- CLIMATE[,61:72]
    RHMINN <- CLIMATE[,73:84]
    RHMAXX <- CLIMATE[,85:96]
    if(soiltype==0){ # simulating rock so turn of soil moisture model and set density equal to bulk density
      BulkDensity<-Density
      cap=0
      runmoist<-0
      PE<-rep(CampNormTbl9_1[1,4],19) #air entry potential J/kg
      KS<-rep(CampNormTbl9_1[1,6],19) #saturated conductivity, kg s/m3
      BB<-rep(CampNormTbl9_1[1,5],19) #soil 'b' parameter
      BD<-rep(BulkDensity/1000,19) # soil bulk density, Mg/m3
    }else{
      if(soiltype<12){ # use soil properties as specified in Campbell and Norman 1998 Table 9.1
      E<-rep(CampNormTbl9_1[soiltype,4],19) #air entry potential J/kg
      KS<-rep(CampNormTbl9_1[soiltype,6],19) #saturated conductivity, kg s/m3
      BB<-rep(CampNormTbl9_1[soiltype,5],19) #soil 'b' parameter
      BD<-rep(BulkDensity/1000,19) # soil bulk density, Mg/m3
      }
    }

    if(runmoist==0){
      # extract soil moisture (this database has longitude from 0-365!)
      longlat1<-longlat
      if(longlat1[1]<0){
        longlat1[1]<-360+longlat1[1]
      }
      if(is.numeric(loc)==FALSE){
        xx<-cbind(longlat1)
      }else{
        xx<-t(cbind(longlat1))
      }
      xx2<-as.numeric(xx)
      cat("extracting soil moisture data", '\n')
      SoilMoist<-raster::extract(soilmoisture,xx)/1000 # this is originally in mm/m
      if(nrow(SoilMoist)>1){SoilMoist<-SoilMoist[1,]}
    }
    if(is.na(max(SoilMoist, ALTT, CLIMATE))==TRUE){
      cat("Sorry, there is no environmental data for this location")
      stop()
    }
    # correct for fact that wind is measured at 10 m height
    # wind shear equation v / vo = (h / ho)^a
    #where
    #v = the velocity at height h (m/s)
    #vo = the velocity at height ho (m/s)
    #a = the wind shear exponent
    #Terrain   Wind Shear Exponent
    #- a -
    #  Open water   0.1
    #Smooth, level, grass-covered   0.15 (or more commonly 1/7)
    #Row crops 	0.2
    #Low bushes with a few trees 	0.2
    #Heavy trees 	0.25
    #Several buildings 	0.25
    #Hilly, mountainous terrain 	0.25
    # source http://www.engineeringtoolbox.com/wind-shear-d_1215.html
    WNMINN<-WNMINN*(1.2/10)^0.15
    WNMAXX<-WNMAXX*(1.2/10)^0.15

    if(timeinterval!=12){ # spline from 12 days to chosen time interval
      TMAXX1 <-suppressWarnings(spline(DOYs12,TMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      TMAXX<-rep(TMAXX1$y,nyears)
      TMINN1 <-suppressWarnings(spline(DOYs12,TMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      TMINN <- rep(TMINN1$y,nyears)
      RHMAXX1 <-suppressWarnings(spline(DOYs12,RHMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      RHMAXX <- rep(RHMAXX1$y,nyears)
      RHMINN1 <-suppressWarnings(spline(DOYs12,RHMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      RHMINN <- rep(RHMINN1$y,nyears)
      CCMAXX1 <-suppressWarnings(spline(DOYs12,CCMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      CCMAXX <- rep(CCMAXX1$y,nyears)
      CCMINN <- CCMAXX
      WNMAXX1 <-suppressWarnings(spline(DOYs12,WNMAXX,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      WNMAXX<-rep(WNMAXX1$y,nyears)
      WNMINN1 <-suppressWarnings(spline(DOYs12,WNMINN,n=timeinterval,xmin=1,xmax=365,method="periodic"))
      WNMINN<-rep(WNMINN1$y,nyears)
      if(runmoist==0){
        SoilMoist1 <-suppressWarnings(spline(DOYs12,SoilMoist,n=timeinterval,xmin=1,xmax=365,method="periodic"))
        SoilMoist<-rep(SoilMoist1$y,nyears)
      }
    }
    if(timeinterval<365){
      TMAXX<-rep(TMAXX,nyears)
      TMINN<-rep(TMINN,nyears)
      RHMAXX<-rep(RHMAXX,nyears)
      RHMINN<-rep(RHMINN,nyears)
      CCMAXX<-rep(CCMAXX,nyears)
      CCMINN<-rep(CCMINN,nyears)
      WNMAXX<-rep(WNMAXX,nyears)
      WNMINN<-rep(WNMINN,nyears)
      if(runmoist==0){
        SoilMoist<-rep(SoilMoist,nyears)
      }
      RAINFALL<-rep(RAINFALL,nyears)
    }
    # get annual mean temp for creating deep soil (2m) boundary condition
    avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
    soilinit<-rep(avetemp,20)
    tannul<-mean(unlist(ALLTEMPS))
    tannulrun<-rep(tannul,DOYnum)

    daymon<-c(31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.) # days in each month

    # if doing daily sims, spread rainfall evenly across days based on mean monthly rainfall and the number of rainy days per month
    if(timeinterval==365){
      RAINFALL1<-1:365
      sort<-matrix(data = 0,nrow = 365,ncol = 2)
      daymon<-c(31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31.) # days in each month
      m<-1
      b<-0
      for (i in 1:12){ #begin loop throught 12 months of year
        ndays=daymon[i]
        for (k in 1:ndays){
          b<-b+1
          sort[m,1]<-i
          sort[m,2]<-b
          if(k<=RAINYDAYS[i] & rainfrac>0){
            if(k==1){
              RAINFALL1[m]=RAINFALL[i]*rainfrac*rainmult # if first day of month, make user-specified fraction of monthly rainfall fall on first day
            }else{
              RAINFALL1[m]=(RAINFALL[i]*(1-rainfrac)*rainmult)/RAINYDAYS[i] # make remaining rain fall evenly over the remaining number of rainy days for the month, starting at the beginning of the month
            }
          }else{
            if(rainfrac==0){
              RAINFALL1[m]=(RAINFALL[i]*rainmult)/RAINYDAYS[i]
            }else{
              RAINFALL1[m]=0.
            }
          }
          m<-m+1
          if(b>RAINYDAYS[i]){
            b<-0
          }
        }
      }
      RAINFALL2<-as.data.frame(cbind(RAINFALL1,sort))
      #RAINFALL2<-RAINFALL2[order(RAINFALL2$V2,RAINFALL2$V3),] # this line scatters the rainy days evenly across each month - snow predictions better if it is commented out so get rainy days all in a row within the month
      RAINFALL<-rep(as.double(RAINFALL2$RAINFALL1),nyears)
      if(TMINN[1]<snowtemp){
        RAINFALL[1]<-0 # this is needed in some cases to allow the integrator to get started
      }
    }else{
      if(timeinterval!=12){
        RAINFALL<-rep(rep(sum(RAINFALL)/timeinterval,timeinterval),nyears) # just spread evenly across every day
      }else{ # running middle day of each month - divide monthly rain by number of days in month
        RAINFALL<-RAINFALL/rep(daymon,nyears)
      }
    }#end check doing daily sims
    dim<-length(RAINFALL)
    if(rungads==1){
      ####### get solar attenuation due to aerosols with program GADS #####################
      relhum<-1.
      season<-0.
      optdep.summer<-as.data.frame(GADS::get.gads(longlat[2],longlat[1],relhum,0))
      optdep.winter<-as.data.frame(GADS::get.gads(longlat[2],longlat[1],relhum,1))
      optdep<-cbind(optdep.winter[,1],rowMeans(cbind(optdep.summer[,2],optdep.winter[,2])))
      optdep<-as.data.frame(optdep)
      colnames(optdep)<-c("LAMBDA","OPTDEPTH")
      a<-lm(OPTDEPTH~poly(LAMBDA, 6, raw=TRUE),data=optdep)
      LAMBDA<-c(290,295,300,305,310,315,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,1080,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,1380,1400,1420,1440,1460,1480,1500,1540,1580,1600,1620,1640,1660,1700,1720,1780,1800,1860,1900,1950,2000,2020,2050,2100,2120,2150,2200,2260,2300,2320,2350,2380,2400,2420,2450,2490,2500,2600,2700,2800,2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000)
      TAI<-predict(a,data.frame(LAMBDA))
      ################ end GADS ##################################################
    }else{ # use the original profile from Elterman, L. 1970. Vertical-attenuation model with eight surface meteorological ranges 2 to 13 kilometers. U. S. Airforce Cambridge Research Laboratory, Bedford, Mass.
      TAI<-c(0.0670358341290886,0.0662612704779235,0.065497075238002,0.0647431301168489,0.0639993178022531,0.0632655219571553,0.0625416272145492,0.0611230843885423,0.0597427855962549,0.0583998423063099,0.0570933810229656,0.0558225431259535,0.0545864847111214,0.0533843764318805,0.0522154033414562,0.0499736739981675,0.047855059159556,0.0458535417401334,0.0439633201842001,0.0421788036108921,0.0404946070106968,0.0389055464934382,0.0374066345877315,0.0359930755919066,0.0346602609764008,0.0334037648376212,0.0322193394032758,0.0311029105891739,0.0300505736074963,0.0290585886265337,0.0281233764818952,0.0272415144391857,0.0264097320081524,0.0256249068083005,0.0248840604859789,0.0241843546829336,0.0235230870563317,0.0228976873502544,0.0223057135186581,0.0217448478998064,0.0212128934421699,0.0207077699817964,0.0202275105711489,0.0197702578594144,0.0193342605242809,0.0189178697551836,0.0177713140039894,0.0174187914242432,0.0170790495503944,0.0167509836728154,0.0164335684174899,0.0161258546410128,0.0158269663770596,0.0155360978343254,0.0152525104459325,0.0149755299703076,0.0147045436435285,0.0144389973831391,0.0141783930434343,0.0134220329447663,0.0131772403830191,0.0129356456025128,0.0126970313213065,0.0124612184223418,0.0122280636204822,0.01199745718102,0.0115436048739351,0.0110993711778668,0.0108808815754663,0.0106648652077878,0.0104513876347606,0.0102405315676965,0.00982708969547694,0.00962473896278535,0.00903679230300494,0.00884767454432418,0.0083031278398166,0.00796072474935954,0.00755817587626185,0.00718610751850881,0.00704629977586921,0.00684663903049612,0.00654155580333479,0.00642947339729728,0.00627223096874308,0.00603955966866779,0.00580920937536261,0.00568506186880564,0.00563167068287251,0.00556222005081865,0.00550522989971023,0.00547395763028062,0.0054478983436216,0.00541823364504573,0.00539532163908382,0.00539239864119488,0.00541690124712384,0.00551525885358836,0.00564825853509463,0.00577220185074264,0.00584222986640171,0.00581645238345584,0.00566088137411449,0.00535516862329704,0.00489914757707667,0.00432017939770409,0.0036813032251836,0.00309019064543606,0.00270890436501562,0.00276446109239711,0.00356019862584603)
    } #end check if running gads

    ################ soil properties  ##################################################

    Intrvls <-(1:dim) # user-supplied last day of year in each time interval sequence
    Numint <- dim  # number of time intervals
    Numtyps <- 2 # number of soil types
    Nodes <- matrix(data = 0, nrow = 10, ncol = dim) # array of all possible soil nodes for max time span of 20 years
    Nodes[1,1:dim]<-3 # deepest node for first substrate type
    Nodes[2,1:dim]<-9 # deepest node for second substrate type
    REFLS<-rep(REFL,dim) # soil reflectances
    PCTWET<-rep(PCTWET,dim) # soil wetness
    Density<-Density/1000 # density of minerals - convert to Mg/m3
    BulkDensity<-BulkDensity/1000 # density of minerals - convert to Mg/m3
    if(runmoist==0){
      moists2<-matrix(nrow= 10, ncol = dim, data=0) # set up an empty vector for soil moisture values through time
      moists2[1,]<-SoilMoist # fill the first row with monthly soil moisture values
      moists2[2,]<-moists2[1,] # make this row same as first row
      moists<-moists2
    }else{
      moists2<-matrix(nrow=10, ncol = dim, data=0) # set up an empty vector for soil moisture values through time
      moists2[1:10,]<-SoilMoist_Init
      moists2[moists2>(1-BulkDensity/Density)]<-(1-BulkDensity/Density)
      moists<-moists2
    }


    # now make the soil properties matrix
    # columns are:
    #1) bulk density (Mg/m3)
    #2) volumetric water content at saturation (0.1 bar matric potential) (m3/m3)
    #3) clay content (%)
    #4) thermal conductivity (W/mK)
    #5) specific heat capacity (J/kg-K)
    #6) mineral density (Mg/m3)
    soilprops<-matrix(data = 0, nrow = 10, ncol = 6) # create an empty soil properties matrix
    soilprops[1,1]<-BulkDensity # insert soil bulk density to profile 1
    soilprops[2,1]<-BulkDensity # insert soil bulk density to profile 2
    soilprops[1,2]<-min(0.26,1-BulkDensity/Density) # insert saturated water content to profile 1
    soilprops[2,2]<-min(0.26,1-BulkDensity/Density) # insert saturated water content to profile 2
    soilprops[1,3]<-Clay     # insert \% clay to profile 1
    soilprops[2,3]<-Clay     # insert\% clay to profile 2
    if(cap==1){ # insert thermal conductivity to profile 1, and see if 'organic cap' added on top
      soilprops[1,4]<-0.2 # mineral thermal conductivity
    }else{
      soilprops[1,4]<-Thcond # mineral thermal conductivity
    }
    soilprops[2,4]<-Thcond # insert thermal conductivity to profile 2
    if(cap==1){ # insert specific heat to profile 1, and see if 'organic cap' added on top
      soilprops[1,5]<-1920 # mineral heat capacity
    }else{
      soilprops[1,5]<-SpecHeat
    }
    soilprops[2,5]<-SpecHeat # insert specific heat to profile 2
    soilprops[1,6]<-Density # insert mineral density to profile 1
    soilprops[2,6]<-Density # insert mineral density to profile 2
    #########################################################################################

    # Next four parameters are segmented velocity profiles due to bushes, rocks etc. on the surface
    #IF NO EXPERIMENTAL WIND PROFILE DATA SET ALL THESE TO ZERO! (then roughness height is based on the parameter RUF)
    Z01 <- 0. # Top (1st) segment roughness height(m)
    Z02 <- 0. # 2nd segment roughness height(m)
    ZH1 <- 0. # Top of (1st) segment, height above surface(m)
    ZH2 <- 0. # 2nd segment, height above surface(m)
    SNOW <- rep(0,timeinterval*nyears) # no snow simulated on surface

    # microclimate input parameters list
    microinput<-c(dim,RUF,ERR,Usrhyt,Numtyps,Numint,Z01,Z02,ZH1,ZH2,idayst,ida,HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,CMH2O,microdaily,tannul,EC,VIEWF,snowtemp,snowdens,snowmelt,undercatch,rainmult,runshade,runmoist,maxpool,evenrain,snowmodel,rainmelt,writecsv,densfun)

    DOY1=matrix(data = 0., nrow = dim, ncol = 1)
    SLES1=matrix(data = 0., nrow = dim, ncol = 1)
    Intrvls1=matrix(data = 0., nrow = dim, ncol = 1)
    MAXSHADES1=matrix(data = 0., nrow = dim, ncol = 1)
    MINSHADES1=matrix(data = 0., nrow = dim, ncol = 1)
    TMAXX1=matrix(data = 0., nrow = dim, ncol = 1)
    TMINN1=matrix(data = 0., nrow = dim, ncol = 1)
    CCMAXX1=matrix(data = 0., nrow = dim, ncol = 1)
    CCMINN1=matrix(data = 0., nrow = dim, ncol = 1)
    RHMAXX1=matrix(data = 0., nrow = dim, ncol = 1)
    RHMINN1=matrix(data = 0., nrow = dim, ncol = 1)
    WNMAXX1=matrix(data = 0., nrow = dim, ncol = 1)
    WNMINN1=matrix(data = 0., nrow = dim, ncol = 1)
    SNOW1=matrix(data = 0., nrow = dim, ncol = 1)
    REFLS1=matrix(data = 0., nrow = dim, ncol = 1)
    PCTWET1=matrix(data = 0., nrow = dim, ncol = 1)
    RAINFALL1=matrix(data = 0, nrow = dim, ncol = 1)
    tannul1=matrix(data = 0, nrow = dim, ncol = 1)
    moists1=matrix(data = 0., nrow = 10, ncol = dim)
    DOY1[1:dim]<-DOY
    SLES1[1:dim]<-SLES
    Intrvls1[1:dim]<-Intrvls
    MAXSHADES1[1:dim]<-MAXSHADES
    MINSHADES1[1:dim]<-MINSHADES
    TMAXX1[1:dim]<-TMAXX
    TMINN1[1:dim]<-TMINN
    CCMAXX1[1:dim]<-CCMAXX
    CCMINN1[1:dim]<-CCMINN
    RHMAXX1[1:dim]<-RHMAXX
    RHMINN1[1:dim]<-RHMINN
    WNMAXX1[1:dim]<-WNMAXX
    WNMINN1[1:dim]<-WNMINN
    SNOW1[1:dim]<-SNOW
    REFLS1[1:dim]<-REFLS
    PCTWET1[1:dim]<-PCTWET
    RAINFALL1[1:dim]<-RAINFALL
    tannul1[1:dim]<-tannul
    moists1[1:10,1:dim]<-moists

    if(shore==0){
      tides<-matrix(data = 0., nrow = 24*dim, ncol = 3) # make an empty matrix
    }
    # all microclimate data input list - all these variables are expected by the input argument of the fortran micro2014 subroutine
    micro<-list(tides=tides,microinput=microinput,DOY=DOY,SLES=SLES1,DEP=DEP,Intrvls=Intrvls1,Nodes=Nodes,MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TIMAXS=TIMAXS,TIMINS=TIMINS,TMAXX=TMAXX1,TMINN=TMINN1,RHMAXX=RHMAXX1,RHMINN=RHMINN1,CCMAXX=CCMAXX1,CCMINN=CCMINN1,WNMAXX=WNMAXX1,WNMINN=WNMINN1,SNOW=SNOW1,REFLS=REFLS1,PCTWET=PCTWET1,soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists1,RAINFALL=RAINFALL1,tannulrun=tannulrun,PE=PE,KS=KS,BB=BB,BD=BD,L=L,LAI=LAI,snowmodel=snowmodel)

    # write all input to csv files in their own folder
    if(write_input==1){
      if(dir.exists("micro csv input")==FALSE){
        dir.create("micro csv input")
      }
      write.table(as.matrix(microinput), file = "micro csv input/microinput.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(DOY, file = "micro csv input/DOY.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(SLES, file = "micro csv input/SLES.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(DEP, file = "micro csv input/DEP.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(Intrvls, file = "micro csv input/Intrvls.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(Nodes, file = "micro csv input/Nodes.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(MAXSHADES, file = "micro csv input/Maxshades.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(MINSHADES, file = "micro csv input/Minshades.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TIMAXS, file = "micro csv input/TIMAXS.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TIMINS, file = "micro csv input/TIMINS.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TMAXX, file = "micro csv input/TMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TMINN, file = "micro csv input/TMINN.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(RHMAXX, file = "micro csv input/RHMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(RHMINN, file = "micro csv input/RHMINN.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(CCMAXX, file = "micro csv input/CCMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(CCMINN, file = "micro csv input/CCMINN.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(WNMAXX, file = "micro csv input/WNMAXX.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(WNMINN, file = "micro csv input/WNMINN.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(SNOW, file = "micro csv input/SNOW.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(REFLS, file = "micro csv input/REFLS.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(PCTWET, file = "micro csv input/PCTWET.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(soilinit, file = "micro csv input/soilinit.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(hori, file = "micro csv input/hori.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(TAI, file = "micro csv input/TAI.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(soilprops, file="micro csv input/soilprop.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(moists,file="micro csv input/moists.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(RAINFALL,file="micro csv input/rain.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(tannulrun,file="micro csv input/tannulrun.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(PE,file="micro csv input/PE.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(BD,file="micro csv input/BD.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(BB,file="micro csv input/BB.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(KS,file="micro csv input/KS.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(L,file="micro csv input/L.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(LAI,file="micro csv input/LAI.csv", sep = ",", col.names = NA, qmethod = "double")
      write.table(tides,file="micro csv input/tides.csv", sep = ",", col.names = NA, qmethod = "double")
    }

    if(is.numeric(loc[1])){
      location<-paste("long",loc[1],"lat",loc[2])
    }else{
      location<-loc
    }
    cat(paste('running microclimate model for',timeinterval,'days by',nyears,'years at site',location,'\n'))
    ptm <- proc.time() # Start timing
    microut<-microclimate(micro)
    print(proc.time() - ptm) # Stop the clock

    metout<-microut$metout # retrieve above ground microclimatic conditions, min shade
    shadmet<-microut$shadmet # retrieve above ground microclimatic conditions, max shade
    soil<-microut$soil # retrieve soil temperatures, minimum shade
    shadsoil<-microut$shadsoil # retrieve soil temperatures, maximum shade
    if(runmoist==1){
      soilmoist<-microut$soilmoist # retrieve soil moisture, minimum shade
      shadmoist<-microut$shadmoist # retrieve soil moisture, maximum shade
      humid<-microut$humid # retrieve soil humidity, minimum shade
      shadhumid<-microut$shadhumid # retrieve soil humidity, maximum shade
      soilpot<-microut$soilpot # retrieve soil water potential, minimum shade
      shadpot<-microut$shadpot # retrieve soil water potential, maximum shade
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
    result<-as.matrix(cbind(longlat,max(soil[,3])))
  } # end error trapping
} # end of NicheMapR_Setup_micro function
