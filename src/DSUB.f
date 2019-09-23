      SUBROUTINE DSUB(TIME,T,DTDT)
      use commondat
      IMPLICIT NONE
      EXTERNAL TAB

C     NicheMapR: software for biophysical mechanistic niche modelling

C     Copyright (C) 2018 Michael R. Kearney and Warren P. Porter

c     This program is free software: you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published by
c     the Free Software Foundation, either version 3 of the License, or (at
c      your option) any later version.

c     This program is distributed in the hope that it will be useful, but
c     WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c     General Public License for more details.

c     You should have received a copy of the GNU General Public License
c     along with this program. If not, see http://www.gnu.org/licenses/.

      double precision A,ALAT,ALONC,ALTT,AMOL,AMULT,ARAD,AZMUTH
      double precision B,BEGP,BP,C,CC,CLEAR,CLOD,CLOUD,CLR,CMH2O,CP,CRAD
     & ,CS,CZ,CZSL
      double precision D,D0,DAS,DENAIR,DENDAY,DEPP,DEP,DP,DTAU,DTDT
      double precision ERRP,E,END,ERR1,ESAT,F,F1,F2,FIN,GRDIN
      double precision H,HC,HD,HEMIS,HGTP,HTRN
      double precision MASS,MAXSHD,MON,OUT
      double precision PATMOS,PCTEVP,PI,PLT,PRESS,PRTP,PSTD,PTWET,PUNSH
      double precision QCOND,QCONV,QEVAP
      double precision QRAD,QRADGR,QRADSK,QRADVG,QSOLAR,RB,RCSP,REFL,RH
     & ,RS
      double precision RUFP,RW,SAB,SABNEW,SHAYD,SHDLIZ,SIGP,SIOUT,SKYIN
     & ,SLEP
      double precision SLOPE,SMET,SOK,SOLR,SPDAY,SRAD,STP,SUN,IRDOWN
      double precision T,TAB,TAIR,TCI,TDS,TGRD,TKDAY,TIM,TIMCOR,TIME
     & ,TIMEF
      double precision TS,TSI,TSKY,TSNHR,TSRHR,TTEST,TVINC,TVIR
      double precision VD,VELR,VV,WB,WC,WTRPOT,X,XXX
      double precision ZENR,ZSLR,ZZ,Z01,Z02,ZH1,ZH2,HRAD,QRADHL,VIEWF,TT
      double precision sle,err,soilprop,moist,Thconduct,Density,Spheat
     & ,snowage,grasshade,ZH
      double precision rainfall,minsnow,inrad,refrad,snowcond,intercept
     & ,prevden
      double precision condep,rainmult,surflux,ep,maxpool,tide
      double precision snownode,maxsnode1,snode,daysincesnow,lastday,
     &undercatch,rainmeltf,snowalbedo,tt_past,densfun
      double precision rww,pc,rl,sp,r1,im
      double precision snowdens,snowmelt,snowtemp,cursnow,cummelted
     & ,qfreze,xtrain,qphase,sumphase,sumphase2,NONP

      INTEGER I,II,I1,I2,IALT,IDA,IDAYST,IEND,IEP,maxsnode2
      INTEGER IOUT,IPINT,IPRINT,ISTART,ITEST
      INTEGER IUV,J,JULNUM,L,methour,trouble,IRmode
      INTEGER M1,DOY,N,N1,NAIR,NOSCAT,NOUT,hourly,rainhourly
      INTEGER Numtyps,hour,runmoist,evenrain,sat,runsnow,maxcount

      CHARACTER(3) INAME,SYMBOL

      DIMENSION T(30),TT(30), DTDT(18), DEPP(30),sumphase2(10)
      DIMENSION DENDAY(30),SPDAY(30),TKDAY(30),densfun(4),qphase(10)
      DIMENSION soilprop(10,5),moist(10),snownode(10),snode(10)
      DIMENSION Thconduct(30),Density(30),Spheat(30),tt_past(30)

      COMMON/AIRRAY/ZZ(10),VV(10)
      COMMON/DMYCRO/Z01,Z02,ZH1,ZH2
      COMMON/CMYCRO/ZH,D0
      COMMON/WDSUB/TSKY,ARAD,CRAD,CLOUD,CLR,SOLR
      COMMON/PAR/SIGP,RCSP,SOK,SAB,HGTP,RUFP,BEGP,MON,PRTP,ERRP,END,
     1 SLEP,DAS,NONP,SUN,PLT,FIN,STP
      COMMON/NONSCR/M1,N1,XXX,TIMEF,DTAU,ERR1,H,NOUT,NAIR,IPRINT
      COMMON/ARRAY/X(30),WC(20),C(20),DEP(30),OUT(101),IOUT(100),
     1 ITEST(23)
      COMMON/CARRAY/INAME(20),SYMBOL(23)
      COMMON/FUN1/SKYIN,GRDIN,MASS,PCTEVP,TCI,TSI,TIM,TGRD,CC,CS,RS,
     *A,B,D,F,HTRN,F1,F2,TS,SMET,RB,SHDLIZ
      COMMON/SIUNIT/SIOUT(10)
      COMMON/WIOCONS/PUNSH,ALAT,AMULT,PRESS,CMH2O,REFL,ALONC,TIMCOR,
     * AZMUTH,SLOPE,TSNHR,TSRHR,Hemis
      COMMON/WIOCONS2/IPINT,NOSCAT,IUV,IALT,IDAYST,IDA,IEP,ISTART,IEND
      COMMON/GROUND/SHAYD,ALTT,MAXSHD,SABNEW,PTWET,rainfall
      COMMON/SNOWPRED/snowtemp,snowdens,snowmelt,snownode,minsnow
     &,maxsnode1,snode,cursnow,daysincesnow,lastday,undercatch,rainmeltf
     &,densfun,snowcond,intercept,snowage,prevden,grasshade
      common/soilmoist/condep,rainmult,maxpool
      common/soilmoist3/runmoist,evenrain,maxcount
      common/soimoist2/rww,pc,rl,sp,r1,im
      common/hour/hourly,rainhourly
      common/IR/IRmode
      COMMON/DAYJUL/JULNUM,DOY
c    Variable soil properties data from Iomet1
      COMMON/SOYVAR1/Numtyps
      COMMON/SOYVAR2/Thconduct,Density,Spheat
      COMMON/SOYFILS/DENDAY,SPDAY,TKDAY
      COMMON/VIEWFACT/VIEWF
      COMMON/NICHEMAPRIO/SLE,ERR,soilprop,surflux
      common/moistcom/moist,ep
      common/snowmod/runsnow,trouble
      COMMON/temps/TT,TT_past,cummelted
      COMMON/WOSUB/DEPP
      COMMON/melt/QFREZE,xtrain,qphase,sumphase,sumphase2

C     NOTATION
C     Key Variables
C     DOY = 'DAY OF YEAR', = simulation day number
C     IDA = TOTAL DAYS OF SIMULATION
C     IDAST = STARTING DAY OF SIMULATION
C     Numtyps = # of substrate types, e.g. snow, rock, sandy soil
C     Intrvls = duration of each vertical soil configuration just defined above.
C     Thconds = substrate conductivities, e.g. snow, soil type(s)
C     Densitys = substrate densities, e.g. snow, soil type(s)
C     Spheats = substrate specific heats, e.g. snow, soil type(s)
C     Nodes = Deepest node for the each substrate type number for each time interval (duration) of vertical arrangement of substrates
C     Nodes(max node depth,subst type) are real numbers. The number to the left of the decimal point is the deepest node for the substrate type, which is to the right of the decimal point.
      j=1
      if(runsnow.eq.1)then
       if(cursnow.lt.minsnow)then
        maxsnode1=0.D0
        DEPP(9)=0.D0
       endif
       do 555 i=1,8
        if(snode(i).lt.1e-8)then
         t(i)=t(1)
         tt(i)=tt(1)
         if((snode(8).lt.1e-8).and.(cursnow.lt.1e-8))then
          t(9)=t(1)
          tt(9)=tt(1)
         endif
        endif
555    continue
      endif

      sat=1
      hour=int(time/60+1)
      if(hour.eq.0)then
       hour=1
      endif
      N=N1
      PI=3.14159
C     CHECK FOR UNSTABLE CONDITIONS OF GROUND SURFACE TEMPERATURE, T(1)
      do 101 i=1,N
       IF(T(i).GE. 81)THEN
        T(i) = 81.
       ELSE
        IF(T(i).LT.-81)THEN
         T(i) = -81
        ENDIF
       ENDIF
101   continue

c     IF(TIME-BEGP) 1,1,100
C**** FIND HEIGHTS ABOVE GROUND FOR DETERMINING TEMP AND VELOCITY AND
C**** PUT DEPTHS IN DEPP(I) AND HEIGHTS IN ZZ(I)
C     HEIGHTS IN ZZ ARE FROM THE GROUND UP, SO ARE VELOCITIES, VV,
C     BUT DO NOT INCLUDE THE SURFACE
c    1 CONTINUE
      DO 2 I=1,10
       VV(I)=0.D0
       ZZ(I)=0
    2 CONTINUE
      NAIR=0
   72 NAIR=NAIR+1
      IF(DEP(NAIR)) 72,73,73
 73   NAIR=NAIR-1
      DO 74 I=1,N
       J=I+NAIR
       if(runsnow.eq.1)then
        if(cursnow.ge.minsnow)then
         if(i.le.8)then
          DEPP(I)=snode(i)
         else
          DEPP(I)=DEP(J-8)+cursnow
         endif
        else
         if(i.le.8)then
          DEPP(I)=snode(i)
         else
          DEPP(I)=DEP(J-8)
         endif
        endif
       else
        DEPP(I)=DEP(J)
       endif
   74 continue
      IF(NAIR.LE.0) GO TO 77
      DO 76 I=1,NAIR
       J=NAIR+1-I
       ZZ(I)=ABS(DEP(J))
   76 CONTINUE
   77 CONTINUE

C     SURFACE ABSORPTIVITY FOR SUB. DSUB CALCULATIONS
      SABNEW = 1.000 - REFLS(doy)
C     SURFACE REFLECTIVITY FOR SOLRAD CALCULATIONS
      REFL = REFLS(doy)
      SLE = SLES(doy)
      moist=moists(1:10,doy)
      if(runsnow.eq.1)then
       methour=(int(SIOUT(1)/60)+1)+24*(doy-1)
       if(methour.gt.1)then
        if(snowhr(methour-1).gt.0)then
         if(daysincesnow.lt.0.3)then
          daysincesnow=0.3
         endif
         snowalbedo=(-9.8740*dlog(daysincesnow) + 78.3434)/100.
         SABNEW = 1-snowalbedo
C        SETTING THIS MONTH'S PERCENT OF SURFACE WITH FREE WATER/SNOW ON IT
         PTWET = 100
        endif
       endif
      endif
      if(condep.gt.10.)then ! use Fresnel's reflection law, p. 212 Gates 1980
        inrad = TAB('ZEN',TIME)*pi/180.
        refrad=asin(sin(inrad)/1.33)
        SABNEW = 1.-0.5*((sin(inrad-refrad)**2/(sin(inrad+refrad)**2))+
     &   (tan(inrad-refrad)**2/(tan(inrad+refrad)**2)))
      endif
C**** get VALUES OF C, WC
c     DENSITY (RHO) TIMES HEAT CAPACITY = RCSP. Assuming unit area of ground surface
c     WC = mass * specific heat (thermal capacitance) = density*volume*specific heat. Assuming unit area.
      if(runsnow.eq.0)then
       TT=T
       call soilprops(TT,ALTT,soilprop,moist)
       TT=T
       DENDAY(1)=Density(1)
       SPDAY(1)=Spheat(1)
       TKDAY(1)=Thconduct(1)
       RCSP = DENDAY(1)*SPDAY(1)
       WC(1)=RCSP*DEPP(2)/2.
       SOK =TKDAY(1)
       C(1)=SOK/DEPP(2)
       DO 51 I=2,N
        DENDAY(I)=Density(I)
        SPDAY(I)=Spheat(I)
        TKDAY(I)=Thconduct(I)
C       mass*specific heat product (per unit area)
        RCSP=DENDAY(I)*SPDAY(I)
        WC(I)=RCSP*(DEPP(I+1)-DEPP(I-1))/2.
        SOK=TKDAY(I)
        C(I)=SOK/(DEPP(I+1)-DEPP(I))
   51  CONTINUE
      else ! snow model
       maxsnode2=int(maxsnode1)
       TT=T
       call soilprops(TT,ALTT,soilprop,moist)
       TT=T
       DO 551 I=1,N
        DENDAY(I)=Density(I)
        SPDAY(I)=Spheat(I)
        TKDAY(I)=Thconduct(I)
        if(I.eq.1)then
         RCSP=DENDAY(1)*SPDAY(1)
         WC(1)=RCSP*DEPP(2)/2.
         SOK=TKDAY(1)
         if(DEPP(2).lt.1e-8)then
          C(1)=0.D0
         else
          C(1)=SOK/DEPP(2)
         endif
        else
         if((depp(i+1).gt.1e-8).or.(i.eq.18))then
C         mass*specific heat product (per unit area)
          RCSP = DENDAY(I)*SPDAY(I)
          if((depp(i).lt.1e-8).and.(maxsnode1.gt.0))then
           WC(I)=RCSP*DEPP(9-maxsnode2)
          else
           WC(I)=RCSP*(DEPP(I+1)-DEPP(I-1))/2.
          endif
          SOK =TKDAY(I)
          C(I)=SOK/(DEPP(I+1)-DEPP(I))
         else
C         mass*specific heat product (per unit area)
c         RCSP = DENDAY(I)*SPDAY(I)
          WC(I)=0.D0
          SOK =TKDAY(1)
          if(maxsnode2.eq.0)then
           C(I)=SOK/DEPP(10)
          else
           C(I)=SOK/DEPP(9-maxsnode2)
          endif
         endif
        endif
551    continue
      ENDIF
c     End of soil properties

      TAIR=TAB('TAR',TIME)
      ZENR=TAB('ZEN',TIME)
      SOLR=TAB('SOL',TIME)
      CLOUD=TAB('CLD',TIME)
      IRDOWN=TAB('IRD',TIME)

C     Modification by M. Kearney for effect of cloud cover on direct solar radiation, using the
C     Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
      IF ((CLOUD .GT. 0.).and.(HOURLY.eq.0)) THEN ! allowing for hourly cloud to be used for longwave calcs if longwave not provided
         SOLR = SOLR*(0.36+0.64*(1.-(CLOUD/100)))
      ENDIF
C     SOLAR RADIATION ON LEVEL GROUND. SURFACE ABSORPTIVITY, SAB, NOW CHANGING EACH MONTH/TIME INTERVAL, SABNEW
      if((int(grasshade).eq.1).and.(runsnow.eq.1))then
       methour=(int(SIOUT(1)/60)+1)+24*(doy-1)
       QSOLAR=SABNEW*SOLR*((100.- SHAYD)/100.)
       if(methour.gt.1)then
        if(snowhr(methour-1).gt.0)then
         QSOLAR=SABNEW*SOLR*((100.- 0.)/100.) ! no shade, snow on veg
        endif
       endif
      else
       QSOLAR=SABNEW*SOLR*((100.- SHAYD)/100.)
      endif
C     SOLAR MODIFICATION FOR SLOPES
      IF(SLOPE .GT. 0) THEN
C       SLOPING SURFACE
        CZ=COS(PI*ZENR/180.)
        ZSLR=TAB('ZSL',TIME)
        CZSL=COS(PI*ZSLR/180.)
        IF(ZENR .LT. 90.) THEN
C        CORRECT FOR SOLAR ON LEVEL GROUND FOR SOLAR ON A SLOPE
          QSOLAR=(QSOLAR/CZ)*CZSL
         ELSE
C        ZENITH ANGLE 90 DEGREES OR GREATER. NO DIRECT SOLAR.
        ENDIF
       ELSE
C      NO SLOPE. QSOLAR ALREADY CALCULATED, NO ADJUSTMENT FOR SLOPE NEEDED
      ENDIF

C     SWINBANK FORMULA FOR SKY TEMPERATURE
c      TSKY=0.0552*(TAIR+273.)**1.5-273.
C     PERCENT CLOUD COVER 2 DEG. LESS THAN TAIR

C     CONVERTING PERCENT CLOUDS TO FRACTION OF SKY CLEAR
      CLR=1.- (CLOUD/100.)

      if(IRDOWN.gt.0)THEN ! hourly IRdown provided
C      NET IR RADIATION: INCOMING FROM SKY + VEGETATION + HILLSHADE - OUTGOING FROM GROUND
       SRAD=SIGP*SLE*(T(1)+273.)**4
       HRAD=SIGP*SLEP*(TAIR+273.)**4
       QRADGR=((100.-SHAYD)/100.)*SRAD+(SHAYD/100.)*HRAD
       QRAD = IRDOWN - QRADGR
c      TSKY=((QRAD+QRADGR)/(SIGP))**(1./4.)-273
       TSKY=(IRDOWN/SIGP)**(1./4.)-273
      else
C      CLEAR SKY RADIANT TEMPERATURE
       if(IRmode.eq.0)then
c       Campbell and Norman 1998 eq. 10.10 to get emissivity of sky
        RH = TAB('REL',TIME)
        if(RH.gt.100.)then
         RH= 100.
        endif
        WB = 0.
        DP = 999.
C       BP CALCULATED FROM ALTITUDE USING THE STANDARD ATMOSPHERE
C       EQUATIONS FROM SUBROUTINE DRYAIR    (TRACY ET AL,1972)
        PSTD=101325.
        PATMOS=PSTD*((1.-(0.0065*ALTT/288.))**(1./0.190284))
        BP = PATMOS
        CALL WETAIR (TAIR,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,
     &      CP,WTRPOT)
      ARAD=1.72*((E/1000.)/(TAIR+273.16))**(1./7.)*0.0000000567*
     &(TAIR+273.16)**4*60./(4.185*10000.)
c     Below is the Gates formula (7.1)
c        ARAD=(1.22*0.00000005673*(TAIR+273.)**4-171)
       else
c       Swinbank, Eq. 10.11 in Campbell and Norman 1998
        ARAD=(0.0000092*(TAIR+273.16)**2)*0.0000000567*(TAIR+273.16)**4
     &  *60./(4.185*10000.)
       endif
C      APPROXIMATING CLOUD RADIANT TEMPERATURE AS REFERENCE SHADE TEMPERATURE - 2 degrees
       CRAD=SIGP*SLEP*(TAIR+271.)**4
c      Hillshade radiant temperature (approximating as air temperature)
       HRAD=SIGP*SLEP*(TAIR+273.)**4
C      GROUND SURFACE RADIATION TEMPERATURE
       SRAD=SIGP*SLE*(T(1)+273.)**4
C      TOTAL SKY IR AVAILABLE/UNIT AREA
       CLEAR = ARAD*CLR
       CLOD = CRAD*(CLOUD/100.)
c      previously SIGP*SLEP*(CLEAR + CLOUD)*((100.- SHAYD)/100.) but changed by MK to
c      allow the formula in Gates to be used
       if((int(grasshade).eq.1).and.(runsnow.eq.1))then
        methour=(int(SIOUT(1)/60)+1)+24*(doy-1)
        QRADSK=(CLEAR + CLOD)*((100.- SHAYD)/100.)
        if(methour.gt.1)then
         if(snowhr(methour-1).gt.0)then
          QRADSK=(CLEAR + CLOD)*((100.- 0.)/100.) ! no shade, snow on veg
         endif
        endif
       else
        QRADSK=(CLEAR + CLOD)*((100.- SHAYD)/100.)
       endif
C      VEGETATION IR AVAILABLE/UNIT AREA
c      previously QRADVG=SIGP*SLEP*(SHAYD/100.)*CRAD but changed by MK to allow formula
c      in Campbell to be used
       if((int(grasshade).eq.1).and.(runsnow.eq.1))then
        methour=(int(SIOUT(1)/60)+1)+24*(doy-1)
        QRADVG=(SHAYD/100.)*HRAD
C       GROUND SURFACE IR UPWARD/UNIT AREA
c       QRADGR=SIGP*SLEP*SRAD MK commented this out and replaced with below
        QRADGR=((100.-SHAYD)/100.)*SRAD+(SHAYD/100.)*HRAD
        if(methour.gt.1)then
         if(snowhr(methour-1).gt.0)then
          QRADVG=(0./100.)*HRAD ! no shade, snow on veg
          QRADGR=((100.-0.)/100.)*SRAD+(0./100.)*HRAD
         endif
        endif
       else
        QRADVG=(SHAYD/100.)*HRAD
        QRADGR=((100.-SHAYD)/100.)*SRAD+(SHAYD/100.)*HRAD
       endif

c      TOTAL HILLSHADE RADIATION
       QRADHL=HRAD
C      NET IR RADIATION: INCOMING FROM SKY + VEGETATION + HILLSHADE - OUTGOING FROM GROUND
       QRAD = (QRADSK + QRADVG)*VIEWF + QRADHL*(1-VIEWF) - QRADGR
c      TSKY=((QRAD+QRADGR)/(SIGP))**(1./4.)-273
       TSKY=(((QRADSK + QRADVG)*VIEWF + QRADHL*(1-VIEWF))/(SIGP))**
     & (1./4.)-273
c      TSKY=((QRADSK + QRADVG)/(SIGP))**(1./4.)-273
      endif

      if(runsnow.eq.1)then
c     change for snow - check what node is to be used for conduction, will be node 10 if no snow
       j=1
       do 1103 i=2,18
        if(density(i-1).lt.1e-8)then
         j=j+1
        endif
1103   continue
       methour=int((int(SIOUT(1)/60)+1)+24*(doy-1))
       if((doy.eq.1).and.(methour.eq.1))then
        QCOND=C(j)*(T(j+1)-T(1))
       else
        QCOND=C(j)*(T(j+1)-T(j))
       endif
      else
       QCOND=C(1)*(T(2)-T(1))
      endif

c     kludgy implementation of tide in - impose strong convective coupling
c     to substrate based on water temp
      VELR=TAB('VEL',TIME)
      methour=0
      tide=0
      if(TIME.ge.1440)then
       methour=(int(1380/60)+1)+24*(doy-1)
      else
       methour=(int(TIME/60)+1)+24*(doy-1)
      endif
      tide=tides(methour,1)
      if(tide.gt.0)then
          TAIR=tides(methour,2)
          TSKY=TAIR
          VELR=50000
      endif
C     COMPUTE VELOCITY AND TEMPERATURE PROFILES
      IF((ZH1.LE.0.000).AND.(ZH2.LE.0.000))THEN
C      NO SEGMENTED VELOCITY PROFILE (SINGLE LOG PROFILE)
        CALL MICRO(HGTP,RUFP,ZH,D0,TAIR,T(1),VELR,QCONV
     &  ,AMOL,NAIR,ZZ,VV,T,ZENR)
       ELSE
C      SEGMENTED VELOCITY PROFILE (VEGETATION OR OTHER OBJECTS MODIFYING VELOCITY PROFILE)
        CALL MICROSEGMT(HGTP,RUFP,TAIR,T(1),VELR,QCONV,AMOL,NAIR,
     &  ZZ,VV,T,ZENR)
c      QCOND=C(1)*(T(2)-T(1))
      ENDIF

C     SOIL TRANSIENTS; FIRST THE NODE AT THE SURFACE.
C     SIGN CONVENTION IN ALL EQUATIONS: POSITIVE TERMS = HEAT INPUT TO THE NODE, NEG. TERMS = HEAT LOSS FROM NODE
C     THIS SURFACE NODE EQUATION IS A HEAT BALANCE ON THE SOIL SURFACE NODE:
C     QIN = QOUT + QSTORED  REARRANGED TO GET THE RATE OF CHANGE OF TEMPERATURE TERM IN QSTORED, m*c*dT/dt
      if(tide.gt.0)then
        DTDT(1)=(QCOND+QCONV)/WC(1)
        QEVAP=0.D0
      else
       IF(PTWET.lt.1e-8)THEN
C       DRY SURFACE
         if(runsnow.eq.1)then
          DTDT(1)=(QSOLAR+QRAD+QCOND+QCONV+QFREZE)/WC(j)
         else
          DTDT(1)=(QSOLAR+QRAD+QCOND+QCONV)/WC(1)
         endif
         QEVAP=0.D0
        ELSE
C       SNOW or WET SURFACE
C       GETTING THE RELATIVE HUMIDITY FOR THIS POINT IN TIME
         RH = TAB('REL',TIME)
         if(RH.gt.100.)then
          RH = 100.
         endif
         WB = 0.D0
         DP = 999.
C        BP CALCULATED FROM ALTITUDE USING THE STANDARD ATMOSPHERE
C        EQUATIONS FROM SUBROUTINE DRYAIR    (TRACY ET AL,1972)
         PSTD=101325.
         PATMOS=PSTD*((1.-(0.0065*ALTT/288.))**(1./0.190284))
         BP = PATMOS
         CALL WETAIR (TAIR,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,
     &      CP,WTRPOT)
C       COMPUTING THE HEAT & MASS TRANSFER COEFFICIENTS, hc & hd
         IF(T(1).LT.-81.)THEN
           T(1) = -81.
         ENDIF
C       CHECK FOR OUTSIZE T(1)
         TTEST = ABS(T(1))
         IF(TTEST.GT. 100)THEN
          T(1) = TAIR + 1.0
         ENDIF
C       CHECKING FOR DIVIDE BY ZERO
         IF(T(1).EQ.TAIR)THEN
           T(1)= T(1)+0.1
         ENDIF
         HC = max(ABS((QCONV*4.184/60.*10000)/(T(1)-TAIR)),0.5D+0)
         HD = (HC/(CP*DENAIR))*(0.71/0.60)**0.666
         CALL EVAP(T(1),TAIR,RH,HD,QEVAP,SAT)
         if(runsnow.eq.1)then
          DTDT(1)=(QSOLAR+QRAD+QCOND+QCONV+QFREZE-QEVAP)/WC(j)
         else
          DTDT(1)=(QSOLAR+QRAD+QCOND+QCONV-QEVAP)/WC(1)
         endif
       ENDIF ! end check for wet surface
      endif ! end check for tide
C     SETTING UP THE DEEP SOIL TEMPERATURE, TDS, FOR SOIL TRANSIENTS.
C     (N=MM+1); N = # OF SOIL NODES SET UP IN 'DEP' ARRAY IN INPUT DATA
C      # OF SOIL NODES IS ASSUMED = 10, UNLESS SPECIFIED OTHERWISE BY A
C      'NON' PARAMETER LINE PLUS DATA LINE IN DATA INPUT. MAX. # NODES = 18.
C      ALL DEFAULT PARAMETER VALUES ARE IN BLKDATA IN PAR ARRAY.  THEIR NAMES
C      ARE IN SYMBOL ARRAY.
      T(N)=TAB('TDS',TIME)
      TDS=T(N)

C     CHECK FOR MAX # OF NODES
      IF (N .GT. 18) THEN
          WRITE(6,*)'MAX # OF NODES (18) EXCEEDED.'
          GO TO 200
         ELSE
          CONTINUE
      ENDIF

C     COMPUTE DERIVATIVES OF THE REST OF THE SOIL NODES
      DO 10 I=2,N
      if(runsnow.eq.1)then
       if((maxsnode1.gt.0).and.(I.eq.int((9-maxsnode1))))then
         DTDT(I)=(C(I-1)*(T(I-1)-T(I))+C(I)*(T(I+1)-T(I)))/WC(I) ! soil
       else
        if(WC(i-1).gt.0)then
         DTDT(I)=(C(I-1)*(T(I-1)-T(I))+C(I)*(T(I+1)-T(I)))/WC(I) ! snow or soil
        else
         dtdt(i)=dtdt(1) ! empty snow node - make it behave like surface
        endif
       endif
      else
       DTDT(I)=(C(I-1)*(T(I-1)-T(I))+C(I)*(T(I+1)-T(I)))/WC(I) ! below the surface
      endif

   10 continue
C     SET UP THE OUTPUT
      OUT(1)=TIME
      OUT(2)=TAIR
      OUT(3)=TSKY
      OUT(4)=T(1)
      OUT(5)=VELR
      OUT(6)=SOLR
      OUT(7)=CLOUD
      OUT(8)=QSOLAR
      OUT(9)=QRAD
      OUT(10)=QCOND
      OUT(11)=QCONV
      OUT(12)=AMOL
      OUT(13)=H
      OUT(53)=TIME/60.
      OUT(55)=TIME*60.
      OUT(101)=QEVAP
      I1=14
      I2=32
      DO 20 I=I1,I2
       if(runsnow.eq.1)then
        II=8 ! only output non snow profile
       else
       II=0
        endif
       OUT(I)=T(I+II-12)
   20 CONTINUE
      L=I1+N
      OUT(L)=TDS
      IF(NAIR.LE.0) RETURN
C     INSERTING AIR TEMPERATURES FROM THE GROUND UP, BUT NOT INCLUDING THE
C     SURFACE  (T(21-30))
      I1=33
      I2=42
      DO 22 I=I1,I2
       OUT(I)=T(I-12)
   22 CONTINUE
C     INSERTING VELOCITIES FROM THE GROUND UP, BUT NOT INCLUDING THE
C     SURFACE (V(1-10))
      I1=43
      I2=52
      DO 24 I=I1,I2
       OUT(I)=VV(I-42)
   24 CONTINUE

  200 RETURN
      END
