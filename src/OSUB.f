      SUBROUTINE OSUB(TIME,Y)

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

C     THIS SUBROUTINE IS CALLED FROM SFODE, THE NUMERICAL INTEGRATOR.  Y MAY BE EITHER
C     THE VALUE BEING ESTIMATED OR ITS DERIVATIVE.  HERE IT IS NOT NEEDED AND SO SERVES
C     ONLY AS A DUMMY VARIABLE FOR WHAT IS REALLY COMPUTED IN DSUB AND TRANSFERRED HERE
C     IN COMMON STATEMENTS.
      use commondat
C     VERSION 2 SEPT. 2000
      IMPLICIT NONE

      DOUBLE PRECISION ALTT,ZH,D0
      DOUBLE PRECISION C,CP,DENAIR,DEP,DEPP,DP,DTAU,E,ERR1,ESAT
      DOUBLE PRECISION FIN,H,PI
      DOUBLE PRECISION LASTIME,MAXSHD,OUT,lastsurf,prevsnow
      DOUBLE PRECISION OUT2,PTWET,SAB,SABNEW,SHAYD,SIOUT
      DOUBLE PRECISION RH,RHLOCL,RW,T,TANNUL,oldmoist
      DOUBLE PRECISION TAB,TD,TI,TIME,TIMEF,TVINC,TVIR
      DOUBLE PRECISION VD,WC,Y,SOK,END,DAS,oldcondep,densfun
      DOUBLE PRECISION MON,X,ZENR,inrad,refrad,snowcond,intercept
      DOUBLE PRECISION TSKY,ARAD,CRAD,CLOD,CLOUD,CLR,SOLR,QRADVG,QRADGR
      DOUBLE PRECISION RCSP,HGTP,RUFP,BEGP,PRTP,ERRP,snowout,curmoist
     & ,soiltemp,SLEP,NONP,SUN,PLT,rainfall,curhumid,surflux,IRDown
      DOUBLE PRECISION CLEAR,QRADSK,SIGP,TAIR,SRAD,QEVAP,frosttest,vel2m
      DOUBLE PRECISION htovpr,water,gwsurf,EP,zz,vv,curroot2
      DOUBLE PRECISION melthresh,time3,err

      DOUBLE PRECISION FROST,curpot,meanT,meanTpast,temp,SLE
      DOUBLE PRECISION bp,hrad,patmos,pstd,qrad,qradhl,viewf,wb,wtrpot
      DOUBLE PRECISION DENDAY,SPDAY,TKDAY,DENDAY2,SPDAY2,TKDAY2,time2
      DOUBLE PRECISION condep,rainmult,soilprop,moist,wccfinal,HTOFN
      DOUBLE PRECISION Z01,Z02,ZH1,ZH2,qconv,ttest,hc,hd,VELR,AMOL,wcc
      DOUBLE PRECISION curmoist2,curhumid2,curpot2,tt,tt_past,snowalbedo
      DOUBLE PRECISION DRLAM,DRRLAM,SRLAM,trans2,leafpot,curroot

      DOUBLE PRECISION PE,KS,BB,BD,maxpool,L,LAI,tide,minsnow,DD
      DOUBLE PRECISION snownode,maxsnode1,snode,daysincesnow,lastday,
     &undercatch,rainmeltf,rainfallb,snowtest,snowpres
      DOUBLE PRECISION PUNSH,ALAT,AMULT,PRESS,CMH2O,REFL,ALONC,TIMCOR
     & ,TSNHR,TSRHR,HEMIS,AZMUTH,SLOPE
      double precision rww,pc,rl,sp,r1,im
      double precision snowdens,snowmelt,snowtemp,cursnow,qphase,melt,
     & sumphase,sumphase2,snowage,prevden,qphase2,sumlayer,densrat,
     & cummelted,melted,snowfall,rainmelt,cpsnow,netsnow,hcap,meltheat,
     & layermass,xtrain,QFREZE,grasshade

      integer maxsnode2,maxsnode3,maxcount,js,numrun,rainhourly,hourly
      INTEGER I,IEND,IFINAL,ILOCT,IOUT,IPRINT,ITEST,trouble
      INTEGER J,JULNUM,MM,DOY,N,NAIR,ND,NOUT,dew,writecsv,runsnow
      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,slipped,sat
      INTEGER I91,I92,I93,I94,I95,I96,runmoist,evenrain,step,timestep
      INTEGER I97,I98,I99,I100,I101,errout,maxerr,errcount
      INTEGER IPINT,NOSCAT,IUV,IALT,IDAYST,IDA,IEP,ISTART

      INTEGER methour,IRmode,microdaily,runshade,k,lamb,cnd

      CHARACTER(3) SYMBOL,INAME,STP
      CHARACTER(6) NAME, HEAD
C      IOUT IS THE NUMBER OF THE OUTPUT TERM DESIRED AS NUMBERED IN NAME ARRAY
C      FILES 6,7,10 & 12 ARE CONSOLE, OUTPUT, METOUT & SOIL RESPECTIVELY
C      FILES 6,I2,I3 & I10 ARE CONSOLE, OUTPUT, METOUT & SOIL RESPECTIVELY

      DIMENSION DEPP(30),OUT2(55),NAME(55),curmoist(10),densfun(4)
      DIMENSION HEAD(55),snownode(10),snode(10),tt(30),y(30),tt_past(30)
      DIMENSION DRLAM(111),DRRLAM(111),SRLAM(111),qphase(10)

      DIMENSION curhumid(10),soiltemp(10),meanT(30),meanTpast(30)
      DIMENSION moist(10),soilprop(10,5),layermass(10)
      DIMENSION DENDAY(10),SPDAY(10),TKDAY(10),DENDAY2(10),qphase2(10),
     &    SPDAY2(10),TKDAY2(10),temp(83),oldmoist(10),curpot(18)
      DIMENSION curmoist2(18),curhumid2(18),curpot2(18),curroot2(18)
      DIMENSION PE(19),KS(19),BD(19),BB(19),L(19),curroot(18),DD(19)
      DIMENSION sumphase2(10)

      COMMON/TABLE/TI(211),TD(211)
      COMMON/TABLE2/ILOCT(21)
      COMMON/ARRAY/T(30),WC(20),C(20),DEP(30),OUT(101),IOUT(100),
     1 ITEST(23)
      COMMON/PAR/SIGP,RCSP,SOK,SAB,HGTP,RUFP,BEGP,MON,PRTP,ERRP,END,
     1 SLEP,DAS,NONP,SUN,PLT,FIN,STP
      COMMON/CARRAY/INAME(20),SYMBOL(23)
      COMMON/NONSCR/MM,N,X,TIMEF,DTAU,ERR1,H,NOUT,NAIR,IPRINT
C    SI UNITS FROM DSUB FOR OUTPUT TO DATA FILE FOR ANIMAL ENERGY
C    BALANCES
      COMMON/SIUNIT/SIOUT(10)
      COMMON/NDAY/ND
      COMMON/WDSUB/TSKY,ARAD,CRAD,CLOUD,CLR,SOLR
      COMMON/WOSUB/DEPP
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101
      COMMON/DAYS/TANNUL
      COMMON/DAYJUL/JULNUM,DOY
      COMMON/SNOWPRED/snowtemp,snowdens,snowmelt,snownode,minsnow
     &,maxsnode1,snode,cursnow,daysincesnow,lastday,undercatch,rainmeltf
     &,densfun,snowcond,intercept,snowage,prevden,grasshade
C     PERCENT GROUND SHADE & ELEVATION (M) TO METOUT
      COMMON/GROUND/SHAYD,ALTT,MAXSHD,SABNEW,PTWET,rainfall
      COMMON/LOCLHUM/RHLOCL
      COMMON/VIEWFACT/VIEWF
      COMMON/SOYFILS/DENDAY,SPDAY,TKDAY
      common/prevtime/lastime,temp,lastsurf
      common/prevtime2/slipped
      common/soilmoist/condep,rainmult,maxpool
      common/soilmoist3/runmoist,evenrain,maxcount
      common/soimoist2/rww,pc,rl,sp,r1,im
      COMMON/DAILY/microdaily
      COMMON/NICHEMAPRIO/SLE,ERR,soilprop,surflux
      common/moistcom/moist,ep
      COMMON/CMYCRO/ZH,D0
      COMMON/DMYCRO/Z01,Z02,ZH1,ZH2
      COMMON/AIRRAY/ZZ(10),VV(10)
      common/campbell/PE,KS,BB,BD,L,LAI,DD
      common/shaderun/runshade
      common/curmoist/curmoist2
      common/write/writecsv
      common/snowmod/runsnow,trouble
      COMMON/temps/TT,TT_past,cummelted
      common/hour/hourly,rainhourly
      COMMON/LAMBDA/DRLAM,DRRLAM,SRLAM,LAMB
      COMMON/IR/IRmode
      COMMON/melt/QFREZE,xtrain,qphase,sumphase,sumphase2
      COMMON/WIOCONS/PUNSH,ALAT,AMULT,PRESS,CMH2O,REFL,ALONC,TIMCOR,
     * AZMUTH,SLOPE,TSNHR,TSRHR,Hemis
      COMMON/WIOCONS2/IPINT,NOSCAT,IUV,IALT,IDAYST,IDA,IEP,ISTART,IEND
      common/errormsg/errout,maxerr,errcount
      COMMON/WICHDAY/NUMRUN

      DATA NAME/'TIME ','TAIR','TSKY','TSURF','VEL','SOL  ','TLIZ',
     1 'QSOLAR','QRAD','QCOND','QCONV','MOL ','STEP','T2','T3','T4'
     1,'T5','T6',
     2 'T7','T8','T9','T10','T11','T12','T13','T14','T15','T16','T17',
     3 'T18','T19','T20','TA1','TA2','TA3','TA4','TA5','TA6','TA7'
     4,'TA8','TA9','TA10','VA1','VA2','VA3','VA4','VA5','VAL','VA7'
     5,'VA8','VA9','VA10','T/60  ','M-E','T*60'/
      DATA IFINAL/0/
      PI=3.14159
      if((time.gt.0).and.(int(time).eq.int(lastime)))then
          trouble=trouble+1
      endif
      
C     INTITIALISE
      maxsnode3=0
      rainmelt=0.D0
      snowout=0.D0
      snowpres=0.D0
      methour=0
      if(runsnow.eq.0)then
       maxsnode1=0.D0
       maxsnode2=0
      endif

c     OUTPUT
C     1        TIME                                              TIME
C     2        TAIR(INPUT AT REF. HT.)                           TAR
C     3        TSKY                                              TSKY
C     4        TSURFACE                                          TSURF
C     5        V AIR (INPUT AT REF. HT.)                         VEL
C     6        SOLAR FLUX                                        SOL
C     7        PERCENT CLOUD COVER                               CLOUD
C     8        ABSORBED SOLAR HEAT FLUX (SOIL)                   QSOLAR
C     9        RADIATION BETWEEN SURFACE AND SKY                 QRAD
C     10       CONDUCTION FROM SURFACE TO NEXT NODE              QCOND
C     11       CONVECTION FROM SURFACE TO AIR                    QCONV
C     12       MONIN-OBUKOV LENGTH                               MOL
C     13       INTEGRATION STEP SIZE                             STEP
C     14 TO 32 SOIL TEMP'S AT DEP(I): FROM 1ST BELOW SURF DOWN   TI
C     33 TO 42 AIR TEMP'S AT DEP(I): FROM 1ST ABOVE SURF UPWARD  TAI
C     43 TO 52 AIR VEL. AT DEP(I): FROM 1ST ABOVE SURF UPWARD    VAI
C     53       TIME/60.                                          T/60
C     54       DEW POINT TEMP.                                   TDP
C     55       TIME*60.                                          T*60
C
C     DAY COUNTER TO DELETE OUTPUT FOR REPLICATE DAYS
C     NEEDED TO ESTABLISH SOIL STEADY PERIODICS
      HTOFN=333.500 !J/g
      melted=0.D0
      melthresh=0.4
      time2=0.D0
      time3=0.D0
      melt=0.D0
      layermass(:)=0.D0
      qphase(:)=0.D0
      qphase2(:)=0.D0
      curhumid(:)=0.D0
      curmoist(:)=0.D0
      curpot(:)=0.D0
      curroot(:)=0.D0
      denday2(:)=0.D0
      meant(:)=0.D0
      meantpast(:)=0.D0
      oldmoist(:)=0.D0
      out2(:)=0.D0
      soiltemp(:)=0.D0
      spday2(:)=0.D0
      tkday2(:)=0.D0
      do 876 i=1,25
       if(time.ge.ti(i+11))then
        time2=ti(i+11)
        time3=time2-60.
       endif
876   continue
      if(time.ne.time2)then
          time=time2
      endif
      if(runsnow.eq.1)then
c     check if no snow and make temp equal top temp
       do 555 i=1,8
        if(snode(i).lt.1e-8)then
         tt(i)=tt(1)
         t(i)=tt(1)
         if((snode(8).lt.1e-8).and.(cursnow.lt.1e-8))then
          tt(9)=tt(1)
          t(9)=t(1)
         endif
        endif
  555  continue
      endif
      IF(TIME .GT. 0.0) GO TO 5
       IFINAL=IFINAL+1
       IF(IFINAL .EQ. 1)THEN
C       SETTING SNOW VARIABLE
C       SETTING THIS MONTH'S SURFACE REFLECTIVITY FROM THE ABSORPTIVITIES
        SABNEW = 1.0 - REFLS(DOY)
C       SETTING THIS MONTH'S PERCENT OF SURFACE WITH FREE WATER/SNOW ON IT
        PTWET = PCTWET(DOY)
        rainfall=RAIN(DOY)
        CONTINUE
       ENDIF
    5 CONTINUE

      if(runsnow.eq.1)then
       PTWET=PCTWET(DOY)
       rainfall=RAIN(DOY)
      endif
      if(int(rainhourly).eq.1)then
       methour=0
       methour=(int(SIOUT(1)/60)+1)+24*(DOY-1)
       rainfall=rainhr(int(TIME/60.+1+25*(DOY-1)))
      endif
      if(microdaily.eq.1)then
       if(DOY.gt.1)then
        ND=1
       endif
      endif


C     NO OUTPUT IF REPEATED DAY,IFINAL, NOT THE FINAL ITERATION
      IF (IFINAL .LT. ND) GO TO 200

C     OUTPUT FOR START OF A DAY
      IF ((TIME-BEGP).GT.0.) GO TO 20
C     SET UP OUTPUT FOR C, WC
      IF (IPRINT.NE.1) GO TO 19
 19   CONTINUE
      DO 40 I=1,NOUT
       J=IABS(IOUT(I))
       HEAD(I)=NAME(J)
   40 CONTINUE

   20 CONTINUE
      DO 50 I=1,NOUT
       J=IABS(IOUT(I))
       OUT2(I)=OUT(J)
   50 CONTINUE

      IF (TIME .LT. 1440.) GO TO 150
      IFINAL=0

  150 CONTINUE
C     PUTTING SOIL DEPTH VALUES AT NODES INTO ARRAY DEPP
      DO 74 I=1,N
       J=I+NAIR
       if(runsnow.eq.1)then
        DEPP(I+8)=DEP(J)
       else
        DEPP(I)=DEP(J)
       endif
74    CONTINUE

      IEND = 14 + N - 2

C     SETTING UP THE 'OUTPUT'
      TAIR=TAB('TAR',TIME)
      ZENR=TAB('ZEN',TIME)
      SOLR=TAB('SOL',TIME)
      CLOUD=TAB('CLD',TIME)
      T(N)=TAB('TDS',TIME)
      VEL2M=TAB('VEL',TIME)/(100.*60.)
      IRDOWN=TAB('IRD',TIME)

c     phase change for freezing moist soil
      methour=int((TIME/60+1)+24*(DOY-1))
      if((methour.eq.1).and.(DOY.eq.1))then
        meanT=tt
        meanTpast=tt_past
      else
       if(runsnow.eq.1)then
        js=9
       else
        js=1
       endif
       do 1131 j=js,js+9 ! loop through soil nodes
        if(j.lt.js+9)then
         meanT(j)=(tt(j)+tt(j+1))/2. ! current temp
         meanTpast(j)=(tt_past(j)+tt_past(j+1))/2. ! last hour's temp
        else
         meanT(j)=tt(j) ! current temp
         meanTpast(j)=tt_past(j) ! last hour's temp
        endif
        if((meanTpast(j).gt.0).and.(meanT(j).le.0))then ! phase change, freezing
         if(j.lt.js+9)then
          layermass(j-js+1)=(dep(j-js+1+4)-dep(j-js+4))*10000*
     &     moist(j-js+1)
         else
          layermass(j-js+1)=(dep(j-js+4)+100-dep(j-js+4))*10000*
     &     moist(j-js+1)! deep soil
         endif
         qphase2(j-js+1)=(meanTpast(j)-meanT(j))*layermass(j-js+1)*4.186
         sumphase2(j-js+1)=qphase2(j-js+1)+sumphase2(j-js+1)
         if(sumphase2(j-js+1).gt.(HTOFN*layermass(j-js+1)))then
          t(j)=-0.1
          t(j+1)=-0.1
          tt(j)=-0.1
          tt(j+1)=-0.1
          y(j)=-0.1
          y(j+1)=-0.1
          sumphase2(j-js+1)=0.
         else
          t(j)=0.1
          t(j+1)=0.1
          tt(j)=0.1
          tt(j+1)=0.1
          y(j)=0.1
          y(j+1)=0.1
         endif
        endif
1131   continue
      endif


C     Modification by M. Kearney for effect of cloud cover on direct solar radiation, using the
C     Angstrom formula (formula 5.33 on P. 177 of "Climate Data and Resources" by Edward Linacre 1992
      IF ((CLOUD .GT. 0.).and.(HOURLY.eq.0)) THEN ! allowing for hourly cloud to be used for longwave calcs if longwave not provided
       SOLR = SOLR*(0.36+0.64*(1.-(CLOUD/100)))
      ENDIF

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
       CLR=1.- (CLOUD/100.)
C      CLEAR SKY RADIANT TEMPERATURE
       if(int(IRmode).eq.0)then
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
     &  (TAIR+273.16)**4*60./(4.185*10000.)
c       Below is the Gates formula (7.1)
c       ARAD=(1.22*0.00000005673*(TAIR+273.)**4-171)
       else
c       Swinbank, Eq. 10.11 in Campbell and Norman 1998
        ARAD=(0.0000092*(TAIR+273.16)**2)*0.0000000567*(TAIR+273.16)**4
     & *60./(4.185*10000.)
       endif

C      APPROXIMATING CLOUD RADIANT TEMPERATURE AS REFERENCE SHADE TEMPERATURE - 2 degrees
       CRAD=SIGP*SLEP*(TAIR+271.)**4
c      Hillshade radiant temperature (approximating as air temperature)
       HRAD=SIGP*SLEP*(TAIR+273.)**4
C      GROUND SURFACE RADIATION TEMPERATURE
       SRAD=SIGP*SLE*(TT(1)+273.)**4
C      TOTAL SKY IR AVAILABLE/UNIT AREA
       CLEAR =  ARAD*CLR
       CLOD =  CRAD*(CLOUD/100.)
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
       TSKY=(((QRADSK + QRADVG)*VIEWF + QRADHL*(1-VIEWF))/(SIGP))**(1./
     & 4.)-273
      ENDIF

c     TSKY=((QRADSK + QRADVG)/(SIGP))**(1./4.)-273
      SIOUT(1) = TIME
C     AIR TEMPERATURE AT ANIMAL HEIGHT (1ST NODE BELOW REFERENCE HEIGHT, 200 CM)
      SIOUT(2) = OUT(34)
C     RELATIVE HUMIDITY AT ANIMAL HEIGHT
      SIOUT(3) = TAB('REL',TIME)
      IF(SIOUT(3).gt.100)then
       SIOUT(3)=100.
      ENDIF
C     VELOCITY AT ANIMAL HEIGHT (OUT(43 - 52) ARE NODES DOWN FROM REFERENCE HEIGHT AIR NODE
C     THE 3 DEFAULT NODE HEIGHTS DEFINED IN SUB. IOMET2. OUT(43) IS REFERENCE HEIGHT.
      SIOUT(4) = OUT(44)/(100.*60.)
C     SUBSTRATE TEMPERATURE
      SIOUT(5) = TT(1)
C     1ST SOIL NODE BELOW THE SURFACE (minsnow)
      SIOUT(6) = TT(2)
C     DEEP SOIL NODE
      SIOUT(7) = TT(N)
C     SOLAR & ZENITH OUTPUT TO METOUT FOR A HORIZONTAL SURFACE
C     FOR ANIMAL CALCULATIONS, ALTHOUGH SLOPE CORRECTIONS DONE
C     IN DSUB FOR THE GROUND
      SIOUT(8) = ZENR
C     CONVERTING FROM CAL/CM2-MIN TO W/M2
      SIOUT(9) = SOLR * 4.185 * 10000. / 60.
C     BLACK BODY EQUIVALENT SKY INFRARED RADIATION (W/M2)
c     SIOUT(10) = ((QRADSK+QRADVG)/SIGP)**0.25 - 273.15
      SIOUT(10) = TSKY
C     FROST

c     convert to W/m2
      QEVAP = OUT(101)*4.185*10000./60.
      if(runsnow.eq.1)then
       methour=(int(SIOUT(1)/60)+1)+24*(DOY-1)
      if(densfun(1).gt.0)then
        if(densfun(2).gt.0)then ! exponential function
         snowdens=(densfun(1)-densfun(2))*(1-EXP(-1*densfun(3)*
     &   cursnow-densfun(4)*snowage))+densfun(2)
        else ! linear function
         snowdens=min(0.9167D+0,densfun(1)*snowage+densfun(2))
        endif
        densrat=snowdens/prevden
        if(methour.gt.1)then
        snowhr(methour-1)=snowhr(methour-1)/densrat
        endif
       endif
      if((OUT(2).le.snowtemp).and.(rainfall.gt.0.0))then
c       compute snow fall using conversion from daily rain to daily snow (disaggregated over 24 hours) and convert from mm rain to cm snow
        if((time.lt.1e-8).or.(int(rainhourly).eq.1))then
c        account for undercatch
         snowfall=((rainfall*rainmult*undercatch*0.1)/snowdens)! snowfall in cm/m2
         if(SHAYD.gt.0)then
           snowfall=snowfall*(1-intercept)! including interception
         endif
         daysincesnow=0.3 !resets daysincesnow and ensures that the equation for snow albedo gives 90% reflectance immediately upon snowfall
        else
         snowfall=0.D0
         daysincesnow=daysincesnow+1./25.
        endif
       else
        snowfall=0.D0
        daysincesnow=daysincesnow+1./25.
       endif
       if(daysincesnow.lt.0.3)then
        daysincesnow=0.3 !ensures that the equation for snow albedo gives 90% reflectance immediately upon snowfall
       endif
      else
       snowfall=0.D0
      endif
      if(out(4).gt.0)then
       HTOVPR=2500.8-2.36*out(4)+0.0016*out(4)**2-0.00006*out(4)**3
      else
       HTOVPR=2834.1-0.29*out(4)-0.004*out(4)**2
      endif
      HTOVPR=HTOVPR*1000 !convert from J/g to J/kg
      WATER = QEVAP/HTOVPR ! kg/s/m2
C     kg/s/m2 TO g/h/m2
      GWSURF  = WATER * 1000. * 3600.
      if(gwsurf.lt.0)then
      gwsurf=0
      endif
c     don't lose water if heat is just going into melting snow
      if((gwsurf.lt.0).or.(out(4).le.0))then
      gwsurf=0
      endif
      methour=0
      methour=(int(SIOUT(1)/60)+1)+24*(DOY-1)
************** begin main code for snow model ********************************
      if(runsnow.eq.1)then
c       write(*,*) "in osub"
c       write(*,*) methour
c      get cm snow lost due to rainfall - from Anderson model
       if((OUT(2).ge.snowtemp).and.(rainfall.gt.0).and.((time.lt.1e-8)
     & .or.(int(rainhourly).eq.1)))then ! mean air temperature warmer than snow temperature and more than 1 mm rain and midnight
        rainmelt=(rainmeltf*rainfall*(OUT(2)))*24./10./snowdens ! melt the snow as a function of the air temperature
        if(int(rainhourly).eq.1)then
         rainmelt=rainmelt/24.
        endif
        if(rainmelt.lt.0)then
         rainmelt=0.D0
        endif
       else
        rainmelt=0.D0
       endif
       netsnow=snowfall-gwsurf*0.0001/snowdens ! convert gwsurf from g/h/m2 (cm3/h/m2) to cm/m2 of evap, and then to snow depth equivalent (1mm/m2 rain/evap = 1000 cm3/m2)
       if(netsnow.lt.0)then
        netsnow=0
       endif
       netsnow = netsnow + xtrain/snowdens ! add any extra snow carried over by xtrain
       xtrain=0.D0 !reset xtrain
       maxsnode2=int(maxsnode1)
       melt=0.D0
       if((methour.eq.1).and.(DOY.eq.1))then
        meanT=tt
       else
        do 98 i=1,8 ! get mean snowpack temperature
         meanT(i)=(tt(i)+tt(i+1))/2. ! current temp
         meanTpast(i)=(tt_past(i)+tt_past(i+1))/2. ! last hour's temp
98      continue
       endif
       if((methour.eq.1).and.(DOY.eq.1))then
        melt=0.D0
       else
        prevsnow=snowhr(methour-1)
        if(prevsnow.ge.minsnow)then ! melt the snow
        WB = 0.D0
        DP = 999.
C       BP CALCULATED FROM ALTITUDE USING THE STANDARD ATMOSPHERE
C       EQUATIONS FROM SUBROUTINE DRYAIR    (TRACY ET AL,1972)
        PSTD=101325.
        PATMOS=PSTD*((1.-(0.0065*ALTT/288.))**(1./0.190284))
        BP = PATMOS
        call WETAIR(0.D+0,WB,100.D+0,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,
     &  DENAIR,CP,WTRPOT) ! get specific heat and mixing ratio of humid air at zero C
        cpsnow = (2.100*snowdens+(1.005+1.82*(RW/1.+RW))* ! based on https://en.wiktionary.org/wiki/humid_heat
     &   (1-snowdens)) ! compute weighted specific heat accounting for ice vs airm SI units
         hcap=cpsnow+HTOFN !J/gC
         do 1100 j=1,(maxsnode2+1) ! loop through snow nodes
          cnd=j+7-maxsnode2
          if(j.ne.maxsnode2+1)then ! not at the bottom yet
           if(meanT(cnd).gt.melthresh)then ! check if greater than threshold above 0 deg C (here 0.4 deg C)
            meltheat=(meanT(cnd)-melthresh)*cpsnow* ! g snow melted (change in temp x heat capacity x g snow (1m2 * height of layer)
     &       (snode(max(1,j+8-maxsnode2))-snode(max(1,cnd)))*10000.D0* ! plus mass times heat of fusion since the layer got above zero
     &        snowdens+HTOFN*(snode(max(1,j+8-maxsnode2))-
     &        snode(max(1,cnd)))*10000.D0*snowdens
            melt=melt+meltheat/HTOFN !(snow depth in cm x snow density in g/cm3 converted units of g/m2) gives J heat input,
           endif                     ! divided by heat of fusion to get g melted
          else
           if(meanT(cnd).gt.melthresh)then ! doing deepest node - depth here is flexible and is the difference between current snow and the
            meltheat=(meanT(cnd)-melthresh)*cpsnow* ! height of the current maximum node used of the 8 possible
     &       (cursnow-snownode(max(1,maxsnode2)))*10000.D0*snowdens+
     &      HTOFN*(cursnow-snownode(max(1,maxsnode2)))*10000.D0*snowdens
            melt=melt+meltheat/HTOFN
           endif
          endif
          if(cursnow.lt.200)then
           if(j.lt.8)then
            if(j.ne.maxsnode2+1)then ! not at bottom yet
             if((meanTpast(cnd).ge.0).and.(meanT(cnd).le.0))then ! phase change, freezing
              layermass(cnd)=max(0.D+0,(snode(max(1,j+8-maxsnode2))-
     &        snode(max(1,cnd)))*10000.D0*snowdens)
              qphase(max(1,cnd))=(meanTpast(max(1,cnd))-meanT(max(1,cnd)
     &        ))*layermass(max(1,cnd))*cpsnow
              t(cnd)=0.D0
              t(cnd+1)=0.D0
              tt(cnd)=0.D0
              tt(cnd+1)=0.D0
              y(cnd)=0.D0
              y(cnd+1)=0.D0
             endif
            else
             if((meanTpast(max(1,cnd)).ge.0.D0).and.
     &(meanT(max(1,cnd)).le.0.))then ! phase change, freezing        
              layermass(cnd)=max(0.D+0,(cursnow-snownode(max(1,maxsnode2
     &        )))*10000*snowdens)
              qphase(max(1,cnd))=(meanTpast(max(1,cnd))-meanT(max(1,cnd
     &        )))*layermass(max(1,cnd))*cpsnow
              t(max(1,cnd))=0.D0
              t(cnd+1)=0.D0
              tt(max(1,cnd))=0.D0
              tt(cnd+1)=0.D0
              y(max(1,cnd))=0.D0
              y(cnd+1)=0.D0
             endif
            endif
           else
            if((meanTpast(8).ge.0.D0).and.(meanT(8).le.0.D0))then ! phase change, freezing
             layermass(8)=max(0.D+0,(cursnow-snownode(8))*10000.D0
     &       *snowdens)
             qphase(8)=(meanTpast(8)-meanT(8))*layermass(8)*
     &       cpsnow
             t(max(1,cnd))=0.D0
             tt(max(1,cnd))=0.D0
             y(max(1,cnd))=0.D0
            endif
           endif
          endif
1100     continue
         sumphase=sum(qphase)+sumphase
         sumlayer=sum(layermass)
         if(cursnow.lt.200)then
          if(sumphase.gt.(HTOFN*sumlayer))then
           do 1121 j=1,(maxsnode2+1) ! loop through snow nodes
            if(j.lt.9)then
             cnd=j+7-maxsnode2
             if((t(max(1,cnd)).lt.-1e-8).and.(qphase(max(1,cnd))
     &        .gt.0.D0))then
              t(max(1,cnd))=-0.5
              t(cnd+1)=-0.5
              tt(max(1,cnd))=-0.5
              tt(cnd+1)=-0.5
              y(max(1,cnd))=-0.5
              y(cnd+1)=-0.5
              qphase(max(1,cnd))=0.
             endif
            endif
1121       continue
           sumphase=0.D0
          endif
         !write(*,*) sumphase
         !write(*,*) snode
         endif
        endif
       endif
       melted=melt*0.0001 ! convert from g/m2 of ice to cm snow
       if(DOY.gt.1)then
        melted=melted*snowmelt
        QFREZE=melted*(1-snowmelt)*79.7/60./10000.D0 ! convert from g snow melted to cal/min/cm2, latent heat of fusion being 79.7 calories per gram, only affecting surface layer!
        if(snowmelt.gt.1)then
         QFREZE=0.D0
        endif
       endif
       cummelted=melted+cummelted
       if((methour.eq.1).and.(DOY.eq.1))then
        melted=0
       endif
       if(methour.gt.1)then
        cursnow=snowhr(methour-1)+netsnow-rainmelt-melted
        snowhr(methour)=cursnow
        if(snowhr(methour).lt.0.D0)then
         snowhr(methour)=0.D0
        endif
        snowpres = snowhr(methour)
        snowtest = snowhr(methour-1)
        if(int(siout(1)).eq.1380)then
         if(cummelted.gt.snowhr(methour-1))then
          cummelted=snowhr(methour-1)
         endif
         cummelted=0.D0
        endif
       endif
       if(snowhr(methour).lt.0.D0)then
        snowhr(methour)=0.D0
       endif
       if(methour.gt.1)then
       if((snowhr(methour).lt.minsnow).and.(snowhr(methour-1).lt.1e-8)
     &)then
        xtrain=snowhr(methour)*snowdens !save snow below min level and add to next hour
        snowhr(methour)=0.D0
       endif
       endif
       cursnow=snowhr(methour)
C      make sure snow is 0 degrees or less, now that melting is done
       do 102 i=1,8
        if(snode(i).gt.0)then
         if(tt(i).gt.0.D0)then
          t(i)=0.D0
          tt(i)=0.D0
          y(i)=0.D0
          if((i.eq.8).and.(tt(9).gt.0.D0))then
           t(9)=0.D0
           tt(9)=0.D0
           y(9)=0.D0
          endif
         endif
        endif
102    continue
c      ensure that recently fallen snow is frozen
       if(methour.gt.1)then
        if((snowhr(methour).gt.0).and.(snowhr(methour-1).lt.1e-8))then
         do 101 i=1,8
          if(snowhr(methour).gt.snownode(i))then
           maxsnode3=i
          endif
101      continue
         if(maxsnode3.eq.0)then
          if(tt(9).gt.0.D0)then
           t(9)=0.D0
           tt(9)=0.D0
           y(9)=0.D0
          endif
         else
          do 989 i=9-maxsnode3,8-maxsnode3
           if(tt(i).gt.0.D0)then
            t(i)=0.D0
            tt(i)=0.D0
            y(i)=0.D0
           endif
989       continue
         endif
        endif
       else
        cursnow=netsnow-rainmelt-melted
        snowhr(methour)=cursnow
        SNOWOUT=snowhr(methour)
        tt_past=tt
        call snowlayer
       endif
       if(methour.gt.1)then
        if(snowhr(methour-1).gt.0.D0)then
         snowalbedo=(-9.8740*dlog(daysincesnow) + 78.3434)/100.
         SABNEW = 1-snowalbedo
C        SETTING THIS MONTH'S PERCENT OF SURFACE WITH FREE WATER/SNOW ON IT
         PTWET = 100.
         SLE = 0.98
        else
         SABNEW = 1.0 - REFLS(DOY)
C        SETTING THIS MONTH'S PERCENT OF SURFACE WITH FREE WATER/SNOW ON IT
         PTWET = PCTWET(DOY)
         SLE = SLES(DOY)
        endif
       else
        SABNEW = 1.0 - REFLS(DOY)
C       SETTING THIS MONTH'S PERCENT OF SURFACE WITH FREE WATER/SNOW ON IT
        PTWET = PCTWET(DOY)
        SLE = SLES(DOY)
       endif
       SNOWOUT=snowhr(methour)
       call snowlayer
       if(methour.gt.1)then
        if((SNOWOUT-snownode(max(1,int(maxsnode1))).lt.0.5) ! avoid getting too close to a given node - save to xtrain and add/subtract later
     &   .and.(snowout.le.snowhr(methour-1)).and.(SNOWOUT.gt.0))then
         xtrain=xtrain+(SNOWOUT-snownode(max(1,int(maxsnode1))))
     &    *snowdens
         SNOWOUT=snownode(max(1,int(maxsnode1)))-.1
         snowhr(methour)=SNOWOUT
        endif
       endif
       call snowlayer
       if(((snowout.lt.minsnow+0.5).or.(snowout.lt.snownode(max(1,int ! avoid getting too close to a given node - save to xtrain and add/subtract later
     &  (maxsnode1)))+0.5)).and.(snowout.gt.0))THEN
        xtrain=0.5*snowdens
        SNOWOUT=snownode(max(1,int(maxsnode1)))+0.5
        snowhr(methour)=SNOWOUT
       endif
       call snowlayer
       tt_past=tt
      else
       xtrain=0.D0
       SNOWOUT=0.D0
      endif
      ! count days of snow in a row
      if((methour.eq.1).and.(DOY.eq.1))then
       snowage=0.D0
      else
       if(snowhr(methour).gt.0.D0)then
        if(snowhr(methour-1).gt.0.D0)then
         snowage=snowage+1./25.
        else
         snowage=0.D0
        endif
       else
        snowage=0.D0
       endif
      endif
      prevden=snowdens
********************* end main code for snow model *******************************


      t=tt
*************** soil moisture model ***********************
      if(runmoist.eq.1)then
       rainfallb=rainfall
       if((DOY.eq.1).and.(methour.eq.1))then
        curmoist=moists(1:10,1)
        j=1
        do 221 i=1,18
         if(mod(i,2).ne.0)then
          curmoist2(i)=curmoist(j)
          j=j+1
         else
          curmoist2(i)=curmoist2(i-1)
         endif
221     continue
       else
        curmoist=soilmoist(methour-1,3:12)
       endif
c      choosing between even rainfall through the day or one event at midnight
       if((snowpres.gt.0).and.(OUT(2).lt.snowtemp))then
        condep=condep+melted
       else
        if(int(rainhourly).eq.0)then
         if(evenrain.eq.0)then
          if((time.gt.1e-8).or.(cursnow.gt.0))then
           rainfallb=0.D0
          endif
          condep=condep+rainfallb*rainmult+melted
         else
          condep=condep+rainfallb/24.*rainmult+melted
         endif
        else
         condep=condep+rainfallb*rainmult+melted
        endif
       endif
       if(condep.lt.0.D0)then
        condep=0.D0
       endif
       soiltemp(1)=OUT(4)
       soiltemp(2:10)=OUT(14:22)
       if(snowout.lt.1e-8)then
c       now compute potential evaporation, EP
        VELR=TAB('VEL',TIME)
C       COMPUTE VELOCITY AND TEMPERATURE PROFILES
        IF((ZH1.LE.0.000).AND.(ZH2.LE.0.000))THEN
C        NO SEGMENTED VELOCITY PROFILE (SINGLE LOG PROFILE)
         CALL MICRO(HGTP,RUFP,ZH,D0,TAIR,soiltemp(1),VELR,QCONV,AMOL,
     &NAIR,ZZ,VV,T,ZENR)
        ELSE
C        SEGMENTED VELOCITY PROFILE (VEGETATION OR OTHER OBJECTS MODIFYING VELOCITY PROFILE)
         CALL MICROSEGMT(HGTP,RUFP,TAIR,soiltemp(1),VELR,QCONV,AMOL,NAIR
     &,ZZ,VV,T,ZENR)
        ENDIF
C       GETTING THE RELATIVE HUMIDITY FOR THIS POINT IN TIME
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
C      COMPUTING THE HEAT & MASS TRANSFER COEFFICIENTS, hc & hd
        IF(soiltemp(1).LT.-81.)THEN
         soiltemp(1)=-81.
        ENDIF
C       CHECKING FOR DIVIDE BY ZERO
        IF(soiltemp(1).EQ.TAIR)THEN
         soiltemp(1)=soiltemp(1)+0.1
        ENDIF
C       CHECKING FOR DIVIDE BY ZERO
        IF(T(1).EQ.TAIR)THEN
         T(1)=T(1)+0.1
        ENDIF
C       CHECK FOR OUTSIZE T(1)
        TTEST=ABS(soiltemp(1))
        IF(TTEST.GT. 100)THEN
         T(1)=TAIR + 1.0
        ENDIF
        HC=max(ABS((QCONV*4.184/60.*10000.D0)/(T(1)-TAIR)),0.5D+0)
        HD=(HC/(CP*DENAIR))*(0.71/0.60)**0.666
        sat=2 ! setting for evap to simulate wet surface and return SI units
        CALL EVAP(soiltemp(1),TAIR,RH,HD,QEVAP,SAT)
        if(soiltemp(1).gt.0)then
         HTOVPR=2500.8-2.36*soiltemp(1)+0.0016*soiltemp(1)**2-0.00006
     &*soiltemp(1)**3
        else
         HTOVPR=2834.1-0.29*soiltemp(1)-0.004*soiltemp(1)**2
        endif
        HTOVPR=HTOVPR*1000.D0
c       evaporation potential, mm/s (kg/s)
        EP = QEVAP/HTOVPR
        if(EP.le.0.D0)then
         EP=0.0000001
        endif
        SAT=1
        CALL EVAP(soiltemp(1),TAIR,RH,HD,QEVAP,SAT)
       else
        EP=0.0000001
       endif ! end check for snow cover - no evap if there is
       if((condep.gt.0.D0).or.((snowout.gt.0.D0).and.(soiltemp(1).gt.
     &  0.D0)))then
        curmoist(1)=1.-BD(1)/DD(1)
        curmoist2(1)=curmoist(1)
       endif

       CALL RELHUMID
       if(RHLOCL.ge.99.)then
        RHLOCL=99.
       endif
       if(rhlocl.lt.0.D0)then
        rhlocl=0.D0
       endif
       oldmoist=curmoist
       oldcondep=condep
       wccfinal=0.D0
       curmoist=oldmoist
       condep=oldcondep
       timestep=3600/10
       step=int(3600/timestep)
       do 222 i=1,step
        call infil(rhlocl/100.,curmoist2,EP,soiltemp,dep(4:13),surflux
     &,wcc,curhumid2,curpot2,timestep,altt,
     &rww,pc,rl,sp,r1,im,maxcount,leafpot,curroot2,trans2)
        if(wcc.lt.0)then
         wcc=0
        endif
        if(surflux.lt.0)then
         surflux=0.D0
        endif
        j=1
        do 228 k=1,18
         if(mod(k,2).ne.0)then
          curmoist(j)=curmoist2(k)
          curhumid(j)=curhumid2(k)
          curpot(j)=curpot2(k)
          curroot(j)=curroot2(k)
          j=j+1
         endif
228     continue
        curmoist(10)=curmoist2(18)
        curhumid(10)=curhumid2(18)
        curpot(10)=curpot2(18)
        wccfinal=wccfinal+wcc
        condep=condep-WCC-surflux
        if(snowout.lt.1e-8)then
         ptwet=surflux/(ep*timestep)*100
        endif
        if(condep.lt.0.D0)then
         condep=0.D0
        endif
        if(condep.gt.0.D0)then
         curmoist(1)=1-BD(1)/DD(1)
         curmoist2(1)=curmoist(1)
        endif
222    continue
       if(condep.lt.0.D0)then
        condep=0.D0
       endif
       if(condep.gt.maxpool)then
        condep=maxpool
       endif
       if(snowout.lt.1e-8)then
        SABNEW = 1.0 - REFLS(DOY)
       endif
c      curmoist(1)=condep/((depp(2)*10)*(1-BD(1)/DD(1)))
       moists(1:10,DOY)=curmoist
       moist(1:10)=curmoist
       if(ptwet.lt.0.D0)then
        ptwet=0.D0
       endif
       if(ptwet.gt.100.)then
        ptwet=100.
       endif
       if(snowout.lt.1e-8)then
        if(condep.gt.0.D0)then
         if(condep.gt.10.)then ! use Fresnel's reflection law, p. 212 Gates 1980
          inrad = TAB('ZEN',TIME)*pi/180.
          refrad=asin(sin(inrad)/1.33)
          SABNEW = 1.-0.5*((sin(inrad-refrad)**2/(sin(inrad+refrad)**2))
     &    +(tan(inrad-refrad)**2/(tan(inrad+refrad)**2)))
         endif
        else
         SABNEW = 1.0 - REFLS(DOY)
        endif
       endif
      else
       moist(1:10)=moists(1:10,DOY)
       curmoist=moist
      endif
c*******************************************************************
c     end check for soil moisture model running

      if(gwsurf.lt.0)then
       gwsurf=0.D0
      endif
      if(int(SIOUT(1)).eq.1440)then
       methour=(int(1380/60)+1)+24*(DOY-1)
      else
       methour=(int(SIOUT(1)/60)+1)+24*(DOY-1)
      endif
      tide=0.D0
      tide=tides(methour,3)
c      get wave splash value and if greater than zero, override prev pctwet value
      if(tide.gt.0.D0)then
       PTWET=tides(methour,3)
      endif

      frosttest=(out(34)+out(4))/2.
      if((QEVAP.lt.-0.0000025).and.(out(34).lt.0))then
c    if((QEVAP.lt.-0.0000025).and.(out(4).lt.0))then
c    if((out(34).lt.0))then
c    if((QEVAP.lt.0).and.(frosttest.lt.0))then
       FROST = 1
      else
       FROST = 0
      endif

      if(QEVAP.lt.0)then
       DEW = 1
      else
       DEW = 0
      endif
C     END OF OUTPUT SETUP



c     SET UP LOCAL RELATIVE HUMIDITY
      CALL RELHUMID

      TKDAY2=(TKDAY/60.)*418.6
      SPDAY2=SPDAY*4185.
      DENDAY2=DENDAY*1.0E+3

      if(slipped.gt.0)then
       if(time.le.1380)then
        IF(SHAYD.LT.MAXSHD)THEN
         methour=int(temp(31)-1)
         metout(methour,1:18)=temp(1:18)
         metout(methour,19)=temp(83)
         soil(methour,1:12)=temp(19:30)
         if(runsnow.eq.1)then
          shdsnow(methour,1)=JULDAY(DOY)
          sunsnow(methour,2)=SIOUT(1)
          sunsnow(methour,3)=TEMP(62)
          sunsnow(methour,4)=TEMP(63)
          sunsnow(methour,5)=TEMP(64)
          sunsnow(methour,6)=TEMP(65)
          sunsnow(methour,7)=TEMP(66)
          sunsnow(methour,8)=TEMP(67)
          sunsnow(methour,9)=TEMP(68)
          sunsnow(methour,10)=TEMP(69)
          sunsnow(methour,11)=TEMP(70)
          soil(methour,3)=TEMP(70) ! make soil surface node equal to temp of bottom snow node
         endif
         if(runmoist.eq.1)then
          soilmoist(methour,1)=temp(19)
          soilmoist(methour,2)=temp(20)
          soilmoist(methour,3:12)=temp(32:41)
          humid(methour,1)=temp(19)
          humid(methour,2)=temp(20)
          humid(methour,3:12)=temp(42:51)
          soilpot(methour,1)=temp(19)
          soilpot(methour,2)=temp(20)
          soilpot(methour,3:12)=temp(52:61)
          plant(methour,1)=temp(19)
          plant(methour,2)=temp(20)
          plant(methour,3)=temp(71)
          plant(methour,4)=temp(72)
          plant(methour,5:14)=temp(73:82)
         endif
        ELSE
         methour=int(temp(31)-1)
         shadmet(methour,1:18)=temp(1:18)
         shadmet(methour,19)=temp(83)
         shadsoil(methour,1:12)=temp(19:30)
         if(runsnow.eq.1)then
          shdsnow(methour,1)=JULDAY(DOY)
          shdsnow(methour,2)=SIOUT(1)
          shdsnow(methour,3)=TEMP(62)
          shdsnow(methour,4)=TEMP(63)
          shdsnow(methour,5)=TEMP(64)
          shdsnow(methour,6)=TEMP(65)
          shdsnow(methour,7)=TEMP(66)
          shdsnow(methour,8)=TEMP(67)
          shdsnow(methour,9)=TEMP(68)
          shdsnow(methour,10)=TEMP(69)
          shdsnow(methour,11)=TEMP(70)
          shadsoil(methour,3)=TEMP(70) ! make soil surface node equal to temp of bottom snow node
         endif
         if(runmoist.eq.1)then
          shadmoist(methour,1)=temp(19)
          shadmoist(methour,2)=temp(20)
          shadmoist(methour,3:12)=temp(32:41)
          shadhumid(methour,1)=temp(19)
          shadhumid(methour,2)=temp(20)
          shadhumid(methour,3:12)=temp(42:51)
          shadpot(methour,1)=temp(19)
          shadpot(methour,2)=temp(20)
          shadpot(methour,3:12)=temp(52:61)
          shadplant(methour,1)=temp(19)
          shadplant(methour,2)=temp(20)
          shadplant(methour,3)=temp(71)
          shadplant(methour,4)=temp(72)
          shadplant(methour,5:14)=temp(73:82)
         endif
        ENDIF
       ENDIF
       slipped=0
c     end check for previous slippage
      endif

      IF (((DOY.eq.1).and.(time.lt.1e-8)).or.(int(TIME).NE.int(LASTIME))
     & .and.(int(max(maxval(abs(out(14:22))),abs(out(4)))).ne.81).or.
     & ((int(lastsurf).eq.81)).or.(int(abs(out(4))).eq.81).or.
     &((int(TIME).eq.int(LASTIME)).and.(int(max(maxval(abs(out(14:22))),
     &abs(out(4)))).ne.81))) THEN
       if(SIOUT(1).le.1380)then
        IF(SHAYD.LT.MAXSHD)THEN
         methour=0
         methour=(int(SIOUT(1)/60)+1)+24*(DOY-1)
         metout(methour,1)=JULDAY(DOY)
         metout(methour,2)=SIOUT(1)
         metout(methour,3)=SIOUT(2)
         metout(methour,4)=OUT(2)
         metout(methour,5)=RHLOCL
         metout(methour,6)=SIOUT(3)
         metout(methour,7)=SIOUT(4)
         metout(methour,8)=VEL2M
         metout(methour,9)=melted/1000.
         metout(methour,10)=condep
         metout(methour,11)=ptwet
         metout(methour,12)=SIOUT(8)
         metout(methour,13)=SIOUT(9)
         metout(methour,14)=SIOUT(10)
         metout(methour,15)=DEW
         metout(methour,16)=FROST
         metout(methour,17)=snowfall
         metout(methour,18)=cursnow
         metout(methour,19)=snowdens
         soil(methour,1)=JULDAY(DOY)
         soil(methour,2)=SIOUT(1)
         soil(methour,3)=OUT(4)
         soil(methour,4:12)=OUT(14:22)
         if(runsnow.eq.1)then
          sunsnow(methour,1)=JULDAY(DOY)
          sunsnow(methour,2)=SIOUT(1)
          sunsnow(methour,3)=TT(1)
          sunsnow(methour,4)=TT(2)
          sunsnow(methour,5)=TT(3)
          sunsnow(methour,6)=TT(4)
          sunsnow(methour,7)=TT(5)
          sunsnow(methour,8)=TT(6)
          sunsnow(methour,9)=TT(7)
          sunsnow(methour,10)=TT(8)
          sunsnow(methour,11)=TT(9)
          soil(methour,3)=TT(9)! make soil surface node equal to temp of bottom snow node
         endif
         if(runmoist.eq.1)then
          soilmoist(methour,1)=JULDAY(DOY)
          soilmoist(methour,2)=SIOUT(1)
          soilmoist(methour,3:12)=curmoist(1:10)
          humid(methour,1)=JULDAY(DOY)
          humid(methour,2)=SIOUT(1)
          humid(methour,3:12)=curhumid(1:10)
          soilpot(methour,1)=JULDAY(DOY)
          soilpot(methour,2)=SIOUT(1)
          soilpot(methour,3:12)=curpot(1:10)
          plant(methour,1)=JULDAY(DOY)
          plant(methour,2)=SIOUT(1)
          plant(methour,3)=trans2 * 1000 * 3600 !g/m2/h
          plant(methour,4)=leafpot
          plant(methour,5:14)=curroot(1:10)
         endif
         if(writecsv.eq.1)then
          WRITE(I3,154) metout(methour,1),",",metout(methour,2),",",
     &metout(methour,3),",",metout(methour,4),",",metout(methour,5),",",
     &metout(methour,6),",",metout(methour,7),",",metout(methour,8),",",
     &metout(methour,9),",",metout(methour,10),",",metout(methour,11),"
     &,",metout(methour,12),",",metout(methour,13),",",metout(methour,14
     &),",",metout(methour,15),",",metout(methour,16),",",
     &metout(methour,17),",",metout(methour,18),",",metout(methour,19)
          WRITE(I10,160) soil(methour,1),",",soil(methour,2),",",
     &soil(methour,3),",",soil(methour,4),",",soil(methour,5),",",
     &soil(methour,6),",",soil(methour,7),",",soil(methour,8),",",
     &soil(methour,9),",",soil(methour,10),",",soil(methour,11),"
     &,",soil(methour,12)
          if(runmoist.eq.1)then
      WRITE(I91,163) soilmoist(methour,1),",",soilmoist(methour,2),","
     &,soilmoist(methour,3),",",soilmoist(methour,4),",",soilmoist(metho
     &ur,5),",",soilmoist(methour,6),",",soilmoist(methour,7),",",soilmo
     &ist(methour,8),",",soilmoist(methour,9),",",soilmoist(methour,10)
     &,",",soilmoist(methour,11),",",soilmoist(methour,12)
      WRITE(I92,160) humid(methour,1),",",humid(methour,2),","
     &,humid(methour,3),",",humid(methour,4),",",humid(metho
     &ur,5),",",humid(methour,6),",",humid(methour,7),",",humid
     &(methour,8),",",humid(methour,9),",",humid(methour,10)
     &,",",humid(methour,11),",",humid(methour,12)
      WRITE(I95,161) soilpot(methour,1),",",soilpot(methour,2),","
     &,soilpot(methour,3),",",soilpot(methour,4),",",soilpot(metho
     &ur,5),",",soilpot(methour,6),",",soilpot(methour,7),",",soilpot
     &(methour,8),",",soilpot(methour,9),",",soilpot(methour,10)
     &,",",soilpot(methour,11),",",soilpot(methour,12)
      WRITE(I100,161) plant(methour,1),",",plant(methour,2),","
     &,plant(methour,3),",",plant(methour,4),",",plant(metho
     &ur,5),",",plant(methour,6),",",plant(methour,7),",",plant
     &(methour,8),",",plant(methour,9),",",plant(methour,10)
     &,",",plant(methour,11),",",plant(methour,12),",",plant(methour,13)
     &,",",plant(methour,14)
        endif
      if(runsnow.eq.1)then
      WRITE(I7,159) JULDAY(DOY),",",SIOUT(1),",",
     &TT(1),",",TT(2),",",TT(3),",",
     &TT(4),",",TT(5),",",TT(6),",",
     &TT(7),",",TT(8),",",TT(9)
      endif
        if(lamb.eq.1)then
           WRITE(I97,164) DRLAMBDA(methour
     &,1),",",DRLAMBDA(methour,2),",",DRLAMBDA(methour,3),",",DRLAMBDA
     &(methour,4),",",DRLAMBDA(methour,5),",",DRLAMBDA(methour,6),","
     &,DRLAMBDA(methour,7),",",DRLAMBDA(methour,8),",",DRLAMBDA(methou
     &r,9),",",DRLAMBDA(methour,10),",",DRLAMBDA(methour,11),",",DRLA
     &MBDA(methour,12),",",DRLAMBDA(methour,13),",",DRLAMBDA(methour,1
     &4),",",DRLAMBDA(methour,15),",",DRLAMBDA(methour,16),",",DRLAMB
     &DA(methour,17),",",DRLAMBDA(methour,18),",",DRLAMBDA(methour,19
     &),",",DRLAMBDA(methour,20),",",DRLAMBDA(methour,21),",",DRLAMBDA
     &(methour,22),",",DRLAMBDA(methour,23),",",DRLAMBDA(methour,24)
     &,",",DRLAMBDA(methour,25),",",DRLAMBDA(methour,26),",",DRLAMBDA
     &(methour,27),",",DRLAMBDA(methour,28),",",DRLAMBDA(methour,29)
     &,",",DRLAMBDA(methour,30),",",DRLAMBDA(methour,31),",",DRLAMBDA
     &(methour,32),",",DRLAMBDA(methour,33),",",DRLAMBDA(methour,34)
     &,",",DRLAMBDA(methour,35),",",DRLAMBDA(methour,36),",",DRLAMBDA
     &(methour,37),",",DRLAMBDA(methour,38),",",DRLAMBDA(methour,39)
     &,",",DRLAMBDA(methour,40),",",DRLAMBDA(methour,41),",",DRLAMBDA
     &(methour,42),",",DRLAMBDA(methour,43),",",DRLAMBDA(methour,44)
     &,",",DRLAMBDA(methour,45),",",DRLAMBDA(methour,46),",",DRLAMBDA
     &(methour,47),",",DRLAMBDA(methour,48),",",DRLAMBDA(methour,49)
     &,",",DRLAMBDA(methour,50),",",DRLAMBDA(methour,51),",",DRLAMBDA
     &(methour,52),",",DRLAMBDA(methour,53),",",DRLAMBDA(methour,54)
     &,",",DRLAMBDA(methour,55),",",DRLAMBDA(methour,56),",",DRLAMBDA
     &(methour,57),",",DRLAMBDA(methour,58),",",DRLAMBDA(methour,59)
     &,",",DRLAMBDA(methour,60),",",DRLAMBDA(methour,61),",",DRLAMBDA
     &(methour,62),",",DRLAMBDA(methour,63),",",DRLAMBDA(methour,64)
     &,",",DRLAMBDA(methour,65),",",DRLAMBDA(methour,66),",",DRLAMBDA
     &(methour,67),",",DRLAMBDA(methour,68),",",DRLAMBDA(methour,69)
     &,",",DRLAMBDA(methour,70),",",DRLAMBDA(methour,71),",",DRLAMBDA
     &(methour,72),",",DRLAMBDA(methour,73),",",DRLAMBDA(methour,74)
     &,",",DRLAMBDA(methour,75),",",DRLAMBDA(methour,76),",",DRLAMBDA
     &(methour,77),",",DRLAMBDA(methour,78),",",DRLAMBDA(methour,79)
     &,",",DRLAMBDA(methour,80),",",DRLAMBDA(methour,81),",",DRLAMBDA
     &(methour,82),",",DRLAMBDA(methour,83),",",DRLAMBDA(methour,84)
     &,",",DRLAMBDA(methour,85),",",DRLAMBDA(methour,86),",",DRLAMBDA
     &(methour,87),",",DRLAMBDA(methour,88),",",DRLAMBDA(methour,89)
     &,",",DRLAMBDA(methour,90),",",DRLAMBDA(methour,91),",",DRLAMBDA
     &(methour,92),",",DRLAMBDA(methour,93),",",DRLAMBDA(methour,94)
     &,",",DRLAMBDA(methour,95),",",DRLAMBDA(methour,96),",",DRLAMBDA
     &(methour,97),",",DRLAMBDA(methour,98),",",DRLAMBDA(methour,99)
     &,",",DRLAMBDA(methour,100),",",DRLAMBDA(methour,101),",",
     &DRLAMBDA(methour,102),",",DRLAMBDA(methour,103),",",DRLAMBDA
     &(methour,104),",",DRLAMBDA(methour,105),",",DRLAMBDA(methour,
     &106),",",DRLAMBDA(methour,107),",",DRLAMBDA(methour,108),",",
     &DRLAMBDA(methour,109),",",DRLAMBDA(methour,110),",",
     &DRLAMBDA(methour,111),",",DRLAMBDA(methour,112),",",
     &DRLAMBDA(methour,113)

      WRITE(I98,164) DRRLAMBDA(methour,
     &1),",",DRRLAMBDA(methour,2),",",DRRLAMBDA(methour,3),",",DRRLAM
     &BDA(methour,4),",",DRRLAMBDA(methour,5),",",DRRLAMBDA(methour,6
     &),",",DRRLAMBDA(methour,7),",",DRRLAMBDA(methour,8),",",DRRLAMBD
     &A(methour,9),",",DRRLAMBDA(methour,10),",",DRRLAMBDA(methour,11
     &),",",DRLAMBDA(methour,12),",",DRRLAMBDA(methour,13),",",DRRLAMB
     &DA(methour,14),",",DRRLAMBDA(methour,15),",",DRRLAMBDA(methour,1
     &6),",",DRLAMBDA(methour,17),",",DRRLAMBDA(methour,18),",",DRRLA
     &MBDA(methour,19),",",DRRLAMBDA(methour,20),",",DRRLAMBDA(methour
     &,21),",",DRRLAMBDA(methour,22),",",DRRLAMBDA(methour,23),",",DR
     &RLAMBDA(methour,24),",",DRRLAMBDA(methour,25),",",DRRLAMBDA(meth
     &our,26),",",DRRLAMBDA(methour,27),",",DRRLAMBDA(methour,28),","
     &,DRRLAMBDA(methour,29),",",DRRLAMBDA(methour,30),",",DRRLAMBDA(m
     &ethour,31),",",DRRLAMBDA(methour,32),",",DRRLAMBDA(methour,33)
     &,",",DRRLAMBDA(methour,34),",",DRRLAMBDA(methour,35),",",DRRLAMB
     &DA(methour,36),",",DRRLAMBDA(methour,37),",",DRRLAMBDA(methour,3
     &8),",",DRRLAMBDA(methour,39),",",DRRLAMBDA(methour,40),",",DRRL
     &AMBDA(methour,41),",",DRRLAMBDA(methour,42),",",DRRLAMBDA(methou
     &r,43),",",DRRLAMBDA(methour,44),",",DRRLAMBDA(methour,45),",",
     &DRRLAMBDA(methour,46),",",DRRLAMBDA(methour,47),",",DRRLAMBDA(me
     &thour,48),",",DRRLAMBDA(methour,49),",",DRRLAMBDA(methour,50)
     &,",",DRRLAMBDA(methour,51),",",DRRLAMBDA(methour,52),",",DRRLAMB
     &DA(methour,53),",",DRRLAMBDA(methour,54),",",DRRLAMBDA(methour
     &,55),",",DRRLAMBDA(methour,56),",",DRRLAMBDA(methour,57),",",DR
     &RLAMBDA(methour,58),",",DRRLAMBDA(methour,59),",",DRRLAMBDA(meth
     &our,60),",",DRRLAMBDA(methour,61),",",DRRLAMBDA(methour,62),","
     &,DRRLAMBDA(methour,63),",",DRRLAMBDA(methour,64),",",DRRLAMBDA(m
     &ethour,65),",",DRRLAMBDA(methour,66),",",DRRLAMBDA(methour,67)
     &,",",DRRLAMBDA(methour,68),",",DRRLAMBDA(methour,69),",",DRRLAMB
     &DA(methour,70),",",DRRLAMBDA(methour,71),",",DRRLAMBDA(methour
     &,72),",",DRRLAMBDA(methour,73),",",DRRLAMBDA(methour,74),",",DR
     &RLAMBDA(methour,75),",",DRRLAMBDA(methour,76),",",DRRLAMBDA(meth
     &our,77),",",DRRLAMBDA(methour,78),",",DRRLAMBDA(methour,79),","
     &,DRRLAMBDA(methour,80),",",DRRLAMBDA(methour,81),",",DRRLAMBDA
     &(methour,82),",",DRRLAMBDA(methour,83),",",DRRLAMBDA(methour,84
     &),",",DRRLAMBDA(methour,85),",",DRRLAMBDA(methour,86),",",DRRLAM
     &BDA(methour,87),",",DRRLAMBDA(methour,88),",",DRRLAMBDA(methour
     &,89),",",DRRLAMBDA(methour,90),",",DRRLAMBDA(methour,91),",",DR
     &RLAMBDA(methour,92),",",DRRLAMBDA(methour,93),",",DRRLAMBDA(meth
     &our,94),",",DRRLAMBDA(methour,95),",",DRRLAMBDA(methour,96),","
     &,DRRLAMBDA(methour,97),",",DRRLAMBDA(methour,98),",",DRRLAMBDA(m
     &ethour,99),",",DRRLAMBDA(methour,100),",",DRRLAMBDA(methour,101
     &),",",DRRLAMBDA(methour,102),",",DRRLAMBDA(methour,103),",",DRRL
     &AMBDA(methour,104),",",DRRLAMBDA(methour,105),",",DRRLAMBDA(meth
     &our,106),",",DRRLAMBDA(methour,107),",",DRRLAMBDA(methour,108)
     &,",",DRRLAMBDA(methour,109),",",DRRLAMBDA(methour,110),",",
     &DRRLAMBDA(methour,111),",",DRRLAMBDA(methour,112),",",
     &DRRLAMBDA(methour,113)

        WRITE(I99,164) SRLAMBDA(methour,1
     &),",",SRLAMBDA(methour,2),",",SRLAMBDA(methour,3),",",SRLAMBDA
     &(methour,4),",",SRLAMBDA(methour,5),",",SRLAMBDA(methour,6),","
     &,SRLAMBDA(methour,7),",",SRLAMBDA(methour,8),",",SRLAMBDA(methou
     &r,9),",",SRLAMBDA(methour,10),",",SRLAMBDA(methour,11),",",SRLA
     &MBDA(methour,12),",",SRLAMBDA(methour,13),",",SRLAMBDA(methour,1
     &4),",",SRLAMBDA(methour,15),",",SRLAMBDA(methour,16),",",SRLAMB
     &DA(methour,17),",",SRLAMBDA(methour,18),",",SRLAMBDA(methour,19
     &),",",SRLAMBDA(methour,20),",",SRLAMBDA(methour,21),",",SRLAMBDA
     &(methour,22),",",SRLAMBDA(methour,23),",",SRLAMBDA(methour,24)
     &,",",SRLAMBDA(methour,25),",",SRLAMBDA(methour,26),",",SRLAMBDA
     &(methour,27),",",SRLAMBDA(methour,28),",",SRLAMBDA(methour,29)
     &,",",SRLAMBDA(methour,30),",",SRLAMBDA(methour,31),",",SRLAMBDA
     &(methour,32),",",SRLAMBDA(methour,33),",",SRLAMBDA(methour,34)
     &,",",SRLAMBDA(methour,35),",",SRLAMBDA(methour,36),",",SRLAMBDA
     &(methour,37),",",SRLAMBDA(methour,38),",",SRLAMBDA(methour,39)
     &,",",SRLAMBDA(methour,40),",",SRLAMBDA(methour,41),",",SRLAMBDA
     &(methour,42),",",SRLAMBDA(methour,43),",",SRLAMBDA(methour,44)
     &,",",SRLAMBDA(methour,45),",",SRLAMBDA(methour,46),",",SRLAMBDA
     &(methour,47),",",SRLAMBDA(methour,48),",",SRLAMBDA(methour,49)
     &,",",SRLAMBDA(methour,50),",",SRLAMBDA(methour,51),",",SRLAMBDA
     &(methour,52),",",SRLAMBDA(methour,53),",",SRLAMBDA(methour,54)
     &,",",SRLAMBDA(methour,55),",",SRLAMBDA(methour,56),",",SRLAMBDA
     &(methour,57),",",SRLAMBDA(methour,58),",",SRLAMBDA(methour,59)
     &,",",SRLAMBDA(methour,60),",",SRLAMBDA(methour,61),",",SRLAMBDA
     &(methour,62),",",SRLAMBDA(methour,63),",",SRLAMBDA(methour,64)
     &,",",SRLAMBDA(methour,65),",",SRLAMBDA(methour,66),",",SRLAMBDA
     &(methour,67),",",SRLAMBDA(methour,68),",",SRLAMBDA(methour,69)
     &,",",SRLAMBDA(methour,70),",",SRLAMBDA(methour,71),",",SRLAMBDA
     &(methour,72),",",SRLAMBDA(methour,73),",",SRLAMBDA(methour,74)
     &,",",SRLAMBDA(methour,75),",",SRLAMBDA(methour,76),",",SRLAMBDA
     &(methour,77),",",SRLAMBDA(methour,78),",",SRLAMBDA(methour,79)
     &,",",SRLAMBDA(methour,80),",",SRLAMBDA(methour,81),",",SRLAMBDA
     &(methour,82),",",SRLAMBDA(methour,83),",",SRLAMBDA(methour,84)
     &,",",SRLAMBDA(methour,85),",",SRLAMBDA(methour,86),",",SRLAMBDA
     &(methour,87),",",SRLAMBDA(methour,88),",",SRLAMBDA(methour,89)
     &,",",SRLAMBDA(methour,90),",",SRLAMBDA(methour,91),",",SRLAMBDA
     &(methour,92),",",SRLAMBDA(methour,93),",",SRLAMBDA(methour,94)
     &,",",SRLAMBDA(methour,95),",",SRLAMBDA(methour,96),",",SRLAMBDA
     &(methour,97),",",SRLAMBDA(methour,98),",",SRLAMBDA(methour,99)
     &,",",SRLAMBDA(methour,100),",",SRLAMBDA(methour,101),",",
     &SRLAMBDA(methour,102),",",SRLAMBDA(methour,103),",",SRLAMBDA
     &(methour,104),",",SRLAMBDA(methour,105),",",SRLAMBDA(methour,
     &106),",",SRLAMBDA(methour,107),",",SRLAMBDA(methour,108),",",
     &SRLAMBDA(methour,109),",",SRLAMBDA(methour,110),",",
     &SRLAMBDA(methour,111),",",SRLAMBDA(methour,112),",",
     &SRLAMBDA(methour,113)
          endif
         endif
        ELSE
         methour=0
         methour=int(SIOUT(1)/60)+1+24*(DOY-1)
         shadmet(methour,1)=JULDAY(DOY)
         shadmet(methour,2)=SIOUT(1)
         shadmet(methour,3)=SIOUT(2)
         shadmet(methour,4)=OUT(2)
         shadmet(methour,5)=RHLOCL
         shadmet(methour,6)=SIOUT(3)
         shadmet(methour,7)=SIOUT(4)
         shadmet(methour,8)=VEL2M
         shadmet(methour,9)=melt/1000.
         shadmet(methour,10)=condep
         shadmet(methour,11)=ptwet
         shadmet(methour,12)=SIOUT(8)
         shadmet(methour,13)=SIOUT(9)
         shadmet(methour,14)=SIOUT(10)
         shadmet(methour,15)=DEW
         shadmet(methour,16)=FROST
         shadmet(methour,17)=snowfall
         shadmet(methour,18)=cursnow
         shadmet(methour,19)=snowdens
         shadsoil(methour,1)=JULDAY(DOY)
         shadsoil(methour,2)=SIOUT(1)
         shadsoil(methour,3)=OUT(4)
         shadsoil(methour,4:12)=OUT(14:22)
         if(runsnow.eq.1)then
          shdsnow(methour,1)=JULDAY(DOY)
          shdsnow(methour,2)=SIOUT(1)
          shdsnow(methour,3)=TT(1)
          shdsnow(methour,4)=TT(2)
          shdsnow(methour,5)=TT(3)
          shdsnow(methour,6)=TT(4)
          shdsnow(methour,7)=TT(5)
          shdsnow(methour,8)=TT(6)
          shdsnow(methour,9)=TT(7)
          shdsnow(methour,10)=TT(8)
          shdsnow(methour,11)=TT(9)
          shadsoil(methour,3)=TT(9) ! make soil surface node equal to temp of bottom snow node
         endif
         if(runmoist.eq.1)then
          shadmoist(methour,1)=JULDAY(DOY)
          shadmoist(methour,2)=SIOUT(1)
          shadmoist(methour,3:12)=curmoist(1:10)
          shadhumid(methour,1)=JULDAY(DOY)
          shadhumid(methour,2)=SIOUT(1)
          shadhumid(methour,3:12)=curhumid(1:10)
          shadpot(methour,1)=JULDAY(DOY)
          shadpot(methour,2)=SIOUT(1)
          shadpot(methour,3:12)=curpot(1:10)
          shadplant(methour,1)=JULDAY(DOY)
          shadplant(methour,2)=SIOUT(1)
          shadplant(methour,3)=trans2 * 1000 * 3600 !g/m2/h
          shadplant(methour,4)=leafpot
          shadplant(methour,5:14)=curroot(1:10)
         endif
         if((writecsv.eq.1).and.(runshade.eq.1))then
          WRITE(I12,154) shadmet(methour,1),",",shadmet(methour,2),",",
     &shadmet(methour,3),",",shadmet(methour,4),",",shadmet(methour,5)
     &,",",shadmet(methour,6),",",shadmet(methour,7),",",shadmet(methour
     &,8),",",shadmet(methour,9),",",shadmet(methour,10),",",shadmet(met
     &hour,11),",",shadmet(methour,12),",",shadmet(methour,13),",",shadm
     &et(methour,14),",",shadmet(methour,15),",",shadmet(methour,16),","
     &,shadmet(methour,17),",",shadmet(methour,18),",",shadmet(methour,1
     &9)
         WRITE(I11,157) shadsoil(methour,1),",",shadsoil(methour,2),",",
     &shadsoil(methour,3),",",shadsoil(methour,4),",",
     &shadsoil(methour,5),",",shadsoil(methour,6),",",
     &shadsoil(methour,7),",",shadsoil(methour,8),",",
     &shadsoil(methour,9),",",shadsoil(methour,10),",",
     &shadsoil(methour,11),",",shadsoil(methour,12)
      if(runsnow.eq.1)then
      WRITE(I8,159) JULDAY(DOY),",",SIOUT(1),",",
     &TT(1),",",TT(2),",",TT(3),",",
     &TT(4),",",TT(5),",",TT(6),",",
     &TT(7),",",TT(8),",",TT(9)
      endif
          if(runmoist.eq.1)then
      WRITE(I93,163) shadmoist(methour,1),",",shadmoist(methour,2),","
     &,shadmoist(methour,3),",",shadmoist(methour,4),",",shadmoist(metho
     &ur,5),",",shadmoist(methour,6),",",shadmoist(methour,7),",",shadmo
     &ist(methour,8),",",shadmoist(methour,9),",",shadmoist(methour,10)
     &,",",shadmoist(methour,11),",",shadmoist(methour,12)
      WRITE(I94,160) shadhumid(methour,1),",",shadhumid(methour,2),","
     &,shadhumid(methour,3),",",shadhumid(methour,4),",",shadhumid(metho
     &ur,5),",",shadhumid(methour,6),",",shadhumid(methour,7),",",shadhu
     &mid(methour,8),",",shadhumid(methour,9),",",shadhumid(methour,10)
     &,",",shadhumid(methour,11),",",shadhumid(methour,12)
      WRITE(I96,161) shadpot(methour,1),",",shadpot(methour,2),","
     &,shadpot(methour,3),",",shadpot(methour,4),",",shadpot(metho
     &ur,5),",",shadpot(methour,6),",",shadpot(methour,7),",",shadpot
     &(methour,8),",",shadpot(methour,9),",",shadpot(methour,10)
     &,",",shadpot(methour,11),",",shadpot(methour,12)
      WRITE(I101,161) shadplant(methour,1),",",shadplant(methour,2),","
     &,shadplant(methour,3),",",shadplant(methour,4),",",shadplant(metho
     &ur,5),",",shadplant(methour,6),",",shadplant(methour,7),","
     &,shadplant(methour,8),",",shadplant(methour,9),","
     &,shadplant(methour,10),",",shadplant(methour,11),","
     &,shadplant(methour,12),",",shadplant(methour,13),","
     &,shadplant(methour,14)
          endif
         endif
        ENDIF
       endif
c     check if duplicate time due to integrator slipping
      else
       if(int(max(maxval(abs(out(14:22))),abs(out(4)))).ne.81)then
        slipped=slipped+1
        errcount=errcount+1
        temp(31)=(int(SIOUT(1)/60)+1)+24*(DOY-1)
        temp(1)=JULDAY(DOY)
        temp(2)=SIOUT(1)
        temp(3)=SIOUT(2)
        temp(4)=OUT(2)
        temp(5)=RHLOCL
        temp(6)=SIOUT(3)
        temp(7)=SIOUT(4)
        temp(8)=VEL2M
        temp(9)=melt/1000.
        temp(10)=condep
        temp(11)=ptwet
        temp(12:14)=SIOUT(8:10)
        temp(15)=DEW
        temp(16)=FROST
        temp(17)=snowfall
        temp(18)=snowout
        temp(19)=JULDAY(DOY)
        temp(20)=SIOUT(1)
        temp(21)=OUT(4)
        temp(22:30)=OUT(14:22)
        if(runsnow.eq.1)then
         temp(62:70)=TT(1:9)
         temp(83)=snowdens
        endif
        if(runmoist.eq.1)then
         temp(32:41)=curmoist(1:10)
         temp(42:51)=curhumid(1:10)
         temp(52:61)=curpot(1:10)
         temp(71)=trans2 * 1000 * 3600 !g/m2/h
         temp(72)=leafpot
         temp(73:82)=curroot(1:10)
        endif
       endif
c     end check for current slippage
      endif

      IF (int(SIOUT(1)) .EQ. 1440) THEN
C      INCREMENT MONTH/TIME INTERVAL COUNTER
       IF (((DOY.eq.1).and.(time.lt.1e-8)).or.(TIME .NE. LASTIME)) THEN
        DOY = DOY + 1
       ENDIF
      ENDIF

      lastime=time
      lastsurf=max(maxval(abs(out(14:22))),abs(out(4)))
  154 FORMAT(1F4.0,A,1F7.2,A,F6.2,A,F6.2,A,F6.2,A,F6.2,A,F7.3,A,F7.3,A,
     &F6.2,A,F6.2,A,F6.2,A,1F6.2,A,1F7.2,A,1F7.2,A,1F2.0,A,1F2.0,A,F7.2,
     &A,F7.2,A,F7.3)
  157 FORMAT(1F4.0,A,1F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A
     & ,F7.2,A,F7.2,A,F7.2,A,F7.2)
  160 FORMAT(1F4.0,A,1F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A
     & ,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2)
  163 FORMAT(1F4.0,A,1F7.2,A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,F7.4,A
     & ,F7.4,A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,F7.4)
  161 FORMAT(1F4.0,A,1F10.2,A,F10.2,A,F10.2,A,F10.2,A,F10.2,A,F10.2
     & ,A,F10.2,A,F10.2,A,F10.2,A,F10.2,A,F10.2,A,F10.2,A,F10.2)
  164 FORMAT(1F4.0,A,1F10.2,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,F10.8,A,
     &F10.8,A,F10.8)
  159 FORMAT(1F4.0,A,1F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A,F7.2,A
     & ,F7.2,A,F7.2,A,F7.2)
  200 RETURN
      END

