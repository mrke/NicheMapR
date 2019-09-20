      SUBROUTINE MICRO(Z,Z0,ZH,D0,T1,T3,V,QC,AMOL,NAIR,ZZ,VV,T,ZEN)
      IMPLICIT NONE

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

C    This subroutine computes a single unsegmented velocity and temperature profile

      double precision A,ADUM,AMOL,AMOLN,D0,DEL,DIFFT,DUM,GAM,PHI,PSI1
      double precision PSI2,QC,RCP,RCPTKG,RHOCP,STB,STO,STS,T,T1,T3,TAVE
     & ,T0,TZO,USTAR
      double precision V,VEL,VV,X,X1,Y,Y1,YY,YY2,Z,Z0,Z01
      double precision Z02,ZEN,ZH,ZH1,ZH2,ZRATIO,ZZ

      INTEGER I,ITER,NAIR,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92
     & ,I93,I94,I95,I96,I97,I98,I99,I100,I101
      DIMENSION ZZ(10),VV(10),T(30)
      COMMON/DMYCRO/Z01,Z02,ZH1,ZH2
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101

C
C**** 1 SEGMENT VELOCITY PROFILE - W. PORTER
C**** VELOCITY PROFILE - Businger, J. A., Wyngaard, J. C., Izumi, Y., & Bradley, E. F. (1971). Flux-Profile Relationships in the Atmospheric Surface Layer. Journal of the Atmospheric Sciences, 28(2), 181–189. doi:10.1175/1520-0469(1971)028<0181:FPRITA>2.0.CO;2
C**** SUBLAYER MODEL - Garratt, J. R., & Hicks, B. B. (1973). Momentum, heat and water vapour transfer to and from natural and artificial surfaces. Quarterly Journal of the Royal Meteorological Society, 99(422), 680–687. doi:10.1002/qj.49709942209
C     Z=REFERENCE HEIGHT
C     Z0=ROUGHNESS HEIGHT
C     T1=TEMPERATURE AT REFERENCE HEIGHT
C     T3=CURRENT ESTIMATE OF GROUND SURFACE TEMP.
C     V=VELOCITY AT REF. HEIGHT
C     QC=COMPUTED (HERE) CONVECTIVE HEAT TRANSFER AT THE SURFACE
C     AMOL=MONIN-OBUKHOV LENGTH
C     NAIR=NO. OF HEIGHTS FOR AIR TEMP'S AND VELOCITIES
C     ZZ=ARRAY OF HEIGHT VALUES
C     VV=ARRAY OF COMPUTED (HERE) VELOCITIES FOR EACH HEIGHT
C     T=ARRAY OF AIR TEMP'S COMPUTED HERE ASSUMING A LOG PROFILE
C
C     THIS SUBROUTINE IS MODIFIED (FEB. 1979)FOR SHEAR OCCURRING ABOVE
C     THE SURFACE DUE TO VEGETATION SPACED OVER THE SURFACE.
C     TEMP. PROFILE REMAINS LOGARITHMIC. VEL. PROFILE LOGARITHMIC IN SEGMENTS
C     SEGMENTS OF VEL. PROFILE= 200-100 CM, 100-30 CM, 30-0 CM.
C     TREF=200 CM, VREF=30 CM FOR SANTA FE, GALAPAGOS
C     Z0 IS PLOTTED FROM 30 CM VEL'S DOWN
C
C     ****WHEN STARTING AT MIDNIGHT ON THE VERY FIRST
C     ITERATION BE SURE  BE SURE  BE SURE
C     INITIAL TSURF GUESS IS LESS THAN TREF
C     SO MICRO WILL GO TO LOWER HALF
C
C     DEFINING Z0'S FROM THE TOP DOWN            Steve's
C     GALAPAGOS  TEXAS, WASHINGTON    NEVADA     Carlsbad NM
C     Z01=16.8    11.16               3.67       8.353
C     Z02= 6.42   10.57               3.29       3.015
C     Z0 = LOWEST (REAL) ROUGHNESS HEIGHT
C     Z0 =       0.021                0.90       0.268
C
C     DEFINING HEIGHTS WHERE Z0 CHANGES FROM THE TOP
C     Z = 'FREE STREAM' REFERENCE HEIGHT = 200 CM originally, now set by user
C     ZH1=100.   84                   80           50
C     ZH2=30.    13.                  60.          25

      RHOCP(TAVE) = 0.08472/TAVE ! note this is a function, internally defined
      PHI(Z)=(1.-GAM*Z/AMOL)**.25
      PSI1(X)=2.*dLOG((1.+X)/2.)+dLOG((1.+X*X)/2.)-2.*ATAN(X)+3.14159/2
      PSI2(X)=2.*dLOG((1.+X*X)/2.)
      GAM=16.
      RCPTKG=6.003E-8 !RHO*CP*T/(K*G) = 6.003E-8 IN CAL-MIN-CM-C UNITS

C     COMPUTING VEL. PROFILE PARAMETERS FROM 200 CM REFERENCE VELOCITY
      ZRATIO = Z/Z0 + 1 ! ratio of reference to roughness height
      DUM=dLOG(ZRATIO)
      VEL = V ! wind speed at reference height
      USTAR = 0.4*V/DUM ! friction velocity
      DIFFT=T1-T3 ! temp at reference height minus ground temp
      TAVE=(T3+T1+546.)/2. ! ave temp in Kelvin
      RCP=RHOCP(TAVE)
      AMOL=-30.0 ! initial Monin-Obukhov length
      ITER=0 ! initialise counter

C	  Paul edit 9/12/19: adding alternative Campbell and Norman 1998 vertical air temperature profile calculation option
      IF(ZH.GT.0)GO TO 1500


C     CHECK FOR FREE CONVECTION (LAPSE) CONDITIONS
      IF(T1.GE.T3)GO TO 1000
      IF(T3.LE.81.)GO TO 1000
      IF(ZEN .GE. 81.)GO TO 1000

C     NEGLECTING FREE CONV. CORRECTION (4%)FOR SEGMENTED PROFILES.

C     ITERATING TO FIND THE MONIN-OBUKHOV LENGTH (AMOL)
   1  X=PHI(Z)
      Y=PSI1(X)
      YY=PSI2(X)
      USTAR=0.4*V/(dLOG(Z/Z0)-Y)
      STS=.62/(Z0*USTAR/12.)**.45
C     BULK STANTON NO.
      STB=(.64/DUM)*(1.-.1*Z/AMOL)
      STO=STB/(1.+STB/STS)
C
      QC=RCP*DIFFT*USTAR*STO
C
      AMOLN=RCPTKG*USTAR**3/QC
      DEL=ABS((AMOLN-AMOL)/AMOL)
      IF (DEL .LT. 1.0E-02) THEN
       GO TO 2
      ENDIF
      AMOL=AMOLN
      ITER=ITER+1
      IF(ITER .GT.30) GO TO 2000
      GO TO 1
C     END OF ITERATION LOOP TO FIND MONIN-OBUKHOV LENGTH

    2 CONTINUE
      IF(NAIR.LE.0) RETURN
      DO 3 I=1,NAIR
       X1=PHI(ZZ(I))
       Y1=PSI1(X1)
       YY2=PSI2(X1)
C      FILL OUT VELOCITY AND TEMP. PROFILES
       ADUM=ZZ(I)/Z0-Y1
       VV(I)=2.5*USTAR*dLOG(ADUM)
       TZO=(T1*STB+T3*STS)/(STB+STS)
       T(I+20)=TZO+(T1-TZO)*dLOG(ZZ(I)/Z0-YY2)/dLOG(Z/Z0-YY)
    3 CONTINUE
      RETURN

C     CALC'S BELOW WHEN NO FREE CONV. ENHANCEMENT OF VEL,TEMP PROFILES
 1000 CONTINUE
      STS=.62/(Z0*USTAR/12.)**.45 !SUBLAYER STANTON NO.
      STB=.64/DUM ! BULK STANTON NO.

      QC=RCP*DIFFT*USTAR*STB/(1+STB/STS) ! convective heat transfer at the surface

      IF(NAIR.LE.0) RETURN
      DO 4 I=1,NAIR
C      FILL OUT VEL. AND TEMP. PROFILES
       VV(I)=2.5*USTAR*dLOG(ZZ(I)/Z0+1)
C      COMPUTING FICTITIOUS TEMP. AT TOP OF SUBLAYER
       TZO=(T1*STB+T3*STS)/(STB+STS)
       T(I+20)=TZO+(T1-TZO)*dLOG(ZZ(I)/Z0+1)/DUM
    4 CONTINUE
      RETURN
    
1500  CONTINUE
      STS=.62/(Z0*USTAR/12.)**.45 !SUBLAYER STANTON NO.
      STB=.64/DUM ! BULK STANTON NO.

      QC=RCP*DIFFT*USTAR*STB/(1+STB/STS) ! convective heat transfer at the surface
C	  Use vertical temperature profile from Campbell and Norman 1998
      IF(NAIR.LE.0) RETURN
      DO 5 I=1,NAIR
C      FILL OUT VEL. AND TEMP. PROFILES
       IF((T1.GE.T3).or.(T3.LE.81.).or.(ZEN .GE. 81.))THEN
        VV(I)=2.5*USTAR*dLOG(ZZ(I)/Z0+1)
       ELSE
        X1=PHI(ZZ(I))
        Y1=PSI1(X1)
        ADUM=ZZ(I)/Z0-Y1
        VV(I)=2.5*USTAR*dLOG(ADUM)
       ENDIF
       A=(T1-T3)/(1-dLOG((Z-D0)/ZH))
       T0=T1+A*dLOG((Z-D0)/ZH)
       T(I+20)=T0-A*dLOG((ZZ(I)-D0)/ZH)
    5 CONTINUE
      RETURN

 2000 continue
      RETURN
      END