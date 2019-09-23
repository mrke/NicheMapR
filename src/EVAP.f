      SUBROUTINE EVAP(TSURF,TAIR,RELHUM,HD,QEVAP,SAT) 

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

C     THIS SUBROUTINE COMPUTES SURFACE WATER LOSS.  
c     warning - this is all in SI units!

      IMPLICIT NONE  
      double precision AIRVD,ALTT,BP,CP,DB,DENAIR,DP,E,EFFSUR,ESAT
      double precision HD,HTOVPR,MAXSHD,PSTD,PTWET,GWSURF
      double precision QEVAP,RH,RELHUM,RW,SABNEW,SHAYD,SURFVD
      double precision TAIR,TSURF,TVIR 
      double precision TVINC,VD,WATER,WB,WTRPOT,rainfall
      integer sat

      COMMON/GROUND/SHAYD,ALTT,MAXSHD,SABNEW,PTWET,rainfall
C     CALCULATING SURFACE SATURATION VAPOR DENSITY
      RH = 100.
C     CHECK FOR TOO LOW A SURFACE TEMPERATURE
      IF (TSURF .LT. -81.) THEN
       DB = -81.   
      ELSE
       DB = TSURF
      ENDIF

C     SETTING 3 PARAMETERS FOR WETAIR, SINCE RH IS KNOWN (SEE WETAIR LISTING)  
      WB=0.D0
      DP=999.  
C     BP CALCULATED FROM ALTITUDE USING THE STANDARD ATMOSPHERE
C     EQUATIONS FROM SUBROUTINE DRYAIR    (TRACY ET AL,1972)
      PSTD=101325.  
      BP = PSTD*((1.-(0.0065*ALTT/288.))**(1./0.190284)) 

      CALL WETAIR(DB,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,  
     * DENAIR,CP,WTRPOT)
      SURFVD = VD 

C     CALCULATING AIR VAPOR DENSITY
      RH = RELHUM
      DB = TAIR
      CALL WETAIR(DB,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,
     * DENAIR,CP,WTRPOT)
      AIRVD = VD   

C     COMPUTE SURFACE WATER LOSS BASED ON EFFECT WET AREA/UNIT SURFACE AREA (M^2)
C     CONVERTING TO THE FRACTION OF A UNIT SURFACE AREA THAT IS WET
      if(sat.eq.2)then
       EFFSUR = 1. ! assume wet surface for water budget calcs
      else
       EFFSUR = PTWET/100.
      endif
C     AMOUNT OF WATER EVAPORATED FROM THE SURFACE (KG/S)
      WATER = EFFSUR * HD *(SURFVD - AIRVD)  

C     FROM DRYAIR: LATENT HEAT OF VAPORIZATION 
c      HTOVPR = 2.5012E+06 - 2.3787E+03 * TAIR
      if(TSURF.gt.0)then
       HTOVPR=2500.8-2.36*TSURF+0.0016*TSURF**2-0.00006*TSURF**3 
      else
       HTOVPR=2834.1-0.29*TSURF-0.004*TSURF**2 
      endif
      HTOVPR=HTOVPR*1000
c     convert qevap to cal/min/cm2 for dsub 
      if(sat.eq.2)then
       QEVAP = WATER * HTOVPR ! keep in SI for water budget calcs  
      else
       QEVAP = WATER * HTOVPR / 4.184 * 60. / 10000.  
      endif
C     KG/S TO G/S 
      GWSURF  = WATER * 1000. 
c     don't lose water if heat is just going into melting snow
      if(TSURF.le.0)then
       gwsurf=0.D0
      endif
      
      RETURN
      END
