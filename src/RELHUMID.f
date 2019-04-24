      Subroutine RELHUMID

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
      implicit none
      double precision ALT,ALTT,BP,CP,DB,DENAIR,DP,E,ESAT,PATMOS,DENSTY
      double precision RH,RHLOCL,THCOND,HTOVPR,TCOEFF,GGROUP,WB,DIFVPR
      double precision SHAYD,SIOUT,MAXSHD,SABNEW,PTWET,VISKIN
      double precision T,WC,C,DEP,OUT,RW,TVINC,TVIR,VAPREF,VD,VISDYN
      double precision WTRPOT,rainfall

      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,ITEST,IOUT
      integer I91,I92,I93,I94,I95,I96
      INTEGER I97,I98,I99,I100,I101
      
      COMMON/GROUND/SHAYD,ALTT,MAXSHD,SABNEW,PTWET,rainfall
      COMMON/SIUNIT/SIOUT(10)
      COMMON/LOCLHUM/RHLOCL
      COMMON/ARRAY/T(30),WC(20),C(20),DEP(30),OUT(101),IOUT(100),
     1 ITEST(23)
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101

c     setting variable values for sub's. dryair, wetair
      WB = 0.
      DP = 999.
c     altitude, alt, known, but not barometric pressure, bp
      ALT = ALTT
      BP = 0.
c     Dry bulb temperature is the 2m Tair
      DB = OUT(2)

c     call dryair to get atmospheric pressure, patmos, from altitude, alt, etc. for reference height
      call DRYAIR(DB,BP,ALT,PATMOS,DENSTY,VISDYN,VISKIN,DIFVPR,
     *THCOND,HTOVPR,TCOEFF,GGROUP)

      RH = SIOUT(3)
      if(RH.gt.100.)then
            RH= 100.
      endif
c     get the vapor pressure, e, at Ta, 2m
      call WETAIR(DB,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,CP,
     &  WTRPOT)
      vapref = e
c     get the saturation vapor pressure, esat, at Tlocal
      RH = 100.
      DB = SIOUT(2)
      call WETAIR(DB,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,CP,
     *  WTRPOT)
c     Definition of relative humidity using the vapor density at reference height
      RHLOCL = (vapref/esat)* 100.

      if(RHLOCL.gt.100.)then
        RHLOCL = 100.
      endif
      if(RHLOCL.lt.0.000)then
        RHLOCL = 0.01
      endif

      return
      end