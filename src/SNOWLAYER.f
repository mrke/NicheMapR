      subroutine snowlayer
      use commondat
      IMPLICIT NONE
      EXTERNAL WETAIR

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

C     Computes snow layer and thermal properties

      DOUBLE PRECISION daysincesnow,lastday,Thconduct,Density,Spheat
      DOUBLE PRECISION densfun,maxsnode1,minsnow,intercept
      DOUBLE PRECISION siout,snownode,snode,undercatch,rainmeltf,TT
      DOUBLE PRECISION DENDAY,SPDAY,TKDAY,T,WC,C,DEP,OUT,snowcond
      double precision snowdens,snowmelt,snowtemp,cursnow,snowage
     & ,prevden,grasshade
     
      INTEGER I,JULNUM,DOY,Numtyps,ITEST,NON,SNON,methour
      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,IOUT,maxsnode
      INTEGER I91,I92,I93,I94,I95,I96,I97,I98,I99,I100,I101
      DIMENSION DENDAY(30),SPDAY(30),TKDAY(30),snode(10),densfun(4)
      DIMENSION snownode(10),Thconduct(30),Density(30),Spheat(30),TT(30)

      COMMON/SOYVAR1/Numtyps
      COMMON/SOYVAR2/Thconduct,Density,Spheat
      COMMON/SOYFILS/DENDAY,SPDAY,TKDAY
      COMMON/SOILND/NON,SNON
      COMMON/DAYJUL/JULNUM,DOY
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101
      COMMON/SNOWPRED/snowtemp,snowdens,snowmelt,snownode,minsnow
     &,maxsnode1,snode,cursnow,daysincesnow,lastday,undercatch,rainmeltf
     &,densfun,snowcond,intercept,snowage,prevden,grasshade
      COMMON/SIUNIT/SIOUT(10)
      COMMON/ARRAY/T(30),WC(20),C(20),DEP(30),OUT(101),IOUT(100),
     1 ITEST(23)

      maxsnode=0
      if(cursnow.gt.300)then
       maxsnode1=0.D0
      endif

      methour=int(SIOUT(1)/60)+1+24*(DOY-1)

      if(cursnow.lt.minsnow)then ! get rid of snow
       maxsnode=0
       maxsnode1=0.D0
       do 34 i=1,8
        snode(i)=0.D0
34     continue
       snowhr(methour)=0.D0
        do 52 i=1,8 ! set temperature of snow nodes to that of soil node 1
         snode(i)=0.
         t(i)=t(1)
         tt(i)=tt(1)
52      continue
       cursnow=0.D0
       goto 900
      endif

      if(snowhr(methour).ge.minsnow)then
       do 1 i=1,8
        if(snowhr(methour).gt.snownode(i))then
         maxsnode=i
        endif
1      continue

       if(maxsnode.gt.7)then
        maxsnode=7
       endif

c     now build up the snow nodes accordingly but start from the bottom
       do 3 i=1,maxsnode
        snode(i+(8-maxsnode))=snownode(i)
3      continue

c     check if snow depth is going down a node
       if(maxsnode.lt.maxsnode1)then
        do 5 i=1,8-maxsnode
         snode(i)=0.D0
         t(i)=t(1)
5       continue
       endif

c     current snow depth in cm
       cursnow=snowhr(methour)
       maxsnode1=real(maxsnode,8)

      else
       do 4 i=1,8
        snode(i)=0.D0
        maxsnode=0
        cursnow=snowhr(methour)
        maxsnode1=real(maxsnode,8)
4      continue
      endif

  900 CONTINUE
      RETURN
      END
