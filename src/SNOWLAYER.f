      subroutine snowlayer
c      use commondat
      IMPLICIT NONE
c      EXTERNAL WETAIR

C     Michael Kearney 2012
C     Computes snow layer and thermal properties

      REAL daysincesnow,lastday,Thconduct,Density,Spheat
      REAL snowtemp,snowdens,snowmelt,densfun,maxsnode1,minsnow
      REAL siout,snownode,snode,cursnow,undercatch,rainmeltf,TT
      REAL DENDAY,SPDAY,TKDAY,T,WC,C,DEP,IOUT,OUT

      INTEGER DAYCT,I,JULNUM,MOY,Numtyps,ITEST,NON,SNON,methour
      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,maxsnode
      INTEGER I91,I92,I93,I94,I95,I96,I97,I98,I99,I100,I101
C    Day's soil properties
      DIMENSION DENDAY(30),SPDAY(30),TKDAY(30),snode(8),densfun(2)
      DIMENSION snownode(8),Thconduct(30),Density(30),Spheat(30),TT(30)

      COMMON/SOYVAR1/Numtyps
      COMMON/SOYVAR2/Thconduct,Density,Spheat
      COMMON/SOYFILS/DENDAY,SPDAY,TKDAY
      COMMON/SOILND/NON,SNON
      COMMON/DAYJUL/JULNUM,MOY
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101
      COMMON/SNOWPRED/snowtemp,snowdens,snowmelt,snownode,minsnow
     &,maxsnode1,snode,cursnow,daysincesnow,lastday,undercatch,rainmeltf
     &,densfun
      COMMON/SIUNIT/SIOUT(10)
      COMMON/ARRAY/T(30),WC(20),C(20),DEP(30),IOUT(100),
     1 OUT(101),ITEST(23)

      DATA DAYCT/1/

      if(cursnow.gt.300)then
       maxsnode1=0.
      endif

      methour=int(SIOUT(1)/60)+1+24*(moy-1)

C     if((cursnow.gt.2).and.(cursnow.lt.2.5))then
C      maxsnode=0
C      maxsnode1=0.
C      do 33 i=1,8
C       snode(i)=0.
C33     continue
C      snode(8)=1.5
C      cursnow=1.5
C      snowhr(methour)=1.5
C       do 51 i=1,8
Cc       snode(i)=0.
C        t(i)=t(1)
C        tt(i)=tt(1)
C51      continue
C      goto 900
C     endif

      if(cursnow.lt.minsnow)then
       maxsnode=0
       maxsnode1=0.
       do 34 i=1,8
        snode(i)=0.
34     continue
       snowhr(methour)=0.
        do 52 i=1,8
         snode(i)=0.
         t(i)=t(1)
         tt(i)=tt(1)
52      continue
       cursnow=0.
       goto 900
      endif

c    check that max snow depth is not less than current snow depth
c    if(snowhr(methour).gt.snownode(8))then
c     snownode(8)=snowhr(methour)
c    endif

      if(snowhr(methour).ge.minsnow)then
       do 1 i=1,8
        if(snowhr(methour).gt.snownode(i))then
         maxsnode=i
        endif
1      continue

       if(maxsnode.gt.7)then
        maxsnode=7
       endif

c    now buildup the snow nodes accordingly but start from the bottom
       do 3 i=1,maxsnode
        snode(i+(8-maxsnode))=snownode(i)
3      continue

c    check if snow depth is going down a node
       if(maxsnode.lt.maxsnode1)then
        do 5 i=1,8-maxsnode
         snode(i)=0.
         t(i)=t(1)
5       continue
       endif

c    current snow depth in cm
       cursnow=snowhr(methour)
       maxsnode1=real(maxsnode)

      else
       do 4 i=1,8
        snode(i)=0.
        maxsnode=0
        cursnow=snowhr(methour)
        maxsnode1=real(maxsnode)
c      maxsnode1=maxsnode
4      continue
      endif

  900 CONTINUE


      RETURN
      END
