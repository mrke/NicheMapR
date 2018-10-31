      SUBROUTINE VSINE(IVAR,VMIN,VMAX,TIMIN,TIMAX)
      use commondat
      implicit none
      
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

C     SUB. SINE PREDICTS VARIABLE ON THE HOUR GIVEN THE MAXIMUM   
C     AND MINIMUM VALUES FOR A SPECIFIED DAY. 

C     ******THIS PROGRAM COMPUTES TIME AND VALUES FOR MICROMET. ***** 

      double precision TIMARY,TIMSR,TIMSS,TIMAX,TIMIN,tmin2,tmax2,TIMTMX
     &  ,time,VMAX,VMIN,XA,YA,itair,icld,iwind,irelhum,vinit,tmin,tmax
      double precision DEPS,vave
      double precision vsmin,vsmax,slope1,slope2,slope3,frmin,timec
      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,IVAR,jsum,itest1,i
      integer microdaily,cnt,non,iday,itest2,itime,jct,I91,I92,I93
     & ,I94,I95,I96 
      INTEGER I97,I98,I99,I100,I101
      
      DIMENSION TIMARY (25),XA(25),YA(25)
      DIMENSION DEPS(21)

      COMMON /WSINE/TIMSR,TIMSS,TIMTMX,TMIN,TMAX,TMIN2,TMAX2
      COMMON /SOILND/NON
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101
      common/init/itair,icld,iwind,irelhum,iday
      COMMON/DAILY/microdaily
      COMMON/dataky/DEPS,CNT
     
C     THIS LOOP READS IN THE LABEL AND APPROPRIATE DATA FOR EACH DAY.   
C     ****  TIME IS READ IN IN HUNDRED HOURS  EG. 2:45 PM = 1445 HOURS.     
C     VMIN - MINIMUM VALUE
C     VMAX - MAXIMUM VALUE     
C     TIMSR - TIME OF SUNRISE (FROM SOLRAD)   
C     TIMSS - TIME OF SUNSET (FROM SOLRAD)
C     TIMIN - TIME OF MINIMUM.  USUALLY ABOUT SUNRISE
C     TIMAX - TIME OF MAXIMUM.  USUALLY 1300 HOURS    
C     IF IN THE SAME TIME ZONE AS THE REFERENCE MERIDIAN
      vinit=0
C     TIMSR,TIMSS,TIMAX PROVIDED BY SOLRAD 
      VAVE=(VMAX+VMIN)/2   
      if((microdaily.eq.1).and.(iday.gt.1))then
       IF (IVAR .EQ. 1) THEN
        vinit=irelhum
       Endif
       IF (IVAR .EQ. 2) THEN
        vinit=icld
       Endif
       IF (IVAR .EQ. 3) THEN
        vinit=iwind
       ENDIF
      endif

C     X = TIME, Y = DEPENDENT VARIABLE, E.G. TEMPERATURE

C     SETTING UP X,Y ARRAYS
      VSMIN = VMIN + 0.01*VAVE
      VSMAX = VMAX - 0.01*VAVE  
      if (timin .lt. timax) then
c      morning minimum, afternoon maximum 
       ITEST1 = int(TIMIN)/100
       ITEST2 = int(TIMAX)/100
C      Slope from midnight to morning minimum       
       SLOPE1 = (VAVE - VSMIN)/(100.- TIMIN)
C      Slope from morning minimum to aft. maximum        
       SLOPE2 = (VSMIN - VSMAX)/(TIMIN - TIMAX) 
C      Slope from aft. max. to midnight ave. for the day        
       SLOPE3 = (VSMAX - VAVE)/(TIMAX - 2400.)
      else
c      morning maximum, afternoon minimum 
       ITEST1 = int(TIMAX)/100
       ITEST2 = int(TIMIN)/100
C      Slope from midnight to morning maximum         
       SLOPE1 = (VAVE - VSMAX)/(100. - TIMAX)
C      Slope from morning maximum to aft. minimum           
       SLOPE2 = (VSMAX - VSMIN)/(TIMAX - TIMIN) 
C      Slope from aft. min. to midnight ave. for the day        
       SLOPE3 = (VSMIN - VAVE)/(TIMIN - 2400.) 
      endif  

      DO 12 I = 1,25
C      FINDING THE INDEPENDENT VALUE, X      
       XA(I) = FLOAT(I)*100. - 100.
C      CONVERTING FROM MILITARY TIME TO BIOME TIME   
       TIME = XA(I)
       ITIME=int(TIME/100.)
       FRMIN=(TIME/100)-ITIME     
       ITIME=ITIME*60
       FRMIN=FRMIN*100.   
       TIMEC=ITIME+FRMIN  
       TIMARY(I) = TIMEC  

C      FINDING THE DEPENDENT VALUE, Y
       IF (I .EQ. 1) THEN
        if((microdaily.eq.1).and.(iday.gt.1))then
         YA(I) = VINIT
        else
         YA(I) = VAVE 
        endif
        GO TO 12
       ENDIF
        
       IF (I .LT. ITEST1) THEN
C       ON FIRST SLOPE: y2 = y1 - m * (x1-x2)
        YA(I) = YA(1) - (SLOPE1 * (XA(1) -XA(I)))
        GO TO 12
       ENDIF

C      AT FIRST MINIMUM OR MAXIMUM?
       IF (I .EQ. ITEST1) THEN
C       AT MINIMUM OR MAXIMUM
        if (timin .lt. timax) then
         YA(I) = VSMIN 
        else
         YA(I) = VSMAX
        endif   
        GO TO 12
       ENDIF
         
       IF (I .LT. ITEST2) THEN 
C       ON SECOND SLOPE 
        if (timin .lt. timax) then
         YA(I) = VSMIN - (SLOPE2 * (TIMIN - XA(I)))
        else 
         YA(I) = VSMAX - (SLOPE2 * (TIMAX - XA(I)))   
C        Correcting for occasional overshoot when times not exactly 
C        on the hour
         if (ya(i) .gt. vsmax) then
          ya(i) = vsmax
         endif              
        endif   
        GO TO 12
       ENDIF

C      AT SECOND MIN OR MAX?
       IF (I .EQ. ITEST2) THEN
        if (timin .lt. timax) then
         YA(I) = VSMAX 
        else
         YA(I) = VSMIN
        endif
        GO TO 12
       ENDIF
        
       IF (I .GT. ITEST2) THEN
C       ON SLOPE 3 
        if (timin .lt. timax) then
         YA(I) = VSMAX - (SLOPE3 * (TIMAX - XA(I)))  
        else 
         YA(I) = VSMIN - (SLOPE3 * (TIMIN - XA(I)))  
        endif  
       ENDIF 
12    CONTINUE 
C     END OF LOOP LOOKING FOR DEPENDENT VALUE

C     SET POSSIBLE NEG. VALUES TO POSITIVE VALUES
      jsum = 0

c     Checking for negative values at inflection points
      DO 100 JCT = 1,25
       if(ya(jct).lt.0.0)then
        jsum = jsum + 1
        if(ya(jct+1).gt.0.0)then
c        correct the neg. value using prior and post pos. values
         ya(jct) = (ya(jct-1) + ya(jct+1))/2.
        else
c        more than 1 negative value sequentially
         ya(jct) = abs(ya(jct))
        endif
       endif
100   CONTINUE

      if(ivar .eq. 1)then
       irelhum=YA(25)
       do 2011 i=1,25
        RELS((CNT-1)*25+i)=YA(i)
2011   continue
      endif
      if(ivar .eq. 2)then
       icld=YA(25)
       do 2012 i=1,25
        CLDS((CNT-1)*25+i)=YA(i)
2012   continue
      endif
      if(ivar .eq. 3)then
       iwind=YA(25)
       do 2013 i=1,25
        VELS((CNT-1)*25+i)=YA(i)
2013   continue
      endif

      RETURN
      END   
