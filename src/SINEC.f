      SUBROUTINE SINEC 
      use commondat
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

C     SUB. SINEC IS CALLED BY SUB. SOLRAD.
C     SINEC PREDICTS TEMPERATURE ON THE HOUR GIVEN THE AVERAGE MAXIMUM   
C     AND AVERAGE MINIMUM TEMPERATURES ** IN DEG. C ** FOR A SPECIFIED DAY. 

C     ******THIS PROGRAM GENERATES TIME AND TEMPERATURE FOR MICROMET. ***** 
      DOUBLE PRECISION TIMSR,TIMSS,TIMTMX,TMAX,TMIN
      DOUBLE PRECISION TANNUL,itair,icld,iwind,irelhum,tmin2,tmax2
      DOUBLE PRECISION DEPS
      DOUBLE PRECISION A,TSR,TREF,SS,SY,ZS,TSS,TAU,T,TI,E,X,Y,Z,FRMIN
      DOUBLE PRECISION TIN,TDS,TIMARY,DUMAIR,TAIRRY,TIMEC,TASUM

      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I,J,TIME,ITIME
      integer microdaily,cnt,IAIR,NON,iday,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101,I102,I103,I104,I105,I106,I107
      
      DIMENSION DEPS(21),TIMARY(25),DUMAIR(25),TAIRRY(25)

      COMMON/WSINE/TIMSR,TIMSS,TIMTMX,TMIN,TMAX,TMIN2,TMAX2
      COMMON/SOILND/NON
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101,I102,I103,I104,I105,I106,I107
      COMMON/DAYS/Tannul
      common/init/itair,icld,iwind,irelhum,iday
      COMMON/DAILY/microdaily
      COMMON/dataky/DEPS,CNT
     
C     THIS LOOP READS IN THE LABEL AND APPROPRIATE DATA FOR EACH DAY.   
C ****TIME IS READ IN IN HUNDRED HOURS  EG. 2:45 PM = 1445 HOURS.     
C     TMIN - MINIMUM TEMPERATURE   (C)     
C     TMAX - MAXIMUM TEMPERATURE   (C)     
C     TIMSR - TIME OF SUNRISE (FROM SOLRAD) = TIMIN   
C     TIMSS - TIME OF SUNSET (FROM SOLRAD)
C     TIMTMX - TIME OF MAXIMUM TEMPERATURE.  USUALLY 1300 HOURS = TIMAX   
C     IF IN THE SAME TIME ZONE AS THE REFERENCE MERIDIAN

C     TIMSR,TIMSS,TIMTMX PROVIDED BY SOLRAD 
      TMIN=TMIN+273.    
      TMAX=TMAX+273.
      TMIN2=TMIN2+273.
      TMAX2=TMAX2+273.    
      A=(TMAX-TMIN)/2.   
      TSR=TMIN
      TREF=(TIMTMX-TIMSR)/2.+TIMSR  
      SS=360.*(TIMSS-TREF)/(2.*(TIMTMX-TIMSR))  
      SY=SS/57.29577    
      ZS=SIN(SY)
      TSS=A*ZS+TMIN+A   
      TAU=3./((2400.-TIMSS)+TIMSR)  

C     THIS LOOP COMPUTES TEMPERATURE EACH HOUR OF THE DAY   
      DO 50 I=1,24  
         J=I+1  
         TIME=I*100
         IF (TIME-TIMSR) 20,21,11   
c          sunrise
   21      T=TMIN 
           T=T-273.   
           GO TO 10   
   11    IF (TIME-TIMSS) 15,15,13 
   20       TI=(2400.-TIMSS)+TIME
         if((microdaily.eq.1).and.(iday.gt.1))then
          T=((TMIN-273-ITAIR)/TIMSR)*TIME+ITAIR
         else
          E=TI*TAU   
          T=(TSS-TSR)/EXP(E)+TSR     
          T=T-273.
         endif       
         GO TO 10   
c        after sunset
   13    TI=(2400.-TIMSS)-(2400.-TIME)  
          E=TI*TAU   
           if(microdaily.eq.1)then
           T=(TSS-TMIN2)/EXP(E)+TMIN2 
           else
            T=(TSS-TSR)/EXP(E)+TSR  
           endif 
           T=T-273.   
           GO TO 10   
c        before or at sunset
   15    X=360.*(TIME-TREF)/(2.*(TIMTMX-TIMSR)) 
         Y=X/57.29577   
         Z=SIN(Y)   
         T=A*Z+TMIN+A   
         T=T-273.   
C         CONVERTING FROM MILITARY TIME TO BIOME TIME   
   10    ITIME=int(TIME)/100
         FRMIN=(TIME/100.)-ITIME     
         ITIME=ITIME*60
         FRMIN=FRMIN*100.   
         TIMEC=ITIME+FRMIN  
         TIMARY(J) = TIMEC  
         TAIRRY(J) = T  
50    CONTINUE  

C     WRITE START "CARD" FOR MICROMET DATA T(0)=T(1440)    
      TIMARY(1) = 0.0 
      if((microdaily.eq.1).and.(iday.gt.1))then 
       TAIRRY(1) = itair 
      else
       TAIRRY(1) = T
      endif     

C     PROCESSING AIR TEMPERATURES FOR SOIL INITIAL TEMPERATURES
      TASUM = 0.000
      DO 2005 IAIR = 1,25
       TASUM = TAIRRY(IAIR) + TASUM
2005  CONTINUE
C     MONTHLY AVERAGE
      TIN = TASUM/25.
C     ANNUAL TEMPERATURE AVERAGE
      TDS = TANNUL
 
c     writing initial soil temperature array from average of air temps
      DO 2010 IAIR = 1,NON
       DUMAIR(IAIR) = TIN
       TINS(IAIR,cnt)=TIN
2010  CONTINUE

      TDSS(CNT)=TDS

      do 2011 i=1,25
       TARS((CNT-1)*25+i)=TAIRRY(i)
2011  continue
      ITAIR=TAIRRY(25)

      RETURN
      END   
