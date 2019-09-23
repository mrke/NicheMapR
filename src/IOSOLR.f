       SUBROUTINE IOSOLR

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

C     THIS SUBROUTINE CALLED BY SOLRAD SETS UP I/O FOR SOLRAD
      use commondat
      IMPLICIT NONE
      double precision ALAT,ALMINT,ALONC,ALONG,ALREF,ALT,ALTT
      double precision AMINUT,AMULT,AZMUTH
      double precision CMH2O,DELONG
      double precision HEMIS,julstnd
      double precision MAXSHD
      double precision PATMOS,PRESS,PSTD,PTWET,PUNSH,REFL
      double precision SABNEW,SHDMAX,SHDMIN,SHAYD,SLOPE,SUMAX,SUMIN
      double precision TANNUL,TIMCOR,TIMAXS,TIMINS
      double precision TSRHR,TSNHR,USRHYT,intercept
      double precision tannul2,rainfall,Thconduct,Density,Spheat
      double precision snownode,minsnow,maxsnode1,snode,snowage,prevden,
     &daysincesnow,lastday,undercatch,rainmeltf,densfun,snowcond
      double precision snowdens,snowmelt,snowtemp,cursnow,grasshade

      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,CONS,IMN,I91,I92
     & ,I93,I94,I95,I96,I97,I98,I99,I100,I101
      INTEGER IALT,IEND,IEP,IPINT,ISTART
      INTEGER IUV,IEMON,ISMON,JULNUM,DOY,NOSCAT
      INTEGER NUMRUN,NUMTYPS,IDA,IDAYST
      INTEGER microdaily

      CHARACTER(1) ANS1,ANS2,ANS3,ANS4,ANS5,ANS6,ANS7,ANS9,ANS10
      CHARACTER(1) ANS11,ANS13,ANS14,ANS15,ANS16,ANS17,ANS18
      CHARACTER(1) SINE,SNSLOP
      CHARACTER(80) LABL1,LABL2,LABL3
      CHARACTER(12) FNAME

      DIMENSION TIMINS(4),TIMAXS(4),snownode(10),snode(10),densfun(4)
      DIMENSION julstnd(2),Thconduct(30),Density(30),Spheat(30)

c     Variable substrate properties, times & locations

      COMMON/LABEL/LABL1,LABL2,LABL3,FNAME,SINE,ANS14,SNSLOP
      COMMON/LABELS/ANS16,ANS17,ANS18
      COMMON/ENDS/JULSTND
      COMMON/DAYS/TANNUL
      COMMON/DAYSS/TIMINS,TIMAXS
      COMMON/SNOWPRED/snowtemp,snowdens,snowmelt,snownode,minsnow
     &,maxsnode1,snode,cursnow,daysincesnow,lastday,undercatch,rainmeltf
     &,densfun,snowcond,intercept,snowage,prevden,grasshade
      COMMON/WIOCONS/PUNSH,ALAT,AMULT,PRESS,CMH2O,REFL,ALONC,TIMCOR,
     * AZMUTH,SLOPE,TSNHR,TSRHR,Hemis
      COMMON/WIOCONS2/IPINT,NOSCAT,IUV,IALT,IDAYST,IDA,IEP,ISTART,IEND
      COMMON/GROUND/SHAYD,ALTT,MAXSHD,SABNEW,PTWET,rainfall
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101
      COMMON/HYTE/USRHYT
      COMMON/DAYJUL/JULNUM,DOY
      COMMON/WICHDAY/NUMRUN
c     Variable soil properties data for Dsub
      COMMON/SOYVAR1/Numtyps
      COMMON/SOYVAR2/Thconduct,Density,Spheat

c     for NicheMapR
      COMMON/LATLONGS/AMINUT,ALONG,ALMINT,ALREF
      COMMON/DAILY/microdaily
      common/deep/tannul2

C     DEFINING DEFAULT PARAMETERS FROM DELETED USER INTERFACE
      CONS = 0
      ISTART = 1
      IEND = 24
C     TOTAL SOLAR INTEGRATED
      ANS1 = 'T'
C     GLOBAL SOLAR OUT (BEAM & DIFFUSE OUT SEPARATELY, NOW)
      ANS2 = 'G'
C     ROUTINE SCATTERING CALCULATION, NO USE OF SUB. GAMMA
      ANS3 = 'R'
C     NO PLOT FORMAT (USED FOR SPECTRAL SOLAR OUTPUT)
      ANS4 = 'N'
C     NORTHERN OR SOUTHERN HEMISPHERE (N/S)
      ANS5 = 'N'
C     USE STANDARD ATM. PRESSURE FOR ELEVATION
      ANS6 = 'Y'
C     12 MIDDAYS (JULIAN) OF EACH MONTH SUPPLIED
      ANS7 = 'E'
C     HORIZONTAL SLOPE?
      ANS9 = 'N'
C     CALCULATE HOURLY TEMPERATURES FROM TMAX, TMIN
      ANS10 = 'Y'
C     DATAFILE ENTRY FOR TMAX, TMIN TAIRS
      ANS11 = 'D'
C     CORRECTNESS OF DATA QUERY FOR AIR TEMP'S FROM KEYBOARD
      ANS13 = 'Y'
C     MICROMET OUTPUT FORMAT
      ANS14 = 'M'
C     DON'T CHOOSE START & END HOURS; 0 - 24 HOURS USED AS DEFAULT
      ANS15 = 'N'
C     RELATIVE HUMIDITY CONSTANT?
      ANS16 = 'N'
C     CLOUD COVER CONSTANT?
      ANS17 = 'N'
C     WIND SPEED CONSTANT?
      ANS18 = 'N'
C     SUN ON HORIZONTAL SURFACE?
      SNSLOP = 'H'
C     HOUR(S)  AFTER SUNRISE WHEN TAIR = TAIR, MIN
      TSRHR = TIMINS(1)!0.0
C     HOUR(S) AFTER SOLAR NOON WHEN TAIR = TAIR, MAX
      TSNHR = TIMAXS(1)!1.0

      IPINT = 0
      NOSCAT = 1
      PUNSH = 1

C     LOCALITY DATA

      IF(IDAYST.GT.1)THEN
C      CORRECT THE VALUES FOR COUNTERS
       ISMON = IDAYST
       IEMON = IDA
       IDAYST = ISMON - (IDAYST - 1)
       IDA = IEMON - ISMON + 1
      ENDIF

      SHDMIN = MINSHADES(1)
      SHDMAX = MAXSHADES(1)

C     ALTT(M) -> ALT(KM)
      ALT = ALTT/1000.

C     FINDING MEAN ANNUAL TEMPERATURE
      SUMAX = 0.0
      SUMIN = 0.0
      DO 399 IMN = 1,IDA
       SUMAX = SUMAX + TMAXX(IMN)
       SUMIN = SUMIN + TMINN(IMN)
399   CONTINUE
c     kearney added this for daily calculations
      if(microdaily.eq.1)then
       TANNUL = TANNUL2
      else
       TANNUL = ((SUMAX+SUMIN)/2.)/IDA
      endif

      ALAT = ABS(ALAT)
      AMINUT = ABS(AMINUT)

C     PUTTING DEGREES AND MINUTES INTO DECIMAL FORM
      ALAT = ALAT + AMINUT/60.

C     PUTTING DEGREES AND MINUTES INTO DECIMAL FORM
      ALONG = ALONG + ALMINT/60.

C     CORRECTING FOR TIME ZONES
      DELONG = (ALONG - ALREF)/15.

C     ALTITUDE

C     CONVERTING TO NEAREST INTEGER
      IALT = NINT(ALT)

C     ATMOSPHERIC PRESSURE
      PSTD=101325.
      PATMOS=PSTD*((1.-(.0065*ALT/288.))**(1./.190284))

C     CONVERTING PASCALS TO MILLIBARS FOR SOLRAD
      PRESS = PATMOS/100.
      SINE = ANS10
C     STORE TIME OF MIN, MAX VALUES RELATIVE TO SUNRISE &
C     SOLAR NOON
      TIMINS(1)=TSRHR
      TIMAXS(1)=TSNHR

      RETURN
      END

