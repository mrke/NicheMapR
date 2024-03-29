      SUBROUTINE PHOTOPERIOD(LAT,DAYCOUNT,IHOUR,DOY,STARTDAY,PHOTOSTART,
     & PHOTOFINISH,PHOTODIRS,PHOTODIRF,DAYLENGTHSTART,DAYLENGTHFINISH,
     & BREEDSTART,BREEDEND,BREED,LAMBDA,PREVDAYLENGTH,LENGTHDAYDIR
     &,LENGTHDAY)

C     NICHEMAPR: SOFTWARE FOR BIOPHYSICAL MECHANISTIC NICHE MODELLING

C     COPYRIGHT (C) 2018 MICHAEL R. KEARNEY AND WARREN P. PORTER

C     THIS PROGRAM IS FREE SOFTWARE: YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C     IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C     THE FREE SOFTWARE FOUNDATION, EITHER VERSION 3 OF THE LICENSE, OR (AT
C      YOUR OPTION) ANY LATER VERSION.

C     THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT
C     WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C     MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. SEE THE GNU
C     GENERAL PUBLIC LICENSE FOR MORE DETAILS.

C     YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C     ALONG WITH THIS PROGRAM. IF NOT, SEE HTTP://WWW.GNU.ORG/LICENSES/.

C     COMPUTATION OF DAY LENGTHAND PHOTOPERIOD DIRECTIONAL SHIFT TO COMPUTE
C     BREEDING PHENOLOGY FOR DEB MODEL

      USE AACOMMONDAT
      
      IMPLICIT NONE

      DOUBLE PRECISION DAYLENGTHFINISH,DAYLENGTHSTART,FLENGTH,LAMBDA,LAT
      DOUBLE PRECISION LENGTHDAY,LENGTHDAYDIR,MLENGTH,PI,PREVDAYLENGTH

      INTEGER BREED,BREEDEND,BREEDSTART,DAYCOUNT,DOY,I,IHOUR,PHOTODIRF
      INTEGER PHOTODIRS,PHOTOFINISH,PHOTOSTART,STARTDAY

      PI=3.14159265

C     SETTING STARTDAY FOR PHOTOPERIOD-SPECIFIED BREEDING PERIODS WHEN STARTDAY SET TO ZERO BY USER

      IF((DAYCOUNT.EQ.1).AND.(IHOUR.LT.2).AND.(PHOTOSTART.EQ.5).AND.
     &  (STARTDAY.EQ.0))THEN
       DO 844 I=1,365
        MLENGTH=1.-TAN(LAT*PI/180.)*TAN(23.439*PI/180.*COS(0.0172*
     &   (DBLE(I)-1)))
        IF(MLENGTH.GT.2.)THEN
         MLENGTH=2.
        ENDIF
        IF(MLENGTH.LT.0.)THEN
         MLENGTH=0.
        ENDIF
        FLENGTH=ACOS(1.-MLENGTH)/(2.*PI)*360./180.
        IF(((FLENGTH*24.).GT.DAYLENGTHSTART).AND.((FLENGTH*24.)
     &      .GT.LENGTHDAY).AND.(PHOTODIRS.EQ.1))THEN
         STARTDAY=I
        ENDIF
        IF(((FLENGTH*24.).GT.DAYLENGTHSTART).AND.((FLENGTH*24.)
     &      .LT.LENGTHDAY).AND.(PHOTODIRS.EQ.0))THEN
         STARTDAY=I
        ENDIF
        LENGTHDAY=FLENGTH*24.
844    CONTINUE
      ENDIF

      IF((DAYCOUNT.EQ.1).AND.(IHOUR.LT.2).AND.(PHOTOSTART.LT.5).AND.
     &  (STARTDAY.EQ.0))THEN
       IF(LAT.LT.0)THEN !SOUTHERN HEMISPHERE
        IF(PHOTOSTART.EQ.1)THEN
         STARTDAY=357
        ENDIF
        IF(PHOTOSTART.EQ.2)THEN
        STARTDAY=80
        ENDIF
        IF(PHOTOSTART.EQ.3)THEN
         STARTDAY=173
        ENDIF
        IF(PHOTOSTART.EQ.4)THEN
         STARTDAY=266
        ENDIF
       ELSE !NORTHERN HEMISPHERE
        IF(PHOTOSTART.EQ.1)THEN
         STARTDAY=173
        ENDIF
        IF(PHOTOSTART.EQ.2)THEN
        STARTDAY=266
        ENDIF
        IF(PHOTOSTART.EQ.3)THEN
         STARTDAY=357
        ENDIF
        IF(PHOTOSTART.EQ.4)THEN
         STARTDAY=80
        ENDIF
       ENDIF
      ENDIF
C     SAVING LAST DAY'S LENGTH FOR CHECK FOR INCREASING OR DECREASING DAY LENGTH
      IF(DAYCOUNT.GT.1)THEN
       PREVDAYLENGTH=LENGTHDAY
      ELSE
       IF(DOY.EQ.1)THEN
        MLENGTH=1.-TAN(LAT*PI/180.)*TAN(23.439*PI/180.*COS(0.0172*
     &   (365.-1.)))
        IF(MLENGTH.GT.2.)THEN
         MLENGTH=2.
        ENDIF
        IF(MLENGTH.LT.0.)THEN
         MLENGTH=0.
        ENDIF
        FLENGTH=ACOS(1.-MLENGTH)/(2.*PI)*360./180.
        PREVDAYLENGTH=FLENGTH*24.
       ELSE
        MLENGTH=1.-TAN(LAT*PI/180.)*TAN(23.439*PI/180.*COS(0.0172*
     &   (REAL(DOY,8)-2.)))
        IF(MLENGTH.GT.2.)THEN
         MLENGTH=2.
        ENDIF
        IF(MLENGTH.LT.0.)THEN
         MLENGTH=0.
        ENDIF
        FLENGTH=ACOS(1.-MLENGTH)/(2.*PI)*360./180.
        PREVDAYLENGTH=FLENGTH*24.
       ENDIF
      ENDIF
C     CALCULATED LENGTH OF CURRENT DAY
      MLENGTH=1.-TAN(LAT*PI/180.)*TAN(23.439*PI/180.*COS(0.0172*
     & (REAL(DOY,8)-1.)))
      IF(MLENGTH.GT.2.)THEN
       MLENGTH=2.
      ENDIF
      IF(MLENGTH.LT.0.)THEN
       MLENGTH=0.
      ENDIF
      FLENGTH=ACOS(1.-MLENGTH)/(2.*PI)*360./180.
      LENGTHDAY=FLENGTH*24.

C     CHECK FOR INCREASING OR DECREASING DAY LENGTH
      IF(DAYCOUNT.GT.1)THEN
       IF(PREVDAYLENGTH.LT.LENGTHDAY)THEN
        LENGTHDAYDIR=1
       ELSE
        LENGTHDAYDIR=0
       ENDIF
      ENDIF

      IF((DAYCOUNT.EQ.1).AND.(IHOUR.LT.2))THEN
       IF(PHOTOSTART.EQ.0)THEN
        LAMBDA=1.
       ENDIF
C      SET LAMBDA FOR PHOTOPERIOD-SPECIFIC BREEDING SEASON
       IF((PHOTOSTART.EQ.5).OR.(PHOTOFINISH.EQ.5))THEN
        LAMBDA=1 ! DEFAULT STARTING VALUE
        BREEDSTART=0
        BREEDEND=0
C       NOW CHECK FOR WHEN BREEDING SEASON WOULD START
        BREED=0
        DO 845 I=1,365
         IF(I.EQ.1)THEN ! GET DAY LENGTH OF LAST DAY OF YEAR
          MLENGTH=1.-TAN(LAT*PI/180.)*TAN(23.439*PI/180.*
     &    COS(0.0172*(365.-1.)))
          IF(MLENGTH.GT.2.)THEN
           MLENGTH=2.
          ENDIF
          IF(MLENGTH.LT.0.)THEN
           MLENGTH=0.
          ENDIF
          FLENGTH=ACOS(1.-MLENGTH)/(2.*PI)*360./180.
          LENGTHDAY=FLENGTH*24.
         ENDIF
         MLENGTH=1.-TAN(LAT*PI/180.)*TAN(23.439*PI/180.*
     &   COS(0.0172*(DBLE(I)-1)))
         IF(MLENGTH.GT.2.)THEN
          MLENGTH=2.
         ENDIF
         IF(MLENGTH.LT.0.)THEN
          MLENGTH=0.
         ENDIF
         FLENGTH=ACOS(1.-MLENGTH)/(2.*PI)*360./180.
         IF(((FLENGTH*24.).GT.DAYLENGTHSTART).AND.(BREED.EQ.0).AND.
     &(LENGTHDAY.LE.DAYLENGTHSTART).AND.((FLENGTH*24.).GT.LENGTHDAY)
     &.AND.(PHOTODIRS.EQ.1))THEN
          BREEDSTART=I
         ENDIF
         IF(((FLENGTH*24.).LT.DAYLENGTHSTART).AND.(BREED.EQ.0).AND.
     &(LENGTHDAY.GE.DAYLENGTHSTART).AND.((FLENGTH*24.).LT.LENGTHDAY)
     &.AND.(PHOTODIRS.EQ.0))THEN
          BREEDSTART=I
         ENDIF
         LENGTHDAY=FLENGTH*24.
845     CONTINUE

C       NOW CHECK FOR WHEN BREEDING SEASON WOULD END
        BREED=1
        DO 846 I=1,365
         IF(I.EQ.1)THEN ! GET DAY LENGTH OF LAST DAY OF YEAR
          MLENGTH=1.-TAN(LAT*PI/180.)*TAN(23.439*PI/180.*
     &COS(0.0172*(365.-1.)))
         IF(MLENGTH.GT.2.)THEN
          MLENGTH=2.
         ENDIF
         IF(MLENGTH.LT.0.)THEN
          MLENGTH=0.
         ENDIF
         FLENGTH=ACOS(1.-MLENGTH)/(2.*PI)*360./180.
         LENGTHDAY=FLENGTH*24.
        ENDIF
         MLENGTH=1-TAN(LAT*PI/180.)*TAN(23.439*PI/180.*
     &COS(0.0172*(DBLE(I)-1)))
         IF(MLENGTH.GT.2.)THEN
          MLENGTH=2.
         ENDIF
         IF(MLENGTH.LT.0.)THEN
          MLENGTH=0.
         ENDIF
         FLENGTH=ACOS(1.-MLENGTH)/(2.*PI)*360./180.
         IF(((FLENGTH*24.).GT.DAYLENGTHFINISH).AND.(BREED.EQ.1).AND.
     &(LENGTHDAY.LE.DAYLENGTHFINISH).AND.((FLENGTH*24.).GT.LENGTHDAY)
     &.AND.(PHOTODIRF.EQ.1))THEN
          BREEDEND=I
         ENDIF
         IF(((FLENGTH*24.).LT.DAYLENGTHFINISH).AND.(BREED.EQ.1).AND.
     &(LENGTHDAY.GE.DAYLENGTHFINISH).AND.((FLENGTH*24.).LT.LENGTHDAY)
     &.AND.(PHOTODIRF.EQ.0))THEN
          BREEDEND=I
         ENDIF
         LENGTHDAY=FLENGTH*24.
846     CONTINUE

C       NOW CHECK FOR SOME COMBINATION OF SEASONAL AND PHOTOPERIOD CUES
        IF(LAT.LT.0)THEN
C        SOUTHERN HEMISPHERE
         IF(PHOTOSTART.EQ.1)THEN
          BREEDSTART=357
         ENDIF
         IF(PHOTOSTART.EQ.2)THEN
          BREEDSTART=80
         ENDIF
         IF(PHOTOSTART.EQ.3)THEN
          BREEDSTART=173
         ENDIF
         IF(PHOTOSTART.EQ.4)THEN
          BREEDSTART=266
         ENDIF
         IF(PHOTOFINISH.EQ.1)THEN
          BREEDEND=357
         ENDIF
         IF(PHOTOFINISH.EQ.2)THEN
          BREEDEND=80
         ENDIF
         IF(PHOTOFINISH.EQ.3)THEN
          BREEDEND=173
         ENDIF
         IF(PHOTOFINISH.EQ.4)THEN
          BREEDEND=266
         ENDIF
        ELSE
C        NORTHERN HEMISPHERE
         IF(PHOTOSTART.EQ.1)THEN
          BREEDSTART=173
         ENDIF
         IF(PHOTOSTART.EQ.2)THEN
          BREEDSTART=266
         ENDIF
         IF(PHOTOSTART.EQ.3)THEN
          BREEDSTART=357
         ENDIF
         IF(PHOTOSTART.EQ.4)THEN
          BREEDSTART=80
         ENDIF
         IF(PHOTOFINISH.EQ.1)THEN
          BREEDEND=173
         ENDIF
         IF(PHOTOFINISH.EQ.2)THEN
          BREEDEND=266
         ENDIF
         IF(PHOTOFINISH.EQ.3)THEN
          BREEDEND=357
         ENDIF
         IF(PHOTOFINISH.EQ.4)THEN
          BREEDEND=80
         ENDIF
        ENDIF
        IF((BREEDSTART.GT.0).AND.(BREEDEND.GT.0))THEN
         IF(BREEDEND-BREEDSTART.LE.0)THEN
          LAMBDA=DBLE((365-BREEDSTART+BREEDEND)/365)
         ELSE
          LAMBDA=DBLE((BREEDEND-BREEDSTART)/365)
         ENDIF
        ENDIF
       ENDIF

C      GETTING LAMBDA FOR SEASON-SPECIFIED BREEDING
       IF((PHOTOSTART.NE.0).AND.(PHOTOSTART.NE.5).AND.
     &   (PHOTOFINISH.NE.5))THEN
        IF((PHOTOSTART.EQ.1).AND.(PHOTOFINISH.EQ.2))THEN
         LAMBDA=3./12.
        ENDIF
        IF((PHOTOSTART.EQ.1).AND.(PHOTOFINISH.EQ.3))THEN
         LAMBDA=6./12.
        ENDIF
        IF((PHOTOSTART.EQ.1).AND.(PHOTOFINISH.EQ.4))THEN
         LAMBDA=9./12.
        ENDIF
        IF((PHOTOSTART.EQ.2).AND.(PHOTOFINISH.EQ.3))THEN
         LAMBDA=3./12.
        ENDIF
        IF((PHOTOSTART.EQ.2).AND.(PHOTOFINISH.EQ.4))THEN
         LAMBDA=6./12.
        ENDIF
        IF((PHOTOSTART.EQ.2).AND.(PHOTOFINISH.EQ.1))THEN
         LAMBDA=9./12.
        ENDIF
        IF((PHOTOSTART.EQ.3).AND.(PHOTOFINISH.EQ.4))THEN
         LAMBDA=3./12.
        ENDIF
        IF((PHOTOSTART.EQ.3).AND.(PHOTOFINISH.EQ.1))THEN
         LAMBDA=6./12.
        ENDIF
        IF((PHOTOSTART.EQ.3).AND.(PHOTOFINISH.EQ.2))THEN
         LAMBDA=9./12.
        ENDIF
        IF((PHOTOSTART.EQ.4).AND.(PHOTOFINISH.EQ.1))THEN
         LAMBDA=3./12.
        ENDIF
        IF((PHOTOSTART.EQ.4).AND.(PHOTOFINISH.EQ.2))THEN
         LAMBDA=6./12.
        ENDIF
        IF((PHOTOSTART.EQ.4).AND.(PHOTOFINISH.EQ.3))THEN
         LAMBDA=9./12.
        ENDIF
       ENDIF
      ENDIF


      RETURN
      END