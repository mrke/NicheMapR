      FUNCTION TAB(NAME,AIND)

C     NicheMapR: software for biophysical mechanistic niche modelling

C     Copyright (C) 2018 Michael R. Kearney and Warren P. Porter

c     This program is free software: you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as published by
c     the FDOYree Software Foundation, either version 3 of the License, or (at
c      your option) any later version.

c     This program is distributed in the hope that it will be useful, but
c     WITHOUT ANY WARRANTY; without even the implied warranty of
c     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c     General Public License for more details.

c     You should have received a copy of the GNU General Public License
c     along with this program. If not, see http://www.gnu.org/licenses/.

c     This function does linear interpolation of intermediate values for
c     intermediate times when the integrator uses intermediate times &
c     values to solve the equations.
      IMPLICIT NONE
      double precision AIND,CONST,TD,TI,TAB
      INTEGER I,I1,IA,IFINAL,ILOCT,ISTART,JULNUM,DOY,errcount,maxerr
      INTEGER errout,numrun
      CHARACTER(3) INAME,NAME,SYMBOL
      COMMON/TABLE/TI(211),TD(211)
      COMMON/TABLE2/ILOCT(21)
      COMMON/CARRAY/INAME(20),SYMBOL(23)
      COMMON/DAYJUL/JULNUM,DOY
      common/errormsg/errout,maxerr,errcount
      COMMON/WICHDAY/NUMRUN

C     NAME = TABLE NAME
C     AIND = TIME
C     TI = INDEPENDENT VARIABLE (TIME)
C     TD = DEPENDENT VARIABLE (e.g. Qsolar, Tair, Wind, etc.)
c     TI(istart) = 0.00  (minutes @ start of day)
c     TI(ifinal) = 1440. (minutes @ end of day)
      TAB=0.D+0
      
C     MATCHING SUPPLIED NAME (NAME) WITH REFERENCE TABLE LIST (INAME)
      DO 10 I=1,20
      IA=I
      IF(INAME(I).EQ.NAME) GOTO 11
   10 CONTINUE
      GOTO 500
   11 CONTINUE
C     ERROR CHECK FOR TIME REQUESTED WITHIN DATA BOUNDS
      ISTART=ILOCT(IA)
      IA=IA+1
      IFINAL=ILOCT(IA)-1
      IF(AIND.LT.TI(ISTART)) GOTO 505
      IF(AIND.GT.TI(IFINAL)) GOTO 510
C
C     PERFORM THE LINEAR INTERPOLATION
C
      IF(AIND.EQ.TI(1)) GOTO 210
      DO 100 I=ISTART,IFINAL
c      Search the times for time greater than current integrator time
       IF(AIND.GT. TI(I)) GOTO 100
c      Set the counter for the time value below current integrator time
       I1=I-1
c      Compute the location in the interval between higher and lower time
       CONST=(TI(I1)-AIND)/(TI(I1)-TI(I))
c      Interpolate linearly the dependent value & return
       TAB=TD(I1)-CONST*(TD(I1)-TD(I))
       GOTO 200
  100 CONTINUE
  200 RETURN

  210 TAB=TD(ISTART)
      RETURN

  500 WRITE(6,501) NAME
  501 FORMAT('0****TABLE',A3,' HAS NEVER BEEN DEFINED. PROGRAM TERMINATE
     1D****')
      errcount=maxerr+1
      numrun=2
      DOY=julnum
      RETURN

  505 CONTINUE
  510 WRITE(6,511) AIND,NAME
  511 FORMAT('0****THE INDEPENDENT VARIABLE OF VALUE',E12.4,
     1  ' OF TABLE ',A3,' IS OUTSIDE THE RANGE GIVEN',
     2 '     PROGRAM TERMINATED')
      errcount=maxerr+1
      numrun=2
      DOY=julnum
      END
