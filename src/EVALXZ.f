      SUBROUTINE EVALXZ
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

C     SUBROUTINE EVALXZ EVALUATES THE VALUE OF THE FUNCTION DTEMP/DTIME
C     FROM THE RUNGE KUTTA ALOGRITHM FOR SOLVING ORDINARY DIFFERENTIAL
C     EQUATIONS, THE ODES BEING CODED INTO DSUB.  THE MAXIMUM
C     NUMBER OF DIFFERENTIAL EQUATIONS ALLOWED IS 40.
      DIMENSION YII(40)

      double precision G,H,X,Y,YC,YII,XMAX,ERR,DX,YP,YI,TNEW,TERR,ERROR
     & ,f,dummy,moist,ep
      INTEGER I,N,NN,NFACTR,NCOND,IPRINT
      COMMON/NONSCR/N,NN,X,XMAX,DX,ERR,H,NCOND,NFACTR,IPRINT
      COMMON/WORK/ERROR(40),TERR(40),TNEW(40),YI(40),YP(40),YC(40),
     1 F(6,40),G(40),Y(40)
     2  ,DUMMY(1160)
      common/moistcom/moist(10),ep
C          SAVE THE INITIAL VALUES READ INTO EVAL
      DO 5 I=1,N
    5  YII(I)=Y(I)
C          FIND THE DERVATIVEA AT XO AND T(I)
      CALL DSUB(X,Y,G)
      DO 10 I=1,N
      YC(I) = Y(I) +(H*G(I))/6.0
   10 Y(I) = YII(I) + (H*G(I))/2.0
      X =X +H/2.0
      CALL DSUB(X,Y,G)
      DO 11 I=1,N
      YC(I)= YC(I)+(H*G(I))/3.0
   11  Y(I) = YII(I) + (H*G(I))/2.0
      CALL DSUB(X,Y,G)
      DO 12 I=1,N
      YC(I) = YC(I)+(H*G(I))/3.0
   12  Y(I) = YII(I) +  H*G(I)
      X =X +H/2.0
      CALL DSUB(X,Y,G)
      DO 13 I=1,N
   13 YC(I) = YC(I) + (H*G(I))/6.0
      CALL DSUB(X,YC,G)
      RETURN
      END