      SUBROUTINE POLINT(XA,YA,N,X,Y,DY)
      IMPLICIT NONE
      INTEGER N,NMAX
      DOUBLE PRECISION DY,X,Y,XA(N),YA(N)
      PARAMETER (NMAX=10)
      INTEGER I,M,NS
      DOUBLE PRECISION DEN,DIF,DIFT,HO,HP,W,C(NMAX),D(NMAX)
      NS=1
      DIF=ABS(X-XA(1))
      DO 11 I=1,N
       DIFT=ABS(X-XA(I))
       IF (DIFT.LT.DIF) THEN
        NS=I
        DIF=DIFT
       ENDIF
       C(I)=YA(I)
       D(I)=YA(I)
11    CONTINUE
      Y=YA(NS)
      NS=NS-1
      DO 13 M=1,N-1
       DO 12 I=1,N-M
        HO=XA(I)-X
        HP=XA(I+M)-X
        W=C(I+1)-D(I)
        DEN=HO-HP
C       IF(DEN.EQ.0.)PAUSE 'FAILURE IN POLINT'
        DEN=W/DEN
        D(I)=HP*DEN
        C(I)=HO*DEN
12     CONTINUE
       IF (2*NS.LT.N-M)THEN
        DY=C(NS+1)
       ELSE
        DY=D(NS)
        NS=NS-1
       ENDIF
       Y=Y+DY
13    CONTINUE
      RETURN
      END
C  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE |2J.012U2.
