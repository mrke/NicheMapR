      SUBROUTINE TRAPZD(A,B,S,N)
      INTEGER N
      DOUBLE PRECISION A,B,S,FUNC
      EXTERNAL FUNC
      INTEGER IT,J
      DOUBLE PRECISION DEL,SUM,TNM,X
      IF (N.EQ.1) THEN
       S=0.5*(B-A)*(FUNC(A)+FUNC(B))
      ELSE
       IT=2**(N-2)
       TNM=IT
       DEL=(B-A)/TNM
       X=A+0.5*DEL
       SUM=0.
       DO 11 J=1,IT
        SUM=SUM+FUNC(X)
        X=X+DEL
11     CONTINUE
       S=0.5*(S+(B-A)*SUM/TNM)
      ENDIF
      RETURN
      END
C  (C) COPR. 1986-92 NUMERICAL RECIPES SOFTWARE |2J.012U2.