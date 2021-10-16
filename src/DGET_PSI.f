      SUBROUTINE DGET_PSI(N,A,MPSI,DMPSI,RPAR,IPAR)

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

C     EQUATIONS TO COMPUTE LIQUID WATER EXCHANGE WITH SOIL FOR EGG OR SKIN

      IMPLICIT NONE
      EXTERNAL ISPLINE
      
      INTEGER IPAR,N
      DOUBLE PRECISION ISPLINE,M,PSI_E,DMPSI,PSI_S,RPAR,D_PSI_E,D_M,PI
      DOUBLE PRECISION M_S,K_S,B,M_A,MPSI,A
      DOUBLE PRECISION A_S,K_E,SPEC_HYD,K_SAT,P_E
      DIMENSION DMPSI(N),MPSI(N),RPAR(34),IPAR(31)
      
      A_S=RPAR(1)
      K_E=RPAR(2)
      SPEC_HYD=RPAR(3)
      K_SAT=RPAR(4)
      P_E=RPAR(5)
      B=RPAR(6)
      M_A=RPAR(7)
      PSI_S=RPAR(8)
      A=MPSI(1)! % S, TIME
      M=MPSI(2)! % KG, EGG MASS
      PSI_E=MPSI(3)! % J/KG, EGG WATER POTENTIAL
      PI=3.14159265

      K_S=K_SAT*(P_E/PSI_S)**(2+3/B)                                ! equation 9.2 from Campbell and Norman 2001
      M_S=(PSI_S-PSI_E)/(1./(A_S*K_E)+1./((2.*PI*A_S)**(1./2.)*K_S))  ! DARCY'S LAW 
      D_M=M_S-M_A                                                     ! KG/S, CHANGE IN MASS
      D_PSI_E=D_M/(SPEC_HYD*M)                                        ! J/KG/S, change in egg water potential

      DMPSI(1)=1.0D+00
      DMPSI(2)=D_M
      DMPSI(3)=D_PSI_E

      RETURN
      END