      SUBROUTINE DCHXY (TAU1,CFA,NCASE,CHX,CHY,NOMITR) 
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
       
C     COMPUTATIONS OF THE X - AND Y - FUNCTIONS OF CHANDRASEKHAR
C     THIS SUBROUTINE CARRIES OUT ALL COMPUTATIONS IN DOUBLE PRECISION..
C     THE MAIN PROGRAM MUST CONTAIN A STATEMENT DECLARING TAU1,CFA(3),  
C     CHX(101),CHY(101) IN DOUBLE PRECISION...  
C     THE INPUT PARAMETERS ARE AS FOLLOWS ...   
C     TAU1   NORMAL OPTICAL THICKNESS OF THE ATMOSPHERE...  
C     CFA(1),CFA(2),CFA(3)   COEFFIECIENTS REPRESENTING THE CHARACTERI- 
C     STIC FUNCTION IN A POLYNOMIAL FORM HAVING THREE TERMS GIVEN   
C     BY ( A SUB J ) * ( MU**2*(J-1)).. J = 1,2,3,...   
C     THIS SUBROUTINE RETURNS THE VALUES OF X FUNCTION AND Y FUNCTION   
C     AT 101 VALUES OF MU GIVEN BY 0.00(0.01)1.00....   
C     THE COMPUTATIONS ARE STARTED WITH THE FOURTH APPROXIMATION
C     GIVEN BY S. CHANDRASEKHAR IN SEC.59 OF THE 'RADIATIVE TRANSFER'   
C     DOVER PUBLICATIONS,INC.,NEW YORK,1960.... 
C     THE VALUES SO OBTAINED ARE THEN ITERATED USING THE PROCEDURE  
C     DESCRIBED IN SEC.60.. 
C     THE ITERATION IS TERMINATED WHEN THE APPROPRIATELY CORRECTED  
C     SUCCESSIVELY ITERATED VALUES OF THE Y FUNCTIONS AGREE IN  
C     FOUR SIGNIFICANT FIGURES....  
C     IF NCASE IS NOT EQUAL TO ZERO, A CONSERVATIVE CASE IS ASSUMED,AND,
C     A STANDARD SOLUTION IS RETURNED.....  
C     THE PROGRAM IS TERMINATED IF TAU1 IS GREATER THAN 2.0, OR,
C     THE CHARACTERISTIC FUNCTION IS NEGATIVE FOR ANY VALUE OF MU,OR,   
C     THE INTEGRAL OF CHARACTERISTIC FUNCTION IS GREATER THAN  0.5....  

    5 FORMAT(//T5,'THE PROGRAM IS TERMINATED AS PSI(I) =',D15.6,T50,    
     1 ' FOR I = ',I5//)
   20 FORMAT(//T10,'NO COMPUTATIONS CAN BE DONE AS THE COEFFIECIENTS ARE
     1 = ',3F15.5//)    
   30 FORMAT(//T5,' THE PROGRAM IS TERMINATED BECAUSE TAU1 = ',D15.5//) 
   40 FORMAT(1H1)   
   52 FORMAT(T25,'OPTICAL THICKNESS TAU SUB 1 =',F15.5//)   
   55 FORMAT(/T5,'MU',T25,'X FN PREVIOUS',T45,'X FN CURRENT',T65,   
     1'Y FN PREVIOUS',T85,'Y FN CURRENT'/)  
   60 FORMAT ( 0PF8.2,8X,1P4D20.5 ) 
   62 FORMAT(//T5,'CFA(1) =',D15.6,T35,'CFA(2) =',D15.6,T65,'CFA(3) =', 
     1D15.6//)  
   63 FORMAT(//T5,'N',T15,'ROOT K',T35,'1.0/ROOT K',T55,'UMA',T75,'ACAP'
     2,T95,'1.0/LAMDA'//)   
c   64 FORMAT(I5,5D20.10)
   75 FORMAT(//T8,'MU',T16,'P(MU)',T30,'P(-MU)',T44,'W(MU)',T58,'C 0 (MU
     1)',T72,'C 0 (-MU)',T86,'C 1 (MU)',T100,'C 1(-MU)',T110,'CHECK'//) 
   80 FORMAT(0PF10.2,1P8D14.6)  
   90 FORMAT(/T5,'ITERATION NUMBER OF THE CASE MARKED CURRENT IS =',I5) 
   95 FORMAT(/T5,'ITERATION NO.=',I5,T35,'X SUB 0 =',D15.5,T70,'Y SUB 0 
     1=',D15.5,T100,' MOMENT CHECK =',D15.5)    
   96 FORMAT(//T5,'THE PROGRAM IS TERMINATED BECAUSE ROOTS CANNOT BE    
     1FOUND TO SATISFY THE CRITERION..'//)  
   97 FORMAT(//T5,' THE TROUBLE IS WITH ROOT NUMBER ',I5//) 
      DOUBLE PRECISION TAU1,CFA(3),PSI(101),AMU(101),XA(101),XB(101),  
     1UMA(5),ACAP(5),TEMX(8),TEMY(8),RTK(5),ALAM(5),FNPP(101),FNPN(101),   
     2FNC0(101),FNC1(101),FNX(101),FNY(101),TEMA,TEMB,TEMC,TEMD,
     3PERA,PERB,FNW(101),FMC0(101),FMC1(101),XC(101),XD(101),XE(101)    
      DOUBLE PRECISION CHXA(101),CHYA(101),CHX(101),CHY(101),DEXPI  
      EQUIVALENCE (FNPP,XC),(FNPN,XD),(FNW,XE)  
      EQUIVALENCE (CHXA,FMC0),(CHYA,FMC1)   
      INTEGER NC0(4,6),NC1(4,8) 
      INTEGER NOMITR,NCASE,J,I,NPRT,KMX,N1,IST,K,N,IC,N2
      DATA NC0/3,4,1,2,2,4,1,3,2,3,1,4,1,4,2,3,1,3,2,4,1,2,3,4/ 
      DATA NC1/2,3,4,1,1,3,4,2,1,2,4,3,1,2,3,4,4,1,2,3,3,1,2,4,2,1,3,4, 
     11,2,3,4/  
C     TERMINATE THE PROGRAM IF ALL THE CONDITIONS ARE NOT SATISFIED,AND,
C     COMPUTE MU,PSI(MU) AND WEIGHTS OF SIMPSON QUADRATURE...   
      perb=0.D0
      IF ( TAU1 .LE. 2.0D0 ) GO TO 100  
   99 WRITE(6,30) TAU1  
      STOP  
  100 IF ( TAU1 .LT. 0.0D0 ) GO TO 99   
      PERA = CFA(1) + CFA(2)/3.0D0 + 0.2D0 * CFA(3) 
      IF ( NCASE .NE. 0 ) PERA = 0.5D0  
      IF ( PERA .GE. 0.0D0 ) GO TO 110  
  105 CONTINUE  
      WRITE(6,20)(CFA(J),J=1,3)     
      STOP  
  110 IF ( PERA .GT. 0.5D0 ) GO TO 105  
      DO 120 I = 1,101  
       AMU(I) = (I-1) * 0.01D0   
       TEMA = AMU(I)**2  
       PSI(I) = CFA(1) + CFA(2)*TEMA + CFA(3)*(TEMA**2)  
       IF ( PSI(I) .GT. -1.0D-15 ) GO TO 120     
       WRITE(6,5)PSI(I),I
       WRITE(6,20)(CFA(J),J=1,3)     
       STOP  
  120 CONTINUE  
      XA(1) = 0.01D0/3.0D0  
      TEMA = 4.0D0 * XA(1)  
      TEMB = 2.0D0 * XA(1)  
      DO 130 I = 2,100,2
      XA(I) = TEMA  
      XA(I+1) = TEMB    
  130 CONTINUE  
      XA(101) = XA(1)   
C     SET NPRT = 0 FOR SUPPRESSING ALL INTERMEDIATE OUTPUT FROM THIS    
C     SUBROUTINE....    
      NPRT = 0  
C     COMPUTE ROOTS OF THE CHARACTERISTIC EQUATION...   
      IF ( NCASE .NE. 0 ) GO TO 140 
      KMX = 4   
      N1 = 0
  131 CONTINUE  
      UMA(1) = 0.96028985649753623D+00  
      UMA(2) = 0.79666647741362674D+00  
      UMA(3) = 0.52553240991632899D+00  
      UMA(4) = 0.18343464249564980D+00  
      ACAP(1) = 0.10122853629037626D+00 
      ACAP(2) = 0.22238103445337447D+00 
      ACAP(3) = 0.31370664587788729D+00 
      ACAP(4) = 0.36268378337836198D+00 
      IF ( N1 .EQ. 0 ) GO TO 220    
      GO TO 400 
  140 KMX = 5   
      UMA(1) = 0.97390652851717172D+00  
      UMA(2) = 0.86506336668898451D+00  
      UMA(3) = 0.67940956829902441D+00  
      UMA(4) = 0.43339539412924719D+00  
      UMA(5) = 0.14887433898163121D+00  
      ACAP(1) = 0.66671344308688138D-01 
      ACAP(2) = 0.14945134915058059D+00 
      ACAP(3) = 0.21908636251598204D+00 
      ACAP(4) = 0.26926671930999636D+00 
      ACAP(5) = 0.29552422471475287D+00 
  220 DO 225 I = 1,KMX  
       TEMX(I) = UMA(I)**2   
       TEMY(I) = CFA(1) + CFA(2)*TEMX(I) + CFA(3)*TEMX(I)**2 
       TEMY(I) = 2.0D0 * ACAP(I) * TEMY(I)   
  225 CONTINUE  
      IF ( NCASE .EQ. 0 ) GO TO 235 
      IST = 2   
      RTK(1) = 0.0D0    
      GO TO 238 
  235 IST = 1   
  238 DO 244 I = IST,KMX
       RTK(I) = (1.0D0 - TEMY(I))/TEMX(I)
       IF ( I .NE. 1 ) GO TO 242     
       TEMA = 1.0D0/UMA(1)**2
       IF ( RTK(1) .LT. TEMA ) GO TO 244 
       RTK(1) = 0.5D0 * TEMA 
       GO TO 244 
  242  TEMA = 1.0D0/UMA(I-1)**2  
       TEMB = 1.0D0/UMA(I)**2
       IF ( RTK(I) .GT. TEMA .AND. RTK(I) .LT. TEMB ) GO TO 244  
       RTK(I) = 0.5D0 * ( TEMA + TEMB )  
  244 CONTINUE  
      J = IST   
  245 IF ( J .NE. 1 ) GO TO 246     
      TEMA = 0.0D0  
      TEMB = 1.0D0/UMA(1)**2
      N1 = 0
      GO TO 250 
  246 TEMA = 1.0D0/UMA(J-1)**2  
      TEMB = 1.0D0/UMA(J)**2
      N1 = 0
  250 TEMC = 1.0D0  
      DO 260 I = 1,KMX  
       TEMC = TEMC - TEMY(I)/(1.0D0 - RTK(J)*TEMX(I))
  260 CONTINUE  
      TEMD = DABS (TEMC)
      IF ( TEMD .LT. 1.0D-14) GO TO 358 
      N1 = N1 + 1   
      IF ( N1 .GT. 50 ) GO TO 350   
      IF ( TEMC .GT. 0.0D0 ) TEMA = RTK(J)  
      IF ( TEMC .LT. 0.0D0 ) TEMB = RTK(J)  
      TEMD = 0.0D0  
      DO 270 I = 1,KMX  
       TEMD = TEMD - (TEMY(I)*TEMX(I))/((1.0D0-RTK(J)*TEMX(I))**2)   
  270 CONTINUE  
      TEMC = RTK(J) - TEMC/TEMD     
      IF ( TEMC .GT. TEMA ) GO TO 280   
  275 RTK(J) = 0.5D0 *  ( TEMA + TEMB ) 
      GO TO 250 
  280 IF ( TEMC .GE. TEMB ) GO TO 275   
      RTK(J) = TEMC     
      GO TO 250 
  350 WRITE(6,96)   
      WRITE(6,62) (CFA(I),I=1,3)    
      WRITE(6,97) J     
      STOP  
  358 J = J + 1 
      IF ( J .LE. KMX) GO TO 245    
      DO 360 I = 1,KMX  
       RTK(I) = DSQRT(RTK(I))
  360 CONTINUE  
      IF ( NCASE .EQ. 0 ) GO TO 400 
      N1 = 11   
      KMX = 4   
      DO 370 J = 1,KMX  
       RTK(J) = RTK(J+1) 
  370 CONTINUE  
      GO TO 131 
C     COMPUTE FUNCTIONS LAMDA, P AND W....  
  400 DO  410 J = 1,KMX 
       ALAM(J) = 1.0D0   
       DO 405 I = 1,KMX  
        ALAM(J) = ALAM(J)*((RTK(J)*UMA(I)+1)/(RTK(J)*UMA(I)-1))   
  405  CONTINUE  
       ALAM(J) = DEXP(-RTK(J)*TAU1)/ALAM(J)  
  410 CONTINUE  
      IF ( NPRT .EQ. 0 ) GO TO 418  
      WRITE(6,62)(CFA(J),J=1,3)     
      WRITE(6,52)TAU1   
      WRITE(6,63)   
      DO 414 J = 1,KMX  
       TEMA = 1.0D0/RTK(J)   
  414 CONTINUE  
  418 DO  420 I = 1,101 
       FNPP(I) = 1.0D0   
       FNPN(I) = 1.0D0   
       FNW(I) = 1.0D0    
       DO  419 J = 1,KMX 
        FNPP(I) = FNPP(I)*( AMU(I)/UMA(J)-1.0D0)  
        FNPN(I) = FNPN(I)*(-AMU(I)/UMA(J)-1.0D0)  
        FNW(I) = FNW(I) * ( 1.0D0 - RTK(J)**2 * AMU(I)**2)    
  419  CONTINUE  
  420 CONTINUE  
C     COMPUTE C SUB 0 AND C SUB 1 ......
      TEMX(1) = 1.0D0   
      TEMX(8) = 1.0D0   
      DO 440 K = 2,7    
       TEMX(K) = 1.0D0   
       DO 430 I = 1,2    
        N1 = NC0(I,K-1)   
        DO 425 J = 1,2    
         N2 = NC0(J+2,K-1) 
         TEMX(K) = TEMX(K) *((RTK(N1) + RTK(N2))/(RTK(N1) - RTK(N2)))  
  425   CONTINUE  
  430  CONTINUE  
       TEMX(K) = -TEMX(K)
  440 CONTINUE  
      DO 460 K = 1,4    
       TEMY(K) = 1.0D0   
       N2 = NC1(4,K)     
       DO 455 I = 1,3    
        N1 = NC1(I,K)     
        TEMY(K) = TEMY(K) *((RTK(N1) + RTK(N2))/(RTK(N1) - RTK(N2)))  
  455  CONTINUE  
  460 CONTINUE  
      DO 470 K = 5,8    
       TEMY(K) = 1.0D0   
       N1 = NC1(1,K)     
       DO 465 J = 1,3    
        N2 = NC1(J+1,K)   
        TEMY(K) = TEMY(K) *((RTK(N1) + RTK(N2))/(RTK(N1) - RTK(N2)))  
  465  CONTINUE  
       TEMY(K) = -TEMY(K)
  470 CONTINUE  
      DO 550 I = 1,101  
       TEMA = 1.0D0  
       TEMB = 1.0D0  
       DO 500 J = 1,4    
        TEMA = TEMA * ( 1.0D0 + RTK(J) *AMU(I))   
        TEMB = TEMB * ( 1.0D0 - RTK(J) *AMU(I))   
  500  CONTINUE  
       FNC0(I) = TEMA    
       FMC0(I) = TEMB    
       TEMA = 1.0D0  
       TEMB = 1.0D0  
       DO 510 J = 1,4    
        TEMA = TEMA *(1.0D0 - RTK(J)*AMU(I))*ALAM(J)  
        TEMB = TEMB *(1.0D0 + RTK(J)*AMU(I))*ALAM(J)  
  510  CONTINUE  
       FNC0(I) = FNC0(I) + TEMA  
       FMC0(I) = FMC0(I) + TEMB  
       IST = 2   
  520  TEMA = 1.0D0  
       TEMB = 1.0D0  
       DO 525 K = 1,2    
        N2 = NC0(K+2,IST-1)   
        TEMA = TEMA * (1.0D0-RTK(N2) * AMU(I))*ALAM(N2)   
        TEMB = TEMB * (1.0D0+RTK(N2) * AMU(I))*ALAM(N2)   
  525  CONTINUE  
       DO 530 J = 1,2    
        N1 = NC0(J,IST-1) 
        TEMA = TEMA * (1.0D0 + RTK(N1)*AMU(I))    
        TEMB = TEMB * (1.0D0 - RTK(N1)*AMU(I))    
  530  CONTINUE  
       FNC0(I) = FNC0(I) + TEMA * TEMX(IST)  
       FMC0(I) = FMC0(I) + TEMB * TEMX(IST)  
       IST = IST + 1     
       IF ( IST .LE. 7 ) GO TO 520   
  550 CONTINUE  
      DO 600 I = 1,101  
       FNC1(I) = 0.0D0   
       FMC1(I) = 0.0D0   
       IST = 1   
  555  N2 = NC1(4,IST)   
       TEMA =      ( 1.0D0 - RTK(N2)*AMU(I))*ALAM(N2)
       TEMB =      ( 1.0D0 + RTK(N2)*AMU(I))*ALAM(N2)
       DO 560 J = 1,3    
        N1 = NC1(J,IST)   
        TEMA = TEMA * ( 1.0D0 + RTK(N1)*AMU(I))   
        TEMB = TEMB * ( 1.0D0 - RTK(N1)*AMU(I))   
  560  CONTINUE  
       FNC1(I) = FNC1(I) + TEMY(IST) * TEMA  
       FMC1(I) = FMC1(I) + TEMY(IST) * TEMB  
       IST = IST + 1     
       IF ( IST .LE. 4 ) GO TO 555   
  575  N1 = NC1(1,IST)   
       TEMA  = 1.0D0 + RTK(N1)*AMU(I)
       TEMB = 1.0D0 - RTK(N1) *AMU(I)
       DO 580 J = 1,3    
        N2 = NC1(J+1,IST) 
        TEMA = TEMA*( 1.0D0 - RTK(N2)*AMU(I))*ALAM(N2)
        TEMB = TEMB*( 1.0D0 + RTK(N2)*AMU(I))*ALAM(N2)
  580  CONTINUE  
       FNC1(I) = FNC1(I) + TEMY(IST) * TEMA  
       FMC1(I) = FMC1(I) + TEMY(IST) * TEMB  
       IST = IST + 1     
       IF ( IST .LE. 8 ) GO TO 575   
       FNC1(I) = -FNC1(I)
       FMC1(I) = -FMC1(I)
  600 CONTINUE  
      IF ( NPRT .EQ. 0 ) GO TO 630  
      WRITE(6,40)   
      WRITE(6,62)(CFA(J),J=1,3)     
      WRITE(6,52)TAU1   
      WRITE(6,75)   
      DO 620 I = 1,101  
       TEMD = FNC0(I) * FMC0(I) - FNC1(I) * FMC1(I) - (FNC0(1)**2 -  
     1  FNC1(1)**2) * FNW(I)  
       WRITE(6,80)AMU(I),FNPP(I),FNPN(I),FNW(I),FNC0(I),FMC0(I), 
     1 FNC1(I),FMC1(I),TEMD  
  620 CONTINUE  
      WRITE(6,40)   
C     COMPUTE THE FOURTH APPROXIMATION OF X AND Y FUNCTIONS..   
  630 XB(1) = 0.0D0     
      IF ( TAU1 .EQ. 0.0D0 ) XB(1) = 1.0D0  
       DO 635 I = 2,101  
       XB(I) = DEXP(-TAU1/AMU(I))    
  635 CONTINUE  
      TEMA = 1.0D0/DSQRT(FNC0(1)**2-FNC1(1)**2) 
      DO 640 I = 1,101  
       TEMC = TEMA / FNW(I)  
       FNX(I) = (FNPN(I)*FMC0(I) - XB(I)*FNPP(I)*FNC1(I))*TEMC   
       FNY(I) = (XB(I)*FNPP(I)*FNC0(I) - FNPN(I)*FMC1(I))*TEMC   
  640 CONTINUE  
C     SET UP ITERATION TO OBTAIN APPROXIMATE SOLUTION...    
C     FNX(I) AND FNY(I)    PREVIOUS SET..   
C     CHXA(I) AND CHYA(I)   APPROXIMATE VALUES OBTAINED AFTER ITERATION.
C     CHX(I) AND CHY(I)   CURRENT SET...
      CHXA(1) = 1.0D0   
      CHYA(1) = XB(1)   
      NOMITR = 1
  650 CONTINUE  
      DO  900 I = 2,101 
       DO 750 IC = 1,101 
        XD(IC) = PSI(IC)*(FNX(I)*FNX(IC)-FNY(I)*FNY(IC))/(AMU(I)
     & +AMU(IC))    
        IF ( I .EQ. IC) GO TO 750     
        XE(IC) =PSI(IC)*(FNY(I)*FNX(IC)-FNX(I)*FNY(IC))/(AMU(I)-AMU(IC))  
  750  CONTINUE  
       IF ( I .GT. 3 ) GO TO 770     
C      LINEAR INTERPOLATION....  
       XE(I) = 0.5D0 * ( XE(I+1) + XE(I-1))  
       GO TO 820 
  770  IF ( I .GT. 5 ) GO TO 780     
C      INTERPOLATION USING EVERETT'S FORMULA TWO POINT ON EITHERSIDE... 
       XE(I)=0.0625D0*(9.0D0*(XE(I+1)+XE(I-1))-XE(I+3)-XE(I-3))  
       GO TO 820 
  780  IF ( I .GT. 96 ) GO TO 800    
C      INTERPOLATION USING EVERETT'S FORMULA THREE POINTS ON
C      EITHERSIDE.....   
       XE(I) = 3.0D0 * (XE(I+5) + XE(I-5)) + 150.0D0*(XE(I+1)+XE(I-1))
     1 - 25.0D0 * (XE(I+3)+XE(I-3)) 
       XE(I) = XE(I)/256.0D0 
       GO TO 820 
  800  XE(I) = 5.0D0*XE(I-1)+10.0D0*XE(I-3)+XE(I-5)-10.0D0*XE(I-2)-5.0D0 
     1 * XE(I-4)
  820  CHXA(I) = 0.0D0   
       CHYA(I) = 0.0D0   
       DO 850 IC = 1,101 
        CHXA(I)=CHXA(I)+XA(IC)*XD(IC) 
        CHYA (I) = CHYA(I) + XA(IC)* XE(IC)   
  850  CONTINUE  
       CHXA(I) = 1.0D0 + AMU(I) * CHXA(I)
       CHYA(I) = XB(I) + AMU(I) * CHYA(I)
  900 CONTINUE  
C     CORRECTION TO THE APPROXIMATION...
      IF ( NOMITR .GT. 1 ) GO TO 908
      IF ( TAU1 .EQ. 0.0D0 ) GO TO 917  
      TEMX(1) = -DEXPI(-TAU1)   
      DO 905 N = 2,7    
       TEMA = N - 1  
       TEMX(N) = (XB(101) - TAU1 * TEMX(N-1))/TEMA   
  905 CONTINUE  
      PERB = CFA(1)*(0.5D0-TEMX(3))+CFA(2)*(0.25D0-TEMX(5))+
     1 CFA(3)*(1.0D0/6.0D0 - TEMX(7))   
      PERB = 2.0D0 * PERB   
  908 TEMA = 0.0D0  
      TEMB = 0.0D0  
      DO 910 I = 1,101  
       TEMA = TEMA + CHXA(I) *PSI(I)*XA(I)   
       TEMB = TEMB + CHYA(I) *PSI(I)*XA(I)   
  910 CONTINUE  
      TEMC = (1.0D0 - 2.0D0 * PERA)/(1.0D0 - TEMA + TEMB )  
      TEMC = (1.0D0 - TEMA - TEMB - TEMC )/PERB 
  917 IF ( TAU1 .EQ. 0.0D0 ) TEMC = 0.0D0   
      DO  920 I = 1,101 
       TEMD = TEMC * AMU(I) * ( 1.0D0 - XB(I))   
       CHX(I) = CHXA(I) + TEMD   
       CHY(I) = CHYA(I) + TEMD   
  920 CONTINUE  
      IF ( NPRT .EQ. 0 ) GO TO 925  
      TEMC = TEMA**2 - 2.0D0*TEMA - TEMB**2 + 2.0D0*PERA    
      WRITE (6,95 ) NOMITR,TEMA,TEMB,TEMC   
  925 IF ( NOMITR .EQ. 1 ) GO TO 950
C      CHECK FOR CONVERGENCE...     
      DO 930 I = 2,101  
       TEMA = ( CHY(I) - FNY(I) ) / CHY(I)   
       TEMA = DABS(TEMA) 
       IF ( TEMA .LE. 2.0D-4 ) GO TO 930 
       GO TO 950 
  930 CONTINUE  
      GO TO 960 
  950 NOMITR = NOMITR + 1   
      IF ( NOMITR .GT. 15 ) GO TO 965   
      DO 955 I=1,101    
       FNX(I) = CHX(I)   
       FNY(I) = CHY(I)   
  955 CONTINUE  
      GO TO 650 
  960 IF ( NPRT .EQ. 0 ) GO TO 975  
  965 WRITE(6,40)   
      WRITE(6,52)TAU1   
      WRITE(6,62)(CFA(J),J=1,3)     
      WRITE(6,55)   
      WRITE(6,60)(AMU(I),FNX(I),CHX(I),FNY(I),CHY(I),I=1,101)   
      WRITE(6,90)NOMITR 
      TEMA = 0.0D0  
      TEMB = 0.0D0  
      DO 970 I = 1,101  
       TEMC = PSI(I) * XA(I) 
       TEMA = TEMA + CHX(I) * TEMC   
       TEMB = TEMB + CHY(I) * TEMC   
  970 CONTINUE  
      TEMC = TEMA**2 -2.0D0*TEMA -TEMB**2 + 2.0D0 *PERA     
      WRITE(6,95) NOMITR,TEMA,TEMB,TEMC 
      IF ( NOMITR .GT. 15 ) STOP    
      CONTINUE  
C     GENERATE STANDARD SOLUTION IF NCASE IS NOT EQUAL TO ZERO..
  975 IF ( NCASE .EQ. 0 ) GO TO 1000
      TEMA = 0.0D0  
      TEMB = 0.0D0  
      TEMC = 0.0D0  
      DO  980 I = 1,101 
       TEMD = PSI(I) * AMU(I) * XA(I)
       TEMA = TEMA + TEMD * CHX(I)   
       TEMB = TEMB + TEMD * CHY(I)   
       TEMC = TEMC + PSI(I) * CHY(I) * XA(I)     
  980 CONTINUE  
      TEMD = TEMC/(TEMA + TEMB)     
      DO 990 I = 1,101  
       TEMA = TEMD*AMU(I)*(CHX(I) + CHY(I))  
       CHX(I) = CHX(I) + TEMA
       CHY(I) = CHY(I) - TEMA
  990 CONTINUE  
 1000 RETURN
      END   
