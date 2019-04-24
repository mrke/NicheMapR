      SUBROUTINE GAMMA(TAU1,GAMR,GAML,SBAR)
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
           
C     COMPUTATIONS OF THE RADIATION SCATTERED BY A PLANE PARALLEL   
C     HOMOGENEOUS ATMOSPHERE USING THE METHOD OF THE X AND Y FUNCTIONS. 
C     ALL COMPUTATI6NS ARE PERFORMED IN DOUBLE PRECISION.....   
C     NPNCH SHOULD BE NON-ZERO IF THE SCATTERING FUNCTIONS ARE  
C     REQUIRED 

      DOUBLE PRECISION CNU1,CNU2,CNU3,CNU4,CU3,CU4,SBAR,
     1TAU1,CHX(101),CHY(101),CFA(3),AMU(101),GAML(101),GAMR(101),   
     2X1(101),Y1(101),X2(101),Y2(101),AIL(101),AI(101),XA(101),XB(101) 
      INTEGER I,NST,NTR
      AMU(1) = 0.0D0    
      DO 500 I = 2,101  
      AMU(I) = 0.01D0 * ( I-1)  
  500 CONTINUE  
C     COMPUTATIONS OF X SUB L,Y SUB L, X SUB R , Y SUB R .......
      CFA(1) = 0.75D0   
      CFA(2) = -0.75D0  
      CFA(3) = 0.0D0    
      NST = 111 
      CALL DCHXY(TAU1,CFA,NST,CHX,CHY,NTR)  
      DO 520 I = 1,101  
      X1(I) = CHX(I)    
      Y1(I) = CHY(I)    
  520 CONTINUE  
      CFA(1) = 0.375D0  
      CFA(2) = -0.375D0 
      NST = 0   
      CALL DCHXY(TAU1,CFA,NST,CHX,CHY,NTR)  
      DO 540 I=1,101    
      X2(I) = CHX(I)    
      Y2(I) = CHY(I)    
  540 CONTINUE  
C     COMPUTATIONS OF THE MOMENTS AND THE CONSTANTS....     
      AIL(1) = 0.01D0/3.0D0 
      CNU1 = 4.0D0 * AIL(1) 
      CNU2 = 2.0D0 * AIL(1) 
      DO 580 I = 2,100,2
      AIL(I) = CNU1     
      AIL(I+1) = CNU2   
  580 CONTINUE  
      AIL(101) = AIL(1) 
      DO 590 I = 1,101  
      XA(I) = 0.0D0     
      XB(I) = 0.0D0     
  590 CONTINUE  
      DO 630 I = 1,101  
      CNU1 = AIL(I) * X1(I) * AMU(I)
      XA(1) = XA(1) + CNU1  
      XA(2) = XA(2) + CNU1 * AMU(I) 
      CNU1 = AIL(I) * Y1(I) * AMU(I)
      XA(3) = XA(3) + CNU1  
      XA(4) = XA(4) + CNU1 * AMU(I) 
      CNU1 = AIL(I) * X2(I) 
      XB(1) = XB(1) + CNU1  
      CNU1 = CNU1 * AMU(I)  
      XB(2) = XB(2) + CNU1  
      CNU1 = CNU1 * AMU(I)  
      XB(3) = XB(3) + CNU1  
      XB(4) = XB(4) + CNU1 * AMU(I) 
      CNU1 = AIL(I) * Y2(I) 
      XB(5) = XB(5) + CNU1  
      CNU1 = CNU1 * AMU(I)  
      XB(6) = XB(6) + CNU1  
      CNU1 = CNU1 * AMU(I)  
      XB(7) = XB(7) + CNU1  
      XB(8) = XB(8) + CNU1 * AMU(I) 
  630 CONTINUE  
      AI(1) = XB(1) + XB(5) - 8.0D0/3.0D0   
      AI(2) = XB(2) + XB(6) 
      AI(3) = XB(3) + XB(7) 
      AI(4) = XB(1) - XB(5) - 8.0D0/3.0D0   
      AI(5) = XB(2) - XB(6) 
      AI(6) = XB(3) - XB(7) 
      AI(7) = XB(4) - XB(8) 
      AI(8) = XA(1) + XA(3) 
      AI(9) = XA(2) + XA(4) 
      AI(10) = XA(1) - XA(3)
      AI(11) = XA(2) - XA(4)
      AI(12) = (AI(1)-AI(3))/((AI(4)-AI(6))*TAU1+2.0D0*(AI(5)-AI(7)))   
      AI(13)=1.0D0/(AI(4)*AI(10) - AI(5) * AI(11))  
      AI(14)=1.0D0/(AI(1)*AI(8)-AI(2)*AI(9)-2.0D0*AI(12)*(AI(5)*AI(8)-  
     1AI(4)*AI(9)))     
      AI(15)=2.0D0 * (AI(8) * AI(10) - AI(9) *AI(11))   
      AI(16) = AI(13) * AI(15)  
      AI(17) = AI(14) * AI(15)  
      CNU1 = 0.5D0*(AI(16) - AI(17))
      CNU2 = 0.5D0*( AI(16) + AI(17))   
      AI(15) = AI(13)*(AI(5)*AI(8)-AI(4)*AI(9)) 
      AI(16)=AI(14)*(AI(2)*AI(10)-AI(1)*AI(11)-2.0D0*AI(12)*(AI(4)*AI(10
     1)-AI(5)*AI(11)))  
      CNU3 = 0.5D0*(AI(15)-AI(16))  
      CNU4 = 0.5D0*(AI(15)+AI(16))  
      AI(15) = AI(13)*(AI(2)*AI(10)-AI(1)*AI(11))   
      AI(16) = AI(14)*(AI(5)*AI(8)-AI(4)*AI(9)) 
      CU3 = 0.5D0*(AI(15)-AI(16))   
      CU4 = 0.5D0*(AI(15)+AI(16))   
      AI(15) = AI(14)*(AI(1)*AI(8)-AI(2)*AI(9)) 
      SBAR=1.0D0-0.375D0 * AI(12) * (AI(4)-AI(6)) * 
     1((CNU2 - CNU1 ) * AI(8) + (CU4 - CU3 ) *AI(2) -AI(15) * AI(6))    
      AI(20)=0.375D0*AI(12)*(CNU2-CNU1)*(AI(4)-AI(6))   
      AI(21) =0.375D0 * AI(12)*(AI(4)-AI(6))    
      AI(22) = AI(21) * ( CU4 - CU3 )   
      AI(23) = AI(21)*AI(15)
      DO 680 I = 1,101  
      GAML(I) = AI(20)*(X1(I) + Y1(I))  
      GAMR(I) =AI(22)*(X2(I) + Y2(I)) -AMU(I)*AI(23)*(X2(I)-Y2(I))  
  680 CONTINUE  
      RETURN
      END   