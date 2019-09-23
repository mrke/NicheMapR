      SUBROUTINE infil(ha,moistt,ET,temp,depth,fl,sw,humid,potent,dt,
     &altt,rw,pc,rl,sp,r1,im,maxcount,leafpot,rootpot,trans)

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

C     Computes water infiltration and redistribution with evaporation, from
c     a bare soil surface, based on program 9.1 of Campbell 1985

      double precision A,B,C,F,P,Z,V,DP,W,WN,K,CP,WS,B1,N,N1
      double precision WD,GR,IM,SE,H,JV,DJ,PP,EP,MW,T,R,DV
     &,VP,KV,ECUR
      double precision RR,RS,PR,BZ,RW,PC,RL,PI,SP,R1,TP,PB,RB
      double precision KS,PE,BB,BD,LAI,L,moistt,HA,ET,temp,depth,humid
     & ,potent,SW,FL
      double precision altt,dpp,pstd,bp,ESAT,VD,RWW,TVIR,TVINC,DENAIR
     & ,CPP,WTRPOT
     & ,RH100,DB,rootpot,leafpot,trans,DD,WB
      double precision SL,PL,FF,E,XP,TR
      Integer M,I,count,maxcount,j,DT

      DIMENSION A(19),B(19),C(19),F(19),P(19),Z(19),V(19),DP(19),W(19)
      DIMENSION WN(19),K(19),CP(19),H(19),JV(19),DJ(19),temp(10)
     &,moistt(18),T(19),depth(10),humid(18),potent(18),PE(19),
     &KS(19),BB(19),PP(19),B1(19),N(19),N1(19),WS(19),DV(19),rootpot(18)
      DIMENSION RR(19),L(19),E(19),RS(19),PR(19),BZ(19),BD(19),DD(19)
      common/campbell/PE,KS,BB,BD,L,LAI,DD

c     P matric potential J/kg
c     Z depth nodes
c     W water content m3/m3
c     WN water content m3/m3
c     K hydraulic conductivity, kg s/m3
c     M number of elements
c     BD soil bulk density, Mg/m3
c     DD soil mineral density, Mg/m3
c     KS saturated conductivity, kg s/m3
c     PE air entry potential J/kg
c     BB soil 'b' parameter
c     RW=2.5E+10 ! resistance per unit length of root, m3 kg-1 s-1
c     PC=-1500. ! critical leaf water potential for stomatal closure, J kg-1
c     RL=2000000. ! resistance per unit length of leaf, m3 kg-1 s-1
c     SP=10. ! stability parameter, -
c     R1=0.001 ! root radius, m
c     IM=0.000001 ! maximum overall mass balance error allowed, kg

      PI=3.14159
      MW=0.018 ! molar mass of water, kg/mol
      R=8.310001 ! gas constant, J/mol/K
      GR=9.8 !gravitational constant m/s/s
     
      A(1:19)=0
      B=A
      C=A
      F=A
      P=A
      Z=A
      V=A
      DP=A
      W=A
      WN=A
      K=A
      CP=A
      H=A
      JV=A
      DJ=A
      B1=A
      N=A
      WS=A
      DV=A
      RR=A
      E=A
      RS=A
      PR=A
      BZ=A
      E=A
      M=18

      PE=ABS(PE)*(-1) !air entry potential J/kg
      PP=PE !initial water potential J/kg
      PP=ABS(PP)*(-1)
      WS=1-BD/DD !saturation water content m3/m3, assuming max density of 2.6 Mg/m3
      Z(M+1)=2 ! depth to lower boundary, m
      WD=1000. ! density of water kg/m3
      DV=0.000024 ! binary diffusion coefficient
  
      B1=1/BB
      N=2+3/BB
      N1=1-N

c     prep for wetair call      
      WB = 0.
      DPP = 999.
      PSTD=101325.
      BP=PSTD*((1.-(0.0065*ALTT/288.))**(1./0.190284))      

      j=2
      do 121 I=3,18
       if(MOD(I, 2).ne.0)then
        Z(I)=depth(j)/100
        j=j+1
       else
        Z(I)=Z(I-1)+(depth(J)/100-Z(I-1))/2
       endif
121   continue
      j=1
      do 120 I=1,19
       if(MOD(I, 2).ne.0)then
        T(I)=temp(j)
        j=j+1
       else
        T(I)=T(I-1)+(temp(j)-T(I-1))/2
       endif
120   continue

      T=T+273 ! temperature from deg C to Kelvin

      Z(2)=0
      Z(1)=0

c     # setting initial water content, m3/m3
      do 2 I=2,M
       WN(I)=moistt(i-1)
       P(I)=PE(I)*(WS(I)/WN(I))**BB(I)
       H(I)=exp(MW*P(I)/(R*T(I-1)))
       K(I)=KS(I)*(PE(I)/P(I))**N(I)
       W(I)=WN(I)
2     continue

      do 22 I=2,M
       V(I)=WD*(Z(I+1)-Z(I-1))/2
22    continue
      
c     # lower boundary condition set to saturated (stays constant)
      P(M+1)=PE(M)*(WS(M+1)/WS(M+1))**BB(M) ! potential
      H(M+1)=1. ! humidity
      W(M+1)=WS(M+1) ! water content
      WN(M+1)=WS(M+1) ! water content
      Z(1)=-1E10 ! depth at node 1, m
      Z(M+1)=1E20 ! depth at deepest node, m
      K(M+1)=KS(M)*(PE(M)/P(M+1))**N(M+1) ! lower boundary conductivity

c     initialize root water uptake variables
      do 98 I=2,M
       if(L(I).gt.0)then
        RR(I)=2*RW/(L(I)*(Z(I+1)-Z(I-1))) ! root resistance
        BZ(I)=(1-M)*LOG(PI*R1*R1*L(I))/(2*PI*L(I)*(Z(I+1)-Z(I-1)))
       else
        RR(I)=1E+20 ! root resistance
        BZ(I)=0.D0
       endif
98    continue

      P(1)=P(2)
      K(1)=0

c     evapotranspiration
      EP=exp(-0.82*LAI)*ET ! partition potential evaporation from potential evapotranspiration
      TP=ET-EP ! now get potenital transpiration

c     plant water uptake
      PB=0.D0 ! weighted mean soil water potential, J/kg
      RB=0.D0 ! weighted mean root soil root resistance, m4 /(s kg)
      PL=0.D0 ! leaf water potential, J/kg
      do 99 i=2,M
       RS(I)=BZ(I)/K(I) ! soil resistance
       PB=PB+P(I)/(RR(I)+RS(I))
       RB=RB+1/(RS(I)+RR(I))
99    continue
      PB=PB/RB
      RB=1/RB
      maxcount=500
      count=0
2080  continue
      IF(PL.gt.PB)then
       PL=PB-TP*(RL+RB)
      ENDIF
      XP=(PL/PC)**SP
      SL=TP*(RL+RB)*SP*XP/(PL*(1+XP)*(1+XP))-1.
      FF=PB-PL-TP*(RL+RB)/(1+XP)
      PL=PL-FF/SL
      count=count+1
      if((ABS(FF).gt.10).and.(count.lt.maxcount))then
       goto 2080
      endif
      TR=TP/(1+XP)
      do 100 I=2,M
       E(I)=(P(I)-PL-RL*TR)/(RR(I)+RS(I)) ! root water uptake
100   continue
      count=0
      
c     start of convergence loop ########################################
11    SE=0
      maxcount=500
      count=count+1
      do 3 I=2,M
       K(I)=KS(I)*(PE(I)/P(I))**N(I) ! conductivities for each node
3     continue
      JV(1)=EP*(H(2)-HA)/(1-HA)
      DJ(1)=EP*MW*H(2)/(R*T(I-1)*(1-HA))
      do 4 I=2,M
       RH100 = 100.
       DB = T(I)-273.15
       CALL WETAIR(DB,WB,RH100,DPP,BP,ECUR,ESAT,VD,RWW,TVIR,
     & TVINC,DENAIR,CPP,WTRPOT)    
       VP = VD ! VP is vapour density
       KV=0.66*DV(I)*VP*(WS(I)-(WN(I)+WN(I+1))/2)/(Z(I+1)-Z(I))
       JV(I)=KV*(H(I+1)-H(I))
       DJ(I)=MW*H(I)*KV/(R*T(I-1))
       CP(I)=-1*V(I)*WN(I)/(BB(I)*P(I)*DT)
c      # Jacobian components
       A(I)=-1*K(I-1)/(Z(I)-Z(I-1))+GR*N(I)*K(I-1)/P(I-1)
       C(I)=-1*K(I+1)/(Z(I+1)-Z(I))
       B(I)=K(I)/(Z(I)-Z(I-1))+K(I)/(Z(I+1)-Z(I))+CP(I)-GR*N(I)*K(I)/
     & P(I)+DJ(I-1)+DJ(I)
c      # mass balance
       F(I)=((P(I)*K(I)-P(I-1)*K(I-1))/(Z(I)-Z(I-1))-(P(I+1)*K(I+1)-P(I)
     & *K(I))/(Z(I+1)-Z(I)))/N1(I)+V(I)*(WN(I)-W(I))/DT-GR*(K(I-
     &1)-K(I))+JV(I-1)-JV(I)+E(I)
       SE=SE+abs(F(I))
4     continue

c     # Thomas algorithm (Gauss elimination)
      do 5 I=2,(M-1)
       C(I)=C(I)/B(I)
       if(C(I).lt.1e-8)then
        C(I)=0.D0
       endif
       F(I)=F(I)/B(I)
       B(I+1)=B(I+1)-A(I+1)*C(I)
       F(I+1)=F(I+1)-A(I+1)*F(I)
5     continue
      DP(M)=F(M)/B(M)
      P(M)=P(M)-DP(M)
      if(P(M).gt.PE(M))then
       P(M)=PE(M)
      endif

      do 6 I=(M-1),2,-1
       DP(I)=F(I)-C(I)*DP(I+1)
       P(I)=P(I)-DP(I) ! matric potential J/kg
c     # check that water potential doesn't become too large
       if(P(I).gt.PE(I))then
        P(I)=(P(I)+DP(I)+PE(I))/2
       endif
6     continue

c     # new water balance at end of the time step
      do 7 I=2,M
       WN(I)=max(WS(I)*(PE(I)/P(I))**B1(I),1.D-7)! cap minimum water content
       P(I)=PE(I)*(WS(I)/WN(I))**BB(I)! recompute water potential after capping water content to prevent instabilities
       H(I)=EXP(MW*P(I)/(R*T(I-1)))
7     continue
      H(M+1)=H(M)

c     loop until convergence
      if((SE.gt.IM).and.(count.lt.maxcount))then
       goto 11
      endif
c     end of convergence loop ##########################################

c     flux into soil, mm/m2/s (kg/m2/s)
      SW=((P(2)*K(2)-P(3)*K(3))/(N1(2)*(Z(3)-Z(2)))+GR*K(2)+TR)*DT

      W=WN

      do 9 I=2,M+1
       moistt(I-1)=WN(I)
9     continue
      
      FL=(EP*(H(2)-HA)/(1-HA))*DT
      humid(1:18)=h(2:19)
      potent(1:18)=P(2:19)
      
c     output transpiration rate, leaf and root water potential
      do 10 I=2,M
       PR(I) = -1 * (TR * RS(I) - P(I)) ! root water potential, J/kg
10    continue
      rootpot(1:18) = PR(2:19)
      leafpot = PL
      trans = TR
      
      RETURN
      END

