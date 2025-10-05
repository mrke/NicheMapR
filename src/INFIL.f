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
C     a bare soil surface, based on program 9.1 of Campbell 1985. Equation and
C     page numbers refer to Campbell, G. S. (1985). Soil Physics with Basic: 
C     Transport Models for Soil-Plant Systems. Elsevier.
 

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
     &KS(19),BB(19),PP(19),B1(19),N(19),N1(19),WS(19),rootpot(18)
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

      PI=3.1415926535897932
      MW=0.01801528 ! molar mass of water, kg/mol
      R=8.3143 ! gas constant, J/mol/K
      GR=9.81 !gravitational constant m/s/s
     
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
      RR=A
      E=A
      RS=A
      PR=A
      BZ=A
      E=A
      M=18 ! 10 user-specified depths, adding an extra depth between each of these, but not including the boundary condition depth at node 19

      PE=ABS(PE)*(-1.) !air entry potential J/kg
      PP=PE !initial water potential J/kg
      PP=ABS(PP)*(-1.)
      WS=1-BD/DD !saturation water content m3/m3
      Z(M+1)=depth(10)/100. ! depth to lower boundary, m
      WD=1000. ! density of water kg/m3
      DV=0.000024 ! binary diffusion coefficient for water vapour, m^2/s
  
      B1=1./BB
      N=2.+3./BB ! vector per soil layer of exponents in equation to obtain hydraulic conductivity from saturated hydraulic conductivity, saturated water content and soil water content
      N1=1.-N

c     prep for wetair call      
      WB = 0.
      DPP = 999.
      PSTD=101325.
      BP=PSTD*((1.-(0.0065*ALTT/288.))**(5.255785959124322))      

      j=2
      do 121 I=3,18 ! add in extra nodes between the 10 soil depths specified by the user for soil temperature calcs
       if(MOD(I, 2).ne.0)then
        Z(I)=depth(j)/100.
        j=j+1
       else
        Z(I)=Z(I-1)+(depth(J)/100.-Z(I-1))/2
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

      T=T+273.15 ! temperature from deg C to Kelvin

      Z(2)=0
      Z(1)=0

c     # setting initial water content, m3/m3
      do 2 I=2,M
       WN(I)=moistt(i-1)
       P(I)=PE(I)*(WS(I)/WN(I))**BB(I) ! matric water potential, EQ5.9 (note thetas=W are inverted so not raised to -BB)
       H(I)=exp(MW*P(I)/(R*T(I-1))) ! fractional humidity, EQ5.14
       K(I)=KS(I)*(PE(I)/P(I))**N(I) ! hydraulic conductivity, EQ6.14
       W(I)=WN(I) ! water content
2     continue

      do 22 I=2,M
       V(I)=WD*(Z(I+1)-Z(I-1))/2 ! bulk density x volume per unit area
22    continue
      
c     # lower boundary condition set to saturated (stays constant)
      P(M+1)=PE(M)*(WS(M+1)/WS(M+1))**BB(M) ! water potential
      H(M+1)=1. ! fractional humidity
      W(M+1)=WS(M+1) ! water content
      WN(M+1)=WS(M+1) ! water content
      Z(1)=-1D10 ! depth at node 1, m
      Z(M+1)=1D20 ! depth at deepest node, m
      K(M+1)=KS(M)*(PE(M)/P(M+1))**N(M+1) ! lower boundary conductivity

c     initialize root water uptake variables
      do 98 I=2,M
       if(L(I).gt.0.)then
        RR(I)=RW/(L(I)*(Z(I+1)-Z(I-1))/2.) ! root resistance
        BZ(I)=N1(I)*LOG(PI*R1*R1*L(I))/
     &   (4.*PI*L(I)*(Z(I+1)-Z(I-1))/2.)
       else
        RR(I)=1D+20 ! root resistance
        BZ(I)=0.D0
       endif
98    continue

      P(1)=P(2)
      K(1)=0.D0

c     evapotranspiration
      EP=exp(-0.82*LAI)*ET ! partition potential evaporation from potential evapotranspiration, EQ12.30
      TP=ET-EP ! now get potential transpiration

c     plant water uptake
      PB=0.D0 ! weighted mean soil water potential, psi_bar, J/kg
      RB=0.D0 ! weighted mean root-soil resistance, R_bar, m4 /(s kg)
      PL=0.D0 ! leaf water potential, J/kg
      do 99 i=2,M
       RS(I)=BZ(I)/K(I) ! soil resistance, simplification of EQ11.14, assuming conductivity constant in the rhizosphere
       PB=PB+P(I)/(RS(I)+RR(I)) ! summing over layers
       RB=RB+1.D0/(RS(I)+RR(I)) ! summing over layers
99    continue
      PB=PB/RB ! final step in evaluating psi_bar, first term on right in EQ11.18
      RB=1.D0/RB ! denominator of first and second terms on right in EQ11.18
      maxcount=500
      count=0

c     begin Newton-Raphson procedure to estimate leaf water potential
2080  continue
      IF(PL.gt.PB)then
       PL=PB-TP*(RB+RL) ! variation on EQ11.18
      ENDIF
      XP=(PL/PC)**SP ! part of EQ12.28 determining stomatal closure
      SL=TP*(RB+RL)*SP*XP/(PL*(1.D0+XP)*(1.D0+XP))-1.D0 ! derivative of stomatal function
      FF=PB-PL-TP*(RB+RL)/(1.D0+XP) ! transpiration mass balance (variation on EQ11.18)
      PL=PL-FF/SL
      count=count+1
      if((ABS(FF).gt.10).and.(count.lt.maxcount))then
       goto 2080
      endif
      TR=TP/(1+XP)
      do 100 I=2,M
       E(I)=(P(I)-PL-RL*TR)/(RR(I)+RS(I)) ! root water uptake, EQ11.15
100   continue
      count=0
      
c     start of convergence loop ########################################
      maxcount=500
      RH100 = 100.
11    SE=0
      count=count+1
      do 3 I=2,M
       K(I)=KS(I)*(PE(I)/P(I))**N(I) ! hydraulic conductivities for each node, EQ6.14
3     continue
      JV(1)=EP*(H(2)-HA)/(1.D0-HA) ! vapour flux at soil surface, EQ9.14
      DJ(1)=EP*MW*H(2)/(R*T(1)*(1.D0-HA)) ! derivative of vapour flux at soil surface, combination of EQ9.14 and EQ5.14
      do 4 I=2,M
       DB = T(I)-273.15 ! back to deg C from Kelvin for call to wetair
       CALL WETAIR(DB,WB,RH100,DPP,BP,ECUR,ESAT,VD,RWW,TVIR,
     & TVINC,DENAIR,CPP,WTRPOT) ! getting saturated vapour density for current temperature    
       VP = VD ! VP is vapour density = c'_v in EQ9.7
       KV=0.66*DV*VP*(WS(I)-(WN(I)+WN(I+1))/2.D0)/(Z(I+1)-Z(I)) ! vapour conductivity, EQ9.7, assuming epsilon(psi_g) = b*psi_g^m (eq. 3.10) where b = 0.66 and m = 1 (p.99)
       JV(I)=KV*(H(I+1)-H(I)) ! fluxes of vapour within soil, EQ9.14
       DJ(I)=MW*H(I)*KV/(R*T(I-1)) ! derivatives of vapour fluxes within soil, combination of EQ9.14 and EQ5.14
       CP(I)=-1.D0*V(I)*WN(I)/(BB(I)*P(I)*DT) ! hydraulic capacity = capacitance, d_theta/d_psi
c      # Jacobian components
       A(I)=-1.D0*K(I-1)/(Z(I)-Z(I-1))+GR*N(I)*K(I-1)/P(I-1)! sub-diagonal element in tridagonal matrix
       C(I)=-1.D0*K(I+1)/(Z(I+1)-Z(I))! super-diagonal element in tridagonal matrix
       B(I)=K(I)/(Z(I)-Z(I-1))+K(I)/(Z(I+1)-Z(I))+CP(I)-GR*N(I)*K(I)/
     & P(I)+DJ(I-1)+DJ(I) ! diagonal element in tridagonal matrix
c      # mass balance including vapour fluxes and root water uptake
       F(I)=((P(I)*K(I)-P(I-1)*K(I-1))/(Z(I)-Z(I-1))-(P(I+1)*K(I+1)-P(I)
     & *K(I))/(Z(I+1)-Z(I)))/N1(I)+V(I)*(WN(I)-W(I))/DT-GR*(K(I-
     &1)-K(I))+JV(I-1)-JV(I)+E(I) ! version of equation 8.28 that additionally conatins vapour fluxes and root water uptake
       SE=SE+abs(F(I)) ! total mass balance error
4     continue

c     # Thomas algorithm (Gauss elimination)
      do 5 I=2,(M-1)
       C(I)=C(I)/B(I)
       if(C(I).lt.1D-8)then
        C(I)=0.D0
       endif
       F(I)=F(I)/B(I)
       B(I+1)=B(I+1)-A(I+1)*C(I)
       F(I+1)=F(I+1)-A(I+1)*F(I)
5     continue
      DP(M)=F(M)/B(M) ! change in potential in an iteration step, J/kg
      P(M)=P(M)-DP(M) ! matric potential J/kg
      if(P(M).gt.PE(M))then
       P(M)=PE(M)
      endif

      do 6 I=(M-1),2,-1 
       DP(I)=F(I)-C(I)*DP(I+1) ! change in potential in an iteration step, J/kg
       P(I)=P(I)-DP(I) ! matric potential J/kg
c     # check that water potential doesn't become too large
       if(P(I).gt.PE(I))then
        P(I)=(P(I)+DP(I)+PE(I))/2.D0
       endif
6     continue

c     # new water balance at end of the time step
      do 7 I=2,M
       WN(I)=max(WS(I)*(PE(I)/P(I))**B1(I),1.D-7)! cap minimum water content
       P(I)=PE(I)*(WS(I)/WN(I))**BB(I)! recompute water potential after capping water content to prevent instabilities
       H(I)=EXP(MW*P(I)/(R*T(I-1))) ! recompute fractional humidity
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
      
      FL=(EP*(H(2)-HA)/(1.D0-HA))*DT
      humid(1:18)=h(2:19)
      potent(1:18)=P(2:19)
      
c     output transpiration rate, leaf and root water potential
      do 10 I=2,M
       PR(I)=-1.D0*(TR*RS(I)-P(I)) ! root water potential, J/kg
10    continue
      rootpot(1:18) = PR(2:19)
      leafpot = PL
      trans = TR
      
      RETURN
      END

