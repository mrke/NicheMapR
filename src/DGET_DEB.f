      SUBROUTINE DGET_DEB(N,A,Y,DDEB,RPAR,IPAR)

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

C     EQUATIONS TO COMPUTE RATES OF CHANGE IN RESERVE, STRUCTURAL VOLUME, MATURITY, 
C     STOMACH ENERGY, AGEING, REPRODUCTION AND BATCH BUFFER FOR ABP DEB MODEL FOR
C     AN INSECT (HEMIMETABOLOUS)

      IMPLICIT NONE

      INTEGER IPAR,N,METAB_MODE

      DOUBLE PRECISION A,ACTHR,B,BATCH,BATCHPREP,BREEDING,D_V,DB,DDEB,DE
      DOUBLE PRECISION DES,DH,DHS,DQ,DR,DS,DV,E,E_HB,E_HJ,E_HP,E_M,E_SC
      DOUBLE PRECISION ES,ESM,F_M,F2,F,G,H,H_A,HS,JX,K,K_J,K_M,KAP
      DOUBLE PRECISION KAP_R,KAP_X,L,L_M,L_T,LAMBDA,M_V,MU_E,MU_V,P_A
      DOUBLE PRECISION P_AM,P_B,P_C,P_J,P_M,P_R,PREGNANT,Q,R,RDOT
      DOUBLE PRECISION RDOT1,RPAR,S,S_G,V,V_M,VDOT,W_V,WAITING,X,Y


      DIMENSION Y(N),DDEB(N),RPAR(34),IPAR(31)

      K_J=RPAR(1)
      P_AM=RPAR(2)
      K_M=RPAR(3)
      F_M=RPAR(4)
      VDOT=RPAR(5)
      E_M=RPAR(6)
      L_M=RPAR(7)
      L_T=RPAR(8)
      KAP=RPAR(9)
      G=RPAR(10)
      M_V=RPAR(11)
      MU_E=RPAR(12)
      MU_V=RPAR(13)
      D_V=RPAR(14)
      W_V=RPAR(15)
      ACTHR=RPAR(16)
      X=RPAR(17)
      K=RPAR(18)
      E_HP=RPAR(19)
      E_HB=RPAR(20)
      S_G=RPAR(21)
      H_A=RPAR(22)
      PREGNANT=RPAR(23)
      BATCH=RPAR(24)
      KAP_R=RPAR(25)
      LAMBDA=RPAR(26)
      BREEDING=RPAR(27)
      KAP_X=RPAR(28)
      WAITING=RPAR(29)
      F=RPAR(30)
      ESM=RPAR(31)
      METAB_MODE=INT(RPAR(32))
      E_HJ=RPAR(33)
      P_M=RPAR(34)
      
      ! UNPACK VARIABLES
      A = Y(1)! TIME
      V = Y(2)! cm3, STRUCTURAL VOLUME
      E = Y(3)! J/cm3, RESERVE DENSITY
      H = Y(4)! J, MATURITY
      ES = Y(5)! J, GUT ENERGY
      S = Y(6)! J, STARVATION ENERGY
      Q = Y(7)! -, DAMAGE ENERGY
      HS = Y(8)! -, AGEING ENERGY
      R = Y(9)! J, REPRODUCTION BUFFER ENERGY
      B = Y(10)! J, EGG BATCH ENERGY

      L = V ** (1./3.) ! CM, STRUCTURAL LENGTH
      V_M = L_M ** 3.
      E_SC = E / E_M                ! -, SCALED RESERVE DENSITY
      RDOT = VDOT * (E_SC / L - (1 + L_T / L) / L_M) / (E_SC + G)
      IF((METAB_MODE.EQ.1).AND.(H.GT.E_HJ))THEN
       RDOT = MIN(0., RDOT)
      ENDIF
      P_C = E * (VDOT / L - RDOT) * V ! J / T, MOBILISATION RATE, EQUATION 2.12 DEB3
      DV = V * RDOT
     
      IF(ES > 0.)THEN
       F2=F
      ELSE
       F2=0
      ENDIF
      
      IF(H.LT.E_HB)THEN ! EMBRYO
       ! STRUCTURE
       IF(WAITING.GT.0.)THEN
        DV = 0.
       ENDIF
       dE = (- 1 *  E * vDOT) / L
       dH = (1 - kap) * p_C - k_J * H ! J/d, change in maturity        
       P_J = K_J * H 
       DH = (1. - KAP) * P_C - P_J
       ! NO AGEING OR STOMACH IN EMBRYO
       DS = 0.
       DES = 0.
       DQ = 0.
       DHS = 0.
       DR = 0.
       DB = 0.
      ELSE ! POST-EMBRYO
      
       ! structure and starvation
       IF(V * RDOT1 < 0.)THEN
        DS = -1. * V * RDOT1 * MU_V * D_V / W_V ! J / T, STARVATION ENERGY TO BE SUBTRACTED FROM REPRODUCTION BUFFER IF NECESSARY
        IF(V .LT. DS)THEN
         DV = V * RDOT1
         DS = 0.
        ENDIF
       ELSE
        DS = 0.
       ENDIF
       
       ! assimilation
       P_A = P_Am * f * L ** 2
       IF(METAB_MODE.EQ.1)THEN
        IF((P_A.GT.P_C).AND.(E.EQ.E_M))THEN
         P_A = P_C
        ENDIF
       ENDIF
       
       ! RESERVE
       IF(ES > 0.)THEN
        DE = P_A / L**3 - (E * VDOT) / L
       ELSE
        DE = (- E * VDOT) / L
       ENDIF
       
       ! MATURATION
       P_J = K_J * H  ! adding starvation costs to P_J so maturation time (immature) or reproduction (mature) can be sacrificed to pay somatic maintenance
       IF(H .LT. E_HP)THEN
        DH = (1. - KAP) * P_C - P_J
       ELSE
        DH = 0.
       ENDIF 
       
       ! FEEDING
       IF(ACTHR .GT. 1.)THEN
        ! REGULATES X DYNAMICS
        JX = F_M * ((X / K) / (1. + X / K)) * V ** (2. / 3.)
        DES = JX * F - 1.* (P_AM / KAP_X) * V ** (2. / 3.)         
       ELSE
        DES = -1.* (P_AM / KAP_X) * V ** (2. / 3.)        
       ENDIF
       IF((METAB_MODE.EQ.1).AND.(H.GT.E_HJ))THEN
        RDOT=0. ! no growth in abp after puberty - not setting this to zero messes up ageing calculation
       ENDIF
       
       ! AGEING
       DQ = (Q * (V / V_M)*S_G + H_A)*E_SC* ((VDOT / L) - RDOT)-RDOT*Q
       DHS = Q - RDOT * HS

       ! REPRODUCTION
       p_R = (1 - KAP) * P_C - P_J
        
       IF((R.LE.0.).AND.(B.LE.0).AND.(S.GT.0.).AND.(P_R.LT.S))THEN
        DV = -1. * ABS(P_R) * W_V / (MU_V * D_V)  ! SUBTRACT FROM STRUCTURE SINCE NOT ENOUGH FLOW TO REPRODUCTION TO PAY FOR SOMATIC MAINTENANCE
        P_R = 0.
       ENDIF
       IF((H .LE. E_HP) .OR. (PREGNANT.GT.0.))THEN
        P_B = 0.
       ELSE
        IF(BATCH .GT. 0.)THEN
         IF(METAB_MODE.eq.0)THEN
          BATCHPREP = (KAP_R / LAMBDA) * ((1. - KAP) * (E_M * (VDOT * 
     &    V ** (2. / 3.) + K_M * V) / (1. + (1. / G))) - P_J)
         ELSE
          IF(METAB_MODE.eq.1)THEN
           BATCHPREP = (KAP_R / LAMBDA) * (E_M * V * VDOT / L - 
     &      KAP * P_C - P_J)
          ELSE
           BATCHPREP = (KAP_R / LAMBDA) * (E_M * (VDOT / L - RDOT) * V - 
     &      KAP * P_C - P_J)
          ENDIF
         ENDIF
         IF(BREEDING .LT. 1)THEN
          P_B = 0.
         ELSE
          !IF THE REPRO BUFFER IS LOWER THAN WHAT P_B WOULD BE (SEE BELOW), P_B IS P_R
          IF(R < BATCHPREP)THEN
           P_B = P_R
          ELSE
           !OTHERWISE IT IS A FASTER RATE, AS SPECIFIED IN PECQUERIE ET. AL JSR 2009 ANCHOVY PAPER,
           !WITH LAMBDA (THE FRACTION OF THE YEAR THE ANIMALS BREED IF FOOD/TEMPERATURE NOT LIMITING) = 0.583 OR 7 MONTHS OF THE YEAR
           P_B = MAX(BATCHPREP, P_R)
          ENDIF
         ENDIF
        ELSE
         P_B = P_R
        ENDIF!END CHECK FOR WHETHER BATCH MODE IS OPERATING
       ENDIF!END CHECK FOR IMMATURE OR MATURE

        ! draw from reproduction and then batch buffers under starvation
        if((dS.GT.0).AND.(p_R.GT.dS))THEN
          p_R = p_R - dS
          dS = 0.
        ENDIF
        if((dS.GT.0.).AND.(p_B.GT.dS))THEN
          p_B = p_B - dS
          dS = 0.
        ENDIF
       !ACCUMULATE ENERGY/MATTER IN REPRODUCTION AND BATCH BUFFERS
       IF(H > E_HP)THEN
        DR = P_R * KAP_R - P_B
        DB = P_B
       ELSE
        DR = 0
        DB = 0
       ENDIF
      ENDIF
      
      DDEB(1)=1.0D+00
      DDEB(2)=DV
      DDEB(3)=DE
      DDEB(4)=DH
      DDEB(5)=DES
      DDEB(6)=DS
      DDEB(7)=DQ
      DDEB(8)=DHS
      DDEB(9)=DR
      DDEB(10)=DB

      RETURN
      END