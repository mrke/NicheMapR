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
      DOUBLE PRECISION E_S,P_XM,F,G,H,H_A,HS,K,K_J,K_M,KAP
      DOUBLE PRECISION KAP_R,KAP_X,L,L_M,L_T,LAMBDA,M_V,MU_E,MU_V,P_A
      DOUBLE PRECISION P_AM,P_B,P_C,P_J,P_R,P_X,PREGNANT,Q,R,RDOT
      DOUBLE PRECISION RPAR,S,S_G,V,V_M,VDOT,W_V,WAITING,X,Y


      DIMENSION Y(N),DDEB(N),RPAR(34),IPAR(31)

      K_J=RPAR(1)
      P_AM=RPAR(2)
      K_M=RPAR(3)
      P_XM=RPAR(4)
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
      METAB_MODE=INT(RPAR(31))
      E_HJ=RPAR(32)
      
      ! UNPACK VARIABLES
      A = Y(1)! TIME
      V = Y(2)! cm3, STRUCTURAL VOLUME
      E = Y(3)! J/cm3, RESERVE DENSITY
      H = Y(4)! J, MATURITY
      E_S = MAX(0.,Y(5))! J, GUT ENERGY
      S = Y(6)! J, STARVATION ENERGY
      Q = Y(7)! -, DAMAGE ENERGY
      HS = Y(8)! -, AGEING ENERGY
      R = Y(9)! J, REPRODUCTION BUFFER ENERGY
      B = Y(10)! J, EGG BATCH ENERGY

      L = V ** (1./3.) ! CM, STRUCTURAL LENGTH
      V_M = L_M ** 3.
      E_SC = E / E_M                ! -, SCALED RESERVE DENSITY
      RDOT = VDOT * (E_SC / L - (1 + L_T / L) / L_M) / (E_SC + G)
      P_C = E * (VDOT / L - RDOT) * V ! J / T, MOBILISATION RATE, EQUATION 2.12 DEB3
      IF((METAB_MODE.EQ.1).AND.(H.GE.E_HJ))THEN
       RDOT = MIN(0.0D0, RDOT)
       P_C = E * V * VDOT / L ! J / T, MOBILISATION RATE, EQUATION 2.12 DEB3
      ENDIF
      DV = V * RDOT
      
      IF(H.LT.E_HB)THEN ! EMBRYO
       ! STRUCTURE
       IF(WAITING.GT.0.)THEN
        DV = 0.
       ENDIF
       dE = (- 1 *  E * VDOT) / L
       dH = (1 - KAP) * P_C - K_J * H ! J/d, change in maturity        
       P_J = K_J * H 
       DH = (1. - KAP) * P_C - P_J
       P_R = (1.-KAP) * P_C-P_J
       P_B = 0.
       ! NO AGEING OR STOMACH IN EMBRYO
       DS = 0.
       DES = 0.
       DQ = 0.
       DHS = 0.
       DR = P_R
       DB = 0.
      ELSE ! POST-EMBRYO
      
       ! structure and starvation
       IF(V * RDOT < 0.)THEN
        DS = -1. * V * RDOT * MU_V * D_V / W_V ! J / T, STARVATION ENERGY TO BE SUBTRACTED FROM REPRODUCTION BUFFER IF NECESSARY
        DV = 0.
        IF((R + B) .LT. DS)THEN !# reproduction and batch buffer has run out so draw from structure
         DV = V * RDOT
         DS = 0.
        ENDIF
       ELSE
        DS = 0.
       ENDIF
       
       ! assimilation
       P_A = P_AM * F * L ** 2.
       
       ! RESERVE
       IF(E_S .GT. P_A)THEN
        DE = P_A / L**3. - (E * VDOT) / L
       ELSE
        DE = MAX(0., E_S / L**3.) - (E * VDOT) / L
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
        P_X = F * P_XM * ((X / K) / (1. + X / K)) * V ** (2. / 3.)! J/TIME, FOOD ENERGY INTAKE RATE
       ELSE
        P_X = 0.
       ENDIF
       DES = P_X - (P_AM / KAP_X) * V**(2. / 3.)! J/TIME, CHANGE IN STOMACH ENERGY        

       IF((METAB_MODE.EQ.1).AND.(H.GE.E_HJ))THEN
        RDOT=0. ! no growth in abp after puberty - not setting this to zero messes up ageing calculation
       ENDIF
       
       ! AGEING
       DQ = (Q * (V / V_M)*S_G + H_A)*E_SC* ((VDOT / L) - RDOT)-RDOT*Q
       DHS = Q - RDOT * HS

       IF((METAB_MODE.EQ.1).AND.(H.GE.E_HJ))THEN
        P_C = P_A - DE * V
       ENDIF
       
       ! REPRODUCTION
       IF((METAB_MODE.EQ.1).AND.(H.GE.E_HJ))THEN
        P_R = (1.-KAP)*P_A-P_J
       ELSE
        P_R = (1.-KAP)*P_C-P_J
       ENDIF
      
       IF((R.LE.0.).AND.(B.LE.0).AND.(S.GT.0.).AND.(P_R.LT.S))THEN
        DV = -1. * ABS(P_R) * W_V / (MU_V * D_V)  ! SUBTRACT FROM STRUCTURE SINCE NOT ENOUGH FLOW TO REPRODUCTION TO PAY FOR SOMATIC MAINTENANCE
        P_R = 0.
       ENDIF
       IF((H .LT. E_HP) .OR. (PREGNANT.EQ.1.))THEN
        P_B = 0.
       ELSE
        IF(BATCH .EQ. 1.)THEN
          BATCHPREP = (KAP_R / LAMBDA) * ((1. - KAP) * (E_M * (VDOT * 
     &    V ** (2. / 3.) + K_M * V) / (1. + (1. / G))) - P_J)
         IF(BREEDING .LT. 1)THEN
          P_B = 0.
         ELSE
          !IF THE REPRO BUFFER IS LOWER THAN WHAT P_B WOULD BE (SEE BELOW), P_B IS P_R
          IF(R.LT.BATCHPREP)THEN
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
       P_R = P_R - P_B ! TAKE FINALISED VALUE OF P_B FROM P_R

       ! draw from reproduction and then batch buffers under starvation
       if((dS.GT.0.).AND.(R.GT.dS))THEN
         p_R = p_R - dS
         dS = 0.
       ENDIF
       if((dS.GT.0.).AND.(B.GT.dS))THEN
         p_B = p_B - dS
         dS = 0.
       ENDIF
        DR = P_R
        DB = P_B * KAP_R
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