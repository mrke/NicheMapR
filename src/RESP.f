      SUBROUTINE RESP

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

C     COMPUTES RESPIRATORY HEAT AND WATER LOSS VIA MASS FLOW THROUGH THE LUNGS,
C     GIVEN GAS CONCENTRATIONS, PRESSURE, RESPIRATION RATE AND HUMIDITY

      IMPLICIT NONE

      DOUBLE PRECISION A1,A2,A3,A4,A4B,A5,A6,AIRML1,AIRML2,AIRVOL
      DOUBLE PRECISION AL,ALT,AMASS,ANDENS,AREF,ASILP,ATOT,BP,BREF
      DOUBLE PRECISION CO2GAS,CO2MOL1,CP,CREF,CUSTOMGEOM,DB,DEBQMET
      DOUBLE PRECISION DEBQMET_INIT,DELTAR,DENAIR,DEPSUB,DP,DRYFOOD,E
      DOUBLE PRECISION EGGSHP,EMISAN,EMISSB,EMISSK,ESAT,EVPMOL,EXTREF
      DOUBLE PRECISION F12,F13,F14,F15,F16,F21,F23,F24,F25,F26,F31,F32
      DOUBLE PRECISION F41,F42,F51,F52,F61,FAECES,FATOSB,FATOSK,FLSHCOND
      DOUBLE PRECISION FLUID,FOODWATERCUR,G,GEVAP,GH2OMET,GH2OMET_INIT
      DOUBLE PRECISION GMASS,H2O_BALPAST,HTOVPR,KGEVAP,MLO2,MLO2_INIT
      DOUBLE PRECISION MR_1,MR_2,MR_3,N2GAS,N2MOL1,N2MOL2,NWASTE,O2GAS
      DOUBLE PRECISION O2MOL1,O2MOL2,O2MOLC,O2STP,PCTCO2,PCTN2,PCTO2
      DOUBLE PRECISION PFEWAT,PHI,PHIMAX,PHIMIN,PO2,PTCOND,PTCOND_ORIG
      DOUBLE PRECISION QCOND,QCONV,QIRIN,QIROUT,QMETAB,QRESP,QSEVAP
      DOUBLE PRECISION QSOLAR,R,REFPO2,RELHUM,RELXIT,RGC,RHO1_3,RHSAT
      DOUBLE PRECISION RPCTCO2,RPCTN2,RPCTO2,RQ,RW,SHP,SIDEX,SIG,SUBTK
      DOUBLE PRECISION T_A,T_AH,T_AL,T_H,T_L,T_REF,TA,TAIR,TBASK,TC
      DOUBLE PRECISION TCORR,TDIGPR,TEMERGE,TLUNG,TMAXPR,TMINPR,TPREF
      DOUBLE PRECISION TQSOL,TR,TRANS1,TSKIN,TSUBST,TVINC,TVIR,TWING,VD
      DOUBLE PRECISION VEL,VO2CON,WB,WEVAP,WMOL1,WMOL2,WQSOL,WTRPOT
      DOUBLE PRECISION XCALC,XTRY,TOTGAS,VAIR,VCO2,CO2MOL2
      DOUBLE PRECISION WEYES,WRESP,WCUT,AEFF,CUTFA,HD,PEYES,SKINW,
     & SKINT,HC,CONVAR,PMOUTH,PANT,PANTMAX,EGGPTCOND,POT

      INTEGER DEB1,GEOMETRY,IHOUR,NODNUM,WINGMOD,WINGCALC

      CHARACTER*1 TRANST
      
      DIMENSION CUSTOMGEOM(8),DEBQMET(24),DRYFOOD(24),FAECES(24)
      DIMENSION GH2OMET(24),MLO2(24),NWASTE(24),SHP(3),EGGSHP(3)

      COMMON/AIRGAS/O2GAS,CO2GAS,N2GAS
      COMMON/ARRHEN/T_A,T_AL,T_AH,T_L,T_H,T_REF
      COMMON/BEHAV2/GEOMETRY,NODNUM,CUSTOMGEOM,SHP,EGGSHP
      COMMON/DEBRESP/MLO2,GH2OMET,DEBQMET,MLO2_INIT,GH2OMET_INIT,
     & DEBQMET_INIT,DRYFOOD,FAECES,NWASTE
      COMMON/EVAP1/WEYES,WRESP,WCUT,AEFF,CUTFA,HD,PEYES,SKINW,
     & SKINT,HC,CONVAR,PMOUTH,PANT,PANTMAX 
      COMMON/FUN1/QSOLAR,QIRIN,QMETAB,QRESP,QSEVAP,QIROUT,QCONV,QCOND
      COMMON/FUN2/AMASS,RELHUM,ATOT,FATOSK,FATOSB,EMISAN,SIG,FLSHCOND
      COMMON/FUN3/AL,TA,VEL,PTCOND,SUBTK,DEPSUB,TSUBST,PTCOND_ORIG,
     & EGGPTCOND,POT
      COMMON/FUN4/TSKIN,R,WEVAP,TR,ALT,BP,H2O_BALPAST
      COMMON/GITRAC/PFEWAT,FOODWATERCUR
      COMMON/GUESS/XTRY
      COMMON/REVAP1/TLUNG,DELTAR,EXTREF,RQ,MR_1,MR_2,MR_3,DEB1
      COMMON/REVAP2/GEVAP,AIRVOL,CO2MOL2
      COMMON/TPREFR/TMAXPR,TMINPR,TDIGPR,TPREF,TBASK,TEMERGE
      COMMON/TREG/TC
      COMMON/USROPT/TRANST
      COMMON/WDSUB1/ANDENS,ASILP,EMISSB,EMISSK,FLUID,G,IHOUR
      COMMON/WINGFUN/RHO1_3,TRANS1,AREF,BREF,CREF,PHI,F21,F31,F41,F51
     &,SIDEX,WQSOL,PHIMIN,PHIMAX,TWING,F12,F32,F42,F52
     &,F61,TQSOL,A1,A2,A3,A4,A4B,A5,A6,F13,F14,F15,F16,F23,F24,F25,F26
     &,WINGCALC,WINGMOD

C     NOTE THAT THERE IS NO RECOVERY OF HEAT OR MOISTURE ASSUMED IN THE NOSE

      TAIR = TA
C     KLUGE FOR THE MOMENT TO DECIDE BODY-AIR TEMPERATURE GRADIENT AND CHECK FOR
C     STABILITY BEFORE DOING TLUNG - TAIR,LOCAL
      DELTAR = 1.0

C     DEFINING VARIABLES
C     BP = BAROMETRIC PRESSURE (PA)
C     EXTREF = EXTRACTION EFFICIENCY (PER CENT)
C     GEVAP = GRAMS OF WATER EVAPORATED FROM RESPIRATORY TRACT/S
C     QRESP = HEAT LOSS DUE TO RESPIRATORY EVAPORATION (W)
C     RGC = UNIVERSAL GAS CONSTANT (PA-M3/MOL-K) = (J/MOL-K)
C     RELHUM = RELATIVE HUMIDITY (PER CENT)
C     RQ = RESPIRATORY QUOTIENT (MOL CO2 / MOL O2)
C     TC = ANIMAL CORE TEMPERATURE(C)
C     TMAXPR = PREFERRED MAX. TCORE
C     TMINPR = PREFERRED MIN. TCORE

C     ASSIGNING REFERENCE VALUES TO VARIABLES
C     AIR FRACTIONS FROM SCHMIDT-NIELSEN, 2ND ED. ANIMAL PHYSIOLOGY CITING
C     OTIS, 1964
      RPCTO2 = 0.2095
      RPCTN2 = 0.7902
      RPCTCO2 = 0.0003
      PCTO2 = RPCTO2
      PCTN2 = RPCTN2
      PCTCO2 = RPCTCO2
C     ALLOWING USER TO MODIFY GAS VALUES FOR BURROW, ETC. CONDITIONS
      IF(PCTO2.NE.O2GAS/100.)THEN
       PCTO2=O2GAS/100.
      ELSE
       PCTO2=RPCTO2
      ENDIF
      IF(PCTN2.NE.N2GAS/100.)THEN
       PCTN2=N2GAS/100.
      ELSE
       PCTN2=RPCTN2
      ENDIF
      IF(PCTCO2.NE.CO2GAS/100.)THEN
       PCTCO2=CO2GAS/100.
      ELSE
       PCTCO2=RPCTCO2
      ENDIF
      
C     ERROR CHECKING for % of each air constituent summing to 1.00000
      TOTGAS = PCTO2 + PCTN2 + PCTCO2
      IF(TOTGAS .GT. 1.00000)THEN
       PCTO2 = 1.00000 - (PCTN2 + PCTCO2)
      ENDIF
      IF(TOTGAS .LT. 1.00000)THEN
       PCTO2 = 1.00000 - (PCTN2 + PCTCO2)
      ENDIF      
      
C     UNIVERSAL GAS CONSTANT (PA - LITERS)/(MOL - K)
      RGC=8314.46
C     INITIALIZING FOR SUB. WETAIR
      WB=0.0
      DP=999.

      PO2=BP*PCTO2
      REFPO2=101325.*RPCTO2
C     OXYGEN CONSUMPTION OF LIZARDS (BENNETT & DAWSON, 1976) (M3/S)
      GMASS=AMASS*1000.

      IF((TRANST.EQ.'Y').OR.(TRANST.EQ.'Y'))THEN
       XCALC=TC
       TLUNG=TC
      ELSE
       XCALC=XTRY
      ENDIF

C     CHECK FOR TOO LARGE OR SMALL CORE TEMP
      IF(TC.GT.50.) THEN
       XCALC=50.
      ELSE
       IF(TC.LT.0.0000) THEN
        XCALC=0.01
       ENDIF
      ENDIF
       
      IF(DEB1.EQ.1)THEN
       TCORR=EXP(T_A*(1./(273.15+T_REF)-1./(273.15+XCALC)))/(1.+EXP(T_AL
     & *(1./(273.15+XCALC)-1/T_L))+EXP(T_AH*(1/T_H-1/(273.15+XCALC))))
C      USING DEB CALCS FOR MET WATER AND O2 CONSUMPTION (CONVERT FROM ML/H TO L/S), SO SKIP THE REGRESSIONS!
       IF(IHOUR.EQ.1)THEN
        O2STP=MLO2_INIT/3600./1000.*TCORR
       ELSE
        O2STP=MLO2(IHOUR-1)/3600./1000.*TCORR
       ENDIF
      ELSE
       O2STP=(1./3.6E+06)*(QMETAB/.0056)
      ENDIF

C     CONVERTING STP -> VOL. OF O2 AT ANIMAL TCORE, ATM. PRESS.
      TLUNG = XCALC
C      VO2CON = (O2STP*REFPO2/273.15)*((TLUNG+273.15)/PO2) ! V2 = (V1*P1/T1) * (T2/P2), VOLUME OF O2 AT TLUNG
      VO2CON = (O2STP*PO2/273.15)*((TLUNG+273.15)/PO2) ! CHANGE FROM WPP 16/10/2021
C     N = PV/RT (IDEAL GAS LAW: NUMBER OF MOLES FROM PRESS,VOL,TEMP)
      O2MOLC = BP*VO2CON/(RGC*(XCALC+273.15)) ! MOL OXYGEN CONSUMED
C     MOLES/S O2,N2, & DRY AIR AT 1: (ENTRANCE) (AIR FLOW = F(O2 CONSUMPTION)
      O2MOL1 = O2MOLC/(EXTREF/100.) ! ACTUAL OXYGEN FLOW IN (MOLES/S), ACCOUNTING FOR EFFICIENCY OF EXTRACTION
      N2MOL1 = O2MOL1*(PCTN2/PCTO2) ! ACTUAL NITROGEN FLOW IN (MOLES/S), ACCOUNTING FOR EFFICIENCY OF EXTRACTION
      VAIR = VO2CON/PCTO2 ! CHANGE FROM WPP 16/10/2021
      VCO2 = PCTCO2*VAIR ! CHANGE FROM WPP 16/10/2021
      CO2MOL1 = BP*VCO2/(RGC*(TLUNG+273.15)) ! CHANGE FROM WPP 16/10/2021       
C     DEMAND FOR AIR = F(%O2 IN THE AIR AND ELEVATION)
C     NOTE THAT AS LONG AS ALL 3 PERCENTAGES ADD TO 100%, NO CHANGE IN AIR FLOW,
C     UNLESS YOU CORRECT FOR CHANGE IN %O2 IN THE AIR AND ELEVATION CHANGES
C     RELATIVE TO SEA LEVEL.
C      AIRATO=(PCTN2+PCTO2+PCTCO2)/PCTO2
C      AIRML1=O2MOL1*AIRATO*(RPCTO2/PCTO2)*(REFPO2/PO2)*PANT
      AIRML1 = (O2MOL1 + N2MOL1 + CO2MOL1) * PANT ! CHANGE FROM WPP 16/10/2021
C     AIR VOLUME @ STP (LITERS/S)
      AIRVOL=(AIRML1*RGC*273.15/101325.)

C     COMPUTING THE VAPOR PRESSURE AT SATURATION FOR THE SUBSEQUENT
C     CALCULATION OF ACTUAL MOLES OF WATER BASED ON ACTUAL RELATIVE
C     HUMIDITY.
      RHSAT=100.
      CALL WETAIR(TAIR,WB,RELHUM,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,
     * CP,WTRPOT)
      WMOL1=AIRML1*(ESAT*(RELHUM/100.))/(BP-ESAT*(RELHUM/100.))

C     MOLES AT 2: (EXIT)
      O2MOL2=O2MOL1-O2MOLC ! REMOVE CONSUMED OXYGEN FROM THE TOTAL
      N2MOL2=N2MOL1
C      CO2MOL2=RQ*O2MOLC
      CO2MOL2 = RQ*O2MOLC + CO2MOL1 ! CHANGE FROM WPP 16/10/2021
C     TOTAL MOLES OF AIR AT 2 (EXIT) WILL BE APPROXIMATELY THE SAME
C     AS AT 1, SINCE THE MOLES OF O2 REMOVED = APPROX. THE # MOLES OF CO2
C     ADDED.  AVOGADRO'S # SPECIFIES THE # MOLECULES/MOLE.
C      AIRML2=(O2MOL2+CO2MOL)*((PCTN2+PCTO2)/PCTO2)*(RPCTO2/PCTO2)*
C     & (REFPO2/PO2)*PANT
      AIRML2 = (O2MOL2 + N2MOL2 + CO2MOL2) * PANT ! CHANGE FROM WPP 16/10/2021

C     SETTING UP CALL TO WETAIR; TEMP. OF EXHALED AIR AT BODY TEMP.
      DB=XCALC
C     ASSUMING SATURATED AIR AT EXHALATION
      RELXIT=100.
      CALL WETAIR(DB,WB,RELXIT,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,
     * CP,WTRPOT)
      WMOL2=AIRML2*(ESAT/(BP-ESAT))
C     ENTHALPY = U2-U1, INTERNAL ENERGY ONLY, I.E. LAT. HEAT OF VAP.
C     ONLY INVOLVED, SINCE ASSUME P,V,T CONSTANT, SO NOT SIGNIFICANT
C     FLOW ENERGY, PV. (H = U + PV)

C     MOLES/S LOST BY BREATHING:
      EVPMOL=WMOL2-WMOL1
C     GRAMS/S LOST BY BREATHING = MOLES LOST * GRAM MOLECULAR WEIGHT OF WATER:
      GEVAP=EVPMOL*18.
      IF(GEVAP.GT.100)THEN
      GEVAP=EVPMOL*18.
      ENDIF

      KGEVAP = GEVAP/1000.
C     LATENT HEAT OF VAPORIZATION FROM SUB. DRYAIR
      HTOVPR = 2.5012E+06 - 2.3787E+03*TLUNG
C     HEAT LOSS BY BREATHING (J/S)=(J/KG)*(KG/S)
      QRESP = HTOVPR*KGEVAP

      RETURN
      END