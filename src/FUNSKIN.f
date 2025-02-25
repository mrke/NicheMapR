      FUNCTION FUNSKIN (X)

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

C     THIS FUNCTION IS USED FOR TRANSIENT HEAT BUDGET CALCULATIONS TO COMPUTE
C     SKIN TEMPERATURE - IT IS CALLED BY ONELUMP AND SOLVED VIA ZBRENT

      IMPLICIT NONE

      DOUBLE PRECISION A,A1,A2,A3,A4,A4B,A5,A6,ABSMAX,ABSMIN,AEFF,AIRVOL
      DOUBLE PRECISION AL,ALT,AMASS,ANDENS,AREA,AREF,ASEMAJR,ASIL,ASILP
      DOUBLE PRECISION ASQ,ATOT,B,BP,BREF,BSEMINR,BSQ,C,CO2MOL,CONVAR,CP
      DOUBLE PRECISION CREF,CSEMINR,CSQ,CUSTOMGEOM,CUTFA,DELTAR,DENAIR
      DOUBLE PRECISION DEPSUB,DP,E,EGGSHP,EMISAN,EMISSB,EMISSK,ENB,ESAT
      DOUBLE PRECISION EXTREF,F12,F13,F14,F15,F16,F21,F23,F24,F25,F26
      DOUBLE PRECISION F31,F32,F41,F42,F51,F52,F61,FATCOND,FATOSB,FATOSK
      DOUBLE PRECISION FLSHCOND,FLTYPE,FLUID,FLYMETAB,FLYSPEED,FLYTIME
      DOUBLE PRECISION FUNSKIN,G,GEVAP,GN,H2O_BALPAST,HC,HD,HR,HTOVPR
      DOUBLE PRECISION MR_1,MR_2,MR_3,PANT,PANTMAX,PEYES,PHI,PHIMAX
      DOUBLE PRECISION PHIMIN,PI,PMOUTH,PSI_BODY,PTCOND,PTCOND_ORIG
      DOUBLE PRECISION QCOND,QCONV,QGENET,QIRIN,QIROUT,QMETAB,QRESP
      DOUBLE PRECISION QSEVAP,QSOLAR,QSOLR,QSWEAT,R,R1,RELHUM,RFLESH
      DOUBLE PRECISION RHO1_3,RHUM,RINSUL,RQ,RSKIN,RW,S1,S2,SHADE,SHP
      DOUBLE PRECISION SIDEX,SIG,SKINT,SKINW,SPHEAT,SUBTK,TA,TAVE,TBASK
      DOUBLE PRECISION TC,TDIGPR,TEMERGE,TLUNG,TMAXPR,TMINPR,TOBJ,TPREF
      DOUBLE PRECISION TQSOL,TR,TRAD,TRANS1,TSKIN,TSKY,TSUBST,TVINC,TVIR
      DOUBLE PRECISION TWING,VD,VDAIR,VDSURF,VEL,VOL,WB,WCUT,WEVAP,WEYES
      DOUBLE PRECISION WQSOL,WRESP,WTRPOT,X,XTRY,EGGPTCOND,POT
      DOUBLE PRECISION CONV_ENHANCE,G_VS_AB,G_VS_AD,MR_4
      
      INTEGER CLIMBING,DEB1,FLIGHT,FLYER,FLYTEST,GEOMETRY,IHOUR,LIVE
      INTEGER MICRO,NODNUM,WINGCALC,WINGMOD,LEAF,FLYHIGH

      DIMENSION CUSTOMGEOM(8),SHP(3),EGGSHP(3)
      
      COMMON/ANPARMS/RINSUL,R1,AREA,VOL,FATCOND
      COMMON/BEHAV2/GEOMETRY,NODNUM,CUSTOMGEOM,SHP,EGGSHP
      COMMON/CLIMB/CLIMBING
      COMMON/ELLIPS/ASEMAJR,BSEMINR,CSEMINR
      COMMON/EVAP1/WEYES,WRESP,WCUT,AEFF,CUTFA,HD,PEYES,SKINW,G_VS_AB,
     & G_VS_AD,SKINT,HC,CONVAR,PMOUTH,PANT,PANTMAX,CONV_ENHANCE,LEAF
      COMMON/FLY/FLYTIME,FLYSPEED,FLYMETAB,FLIGHT,FLYER,FLYTEST,FLYHIGH
      COMMON/FUN1/QSOLAR,QIRIN,QMETAB,QRESP,QSEVAP,QIROUT,QCONV,QCOND
      COMMON/FUN2/AMASS,RELHUM,ATOT,FATOSK,FATOSB,EMISAN,SIG,FLSHCOND
      COMMON/FUN3/AL,TA,VEL,PTCOND,SUBTK,DEPSUB,TSUBST,PTCOND_ORIG,
     & EGGPTCOND,POT
      COMMON/FUN4/TSKIN,R,WEVAP,TR,ALT,BP,H2O_BALPAST
      COMMON/FUN6/SPHEAT,ABSMAX,ABSMIN,LIVE
      COMMON/GUESS/XTRY
      COMMON/REVAP1/TLUNG,DELTAR,EXTREF,RQ,MR_1,MR_2,MR_3,MR_4,DEB1
      COMMON/REVAP2/GEVAP,AIRVOL,CO2MOL
      COMMON/SOLN/ENB
      COMMON/TPREFR/TMAXPR,TMINPR,TDIGPR,TPREF,TBASK,TEMERGE
      COMMON/TREG/TC
      COMMON/WATERPOT/PSI_BODY
      COMMON/WCONV/FLTYPE
      COMMON/WDSUB1/ANDENS,ASILP,EMISSB,EMISSK,FLUID,G,IHOUR
      COMMON/WDSUB2/QSOLR,TOBJ,TSKY,MICRO
      COMMON/WINGFUN/RHO1_3,TRANS1,AREF,BREF,CREF,PHI,F21,F31,F41,F51,
     & SIDEX,WQSOL,PHIMIN,PHIMAX,TWING,F12,F32,F42,F52,F61,TQSOL,A1,A2,
     & A3,A4,A4B,A5,A6,F13,F14,F15,F16,F23,F24,F25,F26,WINGCALC,WINGMOD
      COMMON/WMET/QSWEAT
      COMMON/WSOLAR/ASIL,SHADE


      DATA PI/3.14159265/
      RFLESH=0.
      CSQ=0.
      BSQ=0.
      ASQ=0.
      
C     THE GUESSED VARIABLE, X, IS CORE TEMPERATURE (C); SEE SUB. MET
C     FOR DETAILED EXPLANATION OF CALCULATION OF SURF. TEMP., TSKIN,
C     FROM TC AND MASS
C     THIS ASSUMES UNIFORM BODY TEMPERATURE.

C     CONTROL OF BODY TEMPERATURE GUESSES FOR STABILITY PURPOSES
      IF (X .GT. 80.) THEN
       X = 80.
      ELSE
C      IF (X .LT. -20.0) THEN
C       X = TSKY + 0.1
C      ENDIF
      ENDIF

      TSKIN = X
      XTRY = X

C     GET THE METABOLIC RATE
C     CHECKING FOR INANIMATE OBJECT
      IF (LIVE .EQ. 0) THEN
C      INANIMATE
       QMETAB = 0.0
       TSKIN = X
      ELSE
C      ALIVE, BUT IS IT TOO COLD?
       IF (TC .GE. 0.0) THEN
        CALL MET
       ELSE
C       TOO COLD, SUPER LOW METABOLISM
        QMETAB = 0.0001
        TSKIN = X
       ENDIF
      ENDIF

C     GET THE RESPIRATORY WATER LOSS
C     CHECKING FOR FLUID TYPE
      IF (FLTYPE .EQ. 0.00) THEN
C      AIR
C      CALL FOR RESPIRATORY WATER & ENERGY LOSS
       IF (QMETAB .GE. 0.000) THEN
        CALL RESP
       ELSE
C       NEGATIVE METABOLIC RATE. NO PHYSIOLOGICAL MEANING - DEAD.
        QRESP = 0.00000
        QMETAB = 0.00000
       ENDIF
      ENDIF

C     NET INTERNAL HEAT GENERATION
      QGENET = QMETAB - QRESP
C     NET INTERNAL HEAT GENERATION/UNIT VOLUME. USE FOR ESTIMATING SKIN TEMP.
      GN = QGENET/VOL
      IF (LIVE .EQ. 0) THEN
       GN = 0.
      ENDIF

C     COMPUTING SURFACE TEMPERATURE AS DICTATED BY GEOMETRY

C     FIRST SET AVERAGE BODY TEMPERATURE FOR ESTIMATION OF AVEARAGE LUNG TEMPERATURE
      IF(GEOMETRY.EQ.1)THEN
C      CYLINDER: FROM P. 270 BIRD, STEWART & LIGHTFOOT. 1960. TRANSPORT PHENOMENA.
C      TAVE = (GR**2/(8K)) + TSKIN, WHERE TSKIN = TCORE - GR**2/(4K)
C      NOTE:  THESE SHOULD ALL BE SOLVED SIMULTANEOUSLY.  THIS IS AN APPROXIMATION
C      USING CYLINDER GEOMETRY. SUBCUTANEOUS FAT IS ALLOWED IN CYLINDER & SPHERE
C      CALCULATIONS.
       RFLESH = R1 - RINSUL
C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN
       TLUNG = (GN*RFLESH**2.)/(8.*FLSHCOND) + TSKIN
      ENDIF
      IF(GEOMETRY.EQ.2)THEN
C      ELLIPSOID: DERIVED 24 OCTOBER, 1993  W. PORTER
       A = ASEMAJR
       B = BSEMINR
       C = CSEMINR
       ASQ = A**2.
       BSQ = B**2.
       CSQ = C**2.
C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN
       TLUNG = (GN/(4.*FLSHCOND)) * ((ASQ*BSQ*CSQ)/
     &  (ASQ*BSQ+ASQ*CSQ+BSQ*CSQ)) + TSKIN
      ENDIF
      IF(GEOMETRY.EQ.4)THEN
C      SPHERE:
       RFLESH = R1 - RINSUL
       RSKIN = R1
C      FAT LAYER, IF ANY
       S1 = (QGENET/(4.*PI*FLSHCOND))*((RFLESH - RSKIN)/(RFLESH*RSKIN))
C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN (12 BECAUSE TLUNG IS 1/2 THE TC-TSKIN DIFFERENCE, 6*AK1)
       TLUNG = (GN*RFLESH**2.)/(12.*FLSHCOND) + TSKIN
      ENDIF
      IF((GEOMETRY.EQ.3).OR.(GEOMETRY.EQ.5))THEN
C      MODEL LIZARD AS CYLINDER
C      CYLINDER: FROM P. 270 BIRD, STEWART & LIGHTFOOT. 1960. TRANSPORT PHENOMENA.
C      TAVE = (GR**2/(8K)) + TSKIN, WHERE TSKIN = TCORE - GR**2/(4K)
C      NOTE:  THESE SHOULD ALL BE SOLVED SIMULTANEOUSLY.  THIS IS AN APPROXIMATION
C      USING CYLINDER GEOMETRY. SUBCUTANEOUS FAT IS ALLOWED IN CYLINDER & SPHERE
C      CALCULATIONS.
       RFLESH = R1 - RINSUL
C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN
       TLUNG = (GN*RFLESH**2.)/(8.*FLSHCOND) + TSKIN
      ENDIF
C     LIMITING LUNG TEMPERATURE EXTREMES
      IF (TLUNG .GT. TC) THEN
       TLUNG = TC
      ENDIF
      IF (TLUNG .LT. -3.) THEN
       TLUNG = -3.
      ENDIF

      CALL CONV
      CALL RESP
      CALL SEVAP
      CALL RADOUT
      CALL COND

      IF (FLTYPE .EQ. 1.00) THEN
C      WATER ENVIRONMENT
       QSEVAP = 0.00
       WEVAP = 0.0
       QIRIN = 0.0
       QIROUT = 0.00
       QCOND = 0.00
      ENDIF
      
      TRAD=(QIRIN/(EMISAN*1.*ATOT*SIG))**(1./4.)-273.15
      HTOVPR = 2.5012E+06 - 2.3787E+03 * TA
      TAVE=(TRAD+TSKIN)/2.
      HR=4.*EMISAN*SIG*(TAVE+273.15)**3.
      
C     INITIALIZING FOR SUB. WETAIR
      WB = 0.0
      DP = 999.
      CALL WETAIR(TA,WB,RELHUM,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,
     * DENAIR,CP,WTRPOT)
      VDAIR = VD
      RHUM=100.
      CALL WETAIR(TSKIN,WB,RHUM,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,
     * DENAIR,CP,WTRPOT)
      VDSURF = VD
      
      IF((GEOMETRY.EQ.4).OR.(GEOMETRY.EQ.2))THEN
       S2=((ASQ*BSQ*CSQ)/(ASQ*BSQ+ASQ*CSQ+BSQ*CSQ))
       ENB=HC*CONVAR*(TSKIN-TA)+SIG*EMISAN*CONVAR*((TSKIN+273.15)**4.
     & -(TRAD+273.15)**4.)+HD*CONVAR*SKINW*HTOVPR*(VDSURF-VDAIR)-
     & QSOLAR-(2.*FLSHCOND*VOL*(TC-TSKIN))/S2
      ELSE
       ENB=HC*CONVAR*(TSKIN-TA)+SIG*EMISAN*CONVAR*((TSKIN+273.15)**4.
     & -(TRAD+273.15)**4.)+HD*CONVAR*SKINW*HTOVPR*(VDSURF-VDAIR)-
     & QSOLAR-(4.*FLSHCOND*VOL*(TC-TSKIN))/RFLESH**2.
      ENDIF

      FUNSKIN = ENB
      
      RETURN
      END