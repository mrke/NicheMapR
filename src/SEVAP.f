      SUBROUTINE SEVAP

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

C     THIS SUBROUTINE COMPUTES SKIN EVAPORATION BASED ON THE MASS TRANSFER
C     COEFFICIENT, % OF SURFACE OF THE SKIN ACTING AS A FREE WATER SURFACE
C     AND EXPOSED TO THE AIR, AND THE VAPOR DENSITY GRADIENT BETWEEN THE
C     SURFACE AND THE AIR, EACH AT THEIR OWN TEMPERATURE.

      IMPLICIT NONE

      DOUBLE PRECISION A1,A2,A3,A4,A4B,A5,A6,ABSMAX,ABSMIN,AEFF,AIRVOL
      DOUBLE PRECISION AL,ALT,AMASS,ANDENS,AREF,ASILP,ATOT,BP,BREF
      DOUBLE PRECISION CO2MOL,CONVAR,CP,CREF,CUSTOMGEOM,CUTFA,DB,DELTAR
      DOUBLE PRECISION DENAIR,DEPSUB,DP,E,EGGSHP,EGGPTCOND,EMISAN,EMISSB
      DOUBLE PRECISION EMISSK,ESAT,EXTREF,F12,F13,F14,F15,F16,F21,F23
      DOUBLE PRECISION F24,F25,F26,F31,F32,F41,F42,F51,F52,F61,FATOSB
      DOUBLE PRECISION FATOSK,FLSHCOND,FLUID,G,GEVAP,H2O_BALPAST,HC,HD
      DOUBLE PRECISION HDFORC,HDFREE,HTOVPR,MR_1,MR_2,MR_3,MW,PANT
      DOUBLE PRECISION PANTMAX,PATMOS,PEYES,PHI,PHIMAX,PHIMIN,PMOUTH
      DOUBLE PRECISION PRES,PSI_BODY,PSTD,PTCOND,PTCOND_ORIG,QCOND
      DOUBLE PRECISION QCONV,QIRIN,QIROUT,QMETAB,QRESP,QSEVAP,QSOLAR,R
      DOUBLE PRECISION RAINFALL,RELHUM,RG,RH,RHO1_3,RQ,RW,SHP,SIDEX,SIG
      DOUBLE PRECISION SKINT,SKINW,SPHEAT,SUBTK,TA,TAIR,TBASK,TC,TDIGPR
      DOUBLE PRECISION TEMERGE,TLUNG,TMAXPR,TMINPR,TPREF,TQSOL,TR,TRANS1
      DOUBLE PRECISION TSKIN,TSUBST,TVINC,TVIR,TWING,V,VD,VDAIR,VDSURF
      DOUBLE PRECISION VEL,POT,WATER,WB,WCUT,WEVAP,WEYES,WQSOL,WRESP
      DOUBLE PRECISION WTRPOT,XTRY,CONV_ENHANCE,G_VS_AB,G_VS_AD,V_M
      DOUBLE PRECISION EAIR,ESURF,G_VA,G_V,ESATSURF
      
      INTEGER DEB1,IDAY,IHOUR,GEOMETRY,LIVE,NODNUM,WINGCALC,WINGMOD,LEAF

      DIMENSION CUSTOMGEOM(8),SHP(3),EGGSHP(3),PRES(25)

      COMMON/BEHAV2/GEOMETRY,NODNUM,CUSTOMGEOM,SHP,EGGSHP
      COMMON/DAYITR/IDAY
      COMMON/EVAP1/WEYES,WRESP,WCUT,AEFF,CUTFA,HD,PEYES,SKINW,G_VS_AB,
     & G_VS_AD,SKINT,HC,CONVAR,PMOUTH,PANT,PANTMAX,CONV_ENHANCE,LEAF
      COMMON/EVAP2/HDFREE,HDFORC
      COMMON/FUN1/QSOLAR,QIRIN,QMETAB,QRESP,QSEVAP,QIROUT,QCONV,QCOND
      COMMON/FUN2/AMASS,RELHUM,ATOT,FATOSK,FATOSB,EMISAN,SIG,FLSHCOND
      COMMON/FUN3/AL,TA,VEL,PTCOND,SUBTK,DEPSUB,TSUBST,PTCOND_ORIG,
     & EGGPTCOND,POT
      COMMON/FUN4/TSKIN,R,WEVAP,TR,ALT,BP,H2O_BALPAST
      COMMON/FUN6/SPHEAT,ABSMAX,ABSMIN,LIVE
      COMMON/GUESS/XTRY
      COMMON/PRESSURE/PRES
      COMMON/RAINFALLS/RAINFALL
      COMMON/REVAP1/TLUNG,DELTAR,EXTREF,RQ,MR_1,MR_2,MR_3,DEB1
      COMMON/REVAP2/GEVAP,AIRVOL,CO2MOL
      COMMON/TPREFR/TMAXPR,TMINPR,TDIGPR,TPREF,TBASK,TEMERGE
      COMMON/TREG/TC
      COMMON/WATERPOT/PSI_BODY 
      COMMON/WDSUB1/ANDENS,ASILP,EMISSB,EMISSK,FLUID,G,IHOUR
      COMMON/WINGFUN/RHO1_3,TRANS1,AREF,BREF,CREF,PHI,F21,F31,F41,F51
     &,SIDEX,WQSOL,PHIMIN,PHIMAX,TWING,F12,F32,F42,F52
     &,F61,TQSOL,A1,A2,A3,A4,A4B,A5,A6,F13,F14,F15,F16,F23,F24,F25,F26
     &,WINGCALC,WINGMOD
     
      TAIR=TA
      V=VEL
      XTRY=TC
      MW = 0.018 ! molar mass of water, kg/mol
      RG = 8.314 ! gas constant, J/mol/K
      V_m = 44.6 ! molar volume of air, mol/m3
C     CALCULATING SKIN SURFACE SATURATION VAPOR DENSITY
C      RH = 100.
      RH=exp(PSI_BODY/(RG/MW*(TSKIN+273.15)))*100. ! FROM JULY 2022 MAKING SKIN RH A FUNCTION OF BODY WATER POTENTIAL 
      DB=TSKIN
      IF (LIVE .EQ. 1) THEN
C     CHECK FOR TOO LOW A SURFACE TEMPERATURE
       IF(TSKIN.LT.0.) THEN
        DB=0.
       ENDIF
      ENDIF

C     SETTING 3 PARAMETERS FOR WETAIR, SINCE RH IS KNOWN (SEE WETAIR LISTING)
      WB=0.
      DP=999.
C     BP CALCULATED FROM ALTITUDE USING THE STANDARD ATMOSPHERE
C     EQUATIONS FROM SUBROUTINE DRYAIR    (TRACY ET AL,1972)
      PSTD=101325.
C     PATMOS=PSTD*((1.-(.0065*ALT/288.))**(1./.190284))
      PATMOS=PRES(IHOUR)
      BP=PATMOS

      CALL WETAIR(DB,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,CP,
     * WTRPOT)
      VDSURF=VD
      ESURF=E
      ESATSURF=ESAT
      
C     AIR VAPOR DENSITY
      RH=RELHUM
      DB=TAIR
      CALL WETAIR(DB,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,CP,
     * WTRPOT)
      VDAIR=VD
      EAIR=E
      
C     CHECKING FOR LIVING OBJECTS
      IF (LIVE .EQ. 1) THEN
C      OCULAR WATER LOSS
C      CHECKING FOR OPEN EYES (ACTIVE)
       IF((TC.GE.TBASK).AND.(TC.LE.TMAXPR))THEN
C       EYES OPEN
        WEYES=HD*PEYES*ATOT*(VDSURF-VDAIR)
       ELSE
C       EYES CLOSED AND RESTING
        WEYES = 0.0
       ENDIF
       WRESP = GEVAP/1000.
      ELSE
       WEYES=0.0
       WRESP=0.0
      ENDIF
C     END OF LIVE VS INANIMATE

      IF(LIVE.EQ.0)THEN
C      INANIMATE
       IF(LEAF.EQ.0)THEN
        WCUT=AEFF*HD*(VDSURF-VDAIR)
       ELSE
        G_VA=HD*V_M !BOUNDARY CONDUCTANCE, mol/m2/s
        G_V=(0.5*g_vs_ab*g_va)/(g_vs_ab+g_va)+(0.5*g_vs_ad*g_va)/
     &  (g_vs_ad+g_va) ! vapour conductance, mol/m2/s 
        WCUT=0.018015*G_V*(ESATSURF - EAIR)/BP*CONVAR ! kg/s
       ENDIF
       WATER=WCUT
       GO TO 10
      ELSE
C      ANIMATE, CALCULATE BELOW
      ENDIF

      IF(WEYES.GT.0.)THEN
       WCUT=(AEFF-PEYES*ATOT*SKINW)*HD*(VDSURF-VDAIR)
      ELSE
       WCUT=AEFF*HD*(VDSURF-VDAIR)
      ENDIF
      WATER=WEYES+WRESP+WCUT
C     END OF COMPUTING AEFF FOR SURFACE OR NOT

   10 CONTINUE

C     FROM DRYAIR: LATENT HEAT OF VAPORIZATION
      HTOVPR=2.5012E+06-2.3787E+03*TAIR
      QSEVAP=(WEYES+WCUT)*HTOVPR ! FIXED 16/10/2021 - WAS WATER*HTOVPR BUT THIS INCLUDES RESP AND SO WAS DOUBLE COUNTED

C     KG/S TO G/S
      WEYES=WEYES*1000.
      WRESP=WRESP*1000.
      WCUT=WCUT*1000.
      WEVAP=WATER*1000.

      RETURN
      END