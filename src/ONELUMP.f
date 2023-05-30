      SUBROUTINE ONELUMP (N,T,Y,DDSUB,RPAR,IPAR)

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

      IMPLICIT NONE

      EXTERNAL FUNSKIN,ISPLINE

C     THIS PROGRAM USES THE SOLUTION FOR UNIFORM INTERNAL HEAT GENERATION
C     FOR ALL GEOMETRIES, WHICH YIELDS A PARABOLIC TEMPERATURE PROFILE INTERNALLY.
C     THERE IS NO BREAKING THE ANIMAL INTO INTERNAL NODES ANY MORE.
C     AVERAGE INTERNAL TEMPERATURE IS THE AVERAGE OVER THE INTEGRAL FROM CORE TO SKIN.

C     PROGRAM FOR THE DOPRI INTEGRATOR TO SOLVE THE ONELUMP TRANSIENT HEAT BUDGET.

C     N = NUMBER OF EQUATIONS (NODES) FOR ANIMAL
C     T = TIME (MINUTES)
C     Y = TCORE (C)
C     DDSUB = DTCORE/DTIME

      DOUBLE PRECISION A1,A2,A3,A4,A4B,A5,A6,ABSAN,ABSMAX,ABSMIN,ABSSB
      DOUBLE PRECISION AEFF,AHEIT,AIRVOL,AL,ALENTH,ALT,AMASS,ANDENS
      DOUBLE PRECISION AREA,AREF,ASEMAJR,ASIL,ASILN,ASILP,ASQ,AT,ATOT,AV
      DOUBLE PRECISION AWIDTH,BP,BREF,BSEMINR,BSQ,CO2MOL,CONTDEP
      DOUBLE PRECISION CONTDEPTH,CONTH,CONTHOLE,CONTVOL,CONTW,CONTWET
      DOUBLE PRECISION CONVAR,CP,CREF,CSEMINR,CSQ,CUSTOMGEOM,CUTFA
      DOUBLE PRECISION DEBQMET,DEBQMET_INIT,DELTAR,DENAIR,DEPSUB,DP
      DOUBLE PRECISION DRYFOOD,E,EGGSHP,EMISAN,EMISSB,EMISSK,ENARY1
      DOUBLE PRECISION ENARY10,ENARY11,ENARY12,ENARY13,ENARY14,ENARY15
      DOUBLE PRECISION ENARY16,ENARY17,ENARY18,ENARY19,ENARY2,ENARY20
      DOUBLE PRECISION ENARY21,ENARY22,ENARY23,ENARY24,ENARY25,ENARY26
      DOUBLE PRECISION ENARY27,ENARY28,ENARY29,ENARY3,ENARY30,ENARY31
      DOUBLE PRECISION ENARY32,ENARY33,ENARY34,ENARY35,ENARY36,ENARY37
      DOUBLE PRECISION ENARY38,ENARY39,ENARY4,ENARY40,ENARY41,ENARY42
      DOUBLE PRECISION ENARY43,ENARY44,ENARY45,ENARY46,ENARY47,ENARY48
      DOUBLE PRECISION ENARY5,ENARY6,ENARY7,ENARY8,ENARY9,ENBERR,ESAT
      DOUBLE PRECISION EXTREF,F12,F13,F14,F15,F16,F21,F23,F24,F25,F26
      DOUBLE PRECISION F31,F32,F41,F42,F51,F52,F61,FAECES,FATCOND,FATOBJ
      DOUBLE PRECISION FATOSB,FATOSK,FLSHCOND,FLUID,G,GEVAP,GH2OMET
      DOUBLE PRECISION GH2OMET_INIT,GN,H2O_BALPAST,HC,HD,HR,HRN
      DOUBLE PRECISION HTOVPR,J,KTC,MASSLEFT,MLO2,MLO2_INIT,MR_1,MR_2
      DOUBLE PRECISION MR_3,NWASTE,PANT,PANTMAX,PDIF,PEYES,PHI,PHIMAX
      DOUBLE PRECISION PHIMIN,PI,PMOUTH,PSI_BODY,PTCOND,PTCOND_ORIG
      DOUBLE PRECISION QCOND,QCONV,QGENET,QIN,QIRIN,QIROUT,QMETAB,QOUT
      DOUBLE PRECISION QRESP,QSEVAP,QSOL,QSOLAR,QSOLR,QST,R,R1,R2,REFTOL
      DOUBLE PRECISION RELHUM,RH,RHO1_3,RHREF,RHUM,RINSUL,RQ,RW,S2,SHADE
      DOUBLE PRECISION SHP,SIDEX,SIG,SKINT,SKINW,SPHEAT,SUBTK,T,TA,TALOC
      DOUBLE PRECISION TANNUL,TAVE,TB,TBASK,TC,TCINIT,TCPAST,TDIGPR
      DOUBLE PRECISION TEMERGE,TESTX,TI,TIME,TLUNG,TMAXPR,TMINPR,TOBJ
      DOUBLE PRECISION TOTLEN,TPAST,TPREF,TPRINT,TQSOL,TR,TRAD,TRANS1
      DOUBLE PRECISION TRANSAR,TREF,TSKIN,TSKY,TSKYC,TSUB,TSUBST,TVINC
      DOUBLE PRECISION TVIR,TWING,VD,VDAIR,VDSURF,VEL,VLOC,VOL
      DOUBLE PRECISION VOLUMELEFT,VREF,WB,WC,WCUT,WEVAP,WEYES
      DOUBLE PRECISION WQSOL,WRESP,WTRPOT,X,X1,X2,Y,Z,ZBRENT,ZEN
      DOUBLE PRECISION RAINMULT,RPAR,DDSUB,ISPLINE,EGGPTCOND,POT
      DOUBLE PRECISION CONV_ENHANCE,G_VS_AB,G_VS_AD


      INTEGER CONTONLY,CONTYPE,DEB1,GEOMETRY,IDAY,IHOUR,IPAR,IT,JP
      INTEGER LIVE,MICRO,NM,NODNUM,POND,SCENAR,WETMOD,WINGCALC,LEAF
      INTEGER WINGMOD
      INTEGER RAINHOUR,N
      INTEGER, PARAMETER :: nn=25      ! base points for interpolation
      DOUBLE PRECISION,DIMENSION (NN)::B(NN),C(NN),D(NN)      
      
      LOGICAL SUCCES

      DIMENSION Y(N),DDSUB(N),RPAR(35),IPAR(31)
      
      DIMENSION CUSTOMGEOM(8),DEBQMET(24),DRYFOOD(24),EGGSHP(3)
      DIMENSION ENARY10(25),ENARY11(25),ENARY12(25),ENARY13(25)
      DIMENSION ENARY14(25),ENARY15(25),ENARY16(25),ENARY17(25)
      DIMENSION ENARY18(25),ENARY19(25),ENARY2(25),ENARY20(25)
      DIMENSION ENARY21(25),ENARY22(25),ENARY23(25),ENARY24(25)
      DIMENSION ENARY25(25),ENARY26(25),ENARY27(25),ENARY28(25)
      DIMENSION ENARY29(25),ENARY3(25),ENARY30(25),ENARY31(25)
      DIMENSION ENARY32(25),ENARY33(25),ENARY34(25),ENARY35(25)
      DIMENSION ENARY36(25),ENARY37(25),ENARY38(25),ENARY39(25)
      DIMENSION ENARY4(25),ENARY40(25),ENARY41(25),ENARY42(25)
      DIMENSION ENARY43(25),ENARY44(25),ENARY45(25),ENARY46(25)
      DIMENSION ENARY47(25),ENARY48(25),ENARY5(25),ENARY6(25)
      DIMENSION ENARY7(25),ENARY8(25),ENARY9(25),FAECES(24),ENARY1(25)
      DIMENSION GH2OMET(24),MLO2(24),NWASTE(24),QSOL(25),RH(25),HRN(25)
      DIMENSION RHREF(25),SHP(3),TALOC(25),TI(25),TIME(25),TRANSAR(5,25)
      DIMENSION TREF(25),TSKYC(25),TSUB(25),VLOC(25),VREF(25)
      DIMENSION Z(25)

      COMMON/ANPARMS/RINSUL,R1,AREA,VOL,FATCOND
      COMMON/BEHAV2/GEOMETRY,NODNUM,CUSTOMGEOM,SHP,EGGSHP
      COMMON/CONT/CONTH,CONTW,CONTVOL,CONTDEP,CONTHOLE,CONTWET,RAINMULT,
     & WETMOD,CONTONLY,CONTYPE,RAINHOUR
      COMMON/CONTDEPTH/CONTDEPTH
      COMMON/DAYITR/IDAY
      COMMON/DEBRESP/MLO2,GH2OMET,DEBQMET,MLO2_INIT,GH2OMET_INIT,
     & DEBQMET_INIT,DRYFOOD,FAECES,NWASTE
      COMMON/DIMENS/ALENTH,AWIDTH,AHEIT
      COMMON/DSUB1/ENARY1,ENARY2,ENARY3,ENARY4,ENARY9,ENARY10,ENARY11,
     & ENARY12,ENARY17,ENARY18,ENARY19,ENARY20,ENARY21,ENARY22,ENARY23,
     & ENARY24,ENARY25,ENARY26,ENARY27,ENARY28,ENARY45,ENARY46,ENARY47,
     & ENARY48
      COMMON/DSUB2/ENARY5,ENARY6,ENARY7,ENARY8,ENARY13,ENARY14,ENARY15,
     & ENARY16,ENARY29,ENARY30,ENARY31,ENARY32,ENARY33,ENARY34,ENARY35,
     & ENARY36,ENARY37,ENARY38,ENARY39,ENARY40,ENARY41,ENARY42,ENARY43,
     & ENARY44
      COMMON/ELLIPS/ASEMAJR,BSEMINR,CSEMINR
      COMMON/ENVAR1/QSOL,RH,TSKYC,TIME,TALOC,TREF,RHREF,HRN
      COMMON/ENVAR2/TSUB,VREF,Z,TANNUL,VLOC
      COMMON/EVAP1/WEYES,WRESP,WCUT,AEFF,CUTFA,HD,PEYES,SKINW,G_VS_AB,
     & G_VS_AD,SKINT,HC,CONVAR,PMOUTH,PANT,PANTMAX,CONV_ENHANCE,LEAF
      COMMON/FUN1/QSOLAR,QIRIN,QMETAB,QRESP,QSEVAP,QIROUT,QCONV,QCOND
      COMMON/FUN2/AMASS,RELHUM,ATOT,FATOSK,FATOSB,EMISAN,SIG,FLSHCOND
      COMMON/FUN3/AL,TA,VEL,PTCOND,SUBTK,DEPSUB,TSUBST,PTCOND_ORIG,
     & EGGPTCOND,POT
      COMMON/FUN4/TSKIN,R,WEVAP,TR,ALT,BP,H2O_BALPAST
      COMMON/FUN5/WC,ZEN,PDIF,ABSSB,ABSAN,ASILN,FATOBJ,NM
      COMMON/FUN6/SPHEAT,ABSMAX,ABSMIN,LIVE
      COMMON/OUTSUB/IT
      COMMON/PONDTEST/POND
      COMMON/REVAP1/TLUNG,DELTAR,EXTREF,RQ,MR_1,MR_2,MR_3,DEB1
      COMMON/REVAP2/GEVAP,AIRVOL,CO2MOL
      COMMON/SCENARIO/SCENAR
      COMMON/TPREFR/TMAXPR,TMINPR,TDIGPR,TPREF,TBASK,TEMERGE
      COMMON/TRANS/JP
      COMMON/TRANSIENT1/TCINIT,TRANSAR
      COMMON/TREG/TC
      COMMON/USROP2/ENBERR,TPRINT
      COMMON/WATERPOT/PSI_BODY      
      COMMON/WCOND/TOTLEN,AV,AT
      COMMON/WDSUB1/ANDENS,ASILP,EMISSB,EMISSK,FLUID,G,IHOUR
      COMMON/WDSUB2/QSOLR,TOBJ,TSKY,MICRO
      COMMON/WINGFUN/RHO1_3,TRANS1,AREF,BREF,CREF,PHI,F21,F31,F41,F51,
     & SIDEX,WQSOL,PHIMIN,PHIMAX,TWING,F12,F32,F42,F52,F61,TQSOL,A1,A2,
     & A3,A4,A4B,A5,A6,F13,F14,F15,F16,F23,F24,F25,F26,WINGCALC,WINGMOD
      COMMON/WSOLAR/ASIL,SHADE

      DATA TI/0.,60.,120.,180.,240.,300.,360.,420.,480.,540.,600.,660.,
     &720.,780.,840.,900.,960.,1020.,1080.,1140.,1200.,1260.,1320.,
     &1380.,1440./

      TC=Y(2)
      
C     INITIALISE
      PI = 3.14159265
      VOLUMELEFT=VOL
      MASSLEFT=AMASS
      ASQ=0.
      BSQ=0.
      CSQ=0.
      TPAST=T
      QIN=0.
      
      IF(CONTH.GT.0)THEN ! PROPERTIES OF WATER
       SPHEAT=4186.
       FLSHCOND=0.6
       PEYES=0.
      ENDIF
      
C     SPLINE MICROMET INPUT FOR CURRENT TIME
      CALL SPLINE(TI,ENARY1,B,C,D,NN)
      QSOLR=MAX(0.0D0,ISPLINE(T,TI,ENARY1,B,C,D,NN))
      CALL SPLINE(TI,ENARY2,B,C,D,NN)
      ZEN=MIN(MAX(0.0D0,ISPLINE(T,TI,ENARY2,B,C,D,NN)),90.0D0)*PI/180. ! RADIANS
      CALL SPLINE(TI,ENARY3,B,C,D,NN)
      TA=ISPLINE(T,TI,ENARY3,B,C,D,NN)
      CALL SPLINE(TI,ENARY4,B,C,D,NN)
      VEL=MAX(0.0D0,ISPLINE(T,TI,ENARY4,B,C,D,NN))
      CALL SPLINE(TI,ENARY5,B,C,D,NN)
      RELHUM=MIN(MAX(0.0D0,ISPLINE(T,TI,ENARY5,B,C,D,NN)),100.0D0)
      CALL SPLINE(TI,ENARY6,B,C,D,NN)
      TSUBST=ISPLINE(T,TI,ENARY7,B,C,D,NN)
      CALL SPLINE(TI,ENARY6,B,C,D,NN)
      TSKY=ISPLINE(T,TI,ENARY6,B,C,D,NN)
      CALL SPLINE(TI,ENARY6,B,C,D,NN)
      SHADE=MIN(MAX(0.0D0,ISPLINE(T,TI,ENARY8,B,C,D,NN)),100.0D0)   
  
C     CHECK THAT TEMPERATURE ISN'T GOING HAYWIRE KEARNEY ADDED THIS DURING THE FROG WORKSHOP 17/9/2012
      IF((TC.LT.-100).OR.(TC.GT.100))THEN
       TC=TA
       Y(2)=TC
      ENDIF
      TR = (TSUBST + TSKY)/2.
      TOBJ = TSUBST
C     INFRARED CALCULATIONS?
C     CHECKING TO SEE IF WATER ENVIRONMENT
      IF (FLUID .EQ. 0.0) THEN
C      AIR ENVIRONMENT.  TSKY & TSUBST STORED AND SET ACCORDING TO BEHAVIORS IN STEADY STATE RUN PRIOR TO THIS.
       CALL GEOM
       CALL RADIN
       CALL SOLAR
      ELSE
C      WATER ENVIRONMENT
       FATOBJ = 0.00
       FATOSK = 0.00
       FATOSB = 0.00
       QIRIN = 0.00
      ENDIF

C     QIN
C     GETTING HEAT FLUXES BASED ON CURRENT TEMPERATURES
C     ONE LUMP MODEL
      QIN = QSOLAR + QIRIN
C     M**3 = KG/(KG/M**3)

C     GETTING SURFACE HEAT FLUXES BASED ON CURRENT SKIN TEMP'S
C     AND GETTING RESP. EVAP LOSS BASED ON CURRENT CORE TEMP
      TB=TC
      IF(LIVE.EQ.0)THEN
       QMETAB = 0.
       QRESP = 0.
       QGENET = 0.
       GN = 0.
      ELSE
       IF(TC .GT. 50.)THEN
        TB = 50.
       ENDIF
       QMETAB = 0.0056*10.**(0.038*(TB)-1.771)*(AMASS*1000.)**.82
C      USING INTEGRATED AVERAGE OF BODY TEMPERATURES FOR FIRST ESTIMATE TO GET QRESP
       TLUNG = (TC + TSKIN)/2.
       CALL RESP
       QGENET = QMETAB-QRESP
       GN = QGENET/VOL
      ENDIF

      IF((GEOMETRY.EQ.1).OR.(GEOMETRY.EQ.3).OR.(GEOMETRY.EQ.5))THEN
C      USE CYLINDER: TC - TSKIN = (GENPV*R^2)/(4*FLSHCOND)
       TSKIN = TC - (GN*R1**2.)/(4.*FLSHCOND)
       TLUNG = (GN*R1**2.)/(8.*FLSHCOND) + TSKIN
      ENDIF
      IF((GEOMETRY.EQ.2).OR.(GEOMETRY.EQ.4))THEN
       ASQ = ASEMAJR**2.
       BSQ = BSEMINR**2.
       CSQ = CSEMINR**2.
       TSKIN = TC - (GN/(2.*FLSHCOND)) * ((ASQ*BSQ*CSQ)/
     & (ASQ*BSQ+ASQ*CSQ+BSQ*CSQ))
       TLUNG = (GN/(4.*FLSHCOND)) * ((ASQ*BSQ*CSQ)/
     & (ASQ*BSQ+ASQ*CSQ+BSQ*CSQ)) + TSKIN
      ENDIF

C     CHECKING FOR FLUID TYPE
      IF (FLUID .EQ. 0.00) THEN
C      AIR
       CALL CONV
       IF(LIVE.EQ.1)THEN
        CALL RESP
       ENDIF
       CALL SEVAP
       CALL RADOUT
       CALL COND
      ELSE
C      WATER
       QSEVAP = 0.D0
       WEVAP = 0.D0
       QIROUT = 0.D0
       QCOND = 0.D0
      ENDIF

C     QOUT
C     ONE LUMP MODEL
      QOUT = QCONV + QSEVAP + QIROUT + QRESP + QCOND

C     C/MIN = (J/S)*MIN*(S/MIN)/(KG/M^3*M^3*J/KG-C)
C     INTEGRATOR MULTIPLIES BY NUMBER OF MINUTES:YDOT(N) IN J/MIN
C      DDSUB(2)=((QIN + QGENET - QOUT)*60.)/(AMASS*SPHEAT)

      TRAD=(TSKY+TSUBST)/2.
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
       J=(1./(AMASS*SPHEAT))*((QSOLAR+QGENET+CONVAR*(((GN*S2)/(2.*
     & FLSHCOND))*(HC+HR)+HC*TA+HR*TRAD)+AEFF*(HD*HTOVPR*VDAIR))+
     & ((SUBTK*AV)/R1)*((GN*S2)/(2.*FLSHCOND)+TSUBST))
      ELSE
       R2=R1**2.
       J=(1./(AMASS*SPHEAT))*((QSOLAR+QGENET+CONVAR*(((GN*R2)/(4.*
     & FLSHCOND))*(HC+HR)+HC*TA+HR*TRAD)+AEFF*(HD*HTOVPR*VDAIR))+
     & ((SUBTK*AV)/R1)*((GN*R2)/(4.*FLSHCOND)+TSUBST))
      ENDIF

      KTC=(CONVAR*(TC*HC+TC*HR)+AEFF*HD*HTOVPR*VDSURF+TC*((SUBTK*AV)
     & /R1))/(AMASS*SPHEAT)


      DDSUB(2)=(J-KTC)*60. ! CONVERT RATE OF CHANGE FROM DEG/SECOND TO DEG/MINUTE

C     NOW GET THE SKIN TEMPERATURE WITH ZBRENT
      X1 = -50.
      X2 = 80.

C     GUESSING FOR CORE TEMPERATURE
      X = TA
      CALL ZBRAC(FUNSKIN,X1,X2,SUCCES)
C     INVOKING THE ENERGY BALANCE EQUATION VIA ZBRENT AND FUN

      TESTX = ENBERR
      REFTOL =TESTX
      X = ZBRENT(FUNSKIN,X1,X2,TESTX)
C     OUT COMES THE GUESSED VARIABLE VALUE (X) THAT SATISFIES
C     THE ENERGY BALANCE EQUATION
      TSKIN = X
      TLUNG = (TSKIN+TC)/2.
      CALL SEVAP
C     QSTORED (J/S) = KG/M^3* M^3* J/KG-C* C/(MIN*S/MIN)
      QST = AMASS*SPHEAT*DDSUB(2)/((T - TPAST)*60.)

      TPAST = T
      TCPAST = TC
      
C     CONTAINER/POND CALCULATIONS
      IF(CONTH.GT.0)THEN
       MASSLEFT=AMASS-MAX(WEVAP/1000.*60.,0.0D0) !WEVAP IS G/S SO CONVERT TO KG/MIN
       IF(MASSLEFT.LE.0.)THEN
        MASSLEFT=0.
        CONTDEP=0.
       ELSE
        VOLUMELEFT=(MASSLEFT/ANDENS)
       ENDIF
       CONTDEP=VOLUMELEFT/(PI*(CONTW/2./1000.)**2.)*1000.
       IF(CONTDEP.LT.0.01)THEN
        CONTDEP=0.01
       ENDIF
       IF(CONTDEP.GT.CONTH)THEN
        CONTDEP=CONTH
       ENDIF      
       CALL GEOM
      ENDIF
      
      END