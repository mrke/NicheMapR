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

C	  THIS FUNCTION IS USED FOR TRANSIENT HEAT BUDGET CALCULATIONS TO COMPUTE
C     SKIN TEMPERATURE - IT IS CALLED BY DSUB AND SOLVED VIA ZBRENT

      Implicit None

      DOUBLE PRECISION A,A1,A2,A3,A4,A4b,A5,A6,ABSMAX,ABSMIN,AEFF,AirVol
      DOUBLE PRECISION AL,ALT,AMASS,ANDENS,Area,aref,ASEMAJR,ASIL,ASILP
      DOUBLE PRECISION ASQ,ATOT,B,BP,bref,BSEMINR,BSQ,C,CO2MOL,convar,CP
      DOUBLE PRECISION cref,CSEMINR,CSQ,customgeom,CUTFA,DELTAR,DENAIR
      DOUBLE PRECISION DEPSUB,DP,E,EMISAN,EMISSB,EMISSK,Enb,ESAT,EXTREF
      DOUBLE PRECISION F12,f13,f14,f15,f16,F21,f23,f24,f25,f26,f31,F32
      DOUBLE PRECISION f41,F42,f51,F52,f61,Fatcond,FATOSB,FATOSK
      DOUBLE PRECISION Flshcond,FLTYPE,FLUID,flymetab,flyspeed,flytime
      DOUBLE PRECISION FUNskin,G,GEVAP,GN,H2O_BalPast,HC,HD,HR,HTOVPR
      DOUBLE PRECISION MR_1,MR_2,MR_3,O2MAX,O2MIN,PCTEYE,peyes,phi
      DOUBLE PRECISION phimax,phimin,PI,PTCOND,ptcond_orig,QCOND,QCONV
      DOUBLE PRECISION Qgenet,QIRIN,QIROUT,QMETAB,QRESP,Qsevap,QSOLAR
      DOUBLE PRECISION QSOLR,QSWEAT,R,R1,RELHUM,Rflesh,rho1_3,RHUM
      DOUBLE PRECISION Rinsul,RQ,Rskin,RW,S1,S2,shade,shp,sidex,SIG
      DOUBLE PRECISION SkinT,SkinW,SPHEAT,SUBTK,TA,TAVE,tbask,TC,Tdigpr
      DOUBLE PRECISION temerge,Tlung,Tmaxpr,Tminpr,TOBJ,TPREF,TQSOL,TR
      DOUBLE PRECISION TRAD,trans1,Tskin,TSKY,TSUBST,TVINC,TVIR,TWING,VD
      DOUBLE PRECISION VDAIR,VDSURF,VEL,VOL,WB,WCUT,WEVAP,WEYES,WQSOL
      DOUBLE PRECISION WRESP,WTRPOT,X,Xtry

      INTEGER climbing,DEB1,flight,flyer,flytest,geometry,IHOUR,LIVE
      INTEGER MICRO,nodnum,wingcalc,wingmod

      dimension customgeom(8),shp(3)
      
      COMMON/ANPARMS/Rinsul,R1,Area,VOL,Fatcond
      COMMON/Behav2/geometry,nodnum,customgeom,shp
      common/climb/climbing
      COMMON/ELLIPS/ASEMAJR,BSEMINR,CSEMINR
      COMMON/EVAP1/PCTEYE,WEYES,WRESP,WCUT,AEFF,CUTFA,HD,peyes,SkinW,
     & SkinT,HC,convar
      COMMON/fly/flytime,flyspeed,flymetab,flight,flyer,flytest
      COMMON/FUN1/QSOLAR,QIRIN,QMETAB,QRESP,QSEVAP,QIROUT,QCONV,QCOND
      COMMON/FUN2/AMASS,RELHUM,ATOT,FATOSK,FATOSB,EMISAN,SIG,Flshcond
      COMMON/FUN3/AL,TA,VEL,PTCOND,SUBTK,DEPSUB,TSUBST,ptcond_orig
      COMMON/FUN4/Tskin,R,WEVAP,TR,ALT,BP,H2O_BalPast
      COMMON/FUN6/SPHEAT,ABSMAX,ABSMIN,O2MAX,O2MIN,LIVE
      Common/Guess/Xtry
      COMMON/REVAP1/Tlung,DELTAR,EXTREF,RQ,MR_1,MR_2,MR_3,DEB1
      COMMON/REVAP2/GEVAP,AirVol,CO2MOL
      Common/soln/Enb
      COMMON/TPREFR/TMAXPR,TMINPR,TDIGPR,TPREF,tbask,temerge
      Common/Treg/Tc
      COMMON/WCONV/FLTYPE
      COMMON/WDSUB1/ANDENS,ASILP,EMISSB,EMISSK,FLUID,G,IHOUR
      COMMON/WDSUB2/QSOLR,TOBJ,TSKY,MICRO
      COMMON/WINGFUN/rho1_3,trans1,aref,bref,cref,phi,F21,f31,f41,f51,
     & sidex,WQSOL,phimin,phimax,twing,F12,F32,F42,F52,f61,TQSOL,A1,A2,
     & A3,A4,A4b,A5,A6,f13,f14,f15,f16,f23,f24,f25,f26,wingcalc,wingmod
      COMMON/WMET/QSWEAT
      COMMON/WSOLAR/ASIL,Shade


      DATA PI/3.14159265/

C     THE GUESSED VARIABLE, X, IS CORE TEMPERATURE (C); SEE SUB. MET
C     FOR DETAILED EXPLANATION OF CALCULATION OF SURF. TEMP., Tskin,
C     FROM TC AND MASS
c     This assumes uniform body temperature.

C     Control of body temperature guesses for stability purposes
      If (X .gt. 80.) then
       X = 80.
      else
       If (X .lt. -20.0) then
        X = Tsky + 0.1
       Endif
      Endif

      Tskin = X
      Xtry = X

c     Get the metabolic rate
C     CHECKING FOR INANIMATE OBJECT
      IF (LIVE .EQ. 0) THEN
C      Inanimate
       QMETAB = 0.0
       Tskin = X
      ELSE
c      Alive, but is it too cold?
       If (Tc .ge. 0.0) then
        CALL MET
       else
C       Too cold, super low metabolism
        Qmetab = 0.0001
        Tskin = X
       Endif
      ENDIF

c     Get the respiratory water loss
C     Checking for fluid type
      if (FLTYPE .EQ. 0.00) THEN
C      AIR
C      Call for respiratory water & energy loss
       If (Qmetab .ge. 0.000) then
        CALL RESP
       else
c       Negative metabolic rate. No physiological meaning - dead.
        Qresp = 0.00000
        Qmetab = 0.00000
       endif
      endif

C     Net internal heat generation
      Qgenet = Qmetab - Qresp
C     Net internal heat generation/unit volume. Use for estimating skin temp.
      Gn = Qgenet/Vol
      IF (LIVE .EQ. 0) THEN
       GN = 0.
      ENDIF

C     COMPUTING SURFACE TEMPERATURE AS DICTATED BY GEOMETRY

C     First set average body temperature for estimation of avearage lung temperature
      IF(geometry.eq.1)then
C      Cylinder: From p. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
C      Tave = (gR**2/(8k)) + Tskin, where Tskin = Tcore - gR**2/(4k)
C      Note:  these should all be solved simultaneously.  This is an approximation
C      using cylinder geometry. Subcutaneous fat is allowed in cylinder & sphere
C      calculations.
       Rflesh = R1 - Rinsul
C      Computing average torso temperature from core to skin
       Tlung = (Gn*Rflesh**2)/(8.*Flshcond) + Tskin
      endif
      if(geometry.eq.2)then
C      Ellipsoid: Derived 24 October, 1993  W. Porter
       A = ASEMAJR
       B = BSEMINR
       C = CSEMINR
       ASQ = A**2
       BSQ = B**2
       CSQ = C**2
C      Computing average torso temperature from core to skin
       Tlung = (Gn/(4.*Flshcond)) * ((ASQ*BSQ*CSQ)/
     &  (ASQ*BSQ+ASQ*CSQ+BSQ*CSQ)) + Tskin
      endif
      if(geometry.eq.4)then
C      Sphere:
       RFLESH = R1 - RINSUL
       RSKIN = R1
C      FAT LAYER, IF ANY
       S1 = (QGENET/(4.*PI*Flshcond))*((RFLESH - RSKIN)/(RFLESH*RSKIN))
C      COMPUTING AVERAGE TORSO TEMPERATURE FROM CORE TO SKIN (12 BECAUSE TLUNG IS 1/2 THE TC-TSKIN DIFFERENCE, 6*AK1)
       TLUNG = (GN*RFLESH**2)/(12.*Flshcond) + TSKIN
      endif
      IF((geometry.eq.3).OR.(geometry.EQ.5))then
C      Model lizard as cylinder
C      Cylinder: From p. 270 Bird, Stewart & Lightfoot. 1960. Transport Phenomena.
C      Tave = (gR**2/(8k)) + Tskin, where Tskin = Tcore - gR**2/(4k)
C      Note:  these should all be solved simultaneously.  This is an approximation
C      using cylinder geometry. Subcutaneous fat is allowed in cylinder & sphere
C      calculations.
       Rflesh = R1 - Rinsul
C      Computing average torso temperature from core to skin
       Tlung = (Gn*Rflesh**2)/(8.*Flshcond) + Tskin
      endif
C     Limiting lung temperature extremes
      if (Tlung .gt. Tc) then
       Tlung = Tc
      endif
      if (Tlung .lt. -3.) then
       Tlung = -3.
      endif

      CALL CONV
      CALL RESP
      CALL Sevap
      CALL RADOUT
      CALL COND

      if (FLTYPE .EQ. 1.00) THEN
C      WATER Environment
       Qsevap = 0.00
       WEVAP = 0.0
       QIRIN = 0.0
       QIROUT = 0.00
       QCOND = 0.00
      ENDIF
      
      Trad=(QIRIN/(EMISAN*1*ATOT*SIG))**(1./4.)-273.15
      HTOVPR = 2.5012E+06 - 2.3787E+03 * Ta
      Tave=(Trad+Tskin)/2.
      HR=4*EMISAN*SIG*(Tave+273)**3
      
C     INITIALIZING FOR SUB. WETAIR2
      WB = 0.0
      DP = 999.
      CALL WETAIR2(TA,WB,RELHUM,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,
     * DENAIR,CP,WTRPOT)
      VDAIR = VD
      RHUM=100.
      CALL WETAIR2(TSKIN,WB,RHUM,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,
     * DENAIR,CP,WTRPOT)
      VDSURF = VD
      
      if((geometry.eq.4).or.(geometry.eq.2))then
       S2=((ASQ*BSQ*CSQ)/(ASQ*BSQ+ASQ*CSQ+BSQ*CSQ))
       Enb=hc*CONVAR*(Tskin-Ta)+SIG*emisan*CONVAR*((Tskin+273)**4
     & -(Trad+273)**4)+hD*CONVAR*(skinw/100)*HTOVPR*(VDSURF-VDAIR)-
     & QSOLAR-(2*Flshcond*VOL*(Tc-Tskin))/S2
      else
       Enb=hc*CONVAR*(Tskin-Ta)+SIG*emisan*CONVAR*((Tskin+273)**4
     & -(Trad+273)**4)+hD*CONVAR*(skinw/100)*HTOVPR*(VDSURF-VDAIR)-
     & QSOLAR-(4*Flshcond*VOL*(Tc-Tskin))/Rflesh**2
      endif

      FUNskin = Enb
      
      RETURN
      END
