      SUBROUTINE CLIMBER
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

C	  THIS PROGRAM IS CALLED BY THERMO WHEN ANIMAL IS CLIMBING TO ESCAPE
C     HIGH OR LOW TEMPERATURES

      Implicit None

      DOUBLE PRECISION A1,A2,A3,A4,A4b,A5,A6,Acthr,AL,AMASS,Andens,aref
      DOUBLE PRECISION Asil,Asilp,ATOT,bref,cref,customgeom,Depsel
      DOUBLE PRECISION Depsub,EMISAN,Emissb,Emissk,F12,f13,f14,f15,f16
      DOUBLE PRECISION F21,f23,f24,f25,f26,f31,F32,f41,F42,f51,F52,f61
      DOUBLE PRECISION FATOSB,FATOSK,Flshcond,Fluid,G,HSHSOI,HSOIL
      DOUBLE PRECISION MSHSOI,MSOIL,newdep,phi,phimax,phimin,PSHSOI
      DOUBLE PRECISION PSOIL,Ptcond,ptcond_orig,QCOND,QCONV,QIRIN,QIROUT
      DOUBLE PRECISION QMETAB,QRESP,QSEVAP,Qsol,QSOLAR,Qsolr,RELHUM,RH
      DOUBLE PRECISION rho1_3,rhref,Shade,shp,sidex,SIG,Subtk,Ta,Taloc
      DOUBLE PRECISION Tannul,tbask,Tc,Tcores,Tdigpr,temerge,Time,Tmaxpr
      DOUBLE PRECISION Tminpr,Tobj,TPREF,TQSOL,trans1,Tref,Tshsoi,Tsky
      DOUBLE PRECISION TskyC,Tsoil,Tsub,Tsubst,TWING,Vel,VLOC,Vref,WQSOL
      DOUBLE PRECISION Z,Zsoil
      
      Integer climbing,Ihour,geometry,Micro,Nodnum,wingmod,wingcalc
      
      CHARACTER*1 Burrow,CkGrShad,Climb,Crepus,Dayact,Nocturn
      
      DIMENSION Acthr(25),customgeom(8),Depsel(25),HSHSOI(25),HSOIL(25)
      DIMENSION MSHSOI(25),MSOIL(25),PSHSOI(25),PSOIL(25),QSOL(25)
      DIMENSION RH(25),rhref(25),shp(3),Taloc(25),Tcores(25),Time(25)
      DIMENSION TREF(25),TSHSOI(25),TskyC(25),TSOIL(25),TSUB(25)
      DIMENSION VLOC(25),VREF(25),Z(25),ZSOIL(10)

      COMMON/Behav1/Dayact,Burrow,Climb,CkGrShad,Crepus,Nocturn
      COMMON/Behav2/geometry,nodnum,customgeom,shp
      Common/Behav3/Acthr
      common/climb/climbing
      COMMON/DEPTHS/DEPSEL,Tcores
      COMMON/ENVAR1/QSOL,RH,TskyC,TIME,Taloc,TREF,rhref
      COMMON/ENVAR2/TSUB,VREF,Z,Tannul,VLOC
      COMMON/FUN1/QSOLAR,QIRIN,QMETAB,QRESP,QSEVAP,QIROUT,QCONV,QCOND
      COMMON/FUN2/AMASS,RELHUM,ATOT,FATOSK,FATOSB,EMISAN,SIG,Flshcond
      COMMON/FUN3/AL,TA,VEL,PTCOND,SUBTK,DEPSUB,TSUBST,ptcond_orig
      COMMON/SOIL/TSOIL,TSHSOI,ZSOIL,MSOIL,MSHSOI,PSOIL,PSHSOI,HSOIL,
     & HSHSOI
      COMMON/TPREFR/TMAXPR,TMINPR,TDIGPR,TPREF,tbask,temerge
      Common/Treg/Tc
      COMMON/WDSUB1/ANDENS,ASILP,EMISSB,EMISSK,FLUID,G,IHOUR
      COMMON/WDSUB2/QSOLR,TOBJ,TSKY,MICRO
      COMMON/WINGFUN/rho1_3,trans1,aref,bref,cref,phi,F21,f31,f41,f51
     & ,sidex,WQSOL,phimin,phimax,twing,F12,F32,F42,F52,f61,TQSOL,A1,A2
     & ,A3,A4,A4b,A5,A6,f13,f14,f15,f16,f23,f24,f25,f26,wingcalc,wingmod
      COMMON/WSOLAR/ASIL,Shade
      common/wundrg/newdep

      if ((Tref(Ihour) .ge. Tminpr).and.(Tref(Ihour) .le. Tpref)) THEN
C      It's OK up high
c      if (Tref(Ihour) .LE. TDIGPR)THEN
c        Climb to reference height - usually 1.2 to 2m, reporting as 150 cm height
          DEPSEL(IHOUR) = 150.0
c        Reference height air temperature is the same in sun or shade.
          Ta = Tref(Ihour)
          Vel = Vref(ihour)
          climbing = 1
c        No assumption of increase in shade due to climbing.  May still be in the sun.
c        endif
      endif
      RETURN
      END
