      SUBROUTINE RADIN
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

c     Computes longwave radiation absorbed

      Implicit None

      double precision ABSAN,ABSSB,AL,AMASS,ANDENS,ASIL,ASILN,ASILP,ATOT
      double precision DEPSUB,EMISAN,EMISSB,EMISSK,QSOLR,RELHUM
      double precision FATOBJ,FATOSB,FATOSK,Flshcond,FLUID,FSKY,G
      double precision PTCOND,PCTDIF,QCOND,QCONV,QIRIN,QIROBJ,QIROUT
      double precision QIRSKY,QIRSUB,QMETAB,Qresp,Qsevap,QSOL,QSOLAR
      double precision RH,Shade,SIG,SUBTK
      double precision TA,Tannul,Taloc,Time,TKOBJ,TKSKY,ptcond_orig
      double precision TKSUB,TOBJ,Tref,Tsub,TSKY,TskyC,TSUBST
      double precision VEL,Vref,WC,Z,ZEN,VLOC,sidex,WQSOL
      double precision rho1_3,trans1,aref,bref,cref,phi,F21,f31,f41,f51
     &,phimin,phimax,twing,F12,F32,F42,F52,f23,f24,f25,f26
     &,f61,TQSOL,A1,A2,A3,A4,A4b,A5,A6,f13,f14,f15,f16
      double precision ir1,ir2,ir3,ir4,ir5,ir6,rhref

      INTEGER IHOUR,MICRO,NM,wingmod,wingcalc

      DIMENSION TIME(25),QSOL(25),RH(25),TskyC(25),rhref(25)
      DIMENSION Taloc(25),TREF(25),TSUB(25),VREF(25),Z(25)
      DIMENSION VLOC(25)

      COMMON/ENVAR1/QSOL,RH,TskyC,TIME,Taloc,TREF,rhref
      COMMON/ENVAR2/TSUB,VREF,Z,Tannul,VLOC
      COMMON/FUN1/QSOLAR,QIRIN,QMETAB,QRESP,QSEVAP,QIROUT,QCONV,QCOND
      COMMON/FUN2/AMASS,RELHUM,ATOT,FATOSK,FATOSB,EMISAN,SIG,Flshcond
      COMMON/FUN3/AL,TA,VEL,PTCOND,SUBTK,DEPSUB,TSUBST,ptcond_orig
      COMMON/FUN5/WC,ZEN,PCTDIF,ABSSB,ABSAN,ASILN,FATOBJ,NM
      COMMON/WINGFUN/rho1_3,trans1,aref,bref,cref,phi,F21,f31,f41,f51
     &,sidex,WQSOL,phimin,phimax,twing,F12,F32,F42,F52
     &,f61,TQSOL,A1,A2,A3,A4,A4b,A5,A6,f13,f14,f15,f16,f23,f24,f25,f26
     &,wingcalc,wingmod
      COMMON/WDSUB1/ANDENS,ASILP,EMISSB,EMISSK,FLUID,G,IHOUR
      COMMON/WDSUB2/QSOLR,TOBJ,TSKY,MICRO
      COMMON/WSOLAR/ASIL,Shade

      TKSKY=Tsky+273.15
      TKSUB=TSUBST+273.15

      if(wingmod.eq.2)then
c      currently assuming posture is parallel to the ground and that the 
c      view through surfaces 5 and 6 are half ground and half sky
       TKOBJ = TWING + 273.15
c      top of thorax IR
       IR1=EMISAN*SIG*A2*F12*TKOBJ**4
       IR2=EMISAN*SIG*A2*F32*TKOBJ**4
       IR3=EMISAN*SIG*A2*F42*TKSKY**4
       IR4=EMISAN*SIG*A2*F52*((TKSKY+TKSUB)/2.)**4
       IR5=EMISAN*SIG*A2*F52*((TKSKY+TKSUB)/2.)**4
c      bottom of thorax IR
       IR6=EMISAN*1*ATOT*EMISSB*SIG*TKSUB**4
       QIRIN=IR1+IR2+IR3+IR4+IR5+IR6
      else
       TKOBJ=TSUBST+273.15
       FSKY=FATOSK-FATOBJ
       IF(FSKY.LT.0.000) THEN
        FSKY=0.0
       ENDIF
       QIROBJ=EMISAN*FATOBJ*ATOT*EMISSB*SIG*TKOBJ**4
       QIRSKY=EMISAN*FSKY*ATOT*EMISSK*SIG*TKSKY**4
       QIRSUB=EMISAN*FATOSB*ATOT*EMISSB*SIG*TKSUB**4
       QIRIN=QIRSKY+QIRSUB+QIROBJ
      endif

      RETURN
      END
