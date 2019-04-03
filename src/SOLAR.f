      SUBROUTINE SOLAR
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

c     Computes solar radiation absorbed, either for body or (if butterfly) wings

      Implicit None

      double precision ABSSB,ABSAN,AL,AMASS,ANDENS,ASIL,ASILN,ASILP,ATOT
      double precision DEPSUB,EMISAN,EMISSB,EMISSK,FATOSB,FATOSK
      double precision FATOBJ,Flshcond,FLUID,G,PI,PTCOND,PCTDIF,Qsevap
      double precision QCOND,QCONV,QIRIN,QIROUT,QMETAB,Qnorm,Qresp
      double precision QSDIFF,QSDIR,QSOBJ,QSOLAR,QSOLR,QSRSB,QSSKY
      double precision RELHUM,Shade,SIG,SUBTK,ptcond_orig,sidex,WQSOL
      double precision TA,TOBJ,TSKY,TSUBST,VEL,WC,ZEN,Zenith
      double precision rho1_3,trans1,aref,bref,cref,phi,F21,f31,f41,f51
     &,TQSOL,phimin,phimax,twing,F12,F32,F42,F52,f61,f13,f14,f15,f16
      double precision A1,A2,A3,A4,A4b,A5,A6,f23,f24,f25,f26
      double precision SPHEAT,ABSMAX,ABSMIN,O2MAX,O2MIN

      INTEGER IHOUR,MICRO,NM,wingmod,wingcalc,LIVE

      COMMON/FUN1/QSOLAR,QIRIN,QMETAB,QRESP,QSEVAP,QIROUT,QCONV,QCOND
      COMMON/FUN2/AMASS,RELHUM,ATOT,FATOSK,FATOSB,EMISAN,SIG,Flshcond
      COMMON/FUN3/AL,TA,VEL,PTCOND,SUBTK,DEPSUB,TSUBST,ptcond_orig
      COMMON/FUN5/WC,ZEN,PCTDIF,ABSSB,ABSAN,ASILN,FATOBJ,NM
      COMMON/FUN6/SPHEAT,ABSMAX,ABSMIN,O2MAX,O2MIN,LIVE
      COMMON/WINGFUN/rho1_3,trans1,aref,bref,cref,phi,F21,f31,f41,f51
     &,sidex,WQSOL,phimin,phimax,twing,F12,F32,F42,F52
     &,f61,TQSOL,A1,A2,A3,A4,A4b,A5,A6,f13,f14,f15,f16,f23,f24,f25,f26
     &,wingcalc,wingmod
      COMMON/WSOLAR/ASIL,Shade
      COMMON/WDSUB1/ANDENS,ASILP,EMISSB,EMISSK,FLUID,G,IHOUR
      COMMON/WDSUB2/QSOLR,TOBJ,TSKY,MICRO
      Data PI/3.14159265/

C     Checking for scattered skylight only when sun below horizon
      Zenith=Zen*180./PI
C     DIRECT BEAM COMPONENT
      If(Zenith.lt.90.00)then
C      DIRECT BEAM (normal to the direct beam)
       if(live.eq.0)then ! don't make it adjust posture, it's dead!
        Qnorm=QSOLR
       else
        Qnorm=(QSOLR/COS(ZEN))
       endif
       if(qnorm.gt.1367.)then
c       making sure that low sun angles don't lead to solar values
c       greater than the solar constant
        qnorm=1367.
       endif
       if(zenith.ge.90.)then
        qnorm=0.000
       endif
       QSDIR=ABSAN*ASIL*(1.00-PCTDIF)*Qnorm*(1.0-(shade/100.))
      else
       Qsdir=0.00
       qnorm=0.00
      endif

      if(wingmod.eq.2)then
       call wings(rho1_3,absan,trans1,QSOLR,aref,bref,cref,phi,
     & F21,f31,f41,f51,f61,f12,f32,f42,f52,sidex,WQSOL,TQSOL,A1,A2,
     & A3,A4,A4b,A5,A6,f13,f14,f15,f16,f23,f24,f25,f26,asilp)
       QSRSB=ABSAN*1*ATOT/2*(1.0-ABSSB)*QSOLR*(1.0-(shade/100.))
       QSOLAR=QSRSB+TQSOL
      else
C      DIFFUSE COMPONENTS (SKY AND SUBSTRATE)
       QSOBJ=ABSAN*FATOBJ*ATOT*PCTDIF*Qnorm
       if(wingmod.eq.1)then
        if(phi.lt.90)then
         QSSKY=0.
        else
         QSSKY=ABSAN*FATOSK*ATOT*PCTDIF*Qnorm*(1.0-(shade/100.))
        endif
       else
        QSSKY=ABSAN*FATOSK*ATOT*PCTDIF*Qnorm*(1.0-(shade/100.))
       endif
       QSRSB=ABSAN*FATOSB*ATOT*(1.0-ABSSB)*QSOLR*(1.0-(shade/100.))
       QSDIFF=QSSKY+QSRSB+QSOBJ
       if(wingmod.eq.1)then
        if(phi.lt.90)then
         QSOLAR=QSDIFF
        else
         QSOLAR=QSDIR+QSDIFF
        endif
       else
        QSOLAR=QSDIR+QSDIFF
       endif
      endif

      RETURN
      END
