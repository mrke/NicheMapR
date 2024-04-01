      SUBROUTINE SOLAR
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

C     COMPUTES SOLAR RADIATION ABSORBED, EITHER FOR BODY OR (IF BUTTERFLY) WINGS

      IMPLICIT NONE

      DOUBLE PRECISION A1,A2,A3,A4,A4B,A5,A6,ABSAN,ABSMAX,ABSMIN,ABSSB
      DOUBLE PRECISION AL,AMASS,ANDENS,AREF,ASIL,ASILN,ASILP,ATOT,BREF
      DOUBLE PRECISION CREF,DEPSUB,EMISAN,EMISSB,EMISSK,F12,F13,F14,F15
      DOUBLE PRECISION F16,F21,F23,F24,F25,F26,F31,F32,F41,F42,F51,F52
      DOUBLE PRECISION F61,FATOBJ,FATOSB,FATOSK,FLSHCOND,FLUID,G
      DOUBLE PRECISION PDIF,PHI,PHIMAX,PHIMIN,PI,PTCOND
      DOUBLE PRECISION PTCOND_ORIG,QCOND,QCONV,QIRIN,QIROUT,QMETAB,QNORM
      DOUBLE PRECISION QRESP,QSDIFF,QSDIR,QSEVAP,QSOBJ,QSOLAR,QSOLR
      DOUBLE PRECISION QSRSB,QSSKY,RELHUM,RHO1_3,SHADE,SIDEX,SIG,SPHEAT
      DOUBLE PRECISION SUBTK,TA,TOBJ,TQSOL,TRANS1,TSKY,TSUBST,TWING,VEL
      DOUBLE PRECISION WC,WQSOL,ZEN,ZENITH,EGGPTCOND,POT,TOTLEN,AV,AT
      
      INTEGER IHOUR,LIVE,MICRO,NM,WINGCALC,WINGMOD

      COMMON/FUN1/QSOLAR,QIRIN,QMETAB,QRESP,QSEVAP,QIROUT,QCONV,QCOND
      COMMON/FUN2/AMASS,RELHUM,ATOT,FATOSK,FATOSB,EMISAN,SIG,FLSHCOND
      COMMON/FUN3/AL,TA,VEL,PTCOND,SUBTK,DEPSUB,TSUBST,PTCOND_ORIG,
     & EGGPTCOND,POT
      COMMON/FUN5/WC,ZEN,PDIF,ABSSB,ABSAN,ASILN,FATOBJ,NM
      COMMON/FUN6/SPHEAT,ABSMAX,ABSMIN,LIVE
      COMMON/WCOND/TOTLEN,AV,AT
      COMMON/WDSUB1/ANDENS,ASILP,EMISSB,EMISSK,FLUID,G,IHOUR
      COMMON/WDSUB2/QSOLR,TOBJ,TSKY,MICRO
      COMMON/WINGFUN/RHO1_3,TRANS1,AREF,BREF,CREF,PHI,F21,F31,F41,F51
     &,SIDEX,WQSOL,PHIMIN,PHIMAX,TWING,F12,F32,F42,F52
     &,F61,TQSOL,A1,A2,A3,A4,A4B,A5,A6,F13,F14,F15,F16,F23,F24,F25,F26
     &,WINGCALC,WINGMOD
      COMMON/WSOLAR/ASIL,SHADE
      
      DATA PI/3.14159265/

C     CHECKING FOR SCATTERED SKYLIGHT ONLY WHEN SUN BELOW HORIZON
      ZENITH=ZEN*180./PI
C     DIRECT BEAM COMPONENT
      IF(ZENITH.LT.90.00)THEN
C      DIRECT BEAM (NORMAL TO THE DIRECT BEAM)
       IF(LIVE.NE.1)THEN ! DON'T MAKE IT ADJUST POSTURE, IT'S DEAD!
        QNORM=QSOLR
       ELSE
        QNORM=(QSOLR/COS(ZEN))
       ENDIF
       IF(QNORM.GT.1367.)THEN
C       MAKING SURE THAT LOW SUN ANGLES DON'T LEAD TO SOLAR VALUES
C       GREATER THAN THE SOLAR CONSTANT
        QNORM=1367.
       ENDIF
       IF(ZENITH.GE.90.)THEN
        QNORM=0.000
       ENDIF
       QSDIR=ABSAN*ASIL*(1.00-PDIF)*QNORM*(1.0-(SHADE/100.))
      ELSE
       QSDIR=0.00
       QNORM=0.00
      ENDIF

      IF(WINGMOD.EQ.2)THEN
       CALL WINGS(RHO1_3,ABSAN,TRANS1,QSOLR,AREF,BREF,CREF,PHI,
     & F21,F31,F41,F51,F61,F12,F32,F42,F52,SIDEX,WQSOL,TQSOL,A1,A2,
     & A3,A4,A4B,A5,A6,F13,F14,F15,F16,F23,F24,F25,F26,ASILP)
       QSRSB=ABSAN*1.*ATOT/2.*(1.0-ABSSB)*QSOLR*(1.0-(SHADE/100.))
       QSOLAR=QSRSB+TQSOL
      ELSE
C      DIFFUSE COMPONENTS (SKY AND SUBSTRATE)
       QSOBJ=ABSAN*FATOBJ*(ATOT-AT/2.)*PDIF*QNORM
       IF(WINGMOD.EQ.1)THEN
        IF(PHI.LT.90)THEN
         QSSKY=0.
        ELSE
         QSSKY=ABSAN*FATOSK*(ATOT-AT/2.)*PDIF*QNORM*(1.0-(SHADE/100.))
        ENDIF
       ELSE
        QSSKY=ABSAN*FATOSK*(ATOT-AT/2.)*PDIF*QNORM*(1.0-(SHADE/100.))
       ENDIF
       QSRSB=ABSAN*FATOSB*(ATOT-AV-AT/2.)*(1.0-ABSSB)*QSOLR*
     &  (1.0-(SHADE/100.))
       QSDIFF=QSSKY+QSRSB+QSOBJ
       IF(WINGMOD.EQ.1)THEN
        IF(PHI.LT.90)THEN
         QSOLAR=QSDIFF
        ELSE
         QSOLAR=QSDIR+QSDIFF
        ENDIF
       ELSE
        QSOLAR=QSDIR+QSDIFF
       ENDIF
      ENDIF

      RETURN
      END