      SUBROUTINE RYTREC

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

C     THIS SUBROUTINE COMPUTES THE CONFIGURATION FACTOR BETWEEN TWO
C     FINITE RECTANGLES WITH A COMMON SIDE AT RIGHT ANGLES TO EACH OTHER.
C     EQUATIONS FROM HOTTEL, 1931;
C     HOTTEL, H.C., 1931, "RADIANT HEAT TRANSMISSION BETWEEN SURFACES
C     SEPARATED BY NON-ABSORBING MEDIA," TRANS. ASME, VOL. 53,
C     FSP-53-196, PP. 265-273.
C     AND FROM HAMILTON AND MORGAN
C     HAMILTON, D.C. AND MORGAN, W.R., 1952, "RADIANT-INTERCHANGE
C     CONFIGURATION FACTORS," NASA TN 2836.

C     ONE OF THE CLASSIC COMPILATIONS OF CONFIGURATION FACTORS. HAS A FEW TYPOGRAPHICAL
C     ERRORS [SEE, E.G., FEINGOLD (1966), FEINGOLD AND GUPTA (1970), AND BYRD.] CATALOGS
C     TWELVE DIFFERENT DIFFERENTIAL AREA TO FINITE AREA FACTORS, FIVE DIFFERENTIAL STRIP TO FINITE
C     AREA FACTORS, AND ELEVEN FINITE AREA TO FINITE AREA FACTORS. SOME OF THE FACTORS ARE
C     GENERATED BY CONFIGURATION FACTOR ALGEBRA FROM A SMALLER SET OF CALCULATED OR
C     DERIVED FACTORS. THIS IS A PIONEERING WORK IN CATALOGUING USEFUL INFORMATION.
C     (FROM HOWELL'S WEB PAGE)
C     HTTP://WWW.ME.UTEXAS.EDU/~HOWELL/

C     THE WIDTH (FROM THE COMMON EDGE TO THE FAR EDGE) OF SURFACE 2 IS "A"
C     THE COMMON EDGE IS LENGTH B
C     THE WIDTH OF SURFACE 1 IS "C"
C     F12 = FCA

      IMPLICIT NONE
      DOUBLE PRECISION A,B,C,F,GRP1,GRP2,GRP3,GRP4,GRP5,H,H2,HPW,HSQ
      DOUBLE PRECISION PI,W,W2,WSQ,X,Y
      
      COMMON/RECTNGL/A,B,C,F

      PI=3.14159265

C     WIDTH OF BODY = A
C     LENGTH OF BODY & WING (PARALLEL TO LONG AXIS OF BODY) = B
C     WIDTH OF WING (FROM JUNCTION W BODY TO WING TIP) = C
      X=A/B
      Y=C/B
C     USING THE NOTATION FROM HOWELL
      H=X
      W=Y
C     DEFINING SOME COMMON GROUPS
      HSQ=H**2.
      WSQ=W**2.
      H2=1.+HSQ
      W2=1.+WSQ
      HPW=HSQ+WSQ

C     DEFINING SUBGROUPS OF THE GOVERNING EQUATION
      GRP1=1./(W*PI)
      GRP2=W*ATAN(1./W)+ H*ATAN(1./H)-SQRT(HPW)*ATAN(SQRT(1./HPW))
      GRP3=(W2*H2)/(W2+H**2.)
      GRP4=(W**2. *(W2+H**2.))/(W2*HPW)
      GRP5=(H**2. *(H2+W**2.))/(H2*HPW)

C     GOVERNING EQUATION
      F=GRP1*(GRP2 + 0.25*LOG(GRP3*(GRP4**WSQ)*(GRP5**HSQ)))
      RETURN
      END