      Function FUNC(xi)
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

C     Function Xi calculates Xi as part of the adjacent rectangles
C     configuration/numerical integration calculations.


      Implicit none

      double precision angle,ans,degangl,fa,fb,fb1,fb2,fc,fc1,fc2,func
      double precision pi,x,xi
      
      Common/intgl/x,angle
      data pi/3.14159265/

c     x = lengths ratio, a/b
c     xi = all (sequential) values between integration limits
c     sin**2(angle) = 0.5*(1.-cos(2*angle)
c     debugging info: angle in degrees

      degangl = angle*180./pi
      fa = sqrt(1.0 + xi**2*(0.5*(1.-cos(2*angle))))
      fb1 = (x-xi*cos(angle))
      fb2 = sqrt(1.0+xi**2* (0.5*(1.-cos(2*angle))))
      fb = atan(fb1/fb2)
      fc1 = xi*cos(angle)
      fc2 = sqrt(1.0+xi**2* (0.5*(1.-cos(2*angle))))
      fc=atan(fc1/fc2)
      ans = fa * (fb + fc)
      func = ans
      Return
      End