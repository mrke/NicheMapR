      subroutine wetair(db,wb,rh,dp,bp,e,esat,vd,rw,tvir,tvinc,denair,
     +                  cp,wtrpot)
      implicit none

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

c     Subroutine wetair calculates several properties of humid air.
c     This version was taken from
c     Tracy, C. R., W. R. Welch, B. Pinshow, M. R. Kearney, and W. P.
c     Porter. 2016. Properties of air: A manual for use in biophysical ecology.
c     5th edition. The University of Wisconsin, Madison.
c     which is available as a vignette in the NicheMapR package

      double precision tk,db,wb,wbd,wbsat,esat,vapprs,dltae,wtrpot,cp
      double precision tvir,rw,vd,e,bp,dp,rh,denair,tvinc
      tk  = db + 273.15
      esat = vapprs(db)
      if (dp .lt. 999.0) go to 100
      if (rh .gt. -1.0) go to 200
      wbd = db - wb
      wbsat = vapprs(wb)
      dltae = 0.000660 * (1.0 + 0.00115 * wb) * bp * wbd
      e = wbsat - dltae
      go to 300
100   e = vapprs(dp)
      go to 300
200   e = esat * rh / 100
      go to 400
300   rh = (e / esat) * 100
400   rw = ((0.62197 * 1.0053 * e) / (bp - 1.0053 * e))
      vd = e * 0.018016 / (0.998 * 8.31434 * tk)
      tvir = tk * ((1.0 + rw / (18.016 / 28.966)) / (1.0 + rw))
      tvinc = tvir - tk
      denair = 0.0034838 * bp / (0.999 * tvir)
      cp = (1004.84 + (rw * 1864.40)) / (1.0 + rw)
      if (rh .le. 0.0) go to 500
      wtrpot = 4.615D+5 * tk * dlog(rh / 100.0)
      go to 600
500   wtrpot = -999
600   return
      end
