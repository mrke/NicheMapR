      subroutine wetair(db,wb,rh,dp,bp,e,esat,vd,rw,tvir,tvinc,denair,
     +                  cp,wtrpot)
      implicit none
C    COPYRIGHT 1997  WARREN P. PORTER,  ALL RIGHTS RESERVED

c  subroutine wetair calculates several properties of humid air.  this version

c  was taken from "properties of air: a manual for use in biophysical ecology
c  third edition.
      real tk,db,wb,wbd,wbsat,esat,vapprs,dltae,wtrpot,cp,denair,tvinc
      real tvir,rw,vd,e,bp,dp,rh
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
      cp = (1004.84 + (rw * 1846.40)) / (1.0 + rw)
      if (rh .le. 0.0) go to 500
      wtrpot = 4.615e+5 * tk * alog(rh / 100.0)
      go to 600
500   wtrpot = -999
600   return
      end
