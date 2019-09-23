      SUBROUTINE SOLRAD

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

C    25 June 1998 version incorporates skylight before sunrise and after sunset.  Spectral
C    qualities are not included, but total energy is based on data from
C    V.V. Sharonov. 1948. Illumination of a horizontal plane at twilight and at night.
C    Dokl. Akad. Naukk SSSR. 42: 296-300. 1944;
C    Trudy ykubileinoi nauch. Sess. leningr. gos. Univ. (1948) 47.
C    in the book Twilight by G.V. Rozenberg. 1966. (p.19) Plenum Press. 358 p.

C     PROGRAM SOLRAD GENERATES TERRESTRIAL SOLAR RADIATION SPECTRAL AND
C     INTEGRATED ESTIMATES. THIS PROGRAM FOLLOWS THE DEVELOPMENT IN THE
C     PAPER 'COMPUTING CLEAR DAY SOLAR RADIATION SPECTRA FOR THE TERRESTRIAL
C     ECOLOGICAL ENVIRONMENT' , BY E.C. MCCULLOUGH AND W.P. PORTER
C     TO BE PUBLISHED IN ECOLOGY (1971). THE USER ENTERS APPROPRIATE
C     VALUES IN THE TWO BRACKETED SECTIONS BELOW.   IF USER WISHES TO ENTER
C     A VALUE FOR TOTAL ATMOSPHERIC OZONE, THIS IS DONE AFTER STATEMENT 150.

C     ** 24 DECEMBER '88 VERSION FOR THE IBM (WPP) USING MENUS

C     USUALLY YOU RUN SOLRAD1 FIRST, TO GET SOLAR RADIATION AND
C     ZENITH ANGLES, THEN SINE TO COMPUTE AIR TEMPERATURES,
C     SINCE SOLRAD COMPUTES TIME OF
C     SUNRISE NEEDED IN SINE FOR TMIN (1 HOUR BEFORE SUNRISE).

C     THE FOLLOWING CONTAIN THE INTEGRATED VALUES FOR THE SPECTRAL ESTIMATE
C     DESIGNATED, THAT IS, INTEGRATED OVER WAVELENGTH FROM 290 NM TO THE NTH
C     WAVELENGTH ILAM(N)-------

C     DRINT(N)     DIRECT RADIATION COMPONENT
C     DRRINT(N)    DIRECT RAYLEIGH RADIATION COMPONENT
C     SRINT(N)     SCATTERED RADIATION COMPONENT
C     GRINT(N)     GLOBAL RADIATION COMPONENT (DIRECT + SCATTERED)
C     EPINT(N)     ERYTHEMA PRODUCT
C     REPINT(N)    RECIPROCAL OF EPINT
      use commondat

      implicit none
      DOUBLE PRECISION TLAM1,GAMR(101),GAML(101),SBAR

      double precision a,aalat,abc,ailam,airms,ala,ALAT,allat,altdeg
      double precision altfct,altrad,alon,ALONC,amr,AMULT,AZMUTH,ALTT
      double precision ar2,AZSLP,AZSUN,b,bothalf
      double precision c,cazsun,cdrnt,CMH2O
      double precision cdrrnt,cgrnt,csrnt,cz,czsl,drrlam,dzsl
      double precision d,d0,dazsun,ddec,dec,delt,drint,drlam,drrint
      double precision ec,eom,epint,epl,eplam,er,erlam,elog
      double precision fd,fdav,fdq,fdqdav,FRACT,FRACTR,FRACTS,grint
      double precision global,gpl,grlam,h,ha,had,hemis,hh,HHINT,HRINT
      double precision MAXSHD,HSINT,hzen
      double precision ocgrnt,ODNEW,ODRNT,ODSTOR,OSTORE,OSNEW
      double precision oz,OZNEW,ozone,ozsl,ozstor,ozz
      double precision part1,part2,pi,PRESS,PTWET,PUNSH
      double precision REFL,refr,repint
      double precision s,SABNEW,se,SHAYD,sinde,slam
      double precision SLOPE,SLOPEL,skylum,srint
      double precision srlam,szen,ta,Tannul,TAZSUN,td,tdtl
      double precision testcz,testha,testz,testzz
      double precision TIMCOR,TIMAX,TIMAXS,time,tlam3,tlam4
      double precision TIMIN,TIMINS,timsr,timss,timtmx,tlat,tlam,tlam2
      double precision TMAX,TMIN,TSTOR,tw,tdd,TI
      double precision TNEW,to,tophalf,tr,ts,tsn,tsolmx,TSRHR,TSNHR
      double precision vmax,vmin,w,z,zsl,ZSLNEW,zslsto,zz
      double precision soilprop,rainfall
      double precision tannul2,ahoriz,hori,azi,moist
      double precision altfct1,altfct2,altfct3,altfct4,alt,tai
      double precision itair,icld,iwind,irelhum,tmin2,tmax2,sle,err
      double precision DEPS,surflux,ep,densfun,snowcond,intercept
      double precision snownode,minsnow,maxsnode1,snode,daysincesnow,
     &lastday,undercatch,rainmeltf

      double precision snowdens,snowmelt,snowtemp,cursnow,snowage,
     & prevden,grasshade

      INTEGER I,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96
      integer ia,IALT,icorr,iday,idom,ielam,IDA,IDAYST
      integer iedum,IEND,iep,ILAM,intcz,INTIME
      integer IPINT,isos,ISTART,IT,IUV,IVAR,jct,JJ,JTEST,julnum
      integer K,lamax,llat,M,malt,mm,mmm,mom,DOY
      integer n,NCT,ND,nft,nmax,NOSCAT,NUMRUN,nzct
      integer microdaily,cnt,ILOCT
      INTEGER I97,I98,I99,methour,lamb,I100,I101

      CHARACTER(12) FNAME
      CHARACTER(80) LABL1,LABL2,LABL3
      CHARACTER(1) ANS1,ANS2,ANS3,ANS4,ANS5,ANS6,ANS7,ANS10
      character(1) ANS11,ANS14,ANS15,ANS16,ANS17,ANS18
      CHARACTER(1) SINE,SNSLOP,solout

      DIMENSION DRINT(111),SRINT(111),GRINT(111)
      DIMENSION DRRLAM(111),DRRINT(111)
      DIMENSION SLAM(111),TR(111),TA(111),TO(111),TW(111)
      DIMENSION EOM(13),snownode(10),snode(10)
      DIMENSION OZ(19,12),DRLAM(111),ILAM(111)
      DIMENSION AILAM(111),SRLAM(111),GRLAM(111)
      DIMENSION EPLAM(7),EPINT(7),REPINT(7)
      DIMENSION ALTFCT(21,4)
      DIMENSION S(11),ER(7),FD(11,19),FDQ(11,19)
      DIMENSION A(4),EPL(12),IELAM(12),GPL(11),ERLAM(11)
      DIMENSION TIME(25),OCGRNT(25),ODNEW(25),ODRNT(25),ODSTOR(25)
      DIMENSION OZSTOR(25),OZSL(25),OZZ(25),TSTOR(25),ZSLSTO(25)
      DIMENSION OSNEW(25),OSTORE(25),OZNEW(25),TNEW(25),ZSLNEW(25)
      DIMENSION TIMINS(4),TIMAXS(4)
      DIMENSION hori(24),azi(24),tai(111)
      DIMENSION soilprop(10,5),moist(10)
      DIMENSION DEPS(21),densfun(4)

      COMMON/LABEL/LABL1,LABL2,LABL3,FNAME,SINE,ANS14,SNSLOP
      COMMON/LABELS/ANS16,ANS17,ANS18
      COMMON/DAYS/Tannul
      COMMON/DAYSS/TIMINS,TIMAXS
      COMMON/WIOCONS/PUNSH,ALAT,AMULT,PRESS,CMH2O,REFL,ALONC,TIMCOR,
     * AZMUTH,SLOPE,TSNHR,TSRHR,Hemis
      COMMON/WIOCONS2/IPINT,NOSCAT,IUV,IALT,IDAYST,IDA,IEP,ISTART,IEND
      COMMON/WSINE/TIMSR,TIMSS,TIMTMX,TMIN,TMAX,TMIN2,TMAX2
      COMMON/NDAY/ND
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101
      common/solropt/solout
      Common/wfun/airms,cz
      COMMON/DAYJUL/JULNUM,DOY
      COMMON/GROUND/SHAYD,ALTT,MAXSHD,SABNEW,PTWET,rainfall
      COMMON/SNOWPRED/snowtemp,snowdens,snowmelt,snownode,minsnow
     &,maxsnode1,snode,cursnow,daysincesnow,lastday,undercatch,rainmeltf
     &,densfun,snowcond,intercept,snowage,prevden,grasshade
      COMMON/WICHDAY/NUMRUN
      COMMON/DAILY/microdaily
      common/deep/tannul2
      common/horizon/hori,azi
      common/atten/tai,ec
      common/init/itair,icld,iwind,irelhum,iday
      COMMON/NICHEMAPRIO/SLE,ERR,soilprop,surflux
      common/moistcom/moist,ep
      COMMON/dataky/DEPS,CNT
      COMMON/TABLE/TI(211),TDD(211)
      COMMON/TABLE2/ILOCT(21)
      COMMON/LAMBDA/DRLAM,DRRLAM,SRLAM,LAMB

      DATA GRINT/111*0./
      DATA DRRINT/111*0./
      DATA DRINT/111*0./
      DATA SRINT/111*0./
      DATA EPINT/7*0./
      DATA REPINT/7*0./
      DATA AILAM/111*0./
      DATA GRLAM/111*0./
      DATA OZSTOR/25*90./
      DATA TIME/25*0./
      DATA OCGRNT/25*0./
      DATA ODRNT/25*0./
      DATA OZZ/25*90./
      DATA OZSL/25*90./
      DATA OZNEW/25*90./
      DATA ZSLNEW/25*90./
      DATA ZSLSTO/25*90./
      DATA TSTOR/25*0./
      DATA OSTORE/25*0./
      DATA A/4*0./
      DATA EPL/12*0./
      DATA IELAM/12*0/
      DATA GPL/11*0./
      DATA ERLAM/11*0./
      DATA EPLAM/7*0./


      DATA ILAM/290,295,300,305,310,315,320,330,340,350,360,370,380,390,
     &400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,72
     &0,740,760,780,800,820,840,860,880,900,920,940,960,980,1000,1020,10
     &80,1100,1120,1140,1160,1180,1200,1220,1240,1260,1280,1300,1320,138
     &0,1400,1420,1440,1460,1480,1500,1540,1580,1600,1620,1640,1660,1700
     &,1720,1780,1800,1860,1900,1950,2000,2020,2050,2100,2120,2150,2200,
     &2260,2300,2320,2350,2380,2400,2420,2450,2490,2500,2600,2700,2800,2
     &900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000/

      DATA SLAM/48.2,58.4,51.4,60.2,68.6,75.7,81.9,103.7,105,107.4,105.5
     & ,117.3,111.7,109.9,143.3,175.8,182.3,208,208.5,194.6,183.3,178.3,
     &169.5,170.5,164.6,157.6,151.7,146.8,141.8,136.9,131.4,126,120,115,
     &110.7,105,100,95,91,88,85,83,80,77,75,70,61,59,56,54,51,49,48,45,4
     &2,41,40,39,38,34,33,32,31,30,29,28,26,25,24,24,23,22,21,19,16,15,1
     &2,11,10.7,10.3,10,9.7,9,8.8,8.5,7.9,7.4,6.8,6.7,6.6,6.5,6.4,6.2,5.
     &9,5.5,5.4,4.8,4.3,3.9,3.5,3.1,2.6,2.3,1.9,1.7,1.5,1.4,1.2,1.1,1,1/

      DATA TR/1.41,1.31,1.23,1.14,1.05,0.99,0.92,0.81,0.72,0.63,0.56,0.5
     &,0.45,0.4,0.36,0.3,0.25,0.2,0.17,0.15,0.12,0.1,0.09,0.08,0.07,0.06
     &,0.06,0.05,0.04,0.04,0.03,0.03,0.03,0.02,0.02,0.02,0.02,0.01,0.01,
     &0.01,0.01,0.01,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
     &,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
     &,0,0,0,0,0,0,0,0,0,0/

      DATA TA/0.42,0.415,0.412,0.408,0.404,0.4,0.395,0.388,0.379,0.379,0
     &.379,0.375,0.365,0.345,0.314,0.3,0.288,0.28,0.273,0.264,0.258,0.25
     &3,0.248,0.243,0.236,0.232,0.227,0.223,0.217,0.213,0.21,0.208,0.205
     &,0.202,0.201,0.198,0.195,0.193,0.191,0.19,0.188,0.186,0.184,0.183,
     &0.182,0.181,0.178,0.177,0.176,0.175,0.175,0.174,0.173,0.172,0.171,
     &0.17,0.169,0.168,0.167,0.164,0.163,0.163,0.162,0.161,0.161,0.16,0.
     &159,0.157,0.156,0.156,0.155,0.154,0.153,0.152,0.15,0.149,0.146,0.1
     &45,0.142,0.14,0.139,0.137,0.135,0.135,0.133,0.132,0.131,0.13,0.13,
     &0.129,0.129,0.128,0.128,0.128,0.127,0.127,0.126,0.125,0.124,0.123,
     &0.121,0.118,0.117,0.115,0.113,0.11,0.108,0.107,0.105,0.103,0.1/
      DATA TO/11.5,6.3,3.2,1.62,0.83,0.44,0.26,0.03,0.02,0.01,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/

      DATA TW/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &0,0,0,0,0,0,0,0,0,0.123,0.117,0.1,0.23,0.174,0.058,0,0.024,0.027,0
     &.036,0.215,0.25,0.136,0.058,0.047,0.036,0.042,0.098,0.044,0,0.038,
     &0.83,0,0,0.38,0.289,0.258,0.173,0.008,0,0,0,0,0,0,0,0,0.57,0.76,0,
     &0.185,0.291,0.178,0.196,0.112,0.075,0.074,0.07,0.007,0,0,0,0.086,0
     &.122,0.132,0.14,0.207,0.259,0,0,0,0.549,0.297,0.462,0.52,0.374,0.2
     &22,0.614,0.058,0.038,0.03,0.04,0.16/

      DATA OZ/0.31,0.31,0.32,0.32,0.31,0.3,0.27,0.24,0.23,0.22,0.23,0.24
     &,0.27,0.3,0.32,0.33,0.34,0.34,0.33,0.3,0.31,0.31,0.31,0.3,0.29,0.2
     &8,0.25,0.24,0.22,0.24,0.26,0.28,0.32,0.36,0.39,0.4,0.4,0.39,0.3,0.
     &31,0.31,0.3,0.29,0.28,0.26,0.24,0.24,0.23,0.24,0.26,0.29,0.33,0.38
     &,0.42,0.45,0.46,0.46,0.27,0.28,0.29,0.3,0.3,0.29,0.27,0.25,0.24,0.
     &23,0.25,0.27,0.3,0.34,0.38,0.4,0.42,0.43,0.42,0.34,0.35,0.34,0.33,
     &0.32,0.31,0.28,0.25,0.24,0.24,0.26,0.28,0.3,0.34,0.37,0.39,0.4,0.4
     &,0.39,0.38,0.4,0.39,0.38,0.36,0.33,0.28,0.25,0.24,0.24,0.25,0.27,0
     &.3,0.33,0.35,0.36,0.36,0.36,0.34,0.43,0.44,0.43,0.41,0.39,0.35,0.2
     &9,0.25,0.24,0.24,0.25,0.26,0.29,0.31,0.33,0.34,0.34,0.33,0.32,0.45
     &,0.46,0.45,0.42,0.4,0.37,0.31,0.26,0.24,0.23,0.24,0.26,0.28,0.3,0.
     &31,0.32,0.31,0.3,0.3,0.41,0.42,0.43,0.42,0.4,0.38,0.32,0.26,0.24,0
     &.23,0.24,0.26,0.27,0.28,0.3,0.3,0.29,0.28,0.27,0.37,0.38,0.4,0.4,0
     &.39,0.37,0.32,0.26,0.24,0.22,0.23,0.25,0.26,0.27,0.28,0.28,0.28,0.
     &27,0.26,0.34,0.36,0.38,0.39,0.37,0.34,0.29,0.26,0.24,0.22,0.23,0.2
     &5,0.26,0.28,0.29,0.3,0.29,0.29,0.28,0.31,0.32,0.34,0.35,0.35,0.32,
     &0.29,0.25,0.23,0.22,0.23,0.25,0.27,0.29,0.3,0.31,0.31,0.31,0.3/

      DATA ALTFCT/1,0.887,0.784,0.692,0.609,0.534,0.467,0.406,0.353,0.30
     &4,0.262,0.225,0.192,0.164,0.141,0.12,0.103,0.088,0.075,0.064,0.064
     &,1,0.546,0.353,0.26,0.223,0.199,0.181,0.168,0.155,0.142,0.129,0.11
     &6,0.103,0.092,0.081,0.07,0.059,0.048,0.039,0.031,0.031,1,0.999,0.9
     &81,0.973,0.966,0.959,0.952,0.946,0.939,0.932,0.923,0.911,0.895,0.8
     &73,0.846,0.817,0.787,0.756,0.721,0.682,0.682,1,1,1,1,1,1,1,1,1,1,1
     &,1,1,1,1,1,1,1,1,1,1/

      DATA FD/8.00E-05,6.50E-05,4.00E-05,2.30E-05,1.00E-05,4.50E-06,1.00
     &E-06,1.00E-07,5.50E-09,1.00E-09,3.50E-10,1.60E-10,1.00E-10,1.00E-1
     &0,1.00E-10,1.00E-10,1.00E-10,1.00E-10,1.00E-10,1.00E-03,9.50E-04,9
     &.00E-04,8.00E-04,7.00E-04,6.00E-04,4.50E-04,3.00E-04,1.70E-04,8.00
     &E-05,3.30E-05,1.80E-05,1.00E-05,8.00E-06,6.30E-06,5.00E-06,4.10E-0
     &6,3.50E-06,3.00E-06,3.50E-02,3.30E-02,3.10E-02,2.90E-02,2.50E-02,2
     &.20E-02,1.75E-02,1.30E-02,7.50E-03,4.50E-03,2.50E-03,1.30E-03,5.60
     &E-04,2.70E-04,1.40E-04,7.10E-05,4.20E-05,2.60E-05,1.70E-05,1.55E-0
     &1,1.40E-01,1.40E-01,1.30E-01,1.20E-01,1.10E-01,1.00E-01,9.00E-02,8
     &.20E-02,7.00E-02,5.20E-02,3.50E-02,2.20E-02,1.00E-02,4.00E-03,1.20
     &E-03,4.50E-04,1.90E-04,8.00E-05,3.70E-01,3.70E-01,3.60E-01,3.50E-0
     &1,3.30E-01,3.10E-01,2.90E-01,2.60E-01,2.30E-01,2.00E-01,1.70E-01,1
     &.50E-01,1.00E-01,7.00E-02,3.80E-02,1.60E-02,5.00E-03,1.20E-03,3.00
     &E-04,5.50E-01,5.50E-01,5.50E-01,5.30E-01,5.10E-01,5.00E-01,4.70E-0
     &1,4.40E-01,4.00E-01,3.60E-01,3.10E-01,2.70E-01,2.25E-01,1.80E-01,1
     &.20E-01,6.00E-02,2.50E-02,4.50E-03,7.00E-04,6.50E-01,6.50E-01,6.50
     &E-01,6.50E-01,6.20E-01,6.00E-01,5.70E-01,5.50E-01,5.00E-01,4.50E-0
     &1,4.20E-01,3.80E-01,3.25E-01,2.70E-01,1.90E-01,1.15E-01,5.00E-02,1
     &.10E-02,1.20E-03,7.88E-01,7.80E-01,8.00E-01,8.00E-01,8.00E-01,7.60
     &E-01,7.35E-01,7.10E-01,7.00E-01,6.50E-01,6.00E-01,5.50E-01,5.14E-0
     &1,4.50E-01,3.50E-01,2.68E-01,1.50E-01,6.27E-02,1.20E-02,7.48E-01,7
     &.40E-01,7.40E-01,7.30E-01,7.20E-01,7.10E-01,7.04E-01,6.90E-01,6.70
     &E-01,6.20E-01,5.70E-01,5.30E-01,5.16E-01,4.80E-01,3.90E-01,2.90E-0
     &1,1.70E-01,7.62E-02,3.00E-02,7.00E-01,7.00E-01,7.00E-01,6.90E-01,6
     &.80E-01,6.80E-01,6.60E-01,6.50E-01,6.30E-01,6.00E-01,5.60E-01,5.21
     &E-01,5.00E-01,4.70E-01,3.90E-01,3.00E-01,1.85E-01,9.00E-02,2.60E-0
     &2,6.51E-01,6.50E-01,6.50E-01,6.40E-01,6.30E-01,6.25E-01,6.22E-01,6
     &.00E-01,5.90E-01,5.70E-01,5.50E-01,5.20E-01,4.89E-01,4.60E-01,3.90
     &E-01,3.08E-01,2.00E-01,9.55E-02,2.20E-02/

      DATA FDQ/8.00E-06,7.00E-06,5.20E-06,3.50E-06,1.70E-06,5.50E-07,1.0
     &0E-07,2.50E-08,6.00E-09,1.50E-09,3.00E-10,6.00E-11,1.00E-11,1.00E-
     &11,1.00E-11,1.00E-11,1.00E-11,1.00E-11,1.00E-11,6.10E-04,6.00E-04,
     &5.50E-04,4.50E-04,3.40E-04,2.30E-04,1.20E-04,5.50E-05,2.60E-05,1.2
     &0E-05,6.00E-06,3.50E-06,2.00E-06,1.50E-06,1.00E-06,6.00E-07,4.00E-
     &07,2.30E-07,1.00E-07,2.40E-02,2.30E-02,2.20E-02,2.10E-02,1.80E-02,
     &1.50E-02,1.20E-02,7.50E-03,4.80E-03,2.50E-03,1.20E-03,5.00E-04,2.5
     &0E-04,1.00E-04,4.50E-05,2.00E-05,1.00E-05,5.50E-06,2.50E-06,1.30E-
     &01,1.20E-01,1.10E-01,1.00E-01,9.10E-02,8.50E-02,7.20E-02,6.50E-02,
     &5.20E-02,4.00E-02,2.80E-02,1.70E-02,1.00E-02,3.00E-03,1.00E-03,4.0
     &0E-04,1.70E-04,6.70E-05,2.50E-05,3.40E-01,3.30E-01,3.20E-01,3.00E-
     &01,2.90E-01,2.70E-01,2.50E-01,2.00E-01,1.70E-01,1.50E-01,1.20E-01,
     &8.00E-02,5.80E-02,3.00E-02,1.30E-02,6.20E-03,2.00E-03,6.00E-04,1.9
     &0E-04,5.40E-01,5.30E-01,5.10E-01,5.00E-01,4.80E-01,4.50E-01,4.20E-
     &01,3.70E-01,3.20E-01,2.70E-01,2.20E-01,1.70E-01,1.20E-01,8.00E-02,
     &5.00E-02,2.30E-02,8.00E-03,4.70E-03,9.00E-04,6.50E-01,6.40E-01,6.2
     &0E-01,6.10E-01,6.00E-01,5.60E-01,5.10E-01,4.70E-01,4.20E-01,3.60E-
     &01,3.00E-01,2.50E-01,1.80E-01,1.30E-01,8.00E-02,5.00E-02,2.00E-02,
     &4.60E-03,1.50E-03,8.20E-01,8.10E-01,8.10E-01,8.00E-01,7.90E-01,7.5
     &0E-01,7.00E-01,6.50E-01,6.00E-01,5.40E-01,4.60E-01,3.80E-01,3.20E-
     &01,2.50E-01,1.80E-01,1.20E-01,6.00E-02,2.50E-02,8.00E-03,8.60E-01,
     &8.40E-01,8.20E-01,8.00E-01,7.50E-01,7.00E-01,6.80E-01,6.30E-01,5.8
     &0E-01,5.00E-01,4.40E-01,3.70E-01,3.20E-01,2.40E-01,1.70E-01,1.20E-
     &01,6.00E-02,2.70E-02,1.00E-02,8.00E-01,7.90E-01,7.70E-01,7.60E-01,
     &7.30E-01,6.85E-01,6.40E-01,5.90E-01,5.40E-01,4.75E-01,4.20E-01,3.6
     &5E-01,3.20E-01,2.50E-01,1.80E-01,1.30E-01,7.00E-02,2.90E-02,1.10E-
     &02,7.50E-01,7.40E-01,7.30E-01,7.20E-01,7.00E-01,6.70E-01,6.10E-01,
     &5.50E-01,5.00E-01,4.50E-01,4.00E-01,3.60E-01,3.20E-01,2.60E-01,1.9
     &0E-01,1.40E-01,8.00E-02,3.10E-02,1.20E-02/

      DATA ER/6.19E+00,6.86E+00,1.16E+01,2.51E+01,2.24E+02,5.60E+02,1.16
     &E+03/

      DATA S/0.2,0.255,0.315,0.365,0.394,0.405,0.405,0.395,0.37,0.343,0.
     &32/

      DATA EOM/0,31,59,90,120,151,181,212,243,273,304,334,365/

      NUMRUN=1
      DRRLAM(:)=GRINT
      DRLAM(:)=GRINT
      SRLAM(:)=GRINT
C    DEFINING CONSTANTS THAT USERS USED TO SPECIFY
      NMAX = 111
      IEP = 1
      AMR = 25.
C    DEFINING DEFAULT PARAMETERS FROM DELETED USER INTERFACE
C    TOTAL SOLAR INTEGRATED
      ANS1 = 'T'
C    GLOBAL SOLAR OUT (BEAM & DIFFUSE OUT SEPARATELY, NOW)
      ANS2 = 'G'
C    ROUTINE SCATTERING CALCULATION, NO USE OF SUB. GAMMA
      ANS3 = 'R'
C    NO PLOT FORMAT (USED FOR SPECTRAL SOLAR OUTPUT)
      ANS4 = 'N'
C    NORTHERN (Y) OR SOUTHERN (N) HEMISPHERE?
      ANS5 = 'N'
C    USE STANDARD ATM. PRESSURE FOR ELEVATION
      ANS6 = 'Y'
C    12 MIDDAYS (JULIAN) OF EACH MONTH SUPPLIED
      ANS7 = 'S'
C    HORIZONTAL SLOPE?
c     ANS9 = 'Y'
C    CALCULATE HOURLY TEMPERATURES FROM TMAX, TMIN
      ANS10 = 'Y'
C    DATAFILE ENTRY FOR Tmax, Tmin Tairs
      ANS11 = 'D'
C    MICROMET OUTPUT FORMAT
      ANS14 = 'M'
C    DON'T CHOOSE START & END HOURS; 0 - 24 HOURS USED AS DEFAULT
      ANS15 = 'N'
C    SUN ON HORIZONTAL SURFACE, 'H' or sloping surface, 'S'
      SNSLOP = 'H'
C    HOUR(S)  AFTER SUNRISE WHEN TAIR = TAIR, MIN
      TSRHR = TIMINS(1)!0.0
C    HOUR(S) AFTER SOLAR NOON WHEN TAIR = TAIR, MAX
      TSNHR = TIMAXS(1)!1.0

      TSN=0
      hh=0
      altdeg=0
      ahoriz=0
      CDRNT = 0.D+0
      M=0
      DZSL=0.D0
      dazsun=0.D0
      
C     FORMAT STATEMENTS
  991 FORMAT(1X,2F4.0,1I6,1E12.4)

C     INITIALIZING Output ARRAYS
      TIME(1)=0.D0
      TIME(25)=1440.
      TSTOR(1)=0.D0
      OCGRNT(1)=0.D0
      OCGRNT(25)=0.D0
      OSTORE(1)=0.D0
      OZZ(1)=90.
      OZZ(25)=90.
      OZSTOR(1)=90.
      OZSTOR(25)=90.
      OZSL(1)=90.
      OZSL(25)=90.
      ZSLSTO(1)=90.
      ZSLSTO(25)=90.

c    addition by Kearney to use GADS-calculated TA
      do 911 i=1,111
      TA(I)=TAI(I)
911   continue

C     THE ARRAY ALTFCT(I,J) IS THE RATIO OF THE TOTAL AMOUNT OF A GIVEN
C     ATMOSPHERIC CONSTITUENT (INDEX J) ABOVE THE ALTITUDE OF INTEREST
C     (INDEX I) TO THAT ABOVE SEA LEVEL. FOR MOLECULAR (J=1), AEROSOL (J=2), AND
C     OZONE(J=3) THESE VALUES CAN BE DERIVED FROM 'STANDARD' PROFILES. FOR
C     WATER VAPOR (J=4) NO SUCH PROFILES EXIST AND ONLY ALTFCT (1,4) IS
C     ENTERED AS 1.00. THE ALTITUDE INDEX,I, RUNS FROM 1 TO 21 CORRESPONDING
C     TO ALTITUDES FROM SEA  LEVEL TO 20 KM ABOVE SEA LEVEL IN 1 KM STEPS.

C     CRIPPS AND RAMSAY ( BR. JOURNAL OF DERMATOL.-JUNE 1970) HAVE
C     DETERMINED THE MINIMAL ENERGY TO PRODUCE AN ERYTHEMA (SUNBURN) AS
C     A FUNCTION OF WAVELENGTH. THE ARRAY ER(N) ARE THESE VALUES FOR WAVE-
C     LENGTHS FROM 290 TO 320 NANOMETERS (NM) IN  5 NM STEPS. THE ARRAY
C     ERLAM(N) ARE VALUES FOR WAVELENGTHS FROM 300 TO 320 NM IN 2 NM STEPS.
      ERLAM(1)=12
      ERLAM(2)=15
      ERLAM(3)=20
      ERLAM(4)=34
      ERLAM(5)=88
      ERLAM(6)=224
      ERLAM(7)=300
      ERLAM(8)=430
      ERLAM(9)=600
      ERLAM(10)=830
      ERLAM(11)=1160
      PI=3.1415926
      W=(2.0*PI)/365.
      SE=0.39784993
      D0=80.0

C     DEFINING NEEDED CONSTANTS, INPUT AND OPTIONS

C     THE OPTION NMAX--- WAVELENGTHS BETWEEN 290 NM AND 4000 NM ARE IN THE ARRAY
C     ILAM(N) WHOSE INDEX N RANGES FROM 1 TO 111. THE PROGRAM WILL COMPUTE
C     VALUES OF SOLAR RADIATION SPECTRA FOR WAVELENGTHS 290 NM TO ILAM(NMAX).

C     THE OPTION IPINT--- IF ONLY TOTAL INTEGRATED HORIZONTAL PLANE FLUXES
C     ARE DESIRED, ENTER IPINT = ZERO. FOR SPECTRAL INTENSITIES AND INTEGRATED
C     VALUES AT EACH WAVELENGTH, ENTER IPINT = 1.

C     THE OPTION NOSCAT =0 CAUSES THE PROGRAM TO ADVOID EXPLICITLY COMPUTING
C     DIFFUSE COMPONENT RADIATION. THIS OPTION IS IMPORTANT FOR THE CASES
C     WHERE THE ELEVATION IS EITHER SEA LEVEL, 5,000 FEET(1.54 KM) OR 10,000
C     FEET(3.08KM), SINCE MCCULLOUGH AND PORTER (PUBLISHED IN ECOLOGY, 1971)
C     HAVE COMPUTED DIRECT TO DIFFUSE RATIOS FOR A RAYLEIGH ATMOSPHERE FOR THESE
C     ELEVATIONS AND TYPICAL MEAN SUMMER AND WINTER UNDERLYING SURFACE ALBEDOS.
C     INDEPENDENT OF WHAT PRINT OPTION (IPINT) IS ENTERED THIS PROGRAM ALWAYS
C     RETURNS A VALUE FOR INTERGRATED DIRECT RAYLEIGH COMPONENT RADIATION.
C     IF DIFFUSE RADIATION COMPONENT IS DESIRED ENTER NOSCAT = 1.

C     THE OPTION IUV=0 CAUSES PROGRAM TO USE DAVE-FURAKAWA THEORY FOR
C     THE WAVELENGTH RANGE 290 .GE. LAMBDA .LE. 360 NANOMETERS. IF IUV = 1 AND
C     /OR WAVELENGTH GREATER THAN 360 NM, THEN SCATTERED RADIATION IS COMPUTED
C     USING SUBROUTINE GAMMA (SEE COMMENTS ABOUT THIS SUBROUTINE IN PROGRAM)

C     THE OPTION IEP=1 CAUSES PROGRAM TO SKIP WRITING OUT ERYTHEMA PRODUCT
C     SPECTRUM AND ONLY WRITE OUT TIME TO PRODUCE MINIMAL ERYTHEMA AND MAX.
C     WAVELENGTH OF ERYTHEMA PRODUCT CURVE

C     FILE 10 IS OUTPUT FOR TIME, SOLAR AND
C     IS OUTPUT FOR TIME, ZENITH ANGLES SENT TO OUTPUT FILE 'solrout'

C**** USER SUPPLIES--NMAX, IPINT, NOSCAT, IUV, IEP --VIA MENUS
C      THAT REPLACE FILE USRDAT
C     CALL MENUS VIA SUBROUTINE IOSOLR

      CALL IOSOLR

C     INPUTING PARAMETERS WHICH DETERMINE TERRESTRIAL RADIATION
C     DAY OF YEAR(D),TIME OF DAY(TD) AND LATITUDE(ALAT)
C     TIME OF DAY 0-24 HOURS (LOCAL STANDARD TIME)
C     PRESS(MB), AMR(KM), AND CMH2O ARE SEALEVEL PRESSURE, METEOROLOGICAL
C     RANGE AND TOTAL PRECIPITIBLE WATER VAPOR ABOVE SEA LEVEL RESPECTIVELY
C     ALTITUDE IS ENTERED AS THE INTERGER IALT (0-20 KM) .1 .GE. CMH2O. .LE. 2
C     REFL IS UNDERLYING SURFACE ALBEDO

C**** USER SUPPLIES--IALT, ALAT, PRESS, AMR, CMH2O,REFL,ALONC,PUNSH,IDA

C     LONGITUDE CORRECTION TO SOLAR NOON (IN HOURS), ALONC, IS FIGURED
C     BY ADDING 1/15 FOR EACH DEGREE OF LONGITUDE WEST OF THE STANDARD
C     MERIDIAN, THAT IS, THE PLACE IN THE TIME ZONE WHERE SOLAR NOON AND
C     NOON LOCAL STANDARD TIME ARE THE SAME.  SUBTRACT 1/15 FOR EACH
C     DEGREE EAST OF THE STANDARD MERIDIAN, FOR EXAMPLE, ALONC FOR
C     MILWAUKEE IS -2/15.  ALONC MAY VARY FROM +1/2 TO -1/2.

C     IF PUNSH IS NON ZERO,TIME OF DAY (IN MINUTES), INTEGRATED
C     GLOBAL RADIATION AND ZENITH ANGLE WILL BE WRITTEN TO ELEMENT solrout

C     **** USRDAT ****

C     READ JULIAN DAYS OF YEAR (1-365) INTO 'DAY'. IDA SPECIFIES HOW MANY
C     DAYS YOU ARE READING IN.

      DO 190 IDAY=IDAYST,IDA
        D=julday(IDAY)
C      SURFACE ABSORPTIVITY FOR SUB. DSUB CALCULATIONS
        SABNEW = 1.000 - REFLS(IDAY)
C      SURFACE REFLECTIVITY FOR SOLRAD CALCULATIONS
        REFL = REFLS(IDAY)
        SLE = SLES(IDAY)

C       THE FOLLOWING SECTION GENERATES DAY OF MONTH (IDOM) AND MONTH(M)
C       FROM DAY OF YEAR (D)

        DO 175 K=1,12
          IF(D-EOM(K)) 175,175,171
  171     IF(D-EOM(K+1))172,172,175
  172       M=K
          IF(D-31)173,173,174
  173       IDOM=int(D)
            GO TO 175
  174       IDOM=int(D-EOM(K))
  175   CONTINUE

C       INITIALIZING DATA COMPRESSION COUNTERS FOR OUTPUT TO BIOME.DATA
        NZCT=0
        NCT=1

C      ***********BEGIN PROGRAM CALCULATIONS************************************

      DO 200 IT=ISTART,IEND
C     Time of day (hour)
      TD=IT
C     Latitude (radians)
      AALAT=(ALAT*PI)/180.
C     Time of solar noon
      TSN=12.+ ALONC

C     COMPUTATION OF EARTH TO SUN RADIUS AND SOLAR DECLINATION

      AR2=(1.0+((2.0*EC)*COS(W*D)))
      ALON =(W*(D-D0))+((2.0*EC)*(SIN(W*D)-SIN(W*D0)))
      SINDE=SE*SIN(ALON)
      DEC=ASIN(SINDE)
C     CORRECTION FOR NORTHERN OR SOUTHERN HEMISPHERE; AMULT = 1 FOR N.
C     HEMISPHERE, -1 FOR SOUTHERN HEMISPHERE
      If (Hemis .eq. 1) then
        Amult = 1.0
       else
        If (Hemis .eq. 2) then
          Amult = -1.0
        Endif
      Endif
      DEC = DEC * AMULT
      DDEC=180.*DEC/PI

C     Checking zenith angle for possible skylight before sunrise or after sunset
      TestHA=((2.0*PI)/24.)*(TD-TSN)
      TestCZ=COS(AALAT)*COS(DEC)*COS(TestHA)+(SIN(AALAT)*SIN(DEC))
C     Zenith angle (radians)
      TestZ=ACOS(TestCZ)
C     Zenith angle (degrees)
      TestZZ = TestZ*180./PI
      If (TestZZ .lt. 107.) then
C      Possible skylight.  Check for Zenith angle 89 deg. or greater
        If (TestZZ .gt. 88.) then
C        Compute skylight based on G.V. Rozenberg. 1966. Twilight. Plenum Press.
C        p. 18,19.  First computing lumens: y = b - mx. Regression of data is:
          Elog = 41.34615384 - 0.423076923*TestZZ
          Skylum = 10.0**Elog
C        Converting lux (lumen/m2) to W/m2 on horizontal surface -
C        Twilight - scattered skylight before sunrise or after sunset
C        From p. 239 Documenta Geigy Scientific Tables. 1966. 6th ed. K. Diem, ed.
C        Mech./elect equiv. of light = 1.46*10^-3 kW/lumen
          SRINT(NMAX) = Skylum*1.46E-03
C        Global has no direct, only scattered before sunrise, after sunset
          GRINT(NMAX) = SRINT(NMAX)
         else
        Endif
       else
      Endif

C     THE FUNCTION H IS SOLAR HOUR ANGLE AT SUNRISE OR SUNSET. IF ABSOLUTE
C     VALUE OF H EXCEEDS 1, THEN THERE IS A LONG DAY OR LONG NIGHT.

C     TESTING COS H TO SEE IF EXCEEDS +1 OR -1
C     Equation 7, McCullough & Porter. 1971 Ecology
      TDTL=-(SIN(DEC)*SIN(AALAT))/(COS(AALAT)*COS(DEC))
      ABC=ABS(TDTL)
      IF(ABC-1.0)10,11,11
C      10 = normal day, 11 = long day or night

C     ESTABLISHING WHETHER A LONG DAY OR LONG NIGHT

   11 IF(TDTL-1.)40,41,41
C     40 = long day,  41 = long night
   41 IF(TDTL .EQ. 1) GO TO 42
      GO TO 43
   42 Continue
   43 Continue
C     NO SUN SKIP CALCULATIONS
      If (TestZZ .le. 106.0) then
       CGRNT=GRINT(NMAX)*.01434
      else
       CGRNT=0.D0
      Endif
      ZZ=90.
      DZSL=90.
      GO TO 1889
   40 H=PI
      GO TO 12

C     TESTING TIME OF DAY TO SEE IF BETWEEN SUNRISE AND SUNSET
C     EFFECT OF ATMOS. REFR. IS NOT TAKEN INTO ACCOUNT

   10 H=ACOS(TDTL)
      H=ABS(H)
   12 TS=TD-TSN
C     Compute time (hour angle) at sunrise
      HH=12.*H/PI
C     Check to see if it is before or after solar noon when TS = 0
      IF(TS)13,13,14
C     SUNRISE?
   13 CONTINUE
C     Check to see if the time before solar noon is the same as sunrise
C     time (hour angle, HH)
      IF(ABS(TS)-HH)20,20,30

C     SUNSET?
   14 CONTINUE

      IF(TS-HH)20,32,32

C     BEFORE SUNRISE  SKIP CALCULATIONS
   30 CONTINUE

      IF(TDTL .EQ. 1) GO TO 33
      GO TO 34
C     LONG NIGHT
   33 Continue
   34 Continue
      CDRNT = 0.D+0
      If (TestZZ .le. 106.0) then
       CGRNT=GRINT(NMAX)*.01434
      else
       CGRNT=0.0
      Endif
      ZZ=90.
      DZSL=90.
      GO TO 1889
C     AFTER SUNSET  SKIP CALCULATIONS
   32 IF(TDTL .EQ. 1) GO TO 35
      GO TO 36
   35 Continue
   36 Continue

      CDRNT = 0.D+0
      If (TestZZ .le. 106.0) then
       CGRNT=GRINT(NMAX)*.01434
      else
       CGRNT=0.0
      Endif
      ZZ=90.
      DZSL=90.
      GO TO 1889

C     COMPUTATION OF COSINE OF ZENITH ANGLE  (THE SUN IS UP)
C     HOUR ANGLE, HA, (RADIANS)
   20 HA=((2.0*PI)/24.)*(TD-TSN)

C     COS(ZENITH ANGLE)
      CZ=COS(AALAT)*COS(DEC)*COS(HA)+(SIN(AALAT)*SIN(DEC))
C     Zenith angle (radians)
      Z=ACOS(CZ)
C     Zenith angle (degrees)
      ZZ = Z*180./PI

C     Solar azimuth (degrees); south reference = 0 degrees.
      If (Amult .eq. 1.)then
C      Northern hemisphere
C      Calculating solar azimuth angle from equation 25 in Hosmer
       Altdeg = (90. - zz)
       Altrad = Altdeg*(pi/180.)
C      Cos(solar azimuth)
       tophalf = sin(dec)-sin(aalat)*sin(altrad)
       bothalf = cos(aalat)*cos(altrad)
       Cazsun = tophalf/bothalf
       If (Cazsun .lt. -0.9999999) then
        Azsun = pi
       else
        If (Cazsun .gt. 0.9999999) then
         Azsun = 0.00000
        else
         Azsun = ACOS(Cazsun)
        Endif
       Endif
       If (HA .le. 0.000) then
C       Morning
        DAZSUN = AZSUN*180./PI
       else
C       Afternoon
        DAZSUN = 180. + (180.- AZSUN*180./PI)
       Endif
      Endif
      If (Amult .lt. 0.)then
C      Southern hemisphere
C      Calculating solar azimuth angle from equation 27 in Hosmer
       Altdeg = (90. - zz)
       Altrad = Altdeg*(pi/180.)
C      Cos(solar azimuth)
       Cazsun = (sin(aalat)*sin(altrad)-sin(dec))/
     &    (cos(aalat)*cos(altrad))
C      This solar azimuth equation above works perfectly for the southern or northern hemisphere.
C      Error checking for instability in trig function
       If (Cazsun .lt. -0.9999999) then
        Azsun = pi
       else
        If (Cazsun .gt. 0.9999999) then
         Azsun = 0.00000
        else
         Azsun = ACOS(Cazsun)
        Endif
       Endif
       If (HA .le. 0.000) then
C       Morning
        DAZSUN = AZSUN*180./PI
       else
C       Afternoon
        DAZSUN = 360. - AZSUN*180./PI
       Endif
      Endif

c     kearney added this to account for horizon angle
      ahoriz = hori(minloc(abs(DAZSUN-azi), DIM=1))

C     COMPUTATION OF COSINE OF SLOPE ZENITH ANGLE IF SLOPING GROUND
      IF (SLOPE .GT. 0.000) THEN
       SLOPEL = (SLOPE*PI)/180.
       AZSLP  = (AZMUTH*PI)/180.
C      SUN AZIMUTH ANGLE
C      FROM EQ. 32 in Hosmer, G.L. 1910.  Practical Astronomy.
C      Wiley & Sons. New York (or Hosmer, G.L. & J.M. Robbins. 1948.
C      Practical Astronomy. Wiley & Sons. 355p. This equation is
C      NORTH referenced (N = 0 deg, S = 180 deg.)

C      NOTE: this equation does EXACTLY the same thing as the cosine equation above
C      AFTER all the quadrant corrections below. It is really superfluous, but has
C      been left in for historical purposes. WPPorter, 29 December, 1999.
C      Computing the azimuth when hour angle is known and
C      CORRECTING FOR HEMISPHERE (Amult)
       TAZSUN = SIN(HA)/(COS(AALAT)*TAN(DEC)-SIN(AALAT)*COS(HA))
C      sun azimuth (radians)
       AZSUN = ATAN(TAZSUN)*AMULT
C      Azimuth in degrees
       DAZSUN = AZSUN*180./PI
C      MAKING CORRECT NORTH REFERENCED SOLAR AZIMUTH
       IF (HA .LE. 0.00) THEN
C       MORNING - EAST OF REFERENCE
        IF (DAZSUN .LE. 0.00) THEN
C        1ST QUADRANT (0 - 90 deg.)
         DAZSUN = -1.* DAZSUN
        ELSE
C        2ND QUADRANT (90 - 180 deg.)
         DAZSUN = 180. - DAZSUN
        ENDIF
       else
C       AFTERNOON - WEST OF REFERENCE
        IF (DAZSUN .LT. 0.00) THEN
C        3RD QUADRANT (180 - 270 deg.)
         DAZSUN = 180. - DAZSUN
        ELSE
C        4TH QUADRANT (270 - 360 deg.)
         DAZSUN = 360. - DAZSUN
        ENDIF
       ENDIF
       If (HA .EQ. 0.0000)then
        DAZSUN = 180.
       Endif

C      SOLAR AZIMUTH ANGLE NORTH REFERENCED IN SELLER'S EQ 3.15
C      EQ. 3.15, p. 35; Sellers, W.D. 1965. Physical Climatology. U.
C      Chicago Press
C      KEY EQUATION:  used to correct solar on a slope in Subroutine DSUB
       CZSL = COS(Z)*COS(SLOPEL) +
     *  SIN(Z)*SIN(SLOPEL)*COS((DAZSUN - AZMUTH)*PI/180.)
       ZSL = ACOS(CZSL)
C      SOLAR ZENITH ANGLE FROM NORMAL TO SLOPE = SLOPE SOLAR ZENITH ANGLE
       DZSL = ZSL*180./PI
C      A CHECK FOR SUN OVER THE FLAT HORIZON, BUT NOT OVER THE SLOPE
C      90 DEGREES WANTED FOR COS (DZSL) = 0.0
       IF (DZSL .GT. 90.00000) THEN
        DZSL = 90.000000
       ENDIF
      ELSE
C      NO SLOPE
      ENDIF

C     COMPUTING OPTICAL AIR MASS (AIRMS) FROM KNOWLEDGE OF TRUE ZENITH ANGLE
C     USING THE ROZENBERG FORMULA (SEE P.159 IN BOOK 'TWILIGHT' BY G.V.
C     ROZENBERG, PLENUM PRESS, NEW YORK, 1966)
C     THE DIFFERENCE BETWEEN APPARENT AND TRUE ZENITH ANGLE IS NEGLECTED
C     FOR Z LESS THAN 88 DEGREES
C     VARIATION OF AIRMS WITH ALTITUDE IS IGNORED SINCE IT IS NEGLIGIBLE UP TO
C     AT LEAST 6 KM. ABOVE SEA LEVEL

      IF(Z-1.5358896)90,93,93
   93 REFR=16.+((Z-1.53589)*15.)/(PI/90.)
      REFR=(REFR/60.)*(PI/180.)
      Z=Z-REFR
   90 AIRMS=1./(COS(Z)+(0.025*EXP(-11.*COS(Z))))
      CZ=COS(Z)
      INTCZ=int(100.*CZ+1.0)
      ZZ=180.*Z/PI

C     ASSIGNMENT OF TOTAL ATMOSPHERIC OZONE FROM AVERAGE VALUES GIVEN ON PAGE
C     114 IN ROBINSON'S BOOK 'SOLAR RADIATION'. THE MATRIX OZ(LLAT,M) ARE THE
C     TOTAL ATMOSPHERIC OZONE THICKNESS (IN CM),WHERE INDEX M IS MONTH OF YEAR
C     ( JAN. .EQ. 1) AND LLAT IS TERRESTRIAL LATITUDE TO NEAREST 10 DEGREES
C      LLAT =1 IS 90 DEGREES SOUTH AND LLAT =19 IS 90 DEGREES NORTH

      TLAT=(ALAT+100.)/10.
      LLAT=int(TLAT)
      ALLAT=LLAT
      ALA =ALLAT+0.5
      IF(TLAT-ALA )150,160,160
  160 LLAT=LLAT+1
  150 OZONE=OZ(LLAT,M)

      HAD = 15.*(TD - TSN)

C     COMPUTATION OF DIRECT COMPONENT RADIATION (DRLAM)

      MALT=IALT+1
c     Kearney modification so that the altitude correction factor changes more smoothly. Fitted
c     polynomials to the ALTFCT data originally in SOLAR.DAT and now stored above
      ALT=ALTT/1000+1
      ALTFCT1=0.00007277*ALT**3+0.00507293*ALT**2-0.12482149*ALT+
     &1.11687469
      ALTFCT2=8.35656E-07*ALT**6-6.26384E-05*ALT**5+1.86967E-03*ALT**4-
     &2.82585E-02*ALT**3+2.26739E-01*ALT**2-9.25268E-01*ALT+1.71321E+00
      ALTFCT3=1.07573E-06*ALT**5-5.14511E-05*ALT**4+7.97960E-04*ALT**3-
     &4.90904E-03*ALT**2 + 2.99258E-03*ALT + 1.00238E+00
      ALTFCT4=1

c     mutliplier to correct hourly solar data for horizon angle
      if(altdeg.lt.ahoriz)then
c	   diffuse only - cut down to diffuse fraction      
      TDD(111+IT)=TDD(111+IT)* (0.12 + 0.83 * ((CCMINN(IDAY) + 
     &  CCMAXX(IDAY))/ 2 / 100)) ! from Butt et al. 2010 Agricultural and Forest Meteorology 150 (2010) 361–368
      endif

      DO 301 N=1,NMAX
       TLAM1=(PRESS/1013.)*TR(N)*ALTFCT1
       TLAM2=(25./AMR)*TA(N)*ALTFCT2
       TLAM3=(OZONE/0.34)*TO(N)*ALTFCT3
       TLAM4=TW(N)*SQRT(AIRMS*CMH2O*ALTFCT4)
       TLAM=((real(TLAM1,8)+TLAM2+TLAM3)*AIRMS)+TLAM4
C      MAKING SURE THAT AT LOW SUN ANGLES AIR MASS DOESN'T MAKE
C      TLAM TOO LARGE
       IF (TLAM .GT. 80.) THEN
        TLAM = 80.
       ENDIF
       part1 = SLAM(N)*AR2*CZ
       if(tlam.gt.0.00000000000000000)then
        part2 = EXP(-1.*TLAM)
       else
        part2 = 0.000000000000000
       endif
       if(part2.lt.1.0e-24)then
        drlam(n) = 0.00000000000000000000000000000000
       else
        DRLAM(N)=(part1*part2)/1000.
       endif
C     SO THE INTEGRATOR DOESN'T GET CONFUSED AT VERY LOW SUN ANGLES
       IF (DRLAM(N) .LT. 0.1E-24) THEN
        DRLAM(N) = 0.1E-24
       ENDIF
       DRRLAM(N)=(SLAM(N)*AR2*CZ)*EXP(-1.*real(TLAM1,8)*AIRMS)/1000.

c      kearney added this to account for horizon angle
       if(altdeg.lt.ahoriz)then
        DRLAM(N)=0.1E-24
        DRRLAM(N)=0.1E-24
       endif

C      COMPUTATION OF SKY (SRLAM) AND GLOBAL RADIATION (GRLAM)
       IF(NOSCAT .EQ.0) GO TO 400
       IF(IUV .GT. 0) GO TO 403
       IF(N .GT. 11) GO TO 400
C      THE OPTION IUV = ZERO HAS CAUSED THE PROGRAM TO ENTER THIS SECTION WHICH
C      COMPUTES SCATTERED RADIATION (SRLAM) FOR 290 NM TO 360 NM USING A THEORY
C      OF RADIATION SCATTERED FROM A RAYLEIGH (MOLECULAR) ATMOSPHERE WITH
C      OZONE ABSORPTION. THE FUNCTIONS NEEDED FOR THE COMPUTATION ARE STORED
C      AS FD(N,I) AND FDQ(N,I) WHERE N IS THE WAVELENGTH INDEX AND I IS
C      (ZENITH ANGLE+ 5)/5 ROUNDED OFF TO THE NEAREST INTERGER VALUE.
C      THE ARRAYS FD   AND FDQ    ARE FOR SEA LEVEL ( P = 1013 MB)
C      IF IUV IS ENTERED AS ANY INTERGER GREATER THAN ZERO THEN SRLAM IS
C      COMPUTED BEGINNING AT STATEMENT NUMBER 403, WHICH IS ALSO WHERE
C      COMPUTATION BEGINS FOR WAVELENGTHS GREATER THAN 360 NM (THAT IS,
C      N GREATER THAN 11).
       B=ZZ/5.0
       IA=int(B)
       C=B-IA
       IF(C.GT.0.50) GO TO 220
       I=IA+1
       GO TO 230
  220  I=IA+2
  230  FDAV=FD(N,I)
       FDQDAV=FDQ(N,I)
       SRLAM(N)=(SLAM(N)/PI)*( FDAV    +(FDQDAV    *(REFL/(1.-(REFL*
     X  S(N))))))/1000.
       SRLAM(N)=SRLAM(N)*AR2
       IF (N .GT. 7) GO TO 401
       EPLAM(N)=(SRLAM(N)+DRLAM(N))/ER(N)
       GO TO 401

C      COMPUTING SRLAM USING DAVE-WARTEN SUBROUTINE ( SEE REFERENCE- DAVE
C      AND WARTEN, ' PROGRAM FOR COMPUTING THE STOKES PARAMETERS OF THE
C      RADIATION EMERGING FROM A PLANE-PARRLLEL NON-ABSORBING, RAYLEIGH
C      ATMOSPHERE,  PAPER NO. 320-3248, IBM SCIENTIFIC CENTER, PALO ALTO, CALIF.)
C      TLAM1 IS THE RAYLEIGH OPTICAL DEPTH AND THE FUNCTIONS GAMR, GAML AND SBAR
C      ARE NEEDED TO COMPUTE SRLAM PER EQUATION (15) IN MCCULLOUGH AND PORTER.
C      EVEN  THOUGH THE SUBROUTINE GAMMA AS LISTED AT THE END OF THE DECK HAS
C      BEEN MODIFIED BY US TO BE MORE EFFICIENT, COMPUTATION OF GAMR, GAML, AND
C      SBAR FOR A GIVEN TLAM1 MAY TAKE AS LONG AS 30 SECONDS. MCCCULLOUGH AND
C      PORTER PRESENT RESULTS WHICH MAY PERMIT YOU TO ADVOID USING THE TIME
C      CONSUMING AND COSTLY GAMMA SUBROUTINE.

  403  IF(TLAM1 .LT. 0.03) GO TO 400
       CALL GAMMA (TLAM1,GAMR,GAML,SBAR)
       SRLAM(N)=(((real(GAML(INTCZ),8)+real(GAMR(INTCZ),8))/(2.*(1.-REFL
     &  *real(SBAR,8))))-EXP(-real(TLAM1,8)*AIRMS))*CZ*SLAM(N)*AR2/1000.
       GO TO 401
  400  SRLAM(N)=0.D0
  401  GRLAM(N)=SRLAM(N)+DRLAM(N)

C      INTEGRATION USING THE TRAPIZOIDAL RULE. THE NTH ELEMENT OF DRINT,DRRINT,
C      SRINT AND GRINT IS THE INTEGRAL DIRECT, DIRECT RAYLEIGH,SCATTERED,
C      AND GLOBAL ENERGY BELOW  THE NTH WAVELENGTH ILAM(N).

       IF(N .EQ. 1) GO TO 402
       AILAM(N)=ILAM(N)
       AILAM(N-1)=ILAM(N-1)
       DRINT(N)=DRINT(N-1)+(((AILAM(N)-AILAM(N-1))*DRLAM(N-1))
     X  +(0.5*(AILAM(N)-AILAM(N-1))*(DRLAM(N)-DRLAM(N-1))))
       DRRINT(N)=DRRINT(N-1)+(((AILAM(N)-AILAM(N-1))*DRRLAM(N-1))
     X  +(0.5*(AILAM(N)-AILAM(N-1))*(DRRLAM(N)-DRRLAM(N-1))))
       SRINT(N)=SRINT(N-1)+(((AILAM(N)-AILAM(N-1))*SRLAM(N-1))
     X  +(0.5*(AILAM(N)-AILAM(N-1))*(SRLAM(N)-SRLAM(N-1))))
       GRINT(N)=GRINT(N-1)+(((AILAM(N)-AILAM(N-1))*GRLAM(N-1))
     X  +(0.5*(AILAM(N)-AILAM(N-1))*(GRLAM(N)-GRLAM(N-1))))
       IF(N.GT.7) GO TO 301
       EPINT(N)=EPINT(N-1)+(((AILAM(N)-AILAM(N-1))*EPLAM(N-1))
     X  +(0.5*(AILAM(N)-AILAM(N-1))*(EPLAM(N)-EPLAM(N-1))))
       IF (EPINT(N)) 301,301,302
  302   REPINT(N)=1./(EPINT(N)*60.)
       GO TO 301
  402  SRINT(1)=0.D0
       DRRINT(1)=0.D0
       DRINT(1)=0.D0
       GRINT(1)=0.D0
       EPINT(1)=0.D0
  301 CONTINUE

C     FINDING LAMDA MAX OF ERYTHEMAL PRODUCT CURVE USING LINEAR FIT OF LOGS
C     REFERENCE--- PHYSICS IN MEDICINE AND BIOLOGY, VOLUME 15, PAGE 723(1970)
      IF (IPINT.EQ.0) GO TO 4000
      A(1)=(DLOG(GRLAM(4))-DLOG(GRLAM(3)))/(ILAM(4)-ILAM(3))
      A(2)=(DLOG(GRLAM(5))-DLOG(GRLAM(4)))/(ILAM(5)-ILAM(4))
      A(3)=(DLOG(GRLAM(6))-DLOG(GRLAM(5)))/(ILAM(6)-ILAM(5))
      A(4)=(DLOG(GRLAM(7))-DLOG(GRLAM(6)))/(ILAM(7)-ILAM(6))
      GPL(1)=GRLAM(3)
      GPL(2)=EXP(DLOG(GRLAM(3))+(A(1)*(302.-300.)))
      GPL(3)=EXP(DLOG(GRLAM(3))+(A(1)*(304.-300.)))
      GPL(4)=EXP(DLOG(GRLAM(4))+(A(2)*(306.-305.)))
      GPL(5)=EXP(DLOG(GRLAM(4))+(A(2)*(308.-305.)))
      GPL(6)=GRLAM(5)
      GPL(7)=EXP(DLOG(GRLAM(5))+(A(3)*(312.-310.)))
      GPL(8)=EXP(DLOG(GRLAM(5))+(A(3)*(314.-310.)))
      GPL(9)=EXP(DLOG(GRLAM(6))+(A(4)*(316.-315.)))
      GPL(10)=EXP(DLOG(GRLAM(6))+(A(4)*(318.-315.)))
      GPL(11)=GRLAM(7)
      DO 1510 MOM=1,11
 1510  EPL(MOM)=GPL(MOM)/ERLAM(MOM)
      EPL(12)=EPLAM(4)
      MMM=1
      ISOS=0
      DO 1500 MM=2,12
       IF(EPL(MM)-EPL(MMM))1500,1550,1520
 1520   MMM=MM
        GO TO 1500
 1550   ISOS=1
 1500 CONTINUE
      LAMAX=IELAM(MMM)
 4000 CONTINUE

C     THE FOLLOWING PRINTS OUT EITHER SPECTRAL AND INTEGRAL FLUXES FOR ALL
C     WAVELENGTHS(IPINT NOT = ZERO) OR TOTAL INTEGRATED FLUX (IPINT= ZERO)
C     THE OPTION IEP=1 CAUSES PROGRAM TO SKIP WRITING OUT ERYTHEMA PRODUCT
C     SPECTRUM AND ONLY WRITE OUT TIME TO PRODUCE MINIMAL ERYTHEMA AND MAX.
C     WAVELENGTH OF ERYTHEMA PRODUCT CURVE

C     INFORMATION CONTAINED IN MATRICES AT THIS STAGE OF PROGRAM IS AS FOLLOWS

C     ILAM(N)      NTH WAVELENGTH (NM)
C     DRLAM(N)     DIRECT COMPONENT RADIATION-NTH SPECTRAL ESTIMATE
C     DRRLAM(N)    DIRECT RAYLEIGH COMPONENT RADIATION-NTH SPECTRAL ESTIMATE
C     SRLAM(N)     DIFFUSE COMPONENT RADIATION-NTH SPECTRAL ESTIMATE
C     GRLAM(N)     GLOBAL RADIATION-NTH SPECTRAL ESTIMATE (GRLAM=DRLAM+SRLAM)
C     EPLAM(N)     ERYTHYMAL PRODUCT (GRLAM TIMES ACTION SPECTRUM ER)-NTH SPEC-
C                  TRAL ESTIMATE
C     IELAM(N)     WAVELENGTHS-BETWEEN 300 AND 320 NMS IN 2 NM STEPS
C     GPL(N)       GLOBAL RADIATION BETWEEN 300 AND 320 NM IN 2 NM STEPS
C     EPL(N)       ERYTHEMAL PRODUCT BETWEEN 300 AND 320 NM IN 2 NM STEPS
C     LAMAX        WVLTH AT WHICH MAX OF ERYTHEMA PROD. CURVE OCCURS

C     THE FOLLOWING CONTAIN THE INTEGRATED VALUES FOR THE SPECTRAL ESTIMATE
C     DESIGNATED, THAT IS, INTEGRATED OVER WAVELENGTH FROM 290 NM TO THE NTH
C     WAVELENGTH ILAM(N)-------

C     DRINT(N)     DIRECT RADIATION COMPONENT
C     DRRINT(N)    DIRECT RAYLEIGH RADIATION COMPONENT
C     SRINT(N)     SCATTERED RADIATION COMPONENT
C     GRINT(N)     GLOBAL RADIATION COMPONENT
C     EPINT(N)     ERYTHEMA PRODUCT
C     REPINT(N)    RECIPROCAL OF EPINT


      methour=(int(TD)+1)+24*(IDAY-1)
      DRLAMBDA(methour,1)=JULDAY(IDAY)
      DRLAMBDA(methour,2)=TD
      DRRLAMBDA(methour,1)=JULDAY(IDAY)
      DRRLAMBDA(methour,2)=TD
      SRLAMBDA(methour,1)=JULDAY(IDAY)
      SRLAMBDA(methour,2)=TD
      DRLAMBDA(methour,3:113)=DRLAM(1:111)*10 ! CONVERTING FROM mW/cm2-nm to W/m2-nm
      DRRLAMBDA(methour,3:113)=DRRLAM(1:111)*10 ! CONVERTING FROM mW/cm2-nm to W/m2-nm
      SRLAMBDA(methour,3:113)=SRLAM(1:111)*10 ! CONVERTING FROM mW/cm2-nm to W/m2-nm

      IF(IPINT .EQ. 0) GO TO 350
C     WRITING PLOT FILE
      IF (IT .LE. IEND) THEN
C      CONVERTING FROM mW/cm2-nm to W/m2-nm
       DO 997 K=1,NMAX
        GLOBAL = GRLAM(K)*10.
        WRITE(11,991)julday(IDAY),TD,ILAM(K),GLOBAL
997    CONTINUE
      ELSE
      ENDIF

      IF(ISOS .EQ. 1) GO TO 1700
      GO TO 200
 1700 Continue
      GO TO 200
  350 Continue

C     CONVERSION TO CAL/CM2-MIN.
C     NO SLOPE CORRECTION WANTED HERE BECAUSE WANT TO HAVE NUMBERS AVAILABLE FOR ANIMAL HEAT BALANCE, TOO.
C     THEY ARE EASIEST TO GET JUST FROM THE COS (ZENITH ANGLE). ONLY WANT OUTPUT FOR SLOPE SUN ANGLE, SO CAN DO CORRECTION IN DSUB FOR MICROMET CALC'S.
      CDRNT=DRINT(NMAX)*.01434
      CDRRNT=DRRINT(NMAX)*.01434
      CSRNT=SRINT(NMAX)*.01434
      CGRNT=GRINT(NMAX)*.01434
      IF(NOSCAT .EQ. 0) GO TO 198
      GO TO 1889
  198 Continue

 1889 IF (PUNSH-0.) 200,200,1990
C     FILLING PUNCH ARRAYS
 1990 CONTINUE
C     MICROMET STORAGE (MINUTES, CAL/CM2-MIN)
      TIME(IT)=TD*60.
      OCGRNT(IT)=CGRNT
      ODRNT(IT)=CDRNT
      OZZ(IT)=ZZ
      OZSL(IT)=DZSL

C     Creating Dataky.dat, input file for Micromet Program
C     COUNT A VALUE AND STORE IT
      NCT=NCT+1
      TSTOR(NCT)=TIME(IT)
      OSTORE(NCT)=OCGRNT(IT)
      ODSTOR(NCT)=ODRNT(IT)
      OZSTOR(NCT)=OZZ(IT)
      ZSLSTO(NCT)=OZSL(IT)

  200 CONTINUE

C     END OF 25 HOUR TIME LOOP **********************

C     POST PROCESSING OF COMPRESSED DATA
C     CORRECT TIME ARRAY FOR TIME ZONE, IF NEEDED
      IF (ABS(TIMCOR) .GT. 0.0000) THEN
       DO 201 JCT = 1,NCT
        TSTOR(JCT) = TSTOR(JCT) + TIMCOR*60.
        IF (TSTOR(JCT) .GT. 1440.) THEN
         TSTOR(JCT) = TSTOR(JCT) - 1500.
        ELSE
         IF (TSTOR(JCT) .LT. 0.000) THEN
          TSTOR(JCT) = TSTOR(JCT) + 1500.
         ELSE
         ENDIF
        ENDIF
201    CONTINUE
      ELSE
      ENDIF
C     REORDER ARRAYS FOR ASCENDING TIME?
      INTIME = 0
      DO 205 JCT = 2,NCT
       IF (TSTOR(JCT) .LT. TSTOR(JCT-1))THEN
C       INVERTED ORDER & MUST REARRANGE; STORE DISPLACEMENT
        INTIME = JCT-1
       ELSE
C       NO PROBLEM, CONTINUE
       ENDIF
205   CONTINUE
      IF (INTIME .GT. 0) THEN
C      TIME ZONE CORRECTION NEEDED
C      MUST REORDER ARRAYS TO PUT SMALLEST TIME FIRST
C ******* THERE IS AN ERROR HERE WHEN LONGITUDE .LT. REF LONG (88 vs. 90)
       DO 206 JCT = 1,NCT
        JTEST = JCT + INTIME
C       FIRST AN ARRAY OVERFLOW TEST
        IF (JTEST .GT. 25) THEN
         JTEST = JTEST - 25
        ELSE
         IF (JTEST .LT. 1) THEN
          JTEST = 25 + JTEST
         ELSE
         ENDIF
        ENDIF
        TNEW(JCT)   = TSTOR(JTEST)
        OSNEW(JCT)  = OSTORE(JTEST)
        ODNEW(JCT)  = ODSTOR(JTEST)
        OZNEW(JCT)  = OZSTOR(JTEST)
        ZSLNEW(JCT) = ZSLSTO(JTEST)
206    CONTINUE
C      NOW PUT BACK IN ORIGINAL ARRAY
       DO 207 JCT = 1,NCT
        TSTOR (JCT) =  TNEW(JCT)
        OSTORE(JCT) =  OSNEW(JCT)
        ODSTOR(JCT) =  ODNEW(JCT)
        OZSTOR(JCT) =  OZNEW(JCT)
        ZSLSTO(JCT) =  ZSLNEW(JCT)
207    CONTINUE
      ELSE
C      NO REORDERING NECESSARY
      ENDIF

C     CORRECT START OF DAY
      IF (TSTOR(2) .GT. 60.) THEN
C      NEED TO ENLARGE ARRAY WITH ZEROS ON FRONT END; ALREADY ONE
C      ZERO THERE AT TIME = 0
       NFT = NCT + 1
       DO 202 JCT = 2,NFT
        IF (JCT .EQ. 2)THEN
C        INSERT FIRST TIME BEFORE DATA & NEW 'ZERO' ENTRIES
         TNEW(JCT) = TSTOR(JCT) - 60.
         OSNEW(JCT)= 0.
         ODNEW(JCT)= 0
         OZNEW(JCT)= 90.
         ZSLNEW(JCT)= 90.
        ELSE
C        MOVE CURRENT ENTRIES BACK ONE PLACE IN ARRAY
         TNEW(JCT) = TSTOR (JCT-1)
         OSNEW(JCT) = OSTORE(JCT-1)
         ODNEW(JCT) = ODSTOR(JCT-1)
         OZNEW(JCT) = OZSTOR(JCT-1)
         ZSLNEW(JCT) = ZSLSTO(JCT-1)
        ENDIF
202    CONTINUE
C      SET SIZE OF NEW ARRAY
       NCT = NFT
C      PUT BACK IN ORIGINAL ARRAY IN CASE NO FURTHER CORRECTIONS
       DO 208 JCT = 1,NCT
        TSTOR (JCT) =  TNEW(JCT)
        OSTORE(JCT) =  OSNEW(JCT)
        ODSTOR(JCT) =  ODNEW(JCT)
        OZSTOR(JCT) =  OZNEW(JCT)
        ZSLSTO(JCT) =  ZSLNEW(JCT)
208    CONTINUE
       ELSE
C      NO NEED FOR FRONT END EXTENSION
      ENDIF
C     CHECK FOR ADD ON AT END OF ARRAY
      IF (TSTOR(NCT) .LT. 1380.) THEN
C      NEED TO ENLARGE ARRAY WITH ZEROS ON END
       NFT = NCT + 2
       DO 204 JCT = 1,2
        IF (JCT .EQ. 1) THEN
         TNEW(NCT+1) = TSTOR(NCT) + 60.
         OSNEW(NCT+1) = 0.0
         ODNEW(NCT+1) = 0.0
         OZNEW(NCT+1) = 90.
         ZSLNEW(NCT+1) = 90.
        ELSE
         TNEW(NCT+JCT) = 1440.
         OSNEW(NCT+JCT) = 0.0
         ODNEW(NCT+JCT) = 0.0
         OZNEW(NCT+JCT) = 90.
         ZSLNEW(NCT+JCT) = 90.
        ENDIF
204    CONTINUE
       NCT = NFT
C      END OF AFTERNOON ADD ON'S
C      PUT ADJUSTED VALUES BACK IN ORIGINAL ARRAYS FOR OUTPUT
       DO 710 ICORR = 1,NCT
        TSTOR (ICORR) = TNEW(ICORR)
        OSTORE(ICORR) = OSNEW(ICORR)
        ODSTOR(ICORR) = ODNEW(ICORR)
        OZSTOR(ICORR) = OZNEW(ICORR)
        ZSLSTO(ICORR) = ZSLNEW(ICORR)
710    CONTINUE
      ELSE
C      NO AFTERNOON ADD ON'S NEEDED, SOLAR GOES TO END OF DAY
      ENDIF

      IF (PUNSH - 0.) 1991,1991,1992
 1992 CONTINUE
      IEDUM=IEND+1
C     # OF TABLE ENTRIES FOR BIOME IS NCT

C     COMPUTING TIME OF SOLAR MAXIMUM (MINUTES)
      TSOLMX = TSN*60.

C     Air temperature calculations
C     SUNSET IN MILITARY TIME
      HHINT = AINT(HH)
      FRACT = (HH - HHINT) * 60.
      HSINT = AINT(TSN)
      FRACTS = (TSN - HSINT) * 60.
      TIMSS = (HSINT*100. + FRACTS) + (HHINT*100. + FRACT)
C     SUNRISE IN MILITARY TIME
      DELT = TSN - HH
      HRINT = AINT(DELT)
      FRACTR = (DELT - HRINT) * 60.
C     TIME OF AIR TEMPERATURE MINIMUM (NOTE: A KLUGE OF 200, I.E. 2
C     HOURS IS BEING USED TO MAKE SINEC WORK WHEN MORNING MINIMUM
C     IS NOT AT TRUE SUNRISE, THE ALGORITHM IS 2 HOURS OFF FOR SOME
C     UNKNOWN (6/7/89) REASON)
      IF(TIMINS(1) .GT. 0) THEN
       TSRHR = TIMINS(1) + 2.
      ELSE
       TSRHR = TIMINS(1)
      ENDIF
      TIMSR = (HRINT*100. + FRACTR) + (TSRHR*100.)
C     TIME OF AIR TEMPERATURE MAXIMUM
      TSNHR = TIMAXS(1) + TIMCOR
      TIMTMX =  (HSINT*100. + FRACTS) + (TSNHR*100.)
C     SETTING TMIN, TMAX FROM ARRAYS OBTAINED FROM SUB. IOSOLR
      TMIN = TMINN(IDAY)
      TMAX = TMAXX(IDAY)
      if(iday.lt.julnum)then
       TMIN2 = TMINN(IDAY+1)
       TMAX2 = TMAXX(IDAY+1)
      else
       TMIN2=TMIN
       TMAX2=TMAX
      endif
C     SETTING SNOW PRESENCE/ABSENCE, SURFACE ABSORPTIVITY AND % WET = f(DAY) DONE IN OSUB.

C     SETTING TIME OF MIN & MAX (HOURS BEFORE SUNRISE OR
C     AFTER SOLAR NOON)
      TIMIN = TIMSR
      TIMAX = TIMTMX
c     Kearney added this to allow for varying deep soil boundary
      if(microdaily.eq.1)then
       TANNUL = TANNULRUN(IDAY)
      endif
      CALL SINEC
C     " PUNCHING" AIR TEMPERATURES OUT DONE IN SUBROUTINE SINE

C     INSERT SINUSOIDALLY VARYING RELATIVE HUMIDITY, CLOUD COVER,
C      &/OR WIND SPEED BY CALL TO SINE
C     RELATIVE HUMIDITY
      IF ((ANS16 .EQ. 'N') .OR. (ANS16 .EQ. 'n'))THEN
C      SETTING VMIN, VMAX FROM ARRAYS OBTAINED FROM SUB. IOSOLR
       VMIN = RHMINN(IDAY)
       VMAX = RHMAXX(IDAY)
C      SETTING MAX & MIN TIMES RELATIVE TO SUNRISE & SOLAR NOON
C      TIME OF MINIMUM
       TSRHR = TIMINS(3)
       TIMSR = (HRINT*100. + FRACTR) + (TSRHR*100.)
C      TIME OF MAXIMUM at sunrise for relative humidity
       TSNHR = TIMAXS(3) + TIMCOR
       TIMTMX =  (HSINT*100. + FRACTS) + (TSNHR*100.)
       TIMIN = TIMTMX
       TIMAX = TIMSR
       IVAR = 1
       CALL VSINE(IVAR,VMIN,VMAX,TIMIN,TIMAX)
C         " PUNCHING" RELATIVE HUMIDITY OUT DONE IN SUBROUTINE VSINE
      ELSE
      ENDIF
C     CLOUD COVER
      IF ((ANS17 .EQ. 'N') .OR. (ANS17 .EQ. 'n'))THEN
C      SETTING VMIN, VMAX FROM ARRAYS OBTAINED FROM SUB. IOSOLR
       VMIN = CCMINN(IDAY)
       VMAX = CCMAXX(IDAY)
C      SETTING TIME OF MIN & MAX (INTEGER HOURS BEFORE SUNRISE OR
C      AFTER SOLAR NOON)
C      TIME OF MINIMUM
       TSRHR = TIMINS(4)
       TIMSR = (HRINT*100. + FRACTR) + (TSRHR*100.)
C      TIME OF MAXIMUM cloud cover
       TSNHR = TIMAXS(4) + TIMCOR
       TIMTMX =  (HSINT*100. + FRACTS) + (TSNHR*100.)
       TIMIN = TIMSR
       TIMAX = TIMTMX
       IVAR = 2
       CALL VSINE(IVAR,VMIN,VMAX,TIMIN,TIMAX)
C       " PUNCHING" CLOUD COVER OUT DONE IN SUBROUTINE VSINE
      ENDIF
C     WIND SPEED
      IF ((ANS18 .EQ. 'N') .OR. (ANS18 .EQ. 'n'))THEN
C      SETTING VMIN, VMAX FROM ARRAYS OBTAINED FROM SUB. IOSOLR
C      & CONVERTING FROM M/S TO CM/MIN
       VMAX = WNMAXX(IDAY)*6000.
       VMIN = WNMINN(IDAY)*6000.
C      SETTING TIME OF MIN & MAX (INTEGER HOURS BEFORE SUNRISE OR
C      AFTER SOLAR NOON)
C      TIME OF MINIMUM
       TSRHR = TIMINS(2)
       TIMSR = (HRINT*100. + FRACTR) + (TSRHR*100.)
C      TIME OF MAXIMUM wind speed
       TSNHR = TIMAXS(2) + TIMCOR
       TIMTMX =  (HSINT*100. + FRACTS) + (TSNHR*100.)
       TIMIN = TIMSR
       TIMAX = TIMTMX
       IVAR = 3
       CALL VSINE(IVAR,VMIN,VMAX,TIMIN,TIMAX)
C       " PUNCHING" WIND SPEED OUT DONE IN SUBROUTINE VSINE
      ENDIF

C     OUTPUT OPTION: IF ANS14 = M -> MICROMET; IF ANS14 = P -> PLOT
      IF ((ANS14 .EQ. 'M') .OR. (ANS14 .EQ. 'm')) THEN
C      "PUNCHING" GLOBAL SOLAR OUT FOR MICROMET INPUT
       do 2015 jj=1,25
        SOLS((CNT-1)*25+jj)=OSTORE(JJ)
2015   continue
C       " PUNCHING" ZENITH ANGLES OUT
       do 2013 jj=1,25
        ZENS((CNT-1)*25+jj)=OZSTOR(JJ)
2013   continue
       CNT=CNT+1
       IF (SLOPE .GT. 0.00) THEN
C       " PUNCHING" SLOPE ZENITH ANGLES OUT
        cnt=cnt-1
        do 2014 jj=1,25
         ZSLS((CNT-1)*25+jj)=ZSLSTO(JJ)
2014    continue
        CNT=CNT+1
       ENDIF
      ELSE
       IF (IPINT .EQ. 0) THEN
C      WRITE PLOT DATA: JULIAN DAY, HOUR OF DAY, SOLAR
C      CHECKING FIRST TO SEE IF SUN ON THE SLOPE IS DESIRED
        IF (SLOPE .GT. 0.00) THEN
         DO 730 IT = 1,25
C         SOLAR ON HORIZONTAL/unit area TO SOLAR NORMAL TO SURFACE
          HZEN = OZSTOR(IT)*PI/180.
          OSTORE(IT) = OSTORE(IT)/COS(HZEN)
C         SOLAR NORMAL TO SURFACE TO SOLAR ON SLOPE/unit area
          SZEN = ZSLSTO(IT)*PI/180.
          OSTORE(IT) = OSTORE(IT)* COS(SZEN)
730      CONTINUE
        ELSE
        ENDIF
       ELSE
       ENDIF
      ENDIF
C     ******************************************************

C     RESETTING ARRAYS TO ZERO FOR NEXT DAY'S CALCULATIONS
      DO 2000 JJ=1,IEDUM
      OZZ(JJ)=90.0
      OZSL(JJ)=90.0
 2000 OCGRNT(JJ)=0.D0

 1991 CONTINUE
  190 CONTINUE

      RETURN
      END
