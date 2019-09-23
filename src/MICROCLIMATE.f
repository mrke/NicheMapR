      subroutine microclimate(nn2,microinput1,julday1,SLES1,DEP1,
     &maxshades1,minshades1,Nodes1,timaxs1
     &,timins1,RHMAXX1,RHMINN1,CCMAXX1,CCMINN1,WNMAXX1,WNMINN1,TMAXX1
     &,TMINN1,REFLS1,PCTWET1,soilinit1,hori1,tai1,soilprop1,
     &moists1,rain1,tannulrun1,tides1,PE1,KS1,BB1,BD1,DD1,L1,LAI1
     &,TAIRhr1,RHhr1,WNhr1,CLDhr1,SOLRhr1,RAINhr1,ZENhr1,IRhr1,metout1
     &,soil1,shadmet1,shadsoil1,soilmoist1,shadmoist1,humid1,shadhumid1
     &,soilpot1,shadpot1,sunsnow1,shdsnow1,plant1,shadplant1,DRLAMBDA1
     &,DRRLAMBDA1,SRLAMBDA1)

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

C     This is a modified version of Warren P. Porter's microclimate mode that is now
C     called as a subroutine from R. It has been modified from the original version by
C     Michael Kearney to handle any time period, to have time-varying soil properties as
C     a function of temperature and water, and has an optional full water balance model
C     based on the equations and code of Campbell 1985 Soil Physics with Basic. It includes
c     an optional snow model that simulates snow accumulation and melt and its effect on
c     soil temperature. It no longer writes out temporary data (to dataky.dat).
c     The program and associated package is fully described in the publication
c     Kearney, M. R., and W. P. Porter. 2017. NicheMapR - an R package for biophysical modelling:
c     the microclimate model. Ecography 40:664–674.
c     The full package is available at http://github.com/mrke/NicheMapR/

C     SOLRAD computes clear sky solar radiation for anywhere
C     IOSOLR is for solrad input and min, max values &
c       for input of time dependent variables, air temperature,
C       relative humidity, cloud cover, and wind speed
c     GAMMA computes scattered radiation in the ultraviolet to high precision.
c     GAMMA calls DCHXY and DEXPI during its calculations.
C     SINEC computes time dependent values of air temperature
C     VSINE computes time dependent values of relative humidity,
C      cloud cover & wind speed.
c     RDCTRL controls reading of variables rearranged for input to the
c       meteorology program
c     RDTABL reads table data of time dependent input variables, e.g. air
c       temperature, solar radiation, etc.
c     PTTABL prints out time dependent tables
c     SFODE is the Gear Predictor-Corrector numerical integration program
c       that solves the simultaneous heat balance equations for temperatures.
c     DSUB is the subroutine containing the heat balance equations in
c       derivative form, which is called by SFODE.
c     OSUB outputs the microclimate calculations.

      use commondat
      IMPLICIT NONE

      integer, INTENT(IN) :: nn2

      double precision, DIMENSION(int(nn2)), intent(in) :: CCMAXX1(nn2),
     &CCMINN1(nn2),RHMAXX1(nn2),RHMINN1(nn2),TMINN1(nn2),TMAXX1(nn2),
     &WNMAXX1(nn2),WNMINN1(nn2),REFLS1(nn2),PCTWET1(nn2),LAI1(nn2),
     &SLES1(nn2),MAXSHADES1(nn2),MINSHADES1(nn2),
     &JULDAY1(nn2),rain1(nn2),tannulrun1(nn2),TAIRhr1(nn2*24),
     &RHhr1(nn2*24),WNhr1(nn2*24),CLDhr1(nn2*24),
     &SOLRhr1(nn2*24),RAINhr1(nn2*24),ZENhr1(nn2*24),IRhr1(nn2*24)
      double precision, DIMENSION(10,int(nn2)), intent(in) :: moists1,
     & Nodes1
      double precision, DIMENSION(24*int(nn2),19), intent(inout) ::
     & METOUT1,SHADMET1
      double precision, DIMENSION(24*int(nn2),12), intent(inout) ::
     &SOIL1,SHADSOIL1,soilpot1,shadpot1,humid1,shadhumid1,soilmoist1,
     &shadmoist1
      double precision, DIMENSION(24*int(nn2),14), intent(inout) ::
     &plant1,shadplant1
      double precision, DIMENSION(24*int(nn2),11), intent(inout) ::
     &sunsnow1,shdsnow1
      double precision, DIMENSION(24*int(nn2),113), intent(inout)
     & ::DRLAMBDA1,DRRLAMBDA1,SRLAMBDA1
      double precision, DIMENSION(24*int(nn2),3), intent(in) :: tides1
      double precision DEP1,timaxs1,timins1,hori1,tai1,
     &microinput1,soilprop1,PE1,KS1,BB1,BD1,L1,soilinit1,DD1

      double precision C,DEP,DTAU,ERR1,H,OUT,PAR,PTWET,SABNEW,soildp,
     & airdp,ZH,D0
      double precision T,TD,TI,TIME,TIMEF,WORK,shayd,altt,MAXSHD,WC
     & ,viewf
      double precision itair,icld,iwind,irelhum,rainfall,surflux,PE,KS
     & ,BB,BD,DD
      double precision RUF,SLE,ERR,Usrhyt,Z01,Z02,ZH1,ZH2,curmoist2,
     & soilprop,deps
      double precision ALAT,ALMINT,ALONC,ALONG,ALREF,AMINUT,AMULT
     & ,AZMUTH,CMH2O
      double precision HEMIS,PRESS,PUNSH,REFL,SLOPE
      double precision TANNUL,TIMCOR,TIMAXS,TIMINS,TSRHR,TSNHR,ep,
     & maxpool,L,LAI
      double precision tannul2,hori,azi,tai,ec,moist,minutes,condep,
     & rainmult
      double precision snownode,maxsnode1,snode,daysincesnow,lastday,
     &undercatch,rainmeltf,depp,Thconduct,Density,Spheat,minsnow,densfun
      double precision DRLAM,DRRLAM,SRLAM,snowcond,intercept
      double precision rww,pc,rl,sp,r1,im
      double precision snowdens,snowmelt,snowtemp,cursnow,qphase,
     & sumphase,sumphase2,snowage,prevden,QFREZE,xtrain,grasshade

      INTEGER I,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I20,KK,LL,I100
      INTEGER I91,I92,I93,I94,I95,I96,IPCH,IPRINT,cnt,II,I97,I98,I99
      Integer ENDMON,IDAY,IDMAIN,IFINAL,ILOCT,IOUT,trouble,maxcount,I101
      Integer ISHADE,ITEST,J,JULNUM,NOCON,NOPRNT,NOSUM,NOTRAN,NON,NAN
      Integer M,MM,MONLY,DOY,N,ND,NDEP,NDMAX,NDUM1,NKHT,NLIZ,maxerr
      INTEGER NOUT,NPOS,NUMRUN,NUMTYPS,writecsv,runshade,lamb,errout
      INTEGER IALT,IEND,IEP,IPINT,ISTART,IUV,NOSCAT,IDA,IDAYST,julstnd
      INTEGER microdaily,DOYF,DOYS,DOYF2,DOYS2,runmoist,evenrain,runsnow
      INTEGER errcount,HOURLY,rainhourly,IRmode,solonly

      CHARACTER(80) LABL1,LABL2,LABL3
      CHARACTER(3) IBLK,INAME,SYMBOL
      CHARACTER(1) solout,SINE,ANS14,SNSLOP
      CHARACTER(12) FNAME

      DIMENSION snownode(10),snode(10),qphase(10)
      DIMENSION microinput1(62)
      DIMENSION soilprop(10,5),soilprop1(10,5),moist(10)
      DIMENSION DEPS(21),curmoist2(18)
      DIMENSION TIMINS(4),TIMAXS(4)
      DIMENSION tai1(111),tai(111)
      DIMENSION Dep1(10),minutes(25)
      DIMENSION TIMINS1(4),TIMAXS1(4)
      DIMENSION soilinit1(20),hori1(24),densfun(4),sumphase2(10)
      DIMENSION julstnd(2),depp(30),Thconduct(30),Density(30),Spheat(30)
      DIMENSION hori(24),azi(24),PE(19),KS(19),BD(19),BB(19),DD(19)
      DIMENSION PE1(19),KS1(19),BD1(19),BB1(19),L(19),L1(19),DD1(19)
      DIMENSION SOILDP(20),AIRDP(6)
      DIMENSION DRLAM(111),DRRLAM(111),SRLAM(111)

      COMMON/WORK/WORK(1720)
      COMMON/LABEL/LABL1,LABL2,LABL3,FNAME,SINE,ANS14,SNSLOP
      COMMON/NONSCR/MM,N,TIME,TIMEF,DTAU,ERR1,H,NOUT,NDUM1,IPRINT
      COMMON/ARRAY/T(30),WC(20),C(20),DEP(30),OUT(101),IOUT(100),
     1 ITEST(23)
      COMMON/CARRAY/INAME(20),SYMBOL(23)
      COMMON/TABLE/TI(211),TD(211)
      COMMON/TABLE2/ILOCT(21)
      COMMON/OPSHUN/MONLY,NOPRNT,NOTRAN,NOCON,IPCH,ISHADE,NKHT,
     & NPOS,NLIZ,NOSUM
c      COMMON/WSTEDI/ASVLIZ(6),SHADE(4)
      COMMON/PAR/PAR(18)
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101
      COMMON/NDAY/ND
      common/solropt/solout
      COMMON/GROUND/SHAYD,ALTT,MAXSHD,SABNEW,PTWET,rainfall
      COMMON/WICHDAY/NUMRUN
      COMMON/DAYJUL/JULNUM,DOY
c    Variable soil properties data from Iomet1
      COMMON/SOYVAR1/Numtyps
      COMMON/SOYVAR2/Thconduct,Density,Spheat
      COMMON/DAILY/microdaily
      common/deep/tannul2
      common/atten/tai,ec
      common/VIEWFACT/VIEWF
      common/init/itair,icld,iwind,irelhum,iday
      COMMON/dataky/DEPS,CNT
      common/campbell/PE,KS,BB,BD,L,LAI,DD
      common/shaderun/runshade
      common/curmoist/curmoist2
      COMMON/WIOMT2/RUF
      Common/Hyte/Usrhyt
      COMMON/DMYCRO/Z01,Z02,ZH1,ZH2
      COMMON/CMYCRO/ZH,D0
      COMMON/NICHEMAPRIO/SLE,ERR,soilprop,surflux
      common/moistcom/moist,ep
      COMMON/ENDS/JULSTND
      COMMON/WIOCONS/PUNSH,ALAT,AMULT,PRESS,CMH2O,REFL,ALONC,TIMCOR,
     * AZMUTH,SLOPE,TSNHR,TSRHR,Hemis
      COMMON/WIOCONS2/IPINT,NOSCAT,IUV,IALT,IDAYST,IDA,IEP,ISTART,IEND
      COMMON/DAYS/TANNUL
      COMMON/DAYSS/TIMINS,TIMAXS
      COMMON/LATLONGS/AMINUT,ALONG,ALMINT,ALREF
      COMMON/SNOWPRED/snowtemp,snowdens,snowmelt,snownode,minsnow
     &,maxsnode1,snode,cursnow,daysincesnow,lastday,undercatch,rainmeltf
     &,densfun,snowcond,intercept,snowage,prevden,grasshade
      common/horizon/hori,azi
      common/soilmoist/condep,rainmult,maxpool
      common/soilmoist3/runmoist,evenrain,maxcount
      common/soimoist2/rww,pc,rl,sp,r1,im
      common/write/writecsv
      common/snowmod/runsnow,trouble
      common/dep/depp
      COMMON/SOILND/NON
      common/hour/hourly,rainhourly
      COMMON/LAMBDA/DRLAM,DRRLAM,SRLAM,LAMB
      common/IR/IRmode
      COMMON/melt/QFREZE,xtrain,qphase,sumphase,sumphase2
      common/errormsg/errout,maxerr,errcount
      COMMON/onlysol/solonly

      DATA IBLK/'   '/
      DATA IFINAL/1/

      DATA MINUTES/0,60,120,180,240,300,360,420,480,540,600,660,720,780
     &    ,840,900,960,1020,1080,1140,1200,1260,1320,1380,1440/
C     SIG=.8126E-10     1
C     RCS=.5            2
C     STC=.05
C     SSA=.7            4
C     HGT=40.           5
C     RUF=.05           6
C     BEG=0.            7
C     MON=0             8
C     PRT=60.           9
C     ERR=1.           10
C     END=0            11
C     SLE=1.           12
C     DAS=0            13
C     NON=10.          14
C     SUN=1.           15
C     PLT=0            16
C     FIN=0            17
C     STP=0            18
C     TIME DEPENDENT TABLES NEEDED
C     SOL  HOURLY SOLAR RADIATION
C     TAR  TEMPERATURE OF AIR AT REFERENCE HEIGHT
C     TDS  TEMPERATURE OF DEEP SOIL
C     VEL  AIR VELOCITY AT REFERENCE HEIGHT
C     THE FOLLOWING ARRAYS MUST ALSO BE INCLUDED
C     TIN  INITIAL TEMPERATURES MUST BE SAME NUMBER AS NUMBER OF NODES
C     DEP  DEPTHS IN THE SOIL AND HEIGHTS IN THE AIR
C     DROP OUT OPTION DEP=0,2.5,5,10,15,20,30,40,50,60
C     IOT  THE VARIABLES TO BE PRINTED ACCORDING TO KEY IN MANUAL
C     DROP OUT OPTION.  IOT=1,2,3,4,5,6

      ALLOCATE(SLES(nn2),RAIN(nn2),TIDES(nn2*24,3),metout(nn2*24,19
     &),shadmet(nn2*24,19),soil(nn2*24,12),shadsoil(nn2*24,12),soilmoist
     &(nn2*24,12),shadmoist(nn2*24,12),soilpot(nn2*24,12),shadpot(nn2*24
     &,12),humid(nn2*24,12),shadhumid(nn2*24,12),maxshades(nn2),
     &minshades(nn2),CCMAXX(nn2),CCMINN(nn2),RHMAXX(nn2),RHMINN(nn2)
     &,WNMAXX(nn2),WNMINN(nn2),TMAXX(nn2),TMINN(nn2),TANNULRUN(nn2),
     &REFLS(nn2),moists(10,nn2),intrvls(nn2),snowhr(nn2*25),
     &nodes(10,nn2),TDSS(nn2),TINS(20,nn2),TARS(nn2*25),RELS(nn2*25),
     &CLDS(nn2*25),VELS(nn2*25),SOLS(nn2*25),ZENS(nn2*25),ZSLS(nn2*25),
     &julday(nn2),LAIs(nn2),pctwet(nn2),rainhr(nn2*25),DRLAMBDA(nn2*24,
     &113),DRRLAMBDA(nn2*24,113),SRLAMBDA(nn2*24,113),sunsnow(nn2*24,11)
     &,shdsnow(nn2*24,11),plant(nn2*24,14),shadplant(nn2*24,14))
C    INITIALIZING MONTH OF YEAR COUNTER
      DOY=1
c    INITIALIZING DATAKY COUNTER
      CNT=1      
      minsnow=2.
      QFREZE=0.D0
      xtrain=0.D0
      M=0
      qphase(:)=0.D0
      sumphase=0.D0
      sumphase2(:)=0.D0
      depp(:)=0.D0
      dep(:)=0.D0
      curmoist2(:)=0.D0
      TI(:)=0.D0
      TD(:)=0.D0
      out(:)=0.D0
      density(:)=0.D0
      spheat(:)=0.D0
      thconduct(:)=0.D0
      snowhr(:)=0.D0
      metout1(:,:)=0.D0
      shadmet1(:,:)=0.D0
      soil1(:,:)=0.D0
      shadsoil1(:,:)=0.D0
      soilmoist1(:,:)=0.D0
      shadmoist1(:,:)=0.D0
      humid1(:,:)=0.D0
      shadhumid1(:,:)=0.D0
      soilpot1(:,:)=0.D0
      shadpot1(:,:)=0.D0
      sunsnow1(:,:)=0.D0
      shdsnow1(:,:)=0.D0
      plant1(:,:)=0.D0
      shadplant1(:,:)=0.D0
      metout(:,:)=0.D0
      shadmet(:,:)=0.D0
      soil(:,:)=0.D0
      shadsoil(:,:)=0.D0
      soilmoist(:,:)=0.D0
      shadmoist(:,:)=0.D0
      humid(:,:)=0.D0
      shadhumid(:,:)=0.D0
      soilpot(:,:)=0.D0
      shadpot(:,:)=0.D0
      sunsnow(:,:)=0.D0
      shdsnow(:,:)=0.D0
      plant(:,:)=0.D0
      shadplant(:,:)=0.D0
      DRLAMBDA(:,:)=0.D0
      DRRLAMBDA(:,:)=0.D0
      SRLAMBDA(:,:)=0.D0
c    Unpacking user input from R
      julnum=int(microinput1(1))

      runsnow=int(microinput1(36))
      if(runsnow.eq.1)then
       snownode(1)=minsnow
       snownode(2)=5
       snownode(3)=10
       snownode(4)=20
       snownode(5)=50
       snownode(6)=100
       snownode(7)=200
       snownode(8)=300
       lastday=1
       daysincesnow=0.D0
       snowage=0.D0
      endif
      do 1920 i=1,10
       soilinit1(i)=soilinit1(1)
1920  continue
      if(runsnow.eq.1)then
       N=18
      else
       N=10
      endif

      tides=tides1
      solonly=int(microinput1(60))
      ZH=microinput1(61)*100. ! Converting to cm
      D0=microinput1(62)*100. ! Converting to cm
      
c    do 901 i=1,2
c    julstnd(i)=julstnd1(i)
c901    continue
      RUF=microinput1(2)*100. ! Converting to cm
      SLES=SLES1
      soilprop=soilprop1
      SLE=SLES(1)
      ERR=microinput1(3)
      Usrhyt=microinput1(4)*100
      rain=rain1
      tannulrun=tannulrun1
      do 9191 i=1,25*nn2
       snowhr(i)=0.D0
9191  continue
      moists=max(moists1,0.01D+0)
      moist(1:10)=max(moists1(1:10,1),0.01D+0)
      surflux=0.
      do 919 i=1,24
      hori(i)=hori1(i)
      azi(i)=i*15
919   continue
      do 918 i=1,111
      tai(i)=tai1(i)
918   continue

      do 902 i=1,10
      DEP(i)=DEP1(i)
902   continue
      numtyps=int(microinput1(6))
      microdaily=int(microinput1(23))
      tannul2=microinput1(24)
      viewf=microinput1(26)
      snowtemp=microinput1(27)
      snowdens=microinput1(28)
      snowmelt=microinput1(29)
      undercatch=microinput1(30)
      condep=0.D0
      rainmult=microinput1(31)
      runshade=int(microinput1(32))
      PE=PE1
      KS=KS1
      BB=BB1
      BD=BD1
      DD=DD1
      L=L1
      runmoist=int(microinput1(33))
      maxpool=microinput1(34)
      evenrain=int(microinput1(35))
      rainmeltf=microinput1(37)
      writecsv=int(microinput1(38))
c    WRITE(I2,*)i,' ',j,' ',Thconds(i,j),' ',Thconds1(i,j)
      densfun(1)=microinput1(39)
      densfun(2)=microinput1(40)
      densfun(3)=microinput1(41)
      densfun(4)=microinput1(42)
      prevden=densfun(2)
      hourly=int(microinput1(43))
      rainhourly=int(microinput1(44))
      lamb=int(microinput1(45))
      IUV=int(microinput1(46))
      rww=microinput1(47)
      pc=microinput1(48)
      rl=microinput1(49)
      sp=microinput1(50)
      r1=microinput1(51)
      im=microinput1(52)
      maxcount=int(microinput1(53))
      IRmode=int(microinput1(54))
      errout=int(microinput1(55))
      maxerr=int(microinput1(56))
      snowcond=microinput1(57)/418.6*60
      intercept=microinput1(58)
      grasshade=microinput1(59)
      if((int(hourly).eq.1).or.(int(rainhourly).eq.1))then
      kk=0
      ll=0
      do 310 j=1,nn2
       do 311 i=1,25
        kk=kk+1
        if(i.lt.25)then
         ll=ll+1
         rainhr(kk)=rainhr1(ll)
        else
         rainhr(kk)=rainhr1(ll-1)
        endif
311    continue
310   continue
      endif

      IDA=int(microinput1(12))
      do 904 i=1,IDA
      Intrvls(i)=int(i)
904   continue
      do 905 i=1,numtyps
       do 906 j=1,IDA
      Nodes(i,j)=int(Nodes1(i,j))
906    continue
905   continue
      Z01=microinput1(7)*100.
      Z02=microinput1(8)*100.
      ZH1=microinput1(9)*100.
      ZH2=microinput1(10)*100.


      IDAYST=int(microinput1(11))
      EC=microinput1(25)
      HEMIS=microinput1(13)
      ALAT=microinput1(14)
      AMINUT=microinput1(15)
      ALONG=microinput1(16)
      ALMINT=microinput1(17)
      ALREF=microinput1(18)
      SLOPE=microinput1(19)
      AZMUTH=microinput1(20)
      ALTT=microinput1(21)
      CMH2O=microinput1(22)
      do 907 i=1,4
      timaxs(i)=timaxs1(i)
      timins(i)=timins1(i)
907   continue
      do 908 i=1,IDA
      julday(i)=julday1(i)
      LAIs(i)=LAI1(i)
      MAXSHADES(i)=MAXSHADES1(i)
      MINSHADES(i)=MINSHADES1(i)
      RHMAXX(i)=RHMAXX1(i)
      RHMINN(i)=RHMINN1(i)
      CCMAXX(i)=CCMAXX1(i)
      CCMINN(i)=CCMINN1(i)
      WNMAXX(i)=max(WNMAXX1(i),0.1D+0)
      WNMINN(i)=max(WNMINN1(i),0.1D+0)
      TMAXX(i)=TMAXX1(i)
      TMINN(i)=TMINN1(i)
      REFLS(i)=REFLS1(i)
      PCTWET(i)=PCTWET1(i)
908   continue



C    ****     COMPUTER READ - WRITE SETUP *************

      I1=5
      I2=7
      I3=10
      I4=11
      I5=4
      I6=3
      I7=2
      I8=1
      I9=12
      I10=13
      I11=14
      I12=15
      I20=20
      I91=21
      I92=22
      I93=23
      I94=24
      I95=25
      I96=26
      I97=27
      I98=28
      I99=29
      I100=30
      I101=31
c    setting solrad output option(y/n) for file Solrout
      solout='N'

      if(writecsv.eq.1)then
      OPEN (I3, FILE = 'metout.csv')
      write(I3,111) "JULDAY",",","TIME",",","TALOC",",","TAREF",",","RHL
     &OC",",","RH",",","VLOC",",","VREF",",","SNOWMELT",",","POOLDEP"
     &,",","PCTWET",",","ZEN",",","SOLR",",","TSKYC",",","DEW",","
     &,"FROST",",","SNOWFALL",",","SNOWDEP",",","SNOWDENS"
C     USE UNIT 13 FOR HOUR, SOIL DEPTH & SOIL TEMPERATURE OUTPUT
      OPEN (I10, FILE = 'soil.csv')
      write(I10,112) "JULDAY",",","TIME",",","DEP1",",","DEP2",","
     &,"DEP3",",","DEP4",",","DEP5",",","DEP6",",","DEP7",",","DEP8",","
     &,"DEP9",",","DEP10"
      if(runmoist.eq.1)then
      OPEN (I91, FILE = 'soilmoist.csv')
      write(I91,112) "JULDAY",",","TIME",",","WC1",",","WC2",","
     &,"WC3",",","WC4",",","WC5",",","WC6",",","WC7",",","WC8",","
     &,"WC9",",","WC10"
      OPEN (I92, FILE = 'humid.csv')
      write(I92,112) "JULDAY",",","TIME",",","HUM1",",","HUM2",","
     &,"HUM3",",","HUM4",",","HUM5",",","HUM6",",","HUM7",",","HUM8"
     &,",","HUM9",",","HUM10"
      OPEN (I100, FILE = 'plant.csv')
      write(I100,116) "JULDAY",",","TIME",",","TRANS",",","LEAFPOT",","
     &,"RTPOT1",",","RTPOT2",",","RTPOT3",",","RTPOT4",",","RTPOT5",","
     &,"RTPOT6",",","RTPOT7",",","RTPOT8"
     &,",","RTPOT9",",","RTPOT10"
       OPEN (I101, FILE = 'shdplant.csv')
      write(I101,116) "JULDAY",",","TIME",",","TRANS",",","LEAFPOT",","
     &,"RTPOT1",",","RTPOT2",",","RTPOT3",",","RTPOT4",",","RTPOT5",","
     &,"RTPOT6",",","RTPOT7",",","RTPOT8"
     &,",","RTPOT9",",","RTPOT10"
      OPEN (I95, FILE = 'soilpot.csv')
      write(I95,112) "JULDAY",",","TIME",",","PT1",",","PT2",","
     &,"PT3",",","PT4",",","PT5",",","PT6",",","PT7",",","PT8",","
     &,"PT9",",","PT10"
      OPEN (I1, FILE = 'soil_properties.csv')
      write(I1,113) "JULDAY",",","TIME",",","DEN1",",","DEN2",","
     &,"DEN3",",","DEN4",",","SPH1",",","SPH2",",","SPH3",",","SPH4",","
     &,"COND1",",","COND2",",","COND3",",","COND4"
      endif
      if(runsnow.eq.1)then
      OPEN (I7, FILE = 'sunsnow.csv')
      write(I7,113) "JULDAY",",","TIME",",","DEP1",",","DEP2",",","DEP3"
     &,",","DEP4",",","DEP5",",","DEP6",",","DEP7",",","DEP8",",","DEP9"
      endif
      if(lamb.eq.1)then
      OPEN (I97, FILE = 'DRLAMBDA.csv')
      write(I97,115) "JULDAY",",","TIME",",","290",",","295",",","300"
     &,",","305",",","310",",","315",",","320",",","330",",","340",","
     &,"350",",","360",",","370",",","380",",","390",",","400",",","420"
     &,",","440",",","460",",","480",",","500",",","520",",","540"
     &,",","560",",","580",",","600",",","620",",","640",",","660",","
     &,"680",",","700",",","720",",","740",",","760",",","780",","
     &,"800",",","820",",","840",",","860",",","880",",","900",",","920"
     &,",","940",",","960",",","980",",","1000",",","1020",",","1080"
     &,",","1100",",","1120",",","1140",","
     &,"1160",",","1180",",","1200",",","1220",",","1240",",","1260",","
     &,"1280",",","1300",",","1320",",","1380",",","1400",",","1420",","
     &,"1440",",","1460",",","1480",",","1500",",","1540",",","1580",","
     &,"1600",",","1620",",","1640",",","1660",",","1700",",","1720",","
     &,"1780",",","1800",",","1860",",","1900",",","1950",",","2000",","
     &,"2020",",","2050",",","2100",",","2120",",","2150",",","2200",","
     &,"2260",",","2300",",","2320",",","2350",",","2380",",","2400",","
     &,"2420",",","2450",",","2490",",","2500",",","2600",",","2700",","
     &,"2800",",","2900",",","3000",",","3100",",","3200",",","3300",","
     &,"3400",",","3500",",","3600",",","3700",",","3800",",","3900",","
     &,"4000"
      OPEN (I98, FILE = 'DRRLAMBDA.csv')
      write(I98,115) "JULDAY",",","TIME",",","290",",","295",",","300"
     &,",","305",",","310",",","315",",","320",",","330",",","340",","
     &,"350",",","360",",","370",",","380",",","390",",","400",",","420"
     &,",","440",",","460",",","480",",","500",",","520",",","540"
     &,",","560",",","580",",","600",",","620",",","640",",","660",","
     &,"680",",","700",",","720",",","740",",","760",",","780",","
     &,"800",",","820",",","840",",","860",",","880",",","900",",","920"
     &,",","940",",","960",",","980",",","1000",",","1020",",","1080"
     &,",","1100",",","1120",",","1140",","
     &,"1160",",","1180",",","1200",",","1220",",","1240",",","1260",","
     &,"1280",",","1300",",","1320",",","1380",",","1400",",","1420",","
     &,"1440",",","1460",",","1480",",","1500",",","1540",",","1580",","
     &,"1600",",","1620",",","1640",",","1660",",","1700",",","1720",","
     &,"1780",",","1800",",","1860",",","1900",",","1950",",","2000",","
     &,"2020",",","2050",",","2100",",","2120",",","2150",",","2200",","
     &,"2260",",","2300",",","2320",",","2350",",","2380",",","2400",","
     &,"2420",",","2450",",","2490",",","2500",",","2600",",","2700",","
     &,"2800",",","2900",",","3000",",","3100",",","3200",",","3300",","
     &,"3400",",","3500",",","3600",",","3700",",","3800",",","3900",","
     &,"4000"
      OPEN (I99, FILE = 'SRLAMBDA.csv')
      write(I99,115) "JULDAY",",","TIME",",","290",",","295",",","300"
     &,",","305",",","310",",","315",",","320",",","330",",","340",","
     &,"350",",","360",",","370",",","380",",","390",",","400",",","420"
     &,",","440",",","460",",","480",",","500",",","520",",","540"
     &,",","560",",","580",",","600",",","620",",","640",",","660",","
     &,"680",",","700",",","720",",","740",",","760",",","780",","
     &,"800",",","820",",","840",",","860",",","880",",","900",",","920"
     &,",","940",",","960",",","980",",","1000",",","1020",",","1080"
     &,",","1100",",","1120",",","1140",","
     &,"1160",",","1180",",","1200",",","1220",",","1240",",","1260",","
     &,"1280",",","1300",",","1320",",","1380",",","1400",",","1420",","
     &,"1440",",","1460",",","1480",",","1500",",","1540",",","1580",","
     &,"1600",",","1620",",","1640",",","1660",",","1700",",","1720",","
     &,"1780",",","1800",",","1860",",","1900",",","1950",",","2000",","
     &,"2020",",","2050",",","2100",",","2120",",","2150",",","2200",","
     &,"2260",",","2300",",","2320",",","2350",",","2380",",","2400",","
     &,"2420",",","2450",",","2490",",","2500",",","2600",",","2700",","
     &,"2800",",","2900",",","3000",",","3100",",","3200",",","3300",","
     &,"3400",",","3500",",","3600",",","3700",",","3800",",","3900",","
     &,"4000"
C     OPEN (I100, FILE = 'plant.csv')
C     write(I100,116) "JULDAY",",","TIME",",","TRANS",",","LEAFPOT",","
C    &,"RTPOT1",",","RTPOT2",",","RTPOT3",",","RTPOT4",",","RTPOT5",","
C    &,"RTPOT6",",","RTPOT7",",","RTPOT8",","
C    &,"RTPOT9",",","RTPOT10"
      endif
      if(runshade.eq.1)then
C     USE UNIT I13 FOR above ground micromet OUTPUT when % shade = 100.
      OPEN (I12, FILE = 'shadmet.csv')
      write(I12,111) "JULDAY",",","TIME",",","TALOC",",","TAREF",",","RH
     &LOC",",","RH",",","VLOC",",","VREF",",","SNOWMELT",",","POOLDEP"
     &,",","PCTWET",",","ZEN",",","SOLR",",","TSKYC",",","DEW",","
     &,"FROST",",","SNOWFALL",",","SNOWDEP"",","SNOWDENS"
C     USE UNIT 14 FOR HOUR, SOIL DEPTH & SOIL TEMPERATURE OUTPUT when % shade = 100.
      OPEN (I11, FILE = 'shadsoil.csv')
      write(I11,112) "JULDAY",",","TIME",",","DEP1",",","DEP2",","
     &,"DEP3",",","DEP4",",","DEP5",",","DEP6",",","DEP7",",","DEP8"
     &,",","DEP9",",","DEP10"
      if(runsnow.eq.1)then
      OPEN (I8, FILE = 'shdsnow.csv')
      write(I8,113) "JULDAY",",","TIME",",","DEP1",",","DEP2",",","DEP3"
     &,",","DEP4",",","DEP5",",","DEP6",",","DEP7",",","DEP8",",","DEP9"
      endif
      if(runmoist.eq.1)then
      OPEN (I93, FILE = 'shadmoist.csv')
      write(I93,112) "JULDAY",",","TIME",",","WC1",",","WC2",","
     &,"WC3",",","WC4",",","WC5",",","WC6",",","WC7",",","WC8",","
     &,"WC9",",","WC10"
      OPEN (I94, FILE = 'shadhumid.csv')
      write(I94,112) "JULDAY",",","TIME",",","HUM1",",","HUM2",","
     &,"HUM3",",","HUM4",",","HUM5",",","HUM6",",","HUM7",",","HUM8"
     &,",","HUM9",",","HUM10"
      OPEN (I96, FILE = 'shadpot.csv')
      write(I96,112) "JULDAY",",","TIME",",","PT1",",","PT2",","
     &,"PT3",",","PT4",",","PT5",",","PT6",",","PT7",",","PT8",","
     &,"PT9",",","PT10"
C     OPEN (I101, FILE = 'shadplant.csv')
C     write(I101,116) "JULDAY",",","TIME",",","TRANS",",","LEAFPOT",","
C    &,"RTPOT1",",","RTPOT2",",","RTPOT3",",","RTPOT4",",","RTPOT5",","
C    &,"RTPOT6",",","RTPOT7",",","RTPOT8"
C    &,",","RTPOT9",",","RTPOT10"
      endif
      endif
      endif

111   format(A8,A1,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1
     &,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1,A8,A1)
112   format(A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1
     &,A7,A1,A7,A1)
113   format(A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1
     &,A7,A1,A7,A1)
115   format(A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1
     &,A7,A1,
     &A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,
     &A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,
     &A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,
     &A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,
     &A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,
     &A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,
     &A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,
     &A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,
     &A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,
     &A7,A1,A7,A1,A7,A1)
116   format(A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1,A7,A1
     &,A7,A1,A7,A1,A7,A1,A7,A1)
C    USE UNIT 12 FOR SOLAR OUTPUT - FILE 6 WILL BE THE CONSOLE
      If ((solout .eq. 'Y').or.(solout .eq. 'Y'))then
        OPEN (I9, FILE = 'solrout')
      Endif

C    ******* END MICRO READ - WRITE SETUP *************

c     Defining maximum number of iteration days for the integrator
c     (REPEATS OF A DAY TO GET a STEADY PERIODIC OF THE DAY).
      NDMAX=10

C    DEFINING THE NUMBER OF DAYS TO REPEAT TO GET A STEADY PERIODIC
c    Kearney changed this for daily simulations
      if(microdaily.eq.1)then
          if(doy.eq.1)then
           ND = 3
           else
           ND = 1
          endif
      else
           ND = 3
      endif

      par(6) = RUF
      PAR(12) = SLE
      PAR(10) = ERR
      PAR(5) = microinput1(5)*100.

C     SET UP SOIL NODES, DEPTHS, AIR NODES, HEIGHTS FOR MICROMET
C    SETTING DEFAULT SOIL (NON) AND AIR NODES (NAN) HEIGHTS & DEPTHS
      if(runsnow.eq.1)then
       NON = 18
      else
       NON = 10
      endif
      NAN = 3

c      Assign depth values to soil depth array as specified in BLKDATA.FOR.  NOTE: THE USER NEEDS ACCESS TO NODE SPECIFICATION IN MICROv.DAT 3/12/11
      DO 300 I=1,NON
       SOILDP(I)=DEP(I)
300   CONTINUE

C    DEFINE AIR HEIGHTS (cm)
      AIRDP(1) = 200.
c    AirDP(2) now provided by user  6/27/98
      AIRDP(2) = Usrhyt
      if(RUF.gt.0.5)then
      AIRDP(3) = RUF+.5
      else
      AIRDP(3) = 0.5
      endif

C    SETTING UP THE NODE NUMBERS AND VALUES FOR THE DEP ARRAY
      NDEP = NON + NAN
      DO 52 I = 1,NAN
        AIRDP(I) = (-1.)*AIRDP(I)
52    CONTINUE
C    WRITING HEIGHTS & DEPTHS TO DATA INPUT FILE
      do 53 i=1,3
       DEPS(i)=AIRDP(i)
53    continue
      do 54 i=1,NON
       DEPS(i+3)=SOILDP(i)
54    continue
      if(hourly.eq.1)then
       TD(112:136)=1.
      endif
      CALL SOLRAD
c     Checking on maximum iteration days
      if (ND .gt. NDMAX) then
       ND=NDMAX
      endif
C    DEFINING A DEBUG PRINTOUT OPTION FOR MAIN (0 = NO PRINTOUT)
      IDMAIN = 0

      DO 10 I=1,23
10    ITEST(I)=0
      NOUT=6
C       INSERTING DEFAULT OUTPUT VARIABLES TO BE PRINTED
      DO 12 I=1,6
 12   IOUT(I)=I
      NDEP=10
C     ZEROING WORK, DEPTH AND OUTPUT ARRAYS
      DO 24 I=1,560
 24   WORK(I)=0.D0
      DO 25 I = 1,100
   25   OUT(I) = 0.0
      IPRINT=1

C    ***********************************************************
      errcount=0
200   CONTINUE
      rain=rain1
      LAI=LAIs(DOY)
      TD(10)=TDSS(DOY)
      TD(11)=TDSS(DOY)
      TI(10)=0.D0
      TI(11)=1440.
      if(int(HOURLY).eq.1)then
       DOYS=(DOY)*24-23
       DOYF=DOY*24
       DOYS2=(DOY)*25-24
       DOYF2=DOY*25
       TD(12:35)=TAIRhr1(DOYS:DOYF)
       TD(37:60)=RHhr1(DOYS:DOYF)
       TD(62:85)=CLDhr1(DOYS:DOYF)
       TD(87:110)=max(WNhr1(DOYS:DOYF),0.1D+0)*6000.
       TD(112:135)=SOLRhr1(DOYS:DOYF)/ 4.185 / 10000. * 60.
       TD(36)=TAIRhr1(DOYF)
       TD(61)=RHhr1(DOYF)
       TD(86)=CLDhr1(DOYF)
       TD(111)=max(WNhr1(DOYF),0.1D+0)*6000.
       TD(136)=SOLRhr1(DOYF)/ 4.185 / 10000. * 60.
       if(ZENhr1(1).gt.0)then
        TD(137:160)=ZENhr1(DOYS:DOYF)
        TD(161)=ZENhr1(DOYF)
       else
        TD(137:161)=ZENS(DOYS2:DOYF2)
       endif
       if(IRhr1(1).gt.0)then
        TD(187:210)=IRhr1(DOYS:DOYF)/ 4.185 / 10000. * 60.
        TD(211)=IRhr1(DOYF)/ 4.185 / 10000. * 60.
       else
        TD(187:211)=IRhr1(DOYF)
       endif
       DOYS=(DOY)*25-24
       DOYF=DOY*25
       TD(162:186)=ZSLS(DOYS:DOYF)
      else
       DOYS=(DOY)*25-24
       DOYF=DOY*25
       TD(12:36)=TARS(DOYS:DOYF)
       TD(37:61)=RELS(DOYS:DOYF)
       TD(62:86)=CLDS(DOYS:DOYF)
       TD(87:111)=VELS(DOYS:DOYF)
       TD(112:136)=SOLS(DOYS:DOYF)
       TD(137:161)=ZENS(DOYS:DOYF)
       TD(162:186)=ZSLS(DOYS:DOYF)
      endif
      if(int(HOURLY).eq.2)then ! passing in hourly solar (and possibly zenith angle) only
       DOYS=(DOY)*24-23
       DOYF=DOY*24
       DOYS2=(DOY)*25-24
       DOYF2=DOY*25
       TD(112:135)=SOLRhr1(DOYS:DOYF)/ 4.185 / 10000. * 60.
       TD(136)=SOLRhr1(DOYF)/ 4.185 / 10000. * 60.
       if(ZENhr1(1).gt.0)then
        TD(137:160)=ZENhr1(DOYS:DOYF)
        TD(161)=ZENhr1(DOYF)
       else
        TD(137:161)=ZENS(DOYS2:DOYF2)
       endif
      endif
      TI(12:36)=minutes
      TI(37:61)=minutes
      TI(62:86)=minutes
      TI(87:111)=minutes
      TI(112:136)=minutes
      TI(137:161)=minutes
      TI(162:186)=minutes
      TI(187:211)=minutes
      if(runsnow.eq.1)then
       ndep=21
       ii=18
      else
       ndep=13
       ii=10
      endif

      if(microdaily.eq.1)then
       do 344 i=1,ii
        T(i)=SOILINIT1(i)
344    continue
      else
       if(ifinal.eq.1)then
        T(1:ii)=TINS(1:ii,doy)
       endif
      endif
      DEP(1:ndep)=DEPS(1:ndep)


      if(runsnow.eq.1)then
       N=18
      else
       N=10
      endif
      MM=N-1

      IF(IOUT(1).LT.0) IPRINT=0
      ERR1=ERR
      DTAU=PAR(9)
      TIMEF=PAR(11)
      TIME=PAR(7)
      DO 325 I=1,N
 325  WORK(I+520)=T(I)
      IF(IPRINT.EQ.0  ) THEN
C      SKIP ANY KIND OF OUTPUT
        GO TO 3000
       ELSE
        CONTINUE
      ENDIF

 3000 CONTINUE

C    INTERNAL LOOP IN MICROCLIMATE PROGRAM TO RUN 2 SUCCESSIVE IDENTICAL DAYS
C    EXCEPT FOR A CHANGE IN THE AMOUNT OF SHADE. THE MONTHLY CHANGES IN SHADE
C    ARE DONE IN SOLRAD.  THE VARIABLE SHAYD IS CHANGED MONTHLY AND ALL THE
C    MINIMUMS ARE RUN FIRST FOR MAXIMUM SUN CONDITION, THEN THE WHOLE YEAR IS
C    REPEATED (NUMRUN = 2) FOR THE MINIMUM SUN CONDITION.  10/11/04  W. PORTER

c500   CONTINUE

C    SETTING THIS MONTH'S SHADE VALUES FOR CALCULATIONS IN DSUB
C    THIS DEPENDS ON WHETHER IT IS THE FIRST RUN (MINIMUM SHADE YEAR SIMULATION)
C    OR THE 2ND RUN (MAXIMUM SHADE FOR THE SAME YEAR SIMULATED)
      IF(IFINAL.EQ.1)THEN
        IF(NUMRUN.EQ.1)THEN
          SHAYD = MINSHADES(DOY)
          IF(SHAYD.gt.MAXSHADES(DOY))then
              SHAYD=MAXSHADES(DOY)-0.1
          ENDIF
          MAXSHD = MAXSHADES(DOY)
         ELSE
C        IT'S THE SECOND RUN WHERE THE VARIABLE SHAYD IS THE MAXIMUM VALUE FOR THE
C        MAX. SHADE BOUNDING CONDITION
          SHAYD = MAXSHADES(DOY)
          IF(SHAYD.gt.MAXSHADES(DOY))then
              SHAYD=MAXSHADES(DOY)-0.1
          ENDIF
          MAXSHD = MAXSHADES(DOY)
          maxsnode1=0.D0 ! reset snow settings
          do 2202 i=1,8
           snode(i)=0
2202      continue
          lastday=1
          daysincesnow=0
c         reset reflectivity and pct wet that might have been changed by snow
          do 9008 i=1,IDA
           REFLS(i)=REFLS1(i)
           PCTWET(i)=PCTWET1(i)
9008      continue
        ENDIF
      ENDIF
C     INCREMENT DAY COUNTER TO GET STEADY PERIODIC
      IFINAL = IFINAL + 1

C     CHECK FOR RESET OF IFINAL, IF # OF DAYS TO BE REPEATED, ND,
C    HAS BEEN EXCEEDED, THEN QUIT THIS DAY'S SIMULATION
      IF (IFINAL .GT. ND) THEN
        IFINAL = 1
C      NEED SOME RESET OF NUMRUN VALUE TO EITHER DO REPEAT DAY WITH NEW SHADE VALUE OR QUIT
      ENDIF
      IF(ITEST(1).EQ.-1) STOP
C    CALL THE PREDICTOR-CORRECTOR NUMBERICAL INTEGRATOR TO DO EACH DAY FOR SET VALUES
      CALL SFODE
      if(microdaily.eq.1)then
       ND=1
       do 101 i=1,ii
        soilinit1(i)=work(520+i)
101    continue
      endif

C    LOOPING FOR THE SECOND DAY WITH MAX SHADE
      IF(NUMRUN.EQ.2)THEN
       IF((DOY.LE.JULNUM).and.(runshade.eq.1))THEN
        DO 3006 I=1,N
 3006    T(I)=WORK(I+520)
        rain=rain1! reset rain in case set to zero in osub to escape an instability
        GO TO 200
       ELSE
        do 910 j=1,19
         do 909 i=1,24*julnum
          metout1(i,j)=metout(i,j)
          if(runshade.eq.1)then
           shadmet1(i,j)=shadmet(i,j)
          endif
909      continue
         i=1
910     continue
        do 912 j=1,12
         do 911 i=1,24*julnum
          soil1(i,j)=soil(i,j)
          soilmoist1(i,j)=soilmoist(i,j)
          humid1(i,j)=humid(i,j)
          soilpot1(i,j)=soilpot(i,j)
          if(runshade.eq.1)then
           if(j.le.12)then
            shadsoil1(i,j)=shadsoil(i,j)
            shadmoist1(i,j)=shadmoist(i,j)
            shadhumid1(i,j)=shadhumid(i,j)
            shadpot1(i,j)=shadpot(i,j)
           endif
          endif
911      continue
         i=1
912     continue
        do 916 j=1,11
         do 915 i=1,24*julnum
          sunsnow1(i,j)=sunsnow(i,j)
          if(runshade.eq.1)then
           if(j.le.11)then
            shdsnow1(i,j)=shdsnow(i,j)
           endif
          endif
915      continue
         i=1
916     continue
        do 920 j=1,14
         do 917 i=1,24*julnum
          plant1(i,j)=plant(i,j)
          if(runshade.eq.1)then
           if(j.le.14)then
            shadplant1(i,j)=shadplant(i,j)
           endif
          endif
917      continue
         i=1
920     continue
        do 913 j=1,113
         do 914 i=1,24*julnum
          DRLAMBDA1(i,j)=DRLAMBDA(i,j)
          DRRLAMBDA1(i,j)=DRRLAMBDA(i,j)
          SRLAMBDA1(i,j)=SRLAMBDA(i,j)
914      continue
         i=1
913     continue
        if(writecsv.eq.1)then
         close (i3)
         close (i10)
         close (i7)         
         if(lamb.eq.1)then
          close (i97)
          close (i98)
          close (i99)
         endif
         if(runmoist.eq.1)then
          close (i91)
          close (i92)
          close (i95)
          close (i100)
          close (i101)
         endif
         if(runshade.eq.1)then
          close (i8)
          close (i11)
          close (i12)
          if(runmoist.eq.1)then
           close (i93)
           close (i94)
           close (i96)
          endif
         endif
        endif
        DEALLOCATE(SLES,RAIN,TIDES,metout
     &,shadmet,soil,shadsoil,soilmoist,shadmoist
     &,soilpot,shadpot,humid,shadhumid,plant,shadplant,
     &maxshades,minshades,CCMAXX,CCMINN,RHMAXX,RHMINN
     &,WNMAXX,WNMINN,TMAXX,TMINN,TANNULRUN,sunsnow,shdsnow
     &,REFLS,moists,intrvls,snowhr,nodes,TDSS,
     &TINS,TARS,RELS,CLDS,VELS,SOLS,ZENS,ZSLS,LAIs,
     &PCTWET,julday,rainhr,DRLAMBDA,DRRLAMBDA,SRLAMBDA)
        RETURN
       ENDIF
      ENDIF

C    CHECK FOR THE END OF A YEAR, DAY COUNTER INCREMENTED IN OUTPUT SUBROUTINE OSUB
      ENDMON = JULNUM + 1
      IF (DOY.EQ.ENDMON)THEN
       if(runshade.eq.0)then
        do 9101 j=1,19
         do 9091 i=1,24*julnum
          metout1(i,j)=metout(i,j)
9091     continue
         i=1
9101    continue
        do 9121 j=1,12
         do 9111 i=1,24*julnum
          soil1(i,j)=soil(i,j)
          soilmoist1(i,j)=soilmoist(i,j)
          humid1(i,j)=humid(i,j)
          soilpot1(i,j)=soilpot(i,j)
9111     continue
         i=1
9121    continue
        do 9161 j=1,11
         do 9151 i=1,24*julnum
          sunsnow1(i,j)=sunsnow(i,j)
9151     continue
         i=1
9161    continue
        do 9201 j=1,14
         do 9171 i=1,24*julnum
          plant1(i,j)=plant(i,j)
9171     continue
         i=1
9201    continue
        do 9131 j=1,113
         do 9141 i=1,24*julnum
          DRLAMBDA1(i,j)=DRLAMBDA(i,j)
          DRRLAMBDA1(i,j)=DRRLAMBDA(i,j)
          SRLAMBDA1(i,j)=SRLAMBDA(i,j)
9141     continue
         i=1
9131    continue
        if(writecsv.eq.1)then
         close (i3)
         close (i10)
         close (i7)         
         if(lamb.eq.1)then
          close (i97)
          close (i98)
          close (i99)
         endif
         if(runmoist.eq.1)then
          close (i91)
          close (i92)
          close (i95)
          close (i100)
          close (i101)
         endif
        endif
        DEALLOCATE(SLES,RAIN,TIDES,metout
     &,shadmet,soil,shadsoil,soilmoist,shadmoist
     &,soilpot,shadpot,humid,shadhumid,plant,shadplant,
     &maxshades,minshades,CCMAXX,CCMINN,RHMAXX,RHMINN
     &,WNMAXX,WNMINN,TMAXX,TMINN,TANNULRUN,sunsnow,shdsnow
     &,REFLS,moists,intrvls,snowhr,nodes,TDSS,
     &TINS,TARS,RELS,CLDS,VELS,SOLS,ZENS,ZSLS,LAIs,
     &PCTWET,julday,rainhr,DRLAMBDA,DRRLAMBDA,SRLAMBDA)
        RETURN
       else
        NUMRUN = 2
        DOY = 1
        errcount=0
        rain=rain1! reset rain in case set to zero in osub to escape an instability
        GO TO 200
       endif
      ENDIF

      DO 3002 I=1,N
 3002 T(I)=WORK(I+520)
      errcount=0
      rain=rain1! reset rain in case set to zero in osub to escape an instability
C    DO ANOTHER DAY
      GO TO 200
      RETURN
      END