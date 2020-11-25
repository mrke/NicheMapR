      subroutine gads(lat51,lon51,relhum1,season1,optdep1)
      implicit none
ccccc ------------------------------------------------------------------c
c     create global distributions of microphysical and optical aerosol  c
c     properties on the base of the GADS database.                      c
c                                                                       c
c     version 2.0                                                       c
c     -----------                                                       c
c     11.02.97                                                          c
c                                                                       c
c     version 2.1                                                       c
c     -----------                                                       c
c     19.11.97 - new file format in ../optdat/                          c
c              - dimensions for phase function now 112.                 c
c                                                                       c
c     version 2.2                                                       c
c     -----------                                                       c
c     21.01.98 - files winter.dat, summer.dat updated                   c
c              - output improved                                        c
c              - optical depth corrected                                c
c                                                                       c
c     version 2.2a                                                      c
c     -----------                                                       c
c     06.03.98 - error with calculation of mass corrected.              c
c                                                                       c
c                                                                       c
c     06.03.98                                                 M. Hess  c
ccccc ------------------------------------------------------------------c


ccccc -----------------------------------------------------------------c
c      Calculation of the optical aerosol data from the microphysical  c
c      raw data.                                                       c
c                                                                      c
c      ATLOPT is a modified version of the OPAC and                    c
c      controls the call of the subroutines:                           c
c                                                                      c
c      - HEAD4                                                         c
c      - LOCATE                                                        c
c      - D4RAW                                                         c
c      - OPTCOM                                                        c
c      - OPTPAR                                                        c
c      - OUT4                                                          c
c                                                                      c
c      These subroutines are attached to ATLOPT.                       c
c                                                                      c
c      ab 14.12.92 another format of the new built-in Mie-Rechnungen   c
c      ab 04.11.93 new Parameter scattering, absorption, omega ratio   c
c      ab 13.05.94 OPTCOM: Russ no longer swells                       c
c      ab 27.07.94 opt. Thickness calculated for all wavelengths       c
c                                                                      c
c      18.11.97 GADS 2.1                                               c
c                                                                      c
c      18.11.97                                                M. Hess c
ccccc -----------------------------------------------------------------c

      integer prnr,acnr,njc,rht,mlamb,mhum,nhum,jnopar,nop,ibuf,nih,nlay
      integer niw,njh,nlmal,ntape,ip,i,iwel,nwel,ihum,il,ih,ilat,ilmal
      integer ilon,j,late,lata,loni,lone,lona,lati,nprog,nil,mbuf,nseas
      integer latx,lonx,na,norm,mixnor,season
      real n,numden,abs,acmr,absbuf,alamb,ahum,extbuf,bacbuf,asybuf
      real asy,hstra,hfta,hm,h1,h0,ext,bre,kbuf,brebuf,bimbuf,bim
      real lat5,lon5,relhum,optdep,bac,boundl
      real boundu,mlay,nh,nl
      real nltyp,pha,phabuf
      real sca,scabuf,sis,sisbuf,parlay,ncomp,khum,natyp,mcomp,mopar
      double precision lat51,lon51,relhum1,optdep1,season1
      
      dimension optdep1(25,2),optdep(25,2),jnopar(13),nil(10),hfta(10)
      dimension hstra(10),h0(2,10),h1(2,10),hm(2,10)
      dimension kbuf(20),extbuf(20),scabuf(20),absbuf(20)
      dimension sisbuf(20),asybuf(20),bacbuf(20),phabuf(112,20),
     *brebuf(20),bimbuf(20),rht(2),n(2),njc(2),acnr(5,2),acmr(5,2),
     *nh(2),atn(2),pat(2),ext(2,5),sca(2,5),abs(2,5),sis(2,5),asy(2,5),
     *comnam(20),bac(2,5),pha(2,5,112),bre(2,5),bim(2,5),ncomp(10),
     *nltyp(10),parlay(10,2),boundl(10),boundu(10),
     *opanam(13),optnam(13),alamb(61),khum(8),ahum(8),nhum(8),
     * chum(8)
      character*1 ws,dum
      character*2 chum
      character*3  atn,pat
      character*8  opanam,optnam
      character*4  comnam
      character*7  cseas
      character*11 tseas
      character*20 catyp
      character*30 typnam
      character*50 area
c     CHARACTER(len=255) :: cwd

      common /prog/   nprog
      common /profi/  nil,hfta,hstra,h0,h1,hm
      common /buffer/ ibuf,kbuf,extbuf,scabuf,absbuf,sisbuf,asybuf,
     * bacbuf,phabuf,brebuf,bimbuf,mbuf
      common /mipoi/  latx,lonx,nl,prnr,rht,n,njc,acnr,acmr,nh,atn,pat
      common /oppoi/  ext,sca,abs,sis,asy,bac,pha,bre,bim
      common /atyp/   natyp,mcomp,ncomp,numden,catyp,typnam,comnam
      common /layer/  nlay,mlay,nltyp,parlay,boundl,boundu
      common /season/ nseas,cseas,tseas
      common /opar/   mopar,jnopar,nop,opanam,optnam
      common /wavel/  alamb,mlamb,niw
      common /hum/    khum,ahum,nih,nhum,mhum,chum
      common /geog/   lata,late,lati,lona,lone,loni,na,area
      common /norm/   norm,mixnor
      common /r/      lat5,lon5,relhum,season,optdep
CCCCC -----------------------------------------------------------------C
c     some definitions for this version                                c
CCCCC -----------------------------------------------------------------C
      alamb = (/0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,
     *            0.9,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.2,3.39,3.5,3.75,
     *            4.0,4.5,5.0,5.5,6.0,6.2,6.5,7.2,7.9,8.2,8.5,8.7,9.0,
     *            9.2,9.5,9.8,10.0,10.6,11.0,11.5,12.5,13.0,14.0,14.8,
     *            15.0,16.4,17.2,18.0,18.5,20.0,21.3,22.5,25.0,27.9,30.,
     *            35.0,40.0/)
      mlamb = 61
      ahum = (/0.,50.,70.,80.,90.,95.,98.,99./)
      nhum = (/0,50,70,80,90,95,98,99/)
      mhum = 8
      chum = (/'00','50','70','80','90','95','98','99'/)

      comnam = (/'inso','waso','soot','ssam','sscm','minm','miam',
     *             'micm','mitr','suso','stco','stma','cucc','cucp',
     *             'cuma','fog-','cir1','cir2','cir3','    '/)

      optnam = (/'ext.coef','sca.coef','abs.coef','sisc.alb',
     *             'asym.par','op.depth',
     *             '        ','turb.fac','li.ratio','pha.func',
     *             'ext.rat ','abs.rat ',
     *             '        '/)

      jnopar = (/1,1,1,1,1,1,0,0,1,0,0,0,0/)
      nop = 7
      
      lat5=real(lat51,4)
      lon5=real(lon51)
      season=int(season1)
      relhum=real(relhum1)
      niw=1
      njh=1
      nih=1
      lata=int(lat5)
      late=int(lat5)
      lati=5
      lona=int(lon5)
      lone=int(lon5)
      loni=5
      nlmal=1
      norm=1
      nprog=4
      ntape=22
      ip=0
      do i=1,13
       if (jnopar(i).eq.1) then
        ip=ip+1
        opanam(ip)=optnam(i)
       end if
      end do

ccccc -----------------------------------------------------------------c
c     Query what should be plotted                                     c
ccccc -----------------------------------------------------------------c

      if(season.eq.1)then
       ws='w'
      else
       ws='s'
      endif

ccccc ----------------------------------------------------------------c
c     Input: wavelength                                               c
ccccc ----------------------------------------------------------------c

      nwel=25

      iwel=1
      do 9999 iwel=1,nwel

       if (ws.eq.'w') then
        open(ntape,file='extdata/glodat/winter.dat')
        read (ntape,'(a1)') dum
        cseas='winter '
       else 
        open(ntape,file='extdata/glodat/summer.dat')
        read (ntape,'(a1)') dum
        cseas='summer '
       end if

ccccc ----------------------------------------------------------------c
c     Input: humidity                                                 c
ccccc ----------------------------------------------------------------c
       ihum=int(relhum)

CCCCC -----------------------------------------------------------------C
C     READING THE HEIGHT PROFILES from the file TAPE9                  c
C     						                                           C 
C     HM    : EFFECTIVE LAYER THICKNESS (HOMOGENEOUS DISTRIBUTION)     C
C     HFTA  : FREE TROP LAYER THICKNESS. AEROSOLS IN KM                C
C     HSTRA : LAYER THICKNESS OF THE STRATOSPH. AEROSOLS IN KM         C
C     NIL   : NUMBER OF LAYERS                                         C
CCCCC -----------------------------------------------------------------C

       call prof

ccccc -----------------------------------------------------------------c
c      Labeling of the output files                                    c
ccccc -----------------------------------------------------------------c

       if(iwel.eq.1)then
        call head4 !(iwel,ihum)
       endif

ccccc -----------------------------------------------------------------c
c      Loop over all required wavelengths and moisture classes         c
ccccc -----------------------------------------------------------------c

       !do il=1,niw
        !do ih=1,njh
ccccc -----------------------------------------------------------------c
c     Loop over all required geographic coordinates                    c
ccccc -----------------------------------------------------------------c
         !do ilat=lata,late,-lati
          !do ilmal=1,nlmal
           !do ilon=lona,lone,loni
       il=1
       ih=1
       ilat=lata
       ilmal=1
       ilon=lona

ccccc -----------------------------------------------------------------c
c      Reading the raw data from the files TAPE201, TAPE207:	       c
c     ------------------------------------------------------	       c
C     LAT       : LATITUDE                                             C
C     LON       : LONGITUDE                                            C
C     NL        : NUMBER OF AEROSOL LAYERS                             C
C                 (=2 FOR MARITIME-MINERAL,=1 FOR MARI.)               C
C     PRNR      : PROFIL NUMBER                                        C
C     NT        : NUMBER OF AEROSOL TYPE                               C
C     PAT       : AEROSOL PROFIL TYPE                                  C
C     NH        : NUMBER OF REL. HUMIDITY CLASS                        C
C     N         : TOTAL NUMBER CONCENTRATION                           C
C     NJC       : NUMBER OF AEROSOL COMPONENT                          C
C     ACNR      : AEROSOL COMPONENT                                    C
C     ACMR      : MIXING RATIO                                         C
C                 (PARTIAL NUMBER CONCENTRATION/TOTAL NUMBER CONC.)    C
ccccc -----------------------------------------------------------------c

            call d4raw (ilat,ilon,ntape)

ccccc -----------------------------------------------------------------c
c      Reading of the optical raw data from the files winter.dat and   c
c      summer.dat                                                      c
ccccc -----------------------------------------------------------------c

            call optcom (iwel,ihum)

ccccc -----------------------------------------------------------------c
c      Calculation of the optical parameters at the current grid point c
ccccc -----------------------------------------------------------------c

            call optpar(iwel,ihum)

           !end do
          !end do
         !end do
        !end do
       !end do
       close (ntape)
9999  continue

      do 9998 i=1,25
       do 9997 j=1,2
        optdep1(i,j)=real(optdep(i,j),8)
9997   continue
9998  continue

      return
      end

CCCCC *****************************************************************C
      SUBROUTINE PROF
C     *****************************************************************C
C                                                                      C
C     -----------------------------------------------------------------C
C     READING THE HEIGHT PROFILES from the file profiles.dat and the   C
C     extinction coefficients of the upper atmosphere from extcof.dat  C
CCCCC -----------------------------------------------------------------C
      integer mlamb,nprof,ip,IIP,IL,nil,IWL,niw
      real alamb,hfta,hstra,h0,h1,hm,extfta,extstr,WAVE,EXTFT,EXTST
      
      dimension alamb(61),hfta(10),hstra(10),nil(10),
     *h0(2,10),h1(2,10),hm(2,10),EXTFTA(61),EXTSTR(61)
      
      common /wavel/  alamb,mlamb,niw
      common /profi/ nil,hfta,hstra,h0,h1,hm
      COMMON /FTASTR/ EXTFTA,EXTSTR

CCCCC -----------------------------------------------------------------C
C     THERE ARE 7 PROFILE TYPES. The following data are imported:      C
c                                                                      c
c           iip: Profile type number                                   c
c       nil(ip): Number of layers for type ip (as in tape201 etc.)     c
c      hfta(ip): Height of the layer for the free troposphere. Aerosol c
c     hstra(ip): Height of the stratospheric aerosol layer             c
c     h0(il,ip): Lower limit of layer il                               c
c     h1(il,ip): Upper limit of the layer                              c
c     hm(il,ip): effective thickness of the layer il for the type ip   c
CCCCC -----------------------------------------------------------------C

      open (8,file='extdata/profiles.dat')

      nprof=7
      DO IP=1,nprof
       READ(8,8010) IIP,NIL(IP),HFTA(IP),HSTRA(IP)
       READ(8,8020) (H0(IL,IP),H1(IL,IP), HM(IL,IP),IL=1,NIL(IP) )
 8010  FORMAT(I3,I3,2F8.2)
 8020  FORMAT(2F5.1,F10.3)
      end do

      close (8)

CCCCC -----------------------------------------------------------------C
C     Reading in the extinction coefficients for the upper atmosphere: C
C                                                                      C
C     EXTINCTION COEFFICIENT -	FREE TROPOSPHERIC AEROSOL  +	       C
C     EXTINCTION COEFFICIENT -	  STRATOSPHERIC   AEROSOL	           C
C                                                                      C
C     SKIP THE FIRST TWO ROWS from TAPE9            		           C
CCCCC -----------------------------------------------------------------C

      open (9,file='extdata/extback.dat')
      IL=1
      READ(9,'(/)')
      DO IWL=1,mlamb
       READ(9,*) WAVE,EXTFT,EXTST
       EXTFTA(IWL)=EXTFT
       EXTSTR(IWL)=EXTST
      end do

      close (9)

      RETURN
      END

ccccc *****************************************************************c
      subroutine head4 !(il,ih)
c     *****************************************************************c
c                                                                      c
c     -----------------------------------------------------------------c
c     Labeling of the TAPE10 output file		        	           c
ccccc -----------------------------------------------------------------c

      integer kop,iop,nseas,jnopar,nop,mlamb,niw,nih,nhum,mhum
      integer lata,late,lati,lona,lone,loni,na,norm,mixnor,jnangle
      integer angle,nia
      real mixrat,numden,dens,rmin,khum,ahum,rmod,rmax
      real sigma,ncomp,natyp,mcomp,mopar,alamb
      character*2 chum
      character*4 comnam
      character opanam*8,cseas*7,tseas*11,optnam*8,opnam*8
      character catyp*20,area*50,typnam*30

      dimension khum(8),ahum(8),nhum(8),chum(8),comnam(20)
      dimension rmin(10,8),rmax(10,8),rmod(10,8),opnam(10),
     * mixrat(10),dens(10,8),ncomp(10),jnopar(13),opanam(13),optnam(13)
      dimension jnangle(112),angle(112),alamb(61)
      
      common /season/ nseas,cseas,tseas
      common /geog/   lata,late,lati,lona,lone,loni,na,area
      common /hum/    khum,ahum,nih,nhum,mhum,chum
      common /norm/   norm,mixnor
      common /numdis/ sigma,rmin,rmax,rmod,mixrat,dens
      common /atyp/   natyp,mcomp,ncomp,numden,catyp,typnam,comnam
      common /opar/   mopar,jnopar,nop,opanam,optnam
      common /wavel/  alamb,mlamb,niw
      common /angle/  jnangle,angle,nia

      if (jnopar(10).eq.1) then
       kop=nop-1
      else
       kop=nop
      end if

      do iop=1,kop
       opnam(iop)=opanam(iop)
      end do

      return
      end

CCCCC *****************************************************************C
      subroutine d4raw (lat,lon,ntape)
C     *****************************************************************C
C                                                                      C
C     -----------------------------------------------------------------C
C     READING IN THE DATA FROM THE RAW DATA-FILES TAPE201-TAPE212        C
CCCCC -----------------------------------------------------------------C

      IMPLICIT CHARACTER*3 (Z)

      integer latx,lonx,prnr,acnr,njc,rht,nprog,JC,ntape,l,nl,ic,il,lat
      integer lon
      real n,itest,nh,acmr,sum
      character*3 atn,pat

      dimension rht(2),n(2),njc(2),acnr(5,2),acmr(5,2),nh(2),atn(2)
      dimension pat(2),itest(2),zat(11),zrh(8)

      common /prog/  nprog
      common /mipoi/ latx,lonx,nl,prnr,rht,n,
     *               njc,acnr,acmr,nh,atn,pat
      common /test/  itest
      common /dat/   zat,zrh

  111 READ(NTAPE,1020,end=999) LATX,LONX,NL,PRNR,
     +   ATN(1),PAT(1),RHT(1),N(1),NJC(1),
     +      ( ACNR(JC,1),ACMR(JC,1),JC=1,3 )
 1020 FORMAT(2I4,2I3,(2A3,I3,E10.3,I4,3(I3,E10.4)))

      IF(LATX.NE.LAT .OR.  LONX.NE.LON) THEN
c         print*,' Achtung, falsche Koordinaten: ',latx,lonx
       GOTO 111
      END IF

      IF (NJC(1).GT.3) THEN
       READ(NTAPE,1025,end=999) ( ACNR(JC,1),ACMR(JC,1),JC=4,NJC(1))
 1025  FORMAT(37X,3(I3,E10.4))
      END IF

      IF(NL.NE.1) THEN
       DO 10 L=2,NL
       READ(NTAPE,1021,end=999)  ATN(L),PAT(L),RHT(L),N(L),NJC(L),
     +  ( ACNR(JC,L),ACMR(JC,L),JC=1,3 )
 1021  FORMAT(14X,2A3,I3,E10.3,I4,5(I3,E10.4))

       IF (NJC(L).GT.3) THEN
        READ(NTAPE,1025,end=999) ( ACNR(JC,L),ACMR(JC,L),JC=4,NJC(L))
       END IF

   10  CONTINUE
      END IF

      do l=1,nl
       sum=0.
       do ic=1,njc(l)
        sum=sum+acmr(ic,l)
       end do
      end do

CCCCC -----------------------------------------------------------------C
C     DETERMINATION OF THE AEROSOL TYPE AND HUMIDITY CLASS             C
C     !!! AEROSOL TYPE IS LAYER-DEPENDENT (NOT TAKEN INTO ACCOUNT)     C
CCCCC -----------------------------------------------------------------C

      DO 60 IL=1,NL
      IF(RHT(IL).LE.30) THEN
       NH(IL)=1
      ELSE IF(RHT(IL).GT.30.AND.RHT(IL).LE.65) THEN
       NH(IL)=2
      ELSE IF(RHT(IL).GT.65.AND.RHT(IL).LE.75) THEN
       NH(IL)=3
      ELSE IF(RHT(IL).GT.75.AND.RHT(IL).LE.85) THEN
       NH(IL)=4
      ELSE IF(RHT(IL).GT.85.AND.RHT(IL).LE.92) THEN
       NH(IL)=5
      ELSE IF(RHT(IL).GT.92.AND.RHT(IL).LE.97) THEN
       NH(IL)=6
      ELSE IF(RHT(IL).EQ.98) THEN
       NH(IL)=7
      ELSE IF(RHT(IL).EQ.99) THEN
       NH(IL)=8
      END IF
   60 CONTINUE

  999 RETURN
      END

CCCCC *****************************************************************C
       subroutine optcom (ilamb,ihum)
C     *****************************************************************C
C                                                                      C
C     -----------------------------------------------------------------C
C      Reading the raw optical data into a buffer for                  C
C      20 components                                                   C
c                                                                      c
c     In the original Mie calculations, all coefficients and the phase c
c     function are given in [1/cm]. The new calculations give the      c
c     results in [1/m]. Therefore, the new results must be marked in   c
c     the first line with the addition 'new': z.B. TAPE741, neu        c
c                                                                      c
c     13.05.94 Swelling of soot excluded.                              c
c     17.11.97 new file format in ../optdat/                           c
c                                                                      c
c     Stand: 17.11.97                                         M. Hess  c
CCCCC -----------------------------------------------------------------C

      integer prnr,acnr,njc,rht,ilamb,ihum,ibuf,il,ihu,mhum,ic,jc,iht
      integer nta,ntap,nbuf,ib,mbuf,ios,iline,ila,mlamb,it,nl,ntheta
      integer jnopar,nop,niw,nih,nhum,latx,lonx
      real n,numden,alamb,ahum,acmr,ext,extbuf,asy,bacbuf,abs
      real asybuf,bre,brebuf,kbuf,khum,bim,mopar,sisca
      real ncomp,absbuf,bac,nh,pha,bimbuf,natyp
      real mcomp,phabuf,sca,scabuf,sis,sisbuf,thet
      real rlamb,extco,scaco,absco,asymf,exn,refr,refi

      character*1  dum
      character*2  chum
      character*3  atn,pat
      character*4  comnam
      character*8  opanam,optnam
      character*10 dum2
      character*21 tap
      character*20 catyp
      character*30 typnam

      logical exists,ende
      
      dimension jnopar(13),opanam(13),optnam(13),alamb(61),khum(8)
      dimension ahum(8),nhum(8),chum(8),ncomp(10),comnam(20)
      dimension rht(2),n(2),njc(2),acnr(5,2),acmr(5,2),nh(2),atn(2)
      dimension pat(2), ext(2,5),sca(2,5),abs(2,5),sis(2,5),asy(2,5),
     *bac(2,5),pha(2,5,112),bre(2,5),bim(2,5)
      dimension kbuf(20),extbuf(20),scabuf(20),absbuf(20),
     *sisbuf(20),asybuf(20),bacbuf(20),phabuf(112,20),
     *brebuf(20),bimbuf(20)

      common /opar/   mopar,jnopar,nop,opanam,optnam
      common /wavel/  alamb,mlamb,niw
      common /hum/    khum,ahum,nih,nhum,mhum,chum
      common /atyp/   natyp,mcomp,ncomp,numden,catyp,typnam,comnam
      common /mipoi/  latx,lonx,nl,prnr,rht,n,njc,acnr,acmr,nh,atn,pat
      common /oppoi/  ext,sca,abs,sis,asy,bac,pha,bre,bim
      common /buffer/ ibuf,kbuf,extbuf,scabuf,absbuf,
     *sisbuf,asybuf,bacbuf,phabuf,brebuf,bimbuf,mbuf

ccccc -----------------------------------------------------------------c
c      Loop over all components occurring at the grid point            c
ccccc -----------------------------------------------------------------c

      do il=1,nl

       if (nih.eq.0) then
        do ihu=1,mhum
        if (nh(il).eq.ihu) then
         khum(ihum)=ihu
        end if
        end do
       end if
      
       do ic=1,njc(il)
        jc=acnr(ic,il)
ccccc -----------------------------------------------------------------c
c     Exclusion? swelling at insoluble, soot and                       c
c     mineral components and at clouds                                 c
ccccc -----------------------------------------------------------------c
        if ( jc.eq.1.or.jc.eq.3.or.(jc.ge.6.and.jc.le.9).or.
     *   jc.gt.10 ) then
         iht=1
         nta=700+(jc*10)+1
        else
c        iht=khum(ihum)
         iht=ihum
         nta=700+(jc*10)+iht
        end if
ccccc -----------------------------------------------------------------c
c     Determination of the file name of the sought component 	       c
c     from component number and moisture class                              c
ccccc -----------------------------------------------------------------c
        tap(1:15)='extdata/optdat/'
        tap(16:19)=comnam(jc)
        tap(20:21)=chum(iht)
        ntap=20

ccccc -----------------------------------------------------------------c
c     Determination of the identification number for the buffer over   c
c     the wavelength                                                   c
c                                                                      c
c         nbuf: ID number of the current component for the buffer      c
c     kbuf(20): ID numbers of the stored components                    c
c         mbuf: Position of the current component in the buffer        c
c         ibuf: Position up to which the buffer is occupied            c
ccccc -----------------------------------------------------------------c
        nbuf=ilamb*1000+nta
c 	  print*,'nbuf= ',nbuf
ccccc -----------------------------------------------------------------c
c      Check the buffer for compliance with nbuf                       c
ccccc -----------------------------------------------------------------c
        exists=.false.
         do ib=1,20
          if (nbuf.eq.kbuf(ib)) then
           exists=.true.
           mbuf=ib
           goto 10
          end if
         end do
   10   continue

ccccc -----------------------------------------------------------------c
c      Reading the component data if it is not in the buffer,          c
c      otherwise transferring it from the buffer                       c
ccccc -----------------------------------------------------------------c
        if (exists) then
         ext(il,ic)=extbuf(mbuf)
         sca(il,ic)=scabuf(mbuf)
         abs(il,ic)=absbuf(mbuf)
         sis(il,ic)=sisbuf(mbuf)
         asy(il,ic)=asybuf(mbuf)
         bac(il,ic)=bacbuf(mbuf)
         bre(il,ic)=brebuf(mbuf)
         bim(il,ic)=bimbuf(mbuf)
        else
        if (ibuf.lt.20) then
         ibuf=ibuf+1
        else
         ibuf=1
        end if

        kbuf(ibuf)=nbuf
        open (ntap,file=tap,iostat=ios)
        if (ios.ne.0) then
         print*,' error while opening file ',tap
         print*,'  latitude: ',latx
         print*,' longitude: ',lonx
         stop
        end if
        do iline=1,100
         read (ntap,220) dum2
         if (dum2.eq.'# optical ') then
          goto 2002
         end if
        end do
 2002   continue
        do iline=1,5
         read (ntap,200) dum
        end do

        do ila=1,mlamb
         read (ntap,500) rlamb,extco,scaco,absco,sisca,asymf,
     *    exn,refr,refi
  500    format(2x,7e10.3,2e11.3)

         if (rlamb.eq.alamb(ilamb)) then
          ext(il,ic)=extco
          sca(il,ic)=scaco
          abs(il,ic)=absco
          sis(il,ic)=sisca
          asy(il,ic)=asymf
          bre(il,ic)=refr
          bim(il,ic)=refi
         end if
        end do
        read (ntap,'(7(/))')
        it=1
        ende=.false.
        do while (.not.ende)
         read (ntap,510,end=511)
     *   thet,(pha(il,ic,min(112,it)),ila=1,ilamb)
  510    format(e11.3,1x,70e10.3)
         it=it+1
        end do
  511   ntheta=it-1

c ENDE NEUER INPUT

        bac(il,ic)=pha(il,ic,min(112,ntheta))

        extbuf(ibuf)=ext(il,ic)
        scabuf(ibuf)=sca(il,ic)
        absbuf(ibuf)=abs(il,ic)
        sisbuf(ibuf)=sis(il,ic)
        asybuf(ibuf)=asy(il,ic)
        brebuf(ibuf)=bre(il,ic)
        bimbuf(ibuf)=bim(il,ic)
        bacbuf(ibuf)=bac(il,ic)

        close (ntap)

        end if
       end do
      end do

ccccc -----------------------------------------------------------------c
c     Formate							       c
ccccc -----------------------------------------------------------------c

  200 format(a1)
  220 format(a10)
      return
      end

CCCCC *****************************************************************C
      subroutine optpar (ilamb,ihum)
C     *****************************************************************C
C                                                                      c
C     -----------------------------------------------------------------C
C     Calculation and printout of the desired optical parameters       C
CCCCC -----------------------------------------------------------------C

      integer prnr,ilamb,ihum,it,l,jc,nl,NJC,itp,kc,iop,kop,il,nlay
      integer nil,nprog,rht,acnr,jnopar,nop,mlamb,niw,norm,mixnor
      integer latx,lonx,angle,nia,jnangle,season
      real numden
      REAL EXTN,ABSN,SCAN,PF18N,supf,phafu
      REAL EXTA,ABSA,SCAA,SSA,ASF,PF18A
      real scar,absr,omer
      
      real mlay,nltyp,parlay,boundl,boundu,acmr
      real hfta,hstra,h0,h1,hm
      real mopar
      real alamb
      real n
      real nh
      real ext,sca,abs,sis,asy
      real bac,pha,bre,bim
      real EXTFTA,EXTSTR
      real ITEST,mcomp
      real oparam,phaf,hftae,odepth
      real hu,ho,z,odeptha,turbr,turbf
      real ncomp
      real smas,smag,natyp
      real lat5,lon5,relhum,optdep
      real SUMME,SUMMA,SUMMS,SUMSSA,SUMASF,SUPF18
      
      
      character*3  atn,pat
      character*4  comnam
      character*8  opanam,optnam
      character*20 catyp
      character*30 typnam
      
      dimension optdep(25,2),EXTN(2),ABSN(2),SCAN(2),PF18N(2),supf(112)
      dimension phafu(112,2),EXTA(2),ABSA(2),SCAA(2),SSA(2),ASF(2)
      dimension PF18A(2),scar(2),absr(2),omer(2)
      dimension nltyp(10),parlay(10,2),boundl(10),boundu(10)
      dimension nil(10),hfta(10),hstra(10),h0(2,10),h1(2,10),hm(2,10)
      dimension jnopar(13),opanam(13),optnam(13)
      dimension alamb(61)
      dimension rht(2),n(2)
      dimension njc(2),acnr(5,2),acmr(5,2),nh(2),atn(2),pat(2)
      dimension ext(2,5),sca(2,5),abs(2,5),sis(2,5),asy(2,5)
      dimension bac(2,5),pha(2,5,112),bre(2,5),bim(2,5)
      dimension EXTFTA(61),EXTSTR(61)
      dimension ITEST(2)
      dimension oparam(10,2),phaf(112,2)
      dimension jnangle(112),angle(112)
      dimension smas(10,8),smag(8),ncomp(10),comnam(20)

      
      common /atyp/   natyp,mcomp,ncomp,numden,catyp,typnam,comnam
      common /norm/   norm,mixnor
      common /layer/  nlay,mlay,nltyp,parlay,boundl,boundu
      common /profi/ nil,hfta,hstra,h0,h1,hm
      common /opar/   mopar,jnopar,nop,opanam,optnam
      common /wavel/  alamb,mlamb,niw
      common /mipoi/  latx,lonx,nl,prnr,rht,n,njc,acnr,acmr,nh,atn,pat
      common /oppoi/  ext,sca,abs,sis,asy,bac,pha,bre,bim
      COMMON /FTASTR/ EXTFTA,EXTSTR
      COMMON /TEST/   ITEST
      common /out/    oparam,phaf
      common /prog/   nprog
      common /angle/  jnangle,angle,nia
      common /masse/  smas,smag
      common /r/      lat5,lon5,relhum,season,optdep

CCCCC ------------------------------------------------------------------C
C     MIXING THE AEROSOL TYPE                                           C
C     SUMM(E,A,S) : TOTAL EXTINCTION, ABSORPTION, SCATTERING            C
C     SUPF18      : TOTAL REVERSE COEFFICIENT                           C
C     SUMASF      : SUBTOTAL ASYMMETRY FACTOR (ASF)                     C
C     SUMASF      : SUBTOTAL OF THE SINGLE SCAT. ALB. (SSA)             C
CCCCC ------------------------------------------------------------------C

      DO 10 L=1,NL

       SUMME = 0.
       SUMMA = 0.
       SUMMS = 0.
       SUMSSA = 0.
       SUMASF = 0.
       SUPF18 = 0.
       if (jnopar(10).eq.1) then
        do it=1,112
         supf(it)=0.
        end do
       end if

       DO 20 JC=1,NJC(L)

c      Calculation of the sums'

        SUMME = SUMME + ACMR(JC,L)*EXT(l,jc)
        SUMMA = SUMMA + ACMR(JC,L)*ABS(l,jc)
        SUMMS = SUMMS + ACMR(JC,L)*SCA(l,jc)
        SUMSSA = SUMSSA + ACMR(JC,L)*sis(l,jc)
     +   *EXT(l,jc)
        SUMASF = SUMASF + ACMR(JC,L)*asy(l,jc)
     +   *SCA(l,jc)
        SUPF18 = SUPF18 + ACMR(JC,L)*bac(l,jc)
        if (jnopar(10).eq.1) then
         do it=1,112
          supf(it)=supf(it)+acmr(jc,l)*pha(l,jc,it)
         end do
        end if

   20  CONTINUE

CCCCC -----------------------------------------------------------------C
C     Standardized optical parameters				       c
CCCCC -----------------------------------------------------------------C

c      Calculation of the standardized values'

       EXTN(L) = SUMME
       ABSN(L) = SUMMA
       SCAN(L) = SUMMS
       PF18N(L) = SUPF18
       if (jnopar(10).eq.1) then
        do it=1,112
         phafu(it,l)=supf(it)
        end do
       end if

       SSA(L) = SUMSSA/SUMME
       ASF(L) = SUMASF/SUMMS

CCCCC -----------------------------------------------------------------C
C     ABSOLUTE OPTICAL PARAMETERS                                      C
CCCCC -----------------------------------------------------------------C

c      print*,' Calculation of the absolute values'

       EXTA(L)= EXTN(L) * N(L)
       ABSA(L)= ABSN(L) * N(L)
       SCAA(L)= SCAN(L) * N(L)
       PF18A(L) = PF18N(L)* N(L)
       if (jnopar(10).eq.1.and.norm.eq.1) then
        do it=1,112
         phafu(it,l)=phafu(it,l)*n(l)
        end do
       end if
       if (norm.eq.1) then
        EXTN(L)= EXTA(L)
        ABSN(L)= ABSA(L)
        SCAN(L)= SCAA(L)
        PF18N(L) = PF18A(L)
       end if

       if (jnopar(10).eq.1) then
        itp=1
        do it=1,112
         if (jnangle(it).eq.1) then
          phaf(itp,l)=phafu(it,l)
          itp=itp+1
         end if
        end do
       end if

       if (jnopar(11).eq.1) then
        scar(l)=scaa(l)/smag(ihum)*1000.  ! unit m**2/g
       end if

       if (jnopar(12).eq.1) then
        absr(l)=absa(l)/smag(ihum)*1000.
       end if

       if (jnopar(13).eq.1) then
        kc=0
        do jc=1,njc(l)
         if (ncomp(jc).eq.3) kc=jc
        end do
        if (kc.ne.0) then
         omer(l)=ssa(l)/smas(kc,ihum)
        else
         omer(l)=99.
        end if
       end if

CCCCC -----------------------------------------------------------------C
C     OUTPUT OF DATA                    						       C
CCCCC -----------------------------------------------------------------C

       iop=0
       kop=0
       if (jnopar(1).eq.1) then
        iop=iop+1
        oparam(iop,l)=extn(l)
       end if
       if (jnopar(2).eq.1) then
        iop=iop+1
        oparam(iop,l)=scan(l)
       end if
       if (jnopar(3).eq.1) then
        iop=iop+1
        oparam(iop,l)=absn(l)
       end if
       if (jnopar(4).eq.1) then
        iop=iop+1
        oparam(iop,l)=ssa(l)
       end if
       if (jnopar(5).eq.1) then
        iop=iop+1
        oparam(iop,l)=asf(l)
       end if
       if (jnopar(9).eq.1) then
        iop=iop+1
        if (jnopar(6).eq.1) kop=kop+1
         if (jnopar(7).eq.1) kop=kop+1
          if (jnopar(8).eq.1) kop=kop+1
           kop=kop+iop
           oparam(kop,l)=exta(l)/pf18a(l)
          end if
          if (jnopar(11).eq.1) then
           iop=iop+1
           oparam(iop,l)=scar(l)
          end if
          if (jnopar(12).eq.1) then
           iop=iop+1
           oparam(iop,l)=absr(l)
          end if
          if (jnopar(13).eq.1) then
           iop=iop+1
           oparam(iop,l)=omer(l)
          end if
   10 CONTINUE

CCCCC -----------------------------------------------------------------C
C     OPTICAL THICKNESS                                                C
CCCCC -----------------------------------------------------------------C

      if (jnopar(6).eq.1) then

CCCCC -----------------------------------------------------------------C
c     Determination of HM, HFTA, HSTR, EXTFTA, EXTSTR aus den	       c
c     read values ??in / layer / for RAWOPT		                       c
CCCCC -----------------------------------------------------------------C

       if (nprog.eq.2) then
        do il=1,nl
         if (nltyp(il).eq.1) then
          hm(il,1)=parlay(il,1)
         else if (nltyp(il).eq.2) then
          hm(il,1) = parlay(il,2)*
     *     (exp(-boundl(il)/parlay(il,2))+
     *     exp(-boundu(il)/parlay(il,2)))
         end if
        end do
        hfta(1)=boundu(nlay-2)-boundl(nlay-2)
        hstra(1)=boundu(nlay-1)-boundl(nlay-1)
        extfta(ilamb)=parlay(nlay-2,1)
        extstr(ilamb)=parlay(nlay-1,1)
       end if

       hu = h1(nl,prnr)
       ho = hu + hfta(prnr)
       z = 8.
       hftae = z * ( exp(-hu/z) - exp(-ho/z) )

       ODEPTH = 0.

       DO IL=1,nl
        ODEPTH = ODEPTH + EXTA(IL) * HM(IL,prnr)
       end do

CCCCC -----------------------------------------------------------------C
C                     + FREE TROP. AEROSOL                             C
CCCCC -----------------------------------------------------------------C

       ODEPTH = ODEPTH + EXTFTA(ilamb)*HFTAE
     +  + EXTSTR(ilamb)*HSTRA(prnr)

       odeptha=odepth/alog(10.)

       turbr=0.008569*alamb(ilamb)**(-4)*(1.+0.0113*alamb(ilamb)**
     *         (-2)+0.00013*alamb(ilamb)**(-4))

       turbf=(odepth+turbr)/turbr

       if (jnopar(9).eq.1) then
        kop=iop
       else
        kop=iop+1
c       print*,' exta= ', exta(il),' hm= ',hm(il,1)
       end if

       if (jnopar(6).eq.1) then
        oparam(kop,1)=odepth
        oparam(kop,2)=0.
        kop=kop+1
        iop=iop+1
       end if

       if (jnopar(7).eq.1) then
        oparam(kop,1)=odeptha
        oparam(kop,2)=0.
        kop=kop+1
        iop=iop+1
       end if

       if (jnopar(8).eq.1) then
        oparam(kop,1)=turbf
        oparam(kop,2)=0.
        iop=iop+1
       end if
      end if

      call out4(iop,ilamb,ihum)

      RETURN
      END

ccccc *****************************************************************c
      subroutine out4(iop,ilamb,ihum)
c     *****************************************************************c
c			                                                     c
c     -----------------------------------------------------------------c
c     OUTPUT OF DATA for atlopt          			                 c
ccccc -----------------------------------------------------------------c

      integer prnr,acnr,njc,rht,l,nl,nop,ihum,ilamb,iop,jnopar,mlamb
      integer latx,lonx,angle,nia,jnangle,niw,nih,nhum,mhum,season
      real alamb
      real khum,ahum
      real oparam,phaf
      real mopar
      real n,acmr,nh
      real lat5,lon5,relhum,optdep      
      character*2 chum
      character opanam*8,atn*3,pat*3,optnam*8      
      
      dimension optdep(25,2),jnangle(112),angle(112)
      dimension alamb(61)
      dimension khum(8),ahum(8),nhum(8),chum(8)
      dimension oparam(10,2),phaf(112,2)
      dimension jnopar(13),opanam(13),optnam(13)
      dimension rht(2),n(2),
     * njc(2),acnr(5,2),acmr(5,2),nh(2),atn(2),pat(2)

      common /angle/  jnangle,angle,nia
      common /wavel/  alamb,mlamb,niw
      common /hum/    khum,ahum,nih,nhum,mhum,chum
      common /out/    oparam,phaf
      common /opar/   mopar,jnopar,nop,opanam,optnam
      common /mipoi/  latx,lonx,nl,prnr,rht,n,
     * njc,acnr,acmr,nh,atn,pat
      common /r/      lat5,lon5,relhum,season,optdep

      do l=1,nl
       if (nih.ne.0) rht(l)=int(ahum(ihum))
        if (nop.gt.iop) then
         if (jnopar(9).eq.1.and.jnopar(10).eq.1) then
          oparam(iop,l)=oparam(nop-1,l)
         else if (jnopar(9).eq.1.and.jnopar(10).eq.0) then
          oparam(iop,l)=oparam(nop,l)
         end if
        end if

        if (iop.le.10) then
         if (l.eq.1) then
          optdep(ilamb,1)=alamb(ilamb)
          optdep(ilamb,2)=oparam(6,l)
         end if
        end if

      end do

      return
      end