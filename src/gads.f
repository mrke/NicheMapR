      subroutine gads(lat51,lon51,relhum1,season1,optdep1)

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
c      Berechnung der optischen Aerosoldaten aus den mikrophysika-     c
c      lischen Rohdaten.                                               c
c                                                                      c
c      ATLOPT ist eine modifizierte Version des Programms OPAC und     c
c      steuert den Aufruf der Unterprogramme:                          c
c                                                                      c
c      - HEAD4                                                         c
c      - LOCATE                                                        c
c      - D4RAW                                                         c
c      - OPTCOM                                                        c
c      - OPTPAR                                                        c
c      - OUT4                                                          c
c                                                                      c
c      Diese Unterprogramme sind an ATLOPT angeh?ngt.                  c
c                                                                      c
c      ab 14.12.92 anderes Format der neuen Mie-Rechnungen eingebaut   c
c      ab 04.11.93 neue Parameter scattering, absorption, omega ratio  c
c      ab 13.05.94 OPTCOM: Russ quillt nicht mehr                      c
c      ab 27.07.94 opt. Dicke fuer alle Wellenlaengen berechnet        c
c                                                                      c
c      18.11.97 GADS 2.1                                               c
c                                                                      c
c      18.11.97                                                M. Hess c
ccccc -----------------------------------------------------------------c

      integer   prnr,acnr,njc,rht
      real      n,numden
	real	  lat5,lon5,relhum,optdep(25,2),season
	double precision	  lat51,lon51,relhum1,optdep1(25,2),season1

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
      CHARACTER(len=255) :: cwd

      common /prog/   nprog
      common /profi/  nil(10),hfta(10),hstra(10),
     *                h0(2,10),h1(2,10),hm(2,10)
      common /buffer/ ibuf,kbuf(20),extbuf(20),scabuf(20),absbuf(20),
     *                sisbuf(20),asybuf(20),bacbuf(20),phabuf(112,20),
     *                brebuf(20),bimbuf(20),mbuf
      common /mipoi/  latx,lonx,nl,prnr,rht(2),n(2),
     *                njc(2),acnr(5,2),acmr(5,2),nh(2),atn(2),pat(2)
      common /oppoi/  ext(2,5),sca(2,5),abs(2,5),sis(2,5),asy(2,5),
     *                bac(2,5),pha(2,5,112),bre(2,5),bim(2,5)

      common /atyp/   natyp,mcomp,ncomp(10),numden,
     *                catyp,typnam,comnam(20)
      common /layer/  nlay,mlay,nltyp(10),parlay(10,2),boundl(10),
     *                boundu(10)
      common /param/  nptyp,mpar,par(5)
      common /season/ nseas,cseas,tseas
      common /opar/   mopar,jnopar(13),nop,opanam(13),optnam(13)
      common /wavel/  mlamb,alamb(61),niw
      common /hum/    khum(8),ahum(8),nih,nhum(8),mhum,chum(8)
      common /geog/   lata,late,lati,lona,lone,loni,na,area
      common /norm/   norm,mixnor
	common /r/      lat5,lon5,relhum,season,optdep

      data alamb /0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,
     *            0.9,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.2,3.39,3.5,3.75,
     *            4.0,4.5,5.0,5.5,6.0,6.2,6.5,7.2,7.9,8.2,8.5,8.7,9.0,
     *            9.2,9.5,9.8,10.0,10.6,11.0,11.5,12.5,13.0,14.0,14.8,
     *            15.0,16.4,17.2,18.0,18.5,20.0,21.3,22.5,25.0,27.9,30.,
     *            35.0,40.0/,mlamb/61/

      data ahum /0.,50.,70.,80.,90.,95.,98.,99./
      data nhum /0,50,70,80,90,95,98,99/,mhum/8/
      data chum /'00','50','70','80','90','95','98','99'/

      data comnam /'inso','waso','soot','ssam','sscm','minm','miam',
     *             'micm','mitr','suso','stco','stma','cucc','cucp',
     *             'cuma','fog-','cir1','cir2','cir3','    '/

      data optnam /'ext.coef','sca.coef','abs.coef','sisc.alb',
     *             'asym.par','op.depth',
     *             '        ','turb.fac','li.ratio','pha.func',
     *             'ext.rat ','abs.rat ',
     *             '        '/

      data jnopar/1,1,1,1,1,1,0,0,1,0,0,0,0/,nop/7/

CCCCC -----------------------------------------------------------------C
c     some definitions for this version                                c
CCCCC -----------------------------------------------------------------C
            !CALL getcwd(cwd)
            !WRITE(*,*) TRIM(cwd)
	lat5=lat51
	lon5=lon51
	season=season1
	relhum=relhum1
      niw=1
      njh=1
      nih=1
      lata=lat5
      late=lat5
      lati=5
      lona=lon5
      lone=lon5
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
c     Abfrage, was geplottet werden soll                               c
ccccc -----------------------------------------------------------------c

 1001 continue
c print*,'  '
c      write(*,154)
c  154 format(' (w)inter or (s)ummer? ')
c      read (*,'(a)') ws
	if(season.eq.1)then
	ws='w'
	else
	ws='s'
	endif


ccccc ----------------------------------------------------------------c
c     Input: wavelength                                               c
ccccc ----------------------------------------------------------------c

c      print*,' please select wavelength: '
c      print*,' '

      nwel=25
c      do 11 iwel=1,22
c         if (nwel.ge.(iwel+44)) then
c            write(*,114) iwel,alamb(iwel),(iwel+22),alamb(iwel+22),
c     *                   (iwel+44),alamb(iwel+44)
c         else if (nwel.ge.(iwel+22)) then
c            write(*,113) iwel,alamb(iwel),(iwel+22),alamb(iwel+22)
c         else
c            write(*,111) iwel,alamb(iwel)
c         end if
c   11 continue

  111 format(5x,'(',i2,')',3x,f5.2,1x,'um')
  113 format(5x,'(',i2,')',3x,f5.2,1x,'um',5x,'(',i2,')',
     *       3x,f5.2,1x,'um')
  114 format(5x,'(',i2,')',3x,f5.2,1x,'um',5x,'(',i2,')',
     *       3x,f5.2,1x,'um',
     *       5x,'(',i2,')',3x,f5.2,1x,'um')

  909 continue
c      write (*,*) '?'
c      read (*,*) iwel
	iwel=1
	do 9999 iwel=1,nwel

      if (ws.eq.'w') then
         open(ntape,file='extdata/glodat/winter.dat')
         read (ntape,'(a1)') dum
         cseas='winter '
      else if (ws.eq.'s') then
         open(ntape,file='extdata/glodat/summer.dat')
         read (ntape,'(a1)') dum
         cseas='summer '
      else
         print*,' wrong input! try again!'
         goto 1001
      end if

      if (iwel.lt.1.or.iwel.gt.nwel) then
         print*,' wrong number! try again! '
         goto 909
      end if

ccccc ----------------------------------------------------------------c
c     Input: humidity                                                 c
ccccc ----------------------------------------------------------------c

c      print*,' please select rel. humidity: '
c      print*,' '
c      do ihum=1,mhum
c         write(*,121) ihum,nhum(ihum)
c  121    format (5x,'(',i2,')',3x,i2,' %')
c      end do

c  908 write (*,*) '?'
c      read (*,*) ihum
	ihum=relhum
c      if (ihum.lt.1.or.ihum.gt.mhum) then
c         print*,' wrong number! try again! '
c         goto 908
c      end if

CCCCC -----------------------------------------------------------------C
C     EINLESEN DER HOEHEN-PROFILE vom File TAPE9		       c
C								       C
C     HM    : EFFEKTIVE SCHICHTDICKE (HOMOGENE VERTEILUNG)             C
C     HFTA  : SCHICHTDICKE DES FREIEN TROP. AEROSOLS IN KM             C
C     HSTRA : SCHICHTDICKE DES STRATOSPH. AEROSOLS IN KM               C
C     NIL   : ANZAHL DER SCHICHTEN                                     C
CCCCC -----------------------------------------------------------------C

c      print*,' Anfang prof'
       call prof
c      print*,' Ende prof'

ccccc -----------------------------------------------------------------c
c      Beschriftung des Output-Files                                   c
ccccc -----------------------------------------------------------------c

c      print*,' Anfang head4'
	 if(iwel.eq.1)then
       call head4 (iwel,ihum)
	 endif
c      print*,' Ende head4'

ccccc -----------------------------------------------------------------c
c      Schleife ?ber alle verlangten Wellenl?ngen und Feuchteklassen   c
ccccc -----------------------------------------------------------------c

       do il=1,niw

          do ih=1,njh

ccccc -----------------------------------------------------------------c
c      Schleife ?ber alle verlangten geographischen Koordinaten        c
ccccc -----------------------------------------------------------------c

             do ilat=lata,late,-lati

                do ilmal=1,nlmal
                   do ilon=lona,lone,loni

ccccc -----------------------------------------------------------------c
c      Einlesen der Rohdaten von den Files TAPE201, TAPE207:	       c
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

c       print*,' Anfang d4raw'
                      call d4raw (ilat,ilon,ntape)
c       print*,' Ende d4raw'

ccccc -----------------------------------------------------------------c
c      Einlesen der optischen Rohdaten von den Files winter.dat and    c
c      summer.dat                                                      c
ccccc -----------------------------------------------------------------c

c       print*,' Anfang optcom'
                      call optcom (iwel,ihum)
c       print*,' ende optcom'

ccccc -----------------------------------------------------------------c
c      Berechnung der optischen Parameter am aktuellen Gitterpunkt     c
ccccc -----------------------------------------------------------------c

c       print*,' Anfang optpar',ilat,ilon
                      call optpar(iwel,ihum)
c       print*,' Ende optpar',ilat,ilon

                   end do
                end do
             end do
          end do
       end do
       close (ntape)
9999	continue

c       close (10)

	 do 9998 i=1,25
	  do 9997 j=1,2
	   optdep1(i,j)=optdep(i,j)
9997	  continue
9998	 continue

       return
       end

CCCCC *****************************************************************C
      SUBROUTINE PROF
C     *****************************************************************C
C								       C
C     -----------------------------------------------------------------C
C     EINLESEN DER HOEHEN-PROFILE vom File profiles.dat und der        C
C     Extinktionskoeffizienten der oberen Atmosph?re von extcof.dat    C
CCCCC -----------------------------------------------------------------C

      common /wavel/  mlamb,alamb(61),niw
      common /profi/ nil(10),hfta(10),hstra(10),
     *		      h0(2,10),h1(2,10),hm(2,10)
      COMMON /FTASTR/ EXTFTA(61),EXTSTR(61)

CCCCC -----------------------------------------------------------------C
C     ES GIBT 7 PROFIL-TYPEN. Folgende Daten werden eingelesen:        C
c                                                                      c
c           iip: Nummer des Profiltyps                                 c
c       nil(ip): Zahl der Schichten fuer Typ ip (wie in tape201 usw.)  c
c      hfta(ip): Hoehe der Schicht fuer das freie troposph. Aerosol    c
c     hstra(ip): Hoehe der Schicht des stratosphaerischen Aerosols     c
c     h0(il,ip): Untergrenze der Schicht il                            c
c     h1(il,ip): Obergrenze der Schicht                                c
c     hm(il,ip): effektive Dicke der Schicht il fuer den Typ ip        c
CCCCC -----------------------------------------------------------------C

      open (8,file='extdata/profiles.dat')

      nprof=7
      DO IP=1,nprof
         READ(8,8010) IIP,NIL(IP),HFTA(IP),HSTRA(IP)
         READ(8,8020) (H0(IL,IP),H1(IL,IP), HM(IL,IP),IL=1,NIL(IP) )
 8010    FORMAT(I3,I3,2F8.2)
 8020    FORMAT(2F5.1,F10.3)
      end do

      close (8)

CCCCC -----------------------------------------------------------------C
C     Einlesen der Extinktionskoeffizienten f?r die obere Atmosph?re:  C
C								       C
C     EXTINCTION COEFFICIENT -	FREE TROPOSPHERIC AEROSOL  +	       C
C     EXTINCTION COEFFICIENT -	  STRATOSPHERIC   AEROSOL	       C
C								       C
C     UEBERSPRINGEN DER ERSTEN BEIDEN ZEILEN von TAPE9		       C
CCCCC -----------------------------------------------------------------C

      open (9,file='extdata/extback.dat')
      IL=1
      READ(9,'(/)')
      DO IWL=1,mlamb
         READ(9,*) WAVE,EXTFT,EXTST
c        do ila=1,niw
c           IF (WAVE.EQ.alamb(ila)) THEN
               EXTFTA(IWL)=EXTFT
               EXTSTR(IWL)=EXTST
c              IL=IL+1
c           END IF
c        end do
      end do

      close (9)

      RETURN
      END

ccccc *****************************************************************c
      subroutine head4 (il,ih)
c     *****************************************************************c
c                                                                      c
c     -----------------------------------------------------------------c
c     Beschriftung des Output-Files TAPE10			           c
ccccc -----------------------------------------------------------------c

      real mixrat,numden

      character*2 chum
      character*4 comnam
      character opanam*8,cseas*7,tseas*11,optnam*8,opnam(10)*8
      character catyp*20,area*50,typnam*30

      common /season/ nseas,cseas,tseas
      common /geog/   lata,late,lati,lona,lone,loni,na,area
      common /hum/    khum(8),ahum(8),nih,nhum(8),mhum,chum(8)
      common /norm/   norm,mixnor
      common /numdis/ sigma(10),rmin(10,8),rmax(10,8),rmod(10,8),
     *                mixrat(10),dens(10,8)
      common /atyp/   natyp,mcomp,ncomp(10),numden,
     *                catyp,typnam,comnam(20)
      common /opar/   mopar,jnopar(13),nop,opanam(13),optnam(13)
      common /wavel/  mlamb,alamb(61),niw
      common /angle/  jnangle(112),angle(112),nia

CCCCC -----------------------------------------------------------------C
C     Kopf des output-files                                            C
CCCCC -----------------------------------------------------------------C

c      if (ih.eq.1.and.il.eq.1) then
c         open (10,file='extdata/aererg.txt')

c         write(10,100) cseas
c  100	   format('# Global Aerosol Data Set, Version 2.2a'/,
c     *          '#'/
c     *          '# ',a7,/
c     *          '#')
c  107	 format('====================================================',
c     *		'============================')
c      end if

CCCCC -----------------------------------------------------------------C
C     Beschriftung f?r verschiedene Wellenlaengen und rel. Feuchten    c	       c
CCCCC -----------------------------------------------------------------C

c      if (nih.ne.0) then
c         WRITE(10,4000) alamb(il),ahum(ih)
c 4000    FORMAT('#'/
c     *          '# wavelength: ',f6.3,3x,'relative humidity: ',F3.0,'%'/
c     *          '#')
c      else
c         WRITE(10,4004) alamb(il)
c 4004    FORMAT(/' wavelength: ',f6.3,3x,'relative humidity: -- %')
c      end if

      if (jnopar(10).eq.1) then
         kop=nop-1
      else
         kop=nop
      end if

      do iop=1,kop
         opnam(iop)=opanam(iop)
      end do

c      if (kop.le.10) then
c         write(10,4001) (opnam(in),in=1,kop)
c 4001    format('   LAT  LON NL LAMB    RELHUM ',10(1x,a8,1x))
c      else
c         write(10,4001) (opnam(in),in=1,10)
c         write(10,4002) (opnam(in),in=11,kop)
c 4002    format('                               ',5(1x,a8,1x))
c      end if

c      write(10,4003)
c 4003 format('#',13x,'  [1/km]  ','  [1/km]  ','  [1/km]  ',
c     *       30x,'   [sr]')

      return
      end

CCCCC *****************************************************************C
      subroutine d4raw (lat,lon,ntape)
C     *****************************************************************C
C                                                                      C
C     -----------------------------------------------------------------C
C     EINLESEN DER DATEN VON DEN ROHDATEN-FILES TAPE201-TAPE212        C
CCCCC -----------------------------------------------------------------C

      IMPLICIT CHARACTER*3 (Z)

      integer prnr,acnr,njc,rht
      real n
      character*3 atn,pat

      common /prog/  nprog
      common /mipoi/ latx,lonx,nl,prnr,rht(2),n(2),
     *               njc(2),acnr(5,2),acmr(5,2),nh(2),atn(2),pat(2)
      common /test/  itest(2)
      common /dat/   zat(11),zrh(8)

  111 READ(NTAPE,1020,end=999) LATX,LONX,NL,PRNR,
     +   ATN(1),PAT(1),RHT(1),N(1),NJC(1),
     +      ( ACNR(JC,1),ACMR(JC,1),JC=1,3 )
 1020 FORMAT(2I4,2I3,(2A3,I3,E10.3,I4,3(I3,E10.4)))

c      write(*,1020) LATX,LONX,NL,PRNR,
c     +   ATN(1),PAT(1),RHT(1),N(1),NJC(1),
c     +      ( ACNR(JC,1),ACMR(JC,1),JC=1,3 )

      IF(LATX.NE.LAT .OR.  LONX.NE.LON) THEN
c         print*,' Achtung, falsche Koordinaten: ',latx,lonx
         GOTO 111
      END IF

      IF (NJC(1).GT.3) THEN
          READ(NTAPE,1025,end=999) ( ACNR(JC,1),ACMR(JC,1),JC=4,NJC(1))
 1025     FORMAT(37X,3(I3,E10.4))
c          write(*,1025) ( ACNR(JC,1),ACMR(JC,1),JC=4,NJC(1))
      END IF

      IF(NL.NE.1) THEN
         DO 10 L=2,NL
         READ(NTAPE,1021,end=999)  ATN(L),PAT(L),RHT(L),N(L),NJC(L),
     +                   ( ACNR(JC,L),ACMR(JC,L),JC=1,3 )
 1021    FORMAT(14X,2A3,I3,E10.3,I4,5(I3,E10.4))

         IF (NJC(L).GT.3) THEN
           READ(NTAPE,1025,end=999) ( ACNR(JC,L),ACMR(JC,L),JC=4,NJC(L))
         END IF

   10    CONTINUE
      END IF

      do l=1,nl
         sum=0.
         do ic=1,njc(l)
            sum=sum+acmr(ic,l)
         end do
         if (abs(sum-1.).ge.0.01) then
            print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
            print*,'***!!!    sum of mixing ratios is not 1.     !!!***'
            print*,'***!!! please have a look at errorfile *.err !!!***'
            print*,'***************************************************'
            write (2,1001) latx,lonx,sum
         end if
      end do
1001  format (2i4,3x,1pe10.3)

CCCCC -----------------------------------------------------------------C
C     BESTIMMUNG DER NUMMER DES AEROSOLTYPS UND DER FEUCHTEKLASSE      C
C     !!! AEROSOLTYP IST SCHICHTABHAENGIG (NOCH NICHT BERUECKSICHTIGT) C
CCCCC -----------------------------------------------------------------C

c     DO 50 IT=1,11
c     IF(ATN(1).EQ.ZAT(IT)) THEN
c	 NT=IT
c     END IF
c  50 CONTINUE

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
C      Einlesen der optischen Rohdaten in einen Puffer f?r             C
C      20 Komponenten                                                  C
c                                                                      c
c     Bei den urspruenglichen Mie-Rechnungen sind alle Koeffizienten   c
c     und die Phasenfunktion in [1/cm] angegeben. Die neuen Rechnungen c
c     geben die Ergebnisse in [1/m]. Daher muessen die neuen           c
c     Ergebnisse in der ersten Zeile durch den Zusatz 'neu' gekenn-    c
c     zeichnet werden: z.B. TAPE741, neu                               c
c                                                                      c
c     13.05.94 Quellung von Russ ausgeschlossen.                       c
c     17.11.97 new file format in ../optdat/                           c
c                                                                      c
c     Stand: 17.11.97                                         M. Hess  c
CCCCC -----------------------------------------------------------------C

      integer prnr,acnr,njc,rht
      real n,numden

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

      common /opar/   mopar,jnopar(13),nop,opanam(13),optnam(13)
      common /wavel/  mlamb,alamb(61),niw
      common /hum/    khum(8),ahum(8),nih,nhum(8),mhum,chum(8)
      common /atyp/   natyp,mcomp,ncomp(10),numden,
     *                catyp,typnam,comnam(20)
      common /mipoi/  latx,lonx,nl,prnr,rht(2),n(2),
     *                njc(2),acnr(5,2),acmr(5,2),nh(2),atn(2),pat(2)
      common /oppoi/  ext(2,5),sca(2,5),abs(2,5),sis(2,5),asy(2,5),
     *		    bac(2,5),pha(2,5,112),bre(2,5),bim(2,5)
      common /buffer/ ibuf,kbuf(20),extbuf(20),scabuf(20),absbuf(20),
     *                sisbuf(20),asybuf(20),bacbuf(20),phabuf(112,20),
     *	          brebuf(20),bimbuf(20),mbuf

ccccc -----------------------------------------------------------------c
c      Schleife ?ber alle am Gitterpunkt vorkommenden Komponenten      c
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

c	  print*,'Anfang Komponenten schleife: ',ic,njc(il)

	     jc=acnr(ic,il)

ccccc -----------------------------------------------------------------c
c     Ausschlu? der Quellung bei insoluble, Russ und den               c
c     mineralischen Komponenten und bei den Wolken                     c                          c
ccccc -----------------------------------------------------------------c

            if ( jc.eq.1.or.jc.eq.3.or.(jc.ge.6.and.jc.le.9).or.
     *           jc.gt.10 ) then
               iht=1
               nta=700+(jc*10)+1
            else
c               iht=khum(ihum)
                iht=ihum
                nta=700+(jc*10)+iht
            end if

ccccc -----------------------------------------------------------------c
c     Bestimmung des Filenamens der gesuchten Komponente aus	       c
c     Komponentennummer und Feuchteklasse                              c
ccccc -----------------------------------------------------------------c

          tap(1:15)='extdata/optdat/'
          tap(16:19)=comnam(jc)
          tap(20:21)=chum(iht)
c			write(10,*) tap
            ntap=20

ccccc -----------------------------------------------------------------c
c     Bestimmung der Kennnummer f?r den Puffer ?ber die Wellenl?nge    c
c                                                                      c
c         nbuf: Kennummer der aktuellen Komponente fuer den Puffer     c
c     kbuf(20): Kennummern der gespeicherten Komponenten               c
c         mbuf: Position der aktuellen Komponente im Puffer            c
c         ibuf: Position bis zu der der Puffer belegt ist              c
ccccc -----------------------------------------------------------------c

	     nbuf=ilamb*1000+nta
c 	     print*,'nbuf= ',nbuf

ccccc -----------------------------------------------------------------c
c      ?berpr?fung des Puffers auf ?bereinstimmung mit nbuf            c
ccccc -----------------------------------------------------------------c

	     exists=.false.
	     do ib=1,20
                if (nbuf.eq.kbuf(ib)) then
                   exists=.true.
                   mbuf=ib
                   goto 10
	        end if
	     end do
   10      continue
c           print*,'mbuf= ',mbuf

ccccc -----------------------------------------------------------------c
c      Einlesen der Komponentendaten, falls sie nicht im Puffer        c
c      stehen, sonst ?bernahme aus dem Puffer                          c
ccccc -----------------------------------------------------------------c

c       print*,exists

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
c               print*,' ibuf= ',ibuf
               ibuf=1
            end if
c            print*,ibuf

            kbuf(ibuf)=nbuf
            open (ntap,file=tap,iostat=ios)
c            print*,'opened file ',tap,iostat

            if (ios.ne.0) then
               print*,' error while opening file ',tap
               print*,'  latitude: ',latx
               print*,' longitude: ',lonx
               stop
            end if


c ALTER INPUT
c            read (ntap,200) dum
c            read (ntap,'(22(/))')
c            rlamb=0.
c            do while (rlamb.ne.alamb(ilamb))
c            do ila=1,ilamb
c               read (ntap,500) rlamb,extco,sisca,asymf,exn,refr,refi
c  500          format(2x,8e10.3)
c            end do
c            do ila=ilamb+1,mlamb
c               read (ntap,500) rl
c               print*,rl
c            end do
c
c            read (ntap,'(4(/))')
c
c            ntheta=96
c            do it=1,ntheta
c
c               read (ntap,510,end=511)
c     *               thet,(pha(il,ic,it),ila=1,ilamb)
c  510          format(1x,70e10.3)
c
c               print*,it,thet,pha(il,ic,it)
c
c            end do
c  511       continue
c
c
c  ENDE ALTER INPUT

         do iline=1,100
            read (ntap,220) dum2
            if (dum2.eq.'# optical ') then
               goto 2002
            end if
         end do
 2002    continue
         do iline=1,5
            read (ntap,200) dum
         end do

         do ila=1,mlamb
            read (ntap,500) rlamb,extco,scaco,absco,sisca,asymf,
     *                      exn,refr,refi
  500       format(2x,7e10.3,2e11.3)

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
     *      thet,(pha(il,ic,min(112,it)),ila=1,ilamb)
  510       format(e11.3,1x,70e10.3)
            it=it+1
         end do
  511    ntheta=it-1

c ENDE NEUER INPUT

            bac(il,ic)=pha(il,ic,ntheta)

            extbuf(ibuf)=ext(il,ic)
            scabuf(ibuf)=sca(il,ic)
            absbuf(ibuf)=abs(il,ic)
            sisbuf(ibuf)=sis(il,ic)
            asybuf(ibuf)=asy(il,ic)
            brebuf(ibuf)=bre(il,ic)
            bimbuf(ibuf)=bim(il,ic)
            bacbuf(ibuf)=bac(il,ic)

            close (ntap)

c            print*,'closed file ',ntap

         end if
	 end do
      end do

ccccc -----------------------------------------------------------------c
c     Formate							       c
ccccc -----------------------------------------------------------------c

  100  format(8e10.3)
  200  format(a1)
  220  format(a10)
  300  format(15x,f6.3,/,8x,e10.4,7x,e10.4,7x,e10.4,5x,f7.4,7x,f6.4,/)
  301  format(8x,f6.3//)
  303  format(7x,e10.3,7x,e10.3,7x,e10.3,7x,f7.4/)
  400  format(12(/))
 1010  format(70X,e10.3)

       return
       end

CCCCC *****************************************************************C
      subroutine optpar (ilamb,ihum)
C     *****************************************************************C
C                                                                      c
C     -----------------------------------------------------------------C
C     Berechnung und Ausdruck der gew?nschten optischen Parameter      C
CCCCC -----------------------------------------------------------------C

      integer prnr,acnr,rht
      real n,numden
	real lat5,lon5,relhum,season
	dimension optdep(25,2)
      REAL EXTN(2),ABSN(2),SCAN(2),PF18N(2),supf(112),phafu(112,2)
      REAL EXTA(2),ABSA(2),SCAA(2),SSA(2),ASF(2),PF18A(2)
      real scar(2),absr(2),omer(2)
      character*3  atn,pat
      character*4  comnam
      character*8  opanam,optnam
      character*20 catyp
      character*30 typnam
      common /atyp/   natyp,mcomp,ncomp(10),numden,
     *                catyp,typnam,comnam(20)
      common /norm/   norm,mixnor
      common /layer/  nlay,mlay,nltyp(10),parlay(10,2),boundl(10),
     *                boundu(10)
      common /profi/ nil(10),hfta(10),hstra(10),
     *                h0(2,10),h1(2,10),hm(2,10)
      common /opar/   mopar,jnopar(13),nop,opanam(13),optnam(13)
      common /wavel/  mlamb,alamb(61),niw
      common /mipoi/  latx,lonx,nl,prnr,rht(2),n(2),
     *                njc(2),acnr(5,2),acmr(5,2),nh(2),atn(2),pat(2)
      common /oppoi/  ext(2,5),sca(2,5),abs(2,5),sis(2,5),asy(2,5),
     *                bac(2,5),pha(2,5,112),bre(2,5),bim(2,5)
      COMMON /FTASTR/ EXTFTA(61),EXTSTR(61)
      COMMON /TEST/   ITEST(2)
      common /out/    oparam(10,2),phaf(112,2)
      common /prog/   nprog
      common /angle/  jnangle(112),angle(112),nia
      common /masse/  smas(10,8),smag(8)
	common /r/      lat5,lon5,relhum,season,optdep

CCCCC ------------------------------------------------------------------C
C     MISCHEN DES AEROSOL-TYPS                                          C
C     SUMM(E,A,S) : SUMME EXTINCTION, ABSORPTION, SCATTERING            C
C     SUPF18      : SUMME DES RUECKSTREUKOEFFIZIENTEN                   C
C     SUMASF      : ZWISCHENSUMME DES ASYMMETRIEFAKTORS (ASF)           C
C     SUMASF      : ZWISCHENSUMME DER SINGLE SCAT. ALB. (SSA)           C
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

c      print*,' Berechnung der Summen'

      SUMME = SUMME + ACMR(JC,L)*EXT(l,jc)
      SUMMA = SUMMA + ACMR(JC,L)*ABS(l,jc)
      SUMMS = SUMMS + ACMR(JC,L)*SCA(l,jc)
      SUMSSA = SUMSSA + ACMR(JC,L)*sis(l,jc)
     +                  *EXT(l,jc)
      SUMASF = SUMASF + ACMR(JC,L)*asy(l,jc)
     +                  *SCA(l,jc)
      SUPF18 = SUPF18 + ACMR(JC,L)*bac(l,jc)
      if (jnopar(10).eq.1) then
         do it=1,112
            supf(it)=supf(it)+acmr(jc,l)*pha(l,jc,it)
         end do
      end if

c      print*,jc,l,njc(l),summe,acmr(jc,l),ext(l,jc)

   20 CONTINUE

CCCCC -----------------------------------------------------------------C
C     Normierte  optische Parameter				       c
CCCCC -----------------------------------------------------------------C

c      print*,' Berechnung der normierten Werte'

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
C     ABSOLUTE OPTISCHE PARAMETER                                      C
CCCCC -----------------------------------------------------------------C

c      print*,' Berechnung der absoluten Werte'

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
         scar(l)=scaa(l)/smag(ihum)*1000.  ! Einheit m**2/g
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
C     AUSGABE DER DATEN						       C
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
C     OPTISCHE DICKE                                                   C
CCCCC -----------------------------------------------------------------C

      if (jnopar(6).eq.1) then

CCCCC -----------------------------------------------------------------C
c     Bestimmung von HM, HFTA, HSTR, EXTFTA, EXTSTR aus den	       c
c     eingelesenen Werten in /layer/ fuer RAWOPT		       c
CCCCC -----------------------------------------------------------------C

      if (nprog.eq.2) then
         do il=1,nl
            if (nltyp(il).eq.1) then
               hm(il,1)=parlay(il,1)
            else if (nltyp(il).eq.2) then
               hm(il,1) = parlay(il,2)*
     *                    (exp(-boundl(il)/parlay(il,2))+
     *                     exp(-boundu(il)/parlay(il,2)))
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
     +                   + EXTSTR(ilamb)*HSTRA(prnr)

c        do il=1,nl
c        print*,'   exta= ',exta(il),'    hm(il)= ',hm(il,prnr)
c        end do
c        print*,' ilamb= ' ,ilamb
c        print*,' extfta= ',extfta(ilamb),' hftae= ',hftae
c        print*,' extstr= ',extstr(ilamb),' hstr= ',hstra(prnr)

      odeptha=odepth/alog(10.)

      turbr=0.008569*alamb(ilamb)**(-4)*(1.+0.0113*alamb(ilamb)**
     *         (-2)+0.00013*alamb(ilamb)**(-4))

      turbf=(odepth+turbr)/turbr

      if (jnopar(9).eq.1) then
         kop=iop
      else
         kop=iop+1
c        print*,' exta= ', exta(il),' hm= ',hm(il,1)
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

c      print*,'Aufruf von out4'

      call out4(iop,ilamb,ihum)

c      print*,'Ende out4'

      RETURN
      END

ccccc *****************************************************************c
      subroutine out4(iop,ilamb,ihum)
c     *****************************************************************c
c			                                                     c
c     -----------------------------------------------------------------c
c     AUSGABE DER DATEN	f?r atlopt				                 c
ccccc -----------------------------------------------------------------c

      integer prnr,acnr,njc,rht
      real n
	real lat5,lon5,relhum,season,optdep
	dimension optdep(25,2)

      character*2 chum
      character opanam*8,atn*3,pat*3,optnam*8

      common /angle/  jnangle(112),angle(112),nia
      common /wavel/  mlamb,alamb(61),niw
      common /hum/    khum(8),ahum(8),nih,nhum(8),mhum,chum(8)
      common /out/    oparam(10,2),phaf(112,2)
      common /opar/   mopar,jnopar(13),nop,opanam(13),optnam(13)
      common /mipoi/  latx,lonx,nl,prnr,rht(2),n(2),
     *                njc(2),acnr(5,2),acmr(5,2),nh(2),atn(2),pat(2)
	common /r/      lat5,lon5,relhum,season,optdep

      do l=1,nl
	 if (nih.ne.0) rht(l)=ahum(ihum)
	 if (nop.gt.iop) then
	    if (jnopar(9).eq.1.and.jnopar(10).eq.1) then
	       oparam(iop,l)=oparam(nop-1,l)
	    else if (jnopar(9).eq.1.and.jnopar(10).eq.0) then
	       oparam(iop,l)=oparam(nop,l)
	    end if
	 end if

	 if (iop.le.10) then
	    if (l.eq.1) then
c	       write (10,2020) latx,lonx,nl,alamb(ilamb),ahum(ihum),
c     *			       (oparam(ip,l),ip=1,iop)
c 2020	  FORMAT(2(1x,I4),i3,3x,f6.3,3x,f3.0,1p3e10.3,0p3e10.3,1pe10.3)
           optdep(ilamb,1)=alamb(ilamb)
	     optdep(ilamb,2)=oparam(6,l)
	    else
c	       write (10,3020)
c     *			       (oparam(ip,l),ip=1,iop)
c 3020          FORMAT(13x,1p3e10.3,0p3e10.3,1pe10.3)
	    end if
	 else
c	    if (l.eq.1) then
c	       write (10,2020) latx,lonx,nl,
c     *			       (oparam(ip,l),ip=1,5)
c	       write (10,2030) (oparam(ip,l),ip=6,iop)
c 2030	       FORMAT(11x,1p10e10.3)
c	    else
c	       write (10,3040)
c     *			       (oparam(ip,l),ip=1,5)
c 3040	       FORMAT(13x,10e10.3)
c	       write (10,2030) (oparam(ip,l),ip=6,iop)
c	    end if
	 end if

c	 if(jnopar(10).eq.1) then
c	    write(10,4002)
c 4002	    format('   phase function')
c	    write(10,4010) (phaf(it,l),it=1,nia)
c 4010	    format(8e10.3)
c	    write(10,*) ' '
c	 end if
      end do

      return
      end
