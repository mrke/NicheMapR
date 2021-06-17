#' gads.r
#'
#' R function to extract data from the Fortran Global Aerosol Database software,
#' for a specific location, season and relative humidity
#' see http://www.mpimet.mpg.de/fileadmin/publikationen/Reports/MPI-Report_243.pdf
#' @param lat Latitude in decimal degrees
#' @param lon Longitude in decimal degrees
#' @param relhum Integer specifying relative humidity to use, 1:8 corresponds to 0,50,70,80,90,95,98,99 percent relative humidity
#' @param season Season to obtain data for, 0 = summer, 1 = winter
#' @return optdep A vector of wavelength-specific optical depths
#' @export
gads.r <- function(lat, lon, relhum, season){
  lib.path <- .libPaths()[1]
  lat5s <- seq(-90, 90, 5) #lat range for GADS
  lon5s <- seq(-180, 175, 5) #long range for GADS
  lat5 <- lat5s[which.min(abs(lat5s - lat))] # get nearest latitude square for input location
  lon5 <- lon5s[which.min(abs(lon5s - lon))] # get nearest longitude square for input location
  optdep <- array(0, dim = c(25, 2)) # output array
  niw <- 1
  njh <- 1
  nih <- 1
  lata <- lat5
  late <- lat5
  lati <- 5
  lona <- lon5
  lone <- lon5
  loni <- 5
  nlmal <- 1
  norm <- 1
  nprog <- 4
  ip <- 0
  iwel <- 1
  ihum <- relhum
  il <- 1
  ih <- 1
  ilat <- lata
  imal <- 1
  ilon <- lona
  nwel <- 25


  jnopar <- c(1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0)
  nop <- 7
  kop <- nop
  optnam <- c('ext.coef', 'sca.coef', 'abs.coef', 'sisc.alb', 'asym.par',
              'op.depth','        ', 'turb.fac', 'li.ratio', 'pha.func',
              'ext.rat ','abs.rat ', '        ')
  opanam <- rep("", 13)

  alamb <- c(0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,
             0.9,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.2,3.39,3.5,3.75,
             4.0,4.5,5.0,5.5,6.0,6.2,6.5,7.2,7.9,8.2,8.5,8.7,9.0,
             9.2,9.5,9.8,10.0,10.6,11.0,11.5,12.5,13.0,14.0,14.8,
             15.0,16.4,17.2,18.0,18.5,20.0,21.3,22.5,25.0,27.9,30.,
             35.0,40.0)
  mlamb <- 61

  ahum <- c(0,50,70,80,90,95,98,99)
  nhum <- c(0,50,70,80,90,95,98,99)
  mhum <- 8
  chum <- c('00','50','70','80','90','95','98','99')
  comnam <- c('inso','waso','soot','ssam','sscm','minm','miam',
              'micm','mitr','suso','stco','stma','cucc','cucp',
              'cuma','fog-','cir1','cir2','cir3','    ')
  atn <- c("", "")
  pat <- c("", "")
  rht <- c(0, 0)
  n <- c(0, 0)
  nh <- c(0, 0)
  njc <- c(0, 0)
  acnr <- array(0, dim = c(5, 2))
  acmr <- array(0, dim = c(5, 2))
  khum <- rep(0, 8)
  ext <- array(0, dim = c(2, 5))
  sca <- array(0, dim = c(2, 5))
  abs <- array(0, dim = c(2, 5))
  sis <- array(0, dim = c(2, 5))
  asy <- array(0, dim = c(2, 5))
  bac <- array(0, dim = c(2, 5))
  pha <- array(0, dim = c(112, 2, 5))
  bre <- array(0, dim = c(2, 5))
  bim <- array(0, dim = c(2, 5))
  kbuf <- rep(0, 20)
  extbuf <- rep(0, 20)
  scabuf <- rep(0, 20)
  absbuf <- rep(0, 20)
  sisbuf <- rep(0, 20)
  asybuf <- rep(0, 20)
  bacbuf <- rep(0, 20)
  phabuf <- array(0, dim = c(112, 20))
  brebuf <- rep(0, 20)
  bimbuf <- rep(0, 20)
  extn <- rep(0, 2)
  absn <- rep(0, 2)
  scan <- rep(0, 2)
  pf18n <- rep(0, 2)
  supf <- rep(0, 112)
  phafu <- array(0, dim = c(112, 2))
  exta <- rep(0, 2)
  absa <- rep(0, 2)
  scaa <- rep(0, 2)
  ssa <- rep(0, 2)
  asf <- rep(0, 2)
  pf18a <- rep(0, 2)
  scar <- rep(0, 2)
  absr <- rep(0, 2)
  omer <- rep(0, 2)
  jnangle <- rep(0, 112)
  ncomp <- rep(0, 10)
  oparam <- array(0, dim = c(10, 2))
  nltyp <- rep(0, 10)
  parlay <- array(0, dim = c(10, 2))
  boundl <- rep(0, 10)
  boundu <- rep(0, 10)
  #     READING THE HEIGHT PROFILES from the file TAPE9                  c
  #     						                                           C
  #     HM    : EFFECTIVE LAYER THICKNESS (HOMOGENEOUS DISTRIBUTION)     C
  #     HFTA  : FREE TROP LAYER THICKNESS. AEROSOLS IN KM                C
  #     HSTRA : LAYER THICKNESS OF THE STRATOSPH. AEROSOLS IN KM         C
  #     NIL   : NUMBER OF LAYERS
  iip <- 7
  nil <- c(1, 1, 1, 1, 2, 1, 2, 0, 0)
  hfta <- c(10, 2, 10, 10, 8.5, 6, 8.5, 0, 0, 0)
  hstra <- c(rep(23, 7), 0, 0, 0)
  h0 <- matrix(data = 0, 2, 10)
  h1 <- h0
  hm <- h0
  h0[2, c(5,7)] <- 2
  h1[1, 1:7] <- c(2, 10, 2, 2, 2, 6, 2)
  h1[2, c(5,7)] <- 3.5
  hm[1, 1:7] <- c(2, 5.7, 1.77, 0.86, 0.86, 1.9, 1.77)
  hm[2, c(5,7)] <- 1.5
  extback <- read.table(paste0(lib.path, '/NicheMapR/extdata/extback.dat'), skip = 2)
  extfta <- extback[, 2]
  extstr <- extback[, 3]

  # some definitions for this version

  for(i in 1:13){
    if(jnopar[i] == 1) {
      ip <- ip + 1
      opanam[ip] <- optnam[i]
    }
  }
  # Labeling for different wavelengths and rel. Damp
  opnam <- rep(0, kop)
  opnam[1:kop] <- opanam[1:kop]


  if(season == 1){
    ws <- 'w'
  }else{
    ws <- 's'
  }

  # Input: wavelength

  if(ws == 'w'){
    con <- file(paste0(lib.path, '/NicheMapR/extdata/glodat/winter.dat'), open = "r")
    ntape <- readLines(con) # empty
    close(con)
    cseas <- 'winter '
  }else{
    con <- file(paste0(lib.path, '/NicheMapR/extdata/glodat/summer.dat'), open = "r")
    ntape <- readLines(con, ) # empty
    close(con)
    cseas <- 'summer '
  }
  for(k in 2:length(ntape)){
    line.read <- ntape[k]
    latx <- as.numeric(substr(line.read, 2, 4))
    lonx <- as.numeric(substr(line.read, 5, 8))
    #cat(paste(latx, lonx, '\n'))
    if(!is.na(latx)){
      if(latx == ilat & lonx == ilon){
        break
      }
    }
  }

  for(iwel in 1:nwel){
    if(exists("k2")){
      rm(k2)
    }
    ilamb <- iwel

    # c      Reading the raw data from the files TAPE201, TAPE207:	       c
    # c     ------------------------------------------------------	       c
    # C     LAT       : LATITUDE                                             C
    # C     LON       : LONGITUDE                                            C
    # C     NL        : NUMBER OF AEROSOL LAYERS                             C
    # C                 (=2 FOR MARITIME-MINERAL,=1 FOR MARI.)               C
    # C     PRNR      : PROFIL NUMBER                                        C
    # C     NT        : NUMBER OF AEROSOL TYPE                               C
    # C     PAT       : AEROSOL PROFIL TYPE                                  C
    # C     NH        : NUMBER OF REL. HUMIDITY CLASS                        C
    # C     N         : TOTAL NUMBER CONCENTRATION                           C
    # C     NJC       : NUMBER OF AEROSOL COMPONENT                          C
    # C     ACNR      : AEROSOL COMPONENT                                    C
    # C     ACMR      : MIXING RATIO                                         C
    # C                 (PARTIAL NUMBER CONCENTRATION/TOTAL NUMBER CONC.)    C

    #C     READING IN THE DATA FROM THE RAW DATA-FILES TAPE201-TAPE212        C

    nl <- as.numeric(substr(line.read, 9, 11))
    prnr <- as.numeric(substr(line.read, 12, 14))
    atn[1] <- substr(line.read, 15, 17)
    pat[1] <- substr(line.read, 18, 20)
    rht[1] <- as.numeric(substr(line.read, 21, 23))
    n[1] <- as.numeric(substr(line.read, 24, 33))
    njc[1] <- as.numeric(substr(line.read, 34, 37))
    acnr[1, 1] <- as.numeric(substr(line.read, 38, 40))
    acmr[1, 1] <- as.numeric(substr(line.read, 41, 50))
    acnr[2, 1] <- as.numeric(substr(line.read, 51, 53))
    acmr[2, 1] <- as.numeric(substr(line.read, 54, 63))
    acnr[3, 1] <- as.numeric(substr(line.read, 64, 66))
    acmr[3, 1] <- as.numeric(substr(line.read, 67, 76))
    if(njc[1] > 3){
      k2 <- k
      for(jc in 4:njc[1]){
        k2 <- k2 + 1
        line.read2 <- ntape[k2]
        acnr[jc, 1] <- as.numeric(substr(line.read2, 38, 40))
        acmr[jc, 1] <- as.numeric(substr(line.read2, 41, 50))
      }
    }
    if(nl != 1){
      for(l in 2:nl){
        if(!exists("k2")){
          k2 <- k
        }
        k2 <- k2 + 1
        line.read2 <- ntape[k2]
        atn[l] <- substr(line.read2, 15, 17)
        pat[l] <- substr(line.read2, 18, 20)
        rht[l] <- as.numeric(substr(line.read2, 21, 23))
        n[l] <- as.numeric(substr(line.read2, 24, 33))
        njc[l] <- as.numeric(substr(line.read2, 34, 37))
        acnr[1, l] <- as.numeric(substr(line.read2, 38, 40))
        acmr[1, l] <- as.numeric(substr(line.read2, 41, 50))
        acnr[2, l] <- as.numeric(substr(line.read2, 51, 53))
        acmr[2, l] <- as.numeric(substr(line.read2, 54, 63))
        acnr[3, l] <- as.numeric(substr(line.read2, 64, 66))
        acmr[3, l] <- as.numeric(substr(line.read2, 67, 76))
        if(njc[l] > 3){
          for(jc in 4:njc[l]){
            k2 <- k2 + 1
            line.read2 <- ntape[k2]
            acnr[jc, l] <- as.numeric(substr(line.read2, 38, 40))
            acmr[jc, l] <- as.numeric(substr(line.read2, 41, 50))
          }
        }
      }
    }
    for(l in 1:nl){
      sum <- 0
      for(ic in 1:njc[l]){
        sum <- sum + acmr[ic, l]
      }
      if(abs(sum -1) >= 0.01){
        cat('sum of mixing ratios is not 1. please have a look at errorfile *.err')
      }
    }
    #     DETERMINATION OF THE AEROSOL TYPE AND HUMIDITY CLASS             C
    #     !!! AEROSOL TYPE IS LAYER-DEPENDENT (NOT TAKEN INTO ACCOUNT)     C

    for(il in 1:nl){
      if(rht[il] < 30){
        nh[il] <- 1
      }else{
        if(rht[il] > 30 & rht[il] <= 65){
          nh[il] <- 2
        }else{
          if(rht[il] > 65 & rht[il] <= 75){
            nh[il] <- 3
          }else{
            if(rht[il] > 75 & rht[il] <= 85){
              nh[il] <- 4
            }else{
              if(rht[il] > 85 & rht[il] <= 92){
                nh[il] <- 5
              }else{
                if(rht[il] > 92 & rht[il] <= 97){
                  nh[il] <- 6
                }else{
                  if(rht[il] == 98){
                    nh[il] <- 7
                  }else{
                    if(rht[il] == 99){
                      nh[il] <- 8
                    }}}}}}}}}

    #      Reading of the optical raw data from the files winter.dat and   c
    #      summer.dat

    # subroutine optcom

    #      Calculation of the optical parameters at the current grid point
    ibuf <- 0

    for(il in 1:nl){

      if(nih == 0){
        for(ihu in 1:mhum){
          if(nh[il] == ihu){
            khum(ihum) <- ihu
          }
        }
      }

      for(ic in 1:njc[il]){
        jc <- acnr[ic, il]
        # Exclusion? swelling at insoluble, soot and mineral components and at clouds
        if(jc == 1 | jc == 3 | (jc >= 6 & jc <= 9) | jc > 10){
          iht <- 1
          nta <- 700 + (jc * 10) + 1
        }else{
          iht <- ihum
          nta <- 700 + (jc * 10) + iht
        }
        # Determination of the file name of the sought component from component number and moisture class
        tap <- paste0(lib.path, '/NicheMapR/extdata/optdat/', comnam[jc], chum[iht])
        # c     Determination of the identification number for the buffer over   c
        # c     the wavelength                                                   c
        # c                                                                      c
        # c         nbuf: ID number of the current component for the buffer      c
        # c     kbuf(20): ID numbers of the stored components                    c
        # c         mbuf: Position of the current component in the buffer        c
        # c         ibuf: Position up to which the buffer is occupied            c

        nbuf <- ilamb * 1000 + nta

        # Check the buffer for compliance with nbuf
        exist.chk <- FALSE
        for(ib in 1:20){
          if(nbuf == kbuf[ib]){
            exist.chk <- TRUE
            mbuf <- ib
            break
          }
        }

        # Reading the component data if it is not in the buffer,
        # otherwise transferring it from the buffer

        if(exist.chk){
          ext[il,ic] <- extbuf[mbuf]
          sca[il,ic] <- scabuf[mbuf]
          abs[il,ic] <- absbuf[mbuf]
          sis[il,ic] <- sisbuf[mbuf]
          asy[il,ic] <- asybuf[mbuf]
          bac[il,ic] <- bacbuf[mbuf]
          bre[il,ic] <- brebuf[mbuf]
          bim[il,ic] <- bimbuf[mbuf]
        }else{
          if(ibuf < 20){
            ibuf <- ibuf + 1
          }else{
            ibuf <- 1
          }

          kbuf[ibuf] <- nbuf
          con <- file(tap, open = "r")
          ntap <- readLines(con) # empty
          close(con)
          ntap <- na.omit(ntap)
          for(iline in 1:100){
            if(ntap[iline] == '# optical parameters:'){
              break
            }
          }

          for(ila in (iline + 6):(iline + 5 + mlamb)){
            ntap.out <- ntap[ila]
            rlamb <- as.numeric(substr(ntap.out, 3, 12))
            extco <- as.numeric(substr(ntap.out, 13, 22))
            scaco <- as.numeric(substr(ntap.out, 23, 32))
            absco <- as.numeric(substr(ntap.out, 33, 42))
            sisca <- as.numeric(substr(ntap.out, 43, 52))
            asymf <- as.numeric(substr(ntap.out, 53, 62))
            exn <- as.numeric(substr(ntap.out, 63, 72))
            refr <- as.numeric(substr(ntap.out, 73, 83))
            refi <- as.numeric(substr(ntap.out, 84, 94))

            if(rlamb == alamb[ilamb]){
              ext[il,ic] <- extco
              sca[il,ic] <- scaco
              abs[il,ic] <- absco
              sis[il,ic] <- sisca
              asy[il,ic] <- asymf
              bre[il,ic] <- refr
              bim[il,ic] <- refi
            }
          }
          ntheta <- length(ntap) - (iline + 5 + 8 + mlamb)
          ntap.mat <- matrix(data = NA, nrow = 112, ncol = 62)
          for(ii in 1:ntheta){
            ntap.mat[ii, ] <-  na.omit(as.numeric(unlist(strsplit(ntap[(ila + 8 + ii)], split = " "))))
          }

          pha[1:ntheta, il, ic] <- ntap.mat[, 2]


          bac[il, ic] <- pha[min(ntheta, ntheta), il, ic]

          extbuf[ibuf] <- ext[il,ic]
          scabuf[ibuf] <- sca[il,ic]
          absbuf[ibuf] <- abs[il,ic]
          sisbuf[ibuf] <- sis[il,ic]
          asybuf[ibuf] <- asy[il,ic]
          brebuf[ibuf] <- bre[il,ic]
          bimbuf[ibuf] <- bim[il,ic]
          bacbuf[ibuf] <- bac[il,ic]
        }
      }
    }

    # subroutine optpar

    for(l in 1:nl){
      summe <- 0
      summa <- 0
      summs <- 0
      sumssa <- 0
      sumasf <- 0
      supf18 <- 0
      if(jnopar[10] == 1){
        for(it in 1:112){
          supf[it] <- 0
        }
      }

      for(jc in 1:njc[l]){
        # Calculation of the sums
        summe <- summe + acmr[jc, l] * ext[l,jc]
        summa <- summa + acmr[jc, l] * abs[l,jc]
        summs <- summs + acmr[jc, l] * sca[l,jc]
        sumssa <- sumssa + acmr[jc, l] * sis[l,jc] * ext[l,jc]
        sumasf <- sumasf + acmr[jc, l] * asy[l,jc] * sca[l,jc]
        supf18 <- supf18 + acmr[jc, l] * bac[l,jc]
        if(jnopar[10] == 1){
          for(it in 1:112){
            supf[it] <- supf[it] + acmr[jc, l] * pha[l, jc, it]
          }
        }
      }
      #     Standardized optical parameters				       c
      #     Calculation of the standardized values

      extn[l] <- summe
      absn[l] <- summa
      scan[l] <- summs
      pf18n[l] <- supf18
      if(jnopar[10] == 1){
        for(it in 1:112){
          phafu[it,l] <- supf[it]
        }
      }

      ssa[l] <- sumssa / summe
      asf[l] <- sumasf / summs

      #     ABSOLUTE OPTICAL PARAMETERS                                      C
      exta[l] <- extn[l] * n[l]
      absa[l] <- absn[l] * n[l]
      scaa[l] <- scan[l] * n[l]
      pf18a[l] <- pf18n[l] * n[l]
      if(jnopar[10] == 1 & norm ==1){
        for(it in 1:112){
          phafu[it, l] <- phafu[it, l] * n[l]
        }
      }


      if(norm == 1){
        extn[l] <- exta[l]
        absn[l] <- absa[l]
        scan[l] <- scaa[l]
        pf18n[l] <- pf18a[l]
      }

      if(jnopar[10] == 1){
        itp <- 1
        for(it in 1:112){
          if(jnangle[it] == 1){
            phaf(itp, l) <- phafu(it, l)
            itp <- itp + 1
          }
        }
      }

      if(jnopar[11] == 1){
        scar[l] <- scaa[l] / smag[ihum] * 1000  # unit m**2/g
      }

      if(jnopar[12] == 1){
        absr[l] <- absa[l] / smag[ihum] * 1000
      }

      if(jnopar[13] == 1){
        kc <- 0
        for(jc in 1:njc[l]){
          if(ncomp[jc] == 3){
            kc <- jc
          }
        }
        if(kc != 0){
          omer[l] <- ssa[l] / smas[kc, ihum]
        }else{
          omer[l] <- 99
        }
      }

      iop <- 0
      kop <- 0
      if(jnopar[1] == 1){
        iop <- iop + 1
        oparam[iop, l] <- extn[l]
      }
      if(jnopar[2] == 1){
        iop <- iop + 1
        oparam[iop, l] <- scan[l]
      }
      if(jnopar[3] == 1){
        iop <- iop + 1
        oparam[iop, l] <- absn[l]
      }
      if(jnopar[4] == 1){
        iop <- iop + 1
        oparam[iop, l] <- ssa[l]
      }
      if(jnopar[5] == 1){
        iop <- iop + 1
        oparam[iop, l] <- asf[l]
      }
      if(jnopar[9] == 1){
        iop <- iop + 1
        if(jnopar[6] == 1){kop <- kop + 1}
        if(jnopar[7] == 1){kop <- kop + 1}
        if(jnopar[8] == 1){kop <- kop + 1}
        kop <- kop + iop
        oparam[kop, l] <- exta[l] / pf18a[l]
      }
      if(jnopar[11] == 1){
        iop <- iop + 1
        oparam[iop, l] <- scar[l]
      }
      if(jnopar[12] == 1){
        iop <- iop + 1
        oparam[iop, l] <- absr[l]
      }
      if(jnopar[13] == 1){
        iop <- iop + 1
        oparam[iop, l] <- omer[l]
      }
    }

    #     optical thickness                                                c

    if(jnopar[6] == 1){
      if(nprog == 2){
        for(il in 1:nl){
          if(nltyp[il] == 1){
            hm[il,1] <- parlay[il,1]
          }else{
            if(nltyp[il] == 2){
              hm[il,1] <- parlay[il, 2] + (exp(-1 * boundl[il] / parlay[il, 2]) + exp(-boundu(il)/parlay(il,2)))
            }
          }
        }
        hfta[1] <- boundu[nlay - 2] - boundl[nlay - 2]
        hstra[1] <- boundu[nlay - 1] - boundl[nlay - 1]
        extfta[ilamb] <- parlay[nlay - 2, 1]
        extstr[ilamb] <- parlay[nlay - 1, 1]
      }

      hu <- h1[nl, prnr]
      ho <- hu + hfta[prnr]
      z <- 8
      hftae <- z * (exp(-1 * hu / z) - exp(-1 * ho / z) )

      odepth <- 0

      for(il in 1:nl){
        odepth <- odepth + exta[il] * hm[il, prnr]
      }

      odepth <- odepth + extfta[ilamb] * hftae + extstr[ilamb] * hstra[prnr]
      odeptha <- odepth / log(10)
      turbr <- 0.008569 * alamb[ilamb] ^ (-4) * (1 + 0.0113 * alamb[ilamb] ^ (-2)
                                                 + 0.00013 * alamb[ilamb] ^ (-4))
      turbf <- (odepth + turbr) / turbr

      if(jnopar[9] == 1){
        kop <- iop
      }else{
        kop <- iop + 1
      }

      if(jnopar[6] == 1){
        oparam[kop, 1] <- odepth
        oparam[kop, 2] <- 0
        kop <- kop + 1
        iop <- iop + 1
      }

      if(jnopar[7] == 1){
        oparam[kop,1] <- odeptha
        oparam[kop,2] <- 0
        kop <- kop + 1
        iop <- iop + 1
      }

      if(jnopar[8] == 1){
        oparam[kop, 1] <- turbf
        oparam[kop, 2] <- 0
        iop <- iop + 1
      }
    }

    # subroutine out4
    for(l in 1:nl){
      if(nih != 0){rht[l] <- ahum[ihum]}
      if(nop > iop){
        if(jnopar[9] == 1 & jnopar[10] == 1){
          oparam[iop, l] <- oparam[nop - 1, l]
        }else{
          if(jnopar[9] == 1 & jnopar[10] == 0){
            oparam[iop, l] <- oparam[nop,l]
          }
        }
      }

      if(iop <= 10){
        if(l == 1){
          optdep[ilamb, 1] <- alamb[ilamb]
          optdep[ilamb, 2] <- oparam[6, l]
        }
      }
    }
  }
  optdep[, 1] <- optdep[, 1] * 1000
  return(optdep)
}
