C     NICHEMAPR: SOFTWARE FOR BIOPHYSICAL MECHANISTIC NICHE MODELLING
C     TEMPLATE EQUIVALENT TO endo_devel.R
C     MODIFY THIS ACCORDINGLY AFTER DEVELOPING A ROUTINE FOR THERMOREGULATION
C     USING endo_devel.R

      subroutine SOLVENDO(INPUT,TREG,MORPH,ENBAL,MASBAL)
     
      implicit none
      
      DOUBLE PRECISION ABSAND,ABSANV,ABSSB,AHEIT,AIRML1,AIRML2,AIRVOL
      DOUBLE PRECISION AK1,AK1inc,AK1LAST,AK2,AKMAX,ALENTH,AMASS,ANDENS
      DOUBLE PRECISION ANU,AREACND,AREASKIN,ASEMAJ,ASIL,ASILN,ASILP,ATOT
      DOUBLE PRECISION AWIDTH,BP,BSEMIN,CD,CO2GAS,CONVAR,CONVOUT,CONVSK
      DOUBLE PRECISION CSEMIN,CZ,D,DELTAR,DHAIRD,DHAIRV,DHARA,DIFTOL
      DOUBLE PRECISION DMULT,ELEV,EMISAN,ENBAL,ENVVARS,EVAPGS,EXTREF
      DOUBLE PRECISION FABUSH,FABUSHREF,FAGRD,FAGRDREF,FASKY,FASKYREF
      DOUBLE PRECISION FATPCT,FATTHK,FAVEG,FAVEGREF,FGDREF,FLSHVL,FLTYPE
      DOUBLE PRECISION FLYHR,FSKREF,FURTHRMK,FURTST,FURVARS,FURWET,GEND
      DOUBLE PRECISION GENV,GEOMout,GEOMVARS,GEVAP,GMASS,GR,HC,HCFOR,HD
      DOUBLE PRECISION HCFREE,HDFORC,HDFREE,HTOVPR,INPUT,IPT,IRPROPout
      DOUBLE PRECISION KEFF,KFURCMPRS,KHAIR,LEN,LHAIRD,LHAIRV,MASBAL
      DOUBLE PRECISION MASFAT,MAXPCOND,MAXPTVEN,MORPH,MXWET,N2GAS,N2MOL1
      DOUBLE PRECISION N2MOL2,NESTYP,NTRYD,NTRYV,O2GAS,O2MOL1
      DOUBLE PRECISION O2MOL2,O2STP,ORIENT,PANT,PANTING,PANTLAST,PANTMAX
      DOUBLE PRECISION PANTCOST,PANTMULT,PCOND,PCTBAREVAP,PCTCO2,PCTDIF
      DOUBLE PRECISION PCTEYES,PCTN2,PCTO2,PI,PR,Q10,Q10mult,QBASAL
      DOUBLE PRECISION QBASREF,QCOND,QCONDD,QCONDV,QCONV,QCONVD,QCONVV
      DOUBLE PRECISION QDORSL,QEVAP,QFSEVAPD,QFSEVAPV,QGEN,QGENNETD
      DOUBLE PRECISION QGENNETV,QIR,QIRIN,QIRIND,QIRINV,QIROUT,QIROUTD
      DOUBLE PRECISION QIROUTV,QM1,QM2,QMET,QMIN,QNORM,QRADD,QRADV
      DOUBLE PRECISION QRBSHD,QRBSHV,QRESP,QRGRDD,QRGRDV,QRSKYD,QRSKYV
      DOUBLE PRECISION QRVEGD,QRVEGV,QSDIFF,QSDIR,QSEVAPD,QSEVAPV,QSLR
      DOUBLE PRECISION QSLRD,QSLRV,QSOL,QSOLAR,QSOLR,QSRSB,QSSKY,QSUM
      DOUBLE PRECISION QVENTR,R1,R2,RA,RAISETC,RE,REFLD,REFLV,RELXIT
      DOUBLE PRECISION RESPFN,RESPGEN,RFLESH,RFUR,RH,RHOARA,RHOD,RHOV
      DOUBLE PRECISION RONEST,RP_CO2,RQ,RRAD,RSKIN,SAMODE,SC,SHADE,SHAPE
      DOUBLE PRECISION SHAPEB,SHAPEB_LAST,SHAPEB_MAX,SHAPEB_REF,SHAPEC
      DOUBLE PRECISION sigma,SIMULOUT,SIMULSOLout,SKINW,SKINWLAST
      DOUBLE PRECISION SOLARout,SUBQFAT,SUCCESSD,SUCCESSV,SURFAR,SWEAT
      DOUBLE PRECISION SWEATGS,TA,TAEXIT,TAREF,TBUSH,TC,TCLAST,TCMAX
      DOUBLE PRECISION TCONDSB,TCREF,TENV,TFA,TFAD,TFAV,TGRD,TIMACT
      DOUBLE PRECISION TLOWER,TLUNG,TOL,TRAITS,TREG,TS,TSKCALCAVD
      DOUBLE PRECISION TSKCALCAVV,TSKINMAX,TSKY,TVEG,UNCURL,VEL,VMULT
      DOUBLE PRECISION VOL,VOLFAT,X,XR,Z,ZBRENTin,ZBRENTout,ZEN,ZFUR
      DOUBLE PRECISION ZFURCOMP,ZFURD,ZFURV,ZL

      DOUBLE PRECISION, DIMENSION(3) :: KEFARA,BETARA,B1ARA,DHAR,LHAR,
     & RHOAR,ZZFUR,REFLFR
     
      INTEGER S
    
      DIMENSION IRPROPout(26),GEOMout(25),CONVOUT(14),
     & SOLARout(7),SIMULSOLout(2,15),SIMULOUT(15),FURVARS(9),
     & GEOMVARS(16),TRAITS(9),ENVVARS(17),ZBRENTin(17),ZBRENTout(15),
     & INPUT(87),TREG(15),MORPH(20),ENBAL(10),MASBAL(10)

      PI = ACOS(-1.0d0)
      
      QGEN=input(1)
      QBASAL=input(2)
      TA=input(3)
      SHAPEB_MAX=input(4)
      SHAPEB_REF=input(5)
      SHAPEB=input(6)
      DHAIRD=input(7)
      DHAIRV=input(8)
      LHAIRD=input(9)
      LHAIRV=input(10)
      ZFURD=input(11)
      ZFURV=input(12)
      RHOD=input(13)
      RHOV=input(14)
      REFLD=input(15)
      REFLV=input(16)
      MAXPTVEN=input(17)
      SHAPE=input(18)
      EMISAN=input(19)
      KHAIR=input(20)
      FSKREF=input(21)
      FGDREF=input(22)
      NESTYP=input(23)
      PCTDIF=input(24)
      ABSSB=input(25)
      SAMODE=input(26)
      FLTYPE=input(27)
      ELEV=input(28)
      BP=input(29)
      RP_CO2=input(30)
      SHADE=input(31)
      QSOLR=input(32)
      RoNEST=input(33)
      Z=input(34)
      VEL=input(35)
      TS=input(36)
      TFA=input(37)
      FABUSH=input(38)
      FURTHRMK=input(39)
      RH=input(40)
      TCONDSB=input(41)
      TBUSH=input(42)
      TC=input(43)
      PCTBAREVAP=input(44)
      FLYHR=input(45)
      FURWET=input(46)
      AK1=input(47)
      AK2=input(48)
      PCTEYES=input(49)
      DIFTOL=input(50)
      SKINW=input(51)
      TSKY=input(52)
      TVEG=input(53)
      TAREF=input(54)
      DELTAR=input(55)
      RQ=input(56)
      TIMACT=input(57)
      O2GAS=input(58)
      N2GAS=input(59)
      CO2GAS=input(60)
      RELXIT=input(61)
      PANT=input(62)
      EXTREF=input(63)
      UNCURL=input(64)
      AKMAX=input(65)
      AK1inc=input(66)
      TCMAX=input(67)
      RAISETC=input(68)
      TCREF=input(69)
      Q10=input(70)
      QBASREF=input(71)
      PANTMAX=input(72)
      MXWET=input(73)
      SWEAT=input(74)
      TGRD=input(75)
      AMASS=input(76)
      ANDENS=input(77)
      SUBQFAT=input(78)
      FATPCT=input(79)
      PCOND=input(80)
      MAXPCOND = input(81)
      ZFURCOMP = input(82)
      PANTING=input(83)
      ORIENT=input(84)
      SHAPEC=input(85)
      XR=input(86)
      PANTMULT=input(87)
      
      TSKINMAX=TC ! initialise
      Q10mult=1. ! initialise
      PANTCOST = 0.D0 ! initialise
      
      do while(QGEN < QBASAL)

       !### IRPROP, infrared radiation properties of fur

       !# call the IR properties subroutine
       CALL IRPROP(TA, SHAPEB_MAX, SHAPEB_REF, SHAPEB, DHAIRD, DHAIRV, 
     &  LHAIRD, LHAIRV, ZFURD, ZFURV, RHOD, RHOV, REFLD, REFLV,
     &  MAXPTVEN, ZFURCOMP, KHAIR, IRPROPout)
      
       !# output
       KEFARA = IRPROPout(1:3) !# effective thermal conductivity of fur array, mean, dorsal, ventral (W/mK)
       BETARA = IRPROPout(4:6) !# term involved in computing optical thickess (1/mK2)
       B1ARA = IRPROPout(7:9) !# optical thickness array, mean, dorsal, ventral (m)
       DHAR = IRPROPout(10:12) !# fur diameter array, mean, dorsal, ventral (m)
       LHAR = IRPROPout(13:15) !# fur length array, mean, dorsal, ventral (m)
       RHOAR = IRPROPout(16:18) !# fur density array, mean, dorsal, ventral (1/m2)
       ZZFUR = IRPROPout(19:21) !# fur depth array, mean, dorsal, ventral (m)
       REFLFR = IRPROPout(22:24) !# fur reflectivity array, mean, dorsal, ventral (fractional, 0-1)
       FURTST = IRPROPout(25) !# test of presence of fur (length x diamater x density x depth) (-)
       KFURCMPRS = IRPROPout(26) ! # effictive thermal conductivity of compressed ventral fur (W/mK)
       
       !### GEOM, geometry

       !# input
       DHARA = DHAR(1) !# fur diameter, mean (m) (from IRPROP)
       RHOARA = RHOAR(1) !# hair density, mean (1/m2) (from IRPROP)
       ZFUR = ZZFUR(1) !# fur depth, mean (m) (from IRPROP)

       !# call the subroutine
       CALL GEOM_ENDO(AMASS,ANDENS,FATPCT,SHAPE,ZFUR,SUBQFAT,SHAPEB,
     &  SHAPEB_REF,SHAPEC,DHARA,RHOARA,PCOND,SAMODE,ORIENT,GEOMout)
 
       !# output
       VOL = GEOMout(1) !# volume, m3
       D = GEOMout(2) !# characteristic dimension for convection, m
       MASFAT = GEOMout(3) !# mass body fat, kg
       VOLFAT = GEOMout(4) !# volume body fat, m3
       ALENTH = GEOMout(5) !# length, m
       AWIDTH = GEOMout(6) !# width, m
       AHEIT = GEOMout(7) !# height, m
       ATOT = GEOMout(8) !# total area, m2
       ASIL = GEOMout(9) !# silhouette area, m2
       ASILN = GEOMout(10) !# silhouette area normal to sun, m2
       ASILP = GEOMout(11) !# silhouette area parallel to sun, m2
       GMASS = GEOMout(12) !# mass, g
       AREASKIN = GEOMout(13) !# area of skin, m2
       FLSHVL = GEOMout(14) !# flesh volume, m3
       FATTHK = GEOMout(15) !# fat layer thickness, m
       ASEMAJ = GEOMout(16) !# semimajor axis length, m
       BSEMIN = GEOMout(17) !# b semiminor axis length, m
       CSEMIN = GEOMout(18) !# c semiminor axis length, m (currently only prolate spheroid)
       CONVSK = GEOMout(19) !# area of skin for evaporation (total skin area - hair area), m2
       CONVAR = GEOMout(20) !# area for convection (total area minus ventral area, as determined by PCOND), m2
       R1 = GEOMout(21) !# shape-specific core-skin radius in shortest dimension, m
       R2 = GEOMout(22) !# shape-specific core-fur/feather interface radius in shortest dimension, m
       
       !### SOLAR, solar radiation

       !# solar radiation normal to sun's rays
       ZEN = pi/180.*Z !# convert degrees to radians
       if(Z.lt.90.)then !# compute solar radiation on a surface normal to the direct rays of the sun
        CZ = cos(ZEN)
        QNORM = QSOLR/CZ
       else !# diffuse skylight only
        QNORM = QSOLR
       endif

       ABSAND = 1 - REFLFR(2) !# solar absorptivity of dorsal fur (fractional, 0-1)
       ABSANV = 1 - REFLFR(3) !# solar absorptivity of ventral fur (fractional, 0-1)

C      CORRECT FASKY FOR % VEGETATION SHADE OVERHEAD, ASHADE
       FAVEG = FSKREF*(SHADE/100.)
       FASKY = FSKREF - FAVEG
       FAGRD = FGDREF

       CALL SOLAR_ENDO(ATOT, ABSAND, ABSANV, ABSSB, ASIL, PCTDIF, 
     &  QNORM, SHADE, QSOLR, FASKY, FAVEG, SOLARout)

       QSOLAR = SOLARout(1) !# total (global) solar radiation (W) QSOLAR,QSDIR,QSSKY,QSRSB,QSDIFF,QDORSL,QVENTR
       QSDIR = SOLARout(2) !# direct solar radiaton (W)
       QSSKY = SOLARout(3) !# diffuse solar radiation from sky (W)
       QSRSB = SOLARout(4) !# diffuse solar radiation reflected from substrate (W)
       QSDIFF = SOLARout(5) !# total diffuse solar radiation (W)
       QDORSL = SOLARout(6) !# total dorsal solar radiation (W)
       QVENTR = SOLARout(7) !# total ventral solar radiaton (W)

       !### CONV, convection

       !# input
       SURFAR = CONVAR !# surface area for convection, m2 (from GEOM)
       TENV = TA !# fluid temperature (?C)

       !# run subroutine
       CALL CONV_ENDO(TS, TENV, SHAPE, SURFAR, FLTYPE, FURTST, D, TFA, 
     &  VEL, ZFUR, BP, ELEV, CONVout)

       QCONV = CONVout(1) !# convective heat loss (W)
       HC = CONVout(2) !# combined convection coefficient
       HCFREE = CONVout(3) !# free convection coefficient
       HCFOR = CONVout(4) !# forced convection coefficient
       HD = CONVout(5) !# mass transfer coefficient
       HDFREE = CONVout(6) !# free mass transfer coefficient
       HDFORC = CONVout(7) !# forced mass transfer coefficient
       ANU = CONVout(8) !# Nusselt number (-)
       RE = CONVout(9) !# Reynold's number (-)
       GR = CONVout(10) !# Grasshof number (-)
       PR = CONVout(11) !# Prandlt number (-)
       RA = CONVout(12) !# Rayleigh number (-)
       SC = CONVout(13) !# Schmidt number (-)
       BP = CONVout(14) !# barometric pressure (Pa)
 
       !# reference configuration factors
       FABUSHREF = FABUSH !# nearby bush
       FASKYREF = FASKY !# sky
       FAGRDREF = FAGRD !# ground
       FAVEGREF = FAVEG !# vegetation

       !### SIMULSOL, simultaneous solution of heat balance
       !# repeat for each side, dorsal and ventral, of the animal

       DO 1, S=1,2
        !# set infrared environment
        TVEG = TAREF !# assume vegetation casting shade is at 1.2 m (reference) air temperature (?C)
        TLOWER = TGRD
        !# Calculating solar intensity entering fur. This will depend on whether we are calculating the fur temperature for the dorsal side or the ventral side. The dorsal side will have solar inputs from the direct beam hitting the silhouette area as well as diffuse solar scattered from the sky. The ventral side will have diffuse solar scattered off the substrate.
        !# Resetting config factors and solar depending on whether the dorsal side (S=1) or ventral side (S=2) is being estimated.
        IF(QSOLAR.GT.0.0)THEN
         if(S==1)THEN
          FASKY = FASKYREF/(FASKYREF+FAVEGREF) ! proportion of upward view that is sky
          FAVEG = FAVEGREF/(FASKYREF+FAVEGREF) ! proportion of upward view that is vegetation (shade)
          FAGRD = 0.0
          FABUSH = 0.0
          QSLR = 2.*QSDIR+((QSSKY/FASKYREF)*FASKY) ! direct x 2 because assuming sun in both directions, and unadjusting QSSKY for config factor imposed in SOLAR_ENDO and back to new larger one in both directions 
         else
          FASKY = 0.0
          FAVEG = 0.0
          FAGRD = FAGRDREF/(FAGRDREF+FABUSHREF)
          FABUSH = FABUSHREF/(FAGRDREF+FABUSHREF)
          QSLR = (QVENTR/(1. - FASKYREF - FAVEGREF))*(1.-(2.*PCOND)) ! unadjust by config factor imposed in SOLAR_ENDO to have it coming in both directions, but also cutting down according to fractional area conducting to ground (in both directions)
        ENDIF
        else
         QSLR = 0.0
         if(S==1)then
          FASKY = FASKYREF/(FASKYREF+FAVEGREF)
          FAVEG = FAVEGREF/(FASKYREF+FAVEGREF)
          FAGRD = 0.0
          FABUSH = 0.0
         else
          FASKY = 0.0
          FAVEG = 0.0
          FAGRD = FAGRDREF/(FAGRDREF+FABUSHREF)
          FABUSH = FABUSHREF/(FAGRDREF+FABUSHREF)
         ENDIF
        ENDIF

        !# set fur depth and conductivity
        !# index for KEFARA, the conductivity, is the average (1), front/dorsal (2), back/ventral(3) of the body part
        if((QSOLR.GT.0).OR.(ZZFUR(2).NE.ZZFUR(3)))THEN
         if(S == 1)THEN
          ZL = ZZFUR(2)
          KEFF = KEFARA(2)
         else
          ZL = ZZFUR(3)
          KEFF = KEFARA(3)
         ENDIF
        else
         ZL = ZZFUR(1)
         KEFF = KEFARA(1)
        ENDIF

        RSKIN = R1 !# body radius (including fat), m
        RFLESH = R1 - FATTHK !# body radius flesh only (no fat), m
        RFUR = R1 + ZL !# body radius including fur, m
        D = 2. * RFUR !# diameter, m
        RRAD = RSKIN + (XR * ZL) !# effective radiation radius, m
        LEN = ALENTH !# length, m

        !# Correcting volume to account for subcutaneous fat
        if((SUBQFAT.EQ.1.).AND.(FATTHK.GT.0.0))THEN
         VOL = FLSHVL
        ENDIF

        !# Calculating the "Cd" variable: Qcond = Cd(Tskin-Tsub), where Cd = Conduction area*((kfur/zfur)+(ksub/subdepth))
        IF(S==2)THEN ! doing ventral side, add conduction
         AREACND = ATOT * (PCOND *2)
         CD = AREACND * ((KFURCMPRS/ZFURCOMP))
         CONVAR = CONVAR - AREACND !# Adjust area used for convection to account for PCOND. This is sent in to simulsol & used for CONV, RAD
        ELSE  !# doing dorsal side, no conduction. No need to adjust areas used for convection. 
         AREACND = 0
         CD = AREACND * ((KFURCMPRS/ZFURCOMP))
        ENDIF

        !# package up inputs
        FURVARS = (/LEN,ZFUR,FURTHRMK,KEFF,BETARA,FURTST,ZL/)
        GEOMVARS = (/SHAPE,SUBQFAT,CONVAR,VOL,D,CONVAR,CONVSK,RFUR,
     &   RFLESH,RSKIN,XR,RRAD,ASEMAJ,BSEMIN,CSEMIN,CD/)
        ENVVARS = (/FLTYPE,TA,TS,TBUSH,TVEG,TLOWER,TSKY,TCONDSB,RH,
     &   VEL,BP,ELEV,FASKY,FABUSH,FAVEG,FAGRD,QSLR/)
        TRAITS = (/TC,AK1,AK2,EMISAN,FATTHK,FLYHR,FURWET,PCTBAREVAP,
     &   PCTEYES/)

        !# set IPT, the geometry assumed in SIMULSOL: 1 = cylinder, 2 = sphere, 3 = ellipsoid
        if((SHAPE.eq.1.).or.(SHAPE.eq.3.).or.(SHAPE.eq.5.))THEN
         IPT = 1.
        ENDIF
        if(SHAPE.eq.2.)THEN
         IPT = 2.
        ENDIF
        If(SHAPE.eq.4.)THEN
         IPT = 3.
        ENDIF

        !# call SIMULSOL
        CALL SIMULSOL(DIFTOL, IPT, FURVARS, GEOMVARS, ENVVARS, TRAITS,
     &   TFA, SKINW, TS, SIMULOUT)
        SIMULSOLout(S,:) = SIMULOUT
1      CONTINUE
      
       TSKINMAX=max(SIMULSOLout(1,2), SIMULSOLout(2,2))
      
       !### ZBRENT and RESPFUN

       !# Now compute a weighted mean heat generation for all the parts/components = (dorsal value *(FASKY+FAVEG))+(ventral value*FAGRD)
       GEND = SIMULSOLout(1, 5)
       GENV = SIMULSOLout(2, 5)
       DMULT = FASKYREF + FAVEGREF
       VMULT = 1. - DMULT !# Assume that reflectivity of veg below = ref of soil so VMULT left as 1 - DMULT
       X = GEND * DMULT + GENV * VMULT !# weighted estimate of metabolic heat generation

       !# reset configuration factors
       FABUSH = FABUSHREF !# nearby bush
       FASKY = FASKYREF !# sky
       FAGRD = FAGRDREF !# ground
       FAVEG = FAVEGREF !# vegetation

       !# lung temperature and temperature of exhaled air
       TLUNG =(TC + (SIMULSOLout(1, 2) + SIMULSOLout(2, 2)) * 0.5) * 0.5 !# average of skin and core
       TAEXIT = min(TA + DELTAR, TLUNG) !# temperature of exhaled air, ?C

       !# now guess for metabolic rate that balances the heat budget while allowing metabolic rate
       !# to remain at or above QBASAL, via 'shooting method' ZBRENT
       QMIN = QBASAL
       IF((TA.LT.TC).AND.(TSKINMAX.LT.TC))THEN
        QM1 = QBASAL * 2.* (-1.)
        QM2 = QBASAL * 50.
       ELSE
        QM1 = QBASAL * 50.* (-1.)
        QM2 = QBASAL * 2.
       ENDIF

       QSUM = X
       TOL = AMASS * 0.01
       ZBRENTin = (/TA, O2GAS, N2GAS, CO2GAS, BP, QMIN, RQ, TLUNG,
     &  GMASS, EXTREF, RH, RELXIT, TIMACT, TAEXIT, QSUM, PANT, RP_CO2/)
      
       !# call ZBRENT subroutine which calls RESPFUN
       CALL ZBRENT_ENDO(QM1, QM2, TOL, ZBRENTin, ZBRENTout)
       !colnames(ZBRENTout) = c("RESPFN","QRESP","GEVAP", "PCTO2", "PCTN2", "PCTCO2", "RESPGEN", "O2STP", "O2MOL1", "N2MOL1", "AIRML1", "O2MOL2", "N2MOL2", "AIRML2", "AIRVOL")

       QGEN = ZBRENTout(7) ! Q_GEN,NET
       SHAPEB_LAST = SHAPEB
       AK1LAST = AK1
       TCLAST = TC
       PANTLAST = PANT
       SKINWLAST = SKINW
       if(SHAPEB.lt.SHAPEB_MAX)THEN
        SHAPEB = SHAPEB + UNCURL
       else
        SHAPEB = SHAPEB_MAX
        if(AK1.lt.AKMAX)THEN
         AK1 = AK1 + AK1inc
        else
         AK1 = AKMAX
         if(TC.lt.TCMAX)THEN
          TC = TC + RAISETC
          Q10mult = Q10**((TC - TCREF)/10.)
          QBASAL = QBASREF * Q10mult
         else
          TC = TCMAX
          Q10mult = Q10**((TC - TCREF)/10.)
          if(PANT.lt.PANTMAX)THEN
           PANT = PANT + PANTING
           PANTCOST=((PANT-1D0)/(PANTMAX-1D0)*PANTMULT*QBASREF)
           QBASAL = QBASREF * Q10mult + PANTCOST           
          else
           PANT = PANTMAX
           PANTCOST=((PANT-1D0)/(PANTMAX-1D0)*PANTMULT*QBASREF)
           QBASAL = QBASREF * Q10mult + PANTCOST           
           SKINW = SKINW + SWEAT
           if((SKINW.GT.MXWET).OR.(SWEAT.eq.0.))THEN
            SKINW = MXWET
            RETURN
           ENDIF
          ENDIF
         ENDIF
        ENDIF
       ENDIF
      END DO
      
      ! SIMULSOL output, dorsal
      TFAD=SIMULSOLout(1, 1) ! temperature of feathers/fur-air interface, deg C
      TSKCALCAVD=SIMULSOLout(1, 2) ! averagek skin temperature, deg C
      QCONVD=SIMULSOLout(1, 3) ! convection, W
      QCONDD=SIMULSOLout(1, 4) ! conduction, W
      QGENNETD=SIMULSOLout(1, 5) ! heat generation from flesh, W
      QSEVAPD=SIMULSOLout(1, 6) ! cutaneous evaporative heat loss, W
      QRADD=SIMULSOLout(1, 7) ! radiation lost at fur/feathers/bare skin, W
      QSLRD=SIMULSOLout(1, 8) ! solar radiation, W
      QRSKYD=SIMULSOLout(1, 9) ! sky radiation, W
      QRBSHD=SIMULSOLout(1, 10) ! bush/object radiation, W
      QRVEGD=SIMULSOLout(1, 11) ! overhead vegetation radiation (shade), W
      QRGRDD=SIMULSOLout(1, 12) ! ground radiation, W
      QFSEVAPD=SIMULSOLout(1, 13) ! fur evaporative heat loss, W
      NTRYD=SIMULSOLout(1, 14) ! solution attempts, #
      SUCCESSD=SIMULSOLout(1, 15) ! successful solution found? (0 no, 1 yes)

      ! SIMULSOL output, ventral
      TFAV=SIMULSOLout(2, 1) ! temperature of feathers/fur-air interface, deg C
      TSKCALCAVV=SIMULSOLout(2, 2) ! averagek skin temperature, deg C
      QCONVV=SIMULSOLout(2, 3) ! convection, W
      QCONDV=SIMULSOLout(2, 4) ! conduction, W
      QGENNETV=SIMULSOLout(2, 5) ! heat generation from flesh, W
      QSEVAPV=SIMULSOLout(2, 6) ! cutaneous evaporative heat loss, W
      QRADV=SIMULSOLout(2, 7) ! radiation lost at fur/feathers/bare skin, W
      QSLRV=SIMULSOLout(2, 8) ! solar radiation, W
      QRSKYV=SIMULSOLout(2, 9) ! sky radiation, W
      QRBSHV=SIMULSOLout(2, 10) ! bush/object radiation, W
      QRVEGV=SIMULSOLout(2, 11) ! overhead vegetation radiation (shade), W
      QRGRDV=SIMULSOLout(2, 12) ! ground radiation, W
      QFSEVAPV=SIMULSOLout(2, 13) ! fur evaporative heat loss, W
      NTRYV=SIMULSOLout(2, 14) ! solution attempts, #
      SUCCESSV=SIMULSOLout(2, 15) ! successful solution found? (0 no, 1 yes)

      RESPFN=ZBRENTout(1) ! heat sum (should be near zero), W
      QRESP=ZBRENTout(2) ! respiratory heat loss, W
      GEVAP=ZBRENTout(3) ! respiratory evaporation (g/s)
      PCTO2=ZBRENTout(4) ! O2 concentration (%)
      PCTN2=ZBRENTout(5) ! N2 concentration (%)
      PCTCO2=ZBRENTout(6) ! CO2 concentration (%)
      RESPGEN=ZBRENTout(7) ! metabolic heat (W)
      O2STP=ZBRENTout(8) ! O2 in rate at STP (L/s)
      O2MOL1=ZBRENTout(9) ! O2 in (mol/s)
      N2MOL1=ZBRENTout(10) ! N2 in (mol/s)
      AIRML1=ZBRENTout(11) ! air in (mol/s)
      O2MOL2=ZBRENTout(12) ! O2 out (mol/s)
      N2MOL2=ZBRENTout(13) ! N2 out (mol/s)
      AIRML2=ZBRENTout(14) ! air out (mol/s)
      AIRVOL=ZBRENTout(15) ! air out at STP (L/s)

      HTOVPR = 2.5012E+06 - 2.3787E+03 * TA
      SWEATGS = (SIMULSOLout(1,6) + SIMULSOLout(2,6)) * 0.5 
     &  / HTOVPR * 1000
      EVAPGS = ZBRENTout(3) + SWEATGS

      HTOVPR=2.5012E+06 - 2.3787E+03 * TA ! latent heat of vapourisation, W/kg/C
      SWEATGS=(QSEVAPD + QSEVAPV) * 0.5 / HTOVPR * 1000 ! water lost from skin, g/s
      EVAPGS=GEVAP + SWEATGS ! total evaporative water loss, g/s
      sigma=5.6697E-8
      QIR=QIRIN - QIROUT
      QIROUTD=sigma * EMISAN * AREASKIN * (TSKCALCAVD + 273.15)**4
      QIRIND=QRADD * (-1.) + QIROUTD
      QIROUTV=sigma * EMISAN * AREASKIN * (TSKCALCAVD + 273.15)**4
      QIRINV=QRADV * (-1.) + QIROUTV

      QSOL=QSLRD * DMULT + QSLRV * VMULT ! solar, W
      QIRIN=QIRIND * DMULT + QIRINV * VMULT ! infrared in, W
      QMET=RESPGEN ! metabolism, W
      QEVAP=QSEVAPD * DMULT + QSEVAPV * VMULT + QFSEVAPD * DMULT + 
     & QFSEVAPV * VMULT + QRESP ! evaporation, W
      QIROUT=QIROUTD * DMULT + QIROUTV * VMULT ! infrared out, W
      QCONV=QCONVD * DMULT + QCONVV * VMULT ! convection, W
      QCOND=QCONDD * DMULT + QCONDV * VMULT ! conduction, W

      TREG=(/TC, TLUNG, TSKCALCAVD, TSKCALCAVV, TFAD, TFAV, SHAPEB, 
     & PANT, SKINW, AK1, KEFARA(1), KEFARA(2), KEFARA(3), KFURCMPRS, 
     & Q10mult/)
      !names(treg)=c("TC", "TLUNG", "TSKIN_D", "TSKIN_V", "TFA_D", "TFA_V", "SHAPEB", "PANT", "SKINWET", "K_FLESH", "K_FUR", "K_FUR_D", "K_FUR_V", "K_COMPFUR", "Q10")

      MORPH=(/ATOT, VOL, D, MASFAT, FATTHK, FLSHVL, ALENTH, AWIDTH, 
     & AHEIT, R1, R2, ASIL, ASILN, ASILP, AREASKIN, CONVSK, CONVAR, 
     & AREACND/2., FASKY, FAGRD/)
      !names(morph)=c("AREA", "VOLUME", "CHAR_DIM", "MASS_FAT", "FAT_THICK", "FLESH_VOL", "LENGHT", "WIDTH", "HEIGHT", "DIAM_FLESH", "DIAM_FUR", "AREA_SIL", "AREA_SILN", "AREA_ASILP", "AREA_SKIN", "AREA_SKIN_EVAP", "AREA_CONV", "F_SKY", "F_GROUND")
      
      ENBAL=(/QSOL, QIRIN, QMET, QEVAP, QIROUT, QCONV, QCOND, RESPFN, 
     & max(NTRYD, NTRYV), min(SUCCESSD, SUCCESSV)/)
      !names(enbal)=c("QSOL", "QIRIN", "QMET", "QEVAP", "QIROUT", "QCONV", "QCOND", "ENB", "NTRY", "SUCCESS")

      MASBAL=(/AIRVOL, O2STP, GEVAP, SWEATGS, O2MOL1, 
     & O2MOL2, N2MOL1, N2MOL2, AIRML1, AIRML2/) * 3600
      !names(masbal)=c("AIR_L", "O2_L", "H2OResp_g", "H2OCut_g", "O2_mol_in", "O2_mol_out", "N2_mol_in", "N2_mol_out", "AIR_mol_in", "AIR_mol_out")

      RETURN
      END