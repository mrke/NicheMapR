      subroutine soilprops(TSOI,ALTT,soilprop,moistt)

      use commondat
      IMPLICIT NONE
      EXTERNAL WETAIR
      
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

C     Computes variable thermal conductivity with temperature and water content
c     based on Campbell, G. S., J. D. J. Jungbauer, W. R. Bidlake, and R. D. Hungerford. 1994.
c     Predicting the effect of temperature on soil thermal conductivity. Soil Science 158:307-313.

      double precision TSOI,ALTT,drydensity,bar,moistt,Thconduct,Density
      double precision p_a0,PATMOS,BP,T_K,k_a,k_w,k_m,p_a,hr,D_v0
      double precision rho_hat,lambda,Spheat,rho_hat0,D_v,k_f,g_a,g_c
      double precision WB,RH,DP,ESAT,VD,RW,TVIR,TVINC,DENAIR,CP,WTRPOT,q
      double precision e_a,E,deltax,theta_0,phi_m,theta,phi_g,f_w,k_g
      double precision epsilon_g,epsilon_w,epsilon_m,e_a1,e_a2
      double precision soilprop,dens,spht,HTOFN
      double precision condep,rainmult,maxpool
      double precision DENDAY,SPDAY,TKDAY,minsnow,snowcond2
      double precision snownode,maxsnode1,snode,daysincesnow,
     &lastday,undercatch,rainmeltf,densfun,AZMUTH,SLOPE
      double precision PUNSH,ALAT,AMULT,PRESS,CMH2O,REFL,ALONC,TIMCOR
     & ,TSNHR,TSRHR,HEMIS,snowcond,intercept
      double precision rww,pc,rl,sp,r1,im
      double precision snowdens,snowmelt,snowtemp,cursnow,
     & snowage,prevden,cpsnow,grasshade
     
      INTEGER I,J,II,maxcount
      INTEGER JULNUM,DOY,Numtyps
      INTEGER NON,evenrain,runmoist,runsnow,trouble
      INTEGER I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,ij,I97,I98,I99,I100,I101    
      INTEGER IPINT,NOSCAT,IUV,IALT,IDAYST,IDA,IEP,ISTART,IEND2
      DIMENSION DENDAY(30),SPDAY(30),TKDAY(30)
      dimension soilprop(10,5),TSOI(30),moistt(10)
      DIMENSION Thconduct(30),Density(30),Spheat(30)
      dimension drydensity(10),bar(10),k_m(10)
      dimension spht(10),dens(10),snownode(10),snode(10),densfun(4)

      COMMON/SOYVAR1/Numtyps
      COMMON/SOYVAR2/Thconduct,Density,Spheat
      COMMON/SOYFILS/DENDAY,SPDAY,TKDAY
      COMMON/SOILND/NON
      COMMON/DAYJUL/JULNUM,DOY
      COMMON/WMAIN/I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,I11,I12,I91,I92,I93
     & ,I94,I95,I96,I97,I98,I99,I100,I101
      common/soilmoist/condep,rainmult,maxpool
      common/soilmoist3/runmoist,evenrain,maxcount 
      common/soimoist2/rww,pc,rl,sp,r1,im
      common/snowmod/runsnow,trouble
      COMMON/SNOWPRED/snowtemp,snowdens,snowmelt,snownode,minsnow
     &,maxsnode1,snode,cursnow,daysincesnow,lastday,undercatch,rainmeltf
     &,densfun,snowcond,intercept,snowage,prevden,grasshade
      COMMON/WIOCONS/PUNSH,ALAT,AMULT,PRESS,CMH2O,REFL,ALONC,TIMCOR,
     * AZMUTH,SLOPE,TSNHR,TSRHR,Hemis
      COMMON/WIOCONS2/IPINT,NOSCAT,IUV,IALT,IDAYST,IDA,IEP,ISTART,IEND2
     
      HTOFN=333500. !J/kg

      do 3 j=1,numtyps
      drydensity(j)=soilprop(j,1)
      bar(j)=soilprop(j,2)
      k_m(j)=soilprop(j,3)
      spht(j)=soilprop(j,4)
      dens(j)=soilprop(j,5)
3     continue
      if(runsnow.eq.1)then
       ii=9
       ij=8
      else
       ii=1
       ij=0
      endif
      j=1
      do 2 i=ii,NON ! ii ensures that it starts at soil depth if snow on top 
       if(i.ge.nodes(j,DOY)+ij)then
        j=j+1
        if(j.gt.numtyps)then
         j=numtyps
        endif
       endif
c     don't make it volumetric, but rather mass-specific, so don't multiply by kg/m3 converters (constants from Campbell and Norman 1998, Table 8.2)
       if(runmoist.eq.1)then
        if(runsnow.eq.1)then
         Spheat(i)=drydensity(j)/dens(j)*spht(j)+moistt(i-8)*4184.
        else
         Spheat(i)=drydensity(j)/dens(j)*spht(j)+moistt(i)*4184.
        endif
       else
        Spheat(i)=drydensity(j)/dens(j)*spht(j)+bar(j)*moistt(j)*4184.
       endif
c     constants from Campbell and Norman 1998, Table 8.2
       if(runmoist.eq.1)then
        if(runsnow.eq.1)then
         Density(i)=moistt(i-8)*1000+drydensity(j)/dens(j)*dens(j)*1000
        else
         Density(i)=moistt(i)*1000+drydensity(j)/dens(j)*dens(j)*1000
        endif
       else
        Density(i)=bar(j)*moistt(j)*1000+drydensity(j)/dens(j)*
     &    dens(j)*1000
       endif
c    # standard sea level air pressure, Pa
       p_a0=101325
       PATMOS=p_a0*((1.-(0.0065*ALTT/288.))**(1./0.190284))
       BP = PATMOS

c    # deg K
       T_K=TSOI(i)+273.15
c    # W/mC thermal conductivity of dry air (equation 9 in Campbell et al. 1994)
       k_a=0.024+7.73e-5*TSOI(i)-2.6e-8*TSOI(i)**2
c    # W/mC thermal conductivity of water (equation 8 in Campbell et al. 1994)
       k_w=0.554+2.24e-3*TSOI(i)-9.87e-6*TSOI(i)**2
c    # W/mC thermal conductivity of minerals (value of 2.5 suggested from Campbell and Norman 1998, Table 8.2)
c      k_m= 2.5


c    # air pressure
       p_a=BP
c    # relative humidity
       hr=1.0
c    # vapour diffusivity in air (m2/s), standard value at 0 deg C and sea level pressure (Campbell et al. 1994)
       D_v0=2.12e-5
c    # molar density of air (mol/m3), standard value at 0 deg C and sea level pressure (Campbell et al. 1994)
       rho_hat0=44.65
c    # temperature/pressure-corrected vapour diffusivity in air (m2/s) (p. 309 in Campbell et al. 1994)
       D_v=D_v0*(p_a0/p_a)*(T_K/273.15)**1.75
c    # temperature/pressure-corrected molar density of air (mol/m3) (p. 309 in Campbell et al. 1994)
       rho_hat=rho_hat0*(p_a/p_a0)*(273.15/T_K)
c    # J/mol latent heat of vaporization (Cambell et al. 1994, p. 309)
       lambda=45144-48*TSOI(i)
c      # vapour pressure at three temps, from which slope of vapour pressure function is calculated, assuming near 100% RH

       RH = 99.
       WB = 0.
       DP = 999.

       CALL WETAIR (TSOI(i),WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,
     &      CP,WTRPOT)
       e_a=E
       CALL WETAIR (TSOI(i)-1,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR
     &      ,CP,WTRPOT)

       e_a1=E
       CALL WETAIR (TSOI(i)+1,WB,RH,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR
     &      ,CP,WTRPOT)
      e_a2=E

c    # slope of the vapour pressure function centred at focal temperature
       deltax=(e_a2-e_a1)/2

c     # these could vary with soil texture but the relationship isn't strong
c     # power for liquid recirculation, mean in Table 2 of Campell et al. 1994, excluding peat moss value
c      q_0=4
c      q=q_0*(T_K/303.)**2.
c     using a typical value for 'q' of 4 - program becomes unstable if this is temperature dependent
       q=4.

c    # mean in Table 2 of Campell et al. 1994, excluding peat moss value
       theta_0=0.162

c    # volume fraction of minerals
       phi_m=drydensity(j)/dens(j)
c    # volume fraction of water
       if(runmoist.eq.1)then
        if(runsnow.eq.1)then
         theta=moistt(i-8)
        else
         theta=moistt(i)
        endif
       else
        theta=moistt(j)*bar(j)
       endif
c     # # volume fraction of gas
       phi_g=1-theta-phi_m
       if(phi_g.lt.0)then
        phi_g=0
       endif
c     # eq 8.17 Campbell and Norman 1988
       f_w=1/(1+(theta/theta_0)**(-4.))
c     # eq 8.17 Campbell and Norman 1988, using temperature-specific q
       f_w=1/(1+(theta/theta_0)**(-1.*q))
c    # eq 8.18 Campbell and Norman 1988
       k_g=k_a+lambda*deltax*hr*f_w*rho_hat*D_v/(p_a-e_a)
c     # eq 8.19 Campbell and Norman 1988
       k_f=k_g+f_w*(k_w-k_g)
c    # 0.1 for mineral soils, 0.33 for organic, p 125, Campbell and Norman 1988
       g_a=0.1
c     # p 125, Campbell and Norman
       g_c=1-2*g_a
c     # equation 8.20 in Campbell and Norman 1988
       epsilon_g=2/(3*(1+g_a*(k_g/k_f-1)))+1/(3*(1+g_c*(k_g/k_f-1)))
c    # equation 8.20 in Campbell and Norman 1988
       epsilon_w=2/(3*(1+g_a*(k_w/k_f-1)))+1/(3*(1+g_c*(k_w/k_f-1)))
c    # equation 8.20 in Campbell and Norman 1988
       epsilon_m=2/(3*(1+g_a*(k_m(j)/k_f-1)))+
     &1/(3*(1+g_c*(k_m(j)/k_f-1)))
c    # equation 8.13 in Campbell and Norman 1988
       Thconduct(i)=(theta*epsilon_w*k_w+phi_m*epsilon_m*k_m(j)+phi_g*
     &epsilon_g*k_g)/(theta*epsilon_w+phi_m*epsilon_m+phi_g*epsilon_g)

c     Convert thermal conductivities from W/m-K to cal/min-cm-K for DSUB'S Microclimate calculations
       Thconduct(i)=(Thconduct(i)/418.6)*60.
c     Convert specific heats from J/kg-K to cal/g-K for DSUB'S Microclimate calculations
       Spheat(i)=Spheat(i)/4186.
c     Convert densities from kg/m3 to g/cm3 for DSUB'S Microclimate calculations
C     1 kg/m3 * 1000g/1 kg * 1 m3/1000000.
       Density(i)=Density(i)/1.0E+3
2     continue

      if(runsnow.eq.1)then
       if(densfun(1).gt.0)then ! compute snowdens from linear or exponential function
         if(densfun(2).gt.0)then ! exponential function
           snowdens=(densfun(1)-densfun(2))*(1-EXP(-1*densfun(3)*cursnow
     &      -densfun(4)*snowage))+densfun(2)
         else ! linear function
          snowdens=min(0.9167D+0,densfun(1)*snowage+densfun(2))
         endif
       endif
       if(cursnow.ge.minsnow)then ! snow is present
        call WETAIR(0,WB,100,DP,BP,E,ESAT,VD,RW,TVIR,TVINC,DENAIR,CP,
     &  WTRPOT) ! get specific heat and mixing ratio of humid air at zero C
        cpsnow = (2100*snowdens+(1.005+1.82*(RW/1.+RW))*1000* ! based on https://en.wiktionary.org/wiki/humid_heat
     &   (1-snowdens)) ! compute weighted specific heat accounting for ice vs airm SI units
        snowcond2 = (0.00395+0.00084*(snowdens*1000)-0.0000017756* ! snow thermal conductivity as a function of density (from Aggarwal, R. 2009. Defence Science Journal 59:126–130.)
     &  (snowdens*1000)**2+0.00000000380635*(snowdens*1000)**3)
     & /418.6*60
        do 4 i=1,8 ! from top node down - top only has snow if pack has built up to higher than snow node 7 (2m)
         ! first give all nodes the conductivity and spheat of snow
         if(snowcond.gt.0)then ! check if using user-defined thermal conductivity for snow
          Thconduct(i)=snowcond ! user-defined
         else
          Thconduct(i)=snowcond2 ! Aggarwal function
         endif
         if((TSOI(i).gt.-0.45).and.(TSOI(i).le.0.4))then 
          Spheat(i)=(cpsnow+HTOFN)/4186. ! freezing is occuring, add heat of fusion
         else
          Spheat(i)=cpsnow/4186. ! freezing is not occuring, do not add heat of fusion
         endif
         ! now only give layers with snow the density of snow, otherwise zero (which means they have no influence)
         if((snode(i).gt.0).or.
     &    (((snode(i).lt.1e-8).and.(snode(min(8,i+1)).gt.0))))then
          Density(i)=snowdens
         else
          Density(i)=0.D0
         endif
4       continue
       else
        do 5 i=1,8 ! no snow, give all snow nodes (top 8) the value of the soil layer for conductivity and specific heat and zero for density
         Thconduct(i)=Thconduct(9)
         Spheat(i)=Spheat(9)
         Density(i)=0.D0
5       continue
       endif
      endif
      RETURN
      END
