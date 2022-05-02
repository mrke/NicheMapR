      SUBROUTINE DRYAIR(DB,BP,ALT,PATMOS,DENSTY,VISDYN,VISKIN,DIFVPR,   
     *THCOND,HTOVPR,TCOEFF,GGROUP)   
     
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
     
C***********************************************************************
C   
C SUBROUTINE DRYAIR CALCULATES SEVERAL PROPERTIES OF DRY AIR AND RELATED
C CHARACTERISTICS SHOWN AS OUTPUT VARIABLES BELOW.  THE PROGRAM IS BASED
C ON DATA FROM LIST, R. J. 1971. SMITHSONIAN METEOROLOGICAL TABLES. 
C SMITHSONIAN INSTITUTION PRESS. WASHINGTON, DC.
C   
C THE USER MUST SUPPLY VALUES FOR THE INPUT VARIABLES (DB, BP, AND ALT).
C IF ALT IS KNOWN (-1 000 < ALT < 20 000) BUT NOT BP THEN SET BP=0. 
C   
C*************************** INPUT VARIABLES ***************************
C   
C DB=DRY BULB TEMPERATURE (DEGREE CELSIUS)  
C BP=BAROMETRIC PRESSURE (PASCAL) [BP AT ONE STANDARD ATMOSPHERE IS 
C    101 325 PASCALS (100 PASCALS=1 MILLIBAR)]. 
C ALT=ALTITUDE (METRE) (1 METRE=3.280 839 9 FEET)   
C   
C*************************** OUTPUT VARIABLES **************************
C   
C PATMOS=STANDARD ATMOSPHERIC PRESSURE (PASCAL) 
C DENSTY=DENSITY (KILOGRAM PER CUBIC METRE) 
C VISDYN=DYNAMIC VISCOSITY (KILOGRAM PER METRE SECOND)  
C VISKIN=KINEMATIC VISCOSITY (SQUARE METRE PER SECOND)  
C DIFVPR=DIFFUSIVITY OF WATER VAPOR IN AIR (SQUARE METRE PER SECOND)
C THCOND=THERMAL CONDUCTIVITY (WATT PER METRE KELVIN)   
C HTOVPR=LATENT HEAT OF VAPORIZATION OF WATER (JOULE PER KILOGRAM)  
C TCOEFF=TEMPERATURE COEFFICIENT OF VOLUME EXPANSION (1 PER KELVIN) 
C GGROUP=GROUP OF VARIABLES IN GRASHOF NUMBER (1 PER CUBIC METRE KELVIN)
C   
C***********************************************************************
C   
      IMPLICIT LOGICAL (A-Z)
      DOUBLE PRECISION ALT,BP,C,DB,DENSTY,DIFVPR,GGROUP,HTOVPR,PATMOS
     *,PSTD,TCOEFF,THCOND,TNOT,TSTD,VISDYN,VISKIN,VISNOT  

      TSTD=273.15   
      PSTD=101325.  
      PATMOS=PSTD*((1.-(.0065*ALT/288.))**(1./.190284)) 
      IF(BP .LE. 0.)BP=PATMOS   
      DENSTY=BP/(287.04*(DB+TSTD))  
      VISNOT=1.8325D-5  
      TNOT=296.16   
      C=120.
      VISDYN=(VISNOT*((TNOT+C)/(DB+TSTD+C)))*(((DB+TSTD)/TNOT)**1.5)
      VISKIN=VISDYN/DENSTY  
      DIFVPR=2.26D-5*(((DB+TSTD)/TSTD)**1.81)*(1.D5/BP) 
      THCOND=.02425+(7.038D-5*DB)   
      HTOVPR=2.5012D6-2.3787D3*DB   
      TCOEFF=1./(DB+TSTD)   
      GGROUP=.0980616*TCOEFF/(VISKIN*VISKIN)
      RETURN
      END
