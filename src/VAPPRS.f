      function vapprs(db)

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

      implicit none
      double precision loge,t,db,estar,vapprs
      t=db+273.16 
      if (t .le. 273.16)go to 20 
      loge=-7.90298*(373.16/t-1.)+5.02808*   
     x dlog10(373.16/t)-1.3816D-07 
     x *(10.**(11.344*(1.-t/373.16))-1.)+8.1328D-03 
     x *(10.**(-3.49149* 
     x (373.16/t-1.))-1.)+log10(1013.246)
      estar=(10.**loge)*100  
      go to 30   
20    loge=-9.09718*(273.16/t-1.)-3.56654*   
     x dlog10(273.16/t)+.876793*  
     x (1.-t/273.16)+log10(6.1071)   
      estar=(10.**loge)*100. 
30    continue   
      vapprs=estar   
      return 
      end
