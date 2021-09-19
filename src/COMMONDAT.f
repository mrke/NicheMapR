      module commondat
          
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

      IMPLICIT None
      DOUBLE PRECISION, ALLOCATABLE, public :: SLES(:),RAIN(:)
     &,TIDES(:,:),metout(:,:),snowhr(:),
     &shadmet(:,:),soil(:,:),shadsoil(:,:),soilmoist(:,:),shadmoist
     &(:,:),soilpot(:,:),shadpot(:,:),humid(:,:),shadhumid(:,:),
     &maxshades(:),minshades(:),CCMAXX(:),CCMINN(:),RHMAXX(:),RHMINN(:)
     &,WNMAXX(:),WNMINN(:),TMAXX(:),TMINN(:),TANNULRUN(:)
     &,REFLS(:),moists(:,:),intrvls(:),nodes(:,:),TDSS(:),
     &TINS(:,:),TARS(:),RELS(:),CLDS(:),VELS(:),SOLS(:),ZENS(:),ZSLS(:),
     &PCTWET(:),julday(:),rainhr(:),LAIs(:),DRLAMBDA(:,:),DRRLAMBDA(:,:)
     &,SRLAMBDA(:,:),sunsnow(:,:),shdsnow(:,:),plant(:,:),shadplant(:,:)
     &,tcond(:,:),shadtcond(:,:),specheat(:,:),shadspecheat(:,:)
     &,densit(:,:),shaddensit(:,:)
      end module commondat

