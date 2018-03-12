      module commondat

      IMPLICIT None
      REAL, ALLOCATABLE, public :: SLES(:),RAIN(:),TIDES(:,:),metout(:,:
     &),shadmet(:,:),soil(:,:),shadsoil(:,:),soilmoist(:,:),shadmoist
     &(:,:),soilpot(:,:),shadpot(:,:),humid(:,:),shadhumid(:,:),
     &maxshades(:),minshades(:),CCMAXX(:),CCMINN(:),RHMAXX(:),RHMINN(:)
     &,WNMAXX(:),WNMINN(:),TMAXX(:),TMINN(:),TANNULRUN(:)
     &,REFLS(:),moists(:,:),intrvls(:),snowhr(:),nodes(:,:),TDSS(:),
     &TINS(:,:),TARS(:),RELS(:),CLDS(:),VELS(:),SOLS(:),ZENS(:),ZSLS(:),
     &PCTWET(:),julday(:),rainhr(:),LAIs(:),DRLAMBDA(:,:),DRRLAMBDA(:,:)
     &,SRLAMBDA(:,:),sunsnow(:,:),shdsnow(:,:),plant(:,:),shadplant(:,:)

      end module commondat

