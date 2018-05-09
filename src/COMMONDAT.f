      module commondat

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
      end module commondat

