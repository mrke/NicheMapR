      program microclim

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

c     Main program for non-R testing. Reads in all the files produced by
c     the flag 'write_input' in the micro_global, micro_aust, micro_usa,
c     micro_nz, micro_uk etc. functions

      IMPLICIT NONE

      double precision, DIMENSION(:,:), ALLOCATABLE :: NODES2,metout2
     &,shadmet2,soil2,shadsoil2,soilmoist2,shadmoist2,humid2,shadhumid2,
     &soilpot2,shadpot2,moists2,tides2,DRLAMBDA2,DRRLAMBDA2,SRLAMBDA2,
     &sunsnow2,shdsnow2,plant2,shadplant2,tcond2,shadtcond2,specheat2,
     &shadspecheat2,densit2,shaddensit2
      double precision, DIMENSION(:), ALLOCATABLE :: julday2,LAI2
     &,SLES2,MAXSHADES2,MINSHADES2,RHMAXX2,RHMINN2,CCMAXX2
     &,CCMINN2,WNMAXX2,WNMINN2,REFLS2,PCTWET2
     &,TMAXX2,TMINN2,rain2,tannulrun2,TAIRhr2,RHhr2,WNhr2,CLDhr2,IRhr2
     &,SOLRhr2,RAINhr2,ZENhr2
      double precision soilprop2,microinput2,PE2,BD2,BB2
     &,KS2,L2,hori2,tai2,soilinit2,dep2,DD2

      INTEGER I,J,nn3

      CHARACTER(80) LABEL

      DIMENSION tai2(111),soilinit2(20),hori2(24),
     & Dep2(10),L2(19),microinput2(75),soilprop2(10,5)
      DIMENSION PE2(19),KS2(19),BD2(19),BB2(19),DD2(19)

      OPEN(1,FILE='microinput.csv')
      read(1,*)LABEL
      DO 11 i=1,75
      read(1,*)label,microinput2(i)
11    continue
      close(1)

      nn3=int(microinput2(1))

      allocate (moists2(10,nn3),Nodes2(10,nn3),CCMAXX2(nn3),CCMINN2(nn3)
     &,RHMAXX2(nn3),RHMINN2(nn3),TMINN2(nn3),TMAXX2(nn3),WNMAXX2(nn3),
     &WNMINN2(nn3),REFLS2(nn3),PCTWET2(nn3),SLES2(nn3),
     &MAXSHADES2(nn3),MINSHADES2(nn3),JULDAY2(nn3),LAI2(nn3),
     &rain2(nn3),tannulrun2(nn3),METOUT2(24*nn3,19),SHADMET2(24*nn3,19),
     &SOIL2(24*nn3,12),SHADSOIL2(24*nn3,12),tides2(24*nn3,3),
     &soilpot2(24*nn3,12),shadpot2(24*nn3,12),humid2(24*nn3,12),
     &shadhumid2(24*nn3,12),soilmoist2(24*nn3,12),shadmoist2(24*nn3,12)
     &,sunsnow2(24*nn3,11),shdsnow2(24*nn3,11),plant2(24*nn3,14),
     &shadplant2(24*nn3,14),tcond2(24*nn3,12),shadtcond2(24*nn3,12),
     &specheat2(24*nn3,12),shadspecheat2(24*nn3,12),densit2(24*nn3,12),
     &shaddensit2(24*nn3,14),TAIRhr2(24*nn3),RHhr2(24*nn3),WNhr2(24*nn3)
     &,CLDhr2(24*nn3),IRhr2(24*nn3),SOLRhr2(24*nn3),RAINhr2(24*nn3)
     &,ZENhr2(24*nn3),DRLAMBDA2(24*nn3,113),DRRLAMBDA2(24*nn3,113)
     &,SRLAMBDA2(24*nn3,113))

      OPEN(1,FILE='doy.csv')
      read(1,*)LABEL
      do 13 i=1,nn3
      read(1,*)label,julday2(i)
13    continue
      close(1)

      OPEN(1,FILE='SLES.csv')
      read(1,*)LABEL
      do 14 i=1,nn3
      read(1,*)label,SLES2(i)
14    continue
      close(1)

      OPEN(1,FILE='DEP.csv')
      read(1,*)LABEL
      do 15 i=1,10
      read(1,*)label,DEP2(i)
15    continue
      close(1)

      OPEN(1,FILE='Nodes.csv')
      read(1,*)LABEL
      do 20 i=1,10
      read(1,*)label,(Nodes2(i,j),j=1,nn3)
20    continue
      close(1)

      OPEN(1,FILE='Maxshades.csv')
      read(1,*)LABEL
      do 21 i=1,nn3
      read(1,*)label,Maxshades2(i)
21    continue
      close(1)

      OPEN(1,FILE='Minshades.csv')
      read(1,*)LABEL
      do 22 i=1,nn3
      read(1,*)label,Minshades2(i)
22    continue
      close(1)

C     OPEN(1,FILE='TIMAXS.csv')
C     read(1,*)LABEL
C     do 23 i=1,4
C     read(1,*)label,TIMAXS2(i)
C23    continue
C     close(1)
C
C     OPEN(1,FILE='TIMINS.csv')
C     read(1,*)LABEL
C     do 24 i=1,4
C     read(1,*)label,TIMINS2(i)
C24    continue
C     close(1)

      OPEN(1,FILE='TMAXX.csv')
      read(1,*)LABEL
      do 25 i=1,nn3
      read(1,*)label,TMAXX2(i)
25    continue
      close(1)

      OPEN(1,FILE='TMINN.csv')
      read(1,*)LABEL
      do 26 i=1,nn3
      read(1,*)label,TMINN2(i)
26    continue
      close(1)

      OPEN(1,FILE='RHMAXX.csv')
      read(1,*)LABEL
      do 27 i=1,nn3
      read(1,*)label,RHMAXX2(i)
27    continue
      close(1)

      OPEN(1,FILE='RHMINN.csv')
      read(1,*)LABEL
      do 28 i=1,nn3
      read(1,*)label,RHMINN2(i)
28    continue
      close(1)

      OPEN(1,FILE='CCMAXX.csv')
      read(1,*)LABEL
      do 29 i=1,nn3
      read(1,*)label,CCMAXX2(i)
29    continue
      close(1)

      OPEN(1,FILE='CCMINN.csv')
      read(1,*)LABEL
      do 30 i=1,nn3
      read(1,*)label,CCMINN2(i)
30    continue
      close(1)

      OPEN(1,FILE='WNMAXX.csv')
      read(1,*)LABEL
      do 31 i=1,nn3
      read(1,*)label,WNMAXX2(i)
31    continue
      close(1)

      OPEN(1,FILE='WNMINN.csv')
      read(1,*)LABEL
      do 32 i=1,nn3
      read(1,*)label,WNMINN2(i)
32    continue
      close(1)

      OPEN(1,FILE='REFLS.csv')
      read(1,*)LABEL
      do 34 i=1,nn3
      read(1,*)label,REFLS2(i)
34    continue
      close(1)

      OPEN(1,FILE='PCTWET.csv')
      read(1,*)LABEL
      do 35 i=1,nn3
      read(1,*)label,PCTWET2(i)
35    continue
      close(1)

      OPEN(1,FILE='soilinit.csv')
      read(1,*)LABEL
      do 36 i=1,20
      read(1,*)label,soilinit2(i)
36    continue
      close(1)

      OPEN(1,FILE='hori.csv')
      read(1,*)LABEL
      do 37 i=1,24
      read(1,*)label,hori2(i)
37    continue
      close(1)

      OPEN(1,FILE='TAI.csv')
      read(1,*)LABEL
      do 38 i=1,111
      read(1,*)label,TAI2(i)
38    continue
      close(1)

      OPEN(1,FILE='soilprop.csv')
      read(1,*)LABEL
      do 39 i=1,10
      read(1,*)label,(soilprop2(i,j),j=1,5)
39    continue
      close(1)

      OPEN(1,FILE='moists.csv')
      read(1,*)LABEL
      do 40 i=1,10
      read(1,*)label,(moists2(i,j),j=1,nn3)
40    continue
      close(1)

      OPEN(1,FILE='rain.csv')
      read(1,*)LABEL
      do 41 i=1,nn3
      read(1,*)label,rain2(i)
41    continue
      close(1)

      OPEN(1,FILE='tannulrun.csv')
      read(1,*)LABEL
      do 42 i=1,nn3
      read(1,*)label,tannulrun2(i)
42    continue
      close(1)

      OPEN(1,FILE='PE.csv')
      read(1,*)LABEL
      do 43 i=1,19
      read(1,*)label,PE2(i)
43    continue
      close(1)

      OPEN(1,FILE='KS.csv')
      read(1,*)LABEL
      do 44 i=1,19
      read(1,*)label,KS2(i)
44    continue
      close(1)

      OPEN(1,FILE='BB.csv')
      read(1,*)LABEL
      do 45 i=1,19
      read(1,*)label,BB2(i)
45    continue
      close(1)

      OPEN(1,FILE='BD.csv')
      read(1,*)LABEL
      do 46 i=1,19
      read(1,*)label,BD2(i)
46    continue
      close(1)

      OPEN(1,FILE='DD.csv')
      read(1,*)LABEL
      do 58 i=1,19
      read(1,*)label,DD2(i)
58    continue
      close(1)

      OPEN(1,FILE='L.csv')
      read(1,*)LABEL
      do 47 i=1,19
      read(1,*)label,L2(i)
47    continue
      close(1)

      OPEN(1,FILE='LAI.csv')
      read(1,*)LABEL
      do 48 i=1,(nn3)
      read(1,*)label,LAI2(i)
48    continue
      close(1)

      OPEN(1,FILE='tides.csv')
      read(1,*)LABEL
      do 49 i=1,(nn3*24)
      read(1,*)label,(tides2(i,j),j=1,3)
49    continue
      close(1)

      OPEN(1,FILE='TAIRhr.csv')
      read(1,*)LABEL
      do 50 i=1,nn3*24
      read(1,*)label,TAIRhr2(i)
50    continue
      close(1)

      OPEN(1,FILE='RHhr.csv')
      read(1,*)LABEL
      do 51 i=1,nn3*24
      read(1,*)label,RHhr2(i)
51    continue
      close(1)

      OPEN(1,FILE='WNhr.csv')
      read(1,*)LABEL
      do 52 i=1,nn3*24
      read(1,*)label,WNhr2(i)
52    continue
      close(1)

      OPEN(1,FILE='CLDhr.csv')
      read(1,*)LABEL
      do 53 i=1,nn3*24
      read(1,*)label,CLDhr2(i)
53    continue
      close(1)

      OPEN(1,FILE='SOLRhr.csv')
      read(1,*)LABEL
      do 55 i=1,nn3*24
      read(1,*)label,SOLRhr2(i)
55    continue
      close(1)

      OPEN(1,FILE='RAINhr.csv')
      read(1,*)LABEL
      do 56 i=1,nn3*24
      read(1,*)label,RAINhr2(i)
56    continue
      close(1)

      OPEN(1,FILE='ZENhr.csv')
      read(1,*)LABEL
      do 57 i=1,nn3*24
      read(1,*)label,ZENhr2(i)
57    continue
      close(1)

      OPEN(1,FILE='IRDhr.csv')
      read(1,*)LABEL
      do 59 i=1,nn3*24
      read(1,*)label,IRhr2(i)
59    continue
      close(1)

      call microclimate(nn3,microinput2,julday2,SLES2,DEP2,
     &maxshades2,minshades2,Nodes2,
     &RHMAXX2,RHMINN2,CCMAXX2,CCMINN2,WNMAXX2,WNMINN2,TMAXX2,TMINN2
     &,REFLS2,PCTWET2,soilinit2,hori2,tai2,soilprop2,moists2,
     &rain2,tannulrun2,tides2,PE2,KS2,BB2,BD2,DD2,L2,LAI2,TAIRhr2
     &,RHhr2,WNhr2,CLDhr2,SOLRhr2,RAINhr2,ZENhr2,IRhr2,metout2,soil2
     &,shadmet2,shadsoil2,soilmoist2,shadmoist2,humid2,shadhumid2
     &,soilpot2,shadpot2,sunsnow2,shdsnow2,plant2,shadplant2,tcond2,
     &shadtcond2,specheat2,shadspecheat2,densit2,shaddensit2,DRLAMBDA2
     &,DRRLAMBDA2,SRLAMBDA2)
      end
