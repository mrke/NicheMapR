C     NICHEMAPR: SOFTWARE FOR BIOPHYSICAL MECHANISTIC NICHE MODELLING

C     COPYRIGHT (C) 2020 MICHAEL R. KEARNEY AND WARREN P. PORTER

C     THIS PROGRAM IS FREE SOFTWARE: YOU CAN REDISTRIBUTE IT AND/OR MODIFY
C     IT UNDER THE TERMS OF THE GNU GENERAL PUBLIC LICENSE AS PUBLISHED BY
C     THE FREE SOFTWARE FOUNDATION, EITHER VERSION 3 OF THE LICENSE, OR (AT
C      YOUR OPTION) ANY LATER VERSION.

C     THIS PROGRAM IS DISTRIBUTED IN THE HOPE THAT IT WILL BE USEFUL, BUT
C     WITHOUT ANY WARRANTY; WITHOUT EVEN THE IMPLIED WARRANTY OF
C     MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. SEE THE GNU
C     GENERAL PUBLIC LICENSE FOR MORE DETAILS.

C     YOU SHOULD HAVE RECEIVED A COPY OF THE GNU GENERAL PUBLIC LICENSE
C     ALONG WITH THIS PROGRAM. IF NOT, SEE HTTP://WWW.GNU.ORG/LICENSES/.

C     SUBROUTINE TO COMPUTE SHAPE PROPORTIONS FOR A HUMAN

C     This function finds the shape_b values for each body part as well as
C     the fraction of each body part connected to other body parts, given
C     the known total surface area and total mass, and the fractions of
C     surface area and mass that each body part contributes to the whole body.
C     It essentially breaks the body into separate parts, works out their
C     shapes, then puts them back together to make sure the final sum of
C     areas matches the given area. It reports back the shape_b proportionality
C     constants and the fraction of each body part connected to the rest of the
C     body. The latter fractions are assumed to be undergoing conductive heat
C     exchange to a substrate at core temperature which effectively removes
C     these joining surface areas from any forms of heat exchange

      subroutine get_shapes(MASSs,HEIGHT,AREA,SHAPE_Bs,SHAPE_Bs_min,
     & SHAPE_Bs_max,AREAFRACs,DENSITYs,FATPCTs,SHAPEs,SUBQFATs,DHAIRDs,
     & DHAIRVs,INSDENDs,INSDENVs,tol,maxiter,RESULTS)

      implicit none

      DOUBLE PRECISION MASSs(4),AREA,AREAFRACs(4),DENSITYs(4),FATPCTs(4)
      DOUBLE PRECISION SHAPE_Bs(4),SHAPE_Bs_min(4),SHAPE_Bs_max(4)
      DOUBLE PRECISION SHAPEs(4),SUBQFATs(4),DHAIRDs(4),DHAIRVs(4)
      DOUBLE PRECISION INSDENDs(4),INSDENVs(4),RESULTS(10)
      DOUBLE PRECISION SA_targets(4),SA_currents(4),PJOINs(4)
      DOUBLE PRECISION GEOM_out(22),PART_END(4),PART_TOTAL(4)
      DOUBLE PRECISION AREA_orig,tol,maxiter,stepcount,HEIGHT
      DOUBLE PRECISION HEIGHT_CUR,HEADL,TRUNKL,LEGL,AREA_CUR
      INTEGER part_order(4),i,ii,num_parts
 
      ! Initialize variables
      num_parts = 4
      
      SA_targets = 0.0
      SA_currents = 0.0
      PJOINs = 0.0
      PART_END = 0.0
      PART_TOTAL = 0.0
      HEIGHT_CUR=0.0
      HEADL=0.0
      TRUNKL=0.0
      LEGL=0.0
      
      AREA_orig = AREA
      ! Compute target surface areas
      do i=1,num_parts
      SA_targets(i)=AREA_orig*AREAFRACs(i)
      end do
      DATA part_order/3,4,1,2/ ! arm, leg, head, torso
      stepcount = 0.0
C     call GEOM_ENDO(MASSs(1),DENSITYs(1),DENSITYs(1),FATPCTs(1),SHAPEs(1),
C    &       0D1, SUBQFATs(1),SHAPE_Bs(1),SHAPE_Bs(1),
C    &       (DHAIRDs(1)+DHAIRVs(1))/2D1,(INSDENDs(1)+INSDENVs(1))
C    &       /2D1,PJOINs(1),0D1,0D1,0D1,GEOM_out)
C     open(unit=11, file="debug_output1.txt", status="replace")
C     write(11, *) "Iteration:", stepcount
C     write(11, *) "tol:", tol
C     write(11, *) "maxiter:", maxiter
C     write(11, *) "SHAPE_Bs:", SHAPE_Bs
C     write(11, *) "PJOINs:", PJOINs
C     write(11, *) "SA_targets:", SA_targets
C     write(11, *) "SA_currents:", SA_currents
C     write(11, *) "Sum of Squared Differences:", 
C    * sum((SA_targets - SA_currents)**2.0)
C     write(10, *) "GEOM_out:", GEOM_out
C     close(11)
      do while ((sum((SA_targets-SA_currents)**2.0)) .gt. tol 
C     & .or.((HEIGHT-HEIGHT_CUR)**2.0).GT.tol 
     & .and. stepcount .lt. maxiter)
      if(stepcount .gt. maxiter)then
          print *, 'Warning: Hit maximum iterations on shape adjustment'
      end if

      ! Adjust SHAPE_Bs if needed
      if(stepcount.gt.0)then
          do i=1,num_parts
            !if(((SA_targets(i)-SA_currents(i))**2.0).GT.tol)then
            if(((SA_targets(i)-SA_currents(i))**2.0).GT.tol
C     &       .or.((HEIGHT-HEIGHT_CUR)**2.0).GT.tol
     &        )then
              if(SHAPE_Bs(i)<SHAPE_Bs_max(i).and.
     &         SHAPE_Bs(i)>SHAPE_Bs_min(i)) then
                if((sum(SA_targets)>sum(SA_currents)))then ! current surface area too low
                  SHAPE_Bs(i)=SHAPE_Bs(i)*1.001
                else
                  SHAPE_Bs(i)=SHAPE_Bs(i)*0.999
                end if
              end if
            end if
          end do
      end if
C     if(HEIGHT+0.05>HEIGHT_CUR)then ! current height too low, increase trunk and legs a bit
C       SHAPE_Bs(2)=SHAPE_Bs(2)*1.001
C       SHAPE_Bs(4)=SHAPE_Bs(4)*1.001
C     end if
C     if(HEIGHT-0.05<HEIGHT_CUR)then ! current height too low, increase trunk and legs a bit
C       SHAPE_Bs(2)=SHAPE_Bs(2)*0.999
C       SHAPE_Bs(4)=SHAPE_Bs(4)*0.999
C     end if
      
      ! Compute geometry and update PJOINs
      do i = 1, num_parts
          ii = part_order(i) 

          call GEOM_ENDO(MASSs(ii),DENSITYs(ii),DENSITYs(ii),FATPCTs(ii)
     &       ,SHAPEs(ii),0D1, SUBQFATs(ii),SHAPE_Bs(ii),SHAPE_Bs(ii),
     &       (DHAIRDs(ii)+DHAIRVs(ii))/2D1,(INSDENDs(ii)+INSDENVs(ii))
     &       /2D1,PJOINs(ii),0D1,0D1,0D1,GEOM_out)

          if(i==1.or.i==2)then ! Arm or leg
            part_end(ii)=3.14159*(GEOM_out(6)/2.0)**2.0
            part_total(ii)=2.0*3.14159*(GEOM_out(6)/2.0)**2.0+
     &                    3.14159*GEOM_out(6)*GEOM_out(5)
            PJOINs(ii)=part_end(ii)/part_total(ii)
            if(i==2)then !leg
              LEGL=GEOM_out(5)
            endif
          else if (i==3) then ! Head
           ! assuming that the fraction of the head surface area joining the body 
           ! is in proportion to the fraction of the arms and legs, arbitrarily their average
            PJOINs(ii)=(PJOINs(3)+PJOINs(4))/2.0
            HEADL=GEOM_out(5)
          else if (i==4) then ! Trunk
            part_total(ii)=2.0*3.14159*(GEOM_out(6)/2.0)**2.0+
     &                    3.14159*GEOM_out(6)*GEOM_out(5)
            PJOINs(ii)=(2.0*part_end(3)+2.0*part_end(4))/
     &       part_total(ii)+PJOINs(1)
            TRUNKL=GEOM_out(5)
          end if
      end do

      ! Update current surface areas
      do i = 1,num_parts
          call GEOM_ENDO(MASSs(i),DENSITYs(i),DENSITYs(i),FATPCTs(i),
     &     SHAPEs(i),0D1,SUBQFATs(i),SHAPE_Bs(i),SHAPE_Bs(i),
     &     (DHAIRDs(ii)+DHAIRVs(i))/2D1,(INSDENDs(i)+INSDENVs(i))/
     &     2D1,PJOINs(i),0D1,0D1,0D1,GEOM_out)
          SA_currents(i)=GEOM_out(8)*(1.0-PJOINs(i))
      end do
      stepcount=stepcount+1.0
      HEIGHT_CUR=HEADL+TRUNKL+LEGL
      AREA_CUR=SA_currents(1)+SA_currents(2)+SA_currents(3)*2.0
     & +SA_currents(4)*2.0
C     open(unit=10, file="debug_output.txt", status="replace")
C     write(10, *) "Iteration:", stepcount
C     write(10, *) "SHAPE_Bs:", SHAPE_Bs
C     write(10, *) "PJOINs:", PJOINs
C     write(10, *) "SA_targets:", SA_targets
C     write(10, *) "SA_currents:", SA_currents
C     write(10, *) "Sum of Squared Differences:", 
C    * sum((SA_targets - SA_currents)**2.0)
C     write(10, *) "GEOM_out:", GEOM_out
C     close(10)
      end do

      RESULTS = (/SHAPE_Bs,PJOINs,HEIGHT_CUR,AREA_CUR/)

      end subroutine get_shapes
