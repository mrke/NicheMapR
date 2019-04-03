      SUBROUTINE dget_aELH(N,a,aELH,daELH,RPAR,IPAR)
      IMPLICIT NONE
      integer N,IPAR
      double precision aELH,daELH,RPAR
      double precision a,V,e,r,p_C,dE,dL,E_R,f,k_M,k_E,p_J,p_Am,E_m,g
     &,kap,L,dE_R,dER,e_s,k_J,vdot,H,dE_H
        DIMENSION aELH(N),daELH(N),RPAR(*),IPAR(*)
      f=RPAR(1)
      k_M=RPAR(2)
      vdot=RPAR(3)
      p_J=RPAR(4)
      p_Am=RPAR(5)
      E_m=RPAR(6)
      g=RPAR(7)
      kap=RPAR(8)
      a  = aELH(1)! % d, time since birth
      E  = aELH(2)! % J, reserve
      L  = aELH(3)! % cm, structural length
      H  = aELH(4)! % J, reproduction buffer

      V = L**3D+00                          ! cm^3, structural volume
      e_s = E/ V/ E_m                    ! -, scaled reserve density
      r = (e_s * vdot - g * k_M)/ (e_s + g) ! 1/d, specific growth rate
      p_C = E * (vdot - r)              ! J/d, mobilisation rate

      dE = f * p_Am * V - p_C          ! J/d, change in reserve
      dL = r * L/ 3D+00                   ! cm/d, change in length
      dE_R = (1 - kap) * p_C - k_J * H     ! J/d, change in reprod buffer

      daELH(1)=1.0D+00
      daELH(2)=dE
      daELH(3)=dL
      daELH(4)=dE_H
        RETURN
        END       