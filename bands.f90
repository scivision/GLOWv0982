! Subroutine BANDS

! This software is part of the GLOW model.  Use is governed by the Open Source
! Academic Research License Agreement contained in the file Glowlicense.txt.
! For more information see the file Glow.txt.

! Stan Solomon, 3/2015
! Stan Solomon, 2/2018: changed LBH total band system excitation from aglw(4,3,k)
! to zeta(12,k) so that the fraction of (a,a',w) excitation resulting in (a) state
! emission (i.e., LBH) is accounted for.  This fraction, represented by the branching
! ratio in B(48), is currently esitmated to be 0.7. 

! Version 0.2 is really just a stub that divides total N2 singlet state
! band system excitation rate (a, a', w) into 7 upper states of the N2(a)
! LBH band system (v' = 0-6).  States above v'=6 are presumed to dissociate.
! Frank Condon factors from Ajello and Shemansky, JGR, 90, 9845-9861, 1985

! Input:
!   Use-associated variables from cglow:
!     jmax   number of altitude levels
!     nc     number of emission components (in this case, v' levels)
!     aglw   excited state array (state, species, altitude)
!     zeta   airglow volume emission rate array (emission, species, altitude)
! Output:
!   Use-associated variables from cglow:
!     zlbh   total excitation rate to each LBH v' level (cm-3 s-1)


subroutine bands

  use cglow, only: jmax, nc, aglw, zeta, zlbh

  implicit none

  integer :: k,m
  real :: fcfac(nc)
  data fcfac /0.043,0.114,0.168,0.183,0.160,0.122,0.084,0.0,0.0,0.0/

  do k=1,jmax
    do m=1,nc
      zlbh(m,k)=zeta(12,k)*fcfac(m)
    end do
  end do

  return

end subroutine bands
