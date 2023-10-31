MODULE mod_init
! *** SPHERE at 45 degrees ***
! initializes 
!    uinf 
!    corresponding density on grid
!    exact coefficient of density on grid (since constant)

CONTAINS

SUBROUTINE init(n0)
use globalvars, only: uinf !ra,rb,rc,axmax
use mod_geom !for setgrid
  implicit none
  integer :: n0
  !LOCAL
  real*8 :: rotfact

  rotfact=sqrt(2.d0)/2
  uinf=rotfact*[1, 0, -1];

  call setgrid(n0)
  call setdens
!  call setcoeffdensex

return
END SUBROUTINE init


SUBROUTINE setdens
use globalvars, ONLY: g
  implicit none
  integer :: i,j,l
  real*8 :: rotfact

  rotfact=sqrt(2.d0)/2
  do l=1,2
    do i=-g(l)%noff,g(l)%n+g(l)%noff-1
    do j=0,g(l)%m-1
      ! SLP dens for 45^o flow past sphere uinf=U0*(1,0,-1)/sqrt(2), U0=1
      g(l)%f1(i,j)=-1.5d0 * rotfact
      g(l)%f2(i,j)=0
      g(l)%f3(i,j)= 1.5d0 * rotfact
      g(l)%f1der(i,j,1:9)=0
      g(l)%f2der(i,j,1:9)=0
      g(l)%f3der(i,j,1:9)=0
    enddo
    enddo
  enddo
return
END SUBROUTINE setdens


SUBROUTINE setcoeffdensex
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! coeff for SLP dens for 45^o flow past sphere : f=-1.5Uinf 
! same for both grids
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use globalvars, ONLY : cf1,cf2,cf3
  implicit none
  real*8 :: rotfact

  rotfact=sqrt(2.d0)/2

  cf1(0:3,0:3)=0
  cf2(0:3,0:3)=0
  cf3(0:3,0:3)=0

  cf1(0,0)=-1.5d0*rotfact
  cf3(0,0)=1.5d0*rotfact

  return
END SUBROUTINE setcoeffdensex

END MODULE mod_init
