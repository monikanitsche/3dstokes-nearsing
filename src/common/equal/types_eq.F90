MODULE types
use params
implicit none

  TYPE :: grid
    integer :: n,m,noff
    real*8 :: h
    real*8,dimension(nleft:nright,0:mgrid) :: &
           x,y,z,f1,f2,f3, xnj,ynj,znj,jac
    real*8,dimension(nleft:nright,0:mgrid,1:9) :: f1der,f2der,f3der
    real*8,dimension(nleft:nright) :: alf
    real*8,dimension(0:mgrid) :: bet
  END TYPE grid

  TYPE basept
    integer :: i0,j0,icorr,igrid
    real*8,dimension(3) :: x0,x0b
    real*8 :: d,alfb,betb
    logical :: roundoff
  END TYPE basept

END MODULE types
