MODULE mod_runge
CONTAINS
SUBROUTINE runge(n,xvec)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Solves system using 4th order Runge Kutta
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : nvecmax
use globalvars, ONLY: time,delt 
use mod_geom ! for velo
implicit none
integer :: n,j
real*8, dimension(0:nvecmax,3) :: u,k1,k2,k3,k4,xhold,xvec
real*8 :: umax

call velo(n,xvec,u)
k1(0:n,:) = delt*u(0:n,:)
xhold(0:n,:) = xvec(0:n,:) + k1(0:n,:)/2

call velo(n,xhold,u)
k2(0:n,:) = delt*u(0:n,:)
xhold(0:n,:) = xvec(0:n,:) + k2(0:n,:)/2

call velo(n,xhold,u)
k3(0:n,:) = delt*u(0:n,:)
xhold(0:n,:) = xvec(0:n,:) + k3(0:n,:)

call velo(n,xhold,u)
k4(0:n,:) = delt*u(0:n,:)

xvec(0:n,:) = xvec(0:n,:) + (k1(0:n,:)+2*k2(0:n,:)+2*k3(0:n,:)+k4(0:n,:))/6
time = time+delt

return
END SUBROUTINE runge
END MODULE mod_runge
