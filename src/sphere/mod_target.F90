MODULE mod_target
! contains routines to set target points in different test cases

PUBLIC
CONTAINS

SUBROUTINE target_cross(h,top,right,nvec,xvec)
use params, only : nvecmax
implicit none
integer, intent(out) :: nvec
real*8, intent(in) :: h,top,right
real*8, dimension(0:nvecmax,3), intent(out) :: xvec
!LOCAL
integer :: i,j,k,nx,ny
real*8 :: xx,yy,zz,radsq,eps

!h=0.5d0 top=4; right=5; 

ny=top/h; nx=right/h;
eps=1.d-9

j=-1
do i = -nx,nx
do k = -ny,ny
  xx=i*h
  yy=0
  zz=k*h
  radsq= xx**2 + yy**2 + zz**2
  if (radsq.gt.(1+eps)) then
    j=j+1
    xvec(j,:) = (/ xx , yy , zz /)
  endif
enddo
enddo
nvec=j
print*,'nvec=',nvec
if (nvec.ge.nvecmax) stop 'nvec > nvecmax !'
!print*,xvec(0,:)

write(7,'(i5,a)')nvec,' nvec'
write(7,'(f6.2,2f7.2,a)')h,top,right,' h,top,right'

return
END SUBROUTINE target_cross

END MODULE mod_target
