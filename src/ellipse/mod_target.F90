MODULE mod_target
! contains routines to set target points in different test cases

PUBLIC
CONTAINS

SUBROUTINE target_cross(h,top,right,nvec,xvec)
use params, only : nvecmax
use globalvars, only : ra,rb,rc,uinf
implicit none
integer, intent(out) :: nvec
real*8, intent(in) :: h,top,right
real*8, dimension(0:nvecmax,3), intent(out) :: xvec
!LOCAL
integer :: i,j,k,nx,ny
real*8 :: xx,yy,zz,rsq,eps

!h=0.5d0 top=4; right=5; 

ny=top/h; nx=right/h;
eps=1.d-9

j=-1
do i = -nx,nx
do k = -ny,ny
  xx=i*h
  yy=0
  zz=k*h
  rsq= (xx/ra)**2 + (yy/rb)**2 + (zz/rc)**2
  if (rsq.gt.(1+eps)) then
    j=j+1
    xvec(j,:) = (/ xx , yy , zz /)
  endif
enddo
enddo
nvec=j
print*,'nvec=',nvec
if (nvec.ge.nvecmax) stop 'nvec > nvecmax !'
!print*,xvec(0,:)

write(7,'(3f10.2,a)')ra,rb,rc,' ra,rb,rc'
write(7,'(3f10.2,a)')uinf,' uinf'
write(7,'(i5,a)')nvec,' nvec'
write(7,'(3f10.5,a)')h,top,right,' h,top,right'

return
END SUBROUTINE target_cross

SUBROUTINE target_ray(nvec,xvec)
use params, only : nvecmax
use globalvars, only : ra,rb,rc,uinf
implicit none
integer, intent(out) :: nvec
real*8, dimension(0:nvecmax,3), intent(out) :: xvec
!LOCAL
integer :: l
real*8 :: xx,yy,zz,d,dd(0:20),xn(3),mag,xb(3)

dd(0:5)= (/ 1.d0,0.8d0,0.7d0,0.5d0,0.4d0,0.3d0 /)
do l=0,11
dd(l+6)=0.2d0/2**l
enddo
do l=0,17
  d=dd(l)
!X-DIRECTION
!  xvec(l,:) = (/ ra+d , 0.d0 , 0.d0 /)

!Y-DIRECTION
  xvec(l,:) = (/ 0.d0, rb+d , 0.d0 /)

!OFFDIRECTION
!  xb=(/ -2.1382770687061301d0,  1.4015774873524995d0, 2.9490880088390157D-002 /)
!xn(1)=xb(1)/ra**2; xn(2)=xb(2)/rb**2; xn(3)=xb(3)/rc**2;
!mag=sqrt(xn(1)**2+xn(2)**2+xn(3)**2);
!  xvec(l,:) = xb+d/mag*xn

!  xvec(l,:) = (/ xx , yy , zz /)

enddo
nvec=17

return
END SUBROUTINE target_ray


SUBROUTINE target_patch(d)
use params, only : nvecmax,pi
use globalvars, only : ra,rb,rc
use mod_geom ! for velo
implicit none
real*8, intent(in) :: d
!LOCAL
real*8, dimension(0:nvecmax,3) :: xvec,u
integer :: m,i,j
real*8, dimension(3) :: xb,xn
real*8 :: alf,bet,beta(0:500),h,mag

!use mxm points on (alf,bet) in [0,pi/2]x[0,pi/2]
!finest grid for errors has h=pi2/320
!m=320 !gives grid of pi/2/320 so 4x4 on finest grid 
m=160 !gives grid of pi/2/160 so 2x2 on finest grid 
h=(pi/2)/m
!m=64 !gives grid of pi/2/160 so 2x2 on finest grid 
!h=(pi/5)/m

do i=0,m
!do i=94,94
!do i=160,160
alf=i*h;
!print*,alf
!alf=0.6+i*h;
  do j=0,m
    bet=j*h;
    beta(j)=bet
    xb=(/ ra*cos(alf)*cos(bet),rb*sin(alf)*cos(bet),rc*sin(bet) /)
    xn(1)=xb(1)/ra**2; xn(2)=xb(2)/rb**2; xn(3)=xb(3)/rc**2;
    mag=sqrt(xn(1)**2+xn(2)**2+xn(3)**2);
    xvec(j,:) = xb+d*xn/mag
!print*,alf,bet
  enddo
!print*,beta(79)
  call velo(m,xvec,u)
  do j=0,m
    write(8,*)alf,beta(j),u(j,:)
  enddo
enddo
! alf=0 bet\in [0,pi/2] x=(3,0,0) to (0,0,1) (y=0)
! alf=pi/2 bet\in [0,pi/2] x=(0,2,0) to (0,0,1) (x=0)

return
END SUBROUTINE target_patch

END MODULE mod_target
