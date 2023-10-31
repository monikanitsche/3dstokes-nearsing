MODULE mod_geom
! contains all routines corresponding strictly to sphere geometry
!   setgrid : initializes grid, calls mod_sphere_velo1/setdens
!   velo
!   setbase : sets base of target point
!   compcoeffx : computes coeffx at base

PUBLIC
!PUBLIC :: correct_sphere

CONTAINS

SUBROUTINE setgrid(n0)
!cccccccccccccccccccc
! sets grid of gauss points, assume delalf=delbet=h
!cccccccccccccccccccc
   use params, ONLY : pi,pi2,nwinmax
   use globalvars, ONLY : g,h
   implicit none
   integer, INTENT(IN) :: n0
   ! LOCAL
   integer :: i,j,l,n,m,noff
   real*8 alfa,beta

   n=n0; m=n/2; noff=nwinmax; h=pi2/n;
   do l=1,2
     g(l)%n=n
     g(l)%m=m
     g(l)%h = h
     g(l)%noff = noff
     do i=-noff,n+noff
       g(l)%alf(i)=-pi+h*i
     enddo
     do j=0,g(l)%m
       g(l)%bet(j)=-pi/2+h*j
     enddo
   enddo

   do i=-noff,n+noff
   do j=0,m
     alfa=-pi+h*i
     beta=-pi/2+h*j
     g(1)%x(i,j)=cos(alfa)*cos(beta)
     g(1)%y(i,j)=sin(alfa)*cos(beta)
     g(1)%z(i,j)=sin(beta)

     g(2)%x(i,j)=sin(beta)
     g(2)%y(i,j)=cos(alfa)*cos(beta)
     g(2)%z(i,j)=sin(alfa)*cos(beta)
   enddo
   enddo

   ! nJ: needed for double layer
   do l=1,2
     do i=-g(l)%noff,g(l)%n+g(l)%noff-1
     do j=0,g(l)%m-1
       g(l)%jac(i,j)=cos(g(l)%bet(j))
       g(l)%xnj(i,j)=g(l)%x(i,j)*g(l)%jac(i,j)
       g(l)%ynj(i,j)=g(l)%y(i,j)*g(l)%jac(i,j)
       g(l)%znj(i,j)=g(l)%z(i,j)*g(l)%jac(i,j)
     enddo
     enddo
   enddo

!      call setdens(g)  ! called in initialization instead
   return
   END SUBROUTINE setgrid

SUBROUTINE velo(n,x,u)
! computes velocity at x,u for 45^o flow past sphere centered at origin
! or any other flow, uinf is variable
  use params, ONLY: nvecmax
  use types, ONLY: basept
  use globalvars, ONLY: xfar,uinf
  use mod_stokes ! for compslp
  implicit none
  integer, INTENT(IN) :: n
  real*8, dimension(0:nvecmax,3), INTENT(INOUT) :: x,u
  !LOCAL
  integer :: j
  real*8, dimension(3) :: slp !,uex 
  real*8  :: xdist,xfarsq
  TYPE (basept) :: t

  xfarsq=xfar**2
  do j=0,n
    xdist=x(j,1)**2+x(j,2)**2+x(j,3)**2
    if (xdist.ge.xfarsq) then  ! dont track particles outside domain
      u(j,:)=0
    else
      call setbase(x(j,:),t)
      call compslp(x(j,:),t,slp)
      u(j,:) = slp + uinf
!     call exactvelo45(x(j,:),uex)
    endif
  enddo
  return
END SUBROUTINE velo


SUBROUTINE exactvelo45(x0,uu) 
!%% Exact solution for stokes flow past sphere with Uinf=(1,0,0, taken from ock&ock
! tested in 3d/matlab/testexactsphere/velosphererk
implicit none
real*8, dimension(3) :: x0,uu
! LOCAL
real*8 :: x,y,z,rho,rho3,phi,theta,urho,uphi,u,v,w,rotfact

rotfact=sqrt(2.d0)/2

!x=x0(1)
!y=x0(2)
!z=x0(3)
x=rotfact*(x0(1)-x0(3))
y=x0(2)
z=rotfact*(x0(1)+x0(3))

!find spherical coords rho,phi,theta (follow stewart convention theta=polar coord)
rho=sqrt(x**2+y**2+z**2)
phi=acos(x/rho)
theta=atan2(z,y)
rho3=rho**3; 

urho=(  1 -3/(2*rho) +1/(2*rho3) )*cos(phi)
uphi=( -1 +3/(4*rho) +1/(4*rho3) )*sin(phi)

u=urho*x/rho - uphi*sin(phi)
v=urho*y/rho + uphi*cos(phi)*cos(theta)
w=urho*z/rho + uphi*cos(phi)*sin(theta)

uu(1)=rotfact*(u+w)
uu(2)=v
uu(3)=rotfact*(w-u)
!uu=[u,v,w]
return
END SUBROUTINE exactvelo45


SUBROUTINE setbase(x0,t)
! finds orthogonal projection of x0 onto sphere (center=origin, rad=1), 
! then determines which grid
! and sets xbase: d,x0,x0b,icorr,roundoff, basepoint alfb and betb,i0,j0
use params, ONLY: pi
use types, ONLY: basept
use globalvars, ONLY: correct,h
implicit none
real*8, dimension(3), INTENT(IN)  :: x0
TYPE (basept), INTENT(OUT) :: t
!LOCAL
real*8 :: xmag,disttobase,alf,bet

t%igrid=2
xmag=sqrt(x0(1)**2+x0(2)**2+x0(3)**2)
t%d=abs(xmag-1); t%x0=x0; t%x0b=x0/xmag   

if (abs(t%x0b(3))<1/sqrt(2.d0)) t%igrid=1

t%icorr=0; t%roundoff=.false.
!if correction needed set 
!       icorr=1
!       alfb and betb (depending on grid) 
!       nearest gridpt (i0,j0) 

if (correct) then
if (t%d<6*h) then 
   t%icorr=1
   if (t%igrid.eq.1) then
      t%betb=asin(t%x0b(3))
      t%alfb=atan2(t%x0b(2),t%x0b(1))
   else  ! t%igrid=2, from grid 1 to grid 2 :   x->y,  y->z,  z->x
      t%betb=asin(t%x0b(1))
      t%alfb=atan2(t%x0b(3),t%x0b(2))
   endif

   t%i0=nint((t%alfb+pi)/h); t%j0=nint((t%betb+pi/2)/h)
   !if roundoff correction needed set roundoff=.true.
   if (t%d<h/4) then
      alf=-pi+h*t%i0; bet=-pi/2+h*t%j0
      disttobase=sqrt((t%alfb-alf)**2+ (t%betb-bet)**2)
      if (disttobase<h/4) t%roundoff=.true.
   endif 
   call compcoeffx(t)
endif
endif
return
END SUBROUTINE setbase


SUBROUTINE compcoeffx(t)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! computes expansions of x,y,z,n1,n2,n3,jac about basepoint alfb,betb
! use known exact formulas for the solid sphere
!
! Notation: cjk - coefficient of alf^j bet^k
! for fourth order, double layer, need 
!         x,y,z to h^4 gives x^3 to h^6 and e to h^5 (skipping..)
!         e     to h^5 (but skip O(d) in highest term)
!         ns    to h^3
!         jac   to h^3
! for fourth order sgle layer need a little less
!         x,y,z to h^3 gives x^2 to h^4, e to h^4 (skipping O(d) at end)
!         e     to h^4 
!         ns    to h^2
!         jac   to h^2
! cmax: max derivative in x,y,z
! gmax: max derivative in density (incl ns,jac)
! emax: max derivative in e
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : cmax
use types, ONLY : basept
use globalvars, ONLY : cx,cy,cz,cnjac1,cnjac2,cnjac3,cj,e,calf,cbet,calfbet
use products
implicit none
TYPE (basept), INTENT(IN) :: t
!LOCAL
integer :: j,k,p
real*8 :: alfb,betb,xb,yb,zb,calfsq,cbetsq,dx,dy,dz
real*8                :: sinalf,cosalf,sinbet,cosbet,fact
real*8,dimension(0:4) :: dalfsin,dalfcos,dbetsin,dbetcos,factorial
real*8,dimension(0:cmax,0:cmax) :: cxx,cyy,czz,cnx,cny,cnz

alfb=t%alfb
betb=t%betb

xb=t%x0b(1)
yb=t%x0b(2)
zb=t%x0b(3)


dx=xb-t%x0(1)  !=dx O(d)
dy=yb-t%x0(2)  !=dy O(d)
dz=zb-t%x0(3)  !=dz O(d)

sinalf=sin(alfb)
cosalf=cos(alfb)
sinbet=sin(betb)
cosbet=cos(betb)
dalfsin(0:4)=(/ sinalf, cosalf,-sinalf,-cosalf,sinalf /)
dalfcos(0:4)=(/ cosalf,-sinalf,-cosalf, sinalf,cosalf /)
dbetsin(0:4)=(/ sinbet, cosbet,-sinbet,-cosbet,sinbet /)
dbetcos(0:4)=(/ cosbet,-sinbet,-cosbet, sinbet,cosbet /)
factorial(0:4)=(/ 1,1,2,6,24 /)

if (t%igrid.eq.1) then
  do j=0,4
  do p=0,j
    k=j-p
    fact=factorial(p)*factorial(k)
    cx(p,k)=dalfcos(p)*dbetcos(k)/fact
    cy(p,k)=dalfsin(p)*dbetcos(k)/fact
    cz(p,k)=0
    cj(p,k)=0
  enddo
  enddo
  do k=0,4
    fact=factorial(k)
    cz(0,k)=dbetsin(k)/fact
    cj(0,k)=dbetcos(k)/fact
  enddo

elseif (t%igrid.eq.2) then
  ! from grid 1 -> grid 2
  ! x -> y
  ! y -> z
  ! z -> x
  do j=0,4
  do p=0,j
    k=j-p
    fact=factorial(p)*factorial(k)
    cy(p,k)=dalfcos(p)*dbetcos(k)/fact
    cz(p,k)=dalfsin(p)*dbetcos(k)/fact
    cx(p,k)=0
    cj(p,k)=0
  enddo
  enddo
  do k=0,4
    fact=factorial(k)
    cx(0,k)=dbetsin(k)/fact
    cj(0,k)=dbetcos(k)/fact
  enddo
endif

cx(0,0)=dx
cy(0,0)=dy
cz(0,0)=dz
! only need j for 0:3, so this is a little overkill

calfsq =     cx(1,0)**2+cy(1,0)**2+cz(1,0)**2 &
        + 2*(dx*cx(2,0)+dy*cy(2,0)+dz*cz(2,0))
cbetsq =     cx(0,1)**2+cy(0,1)**2+cz(0,1)**2 &
        + 2*(dx*cx(0,2)+dy*cy(0,2)+dz*cz(0,2))
calf = sqrt(calfsq)
cbet = sqrt(cbetsq)
calfbet=0

call prod445(cx,cx,cxx)
call prod445(cy,cy,cyy)
call prod445(cz,cz,czz)
e=cxx+cyy+czz
! only need e for 3:4, so this is a little overkill

! Expansion of normal: for sphere, n=x
cnx=cx
cny=cy
cnz=cz
cnx(0,0)=xb
cny(0,0)=yb
cnz(0,0)=zb

call prod333(cnx,cj,cnjac1)
call prod333(cny,cj,cnjac2)
call prod333(cnz,cj,cnjac3)

! Note:  while x is to 0:4, only need j and n to 0:3, 
! so above is a little more ops than needed
return
END SUBROUTINE compcoeffx

END MODULE mod_geom
