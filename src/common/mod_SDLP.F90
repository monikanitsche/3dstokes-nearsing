MODULE mod_SDLP

CONTAINS
SUBROUTINE SLPatP(tpt,g,slp) 
! computes Stokes single layer potential (3 components) at targpt x
! where using grid and forces given in grid
use types
use params, ONLY: nleft
implicit none
TYPE (basept), INTENT(IN) :: tpt   
TYPE (grid), INTENT(IN) :: g
real*8, dimension(3),INTENT(OUT) :: slp
!LOCAL
real*8, dimension(3) :: fp2,fp1,fm1,fm2,dtop,dbot,dtop3,dbot3,x
integer :: i,j,n,m,imod,iwt
real*8 :: h,h2,h3,h5
real*8, dimension(3) :: f,xh
real*8 jac,dot,rho,rhosq

!print*,'here'
!print*,gg(1)%f1(nleft, 0)
!print*,gg(1)%f1(320, 100)
!print*,g%f1(320, 100)

x=tpt%x0
h=g%h
n=g%n
m=g%m
h2=h**2
h3=h2*h
h5=h2*h3

slp=0  !initializes all 3 elements of slp to 0
do i=0,n-1
do j=1,m-1

iwt=1
if (tpt%roundoff) then
if ((i.eq.mod(tpt%i0,n)).and.(j.eq.tpt%j0)) then
  iwt=0
endif
endif

!  slp=slp+slpcomp(i,j,x,g)
f(1)=g%f1(i,j)
f(2)=g%f2(i,j)
f(3)=g%f3(i,j)
xh(1)=g%x(i,j)-x(1)
xh(2)=g%y(i,j)-x(2)
xh(3)=g%z(i,j)-x(3)
jac=g%jac(i,j)
dot= f(1)*xh(1) + f(2)*xh(2) + f(3)*xh(3)
rhosq= xh(1)**2+xh(2)**2+xh(3)**2
rho= sqrt(rhosq)
slp=slp + iwt*( f + dot*xh/rhosq ) *jac/rho
enddo
enddo
slp=slp*h2

dtop=0
dbot=0
dtop3=0
dbot3=0
do i=0,n-1
  imod=mod(i+n/2,n)
  fm1= slpcomp( i  ,m-1,x,g)
  fp1=-slpcomp(imod,m-1,x,g)
  fm2= slpcomp( i  ,m-2,x,g)
  fp2=-slpcomp(imod,m-2,x,g)
  dtop = dtop  + (-fp2+8*fp1-8*fm1+fm2)/(12*h)
!  dtop3= dtop3 + (fp2-2*fp1+2*fm1-fm2)/(2*h3)
  fp1= slpcomp( i  , 1 ,x,g)
  fm1=-slpcomp(imod, 1 ,x,g)
  fp2= slpcomp( i  , 2 ,x,g)
  fm2=-slpcomp(imod, 2 ,x,g)
  dbot = dbot + (-fp2+8*fp1-8*fm1+fm2)/(12*h)
!  dbot3= dbot3 + (fp2-2*fp1+2*fm1-fm2)/(2*h3)
enddo
slp= slp - h3/12*(dtop-dbot) !+ h5/720*(dtop3-dbot3)

slp=slp/pi8     !divide here and in corrections
return
END SUBROUTINE SLPatP


FUNCTION slpcomp(i,j,x,g)
use types, ONLY : grid
implicit none
integer :: i,j
!real*8, dimension(3), INTENT(IN) :: x
!real*8, dimension(3), INTENT(OUT) :: slpcomp
!TYPE (grid), INTENT(IN) :: g
real*8, dimension(3) :: x
real*8, dimension(3) :: slpcomp
TYPE (grid) :: g
!LOCAL
real*8, dimension(3) :: f,xh
real*8 jac,dot,rho,rhosq

f(1)=g%f1(i,j)
f(2)=g%f2(i,j)
f(3)=g%f3(i,j)
xh(1)=g%x(i,j)-x(1)
xh(2)=g%y(i,j)-x(2)
xh(3)=g%z(i,j)-x(3)
jac=g%jac(i,j)
dot= f(1)*xh(1) + f(2)*xh(2) + f(3)*xh(3)
rhosq= xh(1)**2+xh(2)**2+xh(3)**2
rho= sqrt(rhosq)
slpcomp = ( f + dot*xh/rhosq ) *jac/rho
!      = ( f3/rho + dot*zh/rho3 ) *jac
return
END FUNCTION slpcomp


SUBROUTINE DLPatP(tpt,g,dlp) 
! computes Stokes single layer potential (3 components) at targpt x
! where using grid and forces given in grid
use types
implicit none
TYPE (basept), INTENT(IN) :: tpt   
TYPE (grid), INTENT(IN) :: g   
real*8, dimension(3),INTENT(OUT) :: dlp
!LOCAL
real*8, dimension(3) :: fp2,fp1,fm1,fm2,dtop,dbot,x
integer :: i,j,n,m,imod,iwt
real*8 :: h,h2,h3,jac,dotf,dotn,rhosq,rho5,fact
real*8, dimension(3) :: f,xh,xn

x=tpt%x0
h=g%h
n=g%n
m=g%m
h2=h**2
h3=h2*h

dlp=0  
do i=0,n-1
do j=1,m-1
   iwt=1
   if (tpt%roundoff) then
      if ((i.eq.mod(tpt%i0,n)).and.(j.eq.tpt%j0)) then
         iwt=0
      endif
   endif
   !  dlp=dlp+dlpcomp(i,j,x,g)
   f(1)=g%f1(i,j)
   f(2)=g%f2(i,j)
   f(3)=g%f3(i,j)
   xh(1)=g%x(i,j)-x(1)
   xh(2)=g%y(i,j)-x(2)
   xh(3)=g%z(i,j)-x(3)
   xn(1)=g%xnj(i,j)
   xn(2)=g%ynj(i,j)
   xn(3)=g%znj(i,j)
!   jac=g%jac(j)
   dotf =  f(1)*xh(1) +  f(2)*xh(2) +  f(3)*xh(3)
   dotn = xn(1)*xh(1) + xn(2)*xh(2) + xn(3)*xh(3)
   rhosq= xh(1)**2+xh(2)**2+xh(3)**2
   rho5 = sqrt(rhosq)**5
!   fact = dotf*dotn*jac/rho5
   fact = dotf*dotn/rho5
!if ((i.ne.49).or.(j.ne.22) ) then
   dlp=dlp+iwt*fact*xh
!else
!print*,fact*xh
!endif
enddo
enddo
dlp=dlp*h2

dtop=0
dbot=0
do i=0,n-1
  imod=mod(i+n/2,n)
  fm1= dlpcomp( i  ,m-1,x,g)
  fp1=-dlpcomp(imod,m-1,x,g)
  fm2= dlpcomp( i  ,m-2,x,g)
  fp2=-dlpcomp(imod,m-2,x,g)
  dtop = dtop  + (-fp2+8*fp1-8*fm1+fm2)/(12*h)
  fp1= dlpcomp( i  , 1 ,x,g)
  fm1=-dlpcomp(imod, 1 ,x,g)
  fp2= dlpcomp( i  , 2 ,x,g)
  fm2=-dlpcomp(imod, 2 ,x,g)
  dbot = dbot + (-fp2+8*fp1-8*fm1+fm2)/(12*h)
enddo
dlp= dlp - h3/12*(dtop-dbot) !+ h5/720*(dtop3-dbot3)

dlp=-3*dlp/pi4     !divide here and in corrections
return
END SUBROUTINE DLPatP


FUNCTION dlpcomp(i,j,x,g)
use types, ONLY : grid
implicit none
integer :: i,j
real*8, dimension(3) :: x,dlpcomp
TYPE (grid) :: g
!LOCAL
real*8, dimension(3) :: f,xh,xn
real*8 jac,dotf,dotn,rhosq,rho5,fact

f(1)=g%f1(i,j)
f(2)=g%f2(i,j)
f(3)=g%f3(i,j)
xh(1)=g%x(i,j)-x(1)
xh(2)=g%y(i,j)-x(2)
xh(3)=g%z(i,j)-x(3)
xn(1)=g%xnj(i,j)
xn(2)=g%ynj(i,j)
xn(3)=g%znj(i,j)
!jac=g%jac(j)
dotf =  f(1)*xh(1) +  f(2)*xh(2) +  f(3)*xh(3)
dotn = xn(1)*xh(1) + xn(2)*xh(2) + xn(3)*xh(3)
rhosq= xh(1)**2+xh(2)**2+xh(3)**2
rho5 = sqrt(rhosq)**5
!fact = dotf*dotn*jac/rho5
fact = dotf*dotn/rho5
dlpcomp = fact*xh
return
END FUNCTION dlpcomp

END MODULE mod_SDLP
