MODULE mod_geom
! contains all routines corresponding strictly to ellip geometry
!   setgrid : initializes grid, calls mod_sphere_velo1/setdens
! velo
!   setbase : sets base of target point
!   compcoeffx : computes coeffx at base

PUBLIC
!PUBLIC :: correct_sphere

CONTAINS

SUBROUTINE setgrid(n0)
   !cccccccccccccccccccc
   ! sets grid of gauss points, assume delalf=delbet=h
   !cccccccccccccccccccc
   use params, only: pi,pi2,nwinmax
   use globalvars, ONLY : g,h,ra,rb,rc
   implicit none
   integer, INTENT(IN) :: n0
   ! LOCAL
   integer :: i,j,l,n,m,noff
   real*8 :: alfa,beta,abc,cosbet,abccosbet,xx,yy,zz

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
     g(1)%x(i,j)=ra*cos(alfa)*cos(beta)
     g(1)%y(i,j)=rb*sin(alfa)*cos(beta)
     g(1)%z(i,j)=rc*sin(beta)

     g(2)%x(i,j)=ra*sin(beta)
     g(2)%y(i,j)=rb*cos(alfa)*cos(beta)
     g(2)%z(i,j)=rc*sin(alfa)*cos(beta)
   enddo
   enddo
!   print*,g(2)%x(10,15)
!   print*,g(2)%z(10,15)

   ! nJ: needed for double layer
   abc=ra*rb*rc
   do l=1,2
     do j=0,g(l)%m-1
       cosbet=cos(g(l)%bet(j))
       abccosbet=abc*cosbet
       do i=-g(l)%noff,g(l)%n+g(l)%noff-1
         xx=g(l)%x(i,j)
         yy=g(l)%y(i,j)
         zz=g(l)%z(i,j)
         g(l)%xnj(i,j)=abccosbet*xx/ra**2 !nj(1)
         g(l)%ynj(i,j)=abccosbet*yy/rb**2 !nj(2)
         g(l)%znj(i,j)=abccosbet*zz/rc**2 !nj(3)
         g(l)%jac(i,j)=abccosbet*sqrt(xx**2/ra**4+yy**2/rb**4+zz**2/rc**4)
       enddo
     enddo
   enddo

   return
   END SUBROUTINE setgrid


SUBROUTINE velo(n,x,u)
! computes velo at x,u for flow past 3-2-1 ellipse, uinf (1,0,0)
use params, ONLY: nvecmax
use types, ONLY: basept
use globalvars, ONLY: xfar,uinf 
use mod_stokes !for compslp, compdlp
implicit none
integer, INTENT(IN) :: n
real*8, dimension(0:nvecmax,3), INTENT(IN) :: x
real*8, dimension(0:nvecmax,3), INTENT(OUT) :: u
!LOCAL
integer :: j,j0
real*8, dimension(3) :: slp,dlp
real*8 :: xdist,xfarsq
TYPE (basept) :: t

!j0=159
!j0=455
j0=79
xfarsq=xfar**2
do j=0,n
!do j=j0,j0
!print*,x(j,:)
   xdist=x(j,1)**2+x(j,2)**2+x(j,3)**2
   if (xdist.ge.xfarsq) then  ! dont track particles outside domain
      u(j,:)=0
   else
      call setbase(x(j,:),t)
      call compsdlp(x(j,:),t,slp,dlp)
      u(j,:) = slp + dlp + uinf
!print*,slp
!print*,t%alfb,t%betb,dlp
   endif
print*,dlp !,t%icorr
!stop
enddo

return
END SUBROUTINE velo


SUBROUTINE setbase(x0,t)
! finds projection of x0 onto ellipse (base)
! then determines which grid
! and sets tpt: d,x0,x0b,icorr,roundoff, alfb, betb,i0,j0
use params, ONLY :pi
use types, ONLY : basept
use globalvars, ONLY :ra,rb,rc,axmax,correct,h
implicit none
real*8, dimension(3), INTENT(IN)  :: x0
TYPE (basept), INTENT(OUT) :: t
!LOCAL
real*8 :: d,alfb,betb,sinalf,sinbet,cosalf,cosbet,xalf,xbet,dels,alf,bet,disttobase
real*8, dimension(3) :: x0b

t%igrid=2
t%x0=x0
!print*,'in setbase'
!print*,x0
call findproj(x0,ra,rb,rc,x0b,d)
t%x0b=x0b
!print*,x0b
!print*,correct
t%d=d
if (abs(x0b(3)/rc)<1/sqrt(2.d0)) t%igrid=1

!print*,'h',h
!print*,'here',t%igrid,d,6*axmax*h

t%icorr=0
t%roundoff=.false.
!if correction possibly needed set 
!       alfb and betb (depending on grid) to get dels
!       nearest gridpt (i0,j0) 
if (correct) then
!print*,d,6*axmax*h,d-6*axmax*h
if (d<6*axmax*h) then
  if (t%igrid.eq.1) then
    betb= asin(x0b(3)/rc)
    alfb=atan2(x0b(2)/rb,x0b(1)/ra)
    sinalf=sin(alfb); cosalf=cos(alfb);
    sinbet=sin(betb); cosbet=cos(betb);
    xalf= sqrt( cosbet**2*( (ra*sinalf)**2 + (rb*cosalf)**2 ) )
    xbet= sqrt( sinbet**2*( (ra*cosalf)**2 + (rb*sinalf)**2 ) + (rc*cosbet)**2)
    dels=h*max(xalf,xbet)
  else  ! t%igrid=2, from grid 1 to grid 2 :   x->y,  y->z,  z->x
    betb= asin(x0b(1)/ra)
    alfb=atan2(x0b(3)/rc,x0b(2)/rb)
    sinalf=sin(alfb); cosalf=cos(alfb);
    sinbet=sin(betb); cosbet=cos(betb);
    xalf= sqrt( cosbet**2*( (rb*sinalf)**2 + (rc*cosalf)**2 ) )
    xbet= sqrt( sinbet**2*( (rb*cosalf)**2 + (rc*sinalf)**2 ) + (ra*cosbet)**2)
    dels=h*max(xalf,xbet)
  endif
!print*,alfb
!print*,betb

  t%alfb=alfb
  t%betb=betb
!print*,dels

  !if correction needed set icorr=1
!  if (t%d<6*dels) t%icorr=1
  if (t%d<6*dels) then 
    t%icorr=1
    t%i0=nint((t%alfb+pi)/h); t%j0=nint((t%betb+pi/2)/h)
    !print*,'h',t%alfb,t%i0,h,t%icorr
    !if roundoff correction needed set roundoff=.true.
    if (t%d<dels/4) then
      alf=-pi+h*t%i0
      bet=-pi/2+h*t%j0
      disttobase=sqrt((t%alfb-alf)**2+ (t%betb-bet)**2)
      if (disttobase<dels/4) t%roundoff=.true.
    endif
  endif
!print*,x0
!print*
!print*,x0b,t%d
!print*,t%roundoff
!print*,t%d
!print*,t%igrid, t%icorr
!print*,t%alfb, t%betb
!print*,t%i0,t%j0
!print*,dels
  call compcoeffx(t)
endif
endif
!print*,t%icorr
!stop

return
END SUBROUTINE setbase


SUBROUTINE findproj(x0,a,b,c,x0b,d)
! Input: target pt x0, ellip parameters a,b,c
! output: basepoint x0b, distance to base d
implicit none
real*8, dimension(3), INTENT(IN) :: x0
real*8, INTENT(IN) :: a,b,c
real*8, dimension(3), INTENT(OUT) :: x0b
real*8, INTENT(OUT) :: d
!LOCAL
integer :: k
real*8 :: a2,b2,c2,x02,y02,z02,f,fp,lamold,lam,err
a2=a**2; b2=b**2; c2=c**2;
x02=x0(1)**2; y02=x0(2)**2; z02=x0(3)**2
lam=0
do k=1,10
  f=a2*x02/(a2-lam)**2+b2*y02/(b2-lam)**2+c2*z02/(c2-lam)**2-1
  fp=2*(a2*x02/(a2-lam)**3+b2*y02/(b2-lam)**3+c2*z02/(c2-lam)**3)
  lamold=lam; lam=lam-f/fp; err=abs(lam-lamold)
!  write(*,'(a,i3,f10.6,e10.2,2f10.6)')'newton',k,lam,err,f,fp
  if (err.lt.1e-11) goto 900
enddo
900 continue
!print*,'newton steps',k,err
if (err.gt.1.e-10) then
  print*,'warning, error in Newtons large',err; !stop
endif

x0b(1)=x0(1)/(1-lam/a2)
x0b(2)=x0(2)/(1-lam/b2)
x0b(3)=x0(3)/(1-lam/c2)
d= sqrt( (x0(1)-x0b(1))**2 + (x0(2)-x0b(2))**2 + (x0(3)-x0b(3))**2 )

return
END SUBROUTINE findproj


SUBROUTINE compcoeffx(t)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! computes expansions of x,y,z,n1,n2,n3,jac about basepoint alfb,betb
! use known exact formulas for the ellipse
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
use globalvars, ONLY : ra,rb,rc,cx,cy,cz,cnjac1,cnjac2,cnjac3,cj,e,calf,cbet,calfbet
use products
implicit none
TYPE (basept), INTENT(IN) :: t
!LOCAL
integer :: j,k,p
real*8 :: alfb,betb,xb,yb,zb,dx,dy,dz,calfsq,cbetsq
real*8 :: sinalf,cosalf,sinbet,cosbet,fact,sin2alf,cos2alf,sin2bet,cos2bet,cosbetsq,cosalfsq,sinalfsq,sinbetsq
real*8,dimension(0:4) :: dalfsin,dalfcos,dbetsin,dbetcos,factorial
real*8,dimension(0:5) :: dcosbetsq,dcosbetsinbet
real*8,dimension(0:cmax,0:cmax) :: cxx,cyy,czz
real*8 :: abc,a2,b2,c2,a4,b4,c4,root,root3,ac,bc,ab, &
         arg,arg10,arg01,arg20,arg11,arg02

alfb=t%alfb
betb=t%betb
!print*,alfb,betb

xb=t%x0b(1)
yb=t%x0b(2)
zb=t%x0b(3)

dx=xb-t%x0(1)  !=dx O(d)
dy=yb-t%x0(2)  !=dy O(d)
dz=zb-t%x0(3)  !=dz O(d)

sinalf=sin(alfb); cosalf=cos(alfb)
sinbet=sin(betb); cosbet=cos(betb)

dalfsin(0:4)=(/ sinalf, cosalf,-sinalf,-cosalf,sinalf /)
dalfcos(0:4)=(/ cosalf,-sinalf,-cosalf, sinalf,cosalf /)
dbetsin(0:4)=(/ sinbet, cosbet,-sinbet,-cosbet,sinbet /)
dbetcos(0:4)=(/ cosbet,-sinbet,-cosbet, sinbet,cosbet /)
factorial(0:4)=(/ 1,1,2,6,24 /)

!derivatives of cosbet^2,cosbet*sinbet : see expjacgrid2.nb
cosbetsq=cosbet**2
sinbetsq=sinbet**2
dcosbetsq(0:5)=(/ cosbetsq , -2*cosbet*sinbet , &
      -2*cosbetsq+2*sinbetsq , 8*cosbet*sinbet , &
      8*cosbetsq-8*sinbetsq , -32*cosbet*sinbet /)
dcosbetsinbet(0:5)=(/ cosbet*sinbet, cosbetsq-sinbetsq, &
      -4*cosbet*sinbet, -4*cosbetsq+4*sinbetsq,  &
      16*cosbet*sinbet, 16*cosbetsq-16*sinbetsq /)

abc=ra*rb*rc
a4=ra**4; b4=rb**4; c4=rc**4
a2=ra**2; b2=rb**2; c2=rc**2
sinalfsq=sinalf**2
cosalfsq=cosalf**2
sin2alf=2*sinalf*cosalf
cos2alf=cos(2*alfb)
sin2bet=2*sinbet*cosbet
cos2bet=cos(2*betb)

if (t%igrid.eq.1) then
  do j=0,4
  do p=0,j
    k=j-p
    fact=factorial(p)*factorial(k)
    cx(p,k)=ra*dalfcos(p)*dbetcos(k)/fact
    cy(p,k)=rb*dalfsin(p)*dbetcos(k)/fact
    cz(p,k)=0
  enddo
  enddo
  do k=0,4
    fact=factorial(k)
    cz(0,k)=rc*dbetsin(k)/fact
  enddo
 
  !For expansion of nJ, see expjacgrid1v
  bc=rb*rc
  cnjac1(0,0)=   bc*cosalf*cosbetsq
  cnjac1(1,0)=  -bc*sinalf*cosbetsq
  cnjac1(0,1)=  -bc*cosalf*sin2bet
  cnjac1(2,0)=  -bc*cosalf*cosbetsq
  cnjac1(1,1)=   bc*sinalf*sin2bet
  cnjac1(0,2)=-2*bc*cosalf*cos2bet
  cnjac1(3,0)=   bc*sinalf*cosbetsq
  cnjac1(2,1)=   bc*cosalf*sin2bet
  cnjac1(1,2)= 2*bc*sinalf*cos2bet
  cnjac1(0,3)= 4*bc*cosalf*sin2bet

  ac=ra*rc
  cnjac2(0,0)=    ac*sinalf*cosbetsq
  cnjac2(1,0)=    ac*cosalf*cosbetsq 
  cnjac2(0,1)=   -ac*sinalf*sin2bet
  cnjac2(2,0)=   -ac*sinalf*cosbetsq 
  cnjac2(1,1)=   -ac*cosalf*sin2bet
  cnjac2(0,2)= -2*ac*sinalf*cos2bet
  cnjac2(3,0)=   -ac*cosalf*cosbetsq
  cnjac2(2,1)=    ac*sinalf*sin2bet
  cnjac2(1,2)= -2*ac*cosalf*cos2bet
  cnjac2(0,3)=  4*ac*sinalf*sin2bet

  ab=ra*rb
  cnjac3(0,0)=    ab*sin2bet/2
  cnjac3(1,0)= 0
  cnjac3(0,1)=    ab*cos2bet
  cnjac3(2,0)= 0  
  cnjac3(1,1)= 0
  cnjac3(0,2)= -2*ab*sin2bet
  cnjac3(3,0)= 0  
  cnjac3(2,1)= 0
  cnjac3(1,2)= 0
  cnjac3(0,3)= -4*ab*cos2bet

  !Expansion of J, see expjacgrid1v
  arg=( xb**2/a4 + yb**2/b4 + zb**2/c4 ); root=sqrt(arg); root3=root**3
  arg10= (a2-b2)*cosbetsq*sin2alf/(a2*b2)
  arg01= 2*cosbet*(1/c2-cosalfsq/a2-sinalfsq/b2)*sinbet
  arg20= 2*(a2-b2)*cos2alf*cosbetsq/(a2*b2)
  arg11= -(a2-b2)*sin2alf*sin2bet/(a2*b2)
  arg02= 2*cos2bet*(1/c2-cosalfsq/a2-sinalfsq/b2)

else if (t%igrid.eq.2) then
  do j=0,4
  do p=0,j
    k=j-p
    fact=factorial(p)*factorial(k)
    cy(p,k)=rb*dalfcos(p)*dbetcos(k)/fact
    cz(p,k)=rc*dalfsin(p)*dbetcos(k)/fact
    cx(p,k)=0
  enddo
  enddo
  do k=0,4
    fact=factorial(k)
    cx(0,k)=ra*dbetsin(k)/fact
  enddo

  bc=rb*rc
  cnjac1(0,0)=    bc*sin2bet/2
  cnjac1(1,0)=  0
  cnjac1(0,1)=    bc*cos2bet
  cnjac1(2,0)=  0
  cnjac1(1,1)=  0
  cnjac1(0,2)= -2*bc*sin2bet 
  cnjac1(3,0)=  0  
  cnjac1(2,1)=  0
  cnjac1(1,2)=  0
  cnjac1(0,3)= -4*bc*cos2bet

  ac=ra*rc
  cnjac2(0,0)=    ac*cosalf*cosbetsq 
  cnjac2(1,0)=   -ac*sinalf*cosbetsq
  cnjac2(0,1)=   -ac*cosalf*sin2bet
  cnjac2(2,0)=   -ac*cosalf*cosbetsq  
  cnjac2(1,1)=    ac*sinalf*sin2bet
  cnjac2(0,2)= -2*ac*cosalf*cos2bet 
  cnjac2(3,0)=    ac*sinalf*cosbetsq
  cnjac2(2,1)=    ac*cosalf*sin2bet
  cnjac2(1,2)=  2*ac*sinalf*cos2bet
  cnjac2(0,3)=  4*ac*cosalf*sin2bet

  ab=ra*rb
  cnjac3(0,0)=    ab*sinalf*cosbetsq
  cnjac3(1,0)=    ab*cosalf*cosbetsq
  cnjac3(0,1)=   -ab*sinalf*sin2bet
  cnjac3(2,0)=   -ab*sinalf*cosbetsq
  cnjac3(1,1)=   -ab*cosalf*sin2bet
  cnjac3(0,2)= -2*ab*sinalf*cos2bet
  cnjac3(3,0)=   -ab*cosalf*cosbetsq
  cnjac3(2,1)=    ab*sinalf*sin2bet
  cnjac3(1,2)= -2*ab*cosalf*cos2bet
  cnjac3(0,3)=  4*ab*sinalf*sin2bet

  arg=( xb**2/a4 + yb**2/b4 + zb**2/c4 ); root=sqrt(arg); root3=root**3
  arg10= (b2-c2)*cosbetsq*sin2alf/(b2*c2)
  arg01= 2*cosbet*(1/a2-cosalfsq/b2-sinalfsq/c2)*sinbet
  arg20= 2*(b2-c2)*cos2alf*cosbetsq/(b2*c2)
  arg11= -(b2-c2)*sin2alf*sin2bet/(b2*c2)
  arg02= 2*cos2bet*(1/a2-cosalfsq/b2-sinalfsq/c2)
endif

cx(0,0)=dx
cy(0,0)=dy
cz(0,0)=dz

!print*,dx*cx(1,1)+dy*cy(1,1)+dz*cz(1,1)
!print*,cx(1,0)*cx(0,1)+cy(1,0)*cy(0,1)+cz(1,0)*cz(0,1)

!print*,dx,dy,dz
calfsq =     cx(1,0)**2+cy(1,0)**2+cz(1,0)**2 &
        + 2*(dx*cx(2,0)+dy*cy(2,0)+dz*cz(2,0))
cbetsq =     cx(0,1)**2+cy(0,1)**2+cz(0,1)**2 &
        + 2*(dx*cx(0,2)+dy*cy(0,2)+dz*cz(0,2))
calf = sqrt(calfsq)
cbet = sqrt(cbetsq)
calfbet= cx(1,0)*cx(0,1)+cy(1,0)*cy(0,1)+cz(1,0)*cz(0,1) &
        + (dx*cx(1,1)+dy*cy(1,1)+dz*cz(1,1))
!calfbet= 0
!print*,'calf,cbet,calfbet'
!print*,calf/cbet,calfbet/(calf*cbet)
!print*,cx(1,0)
!print*,alfb,betb,calf,cbet,calfbet

call prod445(cx,cx,cxx)
call prod445(cy,cy,cyy)
call prod445(cz,cz,czz)
e=cxx+cyy+czz

!Expansion of J, see expjacgrid1v
!do k=0,3
!  cnjac1(k,0)=abc*cosbet*cx(k,0)/a2
!  cnjac2(k,0)=abc*cosbet*cy(k,0)/b2
!  cnjac3(k,0)=abc*cosbet*cz(k,0)/c2
!enddo

! Expansion of jac and njac, see expjac.nb
! coeff of njac
do j=0,3
do p=0,j
  k=j-p
  fact=factorial(p)*factorial(k)
  cnjac1(p,k)=cnjac1(p,k)/fact
  cnjac2(p,k)=cnjac2(p,k)/fact
  cnjac3(p,k)=cnjac3(p,k)/fact
enddo
enddo

! coeff of jac
cj(0,0)= cosbet*root
cj(1,0)= arg10*cosbet/(2*root)
cj(0,1)= (arg01*cosbet-2*arg*sinbet)/(2*root)
cj(2,0)= -(arg10**2-2*arg*arg20)*cosbet/(4*root3)
cj(1,1)= -( (arg01*arg10-2*arg*arg11)*cosbet + 2*arg*arg10*sinbet )/(4*root3)
cj(0,2)= ( -(4*arg**2+arg01**2-2*arg*arg02)*cosbet - 4*arg*arg01*sinbet )/(4*root3)
do j=0,2
do p=0,j
  k=j-p
  fact=factorial(p)*factorial(k)
  cj(p,k)=abc*cj(p,k)/fact
enddo
enddo

return
END SUBROUTINE compcoeffx

END MODULE mod_geom
