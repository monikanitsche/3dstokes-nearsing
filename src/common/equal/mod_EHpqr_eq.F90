#include "flags.h"
MODULE mod_EHpqr
! Computes the corrections E[H_pqr] for all pqr
!
! contains
!    initsubgrid      (sets local subgrid)
!    compallcorrtrap (returns E[Hjks] given basept, grid
!    trapcomponents  (components of basis functions at gridpoints)
!    compallindfact   (computes scaling factors, and indeces of all Hjk)
!    allder_c  (sets 1st,3rd and mixed derivatives of all Hjk)
!    allintFbackrecur2  (computes all integrals of Hjk recursively, c.ne.0)
!    compI005(u1,u2,v1,v2,cc,Iint)
!    intf005dx(a,b,y,c,intfdx)
!    intf005dy(a,b,x,c,intfdy)
!    setgauss(a,b,Ng,xg,wg)


!PRIVATE    !set the default for module
PUBLIC

CONTAINS

SUBROUTINE initsubgrid(t,g,n1,n2,m1,m2,wtalf,wtbet)
! sets endpoints of subgrid and trapezoid weights
use types
implicit none
TYPE (basept),INTENT(IN) :: t
TYPE (grid),INTENT(IN) :: g
integer, INTENT(OUT) :: n1,n2,m1,m2
real*8, dimension(0:*),INTENT(OUT) :: wtalf,wtbet
!LOCAL
integer :: i,j,nwin

!a=-.25d0 !b= .25d0 !c=-.25d0 !d= .25d0
!a=-.5d0 !!      b= .5d0 !!      c=-.5d0 !!      d= .5d0
!i0=nint((t%alfb+pi)/h) !j0=nint((t%betb+pi/2)/h)
!!n1=nint(t%i0+a*pi/h) !!n2=nint(t%i0+b*pi/h)
!!m1=max(nint(t%j0+c*pi/h),0) !!m2=min(nint(t%j0+d*pi/h),g%m)

if (g%n.le.80) then
  nwin=5
elseif (g%n.le.160) then
  nwin=9
elseif (g%n.le.320) then
  nwin=15
else  !if (g%n.le.640) then
  nwin=26
endif
!nwin=g%n/5
!nwin=2*nwin

n1=t%i0-nwin
n2=t%i0+nwin
m1=max(t%j0-nwin,0)
m2=min(t%j0+nwin,g%m)
!print*,t%j0,m1,m2

!if (g%n.eq.320) then
!n1=208-nwin
!n2=208+nwin
!endif
!if (g%n.eq.160) then
!m1=40-nwin
!m2=40+nwin
!endif
!print*,t%i0,t%j0,nwin!,m1,m2,g%bet(m1)
!print*,n1,n2,m1,m2,g%noff
!print*,'in initsub'
!print*,t%i0,t%i0-nwin,g%alf(n1), g%n
!print*,t%i0,t%j0,g%alf(t%i0),g%bet(t%j0),t%alfb,t%betb
!print*,g%alf(n1),g%alf(n2),g%bet(m1),g%bet(m2)
!print*,g%bet(m2),g%bet(m1),t%betb
!print*,(g%bet(m2)+g%bet(m1))/2,t%betb
!print*,(g%bet(m2)+g%bet(m1))/2-t%betb,g%h
!print*,g%bet(t%j0),g%bet(t%j0)-t%betb,g%h

if ((n1.lt.-g%noff).or.(n2.gt.g%n+g%noff)) STOP 'subgrid out of bounds'

do i=n1+1,n2-1
wtalf(i-n1)=1
enddo
do j=m1+1,m2-1
wtbet(j-m1)=1
enddo
wtalf(0)=1/2.d0
wtbet(0)=1/2.d0
wtalf(n2-n1)=1/2.d0
wtbet(m2-m1)=1/2.d0

return
END SUBROUTINE initsubgrid


SUBROUTINE compallcorrtrap(t,g,calf,cbet,calfbet,corr,jind)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Input: alf,bet,h - original grid
!            alfb,betb - coordinates of projection
!            d,calf,cbet - values that determine rho
!            n1,n2,m1,m2 - indeces of corners of window
!     Output: corr(j) = int H(j) - T[H(j)]
!     Note: Hj(alf,bet)=alf^p bet^q / rho^r, rho^2=d^2+calf*alf^2+cbet*bet^2
!           and
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : nbmax,rmax,pmax,nwinmax
use types
use old
implicit none
TYPE (basept), INTENT(IN) :: t
TYPE (grid),       INTENT(IN) :: g
real*8,            INTENT(IN) :: calf,cbet,calfbet
real*8,            INTENT(OUT) :: corr(*)
integer,           INTENT(OUT) :: jind(rmax,0:pmax,0:pmax)
! LOCAL
integer :: i,j,k,nb,n1,n2,m1,m2
integer, dimension(nbmax) :: r,p,q
real*8, dimension(nbmax) :: trap,inth,fact,F1,F2,F3,F4,term, F11, &
            dfu1,dfu2,dfuuu1,dfuuu2,dfv1,dfv2,dfvvv1,dfvvv2,dfuv1,dfuv2,dfuv3,dfuv4, &
            xjunk1,xjunk2,xjunk3,xjunk4
real*8 :: cu,cv,delu,delv,aa,bb,cc,dd,cu3,cv3,cuv,h,h3,h4,h5,weight,uu,vv,c
real*8, dimension(0:nwinmax) :: wtalf,wtbet

!print*
!print*,t%x0
!print*,calf,cbet,calfbet
!stop
!print*

cu=calf/t%d
cv=cbet/t%d
c=calfbet/(calf*cbet)
!print*,cu,cv,c
call compallindfact(t%d,cu,cv,nb,r,p,q,fact,jind)
call initsubgrid(t,g,n1,n2,m1,m2,wtalf,wtbet)

!print*,t%alfb,t%betb
trap(1:nb)=0
do i=n1,n2
do j=m1,m2
  weight=wtalf(i-n1)*wtbet(j-m1)
  if (t%roundoff) then
    if ((i.eq.t%i0).and.(j.eq.t%j0)) then
      weight=0
    endif
  endif
  call trapcomponents(g%alf(i)-t%alfb,g%bet(j)-t%betb,cu,cv,c,nb,r,p,q,trap,weight)
enddo
enddo
delu=g%h*cu
delv=g%h*cv
trap(1:nb)=trap(1:nb)*delu*delv

!print*,t%alfb,g%alf(n1),cu
aa=(g%alf(n1)-t%alfb)*cu
bb=(g%alf(n2)-t%alfb)*cu
cc=(g%bet(m1)-t%betb)*cv
dd=(g%bet(m2)-t%betb)*cv

!print*,t%x0
!print*,t%alfb
!print*,t%betb
!print*,g%bet(m1) 
!print*,t%betb
!print*,g%bet(m2)

!print*,aa,bb,cc,dd,c
!stop

h=g%h
h3=h**3
h4=h3*h
h5=h4*h
cu3=cu**3
cv3=cv**3
cuv=cu*cv
term(1:nb)=0
!print*,'c=',c,cu,cv,h3/12*cu
!print*,m2-m1,n2-n1
do j=m1,m2
  vv=(g%bet(j)-t%betb)*cv
!  call allder(aa,vv,dfu1,dfuuu1,xjunk1,xjunk2,xjunk3,jind)
!  call allder(bb,vv,dfu2,dfuuu2,xjunk1,xjunk2,xjunk3,jind)
  call allder_c(aa,vv,c,dfu1,dfuuu1,xjunk1,xjunk2,xjunk3,jind)
  call allder_c(bb,vv,c,dfu2,dfuuu2,xjunk1,xjunk2,xjunk3,jind)

  term=term + (-h3/12*cu*(dfu2-dfu1) + h5/720*cu3*(dfuuu2-dfuuu1))*wtbet(j-m1)
!  term=term + (-h3/12*cu*(dfu2-dfu1) )*wtbet(j-m1)
enddo
do i=n1,n2
  uu=(g%alf(i)-t%alfb)*cu
!  call allder(uu,cc,xjunk1,xjunk2,dfv1,dfvvv1,xjunk3,jind)
!  call allder(uu,dd,xjunk1,xjunk2,dfv2,dfvvv2,xjunk3,jind)

  call allder_c(uu,cc,c,xjunk1,xjunk2,dfv1,dfvvv1,xjunk3,jind)
  call allder_c(uu,dd,c,xjunk1,xjunk2,dfv2,dfvvv2,xjunk3,jind)
  term=term + (-h3/12*cv*(dfv2-dfv1) + h5/720*cv3*(dfvvv2-dfvvv1))*wtalf(i-n1)
!  term=term + (-h3/12*cv*(dfv2-dfv1) )*wtalf(i-n1)
enddo
!call allder(aa,cc,xjunk1,xjunk2,xjunk3,xjunk4,dfuv1,jind)
!call allder(bb,cc,xjunk1,xjunk2,xjunk3,xjunk4,dfuv2,jind)
!call allder(aa,dd,xjunk1,xjunk2,xjunk3,xjunk4,dfuv3,jind)
!call allder(bb,dd,xjunk1,xjunk2,xjunk3,xjunk4,dfuv4,jind)

call allder_c(aa,cc,c,xjunk1,xjunk2,xjunk3,xjunk4,dfuv1,jind)
call allder_c(bb,cc,c,xjunk1,xjunk2,xjunk3,xjunk4,dfuv2,jind)
call allder_c(aa,dd,c,xjunk1,xjunk2,xjunk3,xjunk4,dfuv3,jind)
call allder_c(bb,dd,c,xjunk1,xjunk2,xjunk3,xjunk4,dfuv4,jind)
term=term + h4/12**2*cuv*(dfuv4-dfuv2-dfuv3+dfuv1)

!call allFs(aa,cc,F1,jind)
!call allFs(bb,cc,F2,jind)
!call allFs(aa,dd,F3,jind)
!call allFs(bb,dd,F4,jind)
!do j=1,nb
!   inth(j) = F4(j) - F3(j) - F2(j) + F1(j)
!enddo

call allintFbackrecur2(aa,bb,cc,dd,c,inth,r,p,q,nb)
do j=1,nb
   corr(j)= fact(j)*( (inth(j)-trap(j))/cuv  - term(j) )  
!print*,j,corr(j),inth(j),trap(j),inth(j)-trap(j),term(j)*cuv ! ALL GOOD!!
enddo
!j=22
!print*,aa,bb,cc,dd
!print*,j,corr(j),inth(j),trap(j),inth(j)-trap(j),term(j)*cuv ! ALL GOOD!!
!print*,p(j),q(j),(r(j)-1)/2
!stop

return
END SUBROUTINE compallcorrtrap


SUBROUTINE trapcomponents(alfhat,bethat,cu,cv,c,nb,r,p,q,trap,wt)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c     adds a component of the gauss sum Hj(alfhat,bethat) for j=1,..,nb
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : pmax
implicit none
integer, intent(IN) :: nb,r(*),p(*),q(*)
real*8, intent(IN) :: alfhat,bethat,cu,cv,c,wt
real*8, intent(INOUT) :: trap(*)
!     LOCAL
integer :: j
real*8 :: u,v,uv,up(0:pmax),vp(0:pmax),rho

u=alfhat*cu
v=bethat*cv
uv=u*v

do j=0,pmax
  up(j)=u**j
  vp(j)=v**j
enddo
rho=sqrt(1+u**2 + 2*c*uv + v**2)

do j=1,nb
   trap(j)=trap(j) + up(p(j)) * vp(q(j)) / rho**r(j) * wt
enddo

return
END SUBROUTINE trapcomponents


SUBROUTINE compallindfact(d,cu,cv,nb,r,p,q,fact,jind)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     Input  : d,calf,dalf
!     Output : number of basis functions nb 
!              for j=1,nb: r(j),p(j),q(j),fact(j)
!     Hj(alf,bet)=alf^p bet^q / rho^r, 
!     rho^2=d^2+calf*alf^2+2*calfbet*alf*bet + cbet*bet^2
!
!     Hj(u,v) = u^p v^q / rho^r, rho^2= 1 + u^2 + 2*c*u*v +  v^2
!
!     and fact st int H(alf,bet)dalf dbet = fact*H(u,v)du,dv/ cu*cv
!     see line 164
!
!     These are all functions needed for O(h^4) method
!     alfj/r: j=0 to 2
!     alfj/r^3: j=0 to 4
!     alfj/r^5: j=0 to 6
!     alfj/r^7: j=3 to 8
!     alfj/r^9: j=6 to 10
!     alfj/r^11: j=9 to 12
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : rmax,pmax
implicit none
!     INPUT
real*8 :: d,cu,cv
!     OUTPUT
integer :: nb,r(*),p(*),q(*),jind(rmax,0:pmax,0:pmax)
real*8 :: fact(*)
!     LOCAL
integer :: ind,j,ipplusq,jj
integer :: rr,expd,expa,expb
 
!     Hrpq = alf^p bet^q / rho^r
!     fact=d**(p+q+2-r)/calf**(p+1)/cbet**(q+1)

ind=0
! H1s
rr=1
!      do ipplusq=0,0 !2nd order
!      do ipplusq=0,1 !3rd order
do ipplusq=0,2 !4th order
!      do ipplusq=0,3 !5th order
do j=0,ipplusq
   ind=ind+1
   r(ind)=rr
   p(ind)=ipplusq-j
   q(ind)=j
   jind(r(ind),p(ind),q(ind))=ind
enddo
enddo
! H3s
rr=3
!      do ipplusq=0,2 !2nd order
!      do ipplusq=0,3 !3rd order
do ipplusq=0,4 !4th order
!      do ipplusq=0,5 !5th order
do j=0,ipplusq
   ind=ind+1
   r(ind)=rr
   p(ind)=ipplusq-j
   q(ind)=j
   jind(r(ind),p(ind),q(ind))=ind
enddo
enddo
! H5s
rr=5
!      do ipplusq=0,4 !2nd order
!      do ipplusq=0,5 !3rd order
do ipplusq=0,6 !4th order
!      do ipplusq=0,7 !5th order
do j=0,ipplusq
   ind=ind+1
   r(ind)=rr
   p(ind)=ipplusq-j
   q(ind)=j
   jind(r(ind),p(ind),q(ind))=ind
enddo
enddo
! H7s
rr=7
!      do ipplusq=3,6 !2nd order
!      do ipplusq=3,7 !3rd order
do ipplusq=3,8 !4th order
!      do ipplusq=3,9 !5th order
do j=0,ipplusq
   ind=ind+1
   r(ind)=rr
   p(ind)=ipplusq-j
   q(ind)=j
   jind(r(ind),p(ind),q(ind))=ind
enddo
enddo
! H9s
rr=9
!      do ipplusq=6,0 !2nd order
!      do ipplusq=6,9 !3rd order
do ipplusq=6,10 !4th order
!      do ipplusq=6,11 !5th order
do j=0,ipplusq
   ind=ind+1
   r(ind)=rr
   p(ind)=ipplusq-j
   q(ind)=j
   jind(r(ind),p(ind),q(ind))=ind
enddo
enddo
! H11s
rr=11
!      do ipplusq=9,0 !3rd order
do ipplusq=9,12 !4th order
!      do ipplusq=9,13 !5th order
do j=0,ipplusq
   ind=ind+1
   r(ind)=rr
   p(ind)=ipplusq-j
   q(ind)=j
   jind(r(ind),p(ind),q(ind))=ind
enddo
enddo
! H13s
rr=13
do ipplusq=12,0 !4th order
!      do ipplusq=12,15 !5th order
do j=0,ipplusq
   ind=ind+1
   r(ind)=rr
   p(ind)=ipplusq-j
   q(ind)=j
   jind(r(ind),p(ind),q(ind))=ind
enddo
enddo


do j=1,ind
   fact(j) =1/(d**r(j)*cu**p(j)*cv**q(j))   
!  Note: cu=calf/d, cv=cbet/d so fact = 1/(d^(r-p-q)*calf^p*cbet^q)
enddo
nb=ind

return
END SUBROUTINE compallindfact


SUBROUTINE allder_c(u,v,c,fu,fuuu,fv,fvvv,fuv,jind)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Computes fu(u,v) for all basis functions 
! (does not output fp,fp3 since right now not needed)
! input : u,v (scalar)
! output : fu (vector)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : rmax,pmax
implicit none
integer, intent(IN) :: jind(rmax,0:pmax,0:pmax)
MODE, intent(IN) :: u,v,c
MODE, intent(OUT) :: fu(*),fuuu(*),fv(*),fvvv(*),fuv(*)
! LOCAL
MODE :: u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12, &
        v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12, &
        rho,rho2,rho3,rho5,rho7,rho9,rho11,rho13,rho15,rho17,uv,cuv,c2

c2=c**2

u2=u**2
u3=u2*u
u4=u2*u2
u5=u2*u3
u6=u2*u4
u7=u3*u4
u8=u4*u4
u9=u5*u4
u10=u5*u5
u11=u6*u5
u12=u6*u6

v2=v**2
v3=v2*v
v4=v2*v2
v5=v2*v3
v6=v2*v4
v7=v3*v4
v8=v4*v4
v9=v5*v4
v10=v5*v5
v11=v6*v5
v12=v6*v6

uv=u*v
cuv=c*uv

rho2 = 1+u2+2*c*uv+ v2
rho  = sqrt(rho2)
rho3 = rho2*rho
rho5 = rho2*rho3
rho7 = rho2*rho5
rho9 = rho2*rho7
rho11 = rho2*rho9
rho13 = rho2*rho11
rho15 = rho2*rho13
rho17 = rho2*rho15

! Hrpq = u^p v^q / rho^r

  fu(jind(1,0,0)) = (-u-c*v)/rho3
fuuu(jind(1,0,0)) = (-3*(u+c*v)*(-3+2*u2+4*c*uv+(-3+5*c2)*v2))/rho7
  fv(jind(1,0,0)) = (-(c*u)-v)/rho3
fvvv(jind(1,0,0)) = (-3*(c*u+v)*(-3+(-3+5*c2)*u2+4*c*uv+2*v2))/rho7
 fuv(jind(1,0,0)) = (3*uv+c2*uv+c*(-1+2*u2+2*v2))/rho5

  fu(jind(1,1,0)) = (1+c*uv+v2)/rho3
fuuu(jind(1,1,0)) = (3*(2*c*u3*v+2*u2*(2+(2+c2)*v2)+c*uv*(5+(5+c2)*v2)+(1+v2)*(-1+(-1+3*c2)*v2)))/rho7
  fv(jind(1,1,0)) = -((u*(c*u+v))/rho3)
fvvv(jind(1,1,0)) = (-3*u*(c*u+v)*(-3+(-3+5*c2)*u2+4*c*uv+2*v2))/rho7
 fuv(jind(1,1,0)) = (-(c2*u2*v)+c*u*(-2+u2-v2)-v*(1-2*u2+v2))/rho5

  fu(jind(1,0,1)) =  -((v*(u+c*v))/rho3)
fuuu(jind(1,0,1)) = (-3*v*(u+c*v)*(-3+2*u2+4*c*uv+(-3+5*c2)*v2))/rho7
  fv(jind(1,0,1)) = (1+u2+c*uv)/rho3
fvvv(jind(1,0,1)) = (3*(-1+(-1+3*c2)*u4+c*(5+c2)*u3*v+4*v2+c*uv*(5+2*v2)+u2*(-2+4*v2+c2*(3+2*v2))))/rho7
 fuv(jind(1,0,1)) = (-u3-c*u2*v+c*v*(-2+v2)-u*(1+(-2+c2)*v2))/rho5

  fu(jind(1,2,0)) = (u*(2+u2+3*c*uv+2*v2))/rho3
fuuu(jind(1,2,0)) = (-3*(2*c*v*(1+v2)**2+2*u*(1+v2)*(2+(2+c2)*v2)+c*u2*v*(5+(5+c2)*v2)+u3*(-1+(-1+3*c2)*v2)))/rho7
  fv(jind(1,2,0)) = -((u2*(c*u+v))/rho3)
fvvv(jind(1,2,0)) = (-3*u2*(c*u+v)*(-3+(-3+5*c2)*u2+4*c*uv+2*v2))/rho7
 fuv(jind(1,2,0)) = -((u*(3*c2*u2*v+v*(2-u2+2*v2)+c*u*(3+4*v2)))/rho5)

  fu(jind(1,1,1)) = (v*(1+c*uv+v2))/rho3
fuuu(jind(1,1,1)) = (3*v*(2*c*u3*v+2*u2*(2+(2+c2)*v2)+c*uv*(5+(5+c2)*v2)+(1+v2)*(-1+(-1+3*c2)*v2)))/rho7
  fv(jind(1,1,1)) = (u*(1+u2+c*uv))/rho3
fvvv(jind(1,1,1)) = (3*u*(-1+(-1+3*c2)*u4+c*(5+c2)*u3*v+4*v2+c*uv*(5+2*v2)+u2*(-2+4*v2+c2*(3+2*v2))))/rho7
 fuv(jind(1,1,1)) = (1+2*c*u3*v+v2+c*uv*(1+2*v2)+u2*(1+(3+c2)*v2))/rho5

  fu(jind(1,0,2)) = -((v2*(u+c*v))/rho3)
fuuu(jind(1,0,2)) = (-3*v2*(u+c*v)*(-3+2*u2+4*c*uv+(-3+5*c2)*v2))/rho7
  fv(jind(1,0,2)) = (v*(2+2*u2+3*c*uv+v2))/rho3
fvvv(jind(1,0,2)) = (-3*(c**3*u3*v2+(1+u2)*v*(4+4*u2-v2)+c2*u2*v*(2+2*u2+3*v2)+c*u*(1+u2)*(2+2*u2+5*v2)))/rho7
 fuv(jind(1,0,2)) = -((v*(2*u3+3*c*v+4*c*u2*v+u*(2+(-1+3*c2)*v2)))/rho5)

  fu(jind(3,0,0)) = (-3*(u+c*v))/rho5
fuuu(jind(3,0,0)) = (-15*(u+c*v)*(-3+4*u2+8*c*uv+(-3+7*c2)*v2))/rho9
  fv(jind(3,0,0)) = (-3*(c*u+v))/rho5
fvvv(jind(3,0,0)) = (-15*(c*u+v)*(-3+(-3+7*c2)*u2+8*c*uv+4*v2))/rho9
 fuv(jind(3,0,0)) = (3*(5*uv+3*c2*uv+c*(-1+4*u2+4*v2)))/rho7

  fu(jind(3,1,0)) = (1-2*u2-c*uv+v2)/rho5
fuuu(jind(3,1,0)) = (-3*(8*u4+12*c*u3*v+12*u2*(-2+(-2+c2)*v2)+c*uv*(-33+(-33+5*c2)*v2)-3*(1+v2)*(-1+(-1+5*c2)*v2)))/rho9
  fv(jind(3,1,0)) = (-3*u*(c*u+v))/rho5
fvvv(jind(3,1,0)) = (-15*u*(c*u+v)*(-3+(-3+7*c2)*u2+8*cuv+4*v2))/rho9
 fuv(jind(3,1,0)) = (3*(c2*u2*v-v*(1-4*u2+v2)+c*u*(-2+3*u2+v2)))/rho7

  fu(jind(3,0,1)) = (-3*v*(u+c*v))/rho5
fuuu(jind(3,0,1)) = (-15*v*(u+c*v)*(-3+4*u2+8*cuv+(-3+7*c2)*v2))/rho9
  fv(jind(3,0,1)) = (1+u2-cuv-2*v2)/rho5
fvvv(jind(3,0,1)) = (-3*(3+(3-15*c2)*u4+c*(-33+5*c2)*u3*v-24*v2+8*v4+3*cuv*(-11+4*v2)+3*u2*(2-8*v2+c2*(-5+4*v2))))/rho9
 fuv(jind(3,0,1)) = (3*(-u3+c*u2*v+c*v*(-2+3*v2)+u*(-1+(4+c2)*v2)))/rho7

  fu(jind(3,2,0)) = (u*(2-u2+cuv+2*v2))/rho5
fuuu(jind(3,2,0)) = (3*(-2*u5+6*c*u4*v-6*c*v*(1+v2)**2+6*u*(1+v2)*(-2+(-2+c2)*v2)+3*u3*(7+(7+c2)*v2)+c*u2*v*(15+(15+c2)*v2)))/rho9
  fv(jind(3,2,0)) = (-3*u2*(c*u+v))/rho5
fvvv(jind(3,2,0)) = (-15*u2*(c*u+v)*(-3+(-3+7*c2)*u2+8*cuv+4*v2))/rho9
 fuv(jind(3,2,0)) = (3*u*(3*u2*v-c2*u2*v+c*u*(-3+2*u2-2*v2)-2*v*(1+v2)))/rho7

  fu(jind(3,1,1)) = (v*(1-2*u2-cuv+v2))/rho5
fuuu(jind(3,1,1)) = (3*v*(-8*u4-12*c*u3*v+cuv*(33+(33-5*c2)*v2)-12*u2*(-2+(-2+c2)*v2)+3*(1+v2)*(-1+(-1+5*c2)*v2)))/rho9
  fv(jind(3,1,1)) = (u*(1+u2-cuv-2*v2))/rho5
fvvv(jind(3,1,1)) = (3*u*(-3+3*(-1+5*c2)*u4+c*(33-5*c2)*u3*v+24*v2-8*v4+3*cuv*(11-4*v2)-3*u2*(2-8*v2+c2*(-5+4*v2))))/rho9
 fuv(jind(3,1,1)) = (1-2*u4+4*c*u3*v-v2-2*v4+cuv*(-5+4*v2)+u2*(-1+(11+c2)*v2))/rho7

  fu(jind(3,0,2)) = (-3*v2*(u+c*v))/rho5
fuuu(jind(3,0,2)) = (-15*v2*(u+c*v)*(-3+4*u2+8*cuv+(-3+7*c2)*v2))/rho9
  fv(jind(3,0,2)) = (v*(2+2*u2+cuv-v2))/rho5
fvvv(jind(3,0,2)) = (3*c**3*u3*v2+9*c2*u2*v*(2+2*u2+v2)-9*c*u*(2+2*u4-5*v2-2*v4+u2*(4-5*v2))-3*v*(12+12*u4-21*v2+2*v4 &
                     -3*u2*(-8+7*v2)))/rho9
 fuv(jind(3,0,2)) = (-3*v*(2*u3+2*c*u2*v+c*v*(3-2*v2)+u*(2+(-3+c2)*v2)))/rho7

! ifdef THIRD
  fu(jind(3,3,0)) = (3*u2*(1+cuv+v2))/rho5
fuuu(jind(3,3,0)) = (-3*(-6*c*u5*v+6*cuv*(1+v2)**2-2*(1+v2)**3+6*u4*(-2+(-2+c2)*v2)+3*u2*(1+v2)*(7+(7+c2)*v2) &
                    +c*u3*v*(15+(15+c2)*v2)))/rho9
  fv(jind(3,3,0)) = (-3*u3*(c*u+v))/rho5
fvvv(jind(3,3,0)) = (-15*u3*(c*u+v)*(-3+(-3+7*c2)*u2+8*cuv+4*v2))/rho9
 fuv(jind(3,3,0)) = (3*u2*(2*u2*v-3*c2*u2*v+c*u*(-4+u2-5*v2)-3*v*(1+v2)))/rho7

  fu(jind(3,2,1)) = (u*v*(2-u2+cuv+2*v2))/rho5
fuuu(jind(3,2,1)) = (-3*v*(2*u5-6*c*u4*v+6*c*v*(1+v2)**2-6*u*(1+v2)*(-2+(-2+c2)*v2)-3*u3*(7+(7+c2)*v2) &
                    -c*u2*v*(15+(15+c2)*v2)))/rho9
  fv(jind(3,2,1)) = (u2*(1+u2-cuv-2*v2))/rho5
fvvv(jind(3,2,1)) = (3*u2*(-3+3*(-1+5*c2)*u4+c*(33-5*c2)*u3*v+24*v2-8*v4+3*cuv*(11-4*v2)-3*u2*(2-8*v2+c2*(-5+4*v2))))/rho9
 fuv(jind(3,2,1)) = -((u*(u4-5*c*u3*v+cuv*(4+v2)+u2*(-1+(-10+c2)*v2)+2*(-1+v2+2*v4)))/rho7)

  fu(jind(3,1,2)) = (v2*(1-2*u2-cuv+v2))/rho5
fuuu(jind(3,1,2)) = (3*v2*(-8*u4-12*c*u3*v+cuv*(33+(33-5*c2)*v2)-12*u2*(-2+(-2+c2)*v2)+3*(1+v2)*(-1+(-1+5*c2)*v2)))/rho9
  fv(jind(3,1,2)) = (u*v*(2+2*u2+cuv-v2))/rho5
fvvv(jind(3,1,2)) = (-3*u*(-(c**3*u3*v2)-3*c2*u2*v*(2+2*u2+v2)+3*c*u*(2+2*u4-5*v2-2*v4+u2*(4-5*v2)) &
                    +v*(12+12*u4-21*v2+2*v4-3*u2*(-8+7*v2))))/rho9
 fuv(jind(3,1,2)) = -((v*(-2+4*u4+c*u3*v-v2+v4+cuv*(4-5*v2)+u2*(2+(-10+c2)*v2)))/rho7)

  fu(jind(3,0,3)) = (-3*v3*(u+c*v))/rho5
fuuu(jind(3,0,3)) = (-15*v3*(u+c*v)*(-3+4*u2+8*cuv+(-3+7*c2)*v2))/rho9
  fv(jind(3,0,3)) = (3*v2*(1+u2+cuv))/rho5
fvvv(jind(3,0,3)) = (-3*(-2-2*u6+6*c*u5*v+21*v2-12*v4+3*u4*(-2+(7+c2)*v2)+c*u3*v*(12+(15+c2)*v2)+3*cuv*(2+5*v2-2*v4) &
                    +3*u2*(-2+(14+c2)*v2+2*(-2+c2)*v4)))/rho9
 fuv(jind(3,0,3)) = (-3*v2*(3*u3+5*c*u2*v-c*v*(-4+v2)+u*(3+(-2+3*c2)*v2)))/rho7

  fu(jind(3,4,0)) = (u3*(4+u2+5*cuv+4*v2))/rho5
fuuu(jind(3,4,0)) = (3*u*(12*cuv*(1+v2)**2+8*(1+v2)**3+u4*(3+(3-15*c2)*v2)+12*u2*(1+v2)*(-2+(-2+c2)*v2) &
                    +c*u3*v*(-33+(-33+5*c2)*v2)))/rho9
  fv(jind(3,4,0)) = (-3*u4*(c*u+v))/rho5
fvvv(jind(3,4,0)) = (-15*u4*(c*u+v)*(-3+(-3+7*c2)*u2+8*cuv+4*v2))/rho9
 fuv(jind(3,4,0)) = (-3*u3*(5*c2*u2*v+v*(4-u2+4*v2)+c*u*(5+8*v2)))/rho7

  fu(jind(3,3,1)) = (3*u2*v*(1+cuv+v2))/rho5
fuuu(jind(3,3,1)) = (3*v*(6*c*u5*v-6*cuv*(1+v2)**2+2*(1+v2)**3-6*u4*(-2+(-2+c2)*v2)-3*u2*(1+v2)*(7+(7+c2)*v2) &
                    -c*u3*v*(15+(15+c2)*v2)))/rho9
  fv(jind(3,3,1)) = (u3*(1+u2-cuv-2*v2))/rho5
fvvv(jind(3,3,1)) = (3*u3*(-3+3*(-1+5*c2)*u4+c*(33-5*c2)*u3*v+24*v2-8*v4+3*cuv*(11-4*v2)-3*u2*(2-8*v2+c2*(-5+4*v2))))/rho9
 fuv(jind(3,3,1)) = (3*u2*(1+2*c*u3*v-v2-2*v4-cuv*(1+2*v2)+u2*(1-(-3+c2)*v2)))/rho7

  fu(jind(3,2,2)) = (u*v2*(2-u2+cuv+2*v2))/rho5
fuuu(jind(3,2,2)) = (-3*v2*(2*u5-6*c*u4*v+6*c*v*(1+v2)**2-6*u*(1+v2)*(-2+(-2+c2)*v2)-3*u3*(7+(7+c2)*v2) &
                    -c*u2*v*(15+(15+c2)*v2)))/rho9
  fv(jind(3,2,2)) = (u2*v*(2+2*u2+cuv-v2))/rho5
fvvv(jind(3,2,2)) = (-3*u2*(-(c**3*u3*v2)-3*c2*u2*v*(2+2*u2+v2)+3*c*u*(2+2*u4-5*v2-2*v4+u2*(4-5*v2)) &
                    +v*(12+12*u4-21*v2+2*v4-3*u2*(-8+7*v2))))/rho9
 fuv(jind(3,2,2)) = (u*v*(4-2*u4+4*c*u3*v+2*v2-2*v4+cuv*(1+4*v2)+u2*(2+(11+c2)*v2)))/rho7

  fu(jind(3,1,3)) = (v3*(1-2*u2-cuv+v2))/rho5
fuuu(jind(3,1,3)) = (3*v3*(-8*u4-12*c*u3*v+cuv*(33+(33-5*c2)*v2)-12*u2*(-2+(-2+c2)*v2)+3*(1+v2)*(-1+(-1+5*c2)*v2)))/rho9
  fv(jind(3,1,3)) = (3*u*v2*(1+u2+cuv))/rho5
fvvv(jind(3,1,3)) = (3*u*(2+2*u6-6*c*u5*v-21*v2+12*v4-3*u4*(-2+(7+c2)*v2)-c*u3*v*(12+(15+c2)*v2)+3*cuv*(-2-5*v2+2*v4) &
                    -3*u2*(-2+(14+c2)*v2+2*(-2+c2)*v4)))/rho9
 fuv(jind(3,1,3)) = (-3*v2*(-1+2*u4+2*c*u3*v-v2+u2*(1+(-3+c2)*v2)+c*u*(v-2*v3)))/rho7

  fu(jind(3,0,4)) = (-3*v4*(u+c*v))/rho5
fuuu(jind(3,0,4)) = (-15*v4*(u+c*v)*(-3+4*u2+8*cuv+(-3+7*c2)*v2))/rho9
  fv(jind(3,0,4)) = (v3*(4+4*u2+5*cuv+v2))/rho5
fvvv(jind(3,0,4)) = (3*v*(8+8*u6+12*c*u5*v-24*v2+3*v4+3*cuv*(4-11*v2)+12*u4*(2+(-2+c2)*v2)+c*u3*v*(24+(-33+5*c2)*v2) &
                    -3*u2*(-8-4*(-4+c2)*v2+(-1+5*c2)*v4)))/rho9
 fuv(jind(3,0,4)) = (-3*v3*(4*u3+5*c*v+8*c*u2*v+u*(4+(-1+5*c2)*v2)))/rho7

  fu(jind(5,0,0)) = (-5*(u+c*v))/rho7
fuuu(jind(5,0,0)) = (-105*(u+c*v)*(-1+2*u2+4*cuv+(-1+3*c2)*v2))/rho11
  fv(jind(5,0,0)) = (-5*(c*u+v))/rho7
fvvv(jind(5,0,0)) = (-105*(c*u+v)*(-1+(-1+3*c2)*u2+4*cuv+2*v2))/rho11
 fuv(jind(5,0,0)) = (5*(7*u*v+5*c2*u*v+c*(-1+6*u2+6*v2)))/rho9

  fu(jind(5,1,0)) = (1-4*u2-3*cuv+v2)/rho7
fuuu(jind(5,1,0)) = (-15*(8*u4+18*c*u3*v+6*u2*(-2+(-2+3*c2)*v2)+cuv*(-17+(-17+7*c2)*v2)-(1+v2)*(-1+(-1+7*c2)*v2)))/rho11
  fv(jind(5,1,0)) = (-5*u*(c*u+v))/rho7
fvvv(jind(5,1,0)) = (-105*u*(c*u+v)*(-1+(-1+3*c2)*u2+4*cuv+2*v2))/rho11
 fuv(jind(5,1,0)) = (5*(3*c2*u2*v-v*(1-6*u2+v2)+c*u*(-2+5*u2+3*v2)))/rho9

  fu(jind(5,0,1)) = (-5*v*(u+c*v))/rho7
fuuu(jind(5,0,1)) = (-105*v*(u+c*v)*(-1+2*u2+4*cuv+(-1+3*c2)*v2))/rho11
  fv(jind(5,0,1)) = (1+u2-3*cuv-4*v2)/rho7
fvvv(jind(5,0,1)) = (-15*(1+(1-7*c2)*u4+c*(-17+7*c2)*u3*v-12*v2+8*v4+cuv*(-17+18*v2)+u2*(2-12*v2+c2*(-7+18*v2))))/rho11
 fuv(jind(5,0,1)) = (-5*(u+u3-3*c*u2*v-3*(2+c2)*u*v2+c*v*(2-5*v2)))/rho9

  fu(jind(5,2,0)) = (u*(2-3*u2-cuv+2*v2))/rho7
fuuu(jind(5,2,0)) = (-15*(4*u5+4*c*u4*v+2*c*v*(1+v2)**2+c*u2*v*(-15+(-15+c2)*v2)+u3*(-13+(-13+3*c2)*v2) &
                    -2*u*(1+v2)*(-2+(-2+3*c2)*v2)))/rho11
  fv(jind(5,2,0)) = (-5*u2*(c*u+v))/rho7
fvvv(jind(5,2,0)) = (-105*u2*(c*u+v)*(-1+(-1+3*c2)*u2+4*cuv+2*v2))/rho11
 fuv(jind(5,2,0)) = (5*u*(c*u*(-3+4*u2)+5*u2*v+c2*u2*v-2*v*(1+v2)))/rho9

  fu(jind(5,1,1)) = (v*(1-4*u2-3*cuv+v2))/rho7
fuuu(jind(5,1,1)) = (15*v*(-8*u4-18*c*u3*v+cuv*(17+(17-7*c2)*v2)-6*u2*(-2+(-2+3*c2)*v2)+(1+v2)*(-1+(-1+7*c2)*v2)))/rho11
  fv(jind(5,1,1)) = (u*(1+u2-3*cuv-4*v2))/rho7
fvvv(jind(5,1,1)) = (15*u*(-1+(-1+7*c2)*u4+c*(17-7*c2)*u3*v+12*v2-8*v4+cuv*(17-18*v2)+u2*(-2+12*v2+c2*(7-18*v2))))/rho11
 fuv(jind(5,1,1)) = (1-4*u4+14*c*u3*v-3*v2-4*v4+cuv*(-11+14*v2)+3*u2*(-1+3*(3+c2)*v2))/rho9

  fu(jind(5,0,2)) = (-5*v2*(u+c*v))/rho7
fuuu(jind(5,0,2)) = (-105*v2*(u+c*v)*(-1+2*u2+4*cuv+(-1+3*c2)*v2))/rho11
  fv(jind(5,0,2)) = (v*(2+2*u2-cuv-3*v2))/rho7
fvvv(jind(5,0,2)) = (-15*(c**3*u3*v2+3*c2*u2*v*(-2-2*u2+v2)+c*u*(2+2*u4-15*v2+4*v4+u2*(4-15*v2)) &
                    +v*(4+4*u4-13*v2+4*v4+u2*(8-13*v2))))/rho11
 fuv(jind(5,0,2)) = (5*v*(-2*u3+c*v*(-3+4*v2)+u*(-2+(5+c2)*v2)))/rho9

  fu(jind(5,3,0)) = (u2*(3-2*u2+cuv+3*v2))/rho7
fuuu(jind(5,3,0)) = (3*(-8*u6+12*c*u5*v-18*cuv*(1+v2)**2+2*(1+v2)**3+4*u4*(14+(14+c2)*v2) &
                    +c*u3*v*(29+(29+c2)*v2)+3*u2*(1+v2)*(-13+(-13+3*c2)*v2)))/rho11
  fv(jind(5,3,0)) = (-5*u3*(c*u+v))/rho7 
fvvv(jind(5,3,0)) = (-105*u3*(c*u+v)*(-1+(-1+3*c2)*u2+4*cuv+2*v2))/rho11 
fuv(jind(5,3,0)) = (5*u2*(4*u2*v-c2*u2*v+c*u*(-4+3*u2-3*v2)-3*v*(1+v2)))/rho9

  fu(jind(5,2,1)) = (u*v*(2-3*u2-cuv+2*v2))/rho7
fuuu(jind(5,2,1)) = (-15*v*(4*u5+4*c*u4*v+2*c*v*(1+v2)**2+c*u2*v*(-15+(-15+c2)*v2) &
                    +u3*(-13+(-13+3*c2)*v2)-2*u*(1+v2)*(-2+(-2+3*c2)*v2)))/rho11
  fv(jind(5,2,1)) = (u2*(1+u2-3*cuv-4*v2))/rho7
fvvv(jind(5,2,1)) = (15*u2*(-1+(-1+7*c2)*u4+c*(17-7*c2)*u3*v+12*v2-8*v4+cuv*(17-18*v2)+u2*(-2+12*v2+c2*(7-18*v2))))/rho11
 fuv(jind(5,2,1)) = -((u*(-2+3*u4-13*c*u3*v+6*v2+8*v4-3*cuv*(-4+v2)+u2*(1-3*(8+c2)*v2)))/rho9)

  fu(jind(5,1,2)) = (v2*(1-4*u2-3*cuv+v2))/rho7
fuuu(jind(5,1,2)) = (15*v2*(-8*u4-18*c*u3*v+cuv*(17+(17-7*c2)*v2)-6*u2*(-2+(-2+3*c2)*v2)+(1+v2)*(-1+(-1+7*c2)*v2)))/rho11
  fv(jind(5,1,2)) = (u*v*(2+2*u2-cuv-3*v2))/rho7
fvvv(jind(5,1,2)) = (-15*u*(c**3*u3*v2+3*c2*u2*v*(-2-2*u2+v2)+c*u*(2+2*u4-15*v2+4*v4+u2*(4-15*v2)) &
                    +v*(4+4*u4-13*v2+4*v4+u2*(8-13*v2))))/rho11
 fuv(jind(5,1,2)) = -((v*(-2+8*u4-3*c*u3*v+v2+3*v4+cuv*(12-13*v2)-3*u2*(-2+(8+c2)*v2)))/rho9)

  fu(jind(5,0,3)) = (-5*v3*(u + c*v))/rho7
fuuu(jind(5,0,3)) = (-105*v3*(u + c*v)*(-1 + 2*u2 + 4*cuv + (-1 + 3*c2)*v2))/rho11
  fv(jind(5,0,3)) = (v2*(3 + 3*u2 + cuv - 2*v2))/rho7
fvvv(jind(5,0,3)) = (3*(2 + 2*u6 - 18*c*u5*v - 39*v2 + 56*v4 - 8*v6 + c*u3*v*(-36 + (29 + c2)*v2) + u4*(6 + (-39 + 9*c2)*v2)  &
                    + cuv*(-18 + 29*v2 + 12*v4) + u2*(6 + (-78 + 9*c2)*v2 + 4*(14 + c2)*v4)))/rho11 
 fuv(jind(5,0,3)) = (-5*v2*(3*u3 + 3*c*u2*v + c*v*(4 - 3*v2) + u*(3 + (-4 + c2)*v2)))/rho9

  fu(jind(5,4,0)) = (u3*(4-u2+3*cuv+4*v2))/rho7
fuuu(jind(5,4,0)) = (-3*u*(2*u6-18*c*u5*v+12*cuv*(1+v2)**2-8*(1+v2)**3+4*u2*(1+v2)*(14+(14+c2)*v2)+c*u3*v*(29+(29+c2)*v2) &
                    +u4*(-39+(-39+9*c2)*v2)))/rho11
  fv(jind(5,4,0)) = (-5*u4*(c*u+v))/rho7
fvvv(jind(5,4,0)) = (-105*u4*(c*u+v)*(-1+(-1+3*c2)*u2+4*cuv+2*v2))/rho11
 fuv(jind(5,4,0)) = (5*u3*(3*u2*v-3*c2*u2*v+c*u*(-5+2*u2-6*v2)-4*v*(1+v2)))/rho9

  fu(jind(5,3,1)) = (u2*v*(3-2*u2+cuv+3*v2))/rho7
fuuu(jind(5,3,1)) = (3*v*(-8*u6+12*c*u5*v-18*cuv*(1+v2)**2+2*(1+v2)**3+4*u4*(14+(14+c2)*v2)+c*u3*v*(29+(29+c2)*v2) &
                    +3*u2*(1+v2)*(-13+(-13+3*c2)*v2)))/rho11
  fv(jind(5,3,1)) = (u3*(1+u2-3*cuv-4*v2))/rho7
fvvv(jind(5,3,1)) = (15*u3*(-1+(-1+7*c2)*u4+c*(17-7*c2)*u3*v+12*v2-8*v4+cuv*(17-18*v2)+u2*(-2+12*v2+c2*(7-18*v2))))/rho11
 fuv(jind(5,3,1)) = -((u2*(-3+2*u4-12*c*u3*v+9*v2+12*v4+cuv*(13+8*v2)+u2*(-1+3*(-7+c2)*v2)))/rho9)

  fu(jind(5,2,2)) = (u*v2*(2-3*u2-cuv+2*v2))/rho7
fuuu(jind(5,2,2)) = (-15*v2*(4*u5+4*c*u4*v+2*c*v*(1+v2)**2+c*u2*v*(-15+(-15+c2)*v2)+u3*(-13+(-13+3*c2)*v2) &
                    -2*u*(1+v2)*(-2+(-2+3*c2)*v2)))/rho11
  fv(jind(5,2,2)) = (u2*v*(2+2*u2-cuv-3*v2))/rho7
fvvv(jind(5,2,2)) = (-15*u2*(c**3*u3*v2+3*c2*u2*v*(-2-2*u2+v2)+c*u*(2+2*u4-15*v2+4*v4+u2*(4-15*v2))+v*(4+4*u4-13*v2+4*v4 &
                    +u2*(8-13*v2))))/rho11
 fuv(jind(5,2,2)) = (u*v*(-6*u4+6*c*u3*v+3*cuv*(-3+2*v2)+u2*(-2+(23+c2)*v2)-2*(-2+v2+3*v4)))/rho9

  fu(jind(5,1,3)) = (v3*(1-4*u2-3*cuv+v2))/rho7
fuuu(jind(5,1,3)) = (15*v3*(-8*u4-18*c*u3*v+cuv*(17+(17-7*c2)*v2)-6*u2*(-2+(-2+3*c2)*v2)+(1+v2)*(-1+(-1+7*c2)*v2)))/rho11
  fv(jind(5,1,3)) = (u*v2*(3+3*u2+cuv-2*v2))/rho7
fvvv(jind(5,1,3)) = (3*u*(2+2*u6-18*c*u5*v-39*v2+56*v4-8*v6+c*u3*v*(-36+(29+c2)*v2)+u4*(6+(-39+9*c2)*v2)+cuv*(-18+29*v2+12*v4) &
                    +u2*(6+(-78+9*c2)*v2+4*(14+c2)*v4)))/rho11
 fuv(jind(5,1,3)) = (v2*(3-12*u4-8*c*u3*v+v2-2*v4+cuv*(-13+12*v2)-3*u2*(3+(-7+c2)*v2)))/rho9

  fu(jind(5,0,4)) = (-5*v4*(u + c*v))/rho7
fuuu(jind(5,0,4)) = (-105*v4*(u + c*v)*(-1 + 2*u2 + 4*cuv + (-1 + 3*c2)*v2))/rho11
  fv(jind(5,0,4)) = (v3*(4+4*u2+3*cuv-v2))/rho7
fvvv(jind(5,0,4)) = (-3*v*(-8-8*u6+12*c*u5*v+56*v2-39*v4+2*v6+4*u4*(-6+(14+c2)*v2)+c*u3*v*(24+(29+c2)*v2)+cuv*(12+29*v2-18*v4) &
                    +u2*(-24+4*(28+c2)*v2+(-39+9*c2)*v4)))/rho11
 fuv(jind(5,0,4)) = (-5*v3*(4*u3+6*c*u2*v+c*v*(5-2*v2)+u*(4+3*(-1+c2)*v2)))/rho9

  fu(jind(5,5,0)) = (5*u4*(1+cuv+v2))/rho7
fuuu(jind(5,5,0)) = (15*u2*(2*c*u5*v+4*cuv*(1+v2)**2+4*(1+v2)**3+u4*(4+(4-6*c2)*v2)+c*u3*v*(-15+(-15+c2)*v2) &
                    +u2*(1+v2)*(-13+(-13+3*c2)*v2)))/rho11
 fv(jind(5,5,0)) = (-5*u5*(c*u+v))/rho7
fvvv(jind(5,5,0)) = (-105*u5*(c*u+v)*(-1+(-1+3*c2)*u2+4*cuv+2*v2))/rho11
 fuv(jind(5,5,0)) = (5*u4*(2*u2*v-5*c2*u2*v+c*u*(-6+u2-9*v2)-5*v*(1+v2)))/rho9

  fu(jind(5,4,1)) = (u3*v*(4 - u2 + 3*cuv + 4*v2))/rho7
fuuu(jind(5,4,1)) = (-3*u*v*(2*u6 - 18*c*u5*v + 12*cuv*(1 + v2)**2 - 8*(1 + v2)**3 + 4*u2*(1 + v2)*(14 + (14 + c2)*v2) &
                     + c*u3*v*(29 + (29 + c2)*v2) + u4*(-39 + (-39 + 9*c2)*v2)))/rho11
  fv(jind(5,4,1)) = (u4*(1+u2-3*cuv-4*v2))/rho7
fvvv(jind(5,4,1)) = (15*u4*(-1+(-1+7*c2)*u4+c*(17-7*c2)*u3*v+12*v2-8*v4+cuv*(17-18*v2)+u2*(-2+12*v2+c2*(7-18*v2))))/rho11
 fuv(jind(5,4,1)) = -((u3*(u4-11*c*u3*v+cuv*(14+19*v2)+3*u2*(-1+3*(-2+c2)*v2)+4*(-1+3*v2+4*v4)))/rho9)

  fu(jind(5,3,2)) = (u2*v2*(3-2*u2+cuv+3*v2))/rho7
fuuu(jind(5,3,2)) = (3*v2*(-8*u6+12*c*u5*v-18*cuv*(1+v2)**2+2*(1+v2)**3+4*u4*(14+(14+c2)*v2)+c*u3*v*(29+(29+c2)*v2) &
                    +3*u2*(1+v2)*(-13+(-13+3*c2)*v2)))/rho11
  fv(jind(5,3,2)) = (u3*v*(2+2*u2-cuv-3*v2))/rho7
fvvv(jind(5,3,2)) = (-15*u3*(c**3*u3*v2+3*c2*u2*v*(-2-2*u2+v2)+c*u*(2+2*u4-15*v2+4*v4+u2*(4-15*v2)) &
                    +v*(4+4*u4-13*v2+4*v4+u2*(8-13*v2))))/rho11
 fuv(jind(5,3,2)) = -((u2*v*(4*u4-9*c*u3*v+cuv*(6+v2)+u2*(-2+(-22+c2)*v2)+3*(-2+v2+3*v4)))/rho9)

  fu(jind(5,2,3)) = (u*v3*(2-3*u2-cuv+2*v2))/rho7
fuuu(jind(5,2,3)) = (-15*v3*(4*u5+4*c*u4*v+2*c*v*(1+v2)**2+c*u2*v*(-15+(-15+c2)*v2)+u3*(-13+(-13+3*c2)*v2) &
                    -2*u*(1+v2)*(-2+(-2+3*c2)*v2)))/rho11
  fv(jind(5,2,3)) = (u2*v2*(3+3*u2+cuv-2*v2))/rho7
fvvv(jind(5,2,3)) = (3*u2*(2+2*u6-18*c*u5*v-39*v2+56*v4-8*v6+c*u3*v*(-36+(29+c2)*v2)+u4*(6+(-39+9*c2)*v2) &
                    +cuv*(-18+29*v2+12*v4)+u2*(6+(-78+9*c2)*v2+4*(14+c2)*v4)))/rho11
 fuv(jind(5,2,3)) = -((u*v2*(-6+9*u4+c*u3*v-2*v2+4*v4+u2*(3+(-22+c2)*v2)+c*u*(6*v-9*v3)))/rho9)

  fu(jind(5,1,4)) = (v4*(1-4*u2-3*cuv+v2))/rho7
fuuu(jind(5,1,4)) = (15*v4*(-8*u4-18*c*u3*v+cuv*(17+(17-7*c2)*v2)-6*u2*(-2+(-2+3*c2)*v2)+(1+v2)*(-1+(-1+7*c2)*v2)))/rho11
  fv(jind(5,1,4)) = (u*v3*(4+4*u2+3*cuv-v2))/rho7
fvvv(jind(5,1,4)) = (-3*u*v*(-8-8*u6+12*c*u5*v+56*v2-39*v4+2*v6+4*u4*(-6+(14+c2)*v2)+c*u3*v*(24+(29+c2)*v2) &
                    +cuv*(12+29*v2-18*v4)+u2*(-24+4*(28+c2)*v2+(-39+9*c2)*v4)))/rho11
 fuv(jind(5,1,4)) = -((v3*(-4+16*u4+19*c*u3*v-3*v2+v4+cuv*(14-11*v2)+3*u2*(4+3*(-2+c2)*v2)))/rho9)

  fu(jind(5,0,5)) = (-5*v5*(u+c*v))/rho7
fuuu(jind(5,0,5)) = (-105*v5*(u+c*v)*(-1+2*u2+4*cuv+(-1+3*c2)*v2))/rho11
  fv(jind(5,0,5)) = (5*v4*(1+u2+cuv))/rho7
fvvv(jind(5,0,5)) = (15*v2*(4+4*u6+4*c*u5*v-13*v2+4*v4+c*u3*v*(8+(-15+c2)*v2)+u4*(12+(-13+3*c2)*v2) &
                    -u2*(-1+2*v2)*(12+(-2+3*c2)*v2)+cuv*(4-15*v2+2*v4)))/rho11
 fuv(jind(5,0,5)) = (-5*v4*(5*u3+9*c*u2*v-c*v*(-6+v2)+u*(5+(-2+5*c2)*v2)))/rho9

! ifdef FOURTH
  fu(jind(5,6,0)) = (u5*(6+u2+7*cuv+6*v2))/rho7
fuuu(jind(5,6,0)) = (15*u3*(18*cuv*(1+v2)**2+8*(1+v2)**3+u4*(1+(1-7*c2)*v2)+6*u2*(1+v2)*(-2+(-2+3*c2)*v2) &
                    +c*u3*v*(-17+(-17+7*c2)*v2)))/rho11
  fv(jind(5,6,0)) = (-5*u6*(c*u+v))/rho7
fvvv(jind(5,6,0)) = (-105*u6*(c*u+v)*(-1+(-1+3*c2)*u2+4*cuv+2*v2))/rho11
 fuv(jind(5,6,0)) = (-5*u5*(7*c2*u2*v+v*(6-u2+6*v2)+c*u*(7+12*v2)))/rho9

  fu(jind(5,5,1)) = (5*u4*v*(1+cuv+v2))/rho7
fuuu(jind(5,5,1)) = (15*u2*v*(2*c*u5*v+4*cuv*(1+v2)**2+4*(1+v2)**3+u4*(4+(4-6*c2)*v2)+c*u3*v*(-15+(-15+c2)*v2) &
                    +u2*(1+v2)*(-13+(-13+3*c2)*v2)))/rho11
  fv(jind(5,5,1)) = (u5*(1+u2-3*cuv-4*v2))/rho7
fvvv(jind(5,5,1)) = (15*u5*(-1+(-1+7*c2)*u4+c*(17-7*c2)*u3*v+12*v2-8*v4+cuv*(17-18*v2)+u2*(-2+12*v2+c2*(7-18*v2))))/rho11
 fuv(jind(5,5,1)) = (5*u4*(1+2*c*u3*v-3*v2-4*v4-3*cuv*(1+2*v2)+u2*(1-3*(-1+c2)*v2)))/rho9


  fu(jind(5,4,2)) = (u3*v2*(4-u2+3*cuv+4*v2))/rho7
fuuu(jind(5,4,2)) = (-3*u*v2*(2*u6-18*c*u5*v+12*cuv*(1+v2)**2-8*(1+v2)**3+4*u2*(1+v2)*(14+(14+c2)*v2) &
                    +c*u3*v*(29+(29+c2)*v2)+u4*(-39+(-39+9*c2)*v2)))/rho11
  fv(jind(5,4,2)) = (u4*v*(2+2*u2-cuv-3*v2))/rho7
fvvv(jind(5,4,2)) = (-15*u4*(c**3*u3*v2+3*c2*u2*v*(-2-2*u2+v2)+c*u*(2+2*u4-15*v2+4*v4+u2*(4-15*v2)) &
                    +v*(4+4*u4-13*v2+4*v4+u2*(8-13*v2))))/rho11
 fuv(jind(5,4,2)) = -((u3*v*(2*u4-12*c*u3*v+cuv*(3+8*v2)+3*u2*(-2+(-7+c2)*v2)+4*(-2+v2+3*v4)))/rho9)

  fu(jind(5,3,3)) = (u2*v3*(3-2*u2+cuv+3*v2))/rho7
fuuu(jind(5,3,3)) = (3*v3*(-8*u6+12*c*u5*v-18*cuv*(1+v2)**2+2*(1+v2)**3+4*u4*(14+(14+c2)*v2)+c*u3*v*(29+(29+c2)*v2) &
                    +3*u2*(1+v2)*(-13+(-13+3*c2)*v2)))/rho11
  fv(jind(5,3,3)) = (u3*v2*(3+3*u2+cuv-2*v2))/rho7
fvvv(jind(5,3,3)) = (3*u3*(2+2*u6-18*c*u5*v-39*v2+56*v4-8*v6+c*u3*v*(-36+(29+c2)*v2)+u4*(6+(-39+9*c2)*v2) &
                    +cuv*(-18+29*v2+12*v4)+u2*(6+(-78+9*c2)*v2+4*(14+c2)*v4)))/rho11
 fuv(jind(5,3,3)) = (u2*v2*(9-6*u4+6*c*u3*v+3*v2-6*v4+cuv*(1+6*v2)+u2*(3+(23+c2)*v2)))/rho9

  fu(jind(5,2,4)) = (u*v4*(2-3*u2-cuv+2*v2))/rho7
fuuu(jind(5,2,4)) = (-15*v4*(4*u5+4*c*u4*v+2*c*v*(1+v2)**2+c*u2*v*(-15+(-15+c2)*v2)+u3*(-13+(-13+3*c2)*v2) &
                    -2*u*(1+v2)*(-2+(-2+3*c2)*v2)))/rho11
  fv(jind(5,2,4)) = (u2*v3*(4+4*u2+3*cuv-v2))/rho7
fvvv(jind(5,2,4)) = (-3*u2*v*(-8-8*u6+12*c*u5*v+56*v2-39*v4+2*v6+4*u4*(-6+(14+c2)*v2)+c*u3*v*(24+(29+c2)*v2) &
                    +cuv*(12+29*v2-18*v4)+u2*(-24+4*(28+c2)*v2+(-39+9*c2)*v4)))/rho11
 fuv(jind(5,2,4)) = -((u*v3*(12*u4+8*c*u3*v+3*cuv*(1-4*v2)+u2*(4+3*(-7+c2)*v2)+2*(-4-3*v2+v4)))/rho9)

  fu(jind(5,1,5)) = (v5*(1 - 4*u2 - 3*cuv + v2))/rho7 
fuuu(jind(5,1,5)) = (15*v5*(-8*u4 - 18*c*u3*v + cuv*(17 + (17 - 7*c2)*v2) - 6*u2*(-2 + (-2 + 3*c2)*v2)  &
                    + (1 + v2)*(-1 + (-1 + 7*c2)*v2)))/rho11
  fv(jind(5,1,5)) = (5*u*v4*(1 + u2 + cuv))/rho7
fvvv(jind(5,1,5)) = (15*u*v2*(4+4*u6+4*c*u5*v-13*v2+4*v4+c*u3*v*(8+(-15+c2)*v2)+u4*(12+(-13+3*c2)*v2) &
                    -u2*(-1+2*v2)*(12+(-2+3*c2)*v2)+cuv*(4-15*v2+2*v4)))/rho11
 fuv(jind(5,1,5)) = (5*v4*(1-4*u4-6*c*u3*v+v2+cuv*(-3+2*v2)-3*u2*(1+(-1+c2)*v2)))/rho9

  fu(jind(5,0,6)) = (-5*v6*(u+c*v))/rho7
fuuu(jind(5,0,6)) = (-105*v6*(u+c*v)*(-1+2*u2+4*cuv+(-1+3*c2)*v2))/rho11
  fv(jind(5,0,6)) = (v5*(6+6*u2+7*cuv+v2))/rho7
fvvv(jind(5,0,6)) = (15*v3*(8+8*u6+18*c*u5*v-12*v2+v4+cuv*(18-17*v2)+6*u4*(4+(-2+3*c2)*v2)+c*u3*v*(36+(-17+7*c2)*v2) &
                    +u2*(24+6*(-4+3*c2)*v2+(1-7*c2)*v4)))/rho11
 fuv(jind(5,0,6)) = (-5*v5*(6*u3+7*c*v+12*c*u2*v+u*(6+(-1+7*c2)*v2)))/rho9

  fu(jind(7,3,0)) = (u2*(3-4*u2-cuv+3*v2))/rho9
fuuu(jind(7,3,0)) = (-3*(40*u6+30*c*u5*v+30*cuv*(1+v2)**2-2*(1+v2)**3+6*u4*(-22+(-22+3*c2)*v2)+c*u3*v*(-129+(-129+5*c2)*v2) &
                    -3*u2*(1+v2)*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,3,0)) = (-7*u3*(c*u+v))/rho9
fvvv(jind(7,3,0)) = (-63*u3*(c*u+v)*(-3+(-3+11*c2)*u2+16*cuv+8*v2))/rho13
 fuv(jind(7,3,0)) = (7*u2*(c2*u2*v+c*u*(-4+5*u2-v2)-3*v*(1-2*u2+v2)))/rho11

  fu(jind(7,2,1)) = (u*v*(2-5*u2-3*cuv+2*v2))/rho9
fuuu(jind(7,2,1)) = (-21*v*(10*u5+18*c*u4*v+2*c*v*(1+v2)**2+5*c*u2*v*(-5+(-5+c2)*v2)-2*u*(1+v2)*(-2+(-2+5*c2)*v2) &
                    +u3*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,2,1)) = (u2*(1+u2-5*cuv-6*v2))/rho9
fvvv(jind(7,2,1)) = (21*u2*(-1+(-1+9*c2)*u4+c*(23-15*c2)*u3*v+16*v2-16*v4+cuv*(23-40*v2)+u2*(-2+16*v2+c2*(9-40*v2))))/rho13
 fuv(jind(7,2,1)) = (u*(-5*u4+29*c*u3*v+5*cuv*(-4+3*v2)+u2*(-3+(46+15*c2)*v2)-2*(-1+5*v2+6*v4)))/rho11


  fu(jind(7,1,2)) = (v2*(1-6*u2-5*cuv+v2))/rho9
fuuu(jind(7,1,2)) = (21*v2*(-16*u4-40*c*u3*v+cuv*(23+(23-15*c2)*v2)-8*u2*(-2+(-2+5*c2)*v2)+(1+v2)*(-1+(-1+9*c2)*v2)))/rho13
  fv(jind(7,1,2)) = (u*v*(2+2*u2-3*cuv-5*v2))/rho9
fvvv(jind(7,1,2)) = (-21*u*(5*c**3*u3*v2-5*c2*u2*v*(2+2*u2-3*v2)+c*u*(2+2*u4-25*v2+18*v4+u2*(4-25*v2)) &
                    +v*(4+4*u4-19*v2+10*v4+u2*(8-19*v2))))/rho13
 fuv(jind(7,1,2)) = (v*(2-12*u4+15*c*u3*v-3*v2-5*v4+cuv*(-20+29*v2)+u2*(-10+(46+15*c2)*v2)))/rho11

  fu(jind(7,0,3)) = (-7*v3*(u+c*v))/rho9
fuuu(jind(7,0,3)) = (-63*v3*(u+c*v)*(-3+8*u2+16*cuv+(-3+11*c2)*v2))/rho13
  fv(jind(7,0,3)) = (v2*(3+3*u2-cuv-4*v2))/rho9
fvvv(jind(7,0,3)) = (-3*(-2-2*u6+30*c*u5*v+57*v2-132*v4+40*v6+u4*(-6+(57-45*c2)*v2)+c*u3*v*(60+(-129+5*c2)*v2) &
                    +3*cuv*(10-43*v2+10*v4)+3*u2*(-2+(38-15*c2)*v2+(-44+6*c2)*v4)))/rho13
 fuv(jind(7,0,3)) = (7*v2*(-3*u3-c*u2*v+c*v*(-4+5*v2)+u*(-3+(6+c2)*v2)))/rho11

  fu(jind(7,4,0)) = (u3*(4-3*u2+cuv+4*v2))/rho9
fuuu(jind(7,4,0)) = (3*u*(-20*u6+20*c*u5*v-36*cuv*(1+v2)**2+8*(1+v2)**3+5*u4*(23+(23+c2)*v2)+c*u3*v*(47+(47+c2)*v2) &
                    +4*u2*(1+v2)*(-22+(-22+3*c2)*v2)))/rho13
  fv(jind(7,4,0)) = (-7*u4*(c*u+v))/rho9
fvvv(jind(7,4,0)) = (-63*u4*(c*u+v)*(-3+(-3+11*c2)*u2+16*cuv+8*v2))/rho13
 fuv(jind(7,4,0)) = (7*u3*(5*u2*v-c2*u2*v+c*u*(-5+4*u2-4*v2)-4*v*(1+v2)))/rho11

  fu(jind(7,3,1)) = (u2*v*(3-4*u2-cuv+3*v2))/rho9
fuuu(jind(7,3,1)) = (3*v*(-40*u6-30*c*u5*v-30*cuv*(1+v2)**2+2*(1+v2)**3+c*u3*v*(129+(129-5*c2)*v2)-6*u4*(-22+(-22+3*c2)*v2) &
                    +3*u2*(1+v2)*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,3,1)) = (u3*(1+u2-5*cuv-6*v2))/rho9
fvvv(jind(7,3,1)) = (21*u3*(-1+(-1+9*c2)*u4+c*(23-15*c2)*u3*v+16*v2-16*v4+cuv*(23-40*v2)+u2*(-2+16*v2+c2*(9-40*v2))))/rho13
 fuv(jind(7,3,1)) = -((u2*(4*u4-26*c*u3*v+cuv*(23+2*v2)+u2*(1-(41+5*c2)*v2)+3*(-1+5*v2+6*v4)))/rho11)

  fu(jind(7,2,2)) = (u*v2*(2-5*u2-3*cuv+2*v2))/rho9
fuuu(jind(7,2,2)) = (-21*v2*(10*u5+18*c*u4*v+2*c*v*(1+v2)**2+5*c*u2*v*(-5+(-5+c2)*v2)-2*u*(1+v2)*(-2+(-2+5*c2)*v2) &
                    +u3*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,2,2)) = (u2*v*(2+2*u2-3*cuv-5*v2))/rho9
fvvv(jind(7,2,2)) = (-21*u2*(5*c**3*u3*v2-5*c2*u2*v*(2+2*u2-3*v2)+c*u*(2+2*u4-25*v2+18*v4+u2*(4-25*v2))+v*(4+4*u4-19*v2 &
                    +10*v4+u2*(8-19*v2))))/rho13
 fuv(jind(7,2,2)) = (u*v*(4-10*u4+16*c*u3*v-6*v2-10*v4+cuv*(-19+16*v2)+u2*(-6+(43+9*c2)*v2)))/rho11

  fu(jind(7,1,3)) = (v3*(1-6*u2-5*cuv+v2))/rho9
fuuu(jind(7,1,3)) = (21*v3*(-16*u4-40*c*u3*v+cuv*(23+(23-15*c2)*v2)-8*u2*(-2+(-2+5*c2)*v2)+(1+v2)*(-1+(-1+9*c2)*v2)))/rho13
  fv(jind(7,1,3)) = (u*v2*(3+3*u2-cuv-4*v2))/rho9
fvvv(jind(7,1,3)) = (3*u*(2+2*u6-30*c*u5*v-57*v2+132*v4-40*v6+c*u3*v*(-60+(129-5*c2)*v2)+u4*(6+(-57+45*c2)*v2) &
                    -3*cuv*(10-43*v2+10*v4)-3*u2*(-2+(38-15*c2)*v2+(-44+6*c2)*v4)))/rho13
 fuv(jind(7,1,3)) = -((v2*(-3+18*u4+2*c*u3*v+v2+4*v4+cuv*(23-26*v2)+u2*(15-(41+5*c2)*v2)))/rho11)

  fu(jind(7,0,4)) = (-7*v4*(u+c*v))/rho9
fuuu(jind(7,0,4)) = (-63*v4*(u+c*v)*(-3+8*u2+16*cuv+(-3+11*c2)*v2))/rho13
  fv(jind(7,0,4)) = (v3*(4+4*u2+cuv-3*v2))/rho9
fvvv(jind(7,0,4)) = (3*v*(8+8*u6-36*c*u5*v-88*v2+115*v4-20*v6+c*u3*v*(-72+(47+c2)*v2)+4*u4*(6+(-22+3*c2)*v2) &
                    +cuv*(-36+47*v2+20*v4)+u2*(24+4*(-44+3*c2)*v2+5*(23+c2)*v4)))/rho13
 fuv(jind(7,0,4)) = (-7*v3*(4*u3+4*c*u2*v+c*v*(5-4*v2)+u*(4+(-5+c2)*v2)))/rho11

  fu(jind(7,5,0)) = (u4*(5-2*u2+3*cuv+5*v2))/rho9
fuuu(jind(7,5,0)) = (-3*u2*(8*u6-36*c*u5*v+20*cuv*(1+v2)**2-20*(1+v2)**3+5*u2*(1+v2)*(23+(23+c2)*v2) &
                    +c*u3*v*(47+(47+c2)*v2)+4*u4*(-22+(-22+3*c2)*v2)))/rho13
  fv(jind(7,5,0)) = (-7*u5*(c*u+v))/rho9
fvvv(jind(7,5,0)) = (-63*u5*(c*u+v)*(-3+(-3+11*c2)*u2+16*cuv+8*v2))/rho13
 fuv(jind(7,5,0)) = (7*u4*(4*u2*v-3*c2*u2*v+c*u*(-6+3*u2-7*v2)-5*v*(1+v2)))/rho11

  fu(jind(7,4,1)) = (u3*v*(4-3*u2+cuv+4*v2))/rho9
fuuu(jind(7,4,1)) = (3*u*v*(-20*u6+20*c*u5*v-36*cuv*(1+v2)**2+8*(1+v2)**3+5*u4*(23+(23+c2)*v2) &
                    +c*u3*v*(47+(47+c2)*v2)+4*u2*(1+v2)*(-22+(-22+3*c2)*v2)))/rho13
  fv(jind(7,4,1)) = (u4*(1+u2-5*cuv-6*v2))/rho9
fvvv(jind(7,4,1)) = (21*u4*(-1+(-1+9*c2)*u4+c*(23-15*c2)*u3*v+16*v2-16*v4+cuv*(23-40*v2)+u2*(-2+16*v2+c2*(9-40*v2))))/rho13
 fuv(jind(7,4,1)) = (u3*(-3*u4+23*c*u3*v-cuv*(26+19*v2)+u2*(1+(36-5*c2)*v2)-4*(-1+5*v2+6*v4)))/rho11

  fu(jind(7,3,2)) = (u2*v2*(3 - 4*u2 - cuv + 3*v2))/(1 + u2 + 2*cuv + v2)**4.5
fuuu(jind(7,3,2)) = (3*v2*(-40*u6-30*c*u5*v-30*cuv*(1+v2)**2+2*(1+v2)**3+c*u3*v*(129+(129-5*c2)*v2) &
                    -6*u4*(-22+(-22+3*c2)*v2)+3*u2*(1+v2)*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,3,2)) = (u3*v*(2+2*u2-3*cuv-5*v2))/rho9
fvvv(jind(7,3,2)) = (-21*u3*(5*c**3*u3*v2-5*c2*u2*v*(2+2*u2-3*v2)+c*u*(2+2*u4-25*v2+18*v4+u2*(4-25*v2)) &
                    +v*(4+4*u4-19*v2+10*v4+u2*(8-19*v2))))/rho13
 fuv(jind(7,3,2)) = (u2*v*(6-8*u4+17*c*u3*v-9*v2-15*v4+3*cuv*(-6+v2)+u2*(-2+(40+3*c2)*v2)))/rho11

fu(jind(7,2,3)) = (u*v3*(2-5*u2-3*cuv+2*v2))/rho9
fuuu(jind(7,2,3)) = (-21*v3*(10*u5+18*c*u4*v+2*c*v*(1+v2)**2+5*c*u2*v*(-5+(-5+c2)*v2)-2*u*(1+v2)*(-2+(-2+5*c2)*v2) &
                    +u3*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,2,3)) = (u2*v2*(3+3*u2-cuv-4*v2))/rho9
fvvv(jind(7,2,3)) = (3*u2*(2+2*u6-30*c*u5*v-57*v2+132*v4-40*v6+c*u3*v*(-60+(129-5*c2)*v2)+u4*(6+(-57+45*c2)*v2) &
                    -3*cuv*(10-43*v2+10*v4)-3*u2*(-2+(38-15*c2)*v2+(-44+6*c2)*v4)))/rho13
 fuv(jind(7,2,3)) = (u*v2*(-15*u4+3*c*u3*v+cuv*(-18+17*v2)+u2*(-9+(40+3*c2)*v2)-2*(-3+v2+4*v4)))/rho11

  fu(jind(7,1,4)) = (v4*(1-6*u2-5*cuv+v2))/rho9
fuuu(jind(7,1,4)) = (21*v4*(-16*u4-40*c*u3*v+cuv*(23+(23-15*c2)*v2)-8*u2*(-2+(-2+5*c2)*v2)+(1+v2)*(-1+(-1+9*c2)*v2)))/rho13
  fv(jind(7,1,4)) = (u*v3*(4+4*u2+cuv-3*v2))/rho9
fvvv(jind(7,1,4)) = (3*u*v*(8+8*u6-36*c*u5*v-88*v2+115*v4-20*v6+c*u3*v*(-72+(47+c2)*v2)+4*u4*(6+(-22+3*c2)*v2) &
                    +cuv*(-36+47*v2+20*v4)+u2*(24+4*(-44+3*c2)*v2+5*(23+c2)*v4)))/rho13
 fuv(jind(7,1,4)) = (v3*(4-24*u4-19*c*u3*v+v2-3*v4+cuv*(-26+23*v2)+u2*(-20+(36-5*c2)*v2)))/rho11

  fu(jind(7,0,5)) = (-7*v5*(u+c*v))/rho9
fuuu(jind(7,0,5)) = (-63*v5*(u+c*v)*(-3+8*u2+16*cuv+(-3+11*c2)*v2))/rho13
  fv(jind(7,0,5)) = (v4*(5+5*u2+3*cuv-2*v2))/rho9
fvvv(jind(7,0,5)) = (-3*v2*(-20-20*u6+20*c*u5*v+115*v2-88*v4+8*v6+5*u4*(-12+(23+c2)*v2)+c*u3*v*(40+(47+c2)*v2) &
                    +cuv*(20+47*v2-36*v4)+u2*(-60+5*(46+c2)*v2+4*(-22+3*c2)*v4)))/rho13
 fuv(jind(7,0,5)) = (-7*v4*(5*u3+7*c*u2*v-3*c*v*(-2+v2)+u*(5+(-4+3*c2)*v2)))/rho11

  fu(jind(7,6,0)) = (u5*(6-u2+5*cuv+6*v2))/rho9
fuuu(jind(7,6,0)) = (3*u3*(-2*u6+30*c*u5*v+30*cuv*(1+v2)**2+40*(1+v2)**3+u4*(57+(57-45*c2)*v2)+6*u2*(1+v2)*(-22+(-22+3*c2)*v2) &
                    +c*u3*v*(-129+(-129+5*c2)*v2)))/rho13
  fv(jind(7,6,0)) = (-7*u6*(c*u+v))/rho9
fvvv(jind(7,6,0)) = (-63*u6*(c*u+v)*(-3+(-3+11*c2)*u2+16*cuv+8*v2))/rho13
 fuv(jind(7,6,0)) = (7*u5*(-5*c2*u2*v+c*u*(-7+2*u2-10*v2)+3*v*(u2-2*(1+v2))))/rho11

  fu(jind(7,5,1)) = (u4*v*(5-2*u2+3*cuv+5*v2))/rho9
fuuu(jind(7,5,1)) = (-3*u2*v*(8*u6-36*c*u5*v+20*cuv*(1+v2)**2-20*(1+v2)**3+5*u2*(1+v2)*(23+(23+c2)*v2)+c*u3*v*(47+(47+c2)*v2) &
                    +4*u4*(-22+(-22+3*c2)*v2)))/rho13
  fv(jind(7,5,1)) = (u5*(1+u2-5*cuv-6*v2))/rho9
fvvv(jind(7,5,1)) = (21*u5*(-1+(-1+9*c2)*u4+c*(23-15*c2)*u3*v+16*v2-16*v4+cuv*(23-40*v2)+u2*(-2+16*v2+c2*(9-40*v2))))/rho13
 fuv(jind(7,5,1)) = -((u4*(2*u4-20*c*u3*v+cuv*(29+36*v2)+u2*(-3+(-31+15*c2)*v2)+5*(-1+5*v2+6*v4)))/rho11)

  fu(jind(7,4,2)) = (u3*v2*(4-3*u2+cuv+4*v2))/rho9
fuuu(jind(7,4,2)) = (3*u*v2*(-20*u6+20*c*u5*v-36*cuv*(1+v2)**2+8*(1+v2)**3+5*u4*(23+(23+c2)*v2)+c*u3*v*(47+(47+c2)*v2) &
                    +4*u2*(1+v2)*(-22+(-22+3*c2)*v2)))/rho13
  fv(jind(7,4,2)) = (u4*v*(2+2*u2-3*cuv-5*v2))/rho9
fvvv(jind(7,4,2)) = (-21*u4*(5*c**3*u3*v2-5*c2*u2*v*(2+2*u2-3*v2)+c*u*(2+2*u4-25*v2+18*v4+u2*(4-25*v2))+v*(4+4*u4-19*v2 &
                    +10*v4+u2*(8-19*v2))))/rho13
 fuv(jind(7,4,2)) = -((u3*v*(6*u4-18*c*u3*v+cuv*(17+10*v2)+u2*(-2+(-37+3*c2)*v2)+4*(-2+3*v2+5*v4)))/rho11)

  fu(jind(7,3,3)) = (u2*v3*(3-4*u2-cuv+3*v2))/rho9
fuuu(jind(7,3,3)) = (3*v3*(-40*u6-30*c*u5*v-30*cuv*(1+v2)**2+2*(1+v2)**3+c*u3*v*(129+(129-5*c2)*v2)-6*u4*(-22+(-22+3*c2)*v2) &
                    +3*u2*(1+v2)*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,3,3)) = (u3*v2*(3+3*u2-cuv-4*v2))/rho9
fvvv(jind(7,3,3)) = (3*u3*(2+2*u6-30*c*u5*v-57*v2+132*v4-40*v6+c*u3*v*(-60+(129-5*c2)*v2)+u4*(6+(-57+45*c2)*v2) &
                    -3*cuv*(10-43*v2+10*v4)-3*u2*(-2+(38-15*c2)*v2+(-44+6*c2)*v4)))/rho13
 fuv(jind(7,3,3)) = (u2*v2*(-12*u4+8*c*u3*v+cuv*(-13+8*v2)+u2*(-3+(39+c2)*v2)-3*(-3+v2+4*v4)))/rho11

  fu(jind(7,2,4)) = (u*v4*(2-5*u2-3*cuv+2*v2))/rho9
fuuu(jind(7,2,4)) = (-21*v4*(10*u5+18*c*u4*v+2*c*v*(1+v2)**2+5*c*u2*v*(-5+(-5+c2)*v2)-2*u*(1+v2)*(-2+(-2+5*c2)*v2) &
                    +u3*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,2,4)) = (u2*v3*(4+4*u2+cuv-3*v2))/rho9
fvvv(jind(7,2,4)) = (3*u2*v*(8+8*u6-36*c*u5*v-88*v2+115*v4-20*v6+c*u3*v*(-72+(47+c2)*v2)+4*u4*(6+(-22+3*c2)*v2) &
                    +cuv*(-36+47*v2+20*v4)+u2*(24+4*(-44+3*c2)*v2+5*(23+c2)*v4)))/rho13
 fuv(jind(7,2,4)) = -((u*v3*(-8+20*u4+10*c*u3*v-2*v2+6*v4+cuv*(17-18*v2)+u2*(12+(-37+3*c2)*v2)))/rho11)

  fu(jind(7,1,5)) = (v5*(1-6*u2-5*cuv+v2))/rho9
fuuu(jind(7,1,5)) = (21*v5*(-16*u4-40*c*u3*v+cuv*(23+(23-15*c2)*v2)-8*u2*(-2+(-2+5*c2)*v2)+(1+v2)*(-1+(-1+9*c2)*v2)))/rho13
  fv(jind(7,1,5)) = (u*v4*(5+5*u2+3*cuv-2*v2))/rho9
fvvv(jind(7,1,5)) = (-3*u*v2*(-20-20*u6+20*c*u5*v+115*v2-88*v4+8*v6+5*u4*(-12+(23+c2)*v2)+c*u3*v*(40+(47+c2)*v2) &
                    +cuv*(20+47*v2-36*v4)+u2*(-60+5*(46+c2)*v2+4*(-22+3*c2)*v4)))/rho13
 fuv(jind(7,1,5)) = -((v4*(-5+30*u4+36*c*u3*v-3*v2+2*v4+cuv*(29-20*v2)+u2*(25+(-31+15*c2)*v2)))/rho11)

  fu(jind(7,0,6)) = (-7*v6*(u+c*v))/rho9
fuuu(jind(7,0,6)) = (-63*v6*(u+c*v)*(-3+8*u2+16*cuv+(-3+11*c2)*v2))/rho13
  fv(jind(7,0,6)) = (v5*(6+6*u2+5*cuv-v2))/rho9
fvvv(jind(7,0,6)) = (3*v3*(40+40*u6+30*c*u5*v-132*v2+57*v4-2*v6+6*u4*(20+(-22+3*c2)*v2)+c*u3*v*(60+(-129+5*c2)*v2) &
                    +3*cuv*(10-43*v2+10*v4)-3*u2*(-40+(88-6*c2)*v2+(-19+15*c2)*v4)))/rho13
 fuv(jind(7,0,6)) = (-7*v5*(6*u3+10*c*u2*v+c*v*(7-2*v2)+u*(6+(-3+5*c2)*v2)))/rho11

  fu(jind(7,7,0)) = (7*u6*(1+cuv+v2))/rho9
fuuu(jind(7,7,0)) = (21*u4*(2*c*u5*v+18*cuv*(1+v2)**2+10*(1+v2)**3+u4*(4+(4-10*c2)*v2)+5*c*u3*v*(-5+(-5+c2)*v2) &
                    +u2*(1+v2)*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,7,0)) = (-7*u7*(c*u+v))/rho9
fvvv(jind(7,7,0)) = (-63*u7*(c*u+v)*(-3+(-3+11*c2)*u2+16*cuv+8*v2))/rho13
 fuv(jind(7,7,0)) = (7*u6*(2*u2*v-7*c2*u2*v+c*u*(-8+u2-13*v2)-7*v*(1+v2)))/rho11

  fu(jind(7,6,1)) = (u5*v*(6 - u2 + 5*cuv + 6*v2))/(1 + u2 + 2*cuv + v2)**4.5
fuuu(jind(7,6,1)) = (3*u3*v*(-2*u6 + 30*c*u5*v + 30*cuv*(1 + v2)**2 + 40*(1 + v2)**3 + u4*(57 + (57 - 45*c2)*v2)  &
                    + 6*u2*(1 + v2)*(-22 + (-22 + 3*c2)*v2) + c*u3*v*(-129 + (-129 + 5*c2)*v2)))/rho13
  fv(jind(7,6,1)) = (u6*(1+u2-5*cuv-6*v2))/rho9
fvvv(jind(7,6,1)) = (21*u6*(-1+(-1+9*c2)*u4+c*(23-15*c2)*u3*v+16*v2-16*v4+cuv*(23-40*v2)+u2*(-2+16*v2+c2*(9-40*v2))))/rho13
 fuv(jind(7,6,1)) = -((u5*(u4-17*c*u3*v+cuv*(32+53*v2)+u2*(-5+(-26+25*c2)*v2)+6*(-1+5*v2+6*v4)))/rho11)

  fu(jind(7,5,2)) = (u4*v2*(5-2*u2+3*cuv+5*v2))/rho9
fuuu(jind(7,5,2)) = (-3*u2*v2*(8*u6-36*c*u5*v+20*cuv*(1+v2)**2-20*(1+v2)**3+5*u2*(1+v2)*(23+(23+c2)*v2)+c*u3*v*(47+(47+c2)*v2) &
                    +4*u4*(-22+(-22+3*c2)*v2)))/rho13
  fv(jind(7,5,2)) = (u5*v*(2+2*u2-3*cuv-5*v2))/rho9
fvvv(jind(7,5,2)) = (-21*u5*(5*c**3*u3*v2-5*c2*u2*v*(2+2*u2-3*v2)+c*u*(2+2*u4-25*v2+18*v4+u2*(4-25*v2))+v*(4+4*u4-19*v2+10*v4 &
                    +u2*(8-19*v2))))/rho13
 fuv(jind(7,5,2)) = -((u4*v*(4*u4-19*c*u3*v+cuv*(16+23*v2)+u2*(-6+(-34+9*c2)*v2)+5*(-2+3*v2+5*v4)))/rho11)

  fu(jind(7,4,3)) = (u3*v3*(4-3*u2+cuv+4*v2))/rho9
fuuu(jind(7,4,3)) = (3*u*v3*(-20*u6+20*c*u5*v-36*cuv*(1+v2)**2+8*(1+v2)**3+5*u4*(23+(23+c2)*v2)+c*u3*v*(47+(47+c2)*v2) &
                    +4*u2*(1+v2)*(-22+(-22+3*c2)*v2)))/rho13
  fv(jind(7,4,3)) = (u4*v2*(3+3*u2-cuv-4*v2))/rho9
fvvv(jind(7,4,3)) = (3*u4*(2+2*u6-30*c*u5*v-57*v2+132*v4-40*v6+c*u3*v*(-60+(129-5*c2)*v2)+u4*(6+(-57+45*c2)*v2) &
                    -3*cuv*(10-43*v2+10*v4)-3*u2*(-2+(38-15*c2)*v2+(-44+6*c2)*v4)))/rho13
 fuv(jind(7,4,3)) = -((u3*v2*(9*u4-13*c*u3*v+cuv*(8+v2)+u2*(-3+(-38+c2)*v2)+4*(-3+v2+4*v4)))/rho11)

  fu(jind(7,3,4)) = (u2*v4*(3-4*u2-cuv+3*v2))/rho9
fuuu(jind(7,3,4)) = (3*v4*(-40*u6-30*c*u5*v-30*cuv*(1+v2)**2+2*(1+v2)**3+c*u3*v*(129+(129-5*c2)*v2)-6*u4*(-22+(-22+3*c2)*v2) &
                    +3*u2*(1+v2)*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,3,4)) = (u3*v3*(4+4*u2+cuv-3*v2))/rho9
fvvv(jind(7,3,4)) = (3*u3*v*(8+8*u6-36*c*u5*v-88*v2+115*v4-20*v6+c*u3*v*(-72+(47+c2)*v2)+4*u4*(6+(-22+3*c2)*v2) &
                    +cuv*(-36+47*v2+20*v4)+u2*(24+4*(-44+3*c2)*v2+5*(23+c2)*v4)))/rho13
 fuv(jind(7,3,4)) = -((u2*v3*(-12+16*u4+c*u3*v-3*v2+9*v4+cuv*(8-13*v2)+u2*(4+(-38+c2)*v2)))/rho11)

  fu(jind(7,2,5)) = (u*v5*(2-5*u2-3*cuv+2*v2))/rho9
fuuu(jind(7,2,5)) = (-21*v5*(10*u5+18*c*u4*v+2*c*v*(1+v2)**2+5*c*u2*v*(-5+(-5+c2)*v2)-2*u*(1+v2)*(-2+(-2+5*c2)*v2) &
                    +u3*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,2,5)) = (u2*v4*(5+5*u2+3*cuv-2*v2))/rho9
fvvv(jind(7,2,5)) = (-3*u2*v2*(-20-20*u6+20*c*u5*v+115*v2-88*v4+8*v6+5*u4*(-12+(23+c2)*v2)+c*u3*v*(40+(47+c2)*v2) &
                    +cuv*(20+47*v2-36*v4)+u2*(-60+5*(46+c2)*v2+4*(-22+3*c2)*v4)))/rho13
 fuv(jind(7,2,5)) = -((u*v4*(-10+25*u4+23*c*u3*v-6*v2+4*v4+cuv*(16-19*v2)+u2*(15+(-34+9*c2)*v2)))/rho11)

  fu(jind(7,1,6)) = (v6*(1-6*u2-5*cuv+v2))/rho9
fuuu(jind(7,1,6)) = (21*v6*(-16*u4-40*c*u3*v+cuv*(23+(23-15*c2)*v2)-8*u2*(-2+(-2+5*c2)*v2)+(1+v2)*(-1+(-1+9*c2)*v2)))/rho13
  fv(jind(7,1,6)) = (u*v5*(6+6*u2+5*cuv-v2))/rho9
fvvv(jind(7,1,6)) = (3*u*v3*(40+40*u6+30*c*u5*v-132*v2+57*v4-2*v6+6*u4*(20+(-22+3*c2)*v2)+c*u3*v*(60+(-129+5*c2)*v2) &
                    +3*cuv*(10-43*v2+10*v4)-3*u2*(-40+(88-6*c2)*v2+(-19+15*c2)*v4)))/rho13
 fuv(jind(7,1,6)) = -((v5*(-6+36*u4+53*c*u3*v-5*v2+v4+cuv*(32-17*v2)+u2*(30+(-26+25*c2)*v2)))/rho11)

  fu(jind(7,0,7)) = (-7*v7*(u+c*v))/rho9
fuuu(jind(7,0,7)) = (-63*v7*(u+c*v)*(-3+8*u2+16*cuv+(-3+11*c2)*v2))/rho13
  fv(jind(7,0,7)) = (7*v6*(1+u2+cuv))/rho9
fvvv(jind(7,0,7)) = (21*v4*(10+10*u6+18*c*u5*v-19*v2+4*v4+c*u3*v*(36+5*(-5+c2)*v2)+u4*(30+(-19+15*c2)*v2) &
                    +cuv*(18-25*v2+2*v4)+u2*(30+(-38+15*c2)*v2+(4-10*c2)*v4)))/rho13
 fuv(jind(7,0,7)) = (-7*v6*(7*u3+13*c*u2*v-c*v*(-8+v2)+u*(7+(-2+7*c2)*v2)))/rho11

  fu(jind(7,8,0)) = (u7*(8 + u2 + 9*cuv + 8*v2))/rho9 
fuuu(jind(7,8,0)) = (21*u5*(40*cuv*(1 + v2)**2 + 16*(1 + v2)**3  &
                    + u4*(1 + (1 - 9*c2)*v2) + 8*u2*(1 + v2)*(-2 + (-2 + 5*c2)*v2) + c*u3*v*(-23 + (-23 + 15*c2)*v2)))/rho13
  fv(jind(7,8,0)) = (-7*u8*(c*u+v))/rho9 
fvvv(jind(7,8,0)) = (-63*u8*(c*u+v)*(-3+(-3+11*c2)*u2+16*cuv+8*v2))/rho13
 fuv(jind(7,8,0)) = (-7*u7*(9*c2*u2*v+v*(8-u2+8*v2)+c*u*(9+16*v2)))/rho11

  fu(jind(7,7,1)) = (7*u6*v*(1+cuv+v2))/rho9
fuuu(jind(7,7,1)) = (21*u4*v*(2*c*u5*v+18*cuv*(1+v2)**2+10*(1+v2)**3+u4*(4+(4-10*c2)*v2)+5*c*u3*v*(-5+(-5+c2)*v2) &
                    +u2*(1+v2)*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,7,1)) = (u7*(1+u2-5*cuv-6*v2))/rho9
fvvv(jind(7,7,1)) = (21*u7*(-1+(-1+9*c2)*u4+c*(23-15*c2)*u3*v+16*v2-16*v4+cuv*(23-40*v2)+u2*(-2+16*v2+c2*(9-40*v2))))/rho13
 fuv(jind(7,7,1)) = (7*u6*(1+2*c*u3*v-5*v2-6*v4-5*cuv*(1+2*v2)+u2*(1+(3-5*c2)*v2)))/rho11

  fu(jind(7,6,2)) = (u5*v2*(6-u2+5*cuv+6*v2))/rho9
fuuu(jind(7,6,2)) = (3*u3*v2*(-2*u6+30*c*u5*v+30*cuv*(1+v2)**2+40*(1+v2)**3+u4*(57+(57-45*c2)*v2) &
                    +6*u2*(1+v2)*(-22+(-22+3*c2)*v2)+c*u3*v*(-129+(-129+5*c2)*v2)))/rho13
  fv(jind(7,6,2)) = (u6*v*(2+2*u2-3*cuv-5*v2))/rho9
fvvv(jind(7,6,2)) = (-21*u6*(5*c**3*u3*v2-5*c2*u2*v*(2+2*u2-3*v2)+c*u*(2+2*u4-25*v2+18*v4+u2*(4-25*v2)) &
                    +v*(4+4*u4-19*v2+10*v4+u2*(8-19*v2))))/rho13
 fuv(jind(7,6,2)) = -((u5*v*(2*u4-20*c*u3*v+3*cuv*(5+12*v2)+u2*(-10+(-31+15*c2)*v2)+6*(-2+3*v2+5*v4)))/rho11)

  fu(jind(7,5,3)) = (u4*v3*(5-2*u2+3*cuv+5*v2))/rho9
fuuu(jind(7,5,3)) = (-3*u2*v3*(8*u6-36*c*u5*v+20*cuv*(1+v2)**2-20*(1+v2)**3+5*u2*(1+v2)*(23+(23+c2)*v2)+c*u3*v*(47+(47+c2)*v2) &
                    +4*u4*(-22+(-22+3*c2)*v2)))/rho13
  fv(jind(7,5,3)) = (u5*v2*(3+3*u2-cuv-4*v2))/rho9
fvvv(jind(7,5,3)) = (3*u5*(2+2*u6-30*c*u5*v-57*v2+132*v4-40*v6+c*u3*v*(-60+(129-5*c2)*v2)+u4*(6+(-57+45*c2)*v2)-3*cuv*(10-43*v2 &
                    +10*v4)-3*u2*(-2+(38-15*c2)*v2+(-44+6*c2)*v4)))/rho13
 fuv(jind(7,5,3)) = -((u4*v2*(6*u4-18*c*u3*v+cuv*(3+10*v2)+u2*(-9+(-37+3*c2)*v2)+5*(-3+v2+4*v4)))/rho11)

  fu(jind(7,4,4)) = (u3*v4*(4-3*u2+cuv+4*v2))/rho9
fuuu(jind(7,4,4)) = (3*u*v4*(-20*u6+20*c*u5*v-36*cuv*(1+v2)**2+8*(1+v2)**3+5*u4*(23+(23+c2)*v2)+c*u3*v*(47+(47+c2)*v2) &
                    +4*u2*(1+v2)*(-22+(-22+3*c2)*v2)))/rho13
  fv(jind(7,4,4)) = (u4*v3*(4+4*u2+cuv-3*v2))/rho9
fvvv(jind(7,4,4)) = (3*u4*v*(8+8*u6-36*c*u5*v-88*v2+115*v4-20*v6+c*u3*v*(-72+(47+c2)*v2)+4*u4*(6+(-22+3*c2)*v2) &
                    +cuv*(-36+47*v2+20*v4)+u2*(24+4*(-44+3*c2)*v2+5*(23+c2)*v4)))/rho13
 fuv(jind(7,4,4)) = (u3*v3*(-12*u4+8*c*u3*v+cuv*(1+8*v2)+u2*(4+(39+c2)*v2)+4*(4+v2-3*v4)))/rho11

  fu(jind(7,3,5)) = (u2*v5*(3-4*u2-cuv+3*v2))/rho9
fuuu(jind(7,3,5)) = (3*v5*(-40*u6-30*c*u5*v-30*cuv*(1+v2)**2+2*(1+v2)**3+c*u3*v*(129+(129-5*c2)*v2)-6*u4*(-22+(-22+3*c2)*v2) &
                    +3*u2*(1+v2)*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,3,5)) = (u3*v4*(5+5*u2+3*cuv-2*v2))/rho9
fvvv(jind(7,3,5)) = (-3*u3*v2*(-20-20*u6+20*c*u5*v+115*v2-88*v4+8*v6+5*u4*(-12+(23+c2)*v2)+c*u3*v*(40+(47+c2)*v2) &
                    +cuv*(20+47*v2-36*v4)+u2*(-60+5*(46+c2)*v2+4*(-22+3*c2)*v4)))/rho13
 fuv(jind(7,3,5)) = -((u2*v4*(-15+20*u4+10*c*u3*v-9*v2+6*v4+3*cuv*(1-6*v2)+u2*(5+(-37+3*c2)*v2)))/rho11)

  fu(jind(7,2,6)) = (u*v6*(2-5*u2-3*cuv+2*v2))/rho9
fuuu(jind(7,2,6)) = (-21*v6*(10*u5+18*c*u4*v+2*c*v*(1+v2)**2+5*c*u2*v*(-5+(-5+c2)*v2)-2*u*(1+v2)*(-2+(-2+5*c2)*v2) &
                    +u3*(-19+(-19+15*c2)*v2)))/rho13
  fv(jind(7,2,6)) = (u2*v5*(6+6*u2+5*cuv-v2))/rho9
fvvv(jind(7,2,6)) = (3*u2*v3*(40+40*u6+30*c*u5*v-132*v2+57*v4-2*v6+6*u4*(20+(-22+3*c2)*v2)+c*u3*v*(60+(-129+5*c2)*v2) &
                    +3*cuv*(10-43*v2+10*v4)-3*u2*(-40+(88-6*c2)*v2+(-19+15*c2)*v4)))/rho13
 fuv(jind(7,2,6)) = -((u*v5*(30*u4+36*c*u3*v-5*cuv*(-3+4*v2)+u2*(18+(-31+15*c2)*v2)+2*(-6-5*v2+v4)))/rho11)

  fu(jind(7,1,7)) = (v7*(1-6*u2-5*cuv+v2))/rho9
fuuu(jind(7,1,7)) = (21*v7*(-16*u4-40*c*u3*v+cuv*(23+(23-15*c2)*v2)-8*u2*(-2+(-2+5*c2)*v2)+(1+v2)*(-1+(-1+9*c2)*v2)))/rho13
  fv(jind(7,1,7)) = (7*u*v6*(1+u2+cuv))/rho9
fvvv(jind(7,1,7)) = (21*u*v4*(10+10*u6+18*c*u5*v-19*v2+4*v4+c*u3*v*(36+5*(-5+c2)*v2)+u4*(30+(-19+15*c2)*v2) &
                    +cuv*(18-25*v2+2*v4)+u2*(30+(-38+15*c2)*v2+(4-10*c2)*v4)))/rho13
 fuv(jind(7,1,7)) = (7*v6*(1-6*u4-10*c*u3*v+v2+cuv*(-5+2*v2)+u2*(-5+(3-5*c2)*v2)))/rho11

  fu(jind(7,0,8)) = (-7*v8*(u+c*v))/rho9
fuuu(jind(7,0,8)) = (-63*v8*(u+c*v)*(-3+8*u2+16*cuv+(-3+11*c2)*v2))/rho13
  fv(jind(7,0,8)) = (v7*(8+8*u2+9*cuv+v2))/rho9
fvvv(jind(7,0,8)) = (21*v5*(16+16*u6+40*c*u5*v-16*v2+v4+cuv*(40-23*v2)+8*u4*(6+(-2+5*c2)*v2)+c*u3*v*(80+(-23+15*c2)*v2) &
                     +u2*(48+8*(-4+5*c2)*v2+(1-9*c2)*v4)))/rho13
 fuv(jind(7,0,8)) = (-7*v7*(8*u3+9*c*v+16*c*u2*v+u*(8+(-1+9*c2)*v2)))/rho11

  fu(jind(9,6,0)) = (3*u5*(2-u2+cuv+2*v2))/rho11
fuuu(jind(9,6,0)) = (-3*u3*(20*u6-60*c*u5*v+30*cuv*(1+v2)**2-40*(1+v2)**3+15*u4*(-11+(-11+c2)*v2) &
                     +6*u2*(1+v2)*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,6,0)) = (-9*u6*(c*u+v))/rho11
fvvv(jind(9,6,0)) = (-99*u6*(c*u+v)*(-3+(-3+13*c2)*u2+20*cuv+10*v2))/rho15
 fuv(jind(9,6,0)) = (9*u5*(5*u2*v-3*c2*u2*v+c*u*(-7+4*u2-8*v2)-6*v*(1+v2)))/rho13

  fu(jind(9,5,1)) = (u4*v*(5-4*u2+cuv+5*v2))/rho11
fuuu(jind(9,5,1)) = (3*u2*v*(-40*u6+30*c*u5*v-60*cuv*(1+v2)**2+20*(1+v2)**3+15*u2*(1+v2)*(-11+(-11+c2)*v2) &
                     +6*u4*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,5,1)) = (u5*(1+u2-7*cuv-8*v2))/rho11
fvvv(jind(9,5,1)) = (9*u5*(-3+(-3+33*c2)*u4+c*(87-77*c2)*u3*v+60*v2-80*v4+3*cuv*(29-70*v2)+u2*(-6+60*v2+c2*(33-210*v2))))/rho15
 fuv(jind(9,5,1)) = (u4*(-4*u4+38*c*u3*v-cuv*(43+34*v2)+u2*(1+(55-7*c2)*v2)-5*(-1+7*v2+8*v4)))/rho13

  fu(jind(9,4,2)) = (u3*v2*(4-5*u2-cuv+4*v2))/rho11
fuuu(jind(9,4,2)) = (-3*u*v2*(70*u6+42*c*u5*v+60*cuv*(1+v2)**2-8*(1+v2)**3+5*c*u3*v*(-39+(-39+c2)*v2) &
                     +21*u4*(-11+(-11+c2)*v2)-60*u2*(1+v2)*(-2+(-2+c2)*v2)))/rho15
  fv(jind(9,4,2)) = (u4*v*(2+2*u2-5*cuv-7*v2))/rho11
fvvv(jind(9,4,2)) = (-9*u4*(35*c**3*u3*v2-21*c2*u2*v*(2+2*u2-5*v2)+v*(12+12*u4-75*v2+56*v4+u2*(24-75*v2)) &
                     +3*c*u*(2+2*u4-35*v2+40*v4+u2*(4-35*v2))))/rho15
 fuv(jind(9,4,2)) = -((u3*v*(10*u4-32*c*u3*v+cuv*(31+4*v2)+u2*(2-(61+5*c2)*v2)+4*(-2+5*v2+7*v4)))/rho13)

  fu(jind(9,3,3)) = (3*u2*v3*(1-2*u2-cuv+v2))/rho11
fuuu(jind(9,3,3)) = (3*v3*(-112*u6-168*c*u5*v-42*cuv*(1+v2)**2+2*(1+v2)**3-120*u4*(-2+(-2+c2)*v2) &
                    -5*c*u3*v*(-57+(-57+7*c2)*v2)+15*u2*(1+v2)*(-5+(-5+7*c2)*v2)))/rho15

  fv(jind(9,3,3)) = (3*u3*v2*(1+u2-cuv-2*v2))/rho11
fvvv(jind(9,3,3)) = (3*u3*(2+2*u6-42*c*u5*v-75*v2+240*v4-112*v6+c*u3*v*(-84-5*(-57+7*c2)*v2)+3*u4*(2+5*(-5+7*c2)*v2) &
                     -3*cuv*(14-95*v2+56*v4)-3*u2*(-2+(50-35*c2)*v2+40*(-2+c2)*v4)))/rho15
 fuv(jind(9,3,3)) = (-9*u2*v2*(-1+u2+2*u4-2*c*u3*v+v2-(7+c2)*u2*v2+2*v4+cuv*(3-2*v2)))/rho13

  fu(jind(9,2,4)) = (u*v4*(2-7*u2-5*cuv+2*v2))/rho11
fuuu(jind(9,2,4)) = (-9*v4*(56*u5+120*c*u4*v+6*c*v*(1+v2)**2+35*c*u2*v*(-3+(-3+c2)*v2)+15*u3*(-5+(-5+7*c2)*v2) &
                     -6*u*(1+v2)*(-2+(-2+7*c2)*v2)))/rho15
  fv(jind(9,2,4)) = (u2*v3*(4+4*u2-cuv-5*v2))/rho11
fvvv(jind(9,2,4)) = (3*u2*v*(8+8*u6-60*c*u5*v-120*v2+231*v4-70*v6-5*c*u3*v*(24+(-39+c2)*v2)+12*u4*(2+5*(-2+c2)*v2) &
                     -3*cuv*(20-65*v2+14*v4)-3*u2*(-8-20*(-4+c2)*v2+7*(-11+c2)*v4)))/rho15
 fuv(jind(9,2,4)) = -((u*v3*(28*u4+4*c*u3*v+cuv*(31-32*v2)+u2*(20-(61+5*c2)*v2)+2*(-4+v2+5*v4)))/rho13)


  fu(jind(9,1,5)) = (v5*(1-8*u2-7*cuv+v2))/rho11
fuuu(jind(9,1,5)) = (9*v5*(-80*u4-210*c*u3*v+cuv*(87+(87-77*c2)*v2)-30*u2*(-2+(-2+7*c2)*v2)+3*(1+v2)*(-1+(-1+11*c2)*v2)))/rho15
  fv(jind(9,1,5)) = (u*v4*(5+5*u2+cuv-4*v2))/rho11
fvvv(jind(9,1,5)) = (3*u*v2*(20+20*u6-60*c*u5*v-165*v2+204*v4-40*v6+15*u4*(4+(-11+c2)*v2)+c*u3*v*(-120+(69+c2)*v2) &
                     +3*cuv*(-20+23*v2+10*v4)+3*u2*(20+5*(-22+c2)*v2+2*(34+c2)*v4)))/rho15
 fuv(jind(9,1,5)) = (v4*(5-40*u4-34*c*u3*v+v2-4*v4+cuv*(-43+38*v2)+u2*(-35+(55-7*c2)*v2)))/rho13

  fu(jind(9,0,6)) = (-9*v6*(u+c*v))/rho11
fuuu(jind(9,0,6)) = (-99*v6*(u+c*v)*(-3+10*u2+20*cuv+(-3+13*c2)*v2))/rho15
  fv(jind(9,0,6)) = (3*v5*(2+2*u2+cuv-v2))/rho11
fvvv(jind(9,0,6)) = (-3*v3*(-40-40*u6+30*c*u5*v+204*v2-165*v4+20*v6+6*u4*(-20+(34+c2)*v2)+c*u3*v*(60+(69+c2)*v2) &
                     +3*cuv*(10+23*v2-20*v4)+3*u2*(-40+2*(68+c2)*v2+5*(-11+c2)*v4)))/rho15
 fuv(jind(9,0,6)) = (-9*v5*(6*u3+8*c*u2*v+c*v*(7-4*v2)+u*(6+(-5+3*c2)*v2)))/rho13

  fu(jind(9,7,0)) = (u6*(7-2*u2+5*cuv+7*v2))/rho11
fuuu(jind(9,7,0)) = (3*u4*(-8*u6+60*c*u5*v+42*cuv*(1+v2)**2+70*(1+v2)**3+5*c*u3*v*(-39+(-39+c2)*v2) &
                     +21*u2*(1+v2)*(-11+(-11+c2)*v2)-60*u4*(-2+(-2+c2)*v2)))/rho15
  fv(jind(9,7,0)) = (-9*u7*(c*u+v))/rho11
fvvv(jind(9,7,0)) = (-99*u7*(c*u+v)*(-3+(-3+13*c2)*u2+20*cuv+10*v2))/rho15
 fuv(jind(9,7,0)) = (9*u6*(4*u2*v-5*c2*u2*v+c*u*(-8+3*u2-11*v2)-7*v*(1+v2)))/rho13

  fu(jind(9,6,1)) = (3*u5*v*(2-u2+cuv+2*v2))/rho11
fuuu(jind(9,6,1)) = (-3*u3*v*(20*u6-60*c*u5*v+30*cuv*(1+v2)**2-40*(1+v2)**3+15*u4*(-11+(-11+c2)*v2) &
                     +6*u2*(1+v2)*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,6,1)) = (u6*(1+u2-7*cuv-8*v2))/rho11
fvvv(jind(9,6,1)) = (9*u6*(-3+(-3+33*c2)*u4+c*(87-77*c2)*u3*v+60*v2-80*v4+3*cuv*(29-70*v2)+u2*(-6+60*v2+c2*(33-210*v2))))/rho15
 fuv(jind(9,6,1)) = (-3*u5*(u4-11*c*u3*v+cuv*(16+19*v2)+u2*(-1+(-16+7*c2)*v2)+2*(-1+7*v2+8*v4)))/rho13

  fu(jind(9,5,2)) = (u4*v2*(5-4*u2+cuv+5*v2))/rho11
fuuu(jind(9,5,2)) = (3*u2*v2*(-40*u6+30*c*u5*v-60*cuv*(1+v2)**2+20*(1+v2)**3+15*u2*(1+v2)*(-11+(-11+c2)*v2) &
                     +6*u4*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,5,2)) = (u5*v*(2+2*u2-5*cuv-7*v2))/rho11
fvvv(jind(9,5,2)) = (-9*u5*(35*c**3*u3*v2-21*c2*u2*v*(2+2*u2-5*v2)+v*(12+12*u4-75*v2+56*v4+u2*(24-75*v2)) &
                     +3*c*u*(2+2*u4-35*v2+40*v4+u2*(4-35*v2))))/rho15
 fuv(jind(9,5,2)) = -((u4*v*(8*u4-31*c*u3*v+cuv*(32+23*v2)+u2*(-2+(-56+5*c2)*v2)+5*(-2+5*v2+7*v4)))/rho13)

  fu(jind(9,4,3)) = (u3*v3*(4-5*u2-cuv+4*v2))/rho11
fuuu(jind(9,4,3)) = (-3*u*v3*(70*u6+42*c*u5*v+60*cuv*(1+v2)**2-8*(1+v2)**3+5*c*u3*v*(-39+(-39+c2)*v2) &
                     +21*u4*(-11+(-11+c2)*v2)-60*u2*(1+v2)*(-2+(-2+c2)*v2)))/rho15
  fv(jind(9,4,3)) = (3*u4*v2*(1+u2-cuv-2*v2))/rho11
fvvv(jind(9,4,3)) = (3*u4*(2+2*u6-42*c*u5*v-75*v2+240*v4-112*v6+c*u3*v*(-84-5*(-57+7*c2)*v2)+3*u4*(2+5*(-5+7*c2)*v2) &
                     -3*cuv*(14-95*v2+56*v4)-3*u2*(-2+(50-35*c2)*v2+40*(-2+c2)*v4)))/rho15
 fuv(jind(9,4,3)) = (3*u3*v2*(-5*u4+7*c*u3*v+cuv*(-8+v2)+u2*(-1+(20+c2)*v2)-4*(-1+v2+2*v4)))/rho13

  fu(jind(9,3,4)) = (3*u2*v4*(1-2*u2-cuv+v2))/rho11
fuuu(jind(9,3,4)) = (3*v4*(-112*u6-168*c*u5*v-42*cuv*(1+v2)**2+2*(1+v2)**3-120*u4*(-2+(-2+c2)*v2) &
                     -5*c*u3*v*(-57+(-57+7*c2)*v2)+15*u2*(1+v2)*(-5+(-5+7*c2)*v2)))/rho15
  fv(jind(9,3,4)) = (u3*v3*(4+4*u2-cuv-5*v2))/rho11
fvvv(jind(9,3,4)) = (3*u3*v*(8+8*u6-60*c*u5*v-120*v2+231*v4-70*v6-5*c*u3*v*(24+(-39+c2)*v2)+12*u4*(2+5*(-2+c2)*v2) &
                     -3*cuv*(20-65*v2+14*v4)-3*u2*(-8-20*(-4+c2)*v2+7*(-11+c2)*v4)))/rho15
 fuv(jind(9,3,4)) = (3*u2*v3*(4-8*u4+c*u3*v-v2-5*v4+cuv*(-8+7*v2)+u2*(-4+(20+c2)*v2)))/rho13

  fu(jind(9,2,5)) = (u*v5*(2-7*u2-5*cuv+2*v2))/rho11
fuuu(jind(9,2,5)) = (-9*v5*(56*u5+120*c*u4*v+6*c*v*(1+v2)**2+35*c*u2*v*(-3+(-3+c2)*v2)+15*u3*(-5+(-5+7*c2)*v2) &
                     -6*u*(1+v2)*(-2+(-2+7*c2)*v2)))/rho15
  fv(jind(9,2,5)) = (u2*v4*(5+5*u2+cuv-4*v2))/rho11
fvvv(jind(9,2,5)) = (3*u2*v2*(20+20*u6-60*c*u5*v-165*v2+204*v4-40*v6+15*u4*(4+(-11+c2)*v2)+c*u3*v*(-120+(69+c2)*v2) &
                     +3*cuv*(-20+23*v2+10*v4)+3*u2*(20+5*(-22+c2)*v2+2*(34+c2)*v4)))/rho15
 fuv(jind(9,2,5)) = -((u*v4*(-10+35*u4+23*c*u3*v-2*v2+8*v4+cuv*(32-31*v2)+u2*(25+(-56+5*c2)*v2)))/rho13)

  fu(jind(9,1,6)) = (v6*(1-8*u2-7*cuv+v2))/rho11
fuuu(jind(9,1,6)) = (9*v6*(-80*u4-210*c*u3*v+cuv*(87+(87-77*c2)*v2)-30*u2*(-2+(-2+7*c2)*v2)+3*(1+v2)*(-1+(-1+11*c2)*v2)))/rho15
  fv(jind(9,1,6)) = (3*u*v5*(2+2*u2+cuv-v2))/rho11
fvvv(jind(9,1,6)) = (-3*u*v3*(-40-40*u6+30*c*u5*v+204*v2-165*v4+20*v6+6*u4*(-20+(34+c2)*v2)+c*u3*v*(60+(69+c2)*v2) &
                     +3*cuv*(10+23*v2-20*v4)+3*u2*(-40+2*(68+c2)*v2+5*(-11+c2)*v4)))/rho15
 fuv(jind(9,1,6)) = (-3*v5*(-2+16*u4+19*c*u3*v-v2+v4+cuv*(16-11*v2)+u2*(14+(-16+7*c2)*v2)))/rho13

  fu(jind(9,0,7)) = (-9*v7*(u+c*v))/rho11
fuuu(jind(9,0,7)) = (-99*v7*(u+c*v)*(-3+10*u2+20*cuv+(-3+13*c2)*v2))/rho15
  fv(jind(9,0,7)) = (v6*(7+7*u2+5*cuv-2*v2))/rho11
fvvv(jind(9,0,7)) = (3*v4*(70+70*u6+42*c*u5*v-231*v2+120*v4-8*v6+c*u3*v*(84+5*(-39+c2)*v2) &
                     +21*u4*(10+(-11+c2)*v2)+3*cuv*(14-65*v2+20*v4)-3*u2*(-70-7*(-22+c2)*v2+20*(-2+c2)*v4)))/rho15
 fuv(jind(9,0,7)) = (-9*v6*(7*u3+11*c*u2*v+c*v*(8-3*v2)+u*(7+(-4+5*c2)*v2)))/rho13

  fu(jind(9,8,0)) = (u7*(8-u2+7*cuv+8*v2))/rho11
fuuu(jind(9,8,0)) = (3*u5*(-2*u6+42*c*u5*v+168*cuv*(1+v2)**2+112*(1+v2)**3+120*u2*(1+v2)*(-2+(-2+c2)*v2) &
                     +5*c*u3*v*(-57+(-57+7*c2)*v2)-15*u4*(-5+(-5+7*c2)*v2)))/rho15
  fv(jind(9,8,0)) = (-9*u8*(c*u+v))/rho11
fvvv(jind(9,8,0)) = (-99*u8*(c*u+v)*(-3+(-3+13*c2)*u2+20*cuv+10*v2))/rho15
 fuv(jind(9,8,0)) = (9*u7*(3*u2*v-7*c2*u2*v+c*u*(-9+2*u2-14*v2)-8*v*(1+v2)))/rho13

  fu(jind(9,7,1)) = (u6*v*(7-2*u2+5*cuv+7*v2))/rho11
fuuu(jind(9,7,1)) = (3*u4*v*(-8*u6+60*c*u5*v+42*cuv*(1+v2)**2+70*(1+v2)**3+5*c*u3*v*(-39+(-39+c2)*v2) &
                     +21*u2*(1+v2)*(-11+(-11+c2)*v2)-60*u4*(-2+(-2+c2)*v2)))/rho15
  fv(jind(9,7,1)) = (u7*(1+u2-7*cuv-8*v2))/rho11
fvvv(jind(9,7,1)) = (9*u7*(-3+(-3+33*c2)*u4+c*(87-77*c2)*u3*v+60*v2-80*v4+3*cuv*(29-70*v2)+u2*(-6+60*v2+c2*(33-210*v2))))/rho15
 fuv(jind(9,7,1)) = -((u6*(2*u4-28*c*u3*v+cuv*(53+80*v2)+u2*(-5+(-41+35*c2)*v2)+7*(-1+7*v2+8*v4)))/rho13)

  fu(jind(9,6,2)) = (3*u5*v2*(2-u2+cuv+2*v2))/rho11
fuuu(jind(9,6,2)) = (-3*u3*v2*(20*u6-60*c*u5*v+30*cuv*(1+v2)**2-40*(1+v2)**3+15*u4*(-11+(-11+c2)*v2) &
                     +6*u2*(1+v2)*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,6,2)) = (u6*v*(2+2*u2-5*cuv-7*v2))/rho11
fvvv(jind(9,6,2)) = (-9*u6*(35*c**3*u3*v2-21*c2*u2*v*(2+2*u2-5*v2)+v*(12+12*u4-75*v2+56*v4+u2*(24-75*v2)) &
                     +3*c*u*(2+2*u4-35*v2+40*v4+u2*(4-35*v2))))/rho15
 fuv(jind(9,6,2)) = (-3*u5*v*(2*u4-10*c*u3*v+cuv*(11+14*v2)+u2*(-2+(-17+5*c2)*v2)+2*(-2+5*v2+7*v4)))/rho13

  fu(jind(9,5,3)) = (u4*v3*(5-4*u2+cuv+5*v2))/rho11
fuuu(jind(9,5,3)) = (3*u2*v3*(-40*u6+30*c*u5*v-60*cuv*(1+v2)**2+20*(1+v2)**3+15*u2*(1+v2)*(-11+(-11+c2)*v2) &
                     +6*u4*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,5,3)) = (3*u5*v2*(1+u2-cuv-2*v2))/rho11
fvvv(jind(9,5,3)) = (3*u5*(2+2*u6-42*c*u5*v-75*v2+240*v4-112*v6+c*u3*v*(-84-5*(-57+7*c2)*v2)+3*u4*(2+5*(-5+7*c2)*v2) &
                     -3*cuv*(14-95*v2+56*v4)-3*u2*(-2+(50-35*c2)*v2+40*(-2+c2)*v4)))/rho15
 fuv(jind(9,5,3)) = (-3*u4*v2*(4*u4-8*c*u3*v+cuv*(7+4*v2)+u2*(-1+(-19+c2)*v2)+5*(-1+v2+2*v4)))/rho13

  fu(jind(9,4,4)) = (u3*v4*(4 - 5*u2 - cuv + 4*v2))/rho11 
fuuu(jind(9,4,4)) = (-3*u*v4*(70*u6 + 42*c*u5*v + 60*cuv*(1 + v2)**2 - 8*(1 + v2)**3 + 5*c*u3*v*(-39 + (-39 + c2)*v2)  &
                     + 21*u4*(-11 + (-11 + c2)*v2) - 60*u2*(1 + v2)*(-2 + (-2 + c2)*v2)))/rho15 
fv(jind(9,4,4)) = (u4*v3*(4 + 4*u2 - cuv - 5*v2))/rho11
fvvv(jind(9,4,4)) = (3*u4*v*(8 + 8*u6 - 60*c*u5*v - 120*v2 + 231*v4 - 70*v6 - 5*c*u3*v*(24 + (-39 + c2)*v2)  &
                     + 12*u4*(2+5*(-2+c2)*v2)-3*cuv*(20-65*v2+14*v4)-3*u2*(-8 - 20*(-4 + c2)*v2 + 7*(-11 + c2)*v4)))/rho15
 fuv(jind(9,4,4)) = (u3*v3*(-20*u4+10*c*u3*v+cuv*(-17+10*v2)+u2*(-4+(59+c2)*v2)-4*(-4+v2+5*v4)))/rho13

  fu(jind(9,3,5)) = (3*u2*v5*(1 - 2*u2 - cuv + v2))/rho11
fuuu(jind(9,3,5)) = (3*v5*(-112*u6 - 168*c*u5*v - 42*cuv*(1 + v2)**2 + 2*(1 + v2)**3 - 120*u4*(-2 + (-2 + c2)*v2)  &
                     - 5*c*u3*v*(-57 + (-57 + 7*c2)*v2) + 15*u2*(1 + v2)*(-5 + (-5 + 7*c2)*v2)))/rho15
  fv(jind(9,3,5)) = (u3*v4*(5+5*u2+cuv-4*v2))/rho11
fvvv(jind(9,3,5)) = (3*u3*v2*(20+20*u6-60*c*u5*v-165*v2+204*v4-40*v6+15*u4*(4+(-11+c2)*v2)+c*u3*v*(-120+(69+c2)*v2) &
                     +3*cuv*(-20+23*v2+10*v4)+3*u2*(20+5*(-22+c2)*v2+2*(34+c2)*v4)))/rho15
 fuv(jind(9,3,5)) = (-3*u2*v4*(-5+10*u4+4*c*u3*v-v2+4*v4+cuv*(7-8*v2)+u2*(5+(-19+c2)*v2)))/rho13

  fu(jind(9,2,6)) = (u*v6*(2-7*u2-5*cuv+2*v2))/rho11
fuuu(jind(9,2,6)) = (-9*v6*(56*u5+120*c*u4*v+6*c*v*(1+v2)**2+35*c*u2*v*(-3+(-3+c2)*v2)+15*u3*(-5+(-5+7*c2)*v2) &
                     -6*u*(1+v2)*(-2+(-2+7*c2)*v2)))/rho15
  fv(jind(9,2,6)) = (3*u2*v5*(2+2*u2+cuv-v2))/rho11
fvvv(jind(9,2,6)) = (-3*u2*v3*(-40-40*u6+30*c*u5*v+204*v2-165*v4+20*v6+6*u4*(-20+(34+c2)*v2)+c*u3*v*(60+(69+c2)*v2) &
                     +3*cuv*(10+23*v2-20*v4)+3*u2*(-40+2*(68+c2)*v2+5*(-11+c2)*v4)))/rho15
 fuv(jind(9,2,6)) = (-3*u*v5*(14*u4+14*c*u3*v+cuv*(11-10*v2)+u2*(10+(-17+5*c2)*v2)+2*(-2-v2+v4)))/rho13

  fu(jind(9,1,7)) = (v7*(1-8*u2-7*cuv+v2))/rho11
fuuu(jind(9,1,7)) = (9*v7*(-80*u4-210*c*u3*v+cuv*(87+(87-77*c2)*v2)-30*u2*(-2+(-2+7*c2)*v2)+3*(1+v2)*(-1+(-1+11*c2)*v2)))/rho15
  fv(jind(9,1,7)) = (u*v6*(7+7*u2+5*cuv-2*v2))/rho11
fvvv(jind(9,1,7)) = (3*u*v4*(70+70*u6+42*c*u5*v-231*v2+120*v4-8*v6+c*u3*v*(84+5*(-39+c2)*v2) &
                     +21*u4*(10+(-11+c2)*v2)+3*cuv*(14-65*v2+20*v4)-3*u2*(-70-7*(-22+c2)*v2+20*(-2+c2)*v4)))/rho15
 fuv(jind(9,1,7)) = -((v6*(-7+56*u4+80*c*u3*v-5*v2+2*v4+cuv*(53-28*v2)+u2*(49+(-41+35*c2)*v2)))/rho13)

  fu(jind(9,0,8)) = (-9*v8*(u+c*v))/rho11
fuuu(jind(9,0,8)) = (-99*v8*(u+c*v)*(-3+10*u2+20*cuv+(-3+13*c2)*v2))/rho15
  fv(jind(9,0,8)) = (v7*(8+8*u2+7*cuv-v2))/rho11
fvvv(jind(9,0,8)) = (3*v5*(112+112*u6+168*c*u5*v-240*v2+75*v4-2*v6+24*u4*(14+5*(-2+c2)*v2)+c*u3*v*(336+5*(-57+7*c2)*v2) &
                     +3*cuv*(56-95*v2+14*v4)+u2*(336+120*(-4+c2)*v2-15*(-5+7*c2)*v4)))/rho15
 fuv(jind(9,0,8)) = (-9*v7*(8*u3+14*c*u2*v+c*v*(9-2*v2)+u*(8+(-3+7*c2)*v2)))/rho13

  fu(jind(9,9,0)) = (9*u8*(1+cuv+v2))/rho11
fuuu(jind(9,9,0)) = (9*u6*(6*c*u5*v+120*cuv*(1+v2)**2+56*(1+v2)**3+35*c*u3*v*(-3+(-3+c2)*v2)+15*u2*(1+v2)*(-5+(-5+7*c2)*v2) &
                     -6*u4*(-2+(-2+7*c2)*v2)))/rho15
  fv(jind(9,9,0)) = (-9*u9*(c*u+v))/rho11
fvvv(jind(9,9,0)) = (-99*u9*(c*u+v)*(-3+(-3+13*c2)*u2+20*cuv+10*v2))/rho15
 fuv(jind(9,9,0)) = (9*u8*(2*u2*v-9*c2*u2*v+c*u*(-10+u2-17*v2)-9*v*(1+v2)))/rho13

  fu(jind(9,8,1)) = (u7*v*(8-u2+7*cuv+8*v2))/rho11
fuuu(jind(9,8,1)) = (3*u5*v*(-2*u6+42*c*u5*v+168*cuv*(1+v2)**2+112*(1+v2)**3+120*u2*(1+v2)*(-2+(-2+c2)*v2) &
                     +5*c*u3*v*(-57+(-57+7*c2)*v2)-15*u4*(-5+(-5+7*c2)*v2)))/rho15
  fv(jind(9,8,1)) = (u8*(1+u2-7*cuv-8*v2))/rho11
fvvv(jind(9,8,1)) = (9*u8*(-3+(-3+33*c2)*u4+c*(87-77*c2)*u3*v+60*v2-80*v4+3*cuv*(29-70*v2)+u2*(-6+60*v2+c2*(33-210*v2))))/rho15
 fuv(jind(9,8,1)) = -((u7*(u4-23*c*u3*v+cuv*(58+103*v2)+u2*(-7+(-34+49*c2)*v2)+8*(-1+7*v2+8*v4)))/rho13)

  fu(jind(9,7,2)) = (u6*v2*(7-2*u2+5*cuv+7*v2))/rho11
fuuu(jind(9,7,2)) = (3*u4*v2*(-8*u6+60*c*u5*v+42*cuv*(1+v2)**2+70*(1+v2)**3+5*c*u3*v*(-39+(-39+c2)*v2) &
                     +21*u2*(1+v2)*(-11+(-11+c2)*v2)-60*u4*(-2+(-2+c2)*v2)))/rho15
  fv(jind(9,7,2)) = (u7*v*(2+2*u2-5*cuv-7*v2))/rho11
fvvv(jind(9,7,2)) = (-9*u7*(35*c**3*u3*v2-21*c2*u2*v*(2+2*u2-5*v2)+v*(12+12*u4-75*v2+56*v4+u2*(24-75*v2)) &
                     +3*c*u*(2+2*u4-35*v2+40*v4+u2*(4-35*v2))))/rho15
 fuv(jind(9,7,2)) = -((u6*v*(4*u4-29*c*u3*v+cuv*(34+61*v2)+u2*(-10+(-46+25*c2)*v2)+7*(-2+5*v2+7*v4)))/rho13)

  fu(jind(9,6,3)) = (3*u5*v3*(2-u2+cuv+2*v2))/rho11
fuuu(jind(9,6,3)) = (-3*u3*v3*(20*u6-60*c*u5*v+30*cuv*(1+v2)**2-40*(1+v2)**3+15*u4*(-11+(-11+c2)*v2) &
                     +6*u2*(1+v2)*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,6,3)) = (3*u6*v2*(1+u2-cuv-2*v2))/rho11
fvvv(jind(9,6,3)) = (3*u6*(2+2*u6-42*c*u5*v-75*v2+240*v4-112*v6+c*u3*v*(-84-5*(-57+7*c2)*v2)+3*u4*(2+5*(-5+7*c2)*v2) &
                     -3*cuv*(14-95*v2+56*v4)-3*u2*(-2+(50-35*c2)*v2+40*(-2+c2)*v4)))/rho15
 fuv(jind(9,6,3)) = (-9*u5*v2*(u4-3*c*u3*v+cuv*(2+3*v2)+u2*(-1+(-6+c2)*v2)+2*(-1+v2+2*v4)))/rho13

  fu(jind(9,5,4)) = (u4*v4*(5-4*u2+cuv+5*v2))/rho11
fuuu(jind(9,5,4)) = (3*u2*v4*(-40*u6+30*c*u5*v-60*cuv*(1+v2)**2+20*(1+v2)**3+15*u2*(1+v2)*(-11+(-11+c2)*v2) &
                     +6*u4*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,5,4)) = (u5*v3*(4+4*u2-cuv-5*v2))/rho11
fvvv(jind(9,5,4)) = (3*u5*v*(8+8*u6-60*c*u5*v-120*v2+231*v4-70*v6-5*c*u3*v*(24+(-39+c2)*v2)+12*u4*(2+5*(-2+c2)*v2) &
                     -3*cuv*(20-65*v2+14*v4)-3*u2*(-8-20*(-4+c2)*v2+7*(-11+c2)*v4)))/rho15
 fuv(jind(9,5,4)) = -((u4*v3*(16*u4-17*c*u3*v+cuv*(10+v2)+u2*(-4+(-58+c2)*v2)+5*(-4+v2+5*v4)))/rho13)

  fu(jind(9,4,5)) = (u3*v5*(4 - 5*u2 - cuv + 4*v2))/rho11
fuuu(jind(9,4,5)) = (-3*u*v5*(70*u6+42*c*u5*v+60*cuv*(1+v2)**2-8*(1+v2)**3+5*c*u3*v*(-39+(-39+c2)*v2) &
                     +21*u4*(-11+(-11+c2)*v2)-60*u2*(1+v2)*(-2+(-2+c2)*v2)))/rho15
  fv(jind(9,4,5)) = (u4*v4*(5+5*u2+cuv-4*v2))/rho11
fvvv(jind(9,4,5)) = (3*u4*v2*(20+20*u6-60*c*u5*v-165*v2+204*v4-40*v6+15*u4*(4+(-11+c2)*v2)+c*u3*v*(-120+(69+c2)*v2) &
                     +3*cuv*(-20+23*v2+10*v4)+3*u2*(20+5*(-22+c2)*v2+2*(34+c2)*v4)))/rho15
 fuv(jind(9,4,5)) = -((u3*v4*(25*u4+c*u3*v+cuv*(10-17*v2)+u2*(5+(-58+c2)*v2)+4*(-5-v2+4*v4)))/rho13)

  fu(jind(9,3,6)) = (3*u2*v6*(1 - 2*u2 - cuv + v2))/rho11
fuuu(jind(9,3,6)) = (3*v6*(-112*u6 - 168*c*u5*v - 42*cuv*(1 + v2)**2 + 2*(1 + v2)**3 - 120*u4*(-2 + (-2 + c2)*v2)  &
                     - 5*c*u3*v*(-57 + (-57 + 7*c2)*v2) + 15*u2*(1 + v2)*(-5 + (-5 + 7*c2)*v2)))/rho15
  fv(jind(9,3,6)) = (3*u3*v5*(2+2*u2+cuv-v2))/rho11
fvvv(jind(9,3,6)) = (-3*u3*v3*(-40-40*u6+30*c*u5*v+204*v2-165*v4+20*v6+6*u4*(-20+(34+c2)*v2)+c*u3*v*(60+(69+c2)*v2) &
                     +3*cuv*(10+23*v2-20*v4)+3*u2*(-40+2*(68+c2)*v2+5*(-11+c2)*v4)))/rho15
 fuv(jind(9,3,6)) = (-9*u2*v5*(-2+4*u4+3*c*u3*v-v2+v4+cuv*(2-3*v2)+u2*(2+(-6+c2)*v2)))/rho13

  fu(jind(9,2,7)) = (u*v7*(2-7*u2-5*cuv+2*v2))/rho11
fuuu(jind(9,2,7)) = (-9*v7*(56*u5+120*c*u4*v+6*c*v*(1+v2)**2+35*c*u2*v*(-3+(-3+c2)*v2)+15*u3*(-5+(-5+7*c2)*v2) &
                     -6*u*(1+v2)*(-2+(-2+7*c2)*v2)))/rho15
  fv(jind(9,2,7)) = (u2*v6*(7+7*u2+5*cuv-2*v2))/rho11
fvvv(jind(9,2,7)) = (3*u2*v4*(70+70*u6+42*c*u5*v-231*v2+120*v4-8*v6+c*u3*v*(84+5*(-39+c2)*v2) &
                     +21*u4*(10+(-11+c2)*v2)+3*cuv*(14-65*v2+20*v4)-3*u2*(-70-7*(-22+c2)*v2+20*(-2+c2)*v4)))/rho15
 fuv(jind(9,2,7)) = -((u*v6*(49*u4+61*c*u3*v+cuv*(34-29*v2)+u2*(35+(-46+25*c2)*v2)+2*(-7-5*v2+2*v4)))/rho13)

  fu(jind(9,1,8)) = (v8*(1-8*u2-7*cuv+v2))/rho11
fuuu(jind(9,1,8)) = (9*v8*(-80*u4-210*c*u3*v+cuv*(87+(87-77*c2)*v2)-30*u2*(-2+(-2+7*c2)*v2)+3*(1+v2)*(-1+(-1+11*c2)*v2)))/rho15
  fv(jind(9,1,8)) = (u*v7*(8+8*u2+7*cuv-v2))/rho11
fvvv(jind(9,1,8)) = (3*u*v5*(112+112*u6+168*c*u5*v-240*v2+75*v4-2*v6+24*u4*(14+5*(-2+c2)*v2)+c*u3*v*(336+5*(-57+7*c2)*v2) &
                     +3*cuv*(56-95*v2+14*v4)+u2*(336+120*(-4+c2)*v2-15*(-5+7*c2)*v4)))/rho15
 fuv(jind(9,1,8)) = -((v7*(-8+64*u4+103*c*u3*v-7*v2+v4+cuv*(58-23*v2)+u2*(56+(-34+49*c2)*v2)))/rho13)

  fu(jind(9,0,9)) = (-9*v9*(u + c*v))/rho11
fuuu(jind(9,0,9)) = (-99*v9*(u + c*v)*(-3 + 10*u2 + 20*cuv + (-3 + 13*c2)*v2))/rho15
  fv(jind(9,0,9)) = (9*v8*(1+u2+cuv))/rho11
fvvv(jind(9,0,9)) = (9*v6*(56+56*u6+120*c*u5*v-75*v2+12*v4+5*c*u3*v*(48+7*(-3+c2)*v2)+3*u4*(56+5*(-5+7*c2)*v2) &
                     +3*cuv*(40-35*v2+2*v4)-3*u2*(-56+(50-35*c2)*v2+2*(-2+7*c2)*v4)))/rho15
 fuv(jind(9,0,9)) = (-9*v8*(9*u3+17*c*u2*v-c*v*(-10+v2)+u*(9+(-2+9*c2)*v2)))/rho13

  fu(jind(9,10,0)) = (u9*(u2+11*cuv+10*(1+v2)))/rho11
fuuu(jind(9,10,0)) = (9*u7*(210*cuv*(1+v2)**2+80*(1+v2)**3+u4*(3+(3-33*c2)*v2)+30*u2*(1+v2)*(-2+(-2+7*c2)*v2) &
                     +c*u3*v*(-87+(-87+77*c2)*v2)))/rho15
  fv(jind(9,10,0)) = (-9*u10*(c*u+v))/rho11
fvvv(jind(9,10,0)) = (-99*u10*(c*u+v)*(-3+(-3+13*c2)*u2+20*cuv+10*v2))/rho15
 fuv(jind(9,10,0)) = (-9*u9*(-(u2*v)+11*c2*u2*v+10*v*(1+v2)+c*u*(11+20*v2)))/rho13

  fu(jind(9,9,1)) = (9*u8*v*(1+cuv+v2))/rho11
fuuu(jind(9,9,1)) = (9*u6*v*(6*c*u5*v+120*cuv*(1+v2)**2+56*(1+v2)**3+35*c*u3*v*(-3+(-3+c2)*v2)+15*u2*(1+v2)*(-5+(-5+7*c2)*v2) &
                     -6*u4*(-2+(-2+7*c2)*v2)))/rho15
  fv(jind(9,9,1)) = (u9*(1+u2-7*cuv-8*v2))/rho11
fvvv(jind(9,9,1)) = (9*u9*(-3+(-3+33*c2)*u4+c*(87-77*c2)*u3*v+60*v2-80*v4+3*cuv*(29-70*v2)+u2*(-6+60*v2+c2*(33-210*v2))))/rho15
 fuv(jind(9,9,1)) = (9*u8*(1+2*c*u3*v-7*v2-8*v4-7*cuv*(1+2*v2)+u2*(1+(3-7*c2)*v2)))/rho13

  fu(jind(9,8,2)) = (u7*v2*(8-u2+7*cuv+8*v2))/rho11
fuuu(jind(9,8,2)) = (3*u5*v2*(-2*u6+42*c*u5*v+168*cuv*(1+v2)**2+112*(1+v2)**3+120*u2*(1+v2)*(-2+(-2+c2)*v2) &
                     +5*c*u3*v*(-57+(-57+7*c2)*v2)-15*u4*(-5+(-5+7*c2)*v2)))/rho15
  fv(jind(9,8,2)) = (u8*v*(2+2*u2-5*cuv-7*v2))/rho11
fvvv(jind(9,8,2)) = (-9*u8*(35*c**3*u3*v2-21*c2*u2*v*(2+2*u2-5*v2)+v*(12+12*u4-75*v2+56*v4+u2*(24-75*v2)) &
                     +3*c*u*(2+2*u4-35*v2+40*v4+u2*(4-35*v2))))/rho15
 fuv(jind(9,8,2)) = -((u7*v*(2*u4-28*c*u3*v+5*cuv*(7+16*v2)+u2*(-14+(-41+35*c2)*v2)+8*(-2+5*v2+7*v4)))/rho13)

  fu(jind(9,7,3)) = (u6*v3*(7-2*u2+5*cuv+7*v2))/rho11
fuuu(jind(9,7,3)) = (3*u4*v3*(-8*u6+60*c*u5*v+42*cuv*(1+v2)**2+70*(1+v2)**3+5*c*u3*v*(-39+(-39+c2)*v2) &
                     +21*u2*(1+v2)*(-11+(-11+c2)*v2)-60*u4*(-2+(-2+c2)*v2)))/rho15
  fv(jind(9,7,3)) = (3*u7*v2*(1+u2-cuv-2*v2))/rho11
fvvv(jind(9,7,3)) = (3*u7*(2+2*u6-42*c*u5*v-75*v2+240*v4-112*v6+c*u3*v*(-84-5*(-57+7*c2)*v2)+3*u4*(2+5*(-5+7*c2)*v2) &
                     -3*cuv*(14-95*v2+56*v4)-3*u2*(-2+(50-35*c2)*v2+40*(-2+c2)*v4)))/rho15
 fuv(jind(9,7,3)) = (-3*u6*v2*(2*u4-10*c*u3*v+cuv*(5+14*v2)+u2*(-5+(-17+5*c2)*v2)+7*(-1+v2+2*v4)))/rho13

  fu(jind(9,6,4)) = (3*u5*v4*(2-u2+cuv+2*v2))/rho11
fuuu(jind(9,6,4)) = (-3*u3*v4*(20*u6-60*c*u5*v+30*cuv*(1+v2)**2-40*(1+v2)**3+15*u4*(-11+(-11+c2)*v2) &
                     +6*u2*(1+v2)*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,6,4)) = (u6*v3*(4+4*u2-cuv-5*v2))/rho11
fvvv(jind(9,6,4)) = (3*u6*v*(8+8*u6-60*c*u5*v-120*v2+231*v4-70*v6-5*c*u3*v*(24+(-39+c2)*v2)+12*u4*(2+5*(-2+c2)*v2) &
                     -3*cuv*(20-65*v2+14*v4)-3*u2*(-8-20*(-4+c2)*v2+7*(-11+c2)*v4)))/rho15
 fuv(jind(9,6,4)) = (-3*u5*v3*(4*u4-8*c*u3*v+cuv*(1+4*v2)+u2*(-4+(-19+c2)*v2)+2*(-4+v2+5*v4)))/rho13

  fu(jind(9,5,5)) = (u4*v5*(5-4*u2+cuv+5*v2))/rho11
fuuu(jind(9,5,5)) = (3*u2*v5*(-40*u6+30*c*u5*v-60*cuv*(1+v2)**2+20*(1+v2)**3+15*u2*(1+v2)*(-11+(-11+c2)*v2) &
                     +6*u4*(34+(34+c2)*v2)+c*u3*v*(69+(69+c2)*v2)))/rho15
  fv(jind(9,5,5)) = (u5*v4*(5+5*u2+cuv-4*v2))/rho11
fvvv(jind(9,5,5)) = (3*u5*v2*(20+20*u6-60*c*u5*v-165*v2+204*v4-40*v6+15*u4*(4+(-11+c2)*v2)+c*u3*v*(-120+(69+c2)*v2) &
                     +3*cuv*(-20+23*v2+10*v4)+3*u2*(20+5*(-22+c2)*v2+2*(34+c2)*v4)))/rho15
 fuv(jind(9,5,5)) = (u4*v4*(-20*u4+10*c*u3*v+cuv*(1+10*v2)+u2*(5+(59+c2)*v2)+5*(5+v2-4*v4)))/rho13

  fu(jind(9,4,6)) = (u3*v6*(4-5*u2-cuv+4*v2))/rho11
fuuu(jind(9,4,6)) = (-3*u*v6*(70*u6+42*c*u5*v+60*cuv*(1+v2)**2-8*(1+v2)**3 &
                     +5*c*u3*v*(-39+(-39+c2)*v2) &
                     +21*u4*(-11+(-11+c2)*v2)-60*u2*(1+v2)*(-2+(-2+c2)*v2)))/rho15
  fv(jind(9,4,6)) = (3*u4*v5*(2+2*u2+cuv-v2))/rho11
fvvv(jind(9,4,6)) = (-3*u4*v3*(-40-40*u6+30*c*u5*v+204*v2-165*v4+20*v6+6*u4*(-20+(34+c2)*v2)+c*u3*v*(60+(69+c2)*v2) &
                     +3*cuv*(10+23*v2-20*v4)+3*u2*(-40+2*(68+c2)*v2+5*(-11+c2)*v4)))/rho15
 fuv(jind(9,4,6)) = (-3*u3*v5*(10*u4+4*c*u3*v+u2*(2+(-19+c2)*v2)+c*u*(v-8*v3)+4*(-2-v2+v4)))/rho13

  fu(jind(9,3,7)) = (3*u2*v7*(1-2*u2-cuv+v2))/rho11
fuuu(jind(9,3,7)) = (3*v7*(-112*u6-168*c*u5*v-42*cuv*(1+v2)**2+2*(1+v2)**3-120*u4*(-2+(-2+c2)*v2) &
                     -5*c*u3*v*(-57+(-57+7*c2)*v2)+15*u2*(1+v2)*(-5+(-5+7*c2)*v2)))/rho15
  fv(jind(9,3,7)) = (u3*v6*(7+7*u2+5*cuv-2*v2))/rho11
fvvv(jind(9,3,7)) = (3*u3*v4*(70+70*u6+42*c*u5*v-231*v2+120*v4-8*v6+c*u3*v*(84+5*(-39+c2)*v2) &
                     +21*u4*(10+(-11+c2)*v2)+3*cuv*(14-65*v2+20*v4)-3*u2*(-70-7*(-22+c2)*v2+20*(-2+c2)*v4)))/rho15
 fuv(jind(9,3,7)) = (-3*u2*v6*(-7+14*u4+14*c*u3*v-5*v2+2*v4+5*cuv*(1-2*v2)+u2*(7+(-17+5*c2)*v2)))/rho13

  fu(jind(9,2,8)) = (u*v8*(2-7*u2-5*cuv+2*v2))/rho11
fuuu(jind(9,2,8)) = (-9*v8*(56*u5+120*c*u4*v+6*c*v*(1+v2)**2+35*c*u2*v*(-3+(-3+c2)*v2)+15*u3*(-5+(-5+7*c2)*v2) &
                     -6*u*(1+v2)*(-2+(-2+7*c2)*v2)))/rho15
  fv(jind(9,2,8)) = (u2*v7*(8+8*u2+7*cuv-v2))/rho11
fvvv(jind(9,2,8)) = (3*u2*v5*(112+112*u6+168*c*u5*v-240*v2+75*v4-2*v6+24*u4*(14+5*(-2+c2)*v2)+c*u3*v*(336+5*(-57+7*c2)*v2) &
                     +3*cuv*(56-95*v2+14*v4)+u2*(336+120*(-4+c2)*v2-15*(-5+7*c2)*v4)))/rho15
 fuv(jind(9,2,8)) = -((u*v7*(56*u4+80*c*u3*v-7*cuv*(-5+4*v2)+u2*(40+(-41+35*c2)*v2)+2*(-8-7*v2+v4)))/rho13)

  fu(jind(9,1,9)) = (v9*(1-8*u2-7*cuv+v2))/rho11
fuuu(jind(9,1,9)) = (9*v9*(-80*u4-210*c*u3*v+cuv*(87+(87-77*c2)*v2)-30*u2*(-2+(-2+7*c2)*v2)+3*(1+v2)*(-1+(-1+11*c2)*v2)))/rho15
  fv(jind(9,1,9)) = (9*u*v8*(1+u2+cuv))/rho11
fvvv(jind(9,1,9)) = (9*u*v6*(56+56*u6+120*c*u5*v-75*v2+12*v4+5*c*u3*v*(48+7*(-3+c2)*v2)+3*u4*(56+5*(-5+7*c2)*v2) &
                     +3*cuv*(40-35*v2+2*v4)-3*u2*(-56+(50-35*c2)*v2+2*(-2+7*c2)*v4)))/rho15
 fuv(jind(9,1,9)) = (9*v8*(1-8*u4-14*c*u3*v+v2+cuv*(-7+2*v2)+u2*(-7+(3-7*c2)*v2)))/rho13

  fu(jind(9,0,10)) = (-9*v10*(u+c*v))/rho11
fuuu(jind(9,0,10)) = (-99*v10*(u+c*v)*(-3+10*u2+20*cuv+(-3+13*c2)*v2))/rho15
  fv(jind(9,0,10)) = (v9*(10+10*u2+11*cuv+v2))/rho11
fvvv(jind(9,0,10)) = (9*v7*(80+80*u6+210*c*u5*v-60*v2+3*v4+3*cuv*(70-29*v2)+30*u4*(8+(-2+7*c2)*v2)+c*u3*v*(420+(-87+77*c2)*v2) &
                     -3*u2*(-80+(40-70*c2)*v2+(-1+11*c2)*v4)))/rho15
 fuv(jind(9,0,10)) = (-9*v9*(10*u3+11*c*v+20*c*u2*v+u*(10+(-1+11*c2)*v2)))/rho13

  fu(jind(11,9,0)) = (u8*(9-2*u2+7*cuv+9*v2))/rho13
fuuu(jind(11,9,0)) = (3*u6*(-8*u6+84*c*u5*v+216*cuv*(1+v2)**2+168*(1+v2)**3+5*c*u3*v*(-83+(-83+7*c2)*v2) &
                     +9*u2*(1+v2)*(-43+(-43+15*c2)*v2)-4*u4*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,9,0)) = (-11*u9*(c*u+v))/rho13
fvvv(jind(11,9,0)) = (-429*u9*(c*u+v)*(-1+(-1+5*c2)*u2+8*cuv+4*v2))/rho17
 fuv(jind(11,9,0)) = (11*u8*(3*c*u3+4*u2*v-7*c2*u2*v-9*v*(1+v2)-5*c*u*(2+3*v2)))/rho15

  fu(jind(11,8,1)) = (u7*v*(8-3*u2+5*cuv+8*v2))/rho13
fuuu(jind(11,8,1)) = (3*u5*v*(-20*u6+100*c*u5*v+56*cuv*(1+v2)**2+112*(1+v2)**3+u4*(215+(215-75*c2)*v2) &
                     +8*u2*(1+v2)*(-46+(-46+3*c2)*v2)+c*u3*v*(-273+(-273+5*c2)*v2)))/rho17
  fv(jind(11,8,1)) = (u8*(1+u2-9*cuv-10*v2))/rho13
fvvv(jind(11,8,1)) = (33*u8*(-1+(-1+13*c2)*u4+c*(35-39*c2)*u3*v+24*v2-40*v4+cuv*(35-108*v2)+u2*(-2+24*v2+c2*(13-108*v2))))/rho17
 fuv(jind(11,8,1)) = -((u7*(-8+3*u4-43*c*u3*v+72*v2+80*v4+3*cuv*(26+37*v2)+5*u2*(-1+3*(-4+3*c2)*v2)))/rho15)

  fu(jind(11,7,2)) = (u6*v2*(7-4*u2+3*cuv+7*v2))/rho13
fuuu(jind(11,7,2)) = (-3*u4*v2*(40*u6-90*c*u5*v+42*cuv*(1+v2)**2-70*(1+v2)**3+7*u2*(1+v2)*(47+(47+c2)*v2) &
                     +c*u3*v*(95+(95+c2)*v2)+6*u4*(-46+(-46+3*c2)*v2)))/rho17
  fv(jind(11,7,2)) = (u7*v*(2+2*u2-7*cuv-9*v2))/rho13
fvvv(jind(11,7,2)) = (-33*u7*(21*c**3*u3*v2-9*c2*u2*v*(2+2*u2-7*v2)+c*u*(2+2*u4-45*v2+70*v4+u2*(4-45*v2)) &
                     +v*(4+4*u4-31*v2+30*v4+u2*(8-31*v2))))/rho17
 fuv(jind(11,7,2)) = -((u6*v*(8*u4-45*c*u3*v+cuv*(54+65*v2)+3*u2*(-2+(-24+7*c2)*v2)+7*(-2+7*v2+9*v4)))/rho15)

  fu(jind(11,6,3)) = (u5*v3*(6-5*u2+cuv+6*v2))/rho13
fuuu(jind(11,6,3)) = (3*u3*v3*(-70*u6+42*c*u5*v-90*cuv*(1+v2)**2+40*(1+v2)**3+7*u4*(47+(47+c2)*v2)+c*u3*v*(95+(95+c2)*v2) &
                     +6*u2*(1+v2)*(-46+(-46+3*c2)*v2)))/rho17
  fv(jind(11,6,3)) = (u6*v2*(3+3*u2-5*cuv-8*v2))/rho13
fvvv(jind(11,6,3)) = (3*u6*(2+2*u6-54*c*u5*v-93*v2+380*v4-240*v6+c*u3*v*(-108-7*(-71+15*c2)*v2)+3*u4*(2+(-31+63*c2)*v2) &
                     +cuv*(-54+497*v2-450*v4)+u2*(6+3*(-62+63*c2)*v2+(380-350*c2)*v4)))/rho17
 fuv(jind(11,6,3)) = -((u5*v2*(15*u4-39*c*u3*v+cuv*(38+27*v2)+u2*(-3+5*(-16+c2)*v2)+6*(-3+5*v2+8*v4)))/rho15)

  fu(jind(11,5,4)) = (u4*v4*(5-6*u2-cuv+5*v2))/rho13
fuuu(jind(11,5,4)) = (3*u2*v4*(-112*u6-56*c*u5*v-100*cuv*(1+v2)**2+20*(1+v2)**3+c*u3*v*(273+(273-5*c2)*v2) &
                     -8*u4*(-46+(-46+3*c2)*v2)+5*u2*(1+v2)*(-43+(-43+15*c2)*v2)))/rho17
  fv(jind(11,5,4)) = (u5*v3*(4+4*u2-3*cuv-7*v2))/rho13
fvvv(jind(11,5,4)) = (3*u5*v*(8+8*u6-84*c*u5*v-152*v2+387*v4-168*v6+c*u3*v*(-168-5*(-83+7*c2)*v2)+4*u4*(6+(-38+35*c2)*v2) &
                     +cuv*(-84+415*v2-216*v4)+u2*(24+4*(-76+35*c2)*v2-9*(-43+15*c2)*v4)))/rho17
 fuv(jind(11,5,4)) = (u4*v3*(-24*u4+25*c*u3*v+3*cuv*(-10+v2)+u2*(-4+3*(28+c2)*v2)-5*(-4+3*v2+7*v4)))/rho15

  fu(jind(11,4,5)) = (u3*v5*(4-7*u2-3*cuv+4*v2))/rho13
fuuu(jind(11,4,5)) = (-3*u*v5*(168*u6+216*c*u5*v+84*cuv*(1+v2)**2-8*(1+v2)**3+5*c*u3*v*(-83+(-83+7*c2)*v2) &
                     +9*u4*(-43+(-43+15*c2)*v2)-4*u2*(1+v2)*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,4,5)) = (u4*v4*(5+5*u2-cuv-6*v2))/rho13
fvvv(jind(11,4,5)) = (3*u4*v2*(20+20*u6-100*c*u5*v-215*v2+368*v4-112*v6+c*u3*v*(-200+(273-5*c2)*v2)+5*u4*(12+(-43+15*c2)*v2) &
                     +cuv*(-100+273*v2-56*v4)+u2*(60+5*(-86+15*c2)*v2-8*(-46+3*c2)*v4)))/rho17
 fuv(jind(11,4,5)) = (u3*v4*(-35*u4+3*c*u3*v+5*cuv*(-6+5*v2)+3*u2*(-5+(28+c2)*v2)-4*(-5+v2+6*v4)))/rho15

  fu(jind(11,3,6)) = (u2*v6*(3-8*u2-5*cuv+3*v2))/rho13
fuuu(jind(11,3,6)) = (3*v6*(-240*u6-450*c*u5*v-54*cuv*(1+v2)**2+2*(1+v2)**3+u4*(380+(380-350*c2)*v2) &
                     -7*c*u3*v*(-71+(-71+15*c2)*v2)+3*u2*(1+v2)*(-31+(-31+63*c2)*v2)))/rho17
  fv(jind(11,3,6)) = (u3*v5*(6+6*u2+cuv-5*v2))/rho13
fvvv(jind(11,3,6)) = (3*u3*v3*(40+40*u6-90*c*u5*v-276*v2+329*v4-70*v6+c*u3*v*(-180+(95+c2)*v2)+6*u4*(20+(-46+3*c2)*v2) &
                     +cuv*(-90+95*v2+42*v4)+u2*(120+6*(-92+3*c2)*v2+7*(47+c2)*v4)))/rho17
 fuv(jind(11,3,6)) = -((u2*v5*(48*u4+27*c*u3*v+cuv*(38-39*v2)+5*u2*(6+(-16+c2)*v2)+3*(-6-v2+5*v4)))/rho15)

  fu(jind(11,2,7)) = (u*v7*(2-9*u2-7*cuv+2*v2))/rho13
fuuu(jind(11,2,7)) = (-33*v7*(30*u5+70*c*u4*v+2*c*v*(1+v2)**2+3*c*u2*v*(-15+(-15+7*c2)*v2)-2*u*(1+v2)*(-2+(-2+9*c2)*v2) &
                     +u3*(-31+(-31+63*c2)*v2)))/rho17
  fv(jind(11,2,7)) = (u2*v6*(7+7*u2+3*cuv-4*v2))/rho13
fvvv(jind(11,2,7)) = (-3*u2*v4*(-70-70*u6+42*c*u5*v+329*v2-276*v4+40*v6+7*u4*(-30+(47+c2)*v2)+c*u3*v*(84+(95+c2)*v2) &
                     +cuv*(42+95*v2-90*v4)+u2*(-210+7*(94+c2)*v2+6*(-46+3*c2)*v4)))/rho17
 fuv(jind(11,2,7)) = -((u*v6*(-14+63*u4+65*c*u3*v-6*v2+8*v4-9*cuv*(-6+5*v2)+u2*(49+3*(-24+7*c2)*v2)))/rho15)

  fu(jind(11,1,8)) = (v8*(1-10*u2-9*cuv+v2))/rho13
fuuu(jind(11,1,8)) = (33*v8*(-40*u4-108*c*u3*v+cuv*(35+(35-39*c2)*v2)-12*u2*(-2+(-2+9*c2)*v2)+(1+v2)*(-1+(-1+13*c2)*v2)))/rho17
  fv(jind(11,1,8)) = (u*v7*(8+8*u2+5*cuv-3*v2))/rho13
fvvv(jind(11,1,8)) = (3*u*v5*(112+112*u6+56*c*u5*v-368*v2+215*v4-20*v6+8*u4*(42+(-46+3*c2)*v2)+c*u3*v*(112+(-273+5*c2)*v2) &
                     +cuv*(56-273*v2+100*v4)+u2*(336+8*(-92+3*c2)*v2+(215-75*c2)*v4)))/rho17
 fuv(jind(11,1,8)) = -((v7*(-8+80*u4+111*c*u3*v-5*v2+3*v4+cuv*(78-43*v2)+u2*(72+15*(-4+3*c2)*v2)))/rho15)

  fu(jind(11,0,9)) = (-11*v9*(u+c*v))/rho13
fuuu(jind(11,0,9)) = (-429*v9*(u+c*v)*(-1+4*u2+8*cuv+(-1+5*c2)*v2))/rho17
  fv(jind(11,0,9)) = (v8*(9+9*u2+7*cuv-2*v2))/rho13
fvvv(jind(11,0,9)) = (3*v6*(168+168*u6+216*c*u5*v-387*v2+152*v4-8*v6+c*u3*v*(432+5*(-83+7*c2)*v2)+9*u4*(56+(-43+15*c2)*v2) &
                     +cuv*(216-415*v2+84*v4)+u2*(504+9*(-86+15*c2)*v2-4*(-38+35*c2)*v4)))/rho17
 fuv(jind(11,0,9)) = (-11*v8*(9*u3+15*c*u2*v+c*v*(10-3*v2)+u*(9+(-4+7*c2)*v2)))/rho15

  fu(jind(11,10,0)) = (u9*(-u2+9*cuv+10*(1+v2)))/rho13
fuuu(jind(11,10,0)) = (3*u7*(-2*u6+54*c*u5*v+450*cuv*(1+v2)**2+240*(1+v2)**3+u4*(93+(93-189*c2)*v2) &
                      +7*c*u3*v*(-71+(-71+15*c2)*v2)+10*u2*(1+v2)*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,10,0)) = (-11*u10*(c*u+v))/rho13
fvvv(jind(11,10,0)) = (-429*u10*(c*u+v)*(-1+(-1+5*c2)*u2+8*cuv+4*v2))/rho17
 fuv(jind(11,10,0)) = (11*u9*(3*u2*v-9*c2*u2*v+c*u*(-11+2*u2-18*v2)-10*v*(1+v2)))/rho15

  fu(jind(11,9,1)) = (u8*v*(9-2*u2+7*cuv+9*v2))/rho13
fuuu(jind(11,9,1)) = (3*u6*v*(-8*u6+84*c*u5*v+216*cuv*(1+v2)**2+168*(1+v2)**3+5*c*u3*v*(-83+(-83+7*c2)*v2) &
                      +9*u2*(1+v2)*(-43+(-43+15*c2)*v2)-4*u4*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,9,1)) = (u9*(1+u2-9*cuv-10*v2))/rho13
fvvv(jind(11,9,1)) = (33*u9*(-1+(-1+13*c2)*u4+c*(35-39*c2)*u3*v+24*v2-40*v4+cuv*(35-108*v2)+u2*(-2+24*v2+c2*(13-108*v2))))/rho17
 fuv(jind(11,9,1)) = -((u8*(-9+2*u4-36*c*u3*v+81*v2+90*v4+5*cuv*(17+28*v2)+u2*(-7+(-51+63*c2)*v2)))/rho15)

  fu(jind(11,8,2)) = (u7*v2*(8-3*u2+5*cuv+8*v2))/rho13
fuuu(jind(11,8,2)) = (3*u5*v2*(-20*u6+100*c*u5*v+56*cuv*(1+v2)**2+112*(1+v2)**3+u4*(215+(215-75*c2)*v2) &
                      +8*u2*(1+v2)*(-46+(-46+3*c2)*v2)+c*u3*v*(-273+(-273+5*c2)*v2)))/rho17
  fv(jind(11,8,2)) = (u8*v*(2+2*u2-7*cuv-9*v2))/rho13
fvvv(jind(11,8,2)) = (-33*u8*(21*c**3*u3*v2-9*c2*u2*v*(2+2*u2-7*v2)+c*u*(2+2*u4-45*v2+70*v4+u2*(4-45*v2)) &
                      +v*(4+4*u4-31*v2+30*v4+u2*(8-31*v2))))/rho17
 fuv(jind(11,8,2)) = -((u7*v*(6*u4-42*c*u3*v+3*cuv*(19+30*v2)+5*u2*(-2+(-13+7*c2)*v2)+8*(-2+7*v2+9*v4)))/rho15)

  fu(jind(11,7,3)) = (u6*v3*(7-4*u2+3*cuv+7*v2))/rho13
fuuu(jind(11,7,3)) = (-3*u4*v3*(40*u6-90*c*u5*v+42*cuv*(1+v2)**2-70*(1+v2)**3+7*u2*(1+v2)*(47+(47+c2)*v2) &
                      +c*u3*v*(95+(95+c2)*v2)+6*u4*(-46+(-46+3*c2)*v2)))/rho17
  fv(jind(11,7,3)) = (u7*v2*(3+3*u2-5*cuv-8*v2))/rho13
fvvv(jind(11,7,3)) = (3*u7*(2+2*u6-54*c*u5*v-93*v2+380*v4-240*v6+c*u3*v*(-108-7*(-71+15*c2)*v2)+3*u4*(2+(-31+63*c2)*v2) &
                      +cuv*(-54+497*v2-450*v4)+u2*(6+3*(-62+63*c2)*v2+(380-350*c2)*v4)))/rho17
 fuv(jind(11,7,3)) = -((u6*v2*(12*u4-40*c*u3*v+cuv*(37+48*v2)+3*u2*(-3+5*(-5+c2)*v2)+7*(-3+5*v2+8*v4)))/rho15)

  fu(jind(11,6,4)) = (u5*v4*(6-5*u2+cuv+6*v2))/rho13
fuuu(jind(11,6,4)) = (3*u3*v4*(-70*u6+42*c*u5*v-90*cuv*(1+v2)**2+40*(1+v2)**3+7*u4*(47+(47+c2)*v2)+c*u3*v*(95+(95+c2)*v2) &
                      +6*u2*(1+v2)*(-46+(-46+3*c2)*v2)))/rho17
  fv(jind(11,6,4)) = (u6*v3*(4+4*u2-3*cuv-7*v2))/rho13
fvvv(jind(11,6,4)) = (3*u6*v*(8+8*u6-84*c*u5*v-152*v2+387*v4-168*v6+c*u3*v*(-168-5*(-83+7*c2)*v2)+4*u4*(6+(-38+35*c2)*v2) &
                      +cuv*(-84+415*v2-216*v4)+u2*(24+4*(-76+35*c2)*v2-9*(-43+15*c2)*v4)))/rho17
 fuv(jind(11,6,4)) = -((u5*v3*(20*u4-30*c*u3*v+cuv*(25+14*v2)+u2*(-4+3*(-27+c2)*v2)+6*(-4+3*v2+7*v4)))/rho15)

  fu(jind(11,5,5)) = (u4*v5*(5-6*u2-cuv+5*v2))/rho13
fuuu(jind(11,5,5)) = (3*u2*v5*(-112*u6-56*c*u5*v-100*cuv*(1+v2)**2+20*(1+v2)**3+c*u3*v*(273+(273-5*c2)*v2) &
                      -8*u4*(-46+(-46+3*c2)*v2)+5*u2*(1+v2)*(-43+(-43+15*c2)*v2)))/rho17
  fv(jind(11,5,5)) = (u5*v4*(5+5*u2-cuv-6*v2))/rho13
fvvv(jind(11,5,5)) = (3*u5*v2*(20+20*u6-100*c*u5*v-215*v2+368*v4-112*v6+c*u3*v*(-200+(273-5*c2)*v2)+5*u4*(12+(-43+15*c2)*v2) &
                      +cuv*(-100+273*v2-56*v4)+u2*(60+5*(-86+15*c2)*v2-8*(-46+3*c2)*v4)))/rho17
 fuv(jind(11,5,5)) = (u4*v4*(-30*u4+12*c*u3*v+3*cuv*(-7+4*v2)+u2*(-5+(83+c2)*v2)-5*(-5+v2+6*v4)))/rho15

  fu(jind(11,4,6)) = (u3*v6*(4-7*u2-3*cuv+4*v2))/rho13
fuuu(jind(11,4,6)) = (-3*u*v6*(168*u6+216*c*u5*v+84*cuv*(1+v2)**2-8*(1+v2)**3+5*c*u3*v*(-83+(-83+7*c2)*v2) &
                      +9*u4*(-43+(-43+15*c2)*v2)-4*u2*(1+v2)*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,4,6)) = (u4*v5*(6+6*u2+cuv-5*v2))/rho13
fvvv(jind(11,4,6)) = (3*u4*v3*(40+40*u6-90*c*u5*v-276*v2+329*v4-70*v6+c*u3*v*(-180+(95+c2)*v2)+6*u4*(20+(-46+3*c2)*v2) &
                      +cuv*(-90+95*v2+42*v4)+u2*(120+6*(-92+3*c2)*v2+7*(47+c2)*v4)))/rho17
 fuv(jind(11,4,6)) = -((u3*v5*(42*u4+14*c*u3*v-5*cuv*(-5+6*v2)+3*u2*(6+(-27+c2)*v2)+4*(-6-v2+5*v4)))/rho15)

  fu(jind(11,3,7)) = (u2*v7*(3-8*u2-5*cuv+3*v2))/rho13
fuuu(jind(11,3,7)) = (3*v7*(-240*u6-450*c*u5*v-54*cuv*(1+v2)**2+2*(1+v2)**3+u4*(380+(380-350*c2)*v2) &
                      -7*c*u3*v*(-71+(-71+15*c2)*v2)+3*u2*(1+v2)*(-31+(-31+63*c2)*v2)))/rho17
  fv(jind(11,3,7)) = (u3*v6*(7+7*u2+3*cuv-4*v2))/rho13
fvvv(jind(11,3,7)) = (-3*u3*v4*(-70-70*u6+42*c*u5*v+329*v2-276*v4+40*v6+7*u4*(-30+(47+c2)*v2)+c*u3*v*(84+(95+c2)*v2) &
                      +cuv*(42+95*v2-90*v4)+u2*(-210+7*(94+c2)*v2+6*(-46+3*c2)*v4)))/rho17
 fuv(jind(11,3,7)) = -((u2*v6*(56*u4+48*c*u3*v+cuv*(37-40*v2)+5*u2*(7+3*(-5+c2)*v2)+3*(-7-3*v2+4*v4)))/rho15)

  fu(jind(11,2,8)) = (u*v8*(2-9*u2-7*cuv+2*v2))/rho13
fuuu(jind(11,2,8)) = (-33*v8*(30*u5+70*c*u4*v+2*c*v*(1+v2)**2+3*c*u2*v*(-15+(-15+7*c2)*v2)-2*u*(1+v2)*(-2+(-2+9*c2)*v2) &
                      +u3*(-31+(-31+63*c2)*v2)))/rho17
  fv(jind(11,2,8)) = (u2*v7*(8+8*u2+5*cuv-3*v2))/rho13
fvvv(jind(11,2,8)) = (3*u2*v5*(112+112*u6+56*c*u5*v-368*v2+215*v4-20*v6+8*u4*(42+(-46+3*c2)*v2)+c*u3*v*(112+(-273+5*c2)*v2) &
                      +cuv*(56-273*v2+100*v4)+u2*(336+8*(-92+3*c2)*v2+(215-75*c2)*v4)))/rho17
 fuv(jind(11,2,8)) = -((u*v7*(72*u4+90*c*u3*v+u2*(56+5*(-13+7*c2)*v2)+c*u*(57*v-42*v3)+2*(-8-5*v2+3*v4)))/rho15)

  fu(jind(11,1,9)) = (v9*(1-10*u2-9*cuv+v2))/rho13
fuuu(jind(11,1,9)) = (33*v9*(-40*u4-108*c*u3*v+cuv*(35+(35-39*c2)*v2)-12*u2*(-2+(-2+9*c2)*v2)+(1+v2)*(-1+(-1+13*c2)*v2)))/rho17
  fv(jind(11,1,9)) = (u*v8*(9+9*u2+7*cuv-2*v2))/rho13
fvvv(jind(11,1,9)) = (3*u*v6*(168+168*u6+216*c*u5*v-387*v2+152*v4-8*v6+c*u3*v*(432+5*(-83+7*c2)*v2)+9*u4*(56+(-43+15*c2)*v2) &
                      +cuv*(216-415*v2+84*v4)+u2*(504+9*(-86+15*c2)*v2-4*(-38+35*c2)*v4)))/rho17
 fuv(jind(11,1,9)) = -((v8*(-9+90*u4+140*c*u3*v-7*v2+2*v4+cuv*(85-36*v2)+u2*(81+(-51+63*c2)*v2)))/rho15)

  fu(jind(11,0,10)) = (-11*v10*(u+c*v))/rho13
fuuu(jind(11,0,10)) = (-429*v10*(u+c*v)*(-1+4*u2+8*cuv+(-1+5*c2)*v2))/rho17
  fv(jind(11,0,10)) = (v9*(10+10*u2+9*cuv-v2))/rho13
fvvv(jind(11,0,10)) = (3*v7*(240+240*u6+450*c*u5*v-380*v2+93*v4-2*v6+c*u3*v*(900+7*(-71+15*c2)*v2)+10*u4*(72+(-38+35*c2)*v2) &
                      +cuv*(450-497*v2+54*v4)+u2*(720+10*(-76+35*c2)*v2+(93-189*c2)*v4)))/rho17
 fuv(jind(11,0,10)) = (-11*v9*(10*u3+18*c*u2*v+c*v*(11-2*v2)+u*(10+(-3+9*c2)*v2)))/rho15

  fu(jind(11,11,0)) = (11*u10*(1+cuv+v2))/rho13
fuuu(jind(11,11,0)) = (33*u8*(2*c*u5*v+70*cuv*(1+v2)**2+30*(1+v2)**3+u4*(4+(4-18*c2)*v2)+3*c*u3*v*(-15+(-15+7*c2)*v2) &
                      +u2*(1+v2)*(-31+(-31+63*c2)*v2)))/rho17
  fv(jind(11,11,0)) = (-11*u11*(c*u+v))/rho13
fvvv(jind(11,11,0)) = (-429*u11*(c*u+v)*(-1+(-1+5*c2)*u2+8*cuv+4*v2))/rho17
 fuv(jind(11,11,0)) = (11*u10*(c*u3+2*u2*v-11*c2*u2*v-11*v*(1+v2)-3*c*u*(4+7*v2)))/rho15

  fu(jind(11,10,1)) = (u9*v*(-u2+9*cuv+10*(1+v2)))/rho13
fuuu(jind(11,10,1)) = (3*u7*v*(-2*u6+54*c*u5*v+450*cuv*(1+v2)**2+240*(1+v2)**3+u4*(93+(93-189*c2)*v2) &
                      +7*c*u3*v*(-71+(-71+15*c2)*v2)+10*u2*(1+v2)*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,10,1)) = (u10*(1+u2-9*cuv-10*v2))/rho13
fvvv(jind(11,10,1)) = (33*u10*(-1+(-1+13*c2)*u4+c*(35-39*c2)*u3*v+24*v2-40*v4+cuv*(35-108*v2)+u2*(-2+24*v2+c2*(13-108*v2))))/rho17
 fuv(jind(11,10,1)) = -((u9*(u4-29*c*u3*v+cuv*(92+169*v2)+u2*(-9+(-42+81*c2)*v2)+10*(-1+9*v2+10*v4)))/rho15)

  fu(jind(11,9,2)) = (u8*v2*(9-2*u2+7*cuv+9*v2))/rho13
fuuu(jind(11,9,2)) = (3*u6*v2*(-8*u6+84*c*u5*v+216*cuv*(1+v2)**2+168*(1+v2)**3+5*c*u3*v*(-83+(-83+7*c2)*v2) &
                      +9*u2*(1+v2)*(-43+(-43+15*c2)*v2)-4*u4*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,9,2)) = (u9*v*(2+2*u2-7*cuv-9*v2))/rho13
fvvv(jind(11,9,2)) = (-33*u9*(21*c**3*u3*v2-9*c2*u2*v*(2+2*u2-7*v2)+c*u*(2+2*u4-45*v2+70*v4+u2*(4-45*v2)) &
                      +v*(4+4*u4-31*v2+30*v4+u2*(8-31*v2))))/rho17
 fuv(jind(11,9,2)) = -((u8*v*(4*u4-39*c*u3*v+5*cuv*(12+23*v2)+u2*(-14+(-58+49*c2)*v2)+9*(-2+7*v2+9*v4)))/rho15)

  fu(jind(11,8,3)) = (u7*v3*(8-3*u2+5*cuv+8*v2))/rho13
fuuu(jind(11,8,3)) = (3*u5*v3*(-20*u6+100*c*u5*v+56*cuv*(1+v2)**2+112*(1+v2)**3+u4*(215+(215-75*c2)*v2) &
                      +8*u2*(1+v2)*(-46+(-46+3*c2)*v2)+c*u3*v*(-273+(-273+5*c2)*v2)))/rho17
  fv(jind(11,8,3)) = (u8*v2*(3+3*u2-5*cuv-8*v2))/rho13
fvvv(jind(11,8,3)) = (3*u8*(2+2*u6-54*c*u5*v-93*v2+380*v4-240*v6+c*u3*v*(-108-7*(-71+15*c2)*v2)+3*u4*(2+(-31+63*c2)*v2) &
                      +cuv*(-54+497*v2-450*v4)+u2*(6+3*(-62+63*c2)*v2+(380-350*c2)*v4)))/rho17
 fuv(jind(11,8,3)) = -((u7*v2*(9*u4-41*c*u3*v+3*cuv*(12+23*v2)+5*u2*(-3+(-14+5*c2)*v2)+8*(-3+5*v2+8*v4)))/rho15)

  fu(jind(11,7,4)) = (u6*v4*(7-4*u2+3*cuv+7*v2))/rho13
fuuu(jind(11,7,4)) = (-3*u4*v4*(40*u6-90*c*u5*v+42*cuv*(1+v2)**2-70*(1+v2)**3+7*u2*(1+v2)*(47+(47+c2)*v2) &
                      +c*u3*v*(95+(95+c2)*v2)+6*u4*(-46+(-46+3*c2)*v2)))/rho17
  fv(jind(11,7,4)) = (u7*v3*(4+4*u2-3*cuv-7*v2))/rho13
fvvv(jind(11,7,4)) = (3*u7*v*(8+8*u6-84*c*u5*v-152*v2+387*v4-168*v6+c*u3*v*(-168-5*(-83+7*c2)*v2)+4*u4*(6+(-38+35*c2)*v2) &
                      +cuv*(-84+415*v2-216*v4)+u2*(24+4*(-76+35*c2)*v2-9*(-43+15*c2)*v4)))/rho17
 fuv(jind(11,7,4)) = -((u6*v3*(16*u4-35*c*u3*v+cuv*(20+31*v2)+3*u2*(-4+(-26+3*c2)*v2)+7*(-4+3*v2+7*v4)))/rho15)

  fu(jind(11,6,5)) = (u5*v5*(6-5*u2+cuv+6*v2))/rho13
fuuu(jind(11,6,5)) = (3*u3*v5*(-70*u6+42*c*u5*v-90*cuv*(1+v2)**2+40*(1+v2)**3+7*u4*(47+(47+c2)*v2) &
                      +c*u3*v*(95+(95+c2)*v2)+6*u2*(1+v2)*(-46+(-46+3*c2)*v2)))/rho17
  fv(jind(11,6,5)) = (u6*v4*(5+5*u2-cuv-6*v2))/rho13
fvvv(jind(11,6,5)) = (3*u6*v2*(20+20*u6-100*c*u5*v-215*v2+368*v4-112*v6+c*u3*v*(-200+(273-5*c2)*v2)+5*u4*(12+(-43+15*c2)*v2) &
                      +cuv*(-100+273*v2-56*v4)+u2*(60+5*(-86+15*c2)*v2-8*(-46+3*c2)*v4)))/rho17
 fuv(jind(11,6,5)) = -((u5*v4*(25*u4-21*c*u3*v+cuv*(12+v2)+u2*(-5+(-82+c2)*v2)+6*(-5+v2+6*v4)))/rho15)

  fu(jind(11,5,6)) = (u4*v6*(5-6*u2-cuv+5*v2))/rho13
fuuu(jind(11,5,6)) = (3*u2*v6*(-112*u6-56*c*u5*v-100*cuv*(1+v2)**2+20*(1+v2)**3+c*u3*v*(273+(273-5*c2)*v2) &
                      -8*u4*(-46+(-46+3*c2)*v2)+5*u2*(1+v2)*(-43+(-43+15*c2)*v2)))/rho17
  fv(jind(11,5,6)) = (u5*v5*(6+6*u2+cuv-5*v2))/rho13
fvvv(jind(11,5,6)) = (3*u5*v3*(40+40*u6-90*c*u5*v-276*v2+329*v4-70*v6+c*u3*v*(-180+(95+c2)*v2)+6*u4*(20+(-46+3*c2)*v2) &
                       +cuv*(-90+95*v2+42*v4)+u2*(120+6*(-92+3*c2)*v2+7*(47+c2)*v4)))/rho17
 fuv(jind(11,5,6)) = -((u4*v5*(36*u4+c*u3*v-3*cuv*(-4+7*v2)+u2*(6+(-82+c2)*v2)+5*(-6-v2+5*v4)))/rho15)

  fu(jind(11,4,7)) = (u3*v7*(4-7*u2-3*cuv+4*v2))/rho13
fuuu(jind(11,4,7)) = (-3*u*v7*(168*u6+216*c*u5*v+84*cuv*(1+v2)**2-8*(1+v2)**3+5*c*u3*v*(-83+(-83+7*c2)*v2) &
                       +9*u4*(-43+(-43+15*c2)*v2)-4*u2*(1+v2)*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,4,7)) = (u4*v6*(7+7*u2+3*cuv-4*v2))/rho13
fvvv(jind(11,4,7)) = (-3*u4*v4*(-70-70*u6+42*c*u5*v+329*v2-276*v4+40*v6+7*u4*(-30+(47+c2)*v2)+c*u3*v*(84+(95+c2)*v2) &
                       +cuv*(42+95*v2-90*v4)+u2*(-210+7*(94+c2)*v2+6*(-46+3*c2)*v4)))/rho17
 fuv(jind(11,4,7)) = -((u3*v6*(49*u4+31*c*u3*v-5*cuv*(-4+7*v2)+3*u2*(7+(-26+3*c2)*v2)+4*(-7-3*v2+4*v4)))/rho15)

  fu(jind(11,3,8)) = (u2*v8*(3-8*u2-5*cuv+3*v2))/rho13
fuuu(jind(11,3,8)) = (3*v8*(-240*u6-450*c*u5*v-54*cuv*(1+v2)**2+2*(1+v2)**3+u4*(380+(380-350*c2)*v2) &
                       -7*c*u3*v*(-71+(-71+15*c2)*v2)+3*u2*(1+v2)*(-31+(-31+63*c2)*v2)))/rho17
  fv(jind(11,3,8)) = (u3*v7*(8+8*u2+5*cuv-3*v2))/rho13
fvvv(jind(11,3,8)) = (3*u3*v5*(112+112*u6+56*c*u5*v-368*v2+215*v4-20*v6+8*u4*(42+(-46+3*c2)*v2)+c*u3*v*(112+(-273+5*c2)*v2) &
                       +cuv*(56-273*v2+100*v4)+u2*(336+8*(-92+3*c2)*v2+(215-75*c2)*v4)))/rho17
 fuv(jind(11,3,8)) = -((u2*v7*(64*u4+69*c*u3*v+cuv*(36-41*v2)+5*u2*(8+(-14+5*c2)*v2)+3*(-8-5*v2+3*v4)))/rho15)

  fu(jind(11,2,9)) = (u*v9*(2-9*u2-7*cuv+2*v2))/rho13
fuuu(jind(11,2,9)) = (-33*v9*(30*u5+70*c*u4*v+2*c*v*(1+v2)**2+3*c*u2*v*(-15+(-15+7*c2)*v2)-2*u*(1+v2)*(-2+(-2+9*c2)*v2) &
                       +u3*(-31+(-31+63*c2)*v2)))/rho17
  fv(jind(11,2,9)) = (u2*v8*(9+9*u2+7*cuv-2*v2))/rho13
fvvv(jind(11,2,9)) = (3*u2*v6*(168+168*u6+216*c*u5*v-387*v2+152*v4-8*v6+c*u3*v*(432+5*(-83+7*c2)*v2)+9*u4*(56+(-43+15*c2)*v2) &
                       +cuv*(216-415*v2+84*v4)+u2*(504+9*(-86+15*c2)*v2-4*(-38+35*c2)*v4)))/rho17
 fuv(jind(11,2,9)) = -((u*v8*(81*u4+115*c*u3*v+u2*(63+(-58+49*c2)*v2)+c*u*(60*v-39*v3)+2*(-9-7*v2+2*v4)))/rho15)

  fu(jind(11,1,10)) = (v10*(1-10*u2-9*cuv+v2))/rho13
fuuu(jind(11,1,10)) = (33*v10*(-40*u4-108*c*u3*v+cuv*(35+(35-39*c2)*v2)-12*u2*(-2+(-2+9*c2)*v2)+(1+v2)*(-1+(-1+13*c2)*v2)))/rho17
  fv(jind(11,1,10)) = (u*v9*(10+10*u2+9*cuv-v2))/rho13
fvvv(jind(11,1,10)) = (3*u*v7*(240+240*u6+450*c*u5*v-380*v2+93*v4-2*v6+c*u3*v*(900+7*(-71+15*c2)*v2)+10*u4*(72+(-38+35*c2)*v2) &
                       +cuv*(450-497*v2+54*v4)+u2*(720+10*(-76+35*c2)*v2+(93-189*c2)*v4)))/rho17
 fuv(jind(11,1,10)) = -((v9*(-10+100*u4+169*c*u3*v-9*v2+v4+cuv*(92-29*v2)+u2*(90+(-42+81*c2)*v2)))/rho15)

  fu(jind(11,0,11)) = (-11*v11*(u+c*v))/rho13
fuuu(jind(11,0,11)) = (-429*v11*(u+c*v)*(-1+4*u2+8*cuv+(-1+5*c2)*v2))/rho17
  fv(jind(11,0,11)) = (11*v10*(1+u2+cuv))/rho13
fvvv(jind(11,0,11)) = (33*v8*(30+30*u6+70*c*u5*v-31*v2+4*v4+c*u3*v*(140+3*(-15+7*c2)*v2)+u4*(90+(-31+63*c2)*v2) &
                       +cuv*(70-45*v2+2*v4)+u2*(90+(-62+63*c2)*v2+(4-18*c2)*v4)))/rho17
 fuv(jind(11,0,11)) = (-11*v10*(11*u3+21*c*u2*v-c*v*(-12+v2)+u*(11+(-2+11*c2)*v2)))/rho15

  fu(jind(11,12,0)) = (u11*(u2+13*cuv+12*(1+v2)))/rho13
fuuu(jind(11,12,0)) = (33*u9*(108*cuv*(1+v2)**2+40*(1+v2)**3+u4*(1+(1-13*c2)*v2)+12*u2*(1+v2)*(-2+(-2+9*c2)*v2) &
                       +c*u3*v*(-35+(-35+39*c2)*v2)))/rho17
  fv(jind(11,12,0)) = (-11*u12*(c*u+v))/rho13
fvvv(jind(11,12,0)) = (-429*u12*(c*u+v)*(-1+(-1+5*c2)*u2+8*cuv+4*v2))/rho17
 fuv(jind(11,12,0)) = (-11*u11*(-(u2*v)+13*c2*u2*v+12*v*(1+v2)+c*u*(13+24*v2)))/rho15

  fu(jind(11,11,1)) = (11*u10*v*(1 + cuv + v2))/rho13
fuuu(jind(11,11,1)) = (33*u8*v*(2*c*u5*v + 70*cuv*(1 + v2)**2 + 30*(1 + v2)**3 + u4*(4 + (4 - 18*c2)*v2)  &
                       + 3*c*u3*v*(-15 + (-15 + 7*c2)*v2) + u2*(1 + v2)*(-31 + (-31 + 63*c2)*v2)))/rho17 
  fv(jind(11,11,1)) = (u11*(1 + u2 - 9*cuv - 10*v2))/rho13
fvvv(jind(11,11,1)) = (33*u11*(-1 + (-1 + 13*c2)*u4 + c*(35 - 39*c2)*u3*v + 24*v2 - 40*v4 + cuv*(35 - 108*v2) &
                        + u2*(-2 + 24*v2 + c2*(13 - 108*v2))))/rho17
 fuv(jind(11,11,1)) = (11*u10*(1 + 2*c*u3*v - 9*v2 - 10*v4 - 9*cuv*(1 + 2*v2) + u2*(1 + (3 - 9*c2)*v2)))/rho15

  fu(jind(11,10,2)) = (u9*v2*(-u2+9*cuv+10*(1+v2)))/rho13
fuuu(jind(11,10,2)) = (3*u7*v2*(-2*u6+54*c*u5*v+450*cuv*(1+v2)**2+240*(1+v2)**3+u4*(93+(93-189*c2)*v2) &
                      +7*c*u3*v*(-71+(-71+15*c2)*v2)+10*u2*(1+v2)*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,10,2)) = (u10*v*(2+2*u2-7*cuv-9*v2))/rho13
fvvv(jind(11,10,2)) = (-33*u10*(21*c**3*u3*v2-9*c2*u2*v*(2+2*u2-7*v2)+c*u*(2+2*u4-45*v2+70*v4+u2*(4-45*v2)) &
                      +v*(4+4*u4-31*v2+30*v4+u2*(8-31*v2))))/rho17
 fuv(jind(11,10,2)) = -((u9*v*(2*u4-36*c*u3*v+7*cuv*(9+20*v2)+3*u2*(-6+(-17+21*c2)*v2)+10*(-2+7*v2+9*v4)))/rho15)

  fu(jind(11,9,3)) = (u8*v3*(9-2*u2+7*cuv+9*v2))/rho13
fuuu(jind(11,9,3)) = (3*u6*v3*(-8*u6+84*c*u5*v+216*cuv*(1+v2)**2+168*(1+v2)**3+5*c*u3*v*(-83+(-83+7*c2)*v2) &
                     +9*u2*(1+v2)*(-43+(-43+15*c2)*v2)-4*u4*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,9,3)) = (u9*v2*(3+3*u2-5*cuv-8*v2))/rho13
fvvv(jind(11,9,3)) = (3*u9*(2+2*u6-54*c*u5*v-93*v2+380*v4-240*v6+c*u3*v*(-108-7*(-71+15*c2)*v2)+3*u4*(2+(-31+63*c2)*v2) &
                     +cuv*(-54+497*v2-450*v4)+u2*(6+3*(-62+63*c2)*v2+(380-350*c2)*v4)))/rho17
 fuv(jind(11,9,3)) = -((u8*v2*(6*u4-42*c*u3*v+5*cuv*(7+18*v2)+u2*(-21+5*(-13+7*c2)*v2)+9*(-3+5*v2+8*v4)))/rho15)

fu(jind(11,8,4)) = (u7*v4*(8-3*u2+5*cuv+8*v2))/rho13
fuuu(jind(11,8,4)) = (3*u5*v4*(-20*u6+100*c*u5*v+56*cuv*(1+v2)**2+112*(1+v2)**3+u4*(215+(215-75*c2)*v2) &
                     +8*u2*(1+v2)*(-46+(-46+3*c2)*v2)+c*u3*v*(-273+(-273+5*c2)*v2)))/rho17
  fv(jind(11,8,4)) = (u8*v3*(4+4*u2-3*cuv-7*v2))/rho13
fvvv(jind(11,8,4)) = (3*u8*v*(8+8*u6-84*c*u5*v-152*v2+387*v4-168*v6+c*u3*v*(-168-5*(-83+7*c2)*v2)+4*u4*(6+(-38+35*c2)*v2) &
                     +cuv*(-84+415*v2-216*v4)+u2*(24+4*(-76+35*c2)*v2-9*(-43+15*c2)*v4)))/rho17
 fuv(jind(11,8,4)) = -((u7*v3*(12*u4-40*c*u3*v+3*cuv*(5+16*v2)+5*u2*(-4+3*(-5+c2)*v2)+8*(-4+3*v2+7*v4)))/rho15)

  fu(jind(11,7,5)) = (u6*v5*(7-4*u2+3*cuv+7*v2))/rho13
fuuu(jind(11,7,5)) = (-3*u4*v5*(40*u6-90*c*u5*v+42*cuv*(1+v2)**2-70*(1+v2)**3+7*u2*(1+v2)*(47+(47+c2)*v2) &
                     +c*u3*v*(95+(95+c2)*v2)+6*u4*(-46+(-46+3*c2)*v2)))/rho17
  fv(jind(11,7,5)) = (u7*v4*(5+5*u2-cuv-6*v2))/rho13
fvvv(jind(11,7,5)) = (3*u7*v2*(20+20*u6-100*c*u5*v-215*v2+368*v4-112*v6+c*u3*v*(-200+(273-5*c2)*v2)+5*u4*(12+(-43+15*c2)*v2) &
                     +cuv*(-100+273*v2-56*v4)+u2*(60+5*(-86+15*c2)*v2-8*(-46+3*c2)*v4)))/rho17
 fuv(jind(11,7,5)) = -((u6*v4*(20*u4-30*c*u3*v+cuv*(3+14*v2)+3*u2*(-5+(-27+c2)*v2)+7*(-5+v2+6*v4)))/rho15)

  fu(jind(11,6,6)) = (u5*v6*(6-5*u2+cuv+6*v2))/rho13
fuuu(jind(11,6,6)) = (3*u3*v6*(-70*u6+42*c*u5*v-90*cuv*(1+v2)**2+40*(1+v2)**3+7*u4*(47+(47+c2)*v2)+c*u3*v*(95+(95+c2)*v2) &
                      +6*u2*(1+v2)*(-46+(-46+3*c2)*v2)))/rho17
  fv(jind(11,6,6)) = (u6*v5*(6+6*u2+cuv-5*v2))/rho13
fvvv(jind(11,6,6)) = (3*u6*v3*(40+40*u6-90*c*u5*v-276*v2+329*v4-70*v6+c*u3*v*(-180+(95+c2)*v2)+6*u4*(20+(-46+3*c2)*v2) &
                      +cuv*(-90+95*v2+42*v4)+u2*(120+6*(-92+3*c2)*v2+7*(47+c2)*v4)))/rho17
 fuv(jind(11,6,6)) = (u5*v5*(-30*u4+12*c*u3*v+cuv*(1+12*v2)+u2*(6+(83+c2)*v2)+6*(6+v2-5*v4)))/rho15

  fu(jind(11,5,7))= (u4*v7*(5 - 6*u2 - cuv + 5*v2))/rho13
fuuu(jind(11,5,7)) = (3*u2*v7*(-112*u6 - 56*c*u5*v - 100*cuv*(1 + v2)**2 + 20*(1 + v2)**3 + c*u3*v*(273 + (273 - 5*c2)*v2) &
                     - 8*u4*(-46 + (-46 + 3*c2)*v2) + 5*u2*(1 + v2)*(-43 + (-43 + 15*c2)*v2)))/rho17
  fv(jind(11,5,7)) = (u5*v6*(7 + 7*u2 + 3*cuv - 4*v2))/rho13
fvvv(jind(11,5,7)) = (-3*u5*v4*(-70-70*u6+42*c*u5*v+329*v2-276*v4+40*v6+7*u4*(-30+(47+c2)*v2)+c*u3*v*(84+(95+c2)*v2) &
                      + cuv*(42 + 95*v2 - 90*v4) + u2*(-210 + 7*(94 + c2)*v2 + 6*(-46 + 3*c2)*v4)))/rho17
 fuv(jind(11,5,7)) = -(u4*v6*(42*u4 + 14*c*u3*v + u2*(7 + 3*(-27 + c2)*v2) + c*u*(3*v - 30*v3) + 5*(-7 - 3*v2 + 4*v4)))/rho15

  fu(jind(11,4,8)) = (u3*v8*(4-7*u2-3*cuv+4*v2))/rho13
fuuu(jind(11,4,8)) = (-3*u*v8*(168*u6+216*c*u5*v+84*cuv*(1+v2)**2-8*(1+v2)**3+5*c*u3*v*(-83+(-83+7*c2)*v2) &
                     +9*u4*(-43+(-43+15*c2)*v2)-4*u2*(1+v2)*(-38+(-38+35*c2)*v2)))/rho17
  fv(jind(11,4,8)) = (u4*v7*(8+8*u2+5*cuv-3*v2))/rho13
fvvv(jind(11,4,8)) = (3*u4*v5*(112+112*u6+56*c*u5*v-368*v2+215*v4-20*v6+8*u4*(42+(-46+3*c2)*v2)+c*u3*v*(112+(-273+5*c2)*v2) &
                     +cuv*(56-273*v2+100*v4)+u2*(336+8*(-92+3*c2)*v2+(215-75*c2)*v4)))/rho17
 fuv(jind(11,4,8)) = -((u3*v7*(56*u4+48*c*u3*v-5*cuv*(-3+8*v2)+3*u2*(8+5*(-5+c2)*v2)+4*(-8-5*v2+3*v4)))/rho15)

  fu(jind(11,3,9)) = (u2*v9*(3-8*u2-5*cuv+3*v2))/rho13
fuuu(jind(11,3,9)) = (3*v9*(-240*u6-450*c*u5*v-54*cuv*(1+v2)**2+2*(1+v2)**3+u4*(380+(380-350*c2)*v2) &
                     -7*c*u3*v*(-71+(-71+15*c2)*v2)+3*u2*(1+v2)*(-31+(-31+63*c2)*v2)))/rho17
  fv(jind(11,3,9)) = (u3*v8*(9+9*u2+7*cuv-2*v2))/rho13
fvvv(jind(11,3,9)) = (3*u3*v6*(168+168*u6+216*c*u5*v-387*v2+152*v4-8*v6+c*u3*v*(432+5*(-83+7*c2)*v2)+9*u4*(56+(-43+15*c2)*v2) &
                      +cuv*(216-415*v2+84*v4)+u2*(504+9*(-86+15*c2)*v2-4*(-38+35*c2)*v4)))/rho17
 fuv(jind(11,3,9)) = -((u2*v8*(72*u4+90*c*u3*v-7*cuv*(-5+6*v2)+5*u2*(9+(-13+7*c2)*v2)+3*(-9-7*v2+2*v4)))/rho15)

  fu(jind(11,2,10)) = (u*v10*(2-9*u2-7*cuv+2*v2))/rho13
fuuu(jind(11,2,10)) = (-33*v10*(30*u5+70*c*u4*v+2*c*v*(1+v2)**2+3*c*u2*v*(-15+(-15+7*c2)*v2)-2*u*(1+v2)*(-2+(-2+9*c2)*v2) &
                      +u3*(-31+(-31+63*c2)*v2)))/rho17
  fv(jind(11,2,10)) = (u2*v9*(10+10*u2+9*cuv-v2))/rho13
fvvv(jind(11,2,10)) = (3*u2*v7*(240+240*u6+450*c*u5*v-380*v2+93*v4-2*v6+c*u3*v*(900+7*(-71+15*c2)*v2)+10*u4*(72+(-38+35*c2)*v2) &
                      +cuv*(450-497*v2+54*v4)+u2*(720+10*(-76+35*c2)*v2+(93-189*c2)*v4)))/rho17
 fuv(jind(11,2,10)) = -((u*v9*(90*u4+140*c*u3*v-9*cuv*(-7+4*v2)+u2*(70+(-51+63*c2)*v2)+2*(-10-9*v2+v4)))/rho15)

  fu(jind(11,1,11)) = (v11*(1-10*u2-9*cuv+v2))/rho13
fuuu(jind(11,1,11)) = (33*v11*(-40*u4-108*c*u3*v+cuv*(35+(35-39*c2)*v2)-12*u2*(-2+(-2+9*c2)*v2)+(1+v2)*(-1+(-1+13*c2)*v2)))/rho17
  fv(jind(11,1,11)) = (11*u*v10*(1+u2+cuv))/rho13
fvvv(jind(11,1,11)) = (33*u*v8*(30+30*u6+70*c*u5*v-31*v2+4*v4+c*u3*v*(140+3*(-15+7*c2)*v2)+u4*(90+(-31+63*c2)*v2)+ &
                       cuv*(70-45*v2+2*v4)+u2*(90+(-62+63*c2)*v2+(4-18*c2)*v4)))/rho17
 fuv(jind(11,1,11)) = (11*v10*(1-10*u4-18*c*u3*v+v2+cuv*(-9+2*v2)+u2*(-9+(3-9*c2)*v2)))/rho15

  fu(jind(11,0,12)) = (-11*v12*(u+c*v))/rho13
fuuu(jind(11,0,12)) = (-429*v12*(u+c*v)*(-1+4*u2+8*cuv+(-1+5*c2)*v2))/rho17
  fv(jind(11,0,12)) = (v11*(12+12*u2+13*cuv+v2))/rho13
fvvv(jind(11,0,12)) = (33*v9*(40+40*u6+108*c*u5*v-24*v2+v4+cuv*(108-35*v2)+12*u4*(10+(-2+9*c2)*v2)+c*u3*v*(216+(-35+39*c2)*v2) &
                      +u2*(120+12*(-4+9*c2)*v2+(1-13*c2)*v4)))/rho17
 fuv(jind(11,0,12)) = (-11*v11*(12*u3+13*c*v+24*c*u2*v+u*(12+(-1+13*c2)*v2)))/rho15

return
END SUBROUTINE allder_c


SUBROUTINE allintFbackrecur2(u1,u2,v1,v2,c,Hint,rexp,p,q,nb)
! use backward vertical for 005, 185,095, then do lateral for each k
use params, ONLY : nbmax
implicit none
real*8, INTENT(in) :: u1,u2,v1,v2,c
integer, INTENT(in), dimension(nbmax) :: rexp,p,q
integer, INTENT(in) :: nb
real*8, dimension(nbmax), INTENT(OUT) :: Hint
!LOCAL
real*8, dimension(0:12,0:5) :: Fintt,Fintb,Gintr,Gintl
real*8, dimension(0:12,0:12,0:5) :: Iint,intex
integer :: i,j,m,n,k,r,rmin
real*8 :: rho11,rho12,rho21,rho22,minsq,mincsq,rtc, &
      f11,f12,f21,f22,Dsq1b,Dsq1t,Dsq2l,Dsq2r,top,bot,right,left, &
      fact(12),f22mf12,f21mf11,f22mf21,f12mf11,rightmleft,topmbot,top1mbot1,top2mbot2,&
top2,bot2,top1,bot1

rho11=sqrt(1+u1**2+2*c*u1*v1+v1**2)
rho12=sqrt(1+u1**2+2*c*u1*v2+v2**2)
rho21=sqrt(1+u2**2+2*c*u2*v1+v1**2)
rho22=sqrt(1+u2**2+2*c*u2*v2+v2**2)
minsq=1-c**2
rtc=sqrt(minsq)

fact(1)=u2-u1
do j=1,11
fact(j+1)=fact(j)*u2+(-u1)**(j+1)
enddo

Iint=0

!Initialize F(0,0),G(0,0)
f11= -log(-u1-c*v1+rho11)
f12= -log(-u1-c*v2+rho12)
f21= -log(-u2-c*v1+rho21)
f22= -log(-u2-c*v2+rho22)
!F(0,0)= -log(-u-c*v+rho)
Fintt(0,0)= f22-f12
Fintb(0,0)= f21-f11

f11= -log(-v1-c*u1+rho11)
f12= -log(-v2-c*u1+rho12)
f21= -log(-v1-c*u2+rho21)
f22= -log(-v2-c*u2+rho22)
Gintr(0,0)= f22-f21
Gintl(0,0)= f12-f11
!!G(0,0)= -log(-v-c*u+rho)

!Set F(m,k),G(m,k) see notes
!----------
Dsq1b=(1-c**2)*v1**2+1
Dsq1t=(1-c**2)*v2**2+1
Dsq2l=(1-c**2)*u1**2+1
Dsq2r=(1-c**2)*u2**2+1
do k=1,5
  !k=1,..,5: rpl=3,5,7,9,11
  r   = 2*k-1
  f11= (u1+c*v1)/rho11**r 
  f12= (u1+c*v2)/rho12**r 
  f21= (u2+c*v1)/rho21**r 
  f22= (u2+c*v2)/rho22**r 

  Fintt(0,k)=( (f22-f12) + 2*(k-1)*Fintt(0,k-1) )/(Dsq1t*r)
  Fintb(0,k)=( (f21-f11) + 2*(k-1)*Fintb(0,k-1) )/(Dsq1b*r)
!  F(0,k)=( (u+c*v)/rho**r + 2*(k-1)*F(0,k-1) )/(Dsq1*r)

  f11= (v1+c*u1)/rho11**r 
  f12= (v2+c*u1)/rho12**r 
  f21= (v1+c*u2)/rho21**r 
  f22= (v2+c*u2)/rho22**r 
  Gintr(0,k)=( (f22-f21) + 2*(k-1)*Gintr(0,k-1) )/(Dsq2r*r)
  Gintl(0,k)=( (f12-f11) + 2*(k-1)*Gintl(0,k-1) )/(Dsq2l*r)
!  G(0,k)=( (v+c*u)/rho**r + 2*(k-1)*G(0,k-1) )/(Dsq2*r)
enddo

!m=1
f11= rho11
f12= rho12
f21= rho21
f22= rho22
Fintt(1,0)=(f22-f12) - c*v2*Fintt(0,0)
Fintb(1,0)=(f21-f11) - c*v1*Fintb(0,0)
!F(1,0)= rho - c*v*F(0,0) 
Gintr(1,0)=(f22-f21) - c*u2*Gintr(0,0)
Gintl(1,0)=(f12-f11) - c*u1*Gintl(0,0)
!G(1,0)= rho - c*u*G(0,0)
do m=2,12
  f11= u1**(m-1)*rho11 
  f12= u1**(m-1)*rho12 
  f21= u2**(m-1)*rho21 
  f22= u2**(m-1)*rho22 
  f22mf12= f22-f12
  f21mf11= f21-f11
  Fintt(m,0)= ( (f22mf12)- (m-1)*(v2**2+1)*Fintt(m-2,0)-(2*m-1)*c*v2*Fintt(m-1,0) )/m
  Fintb(m,0)= ( (f21mf11)- (m-1)*(v1**2+1)*Fintb(m-2,0)-(2*m-1)*c*v1*Fintb(m-1,0) )/m
!  F(m,0)=( u**(m-1)*rho - (m-1)*(v**2+1)*F(m-2,0)-(2*m-1)*c*v*F(m-1,0) )/m

  f11= v1**(m-1)*rho11 
  f12= v2**(m-1)*rho12 
  f21= v1**(m-1)*rho21 
  f22= v2**(m-1)*rho22 
  f22mf21= f22-f21
  f12mf11= f12-f11
  Gintr(m,0)=( (f22mf21) - (m-1)*(u2**2+1)*Gintr(m-2,0)-(2*m-1)*c*u2*Gintr(m-1,0) )/m
  Gintl(m,0)=( (f12mf11) - (m-1)*(u1**2+1)*Gintl(m-2,0)-(2*m-1)*c*u1*Gintl(m-1,0) )/m
!  G(m,0)=( v**(m-1)*rho - (m-1)*(u**2+1)*G(m-2,0)-(2*m-1)*c*u*G(m-1,0) )/m
enddo

!m=1
do k=1,5
  r   = 2*k-1
  f11= -1/(rho11**r*r) 
  f12= -1/(rho12**r*r) 
  f21= -1/(rho21**r*r) 
  f22= -1/(rho22**r*r) 
  Fintt(1,k)= (f22-f12)- c*v2*Fintt(0,k)
  Fintb(1,k)= (f21-f11)- c*v1*Fintb(0,k)
!  F(1,k)= -1/(rho**r*r) - c*v*F(0,k)

  Gintr(1,k)= (f22-f21) - c*u2*Gintr(0,k)
  Gintl(1,k)= (f12-f11) - c*u1*Gintl(0,k)
!  G(1,k)= -1/(rho**r*r) - c*u*G(0,k)
enddo

do m=2,12
do k=1,5
  r   = 2*k-1
  f11= -u1**(m-1)/rho11**r 
  f12= -u1**(m-1)/rho12**r 
  f21= -u2**(m-1)/rho21**r 
  f22= -u2**(m-1)/rho22**r 
  f22mf12= f22-f12
  f21mf11= f21-f11
  Fintt(m,k)=((f22mf12)+ (m-1)*Fintt(m-2,k-1) )/r - c*v2*Fintt(m-1,k)
  Fintb(m,k)=((f21mf11)+ (m-1)*Fintb(m-2,k-1) )/r - c*v1*Fintb(m-1,k)
!  F(m,k)=( -u**(m-1)/rho**r + (m-1)*F(m-2,k-1) )/r - c*v*F(m-1,k)

  f11= -v1**(m-1)/rho11**r 
  f12= -v2**(m-1)/rho12**r 
  f21= -v1**(m-1)/rho21**r 
  f22= -v2**(m-1)/rho22**r 
  f22mf21= f22-f21
  f12mf11= f12-f11
  Gintr(m,k)=( (f22mf21) + (m-1)*Gintr(m-2,k-1) )/r - c*u2*Gintr(m-1,k)
  Gintl(m,k)=( (f12mf11) + (m-1)*Gintl(m-2,k-1) )/r - c*u1*Gintl(m-1,k)
!  G(m,k)=( -v**(m-1)/rho**r + (m-1)*G(m-2,k-1) )/r - c*u*G(m-1,k)
enddo
enddo

!Checked all Gint,Fint for sample case in /bld/ellip/testresults, Gintr has 10 digits
!all others are great
!do m=0,12
!do k=0,5
!write(*,'(2i3,2e21.12,e12.2)')m,k,Gintr(m,k) !,Gintr(m,k)
!enddo
!enddo

! STEP 0: INITIALIZATION 
call compI005(u1,u2,v1,v2,c,Iint(0,0,5))
!print*,Iint(0,0,5)
!call compI095(u1,u2,v1,v2,c,Iint(0,9,5))
!call compI185(u1,u2,v1,v2,c,Iint(1,8,5))

!ONE: ALL I00k.
! Find I00k using Eq (9) backwards. Start with I005 found numerically. 
! Backwards the values are decreasing: stable. Forwards is unstable
do k=4,0,-1
    right =  u2*Gintr(0,k)
    left  =  u1*Gintl(0,k)
    top   =  v2*Fintt(0,k) 
    bot   =  v1*Fintb(0,k) 
    Iint(0,0,k) = ( (2*k+1)*Iint(0,0,k+1) - (right-left) - (top-bot) )/(2*k-1)
enddo

!!TWO: For all k+1, use Eq (9) when m+n+1=2k. Reduces to only Fs and Gs
!do k=1,4
!do m=0,2*k-1
!   n=2*k-m-1
!   right =  u2**(m+1)*Gintr(n,k)
!   left  =  u1**(m+1)*Gintl(n,k)
!   top   =  v2**(n+1)*Fintt(m,k) 
!   bot   =  v1**(n+1)*Fintb(m,k) 
!   Iint(m,n,k+1) = (  (right-left) + (top-bot) )/(2*k+1)
!enddo
!enddo
mincsq=1-c**2


!TWO: ZERO R. USE (10) and (11) to get all 2k-m-n-1=0
!Note: All values of order 1/r^2  where r=|x| or |y| 
!That is, all approx same size, so stable forward (or backward)

!For n=0, all k, m=2k-n-1 use (11) forward in k,m PERFECT!!
k=0; n=0; m=1
top   =  c*Fintt(0,k) 
bot   =  c*Fintb(0,k) 
right =  Gintr(0,k)
left  =  Gintl(0,k)
Iint(m,n,k+1) = ((top-bot)-(right-left))/mincsq  
do k=1,4
  m=2*k-n+1
  top   =  c*v2**n*Fintt(m-1,k) 
  bot   =  c*v1**n*Fintb(m-1,k) 
  right =  u2**(m-1)*Gintr(n,k)
  left  =  u1**(m-1)*Gintl(n,k)
  Iint(m,n,k+1) = ((m-1)*Iint(m-2,n,k)+(top-bot)-(right-left))/(2*k+1)/mincsq
enddo
!print*,Iint(1,0,5)
!stop

!For m=0, all k, n=2k-m-1, use (10) forward in k,n  PERFECT!!
k=0; m=0; n=1; 
top   =  Fintt(0,k) 
bot   =  Fintb(0,k) 
right =  Gintr(0,k)
left  =  Gintl(0,k)
Iint(m,n,k+1) = (c*(right-left)-(top-bot))/mincsq  
do k=1,4
  n=2*k-m+1
  top   =  v2**(n-1)*Fintt(m,k) 
  bot   =  v1**(n-1)*Fintb(m,k) 
  right =  u2**m*Gintr(n-1,k)
  left  =  u1**m*Gintl(n-1,k)
  Iint(m,n,k+1) = ((n-1)*Iint(m,n-2,k)+c*(right-left)-(top-bot))/(2*k+1)/mincsq
enddo

!For m>=1, all k, n=2k-m-1, use (10) forward in k,n PERFECT, see testnew.m
do k=0,4
do m=1,12
  n=2*k-m+1
  if (n.ge.2) then
    top   =  v2**(n-1)*Fintt(m,k) 
    bot   =  v1**(n-1)*Fintb(m,k) 
    right =  u2**m*Gintr(n-1,k)
    left  =  u1**m*Gintl(n-1,k)
    Iint(m,n,k+1) = ((n-1)*Iint(m,n-2,k)-c*m*Iint(m-1,n-1,k)+c*(right-left)-(top-bot))/(2*k+1)/mincsq
  else if (n.eq.1) then
    top   =  v2**(n-1)*Fintt(m,k) 
    bot   =  v1**(n-1)*Fintb(m,k) 
    right =  u2**m*Gintr(n-1,k)
    left  =  u1**m*Gintl(n-1,k)
    Iint(m,n,k+1) = (-c*m*Iint(m-1,n-1,k)+c*(right-left)-(top-bot))/(2*k+1)/mincsq
  endif
enddo
enddo

!THREE: NONZERO R. USE (12) and (13) lateral, increasing!, to get all 2k-m-n-1.ne.0  
!For n=0, m>=1, use (13) 
n=0
do k=0,5
  m=1
  r=2*k-n-m-1
  if (r.ne.0) then
    top   =  c*Fintt(m-1,k) - mincsq*v2*Fintt(m,k)
    bot   =  c*Fintb(m-1,k) - mincsq*v1*Fintb(m,k)
    right =  (1 + mincsq*u2**2)*Gintr(n,k)
    left  =  (1 + mincsq*u1**2)*Gintl(n,k)
    Iint(m,n,k) = ((top-bot)-(right-left))/(r*mincsq)
  endif
  do m=2,12
    r=2*k-n-m-1
    if (r.ne.0) then
      top   =  c*Fintt(m-1,k) - mincsq*v2*Fintt(m,k)
      bot   =  c*Fintb(m-1,k) - mincsq*v1*Fintb(m,k)
      right =  (u2**(m-1) + mincsq*u2**(m+1))*Gintr(n,k)
      left  =  (u1**(m-1) + mincsq*u1**(m+1))*Gintl(n,k)
      Iint(m,n,k) = ((m-1)*Iint(m-2,n,k)+(top-bot)-(right-left))/(mincsq*r)
    endif
  enddo
enddo
!For m=0, n>=1, use (12) 
m=0
do k=0,5
  n=1
  r=2*k-n-m-1
  if (r.ne.0) then
    right =  c*Gintr(n-1,k) - mincsq*u2*Gintr(n,k)
    left  =  c*Gintl(n-1,k) - mincsq*u1*Gintl(n,k)
    top   =  (1 + mincsq*v2**2)*Fintt(m,k)
    bot   =  (1 + mincsq*v1**2)*Fintb(m,k)
    Iint(m,n,k) = (-(top-bot)+(right-left))/(r*mincsq)
  endif
  do n=2,12
    r=2*k-n-m-1
    if (r.ne.0) then
      right=  c*Gintr(n-1,k) - mincsq*u2*Gintr(n,k)
      left =  c*Gintl(n-1,k) - mincsq*u1*Gintl(n,k)
      top  =  (v2**(n-1) + mincsq*v2**(n+1))*Fintt(m,k)
      bot  =  (v1**(n-1) + mincsq*v1**(n+1))*Fintb(m,k)
      Iint(m,n,k) = ((n-1)*Iint(m,n-2,k)-(top-bot)+(right-left))/(mincsq*r)
    endif
  enddo
enddo

!For m>=1,n>=1 use (12) 
do k=0,5
do m=1,12
  n=1
  r=2*k-n-m-1
  if (r.ne.0) then
    right =  c*u2**m*Gintr(n-1,k) - mincsq*u2**(m+1)*Gintr(n,k)
    left  =  c*u1**m*Gintl(n-1,k) - mincsq*u1**(m+1)*Gintl(n,k)
    top   =  (1 + mincsq*v2**2)*Fintt(m,k)
    bot   =  (1 + mincsq*v1**2)*Fintb(m,k)
    Iint(m,n,k) =  (-c*m*Iint(m-1,n-1,k)-(top-bot)+(right-left))/(r*mincsq)
!     if (m.lt.10) then
!      write(*,'(a,i1,i1,i1,e25.17)')'I',m,n,k,Iint(m,n,k)
!     else
!      write(*,'(a,i2,i1,i1,e25.17)')'I',m,n,k,Iint(m,n,k)
!     endif
  endif
  do n=2,12
    r=2*k-n-m-1
    if (r.ne.0) then
      right=  c*u2**m*Gintr(n-1,k) - mincsq*u2**(m+1)*Gintr(n,k)
      left =  c*u1**m*Gintl(n-1,k) - mincsq*u1**(m+1)*Gintl(n,k)
      top  =  (v2**(n-1) + mincsq*v2**(n+1))*Fintt(m,k)
      bot  =  (v1**(n-1) + mincsq*v1**(n+1))*Fintb(m,k)
      Iint(m,n,k) = ((n-1)*Iint(m,n-2,k)-c*m*Iint(m-1,n-1,k)-(top-bot)+(right-left))/(mincsq*r)
!     if ((m.lt.10).and.(n.lt.10)) then
!      write(*,'(a,i1,i1,i1,e25.17)')'I',m,n,k,Iint(m,n,k)
!     else if ((m.lt.10).and.(n.ge.10)) then
!      write(*,'(a,i1,i2,i1,e25.17)')'I',m,n,k,Iint(m,n,k)
!     else if ((m.ge.10).and.(n.lt.10)) then
!      write(*,'(a,i2,i1,i1,e25.17)')'I',m,n,k,Iint(m,n,k)
!     else
!      write(*,'(a,i2,i2,i1,e25.17)')'I',m,n,k,Iint(m,n,k)
!     endif
    endif
  enddo
enddo
enddo

do j=1,nb
  k=(rexp(j)-1)/2
  Hint(j)=Iint(p(j),q(j),k)
!print*,j,p(j),q(j),k
enddo
END SUBROUTINE allintFbackrecur2


SUBROUTINE compI005(u1,u2,v1,v2,cc,Iint)
! taken from i005/findNgI005
! if u1>Rb, use Rb=infty value, done
! compute the error between outer square S_Rb and inner square S_a using Gauss Quad with Ng points. Then determine Ng for which error < 1.e-13
! confirmed that numbers are accurate using samplei005.m
implicit none
real*8, parameter :: pi= 3.14159265358979324d0,pi2=2*pi
real*8,INTENT(IN) :: u1,u2,v1,v2,cc
real*8,INTENT(OUT) :: Iint
!LOCAL
integer :: j,Ng
real*8, dimension(128) :: x,y,w,f
real*8 :: Rb,diagmax,umin,sum1,sum2,sum3,sum4,Iinf,uleft,uright,vbot,vtop

!Rb=20.7/sqrt(1-cc)
Rb=20.7/sqrt(1-abs(cc)**1.4)**0.56
! For now we will use a 1-2-3 ellipse for which c<=0.4. Set Ng for n=40/80
!Ng=40
Ng=48
Iinf=pi2/9/sqrt(1-cc**2)

umin=min(abs(u1),u2,abs(v1),v2)
diagmax=max(sqrt(u1**2+v1**2), sqrt(u2**2+v2**2))
if (umin.ge.Rb) then
  Iint=Iinf
else if (diagmax.le.Rb) then
! S_a completely contained in S_Rb!!!
  !integrate over R1
  sum1=0
  call setgauss(-Rb,v1,Ng,y,w)
  do j=1,Ng
    call intf005dx(u1,Rb,y(j),cc,f(j))
    sum1=sum1+f(j)*w(j)
  enddo
  !integrate over R2
  sum2=0
  call setgauss(u2,Rb,Ng,x,w)
  do j=1,Ng
    call intf005dy(v1,Rb,x(j),cc,f(j))
    sum2=sum2+f(j)*w(j)
  enddo
  !integrate over R3
  sum3=0
  call setgauss(v2,Rb,Ng,y,w)
  do j=1,Ng
    call intf005dx(-Rb,u2,y(j),cc,f(j))
    sum3=sum3+f(j)*w(j)
  enddo
  !integrate over R4
  sum4=0
  call setgauss(-Rb,u1,Ng,x,w)
  do j=1,Ng
    call intf005dy(-Rb,v2,x(j),cc,f(j))
    sum4=sum4+f(j)*w(j)
  enddo
  Iint=Iinf-(sum1+sum2+sum3+sum4)
else 
! stop ' S_a not inside nor outside B_r'
! Find bounds of the big box
  uleft=min(u1,-Rb)
  uright=max(u2,Rb)
  vbot=min(v1,-Rb)
  vtop=max(v2,Rb)
  sum1=0
  sum2=0
  sum3=0
  sum4=0
!  Ng=8
  Ng=16
! Integrate over same 4 rectangles as above, provided they are nonzero
  if ((v1-vbot).ge.Rb/20) then
    !integrate over R1
    !call setgauss(-Rb,v1,Ng,y,w)
     call setgauss(vbot,v1,Ng,y,w)
     do j=1,Ng
       !call intf005dx(u1,Rb,y(j),cc,f(j))
       call intf005dx(u1,uright,y(j),cc,f(j))
       sum1=sum1+f(j)*w(j)
     enddo
  endif
  if ((uright-u2).ge.Rb/20) then
    !integrate over R2
    !call setgauss(u2,Rb,Ng,x,w)
     call setgauss(u2,uright,Ng,x,w)
     do j=1,Ng
       !call intf005dy(v1,Rb,x(j),cc,f(j))
       call intf005dy(v1,vtop,x(j),cc,f(j))
       sum2=sum2+f(j)*w(j)
     enddo
  endif
  if ((vtop-v2).ge.Rb/20) then
    !integrate over R3
    !call setgauss(v2,Rb,Ng,y,w)
    call setgauss(v2,vtop,Ng,y,w)
    do j=1,Ng
      !call intf005dx(-Rb,u2,y(j),cc,f(j))
      call intf005dx(uleft,u2,y(j),cc,f(j))
      sum3=sum3+f(j)*w(j)
    enddo
  endif
  if ((u1-uleft).ge.Rb/20) then
    !integrate over R4
    !call setgauss(-Rb,u1,Ng,x,w)
    call setgauss(uleft,u1,Ng,x,w)
    do j=1,Ng
      !call intf005dy(-Rb,v2,x(j),cc,f(j))
      call intf005dy(vbot,v2,x(j),cc,f(j))
      sum4=sum4+f(j)*w(j)
    enddo
  endif
  Iint=Iinf-(sum1+sum2+sum3+sum4)
endif
END SUBROUTINE compI005


SUBROUTINE intf005dx(a,b,y,c,intfdx)
! int_a^b f_005 dx
implicit none
real*8, INTENT(IN) :: a,b,y,c
real*8, INTENT(OUT) :: intfdx
!LOCAL
real*8 :: x,rhosq,rho,p,q,num,den,fa,fb,q2,q3,q4

q = -1-(1-c**2)*y**2
q2 = q*q
q3 = q2*q
q4 = q3*q

x=b
rhosq=1+x**2+y**2+2*c*x*y; rho=sqrt(rhosq)
p = rhosq-(x+c*y)*rho 
num=5*q4/9 + 5*q3*p + 120*q2*p**2/7 + 80*q*p**3/3 + 16*p**4
den=p**9
fb=num/den

x=a
rhosq=1+x**2+y**2+2*c*x*y; rho=sqrt(rhosq)
p = rhosq-(x+c*y)*rho 
num=5*q4/9 + 5*q3*p + 120*q2*p**2/7 + 80*q*p**3/3 + 16*p**4
den=p**9
fa=num/den

intfdx=(fb-fa)/5
END SUBROUTINE intf005dx


SUBROUTINE intf005dy(a,b,x,c,intfdy)
! int_a^b f_005 dy
implicit none
real*8, INTENT(IN) :: a,b,x,c
real*8, INTENT(OUT) :: intfdy
!LOCAL
real*8 :: y,rhosq,rho,p,q,num,den,fa,fb,q2,q3,q4

q = -1-(1-c**2)*x**2
q2 = q*q
q3 = q2*q
q4 = q3*q

y=b
rhosq=1+x**2+y**2+2*c*x*y; rho=sqrt(rhosq)
p = rhosq-(y+c*x)*rho 
num=q4/9 + q3*p + 24*q2*p**2/7 + 16*q*p**3/3 + 16*p**4/5
den=p**9
fb=num/den

y=a
rhosq=1+x**2+y**2+2*c*x*y; rho=sqrt(rhosq)
p = rhosq-(y+c*x)*rho 
num=q4/9 + q3*p + 24*q2*p**2/7 + 16*q*p**3/3 + 16*p**4/5
!num=5*q4/9 + 5*q3*p + 120*q2*p**2/7 + 80*q*p**3/3 + 16*p**4
den=p**9

fa=num/den
intfdy=fb-fa
END SUBROUTINE intf005dy


SUBROUTINE setgauss(a,b,Ng,xg,wg)
! Assume S_a completely contained in S_Rb!!!
implicit none
integer, INTENT(IN) :: Ng
real*8, INTENT(IN) :: a,b
real*8, INTENT(OUT) :: xg(*),wg(*)

if (Ng.eq.2) then
  xg(1:2) = (/ -0.577350269189626d0,  &
           0.577350269189626d0 /)
  wg(1:2) = (/ 1.000000000000000d0,  &
          1.000000000000000d0 /)
  xg(1:2)=(xg(1:2)+1)/2*(b-a)+a; wg(1:2) = wg(1:2)/2*(b-a)
else if (Ng.eq.4) then
  xg(1:4) = (/ -0.861136311594053d0,  &
          -0.339981043584856d0,  &
           0.339981043584856d0,  &
           0.861136311594053d0 /)
  wg(1:4) = (/ 0.347854845137454d0,  &
          0.652145154862546d0,  &
          0.652145154862546d0,  &
          0.347854845137454d0 /)
  xg(1:4)=(xg(1:4)+1)/2*(b-a)+a; wg(1:4) = wg(1:4)/2*(b-a)
else if (Ng.eq.8) then
  xg(1:8) = (/ -0.960289856497536d0,  &
          -0.796666477413627d0,  &
          -0.525532409916329d0,  &
          -0.183434642495650d0,  &
           0.183434642495650d0,  &
           0.525532409916329d0,  &
           0.796666477413627d0,  &
           0.960289856497536d0 /)
  wg(1:8) = (/ 0.101228536290376d0,  &
          0.222381034453375d0,  &
          0.313706645877887d0,  &
          0.362683783378363d0,  &
          0.362683783378362d0,  &
          0.313706645877887d0,  &
          0.222381034453374d0,  &
          0.101228536290376d0 /)
  xg(1:8)=(xg(1:8)+1)/2*(b-a)+a; wg(1:8) = wg(1:8)/2*(b-a)
else if (Ng.eq.16) then
  xg(1:16) = (/ -0.989400934991650d0,  &
          -0.944575023073233d0,  &
          -0.865631202387832d0,  &
          -0.755404408355003d0,  &
          -0.617876244402644d0,  &
          -0.458016777657227d0,  &
          -0.281603550779258d0,  &
          -0.095012509837637d0,  &
           0.095012509837638d0,  &
           0.281603550779260d0,  &
           0.458016777657227d0,  &
           0.617876244402644d0,  &
           0.755404408355003d0,  &
           0.865631202387832d0,  &
           0.944575023073233d0,  &
           0.989400934991650d0 /)
  wg(1:16) = (/ 0.027152459411754d0,  &
          0.062253523938648d0,  &
          0.095158511682493d0,  &
          0.124628971255534d0,  &
          0.149595988816577d0,  &
          0.169156519395003d0,  &
          0.182603415044923d0,  &
          0.189450610455069d0,  &
          0.189450610455069d0,  &
          0.182603415044924d0,  &
          0.169156519395002d0,  &
          0.149595988816576d0,  &
          0.124628971255534d0,  &
          0.095158511682493d0,  &
          0.062253523938647d0,  &
          0.027152459411755d0 /)
  xg(1:16)=(xg(1:16)+1)/2*(b-a)+a; wg(1:16) = wg(1:16)/2*(b-a)
else if (Ng.eq.24) then
  xg(1:24) = (/ &
  -0.995187219997021d0, &
  -0.974728555971310d0, &
  -0.938274552002732d0, &
  -0.886415527004401d0, &
  -0.820001985973903d0, &
  -0.740124191578554d0, &
  -0.648093651936976d0, &
  -0.545421471388839d0, &
  -0.433793507626045d0, &
  -0.315042679696163d0, &
  -0.191118867473617d0, &
  -0.064056892862606d0, &
   0.064056892862606d0, &
   0.191118867473617d0, &
   0.315042679696163d0, &
   0.433793507626046d0, &
   0.545421471388840d0, &
   0.648093651936976d0, &
   0.740124191578554d0, &
   0.820001985973903d0, &
   0.886415527004401d0, &
   0.938274552002732d0, &
   0.974728555971309d0, &
   0.995187219997021d0 /)
wg(1:24) = (/ &
   0.012341229799987d0, &
   0.028531388628933d0, &
   0.044277438817421d0, &
   0.059298584915436d0, &
   0.073346481411081d0, &
   0.086190161531954d0, &
   0.097618652104114d0, &
   0.107444270115965d0, &
   0.115505668053726d0, &
   0.121670472927803d0, &
   0.125837456346828d0, &
   0.127938195346752d0, &
   0.127938195346753d0, &
   0.125837456346829d0, &
   0.121670472927803d0, &
   0.115505668053726d0, &
   0.107444270115965d0, &
   0.097618652104114d0, &
   0.086190161531953d0, &
   0.073346481411082d0, &
   0.059298584915436d0, &
   0.044277438817420d0, &
   0.028531388628933d0, &
   0.012341229799987d0 /)
  xg(1:24)=(xg(1:24)+1)/2*(b-a)+a; wg(1:24) = wg(1:24)/2*(b-a)
else if (Ng.eq.32) then
  xg(1:32) = (/ -0.997263861849481d0,  &
          -0.985611511545268d0,  &
          -0.964762255587506d0,  &
          -0.934906075937740d0,  &
          -0.896321155766052d0,  &
          -0.849367613732570d0,  &
          -0.794483795967942d0,  &
          -0.732182118740290d0,  &
          -0.663044266930215d0,  &
          -0.587715757240762d0,  &
          -0.506899908932229d0,  &
          -0.421351276130635d0,  &
          -0.331868602282128d0,  &
          -0.239287362252137d0,  &
          -0.144471961582797d0,  &
          -0.048307665687738d0,  &
           0.048307665687738d0,  &
           0.144471961582796d0,  &
           0.239287362252137d0,  &
           0.331868602282128d0,  &
           0.421351276130635d0,  &
           0.506899908932229d0,  &
           0.587715757240762d0,  &
           0.663044266930215d0,  &
           0.732182118740289d0,  &
           0.794483795967942d0,  &
           0.849367613732570d0,  &
           0.896321155766052d0,  &
           0.934906075937739d0,  &
           0.964762255587506d0,  &
           0.985611511545268d0,  &
           0.997263861849481d0 /)
  wg(1:32) = (/ 0.007018610009470d0,  &
   0.016274394730905d0,  &
   0.025392065309263d0,  &
   0.034273862913021d0,  &
   0.042835898022227d0,  &
   0.050998059262377d0,  &
   0.058684093478536d0,  &
   0.065822222776362d0,  &
   0.072345794108848d0,  &
   0.078193895787070d0,  &
   0.083311924226947d0,  &
   0.087652093004404d0,  &
   0.091173878695764d0,  &
   0.093844399080804d0,  &
   0.095638720079275d0,  &
   0.096540088514728d0,  &
   0.096540088514728d0,  &
   0.095638720079275d0,  &
   0.093844399080804d0,  &
   0.091173878695764d0,  &
   0.087652093004404d0,  &
   0.083311924226946d0,  &
   0.078193895787070d0,  &
   0.072345794108849d0,  &
   0.065822222776362d0,  &
   0.058684093478535d0,  &
   0.050998059262377d0,  &
   0.042835898022227d0,  &
   0.034273862913021d0,  &
   0.025392065309262d0,  &
   0.016274394730905d0,  &
   0.007018610009470d0 /)
xg(1:32)=(xg(1:32)+1)/2*(b-a)+a; wg(1:32) = wg(1:32)/2*(b-a)
else if (Ng.eq.40) then
xg(1:40) =(/ &
  -0.998237709710559d0, &
  -0.990726238699457d0, &
  -0.977259949983774d0, &
  -0.957916819213792d0, &
  -0.932812808278676d0, &
  -0.902098806968874d0, &
  -0.865959503212259d0, &
  -0.824612230833312d0, &
  -0.778305651426519d0, &
  -0.727318255189927d0, &
  -0.671956684614179d0, &
  -0.612553889667980d0, &
  -0.549467125095128d0, &
  -0.483075801686178d0, &
  -0.413779204371605d0, &
  -0.341994090825758d0, &
  -0.268152185007254d0, &
  -0.192697580701371d0, &
  -0.116084070675255d0, &
  -0.038772417506051d0, &
   0.038772417506051d0, &
   0.116084070675255d0, &
   0.192697580701371d0, &
   0.268152185007254d0, &
   0.341994090825759d0, &
   0.413779204371605d0, &
   0.483075801686178d0, &
   0.549467125095128d0, &
   0.612553889667980d0, &
   0.671956684614180d0, &
   0.727318255189927d0, &
   0.778305651426519d0, &
   0.824612230833311d0, &
   0.865959503212260d0, &
   0.902098806968874d0, &
   0.932812808278676d0, &
   0.957916819213791d0, &
   0.977259949983774d0, &
   0.990726238699457d0, &
   0.998237709710559d0/)
wg(1:40) = (/ &
   0.004521277098533d0, &
   0.010498284531153d0, &
   0.016421058381908d0, &
   0.022245849194167d0, &
   0.027937006980023d0, &
   0.033460195282548d0, &
   0.038782167974472d0, &
   0.043870908185673d0, &
   0.048695807635072d0, &
   0.053227846983937d0, &
   0.057439769099392d0, &
   0.061306242492929d0, &
   0.064804013456601d0, &
   0.067912045815234d0, &
   0.070611647391287d0, &
   0.072886582395804d0, &
   0.074723169057968d0, &
   0.076110361900626d0, &
   0.077039818164249d0, &
   0.077505947978425d0, &
   0.077505947978425d0, &
   0.077039818164248d0, &
   0.076110361900626d0, &
   0.074723169057969d0, &
   0.072886582395804d0, &
   0.070611647391286d0, &
   0.067912045815234d0, &
   0.064804013456601d0, &
   0.061306242492929d0, &
   0.057439769099392d0, &
   0.053227846983936d0, &
   0.048695807635072d0, &
   0.043870908185673d0, &
   0.038782167974471d0, &
   0.033460195282548d0, &
   0.027937006980023d0, &
   0.022245849194167d0, &
   0.016421058381908d0, &
   0.010498284531153d0, &
   0.004521277098533d0 /)
xg(1:40)=(xg(1:40)+1)/2*(b-a)+a; wg(1:40) = wg(1:40)/2*(b-a)
else if (Ng.eq.48) then
xg(1:48) =(/ &
  -0.998771007252426d0, &
  -0.993530172266351d0, &
  -0.984124583722827d0, &
  -0.970591592546247d0, &
  -0.952987703160431d0, &
  -0.931386690706554d0, &
  -0.905879136715570d0, &
  -0.876572020274248d0, &
  -0.843588261624393d0, &
  -0.807066204029443d0, &
  -0.767159032515740d0, &
  -0.724034130923815d0, &
  -0.677872379632664d0, &
  -0.628867396776513d0, &
  -0.577224726083972d0, &
  -0.523160974722233d0, &
  -0.466902904750958d0, &
  -0.408686481990717d0, &
  -0.348755886292161d0, &
  -0.287362487355456d0, &
  -0.224763790394689d0, &
  -0.161222356068892d0, &
  -0.097004699209463d0, &
  -0.032380170962869d0, &
   0.032380170962869d0, &
   0.097004699209463d0, &
   0.161222356068892d0, &
   0.224763790394689d0, &
   0.287362487355456d0, &
   0.348755886292161d0, &
   0.408686481990717d0, &
   0.466902904750959d0, &
   0.523160974722233d0, &
   0.577224726083972d0, &
   0.628867396776513d0, &
   0.677872379632664d0, &
   0.724034130923815d0, &
   0.767159032515740d0, &
   0.807066204029443d0, &
   0.843588261624393d0, &
   0.876572020274248d0, &
   0.905879136715570d0, &
   0.931386690706554d0, &
   0.952987703160431d0, &
   0.970591592546247d0, &
   0.984124583722827d0, &
   0.993530172266351d0, &
   0.998771007252426d0 /)
wg(1:48) =(/ &
   0.003153346052306d0, &
   0.007327553901276d0, &
   0.011477234579235d0, &
   0.015579315722944d0, &
   0.019616160457355d0, &
   0.023570760839325d0, &
   0.027426509708357d0, &
   0.031167227832798d0, &
   0.034777222564770d0, &
   0.038241351065830d0, &
   0.041545082943465d0, &
   0.044674560856694d0, &
   0.047616658492491d0, &
   0.050359035553855d0, &
   0.052890189485194d0, &
   0.055199503699984d0, &
   0.057277292100403d0, &
   0.059114839698395d0, &
   0.060704439165894d0, &
   0.062039423159893d0, &
   0.063114192286254d0, &
   0.063924238584648d0, &
   0.064466164435950d0, &
   0.064737696812684d0, &
   0.064737696812684d0, &
   0.064466164435950d0, &
   0.063924238584648d0, &
   0.063114192286254d0, &
   0.062039423159893d0, &
   0.060704439165893d0, &
   0.059114839698396d0, &
   0.057277292100403d0, &
   0.055199503699984d0, &
   0.052890189485194d0, &
   0.050359035553855d0, &
   0.047616658492490d0, &
   0.044674560856694d0, &
   0.041545082943465d0, &
   0.038241351065830d0, &
   0.034777222564770d0, &
   0.031167227832798d0, &
   0.027426509708357d0, &
   0.023570760839324d0, &
   0.019616160457356d0, &
   0.015579315722944d0, &
   0.011477234579234d0, &
   0.007327553901277d0, &
   0.003153346052306d0 /)
   xg(1:48)=(xg(1:48)+1)/2*(b-a)+a; wg(1:48) = wg(1:48)/2*(b-a);
else if (Ng.eq.64) then
xg(1:64) =(/ -0.999305041735772d0, &
        -0.996340116771955d0, &
        -0.991013371476744d0, &
        -0.983336253884626d0, &
        -0.973326827789911d0, &
        -0.961008799652054d0, &
        -0.946411374858403d0, &
         -0.92956917213194d0, &
        -0.910522137078502d0, &
        -0.889315445995114d0, &
        -0.865999398154093d0, &
         -0.84062929625258d0, &
        -0.813265315122798d0, &
        -0.783972358943341d0, &
        -0.752819907260532d0, &
        -0.719881850171611d0, &
        -0.685236313054233d0, &
        -0.648965471254657d0, &
        -0.611155355172393d0, &
        -0.571895646202634d0, &
        -0.531279464019894d0, &
        -0.489403145707053d0, &
        -0.446366017253464d0, &
        -0.402270157963992d0, &
        -0.357220158337668d0, &
        -0.311322871990211d0, &
        -0.264687162208768d0, &
        -0.217423643740007d0, &
        -0.169644420423993d0, &
        -0.121462819296121d0, &
       -0.0729931217877992d0, &
       -0.0243502926634247d0, &
        0.0243502926634245d0, &
        0.0729931217877988d0, &
          0.12146281929612d0, &
         0.169644420423993d0, &
         0.217423643740007d0, &
         0.264687162208767d0, &
         0.311322871990211d0, &
         0.357220158337668d0, &
         0.402270157963992d0, &
         0.446366017253464d0, &
         0.489403145707053d0, &
         0.531279464019894d0, &
         0.571895646202634d0, &
         0.611155355172393d0, &
         0.648965471254657d0, &
         0.685236313054233d0, &
         0.719881850171611d0, &
         0.752819907260532d0, &
         0.783972358943341d0, &
         0.813265315122798d0, &
          0.84062929625258d0, &
         0.865999398154093d0, &
         0.889315445995114d0, &
         0.910522137078503d0, &
          0.92956917213194d0, &
         0.946411374858403d0, &
         0.961008799652054d0, &
         0.973326827789911d0, &
         0.983336253884626d0, &
         0.991013371476744d0, &
         0.996340116771955d0, &
         0.999305041735772d0/)
wg(1:64) = (/ 0.00178328072169649d0, &
       0.00414703326056235d0, &
       0.00650445796897812d0, &
       0.00884675982636353d0, &
        0.0111681394601316d0, &
         0.013463047896719d0, &
         0.015726030476025d0, &
         0.017951715775697d0, &
          0.02013482315353d0, &
        0.0222701738083834d0, &
         0.024352702568711d0, &
        0.0263774697150552d0, &
        0.0283396726142598d0, &
        0.0302346570724022d0, &
        0.0320579283548518d0, &
        0.0338051618371415d0, &
        0.0354722132568822d0, &
        0.0370551285402401d0, &
        0.0385501531786157d0, &
        0.0399537411327198d0, &
        0.0412625632426236d0, &
        0.0424735151236539d0, &
        0.0435837245293231d0, &
        0.0445905581637565d0, &
        0.0454916279274186d0, &
        0.0462847965813137d0, &
          0.04696818281621d0, &
        0.0475401657148305d0, &
        0.0479993885964581d0, &
        0.0483447622348029d0, &
        0.0485754674415039d0, &
        0.0486909570091394d0, &
        0.0486909570091398d0, &
        0.0485754674415038d0, &
        0.0483447622348026d0, &
        0.0479993885964582d0, &
        0.0475401657148303d0, &
        0.0469681828162095d0, &
        0.0462847965813141d0, &
        0.0454916279274184d0, &
        0.0445905581637563d0, &
        0.0435837245293234d0, &
        0.0424735151236538d0, &
        0.0412625632426232d0, &
        0.0399537411327201d0, &
        0.0385501531786157d0, &
        0.0370551285402399d0, &
        0.0354722132568831d0, &
         0.033805161837142d0, &
        0.0320579283548514d0, &
        0.0302346570724021d0, &
        0.0283396726142594d0, &
        0.0263774697150548d0, &
         0.024352702568711d0, &
        0.0222701738083837d0, &
        0.0201348231535299d0, &
        0.0179517157756972d0, &
        0.0157260304760251d0, &
        0.0134630478967187d0, &
        0.0111681394601311d0, &
       0.00884675982636387d0, &
       0.00650445796897822d0, &
       0.00414703326056219d0, &
        0.0017832807216965d0/)
  xg(1:64)=(xg(1:64)+1)/2*(b-a)+a; wg(1:64) = wg(1:64)/2*(b-a)
else if (Ng.eq.128) then
xg(1:128) =(/-0.999824887947132d0, &
        -0.999077459977376d0, &
        -0.997733248625514d0, &
        -0.995792758534981d0, &
        -0.993257112900213d0, &
        -0.990127818491734d0, &
        -0.986406742724586d0, &
        -0.982096108435718d0, &
        -0.977198491463907d0, &
        -0.971716818747136d0, &
        -0.965654366431965d0, &
          -0.9590147578537d0, &
        -0.951801961341264d0, &
         -0.94402028783022d0, &
        -0.935674388277916d0, &
        -0.926769250878948d0, &
         -0.91731019808096d0, &
        -0.907302883401757d0, &
        -0.896753288049158d0, &
        -0.885667717345397d0, &
        -0.874052796958032d0, &
        -0.861915468939548d0, &
        -0.849262987577969d0, &
        -0.836102915060907d0, &
        -0.822443116955644d0, &
        -0.808291757507914d0, &
        -0.793657294762193d0, &
        -0.778548475506412d0, &
        -0.762974330044094d0, &
        -0.746944166797062d0, &
        -0.730467566741909d0, &
        -0.713554377683587d0, &
        -0.696214708369514d0, &
        -0.678458922447719d0, &
        -0.660297632272646d0, &
        -0.641741692562307d0, &
        -0.622802193910585d0, &
        -0.603490456158549d0, &
        -0.583818021628763d0, &
        -0.563796648226618d0, &
         -0.54343830241281d0, &
        -0.522755152051175d0, &
        -0.501759559136144d0, &
        -0.480464072404172d0, &
        -0.458881419833552d0, &
        -0.437024501037104d0, &
        -0.414906379552275d0, &
        -0.392540275033267d0, &
        -0.369939555349859d0, &
        -0.347117728597636d0, &
        -0.324088435024414d0, &
        -0.300865438877677d0, &
        -0.277462620177904d0, &
        -0.253893966422694d0, &
         -0.23017356422666d0, &
        -0.206315590902079d0, &
        -0.182334305985337d0, &
        -0.158244042714225d0, &
        -0.134059199461188d0, &
        -0.109794231127644d0, &
       -0.0854636405045156d0, &
       -0.0610819696041397d0, &
       -0.0366637909687337d0, &
       -0.0122236989606159d0, &
        0.0122236989606157d0, &
        0.0366637909687333d0, &
        0.0610819696041394d0, &
        0.0854636405045153d0, &
         0.109794231127644d0, &
         0.134059199461188d0, &
         0.158244042714225d0, &
         0.182334305985337d0, &
         0.206315590902079d0, &
          0.23017356422666d0, &
         0.253893966422694d0, &
         0.277462620177904d0, &
         0.300865438877677d0, &
         0.324088435024413d0, &
         0.347117728597635d0, &
         0.369939555349859d0, &
         0.392540275033267d0, &
         0.414906379552275d0, &
         0.437024501037104d0, &
         0.458881419833552d0, &
         0.480464072404172d0, &
         0.501759559136144d0, &
         0.522755152051175d0, &
          0.54343830241281d0, &
         0.563796648226618d0, &
         0.583818021628763d0, &
         0.603490456158549d0, &
         0.622802193910585d0, &
         0.641741692562308d0, &
         0.660297632272646d0, &
          0.67845892244772d0, &
         0.696214708369514d0, &
         0.713554377683587d0, &
         0.730467566741908d0, &
         0.746944166797062d0, &
         0.762974330044095d0, &
         0.778548475506412d0, &
         0.793657294762193d0, &
         0.808291757507914d0, &
         0.822443116955644d0, &
         0.836102915060907d0, &
         0.849262987577969d0, &
         0.861915468939548d0, &
         0.874052796958032d0, &
         0.885667717345397d0, &
         0.896753288049158d0, &
         0.907302883401757d0, &
         0.917310198080961d0, &
         0.926769250878948d0, &
         0.935674388277916d0, &
          0.94402028783022d0, &
         0.951801961341264d0, &
           0.9590147578537d0, &
         0.965654366431965d0, &
         0.971716818747137d0, &
         0.977198491463907d0, &
         0.982096108435718d0, &
         0.986406742724586d0, &
         0.990127818491734d0, &
         0.993257112900213d0, &
         0.995792758534981d0, &
         0.997733248625514d0, &
         0.999077459977376d0, &
         0.999824887947132d0 /)
wg(1:128) =(/ 0.000449380960292082d0, &
       0.00104581267934044d0, &
       0.00164250301866894d0, &
       0.00223828843096265d0, &
       0.00283275147145784d0, &
       0.00342552604091012d0, &
       0.00401625498373828d0, &
       0.00460458425670286d0, &
       0.00519016183267656d0, &
       0.00577263754286589d0, &
       0.00635166316170743d0, &
       0.00692689256689895d0, &
       0.00749798192563484d0, &
       0.00806458989048619d0, &
        0.0086263777986166d0, &
       0.00918300987166053d0, &
       0.00973415341500667d0, &
         0.010279479015832d0, &
        0.0108186607395032d0, &
        0.0113513763240804d0, &
        0.0118773073727404d0, &
         0.012396139543951d0, &
        0.0129075627392675d0, &
        0.0134112712886167d0, &
        0.0139069641329522d0, &
         0.014394345004167d0, &
        0.0148731226021472d0, &
        0.0153430107688649d0, &
        0.0158037286593993d0, &
        0.0162550009097853d0, &
        0.0166965578015891d0, &
        0.0171281354231113d0, &
        0.0175494758271176d0, &
        0.0179603271850086d0, &
        0.0183604439373313d0, &
        0.0187495869405448d0, &
        0.0191275236099511d0, &
        0.0194940280587068d0, &
        0.0198488812328304d0, &
        0.0201918710421298d0, &
        0.0205227924869599d0, &
        0.0208414477807513d0, &
        0.0211476464682217d0, &
        0.0214412055392084d0, &
        0.0217219495380521d0, &
        0.0219897106684601d0, &
        0.0222443288937999d0, &
        0.0224856520327449d0, &
        0.0227135358502367d0, &
        0.0229278441436872d0, &
        0.0231284488243865d0, &
        0.0233152299940625d0, &
        0.0234880760165359d0, &
        0.0236468835844475d0, &
        0.0237915577810033d0, &
        0.0239220121367038d0, &
        0.0240381686810241d0, &
        0.0241399579890189d0, &
         0.024227319222815d0, &
        0.0243002001679718d0, &
        0.0243585572646908d0, &
          0.02440235563385d0, &
        0.0244315690978502d0, &
        0.0244461801962621d0, &
        0.0244461801962623d0, &
        0.0244315690978504d0, &
        0.0244023556338498d0, &
        0.0243585572646909d0, &
        0.0243002001679716d0, &
        0.0242273192228151d0, &
         0.024139957989019d0, &
        0.0240381686810241d0, &
        0.0239220121367037d0, &
        0.0237915577810031d0, &
        0.0236468835844477d0, &
        0.0234880760165355d0, &
        0.0233152299940628d0, &
        0.0231284488243867d0, &
        0.0229278441436869d0, &
        0.0227135358502366d0, &
        0.0224856520327449d0, &
        0.0222443288937995d0, &
        0.0219897106684606d0, &
         0.021721949538052d0, &
        0.0214412055392084d0, &
        0.0211476464682218d0, &
         0.020841447780751d0, &
        0.0205227924869596d0, &
        0.0201918710421297d0, &
        0.0198488812328309d0, &
        0.0194940280587067d0, &
        0.0191275236099511d0, &
        0.0187495869405448d0, &
         0.018360443937331d0, &
         0.017960327185009d0, &
         0.017549475827118d0, &
        0.0171281354231116d0, &
        0.0166965578015893d0, &
        0.0162550009097851d0, &
        0.0158037286593994d0, &
         0.015343010768865d0, &
        0.0148731226021471d0, &
         0.014394345004167d0, &
        0.0139069641329518d0, &
        0.0134112712886162d0, &
        0.0129075627392677d0, &
        0.0123961395439509d0, &
        0.0118773073727402d0, &
        0.0113513763240807d0, &
        0.0108186607395032d0, &
        0.0102794790158322d0, &
       0.00973415341500671d0, &
       0.00918300987166061d0, &
       0.00862637779861692d0, &
       0.00806458989048616d0, &
       0.00749798192563537d0, &
       0.00692689256689874d0, &
       0.00635166316170716d0, &
       0.00577263754286538d0, &
       0.00519016183267642d0, &
       0.00460458425670292d0, &
       0.00401625498373874d0, &
       0.00342552604091024d0, &
       0.00283275147145767d0, &
       0.00223828843096271d0, &
       0.00164250301866869d0, &
       0.00104581267934019d0, &
      0.000449380960292315d0 /)
  xg(1:128)=(xg(1:128)+1)/2*(b-a)+a; wg(1:128) = wg(1:128)/2*(b-a)
else
stop 'STOP Ng value not allowed'
endif
END SUBROUTINE setgauss

END MODULE mod_EHpqr
