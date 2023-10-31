#include "flags.h"
MODULE mod_stokes
! contains all (13) routines to compute and correct SLP and DLP 
! for Stokes equation at a target near the surface, 
! for a lat-long grid of that surface
!
! The routines in this module work for any surface described by 
! lat-long grid, and target point x0, 
! and are independent of the specific geometry (surface) and 
! specific problem (density + more)
!
! A. Surface specific routines  in mod_geom
!      setgrid,  
!      velo,  
!      setbase
!      compcoeffx
!    are in mod_geom
!
! B. Problem specific routines for SLP application
!      setdens - set density sgle
!      compcoeffdensex - set coeff of density 
!      exactvelo45 - computes velo past sphere at 45^o
!    are in mod_init
!
!XC. Problem specific routines for DLP application
!X      setdensdb - set density sgle
!X      compcoeffdensdb - set coeff of density 
!X  are in mod_sphere_tstdble
!
!
! The only routine that needs to be called from outside: 
! compslp(x0,t,slp)  t:basepoint
! compdlp(x0,t,dlp) 
!
! contains
!    compslp
!    compdlp
!    correctg      (computes Sum(cjk*E[Hjk])
!    correctt      (computes Sum(cjk*E[Hjk] for dble layer
!    comprho (function, corrects 1/rho terms, given numerator)
!    comprho3(function, corrects 1/rho^3 terms, given numerator)
!    comprho5(function, corrects 1/rho^5 terms, given numerator)


!PRIVATE    !set the default for module
PUBLIC

CONTAINS
   SUBROUTINE compslp(x0,t,slp) 
   !compute and correct SLP
   use params,     ONLY: rmax,pmax,nbmax
   use types,      ONLY : basept
   use globalvars, ONLY: g,calf,cbet,calfbet
   use mod_SDLP         ! for SLPatP
   use mod_compcoeffdens
   use mod_EHpqr
   implicit none
   real*8, dimension(3), INTENT(IN) :: x0
   real*8, dimension(3), INTENT(OUT) :: slp
   TYPE (basept), INTENT(IN) :: t  
!LOCAL
   integer :: l
   real*8, dimension(3) :: corrslp
   integer, dimension(rmax,0:pmax,0:pmax) :: jind
   real*8, dimension(nbmax) :: corrbasis

   l=t%igrid
   call SLPatP(t,g(l),slp)
   if (t%icorr.eq.1) then
      call compcoeffdens(g(l),t)  
!      call compcoeffx(t) now done in setbase(mod_geom), if corr 
      call compallcorrtrap(t,g(l),calf,cbet,calfbet,corrbasis,jind)
      call correctg(t%alfb,t%betb,corrbasis,jind,corrslp)
      slp = slp+corrslp
   endif

   return
   END SUBROUTINE compslp


   SUBROUTINE compdlp(x0,t,dlp) 
   !compute and correct SLP
   use params, ONLY: rmax,pmax,nbmax
   use types, ONLY : basept
   use globalvars, ONLY: g,calf,cbet,calfbet
   use mod_SDLP           ! for DLPatP
   use mod_compcoeffdens
   use mod_EHpqr
   implicit none
   real*8, dimension(3), INTENT(IN) :: x0
   real*8, dimension(3), INTENT(OUT) :: dlp
   TYPE (basept), INTENT(IN) :: t  
!LOCAL
   integer :: l
   real*8, dimension(3) :: corrdlp
   integer, dimension(rmax,0:pmax,0:pmax) :: jind
   real*8, dimension(nbmax) :: corrbasis

   l=t%igrid
   call DLPatP(t,g(l),dlp)
   if (t%icorr.eq.1) then
      call compcoeffdens(g(l),t)  
!      call compcoeffx(t) now done in setbase(mod_geom), if corr 
      call compallcorrtrap(t,g(l),calf,cbet,calfbet,corrbasis,jind)
      call correctt(t%alfb,t%betb,corrbasis,jind,corrdlp)
      dlp = dlp+corrdlp
   endif

   return
   END SUBROUTINE compdlp


   SUBROUTINE compsdlp(x0,t,slp,dlp) 
! compute and correct both slp and dlp 
   use params,     ONLY: rmax,pmax,nbmax,nleft
   use types,      ONLY : basept
   use globalvars, ONLY: g,calf,cbet,calfbet
   use mod_SDLP        ! for SLPatP
   use mod_compcoeffdens
   use mod_EHpqr
   implicit none
   real*8, dimension(3), INTENT(IN) :: x0
   TYPE (basept), INTENT(IN) :: t  
   real*8, dimension(3), INTENT(OUT) :: slp,dlp
!LOCAL
   integer :: l
   real*8, dimension(3) :: corrslp,corrdlp
   integer, dimension(rmax,0:pmax,0:pmax) :: jind
   real*8, dimension(nbmax) :: corrbasis

   l=t%igrid
!print*,size(g)
!print*,size(g(l)%f1,1)
!print*,size(g(l)%f1,2)
!print*,g(l)%f1(320, 100)
!print*,g(l)%f1(nleft, 0)
   call SLPatP(t,g(l),slp)
   call DLPatP(t,g(l),dlp)
!print*,slp
   corrslp=0
   corrdlp=0
!   print*,'here'
!   print*,dlp
   if (t%icorr.eq.1) then
      call compcoeffdens(g(l),t)  
!      call compcoeffx(t) now done in setbase(mod_geom), if corr 
      call compallcorrtrap(t,g(l),calf,cbet,calfbet,corrbasis,jind)
      call correctg(t%alfb,t%betb,corrbasis,jind,corrslp)
      call correctt(t%alfb,t%betb,corrbasis,jind,corrdlp)
!print*,t%alfb,t%betb,slp,slp+corrslp
!print*,dlp,dlp+corrdlp
      slp = slp+corrslp
      dlp = dlp+corrdlp
   endif
!      print*,slp,dlp
!print*,dlp
!print*,corrdlp
!stop
   END SUBROUTINE compsdlp


SUBROUTINE correctg(alfb,betb,corrbasis,jind,corrg)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! A testcase where f=sigma n=n=x (on unit sphere) 
!          with coeff cf1=cx, cf2=cy, cf3=cz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : rmax,pmax,cmax,pi8
!use coeff
use globalvars, ONLY : cx,cy,cz,cj,e,cf1,cf2,cf3
use products
implicit none
integer, dimension(rmax,0:pmax,0:pmax), INTENT(IN) :: jind
real*8, INTENT(IN) :: corrbasis(*),alfb,betb
real*8, INTENT(OUT) :: corrg(*)
! LOCAL
integer :: j,k
real*8, dimension(0:cmax,0:cmax) :: fjac1,fjac2,fjac3,fxjac1,fxjac2,fxjac3,dotfx,numer
real*8 :: corr1,corr2

call prod222(cf1,cj,fjac1) 
call prod222(cf2,cj,fjac2)
call prod222(cf3,cj,fjac3)

call prod233(fjac1,cx,fxjac1) 
call prod233(fjac2,cy,fxjac2)
call prod233(fjac3,cz,fxjac3)

do j=0,4
do k=0,4
  dotfx(j,k) = fxjac1(j,k) + fxjac2(j,k) + fxjac3(j,k)
enddo
enddo

! U-Component --------------------------------------------------
call prod334(dotfx,cx,numer)            
corr1=comprho(fjac1,e,corrbasis,jind)   ! 1/RHO 
corr2=comprho3(numer,e,corrbasis,jind)  ! 1/RHO3 
corrg(1)=corr1+corr2
      
! V-Component --------------------------------------------------
call prod334(dotfx,cy,numer) 
corr1=comprho(fjac2,e,corrbasis,jind)  ! 1/RHO 
corr2=comprho3(numer,e,corrbasis,jind) ! 1/RHO3 
corrg(2)=corr1+corr2
      
! W-Component --------------------------------------------------
call prod334(dotfx,cz,numer) 
corr1=comprho(fjac3,e,corrbasis,jind)  ! 1/RHO 
corr2=comprho3(numer,e,corrbasis,jind) ! 1/RHO3 
corrg(3)=corr1+corr2

corrg(1:3)=corrg(1:3)/pi8     !Divide here and in SLPatP
return
END SUBROUTINE correctg


SUBROUTINE correctt(alfb,betb,corrbasis,jind,corr)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! to start: 
! testcase where f=sigma n=n=x (on unit sphere) 
!          with coeff cf1=cx, cf2=cy, cf3=cz
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : rmax,pmax,cmax,pi4
!use coeff
use globalvars, ONLY : cx,cy,cz,cnjac1,cnjac2,cnjac3,e,cf1,cf2,cf3
use products
implicit none
integer, dimension(rmax,0:pmax,0:pmax), INTENT(IN) :: jind
real*8, INTENT(IN) :: corrbasis(*), alfb,betb
real*8, INTENT(OUT) :: corr(*)
! LOCAL
integer :: j,k
real*8, dimension(0:cmax,0:cmax) :: cnxjac1,cnxjac2,cnxjac3, &
                  cfx1,cfx2,cfx3,dotnxjac,dotfx,g,numer
      
!call prod333(cnx,cj,cnjac1)
!call prod333(cny,cj,cnjac2)
!call prod333(cnz,cj,cnjac3)

call prod344(cnjac1,cx,cnxjac1) 
call prod344(cnjac2,cy,cnxjac2) 
call prod344(cnjac3,cz,cnxjac3) 

call prod344(cf1,cx,cfx1) 
call prod344(cf2,cy,cfx2) 
call prod344(cf3,cz,cfx3) 

do j=0,4
do k=0,4
  dotnxjac(j,k)= cnxjac1(j,k)+cnxjac2(j,k)+cnxjac3(j,k)
  dotfx(j,k)  = cfx1(j,k)+cfx2(j,k)+cfx3(j,k)
enddo
enddo
call prod445(dotnxjac,dotfx,g) 

! U-Component --------------------------------------------------
call prod456(cx,g,numer) 
corr(1)=comprho5(numer,e,corrbasis,jind)     !RHO5

! V-Component --------------------------------------------------
call prod456(cy,g,numer) 
corr(2)=comprho5(numer,e,corrbasis,jind)     !RHO5

! W-Component --------------------------------------------------
call prod456(cz,g,numer) 
corr(3)=comprho5(numer,e,corrbasis,jind)     !RHO5

corr(1:3) = -(3/pi4)*corr(1:3)
return
END SUBROUTINE correctt


FUNCTION comprho(f,e,corr,jind) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     see 3dexpinvr for all cs
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : rmax,pmax,cmax,nbmax
implicit none
! INPUT      
integer, INTENT(IN) :: jind(rmax,0:pmax,0:pmax)
MODE, dimension(0:cmax,0:cmax),INTENT(IN) :: f,e
MODE, INTENT(IN) :: corr(*)
!OUTPUT
MODE :: comprho
! LOCAL      
integer :: j,r,p,j0
MODE :: c(nbmax),sumo

! Coefficients of = 1/r terms, need O(alf^0 ... alf^2) on top
! h^2/d terms
c(jind(1,0,0)) = f(0,0) 

! h^2 terms
c(jind(1,1,0)) = f(1,0) 
c(jind(1,0,1)) = f(0,1)    

! h^2d terms
c(jind(1,2,0)) = f(2,0)  
c(jind(1,1,1)) = f(1,1)    
c(jind(1,0,2)) = f(0,2)    
      

! Coefficients of 1/r3 terms, need O(alf^3... alf^4) on top
! h^2 terms
c(jind(3,3,0)) = -e(3,0)*f(0,0)/2 
c(jind(3,2,1)) = -e(2,1)*f(0,0)/2 
c(jind(3,1,2)) = -e(1,2)*f(0,0)/2 
c(jind(3,0,3)) = -e(0,3)*f(0,0)/2 

#ifdef FOURTH
! h^2d terms
c(jind(3,4,0)) = (-e(4,0)*f(0,0) -e(3,0)*f(1,0))/2   
c(jind(3,3,1)) = (-e(3,1)*f(0,0) -e(3,0)*f(0,1) -e(2,1)*f(1,0))/2
c(jind(3,2,2)) = (-e(2,2)*f(0,0) -e(2,1)*f(0,1) -e(1,2)*f(1,0))/2
c(jind(3,1,3)) = (-e(1,3)*f(0,0) -e(1,2)*f(0,1) -e(0,3)*f(1,0))/2
c(jind(3,0,4)) = (-e(0,4)*f(0,0) -e(0,3)*f(0,1)               )/2

! Coefficients of 1/r5 terms, need O(alf^6) on top
! h^2d terms
c(jind(5,6,0)) = 3*   e(3,0)**2                     *f(0,0)/8
c(jind(5,5,1)) = 3*   e(2,1)*e(3,0)                 *f(0,0)/4
c(jind(5,4,2)) = 3*(2*e(1,2)*e(3,0) +e(2,1)**2     )*f(0,0)/8
c(jind(5,3,3)) = 3*(  e(0,3)*e(3,0) +e(1,2)*e(2,1) )*f(0,0)/4
c(jind(5,2,4)) = 3*(2*e(0,3)*e(2,1) +e(1,2)**2     )*f(0,0)/8
c(jind(5,1,5)) = 3*   e(0,3)*e(1,2)                 *f(0,0)/4
c(jind(5,0,6)) = 3*   e(0,3)**2                     *f(0,0)/8
#endif

! Compute correction for f/rho

! 2nd order method (all h^2/d^2, h^2/d terms)
! r=1, p+q=0
sumo =  c(jind(1,0,0))*corr(jind(1,0,0))  &
#ifdef THIRD
! 3rd order method (all h^2 terms)
! r=1, p+q=1
       + c(jind(1,1,0))*corr(jind(1,1,0))  &
       + c(jind(1,0,1))*corr(jind(1,0,1))  &
! r=3, p+q=3
       + c(jind(3,3,0))*corr(jind(3,3,0))  &
       + c(jind(3,2,1))*corr(jind(3,2,1))  &
       + c(jind(3,1,2))*corr(jind(3,1,2))  &
       + c(jind(3,0,3))*corr(jind(3,0,3))  &
#endif
#ifdef FOURTH
! 4th order method (all h^2 d terms)
! r=1, p+q=2
       + c(jind(1,2,0))*corr(jind(1,2,0)) & 
       + c(jind(1,1,1))*corr(jind(1,1,1))  &
       + c(jind(1,0,2))*corr(jind(1,0,2)) 
! r=3, p+q=4
r=3
j=4
do p=0,j
   j0=jind(r,p,j-p)
   sumo = sumo + c(j0)*corr(j0)
enddo
! r=5, p+q=6
r=5
j=6
do p=0,j
   j0=jind(r,p,j-p)
   sumo = sumo + c(j0)*corr(j0)
enddo
#endif
comprho=sumo
return
END FUNCTION comprho


FUNCTION comprho3(f,e,corr,jind) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     see 3dexpinvr3 for all cs
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : rmax,pmax,cmax,nbmax
implicit none
! INPUT      
integer :: jind(rmax,0:pmax,0:pmax)
MODE :: f(0:cmax,0:cmax),e(0:cmax,0:cmax),corr(*)
! OUTPUT
MODE :: comprho3
! LOCAL      
      integer :: j,r,p,j0
      MODE :: c(nbmax),sumo

! Coefficients of) = 1/r3 terms, need O(alf^0 ... alf^4) on top
! f00=d^2, f10=f01=d
! h^2/d terms
      c(jind(3,0,0)) = f(0,0)   
      c(jind(3,1,0)) = f(1,0)  
      c(jind(3,0,1)) = f(0,1) 
      c(jind(3,2,0)) = f(2,0)
      c(jind(3,1,1)) = f(1,1)  
      c(jind(3,0,2)) = f(0,2) 
! h^2 terms
      c(jind(3,3,0)) = f(3,0)
      c(jind(3,2,1)) = f(2,1)   
      c(jind(3,1,2)) = f(1,2)  
      c(jind(3,0,3)) = f(0,3) 
! h^2d terms
      c(jind(3,4,0)) = f(4,0)
      c(jind(3,3,1)) = f(3,1) 
      c(jind(3,2,2)) = f(2,2)
      c(jind(3,1,3)) = f(1,3)    
      c(jind(3,0,4)) = f(0,4)   

! Coefficients of) = 1/r5 terms, need O(alf^3 ... alf^6) on top
! h^2 terms (sometimes + h^2d)
      c(jind(5,3,0)) = -3*e(3,0)*f(0,0)/2  
      c(jind(5,2,1)) = -3*e(2,1)*f(0,0)/2  
      c(jind(5,1,2)) = -3*e(1,2)*f(0,0)/2  
      c(jind(5,0,3)) = -3*e(0,3)*f(0,0)/2  
      c(jind(5,4,0)) = -3*(e(4,0)*f(0,0) +e(3,0)*f(1,0))/2  
      c(jind(5,3,1)) = -3*(e(3,1)*f(0,0)+e(3,0)*f(0,1)+e(2,1)*f(1,0))/2
      c(jind(5,2,2)) = -3*(e(2,2)*f(0,0)+e(2,1)*f(0,1)+e(1,2)*f(1,0))/2
      c(jind(5,1,3)) = -3*(e(1,3)*f(0,0)+e(1,2)*f(0,1)+e(0,3)*f(1,0))/2 
      c(jind(5,0,4)) = -3*(e(0,4)*f(0,0)+e(0,3)*f(0,1)              )/2 

      c(jind(5,5,0)) = -3*(e(4,0)*f(1,0)+e(3,0)*f(2,0))/2
      c(jind(5,4,1)) = -3*(e(4,0)*f(0,1)+e(3,1)*f(1,0)+e(3,0)*f(1,1) &
                  +e(2,1)*f(2,0))/2
      c(jind(5,3,2)) = -3*(e(3,1)*f(0,1)+e(3,0)*f(0,2)+e(2,2)*f(1,0) &
                  +e(2,1)*f(1,1)+e(1,2)*f(2,0))/2
      c(jind(5,2,3)) = -3*(e(2,2)*f(0,1)+e(2,1)*f(0,2)+e(1,3)*f(1,0) &
                  +e(1,2)*f(1,1)+e(0,3)*f(2,0))/2
      c(jind(5,1,4)) = -3*(e(1,3)*f(0,1)+e(1,2)*f(0,2)+e(0,4)*f(1,0) &
                  +e(0,3)*f(1,1))/2
      c(jind(5,0,5)) = -3*(e(0,4)*f(0,1)+e(0,3)*f(0,2))/2

#ifdef FOURTH
! h^2d terms
      c(jind(5,6,0)) = -3*(e(4,0)*f(2,0)+e(3,0)*f(3,0))/2
      c(jind(5,5,1)) = -3*(e(4,0)*f(1,1)+e(3,1)*f(2,0)+e(3,0)*f(2,1) &
                  +e(2,1)*f(3,0))/2
      c(jind(5,4,2)) = -3*(e(4,0)*f(0,2)+e(3,1)*f(1,1)+e(3,0)*f(1,2) &
                  +e(2,2)*f(2,0)+e(2,1)*f(2,1)+e(1,2)*f(3,0))/2
      c(jind(5,3,3)) = -3*(e(3,1)*f(0,2)+e(3,0)*f(0,3)+e(2,2)*f(1,1) &
                  +e(2,1)*f(1,2)+e(1,3)*f(2,0)+e(1,2)*f(2,1) &
                  +e(0,3)*f(3,0))/2
      c(jind(5,2,4)) = -3*(e(2,2)*f(0,2)+e(2,1)*f(0,3)+e(1,3)*f(1,1) &
                  +e(1,2)*f(1,2)+e(0,4)*f(2,0)+e(0,3)*f(2,1))/2
      c(jind(5,1,5)) = -3*(e(1,3)*f(0,2)+e(1,2)*f(0,3)+e(0,4)*f(1,1) &
                  +e(0,3)*f(1,2))/2
      c(jind(5,0,6)) = -3*(e(0,4)*f(0,2)+e(0,3)*f(0,3))/2

! Coefficients of 1/r7 terms, need O(alf^6... alf^8) on top, f00=d^2, f01=f10=d

! h^2d terms
      c(jind(7,6,0)) = 15*e(3,0)**2*f(0,0)/8
      c(jind(7,5,1)) = 15*e(2,1)*e(3,0)*f(0,0)/4
      c(jind(7,4,2)) = 15*(e(2,1)**2+2*e(1,2)*e(3,0))*f(0,0)/8
      c(jind(7,3,3)) = 15*(e(1,2)*e(2,1)+e(0,3)*e(3,0))*f(0,0)/4
      c(jind(7,2,4)) = 15*(e(1,2)**2+2*e(0,3)*e(2,1))*f(0,0)/8
      c(jind(7,1,5)) = 15*e(0,3)*e(1,2)*f(0,0)/4
      c(jind(7,0,6)) = 15*e(0,3)**2*f(0,0)/8

      c(jind(7,7,0)) = 15*e(3,0)**2*f(1,0)/8
      c(jind(7,6,1)) = 15*e(3,0)*(e(3,0)*f(0,1)+2*e(2,1)*f(1,0))/8
      c(jind(7,5,2)) = 15*(2*e(2,1)*e(3,0)*f(0,1)+e(2,1)**2*f(1,0) &
                  +2*e(1,2)*e(3,0)*f(1,0))/8
      c(jind(7,4,3)) = 15*(e(2,1)**2*f(0,1)+2*e(1,2)*e(2,1)*f(1,0) &
                  +2*e(3,0)*(e(1,2)*f(0,1)+e(0,3)*f(1,0)))/8
      c(jind(7,3,4)) = 15*(2*e(1,2)*e(2,1)*f(0,1)+e(1,2)**2*f(1,0) &
                  +2*e(0,3)*(e(3,0)*f(0,1)+e(2,1)*f(1,0)))/8
      c(jind(7,2,5)) = 15*(e(1,2)**2*f(0,1)+2*e(0,3)*e(2,1)*f(0,1) &
                  +2*e(0,3)*e(1,2)*f(1,0))/8
      c(jind(7,1,6)) = 15*e(0,3)*(2*e(1,2)*f(0,1)+e(0,3)*f(1,0))/8
      c(jind(7,0,7)) = 15*e(0,3)**2*f(0,1)/8

      c(jind(7,8,0)) = 15*e(3,0)**2*f(2,0)/8
      c(jind(7,7,1)) = 15*e(3,0)*(e(3,0)*f(1,1)+2*e(2,1)*f(2,0))/8
      c(jind(7,6,2)) = 15*(e(3,0)**2*f(0,2)+e(2,1)**2*f(2,0) &
                  +2*e(3,0)*(e(2,1)*f(1,1)+e(1,2)*f(2,0)))/8
      c(jind(7,5,3)) = 15*(e(2,1)**2*f(1,1)+2*e(3,0)*(e(1,2)*f(1,1) &
                  +e(0,3)*f(2,0))+2*e(2,1)*(e(3,0)*f(0,2) &
                  +e(1,2)*f(2,0)))/8
      c(jind(7,4,4)) = 15*(e(2,1)**2*f(0,2)+2*e(1,2)*e(3,0)*f(0,2) &
                  +2*e(0,3)*e(3,0)*f(1,1)+e(1,2)**2*f(2,0) &
                  +2*e(2,1)*(e(1,2)*f(1,1)+e(0,3)*f(2,0)))/8
      c(jind(7,3,5)) = 15*(e(1,2)**2*f(1,1)+2*e(0,3)*(e(3,0)*f(0,2)  &
                  +e(2,1)*f(1,1))+2*e(1,2)*(e(2,1)*f(0,2)  &
                  +e(0,3)*f(2,0)))/8
      c(jind(7,2,6)) = 15*(e(1,2)**2*f(0,2)+2*e(0,3)*e(1,2)*f(1,1)  &
                  +e(0,3)*(2*e(2,1)*f(0,2)+e(0,3)*f(2,0)))/8
      c(jind(7,1,7)) = 15*e(0,3)*(2*e(1,2)*f(0,2)+e(0,3)*f(1,1))/8
      c(jind(7,0,8)) = 15*e(0,3)**2*f(0,2)/8
#endif


! 2nd order method (all h^2/d^2, h^2/d terms)
! r=3, p+q=0,2
      sumo=0
      r=3
      do j=0,2
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
!          print*,r,p,j-p,c(j0),corr(j0),c(j0)*corr(j0)
        enddo
      enddo
#ifdef THIRD
! 3rd order method (all h^2 terms)
! r=3, p+q=3
! r=5, p+q=3,5
      r=3
      j=3
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      r=5
      do j=3,5
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      enddo
#endif
#ifdef FOURTH
! 4th order method (all h^2 d terms)
! r=3, p+q=4
! r=5, p+q=6
! r=7, p+q=6,8
      r=3
      j=4
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      r=5
      j=6
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      r=7
      do j=6,8
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      enddo
#endif
      comprho3=sumo
! comparison with old (in simpletestFeb22)
!      j0=jind(3,0,0) ! agrees all digits!
!      j0=jind(3,0,1) ! agrees 14 digits
!      j0=jind(3,1,0) ! agrees 14 digits
!      j0=jind(3,2,0) ! c(j0) agrees ONLY 4 digits
!      write(*,'(3e24.14)')c(j0)*corr(j0),c(j0),corr(j0) 
      return
END FUNCTION comprho3


FUNCTION comprho5(f,e,corr,jind) 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     see 3dGExph4 for all expansions
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use params, ONLY : rmax,pmax,cmax,nbmax
implicit none
! INPUT      
integer :: jind(rmax,0:pmax,0:pmax)
MODE :: f(0:cmax,0:cmax),e(0:cmax,0:cmax),corr(*)
! OUTPUT
MODE :: comprho5
! LOCAL      
integer :: j,r,p,j0
MODE :: c(nbmax),sumo

! Coefficients of) = 1/r5 terms, need O(alf^0 ... alf^6) on top
! f(0,0)=d^3, f(1,0)=f(0,1)=d^2,  f(2,0)=f(1,1)=f(0,2)=d

! h^2/d^2 terms
      c(jind(5,0,0)) =f(0,0) 
      c(jind(5,1,0)) =f(1,0) 
      c(jind(5,0,1)) =f(0,1) 
      c(jind(5,2,0)) =f(2,0) 
      c(jind(5,1,1)) =f(1,1) 
      c(jind(5,0,2)) =f(0,2) 
      c(jind(5,3,0)) =f(3,0) 
      c(jind(5,2,1)) =f(2,1) 
      c(jind(5,1,2)) =f(1,2) 
      c(jind(5,0,3)) =f(0,3) 
! h^2/d terms
      c(jind(5,4,0)) =f(4,0) 
      c(jind(5,3,1)) =f(3,1) 
      c(jind(5,2,2)) =f(2,2) 
      c(jind(5,1,3)) =f(1,3) 
      c(jind(5,0,4)) =f(0,4) 
! h^2 terms
      c(jind(5,5,0)) =f(5,0) 
      c(jind(5,4,1)) =f(4,1) 
      c(jind(5,3,2)) =f(3,2) 
      c(jind(5,2,3)) =f(2,3) 
      c(jind(5,1,4)) =f(1,4) 
      c(jind(5,0,5)) =f(0,5) 
! h^2d terms
      c(jind(5,6,0)) =f(6,0) 
      c(jind(5,5,1)) =f(5,1) 
      c(jind(5,4,2)) =f(4,2) 
      c(jind(5,3,3)) =f(3,3) 
      c(jind(5,2,4)) =f(2,4) 
      c(jind(5,1,5)) =f(1,5) 
      c(jind(5,0,6)) =f(0,6) 

! Coefficients of 1/r7 terms, need O(alf^3 ... alf^8) on top
! h^2/d terms 
      c(jind(7,3,0)) = -5*e(3,0)*f(0,0)/2 
      c(jind(7,2,1)) = -5*e(2,1)*f(0,0)/2 
      c(jind(7,1,2)) = -5*e(1,2)*f(0,0)/2 
      c(jind(7,0,3)) = -5*e(0,3)*f(0,0)/2 

! h^2/d terms +h^2
      c(jind(7,4,0)) = -5*(e(4,0)*f(0,0)+e(3,0)*f(1,0))/2 
      c(jind(7,3,1)) = -5*(e(3,1)*f(0,0)+e(3,0)*f(0,1)+e(2,1)*f(1,0))/2 
      c(jind(7,2,2)) = -5*(e(2,2)*f(0,0)+e(2,1)*f(0,1)+e(1,2)*f(1,0))/2 
      c(jind(7,1,3)) = -5*(e(1,3)*f(0,0)+e(1,2)*f(0,1)+e(0,3)*f(1,0))/2 
      c(jind(7,0,4)) = -5*(e(0,4)*f(0,0)+e(0,3)*f(0,1))/2 

! h^2/d terms +h^2+h^2d
      c(jind(7,5,0)) = -5*(e(5,0)*f(0,0)+e(4,0)*f(1,0)+e(3,0)*f(2,0))/2
      c(jind(7,4,1)) = -5*(e(4,1)*f(0,0)+e(4,0)*f(0,1)+e(3,1)*f(1,0)  &
                  +e(3,0)*f(1,1)+e(2,1)*f(2,0))/2 
      c(jind(7,3,2)) = -5*(e(3,2)*f(0,0)+e(3,1)*f(0,1)+e(3,0)*f(0,2)  &
                  +e(2,2)*f(1,0)+e(2,1)*f(1,1)+e(1,2)*f(2,0))/2
      c(jind(7,2,3)) = -5*(e(2,3)*f(0,0)+e(2,2)*f(0,1)+e(2,1)*f(0,2)  &
                  +e(1,3)*f(1,0)+e(1,2)*f(1,1)+e(0,3)*f(2,0))/2 
      c(jind(7,1,4)) = -5*(e(1,4)*f(0,0)+e(1,3)*f(0,1)+e(1,2)*f(0,2)  &
                  +e(0,4)*f(1,0)+e(0,3)*f(1,1))/2 
      c(jind(7,0,5)) = -5*(e(0,5)*f(0,0)+e(0,4)*f(0,1)+e(0,3)*f(0,2))/2 

! h^2/d terms +h^2+h^2d
      c(jind(7,6,0)) = -5*(e(5,0)*f(1,0)+e(4,0)*f(2,0)+e(3,0)*f(3,0))/2 
      c(jind(7,5,1)) = -5*(e(5,0)*f(0,1)+e(4,1)*f(1,0)+e(4,0)*f(1,1)  &
                  +e(3,1)*f(2,0)+e(3,0)*f(2,1)+e(2,1)*f(3,0))/2
      c(jind(7,4,2)) = -5*(e(4,1)*f(0,1)+e(4,0)*f(0,2)+e(3,2)*f(1,0)  &
                  +e(3,1)*f(1,1)+e(3,0)*f(1,2)+e(2,2)*f(2,0)  &
                  +e(2,1)*f(2,1)+e(1,2)*f(3,0))/2
      c(jind(7,3,3)) = -5*(e(3,2)*f(0,1)+e(3,1)*f(0,2)+e(3,0)*f(0,3)  &
                  +e(2,3)*f(1,0)+e(2,2)*f(1,1)+e(2,1)*f(1,2)  &
                  +e(1,3)*f(2,0)+e(1,2)*f(2,1)+e(0,3)*f(3,0))/2  
      c(jind(7,2,4)) = -5*(e(2,3)*f(0,1)+e(2,2)*f(0,2)+e(2,1)*f(0,3)  &
                  +e(1,4)*f(1,0)+e(1,3)*f(1,1)+e(1,2)*f(1,2)  &
                  +e(0,4)*f(2,0)+e(0,3)*f(2,1))/2
      c(jind(7,1,5)) = -5*(e(1,4)*f(0,1)+e(1,3)*f(0,2)+e(1,2)*f(0,3)  &
                  +e(0,5)*f(1,0)+e(0,4)*f(1,1)+e(0,3)*f(1,2))/2  
      c(jind(7,0,6)) = -5*(e(0,5)*f(0,1)+e(0,4)*f(0,2)+e(0,3)*f(0,3))/2

#ifdef THIRD
! h^2 terms +h^2d
      c(jind(7,7,0)) = -5*(e(5,0)*f(2,0)+e(4,0)*f(3,0)+e(3,0)*f(4,0))/2
      c(jind(7,6,1)) = -5*(e(5,0)*f(1,1)+e(4,1)*f(2,0)+e(4,0)*f(2,1)  &
                  +e(3,1)*f(3,0)+e(3,0)*f(3,1)+e(2,1)*f(4,0))/2
      c(jind(7,5,2)) = -5*(e(5,0)*f(0,2)+e(4,1)*f(1,1)+e(4,0)*f(1,2)  &
                  +e(3,2)*f(2,0)+e(3,1)*f(2,1)+e(3,0)*f(2,2)  &
                  +e(2,2)*f(3,0)+e(2,1)*f(3,1)+e(1,2)*f(4,0))/2
      c(jind(7,4,3)) = -5*(e(4,1)*f(0,2)+e(4,0)*f(0,3)+e(3,2)*f(1,1)  &
                  +e(3,1)*f(1,2)+e(3,0)*f(1,3)+e(2,3)*f(2,0)  &
                  +e(2,2)*f(2,1)+e(2,1)*f(2,2)+e(1,3)*f(3,0)  &
                  +e(1,2)*f(3,1)+e(0,3)*f(4,0))/2
      c(jind(7,3,4)) = -5*(e(3,2)*f(0,2)+e(3,1)*f(0,3)+e(3,0)*f(0,4)  &
                  +e(2,3)*f(1,1)+e(2,2)*f(1,2)+e(2,1)*f(1,3)  &
                  +e(1,4)*f(2,0)+e(1,3)*f(2,1)+e(1,2)*f(2,2)  &
                  +e(0,4)*f(3,0)+e(0,3)*f(3,1))/2
      c(jind(7,2,5)) = -5*(e(2,3)*f(0,2)+e(2,2)*f(0,3)+e(2,1)*f(0,4)  &
                  +e(1,4)*f(1,1)+e(1,3)*f(1,2)+e(1,2)*f(1,3)  &
                  +e(0,5)*f(2,0)+e(0,4)*f(2,1)+e(0,3)*f(2,2))/2
      c(jind(7,1,6)) = -5*(e(1,4)*f(0,2)+e(1,3)*f(0,3)+e(1,2)*f(0,4)  &
                  +e(0,5)*f(1,1)+e(0,4)*f(1,2)+e(0,3)*f(1,3))/2
      c(jind(7,0,7)) = -5*(e(0,5)*f(0,2)+e(0,4)*f(0,3)+e(0,3)*f(0,4))/2

! h^2d terms 
      c(jind(7,8,0)) = -5*(e(5,0)*f(3,0)+e(4,0)*f(4,0)+e(3,0)*f(5,0))/2
      c(jind(7,7,1)) = -5*(e(5,0)*f(2,1)+e(4,1)*f(3,0)+e(4,0)*f(3,1)  &
                  +e(3,1)*f(4,0)+e(3,0)*f(4,1)+e(2,1)*f(5,0))/2
      c(jind(7,6,2)) = -5*(e(5,0)*f(1,2)+e(4,1)*f(2,1)+e(4,0)*f(2,2)  &
                  +e(3,2)*f(3,0)+e(3,1)*f(3,1)+e(3,0)*f(3,2)  &
                  +e(2,2)*f(4,0)+e(2,1)*f(4,1)+e(1,2)*f(5,0))/2
      c(jind(7,5,3)) = -5*(e(5,0)*f(0,3)+e(4,1)*f(1,2)+e(4,0)*f(1,3)  &
                  +e(3,2)*f(2,1)+e(3,1)*f(2,2)+e(3,0)*f(2,3)  &
                  +e(2,3)*f(3,0)+e(2,2)*f(3,1)+e(2,1)*f(3,2)  &
                  +e(1,3)*f(4,0)+e(1,2)*f(4,1)+e(0,3)*f(5,0))/2
      c(jind(7,4,4)) = -5*(e(4,1)*f(0,3)+e(4,0)*f(0,4)+e(3,2)*f(1,2)  &
                  +e(3,1)*f(1,3)+e(3,0)*f(1,4)+e(2,3)*f(2,1)  &
                  +e(2,2)*f(2,2)+e(2,1)*f(2,3)+e(1,4)*f(3,0)  &
                  +e(1,3)*f(3,1)+e(1,2)*f(3,2)+e(0,4)*f(4,0)  &
                  +e(0,3)*f(4,1))/2
      c(jind(7,3,5)) = -5*(e(3,2)*f(0,3)+e(3,1)*f(0,4)+e(3,0)*f(0,5)  &
                  +e(2,3)*f(1,2)+e(2,2)*f(1,3)+e(2,1)*f(1,4)  &
                  +e(1,4)*f(2,1)+e(1,3)*f(2,2)+e(1,2)*f(2,3)  &
                  +e(0,5)*f(3,0)+e(0,4)*f(3,1)+e(0,3)*f(3,2))/2
      c(jind(7,2,6)) = -5*(e(2,3)*f(0,3)+e(2,2)*f(0,4)+e(2,1)*f(0,5)  &
                  +e(1,4)*f(1,2)+e(1,3)*f(1,3)+e(1,2)*f(1,4)  &
                  +e(0,5)*f(2,1)+e(0,4)*f(2,2)+e(0,3)*f(2,3))/2
      c(jind(7,1,7)) = -5*(e(1,4)*f(0,3)+e(1,3)*f(0,4)+e(1,2)*f(0,5)  &
                  +e(0,5)*f(1,2)+e(0,4)*f(1,3)+e(0,3)*f(1,4))/2
      c(jind(7,0,8)) = -5*(e(0,5)*f(0,3)+e(0,4)*f(0,4)+e(0,3)*f(0,5))/2

! Coefficients of 1/r9 terms, need O(alf^6 ... alf^10) on top

! h^2 terms 
      c(jind(9,6,0)) = 35*e(3,0)**2*f(0,0)/8
      c(jind(9,5,1)) = 35*e(2,1)*e(3,0)*f(0,0)/4
      c(jind(9,4,2)) = 35*(e(2,1)**2+2*e(1,2)*e(3,0))*f(0,0)/8
      c(jind(9,3,3)) = 35*(e(1,2)*e(2,1)+e(0,3)*e(3,0))*f(0,0)/4
      c(jind(9,2,4)) = 35*(e(1,2)**2+2*e(0,3)*e(2,1))*f(0,0)/8
      c(jind(9,1,5)) = 35*e(0,3)*e(1,2)*f(0,0)/4
      c(jind(9,0,6)) = 35*e(0,3)**2*f(0,0)/8

! h^2 terms + h^2d
      c(jind(9,7,0)) = 35*e(3,0)*(2*e(4,0)*f(0,0)+e(3,0)*f(1,0))/8
      c(jind(9,6,1)) = 35*(2*e(2,1)*e(4,0)*f(0,0)+e(3,0)**2*f(0,1)  &
                  +2*e(3,0)*(e(3,1)*f(0,0)+e(2,1)*f(1,0)))/8
      c(jind(9,5,2)) = 35*(2*(e(2,2)*e(3,0)+e(2,1)*e(3,1)  &
                  +e(1,2)*e(4,0))*f(0,0)+2*e(2,1)*e(3,0)*f(0,1)  &
                  +(e(2,1)**2+2*e(1,2)*e(3,0))*f(1,0))/8
      c(jind(9,4,3)) = 35*((e(2,1)**2+2*e(1,2)*e(3,0))*f(0,1)  &
                  +2*(e(2,1)*e(2,2)*f(0,0)+e(1,3)*e(3,0)*f(0,0)  &
                  +e(1,2)*e(3,1)*f(0,0)+e(0,3)*e(4,0)*f(0,0)  &
                  +e(1,2)*e(2,1)*f(1,0)+e(0,3)*e(3,0)*f(1,0)))/8
      c(jind(9,3,4)) = 35*(2*(e(1,3)*e(2,1)+e(1,2)*e(2,2)+e(0,4)*e(3,0)  &
                  +e(0,3)*e(3,1))*f(0,0)+2*(e(1,2)*e(2,1)  &
                  +e(0,3)*e(3,0))*f(0,1)+(e(1,2)**2  &
                  +2*e(0,3)*e(2,1))*f(1,0))/8
      c(jind(9,2,5)) = 35*(2*(e(1,2)*e(1,3)+e(0,4)*e(2,1)  &
                  +e(0,3)*e(2,2))*f(0,0)+(e(1,2)**2  &
                  +2*e(0,3)*e(2,1))*f(0,1)+2*e(0,3)*e(1,2)*f(1,0))/8
      c(jind(9,1,6)) = 35*(2*e(0,4)*e(1,2)*f(0,0)  &
                  +e(0,3)*(2*e(1,3)*f(0,0)  &
                  +2*e(1,2)*f(0,1)+e(0,3)*f(1,0)))/8
      c(jind(9,0,7)) = 35*e(0,3)*(2*e(0,4)*f(0,0)+e(0,3)*f(0,1))/8

! h^2 terms + h^2d
      c(jind(9,8,0)) = 35*e(3,0)*(2*e(4,0)*f(1,0)+e(3,0)*f(2,0))/8
      c(jind(9,7,1)) = 35*(2*e(2,1)*e(4,0)*f(1,0)+e(3,0)**2*f(1,1)  &
                  +2*e(3,0)*(e(4,0)*f(0,1)+e(3,1)*f(1,0)  &
                  +e(2,1)*f(2,0)))/8
      c(jind(9,6,2)) = 35*(e(3,0)**2*f(0,2)+2*e(1,2)*e(4,0)*f(1,0)  &
                  +2*e(2,1)*(e(4,0)*f(0,1)+e(3,1)*f(1,0))  &
                  +e(2,1)**2*f(2,0)+2*e(3,0)*(e(3,1)*f(0,1)  &
                  +e(2,2)*f(1,0)+e(2,1)*f(1,1)+e(1,2)*f(2,0)))/8
      c(jind(9,5,3)) = 35*(2*e(2,2)*(e(3,0)*f(0,1)+e(2,1)*f(1,0))  &
                  +e(2,1)**2*f(1,1)+2*e(2,1)*(e(3,1)*f(0,1)  &
                  +e(3,0)*f(0,2)+e(1,2)*f(2,0))  &
                  +2*(e(1,2)*e(4,0)*f(0,1)+e(1,3)*e(3,0)*f(1,0)  &
                  +e(1,2)*e(3,1)*f(1,0)+e(0,3)*e(4,0)*f(1,0)  &
                  +e(1,2)*e(3,0)*f(1,1)+e(0,3)*e(3,0)*f(2,0)))/8
      c(jind(9,4,4)) = 35*(2*e(1,3)*e(3,0)*f(0,1)+2*e(1,2)*e(3,1)*f(0,1)  &
                  +2*e(0,3)*e(4,0)*f(0,1)+e(2,1)**2*f(0,2)  &
                  +2*e(1,2)*e(3,0)*f(0,2)+2*e(1,2)*e(2,2)*f(1,0)  &
                  +2*e(0,4)*e(3,0)*f(1,0)+2*e(0,3)*e(3,1)*f(1,0)  &
                  +2*e(0,3)*e(3,0)*f(1,1)+e(1,2)**2*f(2,0)  &
                  +2*e(2,1)*(e(2,2)*f(0,1)+e(1,3)*f(1,0)  &
                  +e(1,2)*f(1,1)+e(0,3)*f(2,0)))/8
      c(jind(9,3,5)) = 35*(2*e(1,3)*(e(2,1)*f(0,1)+e(1,2)*f(1,0))  &
                  +e(1,2)**2*f(1,1)+2*(e(0,4)*e(3,0)*f(0,1)  &
                  +e(0,3)*e(3,1)*f(0,1)+e(0,3)*e(3,0)*f(0,2)  &
                  +e(0,4)*e(2,1)*f(1,0)+e(0,3)*e(2,2)*f(1,0)  &
                  +e(0,3)*e(2,1)*f(1,1))+2*e(1,2)*(e(2,2)*f(0,1)  &
                  +e(2,1)*f(0,2)+e(0,3)*f(2,0)))/8
      c(jind(9,2,6)) = 35*(2*e(0,4)*e(2,1)*f(0,1)+e(1,2)**2*f(0,2)  &
                  +2*e(1,2)*(e(1,3)*f(0,1)+e(0,4)*f(1,0)  &
                  +e(0,3)*f(1,1))+e(0,3)*(2*e(2,2)*f(0,1)  &
                  +2*e(2,1)*f(0,2)+2*e(1,3)*f(1,0)+e(0,3)*f(2,0)))/8
      c(jind(9,1,7)) = 35*(2*e(0,4)*(e(1,2)*f(0,1)+e(0,3)*f(1,0))  &
                  +e(0,3)*(2*e(1,3)*f(0,1)+2*e(1,2)*f(0,2)  &
                  +e(0,3)*f(1,1)))/8
      c(jind(9,0,8)) = 35*e(0,3)*(2*e(0,4)*f(0,1)+e(0,3)*f(0,2))/8

! h^2 terms + h^2d
      c(jind(9,9,0)) = 35*e(3,0)*(2*e(4,0)*f(2,0)+e(3,0)*f(3,0))/8
      c(jind(9,8,1)) = 35*(2*e(2,1)*e(4,0)*f(2,0)+e(3,0)**2*f(2,1)  &
                  +2*e(3,0)*(e(4,0)*f(1,1)+e(3,1)*f(2,0)  &
                  +e(2,1)*f(3,0)))/8
      c(jind(9,7,2)) = 35*(e(3,0)**2*f(1,2)+2*e(1,2)*e(4,0)*f(2,0)  &
                  +2*e(2,1)*(e(4,0)*f(1,1)+e(3,1)*f(2,0))  &
                  +e(2,1)**2*f(3,0)+2*e(3,0)*(e(4,0)*f(0,2)  &
                  +e(3,1)*f(1,1)+e(2,2)*f(2,0)+e(2,1)*f(2,1)  &
                  +e(1,2)*f(3,0)))/8
      c(jind(9,6,3)) = 35*(e(3,0)**2*f(0,3)+2*(e(1,2)*e(4,0)*f(1,1)  &
                  +e(1,2)*e(3,1)*f(2,0)+e(0,3)*e(4,0)*f(2,0))  &
                  +e(2,1)**2*f(2,1)+2*e(3,0)*(e(3,1)*f(0,2)  &
                  +e(2,2)*f(1,1)+e(2,1)*f(1,2)+e(1,3)*f(2,0)  &
                  +e(1,2)*f(2,1)+e(0,3)*f(3,0))+2*e(2,1)*(e(4,0)*f(0,2)  &
                  +e(3,1)*f(1,1)+e(2,2)*f(2,0)+e(1,2)*f(3,0)))/8
      c(jind(9,5,4)) = 35*(2*e(1,2)*e(4,0)*f(0,2)+2*e(1,3)*e(3,0)*f(1,1)  &
                  +2*e(1,2)*e(3,1)*f(1,1)+2*e(0,3)*e(4,0)*f(1,1)  &
                  +e(2,1)**2*f(1,2)+2*e(1,2)*e(3,0)*f(1,2)  &
                  +2*e(0,4)*e(3,0)*f(2,0)+2*e(0,3)*e(3,1)*f(2,0)  &
                  +2*e(2,2)*(e(3,0)*f(0,2)+e(2,1)*f(1,1)+e(1,2)*f(2,0))  &
                  +2*e(0,3)*e(3,0)*f(2,1)+e(1,2)**2*f(3,0)  &
                  +2*e(2,1)*(e(3,1)*f(0,2)+e(3,0)*f(0,3)+e(1,3)*f(2,0)  &
                  +e(1,2)*f(2,1)+e(0,3)*f(3,0)))/8
      c(jind(9,4,5)) = 35*(2*e(1,2)*e(3,1)*f(0,2)+2*e(0,3)*e(4,0)*f(0,2)  &
                  +e(2,1)**2*f(0,3)+2*e(1,2)*e(3,0)*f(0,3)  &
                  +2*e(1,2)*e(2,2)*f(1,1)+2*e(0,4)*e(3,0)*f(1,1)  &
                  +2*e(0,3)*e(3,1)*f(1,1)+2*e(0,3)*e(3,0)*f(1,2)  &
                  +2*e(0,3)*e(2,2)*f(2,0)+2*e(1,3)*(e(3,0)*f(0,2)  &
                  +e(1,2)*f(2,0))+e(1,2)**2*f(2,1)  &
                  +2*e(2,1)*(e(2,2)*f(0,2)+e(1,3)*f(1,1)+e(1,2)*f(1,2)  &
                  +e(0,4)*f(2,0)+e(0,3)*f(2,1))  &
                  +2*e(0,3)*e(1,2)*f(3,0))/8
      c(jind(9,3,6)) = 35*(2*e(0,4)*e(3,0)*f(0,2)+2*e(0,3)*e(3,1)*f(0,2)  &
                  +2*e(0,3)*e(3,0)*f(0,3)+2*e(0,4)*e(2,1)*f(1,1)  &
                  +2*e(0,3)*e(2,2)*f(1,1)+e(1,2)**2*f(1,2)  &
                  +2*e(0,3)*e(2,1)*f(1,2)+2*e(1,3)*(e(2,1)*f(0,2)  &
                  +e(1,2)*f(1,1)+e(0,3)*f(2,0))+2*e(1,2)*(e(2,2)*f(0,2)  &
                  +e(2,1)*f(0,3)+e(0,4)*f(2,0)+e(0,3)*f(2,1))  &
                  +e(0,3)**2*f(3,0))/8
      c(jind(9,2,7)) = 35*(e(1,2)**2*f(0,3)+2*e(1,2)*(e(1,3)*f(0,2)  &
                  +e(0,4)*f(1,1)+e(0,3)*f(1,2))+2*e(0,4)*(e(2,1)*f(0,2)  &
                  +e(0,3)*f(2,0))+e(0,3)*(2*e(2,2)*f(0,2)  &
                  +2*e(2,1)*f(0,3)+2*e(1,3)*f(1,1)+e(0,3)*f(2,1)))/8
      c(jind(9,1,8)) = 35*(2*e(0,4)*(e(1,2)*f(0,2)+e(0,3)*f(1,1))  &
                  +e(0,3)*(2*e(1,3)*f(0,2)+2*e(1,2)*f(0,3)  &
                  +e(0,3)*f(1,2)))/8
      c(jind(9,0,9)) = 35*e(0,3)*(2*e(0,4)*f(0,2)+e(0,3)*f(0,3))/8
#endif

#ifdef FOURTH
! h^2d terms 
      c(jind(9,10,0)) = 35*e(3,0)*(2*e(4,0)*f(3,0)+e(3,0)*f(4,0))/8
      c(jind(9,9,1)) = 35*(2*e(2,1)*e(4,0)*f(3,0)+e(3,0)**2*f(3,1)  &
                  +2*e(3,0)*(e(4,0)*f(2,1)+e(3,1)*f(3,0)  &
                  +e(2,1)*f(4,0)))/8
      c(jind(9,8,2)) = 35*(e(3,0)**2*f(2,2)+2*e(1,2)*e(4,0)*f(3,0)  &
                  +2*e(2,1)*(e(4,0)*f(2,1)+e(3,1)*f(3,0))  &
                  +e(2,1)**2*f(4,0)+2*e(3,0)*(e(4,0)*f(1,2)  &
                  +e(3,1)*f(2,1)+e(2,2)*f(3,0)+e(2,1)*f(3,1)  &
                  +e(1,2)*f(4,0)))/8
      c(jind(9,7,3)) = 35*(e(3,0)**2*f(1,3)+2*(e(1,2)*e(4,0)*f(2,1)  &
                  +e(1,2)*e(3,1)*f(3,0)+e(0,3)*e(4,0)*f(3,0))  &
                  +e(2,1)**2*f(3,1)+2*e(3,0)*(e(4,0)*f(0,3)  &
                  +e(3,1)*f(1,2)+e(2,2)*f(2,1)+e(2,1)*f(2,2)  &
                  +e(1,3)*f(3,0)+e(1,2)*f(3,1)+e(0,3)*f(4,0))  &
                  +2*e(2,1)*(e(4,0)*f(1,2)+e(3,1)*f(2,1)+e(2,2)*f(3,0)  &
                  +e(1,2)*f(4,0)))/8
      c(jind(9,6,4)) = 35*(e(3,0)**2*f(0,4)+2*e(1,2)*e(4,0)*f(1,2)  &
                  +2*e(1,2)*e(3,1)*f(2,1)+2*e(0,3)*e(4,0)*f(2,1)  &
                  +e(2,1)**2*f(2,2)+2*e(1,2)*e(2,2)*f(3,0)  &
                  +2*e(0,3)*e(3,1)*f(3,0)+2*e(3,0)*(e(3,1)*f(0,3)  &
                  +e(2,2)*f(1,2)+e(2,1)*f(1,3)+e(1,3)*f(2,1)  &
                  +e(1,2)*f(2,2)+e(0,4)*f(3,0)+e(0,3)*f(3,1))  &
                  +e(1,2)**2*f(4,0)+2*e(2,1)*(e(4,0)*f(0,3)  &
                  +e(3,1)*f(1,2)+e(2,2)*f(2,1)+e(1,3)*f(3,0)  &
                  +e(1,2)*f(3,1)+e(0,3)*f(4,0)))/8
      c(jind(9,5,5)) = 35*(2*e(1,2)*e(4,0)*f(0,3)+2*e(1,3)*e(3,0)*f(1,2)  &
                  +2*e(1,2)*e(3,1)*f(1,2)+2*e(0,3)*e(4,0)*f(1,2)  &
                  +e(2,1)**2*f(1,3)+2*e(1,2)*e(3,0)*f(1,3)  &
                  +2*e(0,4)*e(3,0)*f(2,1)+2*e(0,3)*e(3,1)*f(2,1)  &
                  +2*e(0,3)*e(3,0)*f(2,2)+2*e(1,2)*e(1,3)*f(3,0)  &
                  +2*e(2,2)*(e(3,0)*f(0,3)+e(2,1)*f(1,2)+e(1,2)*f(2,1)  &
                  +e(0,3)*f(3,0))+e(1,2)**2*f(3,1)  &
                  +2*e(2,1)*(e(3,1)*f(0,3)+e(3,0)*f(0,4)  &
                  +e(1,3)*f(2,1)+e(1,2)*f(2,2)+e(0,4)*f(3,0)  &
                  +e(0,3)*f(3,1))+2*e(0,3)*e(1,2)*f(4,0))/8
      c(jind(9,4,6)) = 35*(2*e(1,2)*e(3,1)*f(0,3)+2*e(0,3)*e(4,0)*f(0,3)  &
                  +e(2,1)**2*f(0,4)+2*e(1,2)*e(3,0)*f(0,4)  &
                  +2*e(1,2)*e(2,2)*f(1,2)+2*e(0,4)*e(3,0)*f(1,2)  &
                  +2*e(0,3)*e(3,1)*f(1,2)+2*e(0,3)*e(3,0)*f(1,3)  &
                  +2*e(0,3)*e(2,2)*f(2,1)+e(1,2)**2*f(2,2)  &
                  +2*e(2,1)*(e(2,2)*f(0,3)+e(1,3)*f(1,2)  &
                  +e(1,2)*f(1,3)+e(0,4)*f(2,1)+e(0,3)*f(2,2))  &
                  +2*e(0,4)*e(1,2)*f(3,0)+2*e(1,3)*(e(3,0)*f(0,3)  &
                  +e(1,2)*f(2,1)+e(0,3)*f(3,0))+2*e(0,3)*e(1,2)*f(3,1)  &
                  +e(0,3)**2*f(4,0))/8
      c(jind(9,3,7)) = 35*(2*e(0,4)*e(3,0)*f(0,3)+2*e(0,3)*e(3,1)*f(0,3)  &
                  +2*e(0,3)*e(3,0)*f(0,4)+2*e(0,4)*e(2,1)*f(1,2)  &
                  +2*e(0,3)*e(2,2)*f(1,2)+e(1,2)**2*f(1,3)  &
                  +2*e(0,3)*e(2,1)*f(1,3)+2*e(1,3)*(e(2,1)*f(0,3)  &
                  +e(1,2)*f(1,2)+e(0,3)*f(2,1))+2*e(1,2)*(e(2,2)*f(0,3)  &
                  +e(2,1)*f(0,4)+e(0,4)*f(2,1)+e(0,3)*f(2,2))  &
                  +2*e(0,3)*e(0,4)*f(3,0)+e(0,3)**2*f(3,1))/8
      c(jind(9,2,8)) = 35*(e(1,2)**2*f(0,4)+2*e(1,2)*(e(1,3)*f(0,3)  &
                  +e(0,4)*f(1,2)+e(0,3)*f(1,3))+2*e(0,4)*(e(2,1)*f(0,3)  &
                  +e(0,3)*f(2,1))+e(0,3)*(2*e(2,2)*f(0,3)  &
                  +2*e(2,1)*f(0,4)+2*e(1,3)*f(1,2)+e(0,3)*f(2,2)))/8
      c(jind(9,1,9)) = 35*(2*e(0,4)*(e(1,2)*f(0,3)+e(0,3)*f(1,2))  &
                  +e(0,3)*(2*e(1,3)*f(0,3)+2*e(1,2)*f(0,4)  &
                  +e(0,3)*f(1,3)))/8
      c(jind(9,0,10)) = 35*e(0,3)*(2*e(0,4)*f(0,3)+e(0,3)*f(0,4))/8


! Coefficients of 1/r11 terms, need O(alf^9 ... alf^12) on top

! h^2d terms 
      c(jind(11,9,0)) = -105*e(3,0)**3*f(0,0)/16
      c(jind(11,8,1)) = -315*e(2,1)*e(3,0)**2*f(0,0)/16
      c(jind(11,7,2)) = -315*e(3,0)*(e(2,1)**2+e(1,2)*e(3,0))*f(0,0)/16
      c(jind(11,6,3)) = -105*(e(2,1)**3+6*e(1,2)*e(2,1)*e(3,0)  &
                  +3*e(0,3)*e(3,0)**2)*f(0,0)/16
      c(jind(11,5,4)) = -315*(e(1,2)*e(2,1)**2+e(1,2)**2*e(3,0)  &
                  +2*e(0,3)*e(2,1)*e(3,0))*f(0,0)/16
      c(jind(11,4,5)) = -315*(e(1,2)**2*e(2,1)+e(0,3)*e(2,1)**2  &
                  +2*e(0,3)*e(1,2)*e(3,0))*f(0,0)/16
      c(jind(11,3,6)) = -105*(e(1,2)**3+6*e(0,3)*e(1,2)*e(2,1)  &
                  +3*e(0,3)**2*e(3,0))*f(0,0)/16
      c(jind(11,2,7)) = -315*e(0,3)*(e(1,2)**2+e(0,3)*e(2,1))*f(0,0)/16
      c(jind(11,1,8)) = -315*e(0,3)**2*e(1,2)*f(0,0)/16
      c(jind(11,0,9)) = -105*e(0,3)**3*f(0,0)/16

! h^2d terms 
      c(jind(11,10,0)) = -105*e(3,0)**3*f(1,0)/16
      c(jind(11,9,1)) = -105*e(3,0)**2*(e(3,0)*f(0,1)  &
                  +3*e(2,1)*f(1,0))/16
      c(jind(11,8,2)) = -315*e(3,0)*(e(2,1)*e(3,0)*f(0,1)  &
                  +e(2,1)**2*f(1,0)  &
                  +e(1,2)*e(3,0)*f(1,0))/16
      c(jind(11,7,3)) = -105*(3*e(2,1)**2*e(3,0)*f(0,1)  &
                  +e(2,1)**3*f(1,0)  &
                  +6*e(1,2)*e(2,1)*e(3,0)*f(1,0)  &
                  +3*e(3,0)**2*(e(1,2)*f(0,1)+e(0,3)*f(1,0)))/16
      c(jind(11,6,4)) = -105*(e(2,1)**3*f(0,1)  &
                  +3*e(1,2)*e(2,1)**2*f(1,0)  &
                  +6*e(2,1)*e(3,0)*(e(1,2)*f(0,1)+e(0,3)*f(1,0))  &
                  +3*e(3,0)*(e(0,3)*e(3,0)*f(0,1)+e(1,2)**2*f(1,0)))/16
      c(jind(11,5,5)) = -315*(e(1,2)**2*(e(3,0)*f(0,1)+e(2,1)*f(1,0))  &
                  +e(0,3)*e(2,1)*(2*e(3,0)*f(0,1)+e(2,1)*f(1,0))  &
                  +e(1,2)*(e(2,1)**2*f(0,1)+2*e(0,3)*e(3,0)*f(1,0)))/16
      c(jind(11,4,6)) = -105*(3*e(1,2)**2*e(2,1)*f(0,1)  &
                  +e(1,2)**3*f(1,0)  &
                  +6*e(0,3)*e(1,2)*(e(3,0)*f(0,1)+e(2,1)*f(1,0))  &
                  +3*e(0,3)*(e(2,1)**2*f(0,1)+e(0,3)*e(3,0)*f(1,0)))/16
      c(jind(11,3,7)) = -105*(e(1,2)**3*f(0,1)  &
                  +6*e(0,3)*e(1,2)*e(2,1)*f(0,1)  &
                  +3*e(0,3)*e(1,2)**2*f(1,0)+3*e(0,3)**2*(e(3,0)*f(0,1)  &
                  +e(2,1)*f(1,0)))/16
      c(jind(11,2,8)) = -315*e(0,3)*(e(1,2)**2*f(0,1)  &
                  +e(0,3)*e(2,1)*f(0,1)  &
                  +e(0,3)*e(1,2)*f(1,0))/16
      c(jind(11,1,9)) = -105*e(0,3)**2*(3*e(1,2)*f(0,1)  &
                  +e(0,3)*f(1,0))/16
      c(jind(11,0,10)) = -105*e(0,3)**3*f(0,1)/16

! h^2d terms 
      c(jind(11,11,0)) = -105*e(3,0)**3*f(2,0)/16
      c(jind(11,10,1)) = -105*e(3,0)**2*(e(3,0)*f(1,1)  &
                  +3*e(2,1)*f(2,0))/16
      c(jind(11,9,2)) = -105*e(3,0)*(e(3,0)**2*f(0,2)  &
                  +3*e(2,1)**2*f(2,0)  &
                  +3*e(3,0)*(e(2,1)*f(1,1)+e(1,2)*f(2,0)))/16
      c(jind(11,8,3)) = -105*(3*e(2,1)**2*e(3,0)*f(1,1)  &
                  +e(2,1)**3*f(2,0)  &
                  +3*e(3,0)**2*(e(1,2)*f(1,1)+e(0,3)*f(2,0))  &
                  +3*e(2,1)*e(3,0)*(e(3,0)*f(0,2)+2*e(1,2)*f(2,0)))/16
      c(jind(11,7,4)) = -105*(e(2,1)**3*f(1,1)  &
                  +6*e(2,1)*e(3,0)*(e(1,2)*f(1,1)  &
                  +e(0,3)*f(2,0))+3*e(2,1)**2*(e(3,0)*f(0,2)  &
                  +e(1,2)*f(2,0))+3*e(3,0)*(e(1,2)*e(3,0)*f(0,2)  &
                  +e(0,3)*e(3,0)*f(1,1)+e(1,2)**2*f(2,0)))/16
      c(jind(11,6,5)) = -105*(e(2,1)**3*f(0,2)  &
                  +3*e(2,1)**2*(e(1,2)*f(1,1)  &
                  +e(0,3)*f(2,0))+3*e(3,0)*(e(0,3)*e(3,0)*f(0,2)  &
                  +e(1,2)**2*f(1,1)+2*e(0,3)*e(1,2)*f(2,0))  &
                  +3*e(2,1)*(2*e(1,2)*e(3,0)*f(0,2)  &
                  +2*e(0,3)*e(3,0)*f(1,1)+e(1,2)**2*f(2,0)))/16
      c(jind(11,5,6)) = -105*(3*e(1,2)**2*(e(3,0)*f(0,2)  &
                  +e(2,1)*f(1,1))  &
                  +e(1,2)**3*f(2,0)+3*e(1,2)*(e(2,1)**2*f(0,2)  &
                  +2*e(0,3)*e(3,0)*f(1,1)+2*e(0,3)*e(2,1)*f(2,0))  &
                  +3*e(0,3)*(2*e(2,1)*e(3,0)*f(0,2)+e(2,1)**2*f(1,1)  &
                  +e(0,3)*e(3,0)*f(2,0)))/16
      c(jind(11,4,7)) = -105*(e(1,2)**3*f(1,1)  &
                  +6*e(0,3)*e(1,2)*(e(3,0)*f(0,2)  &
                  +e(2,1)*f(1,1))+3*e(1,2)**2*(e(2,1)*f(0,2)  &
                  +e(0,3)*f(2,0))+3*e(0,3)*(e(2,1)**2*f(0,2)  &
                  +e(0,3)*e(3,0)*f(1,1)+e(0,3)*e(2,1)*f(2,0)))/16
      c(jind(11,3,8)) = -105*(e(1,2)**3*f(0,2)  &
                  +3*e(0,3)*e(1,2)**2*f(1,1)  &
                  +3*e(0,3)**2*(e(3,0)*f(0,2)+e(2,1)*f(1,1))  &
                  +3*e(0,3)*e(1,2)*(2*e(2,1)*f(0,2)+e(0,3)*f(2,0)))/16
      c(jind(11,2,9)) = -105*e(0,3)*(3*e(1,2)**2*f(0,2)  &
                  +3*e(0,3)*e(1,2)*f(1,1)+e(0,3)*(3*e(2,1)*f(0,2)  &
                  +e(0,3)*f(2,0)))/16
      c(jind(11,1,10)) = -105*e(0,3)**2*(3*e(1,2)*f(0,2)  &
                  +e(0,3)*f(1,1))/16
      c(jind(11,0,11)) = -105*e(0,3)**3*f(0,2)/16

! h^2d terms 
      c(jind(11,12,0)) = -105*e(3,0)**3*f(3,0)/16
      c(jind(11,11,1)) = -105*e(3,0)**2*(e(3,0)*f(2,1)  &
                  +3*e(2,1)*f(3,0))/16
      c(jind(11,10,2)) = -105*e(3,0)*(e(3,0)**2*f(1,2)  &
                  +3*e(2,1)**2*f(3,0)  &
                  +3*e(3,0)*(e(2,1)*f(2,1)+e(1,2)*f(3,0)))/16
      c(jind(11,9,3)) = -105*(e(3,0)**3*f(0,3)+e(2,1)**3*f(3,0)  &
                  +3*e(3,0)**2*(e(2,1)*f(1,2)+e(1,2)*f(2,1)  &
                  +e(0,3)*f(3,0))+3*e(2,1)*e(3,0)*(e(2,1)*f(2,1)  &
                  +2*e(1,2)*f(3,0)))/16
      c(jind(11,8,4)) = -105*(e(2,1)**3*f(2,1)  &
                  +3*e(2,1)*e(3,0)*(e(3,0)*f(0,3)  &
                  +2*e(1,2)*f(2,1)+2*e(0,3)*f(3,0))  &
                  +3*e(2,1)**2*(e(3,0)*f(1,2)+e(1,2)*f(3,0))  &
                  +3*e(3,0)*(e(1,2)*e(3,0)*f(1,2)+e(0,3)*e(3,0)*f(2,1)  &
                  +e(1,2)**2*f(3,0)))/16
      c(jind(11,7,5)) = -105*(e(2,1)**3*f(1,2)  &
                  +3*e(2,1)**2*(e(3,0)*f(0,3)  &
                  +e(1,2)*f(2,1)+e(0,3)*f(3,0))  &
                  +3*e(3,0)*(e(1,2)*e(3,0)*f(0,3)+e(0,3)*e(3,0)*f(1,2)  &
                  +e(1,2)**2*f(2,1)+2*e(0,3)*e(1,2)*f(3,0))  &
                  +3*e(2,1)*(2*e(1,2)*e(3,0)*f(1,2)  &
                  +2*e(0,3)*e(3,0)*f(2,1)+e(1,2)**2*f(3,0)))/16
      c(jind(11,6,6)) = -105*(e(2,1)**3*f(0,3)  &
                  +3*e(2,1)**2*(e(1,2)*f(1,2)  &
                  +e(0,3)*f(2,1))+3*e(0,3)*e(3,0)*(e(3,0)*f(0,3)  &
                  +2*e(1,2)*f(2,1))+3*e(0,3)**2*e(3,0)*f(3,0)  &
                  +e(1,2)**2*(3*e(3,0)*f(1,2)+e(1,2)*f(3,0))  &
                  +3*e(2,1)*(2*e(0,3)*e(3,0)*f(1,2)+e(1,2)**2*f(2,1)  &
                  +2*e(1,2)*(e(3,0)*f(0,3)+e(0,3)*f(3,0))))/16
      c(jind(11,5,7)) = -105*(e(1,2)**3*f(2,1)  &
                  +3*e(1,2)*(e(2,1)**2*f(0,3)  &
                  +2*e(0,3)*e(3,0)*f(1,2)+2*e(0,3)*e(2,1)*f(2,1))  &
                  +3*e(1,2)**2*(e(3,0)*f(0,3)+e(2,1)*f(1,2)  &
                  +e(0,3)*f(3,0))+3*e(0,3)*(2*e(2,1)*e(3,0)*f(0,3)  &
                  +e(2,1)**2*f(1,2)+e(0,3)*e(3,0)*f(2,1)  &
                  +e(0,3)*e(2,1)*f(3,0)))/16
      c(jind(11,4,8)) = -105*(e(1,2)**3*f(1,2)  &
                  +3*e(1,2)**2*(e(2,1)*f(0,3)  &
                  +e(0,3)*f(2,1))+3*e(0,3)*(e(2,1)**2*f(0,3)  &
                  +e(0,3)*e(3,0)*f(1,2)+e(0,3)*e(2,1)*f(2,1))  &
                  +3*e(0,3)*e(1,2)*(2*e(3,0)*f(0,3)+2*e(2,1)*f(1,2)  &
                  +e(0,3)*f(3,0)))/16
      c(jind(11,3,9)) = -105*(e(1,2)**3*f(0,3)  &
                  +3*e(0,3)*e(1,2)**2*f(1,2)  &
                  +3*e(0,3)*e(1,2)*(2*e(2,1)*f(0,3)+e(0,3)*f(2,1))  &
                  +e(0,3)**2*(3*e(3,0)*f(0,3)+3*e(2,1)*f(1,2)  &
                  +e(0,3)*f(3,0)))/16
      c(jind(11,2,10)) = -105*e(0,3)*(3*e(1,2)**2*f(0,3)  &
                  +3*e(0,3)*e(1,2)*f(1,2)+e(0,3)*(3*e(2,1)*f(0,3)  &
                  +e(0,3)*f(2,1)))/16
      c(jind(11,1,11)) = -105*e(0,3)**2*(3*e(1,2)*f(0,3)  &
                  +e(0,3)*f(1,2))/16
      c(jind(11,0,12)) = -105*e(0,3)**3*f(0,3)/16
#endif

! 2nd order method (all h^2/d^2, h^2/d terms)
! r=5, p+q=0,1,2,3,4
! r=7, p+q=3,4,5,6
      sumo=0
      r=5
      do j=0,4
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      enddo
          j=0
          p=0
          j0=jind(r,p,j-p)
!          print*,r,p,j-p,j0,c(j0),corr(j0),c(j0)*corr(j0)
      r=7
      do j=3,6
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
!          print*,r,p,j-p,c(j0),corr(j0),c(j0)*corr(j0)
        enddo
      enddo
#ifdef THIRD
! 3rd order method (all h^2 terms)
! r=5, p+q=5
! r=7, p+q=7
! r=9, p+q=6,7,8,9
      r=5
      j=5
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      r=7
      j=7
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      r=9
      do j=6,9
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      enddo
#endif
#ifdef FOURTH
! 4th order method (all h^2 d terms)
! r=5, p+q=6
! r=7, p+q=8
! r=9, p+q=10
! r=11,p+q=9,10,11,12
      r=5
      j=6
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      r=7
      j=8
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      r=9
      j=10
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      r=11
      do j=9,12
        do p=0,j
          j0=jind(r,p,j-p)
          sumo = sumo + c(j0)*corr(j0)
        enddo
      enddo
#endif
      comprho5=sumo
      return
END FUNCTION comprho5

END MODULE mod_stokes
