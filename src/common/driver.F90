PROGRAM driver
! applies to any geometry. Uses
! MODULES
! mod_geom  (sphere, ellipse, multi-ellipse) replaces mod_sphere
!       setgrid
!       velo
!       setbase
!       compcoeffx
! mod_init  (case specific: sphere, ellipse, multi-ellipse)
!       for ellipse: sets uinf, ra,rb,rc,filein,datout
!       for sphere: sets uinf, datout
!       sets density matching uinf,ra,rb,rc
!       sets coeffdensityex when available
! mod_target (sphere, ellipse, multi-ellipse)
!       targetcross (sets target for crossection)
!       targetpatchatd (sets target at fixed d over patch, for max)
!       targetshadow (sets target for shadow)
!       targetstreak (sets target for streaklines)
!       targetstream (sets target for streamlines)
! mod_SDLP
! mod_stokes
!      compsdlp
!      compallcorrtrap
!      correctg, correctt
!      compcoeffdens (by interpolation)
!
! reads in problem number ipb
!  if 1 : fixed: comp velo at distance d, everywhere, and print
!  if 2 : fixed: comp velo at distance d, on a patch, print all to find max
!  if 3 : fixed: comp dble at distance d, small domain, and print max, use T11, x0 inside
!  if 4 : fixed: comp velo in crossections, and print
!  if 5 : streamlines: evolve with runge and print every x steps 
!  if 6 : streaklines: evolve with runge and print every x steps 
!  if 7 : shadow: evolve with runge, check for end, print interpolant
!
use globalvars, ONLY : correct,label
use mod_init  !for init (uinf, ra,rb,rc),setdens,
!use mod_geom  !for setgrid
implicit none
integer m1,n1,m2,n2,icorrect,ipb

read*,m1,n1;  !print*,'gridsize n=',n0
read*,m2,n2;  !print*,'gridsize n=',n0
read*,icorrect
read*,ipb
correct=.false.; label='orig_'
if (icorrect.eq.1) then 
   correct=.true.; label='corr_'
endif

call init(m1,n1,m2,n2)       !initializes: uinf
                    !ra,rb,rc,axmax (ellipse only)
                    !grids g(1),g(2)
                    !density of grids

if (ipb.eq.4) then 
  call crossection !set target pts, outfile, compute output
else if (ipb.eq.5) then 
  call ray !set target pts, outfile, compute output
!if (ipb.eq.7) then 
!  call shadow
else if (ipb.eq.2) then 
  call patch
!else if (ipb.eq.4) then 
!  call fixeddble
!else if (ipb.eq.5) then 
!  call evolve(5)
!else if (ipb.eq.6) then 
!  call evolve(6)
!else 
!  call fixed(ipb)
endif
!write(7,'(i5,a)')n0,' n0 (size of grid)'
!call runtime(7)

stop
END PROGRAM driver


SUBROUTINE crossection
use params, ONLY : nvecmax
use globalvars, ONLY : label,xfar
use mod_target ! target_cross
use mod_geom   ! for velo
implicit none
!Local
integer :: nvec,j
character*20 :: flend
real*8 :: htarget,top,right
real*8,dimension(0:nvecmax,3) :: xvec,u

read*
read*,flend
read*,htarget
read*,top,right

xfar =sqrt(top**2+right**2)

open (7,file='datout_cross/info'//label//flend)
open (8,file='datout_cross/velofld'//label//flend)
call target_cross(htarget,top,right,nvec,xvec)
call velo(nvec,xvec,u)
do j=0,nvec
  write(8,*)xvec(j,:),u(j,:)
enddo
END SUBROUTINE crossection


SUBROUTINE ray
use params, ONLY : nvecmax
use globalvars, ONLY : label,xfar
use mod_target ! target_cross
use mod_geom   ! for velo
implicit none
!Local
integer :: nvec,j
character*20 :: flend
real*8,dimension(0:nvecmax,3) :: xvec,u

read*
read*,flend

xfar =10

!open (7,file='datout_ray/inf_'//label//flend)
open (8,file='datout_ray/vel_'//label//flend)
call target_ray(nvec,xvec)
call velo(nvec,xvec,u)
do j=0,nvec
  write(8,*)xvec(j,:),u(j,:)
enddo
END SUBROUTINE ray


SUBROUTINE patch
use params, ONLY : nvecmax
use globalvars, ONLY : label,xfar
use mod_target ! target_cross
implicit none
!Local
integer :: nvec,j,i1,i2,l,patchid
real*8 :: d
real*8,dimension(0:nvecmax,3) :: xvec,u
character*20 :: flend
character*3 :: dlab

read*
read*,flend
xfar =10

!xvec(0,:)=(/ -2.2419955818511128d0, 1.5545424394345515d0, 4.2365148119333791D-002 /)
!print*,xvec(0,:)
!call velo(0,xvec,u)
!print*,u(0,:)

!do l=0,11
do l=2,2
!do l=11,11
  d=0.2d0/2**l
!d=1.6d0
  i2 = mod(l,10)
  i1 = (l - i2)/10
  write(dlab,'(a,2i1,a)')'d',i1,i2
!  print*,dlab
  open (7,file='datout_patch/info'//dlab//label//flend)
  open (8,file='datout_patch/vel'//dlab//label//flend)
  write(7,'(f12.9,a,i3,a)')d,'=d'
  call target_patch(d) 
!  call velo(nvec,xvec,u)
!  do j=0,nvec
!    write(8,*)xvec(j,:),u(j,:)
!  enddo
enddo
END SUBROUTINE patch
