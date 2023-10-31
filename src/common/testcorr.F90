PROGRAM testcorr
! Tests corrections computed in 
! compallcorrtrap in mod_EHpqr
! by outputing the corrections for a fixed point, fixed grid!
! and checking that they converge properly as gridsize n increases
!
! The point was chosen to be one which gave bad result in initial velocity 
! error checks along ray in off direction with c.ne.0
!
! 1) First print only corr(j)*t%icorr (so =0 when not corrected)
!    to check errors decreasing  GOOD!!
! 2) Then print ratios. Method is O(h^6) to reduce effect of boundaries
!
! RESULTS for this test case: 
! 1) errors decrease to 10^-12 to 10^-17 in all cases for n=40,80,160
! 2) Ratios should be 64. They are generally much larger than that
!    except for a few cases in the first step (n=40)
!
! CONCLUSION: PASSED!!!

use params, ONLY : nbmax,rmax,pmax,pi
use types, ONLY : basept
use globalvars, ONLY :ra,rb,rc,g,correct,axmax
use mod_geom  ! for setgrid
use mod_EHpqr 
implicit none
TYPE (basept) :: t
integer, dimension(rmax,0:pmax,0:pmax) :: jind
real*8, dimension(nbmax) :: corr,oldcorr
real*8, dimension(3) :: x0
real*8 :: calf,cbet,calfbet,h
integer :: nb,mstep,l,k,n,j

correct=.true.
axmax=3;

! first test case
!ra=3; rb=2; rc=1; 
!x0=(/-2.4182027985955603d0, 1.8144143528329761d0, 6.4237221566660702D-002 /)
!calf=2.7621055544734592d0
!cbet=1.4775068744576145d0     
!calfbet=-7.3684563742005785D-002

! second test case: maxerror for identity on 132ellipse with f=0,1,0
ra=1; rb=3; rc=2; 
x0=(/6.4377368797743702D-017,  2.2544930010583437d0 ,1.5656556257032237D0 /)
calf=  0.86642251357966082d0; cbet=2.6312319867946359d0; calfbet= -2.4488214736414174D-016;

!x0=(/ cos(alf)*cos(bet),sin(alf)*cos(bet),sin(bet) /)
!x0=x0*(1+d)


nb=179
do mstep=0,17
!do mstep=0,1
!do mstep=0,0
  oldcorr=0
  
!  do k=0,4
  do k=0,3 ! in this testcase n=320,640 are not corrected, pt too far
     n=40*2**k
     call setgrid(n)
     call setbase(x0,t)
     h=g(1)%h
     l=t%igrid
     call compallcorrtrap(t,g(l),calf,cbet,calfbet,corr,jind)
!     write(*,'(i4,10f13.3)')mstep,(oldcorr(j)/corr(j),j=1+mstep*10,min(10+mstep*10,nb))
    write(*,'(i4,10e13.3)')mstep,(corr(j)*t%icorr,j=1+mstep*10,min(10+mstep*10,nb))
     oldcorr=corr
  enddo
print*
enddo

stop
end
