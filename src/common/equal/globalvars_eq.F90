MODULE globalvars
use params, ONLY: cmax
use types, ONLY : grid
implicit none

logical :: correct
character*5 :: label
!real*8, dimension(3) :: uinf
real*8 :: ra,rb,rc,axmax ! 3 axes
real*8, dimension(3) :: uinf


TYPE(grid) :: g(2)

real*8 :: tmax,delt,time,xfar,h

real*8 :: calf,cbet,calfbet
real*8,dimension(0:cmax,0:cmax) :: cx,cy,cz,cj,e,cnjac1,cnjac2,cnjac3
real*8,dimension(0:cmax,0:cmax) :: cf1,cf2,cf3
END MODULE globalvars
