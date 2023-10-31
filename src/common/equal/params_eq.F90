MODULE params
! contains all global parameters
implicit none

! constant parameters
real*8, parameter :: pi= 3.14159265358979324d0,pi2=2*pi,pi4=4*pi,pi8=8*pi

! max gridsizes
integer, parameter :: ngrid=800, nwinmax=2*150  ! at least 2*nwin!!!!
integer, parameter :: mgrid=ngrid/2, nleft=-nwinmax, nright=ngrid+nwinmax

! max number of target points
integer, parameter :: nvecmax=60000  !for 75^2 =5625 grid
!integer, parameter :: nvecmax=105000  !for 75^2 =5625 grid

! FOR CORRECTIONS:
! number and degree of basis functions
integer, parameter :: nbmax=185, rmax=11, pmax=12
! degree of polynomial in numerator
integer, parameter :: cmax=6

END MODULE params
