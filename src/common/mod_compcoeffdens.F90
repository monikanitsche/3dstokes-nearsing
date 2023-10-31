MODULE mod_compcoeffdens
! contains all routines needed to compute the 
! coefficients of the density at the basepoint
! assumes that f1,f2,f3 and all derivatives 
!  f1der(i,j,1:9), f2der(i,j,1:9), f3der(i,j,1:9) 
! have been initialized

PUBLIC

CONTAINS

SUBROUTINE compcoeffdens(g,t)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ONLY FOR GRID 1
! TEST CASE f=n
! for single layer density f=sigma*n=(f1,f2,f3)
! for double layer density f=u=(u,v,w)
! in test case single layer: f=n
! in test case double layer: u=(1,0,0), (0,1,0) or (0,0,1)
!
! in general: 
!    input f at gridpoints, interpolate it and derivatives to alfb,betb
!    output up to 3rd derivatives (for double layer)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
use types
use globalvars, ONLY : cf1,cf2,cf3
implicit none
TYPE(grid)::g
TYPE(basept)::t
integer :: j,p,k

  call comp_f_interpo(g,t)

!do j=0,3
!do p=0,j
!  k=j-p
!print*,p,k,cf1(p,k)
!enddo
!enddo
        
  ! change from derivatives to coefficients
  cf1(2,0) = cf1(2,0)/2   
  cf2(2,0) = cf2(2,0)/2
  cf3(2,0) = cf3(2,0)/2
  cf1(0,2) = cf1(0,2)/2
  cf2(0,2) = cf2(0,2)/2
  cf3(0,2) = cf3(0,2)/2
  
  cf1(3,0) =  cf1(3,0)/6
  cf2(3,0) =  cf2(3,0)/6
  cf3(3,0) =  cf3(3,0)/6
        
  cf1(2,1) =  cf1(2,1)/2
  cf2(2,1) =  cf2(2,1)/2
  cf3(2,1) =  cf3(2,1)/2
        
  cf1(1,2) =  cf1(1,2)/2
  cf2(1,2) =  cf2(1,2)/2
  cf3(1,2) =  cf3(1,2)/2
        
  cf1(0,3) =  cf1(0,3)/6
  cf2(0,3) =  cf2(0,3)/6
  cf3(0,3) =  cf3(0,3)/6

!print*,t%alfb,t%betb,cf3(0,0),cf3(1,0),cf3(0,1),cf3(2,0)
!print*,cf3(1,1),cf3(0,2),cf3(3,0),cf3(2,1),cf3(1,2),cf3(0,3)
return
END SUBROUTINE compcoeffdens


SUBROUTINE comp_f_interpo(g,t)
  ! Finds bicubi interpolant of densities
  use params, only: pi
  use globalvars, only:cf1,cf2,cf3
  use types
  implicit none
  TYPE(grid), INTENT(IN) :: g
  TYPE(basept), INTENT(IN) :: t
  integer :: ik,k,i0,j0
  integer, parameter :: n=16
  integer :: ind1(n),ind2(n)
  real*8 :: alfa_value(n), beta_value(n),f_value(n)
  real*8 :: h,alfb,betb
  real*8 :: comp_value(9)

  alfb=t%alfb
  betb=t%betb
  i0=t%i0
  j0=t%j0
  h=g%h
!print*,'in compfinterp'
!print*,alfb,i0, -pi+i0*h
!print*,betb,i0, -pi/2+j0*h
!print*,g%f1(i0,j0-2)
!print*,g%f1(i0,j0-1)
!print*,g%f1(i0,j0)
!print*,g%f1(i0,j0+1)
!print*,g%f1(i0,j0+2)
!print*,g%f1(i0-2,j0)
!print*,g%f1(i0-1,j0)
!!print*,g%f1(i0,j0)
!print*,g%f1(i0+1,j0)
!print*,g%f1(i0+2,j0)
!print*
!print*,( g%f3(i0+1,j0+1) -g%f3(i0+1,j0-1) -g%f3(i0-1,j0+1) +g%f3(i0-1,j0-1) )/h**2/4 !11 derivative
!print*,( g%f2(i0+1,j0) -g%f2(i0-1,j0) )/(2*h) !20 derivative
!print*,( g%f2(i0+1,j0+1) -2*g%f2(i0+1,j0)+g%f2(i0+1,j0-1) &
!- ( g%f2(i0-1,j0+1) -2*g%f2(i0-1,j0)+g%f2(i0-1,j0-1) ) )/(2*h**3) !21
!print*,( g%f2(i0+2,j0) -2*g%f2(i0+1,j0)+g%f2(i0,j0) &
!- ( g%f2(i0,j0) -2*g%f2(i0-1,j0)+g%f2(i0-2,j0) ) )/(2*h**3) !30 derivative
!print*,g%f3der(i0,j0,4)
!print*,g%f1(i0,j0)
!print*,( g%f1(i0,j0+1) -2*g%f1(i0,j0)+g%f1(i0,j0-1) )/(h**2) !02 derivative
!print*,( g%f1(i0+1,j0) -2*g%f1(i0,j0)+g%f1(i0-1,j0) )/(h**2) !20 derivative
!print*
!print*,g%f1der(i0,j0,5)
!print*,g%f1der(i0,j0,3)
!print*,g%f3der(i0,j0,2)
!print*,g%f3der(i0,j0+1,2)
!print*

  h = g%h; 

! set i,j coordinates of 4x4 nbhd of i0,j0 (16 points)
  ind1(1) = i0-2;  ind1(5) = i0-1; ind1(9)  = i0; ind1(13) = i0+1;
  ind1(2) = i0-2;  ind1(6) = i0-1; ind1(10) = i0; ind1(14) = i0+1;
  ind1(3) = i0-2;  ind1(7) = i0-1; ind1(11) = i0; ind1(15) = i0+1;
  ind1(4) = i0-2;  ind1(8) = i0-1; ind1(12) = i0; ind1(16) = i0+1;

  ind2(1) = j0-2;  ind2(5) = j0-2; ind2(9)  = j0-2; ind2(13) = j0-2;
  ind2(2) = j0-1;  ind2(6) = j0-1; ind2(10) = j0-1; ind2(14) = j0-1;
  ind2(3) = j0;    ind2(7) = j0;   ind2(11) = j0;   ind2(15) = j0;
  ind2(4) = j0+1;  ind2(8) = j0+1; ind2(12) = j0+1; ind2(16) = j0+1;

! set alf,bet at those 16 points
  do k=1,16
    alfa_value(k) = -pi   + ind1(k)*h
    beta_value(k) = -pi/2 + ind2(k)*h
  enddo

! INTERPOLATE ALL VALUES OF F1 TO BASEPOINT
! set values of f at those 16 points
  comp_value(:) = 0; f_value(:) = 0
  do k=1,16
    f_value(k) = g%f1(ind1(k),ind2(k))
  enddo
! interpolate f to alfb,betb
  call comp_value_16by16(comp_value(1),n,alfb,betb,alfa_value,beta_value,f_value)
  !values of function at basepoint
  cf1(0,0) = comp_value(1)

! repeat for all 9 derivatives
  do ik = 1,9  ! f1x,f1y,f1xx,f1xy,f1yy,f1xxx,f1xxy,f1xyy,f1yyy
    f_value(:) = 0
    do k=1,16
      f_value(k) =  g%f1der(ind1(k),ind2(k),ik)
    enddo
    call comp_value_16by16(comp_value(ik),n,alfb,betb,alfa_value,beta_value,f_value)
  enddo
  !values of all derivatives at basepoint
  cf1(1,0) = comp_value(1)
  cf1(0,1) = comp_value(2)
  
  cf1(2,0) = comp_value(3)
  cf1(1,1) = comp_value(4)
  cf1(0,2) = comp_value(5)

  cf1(3,0) = comp_value(6)
  cf1(2,1) = comp_value(7)
  cf1(1,2) = comp_value(8)
  cf1(0,3) = comp_value(9)


! INTERPOLATE ALL VALUES OF F2 TO BASEPOINT
  comp_value(:) = 0; f_value(:) = 0
  do k=1,16
    f_value(k) = g%f2(ind1(k),ind2(k))
  enddo
  call comp_value_16by16(comp_value(1),n,alfb,betb,alfa_value,beta_value,f_value)
  cf2(0,0) = comp_value(1)

  do ik = 1,9  ! f1x,f1y,f1xx,f1yy,f1xy,f1xxx,f1xxy,f1xyy,f1yyy
    f_value(:) = 0
    do k=1,16
      f_value(k) = g%f2der(ind1(k),ind2(k),ik)
    enddo
    call comp_value_16by16(comp_value(ik),n,alfb,betb,alfa_value,beta_value,f_value)
  enddo
  !values of all derivatives at basepoint
  cf2(1,0) = comp_value(1)
  cf2(0,1) = comp_value(2)

  cf2(2,0) = comp_value(3)
  cf2(1,1) = comp_value(4)
  cf2(0,2) = comp_value(5)
  
  cf2(3,0) = comp_value(6)
  cf2(2,1) = comp_value(7)
  cf2(1,2) = comp_value(8)
  cf2(0,3) = comp_value(9)


! INTERPOLATE ALL VALUES OF F3 TO BASEPOINT
  comp_value(:) = 0; f_value(:) = 0
  do k=1,16
    f_value(k) = g%f3(ind1(k),ind2(k))
  enddo
  call comp_value_16by16(comp_value(1),n,alfb,betb,alfa_value,beta_value,f_value)
  cf3(0,0) = comp_value(1)

  do ik = 1,9  ! f1x,f1y,f1xx,f1yy,f1xy,f1xxx,f1xxy,f1xyy,f1yyy
    f_value(:) = 0
    do k=1,16
      f_value(k) = g%f3der(ind1(k),ind2(k),ik)
    enddo
    call comp_value_16by16(comp_value(ik),n,alfb,betb,alfa_value,beta_value,f_value)
  enddo
  cf3(1,0) = comp_value(1)
  cf3(0,1) = comp_value(2)

  cf3(2,0) = comp_value(3)
  cf3(1,1) = comp_value(4)
  cf3(0,2) = comp_value(5)

  cf3(3,0) = comp_value(6)
  cf3(2,1) = comp_value(7)
  cf3(1,2) = comp_value(8)
  cf3(0,3) = comp_value(9)

return
END


SUBROUTINE comp_value_16by16(value_comp,n,x0,y0,x,y,value)
  implicit none
  integer :: i,j,k,nk
  ! input
  integer,intent(in) :: n
  real*8,dimension(n),intent(in) :: x,y,value
  real*8,intent(in)   :: x0,y0
  integer :: num
  ! output 
  real*8,intent(out) :: value_comp

  ! integer :: 
  integer,parameter :: NP=40;
  real*8 :: a(NP,NP),b(NP),c(NP),d(NP),xsol(NP)
  logical sing

  if(n.ne.16)then
    write(*,*) 'bicubic n /= 16, STOP'
    stop
  endif

  do nk = 1,n
    num=0
    do j=0,3
    do i=0,3
      num = num+1
      a(nk,num) = x(nk)**real(i,kind=8) * y(nk)**real(j,kind=8)
    enddo
    enddo
  enddo

  ! do qr decomposition	
  call qrdcmp(a,n,NP,c,d,sing)
  if (sing) then
    write(*,*) 'Singularity in QR decomposition.'   
    stop
  endif   
  do i=1,n
    xsol(i) = value(i)
  enddo    
  call qrsolv(a,n,NP,c,d,xsol)      
         
  num=0; value_comp = 0
  do j=0,3
  do i=0,3
    num = num+1
    value_comp = value_comp + xsol(num) * x0**real(i,kind=8) * y0**real(j,kind=8)
  enddo
  enddo
    
return
END SUBROUTINE comp_value_16by16


SUBROUTINE qrdcmp(a,n,np,c,d,sing)
! *********************
! QR = a
! a : input matrix
! np : physical dimension of a
! n : real_used dimension of a
! a : output, upper triangular matrix R is returned to a
! d : diagonal of R
! sing : returns true if singularis encountered during the decomposition
! c : check NR
! *********************
      INTEGER n,np
      REAL*8 a(np,np),c(n),d(n)
      LOGICAL sing
      INTEGER i,j,k
      REAL*8 scale,sigma,sum,tau
      sing=.false.
      do 17 k=1,n-1
        scale=0
        do 11 i=k,n
          scale=max(scale,abs(a(i,k)))
11      continue
        if(scale.eq.0.)then
          sing=.true.
          c(k)=0
          d(k)=0
        else
          do 12 i=k,n
            a(i,k)=a(i,k)/scale
12        continue
          sum=0
          do 13 i=k,n
            sum=sum+a(i,k)**2
13        continue
          sigma=sign(sqrt(sum),a(k,k))
          a(k,k)=a(k,k)+sigma
          c(k)=sigma*a(k,k)
          d(k)=-scale*sigma
          do 16 j=k+1,n
            sum=0
            do 14 i=k,n
              sum=sum+a(i,k)*a(i,j)
14          continue
            tau=sum/c(k)
            do 15 i=k,n
              a(i,j)=a(i,j)-tau*a(i,k)
15          continue
16        continue
        endif
17    continue
      d(n)=a(n,n)
      if(d(n).eq.0.)sing=.true.
      return
      END SUBROUTINE qrdcmp

      SUBROUTINE qrsolv(a,n,np,c,d,b)
! *********************
! a,c,d : output of 'qrdcmp'
! b : right hand side
! np : physical dimension of a
! n : real_used dimension of a
! ******************
      INTEGER n,np
      REAL*8 a(np,np),b(n),c(n),d(n)
!U    USES rsolv
      INTEGER i,j
      REAL*8 sum,tau
      do 13 j=1,n-1
        sum=0.
        do 11 i=j,n
          sum=sum+a(i,j)*b(i)
11      continue
        tau=sum/c(j)
        do 12 i=j,n
          b(i)=b(i)-tau*a(i,j)
12      continue
13    continue
      call rsolv(a,n,np,d,b)
      return
      END SUBROUTINE qrsolv


      SUBROUTINE rsolv(a,n,np,d,b)
      INTEGER n,np
      REAL*8 a(np,np),b(n),d(n)
      INTEGER i,j
      REAL*8 sum
      b(n)=b(n)/d(n)
      do 12 i=n-1,1,-1
        sum=0
        do 11 j=i+1,n
          sum=sum+a(i,j)*b(j)
11      continue
        b(i)=(b(i)-sum)/d(i)
12    continue
      return
      END SUBROUTINE rsolv

END MODULE mod_compcoeffdens
