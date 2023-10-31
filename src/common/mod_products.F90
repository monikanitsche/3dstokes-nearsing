! #include "flags.h"
MODULE products
!    prod222
!    prod244
!    prod333
!    prod344
!    prod445
!    prod456
!    prod466

!PRIVATE    !set the default for module
!PUBLIC :: correctSLPatP
PUBLIC

CONTAINS
SUBROUTINE prod222(x,y,z)
!     computes coefficients of z=x*y where x,y,z are degree 2
use params, ONLY : cmax
implicit none
real*8, dimension(0:cmax,0:cmax) :: x,y,z

      z(0,0)= x(0,0)*y(0,0)
      z(1,0)= x(0,0)*y(1,0)+x(1,0)*y(0,0)
      z(0,1)= x(0,0)*y(0,1)+x(0,1)*y(0,0)
      z(2,0)= x(0,0)*y(2,0)+x(1,0)*y(1,0)+x(2,0)*y(0,0)
      z(0,2)= x(0,0)*y(0,2)+x(0,1)*y(0,1)+x(0,2)*y(0,0)
      z(1,1)= x(1,0)*y(0,1)+x(0,0)*y(1,1)+x(1,1)*y(0,0)+x(0,1)*y(1,0)
      return 
END SUBROUTINE prod222


SUBROUTINE prod233(x,y,z)
!     computes coefficients of z=x*y where x has degree 2,y,z are degree 3
use params, ONLY : cmax
implicit none
real*8, dimension(0:cmax,0:cmax) :: x,y,z

      z(0,0)= x(0,0)*y(0,0)
      z(1,0)= x(0,0)*y(1,0)+x(1,0)*y(0,0)
      z(0,1)= x(0,0)*y(0,1)+x(0,1)*y(0,0)
      z(2,0)= x(0,0)*y(2,0)+x(1,0)*y(1,0)+x(2,0)*y(0,0)
      z(1,1)= x(1,0)*y(0,1)+x(0,0)*y(1,1)+x(1,1)*y(0,0)+x(0,1)*y(1,0)
      z(0,2)= x(0,0)*y(0,2)+x(0,1)*y(0,1)+x(0,2)*y(0,0)
      z(3,0)= x(0,0)*y(3,0)+x(2,0)*y(1,0)+x(1,0)*y(2,0)
      z(2,1)= x(0,0)*y(2,1)+x(1,0)*y(1,1)+x(2,0)*y(0,1)+x(0,1)*y(2,0)+x(1,1)*y(1,0)
      z(1,2)= x(0,0)*y(1,2)+x(0,1)*y(1,1)+x(0,2)*y(1,0)+x(1,0)*y(0,2)+x(1,1)*y(0,1)
      z(0,3)= x(0,0)*y(0,3)+x(0,2)*y(0,1)+x(0,1)*y(0,2)
      return 
END SUBROUTINE prod233
      
SUBROUTINE prod244(x,y,z)
!     computes coefficients of z=x*y where x has degree 2,y,z are degree 4
use params, ONLY : cmax
implicit none
real*8, dimension(0:cmax,0:cmax) :: x,y,z

      z(0,0)= x(0,0)*y(0,0)
      z(1,0)= x(0,0)*y(1,0)+x(1,0)*y(0,0)
      z(0,1)= x(0,0)*y(0,1)+x(0,1)*y(0,0)
      z(2,0)= x(0,0)*y(2,0)+x(1,0)*y(1,0)+x(2,0)*y(0,0)
      z(1,1)= x(1,0)*y(0,1)+x(0,0)*y(1,1)+x(1,1)*y(0,0)+x(0,1)*y(1,0)
      z(0,2)= x(0,0)*y(0,2)+x(0,1)*y(0,1)+x(0,2)*y(0,0)
      z(3,0)= x(0,0)*y(3,0)+x(2,0)*y(1,0)+x(1,0)*y(2,0)
      z(2,1)= x(0,0)*y(2,1)+x(1,0)*y(1,1)+x(2,0)*y(0,1)+x(0,1)*y(2,0)+x(1,1)*y(1,0)
      z(1,2)= x(0,0)*y(1,2)+x(0,1)*y(1,1)+x(0,2)*y(1,0)+x(1,0)*y(0,2)+x(1,1)*y(0,1)
      z(0,3)= x(0,0)*y(0,3)+x(0,2)*y(0,1)+x(0,1)*y(0,2)
      z(4,0)= x(0,0)*y(4,0)+x(1,0)*y(3,0)+x(2,0)*y(2,0)
      z(3,1)= x(0,0)*y(3,1)+x(1,0)*y(2,1)+x(2,0)*y(1,1)  &
             +x(0,1)*y(3,0)+x(1,1)*y(2,0)
      z(2,2)= x(0,0)*y(2,2)+x(1,0)*y(1,2)+x(2,0)*y(0,2)  &
             +x(0,1)*y(2,1)+x(1,1)*y(1,1)  &
             +x(0,2)*y(2,0)
      z(1,3)= x(0,0)*y(1,3)+x(0,1)*y(1,2)+x(0,2)*y(1,1)  &
             +x(1,0)*y(0,3)+x(1,1)*y(0,2)
      z(0,4)= x(0,0)*y(0,4)+x(0,1)*y(0,3)+x(0,2)*y(0,2)
      return 
END SUBROUTINE prod244
      

SUBROUTINE prod333(x,y,z)
!     computes coefficients of z=x*y where x,y are degree 3, z is degree 4
use params, ONLY : cmax
implicit none
real*8, dimension(0:cmax,0:cmax) :: x,y,z
integer :: j,k

      do j=0,3
      do k=0,3
      z(j,k)=0
      enddo
      enddo

      z(0,0)=x(0,0)*y(0,0)
      z(0,1)=x(0,0)*y(0,1)+x(0,1)*y(0,0)
      z(0,2)=x(0,0)*y(0,2)+x(0,1)*y(0,1)+x(0,2)*y(0,0)
      z(0,3)=x(0,0)*y(0,3)+x(0,1)*y(0,2)+x(0,2)*y(0,1)+x(0,3)*y(0,0)
      do k=0,1
        j=1-k
        z(1,0)=z(1,0)+x(k,0)*y(j,0)
        z(1,1)=z(1,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(1,2)=z(1,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      enddo
      do k=0,2
        j=2-k
        z(2,0)=z(2,0)+x(k,0)*y(j,0)
        z(2,1)=z(2,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      enddo
      do k=0,3
        j=3-k
        z(3,0)=z(3,0)+x(k,0)*y(j,0)
      enddo
      return 
END SUBROUTINE prod333


SUBROUTINE prod334(x,y,z)
!     computes coefficients of z=x*y where x,y are degree 3, z is degree 4
use params, ONLY : cmax
implicit none
real*8, dimension(0:cmax,0:cmax) :: x,y,z
integer :: j,k

      do j=0,4
      do k=0,4
      z(j,k)=0
      enddo
      enddo

      z(0,0)=x(0,0)*y(0,0)
      z(0,1)=x(0,0)*y(0,1)+x(0,1)*y(0,0)
      z(0,2)=x(0,0)*y(0,2)+x(0,1)*y(0,1)+x(0,2)*y(0,0)
      z(0,3)=x(0,0)*y(0,3)+x(0,1)*y(0,2)+x(0,2)*y(0,1)+x(0,3)*y(0,0)
      z(0,4)=              x(0,1)*y(0,3)+x(0,2)*y(0,2)+x(0,3)*y(0,1)
      do k=0,1
        j=1-k
        z(1,0)=z(1,0)+x(k,0)*y(j,0)
        z(1,1)=z(1,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(1,2)=z(1,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      enddo
      k=0
      j=1
        z(1,3)=z(1,3)              +x(k,1)*y(j,2)+x(k,2)*y(j,1)+x(k,3)*y(j,0)
      k=1
      j=0
        z(1,3)=z(1,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)

      do k=0,2
        j=2-k
        z(2,0)=z(2,0)+x(k,0)*y(j,0)
        z(2,1)=z(2,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      enddo
      k=0
      j=2
        z(2,2)=z(2,2)              +x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=1
      j=1
        z(2,2)=z(2,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=2
      j=0
        z(2,2)=z(2,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)

      do k=0,3
        j=3-k
        z(3,0)=z(3,0)+x(k,0)*y(j,0)
      enddo
      k=0
      j=3
        z(3,1)=z(3,1)              +x(k,1)*y(j,0)
      k=1
      j=2
        z(3,1)=z(3,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      k=2
      j=1
        z(3,1)=z(3,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      k=3
      j=0
        z(3,1)=z(3,1)+x(k,0)*y(j,1)

      do k=1,3
        j=4-k
        z(4,0)=z(4,0)+x(k,0)*y(j,0)
      enddo
      return 
END SUBROUTINE prod334
      
      

SUBROUTINE prod344(x,y,z)
!     computes coefficients of z=x*y where x,y are degree 3, z is degree 4
use params, ONLY : cmax
implicit none
real*8, dimension(0:cmax,0:cmax) :: x,y,z
integer :: j,k

      do j=0,4
      do k=0,4
      z(j,k)=0
      enddo
      enddo

      z(0,0)=x(0,0)*y(0,0)
      z(0,1)=x(0,0)*y(0,1)+x(0,1)*y(0,0)
      z(0,2)=x(0,0)*y(0,2)+x(0,1)*y(0,1)+x(0,2)*y(0,0)
      z(0,3)=x(0,0)*y(0,3)+x(0,1)*y(0,2)+x(0,2)*y(0,1)+x(0,3)*y(0,0)
      z(0,4)=x(0,0)*y(0,4)+x(0,1)*y(0,3)+x(0,2)*y(0,2)+x(0,3)*y(0,1)
      do k=0,1
        j=1-k
        z(1,0)=z(1,0)+x(k,0)*y(j,0)
        z(1,1)=z(1,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(1,2)=z(1,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      enddo
      k=0
      j=1
        z(1,3)=z(1,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)+x(k,3)*y(j,0)
      k=1
      j=0
        z(1,3)=z(1,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)

      do k=0,2
        j=2-k
        z(2,0)=z(2,0)+x(k,0)*y(j,0)
        z(2,1)=z(2,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      enddo
      k=0
      j=2
        z(2,2)=z(2,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=1
      j=1
        z(2,2)=z(2,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=2
      j=0
        z(2,2)=z(2,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)

      do k=0,3
        j=3-k
        z(3,0)=z(3,0)+x(k,0)*y(j,0)
      enddo
      k=0
      j=3
        z(3,1)=z(3,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      k=1
      j=2
        z(3,1)=z(3,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      k=2
      j=1
        z(3,1)=z(3,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      k=3
      j=0
        z(3,1)=z(3,1)+x(k,0)*y(j,1)

      do k=0,3
        j=4-k
        z(4,0)=z(4,0)+x(k,0)*y(j,0)
      enddo
      return 
END SUBROUTINE prod344
      

SUBROUTINE prod445(x,y,z)
!     computes coefficients of z=x*y where x,y are degree 3, z is degree 4
use params, ONLY : cmax
implicit none
real*8, dimension(0:cmax,0:cmax) :: x,y,z
integer :: j,k

      do j=0,5
      do k=0,5
      z(j,k)=0
      enddo
      enddo

      z(0,0)=x(0,0)*y(0,0)
      z(0,1)=x(0,0)*y(0,1)+x(0,1)*y(0,0)
      z(0,2)=x(0,0)*y(0,2)+x(0,1)*y(0,1)+x(0,2)*y(0,0)
      z(0,3)=x(0,0)*y(0,3)+x(0,1)*y(0,2)+x(0,2)*y(0,1)+x(0,3)*y(0,0)
      z(0,4)=x(0,0)*y(0,4)+x(0,1)*y(0,3)+x(0,2)*y(0,2)+x(0,3)*y(0,1)+x(0,4)*y(0,0)
      z(0,5)=              x(0,1)*y(0,4)+x(0,2)*y(0,3)+x(0,3)*y(0,2)+x(0,4)*y(0,1)
      do k=0,1
        j=1-k
        z(1,0)=z(1,0)+x(k,0)*y(j,0)
        z(1,1)=z(1,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(1,2)=z(1,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
        z(1,3)=z(1,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)+x(k,3)*y(j,0)
      enddo
      k=0
      j=1
        z(1,4)=z(1,4) +x(k,1)*y(j,3)+x(k,2)*y(j,2) +x(k,3)*y(j,1)+x(k,4)*y(j,0)
      k=1
      j=0
        z(1,4)=z(1,4)+x(k,0)*y(j,4)+x(k,1)*y(j,3)+x(k,2)*y(j,2) +x(k,3)*y(j,1)

      do k=0,2
        j=2-k
        z(2,0)=z(2,0)+x(k,0)*y(j,0)
        z(2,1)=z(2,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(2,2)=z(2,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      enddo
      k=0
      j=2
        z(2,3)=z(2,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)+x(k,3)*y(j,0)
      k=1
      j=1
        z(2,3)=z(2,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)+x(k,3)*y(j,0)
      k=2
      j=0
        z(2,3)=z(2,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)

      do k=0,3
        j=3-k
        z(3,0)=z(3,0)+x(k,0)*y(j,0)
        z(3,1)=z(3,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      enddo
      k=0
      j=3
        z(3,2)=z(3,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=1
      j=2
        z(3,2)=z(3,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=2
      j=1
        z(3,2)=z(3,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=3
      j=0
        z(3,2)=z(3,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)

      do k=0,4
        j=4-k
        z(4,0)=z(4,0)+x(k,0)*y(j,0)
      enddo
      k=0
      j=4
        z(4,1)=z(4,1)+x(k,1)*y(j,0)
      k=1
      j=3
        z(4,1)=z(4,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      k=2
      j=2
        z(4,1)=z(4,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      k=3
      j=1
        z(4,1)=z(4,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      k=4
      j=0
        z(4,1)=z(4,1)+x(k,0)*y(j,1)

      do k=1,4
        j=5-k
        z(5,0)=z(5,0)+x(k,0)*y(j,0)
      enddo
      return 
END SUBROUTINE prod445
      

SUBROUTINE prod456(x,y,z)
!     computes coefficients of z=x*y where x deg4,y degree 5, z is degree 6
use params, ONLY : cmax
implicit none
real*8, dimension(0:cmax,0:cmax) :: x,y,z
integer :: j,k

      do j=0,6
      do k=0,6
      z(j,k)=0
      enddo
      enddo

      z(0,0)=x(0,0)*y(0,0)
      z(0,1)=x(0,0)*y(0,1)+x(0,1)*y(0,0)
      z(0,2)=x(0,0)*y(0,2)+x(0,1)*y(0,1)+x(0,2)*y(0,0)
      z(0,3)=x(0,0)*y(0,3)+x(0,1)*y(0,2)+x(0,2)*y(0,1)+x(0,3)*y(0,0)
      z(0,4)=x(0,0)*y(0,4)+x(0,1)*y(0,3)+x(0,2)*y(0,2)+x(0,3)*y(0,1)+x(0,4)*y(0,0)
      z(0,5)=x(0,0)*y(0,5)+x(0,1)*y(0,4)+x(0,2)*y(0,3)+x(0,3)*y(0,2)+x(0,4)*y(0,1)
      z(0,6)=              x(0,1)*y(0,5)+x(0,2)*y(0,4)+x(0,3)*y(0,3)+x(0,4)*y(0,2) 
      do k=0,1
        j=1-k
        z(1,0)=z(1,0)+x(k,0)*y(j,0)
        z(1,1)=z(1,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(1,2)=z(1,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
        z(1,3)=z(1,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)+x(k,3)*y(j,0)
      enddo
      k=0
      j=1
        z(1,4)=z(1,4)+x(k,0)*y(j,4)+x(k,1)*y(j,3)+x(k,2)*y(j,2)  &
               +x(k,3)*y(j,1)+x(k,4)*y(j,0)
        z(1,5)=z(1,5)+x(k,1)*y(j,4)+x(k,2)*y(j,3)  &
               +x(k,3)*y(j,2)+x(k,4)*y(j,1) 
      k=1
      j=0
        z(1,4)=z(1,4)+x(k,0)*y(j,4)+x(k,1)*y(j,3)+x(k,2)*y(j,2) &
              +x(k,3)*y(j,1)
        z(1,5)=z(1,5)+x(k,0)*y(j,5)+x(k,1)*y(j,4)+x(k,2)*y(j,3) &
              +x(k,3)*y(j,2)
      
      do k=0,2
        j=2-k
        z(2,0)=z(2,0)+x(k,0)*y(j,0)
        z(2,1)=z(2,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(2,2)=z(2,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      enddo
      k=0
      j=2
        z(2,3)=z(2,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)+x(k,3)*y(j,0)
        z(2,4)=z(2,4)+x(k,1)*y(j,3)+x(k,2)*y(j,2)+x(k,3)*y(j,1)+x(k,4)*y(j,0)
      k=1
      j=1
        z(2,3)=z(2,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)+x(k,3)*y(j,0)
        z(2,4)=z(2,4)+x(k,0)*y(j,4)+x(k,1)*y(j,3)+x(k,2)*y(j,2)+x(k,3)*y(j,1)
      k=2
      j=0
        z(2,3)=z(2,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)
        z(2,4)=z(2,4)+x(k,0)*y(j,4)+x(k,1)*y(j,3)+x(k,2)*y(j,2)

      do k=0,3
        j=3-k
        z(3,0)=z(3,0)+x(k,0)*y(j,0)
        z(3,1)=z(3,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      enddo
      k=0
      j=3
        z(3,2)=z(3,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1) +x(k,2)*y(j,0)
        z(3,3)=z(3,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1) +x(k,3)*y(j,0)
      k=1
      j=2
        z(3,2)=z(3,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
        z(3,3)=z(3,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1) +x(k,3)*y(j,0)
      k=2
      j=1
        z(3,2)=z(3,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
        z(3,3)=z(3,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)
      k=3
      j=0
        z(3,2)=z(3,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)
        z(3,3)=z(3,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)

      do k=0,4
        j=4-k
        z(4,0)=z(4,0)+x(k,0)*y(j,0)
      enddo
      k=0
      j=4
        z(4,1)=z(4,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(4,2)=z(4,2)              +x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=1
      j=3
        z(4,1)=z(4,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(4,2)=z(4,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=2
      j=2
        z(4,1)=z(4,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(4,2)=z(4,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      k=3
      j=1
        z(4,1)=z(4,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(4,2)=z(4,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)
      k=4
      j=0
        z(4,1)=z(4,1)+x(k,0)*y(j,1)
        z(4,2)=z(4,2)+x(k,0)*y(j,2)

       k=0
       j=5
        z(5,0)=z(5,0)+x(k,0)*y(j,0)
        z(5,1)=z(5,1)              +x(k,1)*y(j,0)
       k=1
       j=4
        z(5,0)=z(5,0)+x(k,0)*y(j,0)
        z(5,1)=z(5,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
       k=2
       j=3
        z(5,0)=z(5,0)+x(k,0)*y(j,0)
        z(5,1)=z(5,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
       k=3
       j=2
        z(5,0)=z(5,0)+x(k,0)*y(j,0)
        z(5,1)=z(5,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
       k=4
       j=1
        z(5,0)=z(5,0)+x(k,0)*y(j,0)
        z(5,1)=z(5,1)+x(k,0)*y(j,1)

      do k=1,4
        j=6-k
        z(6,0)=z(6,0)+x(k,0)*y(j,0)
      enddo
      return 
END SUBROUTINE prod456
      

SUBROUTINE prod666(x,y,z)
!     (not used)
!     computes coefficients of z=x*y where x,y are degree 3, z is degree 4
use params, ONLY : cmax
implicit none
real*8, dimension(0:cmax,0:cmax) :: x,y,z
integer :: j,k

      do j=0,6
      do k=0,6
      z(j,k)=0
      enddo

      enddo
      z(0,0)=x(0,0)*y(0,0)
      z(0,1)=x(0,0)*y(0,1)+x(0,1)*y(0,0)
      z(0,2)=x(0,0)*y(0,2)+x(0,1)*y(0,1)+x(0,2)*y(0,0)
      z(0,3)=x(0,0)*y(0,3)+x(0,1)*y(0,2)+x(0,2)*y(0,1)+x(0,3)*y(0,0)
      z(0,4)=x(0,0)*y(0,4)+x(0,1)*y(0,3)+x(0,2)*y(0,2)+x(0,3)*y(0,1) &
               +x(0,4)*y(0,0)
      z(0,5)=x(0,0)*y(0,5)+x(0,1)*y(0,4)+x(0,2)*y(0,3)+x(0,3)*y(0,2) &
               +x(0,4)*y(0,1)+x(0,5)*y(0,0)
      z(0,6)=x(0,0)*y(0,6)+x(0,1)*y(0,5)+x(0,2)*y(0,4)+x(0,3)*y(0,3) &
               +x(0,4)*y(0,2)+x(0,5)*y(0,1)+x(0,6)*y(0,0)
      do k=0,1
        j=1-k
        z(1,0)=z(1,0)+x(k,0)*y(j,0)
        z(1,1)=z(1,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(1,2)=z(1,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
        z(1,3)=z(1,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)  &
               +x(k,3)*y(j,0)
        z(1,4)=z(1,4)+x(k,0)*y(j,4)+x(k,1)*y(j,3)+x(k,2)*y(j,2)  &
               +x(k,3)*y(j,1)+x(k,4)*y(j,0)
        z(1,5)=z(1,5)+x(k,0)*y(j,5)+x(k,1)*y(j,4)+x(k,2)*y(j,3)  &
               +x(k,3)*y(j,2)+x(k,4)*y(j,1)+x(k,5)*y(j,0)
      enddo
      do k=0,2
        j=2-k
        z(2,0)=z(2,0)+x(k,0)*y(j,0)
        z(2,1)=z(2,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(2,2)=z(2,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
        z(2,3)=z(2,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)  &
               +x(k,3)*y(j,0)
        z(2,4)=z(2,4)+x(k,0)*y(j,4)+x(k,1)*y(j,3)+x(k,2)*y(j,2)  &
               +x(k,3)*y(j,1)+x(k,4)*y(j,0)
      enddo
      do k=0,3
        j=3-k
        z(3,0)=z(3,0)+x(k,0)*y(j,0)
        z(3,1)=z(3,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(3,2)=z(3,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
        z(3,3)=z(3,3)+x(k,0)*y(j,3)+x(k,1)*y(j,2)+x(k,2)*y(j,1)+x(k,3)*y(j,0)
      enddo
      do k=0,4
        j=4-k
        z(4,0)=z(4,0)+x(k,0)*y(j,0)
        z(4,1)=z(4,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
        z(4,2)=z(4,2)+x(k,0)*y(j,2)+x(k,1)*y(j,1)+x(k,2)*y(j,0)
      enddo
      do k=0,5
        j=5-k
        z(5,0)=z(4,0)+x(k,0)*y(j,0)
        z(5,1)=z(4,1)+x(k,0)*y(j,1)+x(k,1)*y(j,0)
      enddo
      do k=0,6
        j=6-k
        z(6,0)=z(6,0)+x(k,0)*y(j,0)
      enddo
      return 
END SUBROUTINE prod666
      
END MODULE products
