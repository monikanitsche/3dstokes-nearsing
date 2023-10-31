MODULE mod_init
! *** ELLIPSE 3-2-1 uinf=(1,0,0) ****
! initializes 
!    uinf
!    ra,rb,rc 
!    corresponding density on grid

CONTAINS

SUBROUTINE init(n0)
use globalvars, ONLY: uinf,ra,rb,rc,axmax
use mod_geom !for setgrid
implicit none
integer :: n0
!LOCAL
integer :: i2,i3
character(len=34) :: flin

!ra=3; rb=1; rc=2; uinf=(/ 1, 0, 0 /); axmax=3; !use with f1=1
ra=3; rb=2; rc=1; uinf=(/ 1, 0, 0 /); axmax=3; !use with f1=1
!ra=1; rb=3; rc=2; uinf=(/ 1, 0, 0 /); axmax=3; !use with f2=1
!ra=2; rb=3; rc=1; uinf=(/ 1, 0, 0 /); axmax=3; !use with f2=1
!ra=1; rb=2; rc=3; uinf=(/ 1, 0, 0 /); axmax=3; !use with f3=1
!ra=2; rb=1; rc=3; uinf=(/ 1, 0, 0 /); axmax=3; !use with f3=1

call setgrid(n0)  ! requires ra,rb,rc

flin = 'data/ellipse-data-321-uinf100-n000'
i2=mod(n0/10,10); i3=mod((n0-i2*10)/100,10)
write(flin(32:33),'(2i1)')i3,i2
!print*,flin
!call setdens(flin)  !flow past sphere
call setdens2  !f=1,0,0

return
END SUBROUTINE init


SUBROUTINE setdens(flin)
use globalvars, ONLY : g
implicit none
character(len=34), INTENT(in) :: flin
!LOCAL
integer :: n,m,nread,ii,i0,j0,i,l
real*8, dimension(12) :: fread
character(len=1) :: numstr

n=g(1)%n
m=g(1)%m
nread = n*(m-1)


do l=1,2
  g(l)%f1=0;g(l)%f2=0;g(l)%f3=0;
  g(l)%f1der=0;g(l)%f2der=0;g(l)%f3der=0;
  write(numstr,'(i1)')l
  !  Line ii=i0+1+(j0-1)*40
  !print*,flin//'/density1_grid'//numstr//'.txt'
  open(12,file=flin//'/density1_grid'//numstr//'.txt')
!   open(unit = 12,status='old',file=flnm,form='formatted', err = 997)
  do ii=1,nread
    read(12,*) fread(1:12)
    i0=mod(ii-1,n); j0=(ii-(i0+1))/n+1
    g(l)%f1(i0,j0) = fread(3)
    g(l)%f1der(i0,j0,1:9) = fread(4:12)  ! f1x,f1y,f1xx,f1xy,f1yy,f1xxx,f1xxy,f1xyy,f1yyy
  enddo
  close(12)

  open(12,file=flin//'/density2_grid'//numstr//'.txt')
  do ii=1,nread
    read(12,*) fread(1:12)
    i0=mod(ii-1,n); j0=(ii-(i0+1))/n+1
    g(l)%f2(i0,j0) = fread(3)
    g(l)%f2der(i0,j0,1:9) = fread(4:12)  ! f2x,f2y,f2xx,f2xy,f2yy,f2xxx,f2xxy,f2xyy,f2yyy
  enddo
  close(12)

  open(12,file=flin//'/density3_grid'//numstr//'.txt')
  do ii=1,nread
    read(12,*) fread(1:12)
    i0=mod(ii-1,n); j0=(ii-(i0+1))/n+1
    g(l)%f3(i0,j0) = fread(3)
    g(l)%f3der(i0,j0,1:9) = fread(4:12)  ! f3x,f3y,f3xx,f3xy,f3yy
  enddo
  close(12)

  do i=-g(l)%noff,-1
    g(l)%f1(i,1:m-1) = g(l)%f1(i+n,1:m-1)
    g(l)%f2(i,1:m-1) = g(l)%f2(i+n,1:m-1)
    g(l)%f3(i,1:m-1) = g(l)%f3(i+n,1:m-1)
    g(l)%f1der(i,1:m-1,1:9) = g(l)%f1der(i+n,1:m-1,1:9) 
    g(l)%f2der(i,1:m-1,1:9) = g(l)%f2der(i+n,1:m-1,1:9) 
    g(l)%f3der(i,1:m-1,1:9) = g(l)%f3der(i+n,1:m-1,1:9) 
  enddo
  do i=0,g(l)%noff
    g(l)%f1(i+n,1:m-1) = g(l)%f1(i,1:m-1)
    g(l)%f2(i+n,1:m-1) = g(l)%f2(i,1:m-1)
    g(l)%f3(i+n,1:m-1) = g(l)%f3(i,1:m-1)
    g(l)%f1der(i+n,1:m-1,1:9) = g(l)%f1der(i,1:m-1,1:9) 
    g(l)%f2der(i+n,1:m-1,1:9) = g(l)%f2der(i,1:m-1,1:9) 
    g(l)%f3der(i+n,1:m-1,1:9) = g(l)%f3der(i,1:m-1,1:9) 
  enddo
enddo
        
!print*,'in init'
!print*,g(1)%f1(320,100)
return
END SUBROUTINE setdens

SUBROUTINE setdens2
use globalvars, ONLY : g
!LOCAL
integer :: l


do l=1,2
  g(l)%f1=1;g(l)%f2=0;g(l)%f3=0;
!  g(l)%f1=0;g(l)%f2=1;g(l)%f3=0;
!  g(l)%f1=0;g(l)%f2=0;g(l)%f3=1;
enddo
        
return
END SUBROUTINE setdens2
   
END MODULE mod_init
