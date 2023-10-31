SUBROUTINE allintFbackrecur(u1,u2,v1,v2,c,Hint,rexp,p,q,nb)
!SUBROUTINE allintFrecurs(u1,u2,v1,v2,c,Gintr,Gintl,rexp,p,q,nb)
! F(1,0,0)=\int\int 1/rho du dv, rho=sqrt(1+u^2+2*c*u*v+v^2)
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
real*8 :: rho11,rho12,rho21,rho22,minsq,rtc, &
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
!write(*,'(2i3,2e17.8,e12.2)')0,k,(f22-f12),2*(k-1)*Fintt(0,k-1),Fintt(0,k)

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

do m=0,12
do k=0,5
write(*,'(2i3,2e17.8,e12.2)')m,k,Fintt(m,k),Gintr(m,k)
enddo
enddo
stop

!----------
k=5
call compI005(u1,u2,v1,v2,c,Iint(0,0,5))
!print*,Iint(0,0,5)
call compI095(u1,u2,v1,v2,c,Iint(0,9,5))
call compI185(u1,u2,v1,v2,c,Iint(1,8,5))
!print*,Iint(0,9,5)
!print*,Iint(1,8,5)
do m=2,9
n=9-(m-2)
  right= -c*u2**m*Gintr(n-1,k) -  u2**(m-1)*Gintr(n,k) 
  left = -c*u1**m*Gintl(n-1,k) -  u1**(m-1)*Gintl(n,k) 
  top  =  c*v2**n*Fintt(m-1,k) +  v2**(n-1)*Fintt(m,k) 
  bot  =  c*v1**n*Fintb(m-1,k) +  v1**(n-1)*Fintb(m,k) 
  Iint(m,n-2,k)=( (m-1)*Iint(m-2,n,k) -c*(n-m)*Iint(m-1,n-1,k) + (right-left) + (top-bot) )/(n-1)
enddo

m=0; n=1; r=2*k-m-n-1
!m=0
  !n=1
  right=-c*Gintr(0,k)+minsq*u2*Gintr(1,k)
  left =-c*Gintl(0,k)+minsq*u1*Gintl(1,k)
  top  =( 1+minsq*v2**2 )*Fintt(0,k)
  bot  =( 1+minsq*v1**2 )*Fintb(0,k)
  Iint(0,1,k)=( (right-left) + (top-bot) )/(-minsq*r)
!  I(0,1,k)=(-c*G(0,k)+minsq*u*G(1,k)+ ( 1+minsq*v**2 )*F(0,k))/(2*minsq)

  right = ( 1+minsq*u2**2 )*Gintr(0,k)
  left =  ( 1+minsq*u1**2 )*Gintl(0,k)
  top = -c*Fintt(0,k)+minsq*v2*Fintt(1,k)
  bot = -c*Fintb(0,k)+minsq*v1*Fintb(1,k)
  Iint(1,0,k)=( (top-bot) + (right-left) )/(-minsq*r)
!  I(1,0,k)=(-c*F(0,k)+minsq*v*F(1,k)+ ( 1+minsq*u**2 )*G(0,k))/(2*minsq)

  !n>1
  do n=2,12
m=0; r=2*k-m-n-1
    if (n.ne.9) then
      right = -c*Gintr(n-1,k)+minsq*u2*Gintr(n,k) 
      left  = -c*Gintl(n-1,k)+minsq*u1*Gintl(n,k) 
      top =  ( v2**(n-1)+minsq*v2**(n+1) )*Fintt(0,k) 
      bot =  ( v1**(n-1)+minsq*v1**(n+1) )*Fintb(0,k) 
      Iint(0,n,k)=( -(n-1)*Iint(0,n-2,k) + (right-left) + (top-bot) )/(-minsq*r)

      right= ( u2**(n-1)+minsq*u2**(n+1) )*Gintr(0,k) 
      left = ( u1**(n-1)+minsq*u1**(n+1) )*Gintl(0,k) 
      rightmleft=right-left
  if (mod(n-1,2).eq.0) then !n-1>=1, if even: so n>=3
    rightmleft= (u2+u1)*(fact(n-2)+minsq*fact(n))*Gintr(0,k) &
             + (u1**(n-1)+minsq*u1**(n+1))*(Gintr(0,k)-Gintl(0,k))
  endif
      top= -c*Fintt(n-1,k)+minsq*v2*Fintt(n,k) 
      bot= -c*Fintb(n-1,k)+minsq*v1*Fintb(n,k) 
      Iint(n,0,k)=( -(n-1)*Iint(n-2,0,k) + (top-bot) + (rightmleft) )/(-minsq*r)
!    I(n,0,k)=( -(n-1)*I(n-2,0,k)-c*F(n-1,k)+minsq*v*F(n,k) &
!              + ( u**(n-1)+minsq*u**(n+1) )*G(0,k) )/((n+1)*minsq)
    endif
  enddo
!m>0
do m=1,12
n=1; r=2*k-m-n-1
if (r.ne.0) then
  right= c*u2**m*Gintr(0,k)-minsq*u2**(m+1)*Gintr(1,k) 
  left = c*u1**m*Gintl(0,k)-minsq*u1**(m+1)*Gintl(1,k) 
  top = -(1+minsq*v2**2)*Fintt(m,k) 
  bot = -(1+minsq*v1**2)*Fintb(m,k) 
  Iint(m,1,k)=(-c*m*Iint(m-1,0,k)+(right-left) + (top-bot) )/(r*minsq)
endif
  do n=2,12
r=2*k-m-n-1
if (r.ne.0) then
    right = + c*u2**m*Gintr(n-1,k)-minsq*u2**(m+1)*Gintr(n,k) 
    left  = + c*u1**m*Gintl(n-1,k)-minsq*u1**(m+1)*Gintl(n,k) 
    top= - ( v2**(n-1)+minsq*v2**(n+1) )*Fintt(m,k) 
    bot= - ( v1**(n-1)+minsq*v1**(n+1) )*Fintb(m,k) 
    Iint(m,n,k)=((n-1)*Iint(m,n-2,k)-c*m*Iint(m-1,n-1,k) + (right-left) + (top-bot) )/(r*minsq)
endif
!    I(m,n,k)=((n-1)*I(m,n-2,k)-c*m*I(m-1,n-1,k) &
!     + c*u**m*G(n-1,k)-minsq*u**(m+1)*G(n,k) &
!     - ( v**(n-1)+minsq*v**(n+1) )*F(m,k) )/(-(m+n+1)*minsq)
  enddo
enddo
!----------
do k=4,0,-1
  do m=0,12
  do n=0,12
r=2*k-m-n-1
  if (r.ne.0) then
    right = + u2**(m+1)*Gintr(n,k)
    left  = + u1**(m+1)*Gintl(n,k)
    top =   +v2**(n+1)*Fintt(m,k) 
    bot =   +v1**(n+1)*Fintb(m,k) 
    rightmleft= right-left
    topmbot  = top-bot
    Iint(m,n,k) = ( (2*k+1)*Iint(m,n,k+1) - (rightmleft) - (topmbot) )/r
  endif
!    I(m,n,k) = ((2*k-m-n-3)*I(m,n,k-1) &
!        + u**(m+1)*G(n,k-1)+v**(n+1)*F(m,k-1) )/rmin
  enddo
  enddo
  do m=0,12
  do n=0,12
    r=2*k-m-n-1
    if (r.eq.0) then
      if (m.eq.0) then
        Iint(m,n,k)=intex(m,n,k)
        right= u2**(m+1)*Gintr(n,k)
        left = u1**(m+1)*Gintl(n,k)
        top = v2**n*Fintt(m,k) 
        bot = v1**n*Fintb(m,k) 
        Iint(m,n,k)=((2*k+1)*minsq*Iint(m+2,n,k+1)+c*n*Iint(m+1,n-1,k)-c*(top-bot)+(right-left))/(m+1)   
      else
        right= u2**m*Gintr(n+1,k)
        left = u1**m*Gintl(n+1,k)
        top = v2**(n+1)*Fintt(m,k) 
        bot = v1**(n+1)*Fintb(m,k) 
        Iint(m,n,k)=((2*k+1)*minsq*Iint(m,n+2,k+1)+c*m*Iint(m-1,n+1,k)-c*(right-left)+(top-bot))/(n+1)   
!  Iint(m,n,k)=((2*k+1)*minsq*Iint(m,n+2,k+1)+c*m*Iint(m-1,n+1,k)-c*u**m*Gint(n+1,k)+v**(n+1)*Fint(m,k))/(n+1)   !uses Eq (10)
      endif
    endif
  enddo
  enddo
enddo

print*,Iint(1,0,0)
stop

do j=1,nb
  k=(rexp(j)-1)/2
  Hint(j)=Iint(p(j),q(j),k)
enddo

END SUBROUTINE allintFbackrecur

