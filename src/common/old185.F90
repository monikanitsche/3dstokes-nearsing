#include "flags.h"
MODULE mod_old185
! Contains old routines to initialize I185, I095, NOT NEEDED
!
! contains
!    compI185(u1,u2,v1,v2,cc,Iint)
!    intf185dx(a,b,y,c,intfdx)
!    intf185dy2(a,b,x,c,intfdy)
!    compI095(u1,u2,v1,v2,cc,Iint)
!    intf095dx(a,b,y,c,intfdx)
!    intf095dy(a,b,x,c,intfdy)


!PRIVATE    !set the default for module
PUBLIC

CONTAINS

SUBROUTINE compI185(u1,u2,v1,v2,cc,Iint)
! taken from i005/findNgI095
! find the smallest square box S_a bigger than D, the use that int_S_a=0
! compute the integral over 2 small rectangles exactly in long direction, and 
! with Gauss Quad with Ng points in short direction. Use Ng st error < 1.e-13
! confirmed that numbers are accurate using samplei005.m
implicit none
real*8,INTENT(IN) :: u1,u2,v1,v2,cc
real*8,INTENT(OUT) :: Iint
!LOCAL
integer :: j,Ng
real*8, dimension(128) :: x,y,w,f
real*8 :: top,bot,c,d,left,right,sum1,sum2

Ng=4
Ng=40
if (v2.ge.abs(v1)) then
  bot=-v2; top=v1; c=-v2; d=v2
else
  bot=v2; top=abs(v1); c=v1; d=-v1
endif
!integrate over thin horizontal rectangle [u1,u2]x[bot,top]
call setgauss(bot,top,Ng,y,w)
sum1=0
do j=1,Ng
  call intf185dx(u1,u2,y(j),cc,f(j))
  sum1=sum1+f(j)*w(j)
enddo

!integrate over thin vertical rectangle [left,right]x[bot,top]
if (u2.ge.abs(u1)) then
  left=-u2; right=u1
else
  left=u2; right=abs(u1)
endif
sum2=0
call setgauss(left,right,Ng,x,w)
do j=1,Ng
  call intf185dy2(c,d,x(j),cc,f(j))
  sum2=sum2+f(j)*w(j)
enddo

Iint= - (sum1+sum2)
END SUBROUTINE compI185


SUBROUTINE intf185dx(a,b,y,c,intfdx)
! int_a^b f_185 dx
implicit none
real*8, INTENT(IN) :: a,b,y,c
real*8, INTENT(OUT) :: intfdx
!LOCAL
real*8 :: x,rhosq,rho,p,q,num,den,fa,fb

q = -1-(1-c**2)*y**2

x=b
rhosq=1+x**2+y**2+2*c*x*y; rho=sqrt(rhosq)
den=rho**9
num=c*y*(x+c*y)*(128*rhosq**4 - 64*rhosq**3*q +48*rhosq**2*q**2 -40*rhosq*q**3) &
+ 35*(1+c*x*y+y**2)*q**4
fb=num/den

x=a
rhosq=1+x**2+y**2+2*c*x*y; rho=sqrt(rhosq)
den=rho**9
num=c*y*(x+c*y)*(128*rhosq**4 - 64*rhosq**3*q +48*rhosq**2*q**2 -40*rhosq*q**3) &
+ 35*(1+c*x*y+y**2)*q**4
fa=num/den

intfdx=y**8*(fb-fa)/(315*q**5)
END SUBROUTINE intf185dx

SUBROUTINE intf185dy2(a,b,x,c,intfdy)
! int_a^b f_185 dx
implicit none
real*8, INTENT(IN) :: a,b,x,c
real*8, INTENT(OUT) :: intfdy
!LOCAL
real*8 :: x2,x4,x6,x8,c2,c4,c6,c8,y,y2,y4,y6,y8,rhosq,rho,den1,den2,q
real*8 :: coeff0,coeff1,coeff2,coeff3,coeff4,coeff5,coeff6,coeff7,coeff8,coeff9,diff

x2=x**2
x4=x2*x2
x6=x4*x2
x8=x6*x2
c2=c**2
c4=c2**2
c6=c4*c2
c8=c6*c2

rhosq=1+x2+b**2+2*c*x*b; rho=sqrt(rhosq); den2=rho**9
rhosq=1+x2+a**2+2*c*x*a; rho=sqrt(rhosq); den1=rho**9

coeff9=-35-140*(1+c2)*x2+70*(-3-6*c2+c4)*x4-28*(5+15*c2-5*c4+c6)*x6+(-35-140*c2+70*c4-28*c6+5*c8)*x8
coeff8= 9*c*x*(-35-140*(1+c2)*x2+70*(-3-6*c2+c4)*x4-28*(5+15*c2-5*c4+c6)*x6+(-35-140*c2+70*c4-28*c6+5*c8)*x8)
coeff7= -72*c2*x2*(1+x2)*(35+35*(3+c2)*x2-7*(-15-10*c2+c4)*x4+(35+35*c2-7*c4+c6)*x6)
coeff6= 168*c*x*(1+x2)**2*(-5-15*(1+3*c2)*x2-15*(1+6*c2+c4)*x4+(-5-45*c2-15*c4+c6)*x6)
coeff5= -1008*c2*x2*(1+x2)**3*(5+10*(1+c2)*x2+(5+10*c2+c4)*x4)
coeff4= -1008*c*x*(1+x2)**4*(1+2*(1+5*c2)*x2+(1+10*c2+5*c4)*x4)
coeff3= -1344*c2*x2*(1+x2)**5*(3+(3+5*c2)*x2)
coeff2= -576*c*x*(1+x2)**6*(1+(1+7*c2)*x2)
coeff1= -1152*c2*x2*(1+x2)**7
coeff0=-128*c*x*(1+x2)**8

diff=(coeff9*(b**9/den2-a**9/den1) &
    +coeff7*(b**7/den2-a**7/den1) &
    +coeff5*(b**5/den2-a**5/den1) &
    +coeff3*(b**3/den2-a**3/den1) &
    +coeff1*(b   /den2-a   /den1)) &
    +(coeff8*(b**8/den2-a**8/den1) &
    +coeff6*(b**6/den2-a**6/den1) &
    +coeff4*(b**4/den2-a**4/den1) &
    +coeff2*(b**2/den2-a**2/den1) &
    +coeff0*(1   /den2-1   /den1))

q = -1-(1-c2)*x2
intfdy=x*diff/(315*q**5)
END SUBROUTINE intf185dy2


SUBROUTINE compI095(u1,u2,v1,v2,cc,Iint)
! taken from i005/findNgI095
! find the smallest square box S_a bigger than D, the use that int_S_a=0
! compute the integral over 2 small rectangles exactly in long direction, and 
! with Gauss Quad with Ng points in short direction. Use Ng st error < 1.e-13
! confirmed that numbers are accurate using samplei005.m
implicit none
real*8,INTENT(IN) :: u1,u2,v1,v2,cc
real*8,INTENT(OUT) :: Iint
!LOCAL
integer :: j,Ng
real*8, dimension(128) :: x,y,w,f
real*8 :: top,bot,c,d,left,right,sum1,sum2

!Ng=4
Ng=40
if (v2.ge.abs(v1)) then
  bot=-v2; top=v1; c=-v2; d=v2
else
  bot=v2; top=abs(v1); c=v1; d=-v1
endif
!integrate over thin horizontal rectangle [u1,u2]x[bot,top]
call setgauss(bot,top,Ng,y,w)
sum1=0
do j=1,Ng
  call intf095dx(u1,u2,y(j),cc,f(j))
  sum1=sum1+f(j)*w(j)
enddo

!integrate over thin vertical rectangle [left,right]x[bot,top]
if (u2.ge.abs(u1)) then
  left=-u2; right=u1
else
  left=u2; right=abs(u1)
endif
sum2=0
call setgauss(left,right,Ng,x,w)
do j=1,Ng
  call intf095dy(c,d,x(j),cc,f(j))
  sum2=sum2+f(j)*w(j)
enddo

Iint= - (sum1+sum2)
END SUBROUTINE compI095


SUBROUTINE intf095dx(a,b,y,c,intfdx)
! int_a^b f_095 dx
implicit none
real*8, INTENT(IN) :: a,b,y,c
real*8, INTENT(OUT) :: intfdx
!LOCAL
real*8 :: x,rhosq,rho,p,q,num,den,fa,fb

x=b
rhosq=1+x**2+y**2+2*c*x*y
rho=sqrt(rhosq)
p = rhosq-(x+c*y)*rho ! term 1
q = -1-(1-c**2)*y**2
den=5*p**9
num=5*q**4/9 + 5*q**3*p + 120*q**2*p**2/7 + 80*q*p**3/3 + 16*p**4
fb=num/den
x=a
rhosq=1+x**2+y**2+2*c*x*y
rho=sqrt(rhosq)
p = rhosq-(x+c*y)*rho ! term 1
q = -1-(1-c**2)*y**2
num=5*q**4/9 + 5*q**3*p + 120*q**2*p**2/7 + 80*q*p**3/3 + 16*p**4
den=5*p**9
fa=num/den
intfdx=y**9*(fb-fa)
END SUBROUTINE intf095dx


SUBROUTINE intf095dy(a,b,x,c,intfdy)
! int_a^b f_095 dx
implicit none
real*8, INTENT(IN) :: a,b,x,c
real*8, INTENT(OUT) :: intfdy
!LOCAL
real*8 :: x2,x4,c2,c4,c6,c8,y,y2,y4,y6,y8,rhosq,rho,den,num,fa,fb,q

x2=x**2
x4=x2*x2
c2=c**2
c4=c2**2
c6=c4*c2
c8=c6*c2

y=b
y2=y**2
y4=y2*y2
y6=y4*y2
y8=y6*y2
rhosq=1+x2+y2+2*c*x*y
rho=sqrt(rhosq)
den=rho**9
num=128+128*x**18+1152*c*x**17*y+576*y2+1008*y4+840*y6+315*y8 &
+ 192*c*x**15*y*(48+7*(3+5*c2)*y2)+576*x**16*(2+(1+7*c2)*y2) &
+ 1008*c*x**13*y*(32+4*(7+10*c2)*y2+(5+10*c2+c4)*y4) &
+ 144*x**14*(32+4*(8+49*c2)*y2+7*(1+10*c2+5*c4)*y4) &
- 168*x**12*(-64-24*(4+21*c2)*y2-6*(7+60*c2+25*c4)*y4+(-5-45*c2-15*c4+c6)*y6) &
+ 72*c*x**11*y*(896+56*(21+25*c2)*y2+28*(15+25*c2+2*c4)*y4+(35+35*c2-7*c4+c6)*y6) &
+ 9*c*x*y*(128+448*y2+560*y4+280*y6+35*y8)&
- 12*c*x**3*y*(-768-112*(21+5*c2)*y2-840*(3+c2)*y4-210*(5+c2)*y6+35*(-3+c2)*y8) &
+ 9*x2*(128+64*(8+7*c2)*y2+112*(7+10*c2)*y4+280*(2+3*c2)*y6+35*(5+4*c2)*y8) &
-18*x4*(-256-224*(4+7*c2)*y2-56*(21+60*c2+5*c4)*y4-140*(5+15*c2+c4)*y6+35*(-5-8*c2+c4)*y8) &
+126*c*x**5*y*(256+32*(21+10*c2)*y2+8*(75+50*c2+c4)*y4-4*(-50-20*c2+c4)*y6+(15-10*c2+3*c4)*y8) &
+42*x**6*(256+96*(8+21*c2)*y2+120*(7+30*c2+5*c4)*y4-4*(-100-450*c2-60*c4+c6)*y6+3*(25+60*c2-15*c4+2*c6)*y8) &
-36*c*x**7*y*(-1792-560*(7+5*c2)*y2-112*(25+25*c2+c4)*y4-2*(350+210*c2-21*c4+c6)*y6+(-35+35*c2-21*c4+5*c6)*y8) &
-9*x**8*(-1792-2240*(2+7*c2)*y2-560*(7+40*c2+10*c4)*y4+56*(-25-150*c2-30*c4+c6)*y6+(-175-560*c2+210*c4-56*c6+5*c8)*y8) &
-9*x**10*(-1792-448*(8+35*c2)*y2-112*(21+150*c2+50*c4)*y4+56*(-10-75*c2-20*c4+c6)*y6+(-35-140*c2+70*c4-28*c6+5*c8)*y8) &
+c*x**9*y*(80640+6720*(21+20*c2)*y2+1008*(75+100*c2+6*c4)*y4+72*(175+140*c2-21*c4+2*c6)*y6+(315-420*c2+378*c4-180*c6+35*c8)*y8)
fb=num/den

y=a
y2=y**2
y4=y2*y2
y6=y4*y2
y8=y6*y2
rhosq=1+x2+y2+2*c*x*y
rho=sqrt(rhosq)
den=rho**9
num=128+128*x**18+1152*c*x**17*y+576*y2+1008*y4+840*y6+315*y8 &
+ 192*c*x**15*y*(48+7*(3+5*c2)*y2)+576*x**16*(2+(1+7*c2)*y2) &
+ 1008*c*x**13*y*(32+4*(7+10*c2)*y2+(5+10*c2+c4)*y4) &
+ 144*x**14*(32+4*(8+49*c2)*y2+7*(1+10*c2+5*c4)*y4) &
- 168*x**12*(-64-24*(4+21*c2)*y2-6*(7+60*c2+25*c4)*y4+(-5-45*c2-15*c4+c6)*y6) &
+ 72*c*x**11*y*(896+56*(21+25*c2)*y2+28*(15+25*c2+2*c4)*y4+(35+35*c2-7*c4+c6)*y6) &
+ 9*c*x*y*(128+448*y2+560*y4+280*y6+35*y8)&
- 12*c*x**3*y*(-768-112*(21+5*c2)*y2-840*(3+c2)*y4-210*(5+c2)*y6+35*(-3+c2)*y8) &
+ 9*x2*(128+64*(8+7*c2)*y2+112*(7+10*c2)*y4+280*(2+3*c2)*y6+35*(5+4*c2)*y8) &
-18*x4*(-256-224*(4+7*c2)*y2-56*(21+60*c2+5*c4)*y4-140*(5+15*c2+c4)*y6+35*(-5-8*c2+c4)*y8) &
+126*c*x**5*y*(256+32*(21+10*c2)*y2+8*(75+50*c2+c4)*y4-4*(-50-20*c2+c4)*y6+(15-10*c2+3*c4)*y8) &
+42*x**6*(256+96*(8+21*c2)*y2+120*(7+30*c2+5*c4)*y4-4*(-100-450*c2-60*c4+c6)*y6+3*(25+60*c2-15*c4+2*c6)*y8) &
-36*c*x**7*y*(-1792-560*(7+5*c2)*y2-112*(25+25*c2+c4)*y4-2*(350+210*c2-21*c4+c6)*y6+(-35+35*c2-21*c4+5*c6)*y8) &
-9*x**8*(-1792-2240*(2+7*c2)*y2-560*(7+40*c2+10*c4)*y4+56*(-25-150*c2-30*c4+c6)*y6+(-175-560*c2+210*c4-56*c6+5*c8)*y8) &
-9*x**10*(-1792-448*(8+35*c2)*y2-112*(21+150*c2+50*c4)*y4+56*(-10-75*c2-20*c4+c6)*y6+(-35-140*c2+70*c4-28*c6+5*c8)*y8) &
+c*x**9*y*(80640+6720*(21+20*c2)*y2+1008*(75+100*c2+6*c4)*y4+72*(175+140*c2-21*c4+2*c6)*y6+(315-420*c2+378*c4-180*c6+35*c8)*y8)
fa=num/den

q = -1-(1-c2)*x2
intfdy=(fb-fa)/(315*q**5)
END SUBROUTINE intf095dy

END MODULE mod_EHpqr
