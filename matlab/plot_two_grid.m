function plot_two_grid


clc

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
a=3; b=2; c=1;
Nthe = 20; Npsi = Nthe/2; 


dthe = 2*pi/Nthe; dpsi=dthe;

the = -pi:dthe:pi; 
psi =   0:dpsi:pi; 

[thegr,psigr]=meshgrid(the,psi);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
figure(1);clf(1); 

xgrid = a.*sin(psigr).*cos(thegr);
ygrid = b.*sin(psigr).*sin(thegr);
zgrid = c.*cos(psigr);

surf(xgrid,ygrid,zgrid); hold on; 

set(gca,'fontsize',20); 
xlabel x; ylabel y
axis equal
grid on


% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
figure(2);clf(2); 

zgrid = c.*sin(psigr).*cos(thegr);
xgrid = a.*sin(psigr).*sin(thegr);
ygrid = b.*cos(psigr);

surf(xgrid,ygrid,zgrid); hold on; 


set(gca,'fontsize',20); 
xlabel x; ylabel y
axis equal
grid on