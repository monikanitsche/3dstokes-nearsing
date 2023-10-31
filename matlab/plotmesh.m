%n=70; m=20; a=3; b=2; c=1;
n=160; m=40; a=3; b=2; c=1;
%n=280; m=80; a=3; b=2; c=1;
%n=560; m=160; a=3; b=2; c=1;
%n=1120; m=320; a=3; b=2; c=1;
%n=40; m=20; a=3; b=2; c=1;
alf=linspace(-pi,pi,n+1);
bet=linspace(-pi/4,pi/4,m/2+1);

x=a*cos(alf)'*cos(bet);
y=b*sin(alf)'*cos(bet);
z=c*ones(size(alf))'*sin(bet);

%n=60; m=20; a=3; b=2; c=1;
n=120; m=40; a=3; b=2; c=1;
%n=240; m=80; a=3; b=2; c=1;
%n=480; m=160; a=3; b=2; c=1;
%n=960; m=320; a=3; b=2; c=1;
%n=40; m=20; a=3; b=2; c=1;
alf=linspace(-pi,pi,n+1);
bet=linspace(-pi/4,pi/4,m/2+1);
%bet2=linspace(pi/4,pi/4,m+1);
alf=linspace(pi/4,3*pi/4,n/4+1);
bet=linspace(-pi/4,pi/4,m+1);
y2=b*cos(alf)'*cos(bet);
z2=c*sin(alf)'*cos(bet);
x2=a*ones(size(alf))'*sin(bet);

figure(1),clf
mesh(x,y,z)
axis equal
hold on
mesh(x2,y2,z2)
mesh(x2,y2,-z2)

xlabel('x')
ylabel('y')
zlabel('z')
set(gca,'fontsize',15)
