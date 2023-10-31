n=70; m=20; a=3; b=2; c=1;
alf=linspace(-pi,pi,n+1);
bet=linspace(-pi/4,pi/4,m/2+1);

x=a*cos(alf)'*cos(bet);
y=b*sin(alf)'*cos(bet);
z=c*ones(size(alf))'*sin(bet);

n=60; m=20; a=3; b=2; c=1;
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

