
% area of trap = (a+b)*h/2
dt=.0005;
omega = 10;
ss= dt*omega; %step size

% w  = .0117;
% lw =.5;
w= 1;
lw = 1;
l = lw*w;

shape=4;
switch shape
case 1%%%%%circle%%%%%%%
ang = 0:ss:2*pi;
r=pi/2;
phi = 0;
t1=r*cos(ang-phi)';
t2=r*sin(ang-phi)';

case 2%%%%%square%%%%%%%
sL = pi/2; %square gait side length
%theta1
t1 = sL/2:-ss:-sL/2;
l1 = -sL/2*ones(1,length(t1));
b1 = -sL/2:ss:sL/2;
r1 = sL/2*ones(1,length(t1));
%       theta2
t2 = sL/2*ones(1,length(t1));
l2 = sL/2:-ss:-sL/2;
b2 = -sL/2*ones(1,length(t1));
r2 = -sL/2:ss:sL/2;
t1=[t1,l1,b1,r1]';
t2=[t2,l2,b2,r2]';
case 3
    t1=[pi/2,pi/4,pi/2]';
    t2=[pi/2,pi/4,-pi/2]';
    
case 4
    %U shape to straight
    ang = pi/2:-ss:0;
    ang2= pi/2:-ss:0;
    t1=ang';
    t2=ang2';
case 5
    %L shape;
    ang = pi/2:-ss:0;
    ang2= ones(1,length(ang))*pi/2;
    t1=ang';
    t2=ang2';
    
case 6
    %Z Shape
    ang = pi/2:-ss:-pi/2;
    ang2= ones(1,length(ang))*pi/2;
    t1=ang';
    t2=ang2';

%%%%%%%%%%%%%%%%%%
end

c1 = ones(length(t1),2)*l;
c2 = [l*cos(t1),l*sin(t1)]+l;
c3 = [w+l*cos(t2),l*sin(t2)]+l;
c4 = [c1(:,1)+w,c1(:,2)];

xCoords=[c1(:,1),c2(:,1),c3(:,1),c4(:,1)];
yCoords=[c1(:,2),c2(:,2),c3(:,2),c4(:,2)];
Area = zeros(length(t1),1);
Area2 = zeros(length(t1),1);
for i=1:length(t1)
    shapeArea = 0;
    for j=1:size(xCoords,2)
        
        if j==size(xCoords,2)
            avgH    = (yCoords(i,1)+yCoords(i,j))/2;
            W       = xCoords(i,1)- xCoords(i,j);
        else
            avgH    = (yCoords(i,j)+yCoords(i,j+1))/2;
            W       = xCoords(i,j+1)-xCoords(i,j);
        end     
        shapeArea = shapeArea + avgH*W;
    end
    Area(i)=abs(shapeArea);
    Area2(i)=(shapeArea);
end
figure(1);
plot(Area);
ylabel('Area [m^2]');
xlabel('position')
title('Abs(Area)')

figure(2);
hold on;
plot(Area2);
ylabel('Area [m^2]');
title('Area')
xlabel('position')