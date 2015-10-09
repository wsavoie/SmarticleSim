
% area of trap = (a+b)*h/2
dt=.005;
omega = 10;
ss= dt*omega; %step size
ss=2*pi/3/120
% w  = .0117;
% lw =.5;
w= 1;
lw = .5;
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
%     ang = 2*pi:-ss:0;
%     ang2= 2*pi:-ss:0;
    ang = 0:ss:pi;
    ang2= 0:ss:pi;
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



if(find(~ang==ang2)) %if ang1 ~= ang2
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
    hold on;
    plot(Area);
    ylabel('Normalized Area');
    xlabel('position')
    title('Abs(Area)')
    
    figure(2);
    hold on;
    plot(Area2);
    ylabel('Normalized Area');
    title('Area')
    xlabel('position')
else
    figure(1);
    hold on;
    x = linspace(0,1.4,29); % l/w value
    y = ang;    
    z1 = zeros(length(y),length(x));
    z2 = zeros(length(y),length(x));
aa = zeros(length(y),length(x));
w=4;
for ii=1:length(x)
    for jj=1:length(y)
    a=(2*cos(y(jj))*w*x(ii)+2*w)*sin(y(jj))*w*x(ii)*1/2;
    aa(jj,ii)= a;
    end
end

hold on;


[mx,ix]=max(aa,[],1);
[my,iy]=max(aa,[],2);
% surf(x,y*180/pi,abs(aa))
% scatter3(x(:),y(ix)*180/pi,abs(max(aa)),'r','MarkerFaceColor','flat');

surf(x,y*180/pi,abs(aa/max(aa(:))))
scatter3(x(:),y(ix)*180/pi,abs(max(aa)/max(aa(:))),'r','MarkerFaceColor','flat');

    figure(1);
    hold on;
dd =[0.25,0.35,0.4,0.7,1,1.2];
% find(y(ix)==dd)
    title('Interior Area(\alpha,l/w)')
    xlabel('l/w')
    ylabel('\alpha (\circ)');
    zlabel('Normalized Area');
end

angplot = 90;
lwplot = 0.5;
figure(11)
xlabel('\alpha')

ylabel('Interior area')
hold on;
plot(y*180/pi,aa(:,find(x==lwplot)));
title(horzcat('interior area vs angle at l/w=',num2str(lwplot)));

figure(12)
xlabel('l/w')
ylabel('interior area (arb)')
hold on;
title(horzcat('interior area vs l/w at \alpha=',num2str(angplot)));
plot(x,aa(find(int32(y*180/pi)==angplot),:))
% figure(12312)
% plot(x,aa(:,find(y==90*pi/180)))
% figure(11)
% tmp = abs(y-pi/2)
% [idx idx] = min(tmp) %index of closest value
% closest = y(idx) %closest value
% xlabel('\alpha')
% ylabel('interior area')
% plot(y*180/pi,aa(:,find(y==.5)))


% c1 = ones(length(t1),2)*l;
% c2 = [l*cos(t1),l*sin(t1)]+l;
% c3 = [w+l*cos(t2),l*sin(t2)]+l;
% c4 = [c1(:,1)+w,c1(:,2)];
% 
% xCoords=[c1(:,1),c2(:,1),c3(:,1),c4(:,1)];
% yCoords=[c1(:,2),c2(:,2),c3(:,2),c4(:,2)];
% Area = zeros(length(t1),1);
% Area2 = zeros(length(t1),1);


% for i=1:length(t1)
%     shapeArea = 0;
%     for j=1:size(xCoords,2)
%         
%         if j==size(xCoords,2)
%             avgH    = (yCoords(i,1)+yCoords(i,j))/2;
%             W       = xCoords(i,1)- xCoords(i,j);
%         else
%             avgH    = (yCoords(i,j)+yCoords(i,j+1))/2;
%             W       = xCoords(i,j+1)-xCoords(i,j);
%         end     
%         shapeArea = shapeArea + avgH*W;
%     end
%     Area(i)=abs(shapeArea);
%     Area2(i)=(shapeArea);
% end
% 
%     figure(1);
%     plot(ang*180/pi,Area);
%     ylabel('Normalized Area');
%     xlabel('\alpha \circ')
%     title('|Interior Area| vs. \alpha')
%     figure(2);
%        figure(2);
%     hold on;
%     plot(ang*180/pi,Area2);
%     ylabel('Normalized Area');
%     title('Interior Area vs. \alpha')
%     xlabel('\alpha (\circ)')

