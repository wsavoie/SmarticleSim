% figure(1);
% clf;
hold on;
axis([-pi/2,pi/2,-pi/2,pi/2]);
axis square
xlabel('\theta_1');
ylabel('\theta_2');
h=imfreehand(gca);
customGait= getPosition(h);
global_theta_1Pos =customGait(:,1);
global_theta_2Pos =customGait(:,2);
hold on;
title('Custom Drawn Gait');
plot(global_theta_1Pos,global_theta_2Pos,'.k');
xlabel('\theta_1');
ylabel('\theta_2');

figText(gcf,15);
% 

x=global_theta_1Pos;
y=global_theta_2Pos;
figure(2);

hold on;
axis([-pi/2,pi/2,-pi/2,pi/2]);
axis square
h = animatedline('Color','r','LineWidth',3,'MarkerFaceColor','black');
for k=1:length(x)
    addpoints(h,x(k),y(k));
    drawnow
end


%%

test1= [0,0;0,2;1,2;1,4;0,4;1,0];
plot(test1(:,1),test1(:,2),'o-');
d=0;
dt1 = diff(test1);
%%
t = 0:pi/100:2*pi;
y = t;
h = plot(t,y,'YDataSource','y');
for k = 1:0.01:10
   y = exp(sin(t.*k));
   refreshdata(h,'caller')
   drawnow
end
