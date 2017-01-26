function [force,angs]=readContactDistrib(fold,fname,binW)
fold='D:\SimResults\Chrono\SmarticleU\tests\blah';
fname='blah2';
simAm=struct;

d=fullfile(fold,fname,'PostProcess','RingContact.txt');
simD=importdata(d);
%time, xContact, yContact, zContact, xForce, yForce, zForce
simD.data(:,2)=-simD.data(:,2);
simD.data(:,2)=-simD.data(:,5);

%use negative in x dir since it is a left handed coordinate system
pos=[simD.data(:,1) simD.data(:,2) simD.data(:,3) ]; %only need time, x, y
force=[simD.data(:,1) simD.data(:,5) simD.data(:,6)];
%get radius of circle
[match, nomatch]=regexp(simD.textdata,'[0-9.]+','match');
r=str2double(match{1});
angs=deg2rad(atan2(pos(:,3),pos(:,2)+180));

binwidth=deg2rad(5); %degrees

%% P(theta)
% figure(1);
% subplot(2,2,1);
% hold on;
%  xlabel('x');ylabel('y');
% title('contact positions');
% plot(pos(:,2),pos(:,3),'o')
% axis equal
% subplot(2,2,2);
% % hold on;
% 
% polarhistogram(angs,2*pi/binwidth);
% hold on;
% title('P(\theta)');
% 
% 
% subplot(2,2,[3 4]);
% hold on;
% title('P(\theta)'); xlabel('\theta'); ylabel('P(\theta)')
% histogram(angs,2*pi/binwidth);
% figText(gcf,15)
% xlim([-180 180])
%% Contact position evolution ///// theta vs. time

% figure(2);
% subplot(1,2,1)
% hold on;
% scatter3(pos(:,2),pos(:,3),pos(:,1),ones(numel(pos(:,2)),1));
% xlabel('x'); ylabel('y'); zlabel('t(s)');
% title('Contact position evolution');
% view(40,35)
% figText(gcf,15)
% % axis equal 
% subplot(1,2,2)
% hold on;
% title('\theta vs. time');
%  xlabel('t (s)'); ylabel('\theta');
% figText(gcf,15)
% plot(pos(:,1),angs,'.');
% ylim([-180 180]);

%% polar histogram
figure(3)
disc=discretize(angs,0:binwidth:2*pi);
ud=unique(disc);
lenuni=length(ud);
fcs=zeros(lenuni,1);
c=[];
for(i=1:lenuni)
fcs(i)=sum(norm(force(disc==ud(i),2:3)));
c=[c; repmat(ud(i)*binwidth,[round(1000*fcs(i)),1])];
end
% polarhistogram(c,2*pi/binwidth)
% title('Force (mN) vs. \theta')

%% sort of movie of point creations
% figure(4);
% [uniT,cols]=unique(pos(:,1));
% xmax=max(pos(:,2)); xmin=min(pos(:,2));
% ymax=max(pos(:,3)); ymin=min(pos(:,3));
% for(i=2:length(uniT))
%     col=(cols(i-1):cols(i)-1);
%     
%     
%     plot(pos(1:cols(i-1),2),pos(1:cols(i-1),3),'k.');
%   hold on;
%     plot(pos(col,2),pos(col,3),'-r.','markersize',20);
%       hold off;
%     p=sprintf('time=%.3f',uniT(i));
%     text(.1,.1,p,'units','normalized','fontsize',15);
%     axis([xmin,xmax,ymin,ymax]);
%     axis square
%     pause();
% end
%     