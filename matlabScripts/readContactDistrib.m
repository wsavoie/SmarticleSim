function [force,angs,fcs]=readContactDistrib(fname,binWid)
% fold='A:\SmarticleRun\Amoeba_COG_dead_1_pos_+x\f_0.2_rob_5_v_9\';
% fname=fullfile(fold,'PostProcess','RingContact.txt');
% binWid=5;

simD=importdata(fname);
%time, xContact, yContact, zContact, xForce, yForce, zForce
simD.data=simD.data(simD.data(:,1)>.2,:);
simD.data(:,2)=-simD.data(:,2);
simD.data(:,5)=-simD.data(:,5);

%use negative in x dir since it is a left handed coordinate system
pos=[simD.data(:,1) simD.data(:,2) simD.data(:,3) ]; %only need time, x, y
force=[simD.data(:,1) simD.data(:,5) simD.data(:,6)];

% %get radius of circle
% % [match, nomatch]=regexp(simD.textdata,'[0-9.]+','match');
% % r=str2double(match{1});

% angs=atan2(pos(:,3),pos(:,2))+deg2rad(180);
angs=atan2(pos(:,3),pos(:,2));
angs=mod(angs,2*pi);
binW=deg2rad(binWid); %degrees

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
% polarhistogram(angs,2*pi/binW);
% hold on;
% title('P(\theta)');
% 
% 
% subplot(2,2,[3 4]);
% hold on;
% title('P(\theta)'); xlabel('\theta'); ylabel('P(\theta)')
% histogram(angs,2*pi/binW);
% figText(gcf,15)

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
% % ylim([-180 180]);

%% polarForce histogram
% figure(3)
% disc=discretize(angs,0:binW:2*pi);
disc=discretize(angs,linspace(0,2*pi,round(2*pi/binW)+1));
ud=unique(disc);
lenuni=length(ud);
fcs=zeros(lenuni,2);
% polarForceHist=[];

% sum(norm(force(disc==ud(:),2:3)));
for(i=1:lenuni)
fcs(i,1)=ud(i)*binW;
fcs(i,2)=sum(norm(force(disc==ud(i),2:3)));

% fcs(i)=sum(norm(force(disc==ud(i),2:3)));
% polarForceHist=[polarForceHist; repmat(ud(i)*binW,[round(1000*fcs(i)),1])];
end

% polarhistogram(polarForceHist,2*pi/binW)
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