% function [out]=GetFileInfo(direc)
% direc= 'A:\SmarticleRun\Active75Pos%\';
% direc='A:\SmarticleRun\lazy_tor.3_smarts30\lazy .01\';
clf



% mainFolder = direc;
mainFolder = uigetdir('A:\SmarticleRun\Active98Pos%');
matName = '\PostProcess\stressData.mat';
clear('x','simParams','smartPos','angleInfo','frameInfo','filename','file');
load(horzcat(mainFolder,matName));
ds = 20;
smartPos = downsample(smartPos,ds);
inactiveSmartsIdx = find(~smartPos{1}(:,1));
boxAng = angleInfo(1);
sf = 1; %startFrame

% 	ChVector<>boxdim(.28/1.5, .55245, 2 * bucket_rad / 8);
% boxSize = [.28/1.5, .55245];
 boxSize = [.28/1.5 *2.5, .55245];
allX = zeros(size(smartPos,2)-sf,length(inactiveSmartsIdx));
allY = zeros(size(smartPos,2)-sf,length(inactiveSmartsIdx));

for i=sf+1:size(smartPos,2)
    allX(i-sf,:)=-smartPos{i}(inactiveSmartsIdx,2); %second column is x
    allY(i-sf,:)=smartPos{i}(inactiveSmartsIdx,4)/cos(boxAng);
end
% hold on;
% axis equal
% for(j=1:size(allX,2))
%     plot(allX(:,j),allY(:,j));
% end
% line([-boxSize(1), boxSize(1)],[-boxSize(2),-boxSize(2)],'color','r'); %bott line
% line([-boxSize(1), boxSize(1)],[boxSize(2), boxSize(2)],'color','r'); %top line
% line([-boxSize(1),-boxSize(1)],[-boxSize(2),boxSize(2)],'color','r'); %left
% line([boxSize(1),  boxSize(1)],[-boxSize(2),boxSize(2)],'color','r'); %right

%%ajsdsa

% idx = 3

%%s
% dt = linspace(0,length(allX))';
clear tracks;
dt=[(sf-1):length(allX)-2+sf]'*.00025*ds;

%%%%%%%%%%%all tracks%%%%%%%%%%%%%%%%
% tracks=cell(1,size(allX,2));
% for j=1:size(allX,2)
%     tracks{j}= [dt,allX(:,j),allY(:,j)];
% end
%%%%%%%%%%special tracks%%%%%%%%%%%
% tracks=cell(1,1);
% idx = 70:79;
% for i=1:length(idx)
% tracks{i}= [dt,allX(:,idx(i)),allY(:,idx(i))];
% end
%%%%%%%%%%single tracks%%%%%%%%%%

tracks=cell(1);

idx = 3;
tt=find(dt>10.067 & dt<11.5);
tt=1:length(dt);
dt=dt(tt);
tracks{1}= [dt,allX(tt,idx),allY(tt,idx)];
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comet(allX(tt,idx),allY(tt,idx),.2)
%%
clf;
figure(1)
ma = msdanalyzer(2,'m','s');
ma = ma.addAll(tracks);
ma.plotTracks
ma.labelPlotTracks
axis equal;
hold on;
line([-boxSize(1), boxSize(1)],[-boxSize(2),-boxSize(2)],'color','r'); %bott line
line([-boxSize(1), boxSize(1)],[boxSize(2), boxSize(2)],'color','r'); %top line
line([-boxSize(1),-boxSize(1)],[-boxSize(2),boxSize(2)],'color','r'); %left
line([boxSize(1),  boxSize(1)],[-boxSize(2),boxSize(2)],'color','r'); %right

%%
figure(2);
ma = ma.computeMSD;
ma.plotMeanMSD(gca, true);

%%
figure(3)
ma = ma.computeVCorr;
ma.plotMeanVCorr;

% meanIdxs = find(tracks{1}(:,1)>10 & tracks{1}(:,1)<10.5);
% M = mean(ma.vcorr{1}(meanIdxs,2));
% line([0 15], [M M],'color','r','linewidth',3);