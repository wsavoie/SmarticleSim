clear all;
% close all;
load('D:\ChronoCode\chronoPkgs\Smarticles\matlabScripts\amoeba\smarticleExpVids\rmv3\movieInfo.mat');

figure(1)
SPACE_UNITS = 'm';
TIME_UNITS = 's';

%find lowest number of frames
% Ns=arrayfun(@(x) [size(x.data{:},1)],movs(:));
% minFrames=min(Ns);
% ma = ma.addAll({movs(i).data{1}(1:minFrames,:)});

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

%define curve params [] for all
spk=[]; smart=[0]; gait=[1]; rob=[5]; v=[];

props={spk smart gait rob v};
inds=1;
for i=1:length(movs)
    
    cond=true;
    for j=1:length(props)
        %if empty accept all values
        if ~isempty(props{j})
            %in case multiple numbers in property
            %if no matches then cond = false
            if(~any(props{j}==movs(i).pars(j)))
                cond = false;
            end
        end
    end
    if(cond)
        ma = ma.addAll(movs(i).data(1));
        usedMovs(inds)=movs(i);
        inds=inds+1;
    end
end

figure(1)
ma.plotTracks
ma.labelPlotTracks
text(0,0+.01,'start')
plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
y=get(gca,'ylim'); x=get(gca,'xlim');
c=max(abs(x)); xlim([-c,c]);
c=max(abs(y)); ylim([-c,c]);
axis equal
figText(gcf,14)

figure(2)
ma = ma.computeMSD;
ma.plotMeanMSD(gca, true);
ma.plotMSD;
pts('fit with drift');
[fo gof]=ma.fitMeanMSD;
figText(gcf,14)
%%
figure(3)
ma = ma.computeVCorr;
ma.plotMeanVCorr
m=ma.vcorr{1};
nansN=sum(isnan(m(:,2)));
M = mean(m(10:end-nansN-2,2));
hold on;
for i= 1:length(ma.vcorr)
    nansN=sum(isnan(ma.vcorr{i}(:,2)));
    M(i)= mean([ma.vcorr{i}(10:end-nansN,2)]);
    
end
M=mean(M);
line([ma.vcorr{1}(10,1) ma.vcorr{1}(end,1)], [M M],'color','r','linewidth',3);
figText(gcf,14);
text(.5,.7,['mean = ',num2str(M,'%2.3f')],'units','normalized','fontsize',16)


% figure(4)
% ma = ma.computeDrift('velocity');
% ma.plotDrift;
% ma.labelPlotTracks;
% %
% figure(5)
% ma = ma.computeMSD;
% ma.plotMeanMSD(gca, true)
% ma.fitMeanMSD
% 
% figure(6)
% ma = ma.computeVCorr;
% ma.plotMeanVCorr
% m=ma.vcorr{1};
% nansN=sum(isnan(m(:,2)));
% M = mean(m(10:end-nansN-2,2));
% hold on;
% for i= 1:length(ma.vcorr)
%     nansN=sum(isnan(ma.vcorr{i}(:,2)));
%     M(i)= mean([ma.vcorr{i}(10:end-nansN,2)]);
%     
% end
% M=mean(M);
% line([ma.vcorr{1}(10,1) ma.vcorr{1}(end,1)], [M M],'color','r','linewidth',3);
% figText(gcf,14);
% text(.5,.7,['mean = ',num2str(M,'%2.3f')],'units','normalized','fontsize',16)

%%
figure(7)
alld=cell2mat(ma.tracks);
scatter3(alld(:,2),alld(:,3),alld(:,1),'.');
xlabel('x (m)'); ylabel('y (m)'); zlabel('t(s)');

%% mean of position at each timestep for all tracks
figure(8)
%get min track size
idx=min(cellfun(@(x) size(x,1), ma.tracks(:)));
idxDat=cellfun(@(x) x(1:idx,:),ma.tracks(:),'uniformoutput',false);

%get x and y data for each timestep 
xdat=cell2mat(cellfun(@(x) x(:,2),idxDat,'uniformoutput',false)');
ydat=cell2mat(cellfun(@(x) x(:,3),idxDat,'uniformoutput',false)');

mxdat=mean(xdat,2);
mydat=mean(ydat,2);

mmxdat=mean(xdat(end,:));
mmydat=mean(ydat(end,:));

mmydat/mmxdat
% plot(mxdat,mydat);
f=fit(mxdat,mydat,'poly1');
hold on;
plot(mxdat,mydat);
plot(f,'-');
legend off;
xlabel('x (m)'); ylabel('y (m)');
title(['mean of position for all runs after ',num2str(ma.tracks{1}(idx,1)),' s'])
 text(.7,.2,['m= ',num2str(f(1),2)],'units','normalized');
figText(gcf,14);