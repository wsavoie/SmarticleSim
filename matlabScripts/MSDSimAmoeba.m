fold=uigetdir('A:\SmarticleRun\')
load(fullfile(fold,'amoebaData.mat'));

% load('A:\SmarticleRun\Amoeba_newsquare_1_dead\amoebaData.mat');
close all
SPACE_UNITS='m';
TIME_UNITS='s';
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

useCOM=0;

f=[.2]; rob=[5]; v=[];
props={f rob v};
for i=1:length(simAm)
    inds=1;
    cond=true;
    for j=1:length(props)
        %if empty accept all values
        if ~isempty(props{j})
            %in case multiple numbers in property
            %if no matches then cond = false
            if(~any(props{j}==simAm(i).pars(j)))
                cond = false;
            end
        end
    end
    if(cond)
        if(useCOM)
            if(isempty(intersect(fieldnames(simAm),'COM')))
                error('trying to use COM when data does not have it');
            end
            ma = ma.addAll(simAm(i).COM);
        else
            ma = ma.addAll(simAm(i).data);
        end
        
%             usedSimAm(inds)=simAm(i);
        inds=inds+1;
    end
end

if(isempty(ma.tracks))
    error('no tracks found with given f,rob,v, specifications');
end
%%
figure(1)
if(useCOM)
title('Ring COM y vs. x');
else
title('Ring COG y vs. x'); 
end
ma.plotTracks
ma.labelPlotTracks
text(0,0+.01,'start')
plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
y=get(gca,'ylim'); x=get(gca,'xlim');
c=max(abs(x)); xlim([-c,c]);
c=max(abs(y)); ylim([-c,c]);

axis equal
figText(gcf,14)
%%
figure(2)
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

figure(3)

ma = ma.computeMSD;
ma.plotMeanMSD(gca, true);
ma.plotMSD;
[a b]=ma.fitMeanMSD;
figText(gcf,14)

%%
figure(4)
ydat=mean(cell2mat(cellfun(@(x) x(:,2),ma.vcorr,'uniformoutput',0)'),2);
% plot(ydat);
y = fft(ydat);
y(1)=[];
n = length(y);
power = abs(y(1:floor(n/2))).^2; % power of first half of transform data
maxfreq = 1/2;                   % maximum frequency
freq = (1:n/2)/(n/2)*maxfreq;    % equally spaced frequency grid
period=[1./freq]';
subplot(1,2,1);
plot(freq,power);
xlabel(' (vcorr/s)')
ylabel('Power')

subplot(1,2,2);
plot(period,power);
xlabel(' (s/vcorr)')
ylabel('Power')

xlim([0,15]);
