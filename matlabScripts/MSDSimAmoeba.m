fold=uigetdir('A:\SmarticleRun\')
load(fullfile(fold,'amoebaData.mat'));

% load('A:\SmarticleRun\Amoeba_newsquare_1_dead\amoebaData.mat');
close all
clear usedSimAm
SPACE_UNITS='m';
TIME_UNITS='s';
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
inds=1;
%************************************************************
%* Fig numbers:
%* 1. displacement yvsx
%* 2. MSD
%* 3. vcorr
%* 4. fourier transform of vcorr
%* 5. contact P(\theta) POLAR
%* 6. contact P(\theta) LINEAR
%* 7. inner force vs. theta POLAR
%* 8. inner force vs. theta LINEAR
%************************************************************
showFigs=[1 5 7];
useCOM=0;
f=[.2]; rob=[4:5]; v=[20];


props={f rob v};
for i=1:length(simAm)
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
        usedSimAm(inds)=simAm(i);
        inds=inds+1;
        if(useCOM)
            if(isempty(intersect(fieldnames(simAm),'COM')))
                error('trying to use COM when data does not have it');
            end
            ma = ma.addAll(simAm(i).COM);
        else
            ma = ma.addAll(simAm(i).data);
        end
        
        
    end
end

if(isempty(ma.tracks))
    error('no tracks found with given f,rob,v, specifications');
end
%% plot displacement yvsx
if(showFigs(showFigs==1))
    figure(1)
    if(useCOM)
        title('Ring COM');
    else
        title('Ring COG');
    end
    ma.plotTracks
    ma.labelPlotTracks
    text(0,0+.01,'start')
    plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
    y=get(gca,'ylim'); x=get(gca,'xlim');
    c=max(abs(x)); xlim([-c,c]);
    c=max(abs(y)); ylim([-c,c]);
    axis equal
    axis([-.25 .25 -.25 .25]);
    x=xlim; y=ylim;
    plot(x,[0,0],'r');
    plot([0,0],y,'r');
    
    figText(gcf,14)
end

%% plot MSD
if(showFigs(showFigs==2))
    figure(2)
    
    ma = ma.computeMSD;
    ma.plotMeanMSD(gca, true);
    ma.plotMSD;
    [a b]=ma.fitMeanMSD;
    figText(gcf,14)
end
%% plot vcorr

if(showFigs(showFigs==3))
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
end

%% fourier transform of vcorr

if(showFigs(showFigs==4))
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
end
%% plot contact probability distribution and force distribution
if(isfield(usedSimAm, 'innerForce'))
    
    if(showFigs(showFigs==5))
        figure(5);
        
        for i=1:length(usedSimAm)
            polarhistogram(usedSimAm(i).contactAngs,2*pi/deg2rad(usedSimAm(i).binW));
            hold on;
        end
        title('contact distribution on ring P(\theta)');
    end
    if(showFigs(showFigs==6))
        figure(6);

        for i=1:length(usedSimAm)
            histogram(usedSimAm(i).contactAngs,2*pi/deg2rad(usedSimAm(i).binW));
            hold on;
        end
        xlim([0 2*pi]);
        title('contact distribution on ring P(\theta)');
        xlabel('\theta');
        ylabel('P(\theta)');
    end
    
    if(showFigs(showFigs==7))
        figure(7);
        for i=1:length(usedSimAm)
            polarhistogram(usedSimAm(i).polarForceHist,2*pi/deg2rad(usedSimAm(i).binW))
            hold on;
        end
         title('Force (mN) vs. \theta on ring')
    end
    if(showFigs(showFigs==8))
        figure(8); hold on; 
        for i=1:length(usedSimAm)
            histogram(usedSimAm(i).polarForceHist,2*pi/deg2rad(usedSimAm(i).binW))
            
        end
        xlabel('\theta');
        ylabel('Force (mN)');
        title('Force distribution on ring')
        xlim([0 2*pi]);
    end

end

