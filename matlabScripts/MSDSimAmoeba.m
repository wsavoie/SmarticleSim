clear all
fold=uigetdir('A:\SmarticleRun\')
load(fullfile(fold,'amoebaData.mat'));

% load('A:\SmarticleRun\Amoeba_newsquare_1_dead\amoebaData.mat');
close all


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
%* 9. plot contact PD (force col) of alive particles (ONLY FOR 1 RUN)
%*10. plot contact PD (force col) of dead particles (ONLY FOR 1 RUN)
%*11. plot rotation of inactive smarticle
%*12. plot positions of inactive smarticle

%*13. plot fluctuation theorem result
%************************************************************
%%
SPACE_UNITS='m';
TIME_UNITS='s';
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
inds=1;
showFigs=[1 11 12];
useCOM=0;
f=[.2]; rob=[4:5]; v=[14];dirs=[0];

props={f rob v dirs};
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
%% 1 plot displacement yvsx
xx=1;
if(showFigs(showFigs==xx))
    figure(xx)
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

%% 2 plot MSD
xx=2;
if(showFigs(showFigs==xx))
    figure(xx)
    
    ma = ma.computeMSD;
    ma.plotMeanMSD(gca, true);
    ma.plotMSD;
    [a b]=ma.fitMeanMSD;
    figText(gcf,14)
end
%% 3 plot vcorr
xx=3;
if(showFigs(showFigs==xx))
    figure(xx)
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

%% 4 fourier transform of vcorr
xx=4;
if(showFigs(showFigs==xx))
    figure(xx)
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
%% 5 contact P(\theta) POLAR
% if(isfield(usedSimAm, 'innerForce'))
xx=5;
if(showFigs(showFigs==xx))
    figure(xx);
    
    for i=1:length(usedSimAm)
        polarhistogram(usedSimAm(i).contactAngs,2*pi/deg2rad(usedSimAm(i).binW),'binLimits',[0,2*pi]);
        hold on;
    end
    title('contact distribution on ring P(\theta)');
end
%% 6 contact P(\theta) linear
xx=6;
if(showFigs(showFigs==xx))
    figure(xx);
    
    for i=1:length(usedSimAm)
        histogram(usedSimAm(i).contactAngs,2*pi/deg2rad(usedSimAm(i).binW));
        hold on;
    end
    xlim([0 2*pi]);
    title('contact distribution on ring P(\theta)');
    xlabel('\theta');
    ylabel('P(\theta)');
end
%% 7 inner force vs. theta POLAR
xx=7;
if(showFigs(showFigs==xx))
    h=figure(xx);
    for i=1:length(usedSimAm)
        binW=deg2rad(usedSimAm(i).binW);
        polarhistogram('binLimits',[0,2*pi],'BinEdges',[0:binW:2*pi],'BinCounts',[usedSimAm(i).fcs(:,2)])
        %             polarhistogram(usedSimAm(i).polarForceHist,2*pi/deg2rad(usedSimAm(i).binW),'binLimits',[0,2*pi])
        hold on;
    end
    title('Force (mN) vs. \theta on ring')
end
%% 8 inner force vs. theta LINEAR
xx=8;
if(showFigs(showFigs==xx))
    figure(xx); hold on;
    for i=1:length(usedSimAm)
        histogram(usedSimAm(i).polarForceHist,2*pi/deg2rad(usedSimAm(i).binW))
    end
    xlabel('\theta');
    ylabel('Force (mN)');
    title('Force distribution on ring')
    xlim([0 2*pi]);
end

% end
%% 9 plot contact probability distribution colored by force for alive particles
xx=9;
if(showFigs(showFigs==xx))
    figure(xx);
    matlab.graphics.internal.prepareCoordinateSystem('polar', []);
    colormap fire; cb=colorbar; hold on;
    ylabel(cb,'Normalized Total Force');
    cmax=max(usedSimAm(:).fcs(:,2));
    caxis([0,cmax]);
    ca=caxis;
    cbins=100;
    cols=fire(cbins);
    for i=1:length(usedSimAm)
        %get binW, Nbins, binVec for ang matrix
        binW=deg2rad(usedSimAm(i).binW); Nbins=2*pi/binW; binVec=linspace(0,2*pi,Nbins+1);
        
        %discretize data by binVec and mod so it is between 0 and 2pi
        [angDisc]=discretize(usedSimAm(i).contactAngs,binVec);
        
        
        %         [frcDisc]=round(usedSimAm(i).deadFrc(:,1)/binW);
        frcDisc=discretize(usedSimAm(i).fcs(:,2)/cmax,linspace(0,1,cbins));
        %get unique bins from data
        unicc=unique(angDisc);
        if(unicc~=(usedSimAm(i).fcs(:,1)/binW))
            error('unicc and deadfcs are diff');
        end
        for j=1:length(unicc)
            polarhistogram(binW*unicc(j)*ones(nnz(angDisc==unicc(j)),1)-binW+.001,Nbins,'binLimits',[0,2*pi],'facecolor',cols(frcDisc(j),:),'binwidth',binW,'binedges',binVec);
        end
        title('contact distribution on ring P(\theta) (active)');
        %         or just do
        %         figure(23);  polarhistogram(usedSimAm(i).deadContactAngs,Nbins,'binLimits',[0,2*pi]);
        %           hold on;  title('contact distribution on ring P(\theta) (active)'); colorbar;
        
    end
    title('contact distribution on ring P(\theta) (active)');
end
%% 10 plot contact probability distribution colored by force of dead particle
xx=10;
if(showFigs(showFigs==xx))
    figure(xx);
    matlab.graphics.internal.prepareCoordinateSystem('polar', []);
    colormap fire; cb=colorbar; hold on;
    ylabel(cb,'Normalized Total Force');
    cmax=max(usedSimAm(:).deadFcs(:,2));
    caxis([0,1]);
    ca=caxis;
    cbins=100;
    cols=fire(cbins);
    for i=1:length(usedSimAm)
        %get binW, Nbins, binVec for ang matrix
        binW=deg2rad(usedSimAm(i).binW); Nbins=2*pi/binW; binVec=linspace(0,2*pi,Nbins+1);
        
        %discretize data by binVec and mod so it is between 0 and 2pi
        [angDisc]=discretize(usedSimAm(i).deadContactAngs,binVec);
        
        
        %         [frcDisc]=round(usedSimAm(i).deadFrc(:,1)/binW);
        [frcDisc]=discretize(usedSimAm(i).deadFcs(:,2)/cmax,linspace(0,1,cbins));
        %get unique bins from data
        unicc=unique(angDisc);
        if(unicc~=(usedSimAm(i).deadFcs(:,1)/binW))
            error('unicc and deadfcs are diff');
        end
        for j=1:length(unicc)
            polarhistogram(binW*unicc(j)*ones(nnz(angDisc==unicc(j)),1)-binW+.001,Nbins,'binLimits',[0,2*pi],'facecolor',cols(frcDisc(j),:),'binwidth',binW,'binedges',binVec);
        end
        title('contact distribution on ring P(\theta) (inactive)');
        %         or just do
        figure(25);  polarhistogram(usedSimAm(i).deadContactAngs,Nbins,'binLimits',[0,2*pi]);
        hold on; title('contact distribution on ring P(\theta) (inactive)'); colorbar;
        
        
    end
    
end
%% 11 plot rotation of inactive smarticle
xx=11;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    for i=1:length(usedSimAm)
        plot(usedSimAm(i).deadInnerForce(:,1),mod(rad2deg(usedSimAm(i).deadRot),180));
        
    end
    title('inactive smarticle rotation about z-axis');
end
%% 12 plot position of inactive smarticle
xx=12;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    for i=1:length(usedSimAm)
        %         plot(usedSimAm(i).deadInnerForce(:,1),mod(rad2deg(usedSimAm(i).deadRot),180));
        plot(usedSimAm(i).deadPos(:,2),usedSimAm(i).deadPos(:,3));
        plot(usedSimAm(i).deadPos(1,2),usedSimAm(i).deadPos(1,3),'k.','markersize',15);
    end
    minR=0.05; 
    plot(minR,0,'g.','markersize',15);
    cols=get(groot,'defaultaxescolororder');
    viscircles([0,0],usedSimAm(1).r,'color',cols(4,:),'linewidth',4);
    axis equal
    title('inactive smarticle positions about z-axis');
end
%% 13 flucuation
xx=13;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    %random min distance can change
    minR=0.05; 
    maxR=0.0653;
    for i=1:length(usedSimAm)
%         get all where |r| from center > minR 
        N = sqrt(sum(abs(usedSimAm(i).deadPos(:,[2,3])).^2,2));
        t=(usedSimAm(i).deadPos(N>minR))
        pp=usedSimAm(i).deadPos(norm(usedSimAm(i).deadPos(:,2,3))
        plot(usedSimAm(i).deadPos(:,2),usedSimAm(i).deadPos(:,3));
        plot(usedSimAm(i).deadPos(1,2),usedSimAm(i).deadPos(1,3),'k.','markersize',15);
    end
    
    cols=get(groot,'defaultaxescolororder');
    viscircles([0,0],usedSimAm(1).r,'color',cols(4,:),'linewidth',4);
    axis equal
    title('inactive smarticle positions about z-axis');
end