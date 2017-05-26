clear all;
close all;
[filez fold]=uigetfile('*.mat','simData','A:\SmarticleRun\SarahAmeobotData');
load(fullfile(fold,filez));

%************************************************************
%* Fig numbers:
%* 1. displacement yvsx
%* 2. MSD
%* 3. v autocorrelation
%* 4. fourier transform of vcorr
%* 5. fit log of mean MSD
%************************************************************


showFigs=[1 2 5];




SPACE_UNITS='unit';
TIME_UNITS='frames';
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

ma = ma.addAll(simTracks);


%% 1 plot displacement yvsx
xx=1;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    ma.plotTracks
    ma.labelPlotTracks
    %     text(0,0+.01,'start')
    plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
% %     legT=cell(1,length(ma.tracks));
%     for i=1:length(ma.tracks)
%         plot(ma.tracks{i}(end,2),ma.tracks{i}(end,3),'ko','markersize',4,'MarkerFaceColor','r');
%         %         leg(i)=h;
%         legT{i}=['v',num2str(usedSimAm(i).pars(3))];
%     end
%     
    
    
    axis tight
    x=get(gca,'xlim');y=get(gca,'ylim');
    c=max(abs(x));
    axis([-c c -c c]);
%     set(gca,'xtick',-c-.05:.1:c-.05,'ytick',-c-.05:.1:c-.05);
    %plot red grid lines
    plot([-c c],[0,0],'r');
    plot([0,0],[-c c],'r');
    figText(gcf,14)
end

%% 2 plot MSD
xx=2;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    ma = ma.computeMSD;
    ma.plotMeanMSD(gca, 1);
        ma.plotMSD(gca);
    [a, b]=ma.fitMeanMSD;
    %     xlim([0 15])
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
        M(i)= mean(ma.vcorr{i}(10:end-nansN,2));
        
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
%% 5. show msd and mean msd, fit log of MSD for POM, and MOP 
xx=5;
if(showFigs(showFigs==xx))
    figure(xx)
    subplot(1,3,1)
    ma.plotMeanMSD(gca,1);
    ma.plotMSD(gca);
    msd=ma.getMeanMSD;
    tend=ma.msd{1}(end,1)*.25;
    tendIdx=find(msd(:,1)<tend,1,'last');
    
    %POM
    
    subplot(1,3,2)
    hold on;
    ma.plotMeanMSD(gca);
%     msd=msd(msd(:,1)<tend,:);
    msd=msd(1:tendIdx,:);
    msd(1,:)=[];
    [POM,gof1]=fit(log(msd(:,1)),log(msd(:,2)),'poly1');
%     plot(POM,log(msd(:,1)),log(msd(:,2)),'o');
    plot(msd(:,1),msd(:,1).^(POM.p1)*exp(POM.p2),'linewidth',2);
    text(.2,.8,['POM=',num2str(mean(POM.p1),3)],'units','normalized','fontsize',20);
    set(gca,'yscale','log','xscale','log')
    
    
    
    %MOP
    subplot(1,3,3)
    hold on;
%     ma.plotMSD(gca,1);
%     MOP=zeros(length(ma),1);
    for i=1:length(ma.tracks)
        msdRun=ma.msd{i}(1:tendIdx,1:2);
        msdRun(1,:)=[];
        [f2,gof2]=fit(log(msdRun(:,1)),log(msdRun(:,2)),'poly1');
        h1=plot(msdRun(:,1),msdRun(:,2),'.');
        plot(msdRun(:,1),msdRun(:,1).^f2.p1*exp(f2.p2),'color',h1.Color,'linewidth',1.5);
        gof2
%         pause
        MOP(i)=f2.p1;
    end
    text(.2,.8,['MOP=',num2str(mean(MOP),3),'\pm',num2str(std(MOP),3)],'units','normalized','fontsize',20);
    set(gca,'yscale','log','xscale','log')
    figText(gcf,20);
end
        
