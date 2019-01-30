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
%* 6. just do mop plot of MSD
%* 7. gamma vs. N
% 8. A vs. P, showing sqrt relationship
%************************************************************


showFigs=[1 2 5 6 7 8];

showFigs=[6];

if(exist('ma')==0)
    SPACE_UNITS='unit';
    TIME_UNITS='frames';
    ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
    ma = ma.addAll(simTracks);
end


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
    if(isempty(ma.msd))
        ma = ma.computeMSD;
    end
    %     ma = ma.computeMSD;
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
    hold on;
    subplot(1,3,1)
    if(isempty(ma.msd))
        ma = ma.computeMSD;
    end
    scale=1/sqrt(N);
    ma.plotMeanMSD(gca,1);
    ma.plotMSD(gca);
    msd=ma.getMeanMSD;
    tend=ma.msd{1}(end,1)*.25*scale;
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
    %     for i=1:length(ma.tracks)
    %         msdRun=ma.msd{i}(1:tendIdx,1:2);
    %         msdRun(1,:)=[];
    %         [f2,gof2]=fit(log(msdRun(:,1)),log(msdRun(:,2)),'poly1');
    %         h1=plot(msdRun(:,1),msdRun(:,2),'.');
    %         plot(msdRun(:,1),msdRun(:,1).^f2.p1*exp(f2.p2),'color',h1.Color,'linewidth',1.5);
    %         gof2
    % %         pause
    %         MOP(i)=f2.p1;
    %     end
    ma=ma.fitLogLogMSD(.25*scale);
    llfit=ma.loglogfit;
    for i=1:length(ma.tracks)
        msdRun=ma.msd{i}(1:tendIdx,1:2);
        msdRun(1,:)=[];
        h1=plot(msdRun(:,1),msdRun(:,2),'.');
        plot(msdRun(:,1),msdRun(:,1).^(llfit.alpha(i)).*(llfit.gamma(i)),'color',h1.Color,'linewidth',1.5);
        %        plot(msdRun(:,1),msdRun(:,1).^f2.p1*exp(f2.p2),'color',h1.Color,'linewidth',1.5);
    end
    set(gca,'yscale','log','xscale','log')
    figText(gcf,20);
    ma.fitMeanMSD;
    MOP=llfit.alpha;
    
    text(.2,.8,['MOP=',num2str(mean(MOP),3),'\pm',num2str(std(MOP),3)],'units','normalized','fontsize',20);
    set(gca,'yscale','log','xscale','log')
    figText(gcf,20);
end

%% 6. just MOP
xx=6;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    if(isempty(ma.msd))
        ma = ma.computeMSD;
    end
    scale=1/(sqrt(N));
    ma.plotMeanMSD(gca,1);
    ma.plotMSD(gca);
    msd=ma.getMeanMSD;
    tend=ma.msd{1}(end,1)*.25*scale;
    tendIdx=find(msd(:,1)<tend,1,'last');
    
    
    ma=ma.fitLogLogMSD(.25*scale);
    llfit=ma.loglogfit;
    for i=1:length(ma.tracks)
        msdRun=ma.msd{i}(1:tendIdx,1:2);
        msdRun(1,:)=[];
        h1=plot(msdRun(:,1),msdRun(:,2),'.');
        plot(msdRun(:,1),msdRun(:,1).^(llfit.alpha(i)).*(llfit.gamma(i)),'color',h1.Color,'linewidth',1.5);
        %        plot(msdRun(:,1),msdRun(:,1).^f2.p1*exp(f2.p2),'color',h1.Color,'linewidth',1.5);
    end
    set(gca,'yscale','log','xscale','log')
    figText(gcf,20);
    ma.fitMeanMSD;
    MOP=llfit.alpha;
    
    text(.2,.8,['MOP=',num2str(mean(MOP),3),'\pm',num2str(std(MOP,1),3)],'units','normalized','fontsize',20);
    set(gca,'yscale','log','xscale','log')
    figText(gcf,20);
end
%% 7. filled in values from data to plot gamma vs. N
xx=7;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    
    x=        [19   37   61   91   127   169   217];
    runs=sqrt([9    10   10   10   9     8     9]);
    
    %with just 1/N
    y=[.981,0,0,0,0,0,0];
    err=[0.0697,0,0,0,0,0,0];

    
    %3.47/(sqrt(N))
%     y=[.967, 1.02,.898,.927,.944,.934,.919];
%     err=[.0967,.134,.0964,.0608,.0447,0.054,.0559];
    %     err=[,.0634,.047,.0255,.0387,.0295,.0257];
    
    %1/(sqrt(N)*3.47)
    y=[.971,.953,.931,.938,.947,.934,.943];
    err=[.0359,.0261,0.0247,.0198,0.0174,0.0168,.0203];
    
    errorbar(x,y,err./sqrt(runs));
    figText(gcf,20);
    xlabel('N');
    ylabel('\gamma');
end
%% 8. A vs. P, showing sqrt relationship
xx=8;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    warning('off','curvefit:fit:noStartPoint')
    R=[1:40];
    A= 3.*R.*(R+1)+1;
    P=6*R;
    
    %     plot(A,P);
    ft= fittype('a*x.^(0.5)+b',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b'});
    f = fit(A',P',ft);
    %     plot(f,'k-','linewidth','2');
    %plot fitted data
    plot(A,f.a*A.^(1/2)+f.b,'k-','linewidth',2);
    plot(A,P,'o','linewidth',2,'markerfacecolor','w');
    plot(A,3.47*sqrt(A));
    
    legend({'fit with 2 free vars','data','fit with 1 free var'},'location','NW')
    text(.2,.2,{['$P(N)=a\sqrt{N}+b$'],...
        ['$a\approx',num2str(f.a,3),'\quad b\approx',num2str(f.b,3),'$']},...
        'units','normalized','Interpreter','latex'); 
    ylabel('Perimeter (N)');
    xlabel('Area (N)');
    figText(gcf,16);
end
%#ok<*NBRAK>