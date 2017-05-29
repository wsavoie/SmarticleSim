clear all
fold=uigetdir('A:\SmarticleRun\')
load(fullfile(fold,'amoebaData.mat'));

% load('A:\SmarticleRun\Amoeba_newsquare_1_dead\amoebaData.mat');
close all


%************************************************************
%* Fig numbers:
%* 1. displacement yvsx
%* 2. MSD
%* 3. v autocorrelation
%* 4. fourier transform of vcorr
%* 5. contact P(\theta) POLAR
%* 6. contact P(\theta) LINEAR
%* 7. inner force vs. theta POLAR
%* 8. inner force vs. theta LINEAR
%* 9. contact PD (force col) of alive particles (ONLY FOR 1 RUN)
%*10. contact PD (force col) of dead particles (ONLY FOR 1 RUN)
%*11. rotation of inactive smarticle
%*12. positions of inactive smarticle
%*13. v.rhat histogram
%*14. v.r histogram
%*15. v dot fixed final r
%*16. v dot fixed final r full data
%*17. v dot normalized direct r full data
%*18. /19 fluctuation theorem based plot must also use 17!!
%*20. plot vectors at each position relating position of deadparticle
%*21. v crosscorr with deadpos
%*22. for each msd traj get linear fit of log MOP AND POM
%*23. histogram of probability(theta_R-theta_r) make vid
%*24. track length plotting
%*25. rotational track length
%*26. smarticle rotation
%*27 v dot rhat dt sampling histogram
%*28 v dot rhat collision event sampling histogram
%*29 v dot rhat dt sampling scatter
%*30 v dot rhat scatter collision event sampling
%*31 v dot rhat heatmap dt sampling
%*32 v dot rhat heatmap collision event sampling
%*33 log(pdf) vs |vel|
%*34 %% 34 log(pdf) vs |vel| fit line
%*45. test
%************************************************************
%%
SPACE_UNITS='m';
TIME_UNITS='s';
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
inds=1;
showFigs=[1 22];
useCOM=0;
f=[]; rob=[]; v=[];dirs=[0];

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
L=length(usedSimAm);
%% 1 plot displacement yvsx
xx=1;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    %     if(useCOM)
    %         title('Ring COM');
    %     else
    %         title('Ring COG');
    %     end
    ma.plotTracks
    ma.labelPlotTracks
    %     text(0,0+.01,'start')
    plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
    legT=cell(1,length(ma.tracks));
    for i=1:length(ma.tracks)
        plot(ma.tracks{i}(end,2),ma.tracks{i}(end,3),'ko','markersize',4,'MarkerFaceColor','r');
        %         leg(i)=h;
        legT{i}=['v',num2str(usedSimAm(i).pars(3))];
    end
    
    
    
    axis tight
    x=get(gca,'xlim');y=get(gca,'ylim');
    c=max(abs(x));
    if c<=.25
        c=.25;
    elseif c<=.35
    c=.35;            
    else
        c=.45;
    end
    axis([-c c -c c]);
    set(gca,'xtick',-c-.05:.1:c-.05,'ytick',-c-.05:.1:c-.05);
    
    %plot red grid lines
    plot([-c c],[0,0],'r');
    plot([0,0],[-c c],'r');
%     legend(legT);
    r=simAm(1).r;
    t = linspace(0,2*pi);plot(r*cos(t),r*sin(t),'-k','linewidth',2);
    figText(gcf,14)
end

%% 2 plot MSD
xx=2;
if(showFigs(showFigs==xx))
    figure(xx)
    
    ma = ma.computeMSD;
    ma.plotMeanMSD(gca, true);
    %     ma.plotMSD;
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
    set(gca,'yscale','log')
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
    
    for i=1:L
        polarhistogram(usedSimAm(i).contactAngs,2*pi/deg2rad(usedSimAm(i).binW),'binLimits',[0,2*pi]);
        hold on;
    end
    title('contact distribution on ring P(\theta)');
end
%% 6 contact P(\theta) linear
xx=6;
if(showFigs(showFigs==xx))
    figure(xx);
    
    for i=1:L
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
    for i=1:L
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
    for i=1:L
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
    ylabel(cb,'Total Force (N)');
    cmax=max(usedSimAm(:).fcs(:,2));
    caxis([0,cmax]);
    ca=caxis;
    cbins=100;
    cols=fire(cbins);
    for i=1:L
        %get binW, Nbins, binVec for ang matrix
        binW=deg2rad(usedSimAm(i).binW); Nbins=2*pi/binW; binVec=linspace(0,2*pi,Nbins+1);
        
        %discretize data by binVec and mod so it is between 0 and 2pi
        [angDisc]=discretize(usedSimAm(i).contactAngs,binVec);
        
        
        %         [frcDisc]=round(usedSimAm(i).deadFrc(:,1)/binW);
        frcDisc=discretize(usedSimAm(i).fcs(:,2),linspace(0,cmax,cbins));
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
    ylabel(cb,'Total Force (N)');
    cmax=max(usedSimAm(:).deadFcs(:,2));
    caxis([0,cmax]);
    ca=caxis;
    cbins=100;
    cols=fire(cbins);
    for i=1:L
        %get binW, Nbins, binVec for ang matrix
        binW=deg2rad(usedSimAm(i).binW); Nbins=2*pi/binW; binVec=linspace(0,2*pi,Nbins+1);
        
        %discretize data by binVec and mod so it is between 0 and 2pi
        [angDisc]=discretize(usedSimAm(i).deadContactAngs,binVec);
        
        
        %         [frcDisc]=round(usedSimAm(i).deadFrc(:,1)/binW);
        [frcDisc]=discretize(usedSimAm(i).deadFcs(:,2),linspace(0,cmax,cbins));
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
        %         figure(25);  polarhistogram(usedSimAm(i).deadContactAngs,Nbins,'binLimits',[0,2*pi]);
        %         hold on; title('contact distribution on ring P(\theta) (inactive)'); colorbar;
        
        
    end
    
end
%% 11 plot rotation of inactive smarticle
xx=11;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    for i=1:L
        plot(usedSimAm(i).deadInnerForce(:,1),mod(rad2deg(usedSimAm(i).deadRot),180));
        
    end
    title('inactive smarticle rotation about z-axis');
end
%% 12 plot position of inactive smarticle
xx=12;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    
    leg=zeros(L,1);
    legT={};
    
    for i=1:L
        %         plot(usedSimAm(i).deadInnerForce(:,1),mod(rad2deg(usedSimAm(i).deadRot),180));
        %         plot(usedSimAm(i).deadPos(:,2),usedSimAm(i).deadPos(:,3));
        %         plot(usedSimAm(i).deadPos(1,2),usedSimAm(i).deadPos(1,3),'k.','markersize',15);
        %         plot(usedSimAm(i).deadPos(end,2),usedSimAm(i).deadPos(end,3),'r.','markersize',15);
        %
        deadPos=usedSimAm(i).fullDeadSmartPos-usedSimAm(i).fullRingPos;
        h=plot(deadPos(:,1),deadPos(:,2));
        
        plot(deadPos(end,1),deadPos(end,2),'r.','markersize',15);
        leg(i)=h;
        x=['v',num2str(usedSimAm(i).pars(3))];
        legT{i}=['v',num2str(usedSimAm(i).pars(3))];
    end
    minR=0.05;
    %     plot(minR,0,'g.','markersize',15);
    cols=get(groot,'defaultaxescolororder');
    viscircles([0,0],usedSimAm(1).r,'color',cols(4,:),'linewidth',4);
    plot(usedSimAm(1).fullDeadSmartPos(1,1),usedSimAm(1).fullDeadSmartPos(1,2),'k.','markersize',15);
    axis equal
    title('inactive smarticle position relative to ring');
    legend(leg,legT);
end
%% 13 v dot rhat
xx=13;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    %random min distance can change
    minR=0.05;
    maxR=0.0653;
    dotProd=[];
    
    for i=1:L
        %         get all where |r| from center > minR
        %         N = sqrt(sum(abs(usedSimAm(i).deadPos(:,[2,3])).^2,2));
        %         t=(usedSimAm(i).deadPos(N>minR));
        
        [deadPos,idx]=unique(usedSimAm(i).deadPos,'rows');
        ringPos=usedSimAm(i).ringpos(idx,:);
        vRing=diff(ringPos)./diff(deadPos(:,1));
        vRing=(vRing);
        
        dotProd=[dotProd; dot(vRing,normr(deadPos(1:end-1,2:3)),2)];
        %         plot(usedSimAm(i).deadPos(1,2),usedSimAm(i).deadPos(1,3),'k.','markersize',15);
    end
    histogram(dotProd);
    %     figure(23)
    %     title('\vec{v_ring}\cdot\vec{r_inactive}');
    title('$\vec{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    xlabel('$\vec{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('counts');
    text(.1,.8,['mean = ',num2str(mean(dotProd))],'units','normalized');
end

%% 14 flucuation v dot r
xx=14;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[];
    for i=1:L
        %         get all where |r| from center > minR
        %         N = sqrt(sum(abs(usedSimAm(i).deadPos(:,[2,3])).^2,2));
        %         t=(usedSimAm(i).deadPos(N>minR));
        [deadPos,idx]=unique(usedSimAm(i).deadPos,'rows');
        ringPos=usedSimAm(i).ringpos(idx,:);
        vRing=diff(ringPos)./diff(deadPos(:,1));
        dotProd=[dotProd; dot(vRing,deadPos(1:end-1,2:3),2)];
    end
    histogram(dotProd);
    dotProd(isnan(dotProd))=[];
    %     figure(23)
    %     plotv(vRing(1:50:1000,:)')
    title('$\vec{v}_{ring} \cdot\vec{r}_{inactive}$','interpreter', 'latex');
    xlabel('$\vec{v}_{ring} \cdot\vec{r}_{inactive}$','interpreter', 'latex');
    ylabel('counts');
    text(.1,.8,['mean = ',num2str(mean(dotProd))],'units','normalized');
end
%% 15  flucuation v dot fixed final r
xx=15;
if(showFigs(showFigs==xx))
    figure(xx)
    %     ti=100;
    %     te=110;
    hold on;
    %random min distance can change
    dotProd=[];
    for i=1:L
        %         get all where |r| from center > minR
        %         N = sqrt(sum(abs(usedSimAm(i).deadPos(:,[2,3])).^2,2));
        %         t=(usedSimAm(i).deadPos(N>minR));
        
        [deadPos,idx]=unique(usedSimAm(i).deadPos,'rows');
        ringPos=usedSimAm(i).ringpos(idx,:);
        
        %         ringPos=ringPos(ti:te,:);
        %         deadPos=deadPos(ti:te,:);
        %
        vRing=diff(ringPos)./diff(deadPos(:,1));
        vRing=(vRing);
        deadPos=repmat([ringPos(end,:)],[size(vRing,1),1]);
        
        %         dotProd=[dotProd; dot(vRing,normr(deadPos(1:end-1,2:3)),2)];
        dotProd=[dotProd; dot(vRing,normr(deadPos(1:end,1:2)),2)];
        %         plot(usedSimAm(i).deadPos(1,2),usedSimAm(i).deadPos(1,3),'k.','markersize',15);
    end
    histogram(dotProd);
    %     figure(23)
    %     title('\vec{v_ring}\cdot\vec{r_inactive}');
    title('$\vec{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    xlabel('$\vec{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('counts');
    text(.1,.8,['mean = ',num2str(mean(dotProd))],'units','normalized');
end

%% 16 v dot fixed final r full data
xx=16;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[];
    for i=1:L
        %         get all where |r| from center > minR
        %         N = sqrt(sum(abs(usedSimAm(i).deadPos(:,[2,3])).^2,2));
        %         t=(usedSimAm(i).deadPos(N>minR));
        
        %         deadPos=usedSimAm(i).fullDeadSmartPos;
        ringPos=usedSimAm(i).fullRingPos;
        vRing=diff(ringPos)./diff(usedSimAm(i).fullT);
        %         vRing=(vRing);
        deadPos=repmat([ringPos(end,:)],[size(vRing,1),1]);
        
        %         dotProd=[dotProd; dot(vRing,normr(deadPos(1:end-1,2:3)),2)];
        dotProd=[dotProd; dot(vRing,normr(deadPos),2)];
        %         plot(usedSimAm(i).deadPos(1,2),usedSimAm(i).deadPos(1,3),'k.','markersize',15);
    end
    histogram(dotProd);
    %     figure(23)
    %     title('\vec{v_ring}\cdot\vec{r_inactive}');
    title('$\vec{v}_{ring} \cdot\hat{r}_{inactive_{final}}$','interpreter', 'latex');
    xlabel('$\vec{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('counts');
    text(.1,.8,['mean = ',num2str(mean(dotProd))],'units','normalized');
    %
    %     pts('final pos=(',deadPos(1,1),',',deadPos(1,2),')');
    %     pts('norm of pos=',norm(deadPos(1,:)),')');
end
%% 17 v dot normalized direct r full data
xx=17;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[];
    dotProd2t=[];
    for i=1:L
        %         get all where |r| from center > minR
        %         N = sqrt(sum(abs(usedSimAm(i).deadPos(:,[2,3])).^2,2));
        %         t=(usedSimAm(i).deadPos(N>minR));
        
        deadPos=usedSimAm(i).fullDeadSmartPos(1:end,:);
        deadPos2t=usedSimAm(i).fullDeadSmartPos(1:2:end,:);
        
        ringPos=usedSimAm(i).fullRingPos(1:end,:);
        vRing=diff(ringPos)./diff(usedSimAm(i).fullT(1:end,:));
        
        
        ringPos2t=usedSimAm(i).fullRingPos(1:2:end,:);
        vRing2t=diff(ringPos2t)./diff(usedSimAm(i).fullT(1:2:end,:));
        
        dotProd=[dotProd; dot(vRing,normr(deadPos(1:end-1,1:2)),2)];
        dotProd2t=[dotProd2t; dot(vRing2t,normr(deadPos2t(1:end-1,1:2)),2)];
        
    end
    histogram(dotProd,75);
    %     figure(23)
    %     title('\vec{v_ring}\cdot\vec{r_inactive}');
    title('$\vec{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    xlabel('$\vec{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('counts');
    text(.1,.8,['mean = ',num2str(mean(dotProd))],'units','normalized');
    
    %     pts('final pos=(',deadPos(1,1),',',deadPos(1,2),')');
    %     pts('norm of pos=',norm(deadPos(1,:)),')');
    set(gca,'yscale','log');
end
%% 18 19 fluctuation theorem based plot must also use 17!!
xx=18;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    subplot(2,2,1)
    
    dp=histogram(dotProd);
    s=dp.Values;
    s(s==0)=.001;
    title('histogram of $\vec{v}_{ring}\cdot\vec{r}$','interpreter','latex');
    
    subplot(2,2,2);
    hold on;
    title('histogram values of $\vec{v}_{ring}\cdot\vec{r}$','interpreter','latex');
    [~,mind]=max(s);
    cleft=s(1:mind);
    cright=s(mind:end);
    cl=log(cleft);
    cr=log(cright);
    plot(1:length(cl),(cl),'g')
    plot(length(cl):length(s),cr,'r')
    
    subplot(2,2,3);
    hold on;
    title('slopes of left and right side of $\vec{v}_{ring}\cdot\vec{r}$','interpreter','latex');
    indL=round(length(cl)*.4);
    plot(cl(indL:end-1))
    x2=(indL:length(cl));
    [a2 b2]=polyfit([x2]',cl(x2)',1);
    plot(x2-x2(1),a2(1)*(x2)+a2(2),'g','linewidth',2);
    indR=round(length(cr)*.6);
    [a b]=polyfit([1:length(cr(2:indR))]',cr(2:indR)',1);
    plot(cr(2:indR))
    x=1:length(cr(2:indR));
    plot(x,a(1)*(x)+a(2),'r','linewidth',2);
    
    subplot(2,2,4);
    hold on;
    plot(x,a(1)*(x)+a(2),'g','linewidth',1);
    plot(x2-x2(1),-a2(1)*(x2-x2(1))+a(2),'r','linewidth',1);
    pts('left + right= (',a2(1),')+ (',a(1),') =',a2(1)+a(1));
    suptitle('window = t');
    
    %2t figure
    figure(xx+1)
    hold on;
    subplot(2,2,1)
    hold on;
    title('histogram of $\vec{v}_{ring}\cdot\vec{r}$','interpreter','latex');
    dp_2t=histogram(dotProd2t);
    s_2t=dp_2t.Values;
    s_2t(s_2t==0)=.001;
    
    subplot(2,2,2);
    hold on
    title('histogram values of $\vec{v}_{ring}\cdot\vec{r}$','interpreter','latex');
    [~,mind_2t]=max(s_2t);
    cleft_2t=s_2t(1:mind_2t);
    cright_2t=s_2t(mind_2t:end);
    cleft_2t(cleft_2t==0)=.01;
    cright_2t(cright_2t==0)=.01;
    cl_2t=log(cleft_2t);
    cr_2t=log(cright_2t);
    plot(1:length(cl_2t),(cl_2t),'g')
    plot(length(cl_2t):length(s_2t),cr_2t,'r')
    
    subplot(2,2,3);
    hold on;
    title('slopes of left and right side of $\vec{v}_{ring}\cdot\vec{r}$','interpreter','latex');
    indL_2t=round(length(cl_2t)*.4);
    plot(cl_2t(indL_2t:end))
    x2_2t=(indL_2t:length(cl_2t));
    [a2_2t, ~]=polyfit([x2_2t]',cl_2t(x2_2t)',1);
    plot(x2_2t-x2_2t(1),a2_2t(1)*(x2_2t)+a2_2t(2),'g','linewidth',2);
    indR_2t=round(length(cr_2t)*.6);
    [a_2t, ~]=polyfit([1:length(cr_2t(2:indR_2t))]',cr_2t(2:indR_2t)',1);
    plot(cr_2t(2:indR_2t))
    x_2t=1:length(cr_2t(2:indR_2t));
    plot(x_2t,a_2t(1)*(x_2t)+a_2t(2),'r','linewidth',2);
    
    subplot(2,2,4);
    hold on;
    plot(x_2t,a_2t(1)*(x_2t)+a_2t(2),'r','linewidth',1);
    plot(x2_2t-x2_2t(1),-a2_2t(1)*(x2_2t-x2_2t(1))+a_2t(2),'g','linewidth',1);
    pts('left2t + right2t= (',a2_2t(1),')+ (',a_2t(1),') =',a2_2t(1)+a_2t(1));
    suptitle('window = 2t');
end

%% 20 plot vectors at each position relating position of deadparticle
xx=20;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    NE=L;
    xp=-1;
    yp=-1;
    dec=500;
    tick=.25;
    leg=[];
    legT={};
    for i=1:L
        
        xp=mod(xp+1,5);
        if(xp==0)
            yp=yp+1;
        end
        %         pts('xp=',xp,' yp=',yp);
        XP=xp*tick;
        YP=yp*tick;
        
        %         deadPos=normr(usedSimAm(i).fullDeadSmartPos-repmat(usedSimAm(i).fullRingPos(1,:),[size(usedSimAm(i).fullRingPos,1),1]));
        %         ringPos=usedSimAm(i).fullRingPos-repmat(usedSimAm(i).fullRingPos(1,:),[size(usedSimAm(i).fullRingPos,1),1]);
        deadPos=normr(usedSimAm(i).fullDeadSmartPos);
        ringPos=usedSimAm(i).fullRingPos;
        h=plot(ringPos(:,1)+XP,ringPos(:,2)+YP,'linewidth',1.5);
        quiver(ringPos(1:dec:end,1)+XP,ringPos(1:dec:end,2)+YP,deadPos(1:dec:end,1),deadPos(1:dec:end,2),'color','k')
        plot(ringPos(1,1)+XP,ringPos(1,2)+YP,'.k','markersize',15);
        plot(ringPos(end,1)+XP,ringPos(end,2)+YP,'.r','markersize',15);
        
        leg=[leg h];
        x=['v',num2str(usedSimAm(i).pars(3))];
        legT{i}=['v',num2str(usedSimAm(i).pars(3))];
    end
    ax=gca;
    ax.XTick=[-.25:tick:max(xticks)];
    ax.YTick=[-.25:tick:max(yticks)];
    grid on;
    legend(leg,legT);
end

%% 21 vcorr
xx=21;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dt=diff(usedSimAm(1).fullT(1:2));
    dec=1;
    pks=[];
    for i=1:L
        %         deadPos=normr(usedSimAm(i).fullDeadSmartPos-repmat(usedSimAm(i).fullRingPos(1,:),[size(usedSimAm(i).fullRingPos,1),1]));
        %         ringPos=usedSimAm(i).fullRingPos-repmat(usedSimAm(i).fullRingPos(1,:),[size(usedSimAm(i).fullRingPos,1),1]);
        deadPos=normr(usedSimAm(i).fullDeadSmartPos(1:end-1,:));
        ringPos=usedSimAm(i).fullRingPos;
        vRing=diff(ringPos)./diff(usedSimAm(i).fullT(1:end,:));
        vRing=normr(vRing);
        deadPos=normr(deadPos);
        
        vTheta=atan2(vRing(:,2),vRing(:,1))+pi;
        dTheta=atan2(deadPos(:,2),deadPos(:,1))+pi;
        vTheta=decimate(vTheta,dec);
        dTheta=decimate(dTheta,dec);
        [Rmm,lags]=xcorr(vTheta,dTheta,55*dec/dt,'unbiased');
        lags=lags*dt/dec;
        %                 Rmm=Rmm(lags>3);
        %                 lags=lags(lags>3);
        %                 Rmm=Rmm(lags<20);
        %                 lags=lags(lags<20);
        Rmm=Rmm/max(abs(Rmm));
        subplot(1,2,1);
        hold on; title('cross correlation (v_{ring}*r_{inactive})');
        xlabel('Lag (s)');ylabel('normalized correlation');
        h=plot(lags,Rmm);
        [pksY,pksX]=findpeaks(abs(Rmm),lags,'MinPeakDistance',10,'minpeakheight',.75,'SortStr','descend');
        %         [pksY,pksX]=findpeaks(abs(Rmm),lags,'MinPeakDistance',10,'SortStr','descend');
        pksY=Rmm(round(lags,5)==round(pksX(1),5));
        pks(i)=pksX(1);
        plot(pksX(1),pksY,'^','color',h.Color,'markerfacecolor',h.Color)
        %         pause;
        %         [Rmm]=xcorr2(vRing);
        
    end
    subplot(1,2,2);
    hold on;
    title('max(abs(crosscorr)) ');
    xlabel('run number');
    ylabel('time lag(s)');
    plot(pks,'-o');
    plot(abs(pks),'.-','markersize',15);
    set(gca,'xtick',v);
    %     [Rmm,lags]=xcorr(vTheta,dTheta);
    %     %         [Rmm]=xcorr2(vRing);
    %     Rmm=Rmm/max(abs(Rmm));
    %     lags=lags*dt*dec;
    %     plot(lags,Rmm);
    %     xlabel('Lag (s)');
    %     [a,b]=findpeaks(abs(Rmm),lags,'MinPeakDistance',10,'minpeakheight',.75)
    
    %         figure(22);
    %          [Rmm,lags]=xcorr(vTheta);
    % %         [Rmm]=xcorr2(vRing);
    %         plot(lags,Rmm/max(Rmm))
    %         xlabel('Lag (s)');
end
%% 22 for each msd traj get linear fit of log
xx=22;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
   
    if(isempty(ma.msd))
        ma = ma.computeMSD;
    end
    
    mmsd=ma.getMeanMSD;
    tend=ma.msd{1}(end,1)*.25;
    tendIdx=find(mmsd(:,1)<tend,1,'last');
    
    
    %first plot with all data
    subplot(1,3,1)
    ma.plotMeanMSD(gca,1);
    ma.plotMSD(gca);
    
    %POM
    subplot(1,3,2)
    hold on;
     ma.plotMeanMSD(gca);
    mmsd=mmsd(1:tendIdx,:);
    mmsd(1,:)=[];
    [POM,gof1]=fit(log(mmsd(:,1)),log(mmsd(:,2)),'poly1');
    plot(mmsd(:,1),mmsd(:,1).^(POM.p1)*exp(POM.p2),'linewidth',2);
    text(.2,.8,['POM=',num2str(mean(POM.p1),3)],'units','normalized','fontsize',20);
    set(gca,'yscale','log','xscale','log')
    
    %MOP
    subplot(1,3,3)
    hold on;
  ma=ma.fitLogLogMSD;
  llfit=ma.loglogfit;
  for i=1:length(ma.tracks)
       msdRun=ma.msd{i}(1:tendIdx,1:2);
        msdRun(1,:)=[];
        h1=plot(msdRun(:,1),msdRun(:,2),'.');
        plot(msdRun(:,1),msdRun(:,1).^(llfit.alpha(i)).*(llfit.gamma(i)),'color',h1.Color,'linewidth',1.5);
  end
%     for i=1:length(ma.tracks)
%         msdRun=ma.msd{i}(1:tendIdx,1:2);
%         msdRun(1,:)=[];
%         [f2,gof2]=fit(log(msdRun(:,1)),log(msdRun(:,2)),'poly1');
%         h1=plot(msdRun(:,1),msdRun(:,2),'.');
%         plot(msdRun(:,1),msdRun(:,1).^f2.p1*exp(f2.p2),'color',h1.Color,'linewidth',1.5);
% %         pause
%         MOP(i)=f2.p1;
%     end
    MOP=llfit.alpha;
    text(.2,.8,['MOP=',num2str(mean(MOP),3),'\pm',num2str(std(MOP),3)],'units','normalized','fontsize',20);
    set(gca,'yscale','log','xscale','log')
    figText(gcf,20);
    ma.fitMeanMSD;

    
end
%% 23 histogram of probability(theta_R-theta_r) make vid
xx=23;
clear ringPos deadPos;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fps=10;
    writeVid=1;
    if(writeVid)
        outputVideo = VideoWriter(fullfile('','prob(theta).avi'));
        outputVideo.FrameRate=fps;
        open(outputVideo);
    end
    for j=[2:1:60]
        for i=1:L
            tt=usedSimAm(i).fullT;
            ringPos(i,:)=normr(usedSimAm(i).fullRingPos(round(tt,5)==j,:));
            deadPos(i,:)=normr(usedSimAm(i).fullDeadSmartPos(round(tt,5)==j,:));
            %             histogra
        end
        clf;
        ringTheta=atan2(ringPos(:,2),ringPos(:,1));
        deadTheta=atan2(deadPos(:,2),deadPos(:,1));
        res=mod(ringTheta-deadTheta+pi,2*pi);
        res=res-pi;
        %         res=atan2(ringPos(:,2)-deadPos(:,2),ringPos(:,1)-deadPos(:,1));
        %         figure(j);
        bins=21;
        edges=linspace(-pi,pi,bins+1);
        histogram(res,'BinEdges',edges,'Normalization','probability');
        
        %axis limits x
        set(gca,'xticklabel',{'-\pi' '-\pi/2' '0' '\pi/2','\pi'},'xtick',[-pi,-pi/2,0,pi/2,pi]);
        xlim([-pi,pi]);
        
        %axis limits y
        set(gca,'yticklabel',{'0' '0.25' '0.5' '0.75','1'},'ytick',[0 0.25 0.5 0.75 1]);
        ylim([0,1]);
        
        
        % title('Probability distribution of \theta_{ring}(t=0) - \theta_{inactive}(t=t_f)');
        %         ylabel('P(\theta_r(t=t_f)-\theta_{is}(t=0))');
        
        xlabel('angle (rads)');
        ylabel('P(\theta_R(t)-\theta_r(0))');
        
        figtxt=text(.1,.9,['t = ',num2str(j),' s'],'units','normalized');
        figText(gcf,15);
        
        
        if(writeVid)
            writeVideo(outputVideo,getframe(gcf));
        else
            pause(1/fps);
        end
        
        
        %        pause(.25);
        
    end
    
    if(writeVid)
        close(outputVideo);
    end
end
%% 24 track length
xx=24;
if(showFigs(showFigs==xx))
    d=75;
    figure(xx)
    a = 1;
    b = 1/d*ones(1,d);
    hold on;
    %     if(useCOM)
    %         title('Ring COM');
    %     else
    %         title('Ring COG');
    %     end
    ma.plotTracks
    ma.labelPlotTracks
    %     text(0,0+.01,'start')
    plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
    leg=zeros(1,length(ma.tracks));
    legT=cell(1,length(ma.tracks));
    CL=zeros(1,length(ma.tracks));
    for i=1:length(ma.tracks)
        
        f=filter(b,a,ma.tracks{i}(:,2:3));
        %add final position to matrix
        f(end,:)=[ma.tracks{i}(end,2),ma.tracks{i}(end,3)];
        
        plot(f(:,1),f(:,2),'linewidth',1.5);
        legT{i}=['v',num2str(usedSimAm(i).pars(3))];
        
        CLF = hypot(diff(f(:,1)), diff(f(:,2)));
        CL(i) = trapz(CLF);                   % Integrate to calculate arc length
    end
    
    for i=1:length(ma.tracks)
        plot(ma.tracks{i}(end,2),ma.tracks{i}(end,3),'ko','markersize',4,'MarkerFaceColor','r');
    end
    
    axis tight
    x=get(gca,'xlim');y=get(gca,'ylim');
    c=max(abs(x));
    if c<=.25
        c=.25;
    else
        c=.45;
    end
    axis([-c c -c c]);
    set(gca,'xtick',-c-.05:.1:c-.05,'ytick',-c-.05:.1:c-.05);
    
    %plot red grid lines
    plot([-c c],[0,0],'r');
    plot([0,0],[-c c],'r');
    legend(legT);
    figText(gcf,14)
    
    figure
    c=errorbar(1,mean(CL),std(CL));
    pts('avg=',c.YData,',    +-=',c.YPositiveDelta);
    
end

%% 25 rotational track length
xx=25;
if(showFigs(showFigs==xx))
    d=75;
    figure(xx)
    a = 1;
    b = 1/d*ones(1,d);
    hold on;
    %     if(useCOM)
    %         title('Ring COM');
    %     else
    %         title('Ring COG');
    %     end
    %     ma.plotTracks
    ma.labelPlotTracks
    %     text(0,0+.01,'start')
    plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
    leg=zeros(1,L);
    legT=cell(1,L);
    CL=zeros(1,L);
    for i=1:L
        
        dp=usedSimAm(i).fullDeadSmartPos-usedSimAm(i).fullRingPos;
        f=filter(b,a,dp);
        %add final position to matrix
        
        f(end,:)=dp(end,:);
        t=usedSimAm(i).deadInnerForce(:,1);
        plot(dp(:,1),dp(:,2),'linewidth',1.5);
        legT{i}=['v',num2str(usedSimAm(i).pars(3))];
        
        CLF = hypot(diff(dp(:,1)), diff(dp(:,2)));
        CL(i) = trapz(CLF);                   % Integrate to calculate arc length
    end
    
    %    for i=1:length(ma.tracks)
    %        plot(ma.tracks{i}(end,2),ma.tracks{i}(end,3),'ko','markersize',4,'MarkerFaceColor','r');
    %    end
    
    %     axis tight
    %  	x=get(gca,'xlim');y=get(gca,'ylim');
    %     c=max(abs(x));
    %     if c<=.25
    %         c=.25;
    %     else
    %         c=.45;
    %     end
    %     axis([-c c -c c]);
    %     set(gca,'xtick',-c-.05:.1:c-.05,'ytick',-c-.05:.1:c-.05);
    
    %     %plot red grid lines
    %     plot([-c c],[0,0],'r');
    %     plot([0,0],[-c c],'r');
    %     legend(legT);
    %     figText(gcf,14)
    %
    figure(1000)
    c=errorbar(1,mean(CL),std(CL));
    %     close(1000);
    pts('avg=',c.YData,',    +-=',c.YPositiveDelta);
    
end
%% 26 smarticle rotation
xx=26;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    numSmarts=size(usedSimAm(1).AllSmartPos{1},1);
    frames=size(usedSimAm(1).AllSmartPos{1},2);
    sDat=cell(numSmarts,1);
    sPt=30;
    t=[usedSimAm(1).AllSmartPos{sPt:end,2}];
    for j=2:numSmarts
        qq=cellfun(@(x) x(j,:),usedSimAm(1).AllSmartPos(sPt:end,1),'UniformOutput',0);
        sDat{j}=vertcat(qq{:});
    end
    %     t=ma.tracks{1}(1:end,1)';
    %plot rotation
    yyaxis left
    %     q=zeros(size(sDat{1},1),numSmarts);
    for(j=2:numSmarts)
        angs(:,j)=sDat{j}(:,6);
    end
    %     angs=double(mod(int32(angs*180/pi),180));
    %     angs=sgolayfilt(angs,3,25);
    
    
    mq=mean(angs,2);
    stda=std(angs,0,2);
    plot(t,stda/max(stda));
    %     errorbar(t,mq,stdq)
    % v=sqrt((diff(ma.tracks{1}(1:1:end,2))).^2+(diff(ma.tracks{1}(1:1:end,3))).^2)/.1;
    v= [diff(ma.tracks{1}(sPt:end,2))/diff(ma.tracks{1}(1:2,1)),... %vx
        diff(ma.tracks{1}(sPt:end,3))/diff(ma.tracks{1}(1:2,1))];%vy
    vnew=sqrt(v(:,1).^2+v(:,2).^2);
    ylabel('std(all smarticle rotation angles) normalized');
    
    yyaxis right
    %     plot(t(3:end),a,'k')
    %     plot(t(3:end),vnew/max(vnew),'-','linewidth',2);
    plot(t(2:end),vnew/max(vnew),'-','linewidth',2);
    
    
    %     ylim([-.2,.2])
    ylabel('velocity normalized');
    xlabel('time (s)');
    figText(gcf,16);
    
    S=vnew;
    Fs =10;            % Sampling frequency
    T = 1/Fs;             % Sampling period
    L = length(S);             % Length of signal
    tfft = (0:L-1)*T;        % Time vector
    f = Fs*(0:(L/2))/L;
    Y = fft(S);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    figure(7899);
    plot(f(3:end),P1(3:end))
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    title('Single-sided FFT velocity');
    figText(gcf,18);
    
    
    [Rmm,lags]=xcorr(vnew/max(vnew),stda/max(stda));
    figure(2222);
    plot(lags/length(vnew)*max(tfft),Rmm/max(Rmm));
    title('cross-correlation (vel, std(ang))');
    xlabel('lag(s)');
    figText(gcf,16);
end

%% 27 v dot rhat histogram for dt sampled system
xx=27;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[];ringVel=[];
    bins=200;
    for i=1:L
        dt=diff(usedSimAm(i).fullT(1:2));
        ringPos=usedSimAm(i).fullRingPos;
        %subtract ring position to get vector of deadpos in ring
        deadPos=usedSimAm(i).fullDeadSmartPos-ringPos;
        deadPos=deadPos(1:end-1,:); %remove 1 index to equal vel vec size
        xxx=deadPos(sqrt(sum(deadPos.^2,2))>usedSimAm(i).r/3,:);
        %get velocity
        vRing=diff(ringPos)./dt;
        vRing2=sqrt(sum(vRing.^2,2));
%         %         %%%%
%         vRing=vRing(sqrt(sum(deadPos.^2,2))>usedSimAm(i).r/3,:);
%         deadPos=deadPos(sqrt(sum(deadPos.^2,2))>usedSimAm(i).r/3,:);
%         m=mean(sqrt(sum(vRing.^2,2))); %%%****
%         vRing=vRing(sqrt(sum(vRing.^2,2))>m,:);%%%********
%         deadPos=deadPos(sqrt(sum(vRing.^2,2))>m,:);%%%********
%         %         %%%%
        
        %get vhat and deadposhat
        vHat=vRing./sqrt(sum(vRing.^2,2));
        dpHat=deadPos./sqrt(sum(deadPos.^2,2));
        
        dprod=dot(vHat,dpHat,2);
        %         if(i==15)
        %             continue;
        %         end
        %add cat
        dotProd=[dotProd; dprod];
        ringVel=[ringVel;vRing2];
    end
    cc=histogram(dotProd,bins);
    xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('counts')
    title('dt sampling');
    axis tight
    figText(gcf,18)
    
    
%     bins=201;
    x=linspace(-1,1,bins);
    disc=discretize(dotProd,bins);
    y=zeros(1,bins);
    for(i=1:length(y))
        y(i)=sum(ringVel(disc==i));
    end
    figure(40);
    bb=bar(x,y);
    
%     xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
%     ylabel('counts');
        xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    zlabel('counts scaled by velocity')
    title('dt sampling');
    axis tight
    figText(gcf,18)
end

%% 28 v dot rhat histogram for collision event sampled system
xx=28;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[];
    for i=1:L
        [deadPos,idx]=unique(usedSimAm(i).deadPos,'rows');
        t=deadPos(:,1);
        deadPos=deadPos(:,2:3);
        ringPos=usedSimAm(i).ringpos(idx,:);
        
        %center deadpos by ring pos first position
        deadPos=bsxfun(@minus, deadPos,ringPos(1,:));
        %center ringPos
        ringPos=bsxfun(@minus, ringPos,ringPos(1,:));
        %subtract ring position to get vector of deadpos in ring
        deadPos=deadPos-ringPos;
        deadPos=deadPos(1:end-1,:);%remove 1 index to equal vel vec size
        vRing=diff(ringPos)./diff(t(:,1));
        %get vhat and deadPoshat
        vHat=vRing./sqrt(sum(vRing.^2,2));
        dpHat=deadPos./sqrt(sum(deadPos.^2,2));
        dprod=dot(vHat,dpHat,2);
        dotProd=[dotProd; dprod];
    end
    
    histogram(dotProd);
    
    xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('counts');
    title('collision event sampling');
    axis tight
    figText(gcf,18)
end


%% 29 v dot rhat scatter dt sampling
xx=29;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[];    ringVel=[];
    for i=1:L
        dt=diff(usedSimAm(i).fullT(1:2));
        ringPos=usedSimAm(i).fullRingPos;
        %subtract ring position to get vector of deadpos in ring
        deadPos=usedSimAm(i).fullDeadSmartPos-ringPos;
        deadPos=deadPos(1:end-1,:); %remove 1 index to equal vel vec size
        %get velocity
        vRing=diff(ringPos)./dt;
        
        %get vhat and deadposhat
        vHat=vRing./sqrt(sum(vRing.^2,2));
        dpHat=deadPos./sqrt(sum(deadPos.^2,2));
        dprod=dot(vHat,dpHat,2);
        
        dotProd=[dotProd; dprod];
        ringVel=[ringVel; sqrt(sum(vRing.^2,2))];
    end
    
    scatter(dotProd, ringVel,'.');
    xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('|v_{ring}|');
    title('dt sampling');
    figText(gcf,18)
end


%% 30 v dot rhat scatter collision event sampling
xx=30;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[];   ringVel=[];
    for i=1:L
        [deadPos,idx]=unique(usedSimAm(i).deadPos,'rows');
        t=deadPos(:,1);
        deadPos=deadPos(:,2:3);
        ringPos=usedSimAm(i).ringpos(idx,:);
        
        %center deadpos by ring pos first position
        deadPos=bsxfun(@minus, deadPos,ringPos(1,:));
        %center ringPos
        ringPos=bsxfun(@minus, ringPos,ringPos(1,:));
        %subtract ring position to get vector of deadpos in ring
        deadPos=deadPos-ringPos;
        deadPos=deadPos(1:end-1,:);%remove 1 index to equal vel vec size
        vRing=diff(ringPos)./diff(t(:,1));
        %get vhat and deadPoshat
        vHat=vRing./sqrt(sum(vRing.^2,2));
        dpHat=deadPos./sqrt(sum(deadPos.^2,2));
        dprod=dot(vHat,dpHat,2);
        
        dotProd=[dotProd; dprod];
        ringVel=[ringVel, sqrt(sum(vRing.^2,2))];
    end
    
    scatter(dotProd, ringVel,'.');
    
    xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('|v_{ring}|');
    title('collision sampling');
    figText(gcf,18)
end


%% 31 v dot rhat heatmap dt sampling
xx=31;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[]; ringVel=[];
    for i=1:L
        dt=diff(usedSimAm(i).fullT(1:2));
        ringPos=usedSimAm(i).fullRingPos;
        %subtract ring position to get vector of deadpos in ring
        deadPos=usedSimAm(i).fullDeadSmartPos-ringPos;
        deadPos=deadPos(1:end-1,:); %remove 1 index to equal vel vec size
        %get velocity
        vRing=diff(ringPos)./dt;
        
        %get vhat and deadposhat
        vHat=vRing./sqrt(sum(vRing.^2,2));
        dpHat=deadPos./sqrt(sum(deadPos.^2,2));
        dprod=dot(vHat,dpHat,2);
        
        dotProd=[dotProd, dprod];
        ringVel=[ringVel, sqrt(sum(vRing.^2,2))];
    end
    %add 1 so index is never zero
    dotProd1=dotProd+1;
    binCountVel = 100;
    binCountDot = 100;
    
    maxVel = max(ringVel);
    minVel = min(ringVel);
    binWVel=(maxVel-minVel)/binCountVel;
    binWDot=2/binCountDot;
    mat = zeros(binCountDot, binCountVel);
    for idx = 1:length(ringVel)
        if(isnan(ceil(dotProd1(idx)/binCountDot)))
            continue;
        end
        mat(ceil(dotProd1(idx)/binWDot), ceil(ringVel(idx)/binWVel)) ...
            = mat(ceil(dotProd1(idx)/binWDot), ceil(ringVel(idx)/binWVel)) + 1;
    end
    subplot(1,2,1);
    colormap fire
    imagesc(log(mat'));
    axis tight
    xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('|v_{ring}|');
    
    subplot(1,2,2);
    colormap fire
    surf(log(mat)');
    xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('|v_{ring}|');
    zlabel('log(P(|v|) counts');
    title('dt sampling');
    figText(gcf,18)
    axis tight
    axis square
end

%% 32 v dot rhat heatmap collision event sampling
xx=32;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[];     ringVel=[];
    for i=1:L
        [deadPos,idx]=unique(usedSimAm(i).deadPos,'rows');
        t=deadPos(:,1);
        deadPos=deadPos(:,2:3);
        ringPos=usedSimAm(i).ringpos(idx,:);
        
        %center deadpos by ring pos first position
        deadPos=bsxfun(@minus, deadPos,ringPos(1,:));
        %center ringPos
        ringPos=bsxfun(@minus, ringPos,ringPos(1,:));
        %subtract ring position to get vector of deadpos in ring
        deadPos=deadPos-ringPos;
        deadPos=deadPos(1:end-1,:);%remove 1 index to equal vel vec size
        vRing=diff(ringPos)./diff(t(:,1));
        %get vhat and deadPoshat
        vHat=vRing./sqrt(sum(vRing.^2,2));
        dpHat=deadPos./sqrt(sum(deadPos.^2,2));
        dprod=dot(vHat,dpHat,2);
        
        dotProd=[dotProd; dprod];
        ringVel=[ringVel, sqrt(sum(vRing.^2,2))];
    end
    %add 1 so index is never zero
    dotProd1=dotProd+1;
    binCountVel = 100;
    binCountDot = 100;
    
    maxVel = max(ringVel);
    minVel = min(ringVel);
    binWVel=(maxVel-minVel)/binCountVel;
    binWDot=2/binCountDot;
    mat = zeros(binCountDot, binCountVel);
    for idx = 1:length(ringVel)
        if(isnan(ceil(dotProd1(idx)/binCountDot)))
            continue;
        end
        mat(ceil(dotProd1(idx)/binWDot), ceil(ringVel(idx)/binWVel)) ...
            = mat(ceil(dotProd1(idx)/binWDot), ceil(ringVel(idx)/binWVel)) + 1;
    end
    subplot(1,2,1);
    colormap fire
    imagesc(log(mat'));
    axis tight
    xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('|v_{ring}|');
    
    subplot(1,2,2);
    surf(log(mat)');
    xlabel('$\hat{v}_{ring} \cdot\hat{r}_{inactive}$','interpreter', 'latex');
    ylabel('|v_{ring}|');
    zlabel('log(P(|v|) counts');
    title('collision event sampling');
    figText(gcf,18)
    axis tight
    axis square
    
end
%% 33 log(pdf) vs |vel|
xx=33;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[]; ringVel=[]; y=[];
    binCountVel = 200;
    ti=.005;
    tf=0.4;
    for i=1:L
        dt=diff(usedSimAm(i).fullT(1:2));
        ringPos=usedSimAm(i).fullRingPos;
        %subtract ring position to get vector of deadpos in ring
        %         deadPos=usedSimAm(i).fullDeadSmartPos-ringPos;
        %         deadPos=deadPos(1:end-1,:); %remove 1 index to equal vel vec size
        %get velocity
        vRing=diff(ringPos)./dt;
        
        %get vhat and deadposhat
        vHat=vRing./sqrt(sum(vRing.^2,2));
        %         dpHat=deadPos./sqrt(sum(deadPos.^2,2));
        %         dprod=dot(vHat,dpHat,2);
        %
        %         dotProd=[dotProd, dprod];
        ringVel=[ringVel, sqrt(sum(vRing.^2,2))];
        y = [y;histcounts(sqrt(sum(vRing.^2,2)),binCountVel)];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     maxVel = max(ringVel);
    %     minVel = min(ringVel);
    %     binWVel=(maxVel-minVel)/binCountVel;
    %     y = histcounts(ringVel,binCountVel);
    %     y=y/max(y);
    %
    %     x=binWVel*(1:length(y));
    %     plot(x,y,'-o');
    %     % errorbar(x,y,err);
    %     xlabel('|v|');
    %     ylabel('P(|v|)');
    %     figText(gcf,18);
    %     set(gca,'yscale','log');
    %%%%%%%%%%%%%%%%%%%%%%%%
    maxVel = max(max(ringVel));
    minVel = min(min(ringVel));
    binWVel=(maxVel-minVel)/binCountVel;
    
    
    y=y./sum(y,2);
    err=std(y,0,1);
    ym=mean(y,1)';
    x=binWVel*(1:length(y))';
    %     plot(x,y,'-o');
    %     errorbar(x,ym,err);
    y2=(ym(x>ti&x<tf));
    x2=x(x>ti&x<tf);
    
    
    xlabel('|v| (m/s)');
    ylabel('P(|v|)');
    figText(gcf,18);
    set(gca,'yscale','log');
    ft= fittype('a*x.^2.*exp(-x.^2.*b)',...
        'dependent',{'y'},'independent',{'x'},...
        'coefficients',{'a','b'});
    
    
    f = fit(x2,y2,ft);
    plot(f,x2,y2)
    %     plot(x,ym);
    %  set(gca,'xscale','log');
    % x1=0.005:0.0001:0.2;
    % x2=0.005:.0001:0.6;
    % x1=0.005:0.0001:0.6;
    % x2=0.005:.0001:0.6;
    % x3=0.005:0.0001:0.2;
    % plot(x1,exp(-27*x1)/(1*10^(1.4)),'linewidth',2);
    % plot(x2,x2.^(-2.4)/(1*10^(5.4)),'linewidth',2)
    % plot(x1,
    % semilogy(x2,10*exp(-5*x2)./(x2.^(3))/(1*10^(8)),'.-')
    
end
%% 34 log(pdf) vs |vel| fit line
xx=34;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dotProd=[]; ringVel=[]; y=[];
    binCountVel = 200;
    ti=.005;
    tf=0.1;
    m=zeros(1,length(L));
    for i=1:L
        
        dt=diff(usedSimAm(i).fullT(1:2));
        ringPos=usedSimAm(i).fullRingPos;
        vRing=diff(ringPos)./dt;
        
        vHat=vRing./sqrt(sum(vRing.^2,2));
        ringVel=[sqrt(sum(vRing.^2,2))];
        y = [histcounts(sqrt(sum(vRing.^2,2)),binCountVel)];
        y=y/sum(y);
        
        maxVel = max(max(ringVel));
        minVel = min(min(ringVel));
        binWVel=(maxVel-minVel)/binCountVel;
        x=binWVel*(1:length(y));
        
        y2=log(y(x>ti&x<tf));
        x2=x(x>ti&x<tf);
        p = polyfit(x2,y2,1);
        %         f=fit(x2',y2','poly1');
        %         plot(x,y,'o');
        h1= plot(x2,y2,'o');
        %     h1= plot(f,x2,y2,'o');
        %                hold on;
        plot(x2,x2*p(1)+p(2),'color',h1.Color);
        %         hold off;
        %         pause
        m(i)=p(1);
    end
    
    pts(mean(m),'+-',std(m))
    xlabel('|v| (m/s)');
    ylabel('P(|v|)');
    figText(gcf,18);
    
    
end
%% 45 test
xx=45;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    xlabel('time(s)')
    ylabel('velocity(mm/s)')
    t=ma.tracks{1}(2:10:end,1);
    x=sqrt((diff(ma.tracks{1}(1:1:end,2))/.1).^2+(diff(ma.tracks{1}(1:1:end,3))/.1).^2);
    vv=diff(x)/.1;
    
    v=vv(1:10:end);
    %     v=diff(x(1:10:end));
    vnew=mean(reshape(diff(x)/.1,[10,length(t)]));
    plot(t,vnew*1000);
    plot(t,v*1000)
    plot([0,600],[0,0],'k','linewidth',2);
    
    
end