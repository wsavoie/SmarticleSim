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
%*21. vcorr
%*22. for each msd traj get linear fit of log
%*23. %% 23 histogram of probability(theta_R-theta_r) make vid
%************************************************************
%%
SPACE_UNITS='m';
TIME_UNITS='s';
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
inds=1;
showFigs=[2 22];
useCOM=0;
f=[.2]; rob=[]; v=[];dirs=[];

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
    hold on;
    if(useCOM)
        title('Ring COM');
    else
        title('Ring COG');
    end
    ma.plotTracks
    ma.labelPlotTracks
    text(0,0+.01,'start')
    plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
    
    for i=1:length(ma.tracks)
        plot(ma.tracks{i}(end,2),ma.tracks{i}(end,3),'ko','markersize',4,'MarkerFaceColor','r');
    end
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
    ylabel(cb,'Total Force (N)');
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
    for i=1:length(usedSimAm)
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
    
    leg=[];
    legT={};
    for i=1:length(usedSimAm)
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
    
    for i=1:length(usedSimAm)
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
    for i=1:length(usedSimAm)
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
    for i=1:length(usedSimAm)
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
    for i=1:length(usedSimAm)
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
    for i=1:length(usedSimAm)
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
    NE=length(usedSimAm);
    xp=-1;
    yp=-1;
    dec=500;
    tick=.25;
    leg=[];
    legT={};
    for i=1:length(usedSimAm)
        
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

%% 21
xx=21;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    dt=diff(usedSimAm(1).fullT(1:2));
    dec=1;
    pks=[];
    for i=1:length(usedSimAm)
        %         deadPos=normr(usedSimAm(i).fullDeadSmartPos-repmat(usedSimAm(i).fullRingPos(1,:),[size(usedSimAm(i).fullRingPos,1),1]));
        %         ringPos=usedSimAm(i).fullRingPos-repmat(usedSimAm(i).fullRingPos(1,:),[size(usedSimAm(i).fullRingPos,1),1]);
        deadPos=normr(usedSimAm(i).fullDeadSmartPos(1:end-1,:));
        ringPos=usedSimAm(i).fullRingPos;
        vRing=diff(ringPos)./diff(usedSimAm(i).fullT(1:end,:));
        vRing=normr(vRing);
        deadPos=normr(deadPos);
        
        vTheta=atan2(vRing(:,2),vRing(:,1));
        dTheta=atan2(deadPos(:,2),deadPos(:,1));
        vTheta=decimate(vTheta,dec);
        dTheta=decimate(dTheta,dec);
        [Rmm,lags]=xcorr(vTheta,dTheta,55*dec/dt,'unbiased');
        lags=lags*dt/dec;
        %         Rmm=Rmm(lags<0);
        %         lags=lags(lags<0);
        %         Rmm=Rmm(lags>-55);
        %         lags=lags(lags>-55);
        Rmm=Rmm/max(Rmm);
        subplot(1,2,1);
        hold on; title('correlation');
        xlabel('Lag (s)');ylabel('normalized correlation');
        h=plot(lags,Rmm);
        [pksY,pksX]=findpeaks(Rmm,lags,'MinPeakDistance',10,'minpeakheight',.75,'SortStr','descend');
        %         [pksY,pksX]=findpeaks(abs(Rmm),lags,'MinPeakDistance',10,'SortStr','descend');
        pksY=Rmm(round(lags,5)==round(pksX(1),5));
        
        pks(i)=pksX(1);
        plot(pksX(1),pksY,'^','color',h.Color,'markerfacecolor',h.Color)
        %         pause;
        %         [Rmm]=xcorr2(vRing);
        
    end
    subplot(1,2,2);
    hold on;
    plot(pks,'-o');
    plot(abs(pks),'.-','markersize',15);
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
    p=ma.getMeanMSD([]);
    x=p(:,1);
    y=p(:,2);
    y=y(x>1.2&x<15);
    x=x(x>1.2&x<15);
    lx=log(x);
    ly=log(y);
    pom=fit(lx,ly,'poly1');
    %mean of powers
    clear x y lx ly
    fs=zeros(length(ma.msd),1);
    for i=1:length(ma.msd)
        a=ma.msd{i}(:,1:2);
        x=a(:,1);
        y=a(:,2);
        y=y(x>1.2&x<15);
        x=x(x>1.2&x<15);
        [lx]=log(x);
        [ly]=log(y);
        [f,gof]=fit(lx,ly,'poly1');
        fs(i)=f.p1;
    end
    plot(fs);
    
     msd=ma.getMeanMSD;
    msd=msd(msd(:,1)<15&msd(:,1)>0.1,:);
    [lx]=log(msd(:,1));
    [ly]=log(msd(:,2));
    [f2]=polyfit(lx,ly,1);
    meanPow=f2(1);
    
    pts('(*)mean of powers=',mean(fs),'  power of mean=',pom.p1,' meanpow=',meanPow);
    % std(fs);
end
%% 23 histogram of probability(theta_R-theta_r) make vid
xx=23;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    fps=10;
    writeVid=true;
    if(writeVid)
        outputVideo = VideoWriter(fullfile('','prob(theta).avi'));
        outputVideo.FrameRate=fps;
        open(outputVideo);
    end
    for j=[2:1:60]
        for i=1:length(usedSimAm)
            tt=usedSimAm(i).fullT;
            ringPos(i,:)=normr(usedSimAm(i).fullRingPos(round(tt,5)==j,:));
            deadPos(i,:)=normr(usedSimAm(i).fullDeadSmartPos(1,:));
            %             histogra
        end
        ringTheta=atan2(ringPos(:,2),ringPos(:,1));
        deadTheta=atan2(deadPos(:,2),deadPos(:,1));
        res=mod(ringTheta-deadTheta+pi,2*pi);
        res=res-pi;
        %         res=atan2(ringPos(:,2)-deadPos(:,2),ringPos(:,1)-deadPos(:,1));
        %         figure(j);
        bins=11;
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
        clf;
        
        %        pause(.25);
        
    end
    
    if(writeVid)
           close(outputVideo);  
    end
end