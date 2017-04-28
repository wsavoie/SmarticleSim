clear all;
close all;
% load('D:\ChronoCode\chronoPkgs\Smarticles\matlabScripts\amoeba\smarticleExpVids\rmv3\movieInfo.mat');

fold=uigetdir('D:\ChronoCode\chronoPkgs\Smarticles\matlabScripts\amoeba\smarticleExpVids\');
load(fullfile(fold,'movieInfo.mat'));
figure(1)
SPACE_UNITS = 'm';
TIME_UNITS = 's';
fold

%************************************************************
%* Fig numbers:
%* 1. displacement yvsx
%* 2. MSD
%* 3. vcorr
%* 4. drift velocity
%* 5. plot scatter data in x vs y where z axis is time
%* 6. mean of position at each timestep for all tracks
%* 7. rotation
%* 8. 1 2 3 into single plot
%* 9. for each msd traj get linear fit of log
%*10. partial msd fit with log
%*11. for each msd traj get linear fit of log
%*12. rotate each each track by the rotation of ring
%*13. polar plot of rotation of ring (used for checking axes)
%************************************************************
showFigs=[7];

ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);

%define curve params [] for all
spk=[0]; smart=[]; gait=[1]; rob=[5]; v=[];

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
if(isempty(ma.tracks))
    error('no tracks found for params given!');
end

%% 1 plot Y vs. X
xx=1;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    hax1=gca;
    ma.plotTracks
    ma.labelPlotTracks
    %     text(0,0+.01,'start')
    plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
    y=get(gca,'ylim'); x=get(gca,'xlim');
    c=max(abs(x)); xlim([-c,c]);
    c=max(abs(y)); ylim([-c,c]);
    axis equal
    
    axis([-.25 .25 -.25 .25]);
    x=xlim; y=ylim;
    set(gca,'xtick',[-.2:.1:.2],'ytick',[-.2:.1:.2]);
    y=get(gca,'ylim'); x=get(gca,'xlim');
    plot(x,[0,0],'r');
    plot([0,0],y,'r');
    
    for i=1:length(ma.tracks)
        plot(ma.tracks{i}(end,2),ma.tracks{i}(end,3),'ko','markersize',4,'MarkerFaceColor','r');
        %         leg(i)=h;
        x=['v',num2str(usedMovs(i).pars(5))];
        legT{i}=['v',num2str(usedMovs(i).pars(5))];
    end
    legend(legT);
    legend off;
    %     title('Tracked displacement of smarticle ring COG');
    figText(gcf,14)
    
end
%% 2 plot MSD
xx=2;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    hax2=gca;
    ma = ma.computeMSD;
    ma.plotMeanMSD(gca, true);
    %     ma.plotMSD;
    [fo, gof]=ma.fitMeanMSD;
    % [a b]=ma.fitMeanMSD;
    D=fo.p1/2/ma.n_dim;
    ci = confint(fo);
    %     interval=D-ci(1)/2/ma.n_dim; %interval for D
    %     str = sprintf(['D = %.3e \\pm %.3e\n Goodness of fit: R^2 = %.3f.' ], ...
    %         D, interval, gof.adjrsquare);
    %     text(0.05,0.5,str,'units','normalized');
    % fo.p1/2/obj.n_dim, ci(1)/2/obj.n_dim, ci(2)/2/obj.n_dim, gof.adjrsquare);
    figText(gcf,14)
    axis tight
    xlim([0 15]);
end
%% 3 plot vcorr
xx=3;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    hax3=gca;
    ma = ma.computeVCorr;
    ma.plotMeanVCorr
    m=ma.vcorr{1};
    nansN=sum(isnan(m(:,2)));
    M = mean(m(2:end-nansN-2,2));
    hold on;
    for i= 1:length(ma.vcorr)
        nansN=sum(isnan(ma.vcorr{i}(:,2)));
        M(i)= mean([ma.vcorr{i}(10:end-nansN-1,2)]);
        
    end
    M=mean(M);
    line([ma.vcorr{1}(10,1) ma.vcorr{1}(end,1)], [M M],'color','r','linewidth',3);
    figText(gcf,14);
    text(.2,.7,['mean = ',num2str(M,'%.3f')],'units','normalized','fontsize',16)
    axis tight
end
%% 4 drift velocity
xx=4;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    ma = ma.computeDrift('velocity');
    ma.plotDrift;
    ma.labelPlotTracks;
end
%% 5 plot scatter data in x vs y where z axis is time
xx=5;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    alld=cell2mat(ma.tracks);
    scatter3(alld(:,2),alld(:,3),alld(:,1),'.');
    xlabel('x (m)'); ylabel('y (m)'); zlabel('t(s)');
end
%% 6 mean of position at each timestep for all tracks
xx=6;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
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
    
    %     mmydat/mmxdat
    % plot(mxdat,mydat);
    
    ff=polyfitZero(mxdat,mydat,1);
    hold on;
    plot(mxdat,mydat,'k');
    axis([-.25 .25 -.25 .25]);
    x=xlim; y=ylim;
    set(gca,'xtick',[-.2:.1:.2],'ytick',[-.2:.1:.2]);
    %     y=get(gca,'ylim'); x=get(gca,'xlim');
    xfitdat=linspace(x(1),x(2),5);
    plot(xfitdat,ff(1)*xfitdat);
    
    legend off;
    xlabel('x (m)'); ylabel('y (m)');
    %     title(['mean of position for all runs after ',num2str(ma.tracks{1}(idx,1)),' s'])
    text(.7,.2,['m= ',num2str(ff(1),3)],'units','normalized');
    figText(gcf,14);
    
    plot(x,[0,0],'r');
    plot([0,0],y,'r');
end

%% 7 plot rotation
xx=7;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    for(i=1:length(usedMovs))
        r=usedMovs(i).rot; t=usedMovs(i).t;
        
        pp=[r,t];
        pp(any(isnan(pp),2),:)=[];
        dr=diff(pp(:,1));
        f=find(abs(dr)>1);
        while(~isempty(f)) %remove jumps in data due to relabelling (?)
            r(f)=[];
            t(f)=[];
            pp(f,:)=[];
            dr=diff(pp(:,1));
            f=find(abs(dr)>1);
        end
        
        plot(pp(:,2),pp(:,1));
    end
end
%% 8 put 3 mains into subplot
xx=8;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    
    s1=subplot(2,2,1);
    pos=get(s1,'Position');
    delete(s1);
    hax8=copyobj(hax3,hf10);
    set(hax8, 'Position', pos);
    
    
    s1=subplot(2,2,2);
    pos=get(s1,'Position');
    delete(s1);
    hax8=copyobj(hax1,hf10);
    set(hax8, 'Position', pos);
    
    s1=subplot(2,2,[3,4]);
    pos=get(s1,'Position');
    delete(s1);
    hax8=copyobj(hax2,hf10);
    set(hax8, 'Position', pos);
    %    hax1=gca;
    %    hf2=figure(2);
    %    s1=subplot(211);
    %    pos=get(s1,'Position');
    %    delete(s1);
    %    hax2=copyobj(hax1,hf2);
    %    set(hax2, 'Position', pos);
end
%% 9 for each msd traj get linear fit of log
xx=9;
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
    
    pts('(*)mean of powers=',mean(fs),'  power of mean=',pom.p1);
    % std(fs);
end
%% 10 partial msd fit with log
xx=10;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
    for i=1:length(ma.msd)
        
        msd=ma.msd{i};
        msd=msd(msd(:,1)<15&msd(:,1)>0.1,:);
        
        plot(msd(:,1),msd(:,2))
        [lx]=log(msd(:,1));
        [ly]=log(msd(:,2));
        [f]=polyfit(lx,ly,1);
        fs(i)=f(1);
        [fo, gof] = fit(msd(:,1),msd(:,2), 'poly1', 'Weights', msd(:,3));
        %     fo.p1/4
    end
    ma.labelPlotMSD
    [ma]=ma.fitMSD;
    ma.fitMeanMSD;
    
    msd=ma.getMeanMSD;
    msd=msd(msd(:,1)<15&msd(:,1)>0.1,:);
    [lx]=log(msd(:,1));
    [ly]=log(msd(:,2));
    [f2]=polyfit(lx,ly,1);
    meanPow=f2(1)
    a=ma.getMeanMSD;
    cc=ma.lfit;
    cc.a;
    
    figText(gcf,15)
    %     msdanalyzer
end
%% 11 for each msd traj get linear fit of log
xx=11;
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

%% 12 rotate each each track by the rotation of ring
xx=12;
if(showFigs(showFigs==xx))
    figure(xx)
    hold on;
     figure(xx)
    hold on;
    hax1=gca;
%     ma.plotTracks
    ma.labelPlotTracks
    %     text(0,0+.01,'start')
    plot(0,0,'ro','markersize',8,'MarkerFaceColor','k');
    y=get(gca,'ylim'); x=get(gca,'xlim');
    c=max(abs(x)); xlim([-c,c]);
    c=max(abs(y)); ylim([-c,c]);
    axis equal
    
    axis([-.25 .25 -.25 .25]);
    x=xlim; y=ylim;
    set(gca,'xtick',[-.2:.1:.2],'ytick',[-.2:.1:.2]);
    y=get(gca,'ylim'); x=get(gca,'xlim');
    plot(x,[0,0],'r');
    plot([0,0],y,'r');
    for i=1:length(usedMovs)
        % dpos=diff(pos);
        % R = [cosd(theta(2:end)) -sind(theta(2:end)); sind(theta(2:end)) cosd(theta(2:end))];
        pos=usedMovs(i).data{1}(:,2:3);
        r=usedMovs(i).rot-pi/2;
%         r=r-r(1);
        dpos=[pos(1,:);diff(pos)];
        newpos=dpos;
        for(j=1:size(dpos,1))
            R = [cos(r(j)) -sin(r(j)); sin(r(j)) cos(r(j))];
            newpos(j,:)=R*dpos(j,:)';

        end
        
        newpos=cumsum(newpos);
        plot(newpos(:,1),newpos(:,2));
        plot(newpos(end,1),newpos(end,2),'ko','markersize',4,'MarkerFaceColor','r');
        figText(gcf,14)
    end
    
end

%% 13 polar plot of rotation of ring (used for checking axes)
xx=13;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    for i=1:length(usedMovs)
        r=linspace(0,1,length(usedMovs(i).rot));
        polarplot(usedMovs(i).rot,r);
    end
    figText(gcf,14)
end