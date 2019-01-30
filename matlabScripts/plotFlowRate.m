directoryout_name = uigetdir('A:\SmarticleRun');
dirsout = dir(directoryout_name);
% directory_name='D:\SimResults\Chrono\SmarticleU\tests\PostProcess';
dirsout = dirsout(strncmpi('hopper', {dirsout.name}, 1));
if isempty(dirsout)
    dirsout = dir(directoryout_name);2
   dirsout = dirsout(strncmpi('Switch', {dirsout.name}, 1)); 
end
gaitRad = zeros(length(dirsout),1);
finalSmarticleAmt = zeros(length(dirsout),2);
clogLen = zeros(length(dirsout),2);
clogAmt = zeros(length(dirsout),2);
clogLenTemp = 0;
clogAmtTemp = 0;
totalEnergy = 0;
partsPerClog=3;
secsPerClog=0;

%%%get dt in a dirty way
directory_name = horzcat(directoryout_name,'\',dirsout(1).name);
dirs=dir(directory_name);
dirs = dirs(~strncmpi('.', {dirs.name}, 1));
x=importdata(horzcat(directory_name,'\',dirs(1).name,'\PostProcess\volumeFraction.txt'));
dt=x(1,1);
final=15;
tt = 1.5:dt:15;
clear clog2;
clog2= cell(length(dirsout),1);
histClog = [];

%TYPE=1:
%TYPE=2: vary OT switch speed
%TYPE=3: vary OT threshold
%TYPE=4: vary OT thresh and switch speed


for k=1:length(dirsout)
    directory_name = horzcat(directoryout_name,'\',dirsout(k).name);
    dirs=dir(directory_name);
    dirs = dirs(~strncmpi('.', {dirs.name}, 1));
    smartsOut = zeros(length(tt),length(dirs));
    maxNan=0;
    
    for i=1:length(dirs)
        volFracData=importdata(horzcat(directory_name,'\',dirs(i).name,'\PostProcess\volumeFraction.txt'));
        fname = horzcat(directory_name,'\',dirs(i).name,'\PostProcess\flowrate.txt');
        if exist(fname, 'file') == 2
            p=importdata(fname);
        else
            p=[0 0 0];
        end
        p(end+1,:) =[final p(end,2), p(end,3)]; %add last to flow rate file which has last time
        [~,ind] = unique(p(:,1));
        %remove only accept unique times
        p=p(ind,:);
        %only accept flow after 1.5 secs, first few through don't count
        pp = p(p(:,1)>=1,:);
        %reset exited to account for particles removed at beginning
        pp(:,2)=pp(:,2)-pp(1,2);
        
        %interpolate smarticles out with timestep
        if size(pp,1)<2
            smartsOut(:,i)=zeros(length(tt),1);
        else
            smartsOut(:,i)=interp1(pp(:,1),pp(:,2),tt,'next');
        end
        %nans will arise, set to zero
        firstNonNan=find(smartsOut(:,i)>0,1,'first');
        if firstNonNan>1
            smartsOut(1:firstNonNan-1,i)=0;
        end
        newFirstNonNan=find(isnan(smartsOut(:,i)),1,'first');
        
        smartsOut(newFirstNonNan:end,i)=round(smartsOut(max(newFirstNonNan-1,1),i));
        %first indices may be NaN
        clogs=movsum(diff(pp(:,1)),partsPerClog);
        clogs2=diff(pp(:,1));
        clog2(k)={[clog2{k}; clogs2(clogs2>secsPerClog)]};
        realClogs=(clogs(clogs>secsPerClog));
        %         histClog = [histClog; realClogs];
        if(realClogs)
            clogLenTemp(k,i)=mean(realClogs);
            clogAmtTemp(k,i)=length(realClogs);
        else
            clogLenTemp(k,i)=0;
            clogAmtTemp(k,i)=0;
        end
        volFracData(isnan(volFracData))=0;
        totalEnergy(k,i) = sum(volFracData(:,6));
        
    end
    clogLenTemp(isnan(clogLenTemp))=0;
    meanSmartsOut = mean(smartsOut,2);
    errSmartsOut= std(smartsOut,0,2);
    
    finalSmarticleAmt(k,1)=meanSmartsOut(end);
    finalSmarticleAmt(k,2)=errSmartsOut(end);

%         expr='[pi]*([\d.]+)';
%         [m, t] =regexp(dirsout(k).name,expr,'tokens','match');
%         gaitRad(k)= str2double(m{1}{1});
end;
TYPE=3;

figure(1);
hold on;
for ii = 1:length(dirsout)
    switch TYPE
        case 1% vary gait radius
            expr='[pi]*([\d.]+)';
            [m, t] =regexp(dirsout(ii).name,expr,'tokens','match');
            gaitRad(ii)= str2double(m{1}{1});
            xlab = 'Gait Radius (rads)';
            xtickz = num2str(gaitRad,'%1.1f');
        case 2% vary OT switch speed
            expr='[Switch]*([\d.]+)';
            [m, t] =regexp(dirsout(ii).name,expr,'tokens','match');
            gaitRad(ii)= str2double(m{1}{1});
            xlab= 'OT switch speed (s)';
            xtickz = num2str(gaitRad,'%1.2f');
        case 3% vary OT threshold
            %%$example "Switch=0.25_t=.04_tt=$V_r=1.0"
            expr='[tt]*([.\d]+)';
            [m, t] =regexp(dirsout(ii).name,expr,'tokens','match');
            gaitRad(ii)= str2double(m{3}{1});
            xlab= 'OT threshold';
            xtickz = num2str(gaitRad,'%1.2f');
        case 4% vary OT thresh and switch speed
    end
    pts(gaitRad);
    
    [f,x]=ecdf(clog2{ii});
    f=1-f;
    h=plot(x,f,'linewidth',1.5);
    [xfit,yfit] = loglinefit(10,x,f);
    h2=loglog(xfit,yfit,'--','linewidth',2,'color',h.Color);
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end
title('comp. CDF of smarticles clogs');
xlabel('Time(s)');
ylabel('P(t)');
set(gca,'xscale','log','yscale','log')

legend(num2str(gaitRad));
set(gca,'xlim',[0 final]);
[gaitRad, sortidxs]=sortrows(gaitRad,1);
finalSmarticleAmt = finalSmarticleAmt(sortidxs,:);


figure(2);
hold on;

errorbar(gaitRad,finalSmarticleAmt(:,1),finalSmarticleAmt(:,2),'linewidth',2);
set(gca,'xtick',gaitRad,'xticklabel',num2str(gaitRad,'%1.1f'),'box','on');
xlabel(xlab);
ylabel('Smarticles Through Hopper');
figText(gcf,16);
%%
figure(3);
clogLen(:,1)=mean(clogLenTemp,2); clogLen(:,2)=std(clogLenTemp,0,2);
clogAmt(:,1)=mean(clogAmtTemp,2); clogAmt(:,2)=std(clogAmtTemp,0,2);

% subplot(1,2,1)
hold on;
errorbar(gaitRad,clogAmt(:,1),clogAmt(:,2),'linewidth',2);
xlabel(xlab);
ylabel('Avg. Times Clogged');
figText(gcf,16);

% subplot(1,2,2)
% errorbar(gaitRad,clogLen(:,1),clogLen(:,2),'linewidth',2);
% xlabel('Gait Radius (rads)');
% ylabel('Avg. Clog Time (s)');
% % text(.5,.5,['Clog is <= (' num2str(partsPerClog) ' smarticles)/(' secsPerClog ' s)']);
%
% figText(gcf,16);
% text(.1 ,.2,['Clog is <= (' num2str(partsPerClog) ' smarticles)/(' num2str(secsPerClog) ' s)'],'units','normalized');
% % ax.XTick= 0:1:max(ax.XTick);

figure(4)
hold on;
xlabel(xlab);
ylabel('\Sigma Torque (Nm)');
mTe = mean(totalEnergy,2);
mTe= mTe(sortidxs);
TeErr = std(totalEnergy,0,2);
TeErr = TeErr(sortidxs);
errorbar(gaitRad, mean(totalEnergy,2),std(totalEnergy,0,2),'linewidth',2);
figText(gcf,16);
clear all
