directoryout_name = uigetdir('A:\SmarticleRun');
dirsout = dir(directoryout_name);
% directory_name='D:\SimResults\Chrono\SmarticleU\tests\PostProcess';
dirsout = dirsout(~strncmpi('.', {dirsout.name}, 1));
gaitRad = zeros(length(dirsout),1);
finalSmarticleAmt = zeros(length(dirsout),2);
gaitLabel=cell(length(dirsout),1);
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
for k=1:length(dirsout)
    directory_name = horzcat(directoryout_name,'\',dirsout(k).name);
    dirs=dir(directory_name);
    dirs = dirs(~strncmpi('.', {dirs.name}, 1));
    smartsOut = zeros(length(tt),length(dirs));
    maxNan=0;
    figure(1);
    
    for i=1:length(dirs)
        volFracData=importdata(horzcat(directory_name,'\',dirs(i).name,'\PostProcess\volumeFraction.txt'));
%         volFracData=str2double(volFracData.textdata);
%         vv=str2num(strvcat(volFracData.data,strrep(volFracData.textdata,'-1.#IND','0')))
        p=importdata(horzcat(directory_name,'\',dirs(i).name,'\PostProcess\flowrate.txt'));
        [t,ind] = unique(p(:,1));
        p=p(ind,:);
        pp = p(p(:,1)>=1.5,:);
        pp(:,2)=pp(:,2)-pp(1,2);
        smartsOut(:,i)=interp1(pp(:,1),pp(:,2),tt,'next');
        firstNonNan=find(smartsOut(:,i)>0,1,'first');
        if firstNonNan>1
            smartsOut(1:firstNonNan-1,i)=0;
        end
        newFirstNonNan=find(isnan(smartsOut(:,i)),1,'first');
        smartsOut(newFirstNonNan:end,i)=round(smartsOut(newFirstNonNan-1,i));
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
%     clf;
%     histogram(histClog,20);
%     histClog=0;
    clogLenTemp(isnan(clogLenTemp))=0;
    rsmartsOut=smartsOut;

    % tt(maxNan:end)=[];
    meanSmartsOut = mean(smartsOut,2);
    errSmartsOut= std(smartsOut,0,2);
    
    finalSmarticleAmt(k,1)=meanSmartsOut(end);
    finalSmarticleAmt(k,2)=errSmartsOut(end);
    expr='[pi]*([\d.]+)'; 
    [m, t] =regexp(dirsout(k).name,expr,'tokens','match');
    gaitRad(k)= str2double(m{1}{1});
%     gaitLabel(k)={['\pi/' num2str(gaitRad(k))]};
    gaitLabel(k)={[num2str(gaitRad(k))]};
    cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7]};
    
    

    %data = [time exited hopper;   total count;    guid]
    % suptitle('Flow of Smarticles in a Hopper')
    % subplot(2,1,1);

    title('comp. CDF of smarticles clogs');
    xlabel('Time(s)');
    ylabel('P(t)');
    figure(1);
    hold on;
    [f,x]=ecdf(clog2{k});
    f=1-f;
    plot(x,f,'linewidth',2);
    set(gca,'xscale','log','yscale','log')
    loglinefit(10,x,f);
    % plot(data{1}(:,1),data{1}(:,2),'linewidth',2); %starts variable for the legend
    % 
    % plot(tt,smartsOut,'linewidth',2);
    col =get(gca,'colororder');
    colidx = get(gca,'colororderindex');
%     shadedErrorBar(tt,meanSmartsOut,errSmartsOut,{'LineWidth',2,'color',col(colidx,:),},1);

    %loop make sure each different line type has proper name in legend, and
    %that each configuration only appears once in legend
    % leg = legend('Smarticles');
    pts(gaitRad(k));
end;

set(gca,'xlim',[0 final]);
[gaitRad, sortidxs]=sortrows(gaitRad,1);
finalSmarticleAmt = finalSmarticleAmt(sortidxs,:);
gaitLabel=gaitLabel(sortidxs);
% totalEnergy=totalEnergy(sortidxs);
figure(2);
hold on;

errorbar(gaitRad,finalSmarticleAmt(:,1),finalSmarticleAmt(:,2),...
    'linewidth',2);
figText(gcf,16);

% set(gca,'xtick',pi./gaitRad,'xticklabel',gaitLabel,'box','on');
set(gca,'xtick',gaitRad,'xticklabel',gaitLabel);
xlabel('Gait Radius (rads)');
ylabel('Smarticles Through Hopper');
%%
figure(3);
clogLen(:,1)=mean(clogLenTemp,2); clogLen(:,2)=std(clogLenTemp,0,2);
clogAmt(:,1)=mean(clogAmtTemp,2); clogAmt(:,2)=std(clogAmtTemp,0,2);

subplot(1,2,1)
hold on;
errorbar(gaitRad,clogAmt(:,1),clogAmt(:,2),'linewidth',2);
xlabel('Gait Radius (rads)');
ylabel('Avg. Times Clogged');
figText(gcf,16);
subplot(1,2,2)
errorbar(gaitRad,clogLen(:,1),clogLen(:,2),'linewidth',2);
xlabel('Gait Radius (rads)');
ylabel('Avg. Clog Time (s)');
% text(.5,.5,['Clog is <= (' num2str(partsPerClog) ' smarticles)/(' secsPerClog ' s)']);

figText(gcf,16);
text(.1 ,.2,['Clog is <= (' num2str(partsPerClog) ' smarticles)/(' num2str(secsPerClog) ' s)'],'units','normalized');
% ax.XTick= 0:1:max(ax.XTick);

figure(4)
hold on;
xlabel('Gait Radius (rads)');
ylabel('\Sigma Torque (Nm)');
mTe = mean(totalEnergy,2);
mTe= mTe(sortidxs);
TeErr = std(totalEnergy,0,2);
TeErr = TeErr(sortidxs);
errorbar(gaitRad, mean(totalEnergy,2),std(totalEnergy,0,2),'linewidth',2);
figText(gcf,16);



% lineVar=data{1}(:,2);
% yAx=max(lineVar);
% plotNames = {'Stress','Gait','U-Shape','Straight','n-Shape','Vib at \circ','Vib Angle'};
% shapeLines=getShapeLines(data{1}(:,1),data{1}(:,3));
% for i=1:size(shapeLines,1)
%     plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) yAx],'color',shapeLines(i,2:4),'LineWidth',5)
%     text(shapeLines(i,1)+.1,yAx*0.98,plotNames(shapeLines(i,5)))
% end


