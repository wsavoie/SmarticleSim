directoryout_name = uigetdir('A:\SmarticleRun');
dirsout = dir(directoryout_name);
% directory_name='D:\SimResults\Chrono\SmarticleU\tests\PostProcess';
dirsout = dirsout(~strncmpi('.', {dirsout.name}, 1));
figure(1);
hold on;
gaitRad = zeros(length(dirsout),1);
finalSmarticleAmt = zeros(length(dirsout),2);
gaitLabel=cell(length(dirsout),1);
clogLen = zeros(length(dirsout),2);
clogAmt = zeros(length(dirsout),2);
clogLenTemp = 0;
clogAmtTemp = 0;
partsPerClog=3;
secsPerClog=.5;
clear clogAmt clogLen clogLenTemp clogAmtTemp 
for k=1:length(dirsout)
    directory_name = horzcat(directoryout_name,'\',dirsout(k).name);
    dirs=dir(directory_name);
    dirs = dirs(~strncmpi('.', {dirs.name}, 1));
    x=importdata(horzcat(directory_name,'\',dirs(1).name,'\PostProcess\volumeFraction.txt'));
    dt=x(1,1);
    final=15;
    tt = 0:dt:15;
    data = cell(length(dirs),1);
    smartsOut = zeros(length(tt),length(dirs));
    maxNan=0;
    figure(1);
    
    for i=1:length(dirs)
        p=importdata(horzcat(directory_name,'\',dirs(i).name,'\PostProcess\flowrate.txt'));
        [t,ind] = unique(p(:,1));
        p=p(ind,:);
        data{i}=p;
        
        smartsOut(:,i)=interp1(p(:,1),p(:,2),tt);

        maxNan = find(isnan(smartsOut(:,i)),1,'first');
        smartsOut(maxNan:end,i)=round(smartsOut(maxNan-1,i));
        % plot(p(:,1),p(:,2));
        
        clogs=movsum(diff(data{i}(:,1)),partsPerClog);
        realClogs=(clogs(clogs>secsPerClog));
        if(realClogs)
            clogLenTemp(k,i)=mean(realClogs); 
            clogAmtTemp(k,i)=length(realClogs);
        else
            clogLenTemp(k,i)=0;
            clogAmtTemp(k,i)=0;
        end
    end
    clogLenTemp(isnan(clogLenTemp))=0;
    rsmartsOut=smartsOut;
   
    % tt(maxNan:end)=[];
    meanSmartsOut = mean(smartsOut,2);
    errSmartsOut= std(smartsOut,0,2);
    
    finalSmarticleAmt(k,1)=meanSmartsOut(end);
    finalSmarticleAmt(k,2)=errSmartsOut(end);
    expr='[pi]*([\d.]+)'; 
    [m t] =regexp(dirsout(k).name,expr,'tokens','match');
    gaitRad(k)= str2double(m{1}{1});
%     gaitLabel(k)={['\pi/' num2str(gaitRad(k))]};
    gaitLabel(k)={[num2str(gaitRad(k))]};
    cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7]};
    
    

    %data = [time exited hopper;   total count;    guid]
    % suptitle('Flow of Smarticles in a Hopper')
    % subplot(2,1,1);

    title('Smarticle passed through hopper');
    xlabel('Time(s)');
    ylabel('Smarticles');
    figure(1);
    hold on;
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
figure(2);
hold on;

errorbar(gaitRad,finalSmarticleAmt(:,1),finalSmarticleAmt(:,2),...
    'linewidth',2);
figText(gcf,16);

% set(gca,'xtick',pi./gaitRad,'xticklabel',gaitLabel,'box','on');
set(gca,'xtick',gaitRad,'xticklabel',gaitLabel,'box','on');
xlabel('Gait Radius');
ylabel('Smarticles Through Hopper');
%%
figure(3);
clogLen(:,1)=mean(clogLenTemp,2); clogLen(:,2)=std(clogLenTemp,0,2);
clogAmt(:,1)=mean(clogAmtTemp,2); clogAmt(:,2)=std(clogAmtTemp,0,2);

subplot(1,2,1)
hold on;
errorbar(gaitRad,clogAmt(:,1),clogAmt(:,2),'linewidth',2);
xlabel('Gait Radius');
ylabel('Times Clogged');

set(gca,'box','on');

subplot(1,2,2)
errorbar(gaitRad,clogLen(:,1),clogLen(:,2),'linewidth',2);
xlabel('Gait Radius');
ylabel('Clog Time');
% text(.5,.5,['Clog is <= (' num2str(partsPerClog) ' smarticles)/(' secsPerClog ' s)']);
set(gca,'box','on');
figText(gcf,16);
text(.5,.5,['Clog is <= (' num2str(partsPerClog) ' smarticles)/(' num2str(secsPerClog) ' s)'],'units','normalized');
% ax.XTick= 0:1:max(ax.XTick);

% lineVar=data{1}(:,2);
% yAx=max(lineVar);
% plotNames = {'Stress','Gait','U-Shape','Straight','n-Shape','Vib at \circ','Vib Angle'};
% shapeLines=getShapeLines(data{1}(:,1),data{1}(:,3));
% for i=1:size(shapeLines,1)
%     plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) yAx],'color',shapeLines(i,2:4),'LineWidth',5)
%     text(shapeLines(i,1)+.1,yAx*0.98,plotNames(shapeLines(i,5)))
% end


%% subplot 2
% figure(1);
% subplot(2,1,2)
% hold on;
% title('Flow Rate of Smarticles through hopper');
% % flowRate = diff(data(:,2))./diff(data(:,1));
% flowRate = diff(data(:,2))./diff(data(:,1)); 
% sp2=[];
% xlabel('Time(s)');
% ylabel('Flow Rate d(Smarticles)/{dt}');
% % //ylabel('Smarticle Flow Rate');
% 
% 
% hold on;
% n = 1; % designs an Nth order low pass
% f = 1/dt;
% % f = 10;
% we = 1*2*pi*dt;
% if we >1
%     we = .005; % .005*1/dt/2 5
% end
% [B,A] = butter(n,we);
% % filters the data in vector X with the filter
% %   described by vectors A and B to create the filtered data Y.  The
% %   filter is described by the difference equation:
% %
% %     a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
% %                           - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
% ydat2 = filtfilt(B,A,flowRate);
% % plot(data(1:end-2,1),flowRate(2:end),'LineWidth',2);
% 
% dt = .00025;
% smartFlow = [(0:.00025:15)',zeros(15./.00025+1,1)];
% 
% curr=0;
% timeIdx = 1;
% for i=1:length(smartFlow)
% 
%     if smartFlow(i,1)>=data(timeIdx,1)
%         curr=curr+1;
%         timeIdx=timeIdx+1;
%     end
%     smartFlow(i,2)=curr;
%     
%     
% end
% x1 =(.00025:.00025:15);
% y1=diff(smartFlow(:,2)/.00025);
% 
% % factor(length(y1))
% % fac=max(factor(length(y1)));
% fac = 50;
% % plot(smartFlow(:,1),smartFlow(:,2));
% x2 = mean(reshape(x1,fac,[]));
% y2 =  mean(reshape(y1,fac,[]));
% plot(x2,y2,'LineWidth',2);
% % plot(x2,y2,'.-','markersize',22,'LineWidth',2);
% legend('Smarticle flow rate')
% lineVar=flowRate;
% yAx=max(lineVar);
% shapeLines=getShapeLines(data(:,1),data(:,3));
% plotNames = {'Stress','Gait','U-Shape','Straight','n-Shape','Vib at \circ','Vib Angle'};
% for i=1:size(shapeLines,1)
%     plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) yAx],'color',shapeLines(i,2:4),'LineWidth',5)
%     text(shapeLines(i,1)+.1,yAx*0.98,plotNames(shapeLines(i,5)))
% end
% 
% % axis([data(1,1) data(end,1) min(flowRate) max(flowRate)]);
%  axis([data(1,1) data(end,1) min(y2) max(y2)]);
% ax = gca;
% ax.XTick= 0:1:max(ax.XTick);
% axis auto
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 

%% comparing %

% 
% % dir_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests\PostProcess');
% % x = dir('D:\SimResults\Chrono\SmarticleU\tests\MixedSmartWithOT\r1');
% % c = {x(3:end).name};
% % sp1 = [];
% % for j=1:length(c)
% %     
% %     directory_name=horzcat(dir_name,'\',c{j},'\PostProcess')
% %     data=importdata(horzcat(directory_name,'\flowrate.txt'));
% %     
% %     dt = .0005;
% %     FS = 1/dt;
% %     cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7]};
% %     plotNames = {'Smarticles','Gait','U-Shape','Straight','Tetris','Vib at \circ','Vib Angle'};
% %     plotTypes = [1];
% %     figure(1);
% %     gc =[0; find(diff(data(:,3)))];
% %     gc=gc+1; %time index will be 1 off from diff
% %     val=[data(gc,3)]+1;%added value can be zero (global) and matrices are 1 started
% %     hold on;
% % 
% %     % suptitle('Flow of Smarticles in a Hopper')
% %     a= gca;
% %     % hold on;
% %     title('Smarticle passed through hopper');
% %     xlabel('Time(s)');
% %     ylabel('Smarticles');
% % 
% %    
% %     gg = [gc, val];
% %     ggOld = gg;
% %     gg= sortrows(gg,2); 
% %     colormap hot
% %     sp1=[sp1 plot(data(:,1),data(:,2))]; %starts variable for the legend
% %     gg = [gc, val];
% %     ggOld = gg;
% %     gg= sortrows(gg,2); 
% %     %loop make sure each different line type has proper name in legend, and
% %     %that each configuration only appears once in legend
% %     for i=1:size(gg,1)
% %         if i~=1
% %             if gg(i,2)==gg(i-1,2) %if true will mean a repeat in color on legend
% %                 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)
% %                 continue;
% %             else
% %                 if j==length(c)
% %                     sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
% %                 end
% %             end
% %         else
% %             if j==length(c)
% %                sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
% %             end
% %         end
% %     plotTypes=[plotTypes, (gg(i,2)+1)];
% %     end
% % end
% % legend([sp1(1:end)],horzcat({'0','.25','.5','.75','1'},plotNames(plotTypes)));
% % axis([data(1,1) data(end,1) data(1,2),data(end,2)]);
% % ax = gca;
% % ax.XTick= 0:1:max(ax.XTick);