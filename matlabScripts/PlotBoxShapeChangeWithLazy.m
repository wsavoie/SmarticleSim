% function [out]=GetFileInfo(direc)
% direc= 'A:\SmarticleRun\un_eq_no_OT_all_plots_4runs\UnequalOTv2\';
direc='A:\SmarticleRun\TestLazy\';
clf
dirs = dir(direc);
%remove non-folders
dirs = {dirs([dirs.isdir]).name};
% remove . and ..
dirs = dirs(1:end-2);

angles = zeros(length(dirs),1);
idle = 1;
if idle
    idleText=['Ant inspired idle distribution G\approx0.56'];
else
    idleText=['Active Smarticles'];
end
for i=1:length(dirs)
%     if older with no lazy designation in filename
%     [boxAngle]=textscan(dirs{i},'%d%d%d');
    [boxAngle]=textscan(dirs{i},'%d_%7.7f_%d_%d');
    angles(i)=boxAngle{1};
end
idle =boxAngle{2};
mainFolder = direc;
matName = '\PostProcess\stressData.mat';
dirPos = [90 0 -90 -1000 -1000 -1000];
bounds = 2;
uniAngles = unique(angles);
j=1;

% for k=1%=1:length(uniAngles)%
    dataz=0;
    k=find(uniAngles==-40);
    idx=find(angles(:)==uniAngles(k));
    clear('x','simParams','smartPos','frameInfo','filename','file');
    for l=1:length(idx)
        horzcat(mainFolder,dirs{idx(l)},matName)
        load(horzcat(mainFolder,dirs{idx(l)},matName));
        x=zeros(size(smartPos,2),1);
        for i=1:size(smartPos,2)
            for j=1:size(smartPos{i},1)
                if(all(smartPos{i}(j,1:2)>dirPos(frameInfo(i,3))-bounds)...
                        && all(smartPos{i}(j,1:2)<dirPos(frameInfo(i,3))+bounds))
                    x(i)=x(i)+1;%add+1 to smarticles within bounds
                end
            end
            x(i)=x(i)/size(smartPos{i},1);
        end
        if(size(dataz,1)>1)
            dataz(:,l)=x;
        else
            dataz=zeros(size(smartPos,2),length(idx));
            dataz(:,l)=x;
        end
        clear('x');
    end
%     err=mean(dataz,2);
%     err=std(dataz,0,2);
    shadedErrorBar(frameInfo(:,1),mean(dataz,2),std(dataz,1,2),'k');
    title({horzcat('Smarticle Frac. Reaching Directed Position (',177, num2str(bounds),'^\circ)',' Box=',num2str(uniAngles(k)),'^\circ')...
        ,horzcat(idleText)});
    xlabel('Time(s)');
    ylabel('Smarticle Fraction');
    set(gca,'xlim',[0,frameInfo(end,1)],'ylim',[0,1])
% end
% for(i=
% x=zeros(size(smartPos,2),1);
% for(i=1:size(smartPos,2))
% %     if(mod(i-1,skipFrames)==0) %-1 to allow first frame to be written
%         for j=1:size(smartPos{i},1)
%             if(all(smartPos{i}(j,1:2)>dirPos(frameInfo(i,3))-bounds)...
%                     && all(smartPos{i}(j,1:2)<dirPos(frameInfo(i,3))+bounds))
%                 x(i)=x(i)+1;
%             end
%             
%         end
%         x(i)=x(i)/size(smartPos{i},1);
% end
%% 
lineVar=mean(dataz,2);
yAx=1;
plotNames = {'Stress','Gait','U-Shape','Straight','n-Shape','Vib at \circ','Vib Angle'};
shapeLines=getShapeLines(frameInfo(:,1),frameInfo(:,3));
for i=1:size(shapeLines,1)
    plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) yAx],'color',shapeLines(i,2:4),'LineWidth',5)
    text(shapeLines(i,1)+.1,yAx*0.98,plotNames(shapeLines(i,5)))
end