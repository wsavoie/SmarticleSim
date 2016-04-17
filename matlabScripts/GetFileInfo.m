% function [out]=GetFileInfo(direc)

dirs = dir(direc);
%remove non-folders
dirs = {dirs([dirs.isdir]).name};
% remove . and ..
dirs = dirs(1:end-2);

angles = zeros(length(dirs),1);
for i=1:length(dirs)
    [boxAngle]=textscan(dirs{i},'%d%d%d');
    angles(i)=boxAngle{1};
    
end
mainFolder = 'D:\SimResults\Chrono\SmarticleU\tests\BoxAngChangeTorPct30v2\';
matName = '\PostProcess\stressData.mat';
dirPos = [0 0 -90 -1000 -1000 -1000];
bounds = 2;
uniAngles = unique(angles);
j=1;

for k=1%=1:length(uniAngles)%
    dataz=0;
    idx=find(angles(:)==uniAngles(k));
    clear('x','simParams','smartPos','frameInfo','filename','file');
    for(l=1:length(idx))
        horzcat(mainFolder,dirs{idx(l)},matName)
        load(horzcat(mainFolder,dirs{idx(l)},matName));
        x=zeros(size(smartPos,2),1);
        for i=1:size(smartPos,2)
            for j=1:size(smartPos{i},1)
                if(all(smartPos{i}(j,1:2)>dirPos(frameInfo(i,3))-bounds)...
                        && all(smartPos{i}(j,1:2)<dirPos(frameInfo(i,3))+bounds))
                    x(i)=x(i)+1;
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
    shadedErrorBar(frameInfo(:,1),mean(dataz,2),std(dataz,1,2),'k')
    title(horzcat('Smarticle % Reaching Directed Position Box=',num2str(uniAngles(k)),177,num2str(bounds),'^\circ'));
    xlabel('Time(s)');
    ylabel('Smarticle %');
    set(gca,'xlim',[0,30],'ylim',[0,1])
end
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