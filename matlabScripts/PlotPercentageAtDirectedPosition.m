

% mainFolder = 'D:\SimResults\Chrono\SmarticleU\tests\BoxAngChangeTorPct30v2\';
mainFolder='A:\SmarticleRun\temp'
% runName = '-65-20160405-040401';
runName = '-20_0.01_20160614_132810';
matName = '\PostProcess\stressData.mat';
ff = horzcat(mainFolder,runName);
filename=horzcat(ff,'\PostProcess\Stress.txt');

if exist(horzcat(ff,matName), 'file') == 2 
    clear('x','simParams','smartPos','frameInfo','filename','file');
    ff =  horzcat(mainFolder,runName);
    load(horzcat(ff,matName));
else
    [smartPos, simParams, frameInfo]= readAllSmarticlesAngles(filename,0);
end
%smartPos % angle1  angle2   movetype   zHeight
%simParams  dt      fps      frameInt   buckRad

%% if data in memory only run this section
%write to video
%colors relating to the moveType and guid
cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7],[.6039,1,0], [0.623, 0 ,1]};

dirPos = [0 0 -90 -1000 -1000 -1000];
bounds = 2;
%%open video writer
x=zeros(size(smartPos,2),1);
for(i=1:size(smartPos,2))
%     if(mod(i-1,skipFrames)==0) %-1 to allow first frame to be written
        for j=1:size(smartPos{i},1)
            if(all(smartPos{i}(j,1:2)>dirPos(frameInfo(i,3))-bounds)...
                    && all(smartPos{i}(j,1:2)<dirPos(frameInfo(i,3))+bounds))
                x(i)=x(i)+1;
            end
            
        end
        x(i)=x(i)/size(smartPos{i},1);
end
%         h1=histogram(smartPos{i}(:,1),60,'BinWidth',1);
%         set(gca,'xlim',[-120 120],'ylim',[0,12]);
%         smartAmt = size(smartPos{i},1);
%         moveType = (smartPos{i}(1,3));
%         ann.String = {horzcat('Smarticle Amt=',num2str(smartAmt)),...
%             horzcat('Current Directed Motion:',num2str(moveType)),...
%             horzcat('Current Time: ',num2str(i*simParams(1)))};
%         xlabel('Angle (\circ)');
%         ylabel('Smarticle Number');
%         h1.FaceColor = cell2mat(cols(moveType+1));
%         figText(gcf,13);
%         set(gcf, 'Position', [100, 100, 1280, 720]);




 %bottom left corner[x,y],[width to right,height up]
% newAxis=axes('pos',[xpos,ypos), posz(2), posz(2)+sizez(2)]);
shapeLines=getShapeLines(frameInfo(:,1),frameInfo(:,3));
lineVar=x;
yAx=max(lineVar)*1;
plotNames = {'Stress','Gait','U-Shape','Straight','n-Shape','Vib at \circ','Vib Angle'};
figure(2);
hold on;
plot(frameInfo(:,1),x,'LineWidth',2);
title(horzcat('Smarticle % Reaching Directed Position ',177,num2str(bounds),'^\circ'));
xlabel('Time(s)');
ylabel('Smarticle %');

for i=1:size(shapeLines,1)
    plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) yAx],'color',shapeLines(i,2:4),'LineWidth',5)
    text(shapeLines(i,1)+.1,yAx*0.98,plotNames(shapeLines(i,5)))
end
figText(gcf,13)
hold off;