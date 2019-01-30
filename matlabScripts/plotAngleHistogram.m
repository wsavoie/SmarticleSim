

mainFolder = 'D:\SimResults\Chrono\SmarticleU\tests\BoxAngChangeTorPct30v2\';
runName = '-56-20160405-031112';
ff = horzcat(mainFolder,runName);
filename=horzcat(ff,'\PostProcess\Stress.txt');
matName = '\PostProcess\stressData.mat';
if exist(horzcat(ff,matName), 'file') == 2 %file exists
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
VID= 1;
pts('Video =',VID);
%colors relating to the moveType and guid
cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7],[.6039,1,0], [0.623, 0 ,1]};

figure(1);
clf;
% ann=annotation('textbox', [0.6,0.8,0.1,0.1],...
%     'String', {horzcat('Smarticle Amt= 00'),'Current Directed Motion: 00'});

ann=annotation('textbox', [0.45,0.86,0.1,0.1],...
    'String', {horzcat('Smarticle Amt= 00'),'Current Directed Motion: 00'});
skipFrames = simParams(3);

%%open video writer
if(VID)
    %     outputVideo = VideoWriter(fullfile('','histOut.avi'));
    %     outputVideo.FrameRate=simParams(2);
    %     open(outputVideo);
    fold= [uigetdir,'\'];
end

for(i=1:size(smartPos,2))
    if(mod(i-1,skipFrames)==0) %-1 to allow first frame to be written
        %         h1=histogram(smartPos{i}(:,1),60,'BinWidth',1);
        %         set(gca,'xlim',[-120 120],'ylim',[0,12]);
        % h1.FaceColor = cell2mat(cols(moveType+1));
        
        %positions are -180 to 180,
        %I want to have left side of map [90 270]->[-90 90] to be be angle 1
        %
        %and right side of map to be 
        
        [touta1,routa1]=rose(smartPos{i}(:,1),linspace(0,357,120));
        [touta2,routa2]=rose(smartPos{i}(:,2),linspace(0,357,120));
        
        
        smartAmt = size(smartPos{i},1);
        moveType = (frameInfo(i,3));
        ann.String = {horzcat('Smarticle Amt=',num2str(smartAmt)),...
            horzcat('Current Directed Motion:',num2str(moveType)),...
            horzcat('Current Time: ',num2str(i*simParams(1)))};
%         xlabel('Angle (\circ)');
%         ylabel('Smarticle Number');
        
        subplot(1,2,1);
        title('\alpha_1');
        ang1=polarplot(touta1,routa1);
%         ang1.Parent.ThetaTick=[0 30 60 90 120 150 180 -150 -120 -90 -60 -30];
        ang1.Parent.ThetaLim=[-180 180];
        ang1.Parent.RLim=[0 smartAmt];
        subplot(1,2,2)
        title('\alpha_2');
        ang2=polarplot(touta2,routa2);
        ang2.Parent.RLim=[0 smartAmt];
%         ang2.Parent.ThetaTick=[0 30 60 90 120 150 180 -150 -120 -90 -60 -30];
        ang2.Parent.ThetaLim=[-180 180];
        figText(gcf,13);
        %         set(gcf, 'Position',ang2 [100, 100, 1280, 720]);
        if(VID)
            saveas(gcf,[fold, num2str(i)],'png')
            %             writeVideo(outputVideo,getframe(gcf));
        else
            pause(0.001);
        end
        
    end
end
if(VID)
    %     close(outputVideo)
    system(['ffmpeg -framerate 4 -i ',fold,'%d.png -c:v mpeg4 -vtag xvid -r 30 -q 0',fold,'out.avi']);
end
%bottom left corner[x,y],[width to right,height up]
% newAxis=axes('pos',[xpos,ypos), posz(2), posz(2)+sizez(2)]);