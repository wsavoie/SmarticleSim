

mainFolder = 'D:\SimResults\Chrono\SmarticleU\tests\BoxAngChangeTorPct30v2\';
runName = '-65-20160405-040401';
ff = horzcat(mainFolder,runName);
filename=horzcat(ff,'\PostProcess\Stress.txt');
matName = '\PostProcess\stressData.mat';
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
VID= 0;
pts('Video =',VID);
%colors relating to the moveType and guid
cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7],[.6039,1,0], [0.623, 0 ,1]};

figure(1);
clf;
ann=annotation('textbox', [0.6,0.8,0.1,0.1],...
           'String', {horzcat('Smarticle Amt= 00'),'Current Directed Motion: 00'});

skipFrames = simParams(3);       

%%open video writer
if(VID)
    outputVideo = VideoWriter(fullfile('','histOut.avi'));
    outputVideo.FrameRate=simParams(2);
    open(outputVideo);
end
for(i=1:size(smartPos,2))
    if(mod(i-1,skipFrames)==0) %-1 to allow first frame to be written
        h1=histogram(smartPos{i}(:,1),60,'BinWidth',1);
        set(gca,'xlim',[-120 120],'ylim',[0,12]);
        smartAmt = size(smartPos{i},1);
        moveType = (smartPos{i}(1,3));
        ann.String = {horzcat('Smarticle Amt=',num2str(smartAmt)),...
            horzcat('Current Directed Motion:',num2str(moveType)),...
            horzcat('Current Time: ',num2str(i*simParams(1)))};
        xlabel('Angle (\circ)');
        ylabel('Smarticle Number');
        h1.FaceColor = cell2mat(cols(moveType+1));
        figText(gcf,13);
        set(gcf, 'Position', [100, 100, 1280, 720]);
        if(VID)
            writeVideo(outputVideo,getframe(gcf));
        else
            pause(0.001);
        end

    end
end
if(VID)
    close(outputVideo)
end
 %bottom left corner[x,y],[width to right,height up]
% newAxis=axes('pos',[xpos,ypos), posz(2), posz(2)+sizez(2)]);