
%write to video
mainFolder = 'D:\SimResults\Chrono\SmarticleU\tests\BoxAngChangeTorPct30\';
runName = '-50-20160404-024257';
ff = horzcat(mainFolder,runName);
filename=horzcat(ff,'\PostProcess\Stress.txt');

if exist(horzcat(ff,'\PostProcess\stressData.mat'), 'file') == 2 
    clear('simParams','smartPos','frameInfo','filename','file');
    ff = 'D:\SimResults\Chrono\SmarticleU\tests\BoxAngChangeTorPct30\-44-20160404-080211';
    load(horzcat(ff,'\PostProcess\stressData.mat'));
else
    [smartPos, simParams, frameInfo]= readAllSmarticlesAngles(filename,0);
end
%smartPos % angle1  angle2   movetype   zHeight
%simParams  dt      fps      frameInt   buckRad


%% if data in memory only run this sectionw
VID= 0;
pts('Video =',VID);
%colors relating to the moveType and guid
cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7],[.6039,1,0], [0.623, 0 ,1]};



figure(1);
clf;
ann=annotation('textbox', [0.6,0.8,0.1,0.1],...
           'String', {horzcat('Smarticle Amt= 00'),'Current Directed Motion: 00','Time: 00'});

skipFrames = simParams(3);       

%%open video writer
if(VID)
    outputVideo = VideoWriter(fullfile('','ConfigSpaceOut.avi'));
    outputVideo.FrameRate=simParams(2);
    open(outputVideo);
end

lineLen=10;
%deal with how to propagate size to all to minimize line size
rowS = 100;
line = zeros(rowS,5,lineLen);
currFrame = 1;
%row = each smarticle %columns angles %depth= previous angles
for(i=1:size(smartPos,2))
    if(mod(i-1,skipFrames)==0) %-1 to allow first frame to be written
%         moveType = (smartPos{i}(1,3));
        moveType = (frameInfo(i,3));
        cla;
        hold on;
       
        set(gca, 'xlim',[-120,120],'ylim',[-120,120]);
        smartAmt = size(smartPos{i},1);
        
        %maybe output smarticle id to track for sure
%         if(size(line,1)<smartAmt)
%             line = [line; zeros(smartAmt-rowS,2,lineLen)]
%         end
        angs =smartPos{i}(:,1:2);
        c = repmat(cols{moveType+1},[size(angs,1),1]);
        out = [angs c];
        
        newEntry=padarray(out,rowS-smartAmt,'post');
        line = cat(3,newEntry,line);
%         line(:,:,lineLen+1:end)=newEntry;
        
        %
        if(size(line,3)>lineLen+1)
            line(:,:,lineLen+1:end)=[];
        end
        
        for(j=1:smartAmt)
            ang1idxs= find(line(j,1,:)==0);
            zeroedIdxs= find(line(j,2,ang1idxs)==0);
            l = squeeze(line(j,:,:))';
            l(zeroedIdxs,:)=[];
            if(~isempty(l))
%                 a = plot(l(:,1),l(:,2),'o-');
%                 a.MarkerFaceColor = cell2mat(cols(moveType+1));
%                 a.MarkerEdgeColor = 'flat';
%                 a.LineWidth=3;

                a = plot(l(:,1),l(:,2),'-');
                
                a.MarkerFaceColor = cell2mat(cols(moveType+1));
                a.MarkerEdgeColor = 'flat';
                a.LineWidth=3;
                b = scatter(l(:,1),l(:,2),length(l),l(:,3:5),'o','filled','SizeData',50);
                
                
            end
        end
                    
        
        
        ann.String = {horzcat('Smarticle Amt: ',num2str(smartAmt)),...
            horzcat('Directed Motion: ',num2str(moveType))...
            horzcat('Current Time: ',num2str(currFrame*simParams(1)))};
        xlabel('\alpha_2 (\circ)');
        ylabel('\alpha_1 (\circ)');
        title('Smarticle Evolution in Configuration Space');
        axis equal
        set(gca, 'xlim',[-120,120],'ylim',[-120,120]);
        figText(gcf,13);
        set(gcf, 'Position', [2600, 100, 800, 800]);
%         set(gcf, 'Position', [100, 100, 400, 400]);
        if(VID)
            writeVideo(outputVideo,getframe(gcf));
        else
            pause(0.0001);
        end

    end
    currFrame=currFrame+1;
end
if(VID)
    close(outputVideo)
end


