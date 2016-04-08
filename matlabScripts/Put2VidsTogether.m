clf
clc
clear all



source1='C:/Users/root/Desktop/60SecsSmarts.avi';
source2='histOut.avi';
source3='ConfigSpaceOut.avi';
% vid1 = vision.VideoFileReader(source1);
% vid2 = vision.VideoFileReader(source2);


outputVideo = VideoWriter(fullfile('','tripleVid.avi'));
outputVideo.FrameRate=30;
open(outputVideo);
    
    
v1 = VideoReader(source1);
v2 = VideoReader(source2);
v3 = VideoReader(source3);
f1=figure(1);

% while ~isDone(vid1)
% 
% %    frame1 = step(vid1);
% %    frame2 = step(vid2);
% %    
% %    frame = horzcat(frame1, frame2);
% %     
% %    step(vidP,frame);
% 
% 
% end
    axis
    axis1=f1.Children;
    axis2=axes('pos',[.1, .1, .3,.3]);
    axis3=axes('pos',[.6, .1, .4,.3]);
while hasFrame(v1)
    h1=figure(1);

    video1 = readFrame(v1);
    video2 = readFrame(v2);
    video3 = readFrame(v3);
    
    h1.CurrentAxes= h1.Children(3);
    imagesc(flip(video1,1));
    hold on;
%         set(gca,'xdir','reverse');
    h1.CurrentAxes =h1.Children(2);
    a =imshow(video2);
%     
%     x= 3*1800/4;
%     a.XData=[x-540 x];
%     a.YData=[720-360 720];
    
     h1.CurrentAxes=h1.Children(1);
    c =imshow(video3);
%     d= 3*1800/4;
%     c.XData=[1 400];
%     c.YData=[720-400 720];
    hold off;
    
 writeVideo(outputVideo,getframe(gcf));   
end

  close(outputVideo)


% source1='histout.avi';
% source2='outVid.avi';
% v1=VideoReader(source1);
% v2=VideoReader(source2);
% 
% % lastFrame = read(v1, inf);
% % numFrames = v1.NumberOfFrames;
% 
% 
% frames1=floor(v1.Duration*v1.FrameRate);
% frames2=floor(v2.Duration*v2.FrameRate);
% figure(1);
% 
% for f1 = 1 : min(frames1, frames2)
%     subplot(2,1,1)
%     thisframe1 = readFrame(v1,'native');
%     imshow(thisframe1);
% %     drawnow();
%     subplot(2,1,2);
%     thisframe2 = readFrame(v2);
%     imshow(thisframe2);
%     drawnow();
%     pause(1);
% end