clc
close all
clear all



source1='histout.avi';
source2='outVid.avi';
vid1 = vision.VideoFileReader(source1);
vid2 = vision.VideoFileReader(source2);
vidP = vision.VideoPlayer;

while ~isDone(vid1)
   frame1 = step(vid1);
   frame2 = step(vid2);

   frame = horzcat(frame1, frame2);

   step(vidP,frame);
end

release(vid1);
release(vid2);
release(vidP);




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