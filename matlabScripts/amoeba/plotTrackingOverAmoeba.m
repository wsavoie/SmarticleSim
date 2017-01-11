inFile='vid.avi';
outFile='trackedOut.avi';
iv=VideoReader(inFile);
ov = VideoWriter(outFile);
worldR=19.2e-2;
[tracks,fps,conv]=trackMovie(worldR,inFile,[],[],'nocenter');
% conv
ov.FrameRate=fps;
open(ov);

kk=1;

tracks=tracks{1};
tracks=tracks/conv;
N=length(tracks);
%     kk=1;
while hasFrame(iv)

    h1=figure(1);
    hold off;
    im = readFrame(iv);
    imshow(im);
    hold on;
    %plot beginning point
    plot(tracks(1,2),tracks(1,3),'ro','markersize',8,'MarkerFaceColor','r');
    
    %plot end point
    plot(tracks(1:kk,2),tracks(1:kk,3),'c','linewidth',2);
    plot(tracks(kk,2),tracks(kk,3),'ro','markersize',8,'linewidth',3);
     viscircles([tracks(kk,2) tracks(kk,3)],worldR/conv,'color','r');
    writeVideo(ov,getframe(gcf));   
    kk=kk+1;
    pause(.05);
end
  close(ov)