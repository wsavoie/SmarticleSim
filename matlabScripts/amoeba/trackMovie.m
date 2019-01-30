function [tracks,fps,conv]=trackMovie(worldR,filename,varargin)
%varargin1 = radius
%varargin2 = pause each frame by amount and stop at frame 590
%varargin3 = center tracks at first frame
%varargin4 = dark or bright default bright
figure(1);

% worldR=19.2e-2;%radius of boundary in meters
V = VideoReader(filename);
% im = rgb2gray(readFrame(V));
im = readFrame(V);
imshow(im);
if(length(varargin)>0)
    if(~isempty(varargin{1}))
        r=varargin{1};
    else
    [x,y]=ginput(2);
    x=abs(x(2)-x(1));
    y=abs(y(2)-y(1));
    r=max(x,y);
    r=round(r/2);
    end  

end
r
if(length(varargin)>3&&~isempty(varargin{4}))
polarity=varargin{4};
else
    polarity='dark'; %bright
end
% r=60;
th=5;
KK=1;
sens=.985;
fps=V.framerate;

rr=round([r-th r+th]);
N = round(V.duration*V.framerate);
cents=zeros(N,2);
rads=zeros(N,1);
[centers, radii] = imfindcircles(im,rr,'Sensitivity',sens,'ObjectPolarity',polarity);
centers=centers(1,:); radii=radii(1);
viscircles(centers,radii);
if(length(varargin)>1)
    if(~isempty(varargin{2}))
        pp=varargin{2};
    end
end
closeWaitbar;
cents(KK,:)=centers;
rads(KK,:)=radii;
h = waitbar(0,'Please wait...');
steps = N;
i=0;

while hasFrame(V)
%     waitbar(i/steps,h,{['Processing frame: ',num2str(i),'/',num2str(steps)]});
pts('frame',i,'/',steps);
    im = readFrame(V);
        imshow(im);
    [centers, radii] = imfindcircles(im,rr,'Sensitivity',sens,'ObjectPolarity',polarity);
    
    
    %the strongest circle is far from previous circle
    %this can happen due to multiple circles found
    if(length(varargin)>1)
        imshow(im);
    end
    if(isempty(centers))
%         error(['lost tracking on frame ',num2str(KK)]);
        pts('didn''t track on ',KK);
        centers=cents(KK-1,:);
        radii=rads(KK-1);
        viscircles(centers,radii,'color','g');
    else
    if(norm(cents(KK,:)-centers(1,:))>15)
        %         figure(2);
        %         imshow(im);
        %         viscircles(centers(1,:),radii(1),'color','b');
        
        %replicate matrix into number of centers and subtract to
        %create a vector between two points
        x=centers-repmat(cents(KK,:),[size(centers,1),1]);
        %get length of each vector
        nm=sqrt(sum(abs(x).^2,2));
        %use smallest distance
        [~,idx]=min(nm);
        centers=centers(idx,:); radii=radii(idx);
        
        viscircles(centers,radii,'color','k');
                
    else
        
        centers=centers(1,:); radii=radii(1);
        viscircles(centers,radii,'color','r');
    end
    end
    if(length(varargin)>1)
        if(~isempty(varargin{2}))
%             pause(pp);
            pts(KK);
            if(KK==varargin{2})
                
                pause;
            end
        end
    end
    
    KK=KK+1;
    cents(KK,:)=centers;
    rads(KK,:)=radii;
    %     if(mod(KK,100)==0)
    %         pts(KK,'/',N);
    %     end
    %     pause(1);
    
    pause(0.05);
    i=i+1;
end
tt=ones(size(cents,1),3);


%center at first value
c=1;
%check for nocenter arg
if(length(varargin)>2)
    if varargin{3}=='nocenter'
        c=0;
    end
end
if(c)
    cents=cents-repmat(cents(1,:),[size(cents,1),1]);
end
%get conversion between pix and meters
conv = worldR/r; %meters/pix
% tt(:,1)=1/fps:(1/fps):V.duration;
% tt(:,1)=1/fps:1/fps:KK*1/fps;
tt(:,1)=linspace(0,KK*1/fps,size(tt,1));
tt(:,2:3)=cents*conv;

tracks=cell(1);
tracks{1}= tt;
closeWaitbar;
% close(V);