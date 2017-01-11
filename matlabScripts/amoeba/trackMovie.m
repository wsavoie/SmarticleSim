function [tracks,fps,conv]=trackMovie(worldR,filename,varargin)
figure(1);

% worldR=19.2e-2;%radius of boundary in meters
V = VideoReader(filename);
% im = rgb2gray(readFrame(V));
im = readFrame(V);
imshow(im);
if(length(varargin)>0)
    if(~isempty(varargin{1}))
        r=varargin{1};
    end  
    [x,y]=ginput(2);
    x=abs(x(2)-x(1));
    y=abs(y(2)-y(1));
    r=max(x,y);
end
% r=60;
th=5;
KK=1;

fps=V.framerate;

rr=round([r-th r+th]);
N = V.duration*V.framerate-40;
cents=zeros(N,2);
[centers, radii] = imfindcircles(im,rr,'Sensitivity',.989,'ObjectPolarity','dark');
centers=centers(1,:); radii=radii(1);
viscircles(centers,radii);
if(length(varargin)>1)
    if(~isempty(varargin{2}))
        pp=varargin{2};
    end
end

cents(KK,:)=centers;

while hasFrame(V)
    im = readFrame(V);
    %     imshow(im);
    [centers, radii] = imfindcircles(im,rr,'Sensitivity',.989,'ObjectPolarity','dark');
    
    
    %the strongest circle is far from previous circle
    %this can happen due to multiple circles found
    
    if(isempty(centers))
        error(['lost tracking on frame ',num2str(KK)]);
        
    end
    
    if(length(varargin)>1)
        imshow(im);
    end
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
        %         pause;
    else
        
        centers=centers(1,:); radii=radii(1);
        viscircles(centers,radii,'color','r');
    end
    
    if(length(varargin)>1)
        if(~isempty(varargin{1}))
            pause(pp);
            pts(KK);
            if(KK>590)
                
                pause;
            end
        end
    end
    
    KK=KK+1;
    cents(KK,:)=centers;
    %     if(mod(KK,100)==0)
    %         pts(KK,'/',N);
    %     end
    %     pause(1);
    
    
    
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
tt(:,1)=1/fps:1/fps:KK*1/fps;
tt(:,2:3)=cents*conv;

tracks=cell(1);
tracks{1}= tt;
