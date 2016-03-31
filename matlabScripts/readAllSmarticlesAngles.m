function [smartPos simParams]=readAllSmarticlesAngles(filename,varargin)
%reads smarticle file into system and saves
%data into cell array
%out{frameNumber} 
%%angle1 angle2, movetype,zHeight
%%%TODO want to determine size of text
%%%to preallocate size of out matrix

i=1;
fid = fopen(filename);
tline = fgets(fid);%first line contains sim info
simParams = sscanf(tline,'%f,')';



tline = fgets(fid)
lineType=0;%-1= ef, %0= volfrac stuff

frame=1;

x=[];
if isempty(varargin{1})
    frameStart=1;
else
    frameStart=varargin{1};
end




while ischar(tline)
%     disp(tline)
    if(strncmpi('#EF',tline,3))
        if(~isempty(x) && frame>frameStart)
%             timestep(:,:,(frame-frameStart))=x;
            smartPos{frame-frameStart}=x;
            x=[];
        end
        x=[];
        i = i+1;
        lineType=-1;
        frame=frame+1;
     
    
    end
    
    d=sscanf(tline,'%f,')';
    if(length(d<5)) %5 is number of params output per smarticle
        x=[x;d];
    end
    
%     if(lineType>0)
%         %angle1 angle2, movetype,zHeight
%        d=sscanf(tline,'%f,')';
%        x=[x;d];
%     end
   
    lineType=lineType+1;
    tline = fgets(fid);

end
