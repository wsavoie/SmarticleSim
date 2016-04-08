function [smartPos simParams frameInfo]=readAllSmarticlesAngles(filename,varargin)
%reads smarticle file into system and saves
%data into cell array
%out{frameNumber} 
%%angle1 angle2, movetype,zHeight
%%%TODO want to determine size of text
%%%to preallocate size of out matrix

%read in all data 
% path = 'D:\SimResults\Chrono\SmarticleU\tests\BoxAngChangeTorPct30v2\';
% a=dir(horzcat(path,'-*'))
% 
% for i=1:length(a)
%     ff= horzcat(path,a(i).name,'\PostProcess\Stress.txt')
%     readAllSmarticlesAngles(ff,0);
% end


fid = fopen(filename','r');     %# Open the file as a binary
lastLine = '';                   %# Initialize to empty
offset = 1;                      %# Offset from the end of file
fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
newChar = fread(fid,1,'*char');  %# Read one character
while (~strcmp(newChar,char(10))) || (offset == 1)
  lastLine = [newChar lastLine];   %# Add the character to a string
  offset = offset+1;
  fseek(fid,-offset,'eof');        %# Seek to the file end, minus the offset
  newChar = fread(fid,1,'*char');  %# Read one character
end
fclose(fid);  %# Close the fil
totFrames = sscanf(lastLine,'#EF%f')+1;
totFrames



i=1;
fid = fopen(filename);
tline = fgets(fid);%first line contains sim info
simParams = sscanf(tline,'%f,')';



tline = fgets(fid);
lineType=0;%-1= ef, %0= volfrac stuff
frameInfo = zeros(totFrames,5);
frame=1;

x=[];
if isempty(varargin{1})
    frameStart=0;
else
    frameStart=varargin{1};
    totFrames= totFrames-frameStart;
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
     
        
        lineType=lineType+1;
        tline = fgets(fid);
        continue;
    end
    d=sscanf(tline,'%f,')';
    if(length(d)==4) %5 is number of params output per smarticle
        x=[x;d];
    else % output will be about each fram
%         tline
         frameInfo(i,:) = d;
    end
    
%     if(lineType>0)
%         %angle1 angle2, movetype,zHeight
%        d=sscanf(tline,'%f,')';
%        x=[x;d];
%     end
    
    lineType=lineType+1;
    tline = fgets(fid);

end
save(horzcat(filename,'/../stressData.mat'),'simParams','smartPos','frameInfo','filename');
fclose(fid);