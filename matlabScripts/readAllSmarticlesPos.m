function [smartInfo]=readAllSmarticlesPos(filename,varargin)
% function [smartPos simParams frameInfo]=readAllSmarticlesPos(filename,varargin)
% filename='D:\SimResults\Chrono\SmarticleU\tests\PostProcess\Stress.txt';


%read in all data
% path = 'D:\SimResults\Chrono\SmarticleU\tests\BoxAngChangeTorPct30v2\';
% a=dir(horzcat(path,'-*'))
%
% for i=1:length(a)
%     ff= horzcat(path,a(i).name,'\PostProcess\Stress.txt')
%     readAllSmarticlesAngles(ff,0);
% end

%stress_of << mySmarticlesVec[i]->active << ", " << mySmarticlesVec[i]->Get_cm().x << ", " << mySmarticlesVec[i]->Get_cm().y << ", " << mySmarticlesVec[i]->Get_cm().z << std::endl;
fid = fopen(filename,'r');     %# Open the file as a binary
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
% totFrames;
smartInfo=cell(totFrames,2);
t=zeros(totFrames,1);



i=1;
fid = fopen(filename);
tline = fgets(fid);%first line contains sim info
simParams = sscanf(tline,'%f,')';
n=simParams(4); %number of smarticles
% warning('change this to n=simParams(4)');

%%%%%%%%%
% tline = fgets(fid);
% inf=sscanf(tline,'%f,')';
% angleInfo = [inf(2) inf(5)];
%%%%%%%%%%%%%%%

frame=1;
vals=10; %number of values per line for smarticle
x=zeros(n,vals);


while ischar(tline)
    if(i==n+1)
        x(:,4)=-x(:,4);%negative x since left-handed coord
        %         q=quaternion([x(:,7),x(:,8),x(:,9),x(:,10)]);
        %         rz=EulerAngles(q(:),'123');
        %         rot=[pi-rz(3,:)]';
        rot=pi-x(:,8);
        alive=x(:,10);
        %
        angs=x(:,1:2);
        smartInfo{frame,1}=[x(:,4:6),angs,rot,alive];
        % smartInfo{frame}=[x(:,4:6),angs,rot];
        
        
        x=zeros(n,vals);
        %         i = 0;
        i = 1;
        d=sscanf(tline,'%s%f,')';
        smartInfo{frame,2}=d(end);
        tline = fgets(fid);
        
        
        frame=frame+1;
        continue;
    end
    
    d=sscanf(tline,'%f,')';
    if(length(d)>=vals) %4 is number of params output per smarticle
        x(i,:)=d;
        i=i+1;
    end
        tline = fgets(fid);
end
save(horzcat(filename,'\..\PosData.mat'),'smartInfo','t');
fclose(fid);
% toc
