function [] = readBallPosDat(fold,startT)
%t,x y z, id, gui, totTorque

ballFiles=dir2(fold,'folders');
rho=7850;
t1=.00127;
t2=.0005;
w=0.0117;

dat=struct;
d=[];
for(i=1:length(ballFiles))
    c=1;
    [parNames,vals]=parseFileNames(ballFiles(i).name);
    versionFolds=dir2(fullfile(fold,ballFiles(i).name),'folders');
    parNames{end+1}='N';
    parNames{end+1}='v';
    vals(end+2)=1;

    for(j=1:length(versionFolds))
        ballOut=[];
        [vNames,vnum]=parseFileNames(versionFolds(j).name);
        vals(end)=vnum;
        %set val for number of smarticles
        vals(end-1)=vals(2)*vals(3);%
        for(k=1:length(parNames))
            eval(['dat(c).',parNames{k},'=',num2str(vals(k)),';']);
        end
        pars(c,:)=vals;
        l=dat(c).lw*w;
        smartSize=[rho,t1,t2,l,w];
        dat(c).fold=fullfile(versionFolds(j).folder,versionFolds(j).name);
        pts(dat(c).fold);
        dat(c).pars=pars(c,:);
        %%%%%%%%%read volume fractiond data%%%%%%%%%%%%%%
        ballChecks=dir2(fullfile(dat(c).fold,'checkPointFiles'));
        for(q=1:length(ballChecks))
            d(q)=str2double((ballChecks(q).name(11:end-4)));
        end
        [inds]=find(d<startT);
        ballChecks(inds)=[];
        d(inds)=[];
        [sc,indz]=sort(d);
        ballChecks=ballChecks(indz);
        t=zeros(1,length(ballChecks));
        for(rr=1:length(ballChecks))
            [~,tz]=parseFileNames(ballChecks(rr).name);
            M = importdata(fullfile(ballChecks(rr).folder,ballChecks(rr).name), ',',10);
            gui(rr)=str2double(M.textdata{9});
            M=M.data;
            %ballOut=[x,y,z,ID,moveType,active,motTorque1,motTorque2]
            
            totTorque(rr)=sum(sum(abs(M(:,36:37))));
            ballOut(:,:,rr)=horzcat(M(:,[1,2,3,20,31,35,36,37]));
            t(rr)=tz;
            smartInfo(:,:,rr)=M(:,[1:7,18,19]); %[x,y,z,e0,e1,e2,e3,ang1,ang2]
        end
        dat(c).ballOut=ballOut;
        dat(c).smartInfo=smartInfo;
        dat(c).t=t;
        dat(c).gui=gui;
        dat(c).totTorque=totTorque;
        dat(c).smartVol=readSimParams(fullfile(dat(c).fold,'PostProcess','simulation_specific_parameters.txt'));
        dat(c).smartSize=smartSize;
%         datafilePath=fullfile(dat(j).fold,'PostProcess','volumeFraction.txt');
%         fdat=importdata(datafilePath);
%         dat(c).gui=M.
%         dat(j).gui=
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %open files
        c=c+1;
        clear ballOut smartInfo t gui totTorque
    end
    %create a struct for params with nice names
    
startTime=(startT)/1000;
save(fullfile(versionFolds(1).folder,'ballDat.mat'),'dat','startTime','smartSize');
c=1;
ballOut=[];
dat=[];
d=[];
end
% ballFiles=dir2(fullfile(fold,'checkPointFiles'));
% %remove all ballfiles < start time
% for(i=1:length(ballFiles))
% c(i)=str2double((ballFiles(i).name(11:end-4)));
% end
% [inds]=find(c<startTime);
% ballFiles(inds)=[];
% c(inds)=[];
% [sc,indz]=sort(c);
% ballFiles=ballFiles(indz);
% t=zeros(1,length(ballFiles));
% for(i=1:length(ballFiles))
%     [~,tz]=parseFileNames(ballFiles(i).name);
%     M = importdata(fullfile(ballFiles(i).folder,ballFiles(i).name), ',',10);
%     M=M.data;
%     %ballOut=[x,y,z,ID,moveType,active]
%     ballOut(:,:,i)=horzcat(M(:,[1,2,3,20,31,35]));
%     t(i)=tz;
% end

