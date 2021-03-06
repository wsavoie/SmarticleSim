function [] = readSimData(fold)
%READSIMDATA reads through files from chrono sim and compiles and saves it
%into a mat file
%  fold is vibVary folder and will contain final mat file
%  folder heirarchy is as follows:
%  vibVary->lw_#1_nl_#2_npl_#3_vib_#4->v_#5->PostProcess->volumeFraction
%  pars= [lw,nl,npl,vib,N,v]

vibFolds=dir2(fold,'folders');
c=1;
dat=struct;
rho=7850;
t1=.00127;
t2=.0005;
w=0.0117;
for(i=1:length(vibFolds))
    [parNames,vals]=parseFileNames(vibFolds(i).name);
    versionFolds=dir2(fullfile(fold,vibFolds(i).name),'folders');
    parNames{end+1}='N';
    parNames{end+1}='v';
    vals(end+2)=1;
    for(j=1:length(versionFolds))
        [vNames,vnum]=parseFileNames(versionFolds(j).name);
        vals(end)=vnum;
        %set val for number of smarticles
        vals(end-1)=vals(2)*vals(3);%
        for(k=1:length(parNames))
            eval(['dat(c).',parNames{k},'=',num2str(vals(k)),';']);
        end
        pars(c,:)=vals;
        
        dat(c).fold=fullfile(versionFolds(j).folder,versionFolds(j).name);
        dat(c).pars=pars(c,:);
        %%%%%%%%%read volume fractiond data%%%%%%%%%%%%%%
        datafilePath=fullfile(dat(c).fold,'PostProcess','volumeFraction.txt');
        fdat=importdata(datafilePath);
        l=dat(c).lw*w;
        smartSize=[rho,t1,t2,l,w];
        %vol fraction file format:
        %time,smartcount,phi,zmax,zcomz,totTorque,guiVal
        dat(c).t=fdat(:,1);
        dat(c).phi=fdat(:,3);
        dat(c).gui=fdat(:,7);
        dat(c).smartSize=smartSize;
        
        %added this param on 1/11/19
        if size(fdat,2)>7
            dat(c).bucketExist=fdat(:,8);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%read stress stick data%%%%%%%%%%%%%%
        stressDataFile=fullfile(dat(c).fold,'PostProcess','stressHook.txt');
        if(exist(stressDataFile,'file'))
            fdat=importdata(stressDataFile);
            D=fdat.data;
            td=fdat.textdata;
            clear fdat;
            
            %stresshook data fraction file format:
            %tstep, reactForceZ, prismaticState, cylinderWallForceR, cylinderWallForceZ
            dat(c).hookF=D(:,2);
            dat(c).prismSt=D(:,3);
            dat(c).wallF=D(:,4:5);
            
            if size(D,2)>5
                dat(c).hookPos=D(:,6);
            end
            clear D;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %open files
        c=c+1;
    end
    %create a struct for params with nice names
    
    
end

save(fullfile(fold,'vibData.mat'),'dat');
end

