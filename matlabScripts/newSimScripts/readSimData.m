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
for(i=1:length(vibFolds))
    [parNames,vals]=parseFileNames(vibFolds(i).name);
    versionFolds=dir2(fullfile(fold,vibFolds(i).name));
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
        
        %vol fraction file format:
        %time,smartcount,phi,zmax,zcomz,totTorque,guiVal
        dat(c).t=fdat(:,1);
        dat(c).phi=fdat(:,3);
        dat(c).gui=fdat(:,7);
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %open files
        c=c+1;
    end
    %create a struct for params with nice names

    
end

save(fullfile(fold,'simData.mat'),'dat','pars');
end

