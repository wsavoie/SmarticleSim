function [] = saveFolderN(fold)
%SAVEFOLDERN Summary of this function goes here
%   Detailed explanation goes here
if(~exist('createPlane','file'))
    setupGeom3d
end
c=dir2(fold,'folders');
nonEnt=0;
for(q=1:length(c))
    ff=fullfile(fold,c(q).name,'Nout.mat');
    wid='MATLAB:inpolygon:ModelingWorldLower';
    warning('off',wid);
    load(fullfile(fold,c(q).name,'ballDat.mat'));
    for(i=1:length(dat))
        tic
        Nout(:,:,i)=generatePackingFromSimDatPLANE_mex(dat(i),nonEnt);
        toc
    end
    usedNdat=dat;
    save(ff,'Nout','usedNdat','nonEnt');
    clear Nout usedNdat dat
end

end