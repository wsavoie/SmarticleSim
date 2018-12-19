function [] = saveFolderN(fold)
%SAVEFOLDERN Summary of this function goes here
%   Detailed explanation goes here
if(~exist('createPlane','file'))
    setupGeom3d
end
c=dir2(fold,'folders');

for(q=1:length(c))
    ff=fullfile(fold,c(q).name,'Nout.mat');
    wid='MATLAB:inpolygon:ModelingWorldLower';
    warning('off',wid);
    load(fullfile(fold,c(q).name,'ballDat.mat'));
    for(i=1:length(dat))
        tic
        Nout(:,:,i)=generatePackingFromSimDatPLANE_mex(dat(i),1);
        toc
    end
    save(ff,'Nout','dat');
    clear Nout
end

end