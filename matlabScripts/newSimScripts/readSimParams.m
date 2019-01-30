function [smartVol] = readSimParams(fpath)
%READSIMPARAMS read params from simparams file 

fid = fopen(fpath);

for(i=1:12)
tline = fgets(fid);%first line contains sim info
simParams = textscan(tline,'%s%s%s')';
end
smartVol=str2double(simParams{3});
fclose(fid);
end

