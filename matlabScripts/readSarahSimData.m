clear all
tracks=9;
N=127;
for(j=1:tracks)
[filez fold]=uigetfile('*.txt','simData','A:\SmarticleRun\SarahAmeobotData\127');
pathz=fullfile(fold,filez);
simD=importdata(pathz);
% remove first line
simD(1,:)=[];
dat=zeros(simD(end,1)-simD(1,1),3);
dat(1,1)=simD(1,1);
dat(2:end,1)=[2:(size(dat,1))]+dat(1,1)-1;

for(i=1:size(simD,1)-1)
frames= simD(i+1,1)-simD(i,1);
idx=simD(i,1)-simD(1,1)+1;
dat(idx:idx+frames,2)=simD(i,2);
dat(idx:idx+frames,3)=simD(i,3);
end

simD=downsample(dat,100);
% blah=simD;
simD={bsxfun(@minus, simD,simD(1,:))};
simTracks(j)=simD;


    
end
SPACE_UNITS='unit';
    TIME_UNITS='frames';
    ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
    ma = ma.addAll(simTracks);
    ma = ma.computeMSD;
save(fullfile(fold,['sarahDat',num2str(N),'.mat']),'simTracks','N','ma');
% clear all;