fold=uigetdir('D:\ChronoCode\chronoPkgs\Smarticles\matlabScripts\amoeba\smarticleExpVids');
f=dir2(fullfile(fold,'*.avi'));
    clf;
movs=struct;
nMovs=length(f);
movs(nMovs).fname='';
r=.192/2;%radius of boundary in meters
fps=zeros(nMovs,1);
conv=zeros(nMovs,1);
for i=1:nMovs
    pts(i,'/',nMovs);
    [movs(i).data,fps,conv]= trackMovie(r,fullfile(fold,f(i).name),80);
    movs(i).fname=f(i).name;
    movs(i).fps=fps;
    movs(i).conv=conv;
    beep;

    [~,vals]=parseFileNames(f(i).name);
    movs(i).pars=vals;
    pause(1);

end

save(fullfile(fold,'movieInfo.mat'),'movs','fold','r','nMovs')
