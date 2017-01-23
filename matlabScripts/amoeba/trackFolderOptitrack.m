fold=uigetdir('D:\ChronoCode\chronoPkgs\Smarticles\matlabScripts\amoeba\smarticleExpVids');
f=dir2(fullfile(fold,'*.csv'));
%     clf;
movs=struct;
nMovs=length(f);
movs(nMovs).fname='';
r=.192/2;%radius of boundary in meters
dec=12; %decimate amount
conv=zeros(nMovs,1);
for i=1:nMovs
    pts(i,'/',nMovs);
%     [t,x,y,tracks]
    [movs(i).t,movs(i).x,movs(i).y,movs(i).data,movs(i).rot]= trackOptitrack(fullfile(fold,f(i).name),dec);
    movs(i).fname=f(i).name;
    movs(i).fps=120/dec;
    movs(i).conv=1;
    [~,vals]=parseFileNames(f(i).name);
    %%
%     spk=[0]; smart=[-90]; gait=[1]; rob=[5]; v=[nMovs];
    vals=[0 0 1 5 i];
    %%
    movs(i).pars=vals;
    
    
%     pause(1);

end

save(fullfile(fold,'movieInfo.mat'),'movs','fold','nMovs')
