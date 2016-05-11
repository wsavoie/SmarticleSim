path = 'A:\SmarticleRun\lazy_tor.3_smarts30\lazy .00125\';
% path = 'D:\SimResults\Chrono\SmarticleU\tests\lazy .05\';
a=dir(horzcat(path,'-*'))
for i=1:length(a)
    ff= horzcat(path,a(i).name,'\PostProcess\Stress.txt');
    pts(ff);
    matFile =  horzcat(path,a(i).name,'\PostProcess\stressData.mat');
    if exist(matFile, 'file') == 2 
        pts('matrix file already created');
        continue
    end
   
    readAllSmarticlesAngles(ff,0);
end
beep;
beep;