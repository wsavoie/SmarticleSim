path = 'A:\SmarticleRun\un_eq_no_OT_all_plots_4runs\UnequalOTv2\';
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