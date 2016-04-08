path = 'D:\SimResults\Chrono\SmarticleU\tests\BoxAngChangeTorPct30v2\';
a=dir(horzcat(path,'-*'))

for i=1:length(a)
    ff= horzcat(path,a(i).name,'\PostProcess\Stress.txt')
    readAllSmarticlesAngles(ff,0);
end
