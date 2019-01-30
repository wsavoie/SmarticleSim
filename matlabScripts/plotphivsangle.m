% directory_name = uigetdir('\\centos\secured\shared_data\diffWidthVolFracFilling');
% directory_name = '\\centos\secured\shared_data\Cyl4.4VolFracFilling\';
directory_name = 'D:\SimResults\Chrono\SmarticleU\tests\angChange2\';
folds = dir(directory_name);
folds = folds(3:end);
dFile = '\PostProcess\volumefraction.txt';
dFile2 = '\PostProcess\simulation_specific_parameters.txt';

%remove date from files and get all unique vals for lws
expression = '(\d+\.?\d*)';
% expression = '(\d+*\.?\d*)*-\d+-\d+';
a =(regexp({folds(:).name},expression,'match'));
lws=[];
angs=[];
for i=1:length(a)
    lws= [lws, str2double(a{i}{1})];
    angs= [angs,str2double(a{i}{3})];
end
% lws = str2double([a{:}]);
[unilw ridx cidx]= unique(angs);
v=2.08026e-08;
rad =.0216;
nlws=[];
nangs=[];
volz = [];
ww= 0.0117;
x=[];y=[];y2=[];
for i= 1:length(folds)
    try
    data = importdata(horzcat(directory_name,folds(i).name,dFile));
    %%%%%%%%%%%%%%%%%%
%     mZ = mean(data(end-0:end,4));%mean zmax
%     c = data(end,2);
%     v2= .0005*.0005*.00127; %extra volume
%     FID = fopen(horzcat(directory_name,folds(i).name,dFile2));
%     formatSpec='%s';
%     Data = textscan(FID,formatSpec,...            
%                 'Delimiter', '\n', ...
%                 'CollectOutput', true);
%     CStr = Data{1};
%     d= textscan(CStr{12},'%s');
%     d=d{1,1};
%     vol=d{3};
%     vol=str2double(d{3});
%     fclose(FID);
%     y2 = [y2 ((vol+2*v2)*c)/(mZ*pi*.022^2)];
%     volz = [volz vol];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     nlws = [nlws (lws(i)+.0005/ww)];
%      nlws = [nlws lws(i)*ww/(ww-.0005*2)];
    nlws = [nlws lws(i)];
    nangs= [nangs angs(i)];
    catch exception2    
        pts('caught error')
        continue;
    end  
    x= [x data(end,1)];
    y= [y mean(data(end-0:end,3))];
%     y= [y (data(end,2)*v/(rad^2*pi*2*data(end,5)))];
end
lws = nlws;
angs = nangs;
hold all;
OutM= errBarCalc(nangs,y);
errorbar(OutM(:,1),OutM(:,2),OutM(:,3),'o-');
axis([0 120 0.12 0.155])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OutM2= errBarCalc(lws,y2);
% errorbar(OutM2(:,1),OutM2(:,2),OutM2(:,3),'yo-');
xlabel('Angle (\circ)');
ylabel('\phi_f');