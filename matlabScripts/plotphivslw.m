% directory_name = uigetdir('\\centos\secured\shared_data\diffWidthVolFracFilling');
directory_name = '\\centos\secured\shared_data\Cyl4.4VolFracFilling\';
folds = dir(directory_name);
folds = folds(3:end);
dFile = '\PostProcess\volumefraction.txt';

%remove date from files and get all unique vals for lws
expression = '(\d+\.?\d*)';
% expression = '(\d+*\.?\d*)*-\d+-\d+';
a =(regexp({folds(:).name},expression,'match'));
lws=[];
for i=1:length(a)
    lws= [lws, str2double(a{i}{1})];
end
% lws = str2double([a{:}]);
[unilw ridx cidx]= unique(lws);
v=2.08026e-08;
rad =.0216;


x=[];y=[];
for i= 1:length(folds)
    data = importdata(horzcat(directory_name,folds(i).name,dFile));
    x= [x data(end,1)];
    y= [y mean(data(end-3:end,3))];
%     y= [y (data(end,2)*v/(rad^2*pi*2*data(end,5)))];
end


OutM= errBarCalc(lws,y);
errorbar(OutM(:,1),OutM(:,2),OutM(:,3),'o-');
axis([0 1.4 0.05 0.5])