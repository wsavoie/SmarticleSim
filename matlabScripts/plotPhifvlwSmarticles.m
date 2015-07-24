clear all
lw=[];
finalphi = [];
lwlegend={};
close all
directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\Results');
directory_name(length(directory_name)+1)='\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\UShape\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\AngledArms constant angle 120\';

getSingle = false;

files = dir(horzcat(directory_name,'\2015*'));
% files = files(3:end);  %first 2 entries are  . and ..

figure(1);
hold on;
ylabel('\phi');
xlabel('Time(s)');

finalphi = zeros(length(files),2);
lwlegend= cell(length(lw),1);
complFound = zeros(length(files),1);
lw = zeros(length(files),1);


for j= 1:length(files)
    fname = files(j).name;
    pattern = 'lw=([0-9.]*)';
    pars= regexp(fname, pattern, 'tokens');  %find amp and frequency from filename
    lw(j) = str2double(cell2mat(pars{1})); %lw info for that fileoar
    [data]=importdata(horzcat(directory_name,files(j).name,'\PostProcess\volumeFraction.txt'));
    
    FID = fopen(horzcat(directory_name,files(j).name,'\PostProcess\simulation_specific_parameters.txt'), 'r');
    if FID == -1, error('Cannot open file'), end
    Data = textscan(FID, '%s', 'delimiter', '\n', 'whitespace', '');
    CStr = Data{1};
    fclose(FID);
    
    shakingIdx = find(data(:,1),1);
    plot(data(:,1)-data(shakingIdx,1),data(:,3));
    finalphi(j,1:2) = [lw(j),mean(data(end-5:end,3))];
    
    if cell2mat(strfind(CStr(end),'completed'))
        complFound(j)= cell2mat(strfind(CStr(end),'completed'));
    else
        complFound(j)=0;
    end
    
    hold on;
    
end


for i=1:length(lw)
    lwlegend{i}=num2str(lw(i));
end
legend(lwlegend);
figText(gcf,13);
figure(3)
hold on;
xlabel('l/w');
ylabel('\phi_f');
[Y,I]=sort(finalphi(:,1));
complFound=complFound(I);
incomplIdx=find(~complFound);
finalphi=finalphi(I,:);

% plot(finalphi(:,1),finalphi(:,2),'o-');
plot(finalphi(incomplIdx,1),finalphi(incomplIdx,2),'*r');
figText(gcf,13);
%% errorbar calculation
[uout, ~, outidx] = unique( finalphi(:,1));
    OutM = zeros(length(uout),3);
    OutMatrix1 = [uout(:), accumarray( outidx, finalphi(:,2), [], @mean ) ];  %mean of G'
    OutMatrix2 = [uout(:), accumarray( outidx, finalphi(:,2), [], @std ) ]; %std of G'
    
    OutM(:,1:2) = OutMatrix1; % [l/w mean phif ]
    OutM(:,3) = OutMatrix2(:,2); % [std ]
    
    
    errorbar(OutM(:,1),OutM(:,2),OutM(:,3),'o-');