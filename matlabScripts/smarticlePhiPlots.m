plotType=1; %1 for constant angle 2 for constant l/w
if plotType == 1 %l/w
    dv='lw';
    dvunits = '';
    dvlabel = 'l/w';
else             %angle 
    dv='ang';
    dvunits = '\circ';
    dvlabel = 'Angle (\circ)';
    
end

close all
directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\Results');
directory_name(length(directory_name)+1)='\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\constant l-w 1\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\UShape\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\AngledArms\';

files = dir(horzcat(directory_name,'2015*'));
% files = files(3:end);  %first 2 entries are  . and ..
figure(1);
hold on;
ylabel('\phi');
xlabel('Time(s)');

%initialize variables
depVar = zeros(length(files),1);
finalphi = zeros(length(files),2);
depVarLegend= cell(length(depVar),1);
complFound = zeros(length(files),1);


for j= 1:length(files)
    
    fname = files(j).name;
    pattern = horzcat(dv,'=([0-9.]*)');
    pars= regexp(fname, pattern, 'tokens');  %find amp and frequency from filename
    depVar(j)= str2double(cell2mat(pars{1})); %lw info for that fileoar
    data=importdata(horzcat(directory_name,files(j).name,'\PostProcess\volumeFraction.txt'));
    
    shakingIdx=(1);
    %plot stuff before shaking only
    %         shakingIdx = find(data(:,1)>4.9);
    %         plot(data(shakingIdx:end,1),data(shakingIdx:end,3));
    %         finalphi(j,1:2) = [ang(j),mean(data(shakingIdx-5:shakingIdx,3))];
    plot(data(:,1),data(:,3));
     %take average of last few points to get better final phi estimate
    finalphi(j,1:2) = [depVar(j),mean(data(end-5:end,3))];
    
    %look for completion line in sim param file to know which values in
    %plot aren't done yet if plotting while sims are still running
    FID = fopen(horzcat(directory_name,files(j).name,'\PostProcess\simulation_specific_parameters.txt'), 'r');
    if FID == -1, error('Cannot open file'), end
    Data = textscan(FID, '%s', 'delimiter', '\n', 'whitespace', '');
    CStr = Data{1};
    fclose(FID);
    if cell2mat(strfind(CStr(end),'completed'))
        complFound(j)= cell2mat(strfind(CStr(end),'completed'));
    else
        complFound(j)=0;
    end
    hold on;
end


for i=1:length(depVar)
    depVarLegend{i}=horzcat(num2str(depVar(i)),dvunits);.
end
legend(depVarLegend);
% figText(gcf,13); personal script which changes font size of all numbers
% in figure to 13
figure(3)
hold on;
xlabel(dvlabel);
ylabel('\phi_f');

%sort files in order
[Y,I]=sort(finalphi(:,1));
finalphi=finalphi(I,:);
complFound=complFound(I);
incomplIdx=find(~complFound);

% plot(finalphi(:,1),finalphi(:,2),'o-');
plot(finalphi(incomplIdx,1),finalphi(incomplIdx,2),'*r');
% figText(gcf,13);
%% errorbar calculation
[uout, ~, outidx] = unique( finalphi(:,1));
OutM = zeros(length(uout),3);
OutMatrix1  = [uout(:), accumarray( outidx, finalphi(:,2), [], @mean ) ];  %mean of phi'
OutMatrix2  = [uout(:), accumarray( outidx, finalphi(:,2), [], @std ) ]; %std of phi'
OutM(:,1:2) = OutMatrix1; % [angle mean phif ]
OutM(:,3)   = OutMatrix2(:,2); % [std ]


errorbar(OutM(:,1),OutM(:,2),OutM(:,3),'o-');