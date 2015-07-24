clear all
anglelegend={};
close all
directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\Results');
directory_name(length(directory_name)+1)='\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\constant l-w 1\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\UShape\';
% directory_name = 'D:\SimResults\Chrono\SmarticleU\Results\AngledArms\';

%added completed string to file on 2015722
%if getFileDate>2015722 check for completion
inProgress = true;
getSingle = false;
files = dir(horzcat(directory_name,'2015*'));
% files = files(3:end);  %first 2 entries are  . and ..
figure(1);
hold on;
ylabel('\phi');
xlabel('Time(s)');

ang = zeros(length(files),1);
finalphi = zeros(length(files),2);
anglelegend= cell(length(ang),1);
complFound = zeros(length(files),1);

if getSingle
    [data]=importdata('volfraction.text');
    shakingIdx = find(data(:,1),1);
    plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
    finalphi = [finalphi; ang(end),mean(data(end-5:end,3))];
else
    for j= 1:length(files)
        
        fname = files(j).name;
        pattern = 'ang=([0-9.]*)';
        pars= regexp(fname, pattern, 'tokens');  %find amp and frequency from filename
        ang(j)= str2double(cell2mat(pars{1})); %lw info for that fileoar
    	data=importdata(horzcat(directory_name,files(j).name,'\PostProcess\volumeFraction.txt'));
        
        shakingIdx=(1);
        %plot stuff before shaking only
%         shakingIdx = find(data(:,1)>4.9);
%         plot(data(shakingIdx:end,1),data(shakingIdx:end,3)); 
%         finalphi(j,1:2) = [ang(j),mean(data(shakingIdx-5:shakingIdx,3))];
        plot(data(:,1),data(:,3)); 
        finalphi(j,1:2) = [ang(j),mean(data(end-5:end,3))];
        
        
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
end

for i=1:length(ang)
    anglelegend{i}=num2str(ang(i));
    anglelegend{i}=horzcat(anglelegend{i},'\circ');
end
legend(anglelegend);
figText(gcf,13);
figure(3)
hold on;
xlabel('angle (\circ)');
ylabel('\phi_f');

%sort files in order
[Y,I]=sort(finalphi(:,1));
finalphi=finalphi(I,:);
complFound=complFound(I);
incomplIdx=find(~complFound);

% plot(finalphi(:,1),finalphi(:,2),'o-');
plot(finalphi(incomplIdx,1),finalphi(incomplIdx,2),'*r');
figText(gcf,13);
%% errorbar calculation
[uout, ~, outidx] = unique( finalphi(:,1));
    OutM = zeros(length(uout),3);
    OutMatrix1  = [uout(:), accumarray( outidx, finalphi(:,2), [], @mean ) ];  %mean of G'
    OutMatrix2  = [uout(:), accumarray( outidx, finalphi(:,2), [], @std ) ]; %std of G'
    OutM(:,1:2) = OutMatrix1; % [l/w mean phif ]
    OutM(:,3)   = OutMatrix2(:,2); % [std ]
    
    
    errorbar(OutM(:,1),OutM(:,2),OutM(:,3),'o-');