plotType=1; %1 for constant angle 2 for constant l/w


switch plotType
    case 1 %l/w
        dv='lw';
        dvunits = '';
        dvlabel = 'l/w';
        dvaxis =[0,1.4,0.05,0.40]
    case 2 %angle 
        dv='ang';
        dvunits = '\circ';
        dvlabel = 'Angle (\circ)';
        dvaxis =[0,120,0.05,.45]
    case 3
        dv='gamma';
        dvunits= ''; %technically should be (kgm/s^2)
        dvlabel = '\Gamma';
        
end

%volumeFraction.text
%time,   SmarticleInContainer,   phi,   height,   COM.z in container   
 
% close all
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
    depVarLegend{i}=horzcat(num2str(depVar(i)),dvunits);
end
legend(depVarLegend);

figText(gcf,13);% personal script which changes font size of all numbers
% in figure to 13
figure(3)
hold on;
xlabel(horzcat(dv,dvunits));
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
figure(10)
OutM=errBarCalc(finalphi(:,1),finalphi(:,2));
errorbar(OutM(:,1),OutM(:,2),OutM(:,3),'o-');
axis(dvaxis)
[uout, ~, outidx] = unique( finalphi(:,1));
xlabel(dvlabel);
ylabel('\phi_f');
figText(gcf,13);
hold all;
% %% phif = c*vp/ve
% figure(22);
% phif = OutM(:,2);
% %dimensions of smarticle
% t1=.00127*3;
% t2=.0005*3;
% w=.0117;
% l= .379*w;
% vp=(t1) * (t2)* (w + 2 * l);
% %ve=p*V; %find this through simulation
% %
%% quadratic fitting of angle
figure(12)
lw = [1.2 1 .7 .379 .25];
%y= a*x^2 -b*x + c
fitVals = [2.45e-5, -.00278, .3161;
           2.94e-5, -.00338, .3316; 
           2.34e-5, -0.003, .3722;
           1.69e-5, -.00228, .4204;
           1.98e-5, .00264 .4628];
suptitle('\phi_f vs. Angle for different l/w fit values (\phi_f=Ax^2+Bx+C)')
subplot(2,2,1);

plot(lw,fitVals(:,1),'.-','MarkerSize',25);
title('A vs. l/w');
xlabel('l/w');
ylabel('A');

subplot(2,2,2)
plot(lw,fitVals(:,2),'.-','MarkerSize',25);
title('B vs. l/w');
xlabel('l/w');
ylabel('B');

subplot(2,1,2)
plot(lw,fitVals(:,3),'.-','MarkerSize',25)

title('C vs. l/w');
xlabel('l/w');
ylabel('C');