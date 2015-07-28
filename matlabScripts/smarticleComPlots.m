dvunits = '';
dvlabel = 'h(t)/h_0';
dvaxis =[0,1.4,0.05,0.40]
dv=' g'
%volumeFraction.text
%time,   SmarticleInContainer,   phi,   height,   COM.z in container   
 
close all
directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\Results');
directory_name(length(directory_name)+1)='\';

files = dir(horzcat(directory_name,'2015*'));
figure(1);
hold on;

depVar = zeros(length(files),1);
depVarLegend= cell(length(depVar),1);
for j= 1:length(files)
    
    fname = files(j).name;
    pattern = horzcat(dv,'=([0-9.]*)');
    pars= regexp(fname, pattern, 'tokens');  %find amp and frequency from filename
    depVar(j)= str2double(cell2mat(pars{1})); %lw info for that fileoar
    data=importdata(horzcat(directory_name,files(j).name,'\PostProcess\volumeFraction.txt'));
    h_0=data(1,5);
    plot(data(:,1),data(:,5)/h_0);
    hold on;
     %take average of last few points to get better fin1al phi estimate
    
    %look for completion line in sim param file to know which values in
    %plot aren't done yet if plotting while sims are still running
%     FID = fopen(horzcat(directory_name,files(j).name,'\PostProcess\simulation_specific_parameters.txt'), 'r');
%     if FID == -1, error('Cannot open file'), end
%     Data = textscan(FID, '%s', 'delimiter', '\n', 'whitespace', '');
%     CStr = Data{1};
%     fclose(FID);
%     if cell2mat(strfind(CStr(end),'completed'))
%         complFound(j)= cell2mat(strfind(CStr(end),'completed'));
%     else
%         complFound(j)=0;
%     end
%     hold on;
end
for i=1:length(depVar)
    depVarLegend{i}=horzcat(num2str(depVar(i)),dvunits);
end
legend(depVarLegend);

xlabel('t(s)');
ylabel(dvlabel);