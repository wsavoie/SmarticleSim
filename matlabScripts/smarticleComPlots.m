dvunits = '';
dvlabel = 'h(t)/h_0';
dvaxis =[0,1.4,0.05,0.40]
dv=' g'
%volumeFraction.text
%time,   SmarticleInContainer,   phi,   height,   COM.z in container   
 
close all
directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests');
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
    t_0idx=find(data(:,1)>.49,1); %shaking starts at .5
    t_0=data(t_0idx,1);
    h_0=data(t_0idx,5);
    
    xdat =data(t_0idx:end,1)-t_0;
    ydat = data(t_0idx:end,5)/h_0;
%     plot(data(t_0idx:end,1)-t_0,data(t_0idx:end,5)*2/h_0-1);
     plot(xdat,ydat);
    hold on;
     
    
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
% legend(depVarLegend);
set(gca, 'XScale', 'log')
xlabel('t(s)');
ylabel(dvlabel);

f= 30;
gamma=[1.23 1.48 1.7 1.96 2.20 2.53];
beta= [.5 .7 .8 .9 1 1.5]; %changes with gamma
delta=14.3; %changes with l/w
t=logspace(-1.5,6,10000);
for i=1:length(gamma)
    tau = 1/(30)*exp(delta/gamma(i));
    htho = exp(-(t./tau)).^beta(i);
%     htho = exp(-t./tau).^beta(i);
    depVarLegend{length(depVarLegend)+1}=horzcat('\beta= ',num2str(beta(i)),' \Gamma= ',num2str(gamma(i)));
    plot(t,htho)
end
legend(depVarLegend);
axis([10^(-1.5) 10^3.5,0,1.0]);



data=importdata(horzcat(directory_name,files(1).name,'\PostProcess\volumeFraction.txt'));
t_0idx=find(data(:,1)>.49,1); %shaking starts at .5
t_0=data(t_0idx,1);
h_0=data(t_0idx,5);

xdat =data(t_0idx:end,1)-t_0;
ydat = data(t_0idx:end,5)/h_0;

delta=15;
beta = 1.48;
parsin = [delta,beta];
PON = 0;
% fitStretched(parsin,c,xdat,ydat,PON)
deltaFit = zeros(length(files),1);
betaFit = zeros(length(files),1);

for i= 1:length(files)
    fname = files(i).name;
    pattern = horzcat(dv,'=([0-9.]*)');
    pars= regexp(fname, pattern, 'tokens');  %find amp and frequency from filename
    depVar(i)= str2double(cell2mat(pars{1})); %lw info for that fileoar
    data=importdata(horzcat(directory_name,files(i).name,'\PostProcess\volumeFraction.txt'));
    t_0idx=find(data(:,1)>.3,1); %shaking starts at .5
    t_0=data(t_0idx,1);
%     h_0=data(t_0idx,5);
    h_0=data(t_0idx,5);
    xdat =data(t_0idx:end,1)-t_0;
%     ydat = data(t_0idx:end,5)*2/h_0-1;
    ydat = data(t_0idx:end,5)/h_0;
    out =fitStretched(parsin,depVar(i),xdat,ydat,PON);
    deltaFit(i)=out(1,1);
    betaFit(i)=out(1,2);
%          plot(xdat,ydat);
%     hold on;
end


figure(44)
plot(depVar,betaFit,'.b-','MarkerSize',25);
xlabel('\Gamma')
ylabel('\beta')

figure(45)
plot(depVar,deltaFit,'.b-','MarkerSize',25);
xlabel('\Gamma')
ylabel('\Delta')


% plot(1./depVar,log(1/30*exp(deltaFit./depVar)),'g.-','MarkerSize',25);


figure(46)
[uout, ~, outidx] = unique( 1./depVar);
OutM = zeros(length(uout),3);
OutMatrix1  = [uout(:), accumarray( outidx,log(1/30*exp(deltaFit./depVar)), [], @mean ) ];  %mean of phi'
OutMatrix2  = [uout(:), accumarray( outidx, log(1/30*exp(deltaFit./depVar)), [], @std ) ]; %std of phi'
OutM(:,1:2) = OutMatrix1; % [angle mean phif ]
OutM(:,3)   = OutMatrix2(:,2); % [std ]
errorbar(OutM(:,1),OutM(:,2),OutM(:,3),'o-');
ylabel('ln(\tau/\tau_0)')
xlabel('\Gamma^{-1}');
