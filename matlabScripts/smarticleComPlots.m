plotType=2;
PON = 0;
switch plotType
    case 1 %gamma
        dvunits = '';
        dvlabel = '\Gamma';
        dvaxis =[0,1.4,0.05,0.40]
        dv=' g'
        ptype= 'gamma';
    case 2 %angle
        dv='ang1';
        dvunits = '\circ';
        dvlabel = 'Angle (\circ)';
        ptype= 'angle';
        dvaxis =[0,120,0.05,.45]
    case 3 %diff angles, NOT IMPLEMENTED
        dv='gamma';
        dvunits= ''; %technically should be (kgm/s^2)
        dvlabel = '\Gamma';
        
end


%volumeFraction.text
%time,   SmarticleInContainer,   phi,   height,   COM.z in container   
 
close all
% directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests');
directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\Results\Shake with Diff Angles');
directory_name(length(directory_name)+1)='\';

files = dir(horzcat(directory_name,'2015*'));
figure(1);
hold on;

depVar = zeros(length(files),1);
depVarLegend= cell(length(depVar),1);
for j= 1:length(files)
    
    fname = files(j).name;
    pattern = horzcat(dv,'=([0-9.]*)');
    pars= regexp(fname, pattern, 'tokens'); 
    try
    depVar(j)= str2double(cell2mat(pars{1})); %lw info for that fileoar
    catch err
        pts('plottype=',plotType);
        error(horzcat('Try switching plotType! You''re plotting with depvar==',ptype));
    end
    data=importdata(horzcat(directory_name,files(j).name,'\PostProcess\volumeFraction.txt'));
    t_0idx=find(data(:,1)>.35,1); %shaking starts at .5
    t_0=data(t_0idx,1);
    h_0=data(t_0idx,5);
    
    xdat =data(t_0idx:end,1)-t_0;
    ydat = data(t_0idx:end,5)/h_0;
%     plot(data(t_0idx:end,1)-t_0,data(t_0idx:end,5)*2/h_0-1);
     plot(xdat,ydat);
    hold on;
end
for i=1:length(depVar)
    depVarLegend{i}=horzcat(num2str(depVar(i)),dvunits);
end
% legend(depVarLegend);
set(gca, 'XScale', 'log')
xlabel('t(s)');
ylabel('h(t)/h_0');
axis([10^(-1.5) 10^3.5,0,1.0]);
%% plot approximate lines from paper
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

%% start of determining beta and delta
delta=15;
beta = 1.48;
parsin = [delta,beta];

% fitStretched(parsin,c,xdat,ydat,PON)
deltaFit = zeros(length(files),1);
betaFit = zeros(length(files),1);

for i= 1:length(files)
    fname = files(i).name;
    pattern = horzcat(dv,'=([0-9.]*)');
    pars= regexp(fname, pattern, 'tokens');  %find amp and frequency from filename
    
    if(plotType ~=1)
        
        pattern2 = horzcat(' g','=([0-9.]*)'); 
        gammaPars= regexp(fname, pattern2, 'tokens');  %find amp and frequency from filename
    else
        gammaPars = pars;
    end
   
    depVar(i)= str2double(cell2mat(pars{1})); %lw info for that fileoar
    gammaVal(i) = str2double(cell2mat(gammaPars{1}));
    data=importdata(horzcat(directory_name,files(i).name,'\PostProcess\volumeFraction.txt'));
    t_0idx=find(data(:,1)>.35,1); %shaking starts at .5
    t_0=data(t_0idx,1);
    h_0=data(t_0idx,5);
    xdat =data(t_0idx:end,1)-t_0;
    ydat = data(t_0idx:end,5)/h_0;
    out =fitStretched(parsin,gammaVal(i),xdat,ydat,PON);
    deltaFit(i)=out(1,1);
    betaFit(i)=out(1,2);
%          plot(xdat,ydat);
%     hold on;
end


figure(44)
betaPlot = errBarCalc(depVar,betaFit);
hold on;
% plot(depVar,betaFit,'.b-','MarkerSize',25);
errorbar(betaPlot(:,1),betaPlot(:,2),betaPlot(:,3),'o-');
xlabel(dvlabel)
ylabel('\beta')

figure(45)
hold on;
deltaPlot = errBarCalc(depVar,deltaFit);
% plot(depVar,deltaFit,'.b-','MarkerSize',25);
errorbar(deltaPlot(:,1),deltaPlot(:,2),deltaPlot(:,3),'o-');
xlabel(dvlabel)
ylabel('\Delta')

if(plotType==1)%only works for gamma varied runs
    
    figure(10)
    lntauvgamma=errBarCalc(1./depVar,log(1/30*exp(deltaFit./depVar)));
    errorbar(lntauvgamma(:,1),lntauvgamma(:,2),lntauvgamma(:,3),'o-');
    hold on;
    ylabel('ln(\tau/\tau_0)')
    xlabel('\Gamma^{-1}');
end
