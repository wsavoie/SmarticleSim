plotType=3;
PON =0;
cutoff = .25;
f= 30; %shaking frequency
switch plotType
    case 1 %gamma 
        dvunits = '';
        dvlabel = '\Gamma';
        dvaxis =[0,1.4,0.05,0.40];
        dvIdx = 3;
        ptype= 'gamma';
    case 2 %angle
        dv='ang1';
        dvIdx = 2;
        dvunits = '\circ';
        dvlabel = 'Angle (\circ)';
        ptype= 'angle';
        dvaxis =[0,120,0.05,.45];
    case 3 %reading in different angles at different gamma
        dv='ang1';
        dvIdx=2;
        dvunits = '\circ';
        dvlabel = 'Angle (\circ)';
        ptype= 'angle';
        dvaxis =[0,120,0.05,.45];
    case 4 %diff angles, NOT IMPLEMENTED
        dv='gamma';
        dvunits= ''; %technically should be (kgm/s^2)
        dvlabel = '\Gamma';
        
end


%volumeFraction.text
%time,   SmarticleInContainer,   phi,   height,   COM.z in container   
 
% close all
% directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests');
directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\Results\ang=90 lw=varied zcom');
directory_name(length(directory_name)+1)='\';

files = dir(horzcat(directory_name,'2015*'));

depVar = zeros(length(files),6);
depVarLegend= cell(length(depVar),1);

%input params 
delta=10;
beta = .5;
parsin = [delta,beta];
for j= 1:length(files)
    figure(1);
    hold on;

    fname = files(j).name;
    lwPatt = horzcat(' lw','=([0-9.]*)');
    lws= regexp(fname, lwPatt, 'tokens');
    
    angPatt = horzcat(' ang1','=([0-9.]*)');
    angs= regexp(fname, angPatt, 'tokens');
    
    gamPatt = horzcat(' g','=([0-9.]*)');
    gams= regexp(fname, gamPatt, 'tokens'); 
    
    depVar(j,1)= str2double(cell2mat(lws{1})); %lw info for that fileoar
    depVar(j,2)= str2double(cell2mat(angs{1}));
    depVar(j,3)= str2double(cell2mat(gams{1})); %lw info for that fileoar
    
    data=importdata(horzcat(directory_name,files(j).name,'\PostProcess\volumeFraction.txt'));
    t_0idx=find(data(:,1)>cutoff,1); %shaking starts at .5
    t_end=find(data(:,1)>data(end,1)*1/exp(1),1);
    t_0=data(t_0idx,1);
    h_0=data(t_0idx,5);
%     xdat =data(t_0idx:t_end,1)-t_0;
%     ydat = data(t_0idx:t_end,5)/h_0;
    xdat =data(t_0idx:end,1)-t_0;
    ydat = data(t_0idx:end,5)/h_0;
%     plot(data(t_0idx:end,1)-t_0,data(t_0idx:end,5)*2/h_0-1);

%%%%%%%%%filter%%%%%%%%%%%%
n = 4; % designs an Nth order low pass
timestep = .00075;
we = f*2*pi*timestep;
if we >1
    we = .05; %forgot what this does
end
[B,A] = butter(n,we);
% filters the data in vector X with the filter
%   described by vectors A and B to create the filtered data Y.  The
%   filter is described by the difference equation:
%
%     a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                           - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
ydat2 = filtfilt(B,A,ydat);

% set(gca,'XScale','log');

if PON
    figure(20);
end
    modelFun = @(a,q) exp(-(q./(1/f*exp(a(1)/depVar(j,3)))).^a(2));
%     out =fitStretched(parsin,depVar(j,3),xdat,ydat2,PON);
    w = [1./log(xdat+1.00001)];
    out1 = fitnlm(xdat,ydat2,modelFun,parsin,'Weight',w);
    out= [out1.Coefficients{1,1},out1.Coefficients{2,1}];
    depVar(j,4:5)=[out(1,1:2)]; %[delta,beta]
    figure(1)
    hold on;
        xlabel('time(s)');
        ylabel('h(t)/h_0');
    if PON
        hold off;
        plot(xdat,ydat2);
%         set(gca, 'XScale', 'log')
        hold on;
        tau = 1/30*exp(depVar(j,4)/depVar(j,3));
    %     plot(xdat,exp(-(xdat./(tau)).^out(2)),'g-o');
        plot(xdat,exp(-(xdat./depVar(j,4)).^out(2)),'g.-');
    %     plot(xdat,-exp(-(xdat./(tau)).^out(2)),'g-o');
        pause(1)
        
    else
        h=plot(xdat,ydat);
        c = h.Color;
        plot(xdat,ydat2,'Color',c);
        set(gca, 'XScale', 'log')
        
    end
       
end


figure(1);
hold on;
%% calculting tau
depVar(:,6)= 1/f*exp(depVar(:,4)./depVar(:,3));
% 
% %remove all tau<16 ?
% [gtauR, gtauC]= find(log(depVar(:,6))<16);
% depVar=depVar(gtauR,:);


figText(gcf,16)
uniLws= unique(depVar(:,1));
uniAngs = unique(depVar(:,2));
uniGams = unique(depVar(:,3));
for i=1:size(depVar,1)
    depVarLegend{i}=horzcat('l/w=',num2str(depVar(i,1)),' ', num2str(depVar(i,2)),'\circ',' \Gamma=',num2str(depVar(i,3)));
end
depVarLegend= depVarLegend(~cellfun(@isempty, depVarLegend));
legend(depVarLegend);
set(gca, 'XScale', 'log')
xlabel('t(s)');
ylabel('h(t)/h_0');
axis([10^(-1.5) 10^1,0,1.0]);


% % plot approximate lines from paper
% 
% gamma=[1.23 1.48 1.7 1.96 2.20 2.53];
% beta= [.7 .6 .7 .8 .9 1]; %changes with gamma
% delta=5; %changes with l/w
% t=logspace(-6,4,10000);
% for i=1:length(gamma)
% 
%     tau = 1/f*exp(delta/gamma(i));
%     htho = exp(-(t./tau).^beta(i));
% %     htho = exp(-t./tau).^beta(i);
%     depVarLegend{length(depVarLegend)+1}=horzcat('\beta= ',num2str(beta(i)),' \Gamma= ',num2str(gamma(i)));
%     plot(t,htho)
% end
% legend(depVarLegend);


%% plotting beta
if plotType~=3 && plotType ~=2
    figure(44)
    betaPlot = errBarCalc(depVar(:,1),depVar(:,5));
    hold on;
    errorbar(betaPlot(:,1),betaPlot(:,2),betaPlot(:,3),'o-');
    xlabel(dvlabel)
    ylabel('\beta')
end

%% plotting delta
figure(45)
hold on;
if plotType==3
    deltaPlot=errBarCalc(depVar(:,1),depVar(:,4));
    errorbar(deltaPlot(:,1),deltaPlot(:,2),deltaPlot(:,3),'o-');
    xlabel('l/w');
    ylabel('\Delta');
end

tauLeg= {};
if (plotType==1 || plotType==3)%tau vs gamma^-1 only works for gamma varied runs
    figure(10)
    hold on;
    for i=1:length(uniLws) %change this to unique(depVar(:,dv))
        [ri,ci] = find(depVar(:,1)==uniLws(i));
        tau = 1/f*exp(depVar(ri,4)./depVar(ri,3));
        lntauvgamma = errBarCalc(1./depVar(ri,3),log(tau));
%         lntauvgamma = errBarCalc(1./depVar(ri,3),1/30*exp(depVar(ri,4)./depVar(ri,3)));
        errorbar(lntauvgamma(:,1),3*i+lntauvgamma(:,2),lntauvgamma(:,3),'o-');
        tauLeg{length(tauLeg)+1}=horzcat('l/w=',num2str(uniLws(i)),'');

    end
        ylabel('ln(\tau/\tau_0)+const')
%         ylabel('\tau/\tau_0+const')
        xlabel('\Gamma^{-1}');
        legend(tauLeg);
        figText(gcf,13);
        title(horzcat('Angle=',num2str(depVar(1,2)), '\circ'))
%         set(gca, 'YScale', 'log')
%         [rft stats] = robustfit(lntauvgamma(:,1),3*i+lntauvgamma(:,2));
end
if(plotType==4)
    
end
