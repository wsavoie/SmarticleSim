plotType=3;
PON = 0;
cutoff = .05;
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
directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\Results\ang=60 lw = varied');
directory_name(length(directory_name)+1)='\';

files = dir(horzcat(directory_name,'2015*'));
figure(1);
hold on;

depVar = zeros(length(files),5);
depVarLegend= cell(length(depVar),1);

%input params 
deltaFit = zeros(length(files),3);
betaFit = zeros(length(files),3);
delta=15;
beta = .68;
parsin = [delta,beta];
for j= 1:length(files)
    
    fname = files(j).name;
    lwPatt = horzcat(' lw','=([0-9.]*)');
    lws= regexp(fname, lwPatt, 'tokens');
    
    angPatt = horzcat(' ang1','=([0-9.]*)');
    angs= regexp(fname, angPatt, 'tokens');
    
    gamPatt = horzcat(' g','=([0-9.]*)');
    gams= regexp(fname, pattern2, 'tokens'); 
    
    depVar(j,1)= str2double(cell2mat(lws{1})); %lw info for that fileoar
    depVar(j,2)= str2double(cell2mat(angs{1}));
    depVar(j,3)= str2double(cell2mat(gams{1})); %lw info for that fileoar
    
    data=importdata(horzcat(directory_name,files(j).name,'\PostProcess\volumeFraction.txt'));
    t_0idx=find(data(:,1)>cutoff,1); %shaking starts at .5
    t_0=data(t_0idx,1);
    h_0=data(t_0idx,5);
    
    xdat =data(t_0idx:end,1)-t_0;
    ydat = data(t_0idx:end,5)/h_0;
%     plot(data(t_0idx:end,1)-t_0,data(t_0idx:end,5)*2/h_0-1);
     plot(xdat,ydat);
    hold on;
    
    
    out =fitStretched(parsin,depVar(j,3),xdat,ydat,PON);
    depVar(j,4:5)=[out(1,1:2)]; %[delta,beta]

%     deltaFit(j,:)=[depVar(j,1:2),out(1,1)];
%     betaFit(j,:)=[depVar(j,1:2),out(1,2)];
end
figText(gcf,16)
uniLws= unique(depVar(:,1));
uniAngs = unique(depVar(:,2));
uniGams = unique(depVar(:,3));
for i=1:size(depVar,1)
    depVarLegend{i}=horzcat('l/w=',num2str(depVar(i,1)),' ', num2str(depVar(i,2)),'\circ',' \Gamma=',num2str(depVar(i,2)));
end
legend(depVarLegend);
set(gca, 'XScale', 'log')
xlabel('t(s)');
ylabel('h(t)/h_0');
axis([10^(-1.5) 10^3.5,0,1.0]);
%% plot approximate lines from paper
% f= 30;
% gamma=[1.23 1.48 1.7 1.96 2.20 2.53];
% beta= [.5 .7 .8 .9 1 1.5]; %changes with gamma
% delta=14.3; %changes with l/w
% t=logspace(-1.5,6,10000);
% for i=1:length(gamma)
% 
%     tau = 1/(30)*exp(delta/gamma(i));
%     htho = exp(-(t./tau)).^beta(i);
% %     htho = exp(-t./tau).^beta(i);
%     depVarLegend{length(depVarLegend)+1}=horzcat('\beta= ',num2str(beta(i)),' \Gamma= ',num2str(gamma(i)));
%     plot(t,htho)
% end
% legend(depVarLegend);


%% plotting beta
% deltaFit=sortrows(deltaFit,1);
% betaFit=sortrows(betaFit,1);

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
% deltaPlot = errBarCalc(depVar(:,1),deltaFit(:,3));
% % plot(depVar,deltaFit,'.b-','MarkerSize',25);
% errorbar(deltaPlot(:,1),deltaPlot(:,2),deltaPlot(:,3),'o-');
% xlabel(dvlabel)
% ylabel('\Delta')
tauLeg= {};
if (plotType==1 || plotType==3)%tau vs gamma^-1 only works for gamma varied runs
    figure(10)
    hold on;
    for i=1:length(uniLws) %change this to unique(depVar(:,dv))
        [ri,ci] = find(depVar(:,1)==uniLws(i));
        
%         lntauvgamma = errBarCalc(1./depVar(ri,3),log(1/30*exp(depVar(ri,4)./depVar(ri,3))));
        lntauvgamma = errBarCalc(1./depVar(ri,3),1/30*exp(depVar(ri,4)./depVar(ri,3)));
        errorbar(lntauvgamma(:,1),3*i+lntauvgamma(:,2),lntauvgamma(:,3),'o-');
        tauLeg{length(tauLeg)+1}=horzcat('l/w=',num2str(uniLws(i)),'');
        
%     for i=1:length(uniDepVar)
%         [rI,cI]=find(deltaFit(:,1)==uniDepVar(i));%find all gammas used for each depVar
%         gams = deltaFit(rI,2);
%         %next line could be written more arbitrarily where make sure
%         %unidepvar is equal to the indice used in deltaPlot but it works
%         %for now
%         lntauvgamma = errBarCalc(1./gams,log(1/30*exp(deltaFit(rI,3)./gams)));
%         errorbar(lntauvgamma(:,1),3*i+lntauvgamma(:,2),lntauvgamma(:,3),'o-');
%         tauLeg{length(tauLeg)+1}=horzcat(num2str(uniDepVar(i)),dvunits);
%     end
    end
%         ylabel('ln(\tau/\tau_0)+const')
        ylabel('\tau/\tau_0+const')
        xlabel('\Gamma^{-1}');
        legend(tauLeg);
        figText(gcf,13);
        title('Angle=60\circ')
        set(gca, 'YScale', 'log')
        
%     lntauvgamma=errBarCalc(1./depVar,log(1/30*exp(deltaFit./depVar)));
%     end
%     lntauvgamma=errBarCalc(1./depVar(:,1),log(1/30*exp(deltaFit./depVar(:,1))));
%     errorbar(lntauvgamma(:,1),lntauvgamma(:,2),lntauvgamma(:,3),'o-');

end
if(plotType==4)
    
end
