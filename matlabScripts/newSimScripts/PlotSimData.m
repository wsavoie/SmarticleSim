fold=uigetdir('B:\SmartSimResults\12-5');
fname='vibData.mat';
if ~exist(fullfile(fold,fname),'file')
    readSimData(fold);
end
load(fullfile(fold,'vibData.mat'));
%************************************************************
%* Fig numbers:
%* 1. plot phi evo vs time for single trial
%* 2. plot phi evo with gui vs time
%* 3. plot phi change vs vibration amp
%************************************************************

showFigs=[3];
inds=1;
lw=[]; nl=[]; npl=[]; vib=[]; N=[]; v=[];
props={lw nl npl vib N v}; 
[lws,nls,npls,vibs,vs]=separateVec(pars,1);
Ns=nls.*npls;

for i=1:length(dat)
    
    cond=true;
    for j=1:length(props)
        %if empty accept all values
        if ~isempty(props{j})
            %in case multiple numbers in property
            %if no matches then cond = false
            if(~any(props{j}==dat(i).pars(j)))
                cond = false;
            end
        end
    end
    if(cond)
        usedD(inds)=dat(i);
        inds=inds+1;
    end
end
N=length(usedD);
I=6;%ind to be used for multiple different plot numbers below
h=figure(1000);
guiCols=get(gca,'colororder');
close(h);

%nick dat before shaking
lwBS=[0, 0.1 0.2 0.3 0.4 0.6 0.7 1.3 1.4];
phiBS=[.227,0.2,0.18 0.16 .1425,.115,.105,0.061,.06];
%nick data after shaking
nickD=[[0:.1:1.4]',[.27 .235 .208 .185 .163 .148 .136 .13 .118 .11 .105 .1025 .095 .09 .085]'];
% xlabel('l/w');
% ylabel('\phi');
% axis([0 1.4 0.05 0.3]);
% figText(gcf,16);

%% 1. plot phi evo vs time for single trial
xx=1;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    
    t=usedD(I).t;       phi=usedD(I).phi;
    downsample(t,100);  downsample(phi,100);
    plot(t,phi,'linewidth',2);
    xlabel('t');
    ylabel('\phi');
    figText(gcf,18);
  
end
%% 2. plot phi evo with gui vs time
xx=2;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    
    t=usedD(I).t;       
    phi=usedD(I).phi;
    gui=usedD(I).gui;
    downsample(t,100);
    downsample(phi,100);
    downsample(gui,100);
    %find changes in guid
    
    plot(t,phi,'linewidth',2);
    %find changes in guid
    gc=findChangesInGui(gui);
    axis tight
    for(i=1:size(gc,1))
        h(i)=plot([t(gc(i,1)),t(gc(i,1))],get(gca,'ylim'),'linewidth',3,...
            'color',guiCols(gc(i,2),:));
    end

    xlabel('t');
    ylabel('\phi');
    figText(gcf,18);
end

%% 3. plot phi change vs vibration amp
xx=3;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    
%     unLW=unique([usedD.lw]);
%     unVib=unique([usedD.vib]);
    
%     lws=[0.7,0.7,0.7,0.6 0.6];
%     vibs=[1,2,4,1,4];
    unLW=unique([lws]);
    unVib=unique([vibs]);
    
    
    for(i = 1:length(unLW))
        vibVals=[];
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
%             vals=length(cond(cond==1));
            idz=find(cond);
            for k=1:length(idz)
                [dphi,phi0,phif]=changeInPhiAmp(usedD(idz(k)),2);
                out(k)=dphi;
                out2(k,:)=[phi0,phif];
            end
            vibVals(j)=unVib(j);
            phiVals(j)={out};
            phiFVals(j)={out2(:,2)};
        end
        phiValsM=cellfun(@(x) mean(x),phiVals);
        phiValsE=cellfun(@(x) std(x),phiVals);
        errorbar(vibVals,phiValsM,phiValsE,'linewidth',2);
    end
    

    xlabel('Vibration Amp (\circ)');
    ylabel('\Delta\phi');
    figText(gcf,18);
end





% phiM=mean(phi);
% ro=.0234;%bucket radius without thickness taken into account in sim
% rn=0.0224;%corrected bucket radius
% %bucket rad was declared wrong in sim this fixes it
% phiM=phiM*(ro/rn)^2;
% errorbar(lw,phiM,-std(phi)/2,std(phi)/2,0.0714*ones(1,length(lw)),0.0714*ones(1,length(lw)),'linewidth',2);
% 
%     plot(lwBS,phiBS,'o-','linewidth',2);
%     plot(nickD(:,1),nickD(:,2),'ro-','linewidth',2);
%    
% legend({'chrono data (\phi_i) l''=l\pmD','Nick mean (\phi_i)','Nick mean \phi_f'})
