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
%* 4. plot phif vs vibration amp
%* 5. plot phi_i, phi_f for nicks data and for vib data
%* 6. plot phi evolution vs vibration amp for many runs
%* 7. plot <N> vs lw
%* 8. plot pdf of motor torque vs <N>
%* 9. plot radial and z frc on cylinder vs time
%*10. plot radial frc on cylinder vs time const l/w diff vib amps
%*11. plot radial frc on cylinder vs time const l/w diff vib amps
%*12. plot radial and z frc on cylinder vs time const vib diff l/w
%*13. plot topHook force vs time
%*14. hook force for constant l/w
%*15. plot sphericity for different l/w
%************************************************************
set(0, 'DefaultLineLineWidth', 2);
showFigs=[15];

lw=[]; nl=[]; npl=[]; vib=[]; N=[]; v=[];
props={lw nl npl vib N v};
I=2;%ind to be used for multiple different plot numbers below

clear usedD;inds=1;
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
[lws,nls,npls,vibs,vs]=separateVec([vertcat(usedD(:).pars)],1);
Ns=nls.*npls;

N=length(usedD);

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

db=[1,2,3,4,5,7,8,14,15];
% plot(nickD(db,1)',[nickD(db,2)-phiBS']','o-')
% xlabel('l/w');
% ylabel('\Delta\phi')
% figText(gcf,16)

%% 1. plot phi evo vs time for single trial
xx=1;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    
    t=usedD(I).t;       phi=usedD(I).phi;
    downsample(t,100);  downsample(phi,100);
    plot(t,phi);
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
    t=downsample(t,100);
    phi=downsample(phi,100);
    gui=downsample(gui,100);
    %find changes in guid
    
    plot(t,phi);
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
    
    clear out legz phiValsE phiValsM;
    for(i = 1:length(unLW))
        vibVals=[];
        unVib=unique(vibs(lws==unLW(i)));
        phiVals={};
        
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            for k=1:length(idz)
                [dphi,phi0,phif]=changeInPhiAmp(usedD(idz(k)),2);
                out(k)=dphi;
            end
            vibVals(j)=unVib(j);
            phiVals(j)={out};
        end
        phiValsM=cellfun(@(x) mean(x),phiVals);
        phiValsE=cellfun(@(x) std(x),phiVals);
        errorbar(vibVals,phiValsM,phiValsE,'linewidth',2);
        legz{i}=['lw=',num2str(unLW(i))];
    end
    
    legend(legz);
    xlabel('Vibration Amp (\circ)');
    ylabel('\Delta\phi');
    figText(gcf,18);
    clear 'legz' 'out' 'unLW' 'unVib' 'phiVals' 'vibVals' 'phiValsE' 'phiValsM' 'phiFVals';
end


%% 4. plot phif vs vibration amp
xx=4;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    
    %     unLW=unique([usedD.lw]);
    %     unVib=unique([usedD.vib]);
    
    %     lws=[0.7,0.7,0.7,0.6 0.6];
    %     vibs=[1,2,4,1,4];
    unLW=unique([lws]);
    %     unVib=unique([vibs]);
    
    clear out legz phiValsM phiValsE phiFVals;
    for(i = 1:length(unLW))
        vibVals=[];
        unVib=unique(vibs(lws==unLW(i)));
        phiFVals={};
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            for k=1:length(idz)
                [dphi,phi0,phif]=changeInPhiAmp(usedD(idz(k)),2);
                out(k,:)=[phi0,phif];
            end
            
            vibVals(j)=unVib(j);
            phiFVals(j)={out(:,2)};
        end
        phiValsM=cellfun(@(x) mean(x),phiFVals);
        phiValsE=cellfun(@(x) std(x),phiFVals);
        errorbar(vibVals,phiValsM,phiValsE,'linewidth',2);
        legz{i}=['lw=',num2str(unLW(i))];
    end
    
    legend(legz);
    xlabel('Vibration Amp (\circ)');
    ylabel('\phi_f');
    figText(gcf,18);
    clear 'legz' 'out' 'unLW' 'unVib' 'phiValsE' 'phiValsM' 'phiFVals';
end

%% 5. plot phi_i, phi_f for nicks data and for vib data
xx=5;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    
    %nicks stuff is defined above:
    %phi_i is dashed
    %phi_f is
    c=1;
    plot(lwBS,phiBS,'--o');
    legz{c}=['exp \phi_i']; c=c+1;
    plot(nickD(:,1),nickD(:,2),'-o');
    legz{c}=['exp \phi_f'];c=c+1;
    
    
    
    %     unLW=unique([usedD.lw]);
    %     unVib=unique([usedD.vib]);
    
    %     lws=[0.7,0.7,0.7,0.6 0.6];
    %     vibs=[1,2,4,1,4];
    unLW=unique([lws]);
    unVib=5; %i'll just use vib=5 for final vibration
    %(phi_i,phi_f,dphi) and error for lwsimE
    [lwsim,lwsimE]=deal(zeros(length(unLW),3));
    for(i = 1:length(unLW))
        vibVals=[];
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            for k=1:length(idz)
                [dphi,phi0,phif]=changeInPhiAmp(usedD(idz(k)),2);
                out(k,:)=[phi0,phif,dphi];
            end
            vibVals(j)=unVib(j);
            phiVals(j,:)={out};
        end
        lwsim(i,:)=cell2mat(cellfun(@(x) mean(x),phiVals,'UniformOutput',0));
        lwsimE(i,:)=cell2mat(cellfun(@(x) std(x),phiVals,'UniformOutput',0));
    end
    errorbar(unLW,lwsim(:,1),lwsimE(:,1),'--o','linewidth',2);
    legz{c}=['sim \phi_i vib=\pm',num2str(unVib),'\circ'];c=c+1;
    errorbar(unLW,lwsim(:,2),lwsimE(:,2),'-o','linewidth',2);
    legz{c}=['sim \phi_f vib=\pm',num2str(unVib),'\circ'];c=c+1;
    
    legend(legz);
    xlabel('l/w');
    ylabel('\phi');
    figText(gcf,18);
    clear 'legz' 'out' 'unLW' 'unVib' 'phiVals' 'vibVals' 'phiValsE' 'phiValsM' 'phiFVals';
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


%% 6. plot phi evolution vs vibration amp for many runs
xx=6;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    
    %     unLW=unique([usedD.lw]);
    %     unVib=unique([usedD.vib]);
    
    %     lws=[0.7,0.7,0.7,0.6 0.6];
    %     vibs=[1,2,4,1,4];
    %     unLW=[0.1];
    unVib=unique([vibs]);
    unLW=[0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
    
    clear out legz phiValsE phiValsM fs;
    fs=struct;
    for(i = 1:length(unLW))
        vibVals=[];
        %         unVib=unique(vibs(lws==unLW(i)));
        unVib=5;
        phiVals={};
        
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            for k=1:length(idz)
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                t=usedD(idz(k)).t;
                phi=usedD(idz(k)).phi;
                gui=usedD(idz(k)).gui;
                t=downsample(t,100);
                phi=downsample(phi,100);
                gui=downsample(gui,100);
                
                phi=phi(t>1);
                gui=gui(t>1);
                t=t(t>1);
                out(:,k)=phi;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            vibVals(j)=unVib(j);
            phiVals(j)={out};
        end
        y=cell2mat(cellfun(@(x) mean(x,2),phiVals,'UniformOutput',0));
        phiValsE=cell2mat(cellfun(@(x) std(x,0,2),phiVals,'UniformOutput',0));
        errorbar(t,y,phiValsE,'linewidth',2);
        legz{i}=['lw=',num2str(unLW(i))];
        %         set(gca,'xscale','log','yscale','log')
        plaw=fit(t,y,'a*x.^b+c');
        explaw=fit(t,y,'a*exp(-b*x)+c');
        fs(i).f=explaw;
        %         plot(t,ff.p1*t+ff.p2);
    end
    
    legend(legz);
    xlabel('t (s)');
    ylabel('\phi_f');
    figText(gcf,18);
    clear 'legz' 'out' 'unLW' 'unVib' 'phiVals' 'vibVals' 'phiValsE' 'phiValsM' 'phiFVals';
end
%% 7. plot <N> vs lw
xx=7;
if(showFigs(showFigs==xx))
    if(~exist('createPlane','file'))
        setupGeom3d
    end
    clear allnM allnE
    figure(xx);
    hold on;
    unVib=unique([vibs]);
    unLW=unique([usedD.lw]);
    for(i = 1:length(unLW))
        %         unVib=unique(vibs(lws==unLW(i)));
        unVib=5;
        phiVals={};
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            for k=1:length(idz)
                c=fileparts(usedD(idz(k)).fold);
                load(fullfile(c,'Nout.mat'));
                nM=mean(mean(Nout,2),3);
                nE=std(mean(Nout,2),0,3);
                
                allnM(i,:)=[nM(1) nM(end)];
                allnE(i,:)=[nE(1) nE(end)];
            end
        end
    end
    errorbar(unLW,allnM(:,1),allnE(:,1),'--','linewidth',2,'DisplayName','\langleN_i\rangle');
    errorbar(unLW,allnM(:,2),allnE(:,2),'linewidth',2,'DisplayName','\langleN_f\rangle');
    legend;
    xlabel('l/w')
    ylabel('\langleN\rangle');
    figText(gcf,18);
    
end
%% 8. plot pdf of motor torque vs <N>
xx=8;
if(showFigs(showFigs==xx))
    clear allnM allnE
    figure(xx);
    hold on;
    %     unVib=unique([vibs]);
    %     unLW=unique([usedD.lw]);
    unLW=[0.6];
    
    for(i = 1:length(unLW))
        % unVib=unique(vibs(lws==unLW(i)));
        unVib=[5];
        phiVals={};
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            idz=find(cond);
            if(idz)
                c=fileparts(usedD(idz(k)).fold);
                load(fullfile(c,'Nout.mat'));
                %                 ballOut=[x,y,z,ID,moveType,active,motTorque1,motTorque2]
                load(fullfile(c,'ballDat.mat'));
                TT=[];
                NN=[];
                for q=1:length(dat)
                    tors=sum(abs(dat(q).ballOut(:,7:8,:)),2);
                    %%%
                    tors=tors(:);
                    TT=[TT; tors];
                    nn=Nout(:,:,q)';
                    NN=[NN; nn(:)];
                    
                    %%%
                end
            else
                error('condition not found');
            end
            tval=1e-4;
            NN=NN(TT<tval);
            TT=TT(TT<tval);
            [nZ,XE,YE]=histcounts2(TT,NN,[100,18],'Normalization','pdf');
            XE=cumsum(diff(XE))-diff(XE)/2;
            YE=cumsum(diff(YE))-diff(YE)/2;
            [xx,yy]=meshgrid(XE,YE);
            
            [TY,TE]=discretize(TT,100);
            [NY,NE]=discretize(NN,[max(NN)-min(NN)]);
            contour(xx,yy,nZ')
        end
        
    end
    
    xlabel('Torque')
    ylabel('<N>');
    zlabel('count')
    figText(gcf,18);
end
%% 9. plot radial and z frc on cylinder vs time
xx=9;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    clear legz
    ds=7; % downsample
    ts= 2.4; %start time
    %                 dat(c).hookF=td(:,2);
    %             dat(c).prismSt=D(:,3);
    %             dat(c).wallF=D(:,4:5);
    t=usedD(I).t;
    wallR=usedD(I).wallF(:,1);
    wallF=usedD(I).wallF(:,2);
    prismSt=usedD(I).prismSt;
    
    prismSt=prismSt(t>ts);
    wallR=wallR(t>ts);
    wallF=wallF(t>ts);
    t=t(t>ts);
    
    t=downsample(t,ds);
    wallR=downsample(wallR,ds);
    wallF=downsample(wallF,ds);
    prismSt=downsample(prismSt,ds);
    
    %find changes in guid
    legz(1)=plot(t,wallR,'linewidth',2,'DisplayName','F_r');
    legz(2)=plot(t,wallF,'linewidth',2,'DisplayName','F_z');
    %find changes in guid
    gc=findChangesInGui(prismSt);
    gc(:,2)=gc(:,2)+4; %starts at -1, set lowest index to 3 for color
    axis tight
    dispName=["Before Lifting","Lift Begin","Lift End"];
    for(i=1:size(gc,1))
        legz(end+1)=plot([t(gc(i,1)),t(gc(i,1))],get(gca,'ylim'),'linewidth',3,...
            'color',guiCols(gc(i,2),:),'DisplayName',dispName(i));
    end
    
    xlabel('t (s)');
    ylabel('Force (N)');
    figText(gcf,18);
    legend(legz);
end

%% 10. plot radial and z frc on cylinder vs time const l/w diff vib amps
xx=10;
if(showFigs(showFigs==xx))
    clear r z
    figure(xx);
    hold on;
    legz=[];
    ds=7; % downsample
    ts= 2.4; %start time

    
    
    %     unVib=unique([vibs]);
    %     unLW=unique([usedD.lw]);
    unLW=[0.8];

    for(i = 1:length(unLW))
        unVib=unique(vibs(lws==unLW(i)));
        phiVals={};
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            idz=find(cond);
            for k=1:length(idz)
                
                t=usedD(idz(k)).t; %runs everytime even though I only need it once
                wallR=usedD(idz(k)).wallF(:,1);
                wallF=usedD(idz(k)).wallF(:,2);
                prismSt=usedD(idz(k)).prismSt;
                
                prismSt=prismSt(t>ts); %runs everytime even though I only need it once
                wallR=wallR(t>ts);
                wallF=wallF(t>ts);
                t=t(t>ts);
                
                t=downsample(t,ds);
                wallR=downsample(wallR,ds);
                wallF=downsample(wallF,ds);
                prismSt=downsample(prismSt,ds);
                
                r(k,:)=[wallR];
                z(k,:)=[wallF];
                
            end
            a=shadedErrorBar(t,mean(r),std(r),{'linewidth',2,'DisplayName',['F_r,v=',num2str(unVib(j))]},.5);
%             b=shadedErrorBar(t,mean(z),std(z),{'linewidth',2,'DisplayName',['F_z,v=',num2str(unVib(j))]},.5);
            legz(end+1)=a.mainLine;
%             legz(end+1)=b.mainLine;
            
            clear r z a
        end
    end
    

    
 %%%%%%find changes in guid%%%%%%%%%%%%%%%
    gc=findChangesInGui(prismSt);
    gc(:,2)=gc(:,2)+4; %starts at -1, set lowest index to 3 for color
    axis tight
    dispName=["Before Lifting","Lift Begin","Lift End"];
    for(i=1:size(gc,1))
        legz(end+1)=plot([t(gc(i,1)),t(gc(i,1))],get(gca,'ylim'),'linewidth',3,...
            'color',guiCols(gc(i,2),:),'DisplayName',dispName(i));
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
text(.2,.8,['l/w=',num2str(unLW)],'units','normalized');
xlabel('t (s)');
ylabel('Force (N)');
figText(gcf,18);
legend(legz)
    
end

%% 11. plot z frc on cylinder vs time const l/w diff vib amps
xx=11;
if(showFigs(showFigs==xx))
    clear r z
    figure(xx);
    hold on;
    legz=[];
    ds=7; % downsample
    ts= 2.4; %start time

    
    
    %     unVib=unique([vibs]);
    %     unLW=unique([usedD.lw]);
    unLW=[0.8];

    for(i = 1:length(unLW))
        unVib=unique(vibs(lws==unLW(i)));
        phiVals={};
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            idz=find(cond);
            for k=1:length(idz)
                
                t=usedD(idz(k)).t; %runs everytime even though I only need it once
                wallR=usedD(idz(k)).wallF(:,1);
                wallF=usedD(idz(k)).wallF(:,2);
                prismSt=usedD(idz(k)).prismSt;
                
                prismSt=prismSt(t>ts); %runs everytime even though I only need it once
                wallR=wallR(t>ts);
                wallF=wallF(t>ts);
                t=t(t>ts);
                
                t=downsample(t,ds);
                wallR=downsample(wallR,ds);
                wallF=downsample(wallF,ds);
                prismSt=downsample(prismSt,ds);
                
                r(k,:)=[wallR];
                z(k,:)=[wallF];
                
            end
%             a=shadedErrorBar(t,mean(r),std(r),{'linewidth',2,'DisplayName',['F_r,v=',num2str(unVib(j))]},.5);
            a=shadedErrorBar(t,mean(z),std(z),{'linewidth',2,'DisplayName',['F_z,v=',num2str(unVib(j))]},.5);
            legz(end+1)=a.mainLine;
%             legz(end+1)=b.mainLine;
            
            clear r z a
        end
    end
    

    
 %%%%%%find changes in guid%%%%%%%%%%%%%%%
    gc=findChangesInGui(prismSt);
    gc(:,2)=gc(:,2)+4; %starts at -1, set lowest index to 3 for color
    axis tight
    dispName=["Before Lifting","Lift Begin","Lift End"];
    for(i=1:size(gc,1))
        legz(end+1)=plot([t(gc(i,1)),t(gc(i,1))],get(gca,'ylim'),'linewidth',3,...
            'color',guiCols(gc(i,2),:),'DisplayName',dispName(i));
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
text(.2,.8,['l/w=',num2str(unLW)],'units','normalized');
xlabel('t (s)');
ylabel('Force (N)');
figText(gcf,18);
legend(legz)
    
end

%% 12. plot radial and z frc on cylinder vs time const vib diff l/w
xx=12;
if(showFigs(showFigs==xx))
    clear r z
    figure(xx);
    hold on;
    legz=[];
    ds=7; % downsample
    ts= 2.4; %start time

    
    
    %     unVib=unique([vibs]);
    %     unLW=unique([usedD.lw]);
    unVib=[2];
    clear allr allz allrM allrE allzM allzE 
    for(i = 1:length(unVib))
        unLW=unique(lws(vibs==unVib(i)));
        phiVals={};
        for(j=1:length(unLW))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            idz=find(cond);
            for k=1:length(idz)
                
                t=usedD(idz(k)).t; %runs everytime even though I only need it once
                wallR=usedD(idz(k)).wallF(:,1);
                wallF=usedD(idz(k)).wallF(:,2);
                prismSt=usedD(idz(k)).prismSt;
                
                prismSt=prismSt(t>ts); %runs everytime even though I only need it once
                wallR=wallR(t>ts);
                wallF=wallF(t>ts);
                t=t(t>ts);
                
                t=downsample(t,ds);
                wallR=downsample(wallR,ds);
                wallF=downsample(wallF,ds);
                prismSt=downsample(prismSt,ds);
                
                r(k,:)=[wallR];
                z(k,:)=[wallF];
                
            end
            legz(end+1)=errorbar(t,mean(r),std(r),'linewidth',2,'DisplayName',['F_r,l/w=',num2str(unLW(j))]);
            legz(end+1)=errorbar(t,mean(z),std(z),'linewidth',2,'DisplayName',['F_z,l/w=',num2str(unLW(j))]);
            clear r z
        end
    end
   
 %%%%%%find changes in guid%%%%%%%%%%%%%%%
    gc=findChangesInGui(prismSt);
    gc(:,2)=gc(:,2)+4; %starts at -1, set lowest index to 3 for color
    axis tight
    dispName=["Before Lifting","Lift Begin","Lift End"];
    for(i=1:size(gc,1))
        legz(end+1)=plot([t(gc(i,1)),t(gc(i,1))],get(gca,'ylim'),'linewidth',3,...
            'color',guiCols(gc(i,2),:),'DisplayName',dispName(i));
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
text(.2,.2,['vib=',num2str(unVib)],'units','normalized');
xlabel('t (s)');
ylabel('Force (N)');
figText(gcf,18);
    
    
end

%% 13. plot topHook force vs time
xx=13;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    clear legz
    ds=5; %dowsample amt
    ts=1; %start time
    %                 dat(c).hookF=td(:,2);
    %             dat(c).prismSt=D(:,3);
    %             dat(c).wallF=D(:,4:5);
    t=usedD(I).t;
    prismSt=usedD(I).prismSt;
    hookF=-usedD(I).hookF;
    
    prismSt=prismSt(t>ts);
    hookF=hookF(t>ts);
    t=t(t>ts);
    
    t=downsample(t,ds);
    hookF=downsample(hookF,ds);
    prismSt=downsample(prismSt,ds);
    
    %find changes in guid
    legz(1)=plot(t,hookF,'linewidth',2,'DisplayName','F_H');
    %find changes in guid
    gc=findChangesInGui(prismSt);
    gc(:,2)=gc(:,2)+4; %starts at -1, set lowest index to 3 for color
    axis tight
    dispName=["Before Lifting","Lift Begin","Lift End"];
    for(i=1:size(gc,1))
        legz(end+1)=plot([t(gc(i,1)),t(gc(i,1))],get(gca,'ylim'),'linewidth',3,...
            'color',guiCols(gc(i,2),:),'DisplayName',dispName(i));
    end
    
    xlabel('t (s)');
    ylabel('Force (N)');
    figText(gcf,18);
    legend(legz);
end
%% 14.hook force for constant l/w 
xx=14;
if(showFigs(showFigs==xx))
    clear F
    figure(xx);
    hold on;
    legz=[];
    ds=5; % downsample
    ts= 1; %start time

    
    
    %     unVib=unique([vibs]);
    %     unLW=unique([usedD.lw]);
    unLW=[0.8];

    for(i = 1:length(unLW))
        unVib=unique(vibs(lws==unLW(i)));
        phiVals={};
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            idz=find(cond);
            for k=1:length(idz)
                
                t=usedD(idz(k)).t; %runs everytime even though I only need it once
                hookF=-usedD(idz(k)).hookF;
                prismSt=usedD(idz(k)).prismSt;
                
                prismSt=prismSt(t>ts); %runs everytime even though I only need it once
                hookF=hookF(t>ts);
                t=t(t>ts);
                
                t=downsample(t,ds);
                hookF=downsample(hookF,ds);
                prismSt=downsample(prismSt,ds);
                
                F(k,:)=[hookF];
                
            end
            a=shadedErrorBar(t,mean(F),std(F),{'linewidth',2,'DisplayName',['F_H,v=',num2str(unVib(j))]},.5);
%             b=shadedErrorBar(t,mean(z),std(z),{'linewidth',2,'DisplayName',['F_z,v=',num2str(unVib(j))]},.5);
            legz(end+1)=a.mainLine;
%             legz(end+1)=b.mainLine;
            
            clear F a
        end
    end
    

    
 %%%%%%find changes in guid%%%%%%%%%%%%%%%
    gc=findChangesInGui(prismSt);
    gc(:,2)=gc(:,2)+4; %starts at -1, set lowest index to 3 for color
    axis tight
    dispName=["Before Lifting","Lift Begin","Lift End"];
    for(i=1:size(gc,1))
        legz(end+1)=plot([t(gc(i,1)),t(gc(i,1))],get(gca,'ylim'),'linewidth',3,...
            'color',guiCols(gc(i,2),:),'DisplayName',dispName(i));
    end
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
text(.2,.8,['l/w=',num2str(unLW)],'units','normalized');
xlabel('t (s)');
ylabel('Force (N)');
figText(gcf,18);
legend(legz)
    
end

%% 15.plot sphericity for different l/w
xx=15;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    unVib=unique([vibs]);
    unLW=unique([usedD.lw]);
    
    
    bucketRemoveTime=19;
    unVib=5;
    
     warning('remember to set bucketREmoveTime to correct value [14,19]')
    clear psiz psiM psiE
    n=15;%remove n markers with highest parwise distance
    for(i = 1:length(unLW))
        %         unVib=unique(vibs(lws==unLW(i)));
        
        phiVals={};
        for(j=1:length(unVib))
            cond=[lws==unLW(i)]&[vibs==(unVib(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            psiz=[];
            for k=1:length(idz)
                c=fileparts(usedD(idz(k)).fold);
                load(fullfile(c,'ballDat.mat'));
                verz=usedD(idz(k)).v;
                ballOut=dat(verz).ballOut;
                
                %%%%%%%%final vib frame data
%                 vibFrame=find(dat(verz).gui==4,1,'last');
                
                vibFrame=size(ballOut,3)-bucketRemoveTime; %2 remove bucket and wait 2s
                [x,y,z]=separateVec(ballOut(:,1:3,vibFrame),1);
                D=squareform(pdist(ballOut(:,1:3,vibFrame))); %square matrix
                dout=sum(D,1);
                [v,bottInds]=maxk(dout,n);
                x(bottInds)=[];y(bottInds)=[];z(bottInds)=[];
                shp=alphaShape(x,y,z,inf);
                vol=volume(shp);
                ar=surfaceArea(shp);
                psiz(k,1)=pi^(1/3)*(6*vol)^(2/3)/ar;
                
                %%%%%%%%%%%final without wall frame data
                [x,y,z]=separateVec(ballOut(:,1:3,end),1);
                D=squareform(pdist(ballOut(:,1:3,end))); %square matrix
                dout=sum(D,1);
                [v,bottInds]=maxk(dout,n);
                x(bottInds)=[];y(bottInds)=[];z(bottInds)=[];
                shp=alphaShape(x,y,z,inf);
                vol=volume(shp);
                ar=surfaceArea(shp);
                psiz(k,2)=pi^(1/3)*(6*vol)^(2/3)/ar;                  
            end
            
        end
        psiM(i,:)=mean(psiz,1);
        psiE(i,:)=std(psiz,0,1);
    end
    errorbar(unLW,psiM(:,1),psiE(:,1),'--','linewidth',2,'DisplayName','with walls');
    errorbar(unLW,psiM(:,2),psiE(:,2),'linewidth',2,'DisplayName','no walls');
    legend;
    xlabel('l/w')
    ylabel('\langleN\rangle');
    figText(gcf,18);
end