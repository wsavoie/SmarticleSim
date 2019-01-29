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
%*10/11. plot z/r frc on cylinder vs time const l/w=1 or const vib=0
%*12. plot topHook force vs time
%*13. hook force for constant l/w=1 or const vib=0
%*14. plot sphericity for different l/w
%*15. plot sphericity for different l/w outerpoints
%*16. plot p(F) vs n for many lw
%*17. plot p(F) vs n for single lw at diff times
%************************************************************
set(0, 'DefaultLineLineWidth', 2);
showFigs=[12];
% showFigs=[1 7 16 17];
unViba=0;
lw=[]; nl=[]; npl=[]; vib=[]; N=[]; v=[];
props={lw nl npl vib N v};
I=1;%ind to be used for multiple different plot numbers below

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
    samp=10;
    t=usedD(I).t;       phi=usedD(I).phi;
    phi=movmean(phi,samp);
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
    legz=[];
    %nicks stuff is defined above:
    %phi_i is dashed
    legz(end+1)=plot(lwBS,phiBS,'--o','DisplayName',['exp \phi_i']);
    legz(end+1)=plot(nickD(:,1),nickD(:,2),'-o','DisplayName',['exp \phi_f']);
    
    varyz=vibs;
    constz=lws;
    unConstz=[unique(constz)]; %constant value you want to use
    constName='l/w';
    constUnits='';
    varyName='v';
    varyUnits='\circ';
    
    
    
    
    
    %     unLW=unique([lws]);
    %     unVib=5; %i'll just use vib=5 for final vibration
    %(phi_i,phi_f,dphi) and error for lwsimE
    [phiSim,phiSimE]=deal(zeros(length(unConstz),3));
    for(i = 1:length(unConstz))
        vibVals=[];
        unVaryz=unique(varyz(constz==unConstz(i)));
        for(j=1:length(unVaryz))
            cond=[constz==unConstz(i)]&[varyz==(unVaryz(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            for k=1:length(idz)
                [dphi,phi0,phif]=changeInPhiAmp(usedD(idz(k)),2);
                out(k,:)=[phi0,phif,dphi];
            end
            vibVals(j)=unVaryz(j);
            phiVals(j,:)=out;
        end
        phiSim(i,:)=mean(phiVals,1);
        phiSimE(i,:)=std(phiVals,0,1);
        %         phiSim(i,:)=cell2mat(cellfun(@(x) mean(x),phiVals,'UniformOutput',0));
        %         phiSimE(i,:)=cell2mat(cellfun(@(x) std(x),phiVals,'UniformOutput',0));
        
    end
    for j=1:length(unVaryz)
        legz(end+1)=errorbar(unConstz,phiSim(:,1),phiSimE(:,1),'--o','linewidth',2,'DisplayName',['sim \phi_i',varyName,'=\pm',num2str(unVaryz),varyUnits]);
        legz(end+1)=errorbar(unConstz,phiSim(:,2),phiSimE(:,2),'-o','linewidth',2,'DisplayName',['sim \phi_f ',varyName,'=\pm',num2str(unVaryz),varyUnits]);
    end
    legend(legz);
    xlabel(constName);
    ylabel('\phi');
    figText(gcf,18);
    clear 'legz' 'out' 'unVary' 'unConstz' 'phiVals' 'vibVals' 'phiValsE' 'phiValsM' 'phiFVals';
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
%         unViba=0;
        phiVals={};
        for(j=1:length(unViba))
            cond=[lws==unLW(i)]&[vibs==(unViba(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            for k=1:length(idz)
                c=fileparts(usedD(idz(k)).fold);
                load(fullfile(c,'Nout.mat'));
                nM=mean(mean(Nout,2),3);
                nE=std(mean(Nout,2),0,3);
                [nf,kk]= max(nM);
                
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
    axis([0 1.3,0 12])
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
    sWeight=9.8*prod(usedD(I).smartSize(1:3))*(2*usedD(I).smartSize(4)+usedD(I).smartSize(5));
    wallR=usedD(I).wallF(:,1)/(usedD(I).N*sWeight);
    wallF=usedD(I).wallF(:,2)/(usedD(I).N*sWeight);
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
    ylabel('Force/(s_mg*N)');
    figText(gcf,18);
    legend(legz);
end

%% 10/11. plot z/r frc on cylinder vs time const l/w=1 or const vib=0
xx=10;
if(showFigs(showFigs==xx))
    clear r z
    figure(xx);
    hold on;
    figure(xx+1);
    hold on;
    legz=[];
    legza=[];
    legzb=[];
    ds=50; % downsample
    ts= 2.4; %start time
    
    %varyType, vib=1, vary l/w=0
    
    varyType=0;
    if(varyType)
        varyz=vibs;
        constz=lws;
        unConstz=[0.7]; %constant value you want to use
        constName='l/w';
        constUnits='';
        varyName='v';
        varyUnits='\circ';
    else
        varyz=lws;
        constz=vibs;
        unConstz=[0]; %constant value you want to use
        constName='v';
        constUnits='\circ';
        varyName='l/w';
        varyUnits='';
    end
    
    for(i = 1:length(unConstz))
        unVaryz=unique(varyz(constz==unConstz(i)));
        for(j=1:length(unVaryz))
            cond=[constz==unConstz(i)]&[varyz==(unVaryz(j))];
            idz=find(cond);
            for k=1:length(idz)
                %                 smartSize=[rho,t1,t2,l,w];
                sWeight=9.8*prod(usedD(idz(k)).smartSize(1:3))*(2*usedD(idz(k)).smartSize(4)+usedD(idz(k)).smartSize(5));
                t=usedD(idz(k)).t; %runs everytime even though I only need it once
                wallR=usedD(idz(k)).wallF(:,1)/(usedD(idz(k)).N*sWeight);
                wallF=usedD(idz(k)).wallF(:,2)/(usedD(idz(k)).N*sWeight);
                prismSt=usedD(idz(k)).prismSt;
                
                prismSt=prismSt(t>ts); %runs everytime even though I only need it once
                wallR=wallR(t>ts);
                wallF=wallF(t>ts);
                t=t(t>ts);
                
                %                 t=downsample(t,ds);
                %                 wallR=downsample(wallR,ds);
                %                 wallF=downsample(wallF,ds);
                %                 prismSt=downsample(prismSt,ds);
                wallR=movmean(wallR,ds);
                wallF=movmean(wallF,ds);
                r(k,:)=[wallR];
                z(k,:)=[wallF];
                
            end
            %             a=shadedErrorBar(t,mean(r),std(r),{'linewidth',2,'DisplayName',['F_r,v=',num2str(unVib(j))]},.5);
            figure(xx)
            a=shadedErrorBar(t,mean(z,1),std(z,0,1),{'linewidth',2,'DisplayName',['F_z, ',varyName,'=',num2str(unVaryz(j)),varyUnits]},.5);
            figure(xx+1)
            b=shadedErrorBar(t,mean(r,1),std(r,0,1),{'linewidth',2,'DisplayName',['F_r, ',varyName,'=',num2str(unVaryz(j)),varyUnits]},.5);
            
            legza(end+1)=a.mainLine;
            legzb(end+1)=b.mainLine;
            clear r z a b
        end
    end
    
    
    for(qq=xx:xx+1)
        figure(qq);
        if(qq==xx)
            legz=legza;
        else
            legz=legzb;
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
        text(.2,.8,[constName,'=',num2str(unConstz),constUnits],'units','normalized')
        xlabel('t (s)');
        ylabel('Force/(s_mg*N)');
        figText(gcf,18);
        legend(legz);
    end
    
end

%% 12. plot topHook force vs time
xx=12;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    clear legz
    ds=10; %dowsample amt
    ts=1; %start time
    legz=[];
    %                 dat(c).hookF=td(:,2);
    %             dat(c).prismSt=D(:,3);
    %             dat(c).wallF=D(:,4:5);
    t=usedD(I).t;
    prismSt=usedD(I).prismSt;
    sWeight=9.8*prod(usedD(I).smartSize(1:3))*(2*usedD(I).smartSize(4)+usedD(I).smartSize(5));
    hookF=-usedD(I).hookF/(usedD(I).N*sWeight);
    
    prismSt=prismSt(t>ts);
    hookF=hookF(t>ts);
    t0=t;
    t=t(t>ts);
    
%     t=downsample(t,ds);
    hookF=movmean(hookF,ds)/usedD(I).N;%scaled force by numb of smarts
%     prismSt=movmean(prismSt,ds);
    
    %find changes in guid
    legz(1)=plot(t,hookF,'linewidth',2,'DisplayName','F_H');
    %find changes in guid
    gc=findChangesInGui(prismSt);
    gc(:,2)=gc(:,2)+4; %starts at -1, set lowest index to 3 for color
    axis tight
    dispName=["Before Lifting","Lift Begin","Lift End"];
    for(i=1:size(gc,1))
        legz(end+1)=plot(t(gc(i,1))*[1,1],get(gca,'ylim'),'linewidth',3,...
            'color',guiCols(gc(i,2),:),'DisplayName',dispName(i));
    end
    
    if isfield(usedD(I),'bucketExist')
        bucketPt=find(usedD(I).bucketExist(t0>ts)==0,1,'first');
        if(bucketPt)
            legz(end+1)=plot(t(bucketPt)*[1,1],get(gca,'ylim'),'linewidth',3,...
                'DisplayName','bucket removal','color','k');
        end
    end
    
    xlabel('t (s)');
    ylabel('Force/(s_mg*N)');
    figText(gcf,18);
    legend(legz);
    
    %%%%
    xlim([4.49,7.5])
end
%% 13.hook force for constant l/w=1 or const vib=0
xx=13;
if(showFigs(showFigs==xx))
    clear F
    figure(xx);
    hold on;
    legz=[];
    ds=50; % downsample
    ts= 1; %start time
    
    %varyType, vib=1, vary l/w=0
    
    varyType=1;
    if(varyType)
        varyz=vibs;
        constz=lws;
        unConstz=[0.6]; %constant value you want to use
        constName='l/w';
        constUnits='';
        varyName='v';
        varyUnits='\circ';
    else
        varyz=lws;
        constz=vibs;
        unConstz=[0]; %constant value you want to use
        constName='v';
        constUnits='\circ';
        varyName='l/w';
        varyUnits='';
    end
    
    for(i = 1:length(unConstz))
        unVaryz=unique(varyz(constz==unConstz(i)));
        phiVals={};
        for(j=1:length(unVaryz))
            cond=[constz==unConstz(i)]&[varyz==(unVaryz(j))];
            idz=find(cond);
            for k=1:length(idz)
                
                sWeight=9.8*prod(usedD(idz(k)).smartSize(1:3))*(2*usedD(idz(k)).smartSize(4)+usedD(idz(k)).smartSize(5));
                t=usedD(idz(k)).t; %runs everytime even though I only need it once
                hookF=-usedD(idz(k)).hookF/(usedD(idz(k)).N*sWeight);
                prismSt=usedD(idz(k)).prismSt;
                
                prismSt=prismSt(t>ts); %runs everytime even though I only need it once
                hookF=hookF(t>ts);
                t=t(t>ts);
                
                %                 t=downsample(t,ds);
                %                 hookF=downsample(hookF,ds);
                %                 prismSt=downsample(prismSt,ds);
                hookF=movmean(hookF,ds);
                F(k,:)=[hookF];
                
            end
            a=shadedErrorBar(t,mean(F,1),std(F,0,1),{'linewidth',2,'DisplayName',['F_H, ',varyName,'=',num2str(unVaryz(j)),varyUnits]},.5);
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
    text(.2,.8,[constName,'=',num2str(unConstz),constUnits],'units','normalized')
    xlabel('t (s)');
    ylabel('Force/(s_mg*N)');
    figText(gcf,18);
    legend(legz)
    
end

%% 14.plot sphericity for different l/w
xx=14;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    
    
    clear psiz psiM psiE
    n=15;%remove n markers with highest parwise distance
    legz=[];
    varyType=1;
    if(varyType)
        varyz=vibs;
        constz=lws;
        unConstz=[unique(constz)]; %constant value you want to use
        constName='l/w';
        constUnits='';
        varyName='v';
        varyUnits='\circ';
    else
        varyz=lws;
        constz=vibs;
        unConstz=[unique(constz)]; %constant value you want to use
        constName='v';
        constUnits='\circ';
        varyName='l/w';
        varyUnits='';
    end
    
    for(i = 1:length(unConstz))
        unVaryz=unique(varyz(constz==unConstz(i)));
        for(j=1:length(unVaryz))
            cond=[constz==unConstz(i)]&[varyz==(unVaryz(j))];
            idz=find(cond);
            
            for k=1:length(idz)
                bucketRad=usedD(idz(k)).smartSize(5)*2; %2*w_smart
                c=fileparts(usedD(idz(k)).fold);
                load(fullfile(c,'ballDat.mat'));
                verz=usedD(idz(k)).v;
                
                verzID=find([dat(:).v]==verz);
                ballOut=dat(verzID).ballOut;
                if isfield(usedD(idz(k)),'bucketExist')
                    %find time bucket was removed
                    tt=usedD(idz(k)).t(find(usedD(idz(k)).bucketExist==0,1,'first'));
                    %value needs to be rounded and multiplied to fit index
                    %type of frame
                    tRem=round(tt,1)*10;
                    %                     shave of indices for startTime
                    BucketRemoveFrame=tRem-startTime*10;
                    %                     error(['first time running, check if value is expected, vibFrame starts at: ',num2str(BucketRemoveFrame)])
                else %using older file
                    warning('remember to set bucketREmoveTime to correct value [14,19]')
                    bucketRemoveTime=19;
                    BucketRemoveFrame=size(ballOut,3)-bucketRemoveTime; %2 remove bucket and wait 2s
                end
                %%%%%%%%last frame before bucket removal
                %                 [x,y,z]=separateVec(ballOut(:,1:3,BucketRemoveFrame),1);
                hookPos=0;
                if isfield(usedD(idz(k)),'hookPos')
                    hookPos=usedD(idz(k)).hookPos(usedD(idz(k)).t==tt);
                end
                [smartPosOut]=RemoveSmartsFromBucket(ballOut(:,1:3,BucketRemoveFrame),n,hookPos,bucketRad);
                psiz(k,1)=calcSphericity(smartPosOut);
                %%%%%%%%%%%final without wall frame data
                %                 [x,y,z]=separateVec(ballOut(:,1:3,end),1);
                hookPos=0;
                if isfield(usedD(idz(k)),hookPos)
                    hookPos=usedD(idz(k)).hookPos(end);
                end
                [smartPosOut]=RemoveSmartsFromBucket(ballOut(:,1:3,end),n,hookPos,bucketRad);
                psiz(k,2)=calcSphericity(smartPosOut);
            end
            
        end
        psiM(i,:)=mean(psiz,1);
        psiE(i,:)=std(psiz,0,1);
    end
    legz(end+1)=errorbar(unConstz,psiM(:,1),psiE(:,1),'--','linewidth',2,'DisplayName',['with walls, ', varyName,'=',num2str(unVaryz(j)),varyUnits]);
    legz(end+1)=errorbar(unConstz,psiM(:,2),psiE(:,2),'linewidth',2,'DisplayName',['no walls, ', varyName,'=',num2str(unVaryz(j)),varyUnits]);
    legend(legz);
    cc=[];
    if(constUnits)
        cc=[' (',constUnits,')'];
    end
    xlabel([constName,cc])
    ylabel('Psi');
    figText(gcf,18);
end
%% 15.plot sphericity for different l/w outerpoints
xx=15;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    
    
    clear psiz psiM psiE
    n=15;%remove n markers with highest parwise distance
    legz=[];
    varyType=1;
    if(varyType)
        varyz=vibs;
        constz=lws;
        unConstz=[unique(constz)]; %constant value you want to use
        constName='l/w';
        constUnits='';
        varyName='v';
        varyUnits='\circ';
    else
        varyz=lws;
        constz=vibs;
        unConstz=[unique(constz)]; %constant value you want to use
        constName='v';
        constUnits='\circ';
        varyName='l/w';
        varyUnits='';
    end
    if ~exist(fullfile(fold,'posOut.mat'),'file')
        getAllSmartPoints(dat);
    end
    load(fullfile(fold,'posOut.mat'))
    for(i = 1:length(unConstz))
        unVaryz=unique(varyz(constz==unConstz(i)));
        for(j=1:length(unVaryz))
            cond=[constz==unConstz(i)]&[varyz==(unVaryz(j))];
            idz=find(cond);
            for k=1:length(idz)
                bucketRad=usedD(idz(k)).smartSize(5)*2; %2*w_smart
                c=fileparts(usedD(idz(k)).fold);
                load(fullfile(c,'ballDat.mat'));
                verz=usedD(idz(k)).v;
                verzID=find([dat(:).v]==verz);
                %                 ballOut=posOut(idz(k)).ballOut;
                ballOut=dat(verzID).ballOut;
                smartPos=posOut(idz(k)).smartPos;
                if isfield(usedD(idz(k)),'bucketExist')
                    %find time bucket was removed
                    tt=usedD(idz(k)).t(find(usedD(idz(k)).bucketExist==0,1,'first'));
                    %value needs to be rounded and multiplied to fit index
                    %type of frame
                    tRem=round(tt,1)*10;
                    %                     shave of indices for startTime
                    BucketRemoveFrame=tRem-startTime*10;
                    %                     error(['first time running, check if value is expected, vibFrame starts at: ',num2str(vibFrame)])
                else %using older file
                    warning('remember to set bucketREmoveTime to correct value [14,19]')
                    bucketRemoveTime=19;
                    BucketRemoveFrame=size(ballOut,3)-bucketRemoveTime; %2 remove bucket and wait 2s
                end
                %%%%%%%%final vib frame data
                hookPos=0;
                if isfield(usedD(idz(k)),'hookPos')
                    %                     hookPos=usedD(idz(k)).hookPos(tt);
                    hookPos=usedD(idz(k)).hookPos(usedD(idz(k)).t==tt);
                end
                
                %resize smartPos
                aa=squeeze(smartPos(BucketRemoveFrame,:,:,:));
                aa=reshape(aa,size(aa,1)*size(aa,2),3);
                [smartPosOut]=RemoveSmartsFromBucket(aa,n,hookPos,bucketRad);
                psiz(k,1)=calcSphericity(smartPosOut);
                %%%%%%%%%%%final without wall frame data
                %                 [x,y,z]=separateVec(ballOut(:,1:3,end),1);
                hookPos=0;
                if isfield(usedD(idz(k)),hookPos)
                    hookPos=usedD(idz(k)).hookPos(end);
                end
                aa=squeeze(smartPos(end,:,:,:));
                aa=reshape(aa,size(aa,1)*size(aa,2),3);
                [smartPosOut]=RemoveSmartsFromBucket(aa,n,hookPos,bucketRad);
                psiz(k,2)=calcSphericity(smartPosOut);
                clear aa
            end
            
        end
        psiM(i,:)=mean(psiz,1);
        psiE(i,:)=std(psiz,0,1);
    end
    legz(end+1)=errorbar(unConstz,psiM(:,1),psiE(:,1),'--','linewidth',2,'DisplayName',['with walls, ', varyName,'=',num2str(unVaryz(j)),varyUnits]);
    legz(end+1)=errorbar(unConstz,psiM(:,2),psiE(:,2),'linewidth',2,'DisplayName',['no walls, ', varyName,'=',num2str(unVaryz(j)),varyUnits]);
    legend(legz);
    cc=[];
    if(constUnits)
        cc=[' (',constUnits,')'];
    end
    xlabel([constName,cc])
    ylabel('\Psi');
    figText(gcf,18);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 16.  plot p(F) vs n
xx=16;
if(showFigs(showFigs==xx))
    if(~exist('createPlane','file'))
        setupGeom3d
    end
    clear nfinal nTot
    figure(xx);
    hold on;
    unVib=unique([vibs]);
    unLW=unique([usedD.lw]);
    edges=[0:1:30];
%     unViba=5;
    nfinal=zeros(length(unLW),length(diff(edges)),length(unViba));
    colz=colormap(brewermap(length(unLW),'RdYlBu'));
    for(i = 1:length(unLW))
        %         unVib=unique(vibs(lws==unLW(i)));
        phiVals={};
        timeI=40;
        
        for(j=1:length(unViba))
            cond=[lws==unLW(i)]&[vibs==(unViba(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            nTot=[];
            for k=1:length(idz)
                c=fileparts(usedD(idz(k)).fold);
                load(fullfile(c,'Nout.mat'));
                
                nTot=horzcat(nTot,Nout(timeI,:));
            end
            nfinal(i,:,j)=histcounts(nTot,'binedges',edges,'Normalization','probability');
        end
        plot(cumsum(diff(edges)),nfinal(i,:,1)+0.1*(length(unLW)-i),'color',colz(i,:));
    end

    xlabel('N')
    ylabel('P(N) + constant');
    ar=annotation('arrow','Color','k','X',[.8 , .8],'Y',[.8,.5]);
    text(ar.X(1)+0.05,.9,['l/w'],'units','Normalized');

    figText(gcf,18);
    
    
end
%% 17. plot p(F) vs n for single lw at diff times
xx=17;
if(showFigs(showFigs==xx))
    if(~exist('createPlane','file'))
        setupGeom3d
    end
    clear nfinal nTot
    figure(xx);
    hold on;
    unVib=unique([vibs]);
    unLW=unique([usedD.lw]);
    edges=[0:1:14];
%     unViba=5;
    unLW=0.4;
    
    for(i = 1:length(unLW))
        for(j=1:length(unViba))
            cond=[lws==unLW(i)]&[vibs==(unViba(j))];
            %             vals=length(cond(cond==1));
            idz=find(cond);
            nTot=[];
            for k=1:length(idz)
                c=fileparts(usedD(idz(k)).fold);
                load(fullfile(c,'Nout.mat'));
                nTot=horzcat(nTot,Nout);
            end
        end
    end
        colz=colormap(brewermap(size(Nout,1),'RdYlBu'));
        maxline=zeros(size(nTot,1),2);
    for(i=1:size(nTot,1))
        nfinal(i,:)=histcounts(nTot(i,:),'binedges',edges,'Normalization','probability');
        plot(cumsum(diff(edges)),nfinal(i,:)+0.1*i,'color',colz(i,:));
        [maxline(i,2),maxline(i,1)]=max(nfinal(i,:)+0.1*i);
    end
    plot(maxline(:,1),maxline(:,2),'k','linewidth',2);
    xlabel('N')
    ylabel('P(N) + constant');
    ar=annotation('arrow','Color','k','X',[.8 , .8],'Y',[.4,.8]);
    text(ar.X(1)+0.05,.3,['time'],'units','Normalized');

    figText(gcf,18);
    
    
end