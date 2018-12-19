clear all
h=figure(1000);
guiCols=get(gca,'colororder');
guiCols(end+1,:)=[0 0 0];
close(h);
startInd=900;
fold=uigetdir('B:\SmartSimResults\12-5');
fname='ballDat.mat';
if ~exist(fullfile(fold,fname),'file')
    warning('reading in data since it does not exist otherwise');
    readBallPosDat(fold,startInd);
    error('run file again now but go into the folder inside the ballVary folder')
end
load(fullfile(fold,fname));

%************************************************************
%* Fig numbers:
%* 1. plot position of every staple at each checkpoint
%* 2. plot convex shape at each checkpoint
%* 3. plot sphericity at each checkpoint
%* 4. plot stress at each checkpoint
%* 5. scatter plot of positions colored by stress
%* 6. pairwise dist coloration
%* 7. plot <N> vs time
%* 8. plot <N> vs lw
%************************************************************
showFigs=[7];


I=3;%dat index
%% 1. plot position of every staple at each checkpoint
xx=1;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    ballOut=dat(I).ballOut;
    N=size(ballOut,3);
    t=dat(I).t;
    ms=40; %marker size
    for i = 1:N
        cc=zeros(size(ballOut,1),3);
        for(j=1:size(ballOut,1))
            cc(j,:)=guiCols(ballOut(j,5,i),:);
        end
        S=ms*ones(size(ballOut,1),1);%markerSize
        scatter3(ballOut(:,1,i),ballOut(:,2,i),ballOut(:,3,i),S,cc,'filled');
        axis([-0.03 0.03,-0.03 0.03,0, 0.05]);
        pause(.1)
        cla;
        
    end
end
%% 2. plot convex shape at each checkpoint
xx=2;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    %     for j=1:length(dat)
    ballOut=dat(I).ballOut;
    N=size(ballOut,3);
    n=15;%remove n markers with highest parwise distance
    for i = 1:N
        
        [x,y,z]=separateVec(ballOut(:,1:3,i),1);
        
        %%%%%%%%%%%%%%%%% remove farthest markers from calc
        D=squareform(pdist(ballOut(:,1:3,i))); %square matrix
        dout=sum(D,1);
        [v,bottInds]=maxk(dout,n);
        x(bottInds)=[];y(bottInds)=[];z(bottInds)=[];
        %%%%%%%%%%%%%%%%%
        
        
        [K,v]=convhulln([x,y,z]);
        shp=alphaShape(x,y,z,inf);
        vol=volume(shp);
        ar=surfaceArea(shp);
        %             trisurf(K,x,y,z,'Facecolor','cyan')
        [guis]=findChangesInGui(dat(I).gui');
        guis(1,:)=[0,1];
        guis(2:end,1)=guis(2:end,1)-2+startTime*10;
        %get gui from mode of
        plot(shp,'Facecolor',guiCols(dat(I).gui(i),:))
        %         axis([-0.03 0.03,-0.03 0.03,-0.01, 0]);
        title(['t=',num2str(dat(I).t(i)/1000),'s']);
        pause
        cla;
    end
    
    figText(gcf,18);
end

%% 3. plot sphericity at each checkpoint
xx=3;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    for j=1:length(dat)
        ballOut=dat(j).ballOut;
        N=size(ballOut,3);
        
        n=15;%remove n markers with highest parwise distance
        for i = 1:N
            [x,y,z]=separateVec(ballOut(:,1:3,i),1);
            %%%%%%%%%%%%%%%%%remove farthest markers from calc
            D=squareform(pdist(ballOut(:,1:3,i))); %square matrix
            dout=sum(D,1);
            [v,bottInds]=maxk(dout,n);
            x(bottInds)=[];y(bottInds)=[];z(bottInds)=[];
            %%%%%%%%%%%%%%%%
            
            shp=alphaShape(x,y,z,inf);
            vol=volume(shp);
            ar=surfaceArea(shp);
            %         trisurf(K,x,y,z,'Facecolor','cyan')
            %         axis([-0.03 0.03,-0.03 0.03,0, 0.05]);
            %         cla;
            psi(i)=pi^(1/3)*(6*vol)^(2/3)/ar;
            cc=guiCols(ballOut(j,5,1),:);
            
        end
        allPsi{j}=psi;
        allT{j}=dat(j).t;
        plot(dat(j).t/1000,psi,'-','linewidth',2);
    end
    
    [guis]=findChangesInGui(dat(1).gui');
    guis(1,:)=[0,1];
    guis(2:end,1)=guis(2:end,1)-2+startTime*10;
    for k=1:size(guis,1)
        plot(guis(k,1)/10*ones(2,1),ylim,'linewidth',2,'color',guiCols(guis(k,2),:))
    end
    xlabel('t (s)');
    ylabel('\Psi');
    figText(gcf,18);
end
%% 4. plot mot torque at each checkpoint
xx=4;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    for j=1:length(dat)
        ballOut=dat(j).ballOut;
        Nfr=size(ballOut,3);
        t=dat(j).t/1000;
        tor(:,j)=dat(j).totTorque/dat(j).N;
        %         plot(dat(j).t/1000,dat(j).totTorque/dat(j).N,'-o');
    end
    errorbar(t,mean(tor,2),std(tor,0,2),'-o','linewidth',2);
    [guis]=findChangesInGui(dat(1).gui');
    guis(1,:)=[0,1];
    guis(2:end,1)=guis(2:end,1)-2+startTime*10;
    axis tight
    for k=1:size(guis,1)
        plot(guis(k,1)/10*ones(2,1),ylim,'linewidth',2,'color',guiCols(guis(k,2),:))
        
    end
    xlabel('t (s)');
    ylabel('\langle\tau\rangle_N (Nm)');
    figText(gcf,18);
    
end
%% 5. scatter plot of positions colored by stress
xx=5;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    ballOut=dat(I).ballOut;
    fr=size(ballOut,3);
    t=dat(I).t;
    colormap(brewermap(size(ballOut,1),'RdYlBu'));
    ms=40; %marker size
    while(1)
        for i = 1:fr
            S=ms*ones(size(ballOut,1),1);%markerSize
            tor=sum(abs(ballOut(:,7:8,21)),2);
            scatter3(ballOut(:,1,i),ballOut(:,2,i),ballOut(:,3,i),S,tor,'filled');
            %         axis([-0.03 0.03,-0.03 0.03,0, 0.05]);
            title(['t=',num2str(dat(I).t(i)/1000),'s']);
            pause(0.1);
            axis([-.04 .04 -.04 .04 0 0.03])
            cla;
            
        end
    end
end

%% 6. pairwise dist coloration
xx=6;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    ballOut=dat(I).ballOut;
    fr=size(ballOut,3);
    colormap(brewermap(size(ballOut,1),'RdYlBu'));
    ms=40; %marker size
    while(1)
        for i = 1:fr
            S=ms*ones(size(ballOut,1),1);%markerSize
            D=squareform(pdist(ballOut(:,1:3,i))); %square matrix
            dout=sum(D,1);
            scatter3(ballOut(:,1,i),ballOut(:,2,i),ballOut(:,3,i),S,dout,'filled');
            title(['t=',num2str(dat(I).t(i)/1000),'s']);
            pause(0.1);
            %             axis([-.04 .04 -.04 .04 0 0.03])
            cla;
        end
    end
end

%% 7. plot <N> vs time
xx=7;
if(showFigs(showFigs==xx))
    if(~exist('createPlane','file'))
        setupGeom3d
    end
    figure(xx);
    hold on;
    ff=fullfile(fold,'Nout.mat');
    if exist(ff,'file')
        load(ff);
    else
        wid='MATLAB:inpolygon:ModelingWorldLower';
        warning('off',wid);
        error('Nout file doesn''t exist');
        for(i=1:length(dat))
            tic
                Nout(:,:,i)=generatePackingFromSimDatPLANE_mex(dat(i),0);
            toc
        end
        save(ff,'Nout','dat');
    end
    nM=mean(mean(Nout,2),3);
    nE=std(mean(Nout,2),0,3);  
   
    t=dat(1).t'/1000;   y=nM;
     
    explaw=fit(t,y,'a*exp(-b*x)+c');
    errorbar(t,y,nE,'linewidth',2);
    plot(explaw)

    xlabel('t (s)')
    ylabel('\langle N\rangle');
    figText(gcf,18);
    set(gca,'yscale','log');
end
%% 8. plot <N> vs lw
xx=8;
if(showFigs(showFigs==xx))
    if(~exist('createPlane','file'))
        setupGeom3d
    end
    figure(xx);
    hold on;
    [c]=fileparts(fold);
    outerFolds=dir2(c,'folders');
    
    ff=fullfile(fold,'Nout.mat');
    folder    
    xlabel('t (s)')
    ylabel('\langle N\rangle');
    figText(gcf,18);
    set(gca,'yscale','log');
end