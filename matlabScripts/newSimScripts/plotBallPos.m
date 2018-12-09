clear all
h=figure(1000);
guiCols=get(gca,'colororder');
guiCols(end+1,:)=[0 0 0];
close(h);
startInd=500;
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
%************************************************************
showFigs=[2 3];


I=1;%dat index
%% 1. plot position of every staple at each checkpoint
xx=1;
if(showFigs(showFigs==xx))
    figure(xx);
    hold on;
    ballOut=dat(I).ballOut;
    N=size(ballOut,3);
    t=dat(I).t;
    for i = 1:N
        cc=zeros(size(ballOut,1),3);
        for(j=1:size(ballOut,1))
            cc(j,:)=guiCols(ballOut(j,5,i),:);
        end
        S=20*ones(size(ballOut,1),1);%markerSize
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
    for i = 1:N
        
        [x,y,z]=separateVec(ballOut(:,1:3,i),1);
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
        axis([-0.03 0.03,-0.03 0.03,0, 0.05]);
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
        for i = 1:N
            [x,y,z]=separateVec(ballOut(:,1:3,i),1);
            
            %remove certain elements based on position
            Rlowest=3;
            [v,bottInds]=mink(z,Rlowest);
            x(bottInds)=[];y(bottInds)=[];z(bottInds)=[];
            
            
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
        
        %             vals=; %put into units of seconds
        
        plot(dat(j).t/1000,psi,'--','linewidth',2);
        
        
        
        %     cc(j,:)=guiCols(ballOut(j,5,1),:);
        %ballOut=[x,y,z,ID,moveType,active]
        
        
        
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
