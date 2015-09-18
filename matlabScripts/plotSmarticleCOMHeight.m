% directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests');
data = importdata('D:\SimResults\Chrono\SmarticleU\tests\9-17-15--9-18-15\com changing shape\r1\PostProcess\volumeFraction.txt');

time        = data(:,1);
smartcount  = data(:,2);
volfrac     = data(:,3);
zmax        = data(:,4);
zcom        = data(:,5);
meanOT      = data(:,6);
guid        = data(:,7);


% Time << CountInside << volumeFraction << zMax << zCom <<GUID
cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7]};
plotNames = {'Smarticles','Gait','U-Shape','Straight','Tetris','Vib at \circ','Vib Angle'};
plotTypes = [1];
hold on;

gc =[0; find(diff(guid))];
gc=gc+1;
plot(time,zcom);
val=[data(gc,6)]+1;%added value can be zero (global) and matrices are 1 started
sp1=[plot(time,zcom)]; %starts variable for the legend
gg = [gc, val];
ggOld = gg;
gg= sortrows(gg,2); 
for i=1:size(gg,1)
    if i~=1
        if gg(i,2)==gg(i-1,2) %if true will mean a repeat in color on legend
            plot([data(gg(i,1),1),data(gg(i,1),1)],[min(zcom),max(zcom)],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)
            text(data(gg(i,1),1),max(zcom),plotNames((gg(i,2)+1)))
            continue;
        else
            sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(zcom),max(zcom)],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
            text(data(gg(i,1),1),max(zcom),plotNames((gg(i,2)+1)))
        end
    else
        sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(zcom),max(zcom)],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
        text(data(gg(i,1),1),max(zcom),plotNames((gg(i,2)+1)))
    end
plotTypes=[plotTypes, (gg(i,2)+1)];
end


data2 = importdata('D:\SimResults\Chrono\SmarticleU\tests\com changing shape\r2\PostProcess\volumeFraction.txt');
data3 = importdata('D:\SimResults\Chrono\SmarticleU\tests\com changing shape\r3\PostProcess\volumeFraction.txt');
% totData= [zcom,data2(:,5)];
totData= [zcom,data2(:,5),data3(:,5)];

% plot(time,mean(totData,2));
err = std(totData,0,2);
shadedErrorBar(time,mean(totData,2),err,{'color',[12/255,100/255,2/255]});

figure(2);
hold on;
plot(time,zcom);
plot(data2(:,1),data2(:,5));
plot(data3(:,1),data3(:,5));

figure(2);
hold on;
plot(time,volfrac);
plot(data2(:,1),data2(:,3));
% plot(data3(:,1),data3(:,5));