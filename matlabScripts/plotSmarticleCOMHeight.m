% directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests');
% data = importdata('D:\SimResults\Chrono\SmarticleU\tests\9-17-15--9-18-15\com changing shape\r1\PostProcess\volumeFraction.txt');
% data = importdata('D:\SimResults\Chrono\SmarticleU\tests\PostProcess\volumeFraction.txt');
% data = importdata('D:\SimResults\Chrono\SmarticleU\tests\9-20\PostProcess\volumeFraction.txt');
% data = importdata('\\centos\secured\shared_data\diffWidth2VolFracFilling\0.9-20150924-095709\PostProcess\volumeFraction.txt');
% data = importdata('\\centos\secured\shared_data\diffWidth3VolFracFilling\0.1-20150925-125854\PostProcess\volumeFraction.txt');
% data = importdata('\\centos\secured\shared_data\diffWidth3VolFracFilling\0.3-20150925-124007\PostProcess\volumeFraction.txt');
% data = importdata('D:\SimResults\Chrono\SmarticleU\tests\fastShapeChange1013\PostProcess\volumeFraction.txt');
% data = importdata('D:\SimResults\Chrono\SmarticleU\tests\10-19 configChangeWithOTMovie\PostProcess\volumeFraction.txt');
% data = importdata('\\centos\secured\shared_data\ConfigChangeWithOTColoring\0.5-40-90-20151019-150116\PostProcess\volumeFraction.txt');
data = importdata('\\centos\secured\shared_data\ConfigChangeWithOTColoring\0.5-40-90-20151020-115542\PostProcess\volumeFraction.txt');
time        = data(:,1);
smartcount  = data(:,2);
volfrac     = data(:,3);
zmax        = data(:,4);
zcom        = data(:,5);
meanOT      = data(:,6);
guid        = data(:,7);
close all;
% volfrac= (smartcount.*v./(rad^2*pi*2*zcom));
% Time << CountInside << volumeFraction << zMax << zCom <<GUID
cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7]};
plotNames = {'Smarticles','Gait','U-Shape','Straight','Tetris','Vib at \circ','Vib Angle'};
plotTypes = [1];
hold on;

shapeLines=getShapeLines(time,guid);

lineVar= zcom;
plot(time,lineVar)
for i=1:size(shapeLines,1)
    plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) max(lineVar)],'color',shapeLines(i,2:4),'LineWidth',2)
    text(shapeLines(i,1),max(lineVar)*1.02,plotNames(shapeLines(i,5)))
end
% title('\phi vs. time')
title('Mean COM of pile vs. time')
xlabel('time [s]');
ylabel('<h> [m]');
% ylabel('<h> [m], \phi');
% ylabel('\phi');
vibTime = find(time>.75,1);
vibrate = false;
if vibrate
    plot(time(vibTime),volfrac(vibTime),'k.','MarkerSize',20);
    text(time(vibTime-2),1.04*lineVar(vibTime),'\phi')
end


set(gca,'Position',...
        get(gca,'OuterPosition') - ...
        get(gca,'TightInset') * ...
        [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1])
%% Plot meanOT
figure(2);
hold on;
plot(time,meanOT);

for(i=1:size(shapeLines,1))
    plot([shapeLines(i,1),shapeLines(i,1)],[min(meanOT) max(meanOT)],'color',shapeLines(i,2:4),'LineWidth',2)
    text(shapeLines(i,1),max(meanOT)*1.02,plotNames(shapeLines(i,5)))
end

title('Mean reaction torque on arms vs. time')
xlabel('time [s]');
ylabel('<|\tau|> [Nm]')

% data2 = importdata('D:\SimResults\Chrono\SmarticleU\tests\com changing shape\r2\PostProcess\volumeFraction.txt');
% data3 = importdata('D:\SimResults\Chrono\SmarticleU\tests\com changing shape\r3\PostProcess\volumeFraction.txt');
% % totData= [zcom,data2(:,5)];
% totData= [zcom,data2(:,5),data3(:,5)];
% 
% % plot(time,mean(totData,2));
% err = std(totData,0,2);
% shadedErrorBar(time,mean(totData,2),err,{'color',[12/255,100/255,2/255]});
% 
% figure(2);
% hold on;
% plot(time,zcom);
% plot(data2(:,1),data2(:,5));
% plot(data3(:,1),data3(:,5));
% 
% figure(2);
% hold on;
% plot(time,volfrac);
% plot(data2(:,1),data2(:,3));
% % plot(data3(:,1),data3(:,5));

%%

% figure(3);
% x=[];
% y=[];
% data = importdata('\\centos\secured\shared_data\diffWidthVolFracFilling\0.1\PostProcess\volumeFraction.txt');
% y= [y data(end,3)];
% x= [x .1];
% data = importdata('\\centos\secured\shared_data\diffWidthVolFracFilling\0.3\PostProcess\volumeFraction.txt');
% y= [y data(end,3)];
% x= [x .3];
% data = importdata('\\centos\secured\shared_data\diffWidthVolFracFilling\0.5\PostProcess\volumeFraction.txt');
% y= [y data(end,3)];
% x= [x .5];
% data = importdata('\\centos\secured\shared_data\diffWidthVolFracFilling\0.7\PostProcess\volumeFraction.txt');
% y= [y data(end,3)];
% x= [x .7];
% data = importdata('\\centos\secured\shared_data\diffWidthVolFracFilling\0.9\PostProcess\volumeFraction.txt');
% y= [y data(end,3)];
% x= [x .9];
% % 
% plot(x,y,'o-');