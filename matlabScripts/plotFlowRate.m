directory_name = uigetdir('A:\SmarticleRun');
% directory_name='D:\SimResults\Chrono\SmarticleU\tests\PostProcess';
data=importdata(horzcat(directory_name,'\flowrate.txt'));
figure(1);
[t,ind] = unique(data(:,1));
data=data(ind,:);
dt = .00025;
FS = 1/dt;  
cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7]};

hold on;

%data = [time exited hopper;   total count;    guid]
% suptitle('Flow of Smarticles in a Hopper')
subplot(2,1,1);
hold on;
a= gca;
hold on;
title('Smarticle passed through hopper');
xlabel('Time(s)');
ylabel('Smarticles');

plot(data(:,1),data(:,2),'linewidth',2); %starts variable for the legend
%loop make sure each different line type has proper name in legend, and
%that each configuration only appears once in legend
leg = legend('Smarticles');
axis([data(1,1) data(end,1) data(1,2),data(end,2)]);
ax = gca;
ax.XTick= 0:1:max(ax.XTick);

lineVar=data(:,2);
yAx=max(lineVar);
plotNames = {'Stress','Gait','U-Shape','Straight','n-Shape','Vib at \circ','Vib Angle'};
shapeLines=getShapeLines(data(:,1),data(:,3));
for i=1:size(shapeLines,1)
    plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) yAx],'color',shapeLines(i,2:4),'LineWidth',5)
    text(shapeLines(i,1)+.1,yAx*0.98,plotNames(shapeLines(i,5)))
end

%% subplot 2
figure(1);
subplot(2,1,2)
hold on;
title('Flow Rate of Smarticles through hopper');
% flowRate = diff(data(:,2))./diff(data(:,1));
flowRate = diff(data(:,2))./diff(data(:,1)); 
sp2=[];
xlabel('Time(s)');
ylabel('Flow Rate d(Smarticles)/{dt}');
% //ylabel('Smarticle Flow Rate');


hold on;
n = 1; % designs an Nth order low pass
f = 1/dt;
% f = 10;
we = 1*2*pi*dt;
if we >1
    we = .005; % .005*1/dt/2 5
end
[B,A] = butter(n,we);
% filters the data in vector X with the filter
%   described by vectors A and B to create the filtered data Y.  The
%   filter is described by the difference equation:
%
%     a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                           - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
ydat2 = filtfilt(B,A,flowRate);
% plot(data(1:end-2,1),flowRate(2:end),'LineWidth',2);


x1=data(6:end,1);
y1=flowRate(5:end);
factor(length(y1))
fac=max(factor(length(y1)));

x2 = mean(reshape(x1,fac,[]));
y2= mean(reshape(y1,fac,[]));
% plot(x1,y1,'LineWidth',2);
plot(x2,y2,'.-','markersize',22,'LineWidth',2);
legend('Smarticle flow rate')
lineVar=flowRate;
yAx=max(lineVar);
shapeLines=getShapeLines(data(:,1),data(:,3));
plotNames = {'Stress','Gait','U-Shape','Straight','n-Shape','Vib at \circ','Vib Angle'};
for i=1:size(shapeLines,1)
    plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) yAx],'color',shapeLines(i,2:4),'LineWidth',5)
    text(shapeLines(i,1)+.1,yAx*0.98,plotNames(shapeLines(i,5)))
end

% axis([data(1,1) data(end,1) min(flowRate) max(flowRate)]);
 axis([data(1,1) data(end,1) min(y2) max(y2)]);
ax = gca;
ax.XTick= 0:1:max(ax.XTick);

% 
% %% comparing %
% 
% 
% dir_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests\PostProcess');
% x = dir('D:\SimResults\Chrono\SmarticleU\tests\MixedSmartWithOT\r1');
% c = {x(3:end).name};
% sp1 = [];
% for j=1:length(c)
%     
%     directory_name=horzcat(dir_name,'\',c{j},'\PostProcess')
%     data=importdata(horzcat(directory_name,'\flowrate.txt'));
%     
%     dt = .0005;
%     FS = 1/dt;
%     cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7]};
%     plotNames = {'Smarticles','Gait','U-Shape','Straight','Tetris','Vib at \circ','Vib Angle'};
%     plotTypes = [1];
%     figure(1);
%     gc =[0; find(diff(data(:,3)))];
%     gc=gc+1; %time index will be 1 off from diff
%     val=[data(gc,3)]+1;%added value can be zero (global) and matrices are 1 started
%     hold on;
% 
%     % suptitle('Flow of Smarticles in a Hopper')
%     a= gca;
%     % hold on;
%     title('Smarticle passed through hopper');
%     xlabel('Time(s)');
%     ylabel('Smarticles');
% 
%    
%     gg = [gc, val];
%     ggOld = gg;
%     gg= sortrows(gg,2); 
%     colormap hot
%     sp1=[sp1 plot(data(:,1),data(:,2))]; %starts variable for the legend
%     gg = [gc, val];
%     ggOld = gg;
%     gg= sortrows(gg,2); 
%     %loop make sure each different line type has proper name in legend, and
%     %that each configuration only appears once in legend
%     for i=1:size(gg,1)
%         if i~=1
%             if gg(i,2)==gg(i-1,2) %if true will mean a repeat in color on legend
%                 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)
%                 continue;
%             else
%                 if j==length(c)
%                     sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
%                 end
%             end
%         else
%             if j==length(c)
%                sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
%             end
%         end
%     plotTypes=[plotTypes, (gg(i,2)+1)];
%     end
% end
% legend([sp1(1:end)],horzcat({'0','.25','.5','.75','1'},plotNames(plotTypes)));
% axis([data(1,1) data(end,1) data(1,2),data(end,2)]);
% ax = gca;
% ax.XTick= 0:1:max(ax.XTick);