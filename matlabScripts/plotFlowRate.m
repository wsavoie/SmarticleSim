directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests\PostProcess');
% directory_name='D:\SimResults\Chrono\SmarticleU\tests\PostProcess';
data=importdata(horzcat(directory_name,'\flowrate.txt'));

dt = .0005;
FS = 1/dt;
cols = {[1,0,0],[1,.5,0],[113/255,188/255, 255/255],[0,0,0],[100,100,30],[30,27,95]};
plotNames = {'Smarticles','Gait','U-Shape','Straight','Tetris','Vib at \circ','Vib Angle'};
plotTypes = [1];
figure(1);
gc =[0; find(diff(data(:,3)))];
gc=gc+1; %time index will be 1 off from diff
val=[data(gc,3)]+1;%added value can be zero (global) and matrices are 1 started
hold on;

% suptitle('Flow of Smarticles in a Hopper')
subplot(2,1,1);
a= gca;
hold on;
title('Smarticle passed through hopper');
xlabel('Time(s)');
ylabel('Smarticles');

sp1=[plot(data(:,1),data(:,2))]; %starts variable for the legend
gg = [gc, val];
ggOld = gg;
gg= sortrows(gg,2); 
%loop make sure each different line type has proper name in legend, and
%that each configuration only appears once in legend
for i=1:size(gg,1)
    if i~=1
        if gg(i,2)==gg(i-1,2) %if true will mean a repeat in color on legend
            plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)
            continue;
        else
            sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
        end
    else
        sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
    end
plotTypes=[plotTypes, (gg(i,2)+1)];
end
legend([sp1(1:end)],plotNames(plotTypes));
axis([data(1,1) data(end,1) data(1,2),data(end,2)]);
ax = gca;
ax.XTick= 0:1:max(ax.XTick)
%% subplot 2
subplot(2,1,2)
hold on;
title('Flow Rate of Smarticles through hopper');
% flowRate = diff(data(:,2))./diff(data(:,1));
flowRate = diff(data(:,2))./diff(data(1:2,1)); 
% flowRate = diff(data(:,2))./diff(data(:,1));
plotNames = {'Flow Rate','Flow Rate Avged','Gait','U-Shape','Straight','Tetris','Vib at \circ','Vib Angle'};
plotTypes2 = [2];
sp2=[];
% sp2= [plot(data(1:end-1,1),flowRate)]; %plotting noisy regular differential
xlabel('Time(s)');
% ylabel('Flow Rate d(Smarticles)/{dt}');
ylabel('Smarticle Flow Rate');


hold on;
n = 1; % designs an Nth order low pass
timestep = .0005;
f = 1/dt;
% f = 10;
we = f*2*pi*timestep;
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
sp2=[sp2 plot(data(1:end-1,1),ydat2,'LineWidth',2)];


for i=1:size(gg,1)
if i~=1
    if gg(i,2)==gg(i-1,2) %if true will mean a repeat in color on legend
        plot([data(gg(i,1),1),data(gg(i,1),1)],[min(flowRate),max(flowRate)],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)
         continue;
    else
        sp2 = [sp2 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(flowRate),max(flowRate)],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
    end
else
    sp2 = [sp2 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(flowRate),max(flowRate)],'color',cell2mat(cols(gg(i,2))),'LineWidth',2)];
end
plotTypes2=[plotTypes2, (gg(i,2)+2)];
end
legend([sp2(1:end)],plotNames(plotTypes2));
% axis([data(1,1) data(end,1) min(flowRate) max(flowRate)]);
 axis([data(1,1) data(end,1) min(ydat2) max(ydat2)]);
ax = gca;
ax.XTick= 0:1:max(ax.XTick)