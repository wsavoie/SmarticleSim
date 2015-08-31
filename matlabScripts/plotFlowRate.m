% directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests\PostProcess');
directory_name='D:\SimResults\Chrono\SmarticleU\tests\PostProcess';
data=importdata(horzcat(directory_name,'\flowrate.txt'));

dt = .0005;
FS = 1/dt;

cols = ['r','g','b','k','c','m'];
plotNames = {'Smarticles','Gait','Staple','Flat','Tetris','Vib at \circ','Vib Angle'};
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
gg= sortrows(gg,2); 
%loop make sure each different line type has proper name in legend, and
%that each configuration only appears once in legend
for i=1:size(gg,1)
    if i~=1
        if gg(i,2)==gg(i-1,2) %if true will mean a repeat in color on legend
            plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],cols(gg(i,2)),'LineWidth',2)
            continue;
        else
            sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],cols(gg(i,2)),'LineWidth',2)];
        end
    else
        sp1 = [sp1 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(data(:,2)),max(data(:,2))],cols(gg(i,2)),'LineWidth',2)];
    end
plotTypes=[plotTypes, (gg(i,2)+1)];
end
legend([sp1(1:end)],plotNames(plotTypes));
axis([data(1,1) data(end,1) data(1,2),data(end,2)]);
plotTypes
%% subplot 2
subplot(2,1,2)
hold on;
title('Flow Rate of Smarticles through hopper');
% flowRate = diff(data(:,2))./diff(data(:,1));
flowRate = diff(data(:,2))./diff(data(1:2,1)); 
% flowRate = diff(data(:,2))./diff(data(:,1));
plotNames = {'Flow Rate','Flow Rate Avged','Gait','Staple','Flat','Tetris','Vib at \circ','Vib Angle'};
plotTypes2 = [1,2];
sp2= [plot(data(1:end-1,1),flowRate)]; %plotting noisy regular differential
window_size = 100;%average over 50 dts (
simple = tsmovavg(flowRate,'s',window_size,1); %moving averaged differential
sp2 =[sp2 plot(data(1:end-1,1),simple)]; %
xlabel('Time(s)');
ylabel('Flow Rate D(Smarticles)/{d(s)}');
for i=1:size(gg,1)
if i~=1
    if gg(i,2)==gg(i-1,2) %if true will mean a repeat in color on legend
        plot([data(gg(i,1),1),data(gg(i,1),1)],[min(flowRate),max(flowRate)],cols(gg(i,2)),'LineWidth',2)
         continue;
    else
        sp2 = [sp2 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(flowRate),max(flowRate)],cols(gg(i,2)),'LineWidth',2)];
    end
else
    sp2 = [sp2 plot([data(gg(i,1),1),data(gg(i,1),1)],[min(flowRate),max(flowRate)],cols(gg(i,2)),'LineWidth',2)];
end
plotTypes2=[plotTypes2, (gg(i,2)+2)];
end
legend([sp2(1:end)],plotNames(plotTypes2));
axis([data(1,1) data(end,1) min(flowRate) max(flowRate)]);