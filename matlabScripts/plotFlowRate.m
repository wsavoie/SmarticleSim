% directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests\PostProcess');
directory_name='D:\SimResults\Chrono\SmarticleU\tests\PostProcess';
data=importdata(horzcat(directory_name,'\flowrate.txt'));

dt = .0005;
FS = 1/dt;

cols = ['r','g','b','k','c'];
figure(1);
gc = find(diff(data(:,3)));
gc=gc+1; %time index will be 1 off from diff
val=data(gc,3)+1;%added value can be zero (global) and matrices are 1 started
hold on;
% suptitle('Flow of Smarticles in a Hopper')
subplot(2,1,1);
hold on;
title('Smarticle passed through hopper');
xlabel('Time(s)');
ylabel('Smarticles');
plot(data(:,1),data(:,2));
for i=1:size(gc,1)
plot([data(gc(i),1),data(gc(i),1)],[min(data(:,2)),max(data(:,2))],cols(val(i)),'LineWidth',2)
end

subplot(2,1,2)
hold on;
title('Flow Rate of Smarticles through hopper');
% flowRate = diff(data(:,2))./diff(data(:,1));
flowRate = diff(data(:,2))/.0005;
% flowRate = diff(data(:,2))./diff(data(:,1));
plot(data(1:end-1,1),flowRate);
for i=1:size(gc,1)
    plot([data(gc(i),1),data(gc(i),1)],[min(flowRate),max(flowRate)],cols(val(i)),'LineWidth',2)
end
xlabel('Time(s)');
ylabel('Flow Rate D(Smarticles)/{d(s)}');

