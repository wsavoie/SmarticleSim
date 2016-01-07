while 1
clf;
%plot stress
% filename = 'D:\SimResults\Chrono\SmarticleU\tests\s5\PostProcess\Stress.txt';
filename = 'D:\SimResults\Chrono\SmarticleU\tests\PostProcess\Stress.txt';
%stressdata
sd=importdata(filename);
time=sd(:,1);
stress=sd(:,2);
guid=sd(:,3);
cylRad =sd(:,4);

plotNames = {'Stress','Gait','U-Shape','Straight','n-Shape','Vib at \circ','Vib Angle'};

lineVar= stress;
dt=time(1);
fs = 1/dt; % Sampling rate (Hz)
%d1 = designfilt('lowpassiir','FilterOrder',12, ...
%    'PassbandFrequency',100,'SampleRate',1/dt,'DesignMethod','butter');
d1 = designfilt('lowpassiir', 'PassbandFrequency', 15, ...
                'StopbandFrequency', 40, 'PassbandRipple', 15, ...
                'StopbandAttenuation', 50, 'SampleRate', 1/dt, ...
                'DesignMethod', 'butter');
y = filtfilt(d1,lineVar);
figure(1);
hold on;




shapeLines=getShapeLines(time,guid);
% figure(12)
subplot(2,1,1);
hold on;
plot(time,lineVar);
line = plot(time,y,'YDataSource','y','LineWidth',4);

yAx=mean(lineVar)*4;
for i=1:size(shapeLines,1)
    plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) yAx],'color',shapeLines(i,2:4),'LineWidth',5)
    text(shapeLines(i,1),yAx*1.02,plotNames(shapeLines(i,5)))
end

axis([0 time(end) 0 yAx]);




ylabel('Stress');
subplot(2,1,2);
plot(time,cylRad,'Color',line.Color);
xlabel('Time (s)');
ylabel('bucket radius (m)');
axis auto
figText(gcf,14);
refreshdata(line,'caller')
drawnow
end