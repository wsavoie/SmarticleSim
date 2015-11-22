%plot stress
filename = 'D:\SimResults\Chrono\SmarticleU\tests\plotStress\knob\PostProcess\Stress.txt';
% filename = 'D:\SimResults\Chrono\SmarticleU\tests\PostProcess\Stress.txt';
%stressdata
sd=importdata(filename);
time=sd(:,1);
stress=sd(:,2);
guid=sd(:,3);

plotNames = {'Stress','Gait','U-Shape','Straight','Tetris','Vib at \circ','Vib Angle'};

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
plot(time,lineVar);
plot(time,y,'LineWidth',4);


hold on;
shapeLines=getShapeLines(time,guid);

plot(time,lineVar)
for i=1:size(shapeLines,1)
    plot([shapeLines(i,1),shapeLines(i,1)],[min(lineVar) max(lineVar)],'color',shapeLines(i,2:4),'LineWidth',2)
    text(shapeLines(i,1),max(lineVar)*1.02,plotNames(shapeLines(i,5)))
end

xlabel('Time (s)');
ylabel('Stress');
figText(gcf,14);

