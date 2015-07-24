mainFolds = 'D:\SimResults\Chrono\SmarticleU\runs\newContact\';
folds = ['run 10', 'run 9', 'run 8', 'run 7', 'run 6', 'run 5' 'run 4'];
lw=[];
finalphi = [];
lwlegend={};
close all

figure(1);
hold on;
ylabel('\phi');
xlabel('Time(s)');


lw=[lw .2];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\halfwidth\run 23\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

% lw=[lw .3];
% [data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\halfwidth\run 21\PostProcess\volumeFraction.txt');
% % shakingIdx = find(data(:,1)>10,1);
% shakingIdx = find(data(:,1),1);
% plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
% finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

lw=[lw 0.3];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\Results\2015720 lw=0.3\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

lw=[lw 0.4];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\Results\2015720 lw=0.4\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

lw=[lw 0.5];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\Results\2015720 lw=0.5\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

lw=[lw 0.6];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\Results\2015720 lw=0.6\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

lw=[lw .7];  %didnt fill well enough fyi
[data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\halfwidth\run 17\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

lw=[lw .8];  %didnt fill well enough fyi
[data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\halfwidth\run 18\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

lw=[lw 0.9];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\halfwidth\run 19\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];


lw=[lw 1.0];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\Results\2015720 lw=1\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];


lw=[lw 1.1];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\halfwidth\run 20\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];


lw=[lw 1.2];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\Results\2015720 lw=1.2\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];


lw=[lw 1.3];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\Results\2015720 lw=1.3\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

lw=[lw 1.4];
[data]=importdata('D:\SimResults\Chrono\SmarticleU\Results\2015720 lw=1.4\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
shakingIdx = find(data(:,1),1);
plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
finalphi = [finalphi; lw(end),mean(data(end-5:end,3))];

% lw=[lw 1.4];
% [data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\halfwidth\run 15\PostProcess\volumeFraction.txt');
% % shakingIdx = find(data(:,1)>10,1);
% shakingIdx = find(data(:,1),1);
% plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3));
% finalphi = [finalphi; lw(end),data(end,3)];


% lw=[lw .4];
% [data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\newContact\run 4b\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
% plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3)*2);
% finalphi = [finalphi; lw(end),data(end,3)*4];
% 
% lw=[lw .5];
% [data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\newContact\run 5b\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>10,1);
% plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3)*2);
% finalphi = [finalphi; lw(end),data(end,3)*4];
% 
% 
% lw=[lw .6];
% [data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\newContact\run 6\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>5,1);
% plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3)*2);
% finalphi = [finalphi; lw(end),data(end,3)*4];
% 
% 
% lw=[lw .7];
% [data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\newContact\run 7\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>5,1);
% plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3)*2);
% finalphi = [finalphi; lw(end),data(end,3)*4];
% 
% 
% lw=[lw .8];
% [data]=importdata('D:\SimResults\Chrono\SmarticleU\runs\newContact\run 8\PostProcess\volumeFraction.txt');
% shakingIdx = find(data(:,1)>5,1);
% plot(data(shakingIdx:end,1)-data(shakingIdx,1),data(shakingIdx:end,3)*2);
% finalphi = [finalphi; lw(end),data(end,3)*4];
% 



for i=1:length(lw)
    lwlegend{i}=num2str(lw(i));
end
legend(lwlegend);

figure(2)
hold on;
xlabel('l/w');
ylabel('\phi_f');
plot(finalphi(:,1),finalphi(:,2),'o-');