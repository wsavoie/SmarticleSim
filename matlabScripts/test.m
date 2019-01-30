%%
t=0:.001:30;

SPACE_UNITS='m';
TIME_UNITS='s';

ma = msdanalyzer(1, SPACE_UNITS, TIME_UNITS);
% t=0:.01:20;
for i=1:20
p1=sin(t+2*pi*rand());
c1={[t',p1']};
ma=ma.addAll(c1);
end



ma=ma.computeVCorr;
ma.plotMeanVCorr
% p4=[t,mean(p1+p2+p3)];
% autocorr(p4);
%%
% hold on;
% % plot(t,forcex,t,forcey,t,forcez);
% R=0.042528;
% RR=R/.8+.5*R/.8;
% marm=0.00355;
% mtot=0.0321;
% T=2*.4290;
% omeg=2*pi/T;
% farm=(R/2)*omeg^2*marm
% ftot=(R/2)*omeg^2*mtot
% plot(t,forcex);
% plot(t,forcey);
% plot(t,forcez);
% 
% title('Penalty');
% xlabel('time (s)');
% ylabel('Force (N)');
% figText(gcf,16);
% legend({'x force','y force', 'z force'})
% % 
% plot(t1,forcex1);
% plot(t1,forcey1);
% plot(t1,forcez1);
%%
% path = 'A:\SmarticleRun\TestActive\';
% % path = 'D:\SimResults\Chrono\SmarticleU\tests\lazy .05\';
% a=dir(horzcat(path,'-*'));
% for i=1:length(a)
%     ff= horzcat(path,a(i).name,'\PostProcess\Stress.txt');
%     pts(ff);
%     matFile =  horzcat(path,a(i).name,'\PostProcess\stressData.mat');
%     if exist(matFile, 'file') == 2 
%         pts('matrix file already created');
%         continue
%     end
%    
%     readAllSmarticlesAngles(ff,0);
% end
% beep;
% beep;
%%%%%%%%
% 
% SPACE_UNITS = 'µm';
% TIME_UNITS = 's';
% 
% N_PARTICLES = 1;
% N_TIME_STEPS = 1000;
% 
% % Diffusion coefficient. Will set the amplitude of the random displacement
% D  = 1e-3; % µm^2/s
% % Time step between acquisition; fast acquisition!
% dT = 0.05; % s,
% 
% % Mean velocity
% vm = 0.05; % µm/s
% 
% % Area size, just used to disperse particles in 2D. Has no impact on
% % analysis.
% SIZE = 2; % µm
% 
% tracks = cell(N_PARTICLES, 1);
% clf;
% k = sqrt(2 * D * dT);
% for i = 1 : N_PARTICLES
% 
%     % Time
%     time = (0 : N_TIME_STEPS-1)' * dT;
% 
%     % Velocity orientation
%     theta = 2 * pi * rand;
%     rx=2 * pi * rand(N_TIME_STEPS,1);
%     ry=2 * pi * rand(N_TIME_STEPS,1);
%     % Mean velocity
%     v = vm * (1 + 1/4*randn);
% 
%     % Initial position
%     X0 = SIZE .* rand(1, 2);
%     pd = makedist('tLocationScale','mu',0,'sigma',1,'nu',2);
% %     pd = makedist('normal');
%     vv=random(pd,N_TIME_STEPS,2);    
% %     dX_levy = [xvals.*cos(rx) yvals.*sin(ry)];
%     mm=max(abs(vv(:)))
%     dX_levy = [vv]/mm;
%     % Instantaneous displacement:
%     dX_brownian = k * randn(N_TIME_STEPS, 2);
%     
%     dX_directed = v * dT * ...
%         [cos(theta)*ones(N_TIME_STEPS,1) sin(theta)*ones(N_TIME_STEPS,1) ];
% 
%     % Integrate uncorrelated displacement
% %     dX = dX_brownian + dX_directed;
%     dX = dX_levy;
%     dX(1, :) = X0;
%     X = cumsum(dX, 1);
% 
%     % Store
%     tracks{i} = [time X];
% 
% end
% clear i X dX time X0
% ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
% ma = ma.addAll(tracks);
% ma.plotTracks
% ma.labelPlotTracks
x=zeros(50,2);
y=zeros(50,2);
for i = 1:50;
    x(i,1)=.8;
    x(i,2)=.2;
    y(i,1)=i*0.01+0.2;
    y(i,2)=i*0.01+0.2;
end

plot(x(:,1),y(:,1),'o',x(:,2),y(:,2),'o');

%% rename vids
fold=uigetdir('A:\SmarticleRun\');
f = dir2(fold,'folders');
for i=1:length(f)
    ff= fullfile(fold,f(i).name,'video_capture');
    if exist(fullfile(ff,'outVid.avi'),'file')==2
     movefile(fullfile(ff,'outVid.avi'),fullfile(ff,[f(i).name,'.avi']));
    else
        pts(fullfile(ff,'outVid.avi'),' exists');
    end
   
end
%% crosscor
% % load relatedsig.mat
% % [Rmm,lags]=xcorr(s1,s3);
% 
% 
% % Rmm = Rmm(lags>0);
% % lags = lags(lags>0);
% figure(55);
% hold on;
% t=0:0.01:2*pi;
% x=cos(t);
% y=sin(t);
% % both=x+y;
% % x=rand(10,1);
% % y=rand(10,1);
% [Rmm,lags]=xcorr(x,y,'unbiased');
% % Rmm = Rmm(lags>0);
% % lags = lags(lags>0);
% % 
% 
% % figure
% plot(lags*.01/pi,Rmm/max(Rmm),'linewidth',2)
% plot(t/pi,x,'r',t/pi,y,'g');
% 
% xlabel('Lag (s)');
% axis tight

%%
% aa=ma.getMeanMSD;
% [x]=aa(:,1);
% [y]=aa(:,2);
% 
% 
% y=y(x<15);
% x=x(x<15);
% y=y(x>1.2);
% x=x(x>1.2);
% 
% [lx]=log(x);
% [ly]=log(y);
% 
% subplot(1,2,1);
% plot(x,y);
% subplot(1,2,2);
% plot(lx,ly);
% 
% [f,gof]=fit(lx,ly,'poly1')

%%
% minn=0; midd=90; maxx=180; angSp=60;
% p1=1500; oldP1=1500;
% p2=1500; oldP2=1500;
% 
% oldP1=p1;
% oldP2=p2;
% p1=maxx * 10 + 600
% p2=minn * 10 + 600
% 
% oldP1=p1;
% oldP2=p2;
% p1=maxx * 10 + 600
% p2=maxx * 10 + 600
% max(abs(p1-oldP1),abs(p2-oldP2))/10.0