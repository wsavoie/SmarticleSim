% function [out]=smarticleMSD(V,tracks)
%V = videoreader
%
load('D:\ChronoCode\chronoPkgs\Smarticles\matlabScripts\amoeba\smarticleExpVids\new room\movieInfo.mat');
figure(1);
% h=figure(2);
ma = msdanalyzer(2,'m','s');
ma = ma.addAll(movs(5).data);
% ma = ma.computeDrift('velocity');
ma.plotTracks;
figText(gcf,14);
plot(ma.tracks{1}(1,2),ma.tracks{1}(1,3),'ro','markersize',8,'MarkerFaceColor','r');
text(ma.tracks{1}(1,2),ma.tracks{1}(1,3)+.01,'start')
plot(ma.tracks{1}(end,2),ma.tracks{1}(end,3),'ko','markersize',8,'MarkerFaceColor','k');
text(ma.tracks{1}(end,2),ma.tracks{1}(end,3)+.01,'end')
ma.labelPlotTracks

figure(5);
plot(ma.tracks{1}(:,1),ma.tracks{1}(:,2))
ylabel('X (m)');
xlabel('time (s)');
figText(gcf,14);
figure(6);
plot(ma.tracks{1}(:,1),ma.tracks{1}(:,3))
ylabel('Y (m)');
xlabel('time (s)');
figText(gcf,14);
%%
figure(2);

ma = ma.computeMSD;
ma.plotMeanMSD(gca, true);
hold on;
t=ma.msd{1}(:,1);
% f = fit(t,ma.msd{1}(:,2),'b*x^m');
% plot(t,f.b.*t.^f.m,'r');
figText(gcf,14);
% text(.5,.8,['\langler^2\rangle\propto \tau^\alpha, \alpha\approx',num2str(f.m,'%2.1f')],'units','normalized','fontsize',16)
% ma = ma.fitMSD;
A = ma.getMeanMSD;
t = A(:, 1); % delay vector
msd = A(:,2); % msd
std_msd = A(:,3); % we will use inverse of the std as weights for the fit
std_msd(1) = std_msd(2); % avoid infinity weight

% ft = fittype('a*x + c*x^2');
% [fo, gof] = fit(t, msd, ft, 'Weights', 1./std_msd, 'StartPoint', [0 0]);
% 
% hold on
% plot(fo)
% legend off
% ma.labelPlotMSD
% 
% Dfit = fo.a / 4;
% Vfit = sqrt(fo.c);
% 
% ci = confint(fo);
% Dci = ci(:,1) / 4;
% Vci = sqrt(ci(:,2));
% 
% fprintf('Parabolic fit of the average MSD curve with 95%% confidence interval:\n')
% 
% fprintf('D = %.3g [ %.3g - %.3g ] %s, real value was %.3g %s\n', ...
%     Dfit, Dci(1), Dci(2), [SPACE_UNITS '²/' TIME_UNITS], D, [SPACE_UNITS '²/' TIME_UNITS]);
% 
% fprintf('V = %.3g [ %.3g - %.3g ] %s, real value was %.3g %s\n', ...
%     Vfit, Vci(1), Vci(2), [SPACE_UNITS '/' TIME_UNITS], vm, [SPACE_UNITS '/' TIME_UNITS]);
%%
figure(3)
ma = ma.computeVCorr;
ma.plotMeanVCorr;

% meanIdxs = find(tracks{1}(:,1)>10 & tracks{1}(:,1)<10.5);
M = mean(ma.vcorr{1}(10:end,2));
line([ma.vcorr{1}(10,1) ma.vcorr{1}(end,1)], [M M],'color','r','linewidth',3);
figText(gcf,14);
text(.5,.7,['mean = ',num2str(M,'%2.3f')],'units','normalized','fontsize',16)