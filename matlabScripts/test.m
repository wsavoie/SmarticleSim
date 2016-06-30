% path = 'A:\SmarticleRun\Active98Pos%\';
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

SPACE_UNITS = 'µm';
TIME_UNITS = 's';

N_PARTICLES = 1;
N_TIME_STEPS = 1000;

% Diffusion coefficient. Will set the amplitude of the random displacement
D  = 1e-3; % µm^2/s
% Time step between acquisition; fast acquisition!
dT = 0.05; % s,

% Mean velocity
vm = 0.05; % µm/s

% Area size, just used to disperse particles in 2D. Has no impact on
% analysis.
SIZE = 2; % µm

tracks = cell(N_PARTICLES, 1);
clf;
k = sqrt(2 * D * dT);
for i = 1 : N_PARTICLES

    % Time
    time = (0 : N_TIME_STEPS-1)' * dT;

    % Velocity orientation
    theta = 2 * pi * rand;
    rx=2 * pi * rand(N_TIME_STEPS,1);
    ry=2 * pi * rand(N_TIME_STEPS,1);
    % Mean velocity
    v = vm * (1 + 1/4*randn);

    % Initial position
    X0 = SIZE .* rand(1, 2);
    pd = makedist('tLocationScale','mu',0,'sigma',1,'nu',2);
%     pd = makedist('normal');
    vv=random(pd,N_TIME_STEPS,2);    
%     dX_levy = [xvals.*cos(rx) yvals.*sin(ry)];
    mm=max(abs(vv(:)))
    dX_levy = [vv]/mm;
    % Instantaneous displacement:
    dX_brownian = k * randn(N_TIME_STEPS, 2);
    
    dX_directed = v * dT * ...
        [cos(theta)*ones(N_TIME_STEPS,1) sin(theta)*ones(N_TIME_STEPS,1) ];

    % Integrate uncorrelated displacement
%     dX = dX_brownian + dX_directed;
    dX = dX_levy;
    dX(1, :) = X0;
    X = cumsum(dX, 1);

    % Store
    tracks{i} = [time X];

end
clear i X dX time X0
ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
ma = ma.addAll(tracks);
ma.plotTracks
ma.labelPlotTracks