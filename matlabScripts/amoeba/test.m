SPACE_UNITS = 'µm';
TIME_UNITS = 's';

N_PARTICLES = 100;
N_TIME_STEPS = 200;
N_DIM = 2; % 2D
SIZE = 10; % µm

D  = 1e-3; % µm^2/s
dT = 0.03; % s
k = sqrt(2 * D * dT);

% Drift parameters: a circle moving forward
drift_time = (0 : N_TIME_STEPS-1)' * dT;
drift_time = msdanalyzer.roundn(drift_time, 3);
drift_pos = [ 1e-1*drift_time + 1*cos(drift_time), 1e-1*drift_time + 1*sin(drift_time) ];
drift = [ drift_time drift_pos];

% Generate tracks
tracks = cell(N_PARTICLES, 1);
for i_spot = 1 : N_PARTICLES

    % Time
    time = (0 : N_TIME_STEPS-1)' * dT;

    % Initial position
    X0 = SIZE .* rand(1, N_DIM);

    % Integrate uncorrelated displacement
    dX = k * randn(N_TIME_STEPS, N_DIM);  %#ok<*UNRCH>
    dX(1, :) = X0;
    X = cumsum(dX, 1);

    % Add drift
    X = X + drift_pos;

    % Store
    tracks{i_spot} = [ time X];

end

%%
ma = msdanalyzer(2, 'um', 's');
ma = ma.addAll(tracks);
[hps, ha] = ma.plotTracks;
ma.labelPlotTracks
%%
ma = ma.computeMSD;
figure
ma.plotMeanMSD(gca, true);
ma.fitMeanMSD

%%
ma = ma.fitLogLogMSD;
r2fits = ma.loglogfit.r2fit;
alphas = ma.loglogfit.alpha;

R2LIMIT = 0.8;

% Remove bad fits
bad_fits = r2fits < R2LIMIT;
fprintf('Keeping %d fits (R2 > %.2f).\n', sum(~bad_fits), R2LIMIT);
alphas(bad_fits) = [];

% Echo results
fprintf('alpha = %.2f ± %.2f (mean ± std, N = %d).\n', ...
    mean(alphas), std(alphas), numel(alphas));

%%
to_erase = rand( size(drift,1), 1) < 0.1;
to_erase(1, end) = false; % ensure we do not erase the last and first ones
drift(to_erase, :) = [];

%%
ma = ma.computeDrift('manual', drift);

%%
ma.plotDrift(ha) % to plot it in the track figure
size(drift) % measured drift
size(ma.drift) % interpolated drift
ma = ma.computeMSD;
ma.fitMeanMSD
ma = ma.fitLogLogMSD;
r2fits = ma.loglogfit.r2fit;
alphas = ma.loglogfit.alpha;

R2LIMIT = 0.8;

% Remove bad fits
bad_fits = r2fits < R2LIMIT;
fprintf('Keeping %d fits (R2 > %.2f).\n', sum(~bad_fits), R2LIMIT);
alphas(bad_fits) = [];

% Echo results
fprintf('alpha = %.2f ± %.2f (mean ± std, N = %d).\n', ...
    mean(alphas), std(alphas), numel(alphas));

ma = ma.computeDrift('centroid'); % note: no extra parameters.
ma = ma.computeMSD;
ma.fitMeanMSD

%%
b=2; a=1;
p=0:.01:12*pi;
x=a*p-b*sin(p);
z=a-b*cos(p);
y=p;
figure(1);

hold on;


x1 = 0:40 ;
y1 =x1;
z1=0;
hold on;
[xx,yy,zz]=meshgrid(x1,y1,z1);
plot3(x(z>0),y(z>0),z(z>0),'linewidth',2)
surf(xx,yy,zz,'edgecolor','none');
plot3(x,y,z,'linewidth',2);
grid on
alpha(0.5)


figure(2);
plot3(x(z>0),y(z>0),z(z>0),'linewidth',2)
