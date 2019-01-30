function [t,ringPos,smartPos,cogPos]=readDeadRing(fname)
simD=importdata(fname);
% tstep, ringx, ringy, ringz, 
% ringcomX, ringcomY, ringcomZ, 
% ringdeadX, ringdeadY, ringdeadZ

% simD.data=simD.data(simD.data(:,1)>.2,:);
simD.data(:,2)=-simD.data(:,2);
simD.data(:,5)=-simD.data(:,5);
simD.data(:,8)=-simD.data(:,8);
%use negative in x dir since it is a left handed coordinate system
ringPos=[simD.data(:,2) simD.data(:,3) ]; %only need time, x, y
cogPos=[simD.data(:,5) simD.data(:,6)];
smartPos=[simD.data(:,8) simD.data(:,9)];

%subtract off initial position
smartPos=bsxfun(@minus, smartPos,ringPos(1,:));
ringPos=bsxfun(@minus, ringPos,ringPos(1,:));


t=simD.data(:,1);