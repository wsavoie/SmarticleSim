function [dphi,phi_0,phi_f] = changeInPhiAmp(dat,vibInd)
%CHANGEINPHIAMP find change in phi amplitude
%   data is usedD(index), vib index is gui change index, "2" would mean
%   dropping staples into system at one config then vibrate

mI=20;%indices to take mean over

%get height before vibration

gc=findChangesInGui(dat.gui);
gci=gc(vibInd,1);
phi_0=mean(dat.phi((gci-mI):gci));

if(size(gc,1)>gci) %measure volume fraction up to next change
    gce=gc(vibInd+1,1);
    phi_f=mean(dat.phi((gce-mI):gce,1));
else
    phi_f=mean(dat.phi((end-20):end));
end

dphi=phi_f-phi_0;
end

