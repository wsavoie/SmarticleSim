function [force,angs,cogpos,rot,fcs,ringpos]=readInactive(fname,binWid)
% fold='A:\SmarticleRun\Amoeba_COG_dead_1_pos_+x\f_0.2_rob_5_v_9\';
% fname=fullfile(fold,'PostProcess','RingContact.txt');
% binWid=5;

simD=importdata(fname);
%time, xContact, yContact, zContact, xForce, yForce, zForce
%,armposx,armposy,armposz, armrote0,armrote1,armrote2,armrote3,
%,ringx, ringy, ringz
if(isfield(simD,'data'))
simD.data=simD.data(simD.data(:,1)>.2,:);
simD.data(:,2)=-simD.data(:,2);
simD.data(:,5)=-simD.data(:,5);
simD.data(:,8)=-simD.data(:,8);
simD.data(:,15)=-simD.data(:,15);
%use negative in x dir since it is a left handed coordinate system
pos=[simD.data(:,1) simD.data(:,2) simD.data(:,3) ]; %only need time, x, y
force=[simD.data(:,1) simD.data(:,5) simD.data(:,6)];
cogpos=[simD.data(:,1) simD.data(:,8) simD.data(:,9)];
ringpos=[simD.data(:,15) simD.data(:,16)];

q=quaternion([simD.data(:,11),simD.data(:,12),simD.data(:,13),simD.data(:,14)]);

rz=EulerAngles(q(:),'123');
rot=[pi-rz(3,:)]';
% %get radius of circle
% % [match, nomatch]=regexp(simD.textdata,'[0-9.]+','match');
% % r=str2double(match{1});

% angs=atan2(pos(:,3),pos(:,2))+deg2rad(180);

angs=atan2(pos(:,3),pos(:,2));
angs=mod(angs,2*pi);


binW=deg2rad(binWid); %degrees
% disc=discretize(angs,0:binW:2*pi);
disc=discretize(angs,linspace(0,2*pi,round(2*pi/binW)+1));
ud=unique(disc);
lenuni=length(ud);
fcs=zeros(lenuni,2);

% sum(norm(force(disc==ud(:),2:3)));
for(i=1:lenuni)
fcs(i,1)=ud(i)*binW;
fcs(i,2)=sum(norm(force(disc==ud(i),2:3)));
end
else
    [force,angs,cogpos,rot,fcs,ringpos]=deal(0);
end