% ff=[175 175 225 125];
% dd=[8 8 8 8];
tic
% ff=[125 125 125 225 450 225 450 225 450 225 450];
ff=[175 225 175 225 175 225];
dd=[5 5 5 5 5];
% ff=[225];
% dd=[15];
sampRate = 8192; %default sample rate of 8192 hertz
for nn=1:length(ff);
tt=(0:sampRate*dd(nn))/sampRate;
soundsc(sin(tt*ff(nn)*2*pi));
pause(dd(nn));
end
toc
