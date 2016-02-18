ff=[125 175 225 125];
dd=[8 8 8 8];

ff=[275 325];
dd=[16 16];
sampRate = 8192; %default sample rate of 8192 hertz
for nn=1:length(dd);
tt=(0:sampRate*dd(nn))/sampRate;
soundsc(sin(tt*ff(nn)*2*pi));
pause(dd(nn));
end