ff=[125 175 225 125];
dd=[8 8 8 8];

ff=[275 325];
dd=[16 16];

for nn=1:length(dd);
tt=(0:8192*dd(nn))/8192;
soundsc(sin(tt*ff(nn)*2*pi));
pause(dd(nn));
end