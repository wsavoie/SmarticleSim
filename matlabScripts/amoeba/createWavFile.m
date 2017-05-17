load handel.mat
sampRate = 8192; %default sample rate of 8192 hertz
ff=[325];
dd=[66];
tt=(0:sampRate*dd)/sampRate;
y=sin(tt*ff*2*pi);





filename = '325hz 66s.wav';
audiowrite(filename,y,sampRate);