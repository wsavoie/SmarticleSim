tic
Nout=generatePackingFromSimDatPLANE(dat(1),0); 
a1=toc;

tic
Nout=generatePackingFromSimDatPLANE_mex(dat(1),0);
a2=toc;

pts('time elapsed=',a2,' speedup =',a1/a2);