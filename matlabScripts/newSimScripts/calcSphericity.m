function [psi] = calcSphericity(smartLocs)
%CALCSPHERICITY calculate sphericity
%smartloc is locations of smarticles
[x,y,z]=separateVec(smartLocs,1);
shp=alphaShape(x,y,z,inf);
vol=volume(shp);
ar=surfaceArea(shp);
psi=pi^(1/3)*(6*vol)^(2/3)/ar;
end

