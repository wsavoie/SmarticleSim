function [out1,out2] = interpolateGait(pts,ss)
% pts = gaitPts;

d=0;%vector length
dt1 = diff(pts);
%square each element, sum along column, sqrt that value, sum all values
%norm of each column of dt1
curveDist=sum(sqrt(sum(dt1.^2,2)));
curvePts=curveDist/(ss);
y=curvspace(pts,curvePts);

out1 = y(:,1);
out2 = y(:,2);
end
