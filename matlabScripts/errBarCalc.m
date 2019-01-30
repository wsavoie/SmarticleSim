function [OutM] = errBarCalc(xdatin,ydatin)

[uout, ~, outidx] = unique(xdatin,'sorted'); %retrieve unique values and indices of xdata
OutM = zeros(length(uout),3);
OutMatrix1  = [uout(:), accumarray( outidx,ydatin, [], @mean ) ];  %mean of ydata
OutMatrix2  = [uout(:), accumarray( outidx, ydatin, [], @std ) ]; %std of ydata
OutM(:,1:2) = OutMatrix1; % [depvar,ydataMean ]
OutM(:,3)   = OutMatrix2(:,2);%[std];
%OutM = [xdat,ydat,err]