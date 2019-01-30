function [smartOut] = RemoveSmartsFromBucket(smartLocs,topN,hookHeight,brad)
%REMOVESMARTSFROMBUCKET remove indices from smart pos where locations are
%unreasonable for sphericity calculation, this includes outside bucket, far
%away from other staples,below hook if applicable
%smartLocs=[x,y,z]
%topN smarts to remove from system
%hookheight=if nonzero cull smarts underneath hook (z<hookHeight)
%brad= if nonzero, cull smarts outside bucket radius

beforeN=size(smartLocs,1);
smartOut=smartLocs;
if(hookHeight)
    smartOut=smartOut(smartOut(:,3)>hookHeight,:);
end
D=squareform(pdist(smartOut)); %square matrix
distz=sum(D,1);
[v,bottInds]=maxk(distz,topN);
smartOut(bottInds,:)=[];
if(brad)
r=sqrt(sum(smartOut(:,1:2).^2,2));
outside=find(r>brad);
    if(any(outside))
%         warning('found smarts outside bucketrad')
    end
    smartOut(outside,:)=[];

end
afterN=size(smartOut,1);

end

