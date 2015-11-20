function [shapeLines]=getShapeLines(time,guid)


cols = {[1,0,0],[1,.5,0],[.4431,.7373, 1],[0,0,0],[.392,.824,.118],[.7,.4,.7]};

hold on;

gc =[0; find(diff(guid))];
gc=gc+1;
val=[guid(gc)]+1;%added value can be zero (global) and matrices are 1 started
% sp1=[plot(time,zcom)]; %starts variable for the legend

gg = [gc, val];
gg= sortrows(gg,2);
figure(1);
shapeLines=zeros(length(unique(guid)),5);
for i=1:size(gg,1)
    if i~=1
        if gg(i,2)==gg(i-1,2) %if true will mean a repeat in color on legend
            shapeLines(i,:)=[time(gg(i,1)) cell2mat(cols(gg(i,2))),gg(i,2)+1];
            continue;
        else
            shapeLines(i,:)=[time(gg(i,1)) cell2mat(cols(gg(i,2))),gg(i,2)+1];
        end
    else
        shapeLines(i,:)=[time(gg(i,1)) cell2mat(cols(gg(i,2))),gg(i,2)+1];
    end
    % plotTypes=[plotTypes, (gg(i,2)+1)];
end

end