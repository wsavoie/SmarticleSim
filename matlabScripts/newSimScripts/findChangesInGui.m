function [guis] = findChangesInGuid(gui)
%FINDCHANGESINGUID list index and gui value of input 
%   guis=[index,val]
    guis=find(diff(gui))+1; %add 1 because of diff
    guis=vertcat(1,guis); %add index 1 for gui value at start of sim
    guiVals=gui(guis);
    guis=[guis,guiVals];
end

