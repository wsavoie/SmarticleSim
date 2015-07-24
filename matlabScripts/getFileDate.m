function [date] = getFileDate(file,varargin)
%takes in folder or filename and pulls off date from first characters
%varargin 1=string,
%varargin 0(or nothing)=double
pars = {1};
if(isstruct(file))
    fname = file.name;
else
    fname = file;
end
pattern = '[0-9]*';
pars= regexp(fname, pattern, 'match');


%date should always be length of 8
if length(pars{1})==8
    if nargin>1 %if true then varargin=true
        if varargin{1}
        date=pars(1);
        return;
        end
    end
    
    date = str2double(pars{1});
    return;
end
date=[];
    
    
