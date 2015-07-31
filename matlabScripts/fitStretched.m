function [out]=fitStretched(parsin,c,xdat,ydat,PON)
%parsin
%consts = (gamma)
%xdat
%ydat
%pon
global gam;
gam = c;
ParsOut =fminsearch('fiterr2',parsin,[],@stretchedExpFit,xdat,ydat,PON); %tau' tau''
err = fiterr2(parsin,@stretchedExpFit,xdat,ydat,0);
    function curv = stretchedExpFit(a,q)    %a(1) is tau' %a(2) is tau''
        %a(1)=delta a(2)=gamma a(3)=beta
        %a(1)=delta a(3)=beta
        curv = exp(-(q./(1/30*exp(a(1)/gam)) ).^a(2));
    end
out=[ParsOut,err];
end