function [out]=fitStretched(parsin,c,xdata,ydata,PON)
%parsin
%consts = (gamma)
%xdat
%ydat
%pon

global gam;
gam = c;
% options = optimset('TolX',1e-5);
ParsOut =fminsearch('fiterr2',parsin,[],@stretchedExpFit,xdata,ydata,PON);
% err = fiterr2(parsin,@stretchedExpFit,xdata,ydata,0);
    function curv = stretchedExpFit(a,q)
%         %a(1)=delta a(2)=beta
%         tau = 1/(30)*exp(a(1)/gam);
%         curv = exp(-(q./tau).^a(2));

           curv = exp(-(q./(1/(30)*exp(a(1)/gam))).^a(2));
%         curv = 0.2 - (0.2 - 1)*exp(-(q./a(1)).^a(2));
    end
if PON
    pts('delta=', ParsOut(1), ' beta=',ParsOut(2));
% [ParsOut err]
end
% out=[ParsOut,err];
out = [ParsOut];
end