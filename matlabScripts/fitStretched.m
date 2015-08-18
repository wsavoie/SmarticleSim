function [out]=fitStretched2(parsin,c,xdata,ydata,PON)
%parsin
%consts = (gamma)
%xdat
%ydat
%pon

global gam;
gam = c;
options = optimset('TolX',1e-5,'TolFun',1e-8,'MaxFunEvals',1e7,'MaxIter',1e5);
ParsOut =fminsearch(@stretchedExpFit,parsin,options);
% err = fiterr2(parsin,@stretchedExpFit,xdata,ydata,0);
    function [sse, FittedCurve] = stretchedExpFit(parsin)
%         %a(1)=delta a(2)=beta
%         tau = 1/(30)*exp(a(1)/gam);
%         curv = exp(-(q./tau).^a(2));
          FittedCurve = exp(-(xdata./(1/(30)*exp(parsin(1)/gam))).^parsin(2));
%             FittedCurve = exp(-(xdata./(1/(30)*exp(parsin(1)/gam)))).^parsin(2);
          ErrorVector = FittedCurve-ydata;
          sse= sum(ErrorVector.^2);
    end
if PON
    pts('delta=', ParsOut(1), ' beta=',ParsOut(2));
% [ParsOut err]
end
% out=[ParsOut,err];
out = ParsOut;
end