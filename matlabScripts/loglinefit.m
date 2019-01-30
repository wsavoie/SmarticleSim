function [xx,ydat]=loglinefit(p,x,y,varargin)
%makes a loglog slope line of power pow and starting at point p from end
%pow = power of line
%p = # from end where line starts
%x = xdata
%y = ydata
%s = size of 
%hin = handle in;
%hout = handle to return

pfit = polyfit(log(x(end-p:end-1)),log(y(end-p:end-1)),1);
m = pfit(1);
b = exp(pfit(2));
xdat= [x(end-p):10.^(log10(x(end-p))):x(end-p)*10];

d=log10(x(end-p));
xx=logspace(log10(x(end-p)),log10(x(end-1)));
c=abs(log10(y(end-p))-log10(y(end-1)));

ylin = logspace(log10(x(end-p))-d,log10(x(end-1))-d);
yy=logspace(0,c,length(xx));

% yline = linspace(1,10,length(xdat));
% ydat = ((xdat).^pow).*y(end-p);
ydat = y(end-p)*ylin.^(m);
% yy
% hout = loglog(xx,ydat,'--','linewidth',2);

if length(varargin)>1
    pts('slope is ',m);
end
end