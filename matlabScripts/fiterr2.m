function err=fiterr2(a,funct,x,data,PLOTTING)
err=norm((feval(funct,a,x)-data)).^2;
if PLOTTING==1
  figure(10);
vv=length(x);
%vv=1000;
plot((x(1:vv)),data(1:vv),'.');
hold on
plot((x(1:vv)),feval(funct,a,x(1:vv)),'g');
hold off
drawnow
%pause
end
% Fit error function, a=function parameters, funct=function
%                     x=function domain, data=real data to fit

%a
%err