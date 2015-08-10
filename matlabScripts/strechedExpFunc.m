% function curv = sinfit(a,q,c)
 function curv = strechedExpFunc(a,q)
%pars are the params are the ones that do not change 
%a is to be fitted ()
%a(1)=delta a(2)=gamma a(3)=beta
curv = exp(-(q./(1/30*exp(a(1)/a(2))) ).^a(3));


% a(1)*sin(a(2).*q+a(3));
