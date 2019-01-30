
% tau = 10;
% K=11;
% syms s K tau u U
% A= [0 K; 0 -tau^-1];
% B= [0; tau^-1];
% C= [1 0 ; 0 1];
% D =[0; 0];
% 
% phi = (s*C-A)^-1
% 
% Y=C*phi*B*U

% Y/u=  [ K/(s*(s*tau + 1))] y(s)
%       [1/(s*tau + 1)]      v(s)
s=tf('s');

P1 = 1/(s^2 + 10*s + 20);
