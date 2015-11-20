%% figure 4a and b
%b
delta = [6.5 9 11 12 15 14 13 11 9.5];
lw = [0 .13 .16 .29 .4 .55 .66 1.02 1.17];
figure(10);
hold on;
plot(lw,delta,'k.','MarkerSize',25,'MarkerEdgeColor','k');
xlabel('l/w');
ylabel('\Delta');
axis([0 1.3 0 18])
figure(12)
hold on;
plot(lw,exp(delta),'k.-','MarkerSize',25,'MarkerEdgeColor','k');
xlabel('l/w');
ylabel('exp(\Delta)');

%a
%if t_o = 1, tt/tt_0 = tt
gamma=[1.23 1.48 1.7 1.96 2.20 2.5];

f=1/30;

delta=15;
figure(4)
tau=log(f)+delta./gamma;
hold on;
plot(1./gamma,tau,'g.-','MarkerSize',25,'MarkerEdgeColor','g');

xlabel('\Gamma^{-1}')
ylabel('ln(\tau/\tau_0)+const')
axis([.1 1.1 0 30])