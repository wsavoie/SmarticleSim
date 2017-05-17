function [hout]=logslopeline(pow,p,x,y)
%makes a loglog slope line of power pow and starting at point p from end
%pow = power of line
%p = # from end where line starts
%x = xdata
%y = ydata
%s = size of 
%hin = handle in;
%hout = handle to return
xdat= [x(end-p):10.^(log10(x(end-p))):x(end-p)*10];
% x
% xdat
yline = linspace(1,10,length(xdat));
% ydat = ((xdat).^pow).*y(end-p);
ydat = y(end-p)*yline.^(pow);
yline
hout = loglog(xdat,ydat,'--','linewidth',2);
% yline = linspace(1,10,length(xdat));
% ydat= y(end-p);
% hout = loglog(xdat,ydat*yline.^(pow));
%%%%%%%%%%%


%    OutM(:,1:2) = OutMatrix1; % [freq,G']
%     OutM(:,3) = OutMatrix2(:,2); %[freq,G',G'']
%     OutM(:,4) = OutMatrix3(:,2); %[freq,G',G'',std(G')]
%     OutM(:,5) = OutMatrix4(:,2); %[freq,G',G'',std(G'),std(G'')]
%     % h = loglog(out(:,2),out(:,3),'ro-',out(:,2),out(:,4),'bo-' );
%     % h = loglog(out(:,2)*rad/gapheight,out(:,3),'ro-',out(:,2)*rad/gapheight,out(:,4),'bo-' );
%     
%     
%     p1 = 4;% point from end to use for log line placement and power
%     p2 = 4;% point from end to use for log line placement and power
%     pow1 = 4/10;
%     pow2 = 1/2;
%     
%     % gg1=[OutM(end-p,1):1e-3:OutM(end-p,1)*10]; % makes a center point of a loglog line point p
%     % gg2=[OutM(end-p2,1):1e-3:OutM(end-p2,1)*10]; % makes a center point of a loglog line at point p
%     gg1=[OutM(end-p1,1):1e-3:OutM(end-p1,1)*10]; % makes a center point of a loglog line at point p
%     gg2=[OutM(end-p2,1):1e-3:OutM(end-p2,1)*10]; % makes a center point of a loglog line at point p
%     
%     gg3 = linspace(1,10,length(gg1));
%     gg4 = linspace(1,10,length(gg2));
%     
%     gpCent= OutM(end-p1,2); % center of last 5 # of G' at correct order of mag
%     gppCent=OutM(end-p2,3); % center of last 5 # of G''at correct order of mag
% loglog(gg1,gpCent*gg3.^(pow1),'r',...
%         gg2,gppCent*gg4.^(pow2),'b');
end