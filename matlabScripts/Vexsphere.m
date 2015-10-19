

N=0;
vr=2;
vSpace= 4/3*pi*vr^3;
vSphere= 4/3*pi*(1)^3;
its=1000000;
% axis equal
% axis([-vr,vr,-vr,vr,-vr,vr]);
% % axis equal
% [sx,sy,sz]=sphere;
 rc=rand(its,3)*2*(vr-1)-(vr-1);
 tic
for i=1:its
%     surf(sx,sy,sz);
%     axis equal
%     axis([-vr,vr,-vr,vr,-vr,vr])
%     
%     hold on;
% 
%     surf(sx+rc(i,1),sy+rc(i,2),sz+rc(i,3));
%         pause(.1)
%     hold off;
    d=norm(rc(i,:));
    if(d<=1)
        N=N+1;
    else
       continue; 
    end
end
toc
p=its/N;
phi=(vSphere*2)/vSpace
p=its/N*1/phi