clear all;
close all;
figure(1);
horz=26.6;
vert=35;
N=400;%number of smarts
box = [-1 -1 1 1 -1; -1 1 1 -1,-1];
box(1,:)=box(1,:)*horz;
box(2,:)=box(2,:)*vert+(vert-horz);
smartBase = [0 0 3 3; 6 0 0 6; 1 1 1 1];


smarts=ones(3,N*4);
for(i=1:N)
hold on;
R=zeros(3,3); R(3,3)=1;
r=rotz(360*rand);
t=[40*rand; 40*rand]-20;
R(1:2,1:2)=r(1:2,1:2);
R(1:2,3)=t;
smart=R*smartBase;

plot(smart(1,:),smart(2,:),'o-')
smarts(:,4*(i-1)+1:4*i)=[smart(1:2,:);i*ones(1,4)];
end
plot(box(1,:)',box(2,:),'k.-')
axis equal

smarts=smarts';
rsmarts=round(smarts,1);
smartsSorted=sortrows(rsmarts,[1,2]);
distMin=1;
m=-horz;
wind=6.1;
idx=2;

h1=figure(1);
h2=figure(2);
objects=allchild(h1);
copyobj(get(h1,'children'),h2);

PON=1;
while(m+distMin<horz)
    

    s=smartsSorted(smartsSorted(:,1)>m & smartsSorted(:,1)<(m+wind),:)';
    if(PON)
   
%     plot([s(1,1),s(2,1)],[s(1,end),s(2,end)],'-g');
    figure(2);
    line([m m],[-20,-20],'color','k')
    line([m+wind m+wind],[20,-20],'color','k')
    end
    
    if(isempty(s))
        m=m+.1;
        continue;
    end
    [mx,ind]=max(s(2,:));
    m=s(1,ind);
    l(1,idx)=s(1,ind);
    l(2,idx)=mx;  
    idx=idx+1;

end
l(:,1)=[-horz;l(2,2)];
l(:,idx)=[horz;l(2,end)];
figure(1);
plot(l(1,:),l(2,:),'-','linewidth',2);

% finding top line of smarts



% 
% R1=zeros(3,3); R1(3,3)=1;
% [R1,R2,R3]=deal(R1);
% r1=rotz(360*rand); r2=rotz(360*rand); r3=rotz(360*rand);
% t1=[3;-1]; t2=[6;5]; t3=[0;0];
% 
% R1(1:2,1:2)=r1(1:2,1:2);
% R2(1:2,1:2)=r2(1:2,1:2);
% R3(1:2,1:2)=r3(1:2,1:2);
% 
% R1(1:2,3)=t1;
% R2(1:2,3)=t2;
% R3(1:2,3)=t3;
% 
% % smart2=smarts0; smarts2(3,)=[3,-1,1];
% 
% 
% 
% 
% 
% smarts1= [R1*smartBase];
% smarts2= [R2*smartBase];
% smarts3= [R3*smartBase];
% 
% 
% figure(1);
% hold on;
% % fill(box(:,1),box(:,2),'r')
% 
% plot(smarts1(1,:),smarts1(2,:),'o-')
% plot(smarts2(1,:),smarts2(2,:),'o-')
% plot(smarts3(1,:),smarts3(2,:),'o-')
% plot(box(1,:)',box(2,:),'k.-')


% A=polyarea([xy(:,1);xy(1,1)],[xy(:,2);xy(1,2)]); 