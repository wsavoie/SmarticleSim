function [smartPos] = singleRunSmartPoints(dat,verz)
%SINGLERUNSMARTPOINTS Summary of this function goes here
%   Detailed explanation goes here
frs=dat(verz).smartInfo(:,:,1:end); %[x,y,z,e0,e1,e2,e3,ang1,ang2]
t1=dat(verz).smartSize(2);
t2=dat(verz).smartSize(3);
l=dat(verz).smartSize(4);
w=dat(verz).smartSize(5);
smartPos=zeros(size(frs,3),size(frs(:,:,1),1),4,3);
for j=1:size(frs,3)
    d=frs(:,:,j); %data for current frame
    
    %final timestep is d
    %     [e0,e1,e2,e3]=separateVec(d(:,4:7),1); %get quaternion data
    e0=d(:,4); e1=d(:,5); e2=d(:,6); e3=d(:,7);
    l2r=[e2,-e3,-e0,-e1];%w,x,y,z->-y,z,w,x
    
    %     q = [d(:,4:7)];
    %     q=quaternion(l2r); %I need to transpose it before setting to quat for matrix
    %     angles = EulerAngles(q,'XYZ');
    %     angles2 = quat2eul(l2r);
    
    
    %     angles=reshape(angles(3,1,:),[size(angles,3),1]);'
    sp=zeros(4,3,size(d,1));
    for k=1:size(d,1)
        
        sm=d(k,:);
        smPos=sm(1:3);
        smPos(1)=-smPos(1);
        lmod=l + t2;
        %         RR=RotationMatrix(q(j));
        RR=quat2rotm(l2r(k,:));
        %scatter3(smPos(1),smPos(2),smPos(3),'o');
        
        %v1   v4
        %|    |
        %v2-0-v3
        
        yy=[0,t1/2,0]*RR;
        v1=([-w/2-lmod/2*cos(sm(:,8)),0,-(l-t2)*sin(sm(:,8))])*RR+smPos;
        v2=([-w/2,0,0])*RR+smPos;
        v3=([w/2,0,0])*RR+smPos;
        v4=([w/2+lmod/2*cos(sm(:,9)),0,-(l-t2)*sin(sm(:,9))])*RR+smPos;
        
        ptz=[v1;v2;v3;v4];
        smartPos(j,k,:,:)=ptz;
    end
    
end