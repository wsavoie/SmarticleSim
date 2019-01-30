function [] = getAllSmartPoints(D)
%GETALLSMARTPOINTS creates matrix with all smarticle end points at all
%times and puts it into a datafile
% dat input is from vibDAta

%from
posOut=struct;
for(i=1:length(D))
    c=fileparts(D(i).fold);
    if exist(fullfile(c,'ballDat.mat'),'file')
        load(fullfile(c,'ballDat.mat')); %called dat
    else
        error('no ballDat file exists');
    end
    posOut(i).startTime=startTime;
    verz=D(i).v;
    verzID=find([dat(:).v]==verz);
%     posOut(i).ballOut=dat(verz).ballOut;
    posOut(i).ballOut=dat(verzID).ballOut;
    posOut(i).smartPos=singleRunSmartPoints(dat,verzID);
% % %     frs=dat(verz).smartInfo(:,:,1:end); %[x,y,z,e0,e1,e2,e3,ang1,ang2]
% % %     % [rho,t1,t2,l,w]=separateVec(dat(1).smartSize,1);
% % %     t1=dat(verz).smartSize(2);
% % %     t2=dat(verz).smartSize(3);
% % %     l=dat(verz).smartSize(4);
% % %     w=dat(verz).smartSize(5);
% % %     posOut(i).smartPos=zeros(size(frs,3),size(frs(:,:,1),1),4,3);
% % %     for j=1:size(frs,3)
% % %         d=frs(:,:,j); %data for current frame
% % %         
% % %         %final timestep is d
% % %         %     [e0,e1,e2,e3]=separateVec(d(:,4:7),1); %get quaternion data
% % %         e0=d(:,4); e1=d(:,5); e2=d(:,6); e3=d(:,7);
% % %         l2r=[e2,-e3,-e0,-e1];%w,x,y,z->-y,z,w,x
% % %         
% % %         %     q = [d(:,4:7)];
% % %         %     q=quaternion(l2r); %I need to transpose it before setting to quat for matrix
% % %         %     angles = EulerAngles(q,'XYZ');
% % %         %     angles2 = quat2eul(l2r);
% % %         
% % %         
% % %         %     angles=reshape(angles(3,1,:),[size(angles,3),1]);'
% % %         sp=zeros(4,3,size(d,1));
% % %         for k=1:size(d,1)
% % %             
% % %             sm=d(k,:);
% % %             smPos=sm(1:3);
% % %             smPos(1)=-smPos(1);
% % %             lmod=l + t2;
% % %             %         RR=RotationMatrix(q(j));
% % %             RR=quat2rotm(l2r(k,:));
% % %             %scatter3(smPos(1),smPos(2),smPos(3),'o');
% % %             
% % %             %v1   v4
% % %             %|    |
% % %             %v2-0-v3
% % %             
% % %             yy=[0,t1/2,0]*RR;
% % %             v1=([-w/2-lmod/2*cos(sm(:,8)),0,-(l-t2)*sin(sm(:,8))])*RR+smPos;
% % %             v2=([-w/2,0,0])*RR+smPos;
% % %             v3=([w/2,0,0])*RR+smPos;
% % %             v4=([w/2+lmod/2*cos(sm(:,9)),0,-(l-t2)*sin(sm(:,9))])*RR+smPos;
% % %             
% % %             ptz=[v1;v2;v3;v4];
% % %             posOut(i).smartPos(j,k,:,:)=ptz;  
% % %         end
% % % 
% % %     end
    posOut(i).smartSize=dat(verzID).smartSize;
end
save(fullfile(fileparts(fileparts(D(1).fold)),'posOut.mat'),'posOut');