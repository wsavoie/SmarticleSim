function [smartCross] = generatePackingFromSimDatPLANE(dat,countDubs)
%[smartCross]=GENERATEPACKINGFROMSIMDAT(dat,startInd,countDubs)
%
%dat is a 1x1 struct (single version
%countDubs= 1 or 0, flag to count double plane crossings as entanglements

TOL=1e-8;
PLOTZ=0;
% if(PLOTZ)
%     figure(129);
%     hold on;
% end
% startInd=5;
% vInd=1;
sinVal=.99; %represents ~85 degs 
%%%%
% warning('use single struct ind as input param');
% c=dat(vInd).smartInfo(:,:,startInd:end); %[x,y,z,e0,e1,e2,e3,ang1,ang2]
c=dat(1).smartInfo(:,:,1:end); %[x,y,z,e0,e1,e2,e3,ang1,ang2]
% [rho,t1,t2,l,w]=separateVec(dat(1).smartSize,1);
t1=dat(1).smartSize(2);
t2=dat(1).smartSize(3);
l=dat(1).smartSize(4);
w=dat(1).smartSize(5);

smartCross=zeros(size(c,3),size(c,1));
for i=1:size(c,3)
    d=c(:,:,i); %data for current frame
    
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
    for j=1:size(d,1)
        
        sm=d(j,:);
        smPos=sm(1:3);
        smPos(1)=-smPos(1);
        lmod=l + t2;
        %         RR=RotationMatrix(q(j));
        RR=quat2rotm(l2r(j,:));
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
        
        %         if(PLOTZ)
        %             patch(ptz(:,1),ptz(:,2),ptz(:,3),[1,0,0]);
        %         end
        %         patch(ptz1(:,1),ptz1(:,2),ptz1(:,3),[1,0,0]);
        %         patch(ptz2(:,1),ptz2(:,2),ptz2(:,3),[1,0,0]);
        
        sp(:,:,j)=ptz;
        %
        %         drawEdge3d(createEdge3d(v2,v3));
        %         drawEdge3d(createEdge3d(v1,v2));
        %         drawEdge3d(createEdge3d(v4,v3));
        %         axis equal
    end
    
    for j=1:size(d,1)
        polyg=sp(:,:,j);
        %         if(PLOTZ)
        %             pH=drawPolygon3d(polyg,'color',[.2,.5,.3],'linewidth',2);
        %         end
        %
        count=0;
        %filter automatically when ang is nearly straight <5deg
        %sin 85deg=.9962
        skipz=pdist([polyg(1,:);polyg(4,:)])>2*l*sinVal+w;
        if(~skipz)
            
            for k=1:size(d,1)
                if(k~=j)
                    armPts=sp(:,:,k);
                    
                    L1=[armPts(2,:);armPts(2,:);armPts(3,:)];
                    L2=[armPts(1,:);armPts(3,:);armPts(4,:)];
                    L=createLine3d(L1,L2);
                    E=createEdge3d(L1,L2);
                    %                 l1=createLine3d(armPts(2,:),armPts(1,:)); %left barb
                    %                 l2=createLine3d(armPts(2,:),armPts(3,:)); %spine
                    %                 l3=createLine3d(armPts(3,:),armPts(4,:)); %right barb
                    %
                    %                 if(PLOTZ)
                    %                     h=drawEdge3d(E,'linewidth',2);
                    %                     g=drawLine3d(L,'linewidth',1,'color',[.8,.2,.3]);
                    %                 end
                    %intersectionPt
                    
                    skipI=pdist([armPts(1,:);armPts(4,:)])>2*l*sinVal+w;
                    if(~skipI)
                        [iPt,inside]=intersectRayPolygon3d(L,polyg);
                        %straight shapes are not accurate
                        
                        if(any(inside))
                            %                     if(sum(inside)>1)
                            %                         warning('multiple inside')
                            %                         pH=drawPolygon3d(polyg,'color',[.2,.5,.3],'linewidth',2);
                            %                         h=drawEdge3d(E,'linewidth',2);
                            %                         g=drawLine3d(L,'linewidth',1,'color',[.8,.2,.3]);
                            %                     end
                            [r,~,~]=find(inside);
                            
                            %point may on actual edge only on line
                            %so check that it exists on the actual edge
                            if((length(r)>1 && countDubs) || length(r)==1)
                                distz=distancePointEdge3d(iPt(r,:),E(r,:));
                                if(distz<TOL)
                                    %                             if(PLOTZ)
                                    %                                 scatter3(iPt(r,1),iPt(r,2),iPt(r,3),60*ones(size(ptArr,1),1),'filled')
                                    %                             end
                                    count=count+length(r);
                                end
                            end
                        end
                    end
                    
                    %                 if(PLOTZ)
                    %                     delete(h);
                    %                     delete(g);
                    %                 end
                end
            end
        end
        smartCross(i,j)=count;
        %         if(PLOTZ)
        %             delete(pH);
        %         end
    end
end
% beep;