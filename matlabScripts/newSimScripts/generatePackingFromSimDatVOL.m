% function [] = generatePackingFromSimDat(dat)
%GENERATEPACKINGFROMSIMDAT
%%%
smartCross=0;
startInd=5;
vInd=1;
hold on;
%%%%
warning('use single struct ind as input param');
c=dat(vInd).smartInfo(:,:,startInd:end); %[x,y,z,e0,e1,e2,e3,ang1,ang2]
[rho,t1,t2,l,w]=separateVec(dat(vInd).smartSize,1);


for(i=size(c,3))
    d=c(:,:,i);
    
    %final timestep is d
    [e0 e1 e2 e3]=separateVec(d(:,4:7),1);
    l2r=[e2,-e3,-e0,-e1];%w,x,y,z->-y,z,w,x
    
    %     q = [d(:,4:7)];
    q=quaternion(l2r); %I need to transpose it before setting to quat for matrix
    angles = EulerAngles(q,'XYZ');
    %     angles=reshape(angles(3,1,:),[size(angles,3),1]);'
    sp=zeros(4,3,size(d,1));
    for(j=1:size(d,1))
        
        %     [e0 e1 e2 e3]=separateVec(d(:,4:7),1);
        %     l2r=[e2,-e3,-e0,-e1];%w,x,y,z->-y,z,w,x
        %     q=quaternion(l2r); %I need to transpose it before setting to quat for matrix
        
        
        sm=d(j,:);
        lmod=l + t2;
        smPos=sm(1:3);
        smPos(1)=-smPos(1);
        
        RR=rotx(angles(1,j))*roty(angles(2,j))*rotz(angles(3,j));
        RR=RotationMatrix(q(j));
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
        
        ptz1=ptz+ones(4,1)*yy;
        ptz2=ptz-ones(4,1)*yy;
        aptz=[ptz1;ptz2];
        
        %         patch(ptz(:,1),ptz(:,2),ptz(:,3),[1,0,0]);
        %     patch(ptz1(:,1),ptz1(:,2),ptz1(:,3),[1,0,0]);
        %         patch(ptz2(:,1),ptz2(:,2),ptz2(:,3),[1,0,0]);
        sp(:,:,j)=ptz;
        %
        %         drawEdge3d(createEdge3d(v2,v3));
        %         drawEdge3d(createEdge3d(v1,v2));
        %         drawEdge3d(createEdge3d(v4,v3));
        %         axis equal
    end
    
    for(j=1:size(d,1))
        
        RR=rotx(angles(1,j))*roty(angles(2,j))*rotz(angles(3,j));
        RR=RotationMatrix(q(j));
        hh=[0,t1/2,0]*RR;
        ww=[t2,0,0]*RR;
        ww2=[0,0,t2/2]*RR;
        
        intShape=[sp(:,:,j)+ones(4,1)*hh;sp(:,:,j)-ones(4,1)*hh];
        %     shp=alphaShape(intShape(:,1),intShape(:,2),intShape(:,3));
        %     plot(shp);
        
        
        
        
        %     shp=alphaShape(np(:,1),np(:,2),np(:,3))
        %
        %     yy=[0,t2/2,0]*RR;
        %
        %
        %     %     plane=createPlane([p2;p3;p1]);
        %
        %     %     drawPlane3d(createPlane(p1,p3,p2),'facecolor',[0,0,.6]);
        % %     sp=
        %     e1=createEdge3d(sp(2,:,j),sp(1,:,j)); %left barb
        %     e2=createEdge3d(sp(2,:,j),sp(3,:,j)); %spine
        %     e3=createEdge3d(sp(3,:,j),sp(4,:,j)); %right barb
        %     %     drawPlane3d(plane,'facecolor',[.2 .4 .1]);
        %     drawEdge3d(e1);
        %     drawEdge3d(e2);
        %     drawEdge3d(e3);
        %
        %     %     plane
        %     createPlane([p0;p1;p2])
        %     drawPlane3d(plane);
        count=0;
        
        %         ishape=alphaShape(intShape(:,1),intShape(:,2),intShape(:,3));
        %         plot(ishape);
        for(k=1:size(d,1))
            %             pts('hi');
            if(k~=j)
                ff=0;%make sure to not allow double counting intersections
                
                RR=rotx(angles(1,k))*roty(angles(2,k))*rotz(angles(3,k));
                RR=RotationMatrix(q(k));
                hh=[0,t1/2,0]*RR;
                ww=[t2,0,0]*RR;
                p1=sp(1,:,k);
                p2=sp(2,:,k);
                p3=sp(3,:,k);
                p4=sp(4,:,k);
                
                
                if(pdist([intShape(1,:);p1],'euclidean')>2*w)
                    continue;
                end
                
                %%%right barb
                rb=[p4+hh;p3+hh;p4-hh;p3-hh];
                rb=[rb;rb+ones(4,1)*ww];
                %             shp1=alphaShape(rb(:,1),rb(:,2),rb(:,3));
                %             plot(shp1);
                
                %%%left barb
                lb=[p1+hh;p2+hh;p1-hh;p2-hh];
                lb=[lb;lb-ones(4,1)*ww];
                %             shp2=alphaShape(lb(:,1),lb(:,2),lb(:,3));
                %             plot(shp2);
                
                %%%center barb
                ww2=[0,0,t2/2]*RR;
                cb=[p2+hh-ww;p3+hh+ww;p2-hh-ww;p3-hh+ww];
                cb=[cb-ones(4,1)*ww2;cb+ones(4,1)*ww2];
                %                 shp3=alphaShape(cb(:,1),cb(:,2),cb(:,3));
                %                 plot(shp3);
                
                Ir=intersectionHull('vert',intShape,'vert',rb);
                if(Ir.vert)
                    count=count+1;
                    ff=1;
                    %                     shp=alphaShape(Ir.vert(:,1),Ir.vert(:,2),Ir.vert(:,3));
                    %                     plot(shp,'facecolor',[0.6,0.2,0.1]);
                    continue;
                end
                
                Il=intersectionHull('vert',intShape,'vert',lb);
                if(Il.vert & ff==0)
                    count=count+1;
                    ff=1;
                    %                     shp=alphaShape(Il.vert(:,1),Il.vert(:,2),Il.vert(:,3));
                    %                     plot(shp,'facecolor',[0.6,0.2,0.1]);
                    continue;
                end
                
                Ic=intersectionHull('vert',intShape,'vert',cb);
                if(Ic.vert & ff==0)
                    ff=1;
                    count=count+1;
                    %                     shp=alphaShape(Ic.vert(:,1),Ic.vert(:,2),Ic.vert(:,3));
                    %                     plot(shp,'facecolor',[0.6,0.2,0.1]);
                    continue;
                end
            end
        end
        smartCross(1,j)=count;
    end
end
pts('end');