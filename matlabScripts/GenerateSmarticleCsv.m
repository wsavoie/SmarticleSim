%function []=GenerateSmarticleCsv(dt,omega, torqueThresh, angLow, angHigh)
% GENERATESMARTICLECSV generates file containing all instructions
% for the different movement types, global, GUI, overTorque, etc.
%   dt          = timestep to be used in simulation
%   omega       = angular velocity of arms to be used in simulation
%   torqueLimit = torque limit for motor
%   angLow      = lowest angle smarticle is allowed to move to
%   angHigh     = highest angle smarticle is allowed to move to
%top of file will include dT, omega, torqueThresh, angLow, angHigh

%https://www.pololu.com/product/1040/specs
%theoretical:
%servo speed - .12s/60deg (strange units) 8.72 rads/sec
%stall torque- 0.8 kgF·cm= .078 nm=
%stall curr-   270 mA

%measured: 
%measured 5.0N with radius of 6.52mm
%T=5*.00652=.0325
%dimensionless parameter for torque and length=
%for robot: T=.008, L=.05 arm weight = 11.3 g
%rads = 8.7267
%running curr- 100mA
%Kt = T/I, .008/.27 = .0296 Nm/A
%T=0.0588399
stapleSize = false;
dt=.001; %.00025
sizeScale=1;
% omega = 4.9244e-5;
%%limit speed in sim (5.3 for 90 deg, 6.3 for 180 deg, takes time to slow so probably 7
%does nothing on torque limit
omegaLim = 6.3; %if you want to change the speed you must change linear interpolate method
omega = 6.3; %distance between points in move list 
% omega = 10;
rho = 7850.0;%/(sizeScale^3);
%(t2_smarticle) * (t_smarticle)* (w_smarticle + 2 * (l_smarticle));
if stapleSize
    t   = .00127*sizeScale;
    t2  = .0005*sizeScale;
    w_s = .0117*sizeScale;
    l_s = w_s;
    rho = 7850.0;%/(sizeScale^3);
    volume =  t2 * t* (w_s + 2 * (l_s));
    mass = volume*rho;
    torqueLimit=2*9.8*mass*w_s;%.00005;  4.6657e-04
else
      
    t= .029982; %t = height of solar panels
    t2 = .02122;
    w_s = 0.05316;
  
    l_s = .8*w_s; 
%     t   = .0079*sizeScale;
%     t2  = .0053*sizeScale;
%     w_s = .0117*sizeScale;
%     l_s = w_s;
%1.5424 = lw
    rho = 443.0;%/(sizeScale^3);
    volume =  t2 * t* (w_s + 2 * (l_s));
    mass = volume*rho;
%     torqueThresh=.0325; %.008 previously cd 
    torqueLimit=.025; %..0325 previously cd 
end





angLow=-120;
angHigh=120;

global_gait= 1;
gui1_gait = 1;
gui2_gait = 1;
gui3_gait = 2;
midt_gait = 2;
extra1_gait=1;
extra2_gait=1;
torqueLimit
PON= 1;

user=getenv('username');
if strcmp(user,'root')
    directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests');
else
    directory_name = uigetdir('D:\GT Coursework\smarticledata'); 
end

fileloc = horzcat(directory_name,'\','smarticleMoves.csv');
fid = fopen(fileloc,'wt');
%add in all initial values to top of file
%square torque so that we can use the square in the program so we can avoid
%having to use sqrts in the code at each timestep!
fprintf(fid,'%s\n%f\n%f\n%f\n%f\n%f\n',dt,omega,torqueLimit,angLow,angHigh,omegaLim);
%add #1 to denote end of first section
fprintf(fid,'#\n');

%% global function position definition
% define some positions in the angular phase space (TO BE CHANGED)
% ss= 2.75*pi/180; %step size save value as in linear interpolate method
radI=(pi/2); %intial radius to base arclength off of
ss=radI*(.4)*(pi/180); %rad*degs*degToRads   this is arc length each step

%want to make sure arclength/step is >> than rads/step
%i.e  ss>>(omega*dt)

switch global_gait
    case 1% circle gait
%         ang = 0:ss:2*pi;
        r=1.6;
        arcDivs=(2*r*pi)/(ss); %circumference divided by arclength to get number of points
        ang=linspace(0,2*pi,arcDivs);
        ang(end)=[]; %remove 2pi index
        phi = 0;
        offset1=0;
        offset2=0;
        global_alpha_1Pos=r*cos(ang-phi)+offset1;
        global_alpha_2Pos=r*sin(ang-phi)+offset2;
        if PON
            hold on;
            figure(1);
            plot(global_alpha_1Pos,global_alpha_2Pos,'.k');
            xlabel('\alpha_1');
            ylabel('\alpha_2');
            
            title('Circle Gait');
            figText(gcf,15);
            axis equal
        end
    case 2% square gait
        %not implemented yet!
        sL = pi/2; %square gait side length
        %theta1
        t1 = sL/2:-ss:-sL/2;
        l1 = -sL/2*ones(1,length(t1));
        b1 = -sL/2:ss:sL/2;
        r1 = sL/2*ones(1,length(t1));
%       theta2
        t2 = sL/2*ones(1,length(t1));
        l2 = sL/2:-ss:-sL/2;
        b2 = -sL/2*ones(1,length(t1));
        r2 = -sL/2:ss:sL/2;
        
        global_alpha_1Pos=[t1,l1,b1,r1];
        global_alpha_2Pos=[t2,l2,b2,r2];
        if PON
            figure(1);
            hold on;
            plot(global_alpha_1Pos,global_alpha_2Pos,'.');
            xlabel('\alpha_1');
            ylabel('\alpha_2');
            axis([-sL sL -sL sL])
            axis square
            title('Square Gait');
            figText(gcf,15);
        end
    case 3% line gait
%         ang = 0;
        ang=0*pi/180;
        r=1.6;

        centerx = 0;
        centery = 0;
        arcDivs=(r*2)/(ss); %circumference divided by arclength to get number of points
%         ang=linspace(0,2*pi,arcDivs);
        lx1=linspace(-r,r,arcDivs)*cos(ang);
        ly1=linspace(-r,r,arcDivs)*sin(ang);
        lx2=linspace(r,-r,arcDivs)*cos(ang);
        ly2=linspace(r,-r,arcDivs)*sin(ang);
        lx2(1)=[]; ly2(1)=[]; %remove doubled first pos
        lx2(end)=[]; ly2(end)=[]; %remove doubled last pos
%         a1_0= centerx+r*cos(ang);
%         a2_0= centery+r*sin(ang);
%         a1_f= centerx-r*cos(ang);
%         a2_f= centery-r*sin(ang);
                
        lx = [lx1 lx2];
        ly = [ly1 ly2];
        global_alpha_1Pos= lx;
        global_alpha_2Pos= ly;
        if PON
            figure(1);
            hold on;
            plot(global_alpha_1Pos,global_alpha_2Pos,'.k');
            xlabel('\alpha_1');
            ylabel('\alpha_2');
            axis([-2*pi/3 2*pi/3 -2*pi/3 2*pi/3])
            axis square
%             title({horzcat('\alpha_1 = ',num2str(a1_0),'\rightarrow ',num2str(a1_f),' \rightarrow ',num2str(a1_0)),...
%            horzcat('\alpha_2 = ',num2str(ly(),'\rightarrow ',num2str(a2_f),'\rightarrow ',num2str(a2_0))});
            title(['line gait length=',num2str(r),' angle=',num2str(ang*180/pi)]);
            figText(gcf,15);
        end
      case 4% rectified square gait
        %not implemented yet!
        sL = pi/2; %square gait side length
        len = length(0:ss:sL); %gait length
        sL2 = pi/2;%*.8;
        %theta1
        amove1 = linspace(0,sL,len);    %0   ->  max
        amove2 = linspace(sL,sL,len);   %max ->  max
        amove3 = linspace(sL,0,len);    %max ->  0
        amove4 = linspace(0,0,len);     %0   ->  0
%       theta2
        bmove1 = linspace(0,0,len);     %0   ->  0
        bmove2 = linspace(0,sL2,len);   %0   ->  max
        bmove3 = linspace(sL2,sL2,len); %max ->  max
        bmove4 = linspace(sL2,0,len);   %max ->  0
        
        global_alpha_1Pos=[amove1,amove2,amove3,amove4];
        global_alpha_2Pos=[bmove1,bmove2,bmove3,bmove4];
        if PON
            figure(1);
            hold on;
            plot(global_alpha_1Pos,global_alpha_2Pos,'.k');
            xlabel('\alpha_1');
            ylabel('\alpha_2');
            axis([-sL 2*sL -sL 2*sL])
            axis equal
            title('Rectified Square Gait');
            figText(gcf,15);
        end  
    case 5
       
        figure(1);
        clf;
        hold on;
        axis([-1.25*pi/2,1.25*pi/2,-1.25*pi/2,1.25*pi/2]);
%         imr=imread('C:\Users\root\Desktop\geom\tgran.PNG');
%         im2=imagesc(imr);
%         image(linspace(-1.25*pi/2,1.25*pi/2,im2.XData(2)),linspace(1.25*pi/2,-1.25*pi/2,im2.YData(2)),imr)
        
        
        
        axis square
        
        xlabel('\alpha_1');
        ylabel('\alpha_2');
        %%%%%%%%
        
        hold on
%         r=pi/4;
%         th = 0:pi/50:2*pi;
%         x=pi/4;
%         y=pi/4;
%         xunit = r * cos(th) + x;
%         yunit = r * sin(th) + y;
%         plot(xunit, yunit,'k');
%         plot(-xunit,-yunit,'k');
%         plot(xunit, -yunit,'k');
%         plot(-xunit, yunit,'k');
%         x=-pi/4;
%         y=-pi/4;
%         x=[-1.75,1.75];
%         y=[-1.75,1.75];
%         plot(x,y,'k');
%         plot(x,-y,'k');
%         plot([-1.75,1.75],[0,0],'k');
%         plot([0,0],[-1.75,1.75],'k');
        

        hold off
        %%%%%%%%%%
        h=imfreehand(gca);
       
        gaitPts= getPosition(h);
        gaitPts(end+1,:)= gaitPts(1,:);
        [global_alpha_1Pos,global_alpha_2Pos] = ...
            interpolateGait(gaitPts,ss);
        
        hold on;
        title('Custom Drawn Gait');
        plot(global_alpha_1Pos,global_alpha_2Pos,'o-k');
        xlabel('\alpha_1');
        ylabel('\alpha_2');
        axis([-1.25*pi/2,1.25*pi/2,-1.25*pi/2,1.25*pi/2]);
        axis square
        figText(gcf,15);
        
end

%convert the distance between the points to the proper amount
%using dt and omega

if(length(global_alpha_1Pos)~=length(global_alpha_2Pos))
    error('global positions have inequal lengths for the arm position arrays!');
end

for i=1:length(global_alpha_1Pos)
    fprintf(fid,'%f, %f',global_alpha_1Pos(i),global_alpha_2Pos(i));
    %denote next area with a # as last char
    if i == length(global_alpha_1Pos)
        fprintf(fid,'#\n',global_alpha_1Pos(i),global_alpha_2Pos(i));
    else
        fprintf(fid,', \n',global_alpha_1Pos(i),global_alpha_2Pos(i));
    end
end


%% GUI function position definitions

guiSize = 5;
for i=1:guiSize
    switch(i)
        case 1 %gui1
            switch(gui1_gait)
                case 1
                    GUI_alpha_1Pos = [pi/2];
                    GUI_alpha_2Pos = [pi/2];
                case 2
                    GUI_alpha_1Pos = [pi/4];
                    GUI_alpha_2Pos = [pi/4];
            end
        case 2 %gui2
            switch(gui2_gait)
                case 1
                    GUI_alpha_1Pos = [0];
                    GUI_alpha_2Pos = [0];
                case 2
                    GUI_alpha_1Pos = [2*pi/3];
                    GUI_alpha_2Pos = [2*pi/3];
            end
        case 3 %gui3
            switch(gui3_gait)
                case 1
                    GUI_alpha_1Pos = [pi/2];
                    GUI_alpha_2Pos = [-pi/2];
                case 2
                    GUI_alpha_1Pos = [-pi/2];
                    GUI_alpha_2Pos = [-pi/2];
                case 3
                    GUI_alpha_1Pos = [pi];
                    GUI_alpha_2Pos = [pi];
            end
        case 4 %extra1
           switch(extra1_gait)
               case 1
                   GUI_alpha_1Pos = [pi/4];
                   GUI_alpha_2Pos = [pi/4];
           end     
        case 5 %extra2
           switch(extra2_gait)
           case 1
                   GUI_alpha_1Pos = [-pi/4];
                   GUI_alpha_2Pos = [-pi/4];
           end     
    end
    
    if(length(GUI_alpha_1Pos)~=length(GUI_alpha_2Pos))
        error('gui positions have inequal lengths for the arm position arrays!');
    end
    
    for i=1:length(GUI_alpha_1Pos)
        fprintf(fid,'%f, %f',GUI_alpha_1Pos(i),GUI_alpha_2Pos(i));
        %denote next area with a # as last char
        if i == length(GUI_alpha_1Pos)
            fprintf(fid,'#\n',GUI_alpha_1Pos(i),GUI_alpha_2Pos(i));
        else
            fprintf(fid,', \n',GUI_alpha_1Pos(i),GUI_alpha_2Pos(i));
        end
    end
end


%% MIDT function position definitions

switch(midt_gait)
    case 1
        sL = pi/2;
        sL2 = pi/2;
        len = length(0:ss:sL); %gait length
        amove1 = linspace(0,sL,len);    %0   ->  max
        amove2 = linspace(sL,sL,len);   %max ->  max
        amove3 = linspace(sL,0,len);    %max ->  0
        amove4 = linspace(0,0,len);     %0   ->  0
%       theta2
        bmove1 = linspace(0,0,len);     %0   ->  0
        bmove2 = linspace(0,sL2,len);   %0   ->  max
        bmove3 = linspace(sL2,sL2,len); %max ->  max
        bmove4 = linspace(sL2,0,len);   %max ->  0
        
        MIDT_alpha_1Pos=[amove1,amove2,amove3,amove4];
        MIDT_alpha_2Pos=[bmove1,bmove2,bmove3,bmove4];
    case 2
        MIDT_alpha_1Pos = [pi/2];
        MIDT_alpha_2Pos = [-pi/2];
end          
   if(length(MIDT_alpha_1Pos)~=length(MIDT_alpha_2Pos))
        error('MIDT positions have inequal lengths for the arm position arrays!');
    end
for i=1:length(MIDT_alpha_1Pos)
    fprintf(fid,'%f, %f',MIDT_alpha_1Pos(i),MIDT_alpha_2Pos(i));
    %denote next area with a # as last char
    if i == length(MIDT_alpha_1Pos)
        fprintf(fid,'#\n',MIDT_alpha_1Pos(i),MIDT_alpha_2Pos(i));
    else
        fprintf(fid,', \n',MIDT_alpha_1Pos(i),MIDT_alpha_2Pos(i));
    end
end

% %% over-torque function position definition
% % define some positions in the angular phase space (TO BE CHANGED)
% OT_alpha_1Pos = sin(dt*omega+20);
% OT_alpha_2Pos = sin(dt*omega+20);
% 
% %convert the distance between the points to the proper amount
% %using dt and omega
% 
% if(length(OT_alpha_1Pos)~=length(OT_alpha_2Pos))
%     error('global positions have inequal lengths for the arm position arrays!');
% end
% 
% for i=1:length(OT_alpha_1Pos)
%     fprintf(fid,'%f\t%f\n',OT_alpha_1Pos(i),OT_alpha_2Pos(i));
% end
% %add #3 to denote end of 3rd section
% fprintf(fid,'#3\n');
% 
fclose('all');