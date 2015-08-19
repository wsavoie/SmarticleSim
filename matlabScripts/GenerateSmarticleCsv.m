function []=GenerateSmarticleCsv(dt,omega, torqueThresh, angLow, angHigh)
% GENERATESMARTICLECSV generates file containing all instructions
% for the different movement types, global, GUI, overTorque, etc.
%   dt          = timestep to be used in simulation
%   omega       = angular velocity of arms to be used in simulation
%   torqeThresh = threshold at which overtorque relief movement begins
%   angLow      = lowest angle smarticle is allowed to move to
%   angHigh     = highest angle smarticle is allowed to move to
%top of file will include dT, omega, torqueThresh, angLow, angHigh
directory_name = uigetdir('D:\SimResults\Chrono\SmarticleU\tests');
fileloc = horzcat(directory_name,'\','smarticleMoves.csv')
fid = fopen(fileloc,'wt');
%add in all initial values to top of file
%square torque so that we can use the square in the program so we can avoid
%having to use sqrts in the code at each timestep!
fprintf(fid,'%s\n%f\n%f\n%f\n%f\n',dt,omega,torqueThresh^2,angLow,angHigh);
%add #1 to denote end of first section
fprintf(fid,'#\n');

%% global function position definition
% define some positions in the angular phase space (TO BE CHANGED)
ss= dt*omega; %step size
gait= 2;
figure(1);
hold on;
switch gait
    case 1% circle gait
        ang = 0:ss:2*pi;
        r=pi/2;
        phi = 0;
        global_theta_1Pos=r*cos(ang-phi);
        global_theta_2Pos=r*sin(ang-phi);
        plot(global_theta_1Pos,global_theta_2Pos,'.k');
        xlabel('\theta_1');
        ylabel('\theta_2');
        axis square
        title('Circle Gait');
        figText(gcf,15);
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
        
        global_theta_1Pos=[t1,l1,b1,r1];
        global_theta_2Pos=[t2,l2,b2,r2];
        plot(global_theta_1Pos,global_theta_2Pos,'.k');
        xlabel('\theta_1');
        ylabel('\theta_2');
        axis([-sL sL -sL sL])
        axis square
        title('Square Gait');
        figText(gcf,15);
end

%convert the distance between the points to the proper amount
%using dt and omega

if(length(global_theta_1Pos)~=length(global_theta_2Pos))
    error('global positions have inequal lengths for the arm position arrays!');
end

for i=1:length(global_theta_1Pos)
    fprintf(fid,'%f, %f, \n',global_theta_1Pos(i),global_theta_2Pos(i));
end
%add #2 to denote end of 2nd section
% fprintf(fid,'#2\n');

% %% over-torque function position definition
% % define some positions in the angular phase space (TO BE CHANGED)
% OT_theta_1Pos = sin(dt*omega+20);
% OT_theta_2Pos = sin(dt*omega+20);
% 
% %convert the distance between the points to the proper amount
% %using dt and omega
% 
% if(length(OT_theta_1Pos)~=length(OT_theta_2Pos))
%     error('global positions have inequal lengths for the arm position arrays!');
% end
% 
% for i=1:length(OT_theta_1Pos)
%     fprintf(fid,'%f\t%f\n',OT_theta_1Pos(i),OT_theta_2Pos(i));
% end
% %add #3 to denote end of 3rd section
% fprintf(fid,'#3\n');
% 
% %% GUI function position definition
% % define some positions in the angular phase space (TO BE CHANGED)
% GUI_theta_1Pos = sin(dt*omega+20);
% GUI_theta_2Pos = sin(dt*omega+20);
% 
% %convert the distance between the points to the proper amount
% %using dt and omega
% 
% if(length(GUI_theta_1Pos)~=length(GUI_theta_2Pos))
%     error('global positions have inequal lengths for the arm position arrays!');
% end
% 
% for i=1:length(GUI_theta_1Pos)
%     fprintf(fid,'%f\t%f\n',GUI_theta_1Pos(i),GUI_theta_2Pos(i));
% end
% %add #4 to denote end of 4th section
% fprintf(fid,'#4\n');

fclose('all');