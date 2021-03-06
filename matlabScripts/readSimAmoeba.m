% simD=importdata('A:\SmarticleRun\Amoeba\5 rob\f_1.0_v_1\PostProcess\RingPos.txt');
% function readSimAmoeba()
% d=simD.data;
% fh = @eig;
closeWaitbar;
clear all
center=1; %1 to center ring pos data at (0,0)
fold=(uigetdir('A:\SmarticleRun'));
f = dir2(fold,'folders');



binW=5; %bin width for polar plots 
simAm=struct;
h = waitbar(0,'Please wait...');
set(h.Children.Title,'interpreter','none');
steps = length(f);
for i =1:length(f)
    waitbar(i/steps,h,{['Processing: ',num2str(i),'/',num2str(length(f))],f(i).name})
    d=fullfile(fold,f(i).name,'PostProcess','RingPos.txt');
    simD=importdata(d);
%     tstep, x, y, z, globalGUI, comX, comY, comZ
        %changing coordinate system
        %x->-x
        %y->y
        %z->z
    simD.data(:,2)=-simD.data(:,2);
    
    data=[simD.data(:,1:3)];
    
    if(size(simD.data,2)>5)
        simD.data(:,6)=-simD.data(:,6);
        COM=[simD.data(:,1) simD.data(:,6:7)];
    end
    if(center)
        data(:,2)=data(:,2) -data(1,2);
        data(:,3)=data(:,3) -data(1,3);
        if(size(simD.data,2)>5)
            COM(:,2)=COM(:,2)-COM(1,2);
            COM(:,3)=COM(:,3)-COM(1,3);
        end
    end
    
    if(size(simD.data,2)>5)
        simAm(i).COM={COM};
    end
%     simAm(i).
    data={data};
    simAm(i).data=data;
    
    simAm(i).name=f(i).name;
    simAm(i).fps=diff(data{1}(1:2,1));
    [~,vals]=parseFileNames(f(i).name);
    simAm(i).pars=vals;
    %get radius of circle
    [match, nomatch]=regexp(simD.textdata,'[0-9.]+','match');
    simAm(i).r=str2double(match{1});
    if(any(simD.data(100:end,5)~=0))
        pts(simAm(i).name);
        error('file with gui non-zero');
    end
    if(all(diff(simD.data(:,1))~=simAm(i).fps))
        pts(simAm(i).name);
        error('file with bad timestep');
    end
    
%     if(exist(fullfile(fold,f(i).name,'PostProcess','RingContact.txt'),'file')==2)
%         simAm(i).binW=binW;
%         [simAm(i).innerForce,simAm(i).contactAngs,simAm(i).fcs]=readContactDistrib(...
%             fullfile(fold,f(i).name,'PostProcess','RingContact.txt'),binW);
%     end

    if(exist(fullfile(fold,f(i).name,'PostProcess','InactivePos.txt'),'file')==2)
        [simAm(i).deadInnerForce,simAm(i).deadContactAngs, simAm(i).deadPos,simAm(i).deadRot,simAm(i).deadFcs,simAm(i).ringpos]=readInactive(...
        fullfile(fold,f(i).name,'PostProcess','InactivePos.txt'),binW);
    end

    if(exist(fullfile(fold,f(i).name,'PostProcess','RingDead.txt'),'file')==2)
        [simAm(i).fullT,simAm(i).fullRingPos,simAm(i).fullDeadSmartPos,~]=readDeadRing(...
        fullfile(fold,f(i).name,'PostProcess','RingDead.txt'));
    end
    
    if(exist(fullfile(fold,f(i).name,'PostProcess','Stress.txt'),'file')==2)
        [simAm(i).AllSmartPos]=readAllSmarticlesPos(...
        fullfile(fold,f(i).name,'PostProcess','Stress.txt'));
        %          readAllSmarticlesPos(...
%         fullfile(fold,f(i).name,'PostProcess','Stress.txt'));
    end
    
    %rename video files to the name of the system
    ff= fullfile(fold,f(i).name,'video_capture');
    if(exist(fullfile(fold,f(i).name,'video_capture','outvid.avi'),'file')==2)
        movefile(fullfile(ff,'outvid.avi'),fullfile(ff,[f(i).name,'.avi']))
    end
end
closeWaitbar;
save(fullfile(fold,'amoebaData.mat'),'fold','simAm','-v7.3')
% end