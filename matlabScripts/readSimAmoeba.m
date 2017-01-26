% simD=importdata('A:\SmarticleRun\Amoeba\5 rob\f_1.0_v_1\PostProcess\RingPos.txt');
% d=simD.data;
% fh = @eig;
clear all
center=1; %1 to center ring pos data at (0,0)
fold=(uigetdir('A:\SmarticleRun'));
f = dir2(fold,'folders');

simAm=struct;

for i =1:length(f);
    d=fullfile(fold,f(i).name,'PostProcess','RingPos.txt');
    simD=importdata(d);
    simD.data(:,3)=-simD.data(:,3);
    
    data=[simD.data(:,1:3)];
    
    if(size(simD.data,2)>5)
        simD.data(:,7)=-simD.data(:,7);
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
    
    if(exist(fullfile(fold,f(i).name,'PostProcess','RingContact.txt'), 'file')==1)
    
    end
end

save(fullfile(fold,'amoebaData.mat'),'fold','simAm')
