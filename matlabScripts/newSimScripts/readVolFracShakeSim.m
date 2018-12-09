% function [out]=readShake(fold)
fold='B:\SmartSimResults\11_26\LWvary';
    c=dir2(fold,'folders');
for i=1:length(c)
    d=regexp(c(i).name,['\d*.\d*'],'match');
    d=abs(str2double(d));
    lw(i)=d(1);
    n(i)=d(2)*d(3);
    inFolds=dir2(fullfile(fold,c(i).name),'folders');
    for j=1:length(inFolds)
        datafilePath=fullfile(inFolds(j).folder,inFolds(j).name,'PostProcess','volumeFraction.txt');
        fdat=importdata(datafilePath);
        phi(j,i)=fdat(end,3);
    end
end
hold on;
phiM=mean(phi);
ro=.0234;%bucket radius without thickness taken into account in sim
rn=0.0224;%corrected bucket radius
%bucket rad was declared wrong in sim this fixes it
phiM=phiM*(ro/rn)^2;
errorbar(lw,phiM,-std(phi)/2,std(phi)/2,0.0714*ones(1,length(lw)),0.0714*ones(1,length(lw)),'linewidth',2);
lwBS=[0, 0.1 0.2 0.3 0.4 0.6 0.7 1.3 1.4];
phiBS=[.227,0.2,0.18 0.16 .1425,.115,.105,0.061,.06];
    nickD=[[0:.1:1.4]',[.27 .235 .208 .185 .163 .148 .136 .13 .118 .11 .105 .1025 .095 .09 .085]'];
    plot(lwBS,phiBS,'o-','linewidth',2);
    plot(nickD(:,1),nickD(:,2),'ro-','linewidth',2);
    xlabel('l/w');
    ylabel('\phi');
    axis([0 1.4 0.05 0.3]);
    figText(gcf,16);
legend({'chrono data (\phi_i) l''=l\pmD','Nick mean (\phi_i)','Nick mean \phi_f'})
