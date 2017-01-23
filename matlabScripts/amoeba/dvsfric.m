load('A:\SmarticleRun\Amoeba_COG_1_dead_90_deg\amoebaData.mat');
cols = get(gca,'colororder');

col=cols(2,:);
filled =0;


figure(1);
hold on;
xx=[.1:.1:.7];
for k=1:length(xx);
    SPACE_UNITS='m';
    TIME_UNITS='s';
    ma = msdanalyzer(2, SPACE_UNITS, TIME_UNITS);
    f=[xx(k)]; rob=[5]; v=[];
    useCOM=1;
props={f rob v};
for i=1:length(simAm)
    inds=1;
    cond=true;
    for j=1:length(props)
        %if empty accept all values
        if ~isempty(props{j})
            %in case multiple numbers in property
            %if no matches then cond = false
%             simAm(i).pars(j)
            if(~any(round(props{j},2)==round(simAm(i).pars(j),2)))
                cond = false;
            end
        end
    end
    if(cond)
        if(useCOM)
            if(isempty(intersect(fieldnames(simAm),'COM')))
                error('trying to use COM when data does not have it');
            end
            ma = ma.addAll(simAm(i).COM);
        else
            ma = ma.addAll(simAm(i).data);
        end
        
%             usedSimAm(inds)=simAm(i);
        inds=inds+1;
    end
end
    
    ma = ma.computeMSD;
    % ma.plotMeanMSD(gca, true);
    % ma.plotMSD;
    
    dmult=1/4;
    [a b]=ma.fitMeanMSD;
    D(k)=dmult*a.p1;
    % E(k)=gof.rmse;
    % rr(k)=gof.adjrsquare;
end
% errorbar(xx,D,E)
if (filled)
    plot(xx,D,'o-','linewidth',2,'color',col,'markerfacecolor',col);
else
    plot(xx,D,'o-','linewidth',2,'color',col);
end
if useCOM
    title('D vs. \mu_s (COM)');
else
    title('D vs. \mu_s (COG)');
end
xlabel('static coefficient of friction \mu_s');
ylabel('measured diffusion constant D (m^2/s)');