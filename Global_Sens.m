% Pairwise Correlations of Parameters for my chosen topology
% Manually in .txtbc change the Ara input

model = IQMmodel('emptymodel21.txtbc'); 
RoundedICs = round(IQMinitialconditions(model),4);
model = IQMinitialconditions(model,RoundedICs);


model1 = IQMparameters(model,{'ARA'},0);
model2 = IQMparameters(model,{'ARA'},0.003);
model3 = IQMparameters(model,{'ARA'},0.008);
model4 = IQMparameters(model,{'ARA'},0.025);
model5 = IQMparameters(model,{'ARA'},0.075);
model6 = IQMparameters(model,{'ARA'},0.226);
model7 = IQMparameters(model,{'ARA'},0.679);
model8 = IQMparameters(model,{'ARA'},2.031);
model9 = IQMparameters(model,{'ARA'},6.088);
model10 = IQMparameters(model,{'ARA'},18.25);
model11 = IQMparameters(model,{'ARA'},54.82);
model12 = IQMparameters(model,{'ARA'},164.524);
model13 = IQMparameters(model,{'ARA'},493.572);
model14 = IQMparameters(model,{'ARA'},1478.718);
model15 = IQMparameters(model,{'ARA'},4442.816);
model16 = IQMparameters(model,{'ARA'},13321.788);

models=[model1 model2 model3 model4 model5 model6 model7 model8 model9 model10 model11 model12 model13 model14 model15 model16];
parameters = IQMparameters(model);

%Since some are fixed (no range) and I don't want to test these because
%global analysis sampling is risky
parameters = setdiff(parameters, {'ARA'},'stable');
parameters = setdiff(parameters, {'kMAPIA'},'stable');
parameters = setdiff(parameters, {'kMC'},'stable');
parameters = setdiff(parameters, {'kMB'},'stable');
parameters = setdiff(parameters, {'CasConcaseF'},'stable');
parameters = setdiff(parameters, {'kMB3'},'stable');
parameters = setdiff(parameters, {'kMC4'},'stable');
parameters = setdiff(parameters, {'kMC5'},'stable');
parameters = setdiff(parameters, {'kMB6'},'stable');


time = [0:2:900];

Ns=[200000];

%% First Order Sobol

tmpx1FO=cell(length(models),1);
tmpx2FO=cell(length(models),1);

parfor nAra=1:length(models)
    
    model=models(nAra);
    
    tabled_sobol_info_fo=cell(length(parameters)+2,length(Ns)+1);
    tabled_sobol_info_fo{1,1}={'Ns'};
    tabled_sobol_info_fo{2,1}={'sums'};

    for i=1:length(parameters)

        tabled_sobol_info_fo{i+2,1} = parameters{i};

    end

    k=2;
    descending_per_Ns_fo = cell2table(cell(1+length(parameters),length(Ns)), 'VariableNames', {'Ns1'});%Change Ns1, Ns2 to include all Ns used


    for m=Ns
        OPTIONS=struct('statenames',{{'Pc'}},'range',1/2,'firstorder', 1, 'Nsim', m);
        %OPTIONS.Nsim=m;
        res_glob_sobolfo=IQMsensglobalsobol(model,time,parameters,OPTIONS);
        empty_cell_fo=cell(length(res_glob_sobolfo.parameters),2);

        for n=1:length(res_glob_sobolfo.parameters)
        empty_cell_fo{n,1}=res_glob_sobolfo.parameters(n);
        end

        for j=1:length(res_glob_sobolfo.parameters)
        empty_cell_fo{j,2}=res_glob_sobolfo.singlemeasure(j); 
        end

        ranked_global_firstorder=sortrows(empty_cell_fo,[2],'descend'); %generalise for all Ns
        ranked_global_firstorder=cell2table(ranked_global_firstorder, 'VariableNames', {'ranked_params','sobol_fo_values'});
        ranked_global_firstorder.Properties.VariableNames=["m","ranked_global_first_order"];
        descending_per_Ns_fo(:,k-1)=[table(m,'VariableNames',{'m'});ranked_global_firstorder(:,1)];

        %%%%

        tabled_sobol_info_fo{1,k}={Ns(k-1)};
        tabled_sobol_info_fo{2,k}={sum([empty_cell_fo{:,2}])};

        for l=3:length(parameters)+2
        tabled_sobol_info_fo{l,k}={empty_cell_fo{l-2,2}};
        end

        k=k+1;

    end

    tabled_sobol_info_fo=cell2table(tabled_sobol_info_fo, 'VariableNames', {'params','Ns1'});
    
    %SAVE tabled_sobol_info_fo and ranked_global_firstorder for the
    %specific Ara
     

tmpx1FO(nAra)={tabled_sobol_info_fo};
tmpx2FO(nAra)={ranked_global_firstorder};

end



%% Total Effect Sobol

tmpx1T=cell(length(models),1);
tmpx2T=cell(length(models),1);

parfor nAra=1:length(models)

    model=models(nAra);

    tabled_sobol_info_te=cell(length(parameters)+2, length(Ns)+1);
    tabled_sobol_info_te{1,1}={'Ns'};
    tabled_sobol_info_te{2,1}={'sums'};
    
    for i=1:length(parameters)
        
        tabled_sobol_info_te{i+2,1} = parameters{i};
    
    end
    
    %tabled_sobol_info_te=tabled_sobol_info_te';
    
    k=2;
    %Change Ns1, Ns2 to include all Ns used
    descending_per_Ns_te = cell2table(cell(1+length(parameters),length(Ns)), 'VariableNames', {'Ns1','Ns2', 'Ns3', 'Ns4'});
    
    for m=Ns
        OPTIONS=struct('statenames',{{'Pc'}},'range',1/2,'firstorder', 0, 'Nsim', m);
        res_glob_sobolte=IQMsensglobalsobol(model,time,parameters,OPTIONS);
        empty_cell_te=cell(length(res_glob_sobolte.parameters),2);
    
        for n=1:length(res_glob_sobolte.parameters)
        empty_cell_te{n,1}=res_glob_sobolte.parameters(n);
        end
        
        for j=1:length(res_glob_sobolte.parameters)
        empty_cell_te{j,2}=res_glob_sobolte.singlemeasure(j); 
        end
    
        ranked_global_totaleffect=sortrows(empty_cell_te,[2],'descend'); %generalise for all Ns
        ranked_global_totaleffect=cell2table(ranked_global_totaleffect, 'VariableNames', {'ranked_params','sobol_fo_values'});
        ranked_global_totaleffect.Properties.VariableNames=["m","ranked_global_totaleffect"];
        descending_per_Ns_te(:,k-1)=[table(m,'VariableNames',{'m'});ranked_global_totaleffect(:,1)];
        
        
        
        tabled_sobol_info_te{1,k}={Ns(k-1)};
        tabled_sobol_info_te{2,k}={sum([empty_cell_te{:,2}])};
        
        for l=3:length(parameters)+2
        tabled_sobol_info_te{l,k}={empty_cell_te{l-2,2}};
        end
          
        k=k+1;
        
    end
    
    tabled_sobol_info_te=cell2table(tabled_sobol_info_te, 'VariableNames', {'params','Ns1','Ns2', 'Ns3', 'Ns4'});
    
    tmpx1T(nAra)={tabled_sobol_info_te};
    tmpx2T(nAra)={ranked_global_totaleffect};
end
