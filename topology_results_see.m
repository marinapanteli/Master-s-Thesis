table_paramSpecs=readtable('parametersTHESIS_after_arabinose.txt'); % Its content is not used, as can be seen below

k=1;
[~,idxcostmin]=min(viablePoints.cost);
[~,idxparams]=ismember(viablePoints.colnames,table_paramSpecs.names);
table_paramSpecs.p0(idxparams)=[10.^viablePoints.rowmat(idxcostmin,:)];
[r1,~]=ismember(1:50,idxparams);
[~,r2]=find(r1==0);
table_paramSpecs.p0(r2)=zeros(1,length(r2));
dokim=table_paramSpecs.p0;
params=dokim;
[~,inputAraPos] = ismember('ARA',DATA.paramNames);
params(inputAraPos) = 0;

% [~,inputkfGCSG3Pos] = ismember('kfGCSG3',DATA.paramNames);
% [~,inputkbGCSG3Pos] = ismember('kbGCSG3',DATA.paramNames);
% [~,inputkMCPICnPos] = ismember('kMCPIC',DATA.paramNames);
% [~,inputkfSGCASsg3Pos] = ismember('kfSGCASsg3',DATA.paramNames);
% [~,inputkbSGCASsg3Pos] = ismember('kbSGCASsg3',DATA.paramNames);
% [~,inputAranPos] = ismember('ARAN',DATA.paramNames);

%[~,inputAranPos] = ismember('kMC',DATA.paramNames);




      
%2% (similarly for 1%, same for 10%)
%%%For 2% relative perturbation of the top 6 params that are the most sensitive
%%%locally from  Ara 0.075 to Ara 6.088 (where the green stripe curves most)
            params_plus=params;
            params_min=params;
            
%            params_plus(inputkfGCSG3Pos) = 1.02*params(inputkfGCSG3Pos);
%            params_plus(inputkbGCSG3Pos) = 1.02*params(inputkbGCSG3Pos);
%             params_plus(inputkMCPICnPos) = 1.02*params(inputkMCPICnPos);
%             params_plus(inputkfSGCASsg3Pos) = 1.02*params(inputkfSGCASsg3Pos);
%             params_plus(inputkbSGCASsg3Pos) = 1.02*params(inputkbSGCASsg3Pos);
             params_plus(inputAranPos) = 1.02*params(inputAranPos);
            
%            params_min(inputkfGCSG3Pos) = 0.98*params(inputkfGCSG3Pos);
%             params_min(inputkbGCSG3Pos) = 0.98*params(inputkbGCSG3Pos);
%             params_min(inputkMCPICnPos) = 0.98*params(inputkMCPICnPos);
%             params_min(inputkfSGCASsg3Pos) = 0.98*params(inputkfSGCASsg3Pos);
%             params_min(inputkbSGCASsg3Pos) = 0.98*params(inputkbSGCASsg3Pos);
             params_min(inputAranPos) = 0.98*params(inputAranPos);
            

            ICs = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params);
            ICs_plus = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params_plus);
            ICs_min = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params_min);

            [~,FluoPos] = ismember({'Pa','Pb','Pc'},ICs.states);     
             for a = 1:16
                     dokim(inputAraPos) = DATA.exp_data{a,1};
                     params_plus(inputAraPos) = DATA.exp_data{a,1};
                     params_min(inputAraPos) = DATA.exp_data{a,1};

                     simulation = IQMPsimulate(DATA.mex_model,[0:10:900], ICs.statevalues(end,:), DATA.paramNames, dokim);
                     simulation_plus = IQMPsimulate(DATA.mex_model,[0:10:900], ICs_plus.statevalues(end,:), DATA.paramNames, params_plus);
                     simulation_min = IQMPsimulate(DATA.mex_model,[0:10:900], ICs_min.statevalues(end,:), DATA.paramNames, params_min);
                     
                     result900(a,:) =  simulation.statevalues(end,FluoPos);
                     result900plus(a,:) =  simulation_plus.statevalues(end,FluoPos);
                     result900minus(a,:) =  simulation_min.statevalues(end,FluoPos);
             end
                    
            
             
            stdval = [zeros(5,1)+0.2;zeros(6,1)+0.1;zeros(5,1)+0.2];
            out = sumsqr((result900./max(result900) - DATA.exp_data{:,2:4})./stdval);
            result900=flip(result900);
            result900minus=flip(result900minus);
            result900plus=flip(result900plus);

             figure()
        plot(result900(:,1)./max(result900(:,1)),'--y')
        hold on
        plot(result900(:,2)./max(result900(:,2)),'--b')
        plot(result900(:,3)./max(result900(:,3)),'--g')
       
        plot(result900plus(:,3)./max(result900plus(:,3)),'--r')
        plot(result900minus(:,3)./max(result900minus(:,3)),'--k')

  
        errorbar([1:16],flip(DATA.exp_data{:,2}),stdval,'-y');
        errorbar([1:16],flip(DATA.exp_data{:,3}),stdval,'-b');
        errorbar([1:16],flip(DATA.exp_data{:,4}),stdval,'-g');
      

%10%
%%%For 10% relative perturbation of the top 6 params that are the most sensitive
%%%locally from  Ara 0.075 to Ara 6.088 (where the green stripe curves most)
            params_plus=params;
            params_min=params;
            
            params_plus(inputkfGCSG3Pos) = 1.1*params(inputkfGCSG3Pos);
%            params_plus(inputkbGCSG3Pos) = 1.1*params(inputkbGCSG3Pos);
%             params_plus(inputkMCPICnPos) = 1.1*params(inputkMCPICnPos);
%             params_plus(inputkfSGCASsg3Pos) = 1.1*params(inputkfSGCASsg3Pos);
%             params_plus(inputkbSGCASsg3Pos) = 1.1*params(inputkbSGCASsg3Pos);
%             params_plus(inputAranPos) = 1.1*params(inputAranPos);
            
            params_min(inputkfGCSG3Pos) = 0.9*params(inputkfGCSG3Pos);
%             params_min(inputkbGCSG3Pos) = 0.9*params(inputkbGCSG3Pos);
%             params_min(inputkMCPICnPos) = 0.9*params(inputkMCPICnPos);
%             params_min(inputkfSGCASsg3Pos) = 0.9*params(inputkfSGCASsg3Pos);
%             params_min(inputkbSGCASsg3Pos) = 0.9*params(inputkbSGCASsg3Pos);
%             params_min(inputAranPos) = 0.9*params(inputAranPos);
            

            ICs = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params);
            ICs_plus = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params_plus);
            ICs_min = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params_min);

            [~,FluoPos] = ismember({'Pa','Pb','Pc'},ICs.states);     
             for a = 1:16
                     dokim(inputAraPos) = DATA.exp_data{a,1};
                     params_plus(inputAraPos) = DATA.exp_data{a,1};
                     params_min(inputAraPos) = DATA.exp_data{a,1};

                     simulation = IQMPsimulate(DATA.mex_model,[0:10:900], ICs.statevalues(end,:), DATA.paramNames, dokim);
                     simulation_plus = IQMPsimulate(DATA.mex_model,[0:10:900], ICs_plus.statevalues(end,:), DATA.paramNames, params_plus);
                     simulation_min = IQMPsimulate(DATA.mex_model,[0:10:900], ICs_min.statevalues(end,:), DATA.paramNames, params_min);
                     
                     result900(a,:) =  simulation.statevalues(end,FluoPos);
                     result900plus(a,:) =  simulation_plus.statevalues(end,FluoPos);
                     result900minus(a,:) =  simulation_min.statevalues(end,FluoPos);
             end
                    
            
             
            stdval = [zeros(5,1)+0.2;zeros(6,1)+0.1;zeros(5,1)+0.2];
            out = sumsqr((result900./max(result900) - DATA.exp_data{:,2:4})./stdval);
            result900=flip(result900);
            result900minus=flip(result900minus);
            result900plus=flip(result900plus);

             figure()
        plot(result900(:,1)./max(result900(:,1)),'--y')
        hold on
        plot(result900(:,2)./max(result900(:,2)),'--b')
        plot(result900(:,3)./max(result900(:,3)),'--g')
       
        plot(result900plus(:,3)./max(result900plus(:,3)),'--r')
        plot(result900minus(:,3)./max(result900minus(:,3)),'--k')

  
        errorbar([1:16],flip(DATA.exp_data{:,2}),stdval,'-y');
        errorbar([1:16],flip(DATA.exp_data{:,3}),stdval,'-b');
        errorbar([1:16],flip(DATA.exp_data{:,4}),stdval,'-g');
      