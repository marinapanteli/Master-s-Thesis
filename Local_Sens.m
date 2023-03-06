%% Initialising DATA
costFunction = @myCost;
viabilityThreshold = chi2inv(0.95,3*16-41);
paramSpecs = 'parametersTHESIS_after_arabinose2.txt'; % its content does not matter, as can be seen below
folder_save = 'Save';
name_save = 'Saved'; %%%%
data_path = 'GeneratedData.csv'; %THIS
model_path = 'emptymodel2.txtbc'; %THIS
n_threads = 0;
par_profile = 'SGE_BSSE';

%Feasibility computation
 posteriorNames = {'ideal' 'uniform'};
 posteriors = {'ideal' 'uniform'};
 %tunableParamSettings = [-1 0 2 1];
tunableParamSettings = 2*ones(1,50);
%% Setup
% Get paramSpecs for futher use
table_paramSpecs = readtable(paramSpecs);
parallelize=1;
% Launch parallel pool if necessary
% parallelize = 0;
% if n_threads>0
%     parpool(par_profile, n_threads);
%     parallelize = 1;
% end
% Initialise cost function
DATA = myCost(0, data_path, model_path, paramSpecs);
%IQMmakeMEXmodel(DATA.model, DATA.mex_model);
%parallelInitPersistent(@myCost, 0, data_path, model_path, paramSpecs);


%% Local Sensitivity Analysis (! change perturbation in this section accordingly)
model = IQMmodel('emptymodel21.txtbc');

model = IQMparameters(model,{'ARA'},6.088); % I am not putting this in a loop to iterate over all ARA yet

OPTIONS.MaxIter=1000;
OPTIONS.Delta=1e-4;
outp=IQMsensdatastat(model,IQMparameters(IQMconvertNonNum2NumIC(model)),1*ones(1,length(IQMparameters(IQMconvertNonNum2NumIC(model)))),1*ones(1,length(IQMparameters(IQMconvertNonNum2NumIC(model)))), OPTIONS);
IQMsensstat(outp)

%% Parameters for my selected topology and the most sensitive for output Pc, which I decide in the section above

%! Change location of viable points file for your system
load('N:\homes\mpanteli\THESIS\Thesis\topodesign-2.0-master\examples\Marina\Save\13_2;)\Saved_viablePoints_sampled_2151768_00000000000000001100000000000000000000000000000111.mat')
table_paramSpecs=readtable('parametersTHESIS_after_arabinose.txt'); % This does not matter in content of values, as can be seen below in the script
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

[~,inputkfGCSG3Pos] = ismember('kfGCSG3',DATA.paramNames);
[~,inputkbGCSG3Pos] = ismember('kbGCSG3',DATA.paramNames);
[~,inputkMCPICPos] = ismember('kMCPIC',DATA.paramNames);
[~,inputkfSGCASsg3Pos] = ismember('kfSGCASsg3',DATA.paramNames);
[~,inputkbSGCASsg3Pos] = ismember('kbSGCASsg3',DATA.paramNames);
[~,inputAranPos] = ismember('ARAN',DATA.paramNames);

[~,inputkMCPos] = ismember('kMC',DATA.paramNames);
[~,inputkMB3Pos] = ismember('kMB3',DATA.paramNames);
[~,inputdmRNAcaseFPos] = ismember('dmRNAcaseF',DATA.paramNames);
[~,inputKmPos] = ismember('Km',DATA.paramNames);
[~,inputkfSGCASsg2Pos] = ismember('kfSGCASsg2',DATA.paramNames);
[~,inputkbSGCASsg2Pos] = ismember('kbSGCASsg2',DATA.paramNames);
[~,inputkfGBSG2Pos] = ismember('kfGBSG2',DATA.paramNames);
[~,inputkbGBSG2Pos] = ismember('kbGBSG2',DATA.paramNames);
[~,inputVMAXPos] = ismember('VMAX',DATA.paramNames);
[~,inputCasConcaseFPos] = ismember('CasConcaseF',DATA.paramNames);
[~,inputkfGCSG1Pos] = ismember('kfGCSG1',DATA.paramNames);
[~,inputkbGCSG1Pos] = ismember('kbGCSG1',DATA.paramNames);
[~,inputkfSGCASsg1Pos] = ismember('kfSGCASsg1',DATA.paramNames);
[~,inputkbSGCASsg1Pos] = ismember('kbSGCASsg1',DATA.paramNames);



%% 1% Visualization
%%%For 1% relative perturbation of the top 6+ params that are the most sensitive
%%%locally from  Ara 0.075 to Ara 6.088 (where the green stripe curves most)
            params_plus=params;
            params_min=params;
            
            params_plus(inputkfGCSG3Pos) = 1.01*params(inputkfGCSG3Pos);
%            params_plus(inputkbGCSG3Pos) = 1.01*params(inputkbGCSG3Pos);
%             params_plus(inputkMCPICPos) = 1.01*params(inputkMCPICPos);
%             params_plus(inputkfSGCASsg3Pos) = 1.01*params(inputkfSGCASsg3Pos);
%             params_plus(inputkbSGCASsg3Pos) = 1.01*params(inputkbSGCASsg3Pos);
%             params_plus(inputAranPos) = 1.01*params(inputAranPos);

%              params_plus(inputkMCPos) = 1.01*params(inputkMCPos);
%              params_plus(inputkMB3Pos) = 1.01*params(inputkMB3Pos);
%              params_plus(inputdmRNAcaseFPos) = 1.01*params(inputdmRNAcaseFPos);
%              params_plus(inputKmPos) = 1.01*params(inputKmPos);
%              params_plus(inputkfSGCASsg2Pos) = 1.01*params(inputkfSGCASsg2Pos);
%              params_plus(inputkbSGCASsg2Pos) = 1.01*params(inputkbSGCASsg2Pos);
%              params_plus(inputkfGBSG2Pos) = 1.01*params(inputkfGBSG2Pos);
%              params_plus(inputkbGBSG2Pos) = 1.01*params(inputkbGBSG2Pos);
%              params_plus(inputVMAXPos) = 1.01*params(inputVMAXPos);
%              params_plus(inputCasConcaseFPos) = 1.01*params(inputCasConcaseFPos);
%              params_plus(inputkfGCSG1Pos) = 1.01*params(inputkfGCSG1Pos);
%              params_plus(inputkbGCSG1Pos) = 1.01*params(inputkbGCSG1Pos);
%              params_plus(inputkfSGCASsg1Pos) = 1.01*params(inputkfSGCASsg1Pos);
%              params_plus(inputkbSGCASsg1Pos) = 1.01*params(inputkbSGCASsg1Pos);


            params_min(inputkfGCSG3Pos) = 0.99*params(inputkfGCSG3Pos);
%             params_min(inputkbGCSG3Pos) = 0.99*params(inputkbGCSG3Pos);
%             params_min(inputkMCPICPos) = 0.99*params(inputkMCPICPos);
%             params_min(inputkfSGCASsg3Pos) = 0.99*params(inputkfSGCASsg3Pos);
%             params_min(inputkbSGCASsg3Pos) = 0.99*params(inputkbSGCASsg3Pos);
%             params_min(inputAranPos) = 0.99*params(inputAranPos);

%              params_min(inputkMCPos) = 0.99*params(inputkMCPos);
%              params_min(inputkMB3Pos) = 0.99*params(inputkMB3Pos);
%              params_min(inputdmRNAcaseFPos) = 0.99*params(inputdmRNAcaseFPos);
%              params_min(inputKmPos) = 0.99*params(inputKmPos);
%              params_min(inputkfSGCASsg2Pos) = 0.99*params(inputkfSGCASsg2Pos);
%              params_min(inputkbSGCASsg2Pos) = 0.99*params(inputkbSGCASsg2Pos);
%              params_min(inputkfGBSG2Pos) = 0.99*params(inputkfGBSG2Pos);
%              params_min(inputkbGBSG2Pos) = 0.99*params(inputkbGBSG2Pos);
%              params_min(inputVMAXPos) = 0.99*params(inputVMAXPos);
%              params_min(inputCasConcaseFPos) = 0.99*params(inputCasConcaseFPos);
%              params_min(inputkfGCSG1Pos) = 0.99*params(inputkfGCSG1Pos);
%              params_min(inputkbGCSG1Pos) = 0.99*params(inputkbGCSG1Pos);
%              params_min(inputkfSGCASsg1Pos) = 0.99*params(inputkfSGCASsg1Pos);
%              params_min(inputkbSGCASsg1Pos) = 0.99*params(inputkbSGCASsg1Pos);
           

            ICs = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params);
            ICs_plus = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params_plus);
            ICs_min = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params_min);

            [~,FluoPos] = ismember({'Pa','Pb','Pc'},ICs.states);     
             for a = 1:16
                     dokim(inputAraPos) = DATA.exp_data{a,1};
                     params_plus(inputAraPos) = DATA.exp_data{a,1};
                     params_min(inputAraPos) = DATA.exp_data{a,1};

                     simulation = IQMPsimulate(DATA.mex_model,[0:1000:1000000], ICs.statevalues(end,:), DATA.paramNames, dokim);
                     simulation_plus = IQMPsimulate(DATA.mex_model,[0:1000:1000000], ICs_plus.statevalues(end,:), DATA.paramNames, params_plus);
                     simulation_min = IQMPsimulate(DATA.mex_model,[0:1000:1000000], ICs_min.statevalues(end,:), DATA.paramNames, params_min);
                     
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
      
%% 2% Visualization
%%%For 2% relative perturbation of the top 6+ params that are the most sensitive
%%%locally from  Ara 0.075 to Ara 6.088 (where the green stripe curves most)
            params_plus=params;
            params_min=params;
            
%            params_plus(inputkfGCSG3Pos) = 1.02*params(inputkfGCSG3Pos);
%            params_plus(inputkbGCSG3Pos) = 1.02*params(inputkbGCSG3Pos);
%             params_plus(inputkMCPICPos) = 1.02*params(inputkMCPICPos);
%             params_plus(inputkfSGCASsg3Pos) = 1.02*params(inputkfSGCASsg3Pos);
%             params_plus(inputkbSGCASsg3Pos) = 1.02*params(inputkbSGCASsg3Pos);
%             params_plus(inputAranPos) = 1.02*params(inputAranPos);
            
%              params_plus(inputkMCPos) = 1.02*params(inputkMCPos);
%              params_plus(inputkMB3Pos) = 1.02*params(inputkMB3Pos);
%              params_plus(inputdmRNAcaseFPos) = 1.02*params(inputdmRNAcaseFPos);
%              params_plus(inputKmPos) = 1.02*params(inputKmPos);
%              params_plus(inputkfSGCASsg2Pos) = 1.02*params(inputkfSGCASsg2Pos);
%              params_plus(inputkbSGCASsg2Pos) = 1.02*params(inputkbSGCASsg2Pos);
%              params_plus(inputkfGBSG2Pos) = 1.02*params(inputkfGBSG2Pos);
%              params_plus(inputkbGBSG2Pos) = 1.02*params(inputkbGBSG2Pos);
%              params_plus(inputVMAXPos) = 1.02*params(inputVMAXPos);
%              params_plus(inputCasConcaseFPos) = 1.02*params(inputCasConcaseFPos);
%              params_plus(inputkfGCSG1Pos) = 1.02*params(inputkfGCSG1Pos);
%              params_plus(inputkbGCSG1Pos) = 1.02*params(inputkbGCSG1Pos);
%              params_plus(inputkfSGCASsg1Pos) = 1.02*params(inputkfSGCASsg1Pos);
              params_plus(inputkbSGCASsg1Pos) = 1.02*params(inputkbSGCASsg1Pos);


%            params_min(inputkfGCSG3Pos) = 0.98*params(inputkfGCSG3Pos);
%             params_min(inputkbGCSG3Pos) = 0.98*params(inputkbGCSG3Pos);
%             params_min(inputkMCPICPos) = 0.98*params(inputkMCPICPos);
%             params_min(inputkfSGCASsg3Pos) = 0.98*params(inputkfSGCASsg3Pos);
%             params_min(inputkbSGCASsg3Pos) = 0.98*params(inputkbSGCASsg3Pos);
%             params_min(inputAranPos) = 0.98*params(inputAranPos);

%              params_min(inputkMCPos) = 0.98*params(inputkMCPos);
%              params_min(inputkMB3Pos) = 0.98*params(inputkMB3Pos);
%              params_min(inputdmRNAcaseFPos) = 0.98*params(inputdmRNAcaseFPos);
%              params_min(inputKmPos) = 0.98*params(inputKmPos);
%              params_min(inputkfSGCASsg2Pos) = 0.98*params(inputkfSGCASsg2Pos);
%              params_min(inputkbSGCASsg2Pos) = 0.98*params(inputkbSGCASsg2Pos);
%              params_min(inputkfGBSG2Pos) = 0.98*params(inputkfGBSG2Pos);
%              params_min(inputkbGBSG2Pos) = 0.98*params(inputkbGBSG2Pos);
%              params_min(inputVMAXPos) = 0.98*params(inputVMAXPos);
%              params_min(inputCasConcaseFPos) = 0.98*params(inputCasConcaseFPos);
%              params_min(inputkfGCSG1Pos) = 0.98*params(inputkfGCSG1Pos);
%              params_min(inputkbGCSG1Pos) = 0.98*params(inputkbGCSG1Pos);
%              params_min(inputkfSGCASsg1Pos) = 0.98*params(inputkfSGCASsg1Pos);
              params_min(inputkbSGCASsg1Pos) = 0.98*params(inputkbSGCASsg1Pos);


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
      

%% 10% Visualization
%%%For 10% relative perturbation of the top 6+ params that are the most sensitive
%%%locally from  Ara 0.075 to Ara 6.088 (where the green stripe curves most)
            params_plus=params;
            params_min=params;
            
%            params_plus(inputkfGCSG3Pos) = 1.1*params(inputkfGCSG3Pos);
%            params_plus(inputkbGCSG3Pos) = 1.1*params(inputkbGCSG3Pos);
%             params_plus(inputkMCPICPos) = 1.1*params(inputkMCPICPos);
%             params_plus(inputkfSGCASsg3Pos) = 1.1*params(inputkfSGCASsg3Pos);
%             params_plus(inputkbSGCASsg3Pos) = 1.1*params(inputkbSGCASsg3Pos);
%             params_plus(inputAranPos) = 1.1*params(inputAranPos);

%              params_plus(inputkMCPos) = 1.1*params(inputkMCPos);
%              params_plus(inputkMB3Pos) = 1.1*params(inputkMB3Pos);
%              params_plus(inputdmRNAcaseFPos) = 1.1*params(inputdmRNAcaseFPos);
%              params_plus(inputKmPos) = 1.1*params(inputKmPos);
%              params_plus(inputkfSGCASsg2Pos) = 1.1*params(inputkfSGCASsg2Pos);
%              params_plus(inputkbSGCASsg2Pos) = 1.1*params(inputkbSGCASsg2Pos);
%              params_plus(inputkfGBSG2Pos) = 1.1*params(inputkfGBSG2Pos);
%              params_plus(inputkbGBSG2Pos) = 1.1*params(inputkbGBSG2Pos);
%              params_plus(inputVMAXPos) = 1.1*params(inputVMAXPos);
%              params_plus(inputCasConcaseFPos) = 1.1*params(inputCasConcaseFPos);
%              params_plus(inputkfGCSG1Pos) = 1.1*params(inputkfGCSG1Pos);
%              params_plus(inputkbGCSG1Pos) = 1.1*params(inputkbGCSG1Pos);
%              params_plus(inputkfSGCASsg1Pos) = 1.1*params(inputkfSGCASsg1Pos);
              params_plus(inputkbSGCASsg1Pos) = 1.1*params(inputkbSGCASsg1Pos);
            

%            params_min(inputkfGCSG3Pos) = 0.9*params(inputkfGCSG3Pos);
%             params_min(inputkbGCSG3Pos) = 0.9*params(inputkbGCSG3Pos);
%             params_min(inputkMCPICPos) = 0.9*params(inputkMCPICPos);
%             params_min(inputkfSGCASsg3Pos) = 0.9*params(inputkfSGCASsg3Pos);
%             params_min(inputkbSGCASsg3Pos) = 0.9*params(inputkbSGCASsg3Pos);
%             params_min(inputAranPos) = 0.9*params(inputAranPos);

%              params_min(inputkMCPos) = 0.9*params(inputkMCPos);
%              params_min(inputkMB3Pos) = 0.9*params(inputkMB3Pos);
%              params_min(inputdmRNAcaseFPos) = 0.9*params(inputdmRNAcaseFPos);
%              params_min(inputKmPos) = 0.9*params(inputKmPos);
%              params_min(inputkfSGCASsg2Pos) = 0.9*params(inputkfSGCASsg2Pos);
%              params_min(inputkbSGCASsg2Pos) = 0.9*params(inputkbSGCASsg2Pos);
%              params_min(inputkfGBSG2Pos) = 0.9*params(inputkfGBSG2Pos);
%              params_min(inputkbGBSG2Pos) = 0.9*params(inputkbGBSG2Pos);
%              params_min(inputVMAXPos) = 0.9*params(inputVMAXPos);
%              params_min(inputCasConcaseFPos) = 0.9*params(inputCasConcaseFPos);
%              params_min(inputkfGCSG1Pos) = 0.9*params(inputkfGCSG1Pos);
%              params_min(inputkbGCSG1Pos) = 0.9*params(inputkbGCSG1Pos);
%              params_min(inputkfSGCASsg1Pos) = 0.9*params(inputkfSGCASsg1Pos);
              params_min(inputkbSGCASsg1Pos) = 0.9*params(inputkbSGCASsg1Pos);


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
      