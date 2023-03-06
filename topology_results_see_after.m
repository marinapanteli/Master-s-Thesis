%% Import
addpath(genpath(['..' filesep '..' filesep 'source']))

%% Settings
costFunction = @myCost_finaleval;
viabilityThreshold = chi2inv(0.95,2*151-41);
%paramSpecs = 'parametersTHESIS_after_arabinose.txt'; %THIS AND COUPLING TOO
paramSpecs = 'Copy_of_parametersTHESIS_finaleval.txt'; %THIS AND COUPLING TOO
folder_save = '24_pos_in';
name_save = 'Saved'; %%%%
%data_path = 'GeneratedData.csv'; %THIS
data_path = 'GeneratedData_ModuleDOK4.csv';
model_path = 'emptymodel.txtbc'; %THIS
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
 parallelize = 0;
% if n_threads>0
%     parpool(par_profile, n_threads);
%     parallelize = 1;
% end
% Initialise cost function
DATA = myCost_finaleval(0, data_path, model_path, paramSpecs);
IQMmakeMEXmodel(DATA.model, DATA.mex_model);

%%

[~,idxcostmin]=min(viablePoints.cost);
[~,idxparams]=ismember(viablePoints.colnames,table_paramSpecs.names);
table_paramSpecs.p0(idxparams)=[10.^viablePoints.rowmat(idxcostmin,:)];
[r1,~]=ismember(1:50,idxparams);
[~,r2]=find(r1==0);
table_paramSpecs.p0(r2)=zeros(1,length(r2));
params=table_paramSpecs.p0;

inputAraPos=5;

params(inputAraPos)=0;

ICs = IQMPsimulate(DATA.mex_model, [0:50:1000000], [], DATA.paramNames, params);
[~,FluoPos] = ismember({'Pa','Pc'},ICs.states);     
            
 dokim = params;
 dokim(inputAraPos) = 6.088;


 simulation = IQMPsimulate(DATA.mex_model,[0:10:1500], ICs.statevalues(end,:), DATA.paramNames, dokim);

 resulting_traj =  simulation.statevalues(:,FluoPos);

            
out = sumsqr(resulting_traj(:,1:2)./max(resulting_traj) - DATA.exp_data{:,2:3}./max(DATA.exp_data{:,2:3}))*1/(0.15)^2;

%  str = '#EDB120';
% color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;

figure()

plot(resulting_traj(:,1)./max(resulting_traj(:,1)),'--b')
hold on;
plot(resulting_traj(:,2)./max(resulting_traj(:,2)),'--r')
hold on;
plot(DATA.exp_data{:,2}./max(DATA.exp_data{:,2}),'+b')
hold on;
plot(DATA.exp_data{:,3}./max(DATA.exp_data{:,3}),'+r')

xlabel("time (min)")
ylabel("Normalized outputs")

