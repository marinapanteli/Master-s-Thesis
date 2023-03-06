%% Settings
clc;clear
% User chosen parameters and files
costFunction = @myCost_finaleval; 


t_to_ICs = [0:50:100000];                  %
Ara_value = 6.088;                         %
t_to_sim = [0:10:1500];                    %
sigma = 0.15;                              %

i_module=4;                                %

model_name = sprintf('emptymodel.txtbc');            %                
paramSpecs = sprintf('parametersTHESIS_after_arabinose.txt');  %                
data_path = sprintf('GeneratedData_ModuleDOK%d.csv',i_module);   %             

mod = IQMmodel(model_name);
Pa_idx = find(strcmp(IQMstates(mod), 'Pa'));
Pb_idx = find(strcmp(IQMstates(mod), 'Pb'));
Pc_idx = find(strcmp(IQMstates(mod), 'Pc'));

states_idx = [];

if ~isempty(Pa_idx)
    states_idx = [states_idx Pa_idx];
end
if ~isempty(Pb_idx)
    states_idx = [states_idx Pb_idx];
end
if ~isempty(Pc_idx)
    states_idx = [states_idx Pc_idx];
end

%data_path = 'dokimazw.csv';     % ym (i.e., without noise) just for me to try see if cost func works

table_paramSpecs = readtable(paramSpecs);

[~,inputAraPos] = ismember('ARA',table_paramSpecs.names);


%% Initialise cost function
DATA = myCost_finaleval(0,data_path, model_name, paramSpecs, t_to_ICs, Ara_value, states_idx,  t_to_sim,sigma);
IQMmakeMEXmodel(DATA.model, DATA.mex_model);

%viabilityThreshold = chi2inv(0.95,length(states_idx)*length(t_to_sim)-size(table_paramSpecs,1)+1); % !! # outputs, non-fixed params in module
viabilityThreshold = chi2inv(0.95,2*151-41); % !! # outputs, non-fixed params in module


%% Run Meigo

npara = length(table_paramSpecs.p0);
problem = [];
opts = [];
problem.f = costFunction;
problem.f=func2str(problem.f);
problem.x_L = table_paramSpecs.bmin;
problem.x_U = table_paramSpecs.bmax;
problem.vtr = viabilityThreshold;
% opts.tolc = 1e-3;
% opts.iterprint = 1;
% opts.plot = 2;
opts.ndiverse = 10*npara;
opts.maxeval = 10000000;
opts.maxtime = 10000;
% opts.inter_save=1;
% opts.local.solver  = 'fminsearch'; 
% opts.local.iterprint = 1;
% opts.local.n1  = 10;
% opts.local.n2  = 10; 
% opts.local.finish  = 0;
% opts.local.bestx  = 1;
% opts.local.tol    = 1;
problem.x_0 =  table_paramSpecs.p0;
%  load('opt_fineval.mat', 'paraOpt')
%  problem.x_0 = paraOpt;
Results = MEIGO(problem,opts,'ESS');
paraOpt = [Results.xbest];   
myCost_finaleval(paraOpt')                                   % THIS POINT-->loglikelihood