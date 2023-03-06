%% Settings
clc;clear
% User chosen parameters and files
costFunction = @myCost_modules; 


t_to_ICs = [0:50:100000];                  %
Ara_value = 6.088;                         %
t_to_sim = [0:10:1500];                    %
sigma = 0.1;                              %

i_module=4;                                %

module_name = sprintf('Module_%d.txtbc',i_module);            %                
%paramSpecs = sprintf('parametersTHESIS_MOD%d.txt',i_module);  %                
%paramSpecs = sprintf('Copy_of_parametersTHESIS_MOD%d.txt',i_module);  %                
paramSpecs = sprintf('For_eval_of_parametersTHESIS_MOD%d.txt',i_module);  %                
data_path = sprintf('GeneratedData_ModuleDOK%d.csv',i_module);   %             

mod = IQMmodel(module_name);
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

% for k = 1:size(table_paramSpecs,1)
%     if k~=inputAraPos
%     table_paramSpecs.bmin(k) = 0.8*table_paramSpecs.p0(k);
%     table_paramSpecs.bmax(k) = 1.2*table_paramSpecs.p0(k);
%     end
% end

%% Initialise cost function
DATA = myCost_modules(0,data_path, module_name, paramSpecs, t_to_ICs, Ara_value, states_idx,  t_to_sim,sigma);
IQMmakeMEXmodel(DATA.model, DATA.mex_model);

%viabilityThreshold = chi2inv(0.95,length(states_idx)*length(t_to_sim)-size(table_paramSpecs,1)+1); % !! # outputs, non-fixed params in module
viabilityThreshold = chi2inv(0.95,length(states_idx)*length(t_to_sim)-size(table_paramSpecs,1)+5); % !! # outputs, non-fixed params in module


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
opts.maxeval = 200000;
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
Results = MEIGO(problem,opts,'ESS');
paraOpt = [Results.xbest];   
myCost_modules(paraOpt')                                   % THIS POINT-->loglikelihood