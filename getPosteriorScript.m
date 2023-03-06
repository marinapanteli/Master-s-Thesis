%% Import
addpath(genpath(['..' filesep '..' filesep 'source']))

%% Settings
costFunction = @myCost;
alpha_level = 0.95;
path_paramSpecs = 'paramSpecs.txt';
name_save = 'Posterior/Save_1';
data_path = 'raw_data.csv';
model_path = 'model.txt';
n_start = 20;
n_res = 1;
n_threads = 0;
max_eval = 10000;
par_profile = 'SGE_BSSE';


%% Setup
%data
data_file = readtable(data_path);
ndata = numel(data_file(:,2:end));
% Launch parallel pool if necessary
parallelize = 0;
if n_threads>0
    parpool(par_profile, n_threads);
    parallelize = 1;
end
% Initialise cost function
DATA = myCost(0, data_path, model_path, path_paramSpecs);
IQMmakeMEXmodel(DATA.model, DATA.mex_model);
parallelInitPersistent(costFunction, 0, data_path, model_path, path_paramSpecs);

%% Computation
posterior = getPosterior(costFunction, alpha_level, path_paramSpecs, ndata, name_save);
