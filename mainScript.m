%% Import
addpath(genpath(['..' filesep '..' filesep 'source']))

%% Settings
costFunction = @myCost_finaleval;
viabilityThreshold = chi2inv(0.95,2*151-41);
%paramSpecs = 'parametersTHESIS_after_arabinose.txt'; %THIS AND COUPLING TOO
paramSpecs = 'Copy_of_parametersTHESIS_finaleval.txt'; %THIS AND COUPLING TOO
folder_save = '2_pos_in';
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
parallelInitPersistent(@myCost_finaleval, 0, data_path, model_path, paramSpecs);

% Test cost function
%  pars = readtable('parametersTHESIS_after_arabinose.txt');
%  pars = pars.p0;
%  test_cost = myCost_finaleval(pars);
%  disp(test_cost) 

% Set folders paths
folder_out = ['METRICS' filesep name_save];
topofilter_run_path = [folder_save filesep name_save];
viable_points_folder = folder_save;

%% Run
%paramSpecs_new = 'parametersTHESIS_after_arabinose2_for_comp_metrics.txt';
% TopoFilter and save viable points
%if ~(exist(topofilter_run_path, 'file') || exist([topofilter_run_path '.mat'], 'file'))
    [runs, viabilityThreshold, settings] = TFmain(costFunction, viabilityThreshold, paramSpecs,...
                                                   'outputFilename', topofilter_run_path,...
                                                   'parallelize', parallelize,...
                                                   'initSearchNres', 100,...
                                                   'saveViablePoints', true,...
                                                'nmont',16667*10,...
                                                'nelip',33333*10,...
                                                'nvmin',20000,...
                                                'nvmax',95000*10,...
                                               'paramCouplingFcn',@Coupling);
%else
%    loaded = load(topofilter_run_path);
% %    runs = loaded.runs;
% %end
% 
% roothpath = '/homes/mpanteli/THESIS/Thesis/topodesign-2.0-master/examples/Marina/Save/28_pos/28_pos_in';
% fpattern = fullfile(roothpath, 'Saved_viablePoints_projected_*.mat');
% files = dir(fpattern);
% topo_costs = [];
% 
% 
% for k= 1 : length(files)
% fName = files(k).name;
% fFolder = files(k).folder;
% fullname = fullfile(fFolder,fName);
% load(fullname, 'viablePoints')
% [minim, minim_idx] = min(viablePoints.cost);
% topo_costs_info = {minim minim_idx fullname};
% M(k,:) = topo_costs_info;
% end
% % 

% Compute relevant metrics (viable volume, robustness, feasibility)
nTopologies = size(runs{1}.viableProjectionsCollection,1);
OutV = cell(1, nTopologies);
folder_out='METRICS\Saved_edited_28';
viable_points_folder = '28_pos_in_edited';
parfor (iTopo = 1:nTopologies, n_threads)
    OutV{iTopo} = computeMetrics(iTopo, topofilter_run_path, viable_points_folder, folder_out, tunableParamSettings,...
                             posteriorNames, posteriors);
end

all_metrics = gatherMetrics(folder_out);
@Volint
% Compute high feasibility region for one circuit
iMod = 2;
alpha = 0.5;
feasibilityThreshold = 1 - (1 - all_metrics.feasibility_ideal{iMod})/alpha;
feasRegion = feasibilityRegionWrapper(iMod, topofilter_run_path, viable_points_folder,...
                                        folder_out, tunableParamSettings,...
                                        posteriorNames(1), posteriors(1),...
                                        feasibilityThreshold);






