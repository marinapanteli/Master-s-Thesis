function  [OUT_ret, paraOpt_PL] = meigo_modules(costFunction, module_name, paramSpecs, table_paramSpecs_altered, data_path, t_to_ICs, Ara_value, t_to_sim, states_idx,sigma)

%% Settings
table_paramSpecs = table_paramSpecs_altered;
% Splot = false;
% 
% %% Initialise cost function
% DATA = myCost_modules(0,data_path, module_name, paramSpecs, t_to_ICs, Ara_value, states_idx,  t_to_sim);
% IQMmakeMEXmodel(DATA.model, DATA.mex_model);

%viabilityThreshold = chi2inv(0.95,length(states_idx)*length(t_to_sim)-size(table_paramSpecs,1)+1); % !! # outputs, non-fixed params in module
viabilityThreshold = chi2inv(0.95,length(states_idx)*length(t_to_sim)-size(table_paramSpecs,1)+2); % !! # outputs, non-fixed params in module


%% Run Meigo

npara = length(table_paramSpecs.p0);
problem = [];
opts = [];
problem.f = costFunction;
problem.f=func2str(problem.f);
problem.x_L = table_paramSpecs.bmin;
problem.x_U = table_paramSpecs.bmax;
problem.vtr = viabilityThreshold;
opts.tolc = 1e-3;
opts.iterprint = 1;
% opts.plot = 2;
opts.ndiverse = 10*npara;
opts.maxeval = 1000;
opts.maxtime = 10000;
% opts.inter_save=1;
opts.local.solver  = 'fminsearch'; 
opts.local.iterprint = 1;
opts.local.n1  = 10;
opts.local.n2  = 10; 
opts.local.finish  = 0;
opts.local.bestx  = 1;
opts.local.tol    = 1;
problem.x_0 =  table_paramSpecs.p0;
Results = MEIGO(problem,opts,'ESS');
paraOpt_PL = [Results.xbest];   
OUT_ret = myCost_modules(paraOpt_PL')  ;                                 % THIS POINT-->loglikelihood
end