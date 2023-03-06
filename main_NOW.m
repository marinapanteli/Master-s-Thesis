%% Preliminary
clear all;
close all;
clc;

TextSizes.DefaultAxesFontSize = 14;
TextSizes.DefaultTextFontSize = 18;
set(0,TextSizes);

% Seed random number generator
rng(0);

%% User defined files and parameters 

% !! Change following to load from your own directory                   
load('U:\THESIS\Thesis\topodesign-2.0-master\examples\Marina\Save\Saved_viablePoints_sampled_2151768_00000000000000001100000000000000000000000000000111.mat', 'viablePoints')
N_cost = 10; % viable points of model12345 that gives e.g. the tenth smallest cost-->It's just some random set of params that are within the min and max bounds for each param
sigma = 0.1;
t_to_ICs = (0:50:100000);
Ara_value = 6.088;
t_to_sim = (0:10:1500);

PL_search_param_values = 15;
perc_left = 0.8;
perc_right = 1.2;

i=4; % Module index

param_from_meigo = sprintf('paraOpt_fromMEIGO_module%dDOK.mat',i);
%param_from_meigo = sprintf('Copy_paraOpt_fromMEIGO_module%dDOK.mat',i);
module_name = sprintf('Module_%d.txtbc',i);
%paramSpecs = sprintf('parametersTHESIS_MOD%d.txt',i);  
paramSpecs = sprintf('Copy_of_parametersTHESIS_MOD%d.txt',i);  %                
%%data_path = sprintf('GeneratedData_Module%d.csv',i);                        
data_path_dok = sprintf('GeneratedData_ModuleDOK%d.csv',i);                        

%% Set outputs and defined parameter names
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

parnames=IQMparameters(mod);
%% Data
% We generate a data set. It consists of a vector of time points t and an 
% observed output, array yobs. 

% ICs
idx_par = zeros(1,length(parnames));
for i = 1:length(parnames)
    exact_match_mask = strcmp(viablePoints.colnames, parnames(i));
    exact_match_locations = find(exact_match_mask);
    idx_par(i) = exact_match_locations;
end

[m,I]=mink(viablePoints.cost,N_cost); 
params = viablePoints.rowmat(I(end),idx_par); 
params=params';
[~,inputAraPos] = ismember('ARA',parnames);
params(inputAraPos) = 0;
ICs = IQMPsimulate(mod, t_to_ICs, [], parnames, params);

% Forward Simulation
[~,inputAraPos] = ismember('ARA',parnames);
params(inputAraPos) = Ara_value;  

simulation_res = IQMPsimulate(mod,t_to_sim, ICs.statevalues(end,:), parnames, params);

% Add noise (w/ rdm seed)
t =  simulation_res.time';
ym = simulation_res.statevalues(:,states_idx); 
mu = ym;
sigmas = sigma*ones(length(states_idx),length(mu))';

%rng('default');  % For reproducibility
yobs = normrnd(mu,sigmas, length(mu),length(states_idx));  

yobs(yobs<0) = 0;

% T = array2table([t_to_sim' yobs]);
% % T.Properties.VariableNames = {string(simulation_res.states(states_idx))};
% writetable(T,data_path_dok)
%% Parameter Estimation (MLE) was done with MEIGO, using the ESS algorithm
% The resulting parameters are generated using to_run_meigo.m and saved in paraOpt_fromMEIGO_moduleX.mat

load(param_from_meigo, 'paraOpt')       


%% Visualization of fit
% The 'measured' data is visualized in a plot, together with fit for the best
% parameter value found with MEIGO.

% Getting names
table_paramSpecs = readtable(paramSpecs);
parnames = table_paramSpecs.names;

% Forward simulation with generated parameters
params = paraOpt';
[~,inputAraPos] = ismember('ARA',parnames);
params(inputAraPos) = 0;
ICs = IQMPsimulate(mod, t_to_ICs, [], parnames, params);

params(inputAraPos) = Ara_value;                          

ysim = IQMPsimulate(mod,t_to_sim, ICs.statevalues(end,:), parnames, params);

% Plot: Fit  
for i = 1:length(states_idx)
    figure();
    plot(t,yobs(:,i),'bo'); hold on;                                
    plot(t_to_sim,ysim.statevalues(:,states_idx(i)),'r-');                       
    xlabel('time (min)');
    sth = simulation_res.states(states_idx(i));
    ylabel(sprintf('%s (nM)',sth{1}));
    %legend('data', 'fit');
end


%% Profile likelihood calculation -- Parameters
% I found https://github.com/nvanriel/profile-likelihood which calculates
% the PL for each parameter. I will be using MEIGO instead, but nvanriel's
% work is helpful as a guide as well, just to see other people's approach.

costFunction = @myCost_modules;                                    
paramSpecs_tabled = readtable(paramSpecs);                          
model_path = module_name;                                 

% Initialise cost function
DATA = myCost_modules(0, data_path_dok, model_path, paramSpecs, t_to_ICs, Ara_value, states_idx, t_to_sim, sigma);
IQMmakeMEXmodel(DATA.model, DATA.mex_model);

[~,AraPos] = ismember('ARA',paramSpecs_tabled.names);
paraOpt(AraPos) = Ara_value;                                        

paramSpecs_tabled.bmin(AraPos) = Ara_value;
paramSpecs_tabled.bmax(AraPos) = Ara_value;
paramSpecs_tabled.p0(AraPos) = Ara_value;

PL_mat = zeros(PL_search_param_values,2+length(paraOpt),length(paraOpt));




for i=1:length(paraOpt)

     % Module 1
%      if i == 3
%          par_vals = linspace(0.00000001, 0.0004, PL_search_param_values);
%      elseif i==12
%          par_vals = linspace(0.00001, 0.001, PL_search_param_values);
%      else
%         par_vals = linspace(perc_left*paraOpt(i), perc_right*paraOpt(i), PL_search_param_values);
%      end

     % Module 2
%       if i == 3
%          par_vals = linspace(0.00000001, 0.0004, PL_search_param_values);
%      elseif i==9
%          par_vals = linspace(0.01,0.079432823, PL_search_param_values);
%      elseif i==10
%          par_vals = linspace(0.01,0.794328235, PL_search_param_values);
%      elseif i==16
%          par_vals = linspace(0.00001,0.001, PL_search_param_values);
%      elseif i==17
%          par_vals = linspace(0.00001,0.001, PL_search_param_values);
%      elseif i==18
%          par_vals = linspace(6000, 7000, PL_search_param_values);
%      elseif i==20
%          par_vals = linspace(0.00001, 0.004, PL_search_param_values);
%      elseif i==22
%          par_vals = linspace(0.0001, 0.03, PL_search_param_values);
%      elseif i==23
%          par_vals = linspace(0.00001, 0.04, PL_search_param_values);
%       else
%         par_vals = linspace(perc_left*paraOpt(i), perc_right*paraOpt(i), PL_search_param_values);
%      end


    % Module 3
%      if i == 3
%          par_vals = linspace(0.00000001, 0.0004, PL_search_param_values);
%      elseif i==7
%          par_vals = linspace(0.001,0.079432823, PL_search_param_values);
%      elseif i==10
%          par_vals = linspace(0.0001,0.794328235, PL_search_param_values);
%      elseif i==8
%          par_vals = linspace(0.01,0.794328235, PL_search_param_values);
%      elseif i==13
%          par_vals = linspace(0.0001,0.1, PL_search_param_values);
%      elseif i==14
%          par_vals = linspace(0.01,0.794328235, PL_search_param_values);
%      elseif i==15
%          par_vals = linspace(0.79, 1.18, PL_search_param_values);
%       else
%         par_vals = linspace(perc_left*paraOpt(i), perc_right*paraOpt(i), PL_search_param_values);
%      end

    % Module 4
%      if i == 3
%          par_vals = linspace(0.00000001, 0.01, PL_search_param_values);
%      elseif i==13
%          par_vals = linspace(0.0001,10, PL_search_param_values);
%      elseif i==16
%          par_vals = linspace(0.00001,0.001, PL_search_param_values);
%      elseif i==17
%          par_vals = linspace(0.00001,0.001, PL_search_param_values);
%       else
%         par_vals = linspace(perc_left*paraOpt(i), perc_right*paraOpt(i), PL_search_param_values);
%      end

    % Module 5
     if i == 8
         par_vals = linspace(0.2, 80, PL_search_param_values);
     elseif i==13
         par_vals = linspace(0.01,200, PL_search_param_values);
     elseif i==10
         par_vals = linspace(1,2, PL_search_param_values);
      else
        par_vals = linspace(perc_left*paraOpt(i), perc_right*paraOpt(i), PL_search_param_values);
     end

        for j=1:PL_search_param_values
       
            table_paramSpecs_altered = paramSpecs_tabled;
        
            table_paramSpecs_altered.bmin(i) = par_vals(j);
            table_paramSpecs_altered.bmax(i) = par_vals(j);
            table_paramSpecs_altered.p0(i) = par_vals(j);
            [OUT_ret,paraOpt_PL] = meigo_modules(costFunction, module_name, paramSpecs, table_paramSpecs_altered, data_path_dok, t_to_ICs, Ara_value, t_to_sim, states_idx,sigma);  % POINT-->loglikelihood
        
            temp = zeros(1,2+length(paraOpt));
            temp(1) = par_vals(j);
            temp(2) = OUT_ret;
            temp(3:3+length(paraOpt)-1) = paraOpt_PL;
            
            PL_mat(j,:,i) = temp;

        end

end

alpha = 0.95;

% Confidence Intervals: To add to optimal Meigo sum of least squares (w/ the sigma^2):
Da = icdf('chi2',alpha,1);
CI_line = myCost_modules(paraOpt) + Da;

opt_point_cost = myCost_modules(paraOpt);

% Profile Likelihood Figures

counter = 0; % Set so as not to plot the Ara which is kept fixed

for i = 1:size(PL_mat,3)

    if i~=AraPos
        counter = counter + 1;

        opt_point_param_val = paraOpt(i);

        subplot(ceil(size(PL_mat,3)/4),4,counter)
        plot(PL_mat(:,1,i),PL_mat(:,2,i),'*'); hold on;
        plot(opt_point_param_val,opt_point_cost,'o');hold on;
        yline(CI_line,'--r',sprintf('CI=%1.1f',CI_line)); hold on;

        max_limy = max([PL_mat(:,2,i);CI_line;opt_point_cost]);
        min_limy = min([PL_mat(:,2,i);CI_line;opt_point_cost]);
        range_limy = max_limy-min_limy;

        ylim([min_limy-0.5*range_limy max_limy+1.5*range_limy])

        xlabel(paramSpecs_tabled.names(i));
        ylabel('Cost');
    end

end

%% Evaluation of chosen module
% 
% i = 4;
% 
% sigma = 0.1;
% t_to_ICs = (0:50:100000);
% Ara_value = 6.088;
% t_to_sim = (0:10:1500);
% 
% module_name = sprintf('Module_%d.txtbc',i);
% paramSpecs = sprintf('For_eval_of_parametersTHESIS_MOD%d.txt',i);  %                
% data_path_dok = sprintf('GeneratedData_ModuleDOK%d.csv',i);                        
% 
% mod = IQMmodel(module_name);
% Pa_idx = find(strcmp(IQMstates(mod), 'Pa'));
% Pb_idx = find(strcmp(IQMstates(mod), 'Pb'));
% Pc_idx = find(strcmp(IQMstates(mod), 'Pc'));
% 
% states_idx = [];
% 
% if ~isempty(Pa_idx)
%     states_idx = [states_idx Pa_idx];
% end
% if ~isempty(Pb_idx)
%     states_idx = [states_idx Pb_idx];
% end
% if ~isempty(Pc_idx)
%     states_idx = [states_idx Pc_idx];
% end
% 
% parnames=IQMparameters(mod);
% 
% costFunction = @myCost_modules;                                    
% paramSpecs_tabled = readtable(paramSpecs);                          
% model_path = module_name;                                 
% 
% % Initialise cost function
% DATA = myCost_modules(0, data_path_dok, model_path, paramSpecs, t_to_ICs, Ara_value, states_idx, t_to_sim, sigma);
% IQMmakeMEXmodel(DATA.model, DATA.mex_model);
% 
% [~,AraPos] = ismember('ARA',paramSpecs_tabled.names);
% paraOpt(AraPos) = Ara_value;   
% 
% table_paramSpecs = readtable(paramSpecs);
% 
% vthreshold = chi2inv(0.95,length(states_idx)*length(t_to_sim)-size(table_paramSpecs,1)+5); % !! # outputs, non-fixed params in module (change depending on module)
% 
% load('x_init.mat', 'paraOpt') % Param set for mod4 whose cost<vthreshold, i.e., is viable
% x0 = paraOpt;
% 
% bmin = table_paramSpecs.bmin;
% bmax = table_paramSpecs.bmax;
% 
% 
% n = 100000;
% 
% OutM = MCexp(costFunction, vthreshold, x0, bmax', bmin', n);
% M_V = OutM.V;
% 
% % tabcheck = zeros(3028,24);
% % 
% % for i = 1:3028
% %     for j = 1:24
% %         if M_V(i,j)<bmin(j)|M_V(i,j)>bmax(j)
% %             tabcheck(i,j)=1;
% %         end
% %     end
% % end
% % 
% % sum(tabcheck,'all')
% 
% % from here on it won't work
% 
% OutE = ELexp(costFunction, vthreshold, M_V, bmax', bmin', 100000);
% E_V = OutE.V;
% 
% tabcheck = zeros(3028,24);
% 
% for i = 1:3028
%     for j = 1:24
%         if E_V(i,j)<bmin(j)|E_V(i,j)>bmax(j)
%             tabcheck(i,j)=1;
%         end
%     end
% end
% 
% sum(tabcheck,'all')
% 
% 
% 
% 
% V_all = vertcat(M_V,E_V);
% 
% OutV = Volint(costFunction, vthreshold, V_all, bmax', bmin', 10000);
% volume_module = OutV.vol;
% 
% 
% 
% bmax1 = bmax;
% bmin1 = bmin;
% bmax1([1,5,15,18,24]) = [];
% bmin1([1,5,15,18,24]) = [];
% 
% 
% spacevol = prod(bmax1 - bmin1);
% 
% Robustness = volume_module/spacevol;
