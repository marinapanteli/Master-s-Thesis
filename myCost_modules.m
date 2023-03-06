function out = myCost_modules(params_m,  varargin)
% Very interesting description of a cost function.
% varargin : empty or, for setup : data_file, model_file
persistent DATA
if ~isempty(varargin)
    
    %Setup persistent DATA
    data_file = varargin{1};
    model_file = varargin{2};
    paramSpecs = readtable(varargin{3});
    DATA.paramNames = paramSpecs.names;
    DATA.exp_data = readtable(data_file);
    DATA.model = IQMmodel(model_file);
    DATA.mex_model = 'my_mex_module';
    DATA.t_to_ICs = varargin{4};
    DATA.Ara_value = varargin{5};
    DATA.states_idx = varargin{6};
    DATA.t_to_sim = varargin{7};
    DATA.sigma = varargin{8};

    out = DATA;

else

    [~,inputAraPos] = ismember('ARA',DATA.paramNames);
    params_m(inputAraPos) = 0;
    ICs_m = IQMPsimulate(DATA.mex_model, DATA.t_to_ICs, [], DATA.paramNames, params_m);
    FluoPos = DATA.states_idx;        
    
    params_m(inputAraPos) = DATA.Ara_value;                     
    simulation = IQMPsimulate(DATA.mex_model,DATA.t_to_sim, ICs_m.statevalues(end,:), DATA.paramNames, params_m);
    result =  simulation.statevalues(:,FluoPos);

    out = sumsqr(result-DATA.exp_data{:,2:end})*1/(DATA.sigma)^2; % since in module same sigma for all data points assumed
   
end
end
