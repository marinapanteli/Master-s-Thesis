function out = myCost_finaleval(params, varargin)
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
    DATA.mex_model = 'my_mex_model';
    out = DATA;
else
%         arabinose = [13321.788;4442.816;1478.718;493.572;164.524;54.82;18.25;6.088;2.031;0.679;0.226;0.075;0.025;0.008;0.003;0];
         arabinose = 6.088;

            [~,inputAraPos] = ismember('ARA',DATA.paramNames);
            params(inputAraPos) = 0;
            ICs = IQMPsimulate(DATA.mex_model, [0:50:100000], [], DATA.paramNames, params);
            [~,FluoPos] = ismember({'Pa','Pc'},ICs.states);     
            
             params(inputAraPos) = arabinose;
             
             simulation = IQMPsimulate(DATA.mex_model,[0:10:1500], ICs.statevalues(end,:), DATA.paramNames, params);
             
             resultsim =  simulation.statevalues(:,FluoPos);
             

             out = sumsqr(resultsim./max(resultsim) - DATA.exp_data{:,2:end}./max(DATA.exp_data{:,2:end}))*1/(0.15)^2;


end
end
