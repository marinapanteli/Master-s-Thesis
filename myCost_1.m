function out = myCost_1(params, varargin)

% varargin : empty or, for setup : data_file, model_file
persistent DATA Splot
if ~isempty(varargin)
    %Setup persistent DATA
    data_file = varargin{1};
    model_file = varargin{2};
    paramSpecs = readtable(varargin{3});
    Splot = varargin{4};
    DATA.paramNames = paramSpecs.names;
    DATA.exp_data = readtable(data_file);
    DATA.model = IQMmodel(model_file);
    DATA.mex_model = 'my_mex_model';
    out = DATA;
else
    %Compute cost for given param vector
%     try
    params = [0.0346736850000000; params(1:3); 1; params(4:10); 0; 0; params(11:20); 1.479108388000000;1.479108388000000; params(21:26); 6.309573445000000e+03; params(27:33); 1.479108388000000; 0; 0; 0; params(34:35); 1.479108388000000; params(36:37); 1.479108388000000];
    [~,inputAraPos] = ismember('ARA',DATA.paramNames);
    params(inputAraPos) = 0;
    ICs = IQMPsimulate(DATA.mex_model, [0:50:360], [], DATA.paramNames, params);
    [~,FluoPos] = ismember({'Pa','Pb','Pc'},ICs.states);     
    for a = 1:16
            params(inputAraPos) = DATA.exp_data{a,1};
            simulation = IQMPsimulate(DATA.mex_model,[0:10:900], ICs.statevalues(end,:), DATA.paramNames, params);
            result900(a,:) =  simulation.statevalues(end,FluoPos);
    end   
    out = sumsqr(result900./max(result900) - DATA.exp_data{:,2:4});
% if Splot
%      figure()
%         plot(result900(:,1)./max(result900(:,1)),'--y')
%         hold on
%         plot(result900(:,2)./max(result900(:,2)),'--b')
%         plot(result900(:,3)./max(result900(:,3)),'--g')
%         plot(DATA.exp_data{:,2},'-y');
%         plot(DATA.exp_data{:,3},'-b');
%         plot(DATA.exp_data{:,4},'-g');
%  end

end
end
