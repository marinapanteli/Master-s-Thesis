function out = myCost(params, varargin)
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
         arabinose = [13321.788;4442.816;1478.718;493.572;164.524;54.82;18.25;6.088;2.031;0.679;0.226;0.075;0.025;0.008;0.003;0];
%           for i=1:length(viablePoints.projected.names)
%               proj = find(strcmp([DATA.paramNames], viablePoints.projected.names(i)));
%               params(proj) = 0;
%           end
            %params = params';
            %T1%params = [0.0346736850000000; params(1:3); 1; params(4:10); 0; 0; params(11:20); 1.479108388000000;1.479108388000000; params(21:26); 6.309573445000000e+03; params(27:33); 1.479108388000000; 0; 0; 0; params(34:35); 1.479108388000000; params(36:37); 1.479108388000000];
            %T2%params = [0.0346736850000000; params(1:3); 1; params(4:10);0;0; params(11:12); 0;0;params(13:18);1.47910838800000;1.47910838800000;params(19:24);6309.57344500000;params(25:31);1.47910838800000;0;0;0;params(32:33);1.47910838800000;0;0;0];    
            %T3
            %params = [0.0346736850000000; params(1:3); 1; params(4:10); 0; 0;0;0;0;0; params(11:16); 1.479108388000000;1.479108388000000; params(17:22); 6.309573445000000e+03; params(23:29); 1.479108388000000; 0; 0; 0; 0; 0; 0; 0; 0; 0];
            %T4%params = [0.0346736850000000; params(1:3); 1; params(4:22); 1.479108388000000;1.479108388000000; params(23:28); 6.309573445000000e+03; params(29:35); 1.479108388000000; params(36:37); 1.47910838800000; params(38:39); 1.47910838800000; params(40:41); 1.47910838800000];
            [~,inputAraPos] = ismember('ARA',DATA.paramNames);
            params(inputAraPos) = 0;
            ICs = IQMPsimulate(DATA.mex_model, [0:50:10000], [], DATA.paramNames, params);
            [~,FluoPos] = ismember({'Pa','Pb','Pc'},ICs.states);     
             for a = 1:16
                     params(inputAraPos) = DATA.exp_data{a,1};
                     simulation = IQMPsimulate(DATA.mex_model,[0:10:900], ICs.statevalues(end,:), DATA.paramNames, params);
                     result900(a,:) =  simulation.statevalues(end,FluoPos);
             end
             stdval = [zeros(5,1)+0.2;zeros(6,1)+0.1;zeros(5,1)+0.2];
             out = sumsqr((result900./max(result900) - DATA.exp_data{:,2:4})./stdval);
%              figure()
%         plot(result900(:,1)./max(result900(:,1)),'--y')
%         hold on
%         plot(result900(:,2)./max(result900(:,2)),'--b')
%         plot(result900(:,3)./max(result900(:,3)),'--g')
%         errorbar([1:16],DATA.exp_data{:,2},stdval,'-y');
%         errorbar([1:16],DATA.exp_data{:,3},stdval,'-b');
%         errorbar([1:16],DATA.exp_data{:,4},stdval,'-g');
%              %out = sumsqr(result900./max(result900) - DATA.exp_data{:,2:4});

end
end
