function [media, addedComponents] = findMinimalMedia(modelCom, varargin)
% Find one or more minimal media for a single-organism model or a community model
%
% USAGE:
%    [media, addedComponents] = findMinMedia(modelCom, K, baseMedia, candidates, members)
%
% INPUT:
%    modelCom          single-organism model or community model with additional fields 
%                      *.indCom and *.infoCom constructed using `createCommModel`
%
% OPTIONAL INPUTS:
%    K                 identify up to K different media (default 1)
%    baseMedia         cell array of exchange reactions representing the media components that must be present
%                      (default the exchange reactions with *.lb < 0 in the input model)
%    candidates        cell array of exchange reactions representing
%                      candidate metabolites that can be added to the media to support growth
%                      (default all exchange reactions)
%    members           cell array of members required to grow (must be subset of infoCom.spAbbr) 
%                      or index vector (subset of 1:numel(infoCom.spAbbr)
%                      (only appliable for community models, default all members)
%
% OUTPUTS:
%    media             cell array of exchange reactions forming the minimal media
%                      If K > 1, a cell array of K cell arrays is returned, each being a medium
%    addedComponents   cell array of exchange reactions in media but not in baseMedia
%                      If K > 1, a cell array of K cell arrays is returned.

%% argument handling modified from `functionTemplate`
optArgDefaultValidator = {...
    'K',             1,                        @(x) isscalar(x) && floor(x) == x
    'baseMedia',     [],                       @(x) ischar(x) | iscellstr(x); ...
    'candidates',    [],                       @(x) ischar(x) | iscellstr(x); ...
    'members',       [],                       @(x) ischar(x) | iscellstr(x) | isnumeric(x); ...
    };


% get all potentially supplied COBRA parameter names
problemTypes = {'MILP'};

[funParams, cobraParams, solverVarargin] = parseCobraVarargin(varargin, ...
    optArgDefaultValidator(:, 1), optArgDefaultValidator(:, 2), optArgDefaultValidator(:, 3), ...
    problemTypes);

[K, baseMedia, candidates, members] = deal(funParams{:});
%%

bigM = 1000;
[m, n] = size(modelCom.S);

iscom = true;
if ~isfield(modelCom, 'indCom') 
    if ~isfield(modelCom, 'infoCom')
        iscom = false;
        %single-organism model
        modelCom.indCom.EXcom = find(sum(modelCom.S ~= 0, 1) == 1)';
        %         neg = any(modelCom.S(:,modelCom.indCom.EXcom) > 0, 1);
        %         %convert to have positive stoichiometry
        %         for j = 1:size(modelCom.indCom.EXcom,1)
        %             if any(modelCom.S(:,modelCom.indCom.EXcom(j,1)) < 0)
        %                 modelCom.S(:,modelCom.indCom.EXcom(j,1)) = abs(modelCom.S(:,modelCom.indCom.EXcom(j,1)));
        %                 [modelCom.lb(modelCom.indCom.EXcom(neg)), modelCom.ub(modelCom.indCom.EXcom(neg))] = ...
        %                     deal(-modelCom.ub(modelCom.indCom.EXcom(neg)), -modelCom.lb(modelCom.indCom.EXcom(neg)));
        %             end
        %         end
    else
        modelCom.indCom = infoCom2indCom(modelCom, modelCom.infoCom);
    end
end
if iscom
    if isempty(members)
        members = 1:numel(modelCom.infoCom.spAbbr);
    elseif ~isnumeric(members)
        members = find(strcmp(modelCom.infoCom.spAbbr, members));
    end
    if size(modelCom.indCom.EXcom, 2) == 2
        % revert the direction of the uptake reactions for convenience
        modelCom.S(:, modelCom.indCom.EXcom) = -modelCom.S(:, modelCom.indCom.EXcom);
        [modelCom.lb(modelCom.indCom.EXcom),  modelCom.ub(modelCom.indCom.EXcom)] = ...
            deal(-modelCom.ub(modelCom.indCom.EXcom), -modelCom.lb(modelCom.indCom.EXcom));
    end
end

nCom = size(modelCom.indCom.EXcom, 1);

if isempty(baseMedia) && ~iscell(baseMedia)
    % use the default uptake conditions as baseMedia if not given
    on = modelCom.lb(modelCom.indCom.EXcom(:, 1)) < 0;
    baseMedia = modelCom.rxns(modelCom.indCom.EXcom(on, 1));
else
    [yn, on] = ismember(baseMedia, modelCom.rxns(modelCom.indCom.EXcom(:, 1)));
    if ~all(yn)
        error(sprintf('The following identifer(s) in `baseMedia` is(are) not uptake/exchange reactions in the model:\n%s\n', ...
            strjoin(baseMedia(~yn), '\n')))
    end
end
lbMedium = zeros(nCom,1);
lbMedium(on) = 1;

if isempty(candidates)
    ubMedium = ones(nCom,1);
else
    % shut down exchange reactions not in candidates
    ubMedium = lbMedium;
    [yn, on] = ismember(candidates, modelCom.rxns(modelCom.indCom.EXcom(:, 1)));
    if ~all(yn)
        error(sprintf('The following identifer(s) in `candidates` is(are) not uptake/exchange reactions in the model:\n%s\n', ...
            strjoin(candidates(~yn), '\n')))
    end
    ubMedium(on) = 1;
end

lb = modelCom.lb;
lb(modelCom.indCom.EXcom(:, 1)) = -1000;
lb(modelCom.indCom.EXcom(ubMedium == 0, 1)) = 0;

% if exist('Cplex','class')
%     LP = Cplex();
%     %add flux v
%     LP.addCols(zeros(n,1), [], modelCom.lb, ub,[],char(modelCom.rxns));
%     %Sv = 0
%     LP.addRows(zeros(m,1), modelCom.S, zeros(m,1),char(modelCom.mets));
%     %add binary variables a_uptake for medium components
%     LP.addCols(ones(nCom,1), [], lbMedium, ubMedium, char('B'*ones(1,nCom)),char(strcat(modelCom.rxns(modelCom.indCom.EXcom(:,1)),'_onoff')));
%     LP = setCplexParam(LP,solverParam);
%     % v_uptake + bigM * a_uptake >= 0
%     LP.addRows(zeros(nCom,1), [sparse(1:nCom, modelCom.indCom.EXcom(:,1), ones(nCom,1), nCom,n), ...
%         bigM * speye(nCom)], inf(nCom,1), char(strcat(modelCom.rxns(modelCom.indCom.EXcom(:,1)),'_ut')));
%     if iscom
%         LP.Model.lb(modelCom.indCom.spBm(members)) = 0.1;
%     else
%         LP.Model.lb(modelCom.c~=0) = 0.02;
%     end
%     LP.solve;
%     if LP.Solution.status == 101
%         metMMIdCom = find(LP.Solution.x(n+1:end) > 0.1);
%         [metCom,~] = find(modelCom.S(:,modelCom.indCom.EXcom(:,1)));
%         metCom = modelCom.mets(metCom);
%         %     metCom = modelCom.mets(modelCom.metSps == 0);
%         metMM = metCom(metMMIdCom);
%     else
%         metMMIdCom = [];
%         metMM = {};
%     end
% else
    LP = struct();
    
    LP.A = [modelCom.S, sparse(m, nCom); ...  % Sv = 0
        sparse(1:nCom, modelCom.indCom.EXcom(:, 1), 1, nCom, n), bigM * speye(nCom)];   % v_uptake + bigM * a_uptake >= 0
    LP.b = [modelCom.b; zeros(nCom, 1)];
    LP.c = [zeros(n, 1); ones(nCom, 1)];
    LP.lb = [lb; lbMedium];
    LP.ub = [modelCom.ub; ubMedium];
    LP.osense = 1;
    LP.csense = char(['E' * ones(m, 1); 'G' * ones(nCom, 1)]);
    LP.vartype = char(['C' * ones(n, 1); 'B' * ones(nCom, 1)]);
    if iscom
        LP.lb(modelCom.indCom.spBm(members)) = 0.1;
    else
        LP.lb(modelCom.c ~= 0) = 0.02;
    end
    k = 0;
    media = {};
    addedComponents = {};
    stop = false;
    while true
        sol = solveCobraMILP(LP, solverVarargin.MILP{:});
    
        if checkSolFeas(LP, sol) <= cobraParams.MILP.feasTol ...
                && all(abs(sol.int) <= 1e-6 | abs(sol.int - 1) <= 1e-6)
            
            media = [media, {modelCom.rxns(modelCom.indCom.EXcom(sol.full((n + 1):end) > 0.1, 1))}];
            addedComponents = [addedComponents, {setdiff(media{end}, baseMedia, 'stable')}];
            k = k + 1;
            if k == K
                stop = true;
            end
        else
            stop = true;
        end
        if stop
            if k == 0
                warning('A medium that supports growth cannot be identified.\n')
            elseif k == 1
                fprintf('A medium that supports growth has been identified.\n')
                if K == 1
                    media = media{1};
                    addedComponents = addedComponents{1};
                end
            else
                fprintf('%d different media that support growth have been identified.\n', k)
            end
            break
        end
        % add integer cut
        idInt1 = find(sol.int > 0.9);
        % sum(j for a_j = 1 in the last solution, 1- a_j) >= 1
        LP.A = [LP.A; ...
            sparse(1, n), sparse(1, idInt1, ones(numel(idInt1), 1), 1, nCom)];
        LP.b = [LP.b; numel(idInt1) - 1];
        LP.csense = [LP.csense; 'L'];
    end
% end
end