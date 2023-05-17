%% input arguments
% spreadsheet data of taxa abundance and additionally the taxonomy of
% different taxa. The example here the spreadsheet is named 'data'
data = load('allTaxaAbundanceData.mat');
data = data.data;

taxaLv = 'Genus';  % taxonomic level abstraction. Need to match the header name in the spreedsheet
nM = 10;  % number of taxa to include (starting from the highest abundant one in average)
pathToStrainModels = 'modelMat';
% this is for site 1
saveModelName = ['Site1Model_' taxaLv '.mat'];
% columns of abundnace data to include
colAbDataInclude = 8:24;


%% process the data
if pathToStrainModels(end) == filesep
    % no file separator is needed at the end of the path
    pathToStrainModels(end) = '';
end

% header in the spreadsheet
header = data(1, :);
data(1, :) = [];
% find the column corresponding to the taxonomic level
colTaxa = find(strcmpi(header, taxaLv));
% get a unique list of taxa
[taxaList, ~, id] = unique(data(:, colTaxa));

% sum the abundance of all the taxa falling under that taxonomy
taxaAbOverTime = zeros(numel(taxaList), numel(colAbDataInclude));
for j = 1:numel(taxaList)
    taxaAbOverTime(j, :) = sum(cell2mat(data(id == j, colAbDataInclude)), 1);
end

taxaAbMean = mean(taxaAbOverTime, 2);

% sort the taxa by the summed abundance
[~, idSortByAbund] = sort(taxaAbMean, 'descend');
% include only the first nM most abundant taxa
taxaInclude = taxaList(idSortByAbund(1:nM));
% their relative abundance
taxaAbInclude = taxaAbMean(idSortByAbund(1:nM));
taxaAbIncludeOT = taxaAbOverTime(idSortByAbund(1:nM), :);


%% Collate strain models into taxonomic group models
% strains ID under each taxonomic group and their collated model 
strainsInTaxa = cell(numel(taxaInclude), 1);
models = cell(numel(taxaInclude), 1);

for j = 1:numel(taxaInclude)
    % get the strain IDs for each included taxonomic group
    strainsInTaxa{j} = data(id == idSortByAbund(j), strcmp(header, '#OTU ID'));
    
    % load the first existing model. (They might not exist)
    k0 = 1;
    while true
        if exist([pathToStrainModels, filesep, strainsInTaxa{j}{k0}, '.mat'], 'file')
            models{j} = load([pathToStrainModels, filesep, strainsInTaxa{j}{k0}, '.mat']);
            models{j} = models{j}.model;
            break
        end
        k0 = k0 + 1;
    end
    for k = (k0 + 1):numel(strainsInTaxa{j})
        if exist([pathToStrainModels, filesep, strainsInTaxa{j}{k}, '.mat'], 'file')
            % then add new reactions from all other models under the same
            % taxonomic group
            modelK = load([pathToStrainModels, filesep, strainsInTaxa{j}{k}, '.mat']);
            modelK = modelK.model;
            rxnsAdd = setdiff(modelK.rxns, models{j}.rxns);
            rxnsAddId = findRxnIDs(modelK, rxnsAdd);
            for i = 1:numel(rxnsAdd)
                models{j} = addReaction(models{j}, rxnsAdd{i}, ...
                    'reactionName', modelK.rxnNames{rxnsAddId(i)},  ...
                    'metaboliteList', modelK.mets(modelK.S(:, rxnsAddId(i)) ~= 0), ...
                    'stoichCoeffList', modelK.S(modelK.S(:, rxnsAddId(i)) ~= 0, rxnsAddId(i)), ...
                    'lowerBound', modelK.lb(rxnsAddId(i)), 'upperBound', modelK.ub(rxnsAddId(i)));
            end
        end
        if mod(k, 10) == 0
            fprintf('taxa #%d:  %d / %d   %04d-%02d-%02d %02d:%02d:%02.0f\n', j, k, numel(strainsInTaxa{j}), clock)
        end
    end
    
end

%% Try to reduce the model. Not really necessary. Seems only small reduction

modelR = cell(size(models));
[size0, sizeR] = deal(zeros(numel(models), 2));
for j = 1:numel(models)
    [rxnActive, N, loopInfo] = findMinNullGeneral(models{j},0,1,'printLevel',2);
    rxnInactive = ~any(rxnActive, 2);
    modelR{j} = removeRxns(models{j}, models{j}.rxns(rxnInactive));
    size0(j, :) = size(models{j}.S);
    sizeR(j, :) = size(modelR{j}.S);
end

%% Create community model
% print overlapping reactions between models of different taxonomic groups
for j = 1:(numel(models)-1)
    for k = (j+1):numel(models)
        fprintf('#overlaps    #rxns in %-30s  #rxns in %-30s:\n', taxaInclude{j}, taxaInclude{k})
        fprintf('%12d %40d %40d\n', numel(intersect(modelR{j}.rxns,modelR{k}.rxns)), ...
            size(modelR{j}.S,2),size(modelR{k}.S,2));
    end
end
% build the model
options = struct();
options.spBm = repmat({'Growth'}, numel(modelR), 1);
options.spAbbr = regexp(taxaInclude, '_([^_]{2})', 'tokens', 'once');
options.spAbbr = cellfun(@(x) x{1}, options.spAbbr, 'UniformOutput', false);
options.spName = taxaInclude;
options.spATPM = repmat({'R_ATPM'}, numel(modelR), 1);
options.sepUtEx = false;
options.metExId = '[C_e]';
model = createCommModel(modelR, options);

save(saveModelName, 'models', 'taxaInclude', 'taxaAbInclude', 'taxaAbIncludeOT', 'strainsInTaxa', ...
    'modelR', 'size0', 'sizeR', 'model')
