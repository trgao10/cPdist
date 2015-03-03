clearvars;
close all;

%%%% classification level
ClassLevel = 'Genus';

%%%% number of neighbors to vote
VoteNN = 1;

%%%% prepare paths
result_path = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/';
data_path = '../DATA/PNAS/';
spreadsheet_path = [data_path 'GroundTruth.xlsx'];

%%%% load taxa_code and species classification
% taxa_file = [data_path 'taxa_codes_teeth.mat'];
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
[~,~,ClTable] = xlsread(spreadsheet_path);
speciesGroupCol = find(strcmpi(ClTable(1,:), ClassLevel));

%%%% different improvements
MethodsType = {'cPComposedLAST', 'cPComposedLAST', 'cPComposedLAST',...
    'cPComposedLAST', 'cPLAST', 'cPLAST', 'cPMST', 'cPViterbi',...
    'cPViterbi', 'cPViterbi'};
Methods = {'cPComposedLAST', 'cPComposedLASTbalance', 'cPComposedLASTmean',...
    'cPComposedLASTmedian', 'cPLAST', 'cPLASTbalance', 'cPMST', 'cPViterbi',...
    'cPViterbiAngle0.5', 'cPViterbiAngle0.25'};
FeatureFix = {'FeatureFixOff', 'FeatureFixOn'};

for MethodsIdx = 1:length(Methods)
    for FeatureFixIdx = 1:length(FeatureFix)
        %%%% load distance matrix
        %%% old
        % DistMatrixPath = '../DATA/PNAS/teethcP.mat';
        % D = load(DistMatrixPath);
        % D = D.Dist;
        %%% new
        % DistMatrixPath = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/cPDist/cPDistMatrix.mat';
        % D = load(DistMatrixPath, 'cPDistMatrix');
        % D = D.cPDistMatrix;
        %%% improved
        DistMatrixPath = [result_path Methods{MethodsIdx} filesep FeatureFix{FeatureFixIdx} filesep];
        D = load([DistMatrixPath MethodsType{MethodsIdx} 'DistMatrix'], 'ImprDistMatrix');
        D = D.ImprDistMatrix;
        
        %%%
        D = D + diag(Inf(GroupSize,1));
        
        %%%% leave-one-out test
        successCount = 0;
        effectiveCount = 0;
        for j=1:GroupSize
            indexInClTable = find(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),'Names')),taxa_code{j}));
            if ~isnumeric(ClTable{indexInClTable, speciesGroupCol})
                continue;
            else
                effectiveCount = effectiveCount+1;
            end
            [~, nnIndices] = sort(D(j,:));
            delInd = [];
            for s=1:GroupSize
                tmpName = taxa_code{nnIndices(s)};
                tmpIndexInClTable = find(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),'Names')),tmpName));
                if ~isnumeric(ClTable{tmpIndexInClTable, speciesGroupCol})
                    delInd = [delInd s];
                end
            end
            nnIndices(delInd) = [];
            Votes = zeros(1,VoteNN);
            for ii=1:VoteNN
                nnIndexInClTable = find(strcmpi(ClTable(1:end,strcmpi(ClTable(1,:),'Names')),taxa_code{nnIndices(ii)}));
                Votes(ii) = ClTable{nnIndexInClTable, speciesGroupCol};
            end
            if (ClTable{indexInClTable, speciesGroupCol} == mode(Votes))
                successCount = successCount+1;
            else
%                 keyboard
            end
        end
        disp('--------------------------------------------------------------');
        disp(['Method = ' Methods{MethodsIdx} ', ' FeatureFix{FeatureFixIdx}]);
        disp(['Classification Level: ' ClassLevel]);
        disp(['successCount = ' num2str(successCount) '/' num2str(effectiveCount)]);
        disp(['Success rate = ' num2str(successCount/effectiveCount)]);
        disp('--------------------------------------------------------------');
    end
end
