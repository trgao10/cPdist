%% Preparation
clear vars;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% Set Parameters
ImprType = 'MST';
FeatureFix = 'Off';
NumFeatPts = 32;
output_path = [pwd '/morphologika/'];
delete_command = 'rm -f ';
outputDistFile = [output_path 'cP' ImprType 'Dist_FeatFix' uplow(FeatureFix) '_TreeEdges.xlsx'];
outputlmkMSEFile = [output_path 'cP' ImprType 'lmkMSE_FeatFix' uplow(FeatureFix) '_TreeEdges.xlsx'];

%% clear up duplicate outputs
command_text = [delete_command outputDistFile];system(command_text);disp(command_text);
command_text = [delete_command outputlmkMSEFile];system(command_text);disp(command_text);

%% Set Path
ResultPath = ['/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/cP' ImprType '/'];
DistPath = [ResultPath 'FeatureFix' uplow(FeatureFix) '/cP' ImprType 'DistMatrix'];
MSEsPath = [ResultPath 'FeatureFix' uplow(FeatureFix) '/cP' ImprType 'lmkMSEMatrix'];
SamplePath = [pwd '/samples/Teeth/'];
DataPath = '~/Work/MATLAB/DATA/PNAS/';
TaxaPath = [DataPath 'teeth_taxa_table.mat'];
cPDistPath = './results/Teeth/cPdist/cPDistMatrix.mat';

%% load results
load(TaxaPath);
load(DistPath);
load(MSEsPath);
load(cPDistPath);

%% Parse parameters
TaxaCode = load(TaxaPath);
TaxaCode = TaxaCode.taxa_code;
GroupSize = length(TaxaCode);

if strcmpi(ImprType,'Viterbi')
elseif strcmpi(ImprType,'ComposedLAST')
else %%% MST or LAST
    [~,PRED] = ConstructGraph(cPDistMatrix,ImprType);
    RootNode = find(PRED==0);
end

Dists = zeros(GroupSize);
lmkMSEs = zeros(GroupSize);
for j=1:GroupSize
    if j~=RootNode
        Dists(j,PRED(j)) = ImprDistMatrix(j,PRED(j));
        lmkMSEs(j,PRED(j)) = lmkMSEMatrix(j,PRED(j));
    end
end

%% write results to spreadsheets
xlswrite(outputDistFile,[0,TaxaCode;TaxaCode',num2cell(Dists)]);
xlswrite(outputlmkMSEFile,[0,TaxaCode;TaxaCode',num2cell(lmkMSEs)]);


