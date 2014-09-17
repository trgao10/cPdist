%% Preparation
clear vars;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% Set Parameters
ImprType = 'ComposedLAST';
subImprType = [ImprType 'balance'];
FeatureFix = 'On';
NumFeatPts = 64;
output_filename = [pwd '\morphologika\cP' subImprType '_FeatureFix' uplow(FeatureFix) '_morphologika_' num2str(NumFeatPts) '.txt'];
delete_command = 'del ';

%% Set Path
ResultPath = ['D:/Work/MATLAB/ArchivedResults/Teeth/cP' subImprType '/'];
ImprDistPath = [ResultPath 'FeatureFix' uplow(FeatureFix) '/cP' subImprType 'DistMatrix'];
MapsPath = [ResultPath 'FeatureFix' uplow(FeatureFix) '/cP' ImprType 'MapsMatrix'];
SamplePath = [pwd '/samples/Teeth/'];
DataPath = 'D:/Work/MATLAB/ArchivedData/PNAS/';
TaxaPath = [DataPath 'teeth_taxa_table.mat'];
cPDistPath = 'D:/Work/MATLAB/ArchivedResults/Teeth/cPDist/cPDistMatrix.mat';

%% load results
load(TaxaPath);
load(DistPath);
load(cPDistPath);
disp(['loading maps from ' MapsPath '...']);
load(MapsPath);
disp('loaded.');

%% Parse parameters
TaxaCode = load(TaxaPath);
TaxaCode = TaxaCode.taxa_code;
GroupSize = length(TaxaCode);


if strcmpi(ImprType,'Viterbi')
    ImprDistMatrix = sparse(ImprDistMatrix);
    TrilDistMatrix = tril(ImprDistMatrix,-1);
    [~, PRED] = graphminspantree(TrilDistMatrix,'Method','Kruskal');
    RootNode = find(PRED==0);
elseif strcmpi(ImprType,'ComposedLAST') || strcmpi(ImprType,'Dist')
    [~,RootNode] = min(sum(ImprDistMatrix));
else %%% MST or LAST
    if strcmpi(subImprType,'cPLASTbalance')
        options.alpha = 1+sqrt(2);
    end
    [~,PRED] = ConstructGraph(cPDistMatrix,ImprType,options);
    RootNode = find(PRED==0);
end

command_text = [delete_command output_filename];
system(command_text);
disp(command_text);

%% Write Information to Output File
fid = fopen(output_filename,'wt');
if(fid==-1)
    error('Can''t open the file.');
end

%%% header
fprintf(fid, '[Individuals]\n');
fprintf(fid, '%d\n', GroupSize);
fprintf(fid, '[landmarks]\n');
fprintf(fid, '%d\n', NumFeatPts);
fprintf(fid, '[dimensions]\n');
fprintf(fid, '3\n');
fprintf(fid, '[names]\n');
for j=1:length(TaxaCode)
    fprintf(fid, '%s\n', TaxaCode{j});
end

%%% feature points
RootG = load([SamplePath TaxaCode{RootNode} '.mat']);
RootG = RootG.G;
RootNodeFeatureInds = RootG.Aux.DensityPnts(1:NumFeatPts);
fprintf(fid, '\n[rawpoints]\n');
for j=1:length(TaxaCode)
    fprintf(fid,'\n%s\n\n', ['''' TaxaCode{j}]);
    if j==RootNode
        fprintf(fid, '%d %d %d\n', RootG.V(:,RootNodeFeatureInds));
    else
        PropNodeFeatureInds = ImprMapsMatrix{RootNode,j}(RootNodeFeatureInds);
        load([SamplePath TaxaCode{j} '.mat']);
        fprintf(fid, '%d %d %d\n', G.V(:,PropNodeFeatureInds));
        clear G;
    end
end

fclose(fid);




