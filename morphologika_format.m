%% Preparation
clear vars;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% Set Parameters
ImprType = 'MST';
FeatureFix = 'On';
NumFeatPts = 256;
output_filename = [pwd '/morphologika/cP' ImprType '_morphologika.txt'];
delete_command = 'rm -f ';

%% Set Path
ResultPath = ['/media/trgao10/Work/MATLAB/ArchivedResults/cP' ImprType '/'];
DistPath = [ResultPath 'FeatureFix' uplow(FeatureFix) '/cP' ImprType 'DistMatrix'];
MapsPath = [ResultPath 'FeatureFix' uplow(FeatureFix) '/cP' ImprType 'MapsMatrix'];
SamplePath = [pwd '/samples/Teeth/'];
DataPath = '~/Work/MATLAB/DATA/PNAS/';
TaxaPath = [DataPath 'teeth_taxa_table.mat'];
cPDistPath = './results/Teeth/cPdist/cPDistMatrix.mat';

%% load results
load(TaxaPath);
load(cPDistPath);
disp('loading maps...');
load(MapsPath);
disp('loaded.');

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




