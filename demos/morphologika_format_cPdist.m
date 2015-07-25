%% Preparation
clear vars;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%% Set Parameters
ImprType = 'Dist';
% subImprType = ImprType;
% subImprType = [ImprType 'median'];
% FeatureFix = 'On';
NumFeatPts = 1024;
output_filename = [pwd '/morphologika/cPDist_morphologika_' num2str(NumFeatPts) '.txt'];
delete_command = 'rm -rf ';

%% Set Path
ResultPath = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/cPDist';
MapsPath = [ResultPath '/cPMapsMatrix'];
SamplePath = [pwd '/samples/Teeth/'];
DataPath = '~/Work/MATLAB/DATA/PNAS/';
TaxaPath = [DataPath 'teeth_taxa_table.mat'];
cPDistPath = '/media/trgao10/Work/MATLAB/ArchivedResults/Teeth/cPDist/cPDistMatrix.mat';

%% load results
load(TaxaPath);
load(cPDistPath);
disp(['loading maps from ' MapsPath '...']);
load(MapsPath);
disp('loaded.');

%% Parse parameters
TaxaCode = load(TaxaPath);
TaxaCode = TaxaCode.taxa_code;
GroupSize = length(TaxaCode);

% RootNode = 1; %% if one needs to fix a RootNode
[~,RootNode] = min(sum(cPDistMatrix.^2)); %% if one needs to fix a RootNode
disp(['RootNode = ' TaxaCode{RootNode}]);

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
        PropNodeFeatureInds = cPMapsMatrix{RootNode,j}(RootNodeFeatureInds);
        load([SamplePath TaxaCode{j} '.mat']);
        fprintf(fid, '%d %d %d\n', G.V(:,PropNodeFeatureInds));
        clear G;
    end
end

fclose(fid);




