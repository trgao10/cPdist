%% Preparation
clear vars;
close all;
path(pathdef);

%% Set Path
DataPath = '../DATA/PNAS/';
TaxaPath = [DataPath 'teeth_taxa_table.mat'];
LandmarkPath = [DataPath 'landmarks_teeth.mat'];

%% Set Parameters
GroupSize = 116;
NumFeatPts = 18;
output_filename = [pwd '/morphologika/PNAS_Observer_Landmarks_18.txt'];

%% load results
TaxaCode = load(TaxaPath);
TaxaCode = TaxaCode.taxa_code;
load(LandmarkPath);

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
fprintf(fid, '\n[rawpoints]\n');
for j=1:length(TaxaCode)
    fprintf(fid,'\n%s\n\n', ['''' TaxaCode{j}]);
    fprintf(fid, '%d %d %d\n', squeeze(PP(j,:,:))');
end

fclose(fid);




